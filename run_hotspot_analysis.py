"""
Script principal da análise de hotspots de desmatamento no Cerrado
(PRODES, 2020-2026), usando Getis-Ord Gi* ao nível de município.

Fluxo da análise:
  1. Lê os polígonos PRODES em blocos, sem carregar todo o arquivo na memória.
  2. Lê os municípios e o limite do Cerrado, recortando a malha municipal
     pelo bioma para incluir todos os municípios com área dentro do Cerrado.
  3. Intersecta a geometria de cada polígono com os municípios recortados e
     soma por município somente a área de interseção dentro do Cerrado.
  4. Constrói a matriz de pesos espaciais (contiguidade Queen) entre os
     municípios e calcula o Getis-Ord Gi* sobre a área desmatada acumulada.
  5. Classifica cada município como hot spot / cold spot / não significante,
     cataloga por estado (UF) e município, e salva os resultados em
     output/.

Uso:
    python run_hotspot_analysis.py
"""

import logging
import subprocess
import sys
import time

from src import aggregation, config, hotspots, io_utils

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger("run_hotspot_analysis")


def main():
    inicio = time.time()
    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Limite do Cerrado e polígonos municipais do Brasil
    cerrado_boundary = io_utils.load_cerrado_boundary()
    municipalities_br = io_utils.load_municipalities()
    municipalities_cerrado = aggregation.clip_municipalities_to_biome(
        municipalities_br,
        cerrado_boundary,
    ).reset_index(drop=True)

    # 2. Intersecta os polígonos PRODES com municípios recortados em blocos.
    #    Cada parcela é atribuída ao município e bioma onde está localizada.
    deforestation_stats = aggregation.aggregate_polygon_intersections(
        io_utils.iter_deforestation_polygon_chunks(),
        municipalities_cerrado,
    )
    result = aggregation.aggregate_by_municipality(
        deforestation_stats,
        municipalities_cerrado,
    )

    logger.info(
        "Malha de análise: %d municípios em %d estados",
        len(result),
        result["SIGLA_UF"].nunique(),
    )

    # 4. Pesos espaciais (contiguidade) + Getis-Ord Gi*
    w = hotspots.build_contiguity_weights(result)
    z_scores, p_values = hotspots.compute_gi_star(result, "area_desmatada_km2", w)
    result["gi_z_score"] = z_scores
    result["gi_p_value"] = p_values
    result["gi_bin"] = hotspots.classify_hotspots(z_scores, p_values)

    # Roda também sobre a taxa (% da área municipal dentro do Cerrado que
    # foi desmatada), útil para comparar municípios de tamanhos distintos.
    z_pct, p_pct = hotspots.compute_gi_star(result, "pct_area_desmatada", w)
    result["gi_z_score_pct"] = z_pct
    result["gi_p_value_pct"] = p_pct
    result["gi_bin_pct"] = hotspots.classify_hotspots(z_pct, p_pct)

    # 5. Cataloga por estado e município, e salva os resultados
    colunas_saida = [
        "CD_MUN",
        "NM_MUN",
        "SIGLA_UF",
        "NM_UF",
        "AREA_KM2",
        "AREA_CERRADO_KM2",
        "area_desmatada_km2",
        "pct_area_desmatada",
        "n_poligonos",
        "gi_z_score",
        "gi_p_value",
        "gi_bin",
        "gi_z_score_pct",
        "gi_p_value_pct",
        "gi_bin_pct",
        "geometry",
    ]
    result = result[colunas_saida].sort_values("gi_z_score", ascending=False)

    result.to_crs(config.INPUT_CRS).to_file(
        config.HOTSPOT_GEOJSON_PATH,
        driver="GeoJSON",
    )
    logger.info("Resultado salvo em %s", config.HOTSPOT_GEOJSON_PATH)

    result.drop(columns="geometry").to_csv(config.HOTSPOT_CSV_PATH, index=False)
    logger.info("Tabela salva em %s", config.HOTSPOT_CSV_PATH)

    resumo = (
        result.melt(
            id_vars=["SIGLA_UF"],
            value_vars=["gi_bin", "gi_bin_pct"],
            var_name="variavel",
            value_name="classificacao",
        )
        .groupby(["SIGLA_UF", "variavel", "classificacao"])
        .size()
        .unstack(fill_value=0)
    )
    resumo.to_csv(config.STATE_SUMMARY_PATH)
    logger.info("Resumo por estado salvo em %s", config.STATE_SUMMARY_PATH)

    visualization_ok = True
    try:
        # A plotagem roda em um subprocesso separado (ver src/make_maps.py):
        # se a etapa gráfica falhar de forma abrupta (crash nativo do
        # matplotlib, por exemplo), o processo principal não é derrubado
        # junto, e os resultados numéricos/geoespaciais salvos acima
        # permanecem intactos.
        subprocess.run(
            [sys.executable, "-m", "src.make_maps"],
            cwd=config.PROJECT_ROOT,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        visualization_ok = False
        logger.warning(
            "Falha ao gerar os mapas PNG (código de saída %s); os resultados "
            "tabulares/geoespaciais já foram salvos normalmente. Rode "
            "'python -m src.make_maps' novamente mais tarde para gerar "
            "apenas os mapas.",
            exc.returncode,
        )

    if visualization_ok:
        logger.info("Mapas PNG gerados com sucesso em %s", config.OUTPUT_DIR)

    logger.info(
        "Análise concluída em %.1f minutos", (time.time() - inicio) / 60
    )


if __name__ == "__main__":
    main()
