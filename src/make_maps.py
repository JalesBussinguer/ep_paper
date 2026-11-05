"""
Script standalone para gerar os mapas PNG a partir do GeoJSON de resultados
já salvo (output/municipios_hotspots_cerrado.geojson).

É mantido separado do run_hotspot_analysis.py e chamado como um subprocesso
porque a geração de gráficos (matplotlib) depende de bibliotecas nativas que
podem falhar de forma abrupta (crash do processo, não uma exceção Python
capturável) em ambientes instáveis; isolando essa etapa em um processo
separado, uma eventual falha na plotagem não afeta os resultados
numéricos/geoespaciais já salvos pela análise principal.

Uso:
    python -m src.make_maps
"""

import logging

import geopandas as gpd

from . import config, visualization

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def main():
    result = gpd.read_file(
        config.HOTSPOT_GEOJSON_PATH,
        engine="pyogrio",
    ).to_crs(config.METRIC_CRS)

    # Usa o rasterizador baseado em Pillow (ver src/visualization.py), que
    # não depende do matplotlib. Se preferir o mapa via matplotlib/geopandas
    # (visualization.plot_hotspot_map) em um ambiente onde ele funcione
    # normalmente, basta trocar a chamada abaixo.
    visualization.render_hotspot_map_pil(
        result,
        "gi_bin",
        "Hotspots de desmatamento no Cerrado (2020-2026)\n"
        "Getis-Ord Gi* sobre área desmatada acumulada por município",
        config.AREA_MAP_PATH,
    )
    visualization.render_hotspot_map_pil(
        result,
        "gi_bin_pct",
        "Hotspots de desmatamento no Cerrado (2020-2026)\n"
        "Getis-Ord Gi* sobre % da área municipal dentro do Cerrado desmatada",
        config.PERCENT_MAP_PATH,
    )


if __name__ == "__main__":
    main()
