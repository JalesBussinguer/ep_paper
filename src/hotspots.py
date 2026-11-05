"""
Cálculo da estatística local de Getis-Ord Gi* (hot spot analysis) sobre os
municípios, seguindo a metodologia descrita em:
https://doc.esri.com/en/arcgis-pro/latest/tool-reference/spatial-statistics/
h-how-hot-spot-analysis-getis-ord-gi-spatial-stati.html

Resumo da metodologia aplicada:
  1. Cada município é uma unidade de análise com um valor observado (área
     desmatada acumulada no período, em km2).
  2. Construímos uma matriz de pesos espaciais de contiguidade (Queen: dois
     municípios são vizinhos se compartilham uma borda ou um vértice), que
     é a conceituação de vizinhança recomendada pela Esri para dados de
     polígonos (feature type "Polygon").
  3. O Gi* de cada município compara a soma do valor do próprio município
     com a soma dos seus vizinhos (star=True) contra o que seria esperado
     sob aleatoriedade espacial completa (CSR). Um Z score alto e positivo
     indica um cluster de valores altos (hot spot); um Z score baixo e
     negativo indica um cluster de valores baixos (cold spot).
  4. A significância estatística é avaliada via pseudo p-valor por
     permutação (Monte Carlo), mais robusto do que a aproximação normal
     quando os dados não são normalmente distribuídos (o que é o caso de
     área desmatada, que tem distribuição assimétrica com muitos zeros).
"""

import logging

import numpy as np
from esda.getisord import G_Local
from libpysal.weights import Queen

from . import config

logger = logging.getLogger(__name__)


def build_contiguity_weights(municipalities):
    """
    Constrói a matriz de pesos espaciais de contiguidade Queen (borda ou
    vértice compartilhado) entre os municípios, e trata ilhas (municípios
    sem nenhum vizinho contíguo, comuns em ilhas/enclaves) atribuindo o
    vizinho mais próximo, para que nenhum município fique fora da análise.
    """
    logger.info("Construindo matriz de pesos espaciais (contiguidade Queen)")
    w = Queen.from_dataframe(municipalities, use_index=False, silence_warnings=True)

    ilhas = w.islands
    if ilhas:
        logger.info(
            "%d município(s) sem vizinho de contiguidade; conectando ao "
            "vizinho mais próximo (centroide a centroide)",
            len(ilhas),
        )
        from libpysal.weights import KNN

        knn = KNN.from_dataframe(municipalities, k=1, use_index=False)
        neighbors = dict(w.neighbors)
        weights = dict(w.weights)
        for ilha in ilhas:
            vizinho = knn.neighbors[ilha][0]
            neighbors[ilha] = [vizinho]
            weights[ilha] = [1.0]
        from libpysal.weights import W

        w = W(neighbors, weights, silence_warnings=True)

    w.transform = "B"  # pesos binários, conforme esperado pelo G_Local(star=True)
    return w


def compute_gi_star(municipalities, value_col, w):
    """
    Executa o Getis-Ord Gi* (com a própria observação incluída na soma
    local, star=True) sobre a coluna `value_col` dos municípios, usando os
    pesos espaciais `w`, e retorna o Z score e o pseudo p-valor (por
    permutação).
    """
    y = municipalities[value_col].to_numpy(dtype=float)

    logger.info(
        "Calculando Gi* para '%s' (%d municípios, %d permutações)",
        value_col,
        len(y),
        config.GI_STAR_PERMUTATIONS,
    )
    gi = G_Local(
        y,
        w,
        transform="B",
        star=True,
        permutations=config.GI_STAR_PERMUTATIONS,
        seed=config.GI_STAR_SEED,
        n_jobs=1,  # evita multiprocessing (mais estavel em ambientes Windows)
    )

    return gi.Zs, gi.p_sim


def classify_hotspots(z_scores, p_values, z_thresholds=config.CONFIDENCE_Z_THRESHOLDS):
    """
    Classifica cada município em hot spot / cold spot / sem significância,
    reproduzindo a lógica do campo "Gi_Bin" do ArcGIS: um município só é
    considerado estatisticamente significante se o Z score ultrapassar o
    limiar bilateral do nível de confiança E o pseudo p-valor (permutação)
    também for menor ou igual ao alfa correspondente.

    Retorna um array de strings com uma das categorias:
    "Hot Spot 99%", "Hot Spot 95%", "Hot Spot 90%",
    "Cold Spot 99%", "Cold Spot 95%", "Cold Spot 90%", "Não significante".
    """
    alphas = {0.90: 0.10, 0.95: 0.05, 0.99: 0.01}
    # Ordena do nível de confiança mais alto para o mais baixo, para que a
    # primeira classificação que "encaixar" seja a mais restritiva possível.
    niveis_ordenados = sorted(z_thresholds.items(), key=lambda kv: -kv[1])

    labels = np.full(len(z_scores), "Não significante", dtype=object)
    resolvido = np.zeros(len(z_scores), dtype=bool)

    for nivel, z_critico in niveis_ordenados:
        alfa = alphas[nivel]
        percentual = int(round(nivel * 100))

        hot = (~resolvido) & (z_scores >= z_critico) & (p_values <= alfa)
        cold = (~resolvido) & (z_scores <= -z_critico) & (p_values <= alfa)

        labels[hot] = f"Hot Spot {percentual}%"
        labels[cold] = f"Cold Spot {percentual}%"
        resolvido |= hot | cold

    return labels
