"""
Funções de carregamento dos dados de entrada (desmatamento e municípios).

O arquivo de desmatamento é muito grande (~3 GB, ~1,95 milhão de polígonos),
então aqui priorizamos leitura seletiva de colunas via pyogrio/geopandas e
evitamos qualquer cópia desnecessária dos dados.
"""

import logging
import time

import geopandas as gpd
import pyogrio
import shapely

from . import config

logger = logging.getLogger(__name__)


def iter_deforestation_polygon_chunks(
    path=config.DEFORESTATION_PATH,
    columns=config.DEFORESTATION_COLUMNS,
    input_crs=config.INPUT_CRS,
    metric_crs=config.METRIC_CRS,
    chunk_size=100_000,
):
    """
    Lê os polígonos PRODES em blocos e os reprojeta para o CRS métrico.

    A leitura em blocos limita o uso de memória durante a interseção exata
    dos polígonos com municípios recortados pelo limite do Cerrado.
    """
    logger.info("Lendo polígonos de desmatamento de %s", path)
    t0 = time.time()

    total = 0
    with pyogrio.open_arrow(
        path,
        columns=columns,
        batch_size=chunk_size,
        use_pyarrow=True,
    ) as (metadata, reader):
        source_crs = metadata.get("crs") or input_crs
        for chunk in reader:
            polygons = gpd.GeoDataFrame.from_arrow(chunk)
            if polygons.crs is None:
                polygons = polygons.set_crs(source_crs)
            polygons = polygons.to_crs(metric_crs)
            yield polygons

            total += len(polygons)
            logger.info(
                "Lidos %d polígonos PRODES (%.1f s)",
                total,
                time.time() - t0,
            )

    if total == 0:
        raise ValueError(f"Nenhum polígono foi lido do arquivo {path}")
    logger.info(
        "Leitura concluída: %d polígonos em %.1f s", total, time.time() - t0
    )


def load_municipalities(
    path=config.MUNICIPALITIES_PATH,
    columns=config.MUNICIPALITY_COLUMNS,
    input_crs=config.INPUT_CRS,
    metric_crs=config.METRIC_CRS,
):
    """
    Carrega os polígonos de município (unidade territorial de análise) com
    as colunas de identificação (município/estado) e a área oficial do
    município (AREA_KM2), reprojetados para o CRS métrico.
    """
    logger.info("Lendo polígonos de municípios de %s", path)
    t0 = time.time()

    gdf = gpd.read_file(path, columns=columns, engine="pyogrio")
    logger.info(
        "Leitura concluída: %d municípios em %.1f s", len(gdf), time.time() - t0
    )

    if gdf.crs is None:
        gdf = gdf.set_crs(input_crs)
    gdf = gdf.to_crs(metric_crs)

    # AREA_KM2 vem como string no GeoJSON original; convertemos para numérico.
    gdf["AREA_KM2"] = gdf["AREA_KM2"].astype(float)
    invalid_geometry = ~shapely.is_valid(gdf.geometry.array)
    if invalid_geometry.any():
        logger.warning(
            "Corrigindo %d geometrias municipais inválidas antes do recorte",
            int(invalid_geometry.sum()),
        )
        gdf.loc[invalid_geometry, "geometry"] = shapely.make_valid(
            gdf.loc[invalid_geometry].geometry.array
        )

    return gdf


def load_cerrado_boundary(
    path=config.CERRADO_BOUNDARY_PATH,
    columns=config.CERRADO_BOUNDARY_COLUMNS,
    input_crs=config.INPUT_CRS,
    metric_crs=config.METRIC_CRS,
):
    """
    Carrega e dissolve o limite do Cerrado, reprojetando-o para o CRS métrico.
    """
    logger.info("Lendo limite do Cerrado de %s", path)
    cerrado = gpd.read_file(path, columns=columns, engine="pyogrio")

    if cerrado.crs is None:
        cerrado = cerrado.set_crs(input_crs)
    cerrado = cerrado.to_crs(metric_crs)

    cerrado_features = cerrado[
        cerrado["bioma"].fillna("").str.casefold() == "cerrado"
    ]
    if cerrado_features.empty:
        raise ValueError(f"Nenhuma geometria do bioma Cerrado encontrada em {path}")

    boundary = cerrado_features.geometry.union_all()
    if boundary.is_empty or boundary.area <= 0:
        raise ValueError(f"Limite do Cerrado inválido ou vazio em {path}")
    if not shapely.is_valid(boundary):
        logger.warning("Corrigindo geometria inválida no limite do Cerrado")
        boundary = shapely.make_valid(boundary)

    logger.info(
        "Limite do Cerrado carregado: %.1f km²",
        boundary.area / 1_000_000,
    )
    return boundary
