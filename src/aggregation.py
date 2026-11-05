"""
Agregação por município da área de interseção dos polígonos PRODES com o
limite municipal recortado pelo Cerrado.

Polígonos que cruzam limites municipais ou a borda do Cerrado contribuem
somente com a área geométrica efetivamente contida em cada unidade.
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

logger = logging.getLogger(__name__)

MUNICIPALITY_ID_COL = "CD_MUN"


def aggregate_by_municipality(deforestation_stats, municipalities):
    """
    Junta as estatísticas de interseção PRODES ao recorte municipal e
    retorna um registro por município, incluindo unidades com valor zero.

    Municípios sem nenhum polígono de desmatamento no período permanecem no
    resultado com area_desmatada_km2 = 0 e n_poligonos = 0, o que é
    importante para a análise de Getis-Ord Gi*: a malha espacial usada no
    cálculo dos vizinhos (contiguidade) precisa ficar completa, sem buracos,
    mesmo que o valor observado em alguns municípios seja zero.
    """
    result = municipalities.merge(
        deforestation_stats,
        on=MUNICIPALITY_ID_COL,
        how="left",
    )
    result["area_desmatada_km2"] = result["area_desmatada_km2"].fillna(0.0)
    result["n_poligonos"] = result["n_poligonos"].fillna(0).astype(int)

    # Usa como denominador somente a porção municipal dentro do Cerrado,
    # calculada após o recorte geométrico pelo limite oficial do bioma.
    result["pct_area_desmatada"] = (
        100 * result["area_desmatada_km2"] / result["AREA_CERRADO_KM2"]
    )

    return result


def aggregate_polygon_intersections(polygon_chunks, municipalities):
    """
    Intersecta cada polígono PRODES com os municípios recortados pelo Cerrado
    e agrega a área efetivamente sobreposta por município.

    Polígonos que cruzam limites municipais são repartidos entre os
    municípios conforme a área de interseção; porções fora do bioma não
    contribuem para a área desmatada.
    """
    municipality_geometries = municipalities.geometry.reset_index(drop=True)
    chunk_summaries = []
    total_candidates = 0
    total_intersections = 0

    for chunk_number, polygons in enumerate(polygon_chunks, start=1):
        polygons = polygons.copy()
        valid_geometry = shapely.is_valid(polygons.geometry.array)
        invalid_count = int((~valid_geometry).sum())
        if invalid_count:
            geometry_column = polygons.geometry.name
            polygons.loc[~valid_geometry, geometry_column] = shapely.make_valid(
                polygons.loc[~valid_geometry].geometry.array
            )
            logger.warning(
                "Bloco %d: corrigidas %d geometrias PRODES inválidas",
                chunk_number,
                invalid_count,
            )

        usable_geometry = (
            ~shapely.is_missing(polygons.geometry.array)
            & ~shapely.is_empty(polygons.geometry.array)
        )
        discarded_count = int((~usable_geometry).sum())
        if discarded_count:
            logger.warning(
                "Bloco %d: descartadas %d geometrias ausentes ou vazias",
                chunk_number,
                discarded_count,
            )
            polygons = polygons.loc[usable_geometry].copy()
        if polygons.empty:
            continue

        joined = gpd.sjoin(
            polygons,
            municipalities[["CD_MUN", "geometry"]],
            how="inner",
            predicate="intersects",
        )
        if joined.empty:
            logger.info("Bloco %d: nenhuma interseção com os municípios", chunk_number)
            continue

        right_positions = joined["index_right"].to_numpy(dtype=int)
        polygon_geometries = joined.geometry.array
        matched_municipalities = municipality_geometries.take(right_positions).array
        fully_covered = shapely.covered_by(
            polygon_geometries,
            matched_municipalities,
        )
        intersection_areas_km2 = shapely.area(polygon_geometries) / 1_000_000
        partial_intersections = ~fully_covered
        if partial_intersections.any():
            clipped_geometries = shapely.intersection(
                polygon_geometries[partial_intersections],
                matched_municipalities[partial_intersections],
            )
            intersection_areas_km2[partial_intersections] = (
                shapely.area(clipped_geometries) / 1_000_000
            )
        positive_area = intersection_areas_km2 > 0
        if not positive_area.any():
            logger.info("Bloco %d: apenas contatos sem área", chunk_number)
            continue

        intersections = pd.DataFrame(
            {
                "CD_MUN": joined["CD_MUN"].to_numpy()[positive_area],
                "fid": joined["fid"].to_numpy()[positive_area],
                "area_intersecao_km2": intersection_areas_km2[positive_area],
            }
        )
        chunk_summary = (
            intersections.groupby("CD_MUN")
            .agg(
                area_desmatada_km2=("area_intersecao_km2", "sum"),
                n_poligonos=("fid", "nunique"),
            )
            .reset_index()
        )
        chunk_summaries.append(chunk_summary)
        total_candidates += len(joined)
        total_intersections += int(positive_area.sum())
        logger.info(
            "Bloco %d: %d candidatos espaciais, %d integralmente contidos, "
            "%d com área positiva",
            chunk_number,
            len(joined),
            int(fully_covered.sum()),
            int(positive_area.sum()),
        )

    if not chunk_summaries:
        raise ValueError("Nenhum polígono PRODES intersectou os municípios do Cerrado")

    all_summaries = pd.concat(chunk_summaries, ignore_index=True)
    stats = (
        all_summaries.groupby("CD_MUN")
        .agg(
            area_desmatada_km2=("area_desmatada_km2", "sum"),
            n_poligonos=("n_poligonos", "sum"),
        )
        .reset_index()
    )
    logger.info(
        "Interseção concluída: %d candidatos e %d recortes com área positiva",
        total_candidates,
        total_intersections,
    )
    return stats


def iter_polygon_intersections(polygon_chunks, municipalities):
    """
    Gera, para cada bloco de polígonos PRODES, as interseções (polígono x
    município recortado pelo Cerrado) com área positiva e a respectiva
    geometria. Segue a mesma lógica de `aggregate_polygon_intersections`.
    """
    municipality_geometries = municipalities.geometry.reset_index(drop=True)
    for polygons in polygon_chunks:
        polygons = polygons.copy()
        geometry_column = polygons.geometry.name
        valid_geometry = shapely.is_valid(polygons.geometry.array)
        if (~valid_geometry).any():
            polygons.loc[~valid_geometry, geometry_column] = shapely.make_valid(
                polygons.loc[~valid_geometry].geometry.array
            )
        usable = (
            ~shapely.is_missing(polygons.geometry.array)
            & ~shapely.is_empty(polygons.geometry.array)
        )
        polygons = polygons.loc[usable]
        if polygons.empty:
            continue

        joined = gpd.sjoin(
            polygons,
            municipalities[["CD_MUN", "geometry"]],
            how="inner",
            predicate="intersects",
        )
        if joined.empty:
            continue

        matched = municipality_geometries.take(
            joined["index_right"].to_numpy(dtype=int)
        ).array
        geometries = np.array(joined.geometry.array, dtype=object)
        covered = shapely.covered_by(geometries, matched)
        if (~covered).any():
            geometries[~covered] = shapely.intersection(
                geometries[~covered], matched[~covered]
            )
        areas_km2 = shapely.area(geometries) / 1_000_000
        positive = areas_km2 > 0
        if not positive.any():
            continue

        yield gpd.GeoDataFrame(
            {
                "CD_MUN": joined["CD_MUN"].to_numpy()[positive],
                "prodes_fid": joined["fid"].to_numpy()[positive],
                "year": joined["year"].to_numpy()[positive].astype(int),
                "area_km2": areas_km2[positive],
            },
            geometry=geometries[positive],
            crs=polygons.crs,
        )


def clip_municipalities_to_biome(municipalities, biome_boundary):
    """
    Recorta os limites municipais pelo Cerrado e mantém todos os municípios
    com área positiva dentro do bioma, inclusive os sem desmatamento.
    """
    clipped = municipalities.copy()
    clipped.geometry = municipalities.geometry.intersection(biome_boundary)
    invalid_geometry = ~shapely.is_valid(clipped.geometry.array)
    if invalid_geometry.any():
        logger.warning(
            "Corrigidas %d geometrias municipais inválidas após o recorte",
            int(invalid_geometry.sum()),
        )
        clipped.loc[invalid_geometry, "geometry"] = shapely.make_valid(
            clipped.loc[invalid_geometry].geometry.array
        )
    clipped["AREA_CERRADO_KM2"] = clipped.geometry.area / 1_000_000
    clipped = clipped[
        (~clipped.geometry.is_empty)
        & (clipped["AREA_CERRADO_KM2"] > 0)
    ].copy()

    if clipped.empty:
        raise ValueError("O recorte dos municípios não intersectou o limite do Cerrado")

    area_bioma_km2 = biome_boundary.area / 1_000_000
    cobertura_pct = 100 * clipped["AREA_CERRADO_KM2"].sum() / area_bioma_km2
    logger.info(
        "Recorte do Cerrado: %d municípios em %d estados; área municipal "
        "recortada %.1f km² (%.4f%% da área do bioma)",
        len(clipped),
        clipped["SIGLA_UF"].nunique(),
        clipped["AREA_CERRADO_KM2"].sum(),
        cobertura_pct,
    )
    return clipped
