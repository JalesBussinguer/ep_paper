"""
Extrai a composição e as transições MapBiomas dentro de fragmentos PRODES
localizados em municípios hotspot.

Uso:
    python -m src.mapbiomas_transitions
"""

import csv
import gc
import gzip
import json
import logging
import math
import os
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
from exactextract import exact_extract

from . import config, io_utils
from .mapbiomas_legend import format_class_name

logger = logging.getLogger(__name__)

CHUNK_SIZE = 100_000
HOTSPOT_COLUMNS = [
    "CD_MUN",
    "NM_MUN",
    "SIGLA_UF",
    "NM_UF",
    "gi_bin",
    "gi_bin_pct",
    "geometry",
]
BASE_OUTPUT_COLUMNS = [
    "fid",
    "year",
    "CD_MUN",
    "NM_MUN",
    "SIGLA_UF",
    "NM_UF",
    "area_intersecao_km2",
]
YEARLY_CLASS_COLUMNS = [
    f"classe_{year}" for year in config.MAPBIOMAS_YEARS
]
OUTPUT_COLUMNS = [
    *BASE_OUTPUT_COLUMNS,
    *YEARLY_CLASS_COLUMNS,
    "classe_pre_evento",
    "classe_ano_deteccao",
    "classe_primeiro_ano_posterior",
    "status_pos_deteccao",
]
INTERSECTION_ATTRIBUTES = [
    "transition_id",
    "prodes_fid",
    "year",
    "CD_MUN",
    "NM_MUN",
    "SIGLA_UF",
    "NM_UF",
    "area_intersecao_km2",
]
TRANSITION_OUTPUT_COLUMNS = [
    "transition_id",
    "fid",
    "year_prodes",
    "CD_MUN",
    "NM_MUN",
    "SIGLA_UF",
    "NM_UF",
    "ano_mapbiomas_origem",
    "ano_mapbiomas_destino",
    "codigo_classe_origem",
    "classe_origem",
    "codigo_classe_destino",
    "classe_destino",
    "area_transicao_km2",
    "percentual_area_intersecao",
    "status_cobertura",
]
COMPOSITION_OUTPUT_COLUMNS = [
    "transition_id",
    "fid",
    "year_prodes",
    "CD_MUN",
    "NM_MUN",
    "SIGLA_UF",
    "NM_UF",
    "ano_mapbiomas",
    "codigo_classe",
    "classe",
    "area_classe_km2",
    "percentual_area_intersecao",
    "status_cobertura",
]


def class_transition_coverage(current_values, coverage_fractions, previous_values):
    """Resume a cobertura por transição entre as classes anuais dos pixels."""
    current = np.ma.asarray(current_values)
    previous = np.ma.asarray(previous_values)
    current_codes = np.asarray(current.data, dtype=np.int64)
    previous_codes = np.asarray(previous.data, dtype=np.int64)
    coverage = np.asarray(coverage_fractions, dtype=np.float64)
    valid = (
        ~np.ma.getmaskarray(current)
        & ~np.ma.getmaskarray(previous)
        & np.isfinite(coverage)
        & (coverage > 0)
    )
    if not valid.any():
        return "{}"

    if np.any(current_codes[valid] > 255) or np.any(previous_codes[valid] > 255):
        raise ValueError("Código raster excede o limite usado para codificar transições")

    transition_codes = previous_codes[valid] * 256 + current_codes[valid]
    unique_codes, inverse = np.unique(transition_codes, return_inverse=True)
    covered_pixels = np.bincount(inverse, weights=coverage[valid])
    return json.dumps(
        {
            str(int(code)): float(amount)
            for code, amount in zip(unique_codes, covered_pixels)
        },
        separators=(",", ":"),
    )


def _get_hotspot_municipalities(path=config.HOTSPOT_GEOJSON_PATH):
    """Seleciona a união de municípios hotspot nas duas variáveis Gi*."""
    hotspots = gpd.read_file(path, engine="pyogrio")
    absolute_hotspot = hotspots["gi_bin"].astype(str).str.startswith("Hot Spot")
    relative_hotspot = hotspots["gi_bin_pct"].astype(str).str.startswith("Hot Spot")
    selected = hotspots.loc[
        absolute_hotspot | relative_hotspot,
        HOTSPOT_COLUMNS,
    ].to_crs(config.METRIC_CRS)

    if selected.empty:
        raise ValueError(f"Nenhum município hotspot encontrado em {path}")

    selected = selected.reset_index(drop=True)
    logger.info(
        "%d municípios hotspot selecionados pela união das duas análises",
        len(selected),
    )
    return selected


def _validate_rasters():
    """Valida a disponibilidade e compatibilidade das sete grades anuais."""
    paths = {
        year: config.MAPBIOMAS_DIR
        / config.MAPBIOMAS_RASTER_PATTERN.format(year=year)
        for year in config.MAPBIOMAS_YEARS
    }
    missing = [path for path in paths.values() if not path.is_file()]
    if missing:
        missing_paths = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(
            "Faltam rasters anuais MapBiomas exigidos:\n" + missing_paths
        )

    reference = None
    for year, path in paths.items():
        with rasterio.open(path) as dataset:
            if dataset.count != 1:
                raise ValueError(f"Esperada uma banda em {path}; obtidas {dataset.count}")
            current = (
                dataset.crs,
                dataset.transform,
                dataset.width,
                dataset.height,
                dataset.dtypes[0],
            )
            if reference is None:
                reference = current
            elif current != reference:
                raise ValueError(
                    f"Grade raster incompatível no ano {year}: {path}"
                )
    return paths


def _prepare_intersection_features(polygons, municipalities):
    """Recorta os eventos por município hotspot e descarta contatos sem área."""
    valid_geometry = shapely.is_valid(polygons.geometry.array)
    invalid_count = int((~valid_geometry).sum())
    if invalid_count:
        geometry_column = polygons.geometry.name
        polygons.loc[~valid_geometry, geometry_column] = shapely.make_valid(
            polygons.loc[~valid_geometry].geometry.array
        )
        logger.warning("Reparadas %d geometrias PRODES inválidas", invalid_count)

    usable_geometry = (
        ~shapely.is_missing(polygons.geometry.array)
        & ~shapely.is_empty(polygons.geometry.array)
    )
    polygons = polygons.loc[usable_geometry].copy()
    if polygons.empty:
        return gpd.GeoDataFrame(
            columns=BASE_OUTPUT_COLUMNS + ["geometry"],
            geometry="geometry",
            crs=municipalities.crs,
        )

    joined = gpd.sjoin(
        polygons,
        municipalities[["CD_MUN", "NM_MUN", "SIGLA_UF", "NM_UF", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    if joined.empty:
        return gpd.GeoDataFrame(
            columns=BASE_OUTPUT_COLUMNS + ["geometry"],
            geometry="geometry",
            crs=municipalities.crs,
        )

    municipality_geometries = municipalities.geometry.reset_index(drop=True)
    right_positions = joined["index_right"].to_numpy(dtype=int)
    polygon_geometries = joined.geometry.array
    matched_municipalities = municipality_geometries.take(right_positions).array
    fully_covered = shapely.covered_by(
        polygon_geometries,
        matched_municipalities,
    )

    intersection_geometries = polygon_geometries.copy()
    partial_intersections = ~fully_covered
    if partial_intersections.any():
        intersection_geometries[partial_intersections] = shapely.intersection(
            polygon_geometries[partial_intersections],
            matched_municipalities[partial_intersections],
        )

    intersection_geometries = np.asarray(
        [_as_multipolygon(geometry) for geometry in intersection_geometries],
        dtype=object,
    )
    intersection_areas_km2 = shapely.area(intersection_geometries) / 1_000_000
    positive_area = intersection_areas_km2 > 0
    if not positive_area.any():
        return gpd.GeoDataFrame(
            columns=BASE_OUTPUT_COLUMNS + ["geometry"],
            geometry="geometry",
            crs=municipalities.crs,
        )

    data = pd.DataFrame(
        {
            "fid": joined["fid"].to_numpy()[positive_area],
            "year": joined["year"].to_numpy()[positive_area],
            "CD_MUN": joined["CD_MUN"].to_numpy()[positive_area],
            "NM_MUN": joined["NM_MUN"].to_numpy()[positive_area],
            "SIGLA_UF": joined["SIGLA_UF"].to_numpy()[positive_area],
            "NM_UF": joined["NM_UF"].to_numpy()[positive_area],
            "area_intersecao_km2": intersection_areas_km2[positive_area],
        }
    )
    geometries = intersection_geometries[positive_area]
    return gpd.GeoDataFrame(
        data,
        geometry=gpd.GeoSeries(geometries, crs=municipalities.crs),
    )


def _as_multipolygon(geometry):
    """Mantém apenas componentes poligonais e padroniza a geometria."""
    if geometry.geom_type == "Polygon":
        return shapely.MultiPolygon([geometry])
    if geometry.geom_type == "MultiPolygon":
        return geometry

    polygon_parts = []
    pending = [geometry]
    while pending:
        part = pending.pop()
        if part.geom_type == "Polygon":
            polygon_parts.append(part)
        elif part.geom_type in {"MultiPolygon", "GeometryCollection"}:
            pending.extend(part.geoms)

    if not polygon_parts:
        return shapely.MultiPolygon()

    merged = shapely.union_all(polygon_parts)
    if merged.geom_type == "Polygon":
        return shapely.MultiPolygon([merged])
    if merged.geom_type == "MultiPolygon":
        return merged
    return _as_multipolygon(merged)


def _extract_yearly_classes(features, raster_paths):
    """Calcula a classe dominante em cada polígono para cada ano raster."""
    output = features[BASE_OUTPUT_COLUMNS].copy()
    with rasterio.open(raster_paths[config.MAPBIOMAS_YEARS[0]]) as dataset:
        raster_features = features.to_crs(dataset.crs)

    for year, path in raster_paths.items():
        extracted = exact_extract(
            str(path),
            raster_features,
            ["majority"],
            include_cols=["fid"],
            output="pandas",
            strategy="feature-sequential",
        )
        class_values = extracted["majority"]
        output[f"classe_{year}"] = class_values.map(format_class_name)
    return output


def _codes_to_labels(class_codes):
    """Converte códigos raster em nomes oficiais, validando a legenda."""
    values = pd.Series(class_codes, dtype="float64").fillna(0).astype(int)
    labels = {
        int(code): format_class_name(int(code))
        for code in values.unique()
    }
    return values.map(labels)


def _add_event_relative_classes(transitions, class_codes):
    """Identifica referência pré-evento, ano da detecção e primeiro ano seguinte."""
    detection_years = transitions["year"].to_numpy(dtype=int)
    feature_ids = transitions["transition_id"].to_numpy(dtype=int)
    year_to_index = {
        year: index for index, year in enumerate(config.MAPBIOMAS_YEARS)
    }

    def select_by_event_offset(offset):
        selected_codes = np.zeros(len(transitions), dtype=np.uint8)
        target_years = detection_years + offset
        for year, year_index in year_to_index.items():
            mask = target_years == year
            selected_codes[mask] = class_codes[year_index, feature_ids[mask]]
        return _codes_to_labels(selected_codes)

    transitions["classe_pre_evento"] = select_by_event_offset(-1)
    transitions["classe_ano_deteccao"] = select_by_event_offset(0)
    post_years = detection_years + 1
    has_post_map = np.isin(post_years, config.MAPBIOMAS_YEARS)
    transitions["classe_primeiro_ano_posterior"] = None
    transitions.loc[has_post_map, "classe_primeiro_ano_posterior"] = (
        select_by_event_offset(1)[has_post_map].to_numpy()
    )
    transitions["status_pos_deteccao"] = np.where(
        has_post_map,
        "observado",
        "sem mapa posterior",
    )
    return transitions


def _update_summary_counts(transitions, counts):
    """Acumula contagens de classe em cada bloco, sem carregar a série inteira."""
    long_format = transitions.melt(
        id_vars=["CD_MUN", "NM_MUN", "SIGLA_UF", "NM_UF", "year"],
        value_vars=YEARLY_CLASS_COLUMNS,
        var_name="ano_mapbiomas",
        value_name="classe",
    )
    long_format["ano_mapbiomas"] = (
        long_format["ano_mapbiomas"].str.removeprefix("classe_").astype(int)
    )
    grouped = long_format.groupby(
        ["SIGLA_UF", "NM_UF", "CD_MUN", "NM_MUN", "year", "ano_mapbiomas", "classe"],
        dropna=False,
    ).size()
    counts.update(grouped.to_dict())


def _write_intersection_geopackage(polygon_chunks, municipalities, path):
    """Salva em GPKG as interseções dos eventos com municípios hotspot."""
    total_features = 0
    first_batch = True
    for chunk_number, polygons in enumerate(polygon_chunks, start=1):
        intersections = _prepare_intersection_features(polygons, municipalities)
        if intersections.empty:
            continue

        intersections = intersections.rename(columns={"fid": "prodes_fid"})
        intersections.insert(
            0,
            "transition_id",
            np.arange(total_features, total_features + len(intersections)),
        )
        intersections = intersections.to_crs("EPSG:4326")
        pyogrio.write_dataframe(
            intersections,
            path,
            layer="hotspot_fragments",
            driver="GPKG",
            geometry_type="MultiPolygon",
            append=not first_batch,
        )
        first_batch = False
        total_features += len(intersections)
        logger.info(
            "Bloco PRODES %d: %d interseções; total %d",
            chunk_number,
            len(intersections),
            total_features,
        )

    if total_features == 0:
        raise ValueError("Nenhum polígono PRODES intersectou municípios hotspot")
    layer_info = pyogrio.read_info(path, layer="hotspot_fragments")
    if layer_info["geometry_type"] != "MultiPolygon":
        raise RuntimeError(
            "O GeoPackage de interseções não foi padronizado como MultiPolygon: "
            f"{layer_info['geometry_type']}"
        )
    return total_features


def _mode_from_marginal(marginal):
    """Retorna a classe com maior cobertura, desempatando pelo menor código."""
    return max(marginal, key=lambda code: (marginal[code], -code))


def _progress_logger(year_from, year_to):
    last_reported_percent = -10

    def report_progress(fraction, _message):
        nonlocal last_reported_percent
        percent = int(fraction * 100)
        if percent >= last_reported_percent + 10 or percent == 100:
            logger.info(
                "Extração %d-%d: %d%%",
                year_from,
                year_to,
                percent,
            )
            last_reported_percent = percent

    return report_progress


def _extract_class_composition_and_transitions(
    raster_paths,
    geopackage_path,
    count,
    composition_csv_path,
    pixel_transitions_csv_path,
):
    """Extrai composição anual e transições de pixels entre anos consecutivos."""
    class_codes = np.zeros((len(config.MAPBIOMAS_YEARS), count), dtype=np.uint8)
    composition_summary = defaultdict(lambda: [0.0, 0])
    transition_summary = defaultdict(lambda: [0.0, 0])
    years = config.MAPBIOMAS_YEARS
    temporary_gpkg = Path(geopackage_path).with_name(
        f"{Path(geopackage_path).stem}_pixel_transitions.gpkg"
    )

    with (
        gzip.open(composition_csv_path, "wt", encoding="utf-8", newline="") as composition_file,
        gzip.open(
            pixel_transitions_csv_path,
            "wt",
            encoding="utf-8",
            newline="",
        ) as transitions_file,
    ):
        composition_writer = csv.writer(composition_file)
        transition_writer = csv.writer(transitions_file)
        composition_writer.writerow(COMPOSITION_OUTPUT_COLUMNS)
        transition_writer.writerow(TRANSITION_OUTPUT_COLUMNS)

        for pair_index, (year_from, year_to) in enumerate(zip(years[:-1], years[1:])):
            logger.info(
                "Calculando transições dos pixels MapBiomas %d → %d",
                year_from,
                year_to,
            )
            if temporary_gpkg.exists():
                temporary_gpkg.unlink()

            exact_extract(
                str(raster_paths[year_to]),
                str(geopackage_path),
                class_transition_coverage,
                weights=str(raster_paths[year_from]),
                include_cols=INTERSECTION_ATTRIBUTES,
                output="gdal",
                output_options={
                    "filename": str(temporary_gpkg),
                    "driver": "GPKG",
                    "layer_name": "pixel_transitions",
                },
                strategy="feature-sequential",
                progress=_progress_logger(year_from, year_to),
            )

            processed = 0
            missing_coverage_count = 0
            with pyogrio.open_arrow(
                temporary_gpkg,
                layer="pixel_transitions",
                columns=[*INTERSECTION_ATTRIBUTES, "class_transition_coverage"],
                read_geometry=False,
                batch_size=20_000,
                use_pyarrow=True,
            ) as (_, reader):
                for arrow_batch in reader:
                    batch = arrow_batch.to_pandas()
                    feature_ids = batch["transition_id"].to_numpy(dtype=int)
                    expected_ids = np.arange(
                        processed,
                        processed + len(batch),
                        dtype=int,
                    )
                    if not np.array_equal(feature_ids, expected_ids):
                        raise RuntimeError(
                            f"IDs de fragmento ausentes ou fora de ordem "
                            f"na transição {year_from}-{year_to}"
                        )

                    transition_rows = []
                    composition_rows = []
                    for row in batch.itertuples(index=False):
                        histogram = (
                            json.loads(row.class_transition_coverage)
                            if row.class_transition_coverage
                            else {}
                        )
                        has_valid_coverage = bool(histogram)
                        if not has_valid_coverage:
                            missing_coverage_count += 1
                        from_marginal = defaultdict(float)
                        to_marginal = defaultdict(float)
                        decoded_transitions = []
                        for encoded_code, covered_pixels in histogram.items():
                            transition_code = int(encoded_code)
                            from_code, to_code = divmod(transition_code, 256)
                            format_class_name(from_code)
                            format_class_name(to_code)
                            covered_pixels = float(covered_pixels)
                            from_marginal[from_code] += covered_pixels
                            to_marginal[to_code] += covered_pixels
                            decoded_transitions.append(
                                (from_code, to_code, covered_pixels)
                            )

                        if not has_valid_coverage:
                            from_marginal[0] = 0.0
                            to_marginal[0] = 0.0
                            decoded_transitions.append((0, 0, 0.0))
                            total_coverage = 0.0
                            coverage_status = "sem pixels válidos"
                        else:
                            total_coverage = math.fsum(histogram.values())
                            if (
                                not math.isfinite(total_coverage)
                                or total_coverage <= 0
                            ):
                                raise RuntimeError(
                                    f"Cobertura inválida no fragmento "
                                    f"{row.transition_id}"
                                )
                            coverage_status = "observado"

                        area_km2 = float(row.area_intersecao_km2)
                        for from_code, to_code, covered_pixels in decoded_transitions:
                            area_fraction = (
                                covered_pixels / total_coverage
                                if has_valid_coverage
                                else 0.0
                            )
                            transition_area_km2 = area_km2 * area_fraction
                            transition_rows.append(
                                [
                                    row.transition_id,
                                    row.prodes_fid,
                                    row.year,
                                    row.CD_MUN,
                                    row.NM_MUN,
                                    row.SIGLA_UF,
                                    row.NM_UF,
                                    year_from,
                                    year_to,
                                    from_code,
                                    format_class_name(from_code),
                                    to_code,
                                    format_class_name(to_code),
                                    transition_area_km2,
                                    area_fraction * 100,
                                    coverage_status,
                                ]
                            )
                            transition_key = (
                                row.SIGLA_UF,
                                row.NM_UF,
                                row.CD_MUN,
                                row.NM_MUN,
                                row.year,
                                year_from,
                                year_to,
                                from_code,
                                to_code,
                            )
                            transition_total = transition_summary[transition_key]
                            transition_total[0] += transition_area_km2
                            transition_total[1] += 1

                        class_codes[pair_index, row.transition_id] = _mode_from_marginal(
                            from_marginal
                        )
                        if pair_index == len(years) - 2:
                            class_codes[pair_index + 1, row.transition_id] = (
                                _mode_from_marginal(to_marginal)
                            )

                        composition_years = []
                        if pair_index == 0:
                            composition_years.append((year_from, from_marginal))
                        composition_years.append((year_to, to_marginal))
                        for map_year, marginal in composition_years:
                            for code, covered_pixels in sorted(marginal.items()):
                                area_fraction = (
                                    covered_pixels / total_coverage
                                    if has_valid_coverage
                                    else 0.0
                                )
                                class_area_km2 = area_km2 * area_fraction
                                composition_rows.append(
                                    [
                                        row.transition_id,
                                        row.prodes_fid,
                                        row.year,
                                        row.CD_MUN,
                                        row.NM_MUN,
                                        row.SIGLA_UF,
                                        row.NM_UF,
                                        map_year,
                                        code,
                                        format_class_name(code),
                                        class_area_km2,
                                        area_fraction * 100,
                                        coverage_status,
                                    ]
                                )
                                composition_key = (
                                    row.SIGLA_UF,
                                    row.NM_UF,
                                    row.CD_MUN,
                                    row.NM_MUN,
                                    row.year,
                                    map_year,
                                    code,
                                )
                                composition_total = composition_summary[composition_key]
                                composition_total[0] += class_area_km2
                                composition_total[1] += 1

                    transition_writer.writerows(transition_rows)
                    composition_writer.writerows(composition_rows)
                    processed += len(batch)

            if processed != count:
                raise RuntimeError(
                    f"Processados {processed} fragmentos na transição "
                    f"{year_from}-{year_to}; esperado {count}"
                )
            logger.info(
                "Transição %d-%d processada para %d fragmentos; "
                "%d sem pixels válidos",
                year_from,
                year_to,
                processed,
                missing_coverage_count,
            )
            for suffix in ("", "-wal", "-shm", "-journal"):
                sidecar = Path(str(temporary_gpkg) + suffix)
                if sidecar.exists():
                    sidecar.unlink()

    composition_columns = [
        "SIGLA_UF",
        "NM_UF",
        "CD_MUN",
        "NM_MUN",
        "year_prodes",
        "ano_mapbiomas",
        "codigo_classe",
        "classe",
        "area_classe_km2",
        "n_intersecoes",
    ]
    composition_summary_frame = pd.DataFrame(
        [
            (*key, format_class_name(key[-1]), values[0], values[1])
            for key, values in composition_summary.items()
        ],
        columns=composition_columns,
    )
    transition_columns = [
        "SIGLA_UF",
        "NM_UF",
        "CD_MUN",
        "NM_MUN",
        "year_prodes",
        "ano_mapbiomas_origem",
        "ano_mapbiomas_destino",
        "codigo_classe_origem",
        "classe_origem",
        "codigo_classe_destino",
        "classe_destino",
        "area_transicao_km2",
        "n_intersecoes",
    ]
    transition_summary_frame = pd.DataFrame(
        [
            (
                *key[:8],
                format_class_name(key[7]),
                key[8],
                format_class_name(key[8]),
                values[0],
                values[1],
            )
            for key, values in transition_summary.items()
        ],
        columns=transition_columns,
    )
    return class_codes, composition_summary_frame, transition_summary_frame


def _export_transition_tables(geopackage_path, count, class_codes, csv_path):
    """Exporta os registros por polígono e acumula a síntese em fluxo Arrow."""
    summary_counts = Counter()
    written = 0
    base_columns = [
        "transition_id",
        "prodes_fid",
        "year",
        "CD_MUN",
        "NM_MUN",
        "SIGLA_UF",
        "NM_UF",
        "area_intersecao_km2",
    ]
    with pyogrio.open_arrow(
        geopackage_path,
        layer="hotspot_fragments",
        columns=base_columns,
        read_geometry=False,
        batch_size=50_000,
        use_pyarrow=True,
    ) as (_, reader):
        for arrow_batch in reader:
            transitions = arrow_batch.to_pandas()
            feature_ids = transitions["transition_id"].to_numpy(dtype=int)
            for year_index, year in enumerate(config.MAPBIOMAS_YEARS):
                transitions[f"classe_{year}"] = _codes_to_labels(
                    class_codes[year_index, feature_ids]
                ).to_numpy()

            transitions = _add_event_relative_classes(transitions, class_codes)
            transitions = transitions.rename(columns={"prodes_fid": "fid"})
            transitions = transitions[OUTPUT_COLUMNS]
            transitions.to_csv(
                csv_path,
                mode="a",
                header=written == 0,
                index=False,
                encoding="utf-8",
            )
            _update_summary_counts(transitions, summary_counts)
            written += len(transitions)

    if written != count:
        raise RuntimeError(
            f"Exportados {written} registros, mas o GeoPackage contém {count}"
        )

    summary = pd.DataFrame(
        [(*key, value) for key, value in summary_counts.items()],
        columns=[
            "SIGLA_UF",
            "NM_UF",
            "CD_MUN",
            "NM_MUN",
            "year",
            "ano_mapbiomas",
            "classe",
            "n_intersecoes",
        ],
    )
    return summary


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    started = time.time()
    raster_paths = _validate_rasters()
    hotspot_municipalities = _get_hotspot_municipalities()

    output_paths = [
        config.MAPBIOMAS_TRANSITIONS_CSV_PATH,
        config.MAPBIOMAS_TRANSITION_SUMMARY_PATH,
        config.MAPBIOMAS_CLASS_COMPOSITION_CSV_PATH,
        config.MAPBIOMAS_PIXEL_TRANSITIONS_CSV_PATH,
        config.MAPBIOMAS_CLASS_COMPOSITION_SUMMARY_PATH,
        config.MAPBIOMAS_PIXEL_TRANSITIONS_SUMMARY_PATH,
    ]
    temporary_paths = {
        output_path: output_path.with_name(
            f"{output_path.name}.{os.getpid()}.tmp"
        )
        for output_path in output_paths
    }
    config.MAPBIOMAS_TRANSITIONS_CSV_PATH.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        prefix="mapbiomas_hotspot_fragments_",
        suffix=".gpkg",
        dir=config.MAPBIOMAS_DIR,
        delete=False,
    ) as temporary_file:
        geopackage_path = temporary_file.name
    os.unlink(geopackage_path)

    try:
        total_features = _write_intersection_geopackage(
            io_utils.iter_deforestation_polygon_chunks(
                columns=config.DEFORESTATION_COLUMNS,
                chunk_size=CHUNK_SIZE,
            ),
            hotspot_municipalities,
            geopackage_path,
        )
        logger.info(
            "%d interseções hotspot preparadas; iniciando extração de classes "
            "e transições anuais",
            total_features,
        )
        (
            class_codes,
            composition_summary,
            pixel_transition_summary,
        ) = _extract_class_composition_and_transitions(
            raster_paths,
            geopackage_path,
            total_features,
            temporary_paths[config.MAPBIOMAS_CLASS_COMPOSITION_CSV_PATH],
            temporary_paths[config.MAPBIOMAS_PIXEL_TRANSITIONS_CSV_PATH],
        )

        summary = _export_transition_tables(
            geopackage_path,
            total_features,
            class_codes,
            temporary_paths[config.MAPBIOMAS_TRANSITIONS_CSV_PATH],
        )
        summary.to_csv(
            temporary_paths[config.MAPBIOMAS_TRANSITION_SUMMARY_PATH],
            index=False,
            encoding="utf-8",
        )
        composition_summary.to_csv(
            temporary_paths[config.MAPBIOMAS_CLASS_COMPOSITION_SUMMARY_PATH],
            index=False,
            encoding="utf-8",
        )
        pixel_transition_summary.to_csv(
            temporary_paths[config.MAPBIOMAS_PIXEL_TRANSITIONS_SUMMARY_PATH],
            index=False,
            encoding="utf-8",
        )
        for output_path, temporary_path in temporary_paths.items():
            os.replace(temporary_path, output_path)
        logger.info(
            "Análise concluída: %d fragmentos em %.1f minutos",
            total_features,
            (time.time() - started) / 60,
        )
    finally:
        gc.collect()
        pixel_transitions_gpkg = (
            Path(geopackage_path).with_name(
                f"{Path(geopackage_path).stem}_pixel_transitions.gpkg"
            )
        )
        for temporary_path in (
            geopackage_path,
            geopackage_path + "-wal",
            geopackage_path + "-shm",
            geopackage_path + "-journal",
            str(pixel_transitions_gpkg),
            str(pixel_transitions_gpkg) + "-wal",
            str(pixel_transitions_gpkg) + "-shm",
            str(pixel_transitions_gpkg) + "-journal",
            *(str(path) for path in temporary_paths.values()),
        ):
            path = Path(temporary_path)
            if path.exists():
                for attempt in range(3):
                    try:
                        path.unlink()
                        break
                    except PermissionError:
                        if attempt == 2:
                            raise
                        gc.collect()
                        time.sleep(0.5)


if __name__ == "__main__":
    main()
