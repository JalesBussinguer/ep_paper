"""Analisa mudanças de classe antes, durante e depois do evento PRODES."""

import logging
import os
from collections import defaultdict

import numpy as np
import pandas as pd

from . import config
from .mapbiomas_legend import format_class_name

logger = logging.getLogger(__name__)

CHUNK_SIZE = 50_000
INPUT_COLUMNS = [
    "year",
    "CD_MUN",
    "NM_MUN",
    "SIGLA_UF",
    "NM_UF",
    "area_intersecao_km2",
    "classe_pre_evento",
    "classe_ano_deteccao",
    "classe_primeiro_ano_posterior",
]
CLASS_GROUP_NAMES = {
    "vegetacao_natural": "Vegetação natural",
    "uso_antrópico": "Uso antrópico",
    "superficie_nao_vegetada": "Superfície não vegetada",
    "agua": "Água",
    "sem_dados": "Sem dados",
}
CLASS_GROUP_BY_CODE = {
    0: "sem_dados",
    **{
        code: "vegetacao_natural"
        for code in (1, 3, 4, 5, 6, 10, 11, 12, 13, 23, 29, 32, 49, 50)
    },
    **{
        code: "uso_antrópico"
        for code in (9, 14, 15, 18, 19, 20, 21, 24, 30, 31, 36, 48, 75)
    },
    22: "superficie_nao_vegetada",
    25: "superficie_nao_vegetada",
    26: "agua",
    27: "sem_dados",
    33: "agua",
    # Código presente nos rasters (~0,002% dos pixels), mas ausente da legenda Col4
    91: "sem_dados",
}
TRANSITION_PAIRS = (
    ("pré-evento → detecção", "classe_pre_evento", "classe_ano_deteccao"),
    (
        "detecção → pós-evento",
        "classe_ano_deteccao",
        "classe_primeiro_ano_posterior",
    ),
    (
        "pré-evento → pós-evento",
        "classe_pre_evento",
        "classe_primeiro_ano_posterior",
    ),
)
SUMMARY_COLUMNS = [
    "intervalo",
    "ano_prodes",
    "n_intersecoes_validas",
    "n_classes_dominantes_alteradas",
    "percentual_intersecoes_classes_alteradas",
    "n_grupos_alterados",
    "percentual_intersecoes_grupos_alterados",
    "area_intersecoes_validas_km2",
    "area_intersecoes_classes_alteradas_km2",
    "percentual_area_intersecoes_classes_alteradas",
    "area_intersecoes_grupos_alterados_km2",
    "percentual_area_intersecoes_grupos_alterados",
]
CLASS_TRANSITION_COLUMNS = [
    "intervalo",
    "ano_prodes",
    "codigo_classe_origem",
    "classe_origem",
    "grupo_origem",
    "codigo_classe_destino",
    "classe_destino",
    "grupo_destino",
    "n_intersecoes",
    "percentual_intersecoes",
    "area_intersecoes_km2",
    "percentual_area_intersecoes",
]
GROUP_TRANSITION_COLUMNS = [
    "intervalo",
    "ano_prodes",
    "grupo_origem",
    "grupo_destino",
    "n_intersecoes",
    "percentual_intersecoes",
    "area_intersecoes_km2",
    "percentual_area_intersecoes",
]
MUNICIPAL_TRANSITION_COLUMNS = [
    "SIGLA_UF",
    "NM_UF",
    "CD_MUN",
    "NM_MUN",
    "ano_prodes",
    "intervalo",
    "grupo_origem",
    "grupo_destino",
    "n_intersecoes",
    "percentual_intersecoes",
    "area_intersecoes_km2",
    "percentual_area_intersecoes",
]
TRAJECTORY_COLUMNS = [
    "ano_prodes",
    "grupo_pre_evento",
    "grupo_ano_deteccao",
    "grupo_pos_evento",
    "trajetoria",
    "n_intersecoes",
    "percentual_intersecoes",
    "area_intersecoes_km2",
    "percentual_area_intersecoes",
]


def _parse_class_codes(classes):
    """Converte rótulos da legenda em códigos, validando valores inesperados."""
    labels = classes.astype("string").str.strip().replace("", pd.NA)
    codes = pd.to_numeric(labels.str.extract(r"^(\d+)\s*-")[0], errors="coerce")
    codes = codes.mask(labels.eq("Sem dados").fillna(False), 0)
    invalid = labels.notna() & codes.isna()
    if invalid.any():
        unexpected = labels.loc[invalid].drop_duplicates().tolist()
        raise ValueError(f"Rótulos de classe não reconhecidos: {unexpected}")
    return codes.astype("Int16")


def _class_groups(codes):
    """Agrupa códigos MapBiomas em categorias para leitura ecológica."""
    groups = codes.map(CLASS_GROUP_BY_CODE)
    unknown = codes.notna() & groups.isna()
    if unknown.any():
        codes_without_group = codes.loc[unknown].drop_duplicates().tolist()
        raise ValueError(
            f"Códigos sem grupo analítico definido: {codes_without_group}"
        )
    return groups


def _add_accumulated(total, grouped):
    """Acumula contagens e áreas agrupadas entre blocos de leitura."""
    for key, values in grouped.iterrows():
        current = total[key]
        current[0] += int(values["n_intersecoes"])
        current[1] += float(values["area_intersecoes_km2"])


def _aggregate_event_transitions(input_path):
    class_transitions = defaultdict(lambda: [0, 0.0])
    group_transitions = defaultdict(lambda: [0, 0.0])
    municipal_group_transitions = defaultdict(lambda: [0, 0.0])
    summary_stats = defaultdict(
        lambda: {
            "valid": 0,
            "class_changed": 0,
            "group_changed": 0,
            "area_valid": 0.0,
            "area_class_changed": 0.0,
            "area_group_changed": 0.0,
        }
    )
    trajectories = defaultdict(lambda: [0, 0.0])
    processed = 0
    missing_post = 0
    detection_years = set()

    for chunk in pd.read_csv(
        input_path,
        usecols=INPUT_COLUMNS,
        encoding="utf-8",
        dtype={"CD_MUN": "string"},
        chunksize=CHUNK_SIZE,
    ):
        processed += len(chunk)
        year = pd.to_numeric(chunk["year"], errors="raise").astype(int)
        detection_years.update(year.unique())
        area = pd.to_numeric(chunk["area_intersecao_km2"], errors="raise")
        if (~np.isfinite(area)).any() or area.lt(0).any():
            raise ValueError("A área das interseções contém valores inválidos")

        codes = {
            column: _parse_class_codes(chunk[column])
            for column in (
                "classe_pre_evento",
                "classe_ano_deteccao",
                "classe_primeiro_ano_posterior",
            )
        }
        groups = {
            column: _class_groups(class_codes)
            for column, class_codes in codes.items()
        }
        missing_post += int(codes["classe_primeiro_ano_posterior"].isna().sum())

        for interval, source_column, destination_column in TRANSITION_PAIRS:
            valid = (
                codes[source_column].notna()
                & codes[destination_column].notna()
            )
            if not valid.any():
                continue

            pair = pd.DataFrame(
                {
                    "ano_prodes": year.loc[valid].to_numpy(),
                    "codigo_classe_origem": codes[source_column]
                    .loc[valid]
                    .astype(int)
                    .to_numpy(),
                    "codigo_classe_destino": codes[destination_column]
                    .loc[valid]
                    .astype(int)
                    .to_numpy(),
                    "grupo_origem": groups[source_column].loc[valid].to_numpy(),
                    "grupo_destino": groups[destination_column]
                    .loc[valid]
                    .to_numpy(),
                    "area_intersecao_km2": area.loc[valid].to_numpy(),
                    "SIGLA_UF": chunk.loc[valid, "SIGLA_UF"].to_numpy(),
                    "NM_UF": chunk.loc[valid, "NM_UF"].to_numpy(),
                    "CD_MUN": chunk.loc[valid, "CD_MUN"].to_numpy(),
                    "NM_MUN": chunk.loc[valid, "NM_MUN"].to_numpy(),
                }
            )
            pair["intervalo"] = interval
            pair["classe_origem"] = pair["codigo_classe_origem"].map(
                format_class_name
            )
            pair["classe_destino"] = pair["codigo_classe_destino"].map(
                format_class_name
            )
            pair["classe_changed"] = (
                pair["codigo_classe_origem"] != pair["codigo_classe_destino"]
            )
            pair["group_changed"] = pair["grupo_origem"] != pair["grupo_destino"]
            pair["area_class_changed"] = (
                pair["area_intersecao_km2"] * pair["classe_changed"]
            )
            pair["area_group_changed"] = (
                pair["area_intersecao_km2"] * pair["group_changed"]
            )

            summary = pair.groupby("ano_prodes").agg(
                n_intersecoes_validas=("ano_prodes", "size"),
                n_classes_dominantes_alteradas=("classe_changed", "sum"),
                n_grupos_alterados=("group_changed", "sum"),
                area_intersecoes_validas_km2=("area_intersecao_km2", "sum"),
                area_intersecoes_classes_alteradas_km2=(
                    "area_class_changed",
                    "sum",
                ),
                area_intersecoes_grupos_alterados_km2=("area_group_changed", "sum"),
            )
            for detection_year, values in summary.iterrows():
                for cohort_year in ("Todos", int(detection_year)):
                    stats = summary_stats[(interval, cohort_year)]
                    stats["valid"] += int(values["n_intersecoes_validas"])
                    stats["class_changed"] += int(
                        values["n_classes_dominantes_alteradas"]
                    )
                    stats["group_changed"] += int(values["n_grupos_alterados"])
                    stats["area_valid"] += float(
                        values["area_intersecoes_validas_km2"]
                    )
                    stats["area_class_changed"] += float(
                        values["area_intersecoes_classes_alteradas_km2"]
                    )
                    stats["area_group_changed"] += float(
                        values["area_intersecoes_grupos_alterados_km2"]
                    )

            exact_keys = [
                "intervalo",
                "codigo_classe_origem",
                "classe_origem",
                "grupo_origem",
                "codigo_classe_destino",
                "classe_destino",
                "grupo_destino",
            ]
            exact_grouped = pair.groupby(exact_keys, dropna=False).agg(
                n_intersecoes=("ano_prodes", "size"),
                area_intersecoes_km2=("area_intersecao_km2", "sum"),
            )
            _add_accumulated(class_transitions, exact_grouped)

            cohort_exact = pair.groupby(
                ["intervalo", "ano_prodes", *exact_keys[1:]],
                dropna=False,
            ).agg(
                n_intersecoes=("ano_prodes", "size"),
                area_intersecoes_km2=("area_intersecao_km2", "sum"),
            )
            _add_accumulated(class_transitions, cohort_exact)

            group_keys = ["intervalo", "grupo_origem", "grupo_destino"]
            group_grouped = pair.groupby(group_keys, dropna=False).agg(
                n_intersecoes=("ano_prodes", "size"),
                area_intersecoes_km2=("area_intersecao_km2", "sum"),
            )
            _add_accumulated(group_transitions, group_grouped)
            cohort_group = pair.groupby(
                ["intervalo", "ano_prodes", "grupo_origem", "grupo_destino"],
                dropna=False,
            ).agg(
                n_intersecoes=("ano_prodes", "size"),
                area_intersecoes_km2=("area_intersecao_km2", "sum"),
            )
            _add_accumulated(group_transitions, cohort_group)

            municipal_keys = [
                "SIGLA_UF",
                "NM_UF",
                "CD_MUN",
                "NM_MUN",
                "ano_prodes",
                "intervalo",
                "grupo_origem",
                "grupo_destino",
            ]
            municipal_grouped = pair.groupby(municipal_keys, dropna=False).agg(
                n_intersecoes=("ano_prodes", "size"),
                area_intersecoes_km2=("area_intersecao_km2", "sum"),
            )
            _add_accumulated(municipal_group_transitions, municipal_grouped)

        complete = codes["classe_primeiro_ano_posterior"].notna()
        if complete.any():
            trajectory = pd.DataFrame(
                {
                    "ano_prodes": year.loc[complete].to_numpy(),
                    "grupo_pre_evento": groups["classe_pre_evento"]
                    .loc[complete]
                    .to_numpy(),
                    "grupo_ano_deteccao": groups["classe_ano_deteccao"]
                    .loc[complete]
                    .to_numpy(),
                    "grupo_pos_evento": groups["classe_primeiro_ano_posterior"]
                    .loc[complete]
                    .to_numpy(),
                    "area_intersecao_km2": area.loc[complete].to_numpy(),
                }
            )
            trajectory["trajetoria"] = (
                trajectory["grupo_pre_evento"].map(CLASS_GROUP_NAMES)
                + " → "
                + trajectory["grupo_ano_deteccao"].map(CLASS_GROUP_NAMES)
                + " → "
                + trajectory["grupo_pos_evento"].map(CLASS_GROUP_NAMES)
            )
            trajectory_keys = [
                "grupo_pre_evento",
                "grupo_ano_deteccao",
                "grupo_pos_evento",
                "trajetoria",
            ]
            all_trajectory = trajectory.groupby(trajectory_keys).agg(
                n_intersecoes=("ano_prodes", "size"),
                area_intersecoes_km2=("area_intersecao_km2", "sum"),
            )
            _add_accumulated(trajectories, all_trajectory)
            by_year_trajectory = trajectory.groupby(
                ["ano_prodes", *trajectory_keys]
            ).agg(
                n_intersecoes=("ano_prodes", "size"),
                area_intersecoes_km2=("area_intersecao_km2", "sum"),
            )
            _add_accumulated(trajectories, by_year_trajectory)

        if processed % (CHUNK_SIZE * 5) == 0:
            logger.info("Analisadas %d linhas da tabela de polígonos", processed)

    if processed == 0:
        raise ValueError(f"Nenhum registro encontrado em {input_path}")
    for interval, _, _ in TRANSITION_PAIRS:
        for detection_year in detection_years:
            summary_stats[(interval, detection_year)]
        summary_stats[(interval, "Todos")]
    logger.info(
        "Análise concluída para %d interseções; %d sem mapa posterior",
        processed,
        missing_post,
    )
    return (
        class_transitions,
        group_transitions,
        municipal_group_transitions,
        summary_stats,
        trajectories,
    )


def _make_summary_frame(summary_stats):
    records = []
    for (interval, detection_year), stats in summary_stats.items():
        valid_count = stats["valid"]
        valid_area = stats["area_valid"]
        records.append(
            {
                "intervalo": interval,
                "ano_prodes": detection_year,
                "n_intersecoes_validas": valid_count,
                "n_classes_dominantes_alteradas": stats["class_changed"],
                "percentual_intersecoes_classes_alteradas": (
                    stats["class_changed"] / valid_count * 100
                    if valid_count
                    else 0.0
                ),
                "n_grupos_alterados": stats["group_changed"],
                "percentual_intersecoes_grupos_alterados": (
                    stats["group_changed"] / valid_count * 100
                    if valid_count
                    else 0.0
                ),
                "area_intersecoes_validas_km2": valid_area,
                "area_intersecoes_classes_alteradas_km2": stats[
                    "area_class_changed"
                ],
                "percentual_area_intersecoes_classes_alteradas": (
                    stats["area_class_changed"] / valid_area * 100
                    if valid_area
                    else 0.0
                ),
                "area_intersecoes_grupos_alterados_km2": stats[
                    "area_group_changed"
                ],
                "percentual_area_intersecoes_grupos_alterados": (
                    stats["area_group_changed"] / valid_area * 100
                    if valid_area
                    else 0.0
                ),
            }
        )
    return pd.DataFrame(records, columns=SUMMARY_COLUMNS)


def _main_outputs(aggregates):
    (
        class_totals,
        group_totals,
        municipal_totals,
        summary_stats,
        trajectory_totals,
    ) = aggregates
    summary = _make_summary_frame(summary_stats)
    denominators = {
        (row.intervalo, row.ano_prodes): (
            row.n_intersecoes_validas,
            row.area_intersecoes_validas_km2,
        )
        for row in summary.itertuples(index=False)
    }

    class_records = []
    for key, values in class_totals.items():
        if len(key) == 7:
            interval, source_code, source_name, source_group, target_code, target_name, target_group = key
            cohort_year = "Todos"
        else:
            (
                interval,
                cohort_year,
                source_code,
                source_name,
                source_group,
                target_code,
                target_name,
                target_group,
            ) = key
        class_records.append(
            {
                "intervalo": interval,
                "ano_prodes": cohort_year,
                "codigo_classe_origem": source_code,
                "classe_origem": source_name,
                "grupo_origem": CLASS_GROUP_NAMES[source_group],
                "codigo_classe_destino": target_code,
                "classe_destino": target_name,
                "grupo_destino": CLASS_GROUP_NAMES[target_group],
                "n_intersecoes": values[0],
                "area_intersecoes_km2": values[1],
            }
        )
    class_frame = pd.DataFrame(class_records)
    class_frame["ano_prodes"] = class_frame["ano_prodes"].astype(str)
    class_frame["percentual_intersecoes"] = 0.0
    class_frame["percentual_area_intersecoes"] = 0.0
    for idx, row in class_frame.iterrows():
        cohort_year = row.ano_prodes
        denominator_year = "Todos" if cohort_year == "Todos" else int(cohort_year)
        denominator = denominators.get((row.intervalo, denominator_year), (0, 0.0))
        class_frame.loc[idx, "percentual_intersecoes"] = (
            row.n_intersecoes / denominator[0] * 100 if denominator[0] else 0.0
        )
        class_frame.loc[idx, "percentual_area_intersecoes"] = (
            row.area_intersecoes_km2 / denominator[1] * 100
            if denominator[1]
            else 0.0
        )
    class_frame = class_frame[CLASS_TRANSITION_COLUMNS].sort_values(
        ["intervalo", "ano_prodes", "n_intersecoes"],
        ascending=[True, True, False],
    )

    group_records = []
    for key, values in group_totals.items():
        if len(key) == 3:
            interval, source_group, target_group = key
            cohort_year = "Todos"
        else:
            interval, cohort_year, source_group, target_group = key
        group_records.append(
            {
                "intervalo": interval,
                "ano_prodes": cohort_year,
                "grupo_origem": CLASS_GROUP_NAMES[source_group],
                "grupo_destino": CLASS_GROUP_NAMES[target_group],
                "n_intersecoes": values[0],
                "area_intersecoes_km2": values[1],
            }
        )
    group_frame = pd.DataFrame(group_records)
    group_frame["ano_prodes"] = group_frame["ano_prodes"].astype(str)
    group_frame["percentual_intersecoes"] = 0.0
    group_frame["percentual_area_intersecoes"] = 0.0
    for idx, row in group_frame.iterrows():
        denominator_year = "Todos" if row.ano_prodes == "Todos" else int(row.ano_prodes)
        denominator = denominators.get((row.intervalo, denominator_year), (0, 0.0))
        group_frame.loc[idx, "percentual_intersecoes"] = (
            row.n_intersecoes / denominator[0] * 100 if denominator[0] else 0.0
        )
        group_frame.loc[idx, "percentual_area_intersecoes"] = (
            row.area_intersecoes_km2 / denominator[1] * 100
            if denominator[1]
            else 0.0
        )
    group_frame = group_frame[GROUP_TRANSITION_COLUMNS].sort_values(
        ["intervalo", "ano_prodes", "n_intersecoes"],
        ascending=[True, True, False],
    )

    local_records = []
    for key, values in municipal_totals.items():
        local_records.append(
            {
                **dict(
                    zip(
                        MUNICIPAL_TRANSITION_COLUMNS[:8],
                        key,
                    )
                ),
                "grupo_origem": CLASS_GROUP_NAMES[key[6]],
                "grupo_destino": CLASS_GROUP_NAMES[key[7]],
                "n_intersecoes": values[0],
                "area_intersecoes_km2": values[1],
            }
        )
    local_frame = pd.DataFrame(local_records, columns=MUNICIPAL_TRANSITION_COLUMNS)
    local_denominators = local_frame.groupby(
        ["SIGLA_UF", "CD_MUN", "ano_prodes", "intervalo"]
    ).agg(
        total_intersecoes=("n_intersecoes", "sum"),
        total_area=("area_intersecoes_km2", "sum"),
    )
    local_frame = local_frame.merge(
        local_denominators,
        on=["SIGLA_UF", "CD_MUN", "ano_prodes", "intervalo"],
        how="left",
        validate="many_to_one",
    )
    local_frame["percentual_intersecoes"] = (
        local_frame["n_intersecoes"] / local_frame["total_intersecoes"] * 100
    )
    local_frame["percentual_area_intersecoes"] = (
        local_frame["area_intersecoes_km2"] / local_frame["total_area"] * 100
    )
    local_frame = local_frame.drop(columns=["total_intersecoes", "total_area"])

    trajectory_records = []
    for key, values in trajectory_totals.items():
        if len(key) == 4:
            pre_group, event_group, post_group, label = key
            cohort_year = "Todos"
        else:
            cohort_year, pre_group, event_group, post_group, label = key
        trajectory_records.append(
            {
                "ano_prodes": cohort_year,
                "grupo_pre_evento": CLASS_GROUP_NAMES[pre_group],
                "grupo_ano_deteccao": CLASS_GROUP_NAMES[event_group],
                "grupo_pos_evento": CLASS_GROUP_NAMES[post_group],
                "trajetoria": label,
                "n_intersecoes": values[0],
                "area_intersecoes_km2": values[1],
            }
        )
    trajectory_frame = pd.DataFrame(trajectory_records)
    trajectory_frame["ano_prodes"] = trajectory_frame["ano_prodes"].astype(str)
    trajectory_frame["percentual_intersecoes"] = 0.0
    trajectory_frame["percentual_area_intersecoes"] = 0.0
    trajectory_denominators = {}
    for detection_year in year_counts_from_trajectory(trajectory_frame):
        subset = trajectory_frame.loc[
            trajectory_frame["ano_prodes"] == detection_year
        ]
        trajectory_denominators[detection_year] = (
            subset["n_intersecoes"].sum(),
            subset["area_intersecoes_km2"].sum(),
        )
    for idx, row in trajectory_frame.iterrows():
        denominator_year = "Todos" if row.ano_prodes == "Todos" else row.ano_prodes
        denominator = trajectory_denominators.get(denominator_year, (0, 0.0))
        trajectory_frame.loc[idx, "percentual_intersecoes"] = (
            row.n_intersecoes / denominator[0] * 100 if denominator[0] else 0.0
        )
        trajectory_frame.loc[idx, "percentual_area_intersecoes"] = (
            row.area_intersecoes_km2 / denominator[1] * 100
            if denominator[1]
            else 0.0
        )
    trajectory_frame = trajectory_frame[TRAJECTORY_COLUMNS].sort_values(
        ["ano_prodes", "n_intersecoes"],
        ascending=[True, False],
    )
    return summary, class_frame, group_frame, local_frame, trajectory_frame


def year_counts_from_trajectory(frame):
    return frame["ano_prodes"].drop_duplicates().tolist()


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    input_path = config.MAPBIOMAS_TRANSITIONS_CSV_PATH
    if not input_path.is_file():
        raise FileNotFoundError(f"Tabela de transições não encontrada: {input_path}")
    aggregates = _aggregate_event_transitions(input_path)
    output_paths = [
        config.MAPBIOMAS_EVENT_TRANSITION_SUMMARY_PATH,
        config.MAPBIOMAS_EVENT_CLASS_TRANSITIONS_PATH,
        config.MAPBIOMAS_EVENT_GROUP_TRANSITIONS_PATH,
        config.MAPBIOMAS_EVENT_MUNICIPAL_TRANSITIONS_PATH,
        config.MAPBIOMAS_EVENT_GROUP_TRAJECTORIES_PATH,
    ]
    output_paths[0].parent.mkdir(parents=True, exist_ok=True)
    temporary_paths = {
        output_path: output_path.with_name(
            f"{output_path.name}.{os.getpid()}.tmp"
        )
        for output_path in output_paths
    }
    try:
        frames = _main_outputs(aggregates)
        for frame, output_path in zip(frames, output_paths):
            frame.to_csv(
                temporary_paths[output_path],
                index=False,
                encoding="utf-8",
                compression="gzip" if output_path.suffix == ".gz" else None,
            )
        for output_path, temporary_path in temporary_paths.items():
            os.replace(temporary_path, output_path)
        logger.info(
            "Tabelas de análise de transições gravadas em %s",
            config.OUTPUT_DIR,
        )
    finally:
        for temporary_path in temporary_paths.values():
            if temporary_path.exists():
                temporary_path.unlink()


if __name__ == "__main__":
    main()
