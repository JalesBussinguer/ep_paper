"""
Análise por pixel das trajetórias de classe MapBiomas nos fragmentos PRODES
amostrados (municípios hotspot e de controle).

Entradas: tabela da amostra (`amostra_fragmentos_prodes.csv`), contagens de
trajetórias por pixel (`amostra_contagem_pixels.csv.gz`) e a classificação Gi*
dos municípios. Todas as áreas são estimadas com pesos de Horvitz-Thompson
(1 / probabilidade de inclusão). Intervalos de confiança vêm de um bootstrap
por conglomerado (município), que respeita a dependência espacial dentro de
cada município.

Uso:
    python -m src.pixel_analysis
"""

import logging

import numpy as np
import pandas as pd
from scipy import stats

from . import config
from .mapbiomas_event_transition_analysis import CLASS_GROUP_BY_CODE, CLASS_GROUP_NAMES
from .mapbiomas_legend import format_class_name
from .pixel_extraction import YEARS, decode_keys

logger = logging.getLogger(__name__)

N_BOOT = 500
N_PERM = 2000
SEED = 7
FIRST_YEAR = YEARS[0]

GROUP_KEYS = list(CLASS_GROUP_NAMES)  # vegetação, antrópico, não vegetada, água, sem dados
NATURAL, ANTHROPIC = 0, 1
TYPE_LABELS = [
    "Conversão no ano da detecção",
    "Conversão posterior à detecção",
    "Vegetação persistente",
    "Já antrópico antes do evento",
    "Regeneração",
    "Outras transições",
    "Sem mapa posterior (vegetação até a detecção)",
]
MAIN_TYPE_LABELS = TYPE_LABELS[:6]

DESTINATION_BY_CODE = {
    15: "Pastagem",
    21: "Mosaico de Agricultura e Pastagem",
    19: "Agricultura",
    18: "Agricultura",
    20: "Agricultura",
    36: "Agricultura",
    48: "Agricultura",
    14: "Agricultura",
}
DESTINATION_LABELS = [
    "Pastagem",
    "Mosaico de Agricultura e Pastagem",
    "Agricultura",
    "Outro uso antrópico",
]


def group_lookup():
    lookup = np.full(256, -1, dtype=np.int8)
    for code, key in CLASS_GROUP_BY_CODE.items():
        lookup[code] = GROUP_KEYS.index(key)
    return lookup


def load_inputs():
    sample = pd.read_csv(config.SAMPLE_TABLE_PATH, dtype={"CD_MUN": "string"})
    counts = pd.read_csv(config.SAMPLE_PIXEL_COUNTS_PATH, dtype={"key": "uint64"})
    municipalities = pd.read_csv(
        config.HOTSPOT_CSV_PATH,
        dtype={"CD_MUN": "string"},
        usecols=["CD_MUN", "NM_MUN", "SIGLA_UF", "NM_UF", "gi_bin", "gi_bin_pct"],
    )
    return sample, counts, municipalities


def classify_municipalities(municipalities):
    """Define hotspot (em qualquer variável/nível), nível máximo e status."""
    result = municipalities.copy()

    def level(label):
        label = label.fillna("")
        return np.where(
            label.str.startswith("Hot Spot"),
            label.str.extract(r"(\d+)%")[0].astype(float).fillna(0),
            0.0,
        )

    result["hot_nivel_area"] = level(result["gi_bin"])
    result["hot_nivel_pct"] = level(result["gi_bin_pct"])
    result["hot_nivel"] = result[["hot_nivel_area", "hot_nivel_pct"]].max(axis=1)
    cold = (
        result["gi_bin"].fillna("").str.startswith("Cold")
        | result["gi_bin_pct"].fillna("").str.startswith("Cold")
    )
    result["grupo"] = np.where(result["hot_nivel"] > 0, "Hotspot", "Controle")
    result["status"] = np.select(
        [result["hot_nivel"] > 0, cold],
        ["Hotspot", "Cold spot"],
        default="Não significante",
    )
    return result


def build_pixel_table(sample, counts, municipalities):
    """Expande contagens em linhas (fragmento, trajetória) com área ponderada."""
    classes = decode_keys(counts["key"].to_numpy())
    lookup = group_lookup()
    unknown = np.unique(classes[lookup[classes] < 0])
    if len(unknown):
        raise ValueError(f"Códigos MapBiomas sem grupo analítico: {unknown.tolist()}")

    table = counts[["sample_id", "n_pixels"]].copy()
    table = table.merge(
        sample[["sample_id", "year", "CD_MUN", "peso", "pixel_area_km2"]],
        on="sample_id",
        validate="many_to_one",
    )
    table = table.merge(municipalities, on="CD_MUN", how="left", validate="many_to_one")
    if table["grupo"].isna().any():
        raise ValueError("Municípios da amostra sem classificação Gi*")
    table["area_km2"] = table["n_pixels"] * table["pixel_area_km2"]
    table["area_w"] = table["area_km2"] * table["peso"]
    return table.reset_index(drop=True), classes


def class_at(classes, years, offset=None, fixed_year=None):
    """Classe de cada trajetória em `ano + offset` (ou em um ano fixo)."""
    n = len(years)
    index = (np.full(n, fixed_year, dtype=int) if fixed_year is not None
             else years + offset) - FIRST_YEAR
    valid = (index >= 0) & (index < len(YEARS))
    out = np.zeros(n, dtype=np.uint8)
    out[valid] = classes[np.flatnonzero(valid), index[valid]]
    return out, valid


def typology(pre, det, post, post_missing):
    """Tipologia da trajetória por pixel, a partir de grupos (índices)."""
    n_, a_ = NATURAL, ANTHROPIC
    conditions = [
        (pre == n_) & (det == a_),
        (pre == n_) & (det == n_) & (post == a_) & ~post_missing,
        (pre == n_) & (det == n_) & (post == n_) & ~post_missing,
        (pre == a_) & (det == a_) & ((post != n_) | post_missing),
        ((pre == a_) & (det == n_)) | ((pre == a_) & (det == a_) & (post == n_) & ~post_missing),
        (pre == n_) & (det == n_) & post_missing,
    ]
    return np.select(conditions, [0, 1, 2, 3, 4, 6], default=5)


def add_reference_columns(table, classes, pre_mode="t-1", post_mode="t+1"):
    lookup = group_lookup()
    years = table["year"].to_numpy(dtype=int)
    if pre_mode == "t-1":
        pre_cls, _ = class_at(classes, years, offset=-1)
    else:
        pre_cls, _ = class_at(classes, years, fixed_year=FIRST_YEAR)
    det_cls, _ = class_at(classes, years, offset=0)
    post_offset = 1 if post_mode == "t+1" else 2
    post_cls, post_valid = class_at(classes, years, offset=post_offset)
    groups = {
        "pre": lookup[pre_cls],
        "det": lookup[det_cls],
        "post": lookup[post_cls],
    }
    labels = typology(groups["pre"], groups["det"], groups["post"], ~post_valid)
    return pre_cls, det_cls, post_cls, post_valid, groups, labels


# ---------------------------------------------------------------- bootstrap

def _cluster_matrix(frame, category_col, categories, cluster_col="CD_MUN"):
    pivot = frame.pivot_table(
        index=cluster_col, columns=category_col, values="area_w",
        aggfunc="sum", fill_value=0.0,
    ).reindex(columns=categories, fill_value=0.0)
    return pivot.to_numpy()


def bootstrap_draws(matrix, rng, n_boot=N_BOOT):
    """Reamostra municípios com reposição; retorna somas por categoria (B x K)."""
    m = matrix.shape[0]
    counts = rng.multinomial(m, np.full(m, 1.0 / m), size=n_boot)
    return counts @ matrix


def summarize_shares(frame, category_col, categories, by, rng):
    """Área, participação (%) e IC95% bootstrap por subgrupo e categoria."""
    rows = []
    groups = frame.groupby(by, observed=True) if by else [((), frame)]
    for key, part in groups:
        key = key if isinstance(key, tuple) else (key,)
        matrix = _cluster_matrix(part, category_col, categories)
        totals = matrix.sum(axis=0)
        grand = totals.sum()
        draws = bootstrap_draws(matrix, rng)
        draw_shares = 100 * draws / draws.sum(axis=1, keepdims=True)
        low, high = np.percentile(draw_shares, [2.5, 97.5], axis=0)
        for index, category in enumerate(categories):
            record = dict(zip(by or [], key))
            record.update(
                categoria=category,
                area_estimada_km2=totals[index],
                participacao_pct=100 * totals[index] / grand if grand else np.nan,
                ic95_inf=low[index],
                ic95_sup=high[index],
                n_fragmentos_amostra=part["sample_id"].nunique(),
                n_municipios=part["CD_MUN"].nunique(),
            )
            rows.append(record)
    return pd.DataFrame(rows)


def cramers_v(table):
    table = np.asarray(table, dtype=float)
    keep_rows = table.sum(axis=1) > 0
    keep_cols = table.sum(axis=0) > 0
    table = table[keep_rows][:, keep_cols]
    if min(table.shape) < 2:
        return np.nan
    total = table.sum()
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / total
    chi2 = ((table - expected) ** 2 / expected).sum()
    return float(np.sqrt(chi2 / (total * (min(table.shape) - 1))))


def compare_groups(frame, category_col, categories, rng):
    """Hotspot x controle: diferenças, V de Cramér, qui-quadrado e permutação."""
    matrices = {
        g: _cluster_matrix(frame[frame["grupo"] == g], category_col, categories)
        for g in ("Hotspot", "Controle")
    }
    draws = {g: bootstrap_draws(m, rng) for g, m in matrices.items()}
    shares = {
        g: 100 * d / d.sum(axis=1, keepdims=True) for g, d in draws.items()
    }
    observed = {g: 100 * m.sum(axis=0) / m.sum() for g, m in matrices.items()}
    rows = []
    for index, category in enumerate(categories):
        diff_draws = shares["Hotspot"][:, index] - shares["Controle"][:, index]
        low, high = np.percentile(diff_draws, [2.5, 97.5])
        rows.append(
            dict(
                categoria=category,
                hotspot_pct=observed["Hotspot"][index],
                controle_pct=observed["Controle"][index],
                diferenca_pp=observed["Hotspot"][index] - observed["Controle"][index],
                ic95_inf=low,
                ic95_sup=high,
            )
        )
    result = pd.DataFrame(rows)

    table = np.vstack([matrices["Hotspot"].sum(axis=0), matrices["Controle"].sum(axis=0)])
    v_observed = cramers_v(table)
    v_draws = [
        cramers_v(np.vstack([draws["Hotspot"][b], draws["Controle"][b]]))
        for b in range(N_BOOT)
    ]
    # Permutação dos rótulos entre municípios (preserva o tamanho dos grupos).
    all_matrix = np.vstack([matrices["Hotspot"], matrices["Controle"]])
    n_hot = matrices["Hotspot"].shape[0]
    permuted = []
    for _ in range(N_PERM):
        order = rng.permutation(all_matrix.shape[0])
        permuted.append(
            cramers_v(np.vstack([
                all_matrix[order[:n_hot]].sum(axis=0),
                all_matrix[order[n_hot:]].sum(axis=0),
            ]))
        )
    permuted = np.asarray(permuted)
    p_perm = (1 + np.sum(permuted >= v_observed)) / (1 + N_PERM)

    weights = frame["peso"].to_numpy()
    n_eff = weights.sum() ** 2 / (weights ** 2).sum()
    scaled = table / table.sum() * n_eff
    keep = scaled.sum(axis=0) > 0
    chi2, p_chi2, dof, _ = stats.chi2_contingency(scaled[:, keep])
    summary = pd.DataFrame(
        [
            dict(
                n_municipios_hotspot=matrices["Hotspot"].shape[0],
                n_municipios_controle=matrices["Controle"].shape[0],
                v_cramer=v_observed,
                v_cramer_ic95_inf=np.nanpercentile(v_draws, 2.5),
                v_cramer_ic95_sup=np.nanpercentile(v_draws, 97.5),
                p_permutacao_municipios=p_perm,
                qui2_kish=chi2,
                gl=dof,
                p_qui2_kish=p_chi2,
                n_efetivo_kish=n_eff,
            )
        ]
    )
    return result, summary


# ---------------------------------------------------------------- análises

def transition_matrices(table, classes, rng):
    """Matrizes de transição por pixel (grupos e classes) em km2 ponderados."""
    lookup = group_lookup()
    years = table["year"].to_numpy(dtype=int)
    pairs = {
        "pré-evento (t-1) → detecção (t)": (-1, 0),
        "detecção (t) → pós-evento (t+1)": (0, 1),
        "pré-evento (t-1) → pós-evento (t+1)": (-1, 1),
    }
    group_rows, class_rows = [], []
    for interval, (offset_from, offset_to) in pairs.items():
        origin, valid_from = class_at(classes, years, offset=offset_from)
        target, valid_to = class_at(classes, years, offset=offset_to)
        valid = valid_from & valid_to
        part = table.loc[valid].copy()
        part["origem"] = origin[valid]
        part["destino"] = target[valid]
        part["g_origem"] = lookup[part["origem"]]
        part["g_destino"] = lookup[part["destino"]]
        for scope, scope_frame in [("Todos", part)] + [
            (g, part[part["grupo"] == g]) for g in ("Hotspot", "Controle")
        ]:
            total = scope_frame["area_w"].sum()
            grouped = scope_frame.groupby(["g_origem", "g_destino"])["area_w"].sum()
            for (go, gd), area in grouped.items():
                group_rows.append(
                    dict(
                        intervalo=interval,
                        escopo=scope,
                        grupo_origem=CLASS_GROUP_NAMES[GROUP_KEYS[go]],
                        grupo_destino=CLASS_GROUP_NAMES[GROUP_KEYS[gd]],
                        area_estimada_km2=area,
                        participacao_pct=100 * area / total,
                    )
                )
            exact = scope_frame.groupby(["origem", "destino"])["area_w"].sum()
            for (co, cd), area in exact.nlargest(25).items():
                class_rows.append(
                    dict(
                        intervalo=interval,
                        escopo=scope,
                        classe_origem=format_class_name(co),
                        classe_destino=format_class_name(cd),
                        area_estimada_km2=area,
                        participacao_pct=100 * area / total,
                    )
                )
    return pd.DataFrame(group_rows), pd.DataFrame(class_rows)


def sensitivity_table(table, classes, rng):
    rows = []
    for pre_mode in ("t-1", "2019"):
        for post_mode in ("t+1", "t+2"):
            *_, labels = add_reference_columns(table, classes, pre_mode, post_mode)
            scenario = table.assign(tipologia=np.array(TYPE_LABELS)[labels])
            years = scenario["year"].to_numpy()
            if post_mode == "t+2":
                scenario = scenario[years <= YEARS[-1] - 2]
            else:
                scenario = scenario[years <= YEARS[-1] - 1]
            if pre_mode == "2019":
                scenario = scenario[scenario["year"] >= FIRST_YEAR + 1]
            summary = summarize_shares(
                scenario, "tipologia", MAIN_TYPE_LABELS, ["grupo"], rng
            )
            summary["referencia_pre"] = pre_mode
            summary["referencia_pos"] = post_mode
            rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def time_to_conversion(table, classes):
    """Defasagem entre a detecção e o 1º ano antrópico, p/ pixels naturais em t-1."""
    lookup = group_lookup()
    years = table["year"].to_numpy(dtype=int)
    pre, _ = class_at(classes, years, offset=-1)
    natural_before = lookup[pre] == NATURAL
    anthropic = lookup[classes] == ANTHROPIC  # (n, anos)
    year_index = years - FIRST_YEAR
    columns = np.arange(len(YEARS))[None, :]
    eligible = anthropic & (columns >= year_index[:, None])
    converted = eligible.any(axis=1)
    first = np.where(converted, eligible.argmax(axis=1) - year_index, -1)
    destination_class = classes[np.arange(len(classes)), np.clip(year_index + first, 0, len(YEARS) - 1)]
    return natural_before, converted, first, destination_class


def lag_table(table, natural_before, converted, first):
    rows = []
    base = table.loc[natural_before]
    for (cohort, group), part in base.groupby(["year", "grupo"]):
        total = part["area_w"].sum()
        max_lag = YEARS[-1] - cohort
        part_first = first[part.index]
        part_conv = converted[part.index]
        for lag in range(max_lag + 1):
            area = part.loc[part_conv & (part_first <= lag), "area_w"].sum()
            rows.append(
                dict(
                    ano_prodes=cohort,
                    grupo=group,
                    defasagem_anos=lag,
                    area_base_km2=total,
                    area_convertida_acumulada_km2=area,
                    taxa_conversao_acumulada_pct=100 * area / total,
                )
            )
    return pd.DataFrame(rows)


def destination_label(codes):
    return np.array([DESTINATION_BY_CODE.get(int(c), "Outro uso antrópico") for c in codes])


def municipal_table(table, natural_before, converted, destination_class, labels, post_valid):
    frame = table.assign(
        tipologia=np.array(TYPE_LABELS)[labels],
        natural_antes=natural_before,
        convertido=converted & natural_before,
    )
    frame["destino"] = np.where(frame["convertido"], destination_label(destination_class), "")
    rows = []
    for code, part in frame.groupby("CD_MUN"):
        natural_area = part.loc[part["natural_antes"], "area_w"].sum()
        converted_area = part.loc[part["convertido"], "area_w"].sum()
        by_destination = part[part["convertido"]].groupby("destino")["area_w"].sum()
        record = dict(
            CD_MUN=code,
            NM_MUN=part["NM_MUN"].iloc[0],
            SIGLA_UF=part["SIGLA_UF"].iloc[0],
            grupo=part["grupo"].iloc[0],
            hot_nivel=part["hot_nivel"].iloc[0],
            n_fragmentos_amostra=part["sample_id"].nunique(),
            area_total_estimada_km2=part["area_w"].sum(),
            area_natural_pre_km2=natural_area,
            area_convertida_km2=converted_area,
            taxa_conversao_pct=100 * converted_area / natural_area if natural_area else np.nan,
            pct_ja_antropico=100
            * part.loc[part["tipologia"] == TYPE_LABELS[3], "area_w"].sum()
            / part["area_w"].sum(),
            destino_predominante=by_destination.idxmax() if len(by_destination) else "",
        )
        for label in DESTINATION_LABELS:
            record[f"pct_conv_{label}"] = (
                100 * by_destination.get(label, 0.0) / converted_area if converted_area else np.nan
            )
        record["estimativa_confiavel"] = record["n_fragmentos_amostra"] >= 30
        rows.append(record)
    return pd.DataFrame(rows).sort_values("area_convertida_km2", ascending=False)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    config.REFINED_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    sample, counts, municipalities = load_inputs()
    municipalities = classify_municipalities(municipalities)
    table, classes = build_pixel_table(sample, counts, municipalities)
    logger.info("Tabela de trajetórias: %d linhas", len(table))

    pre_cls, det_cls, post_cls, post_valid, groups, labels = add_reference_columns(table, classes)
    table["tipologia"] = np.array(TYPE_LABELS)[labels]
    table["anos_pos"] = post_valid

    # Controle de qualidade da amostra.
    fragment_level = sample.merge(municipalities[["CD_MUN", "grupo"]], on="CD_MUN", how="left")
    qa = (
        fragment_level.assign(
            area_pixel_km2=fragment_level["n_pixels_total"] * fragment_level["pixel_area_km2"]
        )
        .groupby("grupo")
        .apply(
            lambda g: pd.Series(
                dict(
                    n_fragmentos_amostra=len(g),
                    area_populacao_estimada_geom_km2=(g["area_km2"] * g["peso"]).sum(),
                    area_populacao_estimada_pixel_km2=(g["area_pixel_km2"] * g["peso"]).sum(),
                    fragmentos_sem_pixel=int((g["n_pixels_total"] == 0).sum()),
                    n_efetivo_kish=g["peso"].sum() ** 2 / (g["peso"] ** 2).sum(),
                )
            ),
            include_groups=False,
        )
        .reset_index()
    )
    qa.to_csv(config.REFINED_DIR / "qa_amostra.csv", index=False, encoding="utf-8")

    main_frame = table[table["year"] <= YEARS[-1] - 1]  # coortes com mapa posterior

    typology_by_year = summarize_shares(
        table[table["year"] <= YEARS[-1] - 1], "tipologia", MAIN_TYPE_LABELS,
        ["grupo", "year"], rng,
    )
    typology_by_group = summarize_shares(main_frame, "tipologia", MAIN_TYPE_LABELS, ["grupo"], rng)
    typology_by_status = summarize_shares(main_frame, "tipologia", MAIN_TYPE_LABELS, ["status"], rng)
    typology_by_level = summarize_shares(main_frame, "tipologia", MAIN_TYPE_LABELS, ["hot_nivel"], rng)
    typology_by_variable = pd.concat(
        [
            summarize_shares(
                main_frame[main_frame[column] > 0], "tipologia", MAIN_TYPE_LABELS, [], rng
            ).assign(variavel=name)
            for column, name in (("hot_nivel_area", "Hotspot (área absoluta)"),
                                 ("hot_nivel_pct", "Hotspot (% do município)"))
        ]
    )
    typology_2025 = summarize_shares(
        table[table["year"] == YEARS[-1]], "tipologia",
        [TYPE_LABELS[0], TYPE_LABELS[3], TYPE_LABELS[4], TYPE_LABELS[5], TYPE_LABELS[6]],
        ["grupo"], rng,
    )
    typology_by_year.to_csv(config.REFINED_DIR / "tipologia_por_grupo_ano.csv", index=False, encoding="utf-8")
    typology_by_group.to_csv(config.REFINED_DIR / "tipologia_por_grupo.csv", index=False, encoding="utf-8")
    typology_by_status.to_csv(config.REFINED_DIR / "tipologia_por_status_gi.csv", index=False, encoding="utf-8")
    typology_by_level.to_csv(config.REFINED_DIR / "tipologia_por_nivel_confianca.csv", index=False, encoding="utf-8")
    typology_by_variable.to_csv(config.REFINED_DIR / "tipologia_por_variavel_gi.csv", index=False, encoding="utf-8")
    typology_2025.to_csv(config.REFINED_DIR / "tipologia_coorte_2025.csv", index=False, encoding="utf-8")

    comparison, comparison_summary = compare_groups(main_frame, "tipologia", MAIN_TYPE_LABELS, rng)
    comparison.to_csv(config.REFINED_DIR / "comparacao_hotspot_controle.csv", index=False, encoding="utf-8")
    comparison_summary.to_csv(config.REFINED_DIR / "comparacao_hotspot_controle_testes.csv", index=False, encoding="utf-8")

    group_matrix, class_matrix = transition_matrices(table, classes, rng)
    group_matrix.to_csv(config.REFINED_DIR / "matriz_transicao_pixel_grupos.csv", index=False, encoding="utf-8")
    class_matrix.to_csv(config.REFINED_DIR / "matriz_transicao_pixel_classes.csv", index=False, encoding="utf-8")

    sensitivity_table(table, classes, rng).to_csv(
        config.REFINED_DIR / "sensibilidade_referencias_temporais.csv", index=False, encoding="utf-8"
    )

    natural_before, converted, first, destination_class = time_to_conversion(table, classes)
    lag_table(table, natural_before, converted, first).to_csv(
        config.REFINED_DIR / "tempo_ate_conversao.csv", index=False, encoding="utf-8"
    )

    municipal = municipal_table(table, natural_before, converted, destination_class, labels, post_valid)
    municipal.to_csv(config.REFINED_DIR / "destinos_conversao_por_municipio.csv", index=False, encoding="utf-8")

    state = (
        municipal.groupby(["SIGLA_UF", "grupo"])
        .agg(
            area_convertida_km2=("area_convertida_km2", "sum"),
            area_natural_pre_km2=("area_natural_pre_km2", "sum"),
            n_municipios=("CD_MUN", "nunique"),
        )
        .reset_index()
    )
    state["taxa_conversao_pct"] = 100 * state["area_convertida_km2"] / state["area_natural_pre_km2"]
    state.to_csv(config.REFINED_DIR / "conversao_por_estado.csv", index=False, encoding="utf-8")

    # Tabela por pixel para validação e figuras
    table.assign(
        pre=pre_cls, det=det_cls, pos=post_cls,
    )[["sample_id", "year", "CD_MUN", "grupo", "tipologia", "area_w", "area_km2", "pre", "det", "pos"]].to_parquet(
        config.REFINED_DIR / "_trajetorias_pixel.parquet", index=False
    )
    logger.info("Análise por pixel concluída")


if __name__ == "__main__":
    main()
