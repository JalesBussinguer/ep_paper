"""
Robustez do Getis-Ord Gi*: autocorrelação global (Moran's I), sensibilidade à
definição de vizinhança, correção de testes múltiplos (FDR de Benjamini-
Hochberg) e estabilidade temporal dos hotspots (períodos e anos).

Uso:
    python -m src.gi_robustness
"""

import logging

import numpy as np
import pandas as pd
from esda.moran import Moran
from libpysal.weights import DistanceBand, KNN, Queen, higher_order, w_union
from libpysal.weights.util import min_threshold_distance
from scipy import stats

from . import aggregation, config, hotspots, io_utils

logger = logging.getLogger(__name__)

FDR_LEVELS = {0.10: 90, 0.05: 95, 0.01: 99}
PERIODS = {"2020-2022": (2020, 2022), "2023-2025": (2023, 2025)}


def load_municipal_frame():
    """Municípios recortados pelo Cerrado + áreas desmatadas por ano."""
    boundary = io_utils.load_cerrado_boundary()
    municipalities = aggregation.clip_municipalities_to_biome(
        io_utils.load_municipalities(), boundary
    ).reset_index(drop=True)
    yearly = pd.read_csv(config.YEARLY_MUNICIPAL_AREA_PATH, dtype={"CD_MUN": "string"})
    return municipalities, yearly


def attach_area(municipalities, yearly, first_year, last_year):
    subset = yearly[(yearly["year"] >= first_year) & (yearly["year"] <= last_year)]
    area = subset.groupby("CD_MUN")["area_km2"].sum()
    frame = municipalities.copy()
    frame["area_desmatada_km2"] = frame["CD_MUN"].map(area).fillna(0.0)
    frame["pct_area_desmatada"] = 100 * frame["area_desmatada_km2"] / frame["AREA_CERRADO_KM2"]
    return frame


def build_weights(municipalities):
    """Variantes de vizinhança, todas com pesos binários e Gi* com star=True."""
    queen = hotspots.build_contiguity_weights(municipalities)
    variants = {"Queen (base)": queen}

    queen_raw = Queen.from_dataframe(municipalities, use_index=False, silence_warnings=True)
    order_two = w_union(queen_raw, higher_order(queen_raw, 2, silence_warnings=True))
    order_two.transform = "B"
    variants["Queen de 2ª ordem"] = order_two

    knn = KNN.from_dataframe(municipalities, k=8, use_index=False)
    knn.transform = "B"
    variants["8 vizinhos mais próximos"] = knn

    points = np.column_stack(
        [municipalities.geometry.centroid.x, municipalities.geometry.centroid.y]
    )
    threshold = max(150_000.0, float(min_threshold_distance(points)) * 1.001)
    band = DistanceBand(points, threshold=threshold, binary=True, silence_warnings=True)
    band.transform = "B"
    variants[f"Banda de distância ({threshold / 1000:.0f} km)"] = band
    return variants


def two_sided_p(z):
    return 2 * stats.norm.sf(np.abs(z))


def benjamini_hochberg(p_values):
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / (np.arange(len(p_values)) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def classify_fdr(z, p_values):
    """Hotspots após correção FDR (p bilateral derivado do Z score)."""
    adjusted = benjamini_hochberg(two_sided_p(z))
    labels = np.full(len(z), "Não significante", dtype=object)
    resolved = np.zeros(len(z), dtype=bool)
    for alpha, level in sorted(FDR_LEVELS.items()):  # 0.01 primeiro
        hot = ~resolved & (z > 0) & (adjusted <= alpha)
        cold = ~resolved & (z < 0) & (adjusted <= alpha)
        labels[hot] = f"Hot Spot {level}%"
        labels[cold] = f"Cold Spot {level}%"
        resolved |= hot | cold
    return labels, adjusted


def is_hot(labels):
    return pd.Series(labels).astype(str).str.startswith("Hot").to_numpy()


def cohen_kappa(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    observed = (a == b).mean()
    expected = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (observed - expected) / (1 - expected) if expected < 1 else np.nan


def run_gi(frame, column, weights):
    z, p = hotspots.compute_gi_star(frame, column, weights)
    labels = hotspots.classify_hotspots(z, p)
    return z, p, labels


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    config.REFINED_DIR.mkdir(parents=True, exist_ok=True)
    municipalities, yearly = load_municipal_frame()
    years = sorted(yearly["year"].unique())
    full = attach_area(municipalities, yearly, years[0], years[-1])
    variants = build_weights(full)

    # 1. Moran's I global
    moran_rows = []
    queen_weights = hotspots.build_contiguity_weights(full)
    queen_weights.transform = "R"
    for column in ("area_desmatada_km2", "pct_area_desmatada"):
        for transform_name, values in (
            ("original", full[column].to_numpy()),
            ("log(1+x)", np.log1p(full[column].to_numpy())),
        ):
            np.random.seed(config.GI_STAR_SEED)
            moran = Moran(values, queen_weights, permutations=999)
            moran_rows.append(
                dict(
                    variavel=column,
                    transformacao=transform_name,
                    morans_i=moran.I,
                    esperado=moran.EI,
                    z_simulado=moran.z_sim,
                    p_simulado=moran.p_sim,
                )
            )
    pd.DataFrame(moran_rows).to_csv(config.REFINED_DIR / "robustez_gi_moran.csv", index=False, encoding="utf-8")

    # 2. Sensibilidade à vizinhança + FDR
    results = full[["CD_MUN", "NM_MUN", "SIGLA_UF", "area_desmatada_km2", "pct_area_desmatada"]].copy()
    summary_rows = []
    base_hot = {}
    for column, short in (("area_desmatada_km2", "area"), ("pct_area_desmatada", "pct")):
        for name, weights in variants.items():
            z, p, labels = run_gi(full, column, weights)
            fdr_labels, fdr_p = classify_fdr(z, p)
            tag = f"{short}|{name}"
            results[f"z|{tag}"] = z
            results[f"classe|{tag}"] = labels
            results[f"classe_fdr|{tag}"] = fdr_labels
            hot, hot_fdr = is_hot(labels), is_hot(fdr_labels)
            if name == "Queen (base)":
                base_hot[short] = hot
            reference = base_hot[short]
            union = (hot | reference).sum()
            summary_rows.append(
                dict(
                    variavel=short,
                    vizinhanca=name,
                    n_hotspots=int(hot.sum()),
                    n_hotspots_fdr=int(hot_fdr.sum()),
                    n_hotspots_pvalor_perm_90=int(hot.sum()),
                    jaccard_vs_base=(hot & reference).sum() / union if union else np.nan,
                    kappa_vs_base=cohen_kappa(hot, reference),
                    n_cold=int(pd.Series(labels).astype(str).str.startswith("Cold").sum()),
                    vizinhos_medios=float(np.mean(list(weights.cardinalities.values()))),
                )
            )
    pd.DataFrame(summary_rows).to_csv(config.REFINED_DIR / "robustez_gi_vizinhanca.csv", index=False, encoding="utf-8")
    results.to_csv(config.REFINED_DIR / "robustez_gi_municipios.csv", index=False, encoding="utf-8")

    # Reprodução do Gi* original (sanidade)
    original = pd.read_csv(config.HOTSPOT_CSV_PATH, dtype={"CD_MUN": "string"}, usecols=["CD_MUN", "gi_bin"])
    merged = results.merge(original, on="CD_MUN")
    agreement = (merged["classe|area|Queen (base)"] == merged["gi_bin"]).mean()
    logger.info("Concordância com o Gi* original (área): %.4f", agreement)

    # 3. Estabilidade temporal
    queen = variants["Queen (base)"]
    period_rows = []
    period_labels = {}
    for period, (start, end) in PERIODS.items():
        frame = attach_area(municipalities, yearly, start, end)
        for column, short in (("area_desmatada_km2", "area"), ("pct_area_desmatada", "pct")):
            z, p, labels = run_gi(frame, column, queen)
            period_labels[(period, short)] = (z, labels)
            period_rows.append(
                dict(
                    periodo=period,
                    variavel=short,
                    n_hotspots=int(is_hot(labels).sum()),
                    area_total_km2=float(frame["area_desmatada_km2"].sum()),
                )
            )
    persistence = full[["CD_MUN", "NM_MUN", "SIGLA_UF"]].copy()
    for short in ("area", "pct"):
        first = is_hot(period_labels[("2020-2022", short)][1])
        second = is_hot(period_labels[("2023-2025", short)][1])
        persistence[f"status_{short}"] = np.select(
            [first & second, second & ~first, first & ~second],
            ["Persistente", "Emergente (só 2023-2025)", "Desaparece (só 2020-2022)"],
            default="Sem hotspot",
        )
        persistence[f"z_2020_2022_{short}"] = period_labels[("2020-2022", short)][0]
        persistence[f"z_2023_2025_{short}"] = period_labels[("2023-2025", short)][0]
        period_rows.append(
            dict(
                periodo="Sobreposição", variavel=short,
                n_hotspots=int((first & second).sum()),
                area_total_km2=(first & second).sum() / max((first | second).sum(), 1),
            )
        )

    yearly_hot = {}
    for year in years:
        frame = attach_area(municipalities, yearly, year, year)
        z, p, labels = run_gi(frame, "area_desmatada_km2", queen)
        yearly_hot[year] = is_hot(labels)
        period_rows.append(
            dict(periodo=str(year), variavel="area", n_hotspots=int(yearly_hot[year].sum()),
                 area_total_km2=float(frame["area_desmatada_km2"].sum()))
        )
    persistence["anos_como_hotspot"] = np.sum([yearly_hot[y] for y in years], axis=0)
    persistence["area_total_km2"] = full["area_desmatada_km2"].to_numpy()
    persistence.sort_values("anos_como_hotspot", ascending=False).to_csv(
        config.REFINED_DIR / "robustez_gi_estabilidade_municipios.csv", index=False, encoding="utf-8"
    )
    pd.DataFrame(period_rows).to_csv(
        config.REFINED_DIR / "robustez_gi_periodos.csv", index=False, encoding="utf-8"
    )

    # Deslocamento do centro de gravidade da área desmatada
    centroids = municipalities.geometry.centroid.to_crs("EPSG:4326")
    shift_rows = []
    for period, (start, end) in PERIODS.items():
        frame = attach_area(municipalities, yearly, start, end)
        weights_area = frame["area_desmatada_km2"].to_numpy()
        point = (np.average(centroids.x, weights=weights_area), np.average(centroids.y, weights=weights_area))
        shift_rows.append(dict(periodo=period, lon=point[0], lat=point[1]))
    shift = pd.DataFrame(shift_rows)
    shift.to_csv(config.REFINED_DIR / "robustez_gi_centro_gravidade.csv", index=False, encoding="utf-8")
    logger.info("Robustez do Gi* concluída")


if __name__ == "__main__":
    main()
