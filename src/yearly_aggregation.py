"""
Passagem única sobre os polígonos PRODES que gera:

1. a área desmatada por município e ano (insumo da robustez do Gi* por período);
2. uma amostra probabilística (Poisson, estratificada por tamanho) das
   interseções polígono x município, com pesos de inclusão, usada na análise
   por pixel dos municípios hotspot e de controle.

Uso:
    python -m src.yearly_aggregation
"""

import logging
import time

import geopandas as gpd
import numpy as np
import pandas as pd

from . import aggregation, config, io_utils

logger = logging.getLogger(__name__)


def inclusion_probability(area_km2):
    """Probabilidade de inclusão na amostra conforme a classe de tamanho."""
    class_index = np.searchsorted(config.SAMPLE_SIZE_BREAKS_KM2, area_km2, side="left")
    return np.asarray(config.SAMPLE_INCLUSION_PROBABILITIES)[class_index]


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    start = time.time()
    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    boundary = io_utils.load_cerrado_boundary()
    municipalities = aggregation.clip_municipalities_to_biome(
        io_utils.load_municipalities(), boundary
    ).reset_index(drop=True)

    yearly_parts = []
    sample_parts = []
    for chunk_number, intersections in enumerate(
        aggregation.iter_polygon_intersections(
            io_utils.iter_deforestation_polygon_chunks(), municipalities
        ),
        start=1,
    ):
        yearly_parts.append(
            intersections.groupby(["CD_MUN", "year"])
            .agg(area_km2=("area_km2", "sum"), n_poligonos=("prodes_fid", "nunique"))
            .reset_index()
        )
        probability = inclusion_probability(intersections["area_km2"].to_numpy())
        rng = np.random.default_rng([config.SAMPLE_SEED, chunk_number])
        selected = rng.random(len(intersections)) < probability
        if selected.any():
            sample = intersections.loc[selected].copy()
            sample["p_incl"] = probability[selected]
            sample_parts.append(sample)
        logger.info(
            "Bloco %d: %d interseções, %d amostradas",
            chunk_number,
            len(intersections),
            int(selected.sum()),
        )

    yearly = (
        pd.concat(yearly_parts)
        .groupby(["CD_MUN", "year"], as_index=False)
        .agg(area_km2=("area_km2", "sum"), n_poligonos=("n_poligonos", "sum"))
    )
    yearly.to_csv(config.YEARLY_MUNICIPAL_AREA_PATH, index=False, encoding="utf-8")

    sample = pd.concat(sample_parts, ignore_index=True)
    sample = gpd.GeoDataFrame(sample, geometry="geometry", crs=sample_parts[0].crs)
    sample.insert(0, "sample_id", np.arange(len(sample)))
    sample["peso"] = 1.0 / sample["p_incl"]
    sample = sample.to_crs("EPSG:4326")
    sample.to_file(config.SAMPLE_FRAGMENTS_PATH, driver="GPKG")
    sample.drop(columns="geometry").to_csv(
        config.SAMPLE_TABLE_PATH, index=False, encoding="utf-8"
    )
    logger.info(
        "Concluído em %.1f min: %d interseções na amostra (soma de pesos %.0f)",
        (time.time() - start) / 60,
        len(sample),
        sample["peso"].sum(),
    )


if __name__ == "__main__":
    main()
