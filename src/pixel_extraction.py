"""
Extração por pixel das classes MapBiomas (2019-2025) dentro de cada fragmento
amostrado. Para cada fragmento, todos os pixels (regra do centro do pixel) são
lidos nos sete rasters anuais, e as combinações distintas de classes ao longo
da série são contadas. A trajetória de cada pixel fica codificada em um inteiro
de 64 bits (um byte por ano, 2019 no byte menos significativo).

Uso:
    python -m src.pixel_extraction
"""

import logging
import math
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import shapely
from joblib import Parallel, delayed
from rasterio.features import geometry_mask
from rasterio.windows import Window

from . import config

logger = logging.getLogger(__name__)

YEARS = tuple(config.MAPBIOMAS_YEARS)
CHUNK_SIZE = 200


def raster_paths():
    return [
        config.MAPBIOMAS_DIR / config.MAPBIOMAS_RASTER_PATTERN.format(year=year)
        for year in YEARS
    ]


def decode_keys(keys):
    """Converte as chaves de 64 bits em uma matriz (n, anos) de códigos de classe."""
    keys = np.asarray(keys, dtype=np.uint64)
    return np.stack(
        [((keys >> np.uint64(8 * i)) & np.uint64(255)).astype(np.uint8)
         for i in range(len(YEARS))],
        axis=1,
    )


def _process_chunk(records):
    datasets = [rasterio.open(path) for path in raster_paths()]
    try:
        reference = datasets[0]
        transform = reference.transform
        origin_x, origin_y = transform.c, transform.f
        resolution = transform.a
        count_rows, meta_rows = [], []
        for sample_id, wkb in records:
            geometry = shapely.from_wkb(wkb)
            min_x, min_y, max_x, max_y = geometry.bounds
            col_start = max(0, math.floor((min_x - origin_x) / resolution))
            col_end = min(reference.width, math.ceil((max_x - origin_x) / resolution))
            row_start = max(0, math.floor((origin_y - max_y) / resolution))
            row_end = min(reference.height, math.ceil((origin_y - min_y) / resolution))
            latitude = geometry.centroid.y
            pixel_area_km2 = (
                resolution * 111.320 * math.cos(math.radians(latitude))
                * resolution * 110.574
            )
            if col_end <= col_start or row_end <= row_start:
                meta_rows.append((sample_id, 0, pixel_area_km2))
                continue
            window = Window(col_start, row_start, col_end - col_start, row_end - row_start)
            mask = geometry_mask(
                [geometry],
                out_shape=(window.height, window.width),
                transform=reference.window_transform(window),
                all_touched=False,
                invert=True,
            )
            n_pixels = int(mask.sum())
            meta_rows.append((sample_id, n_pixels, pixel_area_km2))
            if n_pixels == 0:
                continue
            keys = np.zeros(n_pixels, dtype=np.uint64)
            for index, dataset in enumerate(datasets):
                values = dataset.read(1, window=window)[mask]
                keys |= values.astype(np.uint64) << np.uint64(8 * index)
            unique_keys, counts = np.unique(keys, return_counts=True)
            count_rows.append(
                pd.DataFrame(
                    {
                        "sample_id": sample_id,
                        "key": unique_keys,
                        "n_pixels": counts.astype(np.int64),
                    }
                )
            )
    finally:
        for dataset in datasets:
            dataset.close()
    counts = (
        pd.concat(count_rows, ignore_index=True)
        if count_rows
        else pd.DataFrame(columns=["sample_id", "key", "n_pixels"])
    )
    meta = pd.DataFrame(meta_rows, columns=["sample_id", "n_pixels_total", "pixel_area_km2"])
    return counts, meta


def extract_pixel_counts(fragments, n_jobs=12):
    """Retorna (contagens por trajetória, metadados por fragmento)."""
    fragments = fragments.to_crs("EPSG:4326")
    wkbs = shapely.to_wkb(fragments.geometry.array)
    ids = fragments["sample_id"].to_numpy()
    # Fragmentos grandes primeiro, para equilibrar a carga entre processos.
    order = np.argsort(-fragments["area_km2"].to_numpy())
    records = [(int(ids[i]), wkbs[i]) for i in order]
    chunks = [records[i:i + CHUNK_SIZE] for i in range(0, len(records), CHUNK_SIZE)]
    logger.info("Extraindo %d fragmentos em %d blocos", len(records), len(chunks))
    start = time.time()
    results = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_process_chunk)(chunk) for chunk in chunks
    )
    logger.info("Extração concluída em %.1f min", (time.time() - start) / 60)
    counts = pd.concat([r[0] for r in results], ignore_index=True)
    meta = pd.concat([r[1] for r in results], ignore_index=True)
    return counts, meta


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    fragments = gpd.read_file(config.SAMPLE_FRAGMENTS_PATH)
    counts, meta = extract_pixel_counts(fragments)
    counts.to_csv(config.SAMPLE_PIXEL_COUNTS_PATH, index=False, compression="gzip")
    sample = pd.read_csv(config.SAMPLE_TABLE_PATH, dtype={"CD_MUN": "string"})
    sample = sample.merge(meta, on="sample_id", validate="one_to_one")
    sample.to_csv(config.SAMPLE_TABLE_PATH, index=False, encoding="utf-8")
    logger.info("Contagens: %d linhas; %d fragmentos sem pixel", len(counts),
                int((meta["n_pixels_total"] == 0).sum()))


if __name__ == "__main__":
    main()
