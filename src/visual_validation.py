"""
Validação visual e concordância PRODES × MapBiomas.

Gera (i) a concordância entre a vegetação natural em t-1 (MapBiomas) e os
polígonos PRODES, (ii) pranchas com painéis MapBiomas pré/detecção/pós para uma
amostra estratificada por tipologia e (iii) uma planilha com links para
conferência manual em imagem de satélite.

Uso:
    python -m src.visual_validation
"""

import logging

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from pyproj import CRS
from rasterio.windows import from_bounds

from . import config
from .figures import GROUP_COLORS, GROUP_ORDER, TYPOLOGY_ORDER, group_lookup
from .pixel_extraction import raster_paths

logger = logging.getLogger(__name__)

PER_TYPOLOGY = 20
PER_SHEET = 3
SHEET_COLUMNS = 3
SHEET_ROWS = 3
PANEL_EXTENT_M = 10_000
SCALE_LENGTH_KM = 2
MIN_AREA_KM2, MAX_AREA_KM2 = 0.05, 3.0
SEED = 20260102
PANEL_CMAP = ListedColormap([GROUP_COLORS[group] for group in GROUP_ORDER])


def agreement_table(trajectories):
    """Fração da área PRODES que era vegetação natural em t-1 segundo o MapBiomas."""
    lookup = group_lookup()
    natural = lookup[trajectories["pre"].to_numpy()] == 0
    frame = trajectories.assign(natural_pre=natural * trajectories["area_w"])
    rows = []
    for keys, part in frame.groupby(["grupo", "year"]):
        rows.append(
            dict(
                grupo=keys[0],
                ano_prodes=keys[1],
                area_estimada_km2=part["area_w"].sum(),
                concordancia_vegetacao_natural_pct=100 * part["natural_pre"].sum() / part["area_w"].sum(),
            )
        )
    for group, part in frame.groupby("grupo"):
        rows.append(
            dict(
                grupo=group,
                ano_prodes="Todos",
                area_estimada_km2=part["area_w"].sum(),
                concordancia_vegetacao_natural_pct=100 * part["natural_pre"].sum() / part["area_w"].sum(),
            )
        )
    return pd.DataFrame(rows)


def stratified_sample(trajectories):
    """Uma tipologia dominante por fragmento; sorteio estratificado."""
    eligible = trajectories[(trajectories["grupo"] == "Hotspot") & (trajectories["year"] <= 2024)]
    dominant = (
        eligible.groupby(["sample_id", "tipologia"])["area_w"].sum().reset_index()
        .sort_values("area_w", ascending=False).drop_duplicates("sample_id")
    )
    sample = pd.read_csv(config.SAMPLE_TABLE_PATH)
    dominant = dominant.merge(sample[["sample_id", "year", "CD_MUN", "area_km2"]], on="sample_id")
    dominant = dominant[dominant["area_km2"].between(MIN_AREA_KM2, MAX_AREA_KM2)]
    picked = [
        part.sample(min(PER_TYPOLOGY, len(part)), random_state=SEED)
        for _, part in dominant.groupby("tipologia")
    ]
    return pd.concat(picked, ignore_index=True)


def read_panel(src, bounds):
    window = from_bounds(*bounds, transform=src.transform)
    data = src.read(1, window=window.round_offsets().round_lengths(), boundless=True, fill_value=0)
    return group_lookup()[data]


def panel_bounds(geometry, crs):
    map_crs = CRS.from_user_input(crs)
    if not map_crs.is_geographic:
        raise ValueError(f"Esperado CRS geográfico para os rasters MapBiomas, recebido {map_crs}.")
    center_x, center_y = geometry.centroid.coords[0]
    geod = map_crs.get_geod()
    half_extent = PANEL_EXTENT_M / 2
    west = geod.fwd(center_x, center_y, 270, half_extent)[0]
    east = geod.fwd(center_x, center_y, 90, half_extent)[0]
    south = geod.fwd(center_x, center_y, 180, half_extent)[1]
    north = geod.fwd(center_x, center_y, 0, half_extent)[1]
    return (west, south, east, north), center_y


def add_page_map_elements(fig, crs):
    scale_left, scale_right, scale_y = 0.375, 0.425, 0.05
    fig.add_artist(Line2D([scale_left, scale_right], [scale_y, scale_y],
                          transform=fig.transFigure, color="black", linewidth=1.2))
    for x in (scale_left, scale_right):
        fig.add_artist(Line2D([x, x], [scale_y - 0.006, scale_y + 0.006],
                              transform=fig.transFigure, color="black", linewidth=0.8))
    fig.text((scale_left + scale_right) / 2, scale_y + 0.008,
             f"{SCALE_LENGTH_KM} km", ha="center", va="bottom", fontsize=7)
    fig.text(0.53, 0.043, "↑ N", ha="center", va="center", fontsize=9, fontweight="bold")
    fig.text(0.63, 0.043, f"CRS: {CRS.from_user_input(crs).to_string()}",
             ha="left", va="center", fontsize=7)


def draw_sheet(rows, geometries, path, title):
    paths = raster_paths()
    year_index = {year: i for i, year in enumerate(config.MAPBIOMAS_YEARS)}
    fig, axes = plt.subplots(SHEET_ROWS, SHEET_COLUMNS, figsize=(15, 15))
    sources = {}
    for position, (_, row) in enumerate(rows.iterrows()):
        geometry = geometries[row["sample_id"]]
        years = [2019, int(row["ano_prodes"]), int(row["ano_prodes"]) + 1]
        labels = ["2019", f"{years[1]} (detecção)", f"{years[2]} (pós)"]
        for k, (year, label) in enumerate(zip(years, labels)):
            if year not in sources:
                sources[year] = rasterio.open(paths[year_index[year]])
            bounds, center_lat = panel_bounds(geometry, sources[year].crs)
            ax = axes[position, k]
            ax.imshow(read_panel(sources[year], bounds), cmap=PANEL_CMAP, vmin=0, vmax=len(GROUP_ORDER) - 1,
                      extent=(bounds[0], bounds[2], bounds[1], bounds[3]), interpolation="nearest")
            parts = [geometry] if geometry.geom_type == "Polygon" else list(geometry.geoms)
            for part in parts:
                x, y = part.exterior.xy
                ax.plot(x, y, color="black", linewidth=0.8)
            ax.set_xlim(bounds[0], bounds[2])
            ax.set_ylim(bounds[1], bounds[3])
            ax.set_aspect(np.cos(np.deg2rad(center_lat)))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"#{row['sample_id']} {label}" if k == 0 else label, fontsize=6.5)
    map_crs = next(iter(sources.values())).crs if sources else config.INPUT_CRS
    for source in sources.values():
        source.close()
    fig.suptitle(title, fontsize=10)
    add_page_map_elements(fig, map_crs)
    fig.subplots_adjust(left=0.035, right=0.965, bottom=0.09, top=0.93, wspace=0.1, hspace=0.16)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    config.VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
    trajectories = pd.read_parquet(config.REFINED_DIR / "_trajetorias_pixel.parquet")
    agreement = agreement_table(trajectories)
    agreement.to_csv(config.REFINED_DIR / "concordancia_prodes_mapbiomas.csv", index=False, encoding="utf-8")

    picked = stratified_sample(trajectories)
    fragments = gpd.read_file(config.SAMPLE_FRAGMENTS_PATH)
    fragments = fragments[fragments["sample_id"].isin(picked["sample_id"])].set_index("sample_id")
    geometries = fragments.geometry.to_dict()
    picked["lon"] = [geometries[i].centroid.x for i in picked["sample_id"]]
    picked["lat"] = [geometries[i].centroid.y for i in picked["sample_id"]]
    picked["link_satelite"] = [
        f"https://www.google.com/maps/@{lat:.6f},{lon:.6f},400m/data=!3m1!1e3"
        for lat, lon in zip(picked["lat"], picked["lon"])
    ]
    picked["revisao_manual"] = ""
    picked = picked.rename(
        columns={"year": "ano_prodes", "tipologia": "tipologia_dominante", "area_w": "area_pixel_km2"}
    )
    picked.to_csv(config.VALIDATION_DIR / "amostra_validacao_visual.csv", index=False, encoding="utf-8")

    for number, category in enumerate(TYPOLOGY_ORDER, start=1):
        part = picked[picked["tipologia_dominante"] == category].reset_index(drop=True)
        for page, start in enumerate(range(0, len(part), PER_SHEET), start=1):
            draw_sheet(
                part.iloc[start:start + PER_SHEET], geometries,
                config.VALIDATION_DIR / f"prancha_{number}_{page}.png",
                f"{category} — MapBiomas por grupo (verde: natural; amarelo: antrópico)",
            )
    logger.info("Validação visual salva em %s", config.VALIDATION_DIR)


if __name__ == "__main__":
    main()
