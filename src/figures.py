"""
Figuras finais do artigo (resolução de publicação, rótulos em português).

Uso:
    python -m src.figures
"""

import logging

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Patch
from matplotlib.path import Path
from matplotlib.ticker import FuncFormatter

from . import config, io_utils
from .mapbiomas_event_transition_analysis import CLASS_GROUP_BY_CODE, CLASS_GROUP_NAMES

logger = logging.getLogger(__name__)

DPI = 300
SIMPLIFY_M = 1000
GROUP_ORDER = ["vegetacao_natural", "uso_antrópico", "superficie_nao_vegetada", "agua", "sem_dados"]
GROUP_COLORS = {
    "vegetacao_natural": "#2E7D32",
    "uso_antrópico": "#E6A100",
    "superficie_nao_vegetada": "#9E9E9E",
    "agua": "#1E88E5",
    "sem_dados": "#D0D0D0",
}
TYPOLOGY_ORDER = [
    "Conversão no ano da detecção",
    "Conversão posterior à detecção",
    "Vegetação persistente",
    "Já antrópico antes do evento",
    "Regeneração",
    "Outras transições",
]
TYPOLOGY_WRAPPED = [
    "Conversão no\nano da detecção",
    "Conversão\nposterior",
    "Vegetação\npersistente",
    "Já antrópico\nantes do evento",
    "Regeneração",
    "Outras\ntransições",
]
STAGE_LABELS = ["Pré-evento (t−1)", "Ano da detecção (t)", "Pós-evento (t+1)"]
HOT_COLOR, CONTROL_COLOR = "#B2182B", "#4D6A7A"
GI_COLORS = {
    "Hot Spot 99%": "#B2182B",
    "Hot Spot 95%": "#EF6548",
    "Hot Spot 90%": "#FDBB84",
    "Não significante": "#EEEEEE",
    "Cold Spot 90%": "#C6DBEF",
    "Cold Spot 95%": "#6BAED6",
    "Cold Spot 99%": "#2171B5",
}
GI_LABELS = {
    "Hot Spot 99%": "Hot spot (99%)",
    "Hot Spot 95%": "Hot spot (95%)",
    "Hot Spot 90%": "Hot spot (90%)",
    "Não significante": "Não significante",
    "Cold Spot 90%": "Cold spot (90%)",
    "Cold Spot 95%": "Cold spot (95%)",
    "Cold Spot 99%": "Cold spot (99%)",
}


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.dpi": 100,
            "savefig.dpi": DPI,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
        }
    )


def pt_number(value, decimals=0):
    """Formata números no padrão brasileiro (ponto de milhar e vírgula decimal)."""
    text = f"{value:,.{decimals}f}"
    return text.replace(",", "§").replace(".", ",").replace("§", ".")


def pt_formatter(decimals=0):
    return FuncFormatter(lambda value, _: pt_number(value, decimals))


def panel_label(ax, letter):
    ax.text(0.0, 1.04, f"({letter})", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def save(fig, name):
    path = config.FIGURES_DIR / name
    fig.savefig(path)
    plt.close(fig)
    logger.info("Figura salva: %s", path.name)


def group_lookup():
    table = np.full(256, len(GROUP_ORDER) - 1, dtype=np.int8)
    for code, group in CLASS_GROUP_BY_CODE.items():
        table[code] = GROUP_ORDER.index(group)
    return table


# ----------------------------------------------------------------- mapas base

def load_map_layers():
    """Municípios do Cerrado recortados pelo bioma, contorno dos estados e do bioma."""
    boundary = io_utils.load_cerrado_boundary()
    municipalities = io_utils.load_municipalities()
    municipalities["geometry"] = municipalities.geometry.simplify(SIMPLIFY_M)
    states = municipalities.dissolve(by="SIGLA_UF").reset_index()[["SIGLA_UF", "geometry"]]
    boundary_simple = gpd.GeoSeries([boundary], crs=config.METRIC_CRS).simplify(SIMPLIFY_M).iloc[0]
    clipped = municipalities.copy()
    clipped["geometry"] = clipped.geometry.intersection(boundary_simple)
    clipped = clipped[~clipped.geometry.is_empty]
    states["geometry"] = states.geometry.intersection(boundary_simple)
    outline = gpd.GeoSeries([boundary_simple], crs=config.METRIC_CRS)
    return clipped, states, outline


def finish_map(ax, outline, states, bounds=None):
    states.boundary.plot(ax=ax, color="#777777", linewidth=0.25, zorder=3)
    outline.boundary.plot(ax=ax, color="black", linewidth=0.6, zorder=4)
    if bounds is not None:
        ax.set_xlim(bounds[0], bounds[2])
        ax.set_ylim(bounds[1], bounds[3])
    ax.set_aspect("equal")
    ax.set_axis_off()


def categorical_map(ax, frame, column, colors, labels, order, title, legend_title=None, legend_loc="upper left"):
    frame.plot(ax=ax, color="#F7F7F7", edgecolor="none")
    for category in order:
        part = frame[frame[column] == category]
        if not part.empty:
            part.plot(ax=ax, color=colors[category], edgecolor="white", linewidth=0.05)
    handles = [Patch(facecolor=colors[c], edgecolor="#777777", linewidth=0.4, label=labels[c]) for c in order]
    ax.legend(handles=handles, title=legend_title, loc=legend_loc, frameon=False, fontsize=7, title_fontsize=7.5)
    ax.set_title(title, fontsize=9.5)


# ------------------------------------------------------------------ figuras

def fig_gi_maps(frame, states, outline):
    data = pd.read_csv(config.HOTSPOT_CSV_PATH, dtype={"CD_MUN": "string"})
    merged = frame.merge(data[["CD_MUN", "gi_bin", "gi_bin_pct"]], on="CD_MUN", how="inner")
    bounds = outline.total_bounds
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4))
    order = list(GI_COLORS)
    for ax, column, title, letter in (
        (axes[0], "gi_bin", "Área desmatada (km²)", "a"),
        (axes[1], "gi_bin_pct", "Proporção da área municipal desmatada (%)", "b"),
    ):
        categorical_map(ax, merged, column, GI_COLORS, GI_LABELS, order, title, legend_title="Gi* (confiança)")
        finish_map(ax, outline, states, bounds)
        panel_label(ax, letter)
    fig.suptitle("Hotspots de desmatamento (Getis-Ord Gi*) por município no Cerrado, 2020–2025", y=0.96, fontsize=11)
    fig.tight_layout()
    save(fig, "fig01_hotspots_gi_star.png")


def fig_annual_and_persistence(frame, states, outline):
    yearly = pd.read_csv(config.YEARLY_MUNICIPAL_AREA_PATH)
    by_year = yearly.groupby("year")["area_km2"].sum()
    periods = pd.read_csv(config.REFINED_DIR / "robustez_gi_periodos.csv")
    annual_hot = periods[periods["periodo"].isin([str(y) for y in by_year.index])]
    stability = pd.read_csv(
        config.REFINED_DIR / "robustez_gi_estabilidade_municipios.csv", dtype={"CD_MUN": "string"}
    )
    shift = pd.read_csv(config.REFINED_DIR / "robustez_gi_centro_gravidade.csv")
    status_colors = {
        "Persistente": "#B2182B",
        "Emergente (só 2023-2025)": "#2166AC",
        "Desaparece (só 2020-2022)": "#F4A582",
    }
    status_labels = {
        "Persistente": "Persistente (ambos os períodos)",
        "Emergente (só 2023-2025)": "Emergente (apenas 2023–2025)",
        "Desaparece (só 2020-2022)": "Extinto (apenas 2020–2022)",
    }
    merged = frame.merge(stability[["CD_MUN", "status_area"]], on="CD_MUN", how="left")

    fig = plt.figure(figsize=(12, 6.2))
    grid = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.2], wspace=0.12)
    ax0 = fig.add_subplot(grid[0])
    ax0.bar(by_year.index.astype(str), by_year.to_numpy(), color="#6D8F6A", width=0.65)
    for index, value in enumerate(by_year.to_numpy()):
        ax0.text(index, value + 600, pt_number(value), ha="center", fontsize=8)
    ax0.set_ylabel("Área desmatada (km²)")
    ax0.yaxis.set_major_formatter(pt_formatter())
    ax0.set_xlabel("Ano PRODES")
    ax0.set_ylim(0, by_year.max() * 1.12)
    twin = ax0.twinx()
    twin.spines["right"].set_visible(True)
    twin.plot(annual_hot["periodo"].astype(str), annual_hot["n_hotspots"], color=HOT_COLOR, marker="o", linewidth=1.5)
    twin.set_ylabel("Municípios hotspot por ano (n)", color=HOT_COLOR)
    twin.tick_params(axis="y", colors=HOT_COLOR)
    twin.set_ylim(0, annual_hot["n_hotspots"].max() * 1.25)
    panel_label(ax0, "a")

    ax1 = fig.add_subplot(grid[1])
    merged.plot(ax=ax1, color="#F7F7F7", edgecolor="none")
    for status, color in status_colors.items():
        merged[merged["status_area"] == status].plot(ax=ax1, color=color, edgecolor="white", linewidth=0.05)
    finish_map(ax1, outline, states, outline.total_bounds)
    points = gpd.GeoSeries(gpd.points_from_xy(shift["lon"], shift["lat"]), crs="EPSG:4326").to_crs(config.METRIC_CRS)
    ax1.annotate("", xy=(points.iloc[1].x, points.iloc[1].y), xytext=(points.iloc[0].x, points.iloc[0].y),
                 arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6), zorder=6)
    ax1.scatter([points.iloc[0].x], [points.iloc[0].y], s=28, color="white", edgecolor="black", zorder=6)
    ax1.scatter([points.iloc[1].x], [points.iloc[1].y], s=28, color="black", zorder=6)
    handles = [Patch(facecolor=status_colors[s], label=status_labels[s]) for s in status_colors]
    handles.append(Line2D([0], [0], color="black", marker="o", mfc="white", label="Centro de gravidade (2020–2022 → 2023–2025)"))
    ax1.legend(handles=handles, loc="upper left", frameon=False, fontsize=7.5)
    panel_label(ax1, "b")
    fig.suptitle("Dinâmica temporal da pressão de desmatamento", y=0.98, fontsize=11)
    save(fig, "fig02_dinamica_temporal.png")


def fig_robustness():
    table = pd.read_csv(config.REFINED_DIR / "robustez_gi_vizinhanca.csv")
    names = list(dict.fromkeys(table["vizinhanca"]))
    short = ["Queen\n(base)", "Queen\n2ª ordem", "8 vizinhos\nmais próximos", "Banda\n150 km"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), gridspec_kw={"width_ratios": [1.25, 1.25, 0.9]})
    positions = np.arange(len(names))
    width = 0.38
    for ax, variable, title, letter in (
        (axes[0], "area", "Área desmatada (km²)", "a"),
        (axes[1], "pct", "Proporção da área municipal (%)", "b"),
    ):
        part = table[table["variavel"] == variable].set_index("vizinhanca").loc[names]
        ax.bar(positions - width / 2, part["n_hotspots"], width, color=HOT_COLOR, label="Pseudo p-valor (permutação)")
        ax.bar(positions + width / 2, part["n_hotspots_fdr"], width, color="#F4A582", label="Após correção FDR (q = 0,05)")
        for x, (a, b) in enumerate(zip(part["n_hotspots"], part["n_hotspots_fdr"])):
            ax.text(x - width / 2, a + 6, str(a), ha="center", fontsize=7.5)
            ax.text(x + width / 2, b + 6, str(b), ha="center", fontsize=7.5)
        ax.set_xticks(positions)
        ax.set_xticklabels(short)
        ax.set_title(title, fontsize=9.5)
        ax.set_ylim(0, 560)
        panel_label(ax, letter)
    axes[0].set_ylabel("Municípios classificados como hot spot (n)")
    axes[0].legend(frameon=False, loc="upper left")
    jaccard = (
        table[table["vizinhanca"] != names[0]]
        .pivot(index="vizinhanca", columns="variavel", values="jaccard_vs_base")
        .loc[names[1:]]
    )
    ax = axes[2]
    ax.barh(np.arange(3) - 0.2, jaccard["area"], 0.38, color=CONTROL_COLOR, label="Área")
    ax.barh(np.arange(3) + 0.2, jaccard["pct"], 0.38, color="#9DB4C0", label="Proporção")
    ax.set_yticks(np.arange(3))
    ax.set_yticklabels(short[1:])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Índice de Jaccard em relação ao Queen")
    ax.xaxis.set_major_formatter(pt_formatter(1))
    ax.legend(frameon=False, loc="lower right")
    ax.invert_yaxis()
    panel_label(ax, "c")
    fig.suptitle("Sensibilidade dos hotspots à definição de vizinhança e à correção para testes múltiplos", y=1.02, fontsize=11)
    fig.tight_layout()
    save(fig, "fig03_robustez_gi_star.png")


def sankey_flows(trajectories):
    lookup = group_lookup()
    frame = pd.DataFrame(
        {
            "pre": lookup[trajectories["pre"].to_numpy()],
            "det": lookup[trajectories["det"].to_numpy()],
            "pos": lookup[trajectories["pos"].to_numpy()],
            "area": trajectories["area_w"].to_numpy(),
        }
    )
    return frame.groupby(["pre", "det"])["area"].sum(), frame.groupby(["det", "pos"])["area"].sum()


def draw_sankey(ax, first, second, title):
    total = float(first.sum())
    gap = 0.03 * total
    stage_tables = [first.groupby(level=0).sum(), first.groupby(level=1).sum(), second.groupby(level=1).sum()]
    x_positions = [0.0, 1.0, 2.0]
    bar_width = 0.14
    tops = []
    for stage, table in enumerate(stage_tables):
        cursor, stage_top = 0.0, {}
        for group_index in range(len(GROUP_ORDER)):
            value = float(table.get(group_index, 0.0))
            if value / total < 0.002:
                continue
            stage_top[group_index] = (cursor, cursor + value)
            ax.add_patch(plt.Rectangle((x_positions[stage] - bar_width / 2, cursor), bar_width, value,
                                       color=GROUP_COLORS[GROUP_ORDER[group_index]], zorder=3))
            if value / total > 0.035:
                ax.text(x_positions[stage] + bar_width / 2 + 0.03, cursor + value / 2,
                        pt_number(100 * value / total, 1) + "%", va="center", ha="left", fontsize=8)
            cursor += value + gap
        tops.append(stage_top)

    def ribbons(flows, left_stage):
        out_cursor = {g: tops[left_stage][g][0] for g in tops[left_stage]}
        in_cursor = {g: tops[left_stage + 1][g][0] for g in tops[left_stage + 1]}
        for (src, dst), value in flows.sort_index().items():
            if src not in out_cursor or dst not in in_cursor or value / total < 0.001:
                continue
            x0 = x_positions[left_stage] + bar_width / 2
            x1 = x_positions[left_stage + 1] - bar_width / 2
            y0, y1 = out_cursor[src], in_cursor[dst]
            mid = (x0 + x1) / 2
            verts = [(x0, y0), (mid, y0), (mid, y1), (x1, y1), (x1, y1 + value), (mid, y1 + value),
                     (mid, y0 + value), (x0, y0 + value), (x0, y0)]
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
                     Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
            ax.add_patch(PathPatch(Path(verts, codes), facecolor=GROUP_COLORS[GROUP_ORDER[src]],
                                   edgecolor="none", alpha=0.45, zorder=2))
            out_cursor[src] += value
            in_cursor[dst] += value

    ribbons(first, 0)
    ribbons(second, 1)
    for x, label in zip(x_positions, STAGE_LABELS):
        ax.text(x, -0.05 * total, label, ha="center", va="top", fontsize=9)
    ax.set_xlim(-0.2, 2.45)
    ax.set_ylim(-0.12 * total, max(t[max(t)][1] for t in tops) * 1.03)
    ax.axis("off")
    ax.set_title(title, fontsize=9.5)


def fig_sankey(trajectories):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    for ax, group, title, letter in (
        (axes[0], "Hotspot", "Municípios hotspot", "a"),
        (axes[1], "Controle", "Municípios de controle", "b"),
    ):
        part = trajectories[(trajectories["year"] <= 2024) & (trajectories["grupo"] == group)]
        draw_sankey(ax, *sankey_flows(part), title)
        panel_label(ax, letter)
    handles = [Patch(facecolor=GROUP_COLORS[g], label=CLASS_GROUP_NAMES[g]) for g in GROUP_ORDER[:4]]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Fluxos de cobertura dos pixels desmatados (coortes PRODES 2020–2024; área estimada)", y=1.0, fontsize=11)
    fig.tight_layout()
    save(fig, "fig04_fluxos_cobertura.png")


def fig_typology():
    table = pd.read_csv(config.REFINED_DIR / "tipologia_por_grupo.csv")
    levels = pd.read_csv(config.REFINED_DIR / "tipologia_por_nivel_confianca.csv")
    levels = levels[levels["categoria"] == TYPOLOGY_ORDER[0]].set_index("hot_nivel")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [1.8, 1]})
    ax = axes[0]
    positions = np.arange(len(TYPOLOGY_ORDER))
    width = 0.38
    for offset, group, color in ((-width / 2, "Hotspot", HOT_COLOR), (width / 2, "Controle", CONTROL_COLOR)):
        part = table[table["grupo"] == group].set_index("categoria").reindex(TYPOLOGY_ORDER)
        error = [part["participacao_pct"] - part["ic95_inf"], part["ic95_sup"] - part["participacao_pct"]]
        ax.bar(positions + offset, part["participacao_pct"], width, yerr=error, capsize=2.5,
               color=color, label=group, error_kw=dict(lw=0.9))
        for x, value, high in zip(positions + offset, part["participacao_pct"], part["ic95_sup"]):
            ax.text(x, high + 1.2, pt_number(value, 1), ha="center", fontsize=7)
    ax.set_xticks(positions)
    ax.set_xticklabels(TYPOLOGY_WRAPPED, fontsize=8)
    ax.set_ylabel("Participação na área desmatada (%)")
    ax.set_ylim(0, 62)
    ax.yaxis.set_major_formatter(pt_formatter())
    ax.legend(frameon=False, title="Grupo de municípios")
    panel_label(ax, "a")

    ax = axes[1]
    labels = ["Controle", "Hotspot\n90%", "Hotspot\n95%", "Hotspot\n99%"]
    colors = [CONTROL_COLOR, "#FDBB84", "#EF6548", "#B2182B"]
    levels_used = levels.loc[[0.0, 90.0, 95.0, 99.0]]
    values = levels_used["participacao_pct"].to_numpy()
    low, high = levels_used["ic95_inf"].to_numpy(), levels_used["ic95_sup"].to_numpy()
    ax.bar(np.arange(4), values, 0.6, color=colors, yerr=[values - low, high - values], capsize=3, error_kw=dict(lw=0.9))
    for x, value in enumerate(values):
        ax.text(x, high[x] + 1.2, pt_number(value, 1), ha="center", fontsize=8)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Conversão no ano da detecção (%)")
    ax.set_ylim(0, 62)
    panel_label(ax, "b")
    fig.suptitle("Tipologia de trajetória por pixel: hotspots × controle (IC 95% por bootstrap de municípios)", y=1.03, fontsize=11)
    fig.tight_layout()
    save(fig, "fig05_tipologia_hotspot_controle.png")


def fig_conversion():
    curves = pd.read_csv(config.REFINED_DIR / "tempo_ate_conversao.csv")
    states = pd.read_csv(config.REFINED_DIR / "conversao_por_estado.csv")
    fig = plt.figure(figsize=(13, 4.6))
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.1], wspace=0.28)
    cmap = plt.get_cmap("viridis")
    years = sorted(curves["ano_prodes"].unique())
    axes = []
    for index, (group, letter) in enumerate((("Hotspot", "a"), ("Controle", "b"))):
        ax = fig.add_subplot(grid[index], sharey=axes[0] if axes else None)
        axes.append(ax)
        for k, year in enumerate(years):
            part = curves[(curves["ano_prodes"] == year) & (curves["grupo"] == group)]
            if part.empty:
                continue
            ax.plot(part["defasagem_anos"], part["taxa_conversao_acumulada_pct"], marker="o", ms=3.5,
                    color=cmap(k / max(len(years) - 1, 1)), label=str(year))
        ax.set_title(f"Municípios {'hotspot' if group == 'Hotspot' else 'de controle'}", fontsize=9.5)
        ax.set_xlabel("Anos após o ano de detecção")
        ax.set_xticks(range(6))
        ax.grid(alpha=0.3, linewidth=0.5)
        panel_label(ax, letter)
    axes[0].set_ylabel("Conversão acumulada da área natural em t−1 (%)")
    axes[1].legend(title="Coorte PRODES", frameon=False, fontsize=7.5, loc="center right")
    ax = fig.add_subplot(grid[2])
    order = states[states["grupo"] == "Hotspot"].sort_values("area_convertida_km2", ascending=False)["SIGLA_UF"].tolist()
    height = 0.38
    ypos = np.arange(len(order))
    for offset, group, color in ((-height / 2, "Hotspot", HOT_COLOR), (height / 2, "Controle", CONTROL_COLOR)):
        part = states[states["grupo"] == group].set_index("SIGLA_UF").reindex(order)
        ax.barh(ypos + offset, part["taxa_conversao_pct"], height, color=color, label=group)
    ax.set_yticks(ypos)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Taxa de conversão até 2025 (%)")
    ax.legend(frameon=False, loc="lower right")
    panel_label(ax, "c")
    fig.suptitle("Persistência e velocidade da conversão da vegetação natural", y=1.03, fontsize=11)
    save(fig, "fig06_conversao_acumulada.png")


def fig_conversion_maps(frame, states, outline):
    destinations = pd.read_csv(config.REFINED_DIR / "destinos_conversao_por_municipio.csv", dtype={"CD_MUN": "string"})
    destinations = destinations[destinations["estimativa_confiavel"]]
    merged = frame.merge(destinations, on="CD_MUN", how="inner")
    hot = merged[merged["grupo"] == "Hotspot"]
    bounds = outline.total_bounds
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.4))
    ax = axes[0]
    frame.plot(ax=ax, color="#F7F7F7", edgecolor="none")
    hot.plot(ax=ax, column="taxa_conversao_pct", cmap="YlOrRd", vmin=0, vmax=100, edgecolor="white",
             linewidth=0.05, legend=True, legend_kwds={"shrink": 0.55, "label": "Taxa de conversão (%)"})
    finish_map(ax, outline, states, bounds)
    ax.set_title("Taxa de conversão da vegetação natural (municípios hotspot)", fontsize=9.5)
    panel_label(ax, "a")
    ax = axes[1]
    colors = {"Mosaico de Agricultura e Pastagem": "#8C6BB1", "Pastagem": "#41AB5D"}
    categorical_map(ax, hot, "destino_predominante", colors, {k: k for k in colors}, list(colors),
                    "Classe de destino predominante", legend_title="Uso após a conversão")
    finish_map(ax, outline, states, bounds)
    panel_label(ax, "b")
    fig.suptitle("Conversão e uso da terra pós-desmatamento nos municípios hotspot", y=0.96, fontsize=11)
    fig.tight_layout()
    save(fig, "fig07_conversao_destino.png")


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    style()
    config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for old in config.FIGURES_DIR.glob("*.png"):
        old.unlink()
    trajectories = pd.read_parquet(config.REFINED_DIR / "_trajetorias_pixel.parquet")
    frame, states, outline = load_map_layers()
    frame["CD_MUN"] = frame["CD_MUN"].astype("string")
    fig_gi_maps(frame, states, outline)
    fig_annual_and_persistence(frame, states, outline)
    fig_robustness()
    fig_sankey(trajectories)
    fig_typology()
    fig_conversion()
    fig_conversion_maps(frame, states, outline)


if __name__ == "__main__":
    main()
