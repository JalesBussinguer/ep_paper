"""
Geração de mapa (PNG) com o resultado da análise de hotspots por município.

Existem duas implementações:
  - `plot_hotspot_map`: usa matplotlib/geopandas (abordagem padrão,
    recomendada quando o ambiente Python tem uma instalação saudável do
    matplotlib).
  - `render_hotspot_map_pil`: um rasterizador simples baseado em Pillow, sem
    depender do matplotlib. Foi adicionado como alternativa robusta, pois em
    alguns ambientes o backend nativo do matplotlib (Agg) pode falhar de
    forma abrupta (crash do processo) por problemas de instalação alheios a
    este projeto.
"""

import logging
import math
import unicodedata

from PIL import Image, ImageDraw, ImageFont

logger = logging.getLogger(__name__)

# Cores seguindo a paleta clássica do "Gi_Bin" do ArcGIS (vermelho para hot
# spots, azul para cold spots, cinza para não significante).
GI_BIN_COLORS = {
    "Hot Spot 99%": "#a50026",
    "Hot Spot 95%": "#e31a1c",
    "Hot Spot 90%": "#fb9a99",
    "Não significante": "#d9d9d9",
    "Cold Spot 90%": "#a6bddb",
    "Cold Spot 95%": "#3690c0",
    "Cold Spot 99%": "#023858",
}

# Ordem de plotagem/legenda (do cold spot mais forte ao hot spot mais forte)
GI_BIN_ORDER = [
    "Cold Spot 99%",
    "Cold Spot 95%",
    "Cold Spot 90%",
    "Não significante",
    "Hot Spot 90%",
    "Hot Spot 95%",
    "Hot Spot 99%",
]


def plot_hotspot_map(municipalities, bin_col, title, output_path):
    """
    Plota um mapa coroplético dos municípios coloridos pela classificação de
    hotspot (`bin_col`) e salva em `output_path` (PNG).
    """
    # Import feito aqui (e não no topo do módulo) para que
    # `render_hotspot_map_pil` continue funcionando mesmo se a instalação do
    # matplotlib estiver com problemas no ambiente.
    import matplotlib

    matplotlib.use("Agg")  # backend não interativo, seguro para rodar via script
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 10))
    bounds = municipalities.total_bounds

    for categoria in GI_BIN_ORDER:
        subset = municipalities[municipalities[bin_col] == categoria]
        if subset.empty:
            continue
        subset.plot(
            ax=ax,
            color=GI_BIN_COLORS[categoria],
            edgecolor="white",
            linewidth=0.1,
            label=categoria,
        )

    ax.set_title(title, fontsize=13)
    _add_matplotlib_map_elements(ax, bounds)
    ax.set_axis_off()
    ax.legend(loc="lower left", fontsize=8, frameon=True, title="Classificação (Gi*)")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    logger.info("Mapa salvo em %s", output_path)


def map_scale_length_km(map_width_m):
    raw_length_km = map_width_m * 0.18 / 1000
    if raw_length_km <= 0:
        raise ValueError("A largura do mapa deve ser positiva.")
    magnitude = 10 ** math.floor(math.log10(raw_length_km))
    return max(step * magnitude for step in (1, 2, 5) if step * magnitude <= raw_length_km)


def _add_matplotlib_map_elements(ax, bounds):
    width, height = bounds[2] - bounds[0], bounds[3] - bounds[1]
    scale_length_km = map_scale_length_km(width)
    padding = 0.16
    ax.set_xlim(bounds[0] - width * padding, bounds[2] + width * padding)
    ax.set_ylim(bounds[1] - height * padding, bounds[3] + height * padding)
    x = bounds[0] - width * padding * 0.75
    y = 0.045
    length = scale_length_km * 1000
    transform = ax.get_xaxis_transform()
    ax.plot([x, x + length], [y, y], color="black", linewidth=3,
            transform=transform, clip_on=False, zorder=7)
    ax.plot([x, x + length], [y, y], color="white", linewidth=1.2,
            transform=transform, clip_on=False, zorder=8)
    for tick_x in (x, x + length):
        ax.plot([tick_x, tick_x], [y - 0.012, y + 0.012], color="black",
                linewidth=1, transform=transform, clip_on=False, zorder=9)
    ax.text(
        x + length / 2,
        y + 0.015,
        f"{scale_length_km:g} km",
        transform=transform,
        ha="center",
        va="bottom",
        fontsize=8,
        clip_on=False,
        zorder=9,
    )
    ax.annotate(
        "N",
        xy=(0.95, 0.94),
        xytext=(0.95, 0.86),
        xycoords="axes fraction",
        textcoords="axes fraction",
        ha="center",
        va="bottom",
        fontweight="bold",
        annotation_clip=False,
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 1.2},
        zorder=10,
    )
    ax.text(
        0.99,
        0.015,
        "SIRGAS 2000 / Brazil Polyconic (EPSG:5880)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6,
        clip_on=False,
        zorder=10,
    )


def _remove_accents(text):
    """
    Remove acentos de um texto (ex.: "Não significante" -> "Nao significante").

    Usado apenas no texto desenhado sobre o mapa (titulo/legenda): a fonte do
    Pillow usada em `render_hotspot_map_pil` depende do freetype, cuja
    instalacao neste ambiente demonstrou corromper glifos acentuados. Isso
    nao afeta os comentarios do codigo nem os dados de saida (CSV/GeoJSON),
    apenas o texto rasterizado no PNG.
    """
    normalizado = unicodedata.normalize("NFKD", text)
    return "".join(caractere for caractere in normalizado if not unicodedata.combining(caractere))


def _hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))


def _polygon_to_pixels(polygon, minx, maxy, scale, margin):
    """Converte os anéis (exterior/interiores) de um polígono shapely para
    coordenadas de pixel (origem no canto superior esquerdo da imagem)."""
    exterior = [
        (margin + (x - minx) * scale, margin + (maxy - y) * scale)
        for x, y in polygon.exterior.coords
    ]
    interiors = [
        [
            (margin + (x - minx) * scale, margin + (maxy - y) * scale)
            for x, y in ring.coords
        ]
        for ring in polygon.interiors
    ]
    return exterior, interiors


def render_hotspot_map_pil(
    municipalities,
    bin_col,
    title,
    output_path,
    width=1600,
    background_color=(255, 255, 255),
):
    """
    Rasteriza um mapa coroplético simples (sem depender do matplotlib),
    desenhando cada município (Polygon/MultiPolygon) com Pillow, colorido de
    acordo com a classificação de hotspot (`bin_col`).

    Os furos (interiores) dos polígonos são "recortados" preenchendo-os com
    a cor de fundo por cima do preenchimento externo.
    """
    margin = 40
    horizontal_padding = 0.12
    map_top = margin + 30
    legend_h = 30 * (len(GI_BIN_ORDER) + 1) + 20

    minx, miny, maxx, maxy = municipalities.total_bounds
    map_width_px = int((width - 2 * margin) * (1 - 2 * horizontal_padding))
    scale = map_width_px / (maxx - minx)
    map_left = margin + int((width - 2 * margin) * horizontal_padding)
    map_right = map_left + map_width_px
    map_height_px = int((maxy - miny) * scale)
    map_bottom = map_top + map_height_px
    height = map_bottom + 80 + legend_h

    img = Image.new("RGB", (width, height), background_color)
    draw = ImageDraw.Draw(img)

    for categoria in GI_BIN_ORDER:
        subset = municipalities[municipalities[bin_col] == categoria]
        if subset.empty:
            continue
        cor = _hex_to_rgb(GI_BIN_COLORS[categoria])
        for geom in subset.geometry:
            if geom is None or geom.is_empty:
                continue
            polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
            for poly in polys:
                exterior, interiors = _polygon_to_pixels(poly, minx, maxy, scale, map_left)
                exterior = [(x, y + map_top - margin) for x, y in exterior]
                interiors = [
                    [(x, y + map_top - margin) for x, y in ring]
                    for ring in interiors
                ]
                draw.polygon(exterior, fill=cor, outline=(255, 255, 255))
                for ring in interiors:
                    draw.polygon(ring, fill=background_color)

    # Usa a fonte bitmap embutida do Pillow (não depende do freetype), pois
    # a instalação de freetype disponível neste ambiente demonstrou
    # instabilidade (glifos corrompidos) ao renderizar fontes TrueType.
    font = ImageFont.load_default()
    font_title = font

    draw.text((margin, 5), _remove_accents(title), fill=(0, 0, 0), font=font_title)

    scale_length_km = map_scale_length_km(maxx - minx)
    scale_length_px = scale_length_km * 1000 * scale
    scale_x = margin
    scale_y = map_bottom + 20
    draw.line(
        (scale_x, scale_y, scale_x + scale_length_px, scale_y),
        fill=(0, 0, 0),
        width=6,
    )
    draw.line(
        (scale_x, scale_y, scale_x + scale_length_px, scale_y),
        fill=(255, 255, 255),
        width=2,
    )
    for tick_x in (scale_x, scale_x + scale_length_px):
        draw.line((tick_x, scale_y - 8, tick_x, scale_y + 8), fill=(0, 0, 0), width=2)
    scale_label = f"{scale_length_km:g} km"
    label_box = draw.textbbox((0, 0), scale_label, font=font)
    draw.rectangle(
        (
            scale_x + scale_length_px / 2 - (label_box[2] - label_box[0]) / 2 - 3,
            scale_y - 25,
            scale_x + scale_length_px / 2 + (label_box[2] - label_box[0]) / 2 + 3,
            scale_y - 8,
        ),
        fill=background_color,
    )
    draw.text(
        (scale_x + scale_length_px / 2, scale_y - 24),
        scale_label,
        fill=(0, 0, 0),
        font=font,
        anchor="mt",
    )

    north_x = width - margin / 2
    north_y = map_top + 20
    draw.text((north_x, north_y - 45), "N", fill=(0, 0, 0), font=font, anchor="mm")
    draw.line((north_x, north_y - 30, north_x, north_y + 5), fill=(0, 0, 0), width=3)
    draw.polygon(
        [(north_x, north_y - 37), (north_x - 8, north_y - 22), (north_x + 8, north_y - 22)],
        fill=(0, 0, 0),
    )

    legend_y = height - legend_h + 10
    draw.text(
        (margin, legend_y - 38),
        _remove_accents("SIRGAS 2000 / Brazil Polyconic (EPSG:5880)"),
        fill=(0, 0, 0),
        font=font,
    )
    draw.text((margin, legend_y - 20), _remove_accents("Classificação (Gi*):"), fill=(0, 0, 0), font=font)
    for i, categoria in enumerate(reversed(GI_BIN_ORDER)):
        y = legend_y + i * 26
        cor = _hex_to_rgb(GI_BIN_COLORS[categoria])
        draw.rectangle([margin, y, margin + 20, y + 18], fill=cor, outline=(0, 0, 0))
        draw.text((margin + 30, y), _remove_accents(categoria), fill=(0, 0, 0), font=font)

    img.save(output_path)
    logger.info("Mapa (Pillow) salvo em %s", output_path)
