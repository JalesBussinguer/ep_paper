"""
Configurações gerais do pipeline de análise de hotspots de desmatamento.

Todos os caminhos e parâmetros usados pelos outros módulos ficam
centralizados aqui, para facilitar ajustes sem precisar mexer na lógica da
análise.
"""

from pathlib import Path

# Raiz do projeto (dois níveis acima deste arquivo: src/config.py -> ep_paper/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"
HOTSPOT_GEOJSON_PATH = OUTPUT_DIR / "municipios_hotspots_cerrado.geojson"
HOTSPOT_CSV_PATH = OUTPUT_DIR / "municipios_hotspots_cerrado.csv"
STATE_SUMMARY_PATH = OUTPUT_DIR / "resumo_hotspots_por_estado_cerrado.csv"
AREA_MAP_PATH = OUTPUT_DIR / "mapa_hotspots_area_desmatada_cerrado.png"
PERCENT_MAP_PATH = OUTPUT_DIR / "mapa_hotspots_pct_desmatado_cerrado.png"
MAPBIOMAS_DIR = Path(r"H:\mapbiomas_maps")
MAPBIOMAS_YEARS = tuple(range(2019, 2026))
MAPBIOMAS_RASTER_PATTERN = "brazil_coverage-col4_10m_{year}.tif"
MAPBIOMAS_TRANSITIONS_CSV_PATH = (
    OUTPUT_DIR / "poligonos_hotspots_transicoes_mapbiomas.csv"
)
MAPBIOMAS_TRANSITION_SUMMARY_PATH = (
    OUTPUT_DIR / "resumo_transicoes_mapbiomas_por_ano.csv"
)
MAPBIOMAS_CLASS_COMPOSITION_CSV_PATH = (
    OUTPUT_DIR / "composicao_classes_poligonos_mapbiomas.csv.gz"
)
MAPBIOMAS_PIXEL_TRANSITIONS_CSV_PATH = (
    OUTPUT_DIR / "transicoes_classes_poligonos_mapbiomas.csv.gz"
)
MAPBIOMAS_CLASS_COMPOSITION_SUMMARY_PATH = (
    OUTPUT_DIR / "resumo_composicao_classes_mapbiomas.csv"
)
MAPBIOMAS_PIXEL_TRANSITIONS_SUMMARY_PATH = (
    OUTPUT_DIR / "resumo_transicoes_classes_mapbiomas.csv"
)
MAPBIOMAS_EVENT_TRANSITION_SUMMARY_PATH = (
    OUTPUT_DIR / "resumo_transicoes_evento_mapbiomas.csv"
)
MAPBIOMAS_EVENT_CLASS_TRANSITIONS_PATH = (
    OUTPUT_DIR / "transicoes_classes_evento_mapbiomas.csv"
)
MAPBIOMAS_EVENT_GROUP_TRANSITIONS_PATH = (
    OUTPUT_DIR / "transicoes_grupos_evento_mapbiomas.csv"
)
MAPBIOMAS_EVENT_MUNICIPAL_TRANSITIONS_PATH = (
    OUTPUT_DIR / "transicoes_municipais_evento_mapbiomas.csv.gz"
)
MAPBIOMAS_EVENT_GROUP_TRAJECTORIES_PATH = (
    OUTPUT_DIR / "trajetorias_grupos_evento_mapbiomas.csv"
)

# Etapa de refinamento metodológico (amostragem por pixel, controle, robustez)
YEARLY_MUNICIPAL_AREA_PATH = OUTPUT_DIR / "area_desmatada_municipio_ano.csv"
SAMPLE_FRAGMENTS_PATH = OUTPUT_DIR / "amostra_fragmentos_prodes.gpkg"
SAMPLE_PIXEL_COUNTS_PATH = OUTPUT_DIR / "amostra_contagem_pixels.csv.gz"
SAMPLE_TABLE_PATH = OUTPUT_DIR / "amostra_fragmentos_prodes.csv"
REFINED_DIR = OUTPUT_DIR / "refinamento"
FIGURES_DIR = REFINED_DIR / "figuras"
VALIDATION_DIR = REFINED_DIR / "validacao_visual"

# Amostragem de Poisson por classe de tamanho (km2): limite superior -> probabilidade
# de inclusão. O peso de cada fragmento é 1 / probabilidade (estimador de
# Horvitz-Thompson), o que permite estimar totais de área da população.
SAMPLE_SIZE_BREAKS_KM2 = (0.005, 0.02, 0.1, 0.5, 2.0, 10.0, float("inf"))
SAMPLE_INCLUSION_PROBABILITIES = (0.005, 0.01, 0.03, 0.15, 0.5, 1.0, 1.0)
SAMPLE_SEED = 20260101

DEFORESTATION_PATH = DATA_DIR / "cerrado_desmatamento_2020_2026.geojson"
MUNICIPALITIES_PATH = DATA_DIR / "municipios.geojson"
CERRADO_BOUNDARY_PATH = DATA_DIR / "cerrado_border.geojson"

# CRS geográfico original dos dados de entrada (SIRGAS 2000)
INPUT_CRS = "EPSG:4674"

# CRS projetado (metros) usado para cálculos métricos (centroide, área,
# distância/contiguidade). SIRGAS 2000 / Brazil Polyconic é adequado para
# análises que cobrem todo o território brasileiro.
METRIC_CRS = "EPSG:5880"

# Identificador necessário para contar cada polígono por município.
# A área é calculada a partir das interseções geométricas no CRS métrico.
DEFORESTATION_COLUMNS = ["fid", "year"]

# Colunas do GeoJSON de municípios que serão mantidas para a análise.
MUNICIPALITY_COLUMNS = [
    "CD_MUN",
    "NM_MUN",
    "CD_UF",
    "NM_UF",
    "SIGLA_UF",
    "AREA_KM2",
]

CERRADO_BOUNDARY_COLUMNS = ["bioma", "cd_bioma"]

# Parâmetros da análise de hotspots (Getis-Ord Gi*)
GI_STAR_PERMUTATIONS = 999
GI_STAR_SEED = 42

# Níveis de confiança (bilateral) usados para classificar os municípios em
# hot spots / cold spots, seguindo a mesma lógica do "Gi_Bin" do ArcGIS (ver
# seção "How Hot Spot Analysis works" da documentação da Esri).
CONFIDENCE_Z_THRESHOLDS = {
    0.90: 1.645,
    0.95: 1.960,
    0.99: 2.576,
}
