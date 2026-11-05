"""
Legenda oficial das classes MapBiomas Brasil, Coleção 4.

Fontes oficiais:
- https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/08/Codigos_da_legenda_Colecao_4.pdf
- https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/08/Descricao-detalhada-de-legenda-Cobertura-10m-Col4.pdf
"""

MAPBIOMAS_CLASS_NAMES = {
    1: "Floresta",
    3: "Formação Florestal",
    4: "Formação Savânica",
    5: "Mangue",
    6: "Floresta Alagável",
    9: "Floresta Plantada",
    10: "Formação Natural não Florestal",
    11: "Área Úmida Natural não Florestal",
    12: "Formação Campestre",
    13: "Outra Formação Natural não Florestal",
    14: "Agropecuária",
    15: "Pastagem",
    18: "Agricultura",
    19: "Cultura Anual e Perene",
    20: "Cultura Semi-Perene",
    21: "Mosaico de Agricultura e Pastagem",
    22: "Área não Vegetada",
    23: "Praia e Duna",
    24: "Infraestrutura Urbana",
    25: "Outra Área não Vegetada",
    26: "Corpos d'Água",
    27: "Não observado",
    29: "Afloramento Rochoso",
    30: "Mineração",
    31: "Aquicultura",
    32: "Apicum",
    33: "Rio, Lago e Oceano",
    36: "Lavoura Perene",
    48: "Outras Lavouras Perenes",
    49: "Restinga Arbórea",
    50: "Restinga Herbácea ou Arbustiva",
    75: "Usina Fotovoltaica",
    91: "Não documentada na legenda",
}


def format_class_name(class_code):
    """Retorna o código e o nome oficial de uma classe MapBiomas."""
    if class_code is None:
        return "Sem dados"

    try:
        if class_code != class_code:
            return "Sem dados"
    except TypeError:
        return "Sem dados"

    code = int(class_code)
    if code == 0:
        return "Sem dados"
    try:
        name = MAPBIOMAS_CLASS_NAMES[code]
    except KeyError as exc:
        raise ValueError(
            f"Código MapBiomas {code} ausente na legenda oficial da Coleção 4"
        ) from exc
    return f"{code} - {name}"
