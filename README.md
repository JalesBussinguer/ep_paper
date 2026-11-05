# Análise de Hotspots de Desmatamento no Cerrado (PRODES 2020-2026)

Análise de hotspots (Getis-Ord Gi*) da pressão de desmatamento no bioma
Cerrado, com base nos polígonos de desmatamento do PRODES (2020-2026) e nos
limites municipais do IBGE.

## Metodologia

Segue a metodologia de Hot Spot Analysis (Getis-Ord Gi*) descrita pela Esri:
https://doc.esri.com/en/arcgis-pro/latest/tool-reference/spatial-statistics/h-how-hot-spot-analysis-getis-ord-gi-spatial-stati.html

1. **Objeto de análise**: cada polígono de desmatamento do PRODES
   (`data/cerrado_desmatamento_2020_2026.geojson`). A geometria completa é
   usada para calcular a área de interseção; o atributo pré-calculado
   `area_km` não é usado nesta etapa.
2. **Recorte espacial**: os limites municipais do IBGE
   (`data/municipios.geojson`) são recortados pela geometria do Cerrado
   (`data/cerrado_border.geojson`). A análise inclui todos os municípios
   com área positiva dentro do bioma, mesmo os que não registraram
   desmatamento. Centroides fora do limite do Cerrado não entram na análise;
   os que caem fora de um limite municipal recortado são associados ao
   município mais próximo.
3. **Agregação municipal**: cada polígono PRODES é intersectado com os
   municípios recortados pelo Cerrado. Se um polígono atravessa fronteiras,
   sua área é repartida conforme a área geométrica de cada interseção; as
   porções fora do Cerrado não são contabilizadas. As áreas das geometrias
   são calculadas em projeção métrica (EPSG:5880), e geometrias inválidas
   são reparadas antes das operações espaciais.
4. **Pesos espaciais**: contiguidade Queen (borda ou vértice compartilhado)
   entre os municípios dessa malha.
5. **Getis-Ord Gi\***: calculado sobre duas variáveis, para cada município:
   - `area_desmatada_km2`: soma da área desmatada no período (variável
     pedida originalmente).
   - `pct_area_desmatada`: % da área municipal **dentro do Cerrado**
     desmatada no período. O denominador é a área geométrica do município
     após o recorte pelo limite do bioma (`AREA_CERRADO_KM2`), não a área
     total municipal (`AREA_KM2`).
6. Cada município é classificado em hot spot / cold spot (90/95/99% de
   confiança, com base no Z score e no pseudo p-valor por permutação) ou
   "não significante", e catalogado por estado (UF) e município.

## Estrutura do projeto

```
data/                                  dados de entrada (não alterados)
  cerrado_border.geojson              limite usado para recortar municípios
src/
  config.py                            caminhos e parâmetros da análise
  io_utils.py                          leitura dos GeoJSON de entrada
  aggregation.py                       join espacial + agregação por município
  hotspots.py                          pesos espaciais + Getis-Ord Gi*
  visualization.py                     geração de mapas (matplotlib e Pillow)
  make_maps.py                         script para gerar só os mapas a partir do resultado já salvo
run_hotspot_analysis.py                script principal (roda tudo, ponta a ponta)
output/                                resultados gerados pela análise
```

## Como rodar

```
pip install -r requirements.txt
python run_hotspot_analysis.py
```

## Transições de classe MapBiomas nos municípios hotspot

Para extrair a classe anual dentro de cada polígono PRODES pertencente a um
município hotspot (hot spot em `gi_bin` **ou** `gi_bin_pct`, em qualquer
nível de confiança), execute:

```
python -m src.mapbiomas_transitions
```

O script espera os GeoTIFFs da Coleção 4, 10 m, em
`H:\mapbiomas_maps\brazil_coverage-col4_10m_2019.tif` até o ano de 2025.
As sete grades são verificadas quanto a CRS, resolução, dimensões e
alinhamento antes da extração.

`poligonos_hotspots_transicoes_mapbiomas.csv` contém uma linha para cada par
polígono PRODES–município hotspot intersectado. Inclui `fid`, ano PRODES,
município/estado, área da interseção e `classe_2019` a `classe_2025`. Cada
classe é a classe dominante por área de pixel dentro da porção do polígono
que cai naquele município hotspot, identificada pelo código e nome em
português conforme os [códigos oficiais da Coleção 4](https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/08/Codigos_da_legenda_Colecao_4.pdf)
e a [descrição oficial do produto Cobertura 10 m](https://brasil.mapbiomas.org/wp-content/uploads/sites/3/2026/08/Descricao-detalhada-de-legenda-Cobertura-10m-Col4.pdf).
Os campos extras `classe_pre_evento`, `classe_ano_deteccao` e
`classe_primeiro_ano_posterior` facilitam a comparação temporal. Para
eventos de 2025 e 2026 não existe mapa do primeiro ano posterior disponível
na série (que termina em 2025); isso é indicado em `status_pos_deteccao`.

`resumo_transicoes_mapbiomas_por_ano.csv` contabiliza os fragmentos por
estado, município, ano PRODES, ano MapBiomas e classe dominante. A extração
é feita por blocos para reduzir o uso de memória. Um polígono que cruza
fronteiras pode gerar uma linha por município hotspot intersectado; apenas
a área da interseção com esse município é classificada.

## Análise das transições associadas ao evento PRODES

Para comparar as classes MapBiomas antes, durante e depois do ano de detecção
PRODES, execute:

```
python -m src.mapbiomas_event_transition_analysis
```

O módulo lê `poligonos_hotspots_transicoes_mapbiomas.csv` e analisa as
colunas `classe_pre_evento`, `classe_ano_deteccao` e
`classe_primeiro_ano_posterior`. A unidade é a interseção entre um polígono
PRODES e um município hotspot. São avaliados três intervalos: pré-evento →
detecção, detecção → pós-evento e pré-evento → pós-evento. Além das transições
exatas entre classes da legenda, os códigos são agrupados em vegetação
natural, uso antrópico, superfície não vegetada, água e sem dados.

O processo gera:

- `resumo_transicoes_evento_mapbiomas.csv`: quantidade e percentual de
  interseções com mudança de classe dominante ou de grupo, e as áreas das
  interseções correspondentes, no total e por coorte anual PRODES.
- `transicoes_classes_evento_mapbiomas.csv`: matriz detalhada de transições
  entre códigos/classes, no total e por ano PRODES.
- `transicoes_grupos_evento_mapbiomas.csv`: transições entre grupos de
  classes, no total e por ano PRODES.
- `transicoes_municipais_evento_mapbiomas.csv.gz`: transições entre grupos
  por município, estado, intervalo e ano PRODES.
- `trajetorias_grupos_evento_mapbiomas.csv`: trajetórias resumidas de grupo
  pré-evento → ano de detecção → pós-evento, no total e por coorte.

Os percentuais de interseções usam como denominador todas as interseções
válidas para aquele intervalo e coorte. Os percentuais de área usam a soma da
área geométrica dessas interseções válidas. A área de interseções cuja classe
dominante mudou é apenas uma medida ponderada pela área do fragmento; ela não
mede a área de pixels que efetivamente mudou de classe. Os polígonos podem
gerar múltiplas interseções municipais e, portanto, não devem ser interpretados
como eventos PRODES únicos. A classe dominante é uma síntese de cada
interseção, não a composição completa dos pixels. Eventos sem classe
posterior disponível (por exemplo, detecções de 2025 na série até 2025) são
excluídos apenas dos intervalos/trajetórias que precisam dessa informação.
As mudanças observadas não demonstram, por si só, causalidade do evento
PRODES.

O dataset de desmatamento tem quase 2 milhões de polígonos (~3 GB), então a
leitura é uma das etapas mais demoradas. O recorte espacial e a preparação
dos municípios também podem levar cerca de um minuto; a análise Gi* em si
trabalha ao nível de município (poucos milhares de unidades).

## Refinamento metodológico (composição por pixel, controle e robustez)

Execução (nesta ordem; `output/refinamento/` recebe as saídas):

```powershell
python -m src.yearly_aggregation   # área município×ano + amostra de Poisson (~16 min)
python -m src.pixel_extraction     # trajetórias 2019-2025 por pixel (~50 min, 12 workers)
python -m src.pixel_analysis       # tipologia, hotspot × controle, matrizes, tempo até conversão
python -m src.gi_robustness        # Moran, vizinhanças, FDR, períodos
python -m src.figures              # Sankey, barras, curvas e mapas
python -m src.visual_validation    # concordância PRODES × MapBiomas e pranchas
```

- **Amostragem**: fragmentos (polígono ∩ município ∩ Cerrado) sorteados por
  Poisson com probabilidade dependente da classe de tamanho; estimativas de
  área por Horvitz-Thompson (`peso = 1 / p_incl`). Hotspots e controle usam as
  mesmas probabilidades.
- **Pixels**: regra do centro do pixel; cada pixel guarda a classe MapBiomas de
  2019 a 2025. Tipologias: conversão no ano da detecção, conversão posterior,
  vegetação persistente, já antrópico antes, regeneração e outras (2025 fica
  sem pós-evento e fora das comparações principais).
- **Hotspot × controle**: hotspot = "Hot Spot" em qualquer variável do Gi*;
  controle = demais municípios. IC por bootstrap de municípios (500), V de
  Cramér e p por permutação de rótulos entre municípios (2.000).
- **Robustez do Gi***: Moran's I, Queen (base) × Queen de 2ª ordem × 8 vizinhos
  × banda de 150 km (Jaccard e kappa), FDR de Benjamini-Hochberg, períodos
  2020-2022 × 2023-2025 e por ano.
- **Limitações**: a unidade é a interseção (não o evento); a resolução é anual;
  não há inferência causal; o código 91 do raster (≈0,002% dos pixels) não
  consta da legenda Col4 e é tratado como "sem dados"; a conferência em imagem
  de satélite é manual (`validacao_visual/amostra_validacao_visual.csv`).

## Saídas (pasta `output/`)

- `municipios_hotspots_cerrado.geojson` / `.csv`: um registro por município
  que intersecta o Cerrado, com:
  `CD_MUN`, `NM_MUN`, `SIGLA_UF`, `NM_UF`, área total do município
  (`AREA_KM2`), área dentro do Cerrado (`AREA_CERRADO_KM2`), área
  desmatada, % da área dentro do Cerrado desmatada, número de polígonos
  PRODES com interseção positiva, Z score, p-valor e classificação de
  hotspot (`gi_bin`), tanto para área absoluta quanto para taxa (`_pct`).
  Um mesmo polígono pode ser contado em mais de um município quando cruza
  limites municipais, mas sua área é repartida pelas interseções.
- `resumo_hotspots_por_estado_cerrado.csv`: contagem de municípios por
  estado, variável analisada (área absoluta ou percentual) e classe Gi*.
- `mapa_hotspots_area_desmatada_cerrado.png` /
  `mapa_hotspots_pct_desmatado_cerrado.png`:
  mapas coroplético dos municípios coloridos pela classificação de hotspot.

## Observação sobre a geração dos mapas

A geração dos mapas PNG é feita através de `src/make_maps.py`, chamado como
um subprocesso separado por `run_hotspot_analysis.py`. Isso isola qualquer
falha de renderização (bibliotecas gráficas nativas) do restante da análise,
que é o resultado principal. O renderizador padrão usado
(`visualization.render_hotspot_map_pil`) usa Pillow em vez de matplotlib,
pois em alguns ambientes a instalação de fontes TrueType/freetype usada pelo
matplotlib pode estar corrompida; há também uma implementação alternativa
via matplotlib/geopandas (`visualization.plot_hotspot_map`) disponível para
ambientes onde o matplotlib funcione normalmente.

## Documentação

- [Metodologia](docs/metodologia.md)
- [Resultados e discussão](docs/resultados.md)
