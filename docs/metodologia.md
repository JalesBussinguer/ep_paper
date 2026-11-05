# Metodologia

## 1. Objetivo e questões de pesquisa

O estudo identifica as áreas de maior pressão de desmatamento no bioma Cerrado entre 2020 e 2025 e investiga o destino da cobertura da terra nos polígonos desmatados. A análise é orientada por três questões:

1. **Onde** a pressão de desmatamento é espacialmente concentrada, em termos de área absoluta e de proporção do território municipal?
2. **O que ocorre** com a cobertura da terra após a detecção do desmatamento (permanência, conversão para uso antrópico ou regeneração)?
3. **Como** essa trajetória difere entre os municípios classificados como *hotspots* e os demais municípios do bioma?

A unidade territorial de análise é o município (recortado pelos limites do bioma). O objeto de análise é o polígono de desmatamento, caracterizado pela localização (centroide) e pela área desmatada.

## 2. Dados

| Conjunto | Arquivo / fonte | Descrição |
|---|---|---|
| Desmatamento | `data/cerrado_desmatamento_2020_2026.geojson` | Polígonos de supressão da vegetação nativa detectados pelo PRODES Cerrado. Os dados disponíveis cobrem os anos de 2020 a 2025. |
| Malha municipal | `data/municipios.geojson` | Malha municipal do IBGE (5.573 municípios). |
| Limite do bioma | `data/cerrado_border.geojson` | Limite do Cerrado (1.996.855,6 km²). |
| Cobertura da terra | `H:\mapbiomas_maps\brazil_coverage-col4_10m_{2019..2025}.tif` | Mapas anuais de uso e cobertura da terra do MapBiomas (Coleção 4, 10 m), de 2019 a 2025, em EPSG:4326, tipo `uint8` e grade comum. |

Os dados vetoriais de entrada estão em SIRGAS 2000 (EPSG:4674). As áreas e os comprimentos foram calculados em SIRGAS 2000 / Policônica (EPSG:5880), sistema métrico apropriado para estatísticas de área em escala nacional. Geometrias inválidas foram reparadas antes de qualquer operação espacial.

## 3. Pré-processamento e unidade de análise

### 3.1 Recorte da malha municipal pelo bioma

A malha municipal foi interseccionada com o limite do Cerrado. Cada município passou a ser representado apenas pela parcela contida no bioma, e a área dessa parcela (`AREA_CERRADO_KM2`) foi usada como denominador de todas as proporções. A malha resultante contém 1.435 municípios, distribuídos em 13 unidades da federação, e cobre 99,9984% da área do bioma.

### 3.2 Agregação por interseção polígono–município

Cada polígono do PRODES foi interseccionado com a geometria município ∩ Cerrado. A área desmatada foi repartida entre as interseções, de modo que um polígono que atravessa o limite entre dois municípios contribui para ambos na proporção da área contida em cada um. Em consequência, a **unidade analítica é a interseção polígono–município**, e não o evento original. O procedimento resultou em 1.972.375 interseções e em 81.788 km² de área desmatada no período. Trinta e sete municípios não apresentaram desmatamento.

Duas variáveis foram derivadas por município:

- **Área desmatada** (`area_desmatada_km2`): soma das áreas das interseções, em km², que mede a magnitude absoluta da pressão.
- **Proporção desmatada** (`pct_area_desmatada`): razão entre a área desmatada e `AREA_CERRADO_KM2`, em percentual, que mede a intensidade relativa ao território.

O uso das duas variáveis é deliberado: a área absoluta favorece municípios extensos, ao passo que a proporção favorece municípios pequenos. A análise conjunta reduz o viés de cada medida isolada.

## 4. Análise de hotspots por Getis-Ord Gi\*

### 4.1 Estatística

Para cada município *i* e variável *x*, a estatística Gi\* (Getis e Ord, 1992; Ord e Getis, 1995) é calculada como

$$G_i^* = \frac{\sum_{j} w_{ij}\,x_j - \bar{X}\sum_{j} w_{ij}}{S\sqrt{\dfrac{n\sum_{j} w_{ij}^2 - \left(\sum_{j} w_{ij}\right)^2}{n-1}}},$$

em que *w<sub>ij</sub>* é o peso espacial entre *i* e *j*, com *w<sub>ii</sub>* > 0 (versão "estrela", que inclui o próprio município), *n* é o número de municípios, $\bar{X}$ é a média e *S* é o desvio padrão da variável. O valor resultante é um escore Z: valores positivos elevados indicam agrupamento de valores altos (*hot spots*) e valores negativos elevados indicam agrupamento de valores baixos (*cold spots*). A implementação segue a documentação da Esri para a ferramenta *Hot Spot Analysis (Getis-Ord Gi\*)* e utiliza a classe `G_Local(star=True)` da biblioteca `esda`.

### 4.2 Pesos espaciais

Adotou-se a contiguidade **Queen** com pesos binários. Municípios sem vizinhos (ilhas) foram conectados ao vizinho mais próximo, para que a estatística fosse definida em toda a malha.

### 4.3 Inferência e classificação

A significância foi avaliada por 999 permutações condicionais (semente 42), que fornecem um pseudo valor-*p*. A classificação em níveis de confiança de 90%, 95% e 99% foi feita pelos limiares do escore Z (1,645; 1,960; 2,576) e pelo pseudo valor-*p*, conforme a documentação da Esri. As classes resultantes são *Hot Spot* e *Cold Spot* com 90%, 95% e 99% de confiança, além de *Não significante*.

### 4.4 Definição de hotspot e de controle

Um município foi classificado como **hotspot** se foi *Hot Spot*, em qualquer nível de confiança, na análise da área desmatada **ou** na da proporção desmatada. O critério resultou em 348 municípios, que concentram 54.727 km² (66,9% da área desmatada total). Os demais municípios constituem o grupo de **controle** (884 municípios com amostra disponível na etapa seguinte).

## 5. Estabilidade temporal dos hotspots

O Gi\* foi recalculado (i) para dois períodos (2020–2022 e 2023–2025) e (ii) para cada ano individualmente, com as áreas anuais por município. Os municípios foram classificados quanto à persistência em **persistentes** (hotspot nos dois períodos), **emergentes** (apenas em 2023–2025), **extintos** (apenas em 2020–2022) e **sem hotspot**. Calculou-se ainda o centro de gravidade da área desmatada, ponderado pela área, para cada ano, com o objetivo de caracterizar o deslocamento da fronteira de desmatamento.

## 6. Análises de robustez do Gi\*

1. **Autocorrelação espacial global:** I de Moran (Queen, 999 permutações), para a área e para a proporção, nas escalas original e logarítmica, como condição prévia à interpretação de agrupamentos.
2. **Definição de vizinhança:** comparação da configuração base com (a) Queen de segunda ordem, (b) 8 vizinhos mais próximos e (c) banda de distância de 150 km. A concordância com a configuração base foi medida pelo índice de Jaccard.
3. **Testes múltiplos:** correção por taxa de descobertas falsas (FDR) de Benjamini-Hochberg (*q* = 0,05), aplicada aos valores-*p* bilaterais derivados do escore Z.
4. **Reprodutibilidade interna:** verificou-se a concordância com o Gi\* original (100%).

## 7. Análise de trajetória da cobertura da terra

### 7.1 Amostragem dos polígonos

Como o conjunto contém aproximadamente 2 milhões de interseções, empregou-se amostragem de Poisson estratificada pela área da interseção, com probabilidades de inclusão crescentes com o tamanho:

| Classe de área (km²) | Probabilidade de inclusão |
|---|---|
| < 0,005 | 0,005 |
| 0,005 – 0,02 | 0,01 |
| 0,02 – 0,1 | 0,03 |
| 0,1 – 0,5 | 0,15 |
| 0,5 – 2 | 0,5 |
| 2 – 10 | 1 |
| ≥ 10 | 1 |

A semente é 20260101. Cada fragmento recebe peso *w* = 1/*p*, e as estimativas de totais seguem o estimador de Horvitz-Thompson. A mesma regra foi aplicada a hotspots e ao controle, o que preserva a comparabilidade. A amostra final contém 50.381 fragmentos, cuja soma de pesos é 1.978.383 (cerca de 0,3% acima do total de interseções, por efeito da aleatoriedade da amostragem).

### 7.2 Extração por pixel

Para cada fragmento amostrado, extraiu-se a janela correspondente nos sete rasters (2019–2025). Foram considerados os pixels cujo centro está contido no fragmento (regra do centro do pixel). A trajetória de cada pixel foi codificada em um inteiro de 64 bits (um byte por ano) e as contagens foram agregadas por fragmento e por trajetória. O processamento foi paralelizado (12 *workers*). Quatrocentos e noventa e um fragmentos, muito pequenos, não continham nenhum centro de pixel.

**Controle de qualidade.** A área estimada pela soma de pixels foi comparada à área geométrica: 54.352 km² contra 54.886 km² nos hotspots e 26.883 km² contra 27.061 km² no controle. A razão mediana pixel/geometria por fragmento foi de 0,993.

### 7.3 Agrupamento das classes

As classes do MapBiomas foram agrupadas, segundo a legenda oficial, em: **vegetação natural** (formações florestais, savânicas e campestres, entre outras), **uso antrópico** (pastagem, agricultura, mosaico de agricultura e pastagem, silvicultura e infraestrutura urbana), **superfície não vegetada**, **água** e **sem dados**. O código 91, que representa cerca de 0,002% dos pixels e não consta da legenda da Coleção 4, foi tratado como "sem dados".

### 7.4 Referências temporais e tipologia

Para um polígono detectado no ano *t*, definiram-se: **pré-evento** (*t* − 1), **detecção** (*t*) e **pós-evento** (*t* + 1). A coorte de 2025 não possui ano posterior e foi excluída das comparações principais (`year ≤ 2024`), embora seja descrita separadamente. Cada pixel foi classificado em uma tipologia mutuamente exclusiva, de acordo com a sequência de grupos (N = vegetação natural; A = uso antrópico):

| Tipo | Sequência (pré, detecção, pós) | Interpretação |
|---|---|---|
| Conversão no ano da detecção | N → A → A | A conversão é observada no mesmo ano da detecção. |
| Conversão posterior | N, N → A | A conversão ocorre no ano seguinte ao da detecção. |
| Vegetação persistente | N, N, N | Não há conversão observada. |
| Já antrópico antes do evento | A, A, não N | O pixel já era antrópico antes da detecção. |
| Regeneração | A → N | Retorno à vegetação natural. |
| Outras transições | demais | Inclui superfícies não vegetadas e água. |
| Sem ano posterior | coorte 2025 | Não classificável. |

Análises de sensibilidade substituíram o pré-evento pelo ano de 2019 e o pós-evento por *t* + 2, combinando as quatro opções.

## 8. Inferência estatística da comparação hotspot × controle

As estimativas são ponderadas pelo inverso da probabilidade de inclusão, e a incerteza é avaliada de acordo com a estrutura hierárquica dos dados (fragmentos aninhados em municípios):

1. **Intervalos de confiança:** *bootstrap* por conglomerado de municípios (500 réplicas, reamostragem multinomial dos municípios), com intervalo percentílico de 95%.
2. **Diferenças** entre hotspot e controle, em pontos percentuais (p.p.), com IC de 95%.
3. **Tamanho de efeito:** V de Cramér, com IC pelo mesmo procedimento.
4. **Teste de hipótese principal:** teste de permutação dos rótulos hotspot/controle entre municípios (2.000 permutações), que preserva a dependência dentro de cada município.
5. **Teste secundário:** qui-quadrado com tamanho efetivo de amostra de Kish. Esse teste é anticonservador e é apresentado apenas como referência.

Avaliou-se ainda o gradiente da conversão em função do nível de confiança do Gi\* (90%, 95% e 99%) e da variável que originou a classificação (área ou proporção).

## 9. Conversão acumulada, destinos e padrão por estado

- **Conversão acumulada:** para cada coorte de detecção, calculou-se, para defasagens de 0 a 5 anos, a proporção da área natural em *t* − 1 que passou a uso antrópico até a defasagem considerada.
- **Destino da conversão:** para os pixels convertidos, identificou-se a classe de destino no primeiro ano antrópico, por município (classe predominante e distribuição em área). Estimativas municipais foram consideradas confiáveis quando o número de fragmentos amostrados era suficiente (255 dos 348 municípios hotspot).
- **Conversão por unidade da federação:** a taxa de conversão foi calculada como a razão entre a área natural em *t* − 1 convertida e a área natural em *t* − 1.

## 10. Validação

1. **Concordância PRODES × MapBiomas:** fração da área do polígono PRODES classificada como vegetação natural no MapBiomas em *t* − 1, o que verifica se o desmatamento detectado incide sobre vegetação nativa segundo uma fonte independente.
2. **Inspeção visual:** amostra estratificada de 120 polígonos (20 por tipologia, com áreas entre 0,05 e 3 km², semente 20260102), com onze pranchas de visualização e *links* para imagens de satélite. A revisão manual por imagem de satélite não foi realizada e permanece como etapa pendente (coluna `revisao_manual`). O alinhamento geométrico entre as duas fontes foi confirmado por inspeção visual das pranchas.

## 11. Reprodutibilidade

Todo o código está em Python 3.13 (comentários em português), organizado em `src/`, e os parâmetros estão centralizados em `src/config.py`. As principais bibliotecas são `geopandas`, `shapely`, `libpysal`, `esda`, `rasterio`, `joblib`, `pandas`, `pyarrow`, `matplotlib` e `pytest`. As sementes aleatórias são fixas (42, 20260101 e 20260102, além das usadas no *bootstrap* e nas permutações).

Ordem de execução:

```text
python run_hotspot_analysis.py     # Gi* original
python -m src.mapbiomas_transitions
python -m src.mapbiomas_event_transition_analysis
python -m src.yearly_aggregation   # áreas anuais e amostra de Poisson
python -m src.pixel_extraction     # trajetórias por pixel
python -m src.pixel_analysis       # tipologias e inferência
python -m src.gi_robustness        # robustez do Gi*
python -m src.visual_validation    # validação visual
python -m src.figures              # figuras em inglês; títulos em output/refinamento/figuras/titulos_figuras.txt
pytest tests -q                    # testes (8)
```

## 12. Limitações

1. A unidade analítica é a interseção polígono–município e não o evento de desmatamento; a contagem de eventos não deve ser inferida diretamente.
2. A resolução temporal é anual; conversões ocorridas dentro do mesmo ano não são ordenáveis.
3. A abordagem descreve associações e não permite inferência causal.
4. A coorte de 2025 não possui ano posterior.
5. A regra do centro do pixel gera pequenas diferenças em relação à área geométrica (razão mediana de 0,993).
6. A definição de hotspot pela união de duas variáveis é mais inclusiva do que cada critério isolado.
7. A coorte de 2020 concentra 48% da área e 84% das interseções, com fragmentos pequenos (mediana de 0,011 km²), o que sugere um regime de mapeamento distinto nesse ano. Recomenda-se confirmar a origem desse padrão com o provedor dos dados e testar a sensibilidade excluindo 2020.
8. A amostragem de Poisson introduz variância própria; os intervalos de confiança refletem sobretudo a variabilidade entre municípios.
9. A inspeção visual por satélite ainda não foi concluída.

## Referências

- Getis, A.; Ord, J. K. The analysis of spatial association by use of distance statistics. *Geographical Analysis*, v. 24, n. 3, p. 189–206, 1992.
- Ord, J. K.; Getis, A. Local spatial autocorrelation statistics: distributional issues and an application. *Geographical Analysis*, v. 27, n. 4, p. 286–306, 1995.
- Benjamini, Y.; Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, v. 57, n. 1, p. 289–300, 1995.
- Horvitz, D. G.; Thompson, D. J. A generalization of sampling without replacement from a finite universe. *Journal of the American Statistical Association*, v. 47, n. 260, p. 663–685, 1952.
- Esri. *How Hot Spot Analysis (Getis-Ord Gi\*) works*. ArcGIS Pro, documentação das ferramentas de estatística espacial.
- INPE. PRODES Cerrado: monitoramento do desmatamento do bioma Cerrado por satélite.
- MapBiomas. Coleção de mapas anuais de cobertura e uso da terra do Brasil.
