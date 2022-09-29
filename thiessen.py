import pandas as pd
import numpy as np
import geopandas as gpd

from shapely.ops import unary_union
from shapely.geometry import Point

from calendar import monthrange
from datetime import date, timedelta
from dateutil.relativedelta import relativedelta

from geovoronoi import points_to_coords
from geovoronoi import calculate_polygon_areas
from geovoronoi import voronoi_regions_from_coords

__author__: "Duarte Junior"
__licence__: "GPL"
__version__: "1.2.1"
__maintener__: "Duarte Junior"
__emal__: "duarte.junior@alu.ufc.br"

"""Este arquivo contém o conjunto de funções necessárias para realizar o 
cálculo do Thiessen utilizando dados de precipitação.

A função principal para o cálculo do Thiesse é a de nome `thiessen`.

Arquivos necessãrios para o cálculo do Thiesse:
 - Shapefile (.shp) da área de interesse.
 - Arquivos (.csv) com os dados de precipitação de postos pluviométricos dentro da área de interesse.
 
Estrutura dos arquivos (.csv) de precipitação.
data,latitude,longitude,pr
2020-04-01,-20.1475,-46.2847222222,27.2
2020-04-02,-20.1475,-46.2847222222,0.0
2020-04-03,-20.1475,-46.2847222222,4.1
2020-04-04,-20.1475,-46.2847222222,0.0
2020-04-05,-20.1475,-46.2847222222,0.0
2020-04-06,-20.1475,-46.2847222222,0.0
...

Exemplos de chamada da função `thiessen`

Thiessen Diário:
Este exemplo cálcula a precipitação média diária para o mês de janeiro de 2020
shp = 'diretorio/arquivo/shape/shape.shp'
pr_files = ['dados_precipitacao/arquivo_01.csv',
            'dados_precipitacao/arquivo_02.csv',
            'dados_precipitacao/arquivo_03.csv',
            'dados_precipitacao/arquivo_04.csv',
            'dados_precipittacao/arquivo_05.csv',
            ...]
start_date = date(2020, 1, 1)
start_date = date(2020, 1, 31)

pr_diaria = thiessen(shp, pr_files, start_date, end_date)

Thiessen Mensal:
Este exemplo cálcula a precipitação média mensal para os meses de janeiro a março de 2020
shp = 'diretorio/arquivo/shape/shape.shp'
pr_files = ['dados_precipitacao/arquivo_01.csv',
            'dados_precipitacao/arquivo_02.csv',
            'dados_precipitacao/arquivo_03.csv',
            'dados_precipitacao/arquivo_04.csv',
            'dados_precipittacao/arquivo_05.csv',
            ...]
start_date = date(2020, 1, 1)
start_date = date(2020, 3, 31)

pr_mensal = thiessen(shp, pr_files, start_date, end_date, kind='monthly')


"""

def thiessen(shp: str, pr_files: str, start_date: date, end_date: date, 
             kind='daily') -> pd.DataFrame: 
    """Cálculo do Thiessen

    Args:
        shp (str): Endereço do arquivo shapefile (.shp) da área de interesse.
        pr_files (str): Diretório que contém os arquivos de precipitação.
                        Os arquivos de precipitação devem estar no formato .csv.
                        Os dados nos arquivos devem estar estruturados 4 colunas.
                        As colunas dos arquivos devem ser: data, latitude, longitude, pr
                        A data dos registros deve estar no formato: YYYY-MM-DD
                        Laitude e longitude devem estar em graus.
        start_date (date): Data de início para o cálculo do Thiessen.
        end_date (date):Data de termínio para o cálculo do Thiessen.
        kind (str, optional): Tipo de Thiessen: 'daily' para diário e 'monthly' para mensal. O padrão é 'daily'.

    Returns:
        pd.DataFrame: Dataframe com a precipitação média na área de interesse.
                         O Dataframe de retorno têm duas colunas: data e pr_media
    """
    if kind == 'monthly':
        pr_mean = monthly(start_date, end_date, pr_files, shp) # Cálculo do Thiessen mensal
    
    else:
        pr_mean = daily(start_date, end_date, pr_files, shp) # Cálculo do Thiessen diário
    
    
    return pr_mean

def voronoi(coords: gpd.GeoDataFrame, shp: str) -> np.array:
    """Gera os polígonos de Voronoi. https://en.wikipedia.org/wiki/Voronoi_diagram

    Args:
        coords (gpd.GeoDataFrame): Objeto do Geopandas com as coordenadas (lat, lon) dos postos pluviométricos.
        shp (str): Endereço do arquivo shapefile da área de estudo.

    Returns:
        np.array: Array com a área dos polígonos de Voronoi.
    """
    poly_shp, _ = voronoi_regions_from_coords(coords, shp)
    voronoi_areas = calculate_polygon_areas(poly_shp)
    return np.array(list(voronoi_areas.values()))


def geo_coords(coords: np.array) -> gpd.GeoDataFrame:
    """Transforma array de coordenadas em objeto do Geopandas.

    Args:
        coords (np.array): Array com as coordenas dos postos pluviométricos.

    Returns:
        gpd.GeoDataFrame: Geodataframe das coordenadas dos postos.
    """
    df1 = pd.DataFrame({'X':coords[:,1], 'Y':coords[:,0]})
    df1['geometry'] = list(zip(df1['X'], df1['Y']))
    df1['geometry'] = df1['geometry'].apply(Point)
    gdf1 = gpd.GeoDataFrame(df1, geometry='geometry')
    return gdf1

def read_pr(files: list, dt: date, monthly: str = False) -> pd.pd.DataFrame:
    """Realiza a leituda dos dados de precipitação para uma determinada data.

    Args:
        files (list): Lista de arquivos com os dados de precipitação.
        dt (date): Data de referência para leitura dos dados.
        monthly (str, optional): Se True calcula o acumulado mensal de cada postos. Defaults to False.

    Returns:
        pd.pd.DataFrame: Dados de precipitação de cada um dos postos dentro da área de interesse.
    """
    
    for i, pr_file in enumerate(files):
        pr = pd.read_csv(pr_file, parse_dates=['data'])
        
        if monthly:
            n_days = monthrange(dt.year, dt.month)[1] - 1                   # Número de dias do mês da data de referência
            end_day = dt + timedelta(days=n_days)                           # Último dia do mês da data de referência
            pr = pr[(pr['data'] >= str(dt)) & (pr['data'] <= str(end_day))] # Seleciona todos os registros dentro do mês de referência
            
            # Caso não haja registros de precipitação no período, Ou o nḿero de
            # registros seja inferior a 20 o posto é ignorado e o programa segue
            # com a leitura dos próximos postos
            if len(pr) <= 20:
                continue
            else:
                pr_m = pr['pr'].sum()                                       # Se há registros suficientes calcula o total de precipitação no mês
            
            # Formata o resultado obtido
            pr = pd.DataFrame([pr['data'].values[0],
                            pr['latitude'].values[0],
                            pr['longitude'].values[0],
                            pr_m]).T
        
        else:                                                               # Executado para o cálculo da precipitação diária
            pr = pr[pr['data'] == str(dt)]
            
            if len(pr) == 0:
                continue
        
        if i == 0:                          
            pr_out = pr
        else:
            pr_out = pd.concat([pr_out, pr])
    
    pr_out.columns = ['data', 'latitude', 'longitude', 'pr']
    
    return pr_out

def monthly(start_date: date, end_date: date, files: list, shp: str) -> pd.DataFrame:
    """Cálculo do Thiessen mensal

    Args:
        start_date (date): Data inicial para o cálculo do Thiessen.
        end_date (date): Data final para o cálculo do Thiessen.
        files (list): Lista de arquivos com os dados de precipitação dos postos pluviométricos.
        shp (str): Endereço do arquivo shapefile (.shp) da área de interesse.

    Returns:
        pd.DataFrame: Dados de precipitação média mensal para a área de interesse.
    """
    pr_mean = []
    boundary = gpd.read_file(shp)
    boundary_shape = unary_union(boundary.geometry)
    print('Calculando o Thiessen para:')
    delta = end_date.month -start_date.month
    
    for i in range(delta+1):
        dt = start_date + relativedelta(months=i)
        print(dt)
        df_pr_month = read_pr(files, dt, monthly=True)
        pr_month = df_pr_month['pr'].values
        coords = df_pr_month[['latitude', 'longitude']].values
        not_nan = np.where(~(pr_month == np.nan))
        pr_month = pr_month[not_nan]
        coords = coords[not_nan]
        
        # É necessário que a área de interesse tenha ao menos três postos
        # pluviométricos para o programa consiga traçar os polígonos de Voronoi
        if len(coords) > 2:
            coords = geo_coords(coords).set_crs(boundary.crs)
            coords = points_to_coords(coords.geometry)
            voronoi_areas = voronoi(coords, boundary_shape)
            pr_th = sum(pr_month * voronoi_areas) / voronoi_areas.sum()
            pr_mean.append([dt, pr_th])
        
    return pd.DataFrame(pr_mean, columns=['data', 'pr_media'])


def daily(start_date: date, end_date: date, files: list, shp: str) -> pd.DataFrame:
    """Calcula o Thiessen diário

    Args:
        start_date (date): Data inicial para o cálculo do Thiessen.
        end_date (date): Data final para o cálculo do Thiessen.
        files (list): Lista de arquivos com os dados de precipitação dos postos pluviométricos.
        shp (str): Endereço do arquivo shapefile (.shp) da área de interesse.

    Returns:
        pd.DataFrame: Dados de precipitação média diára para a área de interesse
    """
    pr_mean = []
    boundary = gpd.read_file(shp)
    boundary_shape = unary_union(boundary.geometry)

    print('Calculando o Thiessen para:')
    delta = end_date - start_date
    for i in range(delta.days + 1):
        dt  = start_date + timedelta(days=i)
        print(dt)
        df_pr_day = read_pr(files, dt)
        pr_day = df_pr_day['pr'].values
        coords = df_pr_day[['latitude', 'longitude']].values
        not_nan = np.where(~(pr_day == np.nan))
        pr_day = pr_day[not_nan]
        coords = coords[not_nan]
        
        # É necessário que a área de interesse tenha ao menos três postos
        # pluviométricos para o programa consiga traçar os polígonos de Voronoi
        if len(coords) > 2:
            coords = geo_coords(coords).set_crs(boundary.crs)
            coords = points_to_coords(coords.geometry)
            voronoi_areas = voronoi(coords, boundary_shape)
            pr_th = sum(pr_day * voronoi_areas) / voronoi_areas.sum()
            pr_mean.append([dt, pr_th])

    return pd.DataFrame(pr_mean, columns=['data', 'pr_media'])
