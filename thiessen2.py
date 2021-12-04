'''
Conjunto de funções para processamento do Thiessen
'''
import numpy as np
import pylab as pl
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point, LineString


def thiessen(px, py, bacx, bacy, pr):
    '''
    Realiza o cálculo do Thiessen.
    
    INPUT:
        px (array): Longitude das estações
        py (array): Latitude das estações
        bacx (array): Longitudes do poligono
        bacy (array): Latitudes do poligonoe
        pr (array): Precipitação observada nas estações
    
    OUTPUT:
        pr_med (float): Precipitação média sengundo método de Thiessen
        
    '''
    alpha = voronoi(px, py, bacx, bacy)
    pr_med = np.sum(alpha*pr)
    return pr_med

def areaxy(x, y):
    """
    Calcula  a área de um polígono bidimensional
    formado pelos vértices com vetores de coordenadas x e y.
    O resultado é sensível à direção: a área é
    positiva se o contorno delimitador é antihorário
    e negativa se horário.

    Adaptado da versão em matlab de:
    Copyright (c) 1995 by Kirill K. Pankratov,
    kirill@plume.mit.edu.
    04/20/94, 05/20/95
    """
    # Calcula a integral de contorno Int -y*dx  (mesmo como Int x*dy).
    lx = len(x)-1
    x = x[1:lx+1]-x[0:lx]
    y = y[0:lx]+y[1:lx+1]
    a = -np.dot(x, y/2.)
    return a


def voronoi(px, py, bacx, bacy, fign=False, cont=False):
    """
    Calculo dos pesos de cada ponto.

    verifica os pontos que estão dentro da bacia
    e retorna os coeficientes de influencia de cada um deles.

    :param px: - longitude dos pontos ( ex: lon dos postos)
    :param py: - latitude dos pontos ( ex: lat dos postos)
    :param bacx: - longitude do polígono ( ex: lon das bacias)
    :param bacy: - latitude do polígono ( ex: lat das bacias)
    :param fign: - nome da figura, opcional
    :param cont: - número da figura, opcional
    :return alpha: - retorna os pesos de cada ponto
    """

    # Cálculo do Thiessen
    mnx = np.min(np.append(px, bacx))
    mny = np.min(np.append(py, bacy))

    mxx = np.max(np.append(px, bacx))
    mxy = np.max(np.append(py, bacy))

    lx = mxx - mnx
    ly = mxy - mny

    # Definindo os limites da area de contorno
    mnx = mnx - 10*lx
    mny = mny - 10*ly
    mxx = mxx + 10*lx
    mxy = mxy + 10*ly

    mdx = (mnx + mxx)/2.
    mdy = (mny + mxy)/2.

    # Descomente a linha abaixo para plotar
    # o gráfico da bacias com os pontos
    if fign:
        pl.plot(px, py, 'r^', bacx, bacy, 'b-')

    px = np.append(px, np.array([[mnx], [mdx], [mxx], [mdx]]))
    py = np.append(py, np.array([[mdy], [mxy], [mdy], [mny]]))

    vorpt = np.column_stack([px,py])

    # Gerar voronoi
    vor = Voronoi(vorpt)

    # formatar vetor de polígonos
    bacv = np.column_stack([bacx, bacy])

    area = np.empty(0)

    for region_idx in vor.point_region:
        region = vor.regions[region_idx]

        # Se a região é finita
        if region and -1 not in region:
            coords = np.empty(0)

            # Obtém coordenadas dos vértices
            for vertex_idx in region:

                coords = np.append(coords, vor.vertices[vertex_idx])

            # Intersecção do posto com o polígono
            coords = coords.reshape(len(coords)//2, 2)

            intersc = Polygon(coords).intersection(Polygon(bacv))

            if isinstance(intersc, Polygon):

                # Extrair lat, lon
                lat, lon = np.array(intersc.exterior.xy)

                # Caso do posto selecionado não ter área de intersecção com a bacia
                if lat.shape[0] == 0:
                    area = np.append(area, 0)
                else:
                    area = np.append(area, areaxy(lat, lon))
                    # Descomente a linha abaixo para plotar
                    # o gráfico da bacias com os pontos
                    if fign:
                        pl.plot(lat, lon, 'b-')
            else:
                aux = 0.

                # Identifica Multipolígono
                for i in range(len(intersc.geoms)):

                    if isinstance(intersc[i], Point) or \
                            isinstance(intersc[i], LineString):
                        continue

                    lat, lon = np.array(intersc[i].exterior.xy)

                    if lat.shape[0] == 0:
                        aux = 0

                    else:
                        aux = aux + areaxy(lat, lon)
                        # Descomente a linha abaixo para plotar
                        # o gráfico da bacias com os pontos
                        if fign:
                            pl.plot(lat, lon, 'b-')

                area = np.append(area, aux)

    # Descomente a linha abaixo para plotar
    # o gráfico da bacias com os pontos
    if fign:
        fign = "{0}_{1}".format(cont, fign)
        pl.savefig(fign)
        pl.close()

    alpha = area/np.sum(area)

    return alpha

