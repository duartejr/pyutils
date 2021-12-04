#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 21:33:17 2018

@author: duarte
"""

import numpy as np
import matplotlib.path as mpltpath
import fiona


def isinpoly3(lon, lat, polig):
    """
    Encontra se os pontos com coordenadas lon e lat
    estão dentro ou fora de um poligono.

    :param lon: - array com a longitude dos pontos
    :param lat: - array com a latitude dos pontos
    :param polig: -  matriz com os pontos do polígono
    :return isin: - array booleano com o valor 0 para pontos fora do polígono, 1 para
                    pontos dentro do polígono. Pontos encima do polígono recebem 0.
    """

    coords = np.c_[lon, lat]

    pth = mpltpath.Path(polig,  closed=True)
    isin = pth.contains_points(coords)

    isin = np.array(isin)
#    isin = np.where(isin == True, 1., 0.)

    return isin

def isinshp(lon, lat, polig):
    lon[lon>180] -= 360
    if len(lon.shape) and len(lat.shape) == 1:
        LON, LAT = np.meshgrid(lon, lat)
        LON = np.reshape(LON, len(lon)*len(lat))
        LAT = np.reshape(LAT, len(lon)*len(lat))
        isin = isinpoly3(LON, LAT, polig)
        isin = np.reshape(isin, (len(lat), len(lon)))
    else:
        LON = np.reshape(lon, len(lon)*len(lon[0]))
        LAT = np.reshape(lat, len(lat)*len(lat[0]))
        isin = isinpoly3(LON, LAT, polig)
        isin = np.reshape(isin, lon.shape)
    
    
    
    return isin