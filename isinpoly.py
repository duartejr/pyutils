#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 21:33:17 2018

@author: duarte
"""

import numpy as np
import matplotlib.path as mpltpath


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
    isin = np.where(isin == True, 1., 0.)

    return isin