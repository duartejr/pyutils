#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 08:24:49 2018

@author: duarte
"""
import numpy as np


def ccs(surplus, deficiency, etp):

    Ih = (surplus/etp)*100
    Ia = (deficiency/etp)*100

    Iu = np.nanmean(Ih - 0.6*Ia)

    if 100 <= Iu:
        clim_class = 'A'
        id_class = 0
    elif 80 <= Iu < 100:
        clim_class = 'B4'
        id_class = 1
    elif 60 <= Iu < 80:
        clim_class = 'B3'
        id_class = 2
    elif 40 <= Iu < 60:
        clim_class = 'B2'
        id_class = 3
    elif 20 <= Iu < 40:
        clim_class = 'B1'
        id_class = 4
    elif 0 <= Iu < 20:
        clim_class = 'C2'
        id_class = 5
    elif -20 <= Iu < 0:
        clim_class = 'C1'
        id_class = 6
    elif -40 <= Iu < -20:
        clim_class = 'D'
        id_class = 7
    else:
        clim_class = 'E'
        id_class = 8

    return Iu, clim_class, id_class
