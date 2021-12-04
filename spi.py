# -*- coding: utf-8 -*-

__author__ = ["Paulo Jarbas Camurca"]
__credits__ = ["Fco Vasconcelos", "Marcelo Rodrigues", "Daniel Pinto"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "pjarbas312@gmail.com"


import numpy as np
from scipy.stats import gamma, norm


def spi(data, scale, nseas):

    '''
    Cálculo do SPI
    Adaptado de  Taesam Lee,  Dec.03,2009 INRS-ETE, Quebec, Canada

    param: data  - Vetor de dados mensais
    param: scale - Valor da escala do spi: 1, 3, 4, ou 6
    param: nseas - Número de estações
    return:      - O valor do spi e os parametros da função gama
    '''

    # Data setting to scaled dataset
    # A1 = np.zeros((0, scale))
    A1 = []

    for i in range(scale):
        A1.append(data[i:len(data)-scale+i+1])
    
    A1 = np.array(A1).T
    
    if A1.ndim == 2:
        XS = A1.sum(axis=1)
    else:
        XS = A1

    sz = XS.shape
    Z = np.zeros(sz)
    phat = np.zeros((nseas, 2))

    for i in range(nseas):
        tind = np.arange(i, len(XS), nseas)
        Xn = XS[tind]
        zeroa = Xn[np.where(Xn==0)]
        nnan = np.sum(np.isnan(Xn))
        Xn_nozero = Xn[np.nonzero(Xn)]  # removing zeros
        Xn_nozero = Xn_nozero[~np.isnan(Xn_nozero)] # removing nan        
        q = len(zeroa)/np.float(len(Xn) - nnan)
        parm = gamma.fit(Xn_nozero, floc=0)
        phat[i, 0] = parm[0]
        phat[i, 1] = parm[2]

        Gam_xs = q+(1-q)*gamma.cdf(Xn, parm[0], scale=parm[2])
        Z[tind] = norm.ppf(Gam_xs)

    return Z
