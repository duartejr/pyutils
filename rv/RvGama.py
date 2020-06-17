# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import gamma


def rv_gama(hind, obs, fcst):

    # Ajuste da distribuição gama para os dados observados
    obs = np.sort(obs)
    n_zeros_obs = len(np.where(obs == 0)[0])
    q = n_zeros_obs/float(len(obs))
    obs = obs[obs > 0]
    gamma_obs = gamma.fit(obs, floc=0)

    # Ajuste da distrivbuição gama para os dados do modelo
    hind = np.sort(hind)
    hind = hind[n_zeros_obs:]
    hind = hind[hind > 0]
    gamma_hind = gamma.fit(hind, floc=0)
    # Correção da previsão usando função gama
    if fcst < hind[0]:
        corr = 0
    else:
        prob_mod = gamma.cdf(fcst, *gamma_hind)
        H = q + (1-q)*prob_mod
        corr = gamma.ppf(H, *gamma_obs)
        if str(corr) == 'inf':
            print('hey there')

    return corr
