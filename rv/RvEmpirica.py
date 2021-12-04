import numpy as np


def closest(fcst_target, clim_list):
    aux = []

    for ii in range(0, clim_list.shape[0]):
        aux.append(abs(fcst_target-clim_list[ii]))

    min_value = np.nanmin(aux)
    idx_min = np.where(aux == min_value)[-1][-1]

    return idx_min


def rv_empirica(hind, obs, fcst):

    xs_hind = np.sort(hind)
    xs_obs = np.sort(obs)
    idx_hind = closest(fcst, xs_hind)
    corr = xs_obs[idx_hind]

    return corr
