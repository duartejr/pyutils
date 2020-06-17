import numpy as np

def rv_media(fcst, fcst_mean, mean_obs):

    cor = mean_obs + (fcst - fcst_mean)*(mean_obs/fcst_mean)

    if cor < 0:
        cor = 0.0

    return cor