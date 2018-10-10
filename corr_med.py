import numpy as np

def rv(fcst, fcst_mean, mean_obs):

	cor = mean_obs + (fcst - fcst_mean)*(mean_obs/fcst_mean)
	cor[cor < 0] = 0.0

	return cor