import numpy as np

def rv(fcst, hist, obs):
    hist = np.array(hist)
    obs = np.array(obs)
    fcst = np.array(fcst)
    
    nhist = len(hist)//12
    nobs = len(obs)//12
    nfcst = len(fcst)//12
    
    hist = hist[:nhist*12]
    obs = obs[:nobs*12]
    fcst = fcst[:nfcst*12]
    
    hist = np.reshape(hist, (len(hist)//12, 12))
    obs = np.reshape(obs, (len(obs)//12, 12))
    fcst = np.reshape(fcst, (len(fcst)//12, 12))
    
    
    hist_mean = np.nanmean(hist, 0)
    hist_std = np.nanstd(hist, 0)
    obs_mean = np.nanmean(obs, 0)
    obs_std = np.nanstd(obs, 0)
    
    for i in range(len(fcst)):
        for j in range(len(fcst[0])):
            fcst[i,j] = ((fcst[i,j] - hist_mean[j])/hist_std[j])*obs_std[j] + obs_mean[j]
            
            if fcst[i, j] < 0:
                fcst[i,j] = 0
            if np.isnan(fcst[i, j]):
                fcst[i, j] = 0
    
    fcst = np.reshape(fcst, (len(fcst)*12))
    
    return fcst