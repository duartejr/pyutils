# -*- coding: utf-8 -*-
import numpy as np
from utils import basinsf


def rv_ons(pr_grad, model, basin, subbasin, month):

    try:
        pr = np.mean(pr_grad[:,2:],0)
    except:
        pr = pr_grad[2:]

    p_tot_mod = np.sum(pr)
    params = basinsf.paramsrv(model, basin, subbasin, month)
    
    if not params:
        return None
    
    p_tot_pr_mod = params[0]*p_tot_mod**2 + params[1]*p_tot_mod

    if p_tot_pr_mod > params[2]:
        p_tot_pr_mod = p_tot_mod
    
    if p_tot_mod == 0:
        p_pr = np.array(pr) * p_tot_pr_mod
    else:
        p_pr = np.array(pr) * p_tot_pr_mod / p_tot_mod
    
    beta = p_pr/pr
    beta[np.isnan(beta)] = 0

    try:
        pr_grad[:,2:] = pr_grad[:,2:]*beta
    except:
        pr_grad[2:] = pr_grad[2:]*beta
    
    return pr_grad