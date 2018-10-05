#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 22:26:16 2018

@author: duarte
"""
import numpy as np


def bhthorthwaite(pr, etp, cad=100):

    ALT = [0]
    EXC = [0]
    DEF = [0]
    Nac = []
    ARM = []
    ALT = []
    ETR = []
    DEF = []
    EXC = []
    
    est_umid = True
    
    for i in range(len(pr)):
        
        if est_umid:
            if pr[i] - etp[i] > 0:
                Nac.append(0)
                ARM.append(cad)
                ALT.append(0)
                ETR.append(etp[i])
                DEF.append(0)
                EXC.append(pr[i]-etp[i])
            else:
                est_umid = False
                Nac.append(pr[i]-etp[i])
                ARM.append(cad*np.exp(Nac[i]/cad))
                ALT.append(ARM[i]-ARM[i-1])
                ETR.append(pr[i]+np.abs(ALT[i]))
                DEF.append(etp[i]-ETR[i])
                
                if ARM[i] <= cad:
                    EXC.append(0)
                else:
                    EXC.append(pr[i]-etp[i]-ALT[i])
        else:
            if pr[i]-etp[i]<0:
                Nac.append(pr[i]-etp[i])
                ARM.append(cad*np.exp(Nac[i]/cad))
                ALT.append(ARM[i]-ARM[i-1])
                ETR.append(pr[i]+np.abs(ALT[i]))
                
            else:
                ARM.append(ARM[i-1]+(pr[i]-etp[i]))
                
                if ARM[i] > cad:
                    ARM[i] = cad
                
                Nac.append(cad*np.log(ARM[i]/cad))
                ALT.append(ARM[i]-ARM[i-1])
                ETR.append(etp[i])
            
            DEF.append(etp[i]-ETR[i])
            
            if ARM[i] < cad:
                EXC.append(0)
            else:
                EXC.append(pr[i]-etp[i]-ALT[i])
    
    return np.array(DEF), np.array(EXC)
            
            
            
            