#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 22:26:16 2018

@author: duarte
"""
import numpy as np


def bhthorthwaite(pr, etp, cad=100):

    V = [0]
    EXC = [0]
    DEF = [0]

    for i in range(len(pr)):
        
        if V[i] + pr[i] < etp[i]:
            etr = V[i]+pr[i]
        else:
            etr = etp[i]

        V.append(V[i]+pr[i]-etr)

        if V[i+1] > cad:
            EXC.append(V[i+1]-cad)
            V[i+1] = cad
        else:
            EXC.append(0)

        if etr < etp[i]:
            DEF.append(etp[i]-etr)
        else:
            DEF.append(0)

    return V[1::], EXC[1::], DEF[1::]
