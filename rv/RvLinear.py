from dateutil.relativedelta import relativedelta
import numpy as np
from sklearn import linear_model
from random import randint
import statsmodels.formula.api as smf
import statsmodels.api as sm


def linear(X, y, e):
    X = np.reshape(X,(len(X),1))
    y = np.reshape(y, (len(y),1))
    regr = linear_model.LinearRegression()
    regr = regr.fit(X, y)
    e = np.reshape(e, (1,1))
    return regr.predict(e)


def gaussian(X, y, e):
    regre = sm.GLM(X, y, family=sm.families.Gaussian())
    #regre = sm.GLM(X, y, family=sm.families.Gamma(link=sm.families.links.log()))
    model = regre.fit()
    return model.predict(e)


def rvlinear1(fcst_cor, erros, family, n_prevs):
    f = np.copy(fcst_cor)
    family = eval(family)
    pos_erro = randint(0, len(erros)-1)
    erro1 = float(erros['d1'].iloc[pos_erro])
    f[0] = f[0]+erro1
    X = np.array(erros['d1'].values, dtype=float)
    
    for d in range(2, n_prevs+1):
        y = np.array(erros['d{}'.format(d)].values, dtype=float)
        erro = family(X, y, erro1)
        f[d-1] = f[d-1]+erro
    
    f[np.where(f<0)] = 0.0
    
    return f


def rvlinear2(fcst_cor, erros, family, n_prevs):
    f = np.copy(fcst_cor)
    family = eval(family)
    pos_erro = randint(0, len(erros)-1)
    erro1 = float(erros['d1'].iloc[pos_erro])
    f[0] = f[0]+erro1
    
    for d in range(2, n_prevs+1):
        X = np.array(erros['d{}'.format(d-1)].values, dtype=float)
        y = np.array(erros['d{}'.format(d)].values, dtype=float)
        erro = family(X, y, erro1)
        f[d-1] = f[d-1]+erro[0]
        erro1 = erro
    
    f[np.where(f<0)] = 0.0
    
    return f


def rv(fcst, erros, method_name='rvlinear1', family='linear', n_cens=1,
       n_prevs=10):
    '''
    Função para correção a partir da regressão dos erros.
    Variáveis de entrada:
        fcst: valor previsto
        erros: série de erros das n últimas previsões
        method_name: metodologia para regressão (rvlinear1, rvlinear2)
        family: tipo de regressão (linear, gaussina, gamma)
        n_cens: número de cenários para gerar
        n_prevs: horizonte de previsão em dias
    '''
    method = eval(method_name)
    #erros = np.array(erros.values)
    fcst_cor = np.copy(fcst)
    fcst_out = []
    
    for i in range(n_cens):
        fcst_cor = method(fcst_cor, erros, family, n_prevs) #calcula correção a partir de método escolhido
        fcst_out.append(fcst_cor)
    
    fcst_out = np.array(fcst_out).mean(0)
    
    return fcst_out
    