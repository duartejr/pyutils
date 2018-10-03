import numpy as np
from numpy import arccos, tan, sin, pi, cos

def hargreaves(t_med, t_max, t_min, y):
    """Cálculo da ETP segundo Hargreaves
    Dados:  t_med: temperátura média
            t_max: temperatura máxima
            t_min: temperatura mínima
            y: latitude
    """
    n = np.arange(1, 366)

    # Converte a latitude de graus em radianos
    # Por exemplo, abaixo -5 graus convertidos -0.0873
    # lat = -5*pi/180
    lat = np.radians(y)

    # Declinação Solar
    declination = 0.4093*sin(2*pi*n/365 - 1.405)

    # Sunset out angle (radians)
    w_s = arccos(-tan(lat)*tan(declination))

    # relative distance between Earth and Sun
    d_r = 1 + 0.333*cos(2*pi*n/365)

    #Extraterrestrial solar radiation (mm/dia)
    S_0 = 15.392*d_r*(w_s*sin(lat)*sin(declination)+cos(lat)*cos(declination)*sin(w_s))

    n_dias = [31,28,31,30,31,30,31,31,30,31,30,31]
    n_end = [1]+list(np.cumsum(n_dias))

    # Caculo de S_)m mes a mes (valor medio dentro de cada mes)
    S_0m = []
    for i in range(12):
        ini = n_end[i+1] - n_dias[i] + 1
        fim = n_end[i+1]
        S_0m.append(np.mean(S_0[ini:fim+1])) # Calculo da média diária ao longo do mês

    #Número de anos de dados CRU
    nl = len(t_med)
    n_anos = nl//12
    HG2 = []

    # Cálculo de HG1 (tm/dia) e HG2 (mm/mes) para cada mês de cada ano do histórico

    for i in range(n_anos):
        ini = 12*i
        fim = 12*i+12
        aux = 0.0023*(t_med[ini:fim] + 17.8)*\
              (np.abs(t_max[ini:fim] - t_min[ini:fim]))**.5*S_0m
       
        aux = np.array(aux)
        if len(np.where(aux<0)[0]):
            aux[np.where(aux<0)[0]] = 0
        HG2.append(n_dias*aux)

    return np.array(HG2).reshape(n_anos*12)
