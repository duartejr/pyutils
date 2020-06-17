import numpy as np
from numpy import arccos, tan, sin, pi, cos, log, exp, sqrt, arctan

def hargreaves(t_med, t_max, t_min, y, daily=False):
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
    # Cálculo de HG1 (tm/dia) e HG2 (mm/mes) para cada mês de cada ano do histórico
    HG2 = []
    HG1 = []
    
    if daily:
        while len(S_0) < len(t_med):
            S_0 = np.concatenate((S_0,S_0))
        S_0 = S_0[:len(t_med)]
        HG1 = 0.0023*(t_med+17.8)*np.abs(t_max-t_min)**.5*S_0
        posnegs = np.where(HG1<0)[0]

        if len(posnegs):
            HG1[posnegs[:]] = 0.0
        return np.array(HG1)
    else:
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

def thornthwaite(Ti, n):
    n_d = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * int(len(Ti)/12)
    n_d = np.array(n_d)
    n = n / n_d
    I = []
    for i in range(0, len(Ti), 12):
        I_m = [sum((0.2 * Ti[i:i+12])**1.514)] * 12
        for j in I_m:
            I.append(j)
    I = np.array(I)
    a = 6.75e-7 * I**3 - 7.71e-5 * I**2 + 1.7912e-2 * I + 0.49239
    etpp = 16 * (10 * Ti/I ) ** a
    etp = etpp * n /12 * n_d/30
    return etp

def penman_monteith(Tmed=False, Tmax=False, Tmin=False, Rn=False, lat=False,
                    U2=False, Patm=False, alt=False, UR=False, Rs=False,
                    n=False, albedo=0.23, Uz=False, G=0.0, a=0.25, b=0.5, z=10,
                    date=False):
    '''
    Parameters
    ----------
    Tmed : TYPE, optional
        DESCRIPTION. Temperatura média (ºC).
    Tmax : TYPE, optional
        DESCRIPTION. Temperatura máxima (ºC).
    Tmin : TYPE, optional
        DESCRIPTION. Temperatura mínma (ºC).
    Rn : TYPE, optional
        DESCRIPTION. Saldo de radiação diário (MJm⁻²dia⁻1).
    lat : TYPE, optional
        DESCRIPTION. Latitude.
    U2 : TYPE, optional
        DESCRIPTION. Velocidade do vento a 2 metros de altura (m/s).
    Patm : TYPE, optional
        DESCRIPTION. Pressão atmosférica (kPa).
    alt : TYPE, optional
        DESCRIPTION. Altitude do local (m).
    UR : TYPE, optional
        DESCRIPTION. Umidade relativa média do ar (%).
    Rs : TYPE, optional
        DESCRIPTION. Radiação solar incidente (MJm⁻²dia⁻¹).
    n : TYPE, optional
        DESCRIPTION. Número diário de horas de sol (h).
    albedo : TYPE, optional
        DESCRIPTION. Coeficiente de reglexão vegetal. The default is 0.23.
    Uz : TYPE, optional
        DESCRIPTION. Velocidade do vento a altura z (m/s).
    G : TYPE, optional
        DESCRIPTION. Fluxo de calor no solo. The default is 0.0.
    a : TYPE, optional
        DESCRIPTION. The default is 0.25.
    b : TYPE, optional
        DESCRIPTION. The default is 0.5.
    z : TYPE, optional
        DESCRIPTION. Altura das medições de velocidade de vento (m). The default is 10.

    Returns
    -------
    pet : TYPE
        DESCRIPTION. Evapotranspiração potencial.

    '''
    
    if not Patm:
        Patm = 101.3 * ((293 - 0.0065 * alt) / 293)**5.26
    
    es = 0.6108 * exp((17.25 * Tmed) / (Tmed + 237.3))
    
    if isinstance(UR, np.ndarray):
        ea = (es * UR) / 100
    else:
        ea = 0.61 * exp((17.27 * Tmin) / (Tmin + 237.3))

    if not Rn:
        sigma = 4.903e-9
        J = np.array([x.dayofyear for x in date])
        dr = 1 + 0.033 * cos(J * (2 * pi) / 365)
        Ra = (118.08 / pi) * dr
        lat = np.radians(lat)
        Rso = (0.75 + 2e-5 * alt) * Ra
        ds = 0.409 * sin(((2 * pi) / 365) * J - 1.39)
        X = (1 - (tan(lat))**2 * (tan(ds))**2)
        for p,i in enumerate(X):
            if i <=0:
                X[p] = 0.00001
        ws = pi / 2 - arctan((-tan(lat) * tan(ds)) / X**0.5)
        
        if not Rs:
            N = 24 / pi * ws
            Rs = (a + b * n / N) * Ra
        
        Rns = (1 - albedo) * Rs
        
        Ra = 37.6*dr*(ws*sin(lat)*sin(ds)+cos(lat)*cos(ds)*sin(ws))
        Rnl = sigma * (((Tmax + 273.16)**4 + (Tmin + 273.16)**4) / 2) *\
        (0.34 - 0.14 * sqrt(ea)) * (1.35 * (Rs / Rso) - 0.35)
        Rn = Rns - Rnl
    
    if not U2:
        U2 = Uz * (4.87) / (log(67.8 * z - 5.42))
    
    delta = (4098 * (0.6108 * exp((17.27 * Tmed) / (Tmed + 237.3)))) / (Tmed + 237.3)**2
    gama = 0.665e-3 * Patm
    
    pet = (0.408 * delta * (Rn - G) + (gama * 900 * U2 * (es - ea)) / (Tmed + 273))/\
          (delta + gama * (1 + 0.34 * U2))
    
    return pet
    