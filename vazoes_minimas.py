import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def curva_permanencia(flow):
    flow.sort(reverse=True)
    x = np.linspace(1/len(flow), 100, len(flow))
    plt.plot(x, flow)
    plt.xlabel('Probabilidade (%)')
    plt.ylabel('Vazão (m³/s)')
    plt.title('Curva de permanência')
    plt.show()


def vazao_minima(flow, q=90):
    '''
    
    Parameters
    ----------
    flow : list
        Série de vazão não ordenada.
    q : int, optional
        Especifica vazão mínima. The default is 90.

    Returns
    -------
    float
        Vazão mínima para para o percentil q informado.

    '''
    flow.sort(reverse=True)
    x = np.linspace(1/len(flow), 100, len(flow))
    p = np.where(x >= q)[0][0]
    print(x)
    return flow[p]


def q710(flow, dates):
    '''
    

    Parameters
    ----------
    flow : list
        Série de vazões.
    dates : list
        Datas referentes a cada uma das vazões.

    Returns
    -------
    q : float
        Vazão mínima de 7 dias com 10 anos de retorno.

    '''
    mean = [np.mean(flow[i:i+7]) for i in range(len(flow)-7)] # media movel de 7 dias
    dates = dates[7:]                                         # datas referente as medias
    df = pd.DataFrame(mean, index=dates, columns=['flow'])    # dataframe reunido vazoes e datas
    df.index = pd.to_datetime(df.index)                       # converte datas para datatype
    df['year'] = [x.year for x in df.index]                   # cria coluna com informacao dos anos
    df_min = df.groupby(df.index.year).transform('min').drop_duplicates() # calcula vazoes minimas de cada ano
    x_mean = df_min['flow'].mean()
    s = df_min['flow'].std()
    
    q = x_mean + s * (0.45 + 0.7797 * np.ln( np.ln( 10 / 9)))
    
    return q
    
