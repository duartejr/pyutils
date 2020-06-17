import geopandas as gpd
import geoplot
from geopandas import GeoDataFrame
from shapely.geometry import Point
import numpy as np
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as mpatches
from matplotlib.patches import PathPatch
from matplotlib.colors import LinearSegmentedColormap
import mapclassify as mc


def clip_area(sf, ax):
    
    def coord_lister(geom):
        coords = list(geom.exterior.coords)
        return coords
    
    vertices = sf.geometry.apply(coord_lister)[0]
    codes = [] 
    codes += [Path.MOVETO]
    codes += [Path.LINETO] * (len(vertices) -2)
    codes += [Path.CLOSEPOLY]
    clip = Path(vertices, codes)
    clip = PathPatch(clip, transform=ax.transData)

    return clip


def mycmap(cmap, nbin):
    return LinearSegmentedColormap.from_list('mycmap', cmap, N=nbin+2)


def set_limits(ax, shp):
    env = shp.geometry.bounds
    ax.set_xlim([env.minx[0]-.1, env.maxx[0]+.1])
    ax.set_ylim([env.miny[0]-.1, env.maxy[0]+.1])


def __discret_legend(data, im,  fontsize, ticks_label=0, values=0):
    if not values:
        values = np.sort(np.unique(data.ravel()))
    colors = [ im.cmap(im.norm(value)) for value in values]
    # create a patch (proxy artist) for every color 
    if ticks_label:
        patches = [mpatches.Patch(color=colors[i], label=ticks_label[i]) for i in range(len(values)) ]
    else:
        patches = [ mpatches.Patch(color=colors[i], label="{l}".format(l=values[i]) ) for i in range(len(values)) ]
    # put those patched as legend-handles into the legend
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
               fontsize=fontsize) 


def __colorbar_legend(im, ticks, ticks_label, cbar_extend, vmin, vmax, nbin, cbar_loc,
                      cbar_label, cmap, norm, lbsize, lbpad, fontsize):
    if not ticks:
        ticks = np.linspace(vmin, vmax, nbin+1)
        ticks = list(ticks)
    
    if cbar_extend == 'both':
        ticks.insert(0,ticks[0]-1)
        ticks.append(ticks[-1]+1)
        if not ticks_label:
            ticks_label = ['{0:.2f}'.format(x) for x in ticks[1:-1]]
        ticks2 = ticks[1:-1]
    elif cbar_extend == 'max':
        ticks.append(ticks[-1]+1)
        if not ticks_label:
            ticks_label = ['{0:.2f}'.format(x) for x in ticks[:-1]]
        ticks2 = ticks[:-1]
    elif cbar_extend == 'min':
        ticks.insert(0,ticks[0]-1)
        if not ticks_label:
            ticks_label = ['{0:.2f}'.format(x) for x in ticks[1:]]
        ticks2 = ticks[1:]
    else:
        if not ticks_label:
            ticks_label = ['{0:.2f}'.format(x) for x in ticks]
        ticks2 = ticks
    
    cb = plt.colorbar(im, orientation=cbar_loc, 
              ticks=ticks2,
              extendfrac = 'auto',
              label=cbar_label, extend=cbar_extend,
              boundaries=ticks, fraction=.045,
              cmap=cmap, norm=norm
              )
    cb.ax.tick_params(labelsize=lbsize)
    
    if cbar_loc == 'vertical':
        cb.ax.set_yticklabels(ticks_label)
    else:
        cb.ax.set_xticklabels(ticks_label)
            
    cb.set_label(cbar_label, labelpad=lbpad, size=fontsize)


def __legend(data, im, type_legend, ticks, ticks_label, cbar_extend, vmin, vmax, nbin,
             cbar_loc, cbar_label, cmap, norm, lbsize, lbpad, fontsize):
    
    if type_legend == 'colorbar':
        __colorbar_legend(im, ticks, ticks_label, cbar_extend, vmin, vmax, nbin, cbar_loc, cbar_label, cmap, norm, lbsize, lbpad, fontsize)
    if type_legend == 'discret':
        __discret_legend(data, im, fontsize, ticks_label, ticks)


def imshow(data, lat, lon, shp, clip=None, vmin=0, vmax=100, cmap='seismic_r',
              nbin=10, cbar_loc='vertical', cbar_label='My legend', grid=1,
              cbar_extend='both', lbsize=10, lbpad=3, fontsize=11, aux_shapes=[],
              title='', figname='', show=False, ticks=0, norm=None, figsize=(15,15),
              ticks_label=[], pos_ticks=0, legend='colorbar'):
    '''
    Gera mapas a partir de matrix 2D.

    Parameters
    ----------
    data : array(2D)
        Matriz com dados a serem usado para gerar o mapa.
    lat : array
        Latitudes correspondentes a cada ponto da matriz de dados.
    lon : array
        Longitudes correspondentes a cada ponto da matriz de dados.
    shp : shapefile
        Arquivo shapefile da região onde será plotado o mapa.
    clip : bool, optional
        Se True serão removidas partes da matriz de dados fora do shapefile.
        The default is None.
    vmin : float, optional
        Mínimo valor a ser utilizado para legenda. The default is 0.
    vmax : float, optional
        Máximo valor a ser utilizado para legenda. The default is 100.
    cmap : str or array, optional
        Se str utiliza-se um dos mapas de cores padrões do matplotlib, se array
        define mapa de cor personalizado. The default is 'seismic_r'.
    nbin : int, optional
        Número de divisões/cores para legenda/mapa. The default is 10.
    cbar_loc : str, optional
        Local para posicionar legenda ('vertical' ou 'horizontal'). The default is 'vertical'.
    cbar_label : str, optional
        Nome da legenda. The default is 'My legend'.
    grid : float, optional
        Tamanho da grade a ser utilizada no mapa. The default is 1.
    cbar_extend : str, optional
        Define se colorbar extender valores. The default is 'both'.
    lbsize : int, optional
        Tamanho das labels da lebenda. The default is 10.
    lbpad : int, optional
        Label pad. The default is 3.
    fontsize : int, optional
        Tamanho da fonte padrão do mapa. The default is 11.
    aux_shapes : array(shp), optional
        Lista de array auxiliares com complemento visual. The default is [].
    title : str, optional
        Título do mapa. The default is ''.
    figname : str, optional
        Se fornecido salva mapa no local e com nome especificado. The default is ''.
    show : bool, optional
        Se True mostra figura. The default is False.
    ticks : array, optional
        Ticks da barra de cores. The default is 0.
    norm : bool, optional
        If bool normalize colors. The default is None.

    Returns
    -------
    None.

    '''
    
    gpd_shp = gpd.read_file(shp)
    fig, ax = plt.subplots(figsize=figsize)

    xx,yy = np.meshgrid(lon, lat)
    
    if clip:
        clip = clip_area(gpd_shp, ax)
    
    if not ticks:
        ticks = np.linspace(vmin, vmax, nbin+1)
        ticks = list(ticks)
    
    if type(cmap) != str:
        cmap = mycmap(cmap, nbin)
        if legend == 'colorbar':
            norm = mpl.colors.BoundaryNorm(ticks, cmap.N)
    else:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    
    im = plt.pcolormesh(xx, yy, data, cmap=cmap, clip_path=clip, norm=norm)
                                    
    __legend(data, im, legend, ticks, ticks_label, cbar_extend, vmin, vmax,
             nbin, cbar_loc, cbar_label, cmap, norm, lbsize, lbpad, fontsize)
    
    gpd_shp.plot(ax=ax, facecolor='none', edgecolor='black', lw=2)
    
    if len(aux_shapes):
        for aux_shape in aux_shapes:
            aux_shape = gpd.read_file(aux_shape)
            aux_shape.plot(ax=ax, facecolor='none', edgecolor='black', lw=2)
    
    #draw parallels and meridians
    parallels = np.arange(-90., 90., grid)
    meridians = np.arange(-180., 180., grid)
    ax.hlines(parallels, lon.min(), lon.max(), linewidth=.3)
    ax.vlines(meridians, lat.min(), lat.max(), linewidth=.3)
    xticks = parallels[(parallels>=lon.min()) & (parallels<=lon.max())]
    ax.set_xticks(xticks)
    xtickslabels = [f'{abs(x):0.1f}ºW' if x < 0 else f'{x:0.1f}º' if x == 0 else f'{x:0.1f}ºE' for x in xticks]
    ax.set_xticklabels(xtickslabels, fontsize=fontsize)
    yticks = meridians[(meridians>=lat.min()) & (meridians<=lat.max())]
    ax.set_yticks(yticks)
    yticklabels = [f'{abs(y):0.1f}ºS' if y < 0 else f'{y:0.1f}º' if y == 0 else f'{y:0.1f}ºN' for y in yticks]
    ax.set_yticklabels(yticklabels, fontsize=fontsize)
    set_limits(ax, gpd_shp)
    
    if title:
        plt.title(title, fontsize=fontsize+4)
    
    if figname:
        plt.savefig(figname, dpi=300)
        plt.close()
    
    if show:
        plt.show()
        

def voronoi(data, shp, nbin=10, grid=.1, fontsize=11, title='', vmin=200, vmax=2000, 
            loc_legend='best', figname='', cmap='seismic', ticks=None, 
            legend_title='Legend', show=True):
    '''
    Gera mapa dos poligonos de Voronoi utilizados no método de Thiessen.

    Parameters
    ----------
    data : DataFrame
        Array de dados da variável de interesse.
    shp : shapefile
        Arquivo .shp da área de interesse.
    cmap : str ou array, optional
        Se str utiliza-se um dos mapas de cores padrões do matplotlib, se array
        define mapa de cor personalizado. The default is 'seismic'.
    nbin : int, optional
        Número de divisões/cores para legenda/mapa. The default is 10.
    grid : float, optional
        Tamanho da grade a ser utilizada no mapa. The default is .1.
    fontsize : int, optional
        Tamanho da fonte padrão do mapa. The default is 11.
    title : str, optional
        Título do mapa. The default is ''.
    vmin : float, optional
        Mínimo valor a ser utilizado para legenda. The default is 200.
    vmax : float, optional
        Máximo valor a ser utilizado para legenda. The default is 2000.
    loc_legend : str, optional
        Localização da legenda, valores padrões do matplolib. The default is 'best'.
    figname : str, optional
        Se fornecido salva figura no local e com nome informado. The default is ''.
    cmap : str ou array, optional
        Se str utiliza-se um dos mapas de cores padrões do matplotlib, se array
        define mapa de cor personalizado. The default is 'seismic'.
    ticks : array, optional
        Ticks label. The default is None.
    legend_title : str, optional
        Nome da legenda. The default is 'Legend'.
    show : bool, optional
        Se True exibe figura. The default is True.

    Returns
    -------
    None.

    '''
    gpd_shp = gpd.read_file(shp)
    gpd_geom = gpd_shp['geometry']
    env = gpd_geom.bounds
    geometry = [Point(xy) for xy in zip(data[0], data[1])]
    crs = 'epsg:4674'
    df = data.drop([0,1], axis=1)
    df = df.rename(columns={2:'var'})
    gdf = GeoDataFrame(df, crs=crs, geometry=geometry)
    scheme = mc.Quantiles(df['var'], k=nbin)
    quantiles = np.linspace(vmin, vmax, nbin)
    scheme.bins = quantiles
    
    if not ticks:
        ticks = np.linspace(vmin, vmax, nbin+1)
        ticks = list(ticks)
    
    if type(cmap) != str:
        cmap = mycmap(cmap, nbin) 
    
    fig, ax = plt.subplots(figsize=(15, 15))
    gpd_shp.plot(ax=ax, facecolor='none', edgecolor='black', lw=2)
    geoplot.voronoi(gdf, clip=gpd_geom, hue='var', scheme=scheme,
                         cmap=cmap, legend=True, ax=ax,
                         legend_kwargs={'fontsize':fontsize+4, 'loc': loc_legend,
                                        'title':legend_title, 'title_fontsize':fontsize+6})
    gpd_shp.plot(ax=ax, facecolor='none', edgecolor='black', lw=2)
    set_limits(ax, gpd_shp)
    
    #draw parallels and meridians
    parallels=np.arange(env.miny[0], env.maxy[0], grid)
    ax.hlines([env.miny[0], env.maxy[0]], env.minx[0], env.maxx[0], linewidth=1)
    ax.vlines([env.minx[0], env.maxx[0]], env.miny[0], env.maxy[0], linewidth=1)
    meridians = np.arange(env.minx[0], env.maxx[0], grid)
    ax.hlines(parallels, env.minx[0], env.maxx[0], linewidth=.3)
    ax.vlines(meridians, env.miny[0], env.maxy[0], linewidth=.3)
    xticks = meridians
    yticks = parallels
    yticklabels = [f'{y:0.1f}ºS' if y < 0 else f'{y:0.1}º' if y == 0 else 'f{y:0.1}ºN' for y in yticks]
    xtickslabels = [f'{x:0.1f}ºW' if x < 0 else f'{x:0.1f}º' if x == 0 else f'{x:0.1f}ºE' for x in xticks]
    
    for i in range(len(xticks)):
        ax.text(xticks[i]-0.03, env.miny[0]-0.02, xtickslabels[i], fontsize=fontsize*1.5)
    for i in range(len(yticks)):
        ax.text(env.maxx[0], yticks[i], yticklabels[i], fontsize=fontsize*1.5)
    
    plt.title(title, fontsize=fontsize*2)
    
    if figname:
        plt.savefig(figname, dpi=300)
    
    if show:
        plt.show()
    

def choropleth(data, shps, var='var', cmap='seismic', nbin=10, ticks=None, 
              vmin=200, vmax=2000, cbar_extend='both', cbar_label='Legend',
              cbar_loc='vertical', lbsize=10, lbpad=3, fontsize=21, grid=0.5,
              title='', figname='', show=True, ticks_label=[], figsize=(15,15),
              pos_ticks=0, legend='colorbar', norm=0):
    '''
    Gera mapas coropléticos.

    Parameters
    ----------
    data : array
        Dados da variável.
    shps : list
        Lista de shapefiles do mapa.
    var : str, optional
        Nome da variável de interesse. The default is 'var'.
    cmap : str ou array, optional
        Se str utiliza-se um dos mapas de cores padrões do matplotlib, se array
        define mapa de cor personalizado. The default is 'seismic'.
    nbin : int, optional
        Número de divisões/cores para legenda/mapa. The default is 10.
    ticks : array, optional
        Ticks label. The default is None.
    vmin : float, optional
        Mínimo valor a ser utilizado para legenda. The default is 200.
    vmax : float, optional
        Máximo valor a ser utilizado para legenda. The default is 2000.
    cbar_extend : str, optional
        Define se colorbar extender valores. The default is 'both'.
    cbar_label : str, optional
        Nome da legenda. The default is 'Legend'.
    cbar_loc : str, optional
        Localização da legenda (vertical, horizontal). The default is 'vertical'.
    lbsize : int, optional
        Tamanho das labels da legenda. The default is 10.
    lbpad : int, optional
        Label pad. The default is 3.
    fontsize : int, optional
        Tamanho padrão da fonte do mapa. The default is 21.
    grid : itn, optional
        Tamanho da grade a ser utilizada no mapa.. The default is 0.5.
    title : str, optional
        Título do mapa. The default is ''.
    figname : str, optional
        Se fornecido salva figura no local e com nome informado. The default is ''.
    show : boll, optional
        Se True exibe figura. The default is True.

    Returns
    -------
    None.

    '''
    
    fig, ax = plt.subplots(figsize=(15, 15))
    
    gpd_shps = []
    
    for i, shp in enumerate(shps):
        gpd_shp = gpd.read_file(shp)
        gpd_shp[var] = data[i]
        
        if i == 0:
            gpd_shps = gpd_shp
        else:
            gpd_shps = gpd_shps.append(gpd_shp)
    
    if not ticks:
        ticks = np.linspace(vmin, vmax, nbin+1)
        ticks = list(ticks)
    
    if type(cmap) != str:
        cmap = mycmap(cmap, nbin)  
    
    gpd_shps.plot(column=var, cmap=cmap, ax=ax)
    gpd_shps.plot(ax=ax, facecolor='none', edgecolor='black', lw=2)
    im = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.BoundaryNorm(ticks, cmap.N))
    
    __legend(data, im, legend, ticks, ticks_label, cbar_extend, vmin, vmax,
             nbin, cbar_loc, cbar_label, cmap, norm, lbsize, lbpad, fontsize)
    
    #draw parallels and meridians
    env = gpd_shps.bounds
    parallels=np.arange(env.miny.min()-grid, env.maxy.max()+grid, grid)
    meridians = np.arange(env.minx.min()-grid, env.maxx.max()+grid, grid)    
    ax.hlines(parallels, env.minx.min()-grid, env.maxx.max()+grid, linewidth=.3)
    ax.vlines(meridians, env.miny.min()-grid, env.maxy.max()+grid, linewidth=.3)
    xticks = meridians
    yticks = parallels
    yticklabels = [f'{y:0.1f}ºS' if y < 0 else f'{y:0.1}º' if y == 0 else 'f{y:0.1}ºN' for y in yticks]
    xtickslabels = [f'{x:0.1f}ºW' if x < 0 else f'{x:0.1f}º' if x == 0 else f'{x:0.1f}ºE' for x in xticks]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xtickslabels, fontsize=fontsize)
    ax.set_yticklabels(yticklabels, fontsize=fontsize)
    ax.set_xlim([env.minx.min(), env.maxx.max()])
    ax.set_ylim([env.miny.min(), env.maxy.max()])
    
    plt.title(title, fontsize=fontsize+6)
    
    if figname:
        plt.savefig(figname, dpi=300)
    
    if show:
        plt.show()
        
    
            