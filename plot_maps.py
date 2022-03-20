import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import cartopy, cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader 
import numpy as np 
import geopandas as gpd
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from shapely.geometry import Point
import mapclassify as mc

def get_bounds(shp):
    gpd_shp = gpd.read_file(shp)['geometry']
    bounds = gpd_shp.bounds.values[0]
    bounds = [bounds[0], bounds[2], bounds[1], bounds[3]]
    return bounds

def external_coords(geom):
    return list(geom.exterior.coords)

def clip_area(shp, ax):
    sf = gpd.read_file(shp)
    vertices = sf.geometry.apply(external_coords)[0]

    codes = [] 
    codes += [Path.MOVETO]
    codes += [Path.LINETO] * (len(vertices) -2)
    codes += [Path.CLOSEPOLY]
    clip = Path(vertices, codes)
    clip = PathPatch(clip, transform=ax.transData)
    return clip

def mycmap(colors, extend):

    color_under = colors[0]
    color_over = colors[-1]

    if extend == 'max':
        colors = colors[:-1]
    elif extend == 'min':
        colors = colors[1:]
    elif extend == 'both':
        colors = colors[1:-1]

    cmap = matplotlib.colors.ListedColormap(colors, N=len(colors))
    cmap.set_over(color_over)
    cmap.set_under(color_under)

    return cmap

def imshow(data, lon, lat, bounds=[], shape=False, decorators=[], label='',
           vmin=np.nan, vmax=np.nan, cmap='jet', extend='both', 
           orientation='vertical', pad=0.05, fraction=0.05, nbin=10, ticks=[],
           title='', fontsize=10, fontweight='bold', loc='center', 
           figsize=(10,10), filename=False, show=True, edgecolorshp='black', 
           facecolorshp='none', linewidthshp=.5, drawcoast=False, clip=False,
           gridlat=5, gridlon=5 ):

    """Genenerate maps using imshow function of Cartopy

    Args:
        data (2D array): Matrix with the dada do plot
        lon (1D array): Vector with the longitude values
        lat (1D array): Vector with the latitude values
        bounds (list, optional): Plot area limits. Defaults to False.
        shape (str, optional): Shapefile that demilit the plot areas. Defaults to False.
        decorators (list): A list (strings) of shapefiles to decorate the plot.
        label (str, optional): Label of the legend. Defaults to ''.
        vmin (float, optional): Minimum value to scale. Defaults to np.nan.
        vmax (float, optional): Maximum value to scale. Defaults to np.nan.
        cmap (str or list, optional): Colormap to plot str os list of colors. Defaults to 'jet'.
        extend (str, optional): Extend colorbar argument: 'both', 'max', 'min', 'none'. Defaults to 'both'.
        orientation (str, optional): Colorbar orientation. Defaults to 'vertical'.
        pad (float, optional): Distance between the plot area and the colorbar. Defaults to 0.05.
        fraction (float, optional): Colorbar size relative to the plot area. Defaults to 0.05.
        nbin (int, optional): Segments number of the colorbar. Defaults to 10.
        ticks (list, optional): Colorbar ticks. Defaults to [].
        title (str, optional): Image title. Defaults to ''.
        fontsize (int, optional): Font size. Defaults to 10.
        fontweight (str, optional): Font weight. Defaults to 'bold'.
        loc (str, optional): Colorbar position. Defaults to 'center'.
        figsize (tuple, optional): Figure size. Defaults to (10,10).
        filename (bool, optional): Name of the file to save the plot. Defaults to False.
        show (bool, optional): If True display the plot in the workscape. Defaults to True.
        edgecolorshp (str, optional): Shape's edge color. Defaults to 'black'.
        facecolorshp (str, optional): Shape's face color. Defaults to 'none'.
        linewidthshp (float, optional): Width of the lat-lon grid. Defaults to .5.
        drawcoast (bool, optional): If True draw the coast lines. Defaults to False.
        gridlon (int, optional): resolution of the longitude grid. Defaults to 5.
        gridlat (int, optional): Resolution of the latitude grid. Defaults to 5.
        clip (bool, optional): If True shows just the area inside the shapefile. Defaults to False.
        
        
    """
    
    
    plt.figure(figsize=figsize)                  # Choose the plot size (width x height, in inches)
    ax = plt.axes(projection=ccrs.PlateCarree()) # Use the Cilindrical Equidistant projection in cartopy
    
    xs, ys = np.meshgrid(lon, lat)
    
    if not len(bounds):                               # Define the image bounds [min.lon, max.lon, min.lat, max.lat]
        if shape:
            bounds = get_bounds(shape)
        else:
            bounds = [min(lon), max(lon), min(lat), max(lat)]
    
    ax.set_extent(bounds, crs=ccrs.PlateCarree())
    
    if clip:
        clip = clip_area(shape, ax)
    
    if type(cmap) != str:
        cmap = mycmap(cmap, extend)

    if np.isnan(vmin):
        vmin = data.min()
    if np.isnan(vmax):
        vmax = data.max()
    
    if not np.array(ticks).any():
        ticks = np.linspace(vmin, vmax, nbin)
        print(ticks)
    
    if shape:
        shapefile = list(shpreader.Reader(shape).geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=edgecolorshp,
                          facecolor=facecolorshp, linewidth=linewidthshp)
    
    if len(decorators):
        for decorator in decorators:
            decorator = list(shpreader.Reader(decorator).geometries())
            ax.add_geometries(decorator, ccrs.PlateCarree(), edgecolor=edgecolorshp,
                          facecolor=facecolorshp, linewidth=linewidthshp)
    
    if drawcoast:
        ax.add_feature(cartopy.feature.COASTLINE, color='black', linewidth=0.8)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, 
                      linestyle='--', linewidth=0.25, 
                      xlocs=np.arange(-180, 180, gridlon),
                      ylocs=np.arange(-90, 90, gridlat),
                      draw_labels=True)
    
    gl.top_labels = False
    gl.right_labels = False
    
    # Plot the image
    if clip:
        img = ax.pcolormesh(xs, ys, data, vmin=vmin, vmax=vmax, cmap=cmap, 
                            transform=ccrs.PlateCarree(), clip_path=clip)
    else:
        img = ax.pcolormesh(xs, ys, data, vmin=vmin, vmax=vmax, cmap=cmap,
                            transform=ccrs.PlateCarree())
    
    # Add a colorbar
    plt.colorbar(img, label=label, extend=extend, orientation=orientation,
                 pad=pad, fraction=fraction, ticks=ticks, format='%.1f')
    
    # Add a title
    plt.title(title , fontweight=fontweight, fontsize=fontsize, loc=loc)

    # Save the image
    if filename:
        plt.savefig(filename)
    
    # Show the image
    if show:
        plt.show()

def voronoi(lon, lat, data, shp, nbin=10, grid=.1, fontsize=11, title='', 
            vmin=200, vmax=2000, loc_legend='best', figname='', cmap='seismic',
            ticks=None, legend_title='Legend', show=True):
    """_summary_

    Args:
        lon (_type_): _description_
        lat (_type_): _description_
        data (_type_): _description_
        shp (_type_): _description_
        nbin (int, optional): _description_. Defaults to 10.
        grid (float, optional): _description_. Defaults to .1.
        fontsize (int, optional): _description_. Defaults to 11.
        title (str, optional): _description_. Defaults to ''.
        vmin (int, optional): _description_. Defaults to 200.
        vmax (int, optional): _description_. Defaults to 2000.
        loc_legend (str, optional): _description_. Defaults to 'best'.
        figname (str, optional): _description_. Defaults to ''.
        cmap (str, optional): _description_. Defaults to 'seismic'.
        ticks (_type_, optional): _description_. Defaults to None.
        legend_title (str, optional): _description_. Defaults to 'Legend'.
        show (bool, optional): _description_. Defaults to True.
    """
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
    geometry = [Point(xy) for xy in zip(lon, lat)]
    crs = 'epsg:4674'
    #df = data.drop([0,1], axis=1)
    df = data
    df = df.rename(columns={1:'var'})
    gdf = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
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
    # geoplot.voronoi(gdf, clip=gpd_geom, hue='var', scheme=scheme,
    #                      cmap=cmap, legend=True, ax=ax,
    #                      legend_kwargs={'fontsize':fontsize+4, 'loc': loc_legend,
    #                                     'title':legend_title, 'title_fontsize':fontsize+6})
    gpd_shp.plot(ax=ax, facecolor='none', edgecolor='black', lw=2)
    set_limits(ax, gpd_shp)
    
    #draw parallels and meridians
    parallels=np.arange(env.miny[0]+.1, env.maxy[0], grid)
    # ax.hlines([env.miny[0], env.maxy[0]], env.minx[0], env.maxx[0], linewidth=1)
    # ax.vlines([env.minx[0], env.maxx[0]], env.miny[0], env.maxy[0], linewidth=1)
    meridians = np.arange(env.minx[0]+.5, env.maxx[0], grid)
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
              pos_ticks=0, legend='colorbar', norm=0, fig=None, ax=None,
              yticklabels=True, xtickslabels=True):
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
    if not fig and not ax:
        fig, ax = plt.subplots(figsize=(15, 15))
#%%
    gpd_shps = []
    
    for i, shp in enumerate(shps):
        gpd_shp = gpd.read_file(shp)
        gpd_shp[var] = data[i]
        
        if i == 0:
            gpd_shps = gpd_shp
        else:
            gpd_shps = gpd_shps.append(gpd_shp)
#%% 
    if not np.array(ticks).any():
        ticks = np.linspace(vmin, vmax, nbin+1)
        ticks = list(ticks)
    
    if type(cmap) != str:
        cmap = mycmap(cmap, nbin)  

#%%  
    
    gpd_shps.plot(column=var, cmap=cmap, vmin=vmin, vmax=vmax, ax=ax)
    gpd_shps.plot(ax=ax, facecolor='none', edgecolor='black', lw=2)
    im = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.BoundaryNorm(ticks, cmap.N))

    if legend:
        __legend(data, im, legend, ticks, ticks_label, cbar_extend, vmin, vmax,
                nbin, cbar_loc, cbar_label, cmap, norm, lbsize, lbpad, fontsize, ax)
#%%
    #draw parallels and meridians
    env = gpd_shps.bounds
    parallels=np.arange(env.miny.min()-grid, env.maxy.max()+grid, grid)
    ax.hlines(parallels, env.minx.min()-grid, env.maxx.max()+grid, linewidth=.3)
    yticks = parallels

    if yticklabels:
        ax.set_yticks(yticks)
        yticklabels = [f'{y:0.1f}ºS' if y < 0 else f'{y:0.1}º' if y == 0 else 'f{y:0.1}ºN' for y in yticks]
        ax.set_yticklabels(yticklabels, fontsize=fontsize)
    else:
        ax.set_yticklabels([])
        ax.yaxis.set_ticks_position('none') 
    
    meridians = np.arange(env.minx.min()-grid, env.maxx.max()+grid, grid)
    ax.vlines(meridians, env.miny.min()-grid, env.maxy.max()+grid, linewidth=.3)
    xticks = meridians
    
    if xtickslabels:
        ax.set_xticks(xticks)
        xtickslabels = [f'{x:0.1f}ºW' if x < 0 else f'{x:0.1f}º' if x == 0 else f'{x:0.1f}ºE' for x in xticks]
        ax.set_xticklabels(xtickslabels, fontsize=fontsize)
    else:
        ax.set_xticklabels([])
        ax.xaxis.set_ticks_position('none') 

    ax.set_xlim([env.minx.min(), env.maxx.max()])
    ax.set_ylim([env.miny.min(), env.maxy.max()])
    
    plt.title(title, fontsize=fontsize+6)
    
    if figname:
        plt.savefig(figname, dpi=300)
    
    if show:
        plt.show()
        

def quiver(u, v, lat, lon, m, shp, clip=None, vmin=0, vmax=100, cmap='seismic_r',
              nbin=10, cbar_loc='vertical', cbar_label='My legend', grid=1,
              cbar_extend='both', lbsize=10, lbpad=3, fontsize=11, aux_shapes=[],
              title='', figname='', show=False, ticks=0, norm=None, figsize=(15,15),
              ticks_label=[], pos_ticks=0, legend='colorbar'):
    
    gpd_shp = gpd.read_file(shp)
    fig, ax = plt.subplots(figsize=figsize)

    xx,yy = np.meshgrid(lat, lon)
    
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
    
    im = plt.pcolormesh(yy, xx, m.T, cmap=cmap, clip_path=clip, norm=norm)

    q = plt.quiver(lon, lat, u, v,
                   scale = 1/0.03)
    qk = ax.quiverkey(q, 0.7, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',
                   coordinates='figure', fontproperties={'size':24})
    
                                    
    __legend(m, im, legend, ticks, ticks_label, cbar_extend, vmin, vmax,
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

def __legend(data, im, type_legend, ticks, ticks_label, cbar_extend, vmin, 
             vmax, nbin, cbar_loc, cbar_label, cmap, norm, lbsize, lbpad,
             fontsize, ax):
    
    if type_legend == 'colorbar':
        __colorbar_legend(im, ticks, ticks_label, cbar_extend, vmin, vmax, nbin,
                        cbar_loc, cbar_label, cmap, norm, lbsize, lbpad, fontsize,
                        ax)
    if type_legend == 'discret':
        __discret_legend(data, im, fontsize, ticks_label, ticks, ax)

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


def __colorbar_legend(im, ticks, ticks_label, cbar_extend, vmin, vmax, nbin,
                      cbar_loc, cbar_label, cmap, norm, lbsize, lbpad, 
                      fontsize, ax):

    if not np.array(ticks).any():
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
    
    axins1 = inset_axes(ax,
                    width="5%",  # width = 5% of parent_bbox width
                    height="95%",  # height : 95%
                    bbox_to_anchor=(.1, 0., 1, 1),
                    bbox_transform=ax.transAxes,
                    loc='center right')
    
    cb = plt.colorbar(im, orientation=cbar_loc, 
              ticks=ticks2,
              extendfrac = 'auto',
              label=cbar_label, extend=cbar_extend,
              boundaries=ticks, fraction=0.038,
              cmap=cmap, norm=norm, cax=axins1
              )
    cb.ax.tick_params(labelsize=lbsize)
    
    if cbar_loc == 'vertical':
        cb.ax.set_yticklabels(ticks_label)
    else:
        cb.ax.set_xticklabels(ticks_label)
            
    cb.set_label(cbar_label, labelpad=lbpad, size=fontsize)

def set_limits(ax, shp):
    env = shp.geometry.bounds
    minx = env.minx.min()
    maxx = env.maxx.max()
    miny = env.miny.min()
    maxy = env.maxy.max()
    ax.set_xlim([minx-.1, maxx+.1])
    ax.set_ylim([miny-.1, maxy+.1])