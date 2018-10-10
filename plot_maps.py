import ogr
from mpl_toolkits.basemap import Basemap
import shapefile
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from pyutils.custom_colorbar import colorbar as cbar


def bbox(shape):
	src_ds = ogr.Open(shape)
	src_lyr = src_ds.GetLayer()
	src_lyr.ResetReading()
	
	for feat in src_lyr:
		# get bounding coords in minx, maxx, miny, maxy format
		env = feat.GetGeometryRef().GetEnvelope()

	return env

def clip_area(sf, ax):
    for shape_rec in sf.shapeRecords():
        vertices = []
        codes = []
        pts = shape_rec.shape.points
        prt = list(shape_rec.shape.parts) + [len(pts)]
        for i in range(len(prt) - 1):
            for j in range(prt[i], prt[i+1]):
                vertices.append((pts[j][0], pts[j][1]))
            codes += [Path.MOVETO]
            codes += [Path.LINETO] * (prt[i+1] - prt[i] -2)
            codes += [Path.CLOSEPOLY]
        clip = Path(vertices, codes)
        clip = PathPatch(clip, transform=ax.transData)
    return clip

def mycmap(cmap, nbin):
    return LinearSegmentedColormap.from_list('mycmap', cmap, N=nbin)

def plot_maps(data, lat, lon, dir_shape, name_shape, clip=True, vmin=0,
              vmax=100, cmap='seismic_r', nbin=10, cbar_loc='vertical',
              cbar_label='My legend', cbar_extend='both', lbsize=10, lbpad=3,
              fontsize=11, aux_shapes=[], parallels=np.arange(-90., 90., 2),
              meridians = np.arange(180., 360., 2), title='', figname='',
              show=False):
    '''
    data(array): Matriz com dados a serem plotados. Shape [lat, lon]
    lat(array): Vetor com coordenadas de latitude
    lon(array): Vetor com coordenadas de longitude
    dir_shape(str): Nome do diretório que contêm shape a ser usado como máscara
    name_shape(str): Nome do shapefile a ser usado como máscara
    clip(voolean): Se True aplica máscara a área de desenho
    vmin(float): Mínimo valor para a colorbar
    vmax(float): Máximo valor para a colorbar
    cmap(str or array): Se str aplica um dos colormaps padrões do matplotlib. Se
                        array cria colormap personalizado.
    nbin(int): Número de divisisões para colormap personalizado.
    cbar_loc(str): Localizaçao da colorbar 'vertical' ou 'horizontal'
    cbar_label(str): Legenda para a colorbar
    cbar_extend(str): Extende limites da colorbar: 'both', 'max', 'min' ou None
    lbsize(int): Tamanho da colorbar
    fontsize(int): Tamanho da fonte da colorbar.
    aux_shapes(list): Lista de shapes para personalizar área de desenho.
    parallels(array): Vetor com paralelos a serem marcados no desenho.
    meridians(array): Vetor com meridianos a serem marcados no desenho.
    title(str): Título do desenho.
    figname(str): Nome para salvar a figura.
    show(boolean): Se True mostra prévia do desenho.
    '''    
    env = bbox(dir_shape)
    
    map = Basemap(llcrnrlat=env[2],
                  llcrnrlon=env[0],
                  urcrnrlat=env[3],
                  urcrnrlon=env[1])
    
    sf = shapefile.Reader(dir_shape+'/'+name_shape)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if clip:
        clip = clip_area(sf, ax)
    
    xx,yy = np.meshgrid(lon, lat)
        
    if type(cmap) != str:
        cmap = mycmap(cmap, nbin)
    
    if clip:
        map.pcolormesh(xx, yy, data, cmap=cmap, clip_path=clip, vmin=vmin,
                            vmax=vmax)
    else:
        map.pcolormesh(xx, yy, data, cmap=cmap, vmin=vmin, vmax=vmax)
            
    ticks = np.arange(vmin, vmax+nbin, (vmax-vmin)/nbin)
    ticks_label = ['{0:.2f}'.format(x) for x in ticks]
    cb = plt.colorbar(orientation=cbar_loc, 
                  ticks= ticks,
                  label=cbar_label, extend=cbar_extend)
    cb.ax.tick_params(labelsize=lbsize)
    if cbar_loc == 'vertical':
        cb.ax.set_yticklabels(ticks_label)
    else:
        cb.ax.set_xticklabels(ticks_label)
            
    cb.set_label(cbar_label, labelpad=lbpad, size=fontsize)
    
    if len(aux_shapes):
        for aux_shape in aux_shapes:
            map.readshapefile(aux_shape, 'myshape', drawbounds=True,
                              linewidth=.6, zorder=2, color='k')
            
    map.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10,
                  linewidth=.3)
    
    map.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10,
                      linewidth=.3)
    plt.title(title, fontsize=13)
    
    if figname:
        plt.savefig(figname)
        plt.close()
    
    if show:
        plt.show()
        plt.close()

def plot_maps_shp(shapes, fill_value=0, bbox=[-75, -30, -34, 6],
                  colors = [(1,0,0), (1,1,1), (0,0,1)],
                  nbins = 10,
                  cbar_loc = 'vertical',
                  vmin=0, vmax=9, ticks_label=False, cbar_extend=None,
                  cbar_label=False,
                  label=False,
                   parallels=np.arange(-90., 90., 2),
                   meridians = np.arange(180., 360., 2), title='',
                   figname=False):
   
    map = Basemap(llcrnrlat=bbox[2],
                  llcrnrlon=bbox[0],
                  urcrnrlat=bbox[3],
                  urcrnrlon=bbox[1])
    
    if type(colors) != str:
        cmap = mycmap(colors, nbins)
        
    axes = [0.07, 0.06, 0.77, 0.84]
    fig = plt.figure(figsize=(6, 5.7))
    ax = fig.add_axes(axes)
    bounds = np.linspace(vmin, vmax, nbins+1)
    
    colors =  matplotlib.cm.get_cmap(cmap)
    
    for c, shape in enumerate(shapes):
        map.readshapefile(shape, 'myshape', drawbounds=True,
                              linewidth=.6, zorder=2, color='k')
        patches = []
        
        for shape in map.myshape:
            patches.append(Polygon(np.array(shape), True))
            
        for d, v in enumerate(bounds[1:]):
            if fill_value[c] <= v:
                n = d
                break
        
        ax.add_collection(PatchCollection(patches,
                                          facecolor=colors(n),
                                          edgecolor='k', linewidth=.25))
    plt.title(title, fontsize=13)
    
    map.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10,
                  linewidth=.3)
    map.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10,
                      linewidth=.3)
    axes = [0.85, 0.17, 0.03, 0.61]
    cbar(cmap, axes, np.linspace(vmin, vmax, nbins+1), label, fig, 'vertical',
         ticks_label=cbar_label)
   
    if figname:
        plt.savefig(figname)
        plt.close()
    else:
        plt.show()
    