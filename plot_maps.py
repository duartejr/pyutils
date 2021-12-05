import matplotlib
import matplotlib.pyplot as plt # Plotting library
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import cartopy, cartopy.crs as ccrs # Plot maps
import cartopy.io.shapereader as shpreader # Import shapefiles
import numpy as np # Scientific computing with Python
import geopandas as gpd
from scipy import ndimage

def get_extent(shp):
    gpd_shp = gpd.read_file(shp)['geometry']
    bounds = gpd_shp.bounds
    return bounds.values[0]

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

def mycmap(colors):
    cmap = matplotlib.colors.ListedColormap(colors[1:-1])
    cmap.set_over(colors[0])
    cmap.set_under(colors[-1])
    return cmap


def imshow(data, shape='', figsize=(10,10), extent=[], 
            vmin=np.nan, vmax=np.nan,
            cmap='jet', label='', extend='both', orientation='vertical',
            pad=0.05, fraction=0.05, title='', fontweight='bold',
            fontsize=10, loc='center', filename=False, show=True, gridlat=5,
            gridlon=5, edgecolorshp='black', facecolorshp='none',
            linewidthshp=.5, drawcoast=False, clip=False, ticks=[], nbin=10):
    """Genenerate maps using imshow function of Cartopy

    Args:
        data (array): Matrix with the dada do plot
        shape (str, optional): Shapefile that demilit the plot areas. Defaults to False.
        figsize (tuple, optional): Figure size. Defaults to (10,10).
        extent (list, optional): Plot area limits. Defaults to False.
        vmin (float, optional): Minimum value to scale. Defaults to np.nan.
        vmax (float, optional): Maximum value to scale. Defaults to np.nan.
        cmap (str or list, optional): Colormap to plot str os list of colors. Defaults to 'jet'.
        label (str, optional): Label to legend. Defaults to ''.
        extend (str, optional): Extend colorbar argument: 'both', 'max', 'min', 'none'. Defaults to 'both'.
        orientation (str, optional): Colorbar orientation. Defaults to 'vertical'.
        pad (float, optional): Distance between the plot area and the colorbar. Defaults to 0.05.
        fraction (float, optional): Colorbar size relative to the plot area. Defaults to 0.05.
        title (str, optional): Image title. Defaults to ''.
        fontweight (str, optional): Font weight. Defaults to 'bold'.
        fontsize (int, optional): Font size. Defaults to 10.
        loc (str, optional): Colorbar position. Defaults to 'center'.
        filename (bool, optional): Name of the file to save the plot. Defaults to False.
        show (bool, optional): If True display the plot in the workscape. Defaults to True.
        gridlat (int, optional): Resolution of the latitude grid. Defaults to 5.
        gridlon (int, optional): resolution of the longitude grid. Defaults to 5.
        edgecolorshp (str, optional): Shape's edge color. Defaults to 'black'.
        facecolorshp (str, optional): Shape's face color. Defaults to 'none'.
        linewidthshp (float, optional): Width of the lat-lon grid. Defaults to .5.
        drawcoast (bool, optional): If True draw the coast lines. Defaults to False.
        clip (bool, optional): If True shows just the area inside the shapefile. Defaults to False.
        ticks (list, optional): Colorbar ticks. Defaults to [].
        nbin (int, optional): Segments number of the colorbar. Defaults to 10.
    """
    
    # Choose the plot size (width x height, in inches)
    plt.figure(figsize=figsize)
    # Use the Cilindrical Equidistant projection in cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    if extent:
        # Define the image extent [min. lon, max. lon, min. lat, max. lat]
        img_extent = [extent[0], extent[2], extent[1], extent[3]]
    else:
        extent = get_extent(shape)
        img_extent = [extent[0], extent[2], extent[1], extent[3]]
    
    if clip:
        clip = clip_area(shape, ax)
    
    if type(cmap) != str:
        cmap = mycmap(cmap)

    if np.isnan(vmin):
        vmin = data.min()
    if np.isnan(vmax):
        vmax = data.max()
    
    if not np.array(ticks).any():
        ticks = np.linspace(vmin, vmax, nbin)
    
    shapefile = list(shpreader.Reader(shape).geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=edgecolorshp,
                      facecolor=facecolorshp, linewidth=linewidthshp)
    
    if drawcoast:
        ax.coastlines(resolution='10m', color='black', linewidth=0.8)

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
        img = ax.imshow(data, origin='lower', extent=img_extent, vmin=vmin,
                        vmax=vmax, cmap=cmap, clip_path=clip)
    else:
        img = ax.imshow(data, origin='lower', extent=img_extent, vmin=vmin,
                        vmax=vmax, cmap=cmap)
    # Add a colorbar
    plt.colorbar(img, label=label, extend=extend, orientation=orientation,
                 pad=pad, fraction=fraction, ticks=ticks)
    # Add a title
    plt.title(title , fontweight='bold', fontsize=10, loc='center')

    # Save the image
    if filename:
        plt.savefig(filename)
    # Show the image
    if show:
        plt.show()
