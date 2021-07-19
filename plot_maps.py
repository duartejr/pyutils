import matplotlib
import matplotlib.pyplot as plt # Plotting library
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import cartopy, cartopy.crs as ccrs # Plot maps
import cartopy.io.shapereader as shpreader # Import shapefiles
import numpy as np # Scientific computing with Python
import geopandas as gpd

def get_extent(shp):
    gpd_shp = gpd.read_file(shp)['geometry']
    bounds = gpd_shp.bounds
    return bounds.values[0]

def clip_area(shp, ax):
    
    sf = gpd.read_file(shp)

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

def mycmap(colors):
    cmap = matplotlib.colors.ListedColormap(colors[1:-1])
    cmap.set_over(colors[0])
    cmap.set_under(colors[-1])
    return cmap

def imshow(data, shape=False, figsize=(10,10), extent=False, vmin=0, vmax=100,
            cmap='jet', label='', extend='both', orientation='vertical',
            pad=0.05, fraction=0.05, title='', fontweight='bold',
            fontsize=10, loc='center', filename=False, show=True, gridlat=5,
            gridlon=5, edgecolorshp='black', facecolorshp='none',
            linewidthshp=.5, drawcoast=False, clip=False, ticks=[], nbin=10):
    # Choose the plot size (width x height, in inches)
    plt.figure(figsize=figsize)
    # Use the Cilindrical Equidistant projection in cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    if extent:
        # Define the image extent [min. lon, max. lon, min. lat, max. lat]
        img_extent = [extent[0], extent[2], extent[1], extent[3]]
    if shape and not extent:
        extent = get_extent(shape)
        img_extent = [extent[0], extent[2], extent[1], extent[3]]
    
    if clip:
        clip = clip_area(shape, ax)
    
    if type(cmap) != str:
        cmap = mycmap(cmap)
    
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
    img = ax.imshow(data, origin='lower', extent=img_extent, vmin=vmin,
                    vmax=vmax, cmap=cmap, clip_path=clip)
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