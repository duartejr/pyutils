import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
from geovoronoi import voronoi_regions_from_coords, calculate_polygon_areas
import numpy as np
from calendar import monthrange
from shapely.ops import cascaded_union
from geovoronoi import points_to_coords


def voronoi(coords, shp):
    poly_shp, _ = voronoi_regions_from_coords(coords, shp)
    voronoi_areas = calculate_polygon_areas(poly_shp)
    return np.array(list(voronoi_areas.values()))

def geo_coords(coords):
    df1 = pd.DataFrame({'X':coords[:,1], 'Y':coords[:,0]})
    df1['geometry'] = list(zip(df1['X'], df1['Y']))
    df1['geometry'] = df1['geometry'].apply(Point)
    gdf1 = gpd.GeoDataFrame(df1, geometry='geometry')
    return gdf1


def __monthly(start, end, files, shp):
    pr_mean = []
    print('Calculating Thiessen for:')
    for year in range(start, end+1):
        for month in range(1, 13):
            print(f'{month}/{year}')
            pr_month = pd.DataFrame(columns=['Latitude', 'Longitude',
                                                  'Anos', 'Meses', 'Total'])
            for pr_file in files:
                pr = pd.read_table(pr_file,
                                    usecols=['Latitude', 'Longitude', 'Anos',
                                            'Meses', 'Total'],
                                    sep=';').dropna()
                pr = pr[(pr.Anos == year) & (pr.Meses == month)]
                pr_month = pr_month.append(pr)
            
                coords = pr_month[['Latitude', 'Longitude']].values
            if len(coords) > 2:
                voronoi_areas =voronoi(coords, shp)
                pr_mean.append([year, month, sum(pr_month['Total'].values * voronoi_areas)/
                                voronoi_areas.sum()])
    return pr_mean


def __daily(start, end, files, shp):
    pr_mean = []
    boundary = gpd.read_file(shp)
    boundary_shape = cascaded_union(boundary.geometry)
    
    print('Calculating Thiessen for:')
    for year in range(start, end+1):
        print(year)
        for month in range(1, 13):
            print(month)
            pr_month = pd.DataFrame([])

            for pr_file in files:
                pr = pd.read_csv(pr_file)
                pr = pr[(pr.Anos == year) & (pr.Meses == month)]
                pr_month = pr_month.append(pr)
            
            for day in range(1, monthrange(year, month)[1]+1):
                print(f'{day}/{month}/{year}')
                pr_day = pr_month.loc[pr_month['Dia'] == day]['pr'].values
                coords = pr_month.loc[pr_month['Dia'] == day][['Latitude', 'Longitude']].values
                not_nan = np.where(~np.isnan(pr_day))
                pr_day = pr_day[not_nan]
                coords = coords[not_nan]
                     
                if len(coords) > 2:
                    coords = geo_coords(coords).set_crs(boundary.crs)
                    coords = points_to_coords(coords.geometry)
                    voronoi_areas = voronoi(coords, boundary_shape)
                    pr_th = sum(pr_day * voronoi_areas) / voronoi_areas.sum()
                    pr_mean.append([year, month, day, pr_th])

    return pr_mean


def thiessen(shp, pr_files, start_year, end_year, min_obs=360, buffer=0,
             calc_monthly=True):

    ##################################################
    # Realiza cálculos de thiseen com postos selecionados #
    if calc_monthly:
        pr_mean = __monthly(start_year, end_year, pr_files, shp)
    else:
        pr_mean = __daily(start_year, end_year, pr_files, shp)
    
    ########################################################
    
    return np.array(pr_mean)

