#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 21:30:26 2018

@author: duarte
"""

from netCDF4 import Dataset
from datetime import datetime as dt
from dateutil.relativedelta import relativedelta
import numpy as np
import geopandas as gpd


def __in_poly(shp, lat, lon, buffer=0):
    gpd_shp = gpd.read_file(shp)['geometry']
    gpd_shp = gpd_shp.buffer(buffer)
    bounds = gpd_shp.bounds
    pos_lon = np.where((lon >= bounds. minx[0]) & (lon <= bounds. maxx[0]))[0]
    pos_lat = np.where((lat >= bounds. miny[0]) & (lat <= bounds. maxy[0]))[0]
    buffer += buffer
    
    if buffer > 0:
        while len(pos_lon) == 0 or len(pos_lat) == 0:
            pos_lon, pos_lat = __in_poly(shp, lat, lon, buffer)

    return pos_lon, pos_lat


def __read_var(file, var_name, plat, plon, ptime):
    if len(Dataset(file).variables[var_name].shape) == 3:
        var = Dataset(file).variables[var_name][ptime, plat, plon]
    else:
        var = Dataset(file).variables[var_name][plat, plon]
        #var = Dataset(file).variables[var_name][ptime, :, :, plat, plon]
    return var


def __dates(time_units, times):
    init_date = time_units.split(' ')
    try:
        init_date = dt.strptime(time_units.split(' ')[-1], '%Y-%m-%d')
    except:
        try:
            init_date = dt.strptime(init_date[-2], '%Y-%m-%d')
        except:
            init_date = dt.strptime(init_date[-1].split('T')[0], '%Y-%m-%d')
        
    if time_units.split(' ')[0] == 'months':
        dates = np.array([init_date + relativedelta(months=x) for x in times])
    if time_units.split(' ')[0][:3].lower() == 'day':
        dates = np.array([init_date + relativedelta(days=x) for x in times])
    if time_units.split(' ')[0] == 'hours':
        dates = np.array([init_date + relativedelta(hours=x) for x in times])
    
    return dates


def __dates_in(dates, start_date, end_date):
    pos_time = np.where((dates>=start_date) & (dates<=end_date))[0]
    return pos_time, dates[pos_time]


def __get_lat_lon_time(file, variables):
    time_units, times = 0, 0
    for var in variables:
        if var.lower() in ['lat', 'latitude', 'y', 'nlat']:
            lat = Dataset(file).variables[var][:]
            
        if var.lower() in ['lon', 'longitude', 'x', 'nlon']:
            lon = Dataset(file).variables[var][:]
            lon[np.where(lon>180)] -= 360
        
        if var.lower() in ['time', 'gc_t', 'gc', 's']:
            time_units = Dataset(file).variables[var].units
            times = Dataset(file).variables[var][:]
            
    return lat, lon, time_units, times


def read_nc(file, shape=None, var_name='pre', freq='month', start_date=0,
            end_date=0, buffer=0):
    '''
    Parameters
    ----------
    file : string
        NetCDF filename.
    shape : string, optional
        Shapefile name.
    var_name : string, optional
        variable name. The default is 'pre'.
    freq : string, optional
        Time frequency. The default is 'month'.
    start_date : datetime, optional
        Initial date. The default is 0.
    end_date : datetime, optional
        Final date. The default is 0.
    buffer : float, optional
        Buffer value to be applied to shapefile. The default is 0.

    Returns
    -------
    var : np.array
        An np.array with the values of the selected variable.
    dates: np.array
        An np.array with the dates in the time interval provided.
    lat_in: np.array
        An np.array with the latitudes in the interval provided.
    lon_in: np.array
        An np.array with the longitudes in the interval provided.

    '''

    variables = [x for x in Dataset(file).variables]

    if var_name not in variables:
        print(f"{var_name} isn't in the {variables}.")
        return

    lat, lon, time_units, times = __get_lat_lon_time(file, variables)

    if time_units:
        dates = __dates(time_units, times)
    else:
        dates = [0]

    if start_date and end_date:
        pos_time_in, dates = __dates_in(dates, start_date, end_date)
    elif len(dates) > 1:
        pos_time_in = range(0, len(dates))
    else:
        pos_time_in = [0]

    if shape:
        pos_lon_in, pos_lat_in = __in_poly(shape, lat, lon, buffer)
        lat_in = lat[pos_lat_in]
        lon_in = lon[pos_lon_in]
    else:
        lat_in = lat
        lon_in = lon
        pos_lon_in = range(len(lon))
        pos_lat_in = range(len(lat))

    var = __read_var(file, var_name, pos_lat_in, pos_lon_in, pos_time_in)

    try:
        if freq == 'month' and 'day' in Dataset(file).variables[var_name].units:
            var *= 30
    except AttributeError:
        print(variables)
        if 'GC' in variables:
            var *= 30
        print(AttributeError)

    if time_units:
        return var, dates, lat_in, lon_in
    else:
        return np.array(var), lat_in, lon_in


def save_nc(var, lat, lon, dates, filename, time_units, var_name, var_units,
            var_long_name):
    output = Dataset(filename, 'w', format='NETCDF4')
    output.createDimension('time', len(dates))
    time = output.createVariable('time', 'f4', 'time')
    time.units = '{0}'.format(time_units)
    time.calendar = 'standard'
    time[:] = dates

    output.createDimension('lat', len(lat))
    lat_nc = output.createVariable('lat', 'f4', 'lat')
    lat_nc.units = 'degress-north'
    lat_nc.long_name = 'latitude'
    lat_nc.standard_name = 'latitude'
    lat_nc[:] = lat

    output.createDimension('lon', len(lon))
    lon_nc = output.createVariable('lon', 'f4', 'lon')
    lon_nc.units = 'degress_east'
    lon_nc.long_name = 'longitude'
    lon_nc.standard_name = 'longitude'
    lon_nc[:] = lon

    var_nc = output.createVariable(var_name, 'f4', ['time', 'lat', 'lon'])
    var_nc.units = var_units
    var_nc.long_name = var_long_name
    var_nc.Fillvalue = np.nan
    var_nc[:] = var

    output.close()
    print('Done ----->', filename)
