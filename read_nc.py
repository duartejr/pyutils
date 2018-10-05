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
import shapefile as shp
from pyutils.isinpoly import isinpoly3
from pyutils.buffer import buffer


def read_nc(file, shape, var_name='pre', var_time='month',
            start_date=0, end_date=0, only_in_shape=False,
            buffdist=0):
    variables = Dataset(file).variables
    for v in variables:
        print(v)
    
    time_units = 0
    
    for var in variables:
        if var.lower() == 'lat' or var.lower() == 'latitude':
            lat = Dataset(file).variables[var][:]
        elif var.lower() == 'lon' or var.lower() == 'longitude':
            lon = Dataset(file).variables[var][:]
        elif var.lower() == 'time' or var == 'GC_T' or var == 'GC':
            time_units = Dataset(file).variables[var].units
            times = Dataset(file).variables[var][:]

    pos_lon = np.where(lon>180)[0]
    lon[pos_lon] = lon[pos_lon]-360

    if time_units:
        init_date = time_units.split(' ')
        try:
            init_date = dt.strptime(time_units.split(' ')[-1], '%Y-%m-%d')
        except:
            init_date = dt.strptime(init_date[-2], '%Y-%m-%d')
   
        if time_units.split(' ')[0] == 'months':
            dates = np.array([init_date + relativedelta(months=x) for x in times])
        if time_units.split(' ')[0][:3].lower() == 'day':
            dates = np.array([init_date + relativedelta(days=x) for x in times])

    else:
        init_date = 0
        dates = [0]

    if start_date and end_date:
        pos_time = np.where((dates>=start_date) & (dates<=end_date))[0]
    else:
        pos_time = range(0,len(dates))

    if buffdist:
        new_shape = '/tmp/myshape'
        shape = '/'.join(shape.split('/')[:-1])
        buffer(shape, new_shape, buffdist)
        shape = new_shape+'/'+shape.split('/')[-1]
            
    shape = shp.Reader(shape)
    shape = shape.shapes()
    lat_resol = lat[1] - lat[0]
    lon_resol = lon[1] - lon[0]

    if only_in_shape:
        polig = []
        h = len(lon)*len(lat)
        Lon, Lat = np.meshgrid(lon,lat)
        Lon = np.reshape(Lon, h)
        Lat = np.reshape(Lat, h)

        for point in shape[0].points:
            polig.append(list(point))

        polig = np.array(polig)
        isin = isinpoly3(Lon, Lat, polig)
        isin = np.nonzero(isin)[0]
        Lons_in = Lon[isin]
        Lats_in = Lat[isin]
        lon_in = np.isin(lon, Lons_in)
        lat_in = np.isin(lat, Lats_in)
        
        if time_units:
            try:
                pr = Dataset(file).variables[var_name][pos_time,lat_in,lon_in]
            except:
                pr = Dataset(file).variables[var_name.lower()][pos_time,lat_in,lon_in]
        else:
            try:
                pr = Dataset(file).variables[var_name][lat_in,lon_in]
            except:
                pr = Dataset(file).variables[var_name.lower()][lat_in,lon_in]

        par = [par for par in zip(Lats_in, Lons_in)]
        
        for i,la in enumerate(lat[lat_in]):
            for j,lo in enumerate(lon[lon_in]):
                if (la,lo) not in par:
                    try:
                        pr[:,i,j] = np.nan
                    except:
                        pr[i,j] = np.nan
        
    else:
        bbox = shape[0].bbox
        
        
        if bbox[1] < 0:
            lat_in = np.where((lat>=(bbox[1]-lat_resol)) & (lat<=bbox[-1]+lat_resol))[0]
        else:
            lat_in = np.where((lat>=bbox[1]-2*lat_resol) & (lat<=bbox[-1]+2*lat_resol))[0]
        if bbox[0] < 0:
            lon_in = np.where((lon>=bbox[0]-2*lon_resol) & (lon<=bbox[2]+2*lon_resol))[0]
        else:
            lon_in = np.where((lon>=bbox[0]+2*lon_resol) & (lon<=bbox[2]+2*lon_resol))[0]
                
        try:
            pr = Dataset(file).variables[var_name][pos_time,lat_in,lon_in]
        except:
            pr = Dataset(file).variables[var_name.lower()][pos_time,lat_in,lon_in]
    
    if time_units:
        return pr, dates[pos_time], lon[lon_in], lat[lat_in]
    else:
        return np.array(pr), np.array(lon)[lon_in], np.array(lat)[lat_in]

def save_nc(var,lat,lon,dates,filename,time_units, var_name, var_units, var_long_name):
    output = Dataset(filename, 'w', format='NETCDF4')
    outtime = np.arange(len(dates))
    output.createDimension('time', len(outtime))
    time = output.createVariable('time', 'f4', 'time')
    init_date = dt.strftime(dates[0], '%Y-%m-%d')
    time.units = '{0} since {1}'.format(time_units, init_date)
    time.calendar = 'standard'
    time[:] = outtime

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
