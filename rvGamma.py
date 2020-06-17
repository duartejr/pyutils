#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 18:36:27 2019

@author: duarte
"""

from rv.RvGama import rv_gama
from utils.get_tri_files import get_tri_files
import os
import multiprocessing as mp
import utils.basinsf as basinsf
import numpy as np
from glob import glob
from datetime import datetime
from dateutil.relativedelta import relativedelta
from utils.save_csv import save_csv

##############################################################################
home = os.path.expanduser('~')+'/MEGA/FCPC'
trimestres = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8],
              'SON': [9, 10, 11]}
years = np.arange(2008, 2018)
##############################################################################


def comp_rv_gama(args):
    basin, subbasin, model = args
    fcst_dir = '{0}/io/pr/fcst/{1}/diario/{2}/{3}'.format(home, model,
                                                          basin, subbasin)
    obs_dir = '{0}/io/pr/obs/hidroweb/diario/{1}/{2}'.format(home, basin,
                                                             subbasin)

    base_name = [subbasin+'_'+str(ano) for ano in years]
    fcst_files = [fn for fn in sorted(glob(fcst_dir+'/*.dat')) if
                  any(1+fn.find(str(base)) for base in base_name)]
    obs_files = [fn for fn in sorted(glob(obs_dir+'/*.csv')) if
                 any(1+fn.find(str(base)) for base in base_name)]
    pr_cor = []
    

    for year in years:

        y = years[np.where(years != year)[0]]

        for month in range(1, 13):
            print(year, month)
            base = '{0}_{1}{2:02d}'.format(subbasin, year, month)
            days = [int(fn.split('_')[-3][-2:]) for fn in fcst_files
                    if 1+fn.find(base)]

            for tr, months in trimestres.items():
                if month in months:
                    tri = tr

            files_tri_param = get_tri_files(fcst_files, subbasin,
                                            trimestres[tri], y)

            pr_tri_param = []
            dates = []

            for file_tri in files_tri_param:
                pr = np.loadtxt(file_tri, delimiter=',')

                try:
                    pr = np.mean(pr[:, 2:], 0)
                except:
                    pr = pr[2:]

                pr_tri_param.append(pr)
                dates.append(file_tri.split('_')[-2])

            pr_obs_param = []

            for date in dates:
                date = datetime.strptime(date, '%Y%m%d')
                pr_week = []
                dts = [date + relativedelta(days=i) for i in range(10)]
                bases = ['{0}_{1}'.format(subbasin, datetime.strftime(x,
                         '%Y%m%d')) for x in dts]
                files = [fn for fn in obs_files if any(1+fn.find(base) for
                                                       base in bases)]
                pr_week = []

                for file in files:
                    try:
                        pr_obs = np.loadtxt(file, delimiter=',')
                    except:
                        pr_obs = np.array([np.nan, np.nan, np.nan])

                    try:
                        pr = np.nanmean(pr_obs[:, 2:])
                    except:
                        pr = np.nanmean(pr_obs[2:])

                    pr_week.append(pr)

                if len(pr_week) < 10:
                    for ns in range(10 - len(pr_week)):
                        pr_week.append(np.nan)

                pr_week = np.array(pr_week)
                pr_obs_param.append(pr_week)

            pr_fcst_param = np.array(pr_tri_param)
            pr_obs_param = np.array(pr_obs_param)

            for day in days:

                base = '{0}_{1}{2:02d}{3:02d}'.format(subbasin, year,
                                                      month, day)
                file_cor = [fn for fn in fcst_files if 1+fn.find(base)][0]
                print(file_cor)
                pr = np.loadtxt(file_cor, delimiter=',')

                try:
                    pr = np.mean(pr[:, 2:], 0)
                except:
                    pr = pr[2:]

                aux = []

                for i in range(len(pr)):
                    aux.append(rv_gama(pr_fcst_param[:, i],
                                                     pr_obs_param[:, i],
                                                     pr[i]))

                pr_cor.append([year, month, day] + aux)

    out_dir = '{0}/io/pr/fcst/{1}/cor/gam/{2}/{3}'.format(home, model, basin,
                                                          subbasin)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    start_date = '{0}01'.format(years[0])
    end_date = '{0}12'.format(years[-1])

    if subbasin != basin:
        out_file = '{0}/pr_fcst_cor_gam_{1}_{2}_{3}_{4}_{5}.csv'.format(out_dir,
                                                                        model,
                                                                        basin,
                                                                        subbasin,
                                                                        start_date,
                                                                        end_date)
    else:
        out_file = '{0}/pr_fcst_cor_gam_{1}_{2}_{3}_{4}.csv'.format(home, model,
                                                                    basin,
                                                                    start_date,
                                                                    end_date)

    pr_cor = np.array(pr_cor)
    save_csv(out_file, pr_cor)


def multiprocess(args):
    p = mp.Pool(mp.cpu_count())
    r = p.map(comp_rv_gama, args)
    return 1


if __name__ == "__main__":
    args = []

    for model in basinsf.models()[0:]:
        print('Modelo -----> {}'.format(model))

        for basin in basinsf.basins():
            print('\tBacia -----> {}'.format(basin))

            for subbasin in basinsf.subbasins(basin):
                print('\t\tSub-bacia -----> {}'.format(subbasin))
                args.append([basin, subbasin, model])
                a = args[0]
                comp_rv_gama(a)
                break
            break
        break

                

    multiprocess(args)
