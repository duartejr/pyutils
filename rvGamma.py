#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gamma distribution rainfall correction for ensemble forecasts.

The `comp_rv_gama` function applies gamma-distribution quantile mapping
to correct a single forecast ensemble. The multiprocessing wrapper and
the original basin-looping entry point depended on a private 'utils'
package that is not part of this repository; they have been removed.
"""

import os
import multiprocessing as mp
from glob import glob
from datetime import datetime

import numpy as np
from dateutil.relativedelta import relativedelta

from rv.RvGama import rv_gama


trimestres = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
}


def comp_rv_gama(fcst_files, obs_files, subbasin, years, out_dir):
    """
    Apply gamma-distribution quantile correction to ensemble forecast files.

    Parameters
    ----------
    fcst_files : list of str
        Sorted list of forecast .dat file paths.
    obs_files : list of str
        Sorted list of observed .csv file paths.
    subbasin : str
        Sub-basin identifier used in file names.
    years : array-like of int
        Years to process (leave-one-out cross-validation).
    out_dir : str
        Directory where corrected output files will be written.

    Returns
    -------
    np.ndarray
        Array of shape (n_cases, 3 + n_lead_days) with columns
        [year, month, day, corrected_values...].
    """
    pr_cor = []

    for year in years:
        y = years[np.where(years != year)[0]]

        for month in range(1, 13):
            base = "{0}_{1}{2:02d}".format(subbasin, year, month)
            days = [
                int(fn.split("_")[-3][-2:]) for fn in fcst_files if 1 + fn.find(base)
            ]

            for tr, months in trimestres.items():
                if month in months:
                    tri = tr

            files_tri_param = [
                fn
                for fn in fcst_files
                if any(
                    fn.find("{0}_{1}{2:02d}".format(subbasin, y_i, m)) != -1
                    for y_i in y
                    for m in trimestres[tri]
                )
            ]

            pr_tri_param = []
            dates = []

            for file_tri in files_tri_param:
                pr = np.loadtxt(file_tri, delimiter=",")
                try:
                    pr = np.mean(pr[:, 2:], 0)
                except IndexError:
                    pr = pr[2:]
                pr_tri_param.append(pr)
                dates.append(file_tri.split("_")[-2])

            pr_obs_param = []

            for date in dates:
                dt = datetime.strptime(date, "%Y%m%d")
                dts = [dt + relativedelta(days=i) for i in range(10)]
                bases = [
                    "{0}_{1}".format(subbasin, x.strftime("%Y%m%d")) for x in dts
                ]
                files = [
                    fn
                    for fn in obs_files
                    if any(fn.find(b) != -1 for b in bases)
                ]
                pr_week = []
                for file in files:
                    try:
                        pr_obs = np.loadtxt(file, delimiter=",")
                    except Exception:
                        pr_obs = np.array([np.nan, np.nan, np.nan])
                    try:
                        pr_week.append(np.nanmean(pr_obs[:, 2:]))
                    except IndexError:
                        pr_week.append(np.nanmean(pr_obs[2:]))

                pr_week.extend([np.nan] * (10 - len(pr_week)))
                pr_obs_param.append(np.array(pr_week))

            pr_fcst_param = np.array(pr_tri_param)
            pr_obs_param = np.array(pr_obs_param)

            for day in days:
                base = "{0}_{1}{2:02d}{3:02d}".format(subbasin, year, month, day)
                file_cor = [fn for fn in fcst_files if 1 + fn.find(base)][0]
                pr = np.loadtxt(file_cor, delimiter=",")
                try:
                    pr = np.mean(pr[:, 2:], 0)
                except IndexError:
                    pr = pr[2:]

                corrected = [
                    rv_gama(pr_fcst_param[:, i], pr_obs_param[:, i], pr[i])
                    for i in range(len(pr))
                ]
                pr_cor.append([year, month, day] + corrected)

    os.makedirs(out_dir, exist_ok=True)
    return np.array(pr_cor)
