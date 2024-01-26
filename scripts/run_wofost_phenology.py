"""
Run WOFOST phenology for selected locations, years and genotypes to
predict the date of anthesis. The sowing date is estimated using the
approach proposed by Holzkaemper et al., 2015 using the earliest and
latest sowing dates for each location and year.

To speed up the simulations, the phenology is run in parallel for each
year. The number of parallel processes is defined by the number of
available CPUs minus 1. Thus, if you have 4 CPUs, 3 processes (i.e., 3 years)
will be run in parallel.

MIT License

Copyright (c) 2023 Lukas Valentin Graf, EOA-team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import geopandas as gpd
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
import pcse

from datetime import date, datetime
from pathlib import Path
from pcse.models import Wofost72_Phenology
from pcse.base import ParameterProvider
from pcse.exceptions import WeatherDataProviderError

from utils import get_agromanager, WeatherDataProvider_from_WeatherStation

# parameterization of the algorithm for sowing date estimation based
# on Holzkaemper et al., 2015
TMAX_SOWING = 12                        # mean air temperature in deg C
PRECIP_THRESHOLDS = [20, 16, 12, 8, 4]  # daily rainfall in mm

# setup logging
logger = logging.getLogger(__name__)


def estimate_sowing_date(
        tmean: pd.Series,
        daily_precip: pd.Series,
        choice: int
) -> date:
    """
    Winter Wheat sowing date estimation based on Holzkaemper et al., 2015
    (https://doi.org/10.1007/s10113-014-0627-7) for a single point location.

    :param tmean:
        mean air temperature in deg Celsius
    :param daily_precip:
        daily precipitation in mm
    :param choice:
        in case there are multiple time windows fulfilling the criteria specify
        which window to select as sowing date.
    :returns:
        estimated winter wheat sowing date index for the location
    """
    df = pd.DataFrame({'tmean': tmean, 'daily_precip': daily_precip})
    # restrict time window between October 1 (earliest sowing date October 7th)
    # and November 7th
    df = df.iloc[:37].copy()

    # find period where tmean < 12 deg C for at least 6 consecutive days
    temp_reduced = np.zeros(df.shape[0], dtype=bool)
    for idx in range(6, df.shape[0]):
        if df.tmean.iloc[idx - 6:idx].le(TMAX_SOWING).all():
            temp_reduced[idx] = True
    # and the amount of daily rainfall was <20, <16, <12, <8, <4
    # for 5 consecutive days
    precip_reduced = np.zeros(df.shape[0], dtype=bool)
    for idx in range(5, df.shape[0]):
        if (
            df.daily_precip.iloc[idx - 5] < PRECIP_THRESHOLDS[0] and
            df.daily_precip.iloc[idx - 4] < PRECIP_THRESHOLDS[1] and
            df.daily_precip.iloc[idx - 3] < PRECIP_THRESHOLDS[2] and
            df.daily_precip.iloc[idx - 2] < PRECIP_THRESHOLDS[3] and
            df.daily_precip.iloc[idx - 1] < PRECIP_THRESHOLDS[4]
        ):
            precip_reduced[idx] = True

    matching_idx = np.where(temp_reduced & precip_reduced)[0]
    # take the first matching index and use that date as the sowing date
    dates = df.index

    if len(matching_idx) == 0:
        raise ValueError('No sowing date found')
    if choice > 0:
        # decrease positive choice numbers to be conform with Python
        choice -= 1
    if choice >= len(matching_idx):
        # if there is no choice, set the sowing date index
        # to -1 to indicate that there is no choice left
        return -1
    sowing_date_idx = dates[matching_idx[choice]]
    return sowing_date_idx


def process_year(
        year: int,
        meteo: dict,
        units: gpd.GeoDataFrame,
        output_dir: Path,
        genotypes: "list[str]",
        fpath_tsum1_opt_dir: Path,
        fpath_wheat_calendar: Path,
        sowing_date: "datetime or None"
) -> None:
    """
    For each year and spatial unit, estimate the sowing date and run the
    WOFOST phenology model using optimized TSUM1 parameters for the Swiss
    winter wheat varieties Arina and CH_Claro.

    :param year: year to process
    :param meteo: dictionary with weather data for the year
    :param units: GeoDataFrame with spatial units
    :param output_dir: directory to store the results
    :param genotypes: list of genotypes to run
    :param fpath_tsum1_opt_dir: directory with optimized TSUM1 parameters
    :param fpath_wheat_calendar: path to the wheat calendar file
    :param sowing_date: fixed sowing date with dummy year
    """

    logger.info(f'Working on growing season {year} - {year+1}')
    fpath_out = output_dir.joinpath(f'results_ww_gs_{year}-{year+1}.gpkg')
    if fpath_out.exists():
        logger.info('Results exist already -> skipping')
        return

    date_start = date(year, 10, 1)
    date_end = date(year+1, 9, 30)
    dates = pd.date_range(date_start, date_end)

    # prepare weather data of the year
    meteo_year = {}
    for variable in meteo.keys():
        meteo_year[variable] = meteo[variable][
            (meteo[variable].date >= date_start) &
            (meteo[variable].date <= date_end)
        ].copy()

    # estimate the sowing date per spatial unit
    yearly_results_list = []
    for _, unit in units.iterrows():
        meteo_year_unit = {}
        try:
            for variable in meteo_year.keys():
                meteo_year_unit[variable] = \
                    meteo_year[variable][str(unit.ID)].copy()
        except Exception as e:
            logger.error(f'{unit.NAME}: {e}')
            continue
        # Estimate sowing date if no fixed sowing date is passed. To
        # To estimate the sowing date, the daily mean air temperature is required
        if sowing_date is None:
            tmean = (
                meteo_year_unit['TmaxD'].values +
                meteo_year_unit['TminD'].values) * 0.5
            tmean = pd.Series(tmean)
            daily_precip = meteo_year_unit['RhiresD'].values
            daily_precip = pd.Series(daily_precip)

            # first (earliest) possible sowing date
            try:
                sowing_date_idx_early = estimate_sowing_date(
                    tmean=tmean,
                    daily_precip=daily_precip,
                    choice=1
                )
                sowing_date_early = dates[sowing_date_idx_early]
            except ValueError:
                logger.info(
                    f'{unit.ID} {year}: No sowing date found. ' +
                    'Using Default.')
                sowing_date_early = datetime(year, 10, 7)

            # last (latest) possible sowing date
            try:
                sowing_date_idx_late = estimate_sowing_date(
                    tmean=tmean,
                    daily_precip=daily_precip,
                    choice=-1
                )
                sowing_date_late = dates[sowing_date_idx_late]
            except ValueError:
                logger.info(
                    f'{unit.ID} {year}: No sowing date found. ' +
                    'Using Default.')
                # use default late sowing date (7th ofNovember)
                sowing_date_late = datetime(year, 11, 7)

        # prepare weather data for WOFOST
        weather = pd.DataFrame({
            'date': list(dates),
            'T_min': meteo_year_unit['TminD'].values,
            'T_max': meteo_year_unit['TmaxD'].values
        })

        # get longitude and latitude
        unit_gdf = gpd.GeoDataFrame(
            [unit], crs=units.crs)
        unit_gdf.geometry = unit_gdf.geometry.centroid
        unit_gdf.to_crs(epsg=4326, inplace=True)
        longitude = unit_gdf.geometry.x.values[0]
        latitude = unit_gdf.geometry.y.values[0]
        # use the true elevation if available
        if 'Elevation' in unit_gdf.columns:
            elevation = unit_gdf.Elevation.values[0]
        # else set elevation to dummy value
        else:
            elevation = 450

        # setup agrar manager and weather data provider for WOFOST
        if sowing_date is None:
            agromanager_early = get_agromanager(
                sowing_date=sowing_date_early,
                fpath_wheat_calender=fpath_wheat_calendar
            )['AgroManagement']
            agromanager_late = get_agromanager(
                sowing_date=sowing_date_late,
                fpath_wheat_calender=fpath_wheat_calendar
            )['AgroManagement']
        else:
            sowing_date=datetime(year,sowing_date.month,sowing_date.day)
            agromanager = get_agromanager(
                sowing_date=sowing_date,
                fpath_wheat_calender=fpath_wheat_calendar
            )['AgroManagement']

        wdp = WeatherDataProvider_from_WeatherStation(
            weather_data=weather,
            elevation=elevation,
            lon=longitude,
            lat=latitude
        )

        # loop over genotypes and calculate the date of anthesis
        for genotype in genotypes:
            with open(
                fpath_tsum1_opt_dir.joinpath(f'{genotype.lower()}.txt')
            ) as src:
                tsum1_opt = float(src.read().split()[0])
                tsum2_opt = float(src.read().split()[1])

            cropd = pcse.fileinput.YAMLCropDataProvider()
            cropd.set_active_crop('wheat', 'Winter_wheat_105')
            cropdata = cropd.copy()
            cropdata['TSUM1'] = tsum1_opt
            cropdata['TSUM2'] = tsum2_opt

            # set the parameters
            parameters = ParameterProvider(
                cropdata=cropdata
            )

            if sowing_date is None:
                    # WOFOST72 Phenology model early sowing date
                    wofsim_early = Wofost72_Phenology(
                        parameters, wdp, agromanager_early)
                    # run till terminate
                    try:
                        wofsim_early.run_till_terminate()
                    except WeatherDataProviderError as e:
                        logger.error(f'{unit.NAME}: {e}')
                        continue
                    summary_early = wofsim_early.get_summary_output()

                    # late sowing date
                    wofsim_late = Wofost72_Phenology(
                        parameters, wdp, agromanager_late)
                    # run till terminate
                    try:
                        wofsim_late.run_till_terminate()
                    except WeatherDataProviderError as e:
                        logger.error(f'{unit.NAME}: {e}')
                        continue
                    summary_late = wofsim_late.get_summary_output()
            else:
                # WOFOST72 Phenology model early sowing date
                wofsim = Wofost72_Phenology(
                        parameters, wdp, agromanager)
                # run till terminate
                try:
                    wofsim.run_till_terminate()
                except WeatherDataProviderError as e:
                    logger.error(f'{unit.NAME}: {e}')
                    continue
                summary = wofsim.get_summary_output()

            # save results in compact format
            if sowing_date is None:
                res_early_sowing = {
                    'id': unit.ID,
                    'harvest_year': year + 1,
                    'sowing_date': sowing_date_early.date().strftime(
                        '%Y-%m-%d'),
                    'sowing_date_type': 'earliest',
                    'genotype': genotype,
                    'emergence_date': summary_early[0]['DOE'].strftime(
                        '%Y-%m-%d'),
                    'anthesis_date': summary_early[0]['DOA'].strftime(
                        '%Y-%m-%d'),
                    'maturity_date': summary_early[0]['DOM'].strftime(
                        '%Y-%m-%d'),
                    'crop_land_area_km2': unit.crop_land_km2,
                    'crop_land_area_perc': unit.crop_land_perc,
                    'geometry': unit.geometry
                }
                res_late_sowing = {
                    'id': unit.ID,
                    'harvest_year': year + 1,
                    'sowing_date': sowing_date_late.date().strftime(
                        '%Y-%m-%d'),
                    'sowing_date_type': 'latest',
                    'genotype': genotype,
                    'emergence_date': summary_late[0]['DOE'].strftime(
                        '%Y-%m-%d'),
                    'anthesis_date': summary_late[0]['DOA'].strftime(
                        '%Y-%m-%d'),
                    'maturity_date': summary_late[0]['DOM'].strftime(
                        '%Y-%m-%d'),
                    'crop_land_area_km2': unit.crop_land_km2,
                    'crop_land_area_perc': unit.crop_land_perc,
                    'geometry': unit.geometry
                }
                yearly_results_list.append(res_early_sowing)
                yearly_results_list.append(res_late_sowing)
            else:
                res = {
                    'id': unit.ID,
                    'harvest_year': year + 1,
                    'sowing_date': sowing_date.date().strftime(
                        '%Y-%m-%d'),
                    'sowing_date_type': 'latest',
                    'genotype': genotype,
                    'emergence_date': summary[0]['DOE'].strftime(
                        '%Y-%m-%d'),
                    'anthesis_date': summary[0]['DOA'].strftime(
                        '%Y-%m-%d'),
                    'maturity_date': summary[0]['DOM'].strftime(
                        '%Y-%m-%d'),
                    'crop_land_area_km2': unit.crop_land_km2,
                    'crop_land_area_perc': unit.crop_land_perc,
                    'geometry': unit.geometry
                }
                yearly_results_list.append(res)

    # concatenate yearly results
    yearly_results_df = pd.DataFrame(yearly_results_list)
    yearly_results_gdf = gpd.GeoDataFrame(
        yearly_results_df,
        geometry=yearly_results_df.geometry,
        crs=units.crs
    )
    yearly_results_gdf.to_file(fpath_out, driver='GPKG')

    logger.info(f'Finished growing season {year} - {year+1}')


def run(
        spatial_units: "Path | gpd.GeoDataFrame",
        input_data_dir: Path,
        output_dir: Path,
        years: "list[int]",
        genotypes: "list[str]",
        fpath_wheat_calendar: Path,
        fpath_tsum1_opt_dir: Path
) -> None:
    """
    Run the WOFOST phenology model for the selected locations and
    years.

    :param spatial_units:
        GeoDataFrame with the selected spatial units (locations)
    :param input_data_dir:
        directory with the input meteorological data
    :param output_dir:
        directory where the results should be stored
    :param years:
        list of years to run the model for.
    :param genotypes:
        list of genotypes to run the model for.
    :param fpath_wheat_calendar:
        path to the wheat calendar file
    :param fpath_tsum1_opt_dir:
        path to the directory with the optimized TSUM1 parameters
        per genotype.
    """
    # read unitality data
    if isinstance(spatial_units, Path):
        units = gpd.read_file(spatial_units)
    elif isinstance(spatial_units, gpd.GeoDataFrame):
        units = spatial_units.copy()
    else:
        raise TypeError(
            'spatial units must be either a Path or a GeoDataFrame')
    
    #add latitude longitude to units
    units_new = gpd.GeoDataFrame.copy(units)
    units_new['geometry']=units_new.geometry.centroid
    units_new=units_new.to_crs(epsg=4326)
    units['latitude']=units_new.geometry.y.values
    units['longitude']=units_new.geometry.x.values

    # read the weather data
    meteo = {}
    for fpath_meteo in input_data_dir.glob('*.csv'):
        variable = fpath_meteo.name.split('_')[0]
        meteo_df = pd.read_csv(fpath_meteo)
        if 'date' in meteo_df.keys():
            meteo_df['date'] = pd.to_datetime(meteo_df.date).dt.date
        elif 'time' in meteo_df.keys():
            meteo_df['date'] = pd.to_datetime(meteo_df.time).dt.date
        else:
            raise ValueError('Neither "time" nor "date" found in meteo data')
        meteo[variable] = meteo_df

    #set fixed sowing date (with dummy year that is set automatically during the processing of the year)
    #use sowing_date=None if sowing date is to be estimated with weather data (requires precipitation data)
    sowing_date=datetime(1999,10,15)

    # loop through the years using a multiprocessing pool
    with mp.Pool(processes=mp.cpu_count()-1) as pool:
        pool.starmap(
            process_year,
            zip(years,
                [meteo] * len(years),
                [units] * len(years),
                [output_dir] * len(years),
                [genotypes] * len(years),
                [fpath_tsum1_opt_dir] * len(years),
                [fpath_wheat_calendar] * len(years),
                [sowing_date] * len(years)))


if __name__ == '__main__':

    # fix working directory for all IDEs
    cwd = Path(__file__).absolute().parent

    # set paths
    # default winter wheat calendar required by WOFOST
    fpath_wheat_calendar = cwd.joinpath('wheat_calendar.txt')
    # directory with the optimized TSUM1 parameters
    fpath_tsum1_opt_dir = cwd.joinpath('genotypes_Tsum1_opt')
    # input data directory with meteorological data.
    # You will need git's large file storage (lfs) to download the data
    # from the repository (on git clone)
    input_data_dir = cwd.parent.joinpath('meteo_data/MeteoSwiss')
    # grid cells to run the model for (GeoPackage file)
    spatial_units = cwd.joinpath('spatial_units.gpkg')

    # set years to run the model for
    years = list(range(1971, 2020))
    # set genotypes
    genotypes = ['Arina'] #, 'CH_Claro']

    # set output directory. The results will be stored in a
    # GeoPackage file per year.
    output_dir = cwd.parent.joinpath('results')

    # run the model
    run(
        spatial_units=spatial_units,
        input_data_dir=input_data_dir,
        output_dir=output_dir,
        years=years,
        genotypes=genotypes,
        fpath_wheat_calendar=fpath_wheat_calendar,
        fpath_tsum1_opt_dir=fpath_tsum1_opt_dir
    )
