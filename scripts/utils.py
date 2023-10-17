"""
Utility functions to run the WOFOST phenology model.
"""

import pandas as pd
import pcse
import yaml

from pathlib import Path


# default values for the weather data (not used by the
# WOFOST phenology model but required by the PCSE framework
# as dummy values)
defaults = {"IRRAD": 20e6,
            "VAP": 5,
            "RAIN": 2,
            "E0": 0.5,
            "ES0": 0.5,
            "ET0": 0.5,
            "WIND": 10,
            "SNOWDEPTH": 0}


class WeatherDataProvider_from_WeatherStation(pcse.base.WeatherDataProvider):
    """
    Custom weather data provider (wdp) for the WOFOST phenology model.
    Adopted from Raphael Portmann (Agroscope):
    https://github.com/raphael-portmann/PhenoSwiss/blob/master/calibrate_and_run_phenology_Switzerland.ipynb
    """
    def __init__(
            self,
            weather_data: pd.DataFrame,
            elevation: float,
            lon: float,
            lat: float,
            default_values: dict = defaults):
        """
        Parameters
        ----------
        weather_data : pd.DataFrame
            Weather data.
        elevation : float
            Elevation of the weather station.
        lon : float
            Longitude of the weather station.
        lat : float
            Latitude of the weather station.
        default_values : dict, optional
            Default values for the weather data, by default defaults
        """
        super().__init__()

        self.elevation = elevation
        self.longitude = lon
        self.latitude = lat

        # loop over all days in the weather DataFrame starting from the actual
        # sowing date
        for _, row in weather_data.iterrows():
            _date = row['date'].date()
            t_min = row['T_min']
            t_max = row['T_max']

            wdc = pcse.base.WeatherDataContainer(
                DAY=_date,
                LAT=lat,
                LON=lon,
                ELEV=elevation,
                TMIN=t_min,
                TMAX=t_max,
                **default_values)
            self._store_WeatherDataContainer(wdc, _date)


def get_agromanager(
        sowing_date: pd.Timestamp,
        fpath_wheat_calender: Path
) -> dict:
    """
    Setup the agromanagement for the WOFOST phenology model.

    Parameters
    ----------
    sowing_date : pd.Timestamp
        Sowing date of the wheat.
    fpath_wheat_calender : Path
        Path to the wheat calender.
    """
    # read the wheat calender
    with open(fpath_wheat_calender, 'r') as f:
        yaml_agro = yaml.load(f, Loader=yaml.FullLoader)

    # create start dates
    standard_date = list(yaml_agro['AgroManagement'][0].keys())[0]
    campaign_date = sowing_date.date()
    crop_start_date = sowing_date.date()

    # get campaign and replace dates
    yaml_agro['AgroManagement'][0][campaign_date] = \
        yaml_agro['AgroManagement'][0].pop(standard_date)
    yaml_agro['AgroManagement'][0][campaign_date]['CropCalendar'][
        'crop_start_date'] = crop_start_date
    return yaml_agro
