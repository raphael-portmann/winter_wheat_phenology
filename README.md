# Simulation of winter wheat heading dates using WOFOST-PHENOLOGY for two Swiss varieties

This repository contains code and data to simulate heading dates of Swiss varieties for nearly five decades (1971 to 2020) in the main wheat production areas of Switzerland. We use the phenology model of the [WOFOST crop growth model](https://www.wur.nl/en/show/A-gentle-introduction-to-WOFOST.htm) implemented in the [Python Crop Simulation Environment](https://github.com/ajwdewit/pcse) made available under [EUPL license, Version 1.1](http://ec.europa.eu/idabc/eupl). The implementation of the phenology model is based on tutorials by [Allard de Wit, Wageningen University and Research](https://github.com/ajwdewit) available [here](https://github.com/ajwdewit/pcse_notebooks) and work by [Raphael Portmann, Agroscope Reckenholz](https://github.com/raphael-portmann), who performed the [model calibration](https://github.com/raphael-portmann/PhenoSwiss/blob/master/calibrate_and_run_phenology_Switzerland.ipynb).

## Setup the environment

Clone the repository from GitHub. Since the [meteorological data](/meteo_data/) is rather huge (nearly 5 GB) we have to use [git large file storage](https://git-lfs.com/), which should be installed prior to cloning to also clone these large files. Please note, that cloning might take a while depending on your internet connection.

```bash
git lfs install
git clone git@github.com:EOA-team/winter_wheat_phenology.git
cd wofost_phenology
```

Create a Python environment and install the dependencies (shown for Linux systems).

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run the simulation

To run the simulations, execute the script [run_wofost_phenology.py](/scripts/run_wofost_phenology.py).
Due to the long run time of the simulations, pre-compiled results can be downloaded [here](http://hdl.handle.net/20.500.11850/637092). If you would like to reproduce the results from scratch, continue as described below.

```
cd scripts
python run_wofost_phenology.py
```

It will automatically generate compute the phenology of two Swiss winter wheat varieties (Arina and CH-Claro) for the time period between 1971 and 2020. It uses [gridded meteorological](/meteo_data/) in a 1 by 1 km spatial resolution available on a daily basis from the [Swiss Federal Office of Meteorology and Climatology MeteoSwiss](https://www.meteoswiss.admin.ch/#tab=forecast-map). The simulation run on all 1 by 1 km grid cells that have a share of 20% crop land area according to the [Swiss national crop and grassland layer produced by Pazur et al. (2022)](https://www.dora.lib4ri.ch/wsl/islandora/object/wsl:29612) produced by the [Swiss Federal Institute for Forest, Snow and Landscape Research WSL](https://www.wsl.ch/en).

The [calibrated model parameter TSUM1](/scripts/genotypes_Tsum1_opt/) per variety were provided by [Raphael Portmann, Agroscope Reckenholz](https://github.com/raphael-portmann) using [this Jupyter notebook](https://github.com/raphael-portmann/PhenoSwiss/blob/master/calibrate_and_run_phenology_Switzerland.ipynb) with data made available by Dario Fossati, Agroscope Changings.

---
**_IMPORTANT NOTE ON COMPUTATION TIME AND REQUIREMENTS_** 

The script will run for 5230 grid cells (5230 km2) for the harvest years 1972 to 2020 for two genotypes and the earliest and latest posssible sowing date resulting in a total of 1 025 080 WOFOST runs. To speed up the computation, the code is parallelized by years using `number-of-cpus minus 1` in its default configuration. **Still, the computation takes some time!** (up to 300 minutes on a Fedora workstation with AMD Ryzen Threadripper PRO 3955WX 16-Cores, 126 GB RAM). If your computational capacities are limited we recommend to download the pre-computed results from [here](http://hdl.handle.net/20.500.11850/637092).

---

The results are written to [results](/results/) and stored as GeoPackage file per year. Alternatively, download the [pre-computed results](http://hdl.handle.net/20.500.11850/637092), unzip them and place them in the [results](/results/) directory.

## Results

The results with the simulated heading dates are stored per year (see [here](/sample_data/results_ww_gs_1971-1972.gpkg) for a sample). On the example of the harvest year 1972, the file contents are explained below:

```python
import geopandas as gpd

# read as GeoDataFrame
gdf = gpd.read_file('sample_data/results_ww_gs_1971-1972.gpkg')
gdf.info()
```
outputs
```python
<class 'geopandas.geodataframe.GeoDataFrame'>
RangeIndex: 20964 entries, 0 to 20963
Data columns (total 10 columns):
 #   Column               Non-Null Count  Dtype   
---  ------               --------------  -----   
 0   id                   20964 non-null  int64   
 1   harvest_year         20964 non-null  int64   
 2   sowing_date          20964 non-null  object  
 3   sowing_date_type     20964 non-null  object  
 4   genotype             20964 non-null  object  
 5   emergence_date       20964 non-null  object  
 6   anthesis_date        20964 non-null  object  
 7   crop_land_area_km2   20964 non-null  float64 
 8   crop_land_area_perc  20964 non-null  float64 
 9   geometry             20964 non-null  geometry
dtypes: float64(2), geometry(1), int64(2), object(5)
memory usage: 1.6+ MB
```
The data is georeferenced in Swiss LV95 coordinates (EPSG:2056):

```bash
>>> gdf.crs 
<Projected CRS: EPSG:2056>
Name: CH1903+ / LV95
Axis Info [cartesian]:
- E[east]: Easting (metre)
- N[north]: Northing (metre)
Area of Use:
- name: Liechtenstein; Switzerland.
- bounds: (5.96, 45.82, 10.49, 47.81)
Coordinate Operation:
- name: Swiss Oblique Mercator 1995
- method: Hotine Oblique Mercator (variant B)
Datum: CH1903+
- Ellipsoid: Bessel 1841
- Prime Meridian: Greenwich

```

In this example, `harvest_year` is always 1972 as the simulation produces a single result per year:

```bash
>>> gdf.harvest_year.unique()
array([1972])
```

`sowing_date` is the sowing date of winter wheat from the sowing date algorithm proposed by [Holzkaemper et al. (2015)](https://doi.org/10.1007/s10113-014-0627-7). We use the `earliest` and `latest` possible sowing date (according to the algorithm) as the exact sowing date(s) within a grid cell are unknown:

```bash
>>> gdf.sowing_date.unique()
array(['1971-10-27', '1971-11-06', '1971-10-30', '1971-10-31',
       '1971-11-01', '1971-10-10', '1971-10-20', '1971-11-02',
       '1971-10-07', '1971-10-25', '1971-10-09', '1971-10-21'],
      dtype=object)
```

Which sowing date type was used can be identified using the `sowing_date_type` attribute:

```bash
>>> gdf.sowing_date_type.unique()
array(['earliest', 'latest'], dtype=object)
```

The results were run per genotype (i.e., variety):

```bash
>>> gdf.genotype.unique()
array(['Arina', 'CH_Claro'], dtype=object)
```

The actual variable of interest is the `anthesis_date` which is -- in this case -- the heading date (in the original WOFOST setup, it is the anthesis, i.e., flowering date. Through the optimization using heading dates, we forced the model to predict heading instead of flowering).

```bash
>>> gdf.anthesis_date.unique()
array(['1972-06-05', '1972-06-06', '1972-06-03', '1972-06-04',
       '1972-06-02', '1972-05-31', '1972-06-01', '1972-06-07',
       '1972-06-09', '1972-06-11', '1972-06-08', '1972-06-10',
       '1972-06-12', '1972-06-14', '1972-06-13', '1972-06-15',
       '1972-06-16', '1972-06-18', '1972-06-20', '1972-06-22',
       '1972-06-17', '1972-06-19', '1972-06-23', '1972-06-21',
       '1972-06-24', '1972-05-30', '1972-05-29', '1972-05-28',
       '1972-06-25', '1972-06-26', '1972-06-27', '1972-06-29',
       '1972-06-28', '1972-07-01', '1972-06-30', '1972-07-02',
       '1972-07-05', '1972-07-03', '1972-07-06', '1972-07-04'],
      dtype=object)
```

The heading date can be easily visualized as a geo-referenced map in Swiss LV95 coordinates (requires `matplotlib` to be installed -- `pip install matplotlib`):

```bash
>>> gdf.plot(column='anthesis_date', cmap='viridis')
<Axes: >
```

![Heading Dates 1972](/sample_data/heading_dates_1972.png)