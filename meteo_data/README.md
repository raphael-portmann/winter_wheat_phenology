# Meteorological data for WOFOST simulations

This repository contains the meteorological data per grid cell required to run WOFOST and determine the sowing date.
The data is provided for each grid cell of the simulation on a daily basis, extracted from the [MeteoSwiss](https://www.meteoschweiz.admin.ch/klima/klima-der-schweiz/raeumliche-klimaanalysen.html) gridded data
product (1 by 1 km resolution, same grid cell layout as used for the simulations).

Since a revision of the of the Ordinance on Meteorology and Climatology (MetV; SR 429.11) in March 2023, the data should be made available free-of-charge as open-governmental data. The publication of the MeteoSwiss products as open-governmental data is currently [work in progress](https://github.com/MeteoSwiss/publication-opendata). From 2025 onwards, users should be able to pull the data from MeteoSwiss directly.

To allow users to reproduce our results, the extracted data is made available as CSV files using Git's large file storage (LFS). Many thanks to [MeteoSwiss](https://www.meteoswiss.admin.ch/#tab=forecast-map) for the permission to do so 🙏.

| File Name | Description | Unit | Original Dataset Description
| --------- | ----------- | ---- | ----------------------------- |
| RhiresD_1971-01-01-2020-12-31.csv | Daily precipitation sum (midnight to midnight) | mm | [RhiresD](https://www.meteoschweiz.admin.ch/dam/jcr:4f51f0f1-0fe3-48b5-9de0-15666327e63c/ProdDoc_RhiresD.pdf) |
| TmaxD_1971-01-01-2020-12-31.csv | Daily maximum air temperature 2m above ground | deg C | [TmaxD](https://www.meteoschweiz.admin.ch/dam/jcr:818a4d17-cb0c-4e8b-92c6-1a1bdf5348b7/ProdDoc_TabsD.pdf) |
| TminD_1971-01-01-2020-12-31.csv | Daily minimum air temperature 2m above ground | deg C | [TminD](https://www.meteoschweiz.admin.ch/dam/jcr:818a4d17-cb0c-4e8b-92c6-1a1bdf5348b7/ProdDoc_TabsD.pdf) |
