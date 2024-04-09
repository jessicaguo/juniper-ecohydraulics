# Research compendium associated with 'Dynamic regulation of water potential in *Juniperus osteosperma* mediates ecosystem carbon fluxes'

#### Authors: Jessica S. Guo, Mallory L. Barnes, William K. Smith, William R.L. Anderegg, Steven A. Kannenberg

Manuscript accepted at New Phytologist. DOI forthcoming. 

### Project summary

We instrumented seven juniper trees with stem psychrometers at the Ameriflux site US-CdM during the growing season of 2021 to examine water potential regulation strategies and ecosystem productivity. *J. osteosperma* experiences temporally-varying hydraulic regulation, including extreme anisohydry associated with wet soils and dry air. Extreme anisohydry may serve to prolong water extraction from the soil and extend the GPP response to precipitation during the monsoon season. After accounting for the effect of light and vegetation greenness, predawn water potential corresponded to GPP better than indirect measures of water stress such as VPD and soil moisture. In drylands, dynamic and pulse-driven moisture conditions mean that ecosystem productivity must account for vegetation water stress, especially when temperature and light are not limiting. 

### Compendium contents

 - `daily_app/` contains a simple Shiny app to explore the data
 - `data_cleaned/`
    - `met_daily.Rdata` are the meteorological variables for US-CdM
    - `psy_daily_site_gapfilled.Rdata` are predawn and midday means, standard deviation, and sample size, gapfilled using a Kalman filter
    - `psy_daily_site.Rdata`are predawn and midday means, standard deviation, and sample size
    - `psy_daily.Rdata` are the predawns and middays calculated for each branch and tree
    - `psy_hourly.csv` is output from `scripts/Plot-psy-raw.Rmd` and contains the cleaned half-hourly values for each branch and tree
    - `US-Cdm_NIRv.csv` is the NIRv extracted MODIS by Mallory Barnes
 - `data_raw/` contains an assortment of .csv files used to produce the cleaned data
 - `docs/` contains the Rmarkdown and Quarto versions of the manuscript
   - `hydry-GPP-quarto/` contains the `index.qmd` of the most recent version of the manuscript
   - `hydry-GPP` contains `Hydry-GPP-ms.Rmd`, a older version of the manuscript using Rmarkdown
   
 - `renv/` records the versions of packages used in these analyses
 - `scripts/`
    - `*.R` are numbered scripts to explore, clean, gapfill, and plot the data  
      - `01-explore-psi.R` plots quick comparison of manual and psychrometer measurements
      - `02-calcluate-daily.R` converts half-hourly water potential and meterological data into daily summaries, including `data_cleaned/psy_daily.Rdata` and `data_cleaned/met_daily.Rdata`
      - `03-explore-flux.R` plots fluxes with met data
      - `04-flux-psy-daily.R` plots daily fluxes and psy, and outputs `data_cleaned/psy_daily_site.Rdata`
      - `05-gapfill-psy.R` uses Kalman imputation to gapfill missing site-level psy and outputs `data_cleaned/psy_daily_site_gapfilled.Rdata`
      - `06-map-JAS-ppt.R` creates map of % precip in JAS for the 4 corners area, from PRISM data
      - `07-plot-pulses.R` plots a stack of responses to each of the 3 monsoon pulses
      - `08-test-correl.R` compares correlation of SWP and VWC
      - `09-plot-psy.R` plots between branch correlations
    - `*.Rmd` are notebooks the detail the data diagnistic process
      - `Plot-psy-raw.Rmd` documents the data processing, including skipping maintenance days and periods of poor instrument contact. Outputs `data_cleaned/psy_hourly.csv`
      - `Plot-psy-diurnal.Rmd` provides additional interactive plots that scrutinize the 3-7 am period for signal loss, which can impact the calculation of predawn water potential
      - `Plot-psy-diagnostic.Rmd` examines the relationship between cleaned half-hourly water potentials and meterological variables
      - `Plot-pd-md.Rmd` plots the timeseries of predawn and midday water potential as calculated from `data_cleaned/psy_hourly.csv`
    - `model-*/` are individual folders for Bayesian analyses
        - `model-pd-md/` contains the JAGS model code, data, and initials for the hydraulic regulation model
      - `model-gpp-nirv/` contains the JAGS model code, data, and initials for the two-part GPP model
      - `model-flux-daily` contains the original GPP model, where GPP was regressed against soil moisture and VPD
      - `model-gapfill-pd/` contains an early attempt to use gapfill the water potential data; see `scripts/05-gapfill-psy.R` for replacement approach
      - `model-pd-daily/` contains the original attempt to predict predawns from environmental variables
  
 - `source/` contains a rounding function for aligning half-hourly time stamps
