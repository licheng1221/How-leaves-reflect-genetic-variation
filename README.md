# Evaluating Potential of Leaf Reflectance Spectra to Monitor Plant Genetic Variation

This repository contains the code and data for the publication "Evaluating Potential of Leaf Reflectance Spectra to Monitor Plant Genetic Variation". 

## Data

The spectral measurement data underlying the results presented in this paper are available as a published dataset at SPECCHIO http://sc22.geo.uzh.ch:8080/SPECCHIO_Web_Interface/search, with Keyword: UZH_SG_Nicotiana_attenuata_ASD

The processed data can be found under the folder "Data_rds", which are the input files for the downstream analyses.

## Code

There are five R scripts in this repository:

- `Function.R`
- `Plots_Models.R`
- `Plots_Models.R`
- `Plots.PCA.R`
- `RawDataProcess.R`

We used `RawDataProcess.R` to process the raw data, which can be downloaded from the public dataset mentioned above. Then, we used `Plots_PCA.R` and `PLots_Models.R` to perform PCA, run all the models, and generate the figures in the manuscript.

## Dependencies

The results and figures were generated with R (version 4.3.0 (2023-04-21)) and the following R packages:

- spectrolab (version 0.0.10)
- DescTools(version 0.99.49)
- lme4 (version 1.1-31)
- nlme (version 3.1-160)
- MASS (version 7.3-58.1)
- emmeans (version 1.8.2)
- tidyverse (version 1.3.2)

We used the R package report (version 0.5.5) to generate formatted reports of models.
