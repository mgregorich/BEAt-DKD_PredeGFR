# Prediction model of future eGFR in people with type 2 diabetes mellitus

Code for the development, validation and web implementation of the risk prediction calcuator within the BEAt-DKD project WP1 Task 5.

This repository provides the accompanying code for model building, internal- external validation and external validation. In addition, the code for the web implementation using Shiny is available.

Gregorich, M., Heinzel, A., Kammer, M., Meiselbach, H., Böger, C., Eckhardt, K. U., Mayer, G., Heinze, G., & Oberbauer, R. (2022). Individual-specific prediction of future eGFR in people with type 2 diabetes mellitus: development and external validation. Kidney International. **(in submission)**


## Contents

- `scr` folder: contains the necessary code files for model building and validation without patient data
- `shiny` folder contains the code files for the Shiny web implementation of the prediction model

## Usage

The predictive model code is included to transparently report the model development and validation, however, it cannot be executed due to the lack of the underlying data of the PROVALID, DIACORE and the GCKD study cohorts.

The shiny app can either be started by downloading the shiny folder and executing the code or can be accessed via [here](https://beatdkd.shinyapps.io/shiny/).


## Prerequisites

The code uses the statistical software `R` (>= 4.0) 
