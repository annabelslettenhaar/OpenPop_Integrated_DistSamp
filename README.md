---
editor_options: 
  markdown: 
    wrap: 72
---

[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

# An Integrated Predator-Prey model on gyrfalcons, ptarmigan, rodents and climate interactions

## What is in this repository?

This repository contains code for a workflow analysing data from the TOV-I 
monitoring project. We use data from gyrfalcon nest monitoring, ptarmigan 
line transect surveys, and rodent trapping surveys from 1991 to 2020. The 
model utilizes the age-structured survey data from the ptarmigan surveys 
and auxiliary data from the gyrfalcons and rodents to jointly estimate 
changes in ptarmigan demographic rates (recruitment rate and survival 
probability), gyrfalcon vital rates (productivity and brood initiation rates), 
and how these dynamics are related to rodent densities and temperature 
changes. It is a multi-area model, meaning it simultaneously models processes 
across three study areas.

The original IDSM workflow that estimates ptarmigan recruitment and survival 
was specifically written for data collected through the [Norwegian monitoring 
program for tetraonid birds](https://honsefugl.nina.no/Innsyn/en) (mainly 
Willow Ptarmigan *Lagopus lagopus*), but can be used for other systems that 
collect age-structured distance sampling data.

The model itself is written and implemented in NIMBLE (see
[here](https://r-nimble.org/) for more information about NIMBLE for R).
The NIMBLE model code is written by a function called
**"writeModelCode_Integ.R"** in the "R" folder. 

Additional R functions used for importing and wrangling the data,
simulating data, preparing data in correct format, setting up and
running the model is contained in the R folder. Most functions have 
roxygen documentation explaining the details for each function. Some functions in this folder remain from the original IDSM, but are not used in this integrated version.

The complete workflow for the analysis for the manuscript with the title: "Capturing three decades of alpine food-web and climate dynamics with an integrated predator–prey model", can be found in the masterscript 
with the name **"Analysis_Integrated_OccProd.R"**. All data used for the analysis is available in the "data" folder, 
which will be automatically called and imported by the functions in the 
masterscript. 

Several functions for extracting, analysing and plotting the model output are available under the "miscellaneous" folder. The analyses of time trends in the gyrfalcon data can be found in a script named "PostHoc_brms_gyrs.R". Finally, in the home folder, a full walk-through of the weather analysis in which we select the 'best' suiting weather variable is available in the quarto document named "PostHoc_WeatherAnalysis.qmd". 

NOTE: This repository is derived from a parent repository that remains under active development. New versions and updates of the parent repository will be released as they become available. All code, data, and supporting information required to reproduce the analyses and results presented in the associated manuscript are provided in this repository. Future updates to the parent repository may introduce additional features, improvements, and functionality beyond those required for this manuscript.

## Additional dependencies

There are two dependencies that need to be manually installed to run
the workflow. First, you need to install NIMBLE (follow instructions
given here: <https://r-nimble.org/download>). Second, the analysis uses
code from the nimbleDistance package
(<https://github.com/scrogster/nimbleDistance>). to estimate the half
normal detection distribution. 

Finally, running the workflow requires access to additional data
(radio-telemetry data on ptarmigan, rodent occupancy data, and
shapefiles for municipalities in Norway). Auxiliary data is now bundled
with the repository, while shapefiles can be downloaded from OSF:
<https://osf.io/7326r/>.
