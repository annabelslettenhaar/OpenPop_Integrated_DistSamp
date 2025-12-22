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
line transect surveys, and rodent trapping surveys from 1990 to 2020. The 
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
As of this first version the NIMBLE model code is written by a function called
"writeModelCode.R" in the "R" folder. 

Additional R functions used for downloading and wrangling the data,
simulating data, preparing data in correct format, setting up and
running the model is contained in the R folder. Most functions have 
roxygen documentation explaining the details for each function.

The complete workflow for the analysis for the manuscript with the title: "Three 
decades of alpine prdator-prey dynamics: an integrated model of gyrfalcons, 
ptarmigan, rodents and climate interactions", can be found in the masterscript 
with the name "Analysis_Integrated_OccProd.R"

All data used for the analysis is available in the "data" folder, 
which will be automatically called and imported by the functions in the 
masterscript.

NOTE: this release is done for the initial submission of the manuscript, 
and still contains redundant scripts and information that should later be 
deleted. A cleaned-up version of this repository will be available upon acceptance 
and publication of the manuscript. Despite that, all the code and relevant 
functions required to reproduce the results reported in this study are present, 
and can be obtained by running "Analysis_Integrated_OccProd.R"

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
