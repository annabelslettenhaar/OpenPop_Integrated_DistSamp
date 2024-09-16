---
editor_options: 
  markdown: 
    wrap: 72
---

[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

# An Integrated Open Population Distance Sampling Model

## What is in this repository?

This repository contains code for a workflow analysing line transect
data using an integrated distance sampling model (IDSM). The model
utilizes age-structured survey data and auxiliary data from marked
individuals to jointly estimate changes in population density and
temporal variation in underlying demographic rates (recruitment rate and
survival probability). It is a multi-area model, meaning it
simultaneously models processes across a defined number of areas, and
shares information both across space and time.

The IDSM workflow is set up as a targets pipeline
(<https://books.ropensci.org/targets/>) and was specifically written for
data collected through the [Norwegian monitoring program for tetraonid
birds](https://honsefugl.nina.no/Innsyn/en) (mainly Willow Ptarmigan
*Lagopus lagopus*), but can be used for other systems that collect
age-structured distance sampling data.

The model itself is written and implemented in NIMBLE (see
[here](https://r-nimble.org/) for more information about NIMBLE for R).
As of v.2.0, the NIMBLE model code is written by a function called
"writeModelCode.R" in the "R" folder. Previous versions relied on NIMBLE
code found in the "NIMBLE code" folder.

Additional R functions used for downloading and wrangling the data,
simulating data, preparing data in correct format, setting up and
running the model, as well as post-hoc analyses and plotting and quality
control of results is contained in the R folder. Refer to each
function's roxygen documentation for detailed documentation.

Complete workflows for analysing simulated data, real data from the
Lierne area only, and data from all areas for which public data on
willow ptarmigan in Norway are available are provided in the following
master scripts:

-   "Analysis_SimData.R" for a single simulated dataset
-   "Analysis_SimData_Replicates.R" for multiple simulated datasets
    including replicate runs
-   "Analysis_RealData_LierneVest.R" for real data from Lierna area only
-   "Analysis_RealData.R" for real data from all areas, manually in R.
-   "\_targets.R" for real data from all area, organized in a targets
    pipeline.
-   "Analysis_RealData_GNUparallel_Setup.R", "rypeIDSM_GNUwrapper.R",
    and "Analysis_RealData_GNUparallel_PostProcessing.R" for real data
    from all areas, to be invoked from terminal (within a Nix shell and
    using GNUparallel).

Note that the real-data (multi-area) workflows have been further
developed in the v2.x releases, i.e. use a newer and reparameterized
version of the model (see "R/writeModelCode.R"). Workflows for simulated
data are fully functional, but use the original parameterisation of the
model (see "NIMBLE Code/"). For these latter workflows, we therefore
recommend using v1.5 of the code.

## Additional dependencies

There are a three dependencies that need to be manually installed to run
the workflow. First, you need to install NIMBLE (follow instructions
given here: <https://r-nimble.org/download>). Second, the analysis uses
code from the nimbleDistance package
(<https://github.com/scrogster/nimbleDistance>). to estimate the half
normal detection distribution. Third, we use the LivingNorwayR package
(<https://livingnorway.github.io/LivingNorwayR/>) to download and
wrangle the line transect distance sampling survey data.

Finally, running the workflow requires access to additional data
(radio-telemetry data on ptarmigan, rodent occupancy data, and
shapefiles for municipalities in Norway). Auxiliary data is now bundled
with the repository, while shapefiles can be downloaded from OSF:
<https://osf.io/7326r/>.

## Citation

Releases of this repository are archived and issued DOIs via Zenodo:
[https://zenodo.org/records/10462269](https://zenodo.org/records/10572340)

The citation of the latest version is:

Chloé R. Nater, ErlendNilsen, christofferhohi, Matthew Grainger, Bernardo Brandão Niebuhr, & Francesco Frassinelli. (2024). 
ErlendNilsen/OpenPop_Integrated_DistSamp: Ptarmigan IDSM v2.1 (v2.1). Zenodo. 
<https://doi.org/10.5281/zenodo.13767267>

## Contact

Erlend Nilsen:
[erlend.nilsen\@nina.no](mailto:erlend.nilsen@nina.no){.email}

Chloé Nater: [chloe.nater\@nina.no](mailto:chloe.nater@nina.no){.email}
