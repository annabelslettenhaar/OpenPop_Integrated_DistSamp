install.packages("rix", repos = c(
  "https://b-rodrigues.r-universe.dev",
  "https://cloud.r-project.org"
))

library("rix")

rix(
    r_ver = "4.4.0",
    r_pkgs = c(
        "coda", 
        "codetools",
        "colorspace", 
        "cowplot", 
        "EnvStats", 
        "extraDistr", 
        "ggforce", 
        "ggplot2", 
        "ggpubr", 
        "MCMCvis", 
        "paletteer", 
        "popbio",
        "qs", 
        "Rage",
        "reshape2",
        "RJSONIO", # https://github.com/LivingNorway/LivingNorwayR/issues/75
        "sf", 
        "see", 
        "terra", 
        "tidybayes", 
        "tidyverse",
        "tmap", 
        "viridis"
    ),
    system_pkgs = c(
        "parallel"
    ),
    git_pkgs = list(
        list(
            package_name = "LivingNorwayR",
            repo_url = "https://github.com/LivingNorway/LivingNorwayR",
            commit = "ef83ae20f0be20c922a56bd4fc0053ec30d8b4a3"
        ),
        list(
            package_name = "nimbleDistance",
            repo_url = "https://github.com/scrogster/nimbleDistance",
            commit = "d5ea9c15dc484a7923a60a19e36f5ed5a107ca4b"
        )
    ),
    ide = "rstudio",
    project_path = ".",
    overwrite = TRUE,
    print = TRUE
)