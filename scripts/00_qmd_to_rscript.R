
library("knitr")

purl("processing.qmd", "scripts/script01_processing.R")
purl("modeling.qmd", "scripts/script02_machine_learning.R")
purl("chemometrics.qmd", "scripts/script03_chemometrics.R")
