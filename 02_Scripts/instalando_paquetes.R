# Registro de instalaciones realizadas.

library(BiocManager)

install.packages("BiocManager",
                 lib = "/home/servicio_social/R/library",
                 repos = "https://cloud.r-project.org")

BiocManager::install(version = "3.23", ask = FALSE)
library(BiocManager)
BiocManager::version()

library(devtools)


BiocManager::install(
  c("SpiecEasi", "pulsar"),
  lib = "/home/servicio_social/R/library",
  ask = FALSE,
  update = FALSE
)


devtools::install_github(
  "GraceYoon/SPRING",
  lib = "/home/servicio_social/R/library",
  dependences = TRUE,
  upgrade = "never"
)


install.packages(
  "pulsar",
  lib = "/home/servicio_social/R/library",
  dependencies = TRUE,
  repos = "https://cloud.r-project.org"
)

devtools::install_github(
  "zdk123/pulsar",
  lib = "/home/servicio_social/R/library",
  dependencies = TRUE,
  upgrade = "never"
)

library(pulsar)

BiocManager::install(
  "SpiecEasi",
  lib = "/home/servicio_social/R/library",
  ask = FALSE,
  update = FALSE
)

library(SpiecEasi)


devtools::install_github(
  "GraceYoon/SPRING",
  lib = "/home/servicio_social/R/library",
  dependencies = TRUE,
  upgrade = "never"
)

library(SPRING)


devtools::install_github(
  "stefpeschel/NetCoMi",
  lib = "/home/servicio_social/R/library",
  repos = c(
    "https://cloud.r-project.org/",
    BiocManager::repositories()
  ),
  build_vignettes = FALSE,
  upgrade = "never"
)

library(NetCoMi)
library(phyloseq)
library(dada2)



BiocManager::install(
  "dada2",
  lib = "/home/servicio_social/R/library",
  ask = FALSE,
  update = FALSE
)

library(dada2)
