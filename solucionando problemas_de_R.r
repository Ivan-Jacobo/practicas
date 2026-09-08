# Arreglando R
# Al empezar a usar R prox el 31 de agosto note que no se podian instalar ciertas librerias,
# entre ellas devtools, fs, biostrings, etc., lo cual me parecio raro ya que todas marcaban el  mismo error
# ese error apuntaba a una incompatibilidad de versiones de las loibrerias y sus dependencias con la version de R actual
# 
# sabemos que el porbolmea venia de una rama mayor porque al usar .linParths() me arrojaba las direcciones de los sitios de los cuales cargaba las librerias rstudio
# Estos problemas eran causados porque se estaban cargando las librerias desde una rama superior a la de mi perfil
# por lo que cree un directorio desde mi perfil que contuviera las librerias de R.
#
# Basicamente se encuentra en /home/servicio_social/R/library
#
# y ahi comence a instalar de nuevo fs y usethis con una version compatible al R actual
# despues comnfigure Rstudio para quee cargara las librerias desde ese directorio usando los siguientes comandos
# file.edit(path.expand("~/.Renviron"))
# dentro escribi la direccion de las librerias, lo guarde, reinicie R y continue con mis descargas
#
# fue necesario instalar pkgload
# despues devtools se pudo instalar
#
# inente netcomi pero no se puso por porblemas con spiec easi y pulsar
# instale pulsar mediante devtools desde github
# 
# instale spieceasi desde biocmanager
# INSTALE SPRING DESDE DEVTOOLS GITHUB
#
# ahora estoy instalando netcomi con un codigo extraño que usa dos repositorios:
# devtools::install_github(
# "stefpeschel/NetCoMi",
# lib = "/home/servicio_social/R/library",
# repos = c(
#  "https://cloud.r-project.org/",
#  BiocManager::repositories()
# ),
# build_vignettes = FALSE,
# upgrade = "never"
# )
#