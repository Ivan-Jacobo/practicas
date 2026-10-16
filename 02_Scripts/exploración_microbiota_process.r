# Instalando microbiota processes

## try http:// if https:// URLs are not supported ## the url of mirror
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("MicrobiotaProcess")

library(phyloseq)
library(MicrobiotaProcess)

# FUnciono a la primera, ahora seguire el tutoria para pasar de un phyloseq a un mpse

# convert phyloseq object to mpse
#data(esophagus)
#esophagus
#mpse4 <- esophagus %>% as.MPSE() 
#mpse4
# convert TreeSummarizedExperiment object to mpse
# library(curatedMetagenomicData)
# tse <- curatedMetagenomicData::curatedMetagenomicData("ZhuF_2020.relative_abundance", dryrun=F)
# tse[[1]] %>% as.MPSE() -> mpse5
# mpse5


readRDS("03_Results/objetos_phyloseq/ps_metagenomico_variante.RDS") -> ps_metagenomico_variante
mpse_metagenomico <- ps_metagenomico_variante %>% as.MPSE()
mpse_metagenomico


# primeros pasos -----------------
# Son relevantes para que se cree la columna de rare abundance en el objeto mpse.
# Esta columna es requerida para analisis posteriores.

mpse_metagenomico %<>% mp_rrarefy()


mpse_metagenomico %<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 400
  )

mpse_metagenomico %>% print(width=180)



# Graficos de barras apiladas relativas
# p1 <- mouse.time.mpse %>%
#  mp_plot_abundance(
#    .abundance=RareAbundance,
#    .group=time, 
#    taxa.class = Phylum, 
#    topn = 20,
#    relative = TRUE
#  )




GBA_rel_1 <- mpse_metagenomico %>% 
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Grupo,
    taxa.class = Phylum,
    topn = 20,
    relative = TRUE
  )


GBA_rel_2 <- mpse_metagenomico %>% 
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Grupo,
    taxa.class = Phylum,
    topn = 20,
    relative = FALSE
  )


GBA_rel_1




# Codigo para hacer los boxplot de abundancias diferenciales
# f.box <- mouse.time.mpse %>%
#  mp_plot_diff_boxplot(
#    .group = time,
#  ) %>%
#  set_diff_boxplot_color(
#    values = c("deepskyblue", "orange"),
#    guide = guide_legend(title=NULL)
#  )

abdif_box <- mpse_metagenomico %>% 
  mp_plot_diff_boxplot(
    .group = Grupo,
  ) %>% 
  set_diff_boxplot_color(
    values = c("deepskyblue", "darkorange"),
    guide = guide_legend(title = NULL)
  )



# f.bar <- mouse.time.mpse %>%
#  mp_plot_diff_boxplot(
#    taxa.class = c(Genus, OTU), # select the taxonomy level to display
#    group.abun = TRUE, # display the mean abundance of each group
#    removeUnknown = TRUE, # whether mask the unknown taxa.
#  ) %>%
#  set_diff_boxplot_color(
#    values = c("deepskyblue", "orange"),
#    guide = guide_legend(title=NULL)
#  )


abdif_bar <- mpse_metagenomico %>% 
  mp_plot_diff_boxplot(
    taxa.class = c(Genus, OTU),
    group.abun = TRUE,
    removeUnknown = TRUE
  ) %>% 
  set_diff_boxplot_color(
    values = c("deepskyblue", "darkorange"),
    guide = guide_legend(title= NULL)
  )


aplot::plot_list(abdif_box, abdif_bar)




