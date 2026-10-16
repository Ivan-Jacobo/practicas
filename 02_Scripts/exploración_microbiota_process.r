# Instalando microbiota processes

## try http:// if https:// URLs are not supported ## the url of mirror
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("MicrobiotaProcess")

library(phyloseq)
library(MicrobiotaProcess)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(tidytree)
library(ggstar)
library(forcats)


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




# Tambien se requieren estas lineas de codigo 

# mouse.time.mpse %<>%
#  mp_cal_abundance( # for each samples
#    .abundance = RareAbundance
#  ) %>%
#  mp_cal_abundance( # for each groups 
#    .abundance=RareAbundance,
#    .group=time
 # )
#mouse.time.mpse


mpse_metagenomico %<>%
  mp_cal_abundance(
    .abundance = RareAbundance
  ) %>% 
  mp_cal_abundance(
    .abundance = RareAbundance,
    .group = Grupo
  )
mpse_metagenomico


# ============================ #

#Graficos de diversidad diferemncial de acuerdo a diferentes indices

mpse_metagenomico %<>%
  mp_cal_alpha(.abundance=RareAbundance)

mpse_metagenomico

#f1 <- mouse.time.mpse %>% 
#  mp_plot_alpha(
#    .group=time, 
#    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
#  ) +
#  scale_fill_manual(values=c("#00A087FF", "#3C5488FF"), guide="none") +
# scale_color_manual(values=c("#00A087FF", "#3C5488FF"), guide="none")

div_1 <- mpse_metagenomico %>% 
  mp_plot_alpha(
    .group = Grupo,
    .alpha = c(Observe, Shannon, Simpson, Pielou)
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF"), guide="none") +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF"), guide="none")


#f2 <- mouse.time.mpse %>%
#  mp_plot_alpha(
#    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
#  )
#f1 / f2

div_2 <- mpse_metagenomico %>% 
  mp_plot_alpha(
    .alpha = c(Observe, Shannon, Simpson, Pielou)
  )

div_1
div_2




#==============================================
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


# =====================================================
# Genación de heatmaps

#h1 <- mouse.time.mpse %>%
 # mp_plot_abundance(
  #  .abundance = RareAbundance,
   # .group = time,
    #taxa.class = Phylum,
    #relative = TRUE,
    #topn = 20,
    #geom = 'heatmap',
    #features.dist = 'euclidean',
   # features.hclust = 'average',
  #  sample.dist = 'bray',
 #   sample.hclust = 'average'
#  )

ht_mp_1 <- mpse_metagenomico %>% 
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Grupo,
    taxa.class = Phylum,
    relative = TRUE,
    topn = 20,
    geom = "heatmap",
    feature.dist = "euclidean",
    feature.hclust = "average",
    sample.dist = "bray",
    sample.hclust = "average"
  )
ht_mp_1


#h2 <- mouse.time.mpse %>%
#  mp_plot_abundance(
#    .abundance = RareAbundance,
#    .group = time,
#    taxa.class = Phylum,
#    relative = FALSE,
#    topn = 20,
#    geom = 'heatmap',
#    features.dist = 'euclidean',
#    features.hclust = 'average',
#    sample.dist = 'bray',
#    sample.hclust = 'average'
#  )

ht_mp_2 <- mpse_metagenomico %>% 
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Grupo,
    taxa.class = Phylum,
    relative = FALSE,
    topn = 15,
    geom = "heatmap",
    feature.dist = "euclidean",
    feature.hclust = "average",
    sample.dist = "bray",
    sample.hclust = "average"
  )
ht_mp_2


# the character (scale or theme) of figure can be adjusted by set_scale_theme
# refer to the mp_plot_dist
#aplot::plot_list(gglist=list(h1, h2), tag_levels="A")

# =============================================#

# Matriz de distancias
  # Esta va a se muy buena para los datos de 16s ya que permite integrar multiples filtrados
  # al grafico


mpse_metagenomico %<>%
  mp_decostand(.abundance=Abundance)


mpse_metagenomico

mpse_metagenomico %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")

graf_dist_1 <- mpse_metagenomico %>%  mp_plot_dist(.distmethod = bray)
graf_dist_1

graf_dist_2 <- mpse_metagenomico %>% mp_plot_dist(.distmethod = bray, .group = Grupo)
# The scale or theme of dot heatmap plot can be adjusted using set_scale_theme function.
graf_dist_2 %>% set_scale_theme(
  x = scale_fill_manual(
    values=c("darkorange", "deepskyblue"), 
    guide = guide_legend(
      keywidth = 0.8, 
      keyheight = 0.3, 
      title.theme = element_text(size=4),
      label.theme = element_text(size=2)
    )
  ), 
  aes_var = Grupo # specific the name of variable 
) %>%
  set_scale_theme(
    x = scale_color_gradient(
      guide = guide_legend(keywidth = 0.3, keyheight = 0.3)
    ),
    aes_var = bray
  ) %>%
  set_scale_theme(
    x = scale_size_continuous(
      range = c(0.1, 1.5),
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  )



# ===============================================================

# Prueba que nos muestra diferencias entre los grupos que tenemos de acuerdo a bray

box_grupos_graf <- mpse_metagenomico %>% mp_plot_dist(.distmethod = bray, .group = Grupo, group.test=TRUE, textsize=2)
box_grupos_graf





# =============================================================== 

# Haceindo un PCoA

mpse_metagenomico %<>%
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")

mpse_metagenomico


mpse_metagenomico %<>%
  mp_adonis(.abundance=hellinger, .formula=~Grupo, distmethod="bray", permutations=10000, action="add")

mpse_metagenomico %>% mp_extract_internal_attr(name=adonis)



componente_1 <- mpse_metagenomico %>%
  mp_plot_ord(
    .ord = pcoa, 
    .group = Grupo, 
    .color = Grupo, 
    .size = 1.2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF")) +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF")) 

componente_1_var <- mpse_metagenomico %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Grupo, 
    .color = Grupo, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#00A087FF", "#3C5488FF"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#00A087FF", "#3C5488FF"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )

componente_1

componente_1_var


# =====================================================================================

mpse_metagenomico %<>%
  mp_cal_clust(
    .abundance = hellinger, 
    distmethod = "bray",
    hclustmethod = "average", # (UPGAE)
    action = "add" # action is used to control which result will be returned
  )

mpse_metagenomico


muestra_clusterizacion <- mpse_metagenomico %>% mp_extract_internal_attr(name='SampleClust')
muestra_clusterizacion


arbol_1 <- ggtree(muestra_clusterizacion) + 
  geom_tippoint(aes(color=Grupo)) +
  geom_tiplab(as_ylab = TRUE) +
  ggplot2::scale_x_continuous(expand=c(0, 0.01))

arbol_1
# los grupos no se ordenaron de manera adecuada, 
# existen algunos sujetos que si se alinearon dentrpo de su grupo, pero hay
# muchos otros sujetos que se agrupan de manera aleatoria dentro del cladograma



# =============================================================================

ab_filo <- mpse_metagenomico %>% 
  mp_extract_abundance(taxa.class=Phylum, topn=15)



ab_filo %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phyla="label")
ab_filo




barras_cladograma <- arbol_1 + 
  geom_fruit(
    data=ab_filo,
    geom=geom_col,
    mapping = aes(x = RelRareAbundanceBySample, 
                  y = Sample, 
                  fill = Phyla
    ),
    orientation = "y",
    #offset = 0.4,
    pwidth = 3, 
    axis.params = list(axis = "x", 
                       title = "The relative abundance of phyla (%)",
                       title.size = 4,
                       text.size = 2, 
                       vjust = 1),
    grid.params = list()
  )

barras_cladograma
# Hay un cluster de individuos con cushing que tiene una mayor abundancia de bacillota y sctinomycetota, uroviricota y verrrucomicrobiota que el resto de muestras

# ==============================================================================================================


mpse_metagenomico %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = Grupo,
    first.test.alpha = 0.01
  )


taxa.tree <- mpse_metagenomico  %>% 
  mp_extract_tree(type="taxatree")


taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Grupo, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))
taxa.tree

#saveRDS(mpse_metagenomico, "03_Results/objetos_mpse/mpse_metagenomico.RDS")

readRDS("03_Results/objetos_mpse/mpse_metagenomico.RDS") -> mpse_metagenomico

mpse_metagenomico

# ===============================================================================

# Creando el super arbol
is.na(Sign_time)
# Since taxa.tree is treedata object, it can be visualized by ggtree and ggtreeExtra
p1 <- ggtree(
  taxa.tree,
  layout="radial",
  size = 0.3
) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21
  )

# display the high light of phylum clade.
p2 <- p1 +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label)
  )


# display the relative abundance of features(OTU)
p3 <- p2 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, Grupo, .fun=min),
      size = RelRareAbundanceBySample,
      fill = Grupo,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))


# display the tip labels of taxa tree
p4 <- p3 + geom_tiplab(size=2, offset=7.2)



# display the LDA of significant OTU.
p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(
      x = LDAmean,
      fill = Sign_Grupo,
      subset = !is.na(LDAmean)
    ),
    orientation = "y",
    offset = 0.3,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2)
  )


# display the significant (FDR) taxonomy after kruskal.test (default)
p6 <- p5 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_Grupo)),
    mapping = aes(size = -log10(fdr),
                  fill = Sign_Grupo,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))

# A partir de aqui se traba todo, debe de ser un grafico MUY pesado
#p6 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)


mpse_metagenomico
# ===============================================================



flag <- packageVersion("ggnewscale") >= "0.5.0"
new.fill <- ifelse(flag, "fill_ggnewscale_1", "fill_new_new")
new.fill2 <- ifelse(flag , "fill_ggnewscale_2", "fill_new")

w <- mpse_metagenomico %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.1
  ) 
# if version of `ggnewscale` is >= 0.5.0, you can also use p$ggnewscale to view the renamed scales.
w <- w  +
  scale_fill_manual(values=c("deepskyblue", "orange")) +
  scale_fill_manual(
    aesthetics = new.fill2, # The fill aes was renamed to `new.fill` for the abundance dotplot layer
    values = c("deepskyblue", "orange")
  ) +
  scale_fill_manual(
    aesthetics = new.fill, # The fill aes for hight light layer of tree was renamed to `new.fill2`
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999"
    )
  )
w





# ==============================================

cladograma_taxones_dif <- mpse_metagenomico %>%
  mp_plot_diff_cladogram(
    label.size = 2.5,
    hilight.alpha = .3,
    bg.tree.size = .5,
    bg.point.size = 2,
    bg.point.stroke = .25
  ) +
  scale_fill_diff_cladogram( # set the color of different group.
    values = c('deepskyblue', 'orange')
  ) +
  scale_size_continuous(range = c(1, 4))

cladograma_taxones_dif



mpse_metagenomico %>% 
  mp_plo



#============================================#
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




# =================================================


grafico.mahattan <- mpse_metagenomico %>%
  mp_plot_diff_manhattan(
    .group = Sign_Grupo,
    .y = fdr,
    .size = 2.4,
    taxa.class = c('OTU', 'Genus'),
    anno.taxa.class = Phylum
  )
grafico.mahattan

mpse_metagenomico
