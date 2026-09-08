                  # DADA2 de sujetos sanos #



#                   - Fase 1 -

via_F <- "../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/_1.fastq"
list.files(via_F)

via_R <- "../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/_2.fastq"
list.files(via_R)

# fastq fordwrard -> sort= ordena, vía = ubicación, pattern = patron de busueda.
FqF <- sort(list.files(via_F, pattern="*_1.fastq", full.names = TRUE))

# fastq reverse
FqR <- sort(list.files(via_R, pattern="*_2.fastq", full.names = TRUE))


#                     - Fase 2 -

# Extract sample names, assuming filenames 
# have format: SAMPLENAME_XXX.fastq

#forward
nombres.muestra <- sapply(strsplit(basename(FqF), "_"), `[`, 1)

# reverse
nombres.muestra_2 <- sapply(strsplit(basename(FqR), "_"), `[`, 1)


#                    - Fase 3 -

# Omiti esta parte porque hice las graficas desde la terminal con fastqc

# Graficas de calidad phred para Forward
#pdf("../../../Actividades/Septimo_semestre/Servicio_social/Quality/Quality_Forward_bacterias_ss_27_38.pdf",width=13,height = 8)
#plotQualityProfile(FqF[27:38])
#dev.off()


# Graficas de calidad phred para Revers
#pdf("../../../Actividades/Septimo_semestre/Servicio_social/Quality/Quality_Reverse_bacterias_ss_27_38.pdf",width=13,height = 8)
#plotQualityProfile(FqR[27:38])
#dev.off()


#                     - Fase 4 -

  # Filtrar de acuerdo a la calidad phred, necesito que sea mayor a 30
  # primero hay que generar una carpeta donde se guardarán las secuencias
  # filtradas y recortadas

# Place filtered files in filtered/ subdirectory
filtroFqF <- file.path(via_F, "filtered", paste0(nombres.muestra, "_F_filtro.fastq.gz"))
filtroFqR <- file.path(via_R, "filtered", paste0(nombres.muestra_2, "_R_filtro.fastq.gz"))
names(filtroFqF) <- nombres.muestra
names(filtroFqR) <- nombres.muestra_2

# Filtrando por calidad -> FqF entrada de datos, filtroFqF lugar donde se almacenan 
# los datos recortados. 
salida <- filterAndTrim(FqF, filtroFqF, FqR, filtroFqR, truncLen=c(290,240),
                        maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE) # Forward corte a 280 y reverse a 240

# Nota:
# La calidad de las secuencias reverse es mala, ya que de las 240 pb en adelante la media de la calidad phred baja de 20


# errores forward
errores_F <- learnErrors(filtroFqF, multithread=TRUE)
saveRDS(errores_F,file="../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/_1.fastq/errF.RDS")



# errores en Reverse
errores_R <- learnErrors(filtroFqR, multithread=TRUE)
saveRDS(errores_R,file="../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/_2.fastq/errR.RDS")




#           - Fase 5 -

# Graficas de los errores
png("../../../Actividades/Septimo_semestre/Servicio_social/Errores/errores_F.png")
plotErrors(errores_F, nominalQ=TRUE) # Estan bien, la idea es que los puntos                                    
dev.off()                             # mantengan la dirección de la liena roja
# NO necesariamente tienen que estar acoplados a la perfección.
png("../../../Actividades/Septimo_semestre/Servicio_social/Errores/Errores_R.png")
plotErrors(errores_R, nominalQ=TRUE)
dev.off()

# Usando el algoritmo de inferencia
# filtroFqF -> Les proporciona las secuencias recortadas
# err=errores_F -> Les proporciona los errores correspondientes a las secuencias
# recortadas en Forward.
dadaFs <- dada(filtroFqF, err=errores_F, multithread=TRUE)


dadaRs <- dada(filtroFqR, err=errores_R, multithread=TRUE)



#           - Fase 6 -


# Hay que mezclar ahora "merge".
# denoising = correccion de errores


mergers <- mergePairs(dadaFs, filtroFqF, dadaRs, filtroFqR, verbose=TRUE)
saveRDS(mergers, file="../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/_1.fastq/Forward_reverse_unidos.RDS")
#mergers


tabla_forward_reverse <- makeSequenceTable(mergers)
dim(tabla_forward_reverse)
# Esto no lo correré----
saveRDS(tabla_forward_reverse, file="../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/tabla_forward_reverse.RDS")
#tabla_forward_reverse

#readRDS("../../../Actividades/Septimo_semestre/Servicio_social/tabla_forward_reverse.RDS") -> tabla_forward_reverse
#View(tabla_forward_reverse.nochim)

# Removiendo quimeras
# Nochim = no chimeras = no quimeras
tabla_forward_reverse.nochim <- removeBimeraDenovo(tabla_forward_reverse, method="consensus", multithread=TRUE, verbose=TRUE)
dim(tabla_forward_reverse.nochim)



# falta correr esta linea
sum(tabla_forward_reverse.nochim)/sum(tabla_forward_reverse) # Creo que es para comparar si son diferentes, ya  ue si fueran iguales el resultado sería uno.


# "estadistico de la asignación"
getN <- function(x) sum(getUniques(x))
track <- cbind(salida, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(tabla_forward_reverse.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- nombres.muestra
head(track)


samdf
#           - Fase 7 -

# Asignación taxonomica
taxa <- assignTaxonomy(tabla_forward_reverse.nochim, "../../../silva_nr99_v144_toGenus_trainset.fa.gz", multithread=TRUE)
saveRDS(taxa, file="../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/taxa.RDS")
#writeRDS(file="03_Results/taxa.RDS")

taxa <- addSpecies(taxa, "../../../silva_v144_assignSpecies.fa.gz")
saveRDS(taxa, file="../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/taxa_especies.RDS")
# Muy bonito y todo peor el 100% de las especies salio con NA jajaja
View(taxa)



# Objetivo más simple  (A)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) # Me da una lista con la taxa y su resolución hasta nivel genero
View(taxa.print)
sample_data(taxa.print)


# Aquí incorporare los metadatos
read.csv("01_Raw_data/metadatos_cs_long_sanos.csv") -> metadatos_bacterias
View(metadatos_bacterias)

metadatos_bacterias$Run -> nombres

row.names(metadatos_bacterias) <- nombres
sample_data(metadatos_bacterias)
saveRDS(taxa, file="../../../../../../../media/DiscoDuroExterno/Ivan/Ivan/cushing_longitudinal/Run_descargas/fastq_1_y_2/metadatos_bacterias.RDS")


#       - Fase 8 -



# =======================================================================================
# Forzar para identificar especies (B)

#library(DECIPHER)
#dna <- DNAStringSet(getSequences(tabla_forward_reverse.nochim))
#load("01_RowData/SILVA_SSU_r138_2019.RData")

# Paso más tardado de todos, quiza duro unos 3 min en el proceso
#ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
#ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
#taxid <- t(sapply(ids, function(x) {
#  m <- match(ranks, x$rank)
#  taxa <- x$taxon[m]
#  taxa[startsWith(taxa, "unclassified_")] <- NA
#  taxa
#}))


#colnames(taxid) <- ranks 
#rownames(taxid) <- getSequences(tabla_forward_reverse.nochim)
#head(taxid)
#View(taxid)

# ================================================================================

metadatos_bacterias # metadatos
taxa.print # tax table
tabla_forward_reverse.nochim # otu table

#       -------- Fase 9: Phyloseq ----------------------------


# Data frame con metadatos #


library(phyloseq)
library(Biostrings)
library(ggplot2)


theme_set(theme_bw())
samples.out <- rownames(tabla_forward_reverse.nochim)


subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)


day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


# ------------------------ #
# Phyloseq #
metadatos_bacterias # metadatos
taxa.print # tax table
tabla_forward_reverse.nochim # otu table

class(metadatos_bacterias)
class(taxa.print)
class(tabla_forward_reverse.nochim)

View(tabla_forward_reverse.nochim)

?phyloseq()

rownames(metadatos_bacterias) <- metadatos_bacterias$Run
View(metadatos_bacterias)

ps_cs_lg_sanos <- phyloseq(otu_table(tabla_forward_reverse.nochim, taxa_are_rows=FALSE), 
               sample_data(metadatos_bacterias), 
               tax_table(taxa))

ps_cs_lg_sanos

sample_names(ps_cs_lg_sanos)
taxa_names(ps_cs_lg_sanos)

Biostrings::DNAStringSet(taxa_names(ps_cs_lg_sanos))



ps_cs_lg_sanos <- prune_samples(sample_names(ps_cs_lg_sanos) != "Mock", ps_cs_lg_sanos) # Remove mock sample
dna <- Biostrings::DNAStringSet(taxa_names(ps_cs_lg_sanos))
names(dna) <- taxa_names(ps_cs_lg_sanos)

ps_cs_lg_sanos <- merge_phyloseq(ps_cs_lg_sanos, dna)

taxa_names(ps_cs_lg_sanos) <- paste0("ASV", seq(ntaxa(ps_cs_lg_sanos)))

taxa_names(ps_cs_lg_sanos)

saveRDS(ps_cs_lg_sanos,file="03_Results/ps_cs_lg_sanos.RDS") # modifique el nombre del archivo


#--------------llegue hasta aquí --------#

# RANK abundance ---------------------
# calculo de de la abundancia de cada taxon
ab_tax <- taxa_sums(ps_cs_lg_sanos)

# Abundancias ordenadas de mayor a menor
ab_tax_ord <- sort(ab_tax, decreasing = TRUE)


# gráfica de barras rank-abundance

par(mgp = c(0, 1, -0.5))
barplot(ab_tax_ord, 
        main = "Rank-Abundance", 
        xlab = "Taxones", 
        ylab = "Abundancia total", 
        col = "sienna", 
        las = 2,  
        cex.names = 0.5,
        cex.axis = 0.7)



# =========================================


rel_ab <- transform_sample_counts(ps_cs_lg_sanos, function(x) x / sum(x))


# quitando el borde negro
sample_data(rel_ab)$weight <- as.numeric(sample_data(rel_ab)$weight)

str(sample_data(rel_ab))


plot_bar(rel_ab, x= "gender", fill = "Order") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



df <- psmelt(rel_ab)

ggplot(df, aes(x = weight, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = sort(unique(df$weight))) +  # solo los días presentes
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



df_ejemplo <- psmelt(ps_cs_lg_sanos)  # Convierte a data.frame amigable

# Sumar abundancias por día
library(tidyverse)

df_sum_ej <- df_ejemplo %>%
  group_by(bmi) %>%
  summarise(Total_abundance = sum(Abundance))

# Hacer barplot
ggplot(df_sum_ej, aes(x = factor(bmi), y = Total_abundance, geom_text(0.1))) +
  geom_bar(stat = "identity", fill = "darkorange") +
  theme_bw() +
  labs(x = "Indice de masa corporal", y = "Abundancia total",
       title = "Abundancia total por ims") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) 






df_sum_peso <- df_ejemplo %>%
  group_by(weight) %>%
  summarise(Total_abundance = sum(Abundance))

# Hacer barplot
ggplot(df_sum_peso, aes(x = factor(weight), y = Total_abundance, geom_text(0.1))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  labs(x = "Peso", y = "Abundancia total",
       title = "Abundancia total por peso") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) 



df_sum_altura <- df_ejemplo %>%
  group_by(height) %>%
  summarise(Total_abundance = sum(Abundance))

# Hacer barplot
ggplot(df_sum_altura, aes(x = factor(height), y = Total_abundance, geom_text(0.1))) +
  geom_bar(stat = "identity", fill = "darkred") +
  theme_bw() +
  labs(x = "Estatura", y = "Abundancia total",
       title = "Abundancia total estatura") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) 






