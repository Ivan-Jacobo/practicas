                  # DADA2 de datos de bacterias #

library(dada2)
library(phyloseq)
hh

#                   - Fase 1 -

via_F <- "../../../Actividades/Septimo_semestre/Servicio_social/bacterias/Fordward/"
list.files(via_F)

via_R <- "../../../Actividades/Septimo_semestre/Servicio_social/bacterias/Reverse/"
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

# Graficas de calidad phred para Forward
pdf("../../../Actividades/Septimo_semestre/Servicio_social/Quality/Quality_Forward_bacterias_ss_27_38.pdf",width=13,height = 8)
plotQualityProfile(FqF[27:38])
dev.off()


# Graficas de calidad phred para Revers
pdf("../../../Actividades/Septimo_semestre/Servicio_social/Quality/Quality_Reverse_bacterias_ss_27_38.pdf",width=13,height = 8)
plotQualityProfile(FqR[27:38])
dev.off()


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
salida <- filterAndTrim(FqF, filtroFqF, FqR, filtroFqR, truncLen=c(280,240),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE) # Forward corte a 280 y reverse a 240

# errores forward
errores_F <- learnErrors(filtroFqF, multithread=TRUE)
saveRDS(errores_F,file="../../../Actividades/Septimo_semestre/Servicio_social/errF.RDS")



# errores en Reverse
errores_R <- learnErrors(filtroFqR, multithread=TRUE)
saveRDS(errores_R,file="../../../Actividades/Septimo_semestre/Servicio_social/errR.RDS")




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

#| echo: true
mergers <- mergePairs(dadaFs, filtroFqF, dadaRs, filtroFqR, verbose=TRUE)
saveRDS(mergers, file="../../../Actividades/Septimo_semestre/Servicio_social/Forward_reverse_unidos.RDS")
mergers


tabla_forward_reverse <- makeSequenceTable(mergers)
dim(tabla_forward_reverse)
# Esto no lo correré----
saveRDS(tabla_forward_reverse, file="../../../Actividades/Septimo_semestre/Servicio_social/tabla_forward_reverse.RDS")
tabla_forward_reverse

readRDS("../../../Actividades/Septimo_semestre/Servicio_social/tabla_forward_reverse.RDS") -> tabla_forward_reverse
View(tabla_forward_reverse.nochim)

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
taxa <- assignTaxonomy(tabla_forward_reverse.nochim, "../../Sexto semestre/Robertongo/DADA2/DADA_2/01_RowData/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
saveRDS(taxa, file="../../../Actividades/Septimo_semestre/Servicio_social/taxa.RDS")
#writeRDS(file="03_Results/taxa.RDS")


# Nota: en el script original me habia equivocado en el archivo que proporcione
# para llegar a la asignación de especie. AQUI la corregí.
taxa <- addSpecies(taxa, "../../Sexto semestre/Robertongo/DADA2/DADA_2/01_RowData/silva_v138.2_assignSpecies.fa.gz")
saveRDS(taxa, file="../../../Actividades/Septimo_semestre/Servicio_social/taxa_especies.RDS")
View(taxa)
readRDS("../../../Actividades/Septimo_semestre/Servicio_social/taxa_especies.RDS") -> taxa

# Objetivo más simple  (A)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) # Me da una lista con la taxa y su resolución hasta nivel genero
View(taxa.print)
sample_data(taxa.print)


# Aquí incorporare los metadatos
read.csv("../../../Actividades/Septimo_semestre/Servicio_social/metadatos_bacterias.csv") -> metadatos_bacterias
View(metadatos_bacterias)

metadatos_bacterias$Run -> nombres

row.names(metadatos_bacterias) <- nombres
sample_data(metadatos_bacterias)
saveRDS(taxa, file="../../../Actividades/Septimo_semestre/Servicio_social/metadatos_bacterias.RDS")


#       - Fase 8 -




# Forzar para identificar especies (B)

library(DECIPHER)
dna <- DNAStringSet(getSequences(tabla_forward_reverse.nochim))
load("01_RowData/SILVA_SSU_r138_2019.RData")


# Paso más tardado de todos, quiza duro unos 3 min en el proceso
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)


ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")


taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))


colnames(taxid) <- ranks 
rownames(taxid) <- getSequences(tabla_forward_reverse.nochim)

head(taxid)
View(taxid)


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

View(taxa)


ps <- phyloseq(otu_table(tabla_forward_reverse.nochim, taxa_are_rows=FALSE), 
               sample_data(metadatos_bacterias), 
               tax_table(taxa))

ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))

names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps

save(ps,file="03_Results/phyloseq_result.RDS") # modifique el nombre del archivo

View(ps)
#--------------llegue hasta aquí --------#

View(sample_data(ps))
View(otu_table(ps))
View(tax_table(ps))

plot_bar(ps, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

ps -> phyloseq_bacterias

saveRDS(phyloseq_bacterias, "01_raw_data/Bacterias_hongos_virus/phyloseq_bacterias.RDS")






