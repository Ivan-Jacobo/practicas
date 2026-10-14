# =====================================================================
# CONSTRUCCIÓN DE UN OBJETO PHYLOSEQ A PARTIR DE KRAKEN2 + BRACKEN
# Proyecto: Controles (n = 50) vs. Síndrome de Cushing (n = 50)
# =====================================================================
#
# ¿QUÉ HACE ESTE SCRIPT?
# Un objeto phyloseq necesita TRES tablas que "hablen" entre sí:
#
#   1. otu_table    -> ¿Cuántas lecturas hay de cada taxón en cada muestra?
#                      (viene de los archivos .bracken, columna new_est_reads)
#   2. tax_table    -> ¿A qué linaje pertenece cada taxón? (Dominio, Filo, ...)
#                      (se reconstruye a partir del archivo .KREPORT)
#   3. sample_data  -> ¿Qué es cada muestra? (Control o Cushing)
#
# Las tres se unen mediante el IDENTIFICADOR TAXONÓMICO (taxid de NCBI)
# y el NOMBRE DE LA MUESTRA.
#
# NOTA SOBRE LOS DATOS:
# En el .KREPORT combinado, S1-S50 = Controles 1-50 y S51-S100 = Cushing 1-50.
# Aquí solo usamos el .KREPORT para sacar la TAXONOMÍA (no las abundancias),
# por lo que esa equivalencia no se necesita. Las abundancias vienen de
# Bracken, donde los archivos se llaman Control1...Control50 y CS1...CS50.
# =====================================================================


# ---------------------------------------------------------------------
# 0. PAQUETES Y RUTAS
# ---------------------------------------------------------------------

library(dplyr)     # manipular tablas: select, filter, left_join, ...
library(tidyr)     # cambiar el formato de tablas: pivot_wider
library(tibble)    # column_to_rownames
library(phyloseq)  # (Bioconductor) crea y maneja el objeto phyloseq

# Si phyloseq no está instalado:
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("phyloseq")

# Definir las rutas al inicio: si algo cambia de lugar, solo se edita aquí.
ruta_kreport <- "01_Raw_data/Controles_Cushing.KREPORT"
dir_bracken  <- "01_Raw_data/Bracken"
ruta_salida  <- "03_Results/ps_metagenomico.RDS"


# ---------------------------------------------------------------------
# 1. TABLA DE ABUNDANCIAS (desde los archivos .bracken)
# ---------------------------------------------------------------------
# Cada archivo .bracken tiene una fila por taxón y, entre otras, estas columnas:
#   taxonomy_id   -> identificador NCBI del taxón (lo usamos como "llave")
#   new_est_reads -> lecturas estimadas por Bracken (nuestra abundancia)
# Necesitamos juntar los 100 archivos en UNA matriz: taxones x muestras.

# 1.1 Listar los archivos. El patrón (expresión regular) significa:
#     empieza con "Control" o "CS", seguido de números, y termina en ".bracken"
archivos <- list.files(
  path       = dir_bracken,
  pattern    = "^(Control|CS)[0-9]+\\.bracken$",
  full.names = TRUE
)

# Comprobación de seguridad: deben ser 100 archivos
stopifnot(length(archivos) == 100)

# 1.2 El nombre de cada muestra es el nombre del archivo sin ".bracken"
nombres_muestras <- sub("\\.bracken$", "", basename(archivos))

# 1.3 Leer cada archivo y quedarnos solo con las 2 columnas útiles
lista_bracken <- lapply(archivos, function(f) {
  read.delim(f) %>%
    select(taxonomy_id, new_est_reads)
})
names(lista_bracken) <- nombres_muestras

# 1.4 Apilar las 100 tablas en una sola (formato "largo").
#     .id = "muestra" crea una columna con el nombre de la muestra de origen.
abundancias_larga <- bind_rows(lista_bracken, .id = "muestra")

# 1.5 Pasar a formato "ancho": una fila por taxón, una columna por muestra.
#     values_fill = 0 -> si un taxón no aparece en una muestra, su abundancia
#     es 0 (no "NA"). Esto reemplaza el full_join + replace_na del script
#     original en un solo paso.
abundancias <- abundancias_larga %>%
  pivot_wider(
    names_from  = muestra,
    values_from = new_est_reads,
    values_fill = 0
  ) %>%
  mutate(taxonomy_id = as.character(taxonomy_id))  # como texto, para poder unir luego

dim(abundancias)   # filas = taxones; columnas = 1 (taxonomy_id) + 100 muestras


# ---------------------------------------------------------------------
# 2. TABLA TAXONÓMICA (desde el .KREPORT)
# ---------------------------------------------------------------------
# El .KREPORT es un árbol: cada taxón aparece debajo de su "padre".
# Columnas que usamos:
#   lvl_type -> rango del taxón (D = dominio, K = reino, P = filo, C = clase,
#               O = orden, F = familia, G = género, S = especie)
#   taxid    -> identificador NCBI (equivale a taxonomy_id de Bracken)
#   name     -> nombre, con espacios al inicio que indican la profundidad
# Objetivo: para cada taxón obtener su linaje completo en columnas.

# skip = 102 salta las líneas de encabezado propias de este archivo combinado
kreport <- read.delim(
  ruta_kreport,
  skip        = 102,
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE
)

# Revisar qué códigos de rango existen realmente en TU archivo.
# Kraken2 usa "D" para dominio (Bacteria, Archaea, Eukaryota, Viruses).
# Además existen códigos intermedios como "G1", "S1", "R", "U" (se ignoran).
table(kreport$lvl_type)

# 2.1 Quitar los espacios de sangría de los nombres
kreport <- kreport %>%
  mutate(name = trimws(name))

# 2.2 Rangos que queremos conservar y cómo se llamarán las columnas
rangos  <- c( "K",       "P",      "C",     "O",     "F",      "G",     "S")
niveles <- c( "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# 2.3 Reconstruir el linaje.
# Truco: el .KREPORT está ordenado "de padres a hijos". Por eso basta con ir
# recorriéndolo de arriba a abajo y recordar el último nombre visto en cada
# rango. Cuando aparece un taxón de un rango, se BORRAN los rangos inferiores
# (porque pertenecen a otra rama del árbol).
#
# Ejemplo: al leer "Escherichia" (G), se borra la especie anterior; al leer
# "Escherichia coli" (S), el género "Escherichia" sigue vigente.

linaje_actual <- rep(NA_character_, length(rangos))   # el linaje "en curso"
tax_matriz <- matrix(
  NA_character_,
  nrow = nrow(kreport), ncol = length(rangos),
  dimnames = list(NULL, niveles)
)

for (i in seq_len(nrow(kreport))) {

  pos <- match(kreport$lvl_type[i], rangos)   # posición del rango (NA si no nos interesa)

  if (!is.na(pos)) {
    linaje_actual[pos] <- kreport$name[i]                     # guardar este taxón
    if (pos < length(rangos)) {
      linaje_actual[(pos + 1):length(rangos)] <- NA           # borrar rangos inferiores
    }
    tax_matriz[i, ] <- linaje_actual                          # copiar el linaje a esta fila
  }
}

# 2.4 Armar la tabla taxonómica y quedarnos solo con filas de rangos válidos
tabla_tax <- data.frame(
  taxid    = as.character(kreport$taxid),
  lvl_type = kreport$lvl_type,
  tax_matriz
) %>%
  filter(lvl_type %in% rangos) %>%
  select(-lvl_type) %>%
  distinct(taxid, .keep_all = TRUE)   # por seguridad: un solo renglón por taxid

head(tabla_tax)


# ---------------------------------------------------------------------
# 3. UNIR ABUNDANCIAS Y TAXONOMÍA
# ---------------------------------------------------------------------
# Todo taxón de Bracken debe existir en el .KREPORT. Si esto falla, los
# archivos no provienen de la misma corrida de Kraken2.
stopifnot(all(abundancias$taxonomy_id %in% tabla_tax$taxid))

# Unimos por taxid. Así abundancias y taxonomía quedan alineadas por
# construcción (misma fila = mismo taxón) y no hace falta verificarlo después.
datos <- abundancias %>%
  left_join(tabla_tax, by = c("taxonomy_id" = "taxid"))

dim(datos)
head(datos[, c("taxonomy_id", "Genus", "Species")])


# ---------------------------------------------------------------------
# 4. LAS TRES TABLAS DE PHYLOSEQ
# ---------------------------------------------------------------------

# 4.1 otu_table: matriz numérica (filas = taxones, columnas = muestras)
otu_mat <- datos %>%
  select(taxonomy_id, all_of(nombres_muestras)) %>%
  column_to_rownames("taxonomy_id") %>%    # el taxid pasa a ser el nombre de fila
  as.matrix()

# Comprobaciones: sin NA y sin valores negativos
stopifnot(!anyNA(otu_mat), all(otu_mat >= 0))

# 4.2 tax_table: matriz de texto (filas = taxones, columnas = rangos)
tax_mat <- datos %>%
  select(taxonomy_id, all_of(niveles)) %>%
  column_to_rownames("taxonomy_id") %>%
  as.matrix()

# 4.3 sample_data: una fila por muestra. El grupo se deduce del nombre.
metadatos <- data.frame(
  Grupo = factor(
    ifelse(grepl("^Control", nombres_muestras), "Control", "Cushing"),
    levels = c("Control", "Cushing")   # "Control" = referencia en análisis posteriores
  ),
  row.names = nombres_muestras
)

table(metadatos$Grupo)   # debe dar 50 y 50


# ---------------------------------------------------------------------
# 5. ENSAMBLAR Y GUARDAR EL OBJETO PHYLOSEQ
# ---------------------------------------------------------------------
# phyloseq empareja las tablas por nombre de fila/columna:
#   - taxones: rownames de otu_mat <-> rownames de tax_mat
#   - muestras: colnames de otu_mat <-> rownames de metadatos

ps_metagenomico_variante <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),   # TRUE: los taxones están en las filas
  tax_table(tax_mat),
  sample_data(metadatos)
)

ps_metagenomico
ps_metagenomico_variante   # debe mostrar: N taxa, 100 samples, 8 taxonomic ranks

# Revisión rápida de calidad: lecturas totales por muestra
summary(sample_sums(ps_metagenomico))

saveRDS(ps_metagenomico_variante, "03_Results/ps_metagenomico_variante.RDS")

