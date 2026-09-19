
##############################################################
# COFFEE RNA-seq COMPLETE ANALYSIS
# Cultivar × Treatment
#
# Catuai vs CR95
# saline vs Xylella
#
# Includes:
# 1. Data loading
# 2. DESeq2 factorial analysis
# 3. All biologically relevant contrasts
# 4. Volcano plots
# 5. MA plots
# 6. CAMERA GO enrichment
# 7. Hypergeometric GO enrichment
# 8. Gene Venn diagram
# 9. GO Venn diagram
# 10. Heatmaps with fixed biological column order
# 11. Interaction gene profiles
# 12. WGCNA
# 13. Cytoscape export
##############################################################


############################
# 0. LIBRARIES
############################

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "ggplot2",
  "VennDiagram",
  "gridExtra",
  "pheatmap",
  "dendsort"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "DESeq2",
  "edgeR",
  "limma",
  "GO.db",
  "AnnotationDbi",
  "WGCNA"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(VennDiagram)
library(grid)
library(gridExtra)
library(pheatmap)
library(dendsort)
library(GO.db)
library(AnnotationDbi)
library(WGCNA)

allowWGCNAThreads()


############################
# 1. WORKING DIRECTORY
############################

setwd(choose.dir())

outpath <- "salido_complete"
dir.create(outpath, showWarnings = FALSE)


############################
# 2. LOAD COUNTS
############################

counts <- read.table(
  "Matrix_hisat2.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)


# Clean gene IDs
rownames(counts) <- trimws(rownames(counts))

# Remove duplicated gene IDs
counts <- counts[!duplicated(rownames(counts)), ]

# Remove zero-count genes
counts <- counts[rowSums(counts) > 0, ]


############################
# 3. LOAD METADATA
############################

metadata <- read.table(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$sample

# Clean sample names
colnames(counts) <- trimws(colnames(counts))
metadata$sample <- trimws(metadata$sample)

# Keep only matching samples
common.samples <- intersect(
  colnames(counts),
  metadata$sample
)

counts <- counts[, common.samples, drop = FALSE]

metadata <- metadata[
  match(common.samples, metadata$sample),
  ,
  drop = FALSE
]

rownames(metadata) <- metadata$sample

stopifnot(
  all(colnames(counts) == rownames(metadata))
)


############################
# 4. FACTORS
############################

metadata$cultivar <- factor(
  metadata$cultivar,
  levels = c("Catuai", "CR95")
)

metadata$treatment <- factor(
  metadata$treatment,
  levels = c("saline", "xylella")
)

stopifnot(!any(is.na(metadata$cultivar)))
stopifnot(!any(is.na(metadata$treatment)))


############################
# 5. FIXED BIOLOGICAL SAMPLE ORDER
#
# IMPORTANT:
#
# Catuai saline
# Catuai xylella
# CR95 saline
# CR95 xylella
#
# This prevents hierarchical clustering from rearranging
# the biological groups.
############################

metadata$group <- interaction(
  metadata$cultivar,
  metadata$treatment,
  sep = "_"
)

desired.order <- c(
  "Catuai_saline",
  "Catuai_xylella",
  "CR95_saline",
  "CR95_xylella"
)

sample.order <- rownames(
  metadata[
    match(
      desired.order,
      metadata$group
    ),
    ,
    drop = FALSE
  ]
)

sample.order <- sample.order[
  !is.na(sample.order)
]

# If duplicated/replicate groups exist, order all samples
sample.order <- rownames(
  metadata[
    order(
      factor(
        metadata$group,
        levels = desired.order
      )
    ),
    ,
    drop = FALSE
  ]
)

counts <- counts[
  ,
  sample.order,
  drop = FALSE
]

metadata <- metadata[
  sample.order,
  ,
  drop = FALSE
]


############################
# 6. DESEQ2 OBJECT
############################

coldata <- metadata[
  ,
  c("cultivar", "treatment"),
  drop = FALSE
]

dds.full <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ cultivar + treatment + cultivar:treatment
)


############################
# 7. EXPRESSION FILTERING
############################

cpm <- counts(dds.full, normalized = FALSE)

group <- interaction(
  coldata$cultivar,
  coldata$treatment,
  drop = TRUE
)

keep <- rep(FALSE, nrow(cpm))

for (g in levels(group)) {
  
  idx <- which(group == g)
  
  # At least half of samples in the group
  cutoff <- ceiling(length(idx) / 2)
  
  keep <- keep |
    (
      rowSums(
        cpm[, idx, drop = FALSE] >= 10
      ) >= cutoff
    )
}

dds <- dds.full[keep, ]

cat(
  "\nGenes retained after filtering:",
  nrow(dds),
  "\n"
)


############################
# 8. RUN DESEQ2
############################

dds <- DESeq(dds)

cat("\nDESeq2 coefficients:\n")
print(resultsNames(dds))


############################
# 9. VST
############################

vst.obj <- vst(
  dds,
  blind = TRUE
)

vst.mat <- assay(vst.obj)

# Reorder VST columns
vst.mat <- vst.mat[
  ,
  sample.order,
  drop = FALSE
]


############################
# 10. SAMPLE ANNOTATION
############################

annotation.samples <- data.frame(
  Cultivar = metadata$cultivar,
  Treatment = metadata$treatment
)

rownames(annotation.samples) <- rownames(metadata)


############################
# 11. PCA
############################

pca <- prcomp(
  t(vst.mat),
  center = TRUE,
  scale. = FALSE
)

percentVar <- 100 *
  pca$sdev^2 /
  sum(pca$sdev^2)

pca.df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  cultivar = metadata[
    rownames(pca$x),
    "cultivar"
  ],
  treatment = metadata[
    rownames(pca$x),
    "treatment"
  ]
)

pca.plot <- ggplot(
  pca.df,
  aes(
    x = PC1,
    y = PC2,
    color = treatment,
    shape = cultivar
  )
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  labs(
    title = "PCA: Coffee cultivars × treatment",
    x = paste0(
      "PC1 (",
      round(percentVar[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(percentVar[2], 1),
      "%)"
    )
  ) +
  theme_minimal(base_size = 14)

print(pca.plot)

ggsave(
  file.path(
    outpath,
    "PCA_cultivar_treatment.png"
  ),
  pca.plot,
  width = 8,
  height = 6,
  dpi = 300
)


############################
# 12. SAMPLE DISTANCE HEATMAP
############################

sample.dist <- dist(
  t(vst.mat)
)

pheatmap(
  as.matrix(sample.dist),
  annotation_col = annotation.samples,
  annotation_row = annotation.samples,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "Sample-to-sample distance"
)


############################
# 13. SAMPLE CORRELATION HEATMAP
############################

sample.cor <- cor(vst.mat)

pheatmap(
  sample.cor,
  annotation_col = annotation.samples,
  annotation_row = annotation.samples,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "Sample expression correlation"
)


##############################################################
# 14. ALL DESEQ2 CONTRASTS
##############################################################

res.list <- list()


############################
# A. Catuai treatment
############################

res.list$Catuai <- results(
  dds,
  contrast = c(
    "treatment",
    "xylella",
    "saline"
  ),
  alpha = 0.05
)


############################
# B. CR95 treatment
############################

interaction.coef <- grep(
  "cultivar.*CR95.*treatment.*xylella|treatment.*xylella.*cultivar.*CR95",
  resultsNames(dds),
  value = TRUE
)

if (length(interaction.coef) != 1) {
  
  stop(
    "Could not uniquely identify interaction coefficient. ",
    "Check resultsNames(dds)."
  )
}

res.list$CR95 <- results(
  dds,
  contrast = list(
    c(
      "treatment_xylella_vs_saline",
      interaction.coef
    )
  ),
  alpha = 0.05
)


############################
# C. Cultivar under saline
############################

cultivar.coef <- grep(
  "^cultivar.*CR95.*Catuai",
  resultsNames(dds),
  value = TRUE
)

if (length(cultivar.coef) != 1) {
  
  stop(
    "Could not uniquely identify cultivar coefficient. ",
    "Check resultsNames(dds)."
  )
}

res.list$CR95_vs_Catuai_saline <- results(
  dds,
  contrast = c(
    "cultivar",
    "CR95",
    "Catuai"
  ),
  alpha = 0.05
)


############################
# D. Cultivar under xylella
############################

res.list$CR95_vs_Catuai_xylella <- results(
  dds,
  contrast = list(
    c(
      cultivar.coef,
      interaction.coef
    )
  ),
  alpha = 0.05
)


############################
# E. Main treatment effect
############################

res.list$Treatment_main <- results(
  dds,
  contrast = c(
    "treatment",
    "xylella",
    "saline"
  ),
  alpha = 0.05
)


############################
# F. Main cultivar effect
############################

res.list$Cultivar_main <- results(
  dds,
  contrast = c(
    "cultivar",
    "CR95",
    "Catuai"
  ),
  alpha = 0.05
)


############################
# G. Interaction
############################

dds.LRT <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ cultivar + treatment + cultivar:treatment
)

dds.LRT <- DESeq(
  dds.LRT,
  test = "LRT",
  reduced = ~ cultivar + treatment
)

res.list$Interaction <- results(
  dds.LRT,
  alpha = 0.05
)


##############################################################
# 15. DEG FUNCTION
##############################################################

alpha <- 0.05
lfc.cutoff <- 1

get.degs <- function(res) {
  
  x <- as.data.frame(res)
  
  x$padj[is.na(x$padj)] <- 1
  
  genes <- rownames(x)[
    x$padj < alpha &
      abs(x$log2FoldChange) >= lfc.cutoff
  ]
  
  genes <- trimws(genes)
  
  genes <- unique(genes)
  
  genes
}


############################
# 16. DEG LISTS
############################

DEG.Catuai <- get.degs(
  res.list$Catuai
)

DEG.CR95 <- get.degs(
  res.list$CR95
)

DEG.saline <- get.degs(
  res.list$CR95_vs_Catuai_saline
)

DEG.xylella <- get.degs(
  res.list$CR95_vs_Catuai_xylella
)

interaction.genes <- rownames(
  res.list$Interaction
)[
  !is.na(
    res.list$Interaction$padj
  ) &
    res.list$Interaction$padj < alpha
]

interaction.genes <- unique(
  trimws(interaction.genes)
)


############################
# 17. PRINT DEG COUNTS
############################

cat("\n============================\n")
cat("DEG COUNTS\n")
cat("============================\n")

cat(
  "Catuai:",
  length(DEG.Catuai),
  "\n"
)

cat(
  "CR95:",
  length(DEG.CR95),
  "\n"
)

cat(
  "CR95 vs Catuai saline:",
  length(DEG.saline),
  "\n"
)

cat(
  "CR95 vs Catuai xylella:",
  length(DEG.xylella),
  "\n"
)

cat(
  "Interaction:",
  length(interaction.genes),
  "\n"
)


##############################################################
# 18. IMPORTANT VENN CHECK
#
# This tells us whether the overlap is REALLY zero.
##############################################################

shared.genes <- intersect(
  DEG.Catuai,
  DEG.CR95
)

catuai.only <- setdiff(
  DEG.Catuai,
  DEG.CR95
)

cr95.only <- setdiff(
  DEG.CR95,
  DEG.Catuai
)

cat("\n============================\n")
cat("VENN CHECK\n")
cat("============================\n")

cat(
  "Catuai total:",
  length(DEG.Catuai),
  "\n"
)

cat(
  "CR95 total:",
  length(DEG.CR95),
  "\n"
)

cat(
  "Shared genes:",
  length(shared.genes),
  "\n"
)

if (length(shared.genes) > 0) {
  
  cat("\nShared genes detected:\n")
  
  print(
    head(
      shared.genes,
      30
    )
  )
  
} else {
  
  cat(
    "\nWARNING: intersection is actually ZERO.\n"
  )
  
  cat(
    "This means the gene IDs in the two DEG lists do not match exactly,\n",
    "OR no genes satisfy both DEG criteria.\n"
  )
}


############################
# Save overlap genes
############################

write.table(
  shared.genes,
  file = file.path(
    outpath,
    "Shared_DEGs.txt"
  ),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


##############################################################
# 19. VENN DIAGRAM
#
# FIXED-SIZE PERFECT CIRCLES
#
# scaled = FALSE means the circles do NOT change size
# according to the number of genes.
#
# Alpha is calculated separately from the DEG count.
##############################################################

# Alpha based on number of genes
# More genes = more opaque
# Fewer genes = more transparent

max.genes <- max(
  length(DEG.Catuai),
  length(DEG.CR95),
  1
)

alpha.catuai <- 0.25 +
  0.55 *
  (
    length(DEG.Catuai) /
      max.genes
  )

alpha.cr95 <- 0.25 +
  0.55 *
  (
    length(DEG.CR95) /
      max.genes
  )

fill.catuai <- adjustcolor(
  "#E76F51",
  alpha.f = alpha.catuai
)

fill.cr95 <- adjustcolor(
  "#457B9D",
  alpha.f = alpha.cr95
)


venn.plot <- venn.diagram(
  
  x = list(
    Catuai = DEG.Catuai,
    CR95 = DEG.CR95
  ),
  
  filename = NULL,
  
  # IMPORTANT
  scaled = FALSE,
  
  # Same size circles
  fill = c(
    fill.catuai,
    fill.cr95
  ),
  
  # Keep transparency from fill
  alpha = 1,
  
  lwd = 2,
  
  lty = "blank",
  
  cex = 1.5,
  
  fontface = "bold",
  
  cat.cex = 1.4,
  
  cat.fontface = "bold",
  
  cat.pos = c(
    -20,
    20
  ),
  
  cat.dist = c(
    0.05,
    0.05
  ),
  
  main = "Differentially expressed genes",
  
  main.cex = 1.6,
  
  main.fontface = "bold"
)

grid.newpage()
grid.draw(venn.plot)


##############################################################
# 20. OVERLAP SUMMARY
##############################################################

overlap.table <- data.frame(
  
  category = c(
    "Catuai only",
    "Shared",
    "CR95 only"
  ),
  
  genes = c(
    length(catuai.only),
    length(shared.genes),
    length(cr95.only)
  )
)

write.csv(
  overlap.table,
  file.path(
    outpath,
    "DEG_overlap_summary.csv"
  ),
  row.names = FALSE
)
# =====================================================================
# SCRIPT COMPLETO: FILTRADO, ANOTACIÓN Y CLASIFICACIÓN DE GENES (DEG)
# =====================================================================

# ---------------------------------------------------------------------
# 1. PROCESAR CATUAI
# ---------------------------------------------------------------------

# Cargar la tabla completa de anotaciones desde el directorio
DEG_summary_Catuai <- read.csv("DEG_summary_Catuai.csv")

# Filtrar solo los genes que están en tu lista DEG.Catuai
anotaciones_Catuai <- DEG_summary_Catuai[DEG_summary_Catuai$gene_id %in% DEG.Catuai, ]

# Seleccionar únicamente las columnas que interesan
columnas_utiles <- c("gene_id", "RefseqAnnot_protein", "log2FoldChange")
tabla_limpia_Catuai <- anotaciones_Catuai[, columnas_utiles]

# Crear la nueva columna con la etiqueta de regulación
tabla_limpia_Catuai$Regulacion <- ifelse(tabla_limpia_Catuai$log2FoldChange > 0, "upregulated", "down")

# Organizar de más a menos importante (por magnitud del cambio absoluto)
Catuai_final <- tabla_limpia_Catuai[order(abs(tabla_limpia_Catuai$log2FoldChange), decreasing = TRUE), ]

# Separar los que tienen anotación de los que no (NA)
con_anotacion_Catuai <- Catuai_final[!is.na(Catuai_final$RefseqAnnot_protein), ]
sin_anotacion_Catuai <- Catuai_final[is.na(Catuai_final$RefseqAnnot_protein), ]

# Guardar los archivos finales en tu directorio
write.csv(con_anotacion_Catuai, file = "Catuai_Genes_Con_Etiqueta.csv", row.names = FALSE)
write.csv(sin_anotacion_Catuai, file = "Catuai_Sin_Anotacion_Con_Etiqueta.csv", row.names = FALSE)


# ---------------------------------------------------------------------
# 2. PROCESAR CR95
# ---------------------------------------------------------------------

# Cargar la tabla completa de anotaciones desde el directorio
DEG_summary_CR95 <- read.csv("DEG_summary_CR95.csv")

# Filtrar solo los genes que están en tu lista DEG.CR95
anotaciones_CR95 <- DEG_summary_CR95[DEG_summary_CR95$gene_id %in% DEG.CR95, ]

# Seleccionar únicamente las columnas que interesan
tabla_limpia_CR95 <- anotaciones_CR95[, columnas_utiles]

# Crear la nueva columna con la etiqueta de regulación
tabla_limpia_CR95$Regulacion <- ifelse(tabla_limpia_CR95$log2FoldChange > 0, "upregulated", "down")

# Organizar de más a menos importante (por magnitud del cambio absoluto)
CR95_final <- tabla_limpia_CR95[order(abs(tabla_limpia_CR95$log2FoldChange), decreasing = TRUE), ]

# Separar los que tienen anotación de los que no (NA)
con_anotacion_CR95 <- CR95_final[!is.na(CR95_final$RefseqAnnot_protein), ]
sin_anotacion_CR95 <- CR95_final[is.na(CR95_final$RefseqAnnot_protein), ]

# Guardar los archivos finales en tu directorio
write.csv(con_anotacion_CR95, file = "CR95_Genes_Con_Etiqueta.csv", row.names = FALSE)
write.csv(sin_anotacion_CR95, file = "CR95_Sin_Anotacion_Con_Etiqueta.csv", row.names = FALSE)

# =====================================================================
# FILTRAR GENES RELACIONADOS CON MAP KINASES (MAPK)
# =====================================================================

# 1. Buscar en Catuai (busca tanto en mayúsculas como en minúsculas ignorando diferencias con ignore.case)
mapk_Catuai <- con_anotacion_Catuai[grepl("mapk|map kinase|mitogen-activated|mek", con_anotacion_Catuai$RefseqAnnot_protein, ignore.case = TRUE), ]

# 2. Buscar en CR95
mapk_CR95 <- con_anotacion_CR95[grepl("mapk|map kinase|mitogen-activated|mek", con_anotacion_CR95$RefseqAnnot_protein, ignore.case = TRUE), ]

# 3. Ver los resultados directamente en tu consola de R
print("Genes MAPK en Catuai:")
print(mapk_Catuai)

print("Genes MAPK en CR95:")
print(mapk_CR95)

# 4. (Opcional) Guardarlos en nuevos archivos CSV
write.csv(mapk_Catuai, file = "Catuai_Genes_MAPK.csv", row.names = FALSE)
write.csv(mapk_CR95, file = "CR95_Genes_MAPK.csv", row.names = FALSE)

# Búsqueda más amplia de quinasas en Catuai
kinasas_Catuai <- con_anotacion_Catuai[grepl("kinase", con_anotacion_Catuai$RefseqAnnot_protein, ignore.case = TRUE), ]
kinasas_CR95 <- con_anotacion_CR95[grepl("kinase", con_anotacion_CR95$RefseqAnnot_protein, ignore.case = TRUE), ]

# Ver los primeros resultados
print(kinasas_Catuai)

#BiocManager::install("pathview")
# =====================================================================
# SCRIPT DIRECTO: CARGAR ARCHIVOS "CON ANOTACIÓN" Y MAPEAR EN KEGG
# =====================================================================
# =====================================================================
# SCRIPT MAESTRO: PROCESAMIENTO Y KEGG PARA CATUAI Y CR95 (ANOTADOS Y NO ANOTADOS)
# =====================================================================
library(pathview)

# 1. Cargar tus archivos originales
orig_Catuai <- read.csv("DEG_summary_Catuai.csv")
orig_CR95 <- read.csv("DEG_summary_CR95.csv")

# 2. Función blindada que genera el mapa y lo renombra al instante
procesar_y_mapear_seguro <- function(df_original, lista_genes, nombre_variedad, nombre_salida) {
  
  # A. Filtrar y limpiar códigos K
  datos_filtrados <- df_original[df_original$gene_id %in% lista_genes, ]
  datos_limpios <- datos_filtrados[!is.na(datos_filtrados$KEGG_ko) & datos_filtrados$KEGG_ko != "", ]
  codigos_limpios <- regmatches(datos_limpios$KEGG_ko, regexpr("K[0-9]{5}", datos_limpios$KEGG_ko))
  
  filas_validas <- codigos_limpios != ""
  if (sum(filas_validas) == 0) {
    message(paste("No se encontraron códigos K válidos para:", nombre_variedad))
    return(NULL)
  }
  
  vector_kegg <- datos_limpios$log2FoldChange[filas_validas]
  names(vector_kegg) <- codigos_limpios[filas_validas]
  
  # B. Generar el mapa (siempre crea por defecto ko04016.pathview.png)
  message(paste("Generando mapa KEGG para:", nombre_variedad))
  pathview(gene.data = vector_kegg, 
           pathway.id = "04016", 
           species = "ko", 
           gene.idType = "kegg")
  
  # C. Renombrar inmediatamente el archivo predeterminado para que no se sobrescriba
  archivo_original_generado <- "ko04016.pathview.png"
  nuevo_nombre_archivo <- paste0(nombre_salida, ".png")
  
  if (file.exists(archivo_original_generado)) {
    # Si ya existía uno viejo con ese nombre, lo borramos para evitar conflictos
    if (file.exists(nuevo_nombre_archivo)) file.remove(nuevo_nombre_archivo)
    
    file.rename(archivo_original_generado, nuevo_nombre_archivo)
    message(paste("¡Guardado con éxito como:", nuevo_nombre_archivo))
  }
}

# 3. Ejecutar para Catuai y CR95 de forma separada y segura
procesar_y_mapear_seguro(orig_Catuai, DEG.Catuai, "Catuai", "Mapa_MAPK_Catuai")
procesar_y_mapear_seguro(orig_CR95, DEG.CR95, "CR95", "Mapa_MAPK_CR95")

 ##############################################################
# 21. VOLCANO FUNCTION
##############################################################

plot.volcano <- function(
    res,
    title
) {
  
  x <- as.data.frame(res)
  
  x$padj[is.na(x$padj)] <- 1
  
  x$significance <- "Not significant"
  
  x$significance[
    x$padj < alpha &
      x$log2FoldChange >= lfc.cutoff
  ] <- "Upregulated"
  
  x$significance[
    x$padj < alpha &
      x$log2FoldChange <= -lfc.cutoff
  ] <- "Downregulated"
  
  ggplot(
    x,
    aes(
      x = log2FoldChange,
      y = -log10(
        pmax(
          padj,
          1e-300
        )
      ),
      color = significance
    )
  ) +
    
    geom_point(
      alpha = 0.65,
      size = 1.4
    ) +
    
    geom_vline(
      xintercept = c(
        -lfc.cutoff,
        lfc.cutoff
      ),
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed"
    ) +
    
    scale_color_manual(
      values = c(
        "Upregulated" = "#C62828",
        "Downregulated" = "#1565C0",
        "Not significant" = "grey75"
      )
    ) +
    
    labs(
      title = title,
      x = "log2 Fold Change",
      y = expression(
        -log[10](adjusted~p~value)
      ),
      color = NULL
    ) +
    
    theme_minimal(
      base_size = 13
    )
}


##############################################################
# 22. VOLCANO PLOTS -- ALL MAIN COMPARISONS
##############################################################

volcano.Catuai <- plot.volcano(
  res.list$Catuai,
  "Catuai: Xylella vs saline"
)

volcano.CR95 <- plot.volcano(
  res.list$CR95,
  "CR95: Xylella vs saline"
)

volcano.saline <- plot.volcano(
  res.list$CR95_vs_Catuai_saline,
  "CR95 vs Catuai: saline"
)

volcano.xylella <- plot.volcano(
  res.list$CR95_vs_Catuai_xylella,
  "CR95 vs Catuai: Xylella"
)

volcano.treatment <- plot.volcano(
  res.list$Treatment_main,
  "Main treatment effect"
)

volcano.cultivar <- plot.volcano(
  res.list$Cultivar_main,
  "Main cultivar effect"
)


grid.arrange(
  volcano.Catuai,
  volcano.CR95,
  volcano.saline,
  volcano.xylella,
  volcano.treatment,
  volcano.cultivar,
  ncol = 2
)


##############################################################
# 23. SAVE VOLCANO PLOTS
##############################################################

ggsave(
  file.path(
    outpath,
    "Volcano_all_comparisons.png"
  ),
  grid.arrange(
    volcano.Catuai,
    volcano.CR95,
    volcano.saline,
    volcano.xylella,
    volcano.treatment,
    volcano.cultivar,
    ncol = 2
  ),
  width = 14,
  height = 18,
  dpi = 300
)


##############################################################
# 24. MA PLOT FUNCTION
##############################################################

plot.ma <- function(
    res,
    title
) {
  
  x <- as.data.frame(res)
  
  x$padj[is.na(x$padj)] <- 1
  
  x$sig <- (
    x$padj < alpha &
      abs(x$log2FoldChange) >= lfc.cutoff
  )
  
  ggplot(
    x,
    aes(
      x = log10(
        baseMean + 1
      ),
      y = log2FoldChange,
      color = sig
    )
  ) +
    
    geom_point(
      alpha = 0.5,
      size = 1.2
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    
    scale_color_manual(
      values = c(
        `TRUE` = "#CC3333",
        `FALSE` = "grey70"
      )
    ) +
    
    labs(
      title = title,
      x = "log10(mean normalized expression + 1)",
      y = "log2 Fold Change",
      color = "DEG"
    ) +
    
    theme_minimal(
      base_size = 13
    )
}


##############################################################
#### 25. MA PLOTS ####
##############################################################

ma.Catuai <- plot.ma(
  res.list$Catuai,
  "MA: Catuai"
)

ma.CR95 <- plot.ma(
  res.list$CR95,
  "MA: CR95"
)

ma.saline <- plot.ma(
  res.list$CR95_vs_Catuai_saline,
  "MA: CR95 vs Catuai under saline"
)

ma.xylella <- plot.ma(
  res.list$CR95_vs_Catuai_xylella,
  "MA: CR95 vs Catuai under Xylella"
)

ma.treatment <- plot.ma(
  res.list$Treatment_main,
  "MA: Main treatment effect"
)

ma.cultivar <- plot.ma(
  res.list$Cultivar_main,
  "MA: Main cultivar effect"
)

grid.arrange(
  ma.Catuai,
  ma.CR95,
  ncol = 2
)

grid.arrange(
  ma.saline,
  ma.xylella,
  ncol = 2
)

grid.arrange(
  ma.treatment,
  ma.cultivar,
  ncol = 2
)


##############################################################
# 26. LOAD ANNOTATION
##############################################################

annotation <- read.table(
  "fullAnnotation.txt",
  header = TRUE,
  sep = "\t",
  quote = "\"",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c(
    "-",
    "NA",
    ""
  )
)

annotation$gene_id <- trimws(
  annotation$gene_id
)

annotation.gene <- annotation[
  !duplicated(annotation$gene_id),
  ,
  drop = FALSE
]


##############################################################
# 27. CAMERA ANALYSIS #####
##############################################################

has_go <- (
  !is.na(annotation$GOs) &
    annotation$GOs != ""
)

gene2go <- strsplit(
  annotation$GOs[has_go],
  ","
)

names(gene2go) <- annotation$gene_id[has_go]

term2gene <- split(
  rep(
    names(gene2go),
    lengths(gene2go)
  ),
  trimws(
    unlist(gene2go)
  )
)

min_set_size <- 5

camera_results <- list()

voom_objects <- list()


for (cv in c("Catuai", "CR95")) {
  
  cat(
    "\n=== CAMERA:",
    cv,
    "===\n"
  )
  
  keep.samples <- (
    metadata$cultivar == cv
  )
  
  metadata.cv <- metadata[
    keep.samples,
    ,
    drop = FALSE
  ]
  
  counts.cv <- counts[
    ,
    keep.samples,
    drop = FALSE
  ]
  
  treatment <- factor(
    metadata.cv$treatment,
    levels = c(
      "saline",
      "xylella"
    )
  )
  
  dge <- DGEList(
    counts = counts.cv,
    group = treatment
  )
  
  keep.expr <- filterByExpr(
    dge
  )
  
  dge <- dge[
    keep.expr,
    ,
    keep.lib.sizes = FALSE
  ]
  
  dge <- calcNormFactors(
    dge
  )
  
  design.cv <- model.matrix(
    ~ treatment
  )
  
  v <- voom(
    dge,
    design.cv,
    plot = FALSE
  )
  
  voom_objects[[cv]] <- v
  
  term2gene.cv <- lapply(
    term2gene,
    function(g) {
      intersect(
        g,
        rownames(v$E)
      )
    }
  )
  
  term2gene.cv <- term2gene.cv[
    lengths(term2gene.cv) >=
      min_set_size
  ]
  
  idx <- lapply(
    term2gene.cv,
    function(g) {
      match(
        g,
        rownames(v$E)
      )
    }
  )
  
  camera.res <- camera(
    v,
    index = idx,
    design = design.cv,
    contrast = 2
  )
  
  camera.res$GOID <-
    rownames(camera.res)
  
  term.info <- suppressMessages(
    AnnotationDbi::select(
      GO.db,
      keys = camera.res$GOID,
      columns = c(
        "TERM",
        "ONTOLOGY"
      ),
      keytype = "GOID"
    )
  )
  
  camera.res <- merge(
    camera.res,
    term.info,
    by = "GOID",
    all.x = TRUE
  )
  
  camera.res$TERM[
    is.na(camera.res$TERM)
  ] <-
    camera.res$GOID[
      is.na(camera.res$TERM)
    ]
  
  camera.res$ONTOLOGY[
    is.na(camera.res$ONTOLOGY)
  ] <- "Unknown"
  
  camera.res <- camera.res[
    order(camera.res$PValue),
    ,
    drop = FALSE
  ]
  
  camera_results[[cv]] <-
    camera.res
  
  write.table(
    camera.res,
    file = file.path(
      outpath,
      paste0(
        "CAMERA_",
        cv,
        ".txt"
      )
    ),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )
}


##############################################################
# 28. CAMERA SIGNIFICANT GO TERMS
##############################################################

sig_cutoff <- 0.05

sig_results <- lapply(
  camera_results,
  function(x) {
    x[
      x$FDR < sig_cutoff,
      ,
      drop = FALSE
    ]
  }
)

cat("\nSignificant CAMERA GO terms:\n")

print(
  sapply(
    sig_results,
    nrow
  )
)


##############################################################
# 29. CAMERA BARPLOTS
##############################################################

plot.camera.barplot <- function(
    res,
    title,
    top_n = 50
) {
  
  if (nrow(res) == 0) {
    return(
      ggplot() +
        labs(
          title = paste(
            title,
            "- no significant GO terms"
          )
        ) +
        theme_void()
    )
  }
  
  res.top <- head(
    res[
      order(res$FDR),
      ,
      drop = FALSE
    ],
    top_n
  )
  
  res.top$TERM <- factor(
    res.top$TERM,
    levels = rev(
      res.top$TERM
    )
  )
  
  ggplot(
    res.top,
    aes(
      x = TERM,
      y = -log10(FDR),
      fill = Direction
    )
  ) +
    
    geom_col() +
    
    coord_flip() +
    
    scale_fill_manual(
      values = c(
        Up = "#66A61E",
        Down = "#FFD92F"
      )
    ) +
    
    labs(
      title = title,
      x = "GO term",
      y = "-log10(FDR)"
    ) +
    
    theme_bw()
}


camera.catuai.plot <-
  plot.camera.barplot(
    sig_results$Catuai,
    "Catuai"
  )

camera.cr95.plot <-
  plot.camera.barplot(
    sig_results$CR95,
    "CR95"
  )

grid.arrange(
  camera.catuai.plot,
  camera.cr95.plot,
  ncol = 2
)



##############################################################
# 30. GO VENN####
##############################################################

sig.catuai.go <- unique(
  camera_results$Catuai$GOID[
    camera_results$Catuai$FDR <
      sig_cutoff
  ]
)

sig.cr95.go <- unique(
  camera_results$CR95$GOID[
    camera_results$CR95$FDR <
      sig_cutoff
  ]
)

cat(
  "\nShared significant GO terms:",
  length(
    intersect(
      sig.catuai.go,
      sig.cr95.go
    )
  ),
  "\n"
)

go.venn <- venn.diagram(
  
  x = list(
    Catuai = sig.catuai.go,
    CR95 = sig.cr95.go
  ),
  
  filename = NULL,
  
  scaled = FALSE,
  
  fill = c(
    "#E5C494",
    "#A6761D"
  ),
  
  alpha = 0.6,
  
  lwd = 2,
  
  lty = "blank",
  
  cex = 1.5,
  
  cat.cex = 1.3,
  
  main =
    "Significant GO terms: Catuai vs CR95"
)

grid.newpage()
grid.draw(go.venn)


##############################################################
# 31. HYPERGEOMETRIC GO ENRICHMENT
#
# Complete analysis.
#
# Target:
# significant DEGs
#
# Background:
# all genes that survived DESeq2 filtering
#
# This is important because enrichment should NOT use
# the entire annotation file as background.
##############################################################

run.hypergeometric <- function(
    target.genes,
    background.genes,
    annotation,
    name
) {
  
  target.genes <- unique(
    trimws(target.genes)
  )
  
  background.genes <- unique(
    trimws(background.genes)
  )
  
  # Only genes present in background
  target.genes <- intersect(
    target.genes,
    background.genes
  )
  
  annotation.use <- annotation[
    annotation$gene_id %in%
      background.genes,
    ,
    drop = FALSE
  ]
  
  annotation.use <- annotation.use[
    !duplicated(
      annotation.use$gene_id
    ),
    ,
    drop = FALSE
  ]
  
  valid <- (
    !is.na(annotation.use$GOs) &
      annotation.use$GOs != ""
  )
  
  annotation.use <- annotation.use[
    valid,
    ,
    drop = FALSE
  ]
  
  if (
    length(target.genes) == 0
  ) {
    
    warning(
      paste(
        "No target genes for",
        name
      )
    )
    
    return(NULL)
  }
  
  go.list <- strsplit(
    annotation.use$GOs,
    ","
  )
  
  names(go.list) <-
    annotation.use$gene_id
  
  gene2go <- data.frame(
    
    gene = rep(
      names(go.list),
      lengths(go.list)
    ),
    
    GO = trimws(
      unlist(go.list)
    ),
    
    stringsAsFactors = FALSE
  )
  
  gene2go <- gene2go[
    !is.na(gene2go$GO) &
      gene2go$GO != "" &
      gene2go$GO != "-",
    ,
    drop = FALSE
  ]
  
  # Remove duplicated gene-GO combinations
  gene2go <- unique(
    gene2go
  )
  
  all.terms <- unique(
    gene2go$GO
  )
  
  N <- length(
    unique(
      background.genes
    )
  )
  
  n <- length(
    unique(
      target.genes
    )
  )
  
  results <- lapply(
    all.terms,
    function(term) {
      
      term.genes <- unique(
        gene2go$gene[
          gene2go$GO == term
        ]
      )
      
      K <- length(
        intersect(
          term.genes,
          background.genes
        )
      )
      
      k <- length(
        intersect(
          term.genes,
          target.genes
        )
      )
      
      # Hypergeometric probability:
      #
      # probability of obtaining k or more genes
      # from this GO category by chance
      #
      p <- phyper(
        k - 1,
        K,
        N - K,
        n,
        lower.tail = FALSE
      )
      
      data.frame(
        GOID = term,
        GeneHits = k,
        BackgroundHits = K,
        TargetGenes = n,
        BackgroundGenes = N,
        pvalue = p,
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- do.call(
    rbind,
    results
  )
  
  results$FDR <- p.adjust(
    results$pvalue,
    method = "BH"
  )
  
  # Add GO names
  term.info <- suppressMessages(
    AnnotationDbi::select(
      GO.db,
      keys = results$GOID,
      columns = c(
        "TERM",
        "ONTOLOGY"
      ),
      keytype = "GOID"
    )
  )
  
  results <- merge(
    results,
    term.info,
    by.x = "GOID",
    by.y = "GOID",
    all.x = TRUE
  )
  
  results$TERM[
    is.na(results$TERM)
  ] <-
    results$GOID[
      is.na(results$TERM)
    ]
  
  results$ONTOLOGY[
    is.na(results$ONTOLOGY)
  ] <- "Unknown"
  
  results <- results[
    order(results$FDR),
    ,
    drop = FALSE
  ]
  
  write.csv(
    results,
    file.path(
      outpath,
      paste0(
        name,
        "_Hypergeometric_GO.csv"
      )
    ),
    row.names = FALSE
  )
  
  return(results)
}


##############################################################
# 32. RUN HYPERGEOMETRIC ANALYSIS
##############################################################

background.genes <- rownames(dds)

hyper.Catuai <- run.hypergeometric(
  DEG.Catuai,
  background.genes,
  annotation,
  "Catuai"
)

hyper.CR95 <- run.hypergeometric(
  DEG.CR95,
  background.genes,
  annotation,
  "CR95"
)

hyper.interaction <- run.hypergeometric(
  interaction.genes,
  background.genes,
  annotation,
  "Interaction"
)




##############################################################
# 33. HYPERGEOMETRIC DOTPLOT
##############################################################

plot.hypergeometric <- function(
    res,
    title,
    top_n = 20
) {
  
  if (
    is.null(res) ||
    nrow(res) == 0
  ) {
    
    return(
      ggplot() +
        labs(
          title =
            paste(
              title,
              "- no results"
            )
        ) +
        theme_void()
    )
  }
  
  sig <- res[
    res$FDR < 0.05,
    ,
    drop = FALSE
  ]
  
  if (nrow(sig) == 0) {
    
    return(
      ggplot() +
        labs(
          title =
            paste(
              title,
              "- no significant GO terms"
            )
        ) +
        theme_void()
    )
  }
  
  sig <- head(
    sig[
      order(sig$FDR),
      ,
      drop = FALSE
    ],
    top_n
  )
  
  sig$TERM <- factor(
    sig$TERM,
    levels = rev(
      sig$TERM
    )
  )
  
  ggplot(
    sig,
    aes(
      x = -log10(FDR),
      y = TERM,
      size = GeneHits
    )
  ) +
    
    geom_point() +
    
    labs(
      title = title,
      x = "-log10(FDR)",
      y = "GO term",
      size = "Genes"
    ) +
    
    theme_minimal(
      base_size = 12
    )
}


hyper.plot.catuai <-
  plot.hypergeometric(
    hyper.Catuai,
    "Hypergeometric GO: Catuai"
  )

hyper.plot.cr95 <-
  plot.hypergeometric(
    hyper.CR95,
    "Hypergeometric GO: CR95"
  )

hyper.plot.interaction <-
  plot.hypergeometric(
    hyper.interaction,
    "Hypergeometric GO: Interaction"
  )

grid.arrange(
  hyper.plot.catuai,
  hyper.plot.cr95,
  hyper.plot.interaction,
  ncol = 1
)


##############################################################
# 34. TOP VARIABLE GENES HEATMAP
#
# IMPORTANT:
# Columns are NOT clustered.
#
# Biological order:
#
# Catuai
# saline | xylella
#
# CR95
# saline | xylella
##############################################################

gene.var <- apply(
  vst.mat,
  1,
  var
)

top.variable.genes <- names(
  sort(
    gene.var,
    decreasing = TRUE
  )
)[
  1:min(
    500,
    length(gene.var)
  )
]

pheatmap(
  t(scale(
    t(
      vst.mat[
        top.variable.genes,
        ,
        drop = FALSE
      ]
    )
  )),
  
  cluster_rows = TRUE,
  
  cluster_cols = FALSE,
  
  annotation_col =
    annotation.samples,
  
  show_rownames = FALSE,
  
  main =
    "Top 500 most variable genes",
  
  gaps_col = c(
    sum(
      metadata$cultivar ==
        "Catuai"
    )
  )
)


##############################################################
# 35. HEATMAP OF Catuai + CR95 DEGs
##############################################################

all.degs <- unique(
  c(
    DEG.Catuai,
    DEG.CR95
  )
)

heatmap.degs <- intersect(
  all.degs,
  rownames(vst.mat)
)

if (
  length(heatmap.degs) > 1
) {
  
  heatmap.mat <- t(
    scale(
      t(
        vst.mat[
          heatmap.degs,
          ,
          drop = FALSE
        ]
      )
    )
  )
  
  pheatmap(
    heatmap.mat,
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    annotation_col =
      annotation.samples,
    
    show_rownames =
      length(heatmap.degs) <= 50,
    
    fontsize_row = 7,
    
    gaps_col = sum(
      metadata$cultivar ==
        "Catuai"
    ),
    
    main =
      paste0(
        "DEGs: Catuai + CR95 (",
        length(heatmap.degs),
        " genes)"
      )
  )
}


##############################################################
# 36. INTERACTION HEATMAP
##############################################################

heatmap.interaction <- intersect(
  interaction.genes,
  rownames(vst.mat)
)

if (
  length(heatmap.interaction) > 1
) {
  
  mat.interaction <- t(
    scale(
      t(
        vst.mat[
          heatmap.interaction,
          ,
          drop = FALSE
        ]
      )
    )
  )
  
  pheatmap(
    mat.interaction,
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    annotation_col =
      annotation.samples,
    
    show_rownames =
      length(heatmap.interaction) <= 50,
    
    fontsize_row = 7,
    
    gaps_col = sum(
      metadata$cultivar ==
        "Catuai"
    ),
    
    main =
      paste0(
        "Significant interaction genes (",
        length(
          heatmap.interaction
        ),
        ")"
      )
  )
}


##############################################################
# 37. TOP 50 INTERACTION GENES
##############################################################

interaction.order <- order(
  res.list$Interaction$padj,
  na.last = NA
)

top.interaction <- rownames(
  res.list$Interaction
)[
  head(
    interaction.order,
    50
  )
]

top.interaction <- intersect(
  top.interaction,
  rownames(vst.mat)
)

if (
  length(top.interaction) > 1
) {
  
  pheatmap(
    t(
      scale(
        t(
          vst.mat[
            top.interaction,
            ,
            drop = FALSE
          ]
        )
      )
    ),
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    annotation_col =
      annotation.samples,
    
    show_rownames = TRUE,
    
    fontsize_row = 7,
    
    gaps_col = sum(
      metadata$cultivar ==
        "Catuai"
    ),
    
    main =
      "Top 50 interaction genes"
  )
}


##############################################################
# 38. EXPRESSION PROFILES OF INTERACTION GENES
##############################################################

group.means <- function(
    gene
) {
  
  data.frame(
    
    sample = colnames(vst.mat),
    
    expression =
      vst.mat[
        gene,
      ],
    
    cultivar =
      metadata$cultivar,
    
    treatment =
      metadata$treatment,
    
    group =
      metadata$group
  )
}

top.profile.genes <-
  head(
    top.interaction,
    12
  )

if (
  length(top.profile.genes) > 0
) {
  
  profile.df <- do.call(
    rbind,
    lapply(
      top.profile.genes,
      function(g) {
        
        x <- group.means(g)
        
        x$gene <- g
        
        x
      }
    )
  )
  
  profile.df$group <- factor(
    profile.df$group,
    levels = desired.order
  )
  
  profile.plot <- ggplot(
    profile.df,
    aes(
      x = treatment,
      y = expression,
      group = cultivar,
      color = cultivar
    )
  ) +
    
    stat_summary(
      fun = mean,
      geom = "line",
      linewidth = 1
    ) +
    
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2
    ) +
    
    facet_wrap(
      ~gene,
      scales = "free_y"
    ) +
    
    labs(
      title =
        "Expression profiles of top interaction genes",
      y = "VST expression"
    ) +
    
    theme_minimal(
      base_size = 12
    )
  
  print(profile.plot)
}


##############################################################
# 39. ANNOTATED DESEQ2 TABLES
##############################################################

make.annotated.summary <- function(
    res,
    filename
) {
  
  x <- as.data.frame(res)
  
  x$gene_id <-
    rownames(x)
  
  x$gene_id <-
    trimws(x$gene_id)
  
  x <- merge(
    x,
    annotation.gene,
    by = "gene_id",
    all.x = TRUE
  )
  
  x <- x[
    order(x$padj),
    ,
    drop = FALSE
  ]
  
  write.csv(
    x,
    file.path(
      outpath,
      filename
    ),
    row.names = FALSE
  )
  
  x
}


summary.Catuai <-
  make.annotated.summary(
    res.list$Catuai,
    "DEG_summary_Catuai.csv"
  )

summary.CR95 <-
  make.annotated.summary(
    res.list$CR95,
    "DEG_summary_CR95.csv"
  )


summary.saline <-
  make.annotated.summary(
    res.list$CR95_vs_Catuai_saline,
    "DEG_summary_Cultivar_saline.csv"
  )


summary.xylella <-
  make.annotated.summary(
    res.list$CR95_vs_Catuai_xylella,
    "DEG_summary_Cultivar_xylella.csv"
  )



summary.interaction <-
  make.annotated.summary(
    res.list$Interaction,
    "Interaction_DEG_summary.csv"
  )


##############################################################
# 40. WGCNA #####THIS ONE IS SELECTING THE INCORRECT MODULE BTW
##############################################################

expr.full <- t(vst.mat)

gene.var.wgcna <- apply(
  expr.full,
  2,
  var
)

top.n.genes <- min(
  5000,
  ncol(expr.full)
)

top.genes <- names(
  sort(
    gene.var.wgcna,
    decreasing = TRUE
  )
)[
  1:top.n.genes
]

datExpr <- expr.full[
  ,
  top.genes,
  drop = FALSE
]

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 0
)

if (!gsg$allOK) {
  
  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes,
    drop = FALSE
  ]
}


##############################################################
# 41. WGCNA SAMPLE CLUSTERING
##############################################################

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

plot(
  sampleTree,
  main = "WGCNA sample clustering",
  xlab = "",
  sub = ""
)


##############################################################
# 42. WGCNA SOFT THRESHOLD
##############################################################

powers <- c(
  1:20
)

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 0
)

plot(
  sft$fitIndices[, 1],
  -sign(
    sft$fitIndices[, 3]
  ) *
    sft$fitIndices[, 2],
  
  xlab =
    "Soft threshold power",
  
  ylab =
    "Signed R²",
  
  main =
    "Scale-free topology fit"
)

abline(
  h = 0.8,
  lty = 2
)

soft.power <-
  sft$powerEstimate

if (
  is.na(soft.power)
) {
  soft.power <- 6
}


##############################################################
# 43. BUILD WGCNA NETWORK
##############################################################

net <- blockwiseModules(
  
  datExpr,
  
  power = soft.power,
  
  TOMType = "signed",
  
  networkType = "signed",
  
  minModuleSize = 30,
  
  reassignThreshold = 0,
  
  mergeCutHeight = 0.25,
  
  numericLabels = TRUE,
  
  pamRespectsDendro = FALSE,
  
  verbose = 0
)

moduleColors <-
  labels2colors(
    net$colors
  )

names(moduleColors) <-
  colnames(datExpr)

print(
  table(moduleColors)
)


##############################################################
# 44. WGCNA TRAITS
##############################################################

cultivar.bin <- as.numeric(
  metadata[
    rownames(datExpr),
    "cultivar"
  ] == "CR95"
)

treatment.bin <- as.numeric(
  metadata[
    rownames(datExpr),
    "treatment"
  ] == "xylella"
)

interaction.bin <-
  cultivar.bin *
  treatment.bin

group.factor <- interaction(
  metadata[
    rownames(datExpr),
    "cultivar"
  ],
  metadata[
    rownames(datExpr),
    "treatment"
  ],
  sep = "_"
)

group.matrix <- model.matrix(
  ~0 + group.factor
)

colnames(group.matrix) <-
  sub(
    "^group.factor",
    "",
    colnames(group.matrix)
  )

traits <- data.frame(
  
  cultivar_CR95 =
    cultivar.bin,
  
  treatment_xylella =
    treatment.bin,
  
  interaction_CR95_xylella =
    interaction.bin,
  
  group.matrix
)

rownames(traits) <-
  rownames(datExpr)


##############################################################
# 45. MODULE-TRAIT RELATIONSHIPS
##############################################################

MEs <- orderMEs(
  net$MEs
)

moduleTraitCor <- cor(
  MEs,
  traits,
  use = "p"
)

moduleTraitP <- corPvalueStudent(
  moduleTraitCor,
  nSamples =
    nrow(datExpr)
)

textMatrix <- paste0(
  signif(
    moduleTraitCor,
    2
  ),
  "\n(",
  signif(
    moduleTraitP,
    2
  ),
  ")"
)

dim(textMatrix) <-
  dim(moduleTraitCor)

labeledHeatmap(
  
  Matrix =
    moduleTraitCor,
  
  xLabels =
    colnames(traits),
  
  yLabels =
    colnames(MEs),
  
  ySymbols =
    colnames(MEs),
  
  colorLabels = FALSE,
  
  colors =
    blueWhiteRed(50),
  
  textMatrix =
    textMatrix,
  
  setStdMargins = FALSE,
  
  cex.text = 0.7,
  
  zlim =
    c(-1, 1),
  
  main =
    "Module-trait relationships"
)

write.csv(
  moduleTraitCor,
  file.path(
    outpath,
    "WGCNA_module_trait_correlations.csv"
  )
)

write.csv(
  moduleTraitP,
  file.path(
    outpath,
    "WGCNA_module_trait_pvalues.csv"
  )
)


##############################################################
# 46. SELECT BEST INTERACTION MODULE
##############################################################

trait.of.interest <-
  "interaction_CR95_xylella"

module.ranking <- data.frame(
  
  ME =
    rownames(moduleTraitCor),
  
  correlation =
    moduleTraitCor[
      ,
      trait.of.interest
    ],
  
  pvalue =
    moduleTraitP[
      ,
      trait.of.interest
    ]
)

module.ranking$absCorrelation <-
  abs(
    module.ranking$correlation
  )

module.ranking$module <-
  sub(
    "^ME",
    "",
    module.ranking$ME
  )

module.ranking <-
  module.ranking[
    module.ranking$module != "grey",
    ,
    drop = FALSE
  ]

module.ranking <-
  module.ranking[
    order(
      -module.ranking$absCorrelation
    ),
    ,
    drop = FALSE
  ]

print(
  module.ranking
)

best.module.number <- module.ranking$module[1]

# Convert module number to its WGCNA color name
best.module.name <- labels2colors(
  as.numeric(best.module.number)
)

cat(
  "\nSelected WGCNA module:",
  best.module.name,
  "\nGenes:",
  sum(moduleColors == best.module.name),
  "\n"
)

module.genes <- names(
  moduleColors[
    moduleColors == best.module.name
  ]
)

if (length(module.genes) == 0) {
  
  stop(
    paste(
      "Selected module",
      best.module.name,
      "contains 0 genes."
    )
  )
}


##############################################################
# 47. MODULE MEMBERSHIP + GENE SIGNIFICANCE
##############################################################

MM <- cor(
  datExpr,
  MEs,
  use = "p"
)

GS.interaction <- cor(
  datExpr,
  traits[
    ,
    trait.of.interest
  ],
  use = "p"
)

names(GS.interaction) <- colnames(datExpr)

MM.best <- MM[
  module.genes,
  paste0(
    "ME",
    best.module.number
  )
]

GS.best <- GS.interaction[
  module.genes
]
hub.table <- data.frame(
  
  gene_id =
    module.genes,
  
  module =
    rep(
      best.module.name,
      length(module.genes)
    ),
  
  moduleMembership =
    as.numeric(MM.best),
  
  geneSignificance =
    as.numeric(GS.best),
  
  absMM =
    abs(
      as.numeric(MM.best)
    ),
  
  absGS =
    abs(
      as.numeric(GS.best)
    )
)

# Add whether gene is an interaction DEG
hub.table$is_interaction_DEG <-
  hub.table$gene_id %in%
  interaction.genes

write.csv(
  hub.table,
  file.path(
    outpath,
    "WGCNA_hub_genes.csv"
  ),
  row.names = FALSE
)


##############################################################
# 48. GS vs MM
##############################################################

ggplot(
  hub.table,
  aes(
    x = absMM,
    y = absGS,
    color = is_interaction_DEG
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  labs(
    title =
      paste0(
        "Hub genes: ",
        best.module.name
      ),
    
    x =
      "Absolute module membership",
    
    y =
      "Absolute gene significance",
    
    color =
      "Interaction DEG"
  ) +
  
  theme_minimal(
    base_size = 13
  )


##############################################################
# 49. CYTOSCAPE EXPORT
##############################################################

top.hubs <- head(
  hub.table$gene_id,
  150
)

module.interaction.degs <-
  intersect(
    module.genes,
    interaction.genes
  )

network.genes <- unique(
  c(
    top.hubs,
    module.interaction.degs
  )
)

max.network.genes <- 200

if (
  length(network.genes) >
  max.network.genes
) {
  
  network.genes <-
    head(
      hub.table$gene_id[
        hub.table$gene_id %in%
          network.genes
      ],
      max.network.genes
    )
}

cat(
  "\nGenes exported to Cytoscape:",
  length(network.genes),
  "\n"
)


##############################################################
# 50. TOM
##############################################################

TOM.network <-
  TOMsimilarityFromExpr(
    
    datExpr[
      ,
      network.genes,
      drop = FALSE
    ],
    
    power =
      soft.power,
    
    networkType =
      "signed"
  )

dimnames(TOM.network) <-
  list(
    network.genes,
    network.genes
  )

tom.values <- TOM.network[
  upper.tri(
    TOM.network
  )
]

tom.threshold <-
  as.numeric(
    quantile(
      tom.values,
      probs = 0.95,
      na.rm = TRUE
    )
  )


##############################################################
# 51. CYTOSCAPE FILES
##############################################################

exportNetworkToCytoscape(
  
  TOM.network,
  
  edgeFile =
    file.path(
      outpath,
      paste0(
        "Cytoscape_edges_",
        best.module.name,
        ".txt"
      )
    ),
  
  nodeFile =
    file.path(
      outpath,
      paste0(
        "Cytoscape_nodes_",
        best.module.name,
        ".txt"
      )
    ),
  
  weighted = TRUE,
  
  threshold =
    tom.threshold,
  
  nodeNames =
    network.genes,
  
  nodeAttr =
    rep(
      best.module.name,
      length(
        network.genes
      )
    )
)


##############################################################
# 52. CYTOSCAPE NODE ATTRIBUTES
##############################################################

cyt.node.attrs <-
  hub.table[
    hub.table$gene_id %in%
      network.genes,
    ,
    drop = FALSE
  ]

cyt.node.attrs <-
  cyt.node.attrs[
    match(
      network.genes,
      cyt.node.attrs$gene_id
    ),
    ,
    drop = FALSE
  ]

cyt.node.attrs$nodeName <-
  cyt.node.attrs$gene_id

adjacency <-
  TOM.network >
  tom.threshold

diag(adjacency) <-
  FALSE

cyt.node.attrs$network_degree <-
  rowSums(adjacency)

write.table(
  
  cyt.node.attrs,
  
  file =
    file.path(
      outpath,
      paste0(
        "Cytoscape_node_attributes_",
        best.module.name,
        ".txt"
      )
    ),
  
  sep = "\t",
  
  row.names = FALSE,
  
  quote = FALSE
)


##############################################################
# 53. FINAL SUMMARY
##############################################################

cat("\n")
cat("============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n")

cat(
  "\nCatuai DEGs:",
  length(DEG.Catuai)
)

cat(
  "\nCR95 DEGs:",
  length(DEG.CR95)
)

cat(
  "\nShared DEGs:",
  length(shared.genes)
)

cat(
  "\nInteraction genes:",
  length(interaction.genes)
)

cat(
  "\nCatuai significant GO terms:",
  length(sig.catuai.go)
)

cat(
  "\nCR95 significant GO terms:",
  length(sig.cr95.go)
)

cat(
  "\nShared GO terms:",
  length(
    intersect(
      sig.catuai.go,
      sig.cr95.go
    )
  )
)

cat(
  "\nWGCNA selected module:",
  best.module.name
)

cat(
  "\n\nResults saved in:",
  outpath
)

cat("\n============================================\n")




#cambiar con el LOC####



counts <- round(counts)
expression <- unlist(counts["LOC113701250", ])
barplot(expression,
        names.arg = colnames(counts),
        las = 2,
        ylab = "Counts",
        main = "LOC113701250")


summary.Catuai[summary.Catuai$gene_id=="LOC113701250",]
x=summary.Catuai[summary.Catuai$gene_id=="LOC113701250",]
x$log2FoldChange
###SUMMARY CTUAI#####
summary.CR95[summary.CR95$gene_id=="LOC113699448",]
y=summary.CR95[summary.CR95$gene_id=="LOC113699448",]
y$log2FoldChange


summary.saline[summary.saline$gene_id=="LOC113701250",]
z=summary.saline[summary.saline$gene_id=="LOC113701250",]
z$log2FoldChange

summary.xylella[summary.xylella$gene_id=="LOC113701250",]
k=summary.xylella[summary.xylella$gene_id=="LOC113701250",]
k$log2FoldChange


summary.interaction[summary.interaction$gene_id=="LOC113701250",]
q=summary.interaction[summary.interaction$gene_id=="LOC113701250",]
q$log2FoldChange



#install.packages("patchwork")

library(ggplot2)
library(patchwork)

# =========================
# 1. EXPRESIÓN DEL GEN
# =========================

expression <- as.numeric(unlist(counts["LOC113701250", ]))

df <- data.frame(
  sample = colnames(counts),
  expression = expression
)

metadata <- data.frame(
  sample = c(
    "Ca1x", "Ca2x", "Ca3x", "Ca4x", "Ca5x",
    "Ca6x", "Ca7x", "Ca10x",
    "CaC1", "CaC2", "CaC3",
    "CaC6", "CaC7", "CaC9", "CaC10"
  ),
  cultivar = c(
    "Catuai", "Catuai", "Catuai", "Catuai", "Catuai",
    "CR95", "CR95", "CR95",
    "Catuai", "Catuai", "Catuai",
    "CR95", "CR95", "CR95", "CR95"
  ),
  treatment = c(
    "Xylella", "Xylella", "Xylella", "Xylella", "Xylella",
    "Xylella", "Xylella", "Xylella",
    "Saline", "Saline", "Saline",
    "Saline", "Saline", "Saline", "Saline"
  )
)

df <- merge(df, metadata, by = "sample", sort = FALSE)


df$cultivar <- factor(df$cultivar,
                      levels = c("Catuai", "CR95"))

df$treatment <- factor(df$treatment,
                       levels = c("Saline", "Xylella"))


# =========================
# 2. GRÁFICO DE EXPRESIÓN
# =========================

p1 <- ggplot(df, aes(x = sample, y = expression, fill = treatment)) +
  
  geom_col(width = 0.7) +
  
  facet_grid(. ~ cultivar,
             scales = "free_x",
             space = "free_x") +
  
  labs(
    title = "LOC113701250",
    subtitle = "B-box zinc finger protein 18-like",
    x = NULL,
    y = "Raw read counts",
    fill = "Treatment"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


# =========================
# 3. LOG2 FOLD CHANGE
# =========================

logFC <- data.frame(
  comparison = c(
    "Xylella",
    "Saline",
    "CR95",
    "Catuai",
    "Interaction"
  ),
  
log2FC = c(
    k$log2FoldChange,
    z$log2FoldChange,
    y$log2FoldChange,
    x$log2FoldChange,
    q$log2FoldChange
  )
)


# =========================
# 4. GRÁFICO DE log2FC
# =========================

p2 <- ggplot(logFC, aes(x = comparison, y = log2FC)) +
  
  geom_col(width = 0.65) +
  
  geom_hline(yintercept = 0, linewidth = 0.5) +
  
  geom_text(
    aes(label = round(log2FC, 2)),
    vjust = -0.3,
    size = 4
  ) +
  
  labs(
    title = "Differential expression",
    x = NULL,
    y = "log2 Fold Change"
  ) +
  
  scale_y_continuous(
  expand = expansion(mult = c(0.05, 0.15))
  ) +

  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )


# =========================
# 5. UNIR LOS DOS GRÁFICOS
# =========================

p1 / p2 +
  plot_layout(heights = c(2.5, 2))
###que factores de transcripcion se encienden y por que?####

