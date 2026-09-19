
#WGCNA

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



###########################
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

rownames(annotation.samples) <- rownames



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




##############################################################
# REINICIAR DISPOSITIVOS GRÁFICOS PREVIOS
##############################################################

# Cierra cualquier dispositivo gráfico colgado
while (dev.cur() > 1) dev.off()
graphics.off()


##############################################################
# 40. WGCNA PREPARACIÓN DE DATOS
##############################################################

expr.full <- t(vst.mat)

gene.var.wgcna <- apply(expr.full, 2, var)
top.n.genes <- min(5000, ncol(expr.full))
top.genes <- names(sort(gene.var.wgcna, decreasing = TRUE))[1:top.n.genes]

datExpr <- expr.full[, top.genes, drop = FALSE]

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}


##############################################################
# 41. WGCNA SAMPLE CLUSTERING
##############################################################

sampleTree <- hclust(dist(datExpr), method = "average")

pdf(
  file = file.path(outpath, "WGCNA_sample_clustering.pdf"),
  width = 12,
  height = 8
)
plot(sampleTree, main = "WGCNA sample clustering", xlab = "", sub = "")
dev.off()


##############################################################
# 42. WGCNA SOFT THRESHOLD
##############################################################

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

pdf(
  file = file.path(outpath, "WGCNA_soft_threshold.pdf"),
  width = 10,
  height = 7
)
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft threshold power",
  ylab = "Signed R²",
  main = "Scale-free topology fit"
)
abline(h = 0.8, lty = 2)
dev.off()

soft.power <- sft$powerEstimate
if (is.na(soft.power)) {
  soft.power <- 6
}

cat("\nSelected soft threshold power:", soft.power, "\n")


##############################################################
# 43. BUILD WGCNA NETWORK (numericLabels = FALSE)
##############################################################

net <- blockwiseModules(
  datExpr,
  power = soft.power,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 0
)

moduleColors <- net$colors
names(moduleColors) <- colnames(datExpr)

print(table(moduleColors))


##############################################################
# 43B. WGCNA GENE DENDROGRAM + MODULE COLORS
##############################################################

pdf(
  file = file.path(outpath, "WGCNA_gene_dendrogram_modules.pdf"),
  width = 14,
  height = 8
)

# Dibujar dendrograma del primer bloque (redestilado sin loops propensos a errores)
plotDendroAndColors(
  dendro = net$dendrograms[[1]],
  colors = moduleColors[net$blockGenes[[1]]],
  groupLabels = "Módulo",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "WGCNA Gene Dendrogram - Bloque 1"
)

dev.off()
##############################################################
# 44. WGCNA TRAITS (MATRIZ DE DISEÑO COMPLETA DE 4 GRUPOS)
##############################################################

metadata.wgcna <- metadata[rownames(datExpr), , drop = FALSE]

# 1. Definir matriz binaria para cada combinación Cultivar_Tratamiento
group.factor <- interaction(metadata.wgcna$cultivar, metadata.wgcna$treatment, sep = "_")
traits <- model.matrix(~ 0 + group.factor)
colnames(traits) <- sub("^group.factor", "", colnames(traits))
rownames(traits) <- rownames(datExpr)

# 2. Agregar término continuo de interacción matemática para análisis correlacional
traits <- as.data.frame(traits)
traits$Interaction_Effect <- ifelse(
  metadata.wgcna$cultivar == "CR95" & metadata.wgcna$treatment == "xylella", 1,
  ifelse(metadata.wgcna$cultivar == "CR95" | metadata.wgcna$treatment == "xylella", -0.5, 0)
)


##############################################################
# 45. MODULE-TRAIT RELATIONSHIPS (HEATMAP CORREGIDO)
##############################################################

MEs <- orderMEs(net$MEs)

moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

textMatrix <- paste0(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitP, 1),
  ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

pdf(
  file = file.path(outpath, "WGCNA_module_trait_relationships.pdf"),
  width = 11,
  height = 9
)

par(mfrow = c(1, 1))
par(mar = c(8, 9, 3, 3))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits),
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Relación Módulo-Trait (Interacción Completa Catuai vs CR95)"
)

dev.off()


##############################################################
# 46. MODELO DE INTERACCIÓN GLOBAL (ANOVA DE 2 VÍAS POR MÓDULO)
##############################################################

interaction.results <- data.frame(
  ME = colnames(MEs),
  interaction_F = NA,
  interaction_pvalue = NA
)

for (i in seq_len(ncol(MEs))) {
  model <- lm(MEs[, i] ~ cultivar * treatment, data = metadata.wgcna)
  model.anova <- anova(model)
  
  interaction.results$interaction_F[i]      <- model.anova["cultivar:treatment", "F value"]
  interaction.results$interaction_pvalue[i] <- model.anova["cultivar:treatment", "Pr(>F)"]
}

interaction.results$interaction_FDR <- p.adjust(
  interaction.results$interaction_pvalue,
  method = "BH"
)

interaction.results$module <- sub("^ME", "", interaction.results$ME)
interaction.results <- interaction.results[interaction.results$module != "grey", , drop = FALSE]
interaction.results <- interaction.results[order(interaction.results$interaction_pvalue), , drop = FALSE]

write.csv(
  interaction.results,
  file.path(outpath, "WGCNA_module_interaction_results.csv"),
  row.names = FALSE
)


##############################################################
# 46B. SELECCIÓN DEL MÓDULO OBJETIVO
##############################################################

best.module.name <- "purple" # Fijado directamente a purple o al menor p-valor si cambia
module.genes <- names(moduleColors[moduleColors == best.module.name])

cat(
  "\nMódulo seleccionado para análisis de interacción:", best.module.name,
  "\nTotal de genes en el módulo:", length(module.genes), "\n"
)


##############################################################
# 47. INTEGRACIÓN DE HUB GENES + DATOS DIFERENCIALES DE DESEQ2
##############################################################

# Calcular Module Membership (MM)
MM <- cor(datExpr, MEs, use = "p")
MM.best <- MM[module.genes, paste0("ME", best.module.name)]

# Extraer resultados de DESeq2 para la interacción
res.int.df <- as.data.frame(res.list$Interaction)

# Construcción de la tabla maestra del módulo
hub.table <- data.frame(
  gene_id = module.genes,
  module = rep(best.module.name, length(module.genes)),
  moduleMembership = as.numeric(MM.best),
  absMM = abs(as.numeric(MM.best))
)

# Cruzar con DESeq2
hub.table$deseq2_log2FC_interaction <- res.int.df[hub.table$gene_id, "log2FoldChange"]
hub.table$deseq2_stat_interaction   <- res.int.df[hub.table$gene_id, "stat"]
hub.table$deseq2_pvalue_interaction <- res.int.df[hub.table$gene_id, "pvalue"]
hub.table$deseq2_padj_interaction   <- res.int.df[hub.table$gene_id, "padj"]

# Criterio de DEG de Interacción real (padj DESeq2 < 0.05)
hub.table$is_interaction_DEG <- !is.na(hub.table$deseq2_padj_interaction) & hub.table$deseq2_padj_interaction < 0.05

# Ordenar por conectividad dentro del módulo (absMM)
hub.table <- hub.table[order(-hub.table$absMM), , drop = FALSE]

write.csv(
  hub.table,
  file.path(outpath, "WGCNA_hub_genes.csv"),
  row.names = FALSE
)


##############################################################
# 48. GRÁFICO GS vs MM (USANDO ESTADÍSTICO DE DESEQ2)
##############################################################

pdf(
  file = file.path(outpath, paste0("WGCNA_GS_vs_MM_", best.module.name, ".pdf")),
  width = 9,
  height = 7
)

par(mfrow = c(1, 1))
plot(
  hub.table$absMM,
  abs(hub.table$deseq2_stat_interaction),
  pch = 19,
  col = ifelse(hub.table$is_interaction_DEG, "turquoise3", "coral1"),
  xlab = "Absolute Module Membership (|MM|)",
  ylab = "Absolute DESeq2 Interaction Statistic (|stat|)",
  main = paste("Módulo:", best.module.name, "- Interacción CR95 vs Catuai")
)
legend(
  "topleft",
  legend = c("Interaction DEG (padj < 0.05)", "No significativo"),
  col = c("turquoise3", "coral1"),
  pch = 19,
  title = "Significancia DESeq2"
)

dev.off()


##############################################################
# 49. CYTOSCAPE EXPORT PREPARATION
##############################################################

top.hubs <- head(hub.table$gene_id, 150)
module.interaction.degs <- intersect(module.genes, interaction.genes)

network.genes <- unique(c(top.hubs, module.interaction.degs))
max.network.genes <- 200

if (length(network.genes) > max.network.genes) {
  network.genes <- head(
    hub.table$gene_id[hub.table$gene_id %in% network.genes],
    max.network.genes
  )
}

cat("\nGenes exported to Cytoscape:", length(network.genes), "\n")


##############################################################
# 50. TOM SIMILARITY
##############################################################

TOM.network <- TOMsimilarityFromExpr(
  datExpr[, network.genes, drop = FALSE],
  power = soft.power,
  networkType = "signed"
)

dimnames(TOM.network) <- list(network.genes, network.genes)

tom.values <- TOM.network[upper.tri(TOM.network)]
tom.threshold <- as.numeric(quantile(tom.values, probs = 0.95, na.rm = TRUE))


##############################################################
# 51. CYTOSCAPE FILES EXPORT
##############################################################

exportNetworkToCytoscape(
  TOM.network,
  edgeFile = file.path(outpath, paste0("Cytoscape_edges_", best.module.name, ".txt")),
  nodeFile = file.path(outpath, paste0("Cytoscape_nodes_", best.module.name, ".txt")),
  weighted = TRUE,
  threshold = tom.threshold,
  nodeNames = network.genes,
  nodeAttr = rep(best.module.name, length(network.genes))
)


##############################################################
# 52. CYTOSCAPE NODE ATTRIBUTES
##############################################################

cyt.node.attrs <- hub.table[hub.table$gene_id %in% network.genes, , drop = FALSE]
cyt.node.attrs <- cyt.node.attrs[match(network.genes, cyt.node.attrs$gene_id), , drop = FALSE]

cyt.node.attrs$nodeName <- cyt.node.attrs$gene_id

adjacency <- TOM.network > tom.threshold
diag(adjacency) <- FALSE

cyt.node.attrs$network_degree <- rowSums(adjacency)

write.table(
  cyt.node.attrs,
  file = file.path(outpath, paste0("Cytoscape_node_attributes_", best.module.name, ".txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


##############################################################
# 53. FINAL SUMMARY
##############################################################

# Limpieza final de dispositivos gráficos por seguridad
while (dev.cur() > 1) dev.off()
graphics.off()

cat("\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n")
cat("\nWGCNA selected interaction module:", best.module.name)
cat("\n\nResults saved in:", outpath)
cat("\n============================================\n")
##############################################################
# CONVERSIÓN DIRECTA A CSV (SOLUCIÓN AL ERROR DE SCAN)
##############################################################

# Cargar el archivo usando read.csv / read.delim de forma robusta
edges_file <- file.path(outpath, "Cytoscape_edges_purple.txt")

# Intentar lectura con delimitador de tabulación/espacio automático
edges_df <- read.csv(edges_file, sep = "", header = TRUE, check.names = FALSE)

# Guardar como CSV limpio separado por comas
write.csv(
  edges_df, 
  file = file.path(outpath, "Cytoscape_edges_purple.csv"), 
  row.names = FALSE
)

# Convertir también la tabla de nodos si existe
nodes_file <- file.path(outpath, "Cytoscape_nodes_purple.txt")
if (file.exists(nodes_file)) {
  nodes_df <- read.csv(nodes_file, sep = "", header = TRUE, check.names = FALSE)
  write.csv(
    nodes_df, 
    file = file.path(outpath, "Cytoscape_nodes_purple.csv"), 
    row.names = FALSE
  )
}

cat("¡Proceso completado! Se generaron exitosamente los archivos .csv\n")

library(RCy3)
# 2. Verificar conexión con Cytoscape abierto
cytoscapePing()
# 1. Cargar Datos con el nombre CORRECTO del archivo
edges_df <- read.csv(file.path(outpath, "Cytoscape_edges_purple.csv"))
nodes_df <- read.csv(file.path(outpath, "Cytoscape_nodes_purple.csv"))

# 2. Filtrar Top 100 interacciones más fuertes
top_edges <- edges_df %>%
  arrange(desc(weight)) %>%
  head(100)

connected_nodes <- unique(c(top_edges$fromNode, top_edges$toNode))
filtered_nodes <- nodes_df %>%
  filter(nodeName %in% connected_nodes)

# 3. Reinyectar la tabla completa de atributos a la red actual
loadTableData(
  data = filtered_nodes,
  data.key.column = "nodeName",
  table.key.column = "name"
)

# 4. Crear y Aplicar Estilo Visual de Publicación (Sin tocar coordenadas/layout)
style_name <- "WGCNA_Publication_Style"

if (style_name %in% getVisualStyleNames()) {
  deleteVisualStyle(style_name)
}
createVisualStyle(style_name)

# Etiquetas: Mostrar el nombre de los genes (nodeName)
setNodeLabelMapping("nodeName", style.name = style_name)

# Tamaño: Mapeado a absMM (Hub Status)
if ("absMM" %in% colnames(filtered_nodes)) {
  setNodeSizeMapping(
    table.column = "absMM",
    mapping.type = "c",
    sizes = c(35, 90),
    style.name = style_name
  )
}

# Color: Mapeado a DEGs de Interacción (Rojo = VÁLIDO, Azul = NO DEG)
if ("is_interaction_DEG" %in% colnames(filtered_nodes)) {
  setNodeColorMapping(
    table.column = "is_interaction_DEG",
    table.column.values = c(TRUE, FALSE),
    colors = c("#FF4B4B", "#4A90E2"),
    mapping.type = "d",
    style.name = style_name
  )
}

# Grosor de bordes: Mapeado al peso de correlación
setEdgeLineWidthMapping(
  table.column = "weight",
  mapping.type = "c",
  widths = c(1, 5),
  style.name = style_name
)

# Aplicar el estilo y poner fondo blanco limpio
setVisualStyle(style_name)
setBGColorBypasser("#FFFFFF")

# Versiones de los paquetes
packageVersion("WGCNA")
packageVersion("dynamicTreeCut")
packageVersion("RCy3")
packageVersion("DESeq2") # Ya que es clave en tu flujo

# Si deseas que R te imprima la cita en formato BibTeX o texto plano:
citation("WGCNA")
citation("RCy3")
citation ("dynamicTreeCut")
