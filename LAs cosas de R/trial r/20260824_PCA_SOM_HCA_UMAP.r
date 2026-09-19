#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")
#install.packages(c("kohonen","ggplot2","gplots","VennDiagram","pheatmap","dendsort","DESeq2","uwot"));

library(kohonen);	#This is the library for the SOM
library(ggplot2);	#This library is for transparency in the colors
library(gplots);	#Easy heatmaps
library(VennDiagram);	#self explanatory
library(pheatmap);	#pretty heatmaps
library(dendsort);	#sorting dendrograms
library(DESeq2);	#Normalization and everything related to that
library(uwot);      #umap

##Colors schemes verbatim from Peter Carl-s code on: http://www.r-bloggers.com/the-paul-tol-21-color-salute/
if(TRUE){
# Function for plotting colors side-by-side
pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

# FOCUS PALETTES
# Red as highlight
redfocus = c("#CB181D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
 
# Green as highlight
greenfocus = c("#41AB5D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
 
# Blue as highlight
bluefocus = c("#0033FF", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
 

# EQUAL WEIGHT
# Generated with rainbow(12, s = 0.6, v = 0.75)
rainbow12equal = c("#BF4D4D", "#BF864D", "#BFBF4D", "#86BF4D", "#4DBF4D", "#4DBF86", "#4DBFBF", "#4D86BF", "#4D4DBF", "#864DBF", "#BF4DBF", "#BF4D86")
rainbow10equal = c("#BF4D4D", "#BF914D", "#A8BF4D", "#63BF4D", "#4DBF7A", "#4DBFBF", "#4D7ABF", "#634DBF", "#A84DBF", "#BF4D91")
rainbow8equal = c("#BF4D4D", "#BFA34D", "#86BF4D", "#4DBF69", "#4DBFBF", "#4D69BF", "#864DBF", "#BF4DA3")
rainbow6equal = c("#BF4D4D", "#BFBF4D", "#4DBF4D", "#4DBFBF", "#4D4DBF", "#BF4DBF")
 
# Generated with package "gplots" function rich.colors(12)
rich12equal = c("#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466", "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300")
rich10equal = c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300")
rich8equal = c("#000041", "#0000CB", "#0081FF", "#02DA81", "#80FE1A", "#FDEE02", "#FFAB00", "#FF3300")
rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
 
# Generated with package "fields" function tim.colors(12), which is said to emulate the default matlab colorset
tim12equal = c("#00008F", "#0000EA", "#0047FF", "#00A2FF", "#00FEFF", "#5AFFA5", "#B5FF4A", "#FFED00", "#FF9200", "#FF3700", "#DB0000", "#800000")
tim10equal = c("#00008F", "#0000FF", "#0070FF", "#00DFFF", "#50FFAF", "#BFFF40", "#FFCF00", "#FF6000", "#EF0000", "#800000")
tim8equal = c("#00008F", "#0020FF", "#00AFFF", "#40FFBF", "#CFFF30", "#FF9F00", "#FF1000", "#800000")
tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
 
# Generated with sort(brewer.pal(8,"Dark2")) #Dark2, Set2
dark8equal = c("#1B9E77", "#666666", "#66A61E", "#7570B3", "#A6761D", "#D95F02", "#E6AB02", "#E7298A")
dark6equal = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")
set8equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#B3B3B3", "#E5C494", "#E78AC3", "#FC8D62", "#FFD92F")
set6equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#E78AC3", "#FC8D62", "#FFD92F")
 

# MONOCHROME PALETTES
# sort(brewer.pal(8,"Greens"))
redmono = c("#99000D", "#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1", "#FEE0D2", "#FFF5F0")
greenmono = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0", "#E5F5E0", "#F7FCF5")
bluemono = c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")
grey8mono = c("#000000","#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
grey6mono = c("#242424", "#494949", "#6D6D6D", "#929292", "#B6B6B6", "#DBDBDB")
 

# Qualitative color schemes by Paul Tol
 tol1qualitative=c("#4477AA")
 tol2qualitative=c("#4477AA", "#CC6677")
 tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
 tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
 tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
 tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
 tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
 tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
 tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
 tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
 tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
 tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
 

tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
# ...and finally, the Paul Tol 21-color salute
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, method = "pearson"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
sort_hclust <- function(x,...){
	as.hclust(dendsort(as.dendrogram(x),...));
}

#set.seed(42);

setwd(choose.dir());    #Hope this works also in mac

counts<-read.delim("Saccharomyces.txt", stringsAsFactors=FALSE, row.names=1);
counts<-counts[rowSums(counts)>0,];	#Remove all genes that have no expression

strain<-as.factor(sapply(strsplit(colnames(counts),split="_"),"[[",1));
names(strain)<-colnames(counts);

trt<-as.factor(sapply(strsplit(colnames(counts),split="_"),"[[",2));
names(trt)<-colnames(counts);


coldata<-data.frame(strain=strain, trt=trt);
rownames(coldata)<-colnames(counts);

dds<-DESeqDataSetFromMatrix(countData=as.matrix(round(counts)), colData=coldata, design = ~trt);	
#dds<-DESeqDataSetFromMatrix(countData=as.matrix(round(counts)), colData=coldata, design = ~trt+strain);	
#dds<-DESeqDataSetFromMatrix(countData=as.matrix(round(counts)), colData=coldata, design = ~trt+strain+trt*strain);	

cpm<-fpm(dds, robust = FALSE);	

mat.eval<-cpm>1; 

aux.eval<-paste0(strain,"_",trt);

keep<-array(FALSE,nrow(mat.eval));
cutoff<-ceiling(table(aux.eval)/2);

for (i in unique(aux.eval)){
	keep<-keep|(rowSums(mat.eval[,which(aux.eval==i)])>=cutoff[i]);
}

sum(keep)/length(keep)

dds.1<-dds[keep,];
dds.1 <- estimateSizeFactors(dds.1);	#You have to estimate size factors

vst.1 <- vst(dds.1);

std.mat<-t(scale(t(assay(vst.1))));

pca<-princomp(std.mat);
barras<-100*(pca$sdev)^2/sum((pca$sdev)^2);
barplot(barras, las=2);
x11();

col.trt<-as.character(trt);
col.trt[col.trt=="sc"]<-"blue";
col.trt[col.trt=="sl"]<-"#FFA500";

pch.strain<-as.character(strain);
pch.strain[pch.strain=="wt"]<-"16";
pch.strain[pch.strain=="st"]<-"15";
pch.strain<-as.numeric(pch.strain);

for (i in 1:4){
	for (j in (i+1):5) {
		plot(pca$loadings[,c(i,j)], main = "PCA Loadings", xlab=paste0("PC ",i," (",round(barras[i],1),"%)"), ylab=paste0("PC ",j," (",round(barras[j],1),"%)"), type = "n"); 
		text(pca$loadings[,c(i,j)], as.character(strain), col = col.trt); 
        #points(pca$loadings[,c(i,j)], col = col.trt, pch = pch.strain); 
        legend("bottom", legend = c("sc","sl"), fill = c("blue","red"), bty = "n");
		#x11()
	}
}

i<-1
j<-2
plot(pca$scores[,c(i,j)], main = "PCA Scores", xlab=paste0("PC ",i," (",round(barras[i],1),"%)"), ylab=paste0("PC ",j," (",round(barras[j],1),"%)"), pch="."); 

el.nombre.que.quieran<-pca$scores;
write.csv(el.nombre.que.quieran,file="El_nombre_que_sea.csv")
write.csv(pca$loadings,file="SC_loadings.csv")


##Diff. expression
dds.1 <- DESeq(dds.1, test="LRT", reduced = ~1);
#dds.trt <- DESeq(dds.1, test="LRT", reduced = ~strain); 	#You can control by strain, treatment, both or neither
#dds.strain <- DESeq(dds.1, test="LRT", reduced = ~trt);
#dds.strain.noInter <- DESeq(dds.1, test="LRT", reduced = ~trt+trt*strain); #Only changes in strain, no interaction
#dds.trt.noInter <- DESeq(dds.1, test="LRT", reduced = ~strain+trt*strain); #Only changes in treatment, no interaction
#dds.noInter <- DESeq(dds.1, test="LRT", reduced = ~trt*strain); #Changes in treatment or strain, no interaction
#dds.Inter <- DESeq(dds.1, test="LRT", reduced = ~strain+trt); #Only changes due to interaction
#And so on...

mat.cpm<-fpm(dds.trt);

p.adj<-results(dds.trt)$padj;
names(p.adj)<-rownames(mat.cpm);

sum(is.na(p.adj))
p.adj[is.na(p.adj)]<-1;

sum(p.adj<0.01)
#[1] 4809
sum(p.adj<0.01)/length(p.adj)
#[1] 0.774006

DEG.cpm<-mat.cpm[p.adj<0.01,];

log2<-log(DEG.cpm+1,2);
cs.log2<-t(scale(t(log2)));


som.xp<-som(cs.log2,grid = somgrid(10, 10, 	"hexagonal", toroidal=TRUE));
classif<-som.xp$unit.classif;
names(classif)<-rownames(cs.log2);

plot(som.xp)
plot(som.xp,type="mapping", pch = ".");	
plot(som.xp,type="changes", pch = ".");	
plot(som.xp,type="counts", pch = ".");	
plot(som.xp,type="dist.neighbours", pch = ".");	
plot(som.xp,type="quality", pch = ".");	



ann.colors <- list(	trt = c(sc="blue",sl="red"),
					strain = c(wt="black",st="hotpink"))

clust.sample<-hclust(dist(t(som.xp$codes[[1]]), method="manhattan"), method="complete")

plot(clust.sample, main = "Manhattan")

clust.sample<-hclust(dist(t(som.xp$codes[[1]]), method="euclidean"), method="ward.D2")

plot(clust.sample, main = "Euclidean - Ward")

clust.neuron<-hclust(dist(som.xp$codes[[1]], method="manhattan"), method="ward.D2")

plot(clust.neuron)

pheatmap(t(som.xp$codes[[1]]), border_color = "grey60", scale = "column", 
		show_rownames = T, show_colnames = F,
		cluster_rows = clust.sample, cutree_rows = 2, annotation_row = data.frame(trt=trt, strain=strain),
		cluster_cols = clust.neuron, cutree_cols = 3,
		annotation_colors = ann.colors);
		
       
clust.pat<-cutree(clust.neuron,3);

table(clust.pat)

clust.aux<-paste0("C",clust.pat);
names(clust.aux)<-names(clust.pat)

ann.colors <- list(	trt = c(sc="blue",sl="red"),
                    strain = c(wt="black",st="hotpink"),
                    clust.aux = c(C1="green",C2="orange",C3="cyan"))


pheatmap(t(som.xp$codes[[1]]), border_color = "grey60", scale = "column", 
         show_rownames = T, show_colnames = F,
         cluster_rows = clust.sample, cutree_rows = 2, annotation_row = data.frame(trt=trt, strain=strain),
         cluster_cols = clust.neuron, cutree_cols = 3, annotation_col = data.frame(clust.aux=clust.aux),
         annotation_colors = ann.colors);




type.pattern<-clust.pat[classif];
names(type.pattern)<-names(classif);

table(type.pattern)



############################## UMAP

#We can run a umap to see explore the dataset

UMAP_test<-umap(std.mat,min_dist=0.5, spread=0.5, ret_model=TRUE);

plot(UMAP_test$embedding, pch = 16, main = "TEST - min_dist=0.5, spread=0.5", col = alpha("black",0.05), cex = 0.1, xlab = "UMAP 1", ylab = "UMAP 2");  

#But which parameters are the best ones?
umap.par<-expand.grid(spread=c(1/8,1/4,1/2,1), min_dist=c(1/16,1/8,1/4,1/2));		#We can try different combinations to decide which one gives the best projection
UMAP<-vector("list",nrow(umap.par));
for (i in 1:nrow(umap.par)){
	UMAP[[i]]<-tryCatch(umap(std.mat,min_dist=umap.par[i,"min_dist"], spread=umap.par[i,"spread"], ret_model=TRUE), error = function(e) {
		return(NA)
	})
}

#Now to plot everything
for (i in 1:length(UMAP)){
	if (any(!is.na(UMAP[[i]]))){
		plot(UMAP[[i]]$embedding, pch = 16, main = paste0("min_dist=",umap.par[i,"min_dist"],"; spread=",umap.par[i,"spread"]), col = alpha("black",0.05), cex = 0.1, xlab = "UMAP 1", ylab = "UMAP 2");  
		x11();
	}
}

