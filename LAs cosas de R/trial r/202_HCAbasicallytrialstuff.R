install.packages("dendextend");	#run this line to install dendextend library
library(dendextend);	#this library makes handling dendrograms easier in R

getwd();	#This will show you in which directory you are 

choose.dir();	#this function lets you choose a directory interactively, and returns the path
#it oly works on windows

setwd(choose.dir());	#If you nest it with setwd, you can choose which directory to work in; 
#again, this only works on windows
#in R studio, you can actually set your directory in a user-friendly manner

iris<-read.csv("Iris Dataset.csv");		#this is the file we saved from excel
#IF your excel uses a decimal comma (1/2 = 0,5 instead of 1/2 = 0.5) then you need to use 
#iris<-read.csv2("iris.csv");

#data(iris);			#Only use this if reading csvs doesn't work 

iris

summary(iris)	#R is very useful because it has a lot of integrated functions, and lots of libraries

boxplot(Sepal.length~Species, data = iris);			#iris is a data frame
boxplot(iris$Sepal.length~iris$Species);			#so you can access columns 
boxplot(iris[,"Sepal.length"]~iris[,"Species"]);	#in a lot of different ways
boxplot(iris[,2]~iris[,6]);							#Which will be useful later

plot(iris$Sepal.length, iris$Petal.length);			#You can do a lot of things quickly

color<-rainbow(3)[factor(iris$Species)];							#This only works because Species is a factor; but we are moving forward quickly
plot(iris$Sepal.length, iris$Petal.length, col = color);	#You can use this to find patterns in your data

pairs(iris[,c(-1,-6)], col = color);			#this is a nice dataset to explore everything


#But we are here for hierarchical clustering!!

#First we need to calculate the distances

#So we take just the measurements
iris.meas<-data.matrix(iris[,c(-1,-6)]);

#Calculate distances
iris.dist<-dist(iris.meas);

#and do hierarchical clustering
iris.clust<-hclust(iris.dist);

#Convert to a dendrogram
iris.dend<-as.dendrogram(iris.clust);

plot(iris.dend);

labels(iris.dend)<-iris$Species;		#This is complicated without the dendextend library

plot(iris.dend);
#is it clustering by species?

abline(h=3.5, col = "red");

##OK, so if you remember, there are different ways of clustering
##which one is R using?

?hclust

#Lets try a different one
iris.clust.ward<-hclust(iris.dist, method = "ward.D2");

iris.dend.ward<-as.dendrogram(iris.clust.ward);
labels(iris.dend.ward)<-iris$Species;		#This is complicated withowt the dendextend library
plot(iris.dend.ward);
abline(h=10, col = "red");					#OK, so it's getting better

#There are also different dissimilarity measurements

?dist

iris.dist<-dist(iris.meas, method = "manhattan");
iris.clust.ward<-hclust(iris.dist, method = "ward.D2");

iris.dend.ward<-as.dendrogram(iris.clust.ward);
labels(iris.dend.ward)<-iris$Species;		#This is complicated without the dendextend library
plot(iris.dend.ward);		#So it is not getting better

#We can also try transforming
log.iris<-log(iris.meas,2)
iris.dist<-dist(log.iris);
iris.clust.ward<-hclust(iris.dist, method = "ward.D2");

iris.dend.ward<-as.dendrogram(iris.clust.ward);
labels(iris.dend.ward)<-iris$Species;		
plot(iris.dend.ward);		#Much better
abline(h=5.1,col="red");


iris.grp<-cutree(iris.dend.ward, k=3); 

iris$Species[iris.grp==1]
iris$Species[iris.grp==2]
iris$Species[iris.grp==3]

table(iris$Species[iris.grp==2])

iris.grp<-cutree(iris.dend.ward, h=5.1); 

table(iris$Species[iris.grp==2])

#but why?
summary(iris)

##still, not the best way to see it
#Let's see the data

install.packages("pheatmap")
library(pheatmap);	#This is a library to make pretty heatmaps

pheatmap(log.iris);

pheatmap(t(log.iris)); #It does look pretty, but doesn't tell me much

?pheatmap	#Let's see what it can do

pheatmap(t(log.iris),scale = "row",clustering_distance_cols = "euclidean", clustering_method = "ward.D2")

pheatmap(t(log.iris),scale = "row",clustering_distance_cols = "euclidean", clustering_method = "ward.D2", clustering_distance_rows = "correlation", cutree_cols = 3)

pheatmap(t(log.iris),scale = "row",clustering_distance_cols = "euclidean", clustering_method = "ward.D2", clustering_distance_rows = "correlation", cutree_cols = 3, annotation_col = iris[,"Species"])

rownames(log.iris)<-rownames(iris)<-paste0("sample",1:nrow(log.iris));

pheatmap(t(log.iris),scale = "row",clustering_distance_cols = "euclidean", clustering_method = "ward.D2", clustering_distance_rows = "correlation", cutree_cols = 3, annotation_col = iris[,"Species"])

pheatmap(t(log.iris),scale = "row",clustering_distance_cols = "euclidean", clustering_method = "ward.D2", clustering_distance_rows = "correlation", cutree_cols = 3, annotation_col = iris[,"Species",drop=FALSE])

pheatmap(t(log.iris),scale = "none",clustering_distance_cols = "euclidean", clustering_method = "ward.D2", clustering_distance_rows = "correlation", cutree_cols = 3, annotation_col = iris[,"Species",drop=FALSE])

pheatmap(t(log.iris[,-2]),scale = "none",clustering_distance_cols = "euclidean", clustering_method = "ward.D2", clustering_distance_rows = "correlation", cutree_cols = 3, annotation_col = iris[,"Species",drop=FALSE])

pheatmap(t(log.iris[,-c(1,2)]),scale = "none",clustering_distance_cols = "euclidean", clustering_method = "ward.D2", clustering_distance_rows = "correlation", cutree_cols = 3, annotation_col = iris[,"Species",drop=FALSE])

