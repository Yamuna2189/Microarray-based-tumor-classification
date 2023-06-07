## 1. install packages 'Bioconductor' and some libraries(R 4.0.3 version)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c("affy","affyPLM","sva",'AnnotationDbi','hgu133plus2.db'))

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)

## 2. Read CEL files and Normalize
datdir <- ReadAffy(celfile.path = "/projectnb/bf528/users/glass_bottom_boat_2022/project_1/samples")
hist(datdir, names = NULL, xlab = "Samples", ylab = "log-intensity", main = "Raw Data")
dat_norm <- rma(datdir) # rma algorithm

## 3. Plot RLE and NUSE
dat_plm <- fitPLM(datdir, normalize = T, background = T)
result_rle <- RLE(dat_plm, type = "stats")
df_rle <- as.data.frame(t(result_rle))
plot_rle <- ggplot(df_rle, aes(x=median)) +geom_histogram(color="white", fill="navy blue",bins = 35)+labs(title="Median RLE Scores",x="RLE Score", y = "Count")+theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=-1))
plot_rle

result_nuse <- NUSE(dat_plm, type = "stats")
df_nuse <- as.data.frame(t(result_nuse))
plot_nuse <- ggplot(df_nuse, aes(x=median)) + geom_histogram(color="black", fill="maroon", bins=35)+labs(title="Median NUSE Scores",x="NUSE Score",y = "Count")+theme(plot.title=element_text(size=18, hjust=0.5, face="bold", colour="black", vjust=-1))
plot_nuse

## 4. csv file for RMA normalized
eset <- exprs(dat_norm)
write.csv(eset, file= 'normalised_data.csv')

## 5. Batch effects Correction
metadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
mod <- model.matrix(~normalizationcombatmod, data = metadata)
combat_data <- ComBat(dat = eset, batch = metadata$normalizationcombatbatch, mod = mod)
write.csv(combat_data, file= 'combatted_data.csv')
getwd()

## 6. PCA 
combat_adj_t <-t(combat_data)
combat_scaled <-t(scale(combat_adj_t, center = TRUE, scale = TRUE))
pca <-prcomp(combat_scaled, scale = FALSE, center = FALSE)
pca$rotation
pca_value <- as.data.frame(pca$rotation)

# examine the percent variability explained by each principal component
pca_summary <- summary(pca)
(pca_summary$importance)

ggplot(data = pca_value, mapping = aes(x = PC1, y = PC2)) +geom_point()+theme_bw() +labs(title = 'PCA plot', x= 'PC1=11.4%', y='PC2=8.4%')

## 7. Plotting Outliers in PC2
ggplot(data = pca_value, mapping = aes(y = PC1)) +geom_boxplot()+theme_bw() + labs(title = 'PCA1 Boxplot for detecting outliers')
ggplot(data = pca_value, mapping = aes(y = PC2)) +geom_boxplot()+theme_bw() +labs(title = 'PCA2 Boxplot for detecting outliers')
system.time(source('Programmer.R'))

