#Loading the necessary libraries.
library(dplyr)
library(gplots) 

#Set the working directory for the analysis.
setwd('/usr4/bf528/yamuna21/R/Project_1/')

#Reading the data frame for the analysis.
df <- read.csv('combatted_data.csv',sep=",",row.names=1)
message("Data frame loaded.")
message("Total genes : ",dim(df)[1])


#Filter 1 : Filter to retain intensity values greater than log2(15) and expression in at least 20% of the samples.
filter_1 <- filter(df,rowSums(df > log2(15)) >= (0.2*ncol(df)))
message("Applying first filter.")
message("Total genes remaining after 1st filter : ",dim(filter_1)[1])

#Calculating degrees of freedom.
dof <- ncol(filter_1) - 1
message("Calculating degrees of freedom.")
message("Degrees of Freedom : ",dof)
#Calculating variance for every gene.
message("Calculating variance for every gene.")
filter_1$variance <- apply(filter_1,1,var)

#Computing median variance.
message("Computing median variance.")
median_variance <- median(filter_1$variance)
message("Median Variance : ",median_variance)

#Calculating variance test statistic.
tstat_vec <- dof * (var(as.double(filter_1[1,1:(ncol(filter_1)-1)]))/median_variance)
message("Calculating variance test statistic.")

for (x in 2:(nrow(filter_1)))
{
  tstat <- dof * (var(as.double(filter_1[x,1:(ncol(filter_1)-1)]))/median_variance)
  tstat_vec <- c(tstat_vec,tstat)
}

filter_1$t_stat <- tstat_vec

#Calculating chi-square distribution.
message("Calculating chi-square distribution.")
chilower <- qchisq((0.01)/2, dof)
chiupper <- qchisq((1-0.01)/2, dof, lower.tail = FALSE)
message("Chi Lower : ",chilower)
message("Chi Upper : ",chiupper)

#Applying second filter.
message("Applying second filter.")
filter_2 <- filter(filter_1,filter_1$t_stat > chiupper | filter_1$t_stat < chilower)
message("Total genes remaining after 2nd filter : ",dim(filter_2)[1])

filter_2_t_test <- filter_2

#Writing .csv file after the second filter.
#write.csv(filter_2[,1:(ncol(filter_2)-2)],file="4_2_2Filters")
write.csv(filter_2[,1:(ncol(filter_2)-2)],file="4_2_2Filters")
message("Writing CSV file.")

#Calculating the coefficient of variation for every gene.
message("Calculating the coefficient of variation for every gene.")
coeff_vec <- (sqrt(filter_2[1,]$variance)) / (mean(as.double(filter_2[1,1:(ncol(filter_2)-2)])))

for (x in 2:(nrow(filter_2)))
{
  coeff <- (sqrt(filter_2[x,]$variance)) / (mean(as.double(filter_2[x,1:(ncol(filter_2)-2)])))
  coeff_vec <- c(coeff_vec,coeff)
}

filter_2$cv <- coeff_vec

#Applying third filter.
message("Applying third filter.")
filter_3 <- filter(filter_2,filter_2$cv>0.186)
message("Total genes remaining after 3rd filter : ",dim(filter_3)[1])

#Writing .csv file after the third filter.
write.csv(filter_3[,1:(ncol(filter_3)-3)],file="4_4_3Filters")
message("Writing CSV file.")

#Performing hierarchical clustering.
message("Performing hierarchical clustering.")
cluster <- hclust(dist(t(filter_3[,1:(ncol(filter_3)-3)])))

#Separating clusters in two.
message("Separating clusters.")
cluster_separated <- cutree(cluster, 2)

#Loading metadata for annotation.
message("Loading metadata for annotation.")
metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')

#Splitting columns at "_" to rename.
col_split <- do.call("rbind",strsplit(colnames(filter_3), split="_"))

#Renaming columns.
colnames(filter_3) <- col_split[,1]

metadata <- subset(metadata,select = c(geo_accession,cit.coloncancermolecularsubtype))

#Checking for sample names in metadata file for subtype annotation.
temp <- metadata[metadata$geo_accession %in% colnames(filter_3),]

#Color coding sample names according to molecular subtype.
colors <- ifelse(temp$cit.coloncancermolecularsubtype == 'C3','red','blue')

#Plotting heatmap.
message("Plotting heatmap.")
heatmap.2(as.matrix(filter_3[,1:(ncol(filter_3)-3)]), xlab='Sample', ylab='Gene', 
          main='Gene Expression',
          trace='none', density.info = 'none',ColSideColors = colors,
          key.xlab='Gene Expression Level', scale='row', margins=c(7,7))

coords <- data.frame(x=130.5953,y=79.84456)

legend(coords,legend = c("C3", "Others"),col = c("Red", "Blue"),lty= 1,lwd = 3,xpd=TRUE,cex=0.50)

#Separating data into 2 groups.
message("Separating data into two groups.")
group_1 <- filter_3[,cluster_separated==1]
group_2 <- filter_3[,cluster_separated==2]

message("Number of samples in group 1 : ",dim(group_1)[2]-2)
message("Number of samples in group 2 : ",dim(group_2)[2]-1)
message("Performing Welch's t-test.")
#Performing Welch's t-test.
p_values <- t.test(group_1[1,1:(ncol(group_1)-2)],group_2[1,1:(ncol(group_2)-1)])$p.value
t_test_values <- as.double(t.test(group_1[1,1:(ncol(group_1)-2)],group_2[1,1:(ncol(group_2)-1)])$statistic)

for (x in 2:(nrow(group_1)))
{
  p_value <- t.test(group_1[x,1:(ncol(group_1)-2)],group_2[x,1:(ncol(group_2)-1)])$p.value
  p_values <- c(p_values,p_value)
  t_test_value <- as.double(t.test(group_1[x,1:(ncol(group_1)-2)],group_2[x,1:(ncol(group_2)-1)])$statistic)
  t_test_values <- c(t_test_values,t_test_value)
}

group_1$p_values <- p_values
group_1$t_test_statistic <- t_test_values

#Adjusting p-values using FDR.
message("Adjusting p-values using FDR.")
group_1$p_adj <- p.adjust(group_1$p_values,method="fdr")

#Writing .csv after adjusting p-values.
message("Writing CSV.")
write.csv(select(group_1,t_test_statistic,p_values,p_adj),file="5_4_Welch")

#Performing same analysis for t-test in the dataset after 2nd filter.
filter_2_clust <- hclust(dist(t(filter_2_t_test[,1:(ncol(filter_2_t_test)-2)])))

filter_2_clust_sep <- cutree(filter_2_clust, 2)

filter_2_group_1 <- filter_2_t_test[,filter_2_clust_sep==1]
filter_2_group_2 <- filter_2_t_test[,filter_2_clust_sep==2]

p_vals <- t.test(filter_2_group_1[1,1:(ncol(filter_2_group_1)-1)],filter_2_group_2[1,1:(ncol(filter_2_group_2)-1)])$p.value
t_test_vals <- as.double(t.test(filter_2_group_1[1,1:(ncol(filter_2_group_1)-1)],filter_2_group_2[1,1:(ncol(filter_2_group_2)-1)])$statistic)

for (x in 2:(nrow(filter_2_group_1)))
{
  p_val <- t.test(filter_2_group_1[x,1:(ncol(filter_2_group_1)-1)],filter_2_group_2[x,1:(ncol(filter_2_group_2)-1)])$p.value
  p_vals <- c(p_vals,p_val)
  t_test_val <- as.double(t.test(filter_2_group_1[x,1:(ncol(filter_2_group_1)-1)],filter_2_group_2[x,1:(ncol(filter_2_group_2)-1)])$statistic)
  t_test_vals <- c(t_test_vals,t_test_val)
}

filter_2_group_1$p_values <- p_vals
filter_2_group_1$t_test_statistic <- t_test_vals

filter_2_group_1$p_adj <- p.adjust(filter_2_group_1$p_values,method="fdr")

#Writing .csv for the t-test after second filter.
message("Writing CSV.")
write.csv(select(filter_2_group_1,t_test_statistic,p_values,p_adj),file="4_5_Welch")

#Displaying genes after filtering adjusted p-values.
adj_p_filter <- filter(group_1,group_1$p_adj<0.05)
message("Number of genes left after adjusted p-value filter : ",dim(adj_p_filter)[1])

getwd()
setwd()
