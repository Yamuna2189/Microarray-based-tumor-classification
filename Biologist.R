# Loading all required packages
library(affy)
library(hgu133plus2.db)
library(sva)
library(AnnotationDbi)
library(GSEABase)
library(dplyr)
library(tidyverse)
library(qusage)

#set directory 
setwd("/usr4/bf527/ryamani/Documents/bf528/project1/Biologist/")

#uplaod the data from the analyst part
sample_data <- read.csv("/usr4/bf527/ryamani/Documents/bf528/project1/Biologist/5_4_Welch", row.names = 1, header = TRUE)

#Order t-test by descending 
desc_t <- sample_data %>% arrange(desc(t_test_statistic))

#Map each probeset ID to gene symbol using AnnotationDbi
sample_keys <- AnnotationDbi::select(hgu133plus2.db, keys = (row.names(desc_t)), columns = ("SYMBOL"))


#Remove duplicate probe IDs and join columns differential expression table and table with mapped probesetID:Symbol 
remdup <- sample_keys[!duplicated(sample_keys[1]),]
adata <- cbind(remdup, desc_t)


#get smallest p values and filter by that
newsam <- adata %>%
  group_by(SYMBOL) %>%
  filter(p_adj == min(p_adj)) %>%
  ungroup(SYMBOL)


#Removes null values 
newsam <- newsam[!is.na(newsam$SYMBOL),]

#first answer
#top 1000 and then top 10 up and down regulated genes
top1000_up <- head(newsam, 1000)
top10_up <- head(newsam, 10)
top1000_down <- head(newsam %>% arrange(t_test_statistic), 1000)
top10_down <- head(newsam %>% arrange(t_test_statistic), 10)
write.csv(top10_up, file = "top10_up_results.csv", row.names = FALSE)
write.csv(top10_down, file = "top10_down_results.csv", row.names = FALSE)

#Gene sets and loading data
hallmarks <- getGmt("/usr4/bf527/ryamani/Documents/bf528/project1/h.all.v7.5.1.symbols.gmt")
kegg <- getGmt("/usr4/bf527/ryamani/Documents/bf528/project1/c2.cp.kegg.v7.5.1.symbols.gmt")
go <- getGmt("/usr4/bf527/ryamani/Documents/bf528/project1/c5.go.v7.5.1.symbols.gmt")

#second answer
#calculate the number of genesets in each databases
hallmark_len <- length(names(hallmarks)) #50
go_len <- length(names(go)) #10402
kegg_len <- length(names(kegg)) #186

# Obtain genes which were not expressed
notexp_up <- subset(newsam, !newsam$SYMBOL %in% top1000_up$SYMBOL)
notexp_down <- subset(newsam, !newsam$SYMBOL %in% top1000_down$SYMBOL)

# Defining fisher test function with genelist, geneset, nde= not differentially expressed
fishertest <- function(genelist, geneset, nde)           
{ de_in_gs <- length(intersect(genelist,geneset))    #de_in_gs: differentially expressed genes that are in geneset
de_notin_gs <- length(genelist) - de_in_gs           #de_notin_gs: differentially expressed genes that are not in geneset 
notde_in_gs <- length(intersect(nde,geneset))        #notde_in_gs: not expressed but in geneset
notde_notin_gs <- length(nde) - notde_in_gs          #notde_notin_gs: not differentially expressed and not in geneset
return(c(de_in_gs,de_notin_gs,notde_in_gs,notde_notin_gs))}   #return fishertest values
values

# data frame for results after fisher test
kegg_fisher <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)
go_fisher <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)
hallmark_fisher <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)


# fisher test values
# kegg
for (i in 1:length(kegg))
{
  geneid <- geneIds(kegg[i])
  fisher_up <- fishertest(top1000_up$SYMBOL, geneid[[names(geneid)]], notexp_up$SYMBOL)
  fisher_down <- fishertest(top1000_down$SYMBOL, geneid[[names(geneid)]], notexp_down$SYMBOL)
  upregulated <- fisher.test(matrix(fisher_up,nrow=2))
  downregulated <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_fisher[nrow(kegg_fisher) +1, ] <- c(names(geneid), upregulated$p.value, upregulated$estimate, 'Up')
  kegg_fisher[nrow(kegg_fisher) +1, ] <- c(names(geneid), downregulated$p.value, downregulated$estimate, 'Down')}

kegg_fisher <- kegg_fisher %>% mutate(pvalue = as.numeric(pvalue), est = as.numeric(estimate)) #7


# go
for (i in 1:length(go))
{
  geneid <- geneIds(go[i])
  fisher_up <- fishertest(top1000_up$SYMBOL, geneid[[names(geneid)]], notexp_up$SYMBOL)
  fisher_down <- fishertest(top1000_down$SYMBOL, geneid[[names(geneid)]], notexp_down$SYMBOL)
  upregulated <- fisher.test(matrix(fisher_up,nrow=2))
  downregulated <- fisher.test(matrix(fisher_down, nrow=2))
  go_fisher[nrow(go_fisher) +1, ] <- c(names(geneid), upregulated$p.value, upregulated$estimate, 'Up')
  go_fisher[nrow(go_fisher) +1, ] <- c(names(geneid), downregulated$p.value, downregulated$estimate, 'Down')}

go_fisher <- go_fisher %>% mutate(pvalue = as.numeric(pvalue), est = as.numeric(estimate))


# hallmarks
for (i in 1:length(hallmarks))
{
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(top1000_up$SYMBOL, geneid[[names(geneid)]], notexp_up$SYMBOL)
  fisher_down <- fishertest(top1000_down$SYMBOL, geneid[[names(geneid)]], notexp_down$SYMBOL)
  upregulated <- fisher.test(matrix(fisher_up,nrow=2))
  downregulated <- fisher.test(matrix(fisher_down, nrow=2))
  hallmark_fisher[nrow(hallmark_fisher) +1, ] <- c(names(geneid), upregulated$p.value, upregulated$estimate, 'Up')
  hallmark_fisher[nrow(hallmark_fisher) +1, ] <- c(names(geneid), downregulated$p.value, downregulated$estimate, 'Down')}

hallmark_fisher <- hallmark_fisher %>% mutate(pvalue = as.numeric(pvalue), est = as.numeric(estimate))

# Using the FDR method to adjust the pvalue in each fisher result
kegg_fisher$FDR <- p.adjust(kegg_fisher$pvalue, method = "BH", n = length(kegg_fisher$pvalue))
write.csv(kegg_fisher, "kegg_FDR.csv")

go_fisher$FDR <- p.adjust(go_fisher$pvalue, method = "BH", n = length(go_fisher$pvalue))
write.csv(go_fisher, "go_FDR.csv")   

hallmark_fisher$FDR <- p.adjust(hallmark_fisher$pvalue, method = "BH", n = length(hallmark_fisher$pvalue))
write.csv(hallmark_fisher, "hallmark_FDR.csv")

#answer 3
#answer 4
#Finding enriched genesets which are significant and calcuate them

# kegg

enriched_kegg <- kegg_fisher[kegg_fisher$pvalue<0.05,]
enriched_kegg_len <- length(enriched_kegg$setname) #15 #counts number of sig enriched gs
enriched_kegg_desc <- enriched_kegg %>% arrange(pvalue)
enriched_kegg3 <- head(enriched_kegg_desc, 3) #top 3 sig enriched gs
view(enriched_kegg_len) #15
write.csv(enriched_kegg3, file="Kegg_enriched_top3.csv")

# go

enriched_go <- go_fisher[go_fisher$pvalue<0.05,]
enriched_go_len <- length(enriched_go$setname) #607
enriched_go_desc <- enriched_go %>% arrange(pvalue)
enriched_go3 <- head(enriched_go_desc, 3)
view(enriched_go_len) #607
write.csv(enriched_go3, file="go_enriched_top3.csv")

# hallmarks

enriched_hallmarks <- hallmark_fisher[hallmark_fisher$pvalue<0.05,]
enriched_hallmarks_len <- length(enriched_hallmarks$setname) #6 
enriched_hallmarks_desc <- enriched_hallmarks %>% arrange(pvalue)
enriched_hallmarks3 <- head(enriched_hallmarks_desc, 3) 
view(enriched_hallmarks_len) #6
write.csv(enriched_hallmarks3, file="hallmarks_enriched_top3.csv")

