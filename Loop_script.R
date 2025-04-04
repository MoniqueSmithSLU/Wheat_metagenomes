##############################################################################
#
# Reference to ARTICLE
# Monique Smith, Sweden's Agricultural University (SLU)
# monique.smith@slu.se
# revisions 1st April 2025

# additional files needed: 
# K01202.biom,K02532.biom,K02847.biom,K03477.biom,K05874.biom,K05876.biom,K08325.biom,
# K08961.biom,K10235.biom,K10857.biom,K14266.biom,K16088.biom,K16210.biom,K17060.biom,
# K18576.biom,K18578.biom,K18651.biom,K18786.biom,K19731.biom,K20455.biom,K00211.biom,
# K00608.biom,K00998.biom,K01085.biom,K01352.biom,K01355.biom,K01819.biom,K02383.biom,
# K03181.biom,K03202.biom,K03397.biom,K04338.biom,K04643.biom,K05372.biom,K05660.biom,
# K06080.biom,K07345.biom,K07786.biom,K07991.biom,K08156.biom,K08566.biom,K10094.biom,
# K10926.biom,K10975.biom,K11009.biom,K11017.biom,K11734.biom,K11743.biom,K11889.biom,
# K11909.biom,K11934.biom,K11964.biom,K12059.biom,K12069.biom,K12083.biom,K12228.biom,
# K12281.biom,K12284.biom,K12285.biom,K12286.biom,K12455.biom,K12688.biom,K12820.biom,
# K13448.biom,K13454.biom,K13585.biom,K13586.biom,K13831.biom,K13964.biom,K14368.biom,
# K14627.biom,K14781.biom,K15968.biom,K16215.biom,K16347.biom,K16348.biom,K16397.biom,
# K16403.biom,K16448.biom,K16552.biom,K16553.biom,K16696.biom,K18373.biom,K18380.biom,
# K18642.biom,K18793.biom,K19033.biom,K19060.biom,K19101.biom,K19213.biom,K19611.biom,
# K19885.biom,K20086.biom,K20088.biom,K20090.biom,K20267.biom,K20268.biom,K20272.biom,
# K20275.biom,K20326.biom,K20555.biom,K20961.biom,K21212.biom,K21239.biom,K21280.biom,
# K21346.biom,K21825.biom,K22335.biom,K20966.biom,K19216.biom,K02399.biom,K10927.biom,
# K12280.biom
#sizeFactorsTAXA.csv - this is from the All_R_code.R script


# Purpose is to analysis the taxa behind each of the differentially abundant genes between Modern and Heritage wheat
# This script is for the supporting data and statistical analysis for 
# Table 3 and Table S2 - other results are in the script All_R_code.R

# David has extracted the reads and realigned the taxonomy and sent me the biom files from Megan
# Data has been extracted from MEGAN6 from the merged samples (large samples were split to meganise and then merged)
# Samples were aligned using AnnoTree so only contain Bacteria and Archaea 
# In Megan - I uncollapse all nodes, select all leaves and go file->export->biom

##############################################################################

# instal specialty packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
#### Load packages and data ####

library(phyloseq)
library(ggplot2)
library(DESeq2)
library("data.table")
library(microbiome)
#library("vsn")
library(psadd)
library("pheatmap")
library("RColorBrewer")
library(dplyr)
library(devtools)
library("ashr")
library(vegan)
library(plyr)
library(car)
library(ranacapa)
library(remotes)
library(psadd)
library(agricolae)
library(emmeans)
library("labdsv")
library(ggpubr)
library(microViz)
library(rio)

# The taxa table from Megan has some mixed classification. There are taxa names in the wrong column. 
# I don't think this needs to be fixed for this analysis, it was more important for plotting with the other data sets. 
#  Matt made a script to move everything into the right columns but I wont use that here

# First I ran the analysis on the top 20 most abundant DA genes but now I am rerunning with all 113 genes together so that p adjustments can be made on the full dataset
# and so the script is reproducible 


# Import the biom tables exactly as exported from MEGAN6 
library(biomformat)

#Create a list of the file paths to the BIOM files to be merged
file_paths <- c("K01202.biom","K02532.biom","K02847.biom","K03477.biom","K05874.biom","K05876.biom","K08325.biom",
                "K08961.biom","K10235.biom","K10857.biom","K14266.biom","K16088.biom","K16210.biom","K17060.biom",
                "K18576.biom","K18578.biom","K18651.biom","K18786.biom","K19731.biom","K20455.biom","K00211.biom",
                "K00608.biom","K00998.biom","K01085.biom","K01352.biom","K01355.biom","K01819.biom","K02383.biom",
                "K03181.biom","K03202.biom","K03397.biom","K04338.biom","K04643.biom","K05372.biom","K05660.biom",
                "K06080.biom","K07345.biom","K07786.biom","K07991.biom","K08156.biom","K08566.biom","K10094.biom",
                "K10926.biom","K10975.biom","K11009.biom","K11017.biom","K11734.biom","K11743.biom","K11889.biom",
                "K11909.biom","K11934.biom","K11964.biom","K12059.biom","K12069.biom","K12083.biom","K12228.biom",
                "K12281.biom","K12284.biom","K12285.biom","K12286.biom","K12455.biom","K12688.biom","K12820.biom",
                "K13448.biom","K13454.biom","K13585.biom","K13586.biom","K13831.biom","K13964.biom","K14368.biom",
                "K14627.biom","K14781.biom","K15968.biom","K16215.biom","K16347.biom","K16348.biom","K16397.biom",
                "K16403.biom","K16448.biom","K16552.biom","K16553.biom","K16696.biom","K18373.biom","K18380.biom",
                "K18642.biom","K18793.biom","K19033.biom","K19060.biom","K19101.biom","K19213.biom","K19611.biom",
                "K19885.biom","K20086.biom","K20088.biom","K20090.biom","K20267.biom","K20268.biom","K20272.biom",
                "K20275.biom","K20326.biom","K20555.biom","K20961.biom","K21212.biom","K21239.biom","K21280.biom",
                "K21346.biom","K21825.biom","K22335.biom","K20966.biom","K19216.biom","K02399.biom","K10927.biom",
                "K12280.biom" 
               )




#Load the BIOM files into a list of biom objects
biom_list <- lapply(file_paths, import_biom)

biom_list # 113 phyloseq objects
head(otu_table(biom_list[[3]]))
head(tax_table(biom_list[[90]]))


# remove the KEGG ID from the sample names so they are just S1, S2...
# first make an empty list
biom_list2=list(1:113)

for (i in c(1:113)) #The loop code will be repeated 20 times, based on the values that “I” refers to
{
  df=biom_list[[i]] #Select the ith dataset in your list
  # change the names of the samples
  old_names <- as.character(sample_names(df)) 
  #old_names2 <- gsub('Red_Llamas','Red.Llamas',old_names)
  #old_names3 <- gsub('Bulk_Soil','Bulk.Soil',old_names2)
  new_names <- strsplit(old_names,split ="_")
  new_names_mat <- matrix(unlist(new_names),ncol=2,byrow=T)
  colnames(new_names_mat) <- c("C1","C2") 
  col_vec <- "C2"
  col_mat <- new_names_mat[,col_vec] 
  sample_names(df) <- col_mat
  biom_list2[[i]]=df

}

# Check
sample_names(biom_list2[[100]])

# get metadata
md <- read.delim("Metadata_merged2.txt", header = TRUE, sep = "\t", dec = ".")

# add a new factor so figures don't refer to Heritage and modern as per reviewers comments
md$Sample_type <- "Bulk_Soil"
md[md$Sample_type %in% c('Modern'),]$Sample_type <- "Modern"
md[md$Sample_type %in% c('Heritage'),]$Sample_type <- "Heritage"

# meta data needs to have the samples as row names
rownames(md) <- md[,1]
md[,1] <- NULL
md

# add size factors to sample data
sizeFactors <- read.csv("sizeFactorsTAXA.csv", header = TRUE, row.names = 1)

#merge by row number
md_sf_merge <- merge(md, sizeFactors, 
                          by = 'row.names', all = TRUE) 
rownames(md_sf_merge) <- md_sf_merge[,1]
md_sf_merge[,1] <- NULL
md_sf_merge
# this is needed to tell R that md is sample data so it can be used in phyloseq
md2 <- sample_data(md_sf_merge)




# now merge the metadata into each biom

# first make an empty list
ps_list=list(1:113)

for (i in c(1:113)) #The loop code will be repeated 113 times, based on the values that “I” refers to
{
  df=biom_list2[[i]] #Select the ith dataset in your list
  # make phyloseq object
  
  ps_list[[i]]= merge_phyloseq(df, md2)
  
}

# check
ps_list

otu_table(ps_list[[65]])
head(tax_table(ps_list[[2]]))
sample_data(ps_list[[10]])

# some KEGGS have missing samples for the zero abundance
# It seems to work out that the missing samples are not needed even though they are in
# the metadata file and matches everything up by sample number 

# Make a list of names for later
names(ps_list)=c("K01202","K02532","K02847","K03477","K05874","K05876","K08325",
                 "K08961","K10235","K10857","K14266","K16088","K16210","K17060",
                 "K18576","K18578","K18651","K18786","K19731","K20455","K00211",
                 "K00608","K00998","K01085","K01352","K01355","K01819","K02383",
                 "K03181","K03202","K03397","K04338","K04643","K05372","K05660",
                 "K06080","K07345","K07786","K07991","K08156","K08566","K10094",
                 "K10926","K10975","K11009","K11017","K11734","K11743","K11889",
                 "K11909","K11934","K11964","K12059","K12069","K12083","K12228",
                 "K12281","K12284","K12285","K12286","K12455","K12688","K12820",
                 "K13448","K13454","K13585","K13586","K13831","K13964","K14368",
                 "K14627","K14781","K15968","K16215","K16347","K16348","K16397",
                 "K16403","K16448","K16552","K16553","K16696","K18373","K18380",
                 "K18642","K18793","K19033","K19060","K19101","K19213","K19611",
                 "K19885","K20086","K20088","K20090","K20267","K20268","K20272",
                 "K20275","K20326","K20555","K20961","K21212","K21239","K21280",
                 "K21346","K21825","K22335","K20966","K19216","K02399","K10927",
                 "K12280" 
)


# Prefiltering ####

# look at the phyla present 
phyla_list_all=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list[[i]] #Select the ith dataset in your list
  
  phyla_list_all[[i]] <- unique(df@tax_table@.Data[,"Rank2"])
  
}
phyla_list_all

# Phyla that stand out that should be renamed to unclassified
# NA, s__bacterium, s__Terrabacteria group bacterium ANGP1, s__bacterium M00.F.Ca.ET.205.01.1.1, s__bacterium JKG1

ps_list2=list(1:113)
for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list[[i]] #Select the ith dataset in your list
  ps_list2[[i]] <- tax_fix(df)
}
ps_list2

# check phyla names
phyla_list_all2=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list2[[i]] #Select the ith dataset in your list
  
  phyla_list_all2[[i]] <- unique(df@tax_table@.Data[,"Rank2"])
  
}
phyla_list_all2





#### Normalisation ####

ps_list_norm=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list2[[i]] #Select the ith dataset in your list
  sizeFactorsvector <- as.list(sample_data(df)[,4])
  # Normalise
  mat <- sweep(otu_table(df), 1 + taxa_are_rows(df), sizeFactorsvector$SizeFactors, FUN = `/`)
  ps_list_norm[[i]] <- df
  otu_table(ps_list_norm[[i]]) <- otu_table(mat, taxa_are_rows = taxa_are_rows(df))

}
  
ps_list_norm

#check 
#head of original

head_list=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list2[[i]] #Select the ith dataset in your list
  
  head_list[[i]] <- head(otu_table(df))
  
}
head_list[[1]]

# head of normalised data
head_list_norm=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list_norm[[i]] #Select the ith dataset in your list
  
  head_list_norm[[i]] <- head(otu_table(df))
  
}
head_list_norm[[1]]

# looks good but I need to round integers for some of the analysis


ps_list_norm_int=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list_norm[[i]] #Select the ith dataset in your list
  ps_list_norm_int[[i]] <- df
  otu_table(ps_list_norm_int[[i]]) <- round(otu_table(df))
  
}

ps_list_norm_int



# graph sequencing depths
# first make an empty list
graph_list2=list(1:113)
TotalReads_list=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list_norm_int[[i]] #Select the ith dataset in your list
  
  sdt = data.table(as(sample_data(df), "data.frame"),
                   TotalReads = sample_sums(df), keep.rownames = TRUE)
  setnames(sdt, "rn", "SampleID")
  graph_list2[[i]] = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle(names(ps_list)[i]) + facet_wrap(~Sample_type)
  
  TotalReads_list[[i]] <- aggregate(sdt$TotalReads, by=list(Sample_type=sdt$Sample_type), FUN=sum)
}

#ggarrange(plotlist=graph_list2) #Plots all graphs in the list
TotalReads_list

# Make a list of names for later
names(TotalReads_list)=names(ps_list)




library(dplyr)
library(tidyr)
# Combine the dataframes and reshape to wide format
combined_df <- TotalReads_list %>% 
  bind_rows(.id = "KEGG_ID") %>%
  pivot_wider(names_from = Sample_type, values_from = x)


write.csv(combined_df, "TotalReads_list_normalised_integer_data_113.csv", row.names=TRUE)





#### Sample differences ####

# Relative abundance ####

ps_list_rel=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list_norm_int[[i]] #Select the ith dataset in your list
  ps_list_rel[[i]] = transform_sample_counts(df, function(x) x / sum(x) )
}

ps_list_rel

# check
otu_table(ps_list_rel[[65]])
# there are NaN's for the KEGG IDs with samples with 0 values


# Change the NAN's to 0's
# Initialize an empty list to store the updated phyloseq objects
updated_ps_list <- list(1:113)

# Loop over each phyloseq object in the list
for (i in seq_along(ps_list_rel)) {
  # Extract the current phyloseq object
  ps <- ps_list_rel[[i]]
  
  # Extract the OTU table from the phyloseq object
  otu_table_matrix <- as(otu_table(ps), "matrix")
  
  # Replace NaN values with 0
  otu_table_matrix[is.nan(otu_table_matrix)] <- 0
  
  # Create a new OTU table
  new_otu_table <- otu_table(otu_table_matrix, taxa_are_rows = TRUE)
  
  # Replace the OTU table in the phyloseq object with the new OTU table
  ps@otu_table <- new_otu_table
  
  # Store the updated phyloseq object in the updated_ps_list
  updated_ps_list[[i]] <- ps
}


updated_ps_list

otu_table(updated_ps_list[[65]])
names(updated_ps_list) <- names(TotalReads_list)

# this seems to have fixed the issue. 

# PERMANOVA on relative abundances ####
# library(vegan)
#library("plyr")
sample_data(updated_ps_list[[65]])

# first make distance matrices
dist_matrices <- vector("list", 113)

for (i in 1:113) {
  
  ps_obj <- updated_ps_list[[i]]
  
  dist_matrices[[i]] <- phyloseq::distance(ps_obj, method="bray")
  # Check for NA values and replace them with 0
  dist_matrices[[i]][is.na(dist_matrices[[i]])] <- 0
}

dist_matrices

# now for permanova but accounting for datasets that are missing samples from one level of Sample_type
rel_perm_results <- vector("list", 113)  # Pre-allocate a list of 113 elements
rel_perm_pvals <- vector("list", 113)
metadata_data <- vector("list", 113)

for (i in 1:113) {
  
  ps_obj <- updated_ps_list[[i]]
  metadata_data[[i]] <- as(sample_data(ps_obj), "data.frame")
  # Check if there is only one unique level of Sample_type
  if (length(unique(metadata_data[[i]]$Sample_type)) == 1) {
    message(paste("Only one level of Sample_type for object", i))
    rel_perm_results[[i]] <- NA
    rel_perm_pvals[[i]] <- NA
  } else {
    # Run permanova analysis
    rel_perm_results[[i]] <- adonis2(dist_matrices[[i]] ~ Sample_type, data = metadata_data[[i]], permutations = 9999)
    
    # Get p-value from the permanova result and store it in the list
    rel_perm_pvals[[i]] <- rel_perm_results[[i]]$'Pr(>F)'[1]
    
  }
}

rel_perm_results
rel_perm_pvals

names(rel_perm_pvals) <- names(TotalReads_list)


# multiple testing correstion
# Unlist the p-value lists
rel_perm_pvals_unlisted <- unlist(rel_perm_pvals)

# FDR correction for obs_ttest_pval_list
rel_perm_pvals_fdr <- p.adjust(rel_perm_pvals_unlisted , method = "fdr")
rel_perm_pvals_fdr

write.csv(rel_perm_pvals_fdr, "RA_permanova_pvalues113.csv", row.names = TRUE)




#### Table S2 ####

#### get relative abundance means of the phyla #### 

phyla_abund_list=list(1:113)
phyla_glom_list=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=updated_ps_list[[i]] #Select the ith dataset in your list
  
  phyla_glom_list[[i]] <- tax_glom(df, taxrank = 'Rank2', NArm=FALSE)
  
  phyla_abund_list[[i]] = psmelt(phyla_glom_list[[i]])
}

phyla_abund_list
phyla_glom_list



# get the means for every Phylum for sample type
library(plyr)


phyla_means_list=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=phyla_abund_list[[i]] #Select the ith dataset in your list
  
  #first I need to sum the Unclassified and low abundances
  sumPhylum.df<- ddply(df, c("Sample", "Rank2", "Sample_type"), summarise,
                       sum.Abund=sum(Abundance))
  #now get the means
  phyla_means_list[[i]] <- ddply(sumPhylum.df, c("Sample_type", "Rank2"), summarise,
                         mean.Abund=mean(sum.Abund))
}

phyla_means_list



names(phyla_means_list)=c("K00211","K00608",
                          "K00998","K01085","K01352","K01355","K01819","K02383","K03181",
                          "K03202","K03397","K04338","K04643","K05372","K05660","K06080",
                          "K07345","K07786","K07991","K08156","K08566","K10094","K10926",
                          "K10975","K11009","K11017","K11734","K11743","K11889","K11909",
                          "K11934","K11964","K12059","K12069","K12083","K12228","K12281",
                          "K12284","K12285","K12286","K12455","K12688","K12820","K13448",
                          "K13454","K13585","K13586","K13831","K13964","K14368","K14627",
                          "K14781","K15968","K16215","K16347","K16348","K16397","K16403",
                          "K16448","K16552","K16553","K16696","K18373","K18380","K18642",
                          "K18793","K19033","K19060","K19101","K19213","K19611","K19885",
                          "K20086","K20088","K20090","K20267","K20268","K20272","K20275",
                          "K20326","K20555","K20961","K21212","K21239","K21280","K21346",
                          "K21825","K22335","K20966","K19216","K02399","K10927","K12280"
)

# subset the unique Sample_type levels
phyla_means_list_Heritage=list(1:113)
phyla_means_list_Modern=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  #df=alpha_data_list[[i]] #Select the ith dataset in your list
  phyla_means_list_Heritage[[i]] <- subset(phyla_means_list[[i]], Sample_type == "Heritage")
  phyla_means_list_Modern[[i]] <- subset(phyla_means_list[[i]], Sample_type == "Modern")
}

names(phyla_means_list_Heritage) <- names(phyla_means_list)
names(phyla_means_list_Modern) <- names(phyla_means_list)

# Get unique Rank2 levels across both lists
all_rank2_levels_Heritage <- sort(unique(unlist(lapply(phyla_means_list_Heritage, function(x) x[["Rank2"]]))))
all_rank2_levels_Modern <- sort(unique(unlist(lapply(phyla_means_list_Modern, function(x) x[["Rank2"]]))))

# Create a data frame to hold the final result
final_df_Heritage <- data.frame(matrix(nrow = 113, ncol = length(all_rank2_levels_Heritage))) 
final_df_Modern <- data.frame(matrix(nrow = 113, ncol = length(all_rank2_levels_Modern))) 
# assign column names
colnames(final_df_Heritage) = all_rank2_levels_Heritage
colnames(final_df_Modern) = all_rank2_levels_Modern
# display
print(final_df_Heritage)
print(final_df_Modern)

final_df_Heritage$KEGG_ID <- names(phyla_means_list_Heritage)
final_df_Modern$KEGG_ID <- names(phyla_means_list_Modern)

# Populate the data frame with data from each list
for (i in 1:113) {
  df_Heritage <- phyla_means_list_Heritage[[i]]
  df_Modern <- phyla_means_list_Modern[[i]]
  
  final_df_Heritage[i, "KEGG_ID"] <- names(phyla_means_list_Heritage)[i]
  final_df_Modern[i, "KEGG_ID"] <- names(phyla_means_list_Modern)[i]
  
  # Populate data from the Heritage_df
  for (rank2 in all_rank2_levels_Heritage) {
    if (rank2 %in% df_Heritage[["Rank2"]]) {
      final_df_Heritage[i, rank2] <- df_Heritage[["mean.Abund"]][df_Heritage[["Rank2"]] == rank2]
    }
  }
  
  # Populate data from the Modern_df
  for (rank2 in all_rank2_levels_Modern) {
    if (rank2 %in% df_Modern[["Rank2"]]) {
      final_df_Modern[i, rank2] <- df_Modern[["mean.Abund"]][df_Modern[["Rank2"]] == rank2]
    }
  }
}

final_df_Modern
final_df_Heritage

write.csv(final_df_Modern, "Phyla_means_113_df_Modern.csv")
write.csv(final_df_Heritage, "Phyla_means_113_df_Heritage.csv")


# same as individual files.... 
# library(rio)
# #export_list(phyla_means_list, file = "%s.csv")
# #added file path so they were kept separate. 
# export_list(phyla_means_list, file = file.path('phyla_means',"%s.csv"))



#### get means of the class #### 

#"K05660"
head(tax_table(updated_ps_list[["K05660"]]))
# there is no Rank3 for some ps objects in the list. 

class_abund_list=list(1:113)

for (i in 1:113) {
  # Try to process each dataset and skip if there's an error
  tryCatch({
    df = updated_ps_list[[i]]
    
    # Check if Rank3 exists in the tax_table
    if (!("Rank3" %in% colnames(tax_table(df)))) {
      cat("Skipping dataset", i, "- Rank3 column does not exist in tax_table\n")
      next  # Skip to the next iteration
    }
    
    # If we made it here, try to glom and melt
    glom <- tax_glom(df, taxrank = 'Rank3', NArm = FALSE)
    class_abund_list[[i]] = psmelt(glom)
    cat("Successfully processed dataset", i, "\n")
    
  }, error = function(e) {
    # This will catch any errors that occur during processing
    cat("Error processing dataset", i, ":", conditionMessage(e), "\n")
  })
}

class_abund_list





# get the means for every class for sample type
library(plyr)


class_means_list=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=class_abund_list[[i]] #Select the ith dataset in your list
  
  #first I need to sum the Unclassified and low abundances
  sumClass.df<- ddply(df, c("Sample", "Rank3", "Sample_type"), summarise,
                       sum.Abund=sum(Abundance))
  #now get the means
  class_means_list[[i]] <- ddply(sumClass.df, c("Sample_type", "Rank3"), summarise,
                                 mean.Abund=mean(sum.Abund),
                                 sd.Abund=sd(sum.Abund))
}

class_means_list



names(class_means_list)=names(ps_list)

# subset the unique Sample_type levels
class_means_list_Heritage=list(1:113)
class_means_list_Modern=list(1:113)

for (i in 1:113) {
  # Check if the current element exists and isn't NULL
  if (!is.null(class_means_list[[i]])) {
    # Check if the dataframe has the Sample_type column
    if ("Sample_type" %in% colnames(class_means_list[[i]])) {
      # Subset the data
      class_means_list_Heritage[[i]] <- subset(class_means_list[[i]], Sample_type == "Heritage")
      class_means_list_Modern[[i]] <- subset(class_means_list[[i]], Sample_type == "Modern")
      cat("Successfully processed dataset", i, "\n")
    } else {
      cat("Dataset", i, "doesn't have 'Sample_type' column\n")
    }
  } else {
    cat("Dataset", i, "is NULL or missing\n")
  }
}

names(class_means_list_Heritage) <- names(class_means_list)
names(class_means_list_Modern) <- names(class_means_list)

# Get unique Rank3 levels across both lists
all_rank3_levels_Heritage <- sort(unique(unlist(lapply(class_means_list_Heritage, function(x) x[["Rank3"]]))))
all_rank3_levels_Modern <- sort(unique(unlist(lapply(class_means_list_Modern, function(x) x[["Rank3"]]))))

# Create a data frame to hold the final result
final_df_Heritage <- data.frame(matrix(nrow = 113, ncol = length(all_rank3_levels_Heritage))) 
final_df_Modern <- data.frame(matrix(nrow = 113, ncol = length(all_rank3_levels_Modern))) 
# assign column names
colnames(final_df_Heritage) = all_rank3_levels_Heritage
colnames(final_df_Modern) = all_rank3_levels_Modern
# display
print(final_df_Heritage)
print(final_df_Modern)

final_df_Heritage$KEGG_ID <- names(class_means_list_Heritage)
final_df_Modern$KEGG_ID <- names(class_means_list_Modern)

# Populate the data frame with data from each list
for (i in 1:113) {
  df_Heritage <- class_means_list_Heritage[[i]]
  df_Modern <- class_means_list_Modern[[i]]
  
  final_df_Heritage[i, "KEGG_ID"] <- names(class_means_list_Heritage)[i]
  final_df_Modern[i, "KEGG_ID"] <- names(class_means_list_Modern)[i]
  
  # Populate data from the Heritage_df
  for (Rank3 in all_rank3_levels_Heritage) {
    if (Rank3 %in% df_Heritage[["Rank3"]]) {
      final_df_Heritage[i, Rank3] <- df_Heritage[["mean.Abund"]][df_Heritage[["Rank3"]] == Rank3]
    }
  }
  
  # Populate data from the Modern_df
  for (Rank3 in all_rank3_levels_Modern) {
    if (Rank3 %in% df_Modern[["Rank3"]]) {
      final_df_Modern[i, Rank3] <- df_Modern[["mean.Abund"]][df_Modern[["Rank3"]] == Rank3]
    }
  }
}

final_df_Modern
final_df_Heritage

write.csv(final_df_Modern, "class_means_df_Modern.csv")
write.csv(final_df_Heritage, "class_means_df_Heritage.csv")



# PERMANOVA on relative abundances PHYLA ####
# library(vegan)
#library("plyr")
sample_data(phyla_glom_list[[65]])
otu_table(phyla_glom_list[[57]])
# first make distance matrices
phyla_dist_matrices <- vector("list", 113)

for (i in 1:113) {
  
  ps_obj <- phyla_glom_list[[i]]
  
  phyla_dist_matrices[[i]] <- phyloseq::distance(ps_obj, method="bray")
  # Check for NA values and replace them with 0
  phyla_dist_matrices[[i]][is.na(phyla_dist_matrices[[i]])] <- 0
}

phyla_dist_matrices



# now for permanova but accounting for datasets that are missing samples from one level of Sample_type
phyla_rel_perm_results <- vector("list", 113)  # Pre-allocate a list of 113 elements
phyla_rel_perm_pvals <- vector("list", 113)
phyla_metadata_data <- vector("list", 113)

for (i in 1:113) {
  
  ps_obj <- phyla_glom_list[[i]]
  phyla_metadata_data[[i]] <- as(sample_data(ps_obj), "data.frame")
  # Check if there is only one unique level of Sample_type
  if (length(unique(phyla_metadata_data[[i]]$Sample_type)) == 1) {
    message(paste("Only one level of Sample_type for object", i))
    phyla_rel_perm_results[[i]] <- NA
    phyla_rel_perm_pvals[[i]] <- NA
  } else {
    # Run permanova analysis
    phyla_rel_perm_results[[i]] <- adonis2(phyla_dist_matrices[[i]] ~ Sample_type, data = phyla_metadata_data[[i]], permutations = 9999)
    
    # Get p-value from the permanova result and store it in the list
    phyla_rel_perm_pvals[[i]] <- phyla_rel_perm_results[[i]]$'Pr(>F)'[1]
    
  }
}

phyla_rel_perm_results
phyla_rel_perm_pvals

names(phyla_rel_perm_pvals) <- names(TotalReads_list)


# multiple testing correction
# Unlist the p-value lists
phyla_rel_perm_pvals_unlisted <- unlist(phyla_rel_perm_pvals)

# FDR correction for obs_ttest_pval_list
phyla_rel_perm_pvals_fdr <- p.adjust(phyla_rel_perm_pvals_unlisted , method = "fdr")
phyla_rel_perm_pvals_fdr

write.csv(phyla_rel_perm_pvals_fdr, "RA_phyla_permanova_pvalues113.csv", row.names = TRUE)





# Table 3 ####
#### Alpha diversity ####

alpha_data_list=list(1:113)
alpha_means_list=list(1:113)
alpha_sd_list=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  df=ps_list_norm_int[[i]] #go back to the normalised integer data
  alpha_div_norm <- estimate_richness(df, measures=c("Observed", "Shannon"))
  alpha_data_list[[i]] <- cbind(as(alpha_div_norm, "data.frame"), as(sample_data(df)[rownames(alpha_div_norm), ], "matrix"))
  alpha_means_list[[i]] <- ddply(alpha_data_list[[i]], .(Sample_type), colwise(mean))
  alpha_sd_list[[i]] <- ddply(alpha_data_list[[i]], .(Sample_type), colwise(sd))
}

alpha_data_list
alpha_means_list
alpha_sd_list

names(alpha_data_list)=names(ps_list)
names(alpha_means_list)=names(ps_list)
names(alpha_sd_list)=names(ps_list)

# Combine the lists into dataframes
alpha_data_df <- do.call(rbind, alpha_data_list)
alpha_means_df <- do.call(rbind, alpha_means_list)
alpha_sd_df <- do.call(rbind, alpha_sd_list)

# Extract the base names (everything before the period)
kegg_ids <- gsub("\\.[0-9]+$", "", rownames(alpha_means_df))
# Add the KEGG_ID column to your dataframe
alpha_means_df$KEGG_ID <- kegg_ids

# Extract the base names (everything before the period)
kegg_ids <- gsub("\\.[0-9]+$", "", rownames(alpha_sd_df))
# Add the KEGG_ID column to your dataframe
alpha_sd_df$KEGG_ID <- kegg_ids

# If you want to remove the original rownames and use 1:nrow instead
rownames(alpha_means_df) <- NULL
rownames(alpha_sd_df) <- NULL


# Transform the merged dataframe to wide format
alpha_means_df_wide <- alpha_means_df %>%
  pivot_wider(id_cols = KEGG_ID,
              names_from = Sample_type,
              values_from = c(Observed, Shannon),
              names_glue = "{Sample_type}.{.value}")

alpha_sd_df_wide <- alpha_sd_df %>%
  pivot_wider(id_cols = KEGG_ID,
              names_from = Sample_type,
              values_from = c(Observed, Shannon),
              names_glue = "{Sample_type}.{.value}")

# # Save the combined dataframes to separate CSV files
# write.csv(alpha_data_df, "alpha_data_list113.csv", row.names = TRUE)
# write.csv(alpha_means_df, "alpha_means_list113.csv", row.names = TRUE)
# write.csv(alpha_sd_df, "alpha_sd_list113.csv", row.names = TRUE)

write.csv(alpha_means_df_wide, "alpha_means_list113_wide.csv", row.names = TRUE)
write.csv(alpha_sd_df_wide, "alpha_sd_list113_wide.csv", row.names = TRUE)




# run some tests on the alpha diversity ####

# 1) remove bulk soil
# 2) check assumptions
# 3) simple t-test

# 1) remove bulk soil
alpha_data_list_plants=list(1:113)

for (i in c(1:113)) #The loop code will be repeated x times, based on the values that “I” refers to
{
  #df=alpha_data_list[[i]] #Select the ith dataset in your list
  alpha_data_list_plants[[i]] <- subset(alpha_data_list[[i]], Sample_type != "Bulk_Soil")
}

alpha_data_list_plants


# 2) check assumptions
obs_norm_list=list(1:113)
obs_var_list=list(1:113)
obs_histo_list=list(1:113)
shan_norm_list=list(1:113)
shan_var_list=list(1:113)
shan_histo_list=list(1:113)

for (i in c(1:113)) {
  # Use tryCatch block to handle errors
  tryCatch({
    df <- alpha_data_list_plants[[i]]
    
    obs_norm_list[[i]] <- shapiro.test(df$Observed)
    obs_var_list[[i]] <- leveneTest(Observed ~ Sample_type, data = df)
    shan_norm_list[[i]] <- shapiro.test(df$Shannon)
    shan_var_list[[i]] <- leveneTest(Shannon ~ Sample_type, data = df)
    
    shan_histo_list[[i]] <- ggplot(df, aes(Shannon)) + geom_histogram() + ggtitle(names(alpha_data_list_plants)[i])
    obs_histo_list[[i]] <- ggplot(df, aes(Observed)) + geom_histogram() + ggtitle(names(alpha_data_list_plants)[i])
  }, error = function(cond) {
    # Return dataframe with NAs
    message("Error occurred in dataframe ", names(alpha_data_list_plants)[i])
    obs_norm_list[[i]] <- NA
    obs_var_list[[i]] <- NA
    shan_norm_list[[i]] <- NA
    shan_var_list[[i]] <- NA
    shan_histo_list[[i]] <- NA
    obs_histo_list[[i]] <- NA
  })
}


ggarrange(plotlist=shan_histo_list) #Plots all graphs in the list
ggarrange(plotlist=obs_histo_list)

names(obs_norm_list)=names(ps_list)
names(obs_var_list)=names(ps_list)
names(shan_norm_list)=names(ps_list)
names(shan_var_list)=names(ps_list)


# results
obs_norm_list # Significant = K02532 (only 0.037) - K00608 (0.0001),K02383 (2.243e-08), K03397 (0.0027), K05660 (4.538e-06), 
                            # K06080 (0.0015),K07991 (0.0001), K10975 (0.01), K11009 (0.016), K11743 (0.001), K12059 (0.016),
                            # K12083 (7.533e-05), K19216 (0.046), K12280 (0.018)
obs_var_list # Significant = K02532 (0.018) - K00608 (0.002), K05372 (0.048), K07991 (0.002), K11743 (0.013), K12083 (0.008), 
                            # K12228 (0.036), K12281 (0.013), K16347 (0.046), K18380 (0.018), K20088 (0.0004), K21346 (0.017), K12280 (0.018)
# significant in both variance and normality deviations = K02532, K00608, K07991, K11743, K12083

shan_norm_list # Significant = K02532 (0.005) - K00608 (0.0001), K01819 (0.037), K02383 (2.243e-08), K03397 (4.17e-06), K05372 (0.017), 
                            # K05660 (9.922e-07), K06080 (1.819e-07), K07991 (2.836e-06), K08156 (0.041), K08566 (0.001), K10926 (0.031), 
                            # K10975 (3.894e-06), K11017 (0.0036), K11743 (0.005), K11889 (0.0006), K11964 (0.032), K12083( 0.048), K12228 (0.001),
                            # K12820 (0.0002), K13831 (0.0006), K15968 (0.0045), K16347 (0.014), K16397 (0.028), K16448 (0.007), 
shan_var_list # Significant = K18651 (0.022), K08325 (0.013), K05874 (0.043)

# I think this is ok to go ahead with the t-tests. They are pretty robust. 

# 3) simple t-test

obs_ttest_list=list(1:113)
obs_ttest_pval_list=list(1:113)
shan_ttest_list=list(1:113)
shan_ttest_pval_list=list(1:113)

for (i in c(1:113)) {
  tryCatch({
    df <- alpha_data_list_plants[[i]]
    
    obs_ttest_list[[i]] <- t.test(Observed ~ Sample_type, data = df)
    obs_ttest_pval_list[[i]] <- obs_ttest_list[[i]]$p.value
    
    shan_ttest_list[[i]] <- t.test(Shannon ~ Sample_type, data = df)
    shan_ttest_pval_list[[i]] <- shan_ttest_list[[i]]$p.value
  }, error = function(cond) {
    # Return dataframe with NAs
    message("Error occurred in dataframe ", names(alpha_data_list_plants)[i])
    obs_ttest_list[[i]] <- NA
    shan_ttest_list[[i]] <- NA
    obs_ttest_pval_list[[i]] <- NA
    shan_ttest_pval_list[[i]] <- NA
  })
}

names(obs_ttest_pval_list)=names(ps_list)
names(shan_ttest_pval_list)=names(ps_list)
obs_ttest_pval_list
shan_ttest_pval_list

# Unlist the p-value lists
obs_ttest_pval_list_unlisted <- unlist(obs_ttest_pval_list)
shan_ttest_pval_list_unlisted <- unlist(shan_ttest_pval_list)

# FDR correction for obs_ttest_pval_list
obs_ttest_pval_list_fdr <- p.adjust(obs_ttest_pval_list_unlisted , method = "fdr")
obs_ttest_pval_list_fdr
# FDR correction for shan_ttest_pval_list
shan_ttest_pval_list_fdr <- p.adjust(shan_ttest_pval_list_unlisted, method = "fdr")
shan_ttest_pval_list_fdr


write.csv(obs_ttest_pval_list_fdr, "Observed_pvalues113.csv", row.names = TRUE)
write.csv(shan_ttest_pval_list_fdr, "Shannon_pvalues113.csv", row.names = TRUE)


#### END ####