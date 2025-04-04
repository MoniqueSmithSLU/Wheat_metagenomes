##############################################################################
#
# Reference to ARTICLE
# Monique Smith, Sweden's Agricultural University (SLU)
# monique.smith@slu.se
# revisions 1st April 2025

# Script to support statistical analysis and figures for shotgun metagenomic data aligned to either function (KEGG database) or taxonomy (annotree)

# additional files needed: 
# annomeg_KEGG2.biom
# annomeg_Taxa.biom
# Metadata_merged2.txt
# tern_e 2.R
# sig_and_core_log2fold_HM_KEGG_annotree_function.csv
# tt_mixedRanks_annomeg_sorted_ranks.txt


# Data has been extracted from MEGAN6 from the merged samples (large samples were split to meganise and then merged)
# Samples were aligned using AnnoTree so only contain Bacteria and Archaea 
# In Megan - I uncollapse all nodes, select all leaves and go file->export->biom

##############################################################################

# install specialty packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocGenerics",force = TRUE)
BiocManager::install("DESeq2",force = TRUE)

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")

library(remotes)
remotes::install_github("gauravsk/ranacapa")

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#### Load packages and data ####

library(phyloseq)
library(ggplot2)
library(DESeq2)
library("data.table")
library(microbiome)
library(psadd)
library("pheatmap")
library("RColorBrewer")
library(devtools)
library("ashr")
library(vegan)
library(plyr)
library(car)
library(ranacapa)
library(remotes)
library(emmeans)
library(grid)
library(tidyverse)
library(dplyr)
library(metagMisc)
library(ggpubr)
library(ranacapa)
library(pairwiseAdonis)
library(ggpubr)
library(patchwork)
library(ggh4x)

# Import the biom table
biom.KEGG = import_biom("annomeg_KEGG2.biom") 
biom.TAXA = import_biom("annomeg_Taxa.biom")

# meta data
md <- read.delim("Metadata_merged2.txt", header = TRUE, sep = "\t", dec = ".")
# add a new factor so figures don't refer to tall and semi-dwarf as per reviewers comments
md$Sample_type <- "Bulk_Soil"
md[md$Height %in% c('Semi_dwarf'),]$Sample_type <- "Modern"
md[md$Height %in% c('Tall'),]$Sample_type <- "Heritage"
# meta data needs to have the samples as row names
rownames(md) <- md[,1]
md[,1] <- NULL
md
# this is needed to tell R that md is sample data
md2 <- sample_data(md)

# now merge into phyloseq object
psKEGG <- merge_phyloseq(biom.KEGG, md2)
psTAXA <- merge_phyloseq(biom.TAXA, md2)

# The taxa table from Megan has some mixed classification. There are taxa names in the wrong column. 
# This needs to be fixed in order to carry out some of the analysis, certainly for plotting. 
# To do this we made a script to move everything into the right columns
# export
#write.csv(tax_table(alltaxa), file = "tt_mixedRanks_annomeg.csv",row.names = TRUE)
# Use script reformat_soil_prokaryote_taxa_df.pl
# Now upload the sorted file
newtaxa <- read.delim("tt_mixedRanks_annomeg_sorted_ranks.txt", header = TRUE, sep = "\t", dec = ".")
newtaxa.df <- as.data.frame(unclass(newtaxa), stringsAsFactors = TRUE)

# set the IDs as row names
newtaxa <- newtaxa %>% 
  tibble::column_to_rownames("taxon_ID")

#Transform into matrixes otu and tax tables (sample table can be left as data frame)
newtaxa <- as.matrix(newtaxa)
#Transform to phyloseq object
TAX <- tax_table(newtaxa)
# combine new phyloseq
psTAXA <- phyloseq(otu_table(psTAXA), TAX, md2)



         
# Prefiltering ####
# removing unclassified reads

psKEGG <- subset_taxa(psKEGG, Rank3 != is.na(Rank3))
print(psKEGG)
#10517 genes, removed 8 

psTAXA <- subset_taxa(psTAXA, Phylum != is.na(Phylum))
print(psTAXA)
#23947 taxa, removed 400


# set the minimum reads to 10 for each gene
psKEGG = prune_taxa(taxa_sums(psKEGG) > 10, psKEGG)
psTAXA = prune_taxa(taxa_sums(psTAXA) > 10, psTAXA)



#### Data transformations ####

# relative abundances using the microbiome package
#compositional' abundances are returned as relative abundances in [0, 1],(convert to percentages by multiplying with a factor of 100).
psKEGG.rel <- microbiome::transform(psKEGG, "compositional")
psTAXA.rel <- microbiome::transform(psTAXA, "compositional")

# rarefy to the minimum number of reads
# find the minimum reads
microbiome::summarize_phyloseq(psKEGG)
microbiome::summarize_phyloseq(psTAXA)
#Min. number of reads = 118292681
# rarefy without replacement to minimum depth
psKEGG.rar <- rarefy_even_depth(psKEGG, rngseed = 118292681)
psTAXA.rar <- rarefy_even_depth(psTAXA, rngseed = 159120960)

#### Fig. 1a ####
#Alpha diversity plots 
# create data for plots
alpha_div_rar <- estimate_richness(psTAXA.rar)
alpha_div_rar_taxa <- cbind(as(alpha_div_rar, "data.frame"), as(sample_data(psTAXA.rar)[rownames(alpha_div_rar), ], "matrix"))
ddply(alpha_div_rar_taxa, .(Sample_type), colwise(sd))
alpha_div_rar_taxa_trim <- alpha_div_rar_taxa[,c("Observed","Chao1","Shannon","Variety","Sample_type")]
alpha_div_rar_taxa_long <- melt(setDT(alpha_div_rar_taxa_trim), id.vars = c("Variety","Sample_type"), variable.name = "alpha_diversity")
#write.table(alpha_div_rar_taxa_long, file = 'dataFig1a.txt', sep = ',')

Fig1a <- ggplot(alpha_div_rar_taxa_long, aes(x=Sample_type, y=value)) + facet_wrap(~alpha_diversity,scales = "free_y") +
  geom_boxplot() +
  geom_jitter(aes(color = Variety, shape = Variety), size = 3, position=position_jitter(0.2), show.legend = FALSE) +
  scale_color_manual(#name="Sample",
    breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
    values=c("grey","#B12122","#F58c8c", "#4069E1" ,"#ACD8E5"),
    #labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")
  ) +
  scale_shape_manual(#name="Sample",
    breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
    values = c(19,17,17,15,15),
    #labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")
  ) +
  ylab(paste("Alpha diversity value")) +
  xlab(paste("Sample type")) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("Bulk_Soil" = "BS", "Modern" = "M",
                            "Heritage" = "H")) +
  ggtitle("Taxonomy") +
  theme_bw()
Fig1a


# supporting statistics for Fig1a ####
# 1) Observed genes
# test for normality- Shapiro-Wilk statistic
shapiro.test(alpha_div_rar_taxa$Observed) # not significant
hist(alpha_div_rar_taxa$Observed)

#homogeneity of variances - levene's test
#library(car)
leveneTest(Observed ~ Sample_type, data = alpha_div_rar_taxa) # not significant

#Analysis of Variance
# library(agricolae)
boxplot(Observed ~ Sample_type, data = alpha_div_rar_taxa)
#anova <- oneway.test(Observed ~ Sample_type, data = alpha_div_rar_taxa)
mod <- lm(Observed ~ Sample_type, data = alpha_div_rar_taxa)
summary(mod)
anova <- Anova(mod, type="III") # accounts for unbalanced design
anova
# contrasts
pairs(emmeans(mod, "Sample_type"))

# 2) Chao1
# test for normality- Shapiro-Wilk statistic
shapiro.test(alpha_div_rar_taxa$Chao1) # not significant
hist(alpha_div_rar_taxa$Chao1)

#homogeneity of variances - levene's test
#library(car)
leveneTest(Chao1 ~ Sample_type, data = alpha_div_rar_taxa) # not significant

#Analysis of Variance
boxplot(Chao1 ~ Sample_type, data = alpha_div_rar_taxa)
mod <- lm(Chao1 ~ Sample_type, data = alpha_div_rar_taxa)
summary(mod)
anova <- Anova(mod, type="III") # accounts for unbalanced design
anova
# contrasts
pairs(emmeans(mod, "Sample_type"))


# 3) Shannon
# test for normality- Shapiro-Wilk statistic
shapiro.test(alpha_div_rar_taxa$Shannon) # 0.04 significant
hist(alpha_div_rar_taxa$Shannon)

#homogeneity of variances - levene's test
#library(car)
leveneTest(Shannon ~ Sample_type, data = alpha_div_rar_taxa) # not significant

#Analysis of Variance
boxplot(Shannon ~ Sample_type, data = alpha_div_rar_taxa)
mod <- lm(Shannon ~ Sample_type, data = alpha_div_rar_taxa)
summary(mod)
anova <- Anova(mod, type="III") # accounts for unbalanced design
anova
# contrasts
#library(emmeans)
pairs(emmeans(mod, "Sample_type"))



#### Fig. 1c #####
#Alpha diversity plots 
# create data for plots
alpha_div_rar <- estimate_richness(psKEGG.rar)
alpha_div_rar_kegg <- cbind(as(alpha_div_rar, "data.frame"), as(sample_data(psKEGG.rar)[rownames(alpha_div_rar), ], "matrix"))
ddply(alpha_div_rar_kegg, .(Sample_type), colwise(sd))

alpha_div_rar_kegg_trim <- alpha_div_rar_kegg[,c("Observed","Chao1","Shannon","Variety","Sample_type")]
alpha_div_rar_kegg_long <- melt(setDT(alpha_div_rar_kegg_trim), id.vars = c("Variety","Sample_type"), variable.name = "alpha_diversity")

#write.table(alpha_div_rar_kegg_long, file = 'dataFig1c.txt', sep = ',')

Fig1c <- ggplot(alpha_div_rar_kegg_long, aes(x=Sample_type, y=value)) + facet_wrap(~alpha_diversity,scales = "free_y") +
  geom_boxplot() +
  geom_jitter(aes(color = Variety, shape = Variety), size = 3, position=position_jitter(0.2), show.legend = FALSE) +
  scale_color_manual(#name="Sample",
                     breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
                     values=c("grey","#B12122","#F58c8c", "#4069E1" ,"#ACD8E5"),
                     #labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")
                     ) +
  scale_shape_manual(#name="Sample",
                     breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
                     values = c(19,17,17,15,15),
                     #labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")
                     ) +
  ylab(paste("Alpha diversity value")) +
  xlab(paste("Sample type")) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("Bulk_Soil" = "BS", "Modern" = "M",
                              "Heritage" = "H")) +
  ggtitle("Functional genes") +
  theme_bw()
Fig1c



# supporting statistics for Fig1c ####
# 1) Observed genes
# test for normality- Shapiro-Wilk statistic
shapiro.test(alpha_div_rar_kegg$Observed) # not significant
hist(alpha_div_rar_kegg$Observed)

#homogeneity of variances - levene's test
#library(car)
leveneTest(Observed ~ Sample_type, data = alpha_div_rar_kegg) # significant

#Analysis of Variance
# library(agricolae)
boxplot(Observed ~ Sample_type, data = alpha_div_rar_kegg)
anova <- oneway.test(Observed ~ Sample_type, data = alpha_div_rar_kegg)
mod <- lm(Observed ~ Sample_type, data = alpha_div_rar_kegg)
summary(mod)
anova <- Anova(mod, type="III") # accounts for unbalanced and unequal variance
anova
# contrasts
#library(emmeans)
pairs(emmeans(mod, "Sample_type"))



# 2) Chao1
# test for normality- Shapiro-Wilk statistic
shapiro.test(alpha_div_rar_kegg$Chao1) # not significant
hist(alpha_div_rar_kegg$Chao1)

#homogeneity of variances - levene's test
#library(car)
leveneTest(Chao1 ~ Sample_type, data = alpha_div_rar_kegg) # not significant

#Analysis of Variance
boxplot(Chao1 ~ Sample_type, data = alpha_div_rar_kegg)
mod <- lm(Chao1 ~ Sample_type, data = alpha_div_rar_kegg)
summary(mod)
anova <- Anova(mod, type="III") # accounts for unbalanced design
anova
# contrasts
#library(emmeans)
pairs(emmeans(mod, "Sample_type"))


# 3) Shannon
# test for normality- Shapiro-Wilk statistic
shapiro.test(alpha_div_rar_kegg$Shannon) # not significant
hist(alpha_div_rar_kegg$Shannon)

#homogeneity of variances - levene's test
#library(car)
leveneTest(Shannon ~ Sample_type, data = alpha_div_rar_kegg) # not significant

#Analysis of Variance
boxplot(Shannon ~ Sample_type, data = alpha_div_rar_kegg)
mod <- lm(Shannon ~ Sample_type, data = alpha_div_rar_kegg)
summary(mod)
anova <- Anova(mod, type="III") # accounts for unbalanced 
anova
# contrasts
#library(emmeans)
pairs(emmeans(mod, "Sample_type"))

# Fig. 1b ####
# FIRST get the data for ggplot from phyloseq object
TAXA.ord <- ordinate(psTAXA.rar, "PCoA", "bray")
# Ordination object plus all metadata, note: we use justDF=T. This will not return a plot but a data.frame
ordip.taxa <- plot_ordination(psTAXA.rar, TAXA.ord, justDF = T)
# Get axis 1 and 2 variation
evals1 <- round(TAXA.ord$values$Eigenvalues[1] / sum(TAXA.ord$values$Eigenvalues) * 100, 2)
evals2 <- round(TAXA.ord$values$Eigenvalues[2] / sum(TAXA.ord$values$Eigenvalues) * 100, 2)
# theme_set(theme_bw(14))

#write.table(ordip.taxa, file = 'dataFig1b.txt', sep = ',')

# PLOT
# use the ordip
# first step get blank plot
p.taxa <- ggplot(ordip.taxa, aes(x = Axis.1, y = Axis.2))

# NEED TO RUN A FEW LINES ABOVE TO MAKE SURE YOURE GRAPHING GENES NOT TAXA
Fig1b <- p.taxa +
  geom_hline(yintercept=0, colour = "grey20", linewidth = 0.4) +
  geom_vline(xintercept=0, colour = "grey20", linewidth = 0.4) +
  # now add a layer of points 
  geom_point(aes(color = Variety, shape = Variety), size = 3) +
  scale_color_manual(name="Sample",
                     breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
                     values=c("grey","#B12122","#F58c8c", "#4069E1" ,"#ACD8E5"),
                     labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")) +
  scale_shape_manual(name="Sample",
                     breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
                     values = c(19,17,17,15,15),
                     labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")) +
  xlab(paste("PCoA 1 (", evals1, "%)", sep = "")) +
  ylab(paste("PCoA 2 (", evals2, "%)", sep = "")) +
  #ggtitle("a) Taxonomy") +
  theme_bw() +
  theme(plot.margin = unit(c(22,5.5,5.5,5.5), 'points'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) # + coord_fixed() # optional to demonstrate that the majority of variation is from PCoA 1

# Print figure
print(Fig1b)






# Fig. 1d ####
# prepare the data
# FIRST get the data for ggplot from phyloseq object
KEGG.ord <- ordinate(psKEGG.rar, "PCoA", "bray")
# Ordination object plus all metadata, note: we use justDF=T. This will not return a plot but a data.frame
ordip <- plot_ordination(psKEGG.rar, KEGG.ord, justDF = T)
# Get axis 1 and 2 variation
evals1 <- round(KEGG.ord$values$Eigenvalues[1] / sum(KEGG.ord$values$Eigenvalues) * 100, 2)
evals2 <- round(KEGG.ord$values$Eigenvalues[2] / sum(KEGG.ord$values$Eigenvalues) * 100, 2)
# theme_set(theme_bw(14))

#write.table(ordip, file = 'dataFig1d.txt', sep = ',')

# PLOT
# first step get blank plot
p <- ggplot(ordip, aes(x = Axis.1, y = Axis.2))

# NEED TO RUN A FEW LINES ABOVE TO MAKE SURE YOURE GRAPHING GENES NOT TAXA
Fig1d <- p +
  geom_hline(yintercept=0, colour = "grey20", linewidth = 0.4) +
  geom_vline(xintercept=0, colour = "grey20", linewidth = 0.4) +
  # now add a layer of points 
  geom_point(aes(color = Variety, shape = Variety), size = 3) +
  scale_color_manual(name="Sample",
    breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
                     values=c("grey","#B12122","#F58c8c", "#4069E1" ,"#ACD8E5"),
                     labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")) +
  scale_shape_manual(name="Sample",
    breaks = c("Bulk_Soil", "Chidham", "Red_Llamas", "Gallant", "Malacca"),
                     values = c(19,17,17,15,15),
                     labels = c("Bulk soil", "Chidham White Chaff", "Red Lammas","Gallant", "Malacca")) +
  xlab(paste("PCoA 1 (", evals1, "%)", sep = "")) +
  ylab(paste("PCoA 2 (", evals2, "%)", sep = "")) +
  #ggtitle("b) Functional genes") +
  theme_bw() +
  theme(plot.margin = unit(c(22,5.5,5.5,5.5), 'points'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) # + coord_fixed() # optional to demonstrate that the majority of variation is from PCoA 1

# Print figure
print(Fig1d)

# Use this to check the colours of elements in the black and white theme
#calc_element("panel.border", theme_bw())

# Fig. 1 combined ####
# combining Figure one is in script 
# Fig1c, Fig1d
Fig1.combined2 <- ggarrange(Fig1a, Fig1b,
                            Fig1c, Fig1d,
                            common.legend = TRUE, 
                            legend = "top",
                            labels = c("a)", "b)", "c)", "d)"))

Fig1.combined2




#### Table 2 ####
#  Permanova of taxonomy
metadata <- as(sample_data(psTAXA), "data.frame")
dist <- phyloseq::distance(psTAXA.rar, method="bray") 

adonis2(dist ~ Sample_type,
        data = metadata, permutations = 9999)

# check dispersion
dispersion <- betadisper(dist, group=metadata$Sample_type)
anova(dispersion)
#permutest(dispersion, pairwise = TRUE, permutations = 9999)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
(mod.HSD <- TukeyHSD(dispersion))
plot(mod.HSD)
# everything is significant so need to be careful when interpretting results
pairwise.adonis(dist, metadata$Sample_type, perm = 9999, p.adjust.m="fdr")
#                     pairs Df SumsOfSqs F.Model     R2 p.value p.adjusted sig
#   1   Bulk_Soil vs Modern  1   0.09946   13.71 0.5549  0.0044    0.00440   *
#   2 Bulk_Soil vs Heritage  1   0.32137   21.52 0.6618  0.0027    0.00405   *
#   3    Modern vs Heritage  1   0.21004   15.62 0.4646  0.0001    0.00030  **

#  Permanova of functional genes
dist <- phyloseq::distance(psKEGG.rar, method="bray") 

adonis2(dist ~ Sample_type,
        data = metadata, permutations = 9999)

adonis <- adonis2(dist ~ Sample_type / Variety, data = metadata)

#           Df SumOfSqs      R2     F Pr(>F)    
# Sample_type    2 0.039362 0.71886 25.57  1e-04 ***
# Residual 20 0.015394 0.28114                 
# Total    22 0.054756 1.00000      


# check dispersion
dispersion <- betadisper(dist, group=metadata$Sample_type)
anova(dispersion)
#permutest(dispersion, pairwise = TRUE, permutations = 9999)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
(mod.HSD <- TukeyHSD(dispersion))
plot(mod.HSD)
# Significant between bulk soils and the two cultivars

# Pairwise
pairwise.adonis(dist, metadata$Sample_type, perm = 9999, p.adjust.m="fdr")




# Differential abundance analysis ####
# statistic analysis for Fig2, Fig3, Table 3, FigS3, Table S1, Table S2


# relevel(dds$Sample_type, ref = "Heritage")

# log2 fold change (MLE): Sample_type Heritage vs Bulk Soil 
# Wald test p-value: Sample_type Heritage vs Bulk Soil 
# What about the other comparisons?
resultsNames(ddsTAXA2)

ddsTAXA2_res.cols <- mcols(ddsTAXA2_res, use.names = TRUE)
summary(ddsTAXA2_res)

# Functional data
# create deseq object
ddsKEGG <- phyloseq_to_deseq2(psKEGG, ~ Sample_type)
ddsTAXA <- phyloseq_to_deseq2(psTAXA, ~ Sample_type)
# run deseq
ddsKEGG2 <- DESeq(ddsKEGG)
ddsTAXA2 <- DESeq(ddsTAXA)
# get size factors for normalization
sizeFactorsTAXA <- sizeFactors(ddsTAXA2)
write.csv(sizeFactorsTAXA, file = "sizeFactorsTAXA.csv")

sizeFactorsKEGG <- sizeFactors(ddsKEGG2)
write.csv(sizeFactorsKEGG, file = "sizeFactorsKEGG.csv")

# save results table
ddsTAXA2_res <- results(ddsTAXA2)
dds_res <- results(ddsKEGG2)
dds_res.cols <- mcols(dds_res, use.names = TRUE)
ddsTAXA2_res.cols <- mcols(ddsTAXA2_res, use.names = TRUE)


# Perform contrasts between groups
ddsTAXA2_HvM <- results(ddsTAXA2, contrast = c("Sample_type", "Heritage", "Modern"), alpha=0.05)
ddsTAXA2_HvB <- results(ddsTAXA2, contrast = c("Sample_type", "Heritage", "Bulk_Soil"), alpha=0.05)
ddsTAXA2_MvB <- results(ddsTAXA2, contrast = c("Sample_type", "Modern", "Bulk_Soil"), alpha=0.05)

# write.table(ddsTAXA2_HvM, file = 'data_DA_TAXA2results_HvM.txt', sep = ',')
# write.table(ddsTAXA2_HvB, file = 'data_DA_TAXA2results_HvB.txt', sep = ',')
# write.table(ddsTAXA2_MvB, file = 'data_DA_TAXA2results_MvB.txt', sep = ',')

ddsKEGG2_HvM <- results(ddsKEGG2, contrast = c("Sample_type", "Heritage", "Modern"), alpha=0.05)
ddsKEGG2_HvB <- results(ddsKEGG2, contrast = c("Sample_type", "Heritage", "Bulk_Soil"), alpha=0.05)
ddsKEGG2_MvB <- results(ddsKEGG2, contrast = c("Sample_type", "Modern", "Bulk_Soil"), alpha=0.05)

#write.table(ddsKEGG2_HvM, file = 'data_DA_results_HvM.txt', sep = ',')
#write.table(ddsKEGG2_HvB, file = 'data_DA_results_HvB.txt', sep = ',')
#write.table(ddsKEGG2_MvB, file = 'data_DA_results_MvB.txt', sep = ',')


#Apply shrinkage to account for rare or highly variable genes
#library("ashr")
#The R package 'ashr' implements an Empirical Bayes approach for large-scale hypothesis testing 
#and false discovery rate (FDR) estimation
# This procedure to moderate (or “shrink”) log2 fold changes from genes with very low counts and highly variable counts, 
# as can be seen by the narrowing of the vertical spread of points on the left side of the MA plot

#Taxa
TAXA2srk_HvM <- lfcShrink(ddsTAXA2, contrast=c("Sample_type", "Heritage", "Modern"), res=ddsTAXA2_HvM , type="ashr", lfcThreshold=1)
head(TAXA2srk_HvM)
#plot
DESeq2::plotMA(TAXA2srk_HvM, alpha = 0.05, ylim=c(-5,5), cex=.8)
abline(h=c(-1,1), col = "black", lwd=2)
table(TAXA2srk_HvM$svalue < 0.05) # 1414

TAXA2srk_HvB <- lfcShrink(ddsTAXA2, contrast=c("Sample_type", "Heritage", "Bulk_Soil"), res=ddsTAXA2_HvB, type="ashr", lfcThreshold=1)
#plot
DESeq2::plotMA(TAXA2srk_HvB, alpha = 0.05, ylim=c(-5,5), cex=.8)
abline(h=c(-1,1), col = "black", lwd=2)
table(TAXA2srk_HvB$svalue < 0.05) # 4975

TAXA2srk_MvB <- lfcShrink(ddsTAXA2, contrast=c("Sample_type", "Modern", "Bulk_Soil"), res=ddsTAXA2_MvB, type="ashr", lfcThreshold=1)
#plot
DESeq2::plotMA(TAXA2srk_MvB, alpha = 0.05, ylim=c(-5,5), cex=.8)
abline(h=c(-1,1), col="black", lwd=2)
table(TAXA2srk_MvB$svalue < 0.05) # 1114


# Function
KEGG2srk_HvM <- lfcShrink(ddsKEGG2, contrast=c("Sample_type", "Heritage", "Modern"), res=ddsKEGG2_HvM , type="ashr", lfcThreshold=1)
#plot
DESeq2::plotMA(KEGG2srk_HvM, alpha = 0.05, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col = "black", lwd=2)
table(KEGG2srk_HvM$svalue < 0.05) # 113

KEGG2srk_HvB <- lfcShrink(ddsKEGG2, contrast=c("Sample_type", "Heritage", "Bulk_Soil"), res=ddsKEGG2_HvB, type="ashr", lfcThreshold=1)
#plot
DESeq2::plotMA(KEGG2srk_HvB, alpha = 0.05, ylim=c(-5,5), cex=.8)
abline(h=c(-1,1), col = "black", lwd=2)
table(KEGG2srk_HvB$svalue < 0.05) # 1712

KEGG2srk_MvB <- lfcShrink(ddsKEGG2, contrast=c("Sample_type", "Modern", "Bulk_Soil"), res=ddsKEGG2_MvB, type="ashr", lfcThreshold=1)
#plot
DESeq2::plotMA(KEGG2srk_MvB, alpha = 0.05, ylim=c(-5,5), cex=.8)
abline(h=c(-1,1), col="black", lwd=2)
table(KEGG2srk_MvB$svalue < 0.05) # 194


# Join tables of significance with gene information #### 
#Taxon comparisons
HM.datTAXA2 = cbind(as(TAXA2srk_HvM, "data.frame"), as(tax_table(psTAXA)[rownames(TAXA2srk_HvM), ], "matrix"))
HB.datTAXA2 = cbind(as(TAXA2srk_HvB, "data.frame"), as(tax_table(psTAXA)[rownames(TAXA2srk_HvB), ], "matrix"))
MB.datTAXA2 = cbind(as(TAXA2srk_MvB, "data.frame"), as(tax_table(psTAXA)[rownames(TAXA2srk_MvB), ], "matrix"))

# write.table(HM.datTAXA2, file = 'data_DA_results_TAXA_srk_HvM.txt', sep = ',')
# write.table(HB.datTAXA2, file = 'data_DA_results_TAXA_srk_HvB.txt', sep = ',')
# write.table(MB.datTAXA2, file = 'data_DA_results_TAXA_srk_MvB.txt', sep = ',')

# Function comparions
HM.datKEGG2 = cbind(as(KEGG2srk_HvM, "data.frame"), as(tax_table(psKEGG)[rownames(KEGG2srk_HvM), ], "matrix"))
HB.datKEGG2 = cbind(as(KEGG2srk_HvB, "data.frame"), as(tax_table(psKEGG)[rownames(KEGG2srk_HvB), ], "matrix"))
MB.datKEGG2 = cbind(as(KEGG2srk_MvB, "data.frame"), as(tax_table(psKEGG)[rownames(KEGG2srk_MvB), ], "matrix"))

# write.table(HM.datKEGG2, file = 'data_DA_results_srk_HvM.txt', sep = ',')
# write.table(HB.datKEGG2, file = 'data_DA_results_srk_HvB.txt', sep = ',')
# write.table(MB.datKEGG2, file = 'data_DA_results_srk_MvB.txt', sep = ',')



#### Make a phyloseq object of significant genes between wheat ####
# the aim here is to subset the full data set and leave only the genes that are significant between
# heritage and modern wheat. Using the subsetted phyloseq object I can sum up the read counts

# First get a list of the significant genes
HM.sigTaxaList <- unique(
  rownames(subset(TAXA2srk_HvM, svalue<=0.05))) # up- and down-regulated in heritage
length(HM.sigTaxaList) #1414


sig.HM.datKEGG2 <- subset(HM.datKEGG2,svalue < alpha)
capture.output(sig.HM.datKEGG2, file = "sig.HM.datKEGG2.txt")
dim(sig.HM.datKEGG2)
list.sig.HM <- row.names(sig.HM.datKEGG2)

#then prune taxa
my_subset <- subset(otu_table(psTAXA), rownames(otu_table(psTAXA)) %in% HM.sigTaxaList)
dim(my_subset)
psTAXA.sigHM <- merge_phyloseq(my_subset, tax_table(psTAXA), sample_data(psTAXA))
print(psTAXA.sigHM)

psKEGG.sigHM = prune_taxa(list.sig.HM, psKEGG)
microbiome::summarize_phyloseq(psKEGG.sigHM)
# a crude look at the total number of reads for each taxa
summary(taxa_sums(psTAXA.sigHM))
summary(taxa_sums(psKEGG.sigHM))

# look at gene counts (Function only)
tdt.sigHM = data.table(
                 TotalCounts = taxa_sums(psKEGG.sigHM),
                 genes = taxa_names(psKEGG.sigHM))
tdt.sigHM

library(tidyverse)
tdt.sigHM2 <- tdt.sigHM %>% remove_rownames %>% column_to_rownames(var="genes")

# merge the results table with the total counts of each gene into one dataframe
sig.HM.sums.datKEGG2 <- merge(sig.HM.datKEGG2, tdt.sigHM2,
                          by = 'row.names', all = FALSE)
head(sig.HM.sums.datKEGG2)

# write to csv - final results for differential abundance between the two wheat cultivar types
write.csv(sig.HM.sums.datKEGG2,"sig_log2fold1_sums_HM_KEGG_annotree.csv", row.names = FALSE)


#### Table of results for each wheat vs bulk soil - Function only ####
# First get a list of the significant genes
sig.HB.datKEGG2 <- subset(HB.datKEGG2,svalue < alpha)
sig.MB.datKEGG2 <- subset(MB.datKEGG2,svalue < alpha)
dim(sig.HB.datKEGG2)
dim(sig.MB.datKEGG2)
list.sig.HB <- row.names(sig.HB.datKEGG2)
list.sig.MB <- row.names(sig.MB.datKEGG2)
#then prune taxa
psKEGG.sigHB = prune_taxa(list.sig.HB, psKEGG)
psKEGG.sigMB = prune_taxa(list.sig.MB, psKEGG)
# a crude look at the total number of reads for each gene
summary(taxa_sums(psKEGG.sigHB))
summary(taxa_sums(psKEGG.sigMB))
# look at gene counts
tdt.sigHB = data.table(
  TotalCounts = taxa_sums(psKEGG.sigHB),
  genes = taxa_names(psKEGG.sigHB))
tdt.sigHB

tdt.sigMB = data.table(
  TotalCounts = taxa_sums(psKEGG.sigMB),
  genes = taxa_names(psKEGG.sigMB))
tdt.sigMB

#library(tidyverse)
tdt.sigHB2 <- tdt.sigHB %>% remove_rownames %>% column_to_rownames(var="genes")
tdt.sigMB2 <- tdt.sigMB %>% remove_rownames %>% column_to_rownames(var="genes")
# merge the results table with the total counts of each gene into one dataframe
sig.HB.sums.datKEGG2 <- merge(sig.HB.datKEGG2, tdt.sigHB2,
                              by = 'row.names', all = FALSE)
sig.MB.sums.datKEGG2 <- merge(sig.MB.datKEGG2, tdt.sigMB2,
                              by = 'row.names', all = FALSE)


#get normalised total counts for each Sample_type
KEGG2SumPerSample_type <- sapply( levels(ddsKEGG2$Sample_type), function(lvl) rowSums(counts(ddsKEGG2,normalized=TRUE)[,ddsKEGG2$Sample_type == lvl, drop=F] ) )
head(KEGG2SumPerSample_type)

KEGG2SumPerSample_typeSigHB <- KEGG2SumPerSample_type[rownames(KEGG2SumPerSample_type) %in% list.sig.HB, ] 
KEGG2SumPerSample_typeSigMB <- KEGG2SumPerSample_type[rownames(KEGG2SumPerSample_type) %in% list.sig.MB, ] 

# merge the results table with the total counts of each gene into one dataframe
sig.HB.sums.datKEGG3 <- sig.HB.sums.datKEGG2 %>% remove_rownames %>% column_to_rownames(var="Row.names")
sig.MB.sums.datKEGG3 <- sig.MB.sums.datKEGG2 %>% remove_rownames %>% column_to_rownames(var="Row.names")

sig.HB.sums.datKEGG2.Sample_types <- merge(sig.HB.sums.datKEGG3, KEGG2SumPerSample_typeSigHB,
                              by = 'row.names', all = FALSE)
sig.MB.sums.datKEGG2.Sample_types <- merge(sig.MB.sums.datKEGG3, KEGG2SumPerSample_typeSigMB,
                              by = 'row.names', all = FALSE)


# write to csv to do more research on these genes
#write.csv(sig.HB.sums.datKEGG2.Sample_types,"sig_genes_heritage_vs_bulksoil.csv", row.names = FALSE)
#write.csv(sig.MB.sums.datKEGG2.Sample_types,"sig_genes_modern_vs_bulksoil.csv", row.names = FALSE)



# core genes for plots ####
# get 10 housekeeping (core) prokaryote genes to add to graphs
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6984700/ for justification of genes (rho, rpoD and gyrA
head(tax_table(psKEGG))
# need to subset by Rank5
# "K03060 DNA-directed RNA polymerase subunit omega [EC:2.7.7.6]" #(rpoZ)
# "K02470 DNA gyrase subunit B [EC:5.6.2.2]" #(gyrB)
# "K03106 signal recognition particle subunit SRP54 [EC:3.6.5.4]" #(ffh)
# "K01915 glutamine synthetase [EC:6.3.1.2]" (glnA)
# "K03628 transcription termination factor Rho" (rho)
# "K02939 large subunit ribosomal protein L9" (rplI)
# "K02879 large subunit ribosomal protein L17" (rplQ)
# "K03168 DNA topoisomerase I [EC:5.6.2.1]" (topA)
# "K01696 tryptophan synthase beta chain [EC:4.2.1.20]"
# "K02112 F-type H+/Na+-transporting ATPase subunit beta [EC:7.1.2.2 7.2.2.1]" (atpD)


ps_coreGenes <- subset_taxa(psKEGG, Rank5=="K03060 DNA-directed RNA polymerase subunit omega [EC:2.7.7.6]" |
                              Rank5== "K02470 DNA gyrase subunit B [EC:5.6.2.2]"|
                              Rank5== "K03106 signal recognition particle subunit SRP54 [EC:3.6.5.4]"|
                              Rank5== "K01915 glutamine synthetase [EC:6.3.1.2]"|
                              Rank5== "K03628 transcription termination factor Rho"|
                              Rank5== "K02939 large subunit ribosomal protein L9"|
                              Rank5== "K02879 large subunit ribosomal protein L17"|
                              Rank5== "K03168 DNA topoisomerase I [EC:5.6.2.1]" |
                              Rank5== "K02112 F-type H+/Na+-transporting ATPase subunit beta [EC:7.1.2.2 7.2.2.1]"|
                              Rank5== "K01696 tryptophan synthase beta chain [EC:4.2.1.20]"
                            )
tax_table(ps_coreGenes)


# getting row names 
rowsCore <- rownames(tax_table(ps_coreGenes))

# extracting data frame rows
HM.KEGG2core <- HM.datKEGG2[rownames(HM.datKEGG2) %in% rowsCore, ] 
print(HM.KEGG2core)

# Fig. 2a ####
# Ternary plot - NOTE run function script in tern_e 2.R
# need to get the means of each compartment for the ternary plots 
TAXA2baseMeanPerLvl <- sapply( levels(ddsTAXA2$Sample_type), function(lvl) rowMeans( counts(ddsTAXA2,normalized=TRUE)[,ddsTAXA2$Sample_type == lvl, drop=F] ) )
head(TAXA2baseMeanPerLvl)
sort(colSums(TAXA2baseMeanPerLvl), dec=T)
# These are normalised means

# define colorcoding for ternary plot
# make a list of taxa that are higher in abundance in each group
# thresholds for the svalue after shrinkage should be lower than the normal 0.05 for p adjusted values #https://support.bioconductor.org/p/133091/
# but I already but a threshold of 1 for the log fold change so I think this is overkill
# from help ?lfcShrink
# lfcThreshold- a non-negative value which specifies a log2 fold change threshold (as in results). This can be used with any shrinkage type. It will produce new p-values or s-values testing whether the LFC is greater in absolute value than the threshold. The s-values returned in combination with apeglm or ashr provide the probability of FSOS events, "false sign or small", among the tests with equal or smaller s-value than a given gene's s-value, where "small" is specified by lfcThreshold.
# svalue- logical, should p-values and adjusted p-values be replaced with s-values when using apeglm or ashr. s-values provide the probability of false signs among the tests with equal or smaller s-value than a given given's s-value. See Stephens (2016) reference on s-values.

log2cutoff <- 1
svaluecutoff <- 0.05

plantTaxa <- unique(c(
  rownames(subset(TAXA2srk_HvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)), # upregulated in heritage
  rownames(subset(TAXA2srk_MvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)))) # upregulated in modern
length(plantTaxa) #3622


heritageTaxa <- unique(
  #rownames(subset(TAXA2srk_HvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)), # upregulated in heritage
  rownames(subset(TAXA2srk_HvM, svalue<=svaluecutoff & log2FoldChange>=log2cutoff))) # upregulated in heritage
length(heritageTaxa) #1344


modernTaxa <- unique(
  rownames(subset(TAXA2srk_HvM, svalue<=svaluecutoff & log2FoldChange<=-1))) # downregulated in heritage
#rownames(subset(TAXA2srk_MvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)))) # upregulated in modern
length(modernTaxa) #70


soilTaxa <- unique(c(
  rownames(subset(TAXA2srk_HvB, svalue<=svaluecutoff & log2FoldChange<=-1)), # downregulated in heritage
  rownames(subset(TAXA2srk_MvB, svalue<=svaluecutoff & log2FoldChange<=-1)))) # upregulated in modern
length(soilTaxa) #1382

# It was too messy to have the non significant taxa on the graph so 
#I am subsetting the data to be only significant and changing the colours
sigTaxa <- c(plantTaxa,soilTaxa,heritageTaxa,modernTaxa)
length(unique(sigTaxa)) #6418 in total but 5072 after adding unique only

TAXA2baseMeanPerLvlSig <- TAXA2baseMeanPerLvl[rownames(TAXA2baseMeanPerLvl) %in% sigTaxa, ] 
dim(TAXA2baseMeanPerLvlSig) #5072

# make a phyloseq object with all significant genes 

#then prune taxa
psTAXA.sigALL = prune_taxa(sigTaxa, psTAXA)

# set colours based on significant Taxa
# # The colour-blind friendly palette with grey:

TAXA2ternary_colors <- ifelse(rownames(TAXA2baseMeanPerLvlSig) %in% plantTaxa, "#004d40","#ffc107")
names(TAXA2ternary_colors) <- rownames(TAXA2baseMeanPerLvlSig)
TAXA2ternary_colors[heritageTaxa] <- "#d81b60"
TAXA2ternary_colors[modernTaxa] <- "#1e88e5"

### Plotting ternary with colours
tern_e(TAXA2baseMeanPerLvlSig, prop=F, col=TAXA2ternary_colors, main = "Taxa",cex = .4)



# Fig. 2b ####
# Ternary plot - NOTE run function script in tern_e 2.R
# need to get the means of each compartment for the ternary plots 
KEGG2baseMeanPerLvl <- sapply( levels(ddsKEGG2$Sample_type), function(lvl) rowMeans( counts(ddsKEGG2,normalized=TRUE)[,ddsKEGG2$Sample_type == lvl, drop=F] ) )
head(KEGG2baseMeanPerLvl)
dim(KEGG2baseMeanPerLvl)
sort(colSums(KEGG2baseMeanPerLvl), dec=T)
# These are normalised means

logKEGG2baseMeanPerLvl <- log(KEGG2baseMeanPerLvl+1)
summary(logKEGG2baseMeanPerLvl)

ddsKEGG2lognorm <- normTransform(ddsKEGG2, f = log2, pc = 1)
ddsKEGG2lognormMeans <- sapply( levels(ddsKEGG2lognorm$Sample_type), function(lvl) rowMeans( assay(ddsKEGG2lognorm,normalized=TRUE)[,ddsKEGG2lognorm$Sample_type == lvl, drop=F] ) )
summary(ddsKEGG2lognormMeans)

# plot ternary without colours for significant genes
tern_e(KEGG2baseMeanPerLvl, prop=T, main = "Genes")
tern_e(logKEGG2baseMeanPerLvl, prop=T, main = "Genes")
tern_e(ddsKEGG2lognormMeans, prop=T, main = "Genes")


# define colourcoding for ternary plot
# make a list of genes that are higher in abundance in each group
# thresholds for the svalue after shrinkage should be lower than the normal 0.05 for p adjusted values #https://support.bioconductor.org/p/133091/
# but I already but a threshold of 1 for the log fold change so I think this is overkill
# from help ?lfcShrink
# lfcThreshold- a non-negative value which specifies a log2 fold change threshold (as in results). This can be used with any shrinkage type. It will produce new p-values or s-values testing whether the LFC is greater in absolute value than the threshold. The s-values returned in combination with apeglm or ashr provide the probability of FSOS events, "false sign or small", among the tests with equal or smaller s-value than a given gene's s-value, where "small" is specified by lfcThreshold.
# svalue- logical, should p-values and adjusted p-values be replaced with s-values when using apeglm or ashr. s-values provide the probability of false signs among the tests with equal or smaller s-value than a given given's s-value. See Stephens (2016) reference on s-values.


log2cutoff <- 1
svaluecutoff <- 0.05

plantGenes <- unique(c(
  rownames(subset(KEGG2srk_HvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)), # upregulated in tall
  rownames(subset(KEGG2srk_MvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)))) # upregulated in semi dwarf
length(plantGenes) #1134


heritageGenes <- unique(
  #rownames(subset(KEGG2srk_HvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)), # upregulated in tall
  rownames(subset(KEGG2srk_HvM, svalue<=svaluecutoff & log2FoldChange>=log2cutoff))) # upregulated in tall
length(heritageGenes) #107


modernGenes <- unique(
  rownames(subset(KEGG2srk_HvM, svalue<=svaluecutoff & log2FoldChange<=-1))) # downregulated in tall
  #rownames(subset(KEGG2srk_MvB, svalue<=svaluecutoff & log2FoldChange>=log2cutoff)))) # upregulated in semi dwarf
length(modernGenes) #6


soilGenes <- unique(c(
  rownames(subset(KEGG2srk_HvB, svalue<=svaluecutoff & log2FoldChange<=-1)), # downregulated in tall
  rownames(subset(KEGG2srk_MvB, svalue<=svaluecutoff & log2FoldChange<=-1)))) # upregulated in semi dwarf
length(soilGenes) #580


# It was too messy to have the non significant genes on the graph so 
#I am subsetting the data to be only significant and changing the colours
sigGenes <- c(plantGenes,soilGenes,heritageGenes,modernGenes)
length(unique(sigGenes)) #1719 when added unique


# add the core genes
sigcoreGenes <- c(sigGenes,rowsCore)

KEGG2baseMeanPerLvlSig <- KEGG2baseMeanPerLvl[rownames(KEGG2baseMeanPerLvl) %in% sigcoreGenes, ] 
dim(KEGG2baseMeanPerLvlSig) #1729 10 core genes added



# set colours based on significant genes
#KEGG2ternary_colors <- ifelse(rownames(KEGG2baseMeanPerLvl) %in% plantGenes, "#f0e442","darkgrey")
KEGG2ternary_colors <- ifelse(rownames(KEGG2baseMeanPerLvlSig) %in% plantGenes, "#004d40","#ffc107")
#names(KEGG2ternary_colors) <- rownames(KEGG2baseMeanPerLvl)
names(KEGG2ternary_colors) <- rownames(KEGG2baseMeanPerLvlSig)
#KEGG2ternary_colors[soilGenes] <- "#009e73"
KEGG2ternary_colors[heritageGenes] <- "#d81b60"
KEGG2ternary_colors[modernGenes] <- "#1e88e5"
KEGG2ternary_colors[rowsCore] <- "black" # these are the core genes

### Plotting ternary with colours
tern_e(KEGG2baseMeanPerLvlSig, prop=F, col=KEGG2ternary_colors, main = "Genes",cex = .4)


### Fig. 2 combined ####
# legend in between the plots
library(grid)
library(stringr)
library(tidyverse)

# Define your legend elements
legend_names <- c("Heritage", "Modern", "Bulk soil", "All wheat", "Housekeeping genes")
legend_colors <- c("#d81b60", "#1e88e5", "#ffc107", "#004d40", "black")

# Open a PDF file
pdf("combined_ternary_plots.pdf", width = 10, height = 5)

# Create a viewport for the entire plot area
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

# Plot 1 (left)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
tern_e(TAXA2baseMeanPerLvlSig, prop_size = FALSE, col = TAXA2ternary_colors,
       main = "a) Taxonomy", cex = 0.4, newpage = FALSE, pop = FALSE)
upViewport()

# Plot 2 (right)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
tern_e(KEGG2baseMeanPerLvlSig, prop_size = FALSE, col = KEGG2ternary_colors,
       main = "b) Function", cex = 0.4, newpage = FALSE, pop = FALSE)
upViewport()

# Create a new viewport that spans both columns for the legend in the middle
pushViewport(viewport(x = 0.5, y = 0.55, width = 0.15, height = 0.6, just = c("center", "center")))

# Add a semi-transparent white background for better readability
grid.rect(gp = gpar(fill = "white", alpha = 0.7,col = NA))

# Add title
grid.text("Enriched in:", x = 0.5, y = 0.95, just = c("center", "top"),
          gp = gpar(fontface = "bold", fontsize = 9))

# Add legend items vertically
for (i in 1:length(legend_names)) {
  y_pos = 0.85 - (i-1) * 0.1
  grid.points(x = 0.2, y = y_pos, pch = 19, size = unit(0.8, "char"),
              gp = gpar(col = legend_colors[i]))
  grid.text(legend_names[i], x = 0.3, y = y_pos, just = "left",
            gp = gpar(fontsize = 8))
}

# Close viewports
popViewport(2)

# Close the PDF file
dev.off()

### Fig. 3a ####
# graph DESeq logfold results
# from a tutorial https://joey711.github.io/phyloseq-extensions/DESeq2.html

#for this plot we are only interested in the differences between Modern and heritage

# thresholds for the svalue after shrinkage should be lower than the normal 0.05 for p adjusted values #https://support.bioconductor.org/p/133091/
# but I already but a threshold of 1 for the log fold change so I think this is overkill
# from help ?lfcShrink
# lfcThreshold- a non-negative value which specifies a log2 fold change threshold (as in results). This can be used with any shrinkage type. It will produce new p-values or s-values testing whether the LFC is greater in absolute value than the threshold. The s-values returned in combination with apeglm or ashr provide the probability of FSOS events, "false sign or small", among the tests with equal or smaller s-value than a given gene's s-value, where "small" is specified by lfcThreshold.
# svalue- logical, should p-values and adjusted p-values be replaced with s-values when using apeglm or ashr. s-values provide the probability of false signs among the tests with equal or smaller s-value than a given given's s-value. See Stephens (2016) reference on s-values.
alpha = 0.05
sigtabHM = HM.datTAXA2[which(HM.datTAXA2$svalue < alpha), ]
sigtabHM$Sample_type <- ifelse(sigtabHM$log2FoldChange > 1, "Heritage", "Modern")
dim(sigtabHM)

# write to csv to do more research on these Taxa
write.csv(sigtabHM,"sig_log2fold_HM_TAXA_annotree.csv", row.names = FALSE)

#### All Phyla to class level ####

levels(as.factor(sigtabHM$Phylum))
sigtabHM$Phylum3 <- sigtabHM$Phylum
# I want to change the name of all the Candidatus taxa so they group together in the plot
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Harrisonbacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "candidate division FCPU426", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Adlerbacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Blackallbacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Magasanikbacteria" , "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Omnitrophica", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Tagabacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Uhrbacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Zambryskibacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Kaiserbacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Nomurabacteria" , "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Saccharibacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Thermoplasmatota", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Yanofskybacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Moranbacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Peregrinibacteria", "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Taylorbacteria" , "Candidate phyla", as.character(Phylum3)))
sigtabHM$Phylum3 <- with(sigtabHM, ifelse(Phylum == "Candidatus Woesebacteria", "Candidate phyla", as.character(Phylum3)))

levels(as.factor(sigtabHM$Phylum3))


# I want to have all phylum at class level except for candidate phyla
sigtabHM$Phylum3.Class <- ifelse(is.na(sigtabHM$Class), paste(sigtabHM$Phylum3, "-","Unclassified"), paste(sigtabHM$Phylum3, "-", sigtabHM$Class))
levels(as.factor(sigtabHM$Phylum3.Class))

# ifelse to get only candidate phyla back to phyla level level
levels(as.factor(sigtabHM$Phylum3))
sigtabHM$Phylum3.Class <- with(sigtabHM, ifelse(Phylum3 == "Candidate phyla", as.character(Phylum3), as.character(Phylum3.Class)))
levels(as.factor(sigtabHM$Phylum3.Class))

# we want to know how many rows (taxa) in each class
detach(package:plyr)
library(dplyr)
head(sigtabHM)
str(sigtabHM)

n_Phylum3.Class <- sigtabHM %>% 
  group_by(Phylum3.Class) %>%
  summarise(OTUs = length(Phylum3.Class))

n_Phylum3.Class$xaxis <- "xaxis"

# Class order highest to lowest abund
x = tapply(n_Phylum3.Class$OTUs, n_Phylum3.Class$Phylum3.Class, function(x) max(x))
x = sort(x, FALSE) # decreasing = false
n_Phylum3.Class$Phylum3.Class = factor(as.character(n_Phylum3.Class$Phylum3.Class), levels=names(x))
sigtabHM$Phylum3.Class = factor(as.character(sigtabHM$Phylum3.Class), levels=names(x))

# put unclassified at the bottom and group phyla together (still based on the number of taxa)
levels(sigtabHM$Phylum3.Class)
# do for both datasets
sigtabHM$Phylum3.Class <- factor(sigtabHM$Phylum3.Class, 
                               levels = c("Candidate phyla",
                                          "Acidobacteria - Unclassified",
                                          "Euryarchaeota - Unclassified","Euryarchaeota - Halobacteria","Euryarchaeota - Methanomicrobia",
                                          "Ignavibacteriae - Ignavibacteria",
                                          "Planctomycetes - Unclassified","Planctomycetes - Phycisphaerae","Planctomycetes - Planctomycetia",
                                          "Thermotogae - Thermotogae",
                                          "Cyanobacteria - Unclassified",
                                          "Thaumarchaeota - Unclassified",
                                          "Elusimicrobia - Unclassified",
                                          "Verrucomicrobia - Verrucomicrobiae",
                                          "Spirochaetes - Unclassified","Spirochaetes - Spirochaetia",
                                          "Firmicutes - Unclassified","Firmicutes - Erysipelotrichia",
                                          "Firmicutes - Bacilli","Firmicutes - Clostridia",
                                          "Bacteroidetes - Unclassified",
                                          "Bacteroidetes - Bacteroidia","Bacteroidetes - Chitinophagia","Bacteroidetes - Cytophagia",
                                          "Bacteroidetes - Flavobacteriia","Bacteroidetes - Sphingobacteriia",
                                          "Actinobacteria - Coriobacteriia",
                                          "Actinobacteria - Actinomycetia",
                                          "Proteobacteria - Unclassified",
                                          "Proteobacteria - Hydrogenophilalia",
                                          "Proteobacteria - Oligoflexia",
                                          "Proteobacteria - Gammaproteobacteria", 
                                          "Proteobacteria - Betaproteobacteria",
                                          "Proteobacteria - Alphaproteobacteria"
                               ))

n_Phylum3.Class$Phylum3.Class <- factor(n_Phylum3.Class$Phylum3.Class, 
                                        levels = c("Candidate phyla",
                                                   "Acidobacteria - Unclassified",
                                                   "Euryarchaeota - Unclassified","Euryarchaeota - Halobacteria","Euryarchaeota - Methanomicrobia",
                                                   "Ignavibacteriae - Ignavibacteria",
                                                   "Planctomycetes - Unclassified","Planctomycetes - Phycisphaerae","Planctomycetes - Planctomycetia",
                                                   "Thermotogae - Thermotogae",
                                                   "Cyanobacteria - Unclassified",
                                                   "Thaumarchaeota - Unclassified",
                                                   "Elusimicrobia - Unclassified",
                                                   "Verrucomicrobia - Verrucomicrobiae",
                                                   "Spirochaetes - Unclassified","Spirochaetes - Spirochaetia",
                                                   "Firmicutes - Unclassified","Firmicutes - Erysipelotrichia",
                                                   "Firmicutes - Bacilli","Firmicutes - Clostridia",
                                                   "Bacteroidetes - Unclassified",
                                                   "Bacteroidetes - Bacteroidia","Bacteroidetes - Chitinophagia","Bacteroidetes - Cytophagia",
                                                   "Bacteroidetes - Flavobacteriia","Bacteroidetes - Sphingobacteriia",
                                                   "Actinobacteria - Coriobacteriia",
                                                   "Actinobacteria - Actinomycetia",
                                                   "Proteobacteria - Unclassified",
                                                   "Proteobacteria - Hydrogenophilalia",
                                                   "Proteobacteria - Oligoflexia",
                                                   "Proteobacteria - Gammaproteobacteria", 
                                                   "Proteobacteria - Betaproteobacteria",
                                                   "Proteobacteria - Alphaproteobacteria"
                                        ))



taxa.abund.plot <- ggplot(data=n_Phylum3.Class, aes(x=xaxis,y=Phylum3.Class)) +
  geom_text(data=n_Phylum3.Class, aes(label=OTUs)) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14),
        #panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(colour = "grey92") ,
        plot.margin = unit(c(6.05,0.2,2,-1), 'lines'), #unit(c(top, right, bottom, left), units) The default value for mar is c(5.1, 4.1, 4.1, 2.1).
        legend.position="top"
  ) +
  xlab("") +
  ylab("")
taxa.abund.plot



library(ggbreak) # break the x axis

taxa.logplot<-ggplot(sigtabHM, aes(x=log2FoldChange, y=Phylum3.Class, fill = log2FoldChange)) + 
  geom_point(size=3,shape=21) + #, position = "jitter") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        legend.position="top",
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        plot.margin = unit(c(0.2,-1,0,0.2), 'lines'),
        text = element_text(size = 14),
        axis.title.x = element_text(hjust=0.9),
        panel.grid.major.y = element_line(colour = "grey92") ,
        panel.grid.major.x = element_line(colour = "grey92") ,
        panel.background = element_blank(),
        axis.line.x.top = element_blank()
  ) +
  scale_x_break(c(-6, -17.5), ticklabels=c(-6,-3,0,3,6)) + 
  scale_x_continuous(limits = c(-19, 6.5),breaks=-18)+
  scale_fill_gradientn(colours = c("#1e88e5", "white", "#d81b60"),
                       name = "Enriched in",
                       values = scales::rescale(c(-15,-2.5, 0, 2.5, 5)),
                       #breaks = unlist(list_val2),
                       breaks = c(-17, 0, 5),
                       labels = c("Modern", "","Heritage")
  ) +
  ylab("") +
  xlab("log2 fold change") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

taxa.logplot



# ggarrange doesn't work with splitting the x axis
# patchwork is the only package I found to be compatible but it wouldn't line up the y axis
# I have had to adjust the margins on the abund plot to squash it down
# then combine
library(patchwork)


#Combine plots
Fig3a <- taxa.logplot +
  taxa.abund.plot +     # had to play around with the margins
  plot_layout(ncol = 2,widths = c(1, .1))

Fig3a






### Fig 3b  ####
# graph DESeq logfold results
# from a tutorial https://joey711.github.io/phyloseq-extensions/DESeq2.html

#for this plot we are only interested in the differences between Modern and heritage

# thresholds for the svalue after shrinkage should be lower than the normal 0.05 for p adjusted values #https://support.bioconductor.org/p/133091/
# but I already but a threshold of 1 for the log fold change so I think this is overkill
# from help ?lfcShrink
# lfcThreshold- a non-negative value which specifies a log2 fold change threshold (as in results). This can be used with any shrinkage type. It will produce new p-values or s-values testing whether the LFC is greater in absolute value than the threshold. The s-values returned in combination with apeglm or ashr provide the probability of FSOS events, "false sign or small", among the tests with equal or smaller s-value than a given gene's s-value, where "small" is specified by lfcThreshold.
# svalue- logical, should p-values and adjusted p-values be replaced with s-values when using apeglm or ashr. s-values provide the probability of false signs among the tests with equal or smaller s-value than a given given's s-value. See Stephens (2016) reference on s-values.


alpha = 0.05
sigtabKEGG2 = HM.datKEGG2[which(HM.datKEGG2$svalue < alpha), ]
dim(sigtabKEGG2)

# add the core genes to the sig genes
sig_and_coreKEGG2 <- rbind(sigtabKEGG2, HM.KEGG2core)
dim(sig_and_coreKEGG2)
str(sig_and_coreKEGG2)

write.csv(sig_and_coreKEGG2,"sig_and_core_log2fold_HM_KEGG_annotree.csv", row.names = T)

# We have added a function column to split up the broad functions from KEGG like "enzymes" and "metabolism"
sig_and_coreKEGG22 <- read.csv("sig_and_core_log2fold_HM_KEGG_annotree_function.csv")


# we want to know how many rows (genes) in each pathway
detach(package:plyr)
library(dplyr)

n_Path <- sig_and_coreKEGG22 %>% 
  group_by(Function) %>%
  summarise(genes = length(Function))

n_Path$xaxis <- "xaxis"

# Function order highest to lowest
x = tapply(n_Path$genes, n_Path$Function, function(x) max(x))
x = sort(x, FALSE)
n_Path$Function = factor(as.character(n_Path$Function), levels=names(x))
sig_and_coreKEGG22$Function = factor(as.character(sig_and_coreKEGG22$Function), levels=names(x))

# manually move housekeeping genes to the bottom but keep teh rest in the order of genes
levels(sig_and_coreKEGG22$Function)
# do for both sig_and_coreKEGG22 and n_Path
n_Path$Function <- factor(n_Path$Function, 
                                               levels = c("Housekeeping",
                                                          "Calcium sensors", "Secondary messaging","Cell wall degradation","Defence",                   
                                                          "Motility","Antibiotic resistance","Cell cycle","Transcriptional regulation","Quorum sensing",           
                                                          "Sensor regulator","Biofilm formation","Secondary metabolism","Secretion systems","Membrane transport",         
                                                          "Primary metabolism"))

sig_and_coreKEGG22$Function <- factor(sig_and_coreKEGG22$Function, 
                          levels = c("Housekeeping",
                                     "Calcium sensors", "Secondary messaging","Cell wall degradation","Defence",                   
                                     "Motility","Antibiotic resistance","Cell cycle","Transcriptional regulation","Quorum sensing",           
                                     "Sensor regulator","Biofilm formation","Secondary metabolism","Secretion systems","Membrane transport",         
                                     "Primary metabolism"))


# we decided plotting the number would be more easily interpreted. 
kegg.abund.plot <- ggplot(data=n_Path, aes(x=xaxis,y=Function)) +
  geom_text(data=n_Path, aes(label=genes)) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92") ,
        panel.background = element_blank(),
        text = element_text(size = 14),
        panel.grid.major.x = element_blank() ,
        plot.margin = unit(c(0.2,0.2,0,0), 'lines'), #unit(c(top, right, bottom, left), units) The default value for mar is c(5.1, 4.1, 4.1, 2.1).
        legend.position="top"
  ) +
  xlab("") +
  ylab("")
kegg.abund.plot

kegg.logplot<-ggplot(sig_and_coreKEGG22, aes(x=log2FoldChange, y=Function, fill = log2FoldChange)) + 
  #theme_bw()+
  geom_point(size=3,shape=21) + #, position = "jitter") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        plot.margin = unit(c(0.2,0,0,0.2), 'lines'),
        text = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92") ,
        panel.grid.major.x = element_line(colour = "grey92") ,
        legend.position="top"#,
        #legend.title=element_blank()
        ) +
  #scale_fill_gradient2(low = "#d81b60", mid = "white", high = "#1e88e5", midpoint = 0) +
  scale_fill_gradientn(colours = c("#1e88e5", "white", "#d81b60"),
                       name = "Enriched in",
                       #values = scales::rescale(c(-15,-7, 0, 2.5, 7)),
                       breaks = c(-3,0,3),
                       #                      #breaks = c(-15, 0, 5),
                       labels = c("Modern", "","Heritage")
  ) +
ylab("") +
  xlab("log2 fold change") +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))#+
  # scale_x_continuous(limits = c(-3.5, 5.5))+
  # geom_vline(xintercept = 4, linetype="solid", 
  #            color = "black", linewidth=0.5) +
  # geom_point(data=n_Path, aes(x=log2FoldChange,y=Function,size = genes))

kegg.logplot



#Combine plots
Fig3b <- kegg.logplot +
  kegg.abund.plot +     # had to play around with the margins
  plot_layout(ncol = 2,widths = c(1, .1))

Fig3b

#### Fig. 3 combined ####
library(patchwork)
library(gridExtra)
library(grid)

# Create plots 3 and 4 as a patchwork
kegg_combined <- kegg.logplot + kegg.abund.plot + 
  plot_layout(ncol = 2, widths = c(1, 0.1))
kegg_grob <- patchworkGrob(kegg_combined)

# Create special grobs for plots 1 and 2
g1 <- create_special_grob(taxa.logplot)
g2 <- create_special_grob(taxa.abund.plot)

# Create the layout as before
pdf("Figure3.pdf", width = 12, height = 7)
grid.newpage()

lay <- grid.layout(nrow = 1, ncol = 3, 
                   widths = unit(c(1.3, 0.1, 1.1), "null"))
pushViewport(viewport(layout = lay))

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
grid.draw(g1)
popViewport()

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
grid.draw(g2)
popViewport()

pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 1))
grid.draw(kegg_grob)
popViewport()

popViewport()
dev.off()

#### Table 3 ####
#get normalised total counts for each Sample_type
# Use the deseq object so that we have the same normalisation
# first do on all genes and then we subset for the differentially abundant between tall and short
 
KEGG2SumPerSample_type <- sapply(levels(ddsKEGG2$Sample_type), function(lvl) rowSums(counts(ddsKEGG2,normalized=TRUE)[,ddsKEGG2$Sample_type == lvl, drop=F] ) )
head(KEGG2SumPerSample_type)

KEGG2SumPerSample_typeSig <- KEGG2SumPerSample_type[rownames(KEGG2SumPerSample_type) %in% list.sig.HM, ] 
dim(KEGG2SumPerSample_typeSig) 


write.csv(KEGG2SumPerSample_typeSig,"KEGG2SumPerSample_typeSig.csv", row.names = TRUE)


# This was used to calculate the total read count and proportion of reads belonging to Heritage or Modern rhizosphere samples
# The other columns were calculated from Loop_script.R

 #### Supplementary figures and Tables ####
# Rarefaction curves Fig S1
# stop R printing scientific numbers
options(scipen = 100, digits = 4)

# Fig. S1a ####
taxa.rar.plot <- ggrare(psTAXA, step = 100000, color = "Sample_type", se = FALSE)

taxa.rar.plot2 <- taxa.rar.plot  + theme_bw() + 
  scale_color_manual(name = "Sample type",
                     breaks = c("Bulk_Soil", "Modern","Heritage"),
                     values = c("#004d40","#1e88e5","#d81b60"),
                     labels = c("Bulk soil", "Modern","Heritage")) +
  theme(axis.text.x = element_text(size = 8,angle = 35, hjust = 1)) +
  ylab("Taxa richness") +
  ggtitle("a) Taxonomy")
taxa.rar.plot2

  # Fig. S1b ####
 kegg.rar.plot <- ggrare(psKEGG, step = 100000, color = "Sample_type", se = FALSE)
 
kegg.rar.plot2 <- kegg.rar.plot  + theme_bw() + 
   scale_color_manual(name = "Sample type",
                      breaks = c("Bulk_Soil", "Modern","Heritage"),
                      values = c("#004d40","#1e88e5","#d81b60"),
                      labels = c("Bulk soil", "Modern","Heritage")) +
   theme(axis.text.x = element_text(size = 8,angle = 35, hjust = 1)) +
   ylab("Gene richness")+
   ggtitle("b) Function")
 kegg.rar.plot2
 
 # Fig S1 final
 FigS1 <- ggarrange(taxa.rar.plot2+ rremove("xlab"),kegg.rar.plot2+ rremove("xlab"),
           common.legend = T, 
           legend = "top",
           align = "h",
           ncol = 2, nrow=1)
 
 annotate_figure(FigS1,
                 bottom = textGrob("Sequence sample size", gp = gpar(cex = 1)))
 
# Fig. S2 ####
 # Relative abundance plots at Phylum level
 glom <- tax_glom(psTAXA.rel, taxrank = 'Phylum', NArm=FALSE)
 data_glom<- psmelt(glom) # create dataframe from phyloseq object
 data_glom$Phylum <- as.character(data_glom$Phylum) #convert to character
 
 #simple way to rename phylum with < 1% abundance
 data_glom$Phylum[data_glom$Abundance < 0.01] <- "< 1% abund."
 
 # change NAs to Unclassified
 data_glom <- data_glom %>%
   mutate(Phylum = ifelse(is.na(Phylum), 'Unclassified', Phylum))
 
 # get the means for every Phylum 
 #first I need to sum the Unclassified and low abundances
 sumPhylum.df<- ddply(data_glom, c("Sample", "Phylum", "Sample_type"), summarise,
                      sum.Abund=sum(Abundance))
 #now get the means
 meanPhylum.df <- ddply(sumPhylum.df, c("Sample_type", "Phylum"), summarise,
                        mean.Abund=mean(sum.Abund),
                        sd.Abund=sd(sum.Abund))
 
 # Phylum order
 x = tapply(meanPhylum.df$mean.Abund, meanPhylum.df$Phylum, function(x) max(x))
 x = sort(x, TRUE)
 meanPhylum.df$Phylum = factor(as.character(meanPhylum.df$Phylum), levels=names(x))
 levels(meanPhylum.df$Phylum)
 
 #Count # phyla to set color palette
 Count = length(unique(data_glom$Phylum))
 Count

 #To change order of bacterial phylum, use factor() to change order of levels (for example, maybe you want to put the most abundant taxa first)... Use unique(data_glom$phylum) to see exact spelling of each phylum
 unique(data_glom$Phylum)
 data_glom$Phylum <- factor(data_glom$Phylum, levels = c("Proteobacteria", "Actinobacteria", "Acidobacteria", "Bacteroidetes",          
                                                         "Verrucomicrobia", "Chloroflexi", "Planctomycetes",         
                                                         "Gemmatimonadetes", "Firmicutes", "Candidatus Rokubacteria", "< 1% abund."))
 
 palette <- c("Proteobacteria"="#41AB5D",
              "Acidobacteria" = "darkorchid",
              "Actinobacteria"="#004529",
              "Bacteroidetes" = "#D94801",
              "Verrucomicrobia"="gold1", 
              "Chloroflexi" ="royalblue4",
              "Planctomycetes" ="deeppink",         
              "Gemmatimonadetes"="cyan2", 
              "Firmicutes"= "#67A9CF",
              "Candidatus Rokubacteria" = "#cc79a7",
              "< 1% abund."="#F7F7F7"
 )
 
 # New facet label names for Sample_type variable
 height.labs <- c("Bulk soil", "Modern", "Heritage")
 names(height.labs) <- c("Bulk_Soil", "Modern","Heritage")
 
 # set order of x axis
 unique(data_glom$Sample)
 sample_data(psTAXA)
 data_glom$Sample <- factor(data_glom$Sample, levels = c("S21", "S22", "S23", # bulk soil
                                                         "S6", "S7",  "S8",  "S9", "S10", # gallant
                                                         "S11", "S12", "S13", "S14", "S15", # Malacca
                                                         "S1", "S2", "S3", "S4", "S5", # Chidham
                                                         "S16", "S17", "S18", "S19", "S20"  # Red llamas
 ))
 
 # Create the plot
 spatial_plot <- ggplot(data=data_glom, aes(x=interaction(Sample,Variety, sep = "!"), y=Abundance, fill=Phylum)) + 
   facet_grid(~Sample_type, scales = "free",labeller = labeller(Sample_type = height.labs))
 spatial_plot + geom_bar(aes(), stat="identity", position="stack") +
   scale_fill_manual(name = "Phylum",
                     values = palette) +
   theme(legend.position="bottom",
         axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0),
         #axis.text.x = element_blank(),
         panel.background = element_blank()) + 
   guides(fill=guide_legend(nrow=3)) + 
   scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "") +
   ylab("Relative abundance") 

 
 
 
 # calculate means
 
 ddply(data_glom, c("Sample_type", "Phylum"), summarise,
       mean.Abund=mean(Abundance),
       sd.Abund=sd(Abundance))
 
 
 # Table S1 ###
 #### Top 20 most abundant genes in heritage vs bulk soil comparison
 # get the normalised counts for all samples from the original deseq object
 ddsKEGG_esf <- estimateSizeFactors(ddsKEGG)
 ddsKEGG_normcounts <- counts(ddsKEGG_esf, normalized=TRUE)
 OTU.norm = otu_table(ddsKEGG_normcounts, taxa_are_rows = TRUE)
 psKEGG_HvsBS <- psKEGG
 otu_table(psKEGG_HvsBS) <- OTU.norm
 
 # use the list of significant genes list.sig.HB  to prune the significant taxa
 psKEGG_HvsBS = prune_taxa(list.sig.HB, psKEGG_HvsBS)
 ps0<-names(sort(taxa_sums(psKEGG_HvsBS), TRUE)[1:20]) #get most abundant ones
 psKEGG_HvsBS_20<-prune_taxa(ps0, psKEGG_HvsBS)
 df_psKEGG_HvsBS_20<- psmelt(psKEGG_HvsBS_20) 
 
 
 
 # split into a list of df, one for each gene
 HvsBS_20List <- split(df_psKEGG_HvsBS_20,df_psKEGG_HvsBS_20$OTU)
 
 HBgenes_names <- c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                    "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 names(HvsBS_20List) <- HBgenes_names
 
 
 # run some tests on the gene abundances 
 
 
 # 1) check assumptions
 # 2) simple anova
 
 
 # 1) check assumptions
 
 #set up objects for the loop
 obs_norm_list=list(1:20)
 obs_var_list=list(1:20)
 obs_histo_list=list(1:20)
 
 # the loop
 for (i in c(1:20)) #The loop code will be repeated three times, based on the values that “I” refers to
 {
   df=HvsBS_20List[[i]] #Select the ith dataset in your list
   
   obs_norm_list[[i]] <-  shapiro.test(df$Abundance)
   obs_var_list[[i]] <-  leveneTest(Abundance ~ Sample_type, data = df)
   obs_histo_list[[i]] <- ggplot(df, aes(Abundance)) + geom_histogram() + ggtitle(names(HvsBS_20List)[i])
 }
 
 # check plots
 ggarrange(plotlist=obs_histo_list)#Plots all graphs in the list
 
 names(obs_norm_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                        "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 names(obs_var_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                       "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 
 # check results
 obs_norm_list # Significant = `12373` (0.027) - same by chance
 obs_var_list # Significant = `7552` (0.046), `3406` (0.048), `2529` (0.038), `15270` (0.037), `14266` (0.0311), `1190` (0.031)
 # variance can be handled 
 
 # I think this is ok to go ahead with the anova. They are pretty robust. 
 
 
 
 # 2) anova
 # empty list
 obs_lm_list = list(1:20)
 obs_anova_list = list(1:20)
 obs_emmeans_list = list(1:20)
 #loop
 for (i in c(1:20)) #The loop code will be repeated three times, based on the values that “I” refers to
 {
   df=HvsBS_20List[[i]] #Select the ith dataset in your list
   obs_lm_list[[i]] <- lm(Abundance ~ Sample_type, data = df)
   obs_anova_list[[i]] <- Anova(obs_lm_list[[i]], type="III")
   obs_emmeans_list[[i]] <- pairs(emmeans(obs_lm_list[[i]], "Sample_type"))
 }
 
 #check results
 names(obs_anova_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                         "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 names(obs_emmeans_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                           "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 obs_anova_list
 obs_emmeans_list
 
 
 # control for multiple testing
 # only need to correct for bulk soil vs sd because so many are boarder line.
 # The bulk soil vs tall are all significant as we would expect since they are the differentially abundant genes 
 # SD and tall only have one gene that is not significant and it is the one that makes sense looking at the data `7305` or K07305 peptide-methionine (R)-S-oxide reductase [EC:1.8.4.12]
 
 # make a list of p values for each contrast in the order of the names above
 BvM.pvalues <- c(0.0546,0.0125,0.0981,0.3092,0.0555,0.0134,0.0184,0.0260,0.0599,0.0223,
                  0.0186,0.0489,0.0332,0.0198,0.0455,0.0138,0.0014,0.0054,0.0178,0.0013)
 
 p.adjust(BvM.pvalues,"fdr")
 
 # [1] 0.06529412 0.03960000 0.10326316 0.30920000 0.06529412 0.03960000 0.03960000 0.04333333 0.06655556 0.04054545 
 # 0.03960000 0.06520000 0.05107692 0.03960000 0.06500000 0.03960000 0.01400000 0.03600000 0.03960000 0.01400000
 

#### Fig. S3 ####
 # first get the data in order
 library(plyr)
 # start with the contrast between Heritage and modern
 alpha = 0.05
 sigtabHM = HM.datTAXA2[which(HM.datTAXA2$svalue < alpha), ]
 sigtabHM$Sample_type <- ifelse(sigtabHM$log2FoldChange > 1, "Heritage", "Modern")
 
 #Get a count of each phylum
 sum.Phylum.dfHM<- ddply(sigtabHM, c("Phylum", "Sample_type"), summarise,
                         count=length(baseMean))
 
 sum.Phylum.dfHM$percentage <- (sum.Phylum.dfHM$count/1414)*100
 
 total <- ddply(sum.Phylum.dfHM, "Sample_type", summarise,
                total=sum(count))
 
 sum.Phylum.dfHM$percentSample_type <- ifelse(sum.Phylum.dfHM$Sample_type == "Heritage", (sum.Phylum.dfHM$count/1344)*100, (sum.Phylum.dfHM$count/70)*100)
 sum.Phylum.dfHM$Phylum <- as.character(sum.Phylum.dfHM$Phylum) #convert to character
 sum.Phylum.dfHM$Phylum[sum.Phylum.dfHM$percentSample_type < 5] <- "< 5% abund."
 sum.Phylum.dfHM2 <- aggregate(percentSample_type ~ Phylum + Sample_type, data=sum.Phylum.dfHM, FUN=sum)
 
 is.num <- sapply(sum.Phylum.dfHM2, is.numeric)
 sum.Phylum.dfHM2[is.num] <- lapply(sum.Phylum.dfHM2[is.num], round, 1)
 
 # Phylum order
 x = tapply(sum.Phylum.dfHM2$percentSample_type, sum.Phylum.dfHM2$Phylum, function(x) max(x))
 x = sort(x, TRUE)
 sum.Phylum.dfHM2$Phylum = factor(as.character(sum.Phylum.dfHM2$Phylum), levels=names(x))
 
 Count = length(unique(sum.Phylum.dfHM2$Phylum))
 Count
 
 # put the rare phyla at the bottom
 sum.Phylum.dfHM2$Phylum <- factor(sum.Phylum.dfHM2$Phylum, levels = c("Proteobacteria",  "Bacteroidetes", "Actinobacteria",
                                                                       "Firmicutes", "Spirochaetes", "Candidatus Thermoplasmatota", "< 5% abund."))
 
 
 #### look at bulk soil comparisons 
# Bulk soil vs heritage AND Bulk soil vs Modern
 sigtabHB = HB.datTAXA2[which(HB.datTAXA2$svalue < alpha), ]
 dim(sigtabHB)
 sigtabMB = MB.datTAXA2[which(MB.datTAXA2$svalue < alpha), ]
 dim(sigtabMB)

 # in both HB and MB its wheat vs bulk soil so positive is enriched in the plant compared with bulk soil
 sigtabHB$Sample_type <- ifelse(sigtabHB$log2FoldChange > 1, "Heritage", "Bulk_soil")
 sigtabMB$Sample_type <- ifelse(sigtabMB$log2FoldChange > 1, "Modern", "Bulk_soil")
 
 # count
 count(sigtabHB$log2FoldChange > 1)
 #  FALSE 1378   TRUE 3597
 count(sigtabMB$log2FoldChange > 1)
 #  FALSE  117  TRUE  997
 
 #first I need to get a count of each phylum
 sum.Phylum.dfHB<- ddply(sigtabHB, c("Phylum", "Sample_type"), summarise, count=length(baseMean))
 sum.Phylum.dfMB<- ddply(sigtabMB, c("Phylum", "Sample_type"), summarise, count=length(baseMean))
 
 #total number of taxa for calculating percentages
 rowsHB <- nrow(sigtabHB)
 rowsMB <- nrow(sigtabMB)
 
 # calculate percentages within sample type (Sample_type)
 totalHB <- ddply(sum.Phylum.dfHB, "Sample_type", summarise, total=sum(count))
 totalMB <- ddply(sum.Phylum.dfMB, "Sample_type", summarise, total=sum(count))
 
 sum.Phylum.dfHB$percentSample_type <- ifelse(sum.Phylum.dfHB$Sample_type == "Heritage", (sum.Phylum.dfHB$count/3597)*100, (sum.Phylum.dfHB$count/1378)*100)
 sum.Phylum.dfMB$percentSample_type <- ifelse(sum.Phylum.dfMB$Sample_type == "Modern", (sum.Phylum.dfMB$count/997)*100, (sum.Phylum.dfMB$count/117)*100)
 
 #convert to character
 sum.Phylum.dfHB$Phylum <- as.character(sum.Phylum.dfHB$Phylum) 
 sum.Phylum.dfMB$Phylum <- as.character(sum.Phylum.dfMB$Phylum)
 
 # group low abundant phyla
 sum.Phylum.dfHB$Phylum[sum.Phylum.dfHB$percentSample_type < 5] <- "< 5% abund."
 sum.Phylum.dfMB$Phylum[sum.Phylum.dfMB$percentSample_type < 5] <- "< 5% abund."
 
 sum.Phylum.df2HB <- aggregate(percentSample_type ~ Phylum + Sample_type, data=sum.Phylum.dfHB, FUN=sum)
 sum.Phylum.df2MB <- aggregate(percentSample_type ~ Phylum + Sample_type, data=sum.Phylum.dfMB, FUN=sum)
 
 is.numHB <- sapply(sum.Phylum.df2HB, is.numeric)
 is.numMB <- sapply(sum.Phylum.df2MB, is.numeric)
 
 sum.Phylum.df2HB[is.numHB] <- lapply(sum.Phylum.df2HB[is.numHB], round, 1)
 sum.Phylum.df2MB[is.numMB] <- lapply(sum.Phylum.df2MB[is.numMB], round, 1)
 
 # Phylum order HB
 x = tapply(sum.Phylum.df2HB$percentSample_type, sum.Phylum.df2HB$Phylum, function(x) max(x))
 x = sort(x, TRUE)
 sum.Phylum.df2HB$Phylum = factor(as.character(sum.Phylum.df2HB$Phylum), levels=names(x))
 Count = length(unique(sum.Phylum.df2HB$Phylum))
 Count #5
 
 # Phylum order MB
 x = tapply(sum.Phylum.df2MB$percentSample_type, sum.Phylum.df2MB$Phylum, function(x) max(x))
 x = sort(x, TRUE)
 sum.Phylum.df2MB$Phylum = factor(as.character(sum.Phylum.df2MB$Phylum), levels=names(x))
 Count = length(unique(sum.Phylum.df2MB$Phylum))
 Count #8
 
 
 # put the rare phyla at the bottom
 sum.Phylum.df2HB$Phylum <- factor(sum.Phylum.df2HB$Phylum, levels = c("Proteobacteria", "Firmicutes",  "Bacteroidetes", 
                                                                       "Actinobacteria", "< 5% abund."))
 
 sum.Phylum.df2MB$Phylum <- factor(sum.Phylum.df2MB$Phylum, levels = c("Proteobacteria",  "Bacteroidetes", "Firmicutes",
                                                                       "Thaumarchaeota", "Chlamydiae", "Candidatus Levybacteria",
                                                                       "Actinobacteria", "< 5% abund."))
 #### Now Plot
 # set colours
 phylumHB <- as.list(levels(sum.Phylum.df2HB$Phylum))
 phylumMB <- as.list(levels(sum.Phylum.df2MB$Phylum))
 phylumHM <- as.list(levels(sum.Phylum.dfHM2$Phylum))
 
 
 # create new list
 uniquephylum <- c(phylumHB,phylumMB,phylumHM)
 unique(uniquephylum)
 # set colours for all so that they are the same in each plot
 palette <- c("Firmicutes"= "#67A9CF",
              "Candidatus Levybacteria"="#ABD9E9",
              "Chlamydiae"= "#E6AB02",
              "Actinobacteria"="#004529",
              "< 5% abund."="#F7F7F7",
              "Thaumarchaeota"="#88419D",
              "Proteobacteria"="#41AB5D",
              "Bacteroidetes" = "#D94801",
              "Spirochaetes"="#FDE0DD",
              "Candidatus Thermoplasmatota"="#B2DF8A")
 
 # First, ensure all phyla are included in each dataset by adding rows with zero values
 # Get a complete list of all phyla across all datasets
 all_phyla <- unique(c(
   levels(sum.Phylum.df2HB$Phylum),
   levels(sum.Phylum.df2MB$Phylum),
   levels(sum.Phylum.dfHM2$Phylum)
 ))
 
 # Function to add missing phyla with zero values to each dataset
 add_missing_phyla <- function(df, all_phyla) {
   existing_phyla <- levels(df$Phylum)
   missing_phyla <- setdiff(all_phyla, existing_phyla)
   
   if (length(missing_phyla) > 0) {
     # Get sample types
     sample_types <- unique(df$Sample_type)
     
     # Create rows for missing phyla
     new_rows <- data.frame()
     for (st in sample_types) {
       for (ph in missing_phyla) {
         new_row <- data.frame(
           Sample_type = st,
           Phylum = ph,
           percentSample_type = 0  # Zero abundance
         )
         new_rows <- rbind(new_rows, new_row)
       }
     }
     
     # Make sure Phylum is a factor with the same levels
     new_rows$Phylum <- factor(new_rows$Phylum, levels = all_phyla)
     df$Phylum <- factor(df$Phylum, levels = all_phyla)
     
     # Combine with original data
     result <- rbind(df, new_rows)
     return(result)
   }
   
   # If no missing phyla, just update the factor levels
   df$Phylum <- factor(df$Phylum, levels = all_phyla)
   return(df)
 }
 
 # Apply the function to each dataset
 sum.Phylum.dfHM2_complete <- add_missing_phyla(sum.Phylum.dfHM2, all_phyla)
 sum.Phylum.df2HB_complete <- add_missing_phyla(sum.Phylum.df2HB, all_phyla)
 sum.Phylum.df2MB_complete <- add_missing_phyla(sum.Phylum.df2MB, all_phyla)
 
 # Now remake the plots with the complete datasets
 
 # plot HM
 phyla_plotHM <- ggplot(data=sum.Phylum.dfHM2_complete, aes(x=Sample_type, y=percentSample_type, fill = factor(Phylum, levels = c("Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", "Spirochaetes","Thaumarchaeota",
                                                                                                                                  "Chlamydiae","Candidatus Levybacteria","Candidatus Thermoplasmatota","< 5% abund.")), 
                                                            label = percentSample_type))
 
 phyla_plotHM2 <- phyla_plotHM + geom_bar(aes(), stat="identity", position="stack") +
   scale_fill_manual(name = "Phylum",
                     values = palette,
                     drop = FALSE) +
   theme(legend.position="right",
         #axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0),
         #axis.text.x = element_blank(),
         axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
         panel.background = element_blank()) + 
   guides(fill=guide_legend(ncol=4)) + 
   geom_text(size = 3, position = position_stack(vjust = 0.5),check_overlap = TRUE) +
   #scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "") +
   ylab("Proportion of differentially abundant taxa") 
 
 # then plot HB
 
 phyla_plotHB <- ggplot(data=sum.Phylum.df2HB_complete, aes(x=Sample_type, y=percentSample_type, fill = factor(Phylum, levels = c("Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", "Spirochaetes","Thaumarchaeota",
                                                                                                                                  "Chlamydiae","Candidatus Levybacteria","Candidatus Thermoplasmatota","< 5% abund.")), 
                                                            label = percentSample_type))
 
 phyla_plotHB2 <- phyla_plotHB + geom_bar(aes(), stat="identity", position="stack") +
   scale_fill_manual(name = "Phylum",
                     values = palette,
                     drop = FALSE) +
   theme(legend.position="right",
         #axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0),
         #axis.text.x = element_blank(),
         axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
         panel.background = element_blank()) + 
   guides(fill=guide_legend(ncol=4)) + 
   geom_text(size = 3, position = position_stack(vjust = 0.5),check_overlap = TRUE) +
   #scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "") +
   ylab("Proportion of differentially abundant taxa") 
 
 
 # then plot MB
 phyla_plotMB <- ggplot(data=sum.Phylum.df2MB_complete, aes(x=Sample_type, y=percentSample_type, fill = factor(Phylum, levels = c("Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", "Spirochaetes","Thaumarchaeota",
                                                                                                                                  "Chlamydiae","Candidatus Levybacteria","Candidatus Thermoplasmatota","< 5% abund.")), 
                                                            label = percentSample_type))
 
 phyla_plotMB2 <- phyla_plotMB + geom_bar(aes(), stat="identity", position="stack") +
   scale_fill_manual(name = "Phylum",
                     values = palette,
                     drop = FALSE) +
   theme(legend.position="right",
         #axis.text.x = element_text(angle = 300, vjust = 0.5, hjust=0),
         #axis.text.x = element_blank(),
         axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
         panel.background = element_blank()) + 
   guides(fill=guide_legend(ncol=4)) + 
   geom_text(size = 3, position = position_stack(vjust = 0.5),check_overlap = TRUE) +
   #scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "") +
   ylab("Proportion of differentially abundant taxa") 
 
 
 library(ggpubr)
 ggarrange(phyla_plotHM2, phyla_plotHB2+ rremove("ylab"), phyla_plotMB2+ rremove("ylab"),
           align = "h",
           common.legend = TRUE,
           legend = "bottom",
           ncol = 3, nrow=1)
 
 #### Table S1 ####
 # Top 20 most abundant genes in tall vs bulk soil comparison 
 # get the normalised counts for all samples from the original deseq object
 ddsKEGG_esf <- estimateSizeFactors(ddsKEGG)
 ddsKEGG_normcounts <- counts(ddsKEGG_esf, normalized=TRUE)
 OTU.norm = otu_table(ddsKEGG_normcounts, taxa_are_rows = TRUE)
 psKEGG_HvsBS <- psKEGG
 otu_table(psKEGG_HvsBS) <- OTU.norm
 # use the list of significant genes list.sig.HB  to prune the significant taxa
 psKEGG_HvsBS = prune_taxa(list.sig.HB, psKEGG_HvsBS)
 ps0<-names(sort(taxa_sums(psKEGG_HvsBS), TRUE)[1:20]) #get most abundant ones
 psKEGG_HvsBS_20<-prune_taxa(ps0, psKEGG_HvsBS)
 psKEGG_HvsBS_20DF<- psmelt(psKEGG_HvsBS_20) 
 
 # split into a list of df, one for each gene
 HvsBS_20List <- split(psKEGG_HvsBS_20DF,psKEGG_HvsBS_20DF$OTU)
 length(HvsBS_20List)
 nrow(HvsBS_20List[[2]])
 head(HvsBS_20List[[20]])
 HBgenes_names <- c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                    "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 names(HvsBS_20List) <- HBgenes_names
 
 # run some tests on the gene abundances 
 # 1) check assumptions
 #set up objects for the loop
 obs_norm_list=list(1:20)
 obs_var_list=list(1:20)
 obs_histo_list=list(1:20)
 
 # the loop
 for (i in c(1:20)) #The loop code will be repeated three times, based on the values that “I” refers to
 {
   df=HvsBS_20List[[i]] #Select the ith dataset in your list
   
   obs_norm_list[[i]] <-  shapiro.test(df$Abundance)
   obs_var_list[[i]] <-  leveneTest(Abundance ~ Sample_type, data = df)
   obs_histo_list[[i]] <- ggplot(df, aes(Abundance)) + geom_histogram() + ggtitle(names(HvsBS_20List)[i])
 }
 
 # check plots
 ggarrange(plotlist=obs_histo_list)#Plots all graphs in the list
 
 names(obs_norm_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                        "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 names(obs_var_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                       "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 
 # check results
 obs_norm_list # Significant = `12373` (0.027) - same by chance
 obs_var_list # Significant = `7552` (0.046), `3406` (0.048), `2529` (0.038), `15270` (0.037), `14266` (0.0311), `1190` (0.031)
 # variance can be handled 
 
 # I think this is ok to go ahead with the anova. They are pretty robust. 
  # 2) anova
 # empty list
 obs_lm_list = list(1:20)
 obs_anova_list = list(1:20)
 obs_emmeans_list = list(1:20)
 #loop
 for (i in c(1:20)) #The loop code will be repeated three times, based on the values that “I” refers to
 {
   df=HvsBS_20List[[i]] #Select the ith dataset in your list
   obs_lm_list[[i]] <- lm(Abundance ~ Sample_type, data = df)
   obs_anova_list[[i]] <- Anova(obs_lm_list[[i]], type="III")
   obs_emmeans_list[[i]] <- pairs(emmeans(obs_lm_list[[i]], "Sample_type"))
 }
 
 #check results
 names(obs_anova_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                         "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 names(obs_emmeans_list)=c("1190", "1198", "12308", "12373", "14266","15270", "15923", "16135", "1811","18900",
                           "2014", "21572","21573","2529", "3406","3566","362","5970","7305", "7552")
 obs_anova_list
 obs_emmeans_list
 
 # control for multiple testing
 # only need to correct for bulk soil vs sd because so many are boarder line.
 # The bulk soil vs tall are all significant as we would expect since they are the differentially abundant genes 
 # SD and tall only have one gene that is not significant and it is the one that makes sense looking at the data `7305` or K07305 peptide-methionine (R)-S-oxide reductase [EC:1.8.4.12]
 
 # make a list of p values for each contrast in the order of the names above
 BvsD.pvalues <- c(0.0546,0.0125,0.0981,0.3092,0.0555,0.0134,0.0184,0.0260,0.0599,0.0223,
                   0.0186,0.0489,0.0332,0.0198,0.0455,0.0138,0.0014,0.0054,0.0178,0.0013)
 
 p.adjust(BvsD.pvalues,"fdr")

#### Table S2 #### 
 # see Loop_script.R

#### END ####