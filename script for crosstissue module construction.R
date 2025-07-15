setwd('')
#install.packages('BiocManager')
#BiocManager::install('bnstruct')
library(bnstruct)
library(WGCNA)
allowWGCNAThreads()
library(reshape2)
library(qgraph)
library(reshape2)
library(dplyr)
adip = read.delim('HMDP_chow_trx_adipose.txt')
liver = read.delim('HMDP_chow_trx_liver.txt')
aor = read.delim('HMDP_chow_trx_aorta.txt')
heart = read.delim('HMDP_chow_trx_heart.txt')
bone = read.delim('HMDP_chow_trx_bone.txt')
traits = read.delim('HMDP_chow_traits.txt')

adip = adip[adip$Strain %in% heart$Strain,]
liver = liver[liver$Strain %in% heart$Strain,]
aor = aor[aor$Strain %in% heart$Strain,]
heart = heart[heart$Strain %in% adip$Strain,]
bone = bone[bone$Strain %in% heart$Strain,]

adip$gene = paste0(adip$gene_symbol, '_adipose')
liver$gene = paste0(liver$gene_symbol, '_liver')
aor$gene = paste0(aor$gene_symbol, '_aorta')
heart$gene = paste0(heart$gene_symbol, '_heart')
bone$gene = paste0(bone$Symbol, '_bone')

a1 = adip %>% select(Strain, gene, expression_value)
l1 = liver %>% select(Strain, gene, expression_value)
ao1 = aor %>% select(Strain, gene, expression_value)
h1 = heart %>% select(Strain, gene, expression_value)
b1 = bone %>% select(Strain, gene, expression_value)

all_tissues = as.data.frame(rbind(a1, l1, ao1, h1, b1))
bin_mat = dcast(all_tissues, Strain ~ gene, value.var = 'expression_value', fun.aggregate = mean)

row.names(bin_mat) = bin_mat$Strain
bin_mat$Strain = NULL


trait = traits[traits$Strain %in% heart$Strain,]
t1 = dcast(trait, Strain ~ trait_name, value.var = 'value', fun.aggregate = mean )
row.names(t1) = t1$Strain
trait_match = t1 %>% select(Strain, 'Percent Fat Mass' = NMR_BF_percentage, 'Body Weight' = 'BW', 'Glucose' = 'Glucose', 'Plasma TG' = 'triglycerides', 'Plasma Cholesterol' = 'total cholesterol', 'Plasma FFA' = 'free fatty acids', 'Insulin' = 'Insulin')


row.names(trait_match) = trait_match$Strain
trait_match$Strain = NULL


femData = bin_mat
datExpr0 = femData
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#FALSE
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 111, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 111, minSize = 4)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

ms_screen = as.data.frame(row.names(datExpr0))
ms_screen$kept = paste0(clust)
ms_screen = ms_screen[grepl("1", ms_screen$kept),]

final_chow_traits = trait_match
final_chow_traits$m = match(row.names(final_chow_traits), ms_screen$`row.names(datExpr0)`, nomatch=0)
final_chow_traits$m = final_chow_traits$m > 0
final_chow_traits = final_chow_traits[!grepl('FALSE', final_chow_traits$m),]
final_chow_traits$m = NULL

datTraits = final_chow_traits
collectGarbage()
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry

traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

#didnt work


powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
datExpr
#net = blockwiseModules(datExpr, power = 3, maxBlockSize = 5000,
# TOMType = "unsigned", minModuleSize = 100,
# reassignThreshold = 0, mergeCutHeight = 0.35,
# numericLabels = TRUE, pamRespectsDendro = FALSE,
# saveTOMs = TRUE,
# saveTOMFileBase = "femaleMouseTOM",
# verbose = 3)

net = blockwiseModules(datExpr, power = 2, maxBlockSize = 5000,
                       TOMType = "unsigned", minModuleSize = 125,
                       reassignThreshold = 1e-3, mergeCutHeight = 0.42,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3)

table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

length(net$colors)

test1 = as.data.frame(row.names(datExpr0))
test1$keep = keepSamples
test1 = test1[!grepl("FALSE", test1$keep),]

all_tog = MEs
all_tog$Strain = row.names(all_tog)
final_chow_traits$Strain = row.names(final_chow_traits)
all_tog1 = inner_join(all_tog, final_chow_traits, by = 'Strain')

module_membership = as.data.frame(net$colors)
colnames(module_membership) = 'module'
module_membership$gene = row.names(module_membership)
#all_tog1$ME0 = NULL
#install.packages('bnstruct')
#library("devtools")
#install_github("sambofra/bnstruct")

row.names(all_tog1) = all_tog1$Strain
all_tog1$Strain = NULL
all_tog1$ME0 = NULL
bics_map = bicorAndPvalue(all_tog1, all_tog1, use = 'p')

map1 = as.data.frame(bics_map$bicor)
map1 = melt(as.matrix(map1))
map1$value[map1$value > 0.999999] <- 0
map2 = dcast(map1, Var1 ~ Var2, value.var = 'value')
row.names(map2) = map2$Var1
map2$Var1 = NULL
#qgraph(map2, minimum = 0.3, cut = 0.5, vsize = 3, legend = TRUE, borders = FALSE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3")

qgraph(map2, minimum = 0.25, cut = 0.5, vsize = 4, legend = F, borders = FALSE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=1, directed=F, labels = colnames(map2))

module_membership$tissue = gsub(".*_","",module_membership$gene)
module_membership$gene_symbol = gsub("_.*","",module_membership$gene)


#example for writing lists of genes to tables
mod_table = module_membership[module_membership$module==1,]
write.table(mod_table$gene_symbol, file = 'ME1 module members.txt', row.names = F, quote = F)

mod_table = module_membership[module_membership$module==2,]
write.table(mod_table$gene_symbol, file = 'M2 module members.txt', row.names = F, quote = F)

