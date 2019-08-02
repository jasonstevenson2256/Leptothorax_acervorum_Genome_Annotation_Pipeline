#### Testing for functional annotation enrichment in the social region ####

# This script is adapted from one written by Alun Junes at UoL (GGB)
# The original script was designed to look for GO term enrichment
# in differentially expressed genes.  This new script looks for
# GO term enrichment in functional annotations of predicted genes
# in the LG2 social region of L. acervorum versus the complete set of
# functional annotations for all genes in the annotated genomes.
# There are two annotated genomes, one for each social phenotype
# of L. acervorum.

# 1. Load packages

#Only run this line once!
install.packages("BiocManager")
BiocManager::install(c("GOstats", "GSEABase", "treemap"))

#Run this block every time you run this script
library(GOstats)
library(GSEABase)
library(treemap)


# 2. Setwd and read in data.  
# The GO_annotations_OT and GO_annotations_PF are two-column tsv
# lists of gene IDs (left column) and GO terms (right column).
# All entries are in quotes and only one GO term can be
# associated with one ID in one row (i.e. gene IDs are repeated
# as many times as there are associated GO terms).  The GO terms
# were assigned using the GOfeat web server.  
# The "evi" column contains the GO evidence code indicating
# how GO terms were assigned.  "IEA" stands for "Inferred from
# Electronic Annotation" 
# http://wiki.geneontology.org/index.php/Inferred_from_Electronic_Annotation_(IEA)

setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\GO FEAT")

GO_annotations_OT <- read.table("OT3b_IDs_and_terms.txt", 
                                header = FALSE, 
                                stringsAsFactors = FALSE)
GO_annotations_OT[,3] <- "IEA"
names(GO_annotations_OT) <- c("genes","GOIds","evi")
GO_annotations_OT <- GO_annotations_OT[c(2,3,1)]

GO_annotations_PF <- read.table("PF_IDs_and_terms.txt", 
                                header = FALSE, 
                                stringsAsFactors = FALSE)
GO_annotations_PF[,3] <- "IEA"
names(GO_annotations_PF) <- c("genes","GOIds","evi")
GO_annotations_PF <- GO_annotations_PF[c(2,3,1)]


# 3. For each set of annotations, create a GOFrame object,
# a gene set collection, and a gene universe
# some warnings might be thrown here, it's fine
GO_frame_OT <- GOFrame(GO_annotations_OT, organism = "Leptothorax acervorum")
goAllFrame_OT <- GOAllFrame(GO_frame_OT)
gsc_OT <- GeneSetCollection(goAllFrame_OT, setType = GOCollection())
universe_OT <- as.vector(unique(GO_annotations_OT[,3]))

GO_frame_PF <- GOFrame(GO_annotations_PF, organism = "Leptothorax acervorum")
goAllFrame_PF <- GOAllFrame(GO_frame_PF)
gsc_PF <- GeneSetCollection(goAllFrame_PF, setType = GOCollection())
universe_PF <- as.vector(unique(GO_annotations_PF[,3]))


# 4. Read in an wrangle the lists of genes we are interested in.
# For each assembly (OT and PF) there are three sets of genes
# which are of interest to compare to the genome as a whole.
# These are:
# a) all genes on scaffolds containing social markers
# b) all genes in scaffolds assigned to LG2 
# c) all genes in scaffolds with social markers which are also assigned to LG2
# These lists are read in as single-columns of character strings

#Read OT gene lists
OT_genes_a <- read.table(file = "Gene_ID_lists\\OT3b_Social_Marker_gene_IDs.txt", 
                         header = FALSE, stringsAsFactors = FALSE)
OT_genes_b <- read.table(file = "Gene_ID_lists\\OT3b_LG2_gene_IDs.txt", 
                         header = FALSE, stringsAsFactors = FALSE)
OT_genes_c <- read.table(file = "Gene_ID_lists\\OT3b_LG2_Social_Marker_gene_IDs.txt", 
                         header = FALSE, stringsAsFactors = FALSE)
#Read PF gene lists
PF_genes_a <- read.table(file = "Gene_ID_lists\\PF_Social_Marker_gene_IDs.txt", 
                         header = FALSE, stringsAsFactors = FALSE)
PF_genes_b <- read.table(file = "Gene_ID_lists\\PF_LG2_gene_IDs.txt", 
                         header = FALSE, stringsAsFactors = FALSE)
PF_genes_c <- read.table(file = "Gene_ID_lists\\PF_LG2_Social_Marker_gene_IDs.txt", 
                         header = FALSE, stringsAsFactors = FALSE)

#Coerce all gene lists to vectors
OT_genes_a <- as.vector(OT_genes_a[,1])
OT_genes_b <- as.vector(OT_genes_b[,1])
OT_genes_c <- as.vector(OT_genes_c[,1])
PF_genes_a <- as.vector(PF_genes_a[,1])
PF_genes_b <- as.vector(PF_genes_b[,1])
PF_genes_c <- as.vector(PF_genes_c[,1])

#Keep only genes which have been annotated
OT_genes_a <- OT_genes_a[OT_genes_a %in% universe_OT]
OT_genes_b <- OT_genes_b[OT_genes_b %in% universe_OT]
OT_genes_c <- OT_genes_c[OT_genes_c %in% universe_OT]
PF_genes_a <- PF_genes_a[PF_genes_a %in% universe_PF]
PF_genes_b <- PF_genes_b[PF_genes_b %in% universe_PF]
PF_genes_c <- PF_genes_c[PF_genes_c %in% universe_PF]


# 5. Generate parameters for the hypergeometric test

#First, a function
#The arguments are as follows:
#genes_of_i = one of the gene lists from part 4
#universe = one of the two universes of genes from part 3
#pvalue_cut = p-value cutoff (e.g. 0.05)
#gscollection = one of the two gsc objects from part 3
Get_GO_params_all <- function(genes_of_i, universe, pvalue_cut, gscollection){
  onto_terms <- c("BP", "CC", "MF")
  directions <- c("over", "under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1, paste(onto_terms[i], directions[j], sep = "_"))
      parameters <- GSEAGOHyperGParams(name = "Leptothorax acervorum Hygeo params", 
                                       geneSetCollection = gscollection, 
                                       universeGeneIds = universe, 
                                       geneIds = genes_of_i, 
                                       ontology = paste(onto_terms[i]), 
                                       pvalueCutoff = pvalue_cut, 
                                       conditional = TRUE, 
                                       testDirection = paste(directions[j]))
      param_list <- c(param_list, parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

#Then call the function for OT
param_list_OT_a <- Get_GO_params_all(genes_of_i = OT_genes_a, 
                                     universe = universe_OT, 
                                     pvalue_cut = 0.05, 
                                     gscollection = gsc_OT)
param_list_OT_b <- Get_GO_params_all(genes_of_i = OT_genes_b, 
                                     universe = universe_OT, 
                                     pvalue_cut = 0.05, 
                                     gscollection = gsc_OT)
param_list_OT_c <- Get_GO_params_all(genes_of_i = OT_genes_c, 
                                     universe = universe_OT, 
                                     pvalue_cut = 0.05, 
                                     gscollection = gsc_OT)
#Then for PF
param_list_PF_a <- Get_GO_params_all(genes_of_i = PF_genes_a, 
                                     universe = universe_PF, 
                                     pvalue_cut = 0.05, 
                                     gscollection = gsc_PF)
param_list_PF_b <- Get_GO_params_all(genes_of_i = PF_genes_b, 
                                     universe = universe_PF, 
                                     pvalue_cut = 0.05, 
                                     gscollection = gsc_PF)
param_list_PF_c <- Get_GO_params_all(genes_of_i = PF_genes_c, 
                                     universe = universe_PF, 
                                     pvalue_cut = 0.05, 
                                     gscollection = gsc_PF)


# 6. Conduct hypergeometric tests using the
# parameter lists created in part 5

#First a function
Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list, res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}

#Then call the function
GO_enrichment_OT_a <- Hyper_G_test(param_list = param_list_OT_a)
GO_enrichment_OT_b <- Hyper_G_test(param_list = param_list_OT_b)
GO_enrichment_OT_c <- Hyper_G_test(param_list = param_list_OT_c)
GO_enrichment_PF_a <- Hyper_G_test(param_list = param_list_PF_a)
GO_enrichment_PF_b <- Hyper_G_test(param_list = param_list_PF_b)
GO_enrichment_PF_c <- Hyper_G_test(param_list = param_list_PF_c)


# 7. Pull out results using summary()

#Over-represented GO terms
#First pull out summaries of the outputs from part 6 for OT
Result_OT_a_MF_over <- summary(GO_enrichment_OT_a[["MF_over"]])
Result_OT_a_CC_over <- summary(GO_enrichment_OT_a[["CC_over"]])
Result_OT_a_BP_over <- summary(GO_enrichment_OT_a[["BP_over"]])
Result_OT_b_MF_over <- summary(GO_enrichment_OT_b[["MF_over"]])
Result_OT_b_CC_over <- summary(GO_enrichment_OT_b[["CC_over"]])
Result_OT_b_BP_over <- summary(GO_enrichment_OT_b[["BP_over"]])
Result_OT_c_MF_over <- summary(GO_enrichment_OT_c[["MF_over"]])
Result_OT_c_CC_over <- summary(GO_enrichment_OT_c[["CC_over"]])
Result_OT_c_BP_over <- summary(GO_enrichment_OT_c[["BP_over"]])
#Then for PF
Result_PF_a_MF_over <- summary(GO_enrichment_PF_a[["MF_over"]])
Result_PF_a_CC_over <- summary(GO_enrichment_PF_a[["CC_over"]])
Result_PF_a_BP_over <- summary(GO_enrichment_PF_a[["BP_over"]])
Result_PF_b_MF_over <- summary(GO_enrichment_PF_b[["MF_over"]])
Result_PF_b_CC_over <- summary(GO_enrichment_PF_b[["CC_over"]])
Result_PF_b_BP_over <- summary(GO_enrichment_PF_b[["BP_over"]])
Result_PF_c_MF_over <- summary(GO_enrichment_PF_c[["MF_over"]])
Result_PF_c_CC_over <- summary(GO_enrichment_PF_c[["CC_over"]])
Result_PF_c_BP_over <- summary(GO_enrichment_PF_c[["BP_over"]])

#Under-represented GO terms
#First pull out summaries of the outputs from part 6 for OT
Result_OT_a_MF_under <- summary(GO_enrichment_OT_a[["MF_under"]])
Result_OT_a_CC_under <- summary(GO_enrichment_OT_a[["CC_under"]])
Result_OT_a_BP_under <- summary(GO_enrichment_OT_a[["BP_under"]])
Result_OT_b_MF_under <- summary(GO_enrichment_OT_b[["MF_under"]])
Result_OT_b_CC_under <- summary(GO_enrichment_OT_b[["CC_under"]])
Result_OT_b_BP_under <- summary(GO_enrichment_OT_b[["BP_under"]])
Result_OT_c_MF_under <- summary(GO_enrichment_OT_c[["MF_under"]])
Result_OT_c_CC_under <- summary(GO_enrichment_OT_c[["CC_under"]])
Result_OT_c_BP_under <- summary(GO_enrichment_OT_c[["BP_under"]])
#Then for PF
Result_PF_a_MF_under <- summary(GO_enrichment_PF_a[["MF_under"]])
Result_PF_a_CC_under <- summary(GO_enrichment_PF_a[["CC_under"]])
Result_PF_a_BP_under <- summary(GO_enrichment_PF_a[["BP_under"]])
Result_PF_b_MF_under <- summary(GO_enrichment_PF_b[["MF_under"]])
Result_PF_b_CC_under <- summary(GO_enrichment_PF_b[["CC_under"]])
Result_PF_b_BP_under <- summary(GO_enrichment_PF_b[["BP_under"]])
Result_PF_c_MF_under <- summary(GO_enrichment_PF_c[["MF_under"]])
Result_PF_c_CC_under <- summary(GO_enrichment_PF_c[["CC_under"]])
Result_PF_c_BP_under <- summary(GO_enrichment_PF_c[["BP_under"]])


# 8. P-value correction
# False-discovery rate p-value correction is used because
# alternatives (e.g. Bonferroni) are typically too strict
# with exploratory analyses such as this, and are liable
# to produce false negatives.

#Over-represented GO terms
#OT
Result_OT_a_MF_over <- Result_OT_a_MF_over[p.adjust(Result_OT_a_MF_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_a_CC_over <- Result_OT_a_CC_over[p.adjust(Result_OT_a_CC_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_a_BP_over <- Result_OT_a_BP_over[p.adjust(Result_OT_a_BP_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_b_MF_over <- Result_OT_b_MF_over[p.adjust(Result_OT_b_MF_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_b_CC_over <- Result_OT_b_CC_over[p.adjust(Result_OT_b_CC_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_b_BP_over <- Result_OT_b_BP_over[p.adjust(Result_OT_b_BP_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_c_MF_over <- Result_OT_c_MF_over[p.adjust(Result_OT_c_MF_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_c_CC_over <- Result_OT_c_CC_over[p.adjust(Result_OT_c_CC_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_c_BP_over <- Result_OT_c_BP_over[p.adjust(Result_OT_c_BP_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
#PF
Result_PF_a_MF_over <- Result_PF_a_MF_over[p.adjust(Result_PF_a_MF_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_a_CC_over <- Result_PF_a_CC_over[p.adjust(Result_PF_a_CC_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_a_BP_over <- Result_PF_a_BP_over[p.adjust(Result_PF_a_BP_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_b_MF_over <- Result_PF_b_MF_over[p.adjust(Result_PF_b_MF_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_b_CC_over <- Result_PF_b_CC_over[p.adjust(Result_PF_b_CC_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_b_BP_over <- Result_PF_b_BP_over[p.adjust(Result_PF_b_BP_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_c_MF_over <- Result_PF_c_MF_over[p.adjust(Result_PF_c_MF_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_c_CC_over <- Result_PF_c_CC_over[p.adjust(Result_PF_c_CC_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_c_BP_over <- Result_PF_c_BP_over[p.adjust(Result_PF_c_BP_over$Pvalue, 
                                                    method = "fdr") < 0.05,]
#Under-represented GO terms
#OT
Result_OT_a_MF_under <- Result_OT_a_MF_under[p.adjust(Result_OT_a_MF_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_a_CC_under <- Result_OT_a_CC_under[p.adjust(Result_OT_a_CC_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_a_BP_under <- Result_OT_a_BP_under[p.adjust(Result_OT_a_BP_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_b_MF_under <- Result_OT_b_MF_under[p.adjust(Result_OT_b_MF_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_b_CC_under <- Result_OT_b_CC_under[p.adjust(Result_OT_b_CC_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_b_BP_under <- Result_OT_b_BP_under[p.adjust(Result_OT_b_BP_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_c_MF_under <- Result_OT_c_MF_under[p.adjust(Result_OT_c_MF_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_c_CC_under <- Result_OT_c_CC_under[p.adjust(Result_OT_c_CC_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_OT_c_BP_under <- Result_OT_c_BP_under[p.adjust(Result_OT_c_BP_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
#PF
Result_PF_a_MF_under <- Result_PF_a_MF_under[p.adjust(Result_PF_a_MF_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_a_CC_under <- Result_PF_a_CC_under[p.adjust(Result_PF_a_CC_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_a_BP_under <- Result_PF_a_BP_under[p.adjust(Result_PF_a_BP_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_b_MF_under <- Result_PF_b_MF_under[p.adjust(Result_PF_b_MF_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_b_CC_under <- Result_PF_b_CC_under[p.adjust(Result_PF_b_CC_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_b_BP_under <- Result_PF_b_BP_under[p.adjust(Result_PF_b_BP_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_c_MF_under <- Result_PF_c_MF_under[p.adjust(Result_PF_c_MF_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_c_CC_under <- Result_PF_c_CC_under[p.adjust(Result_PF_c_CC_under$Pvalue, 
                                                    method = "fdr") < 0.05,]
Result_PF_c_BP_under <- Result_PF_c_BP_under[p.adjust(Result_PF_c_BP_under$Pvalue, 
                                                    method = "fdr") < 0.05,]


# 9. Write files for uploading to revigo website
#You may also want to change quote = TRUE to quote = FALSE if 
#revigo doesn't like the files.

#setwd
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\GO FEAT\\Funtional_Enrichment")

#OT
write.table(Result_OT_a_MF_over[,1:2], file="OT_a_MF_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_a_CC_over[,1:2], file="OT_a_CC_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_a_BP_over[,1:2], file="OT_a_BP_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_b_MF_over[,1:2], file="OT_b_MF_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_b_CC_over[,1:2], file="OT_b_CC_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_b_BP_over[,1:2], file="OT_b_BP_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_c_MF_over[,1:2], file="OT_c_MF_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_c_CC_over[,1:2], file="OT_c_CC_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_c_BP_over[,1:2], file="OT_c_BP_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(Result_OT_a_MF_under[,1:2], file="OT_a_MF_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_a_CC_under[,1:2], file="OT_a_CC_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_a_BP_under[,1:2], file="OT_a_BP_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_b_MF_under[,1:2], file="OT_b_MF_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_b_CC_under[,1:2], file="OT_b_CC_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_b_BP_under[,1:2], file="OT_b_BP_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_c_MF_under[,1:2], file="OT_c_MF_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_c_CC_under[,1:2], file="OT_c_CC_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_OT_c_BP_under[,1:2], file="OT_c_BP_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#PF
write.table(Result_PF_a_MF_over[,1:2], file="PF_a_MF_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_a_CC_over[,1:2], file="PF_a_CC_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_a_BP_over[,1:2], file="PF_a_BP_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_b_MF_over[,1:2], file="PF_b_MF_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_b_CC_over[,1:2], file="PF_b_CC_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_b_BP_over[,1:2], file="PF_b_BP_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_c_MF_over[,1:2], file="PF_c_MF_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_c_CC_over[,1:2], file="PF_c_CC_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_c_BP_over[,1:2], file="PF_c_BP_Over.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(Result_PF_a_MF_under[,1:2], file="PF_a_MF_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_a_CC_under[,1:2], file="PF_a_CC_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_a_BP_under[,1:2], file="PF_a_BP_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_b_MF_under[,1:2], file="PF_b_MF_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_b_CC_under[,1:2], file="PF_b_CC_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_b_BP_under[,1:2], file="PF_b_BP_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_c_MF_under[,1:2], file="PF_c_MF_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_c_CC_under[,1:2], file="PF_c_CC_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Result_PF_c_BP_under[,1:2], file="PF_c_BP_Under.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)