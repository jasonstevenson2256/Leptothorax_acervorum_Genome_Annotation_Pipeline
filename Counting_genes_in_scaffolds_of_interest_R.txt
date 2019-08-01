## Counting genes in scaffolds of interest

# 2x annotations
# 3x scaffold lists of interest
# Lists in single column format
# Annotations in gff3 format

# Per annotation, split gff3 into 3 sets (1 per list of IDs)
# Subset down to records saying "gene"
# Count rows


# 1. Setwd and read gff3 file
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\BRAKER\\FM")
OT3b_GFF <- read.table("augustus.good.genes.gff3.txt", header=FALSE, sep = "\t")
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\BRAKER\\P")
PF_GFF <- read.table("augustus.good.genes.gff3.txt", header=FALSE, sep = "\t")

# 2. Assign meaningful column names
column_names <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
colnames(OT3b_GFF) <- column_names
colnames(PF_GFF) <- column_names

# 3. Read lists of scaffold IDs
# _all = all contigs with markers
# _LG2 = all contigs in LG2
# _LG2_markers = all LG2 contigs with markers
# _all_contigs = every single contig in assembly
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\10x Genomes\\L_acer_FM")
OT3b_all <- read.table("IDs_OT18_12_MP4_re-do_S2_assembly_3b_all_contigs_containing_social_markers.txt", header = FALSE)
OT3b_LG2 <- read.table("IDs_OT18_12_MP4_re-do_S2_assembly_3b_all_contigs_in_LG2.txt", header = FALSE)
OT3b_LG2_markers <- read.table("IDs_OT18_12_MP4_re-do_S2_assembly_3b_all_contigs_in_LG2_containing_social_markers.txt", header = FALSE)
OT3b_all_conts <- read.table("OT_conts.txt", stringsAsFactors = FALSE)

setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\10x Genomes\\L_acer_P")
PF_all <- read.table("IDs_PF_15_M1_S3_all_contigs_containing_social_markers.txt", header = FALSE)
PF_LG2 <- read.table("IDs_PF_15_M1_S3_all_contigs_in_LG2.txt", header = FALSE)
PF_LG2_markers <- read.table("IDs_PF_15_M1_S3_all_contigs_in_LG2_containing_social_markers.txt", header = FALSE)
PF_all_conts <- read.table("PF_conts.txt", stringsAsFactors = FALSE)

# 4. Split GFF3 into records with "gene" in "type" column
OT3b_GFF <- OT3b_GFF[OT3b_GFF$type == "gene", ]
PF_GFF <- PF_GFF[PF_GFF$type == "gene", ]

# 5. Split GFF3 into 3 subsets matching the 3 scaffold lists
OT3b_GFF_all <- OT3b_GFF[OT3b_GFF$seqid %in% OT3b_all$V1, ]
OT3b_GFF_LG2 <- OT3b_GFF[OT3b_GFF$seqid %in% OT3b_LG2$V1, ]
OT3b_GFF_LG2_markers <- OT3b_GFF[OT3b_GFF$seqid %in% OT3b_LG2_markers$V1, ]
OT3b_GFF_all_conts <- OT3b_GFF[OT3b_GFF$seqid %in% OT3b_all_conts$V1, ]

PF_GFF_all <- PF_GFF[PF_GFF$seqid %in% PF_all$V1, ]
PF_GFF_LG2 <- PF_GFF[PF_GFF$seqid %in% PF_LG2$V1, ]
PF_GFF_LG2_markers <- PF_GFF[PF_GFF$seqid %in% PF_LG2_markers$V1, ]
PF_GFF_all_conts <- PF_GFF[PF_GFF$seqid %in% PF_all_conts$V1, ]

# 6. Subset each split to records matching each ID in the corresponding list, i.e. list scaffolds and their genes
# Funtion to splitting and counting
split_and_count <- function(GFF, IDlist){
  # Define a count vector to count hits to IDs
  count_vector <- c()
  # Define range
  list_range <- 1:nrow(IDlist)
  for(i in list_range){
    # Get IDs
    ID <- IDlist$V1[i]
    # Subset GFF to records matching ID
    subGFF <- GFF[GFF$seqid == ID, ]
    # Add count to count_vector
    count_vector <- c(count_vector, nrow(subGFF))
  }
  # Add count_vector as new column to IDlist, save as new df
  result <- IDlist
  result$count <- count_vector
  # Give nice names
  nice_names <- c("ID", "Gene_Count")
  colnames(result) <- nice_names
  #Return the resut
  return(result)
}

# Call split_and_count
count_OT3b_all <- split_and_count(OT3b_GFF_all, OT3b_all)
count_OT3b_LG2 <- split_and_count(OT3b_GFF_LG2, OT3b_LG2)
count_OT3b_LG2_markers <- split_and_count(OT3b_GFF_LG2_markers, OT3b_LG2_markers)
count_OT3b_all_conts <- split_and_count(OT3b_GFF_all_conts, OT3b_all_conts)

count_PF_all <- split_and_count(PF_GFF_all, PF_all)
count_PF_LG2 <- split_and_count(PF_GFF_LG2, PF_LG2)
count_PF_LG2_markers <- split_and_count(PF_GFF_LG2_markers, PF_LG2_markers)
count_PF_all_conts <- split_and_count(PF_GFF_all_conts, PF_all_conts)

# 7. Read in contig lengths
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\10x Genomes")
len_OT_all <- read.table("len_OT_all_contigs_with_social_markers.fa")
len_OT_LG2 <- read.table("len_OT_LG2_contigs.fa")
len_OT_LG2_markers <- read.table("len_OT_LG2_contigs_with_social_markers.fa")
len_OT_all_conts <- read.table("")

len_PF_all <- read.table("len_PF_all_contigs_with_social_markers.fa")
len_PF_LG2 <- read.table("len_PF_LG2_contigs.fa")
len_PF_LG2_markers <- read.table("len_PF_LG2_contigs_with_social_markers.fa")
len_PF_all_conts <- read.table("")

len_all_contigs_OT <- read.table("len_OT_all_contigs_of_interest.fa")
len_all_contigs_PF <- read.table("len_PF_all_contigs_of_interest.fa")

# 8. Order by contig ID
len_OT_all <- len_OT_all[order(len_OT_all$V1),]
len_OT_LG2 <- len_OT_LG2[order(len_OT_LG2$V1),]
len_OT_LG2_markers <- len_OT_LG2_markers[order(len_OT_LG2_markers$V1),]
len_PF_all <- len_PF_all[order(len_PF_all$V1),]
len_PF_LG2 <- len_PF_LG2[order(len_PF_LG2$V1),]
len_PF_LG2_markers <- len_PF_LG2_markers[order(len_PF_LG2_markers$V1),]

count_OT3b_all <- count_OT3b_all[order(count_OT3b_all$ID),]
count_OT3b_LG2 <- count_OT3b_LG2[order(count_OT3b_LG2$ID),]
count_OT3b_LG2_markers <- count_OT3b_LG2_markers[order(count_OT3b_LG2_markers$ID),]
count_PF_all <- count_PF_all[order(count_PF_all$ID),]
count_PF_LG2 <- count_PF_LG2[order(count_PF_LG2$ID),]
count_PF_LG2_markers <- count_PF_LG2_markers[order(count_PF_LG2_markers$ID),]

len_all_contigs_OT <- len_all_contigs_OT[order(len_all_contigs_OT$V1),]
len_all_contigs_PF <- len_all_contigs_PF[order(len_all_contigs_PF$V1),]


# 9. Scatter Plots

#Scatter Plots, green, red, blue. Three scatterplot per assembly, coloured by dataset.
par(mfrow=c(2,3))
plot(count_OT3b_all$Gene_Count, len_OT_all$V2)
plot(count_OT3b_LG2$Gene_Count, len_OT_LG2$V2)
plot(count_OT3b_LG2_markers$Gene_Count, len_OT_LG2_markers$V2)

plot(count_PF_all$Gene_Count, len_PF_all$V2)
plot(count_PF_LG2$Gene_Count, len_PF_LG2$V2)
plot(count_PF_LG2_markers$Gene_Count, len_PF_LG2_markers$V2)

# Three graphs, one per dataset, dataset for both assemblies plotted per graph.
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Images")
png(filename="OT3b_and_PF_contigs_of_interest5.png",
    width = 800,
    height = 800,
    pointsize = 20
    )
par(mfrow=c(2,3))
#par(mfrow=c(1,3))
#Dataset A
plot(count_PF_all$Gene_Count, 
     len_PF_all$V2, 
     #ylim = c(0, 11000000), 
     #xlim = c(0, 1000),
     main = "All_Social_Markers",
     xlab = "",
     ylab = "",
     las = 1,
     cex.axis = 1, cex.lab = 1.2, cex.main = 1.5, #cex = 1.2,
     col="red", pch=16)
abline(lm(len_PF_all$V2~count_PF_all$Gene_Count), col = "red")
points(count_OT3b_all$Gene_Count, 
       len_OT_all$V2, 
       ylim = c(0, 11000000), 
       xlim = c(0, 1000),
       col = "deepskyblue2", pch=16)
abline(lm(len_OT_all$V2~count_OT3b_all$Gene_Count), col = "deepskyblue2")

#Dataset B
plot(count_PF_LG2$Gene_Count,
     len_PF_LG2$V2, 
     ylim = c(0, 11000000), 
     xlim = c(0, 1000),
     main = "LG2",
     xlab = "Gene Count",
     ylab = "",
     las = 1,
     cex.axis = 1, cex.lab = 1.2, cex.main = 1.5, #cex = 1.2,
     col="red", pch=16)
abline(lm(len_PF_LG2$V2~count_PF_LG2$Gene_Count), col = "red")
points(count_OT3b_LG2$Gene_Count, 
     len_OT_LG2$V2,
     col="deepskyblue2", pch=16)
abline(lm(len_OT_LG2$V2~count_OT3b_LG2$Gene_Count), col = "deepskyblue2")

#Dataset C
plot(count_PF_LG2_markers$Gene_Count, 
     len_PF_LG2_markers$V2,
     ylim = c(0, 11000000), 
     xlim = c(0, 1000),
     main = "Social_Region",
     xlab = "",
     ylab = "",
     las = 1,
     cex.axis = 1, cex.lab = 1.2, cex.main = 1.5, #cex = 1.2,
     col="red", pch=16)
abline(lm(len_PF_LG2_markers$V2~count_PF_LG2_markers$Gene_Count), col = "red")
points(count_OT3b_LG2_markers$Gene_Count, 
       len_OT_LG2_markers$V2, col = "deepskyblue2", pch=16)
abline(lm(len_OT_LG2_markers$V2~count_OT3b_LG2_markers$Gene_Count), col = "deepskyblue2")
dev.off()