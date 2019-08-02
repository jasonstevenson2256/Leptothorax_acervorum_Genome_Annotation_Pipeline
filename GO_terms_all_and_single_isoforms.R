#### Find number of GO terms for each gene and each isoform ####

library(data.table)

#1. Setwd and read in all isoform data
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\GO FEAT")
GOfeat_OT <- fread(file="GO_FEAT_OT3b_all.csv")
GOfeat_PF <- fread(file="GO_FEAT_PF_all.csv")

# Read in single isoform data
GOfeat_OT_t1 <- fread(file="GO_FEAT_OT_single_isoforms.txt")
GOfeat_PF_t1 <- fread(file="GO_FEAT_PF_single_isoforms.txt")

# Give column names
column_names <- c("#", "Locus tag", "Length", "Product", "Completeness", "Gene onthology", 
                  "SEED", "Protein databases", "Genome annotation databases", "Family and domain databases", 
                  "Crossreferences", "V12")
colnames(GOfeat_OT) <- column_names
colnames(GOfeat_PF) <- column_names
colnames(GOfeat_OT_t1) <- column_names
colnames(GOfeat_PF_t1) <- column_names

#2. Create paired vectors matching individual gene IDs with individual GO terms
#   - equal lengths, gene IDs are repeated once per associated GO term.  
#   Then merge together into a to dataframe.
#
#First a function
make.result.df <- function(df){
  #Create empty output vectors
  out.ID <- c()
  out.GO <- c()
  #Iterate over gene IDs
  for(i in 1:nrow(df)){
    #For each ID, get the associated GO terms column
    tosplit <- df$`Gene onthology`[i]
    #Split on pipe ("|") characters
    split1 <- strsplit(tosplit, "\\|")
    #Coerce to vector
    split1 <- unlist(split1)
    #Create empty vector to hold split elements of split1
    vec <- c()
    #Iterate over split1 vector
    for(j in 1:length(split1)){
      #Populate vec with a second split of each element in split1
      #to recover just the GO terms
      term <- substr(split1[j], start = 1, stop = 10)
      vec <- c(vec, term)
    }
    #Populate out.ID
    out.ID <- c(out.ID, rep(df$`Locus tag`[i], each=length(vec)))
    #Populate out.GO
    out.GO <- c(out.GO, vec)
  }
  #Stick out.ID and out.GO together in one df
  result <- data.frame(out.ID, out.GO)
  #Provide appropriate column names
  names(result) <- c("ID", "GO_term")
  #Remove rows without GO terms in them
  #result <- result[ grep("GO", result$GO_term) , ]
  #Return result
  return(result)
}
#
#Then call the function
result_OT <- make.result.df(GOfeat_OT)
result_PF <- make.result.df(GOfeat_PF)
result_OT_t1 <- make.result.df(GOfeat_OT_t1)
result_PF_t1 <- make.result.df(GOfeat_PF_t1)

#3. find unique GO_term values for each gene
# For OT
result_OT_stripped <- result_OT
# subset to remive .ts (isoform label)
result_OT_stripped$ID <- sub("\\..*", "", result_OT$ID)
#uniq_OT_stripped <- data.frame(c(table(result_OT_stripped$GO_term)))
# extract unique possible ID tags
OT_IDs <- unique(result_OT_stripped$ID)
# store total go terms per gene
OT_totals<-c()

for(ID in OT_IDs){
  OT_section<-result_OT_stripped[result_OT_stripped$ID==ID,]
  OT_section<-na.omit(OT_section)
  count<-length(unique(OT_section$GO_term))
  OT_totals<-c(OT_totals,count)
}

# For PF
result_PF_stripped <- result_PF
result_PF_stripped$ID <- sub("\\..*", "", result_PF$ID)
#uniq_PF_stripped <- data.frame(c(table(result_PF_stripped$GO_term)))
PF_IDs <- unique(result_PF_stripped$ID)
PF_totals<-c()

for(ID in PF_IDs){
  PF_section<-result_PF_stripped[result_PF_stripped$ID==ID,]
  PF_section<-na.omit(PF_section)
  count<-length(unique(PF_section$GO_term))
  PF_totals<-c(PF_totals,count)
}

# For OT_t1
result_OT_t1_stripped <- result_OT_t1
result_OT_t1_stripped$ID <- sub("\\..*", "", result_OT_t1$ID)
#uniq_OT_t1_stripped <- data.frame(c(table(result_OT_t1_stripped$GO_term)))
OT_t1_IDs <- unique(result_OT_t1_stripped$ID)
OT_t1_totals<-c()

for(ID in OT_t1_IDs){
  OT_t1_section<-result_OT_t1_stripped[result_OT_t1_stripped$ID==ID,]
  OT_t1_section<-na.omit(OT_t1_section)
  count<-length(unique(OT_t1_section$GO_term))
  OT_t1_totals<-c(OT_t1_totals,count)
}

# For PF_t1
result_PF_t1_stripped <- result_PF_t1
result_PF_t1_stripped$ID <- sub("\\..*", "", result_PF_t1$ID)
#uniq_PF_t1_stripped <- data.frame(c(table(result_PF_t1_stripped$GO_term)))
PF_t1_IDs <- unique(result_PF_t1_stripped$ID)
PF_t1_totals<-c()

for(ID in PF_t1_IDs){
  PF_t1_section<-result_PF_t1_stripped[result_PF_t1_stripped$ID==ID,]
  PF_t1_section<-na.omit(PF_t1_section)
  count<-length(unique(PF_t1_section$GO_term))
  PF_t1_totals<-c(PF_t1_totals,count)
}

# concatenate into data frame

stripped_ID_totals_OT<-data.frame(OT_IDs,OT_totals)
names(stripped_ID_totals_OT) <- c('ID','GO_count')
stripped_ID_totals_PF<-data.frame(PF_IDs,PF_totals)
names(stripped_ID_totals_PF) <- c('ID','GO_count')
stripped_ID_totals_OT_t1<-data.frame(OT_t1_IDs,OT_t1_totals)
names(stripped_ID_totals_OT_t1) <- c('ID','GO_count')
stripped_ID_totals_PF_t1<-data.frame(PF_t1_IDs,PF_t1_totals)
names(stripped_ID_totals_PF_t1) <- c('ID','GO_count')

#4. Make graphs
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Images")

# Scatter plots of numer of terms in all vs single isoforms
# Make new dataframes of results ordered by GO count
OT_results_ordered <- stripped_ID_totals_OT[order(as.numeric(as.character(stripped_ID_totals_OT$GO_count))), ]
PF_results_ordered <- stripped_ID_totals_PF[order(as.numeric(as.character(stripped_ID_totals_PF$GO_count))), ]
# Order t1s in same order of IDs from above data frames
OT_t1_results_ordered <- stripped_ID_totals_OT_t1[order(match(stripped_ID_totals_OT_t1$ID, OT_results_ordered$ID)),]
PF_t1_results_ordered <- stripped_ID_totals_PF_t1[order(match(stripped_ID_totals_PF_t1$ID, PF_results_ordered$ID)),]

# Plot OT all vs single isoforms
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Images")

png(filename="OT3b_all_vs_single_isoforms_scatter.png")
par(mfrow=c(1,1))
plot(OT_t1_results_ordered$GO_count,
     main = "OT3b",
     xlab = "Genes (ranked by number of GO terms)",
     ylab = "GO Count",
     ylim = c(0, 40),
     las = 1,
     cex.axis = 0.8, cex.lab = 1, cex.main = 1.2, cex = 0.8,
     col="blue", pch=16)
points(OT_results_ordered$GO_count,
       col="deepskyblue2", pch=16, cex = 0.8)
dev.off()

# Plot PF all vs single isoforms
png(filename="PF_all_vs_single_isoforms_scatter.png")
par(mfrow=c(1,1))
plot(PF_t1_results_ordered$GO_count,
     main = "PF",
     xlab = "Genes (ranked by number of GO terms)",
     ylab = "GO Count",
     las = 1,
     cex.axis = 0.8, cex.lab = 1, cex.main = 1.2, cex = 0.8,
     col="orange", pch=16)
points(PF_results_ordered$GO_count,
       col="red", pch=16, cex = 0.8)
dev.off()

# OT vs PF all isoforms
par(mfrow=c(1,1))
plot(OT_results_ordered$GO_count,
     main = "OT vs PF All Isoforms",
     xlab = "Isoforms",
     ylab = "GO Count",
     las = 1,
     cex.axis = 1, cex.lab = 1.2, cex.main = 1.5, #cex = 1.2,
     col="deepskyblue2", pch=16)
points(PF_results_ordered$GO_count,
       col="red", pch=16)

