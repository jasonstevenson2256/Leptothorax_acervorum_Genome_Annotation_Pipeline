#### Wrangling GOfeat results files into correct format for analysis ####

library(data.table)

#1. Setwd and read in data
setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\GO FEAT")
GOfeat_OT <- fread(file="GO_FEAT_OT3b_all.csv")
GOfeat_PF <- fread(file="GO_FEAT_PF_all.csv")

##1.1 Optional - read in single isoform data instead
GOfeat_OT <- fread(file="GO_FEAT_OT_single_isoforms.txt")
GOfeat_PF <- fread(file="GO_FEAT_PF_single_isoforms.txt")

column_names <- c("#", "Locus tag", "Length", "Product", "Completeness", "Gene onthology", 
                  "SEED", "Protein databases", "Genome annotation databases", "Family and domain databases", 
                  "Crossreferences", "V12")
colnames(GOfeat_OT) <- column_names
colnames(GOfeat_PF) <- column_names

#2. Subset data down to just records with GO terms
GOfeat_OT <- GOfeat_OT[!GOfeat_OT$`Gene onthology`=="", ]
GOfeat_PF <- GOfeat_PF[!GOfeat_PF$`Gene onthology`=="", ]

#3. Create paired vectors matching individual gene IDs with individual GO terms
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
  result <- result[ grep("GO", result$GO_term) , ]
  #Return result
  return(result)
}
#
#Then call the function
result_OT <- make.result.df(GOfeat_OT)
result_PF <- make.result.df(GOfeat_PF)

#4. Write result dfs to file
write.table(result_OT, file="OT3b_IDs_and_terms.txt", sep="\t", 
            row.names=FALSE, col.names=FALSE, quote=TRUE)
write.table(result_PF, file="PF_IDs_and_terms.txt", sep="\t", 
            row.names=FALSE, col.names=FALSE, quote=TRUE)


