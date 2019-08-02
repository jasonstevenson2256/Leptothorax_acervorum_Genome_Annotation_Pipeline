setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Data\\10x Assembly Results\\BLAST")

OvP <- read.table("OTvPF_alignment_distribution.txt")
PvO <- read.table("PFvOT_alignment_distribution.txt")

setwd("\\Users\\Jason\\Documents\\University work (updated)\\MSc\\Project\\Images")

png(filename="OTvPF_align_hist.png")
par(mfrow=c(1,1))
hist(OvP$V1,
     breaks = c(20),
     xlim = c(-500, 50000),
     main = "Q=OT, S=PF",
     xlab="Missing Alignment (bp)",
     border = "black",
     col = "deepskyblue2")
dev.off()

png(filename="PFvOT_align_hist.png")
par(mfrow=c(1,1))
hist(PvO$V1, 
     breaks = c(50),
     xlim = c(-500, 50000),
     main = "Q=PF, S=OT",
     xlab="Missing Alignment (bp)",
     border = "black",
     col = "red")
dev.off()