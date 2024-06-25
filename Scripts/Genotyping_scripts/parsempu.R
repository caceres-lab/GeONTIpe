
#repos <- "https://cran.rstudio.com"
#
#packages <- c("stringr", "plot.matrix", "tidyr")
#
#for (pkg in packages) {
#  if (!requireNamespace(pkg, quietly = TRUE)) {
#    install.packages(pkg, repos = repos, dependencies = TRUE)
#    library(pkg, character.only = TRUE)
#  } else {
#    library(pkg, character.only = TRUE)
#  }
#}

library(stringr)
library(dplyr)
library(plot.matrix)
library(tidyr)

input <- commandArgs(T)

## For testing
#input <- c("varianteHsInv0410",
#           "HsInv0410_Genotype.txt",
#           "snpsHsInv0410.txt",
#           "SDs",
#           "HsInv0410")

setwd(paste(input[9]))

mpu <- read.table(input[1], sep = "\t", header = F)

min_prop_snps <- as.numeric(input[6])
min_cov <- as.numeric(input[7])
tree_distance <- as.numeric(input[8])

infoReadsT <- read.table(input[2], sep = "\t", header = T)
rownames(infoReadsT) <- infoReadsT$Read
infoReads <- infoReadsT[infoReadsT$Read != "FinalGenotype", "Read"]

vcf <- read.table(input[3], sep = "\t", header = F)
'%!in%' <- Negate("%in%")
dups <- vcf$V2[duplicated(vcf$V2)]

vcfaux <- vcf[vcf$V2 %!in% dups,]

for(i in dups){
  ss <- vcf[vcf$V2 == i,]
  vcfaux <- rbind(vcfaux, ss[1,])
}

vcf <- vcfaux[order(vcfaux$V2),]
rownames(vcf) <- vcf$V2

mpu <- mpu %>% distinct(V1, .keep_all = TRUE)
rownames(mpu) <- mpu$V1

mpu$nreads <- sapply(mpu$V4, function(x){
  length(str_split(x, ",", simplify = T))
})

SDs <- read.table(input[4], sep = "\t", header = F)
inv <- input[5]

colnames(SDs) <- c("Inv", "Chr", "Init", "End")
rownames(SDs) <- SDs$Inv

sp <- mpu[,c(1,3)]
colnames(sp) <- c("Pos", "Or")
rownames(sp) <- sp$Pos

sp$Mod <- sapply(sp$Or, function(x){
  ## Make commas and dots uniform
  y <- str_replace_all(x, ",", ".")
  
  ## Lowercase nucleotides
  y <- str_replace_all(y, "A", "a")
  y <- str_replace_all(y, "C", "c")
  y <- str_replace_all(y, "G", "g")
  y <- str_replace_all(y, "T", "t")
  
  ## Starts and ends
  y <- str_replace_all(y, "\\.\\$", ".")
  y <- str_replace_all(y, "\\.\\^.", ".")
  
  ## Deletions
  y <- str_replace_all(y, "\\*", "N")
  
  ## Next base deletions or insertions
  y <- str_replace_all(y, "\\.\\+", "I")
  y <- str_replace_all(y, "\\.-", "D")
  
  ## End substituting the matches
  y <- str_replace_all(y, "\\.", "M")
  
  ## Add and remove commas to separate the variants
  y <- str_replace_all(y, "M", ",M")
  y <- str_replace_all(y, "N", ",N")
  y <- str_replace_all(y, "I", ",I")
  y <- str_replace_all(y, "D", ",D")
  
  y <- str_replace_all(y, "^,", "")
  
  ## Add commas between alternative nucleotides
  ## Also checking if the deletions/insertions have the exact length
  
  if(grepl("a|c|g|t", y)){
    string <- y
    spl_st <- as.character(str_split(string, ",", simplify = T))
    
    nst <- c()
    
    for(row in spl_st){
      if(grepl("D|I", row)){
        num <- as.numeric(str_match(row, "\\d+"))
        if(num > 1){
          nrow <- gsub(pattern = paste0("(",as.character(num),"[a-z]{",as.character(num),"})"),
                       x = row,
                       replacement = "\\1,", perl = T)
          
        } else {
          nrow <- gsub(pattern = paste0("([1-9]+[a-z])"),
                       x = row,
                       replacement = "\\1,")
        }
        
      } else {
        if(nchar(row > 1)){
          aux <- as.character(str_split(row, "", simplify = T))
          nrow <- paste0(paste0(aux, collapse = ","),",")
        } else {
          nrow <- row
        }
      }
      
      nst <- c(nst, nrow)
    }
    
    y <- paste0(nst, collapse = "")
    y <- str_replace(y, ",$", "")
  }
  
  ## Check again the not introduced commas after filtering the I/D rearrangements
  patt2det <- grepl(",[a-z]{2}", y)
  
  while(patt2det){
    y <- gsub(pattern = "(,[a-z]{1})([a-z]{1})",
              x = y,
              replacement = "\\1,\\2")
    patt2det <- grepl(",[a-z]{2}", y)
  }
  
  y <- gsub(pattern = "([a-z]{1})([A-Z]{1})",
            x = y,
            replacement = "\\1,\\2")
  y <- str_replace(y, ",$", "")
  
  y
})

multi <- sp[sapply(sp$Mod, function(x){
  if(length(unique(as.character(str_split(x, ",", simplify = T)))) > 2 || length(unique(as.character(str_split(x, ",", simplify = T)))) == 1){
    T
  } else {F}
  
}),]

## Filter multiallelics
sp <- sp[!rownames(sp) %in% rownames(multi), ]

sp$type <- sapply(sp$Mod, function(x){
  spx <- str_split(x, ",", simplify = T)
  if(!any(spx %in% c("a", "c", "t", "g"))){
    "INDEL"
  } else {
    "SNP"
  }
  
})

sp <- sp[sp$type == "SNP",]

sp$maf <- sapply(sp$Mod, function(x){
  string <- str_split(x, "," , simplify = T)
  tot <- length(string)
  numOcc <- c()
  for (st in unique(as.character(string))){
    numOcc <- c(numOcc, length(string[string == st]))
  }
  nn <- sort(numOcc, TRUE)[2]
  freq <- nn/tot
})

sp$cov <- sapply(sp$Mod, function(x){
  length(str_split(x, ",", simplify = T))
})

sp <- sp[sp$maf > min_prop_snps, ]
sp <- sp[sp$cov > min_cov,]

if (dim(sp)[1] == 0) {
  message_error <- "Not biallelic variants"
  write.table(message_error, paste0(inv,"/SNP/",inv, "_error1.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  stop(message_error)
}

sp$reads <- mpu[as.character(sp$Pos),"V4"]

sp$First <- sapply(sp$Mod, function(x){
  string <- str_split(x, "," , simplify = T)
  cuenta <- as.data.frame(matrix(nrow = 0, ncol = 2))
  for (a in unique(as.character(string))){
    cuenta <- rbind(cuenta, c(a, length(string[string == a])))
  }
  colnames(cuenta) <- c("Nuc", "Freq")
  cuenta <- cuenta[order(cuenta$Freq, decreasing = T),]
  
  cuenta[1,"Nuc"]
  
})

sp$Second <- sapply(sp$Mod, function(x){
  string <- str_split(x, "," , simplify = T)
  cuenta <- as.data.frame(matrix(nrow = 0, ncol = 2))
  for (a in unique(as.character(string))){
    cuenta <- rbind(cuenta, c(a, length(string[string == a])))
  }
  colnames(cuenta) <- c("Nuc", "Freq")
  cuenta <- cuenta[order(cuenta$Freq, decreasing = T),]
  
  cuenta[2,"Nuc"]
  
})

VarTable <- as.data.frame(sp$Pos)
colnames(VarTable) <- "Var"
rownames(VarTable) <- VarTable$Var
VarTable$First <- sp$First
VarTable$Second <- sp$Second

allReads <- unique(as.character(str_split(paste0(sp$reads, collapse = ","), ",", simplify = T)))
allReads <- allReads[allReads != "*"]

TagMtx <- as.data.frame(matrix(nrow = length(allReads), ncol = dim(sp)[1]))
rownames(TagMtx) <- allReads
colnames(TagMtx) <- rownames(sp)

for(pos in sp$Pos){
  vars <- as.character(str_split(sp[as.character(pos),"Mod"], ",", simplify = T))
  reads <- as.character(str_split(sp[as.character(pos),"reads"], ",", simplify = T))
  
  for(i in seq(length(reads))){
    if(reads[i] != "*"){
      read <- reads[i]
      var <- vars[i]
      
      if(var == VarTable[as.character(pos),"First"]){
        TagMtx[read, as.character(pos)] <- 1
        
      } else {
        TagMtx[read, as.character(pos)] <- -1
      }
      
    }
    
  }
  
}

cutTag <- TagMtx
cutTag <- cutTag[infoReads,colnames(cutTag), drop = F]
cutTag <- cutTag[,colSums(is.na(cutTag))<nrow(cutTag), drop =F]

### Filter out by real SNPs checking

SNPChecking <- as.data.frame(VarTable$Var)
colnames(SNPChecking) <- "Position"
rownames(SNPChecking) <- SNPChecking$Position

rownames(vcf) <- vcf$V2
SNPChecking$REF <- vcf[rownames(SNPChecking),"V4"]
SNPChecking$ALT <- vcf[rownames(SNPChecking),"V5"]
SNPChecking$ReadALT <- sapply(rownames(SNPChecking), function(x){
  vars <- toupper(VarTable[x, c(2,3)])
  vars[vars != "M"]
})

SNPChecking$RightAlt <- (SNPChecking$ALT == SNPChecking$ReadALT)

InfoPos <- SNPChecking[SNPChecking$RightAlt, "Position"]
cutTag <- cutTag[,colnames(cutTag) %in% InfoPos, drop = F]

Side <- as.data.frame(matrix(ncol = 0, nrow = 1))
for(c in colnames(cutTag)){
  num <- as.numeric(c)
  if(num > SDs[SDs$Inv == paste0(inv,"_A"),"Init"] && num < SDs[SDs$Inv == paste0(inv,"_A"),"End"]){
    Side[[c]] <- "1"
  } else if(num > SDs[SDs$Inv == paste0(inv,"_BC"),"Init"] && num < SDs[SDs$Inv == paste0(inv,"_BC"),"End"]){
    Side[[c]] <- "2"
  } else if(num > SDs[SDs$Inv == paste0(inv,"_D"),"Init"] && num < SDs[SDs$Inv == paste0(inv,"_D"),"End"]){
    Side[[c]] <- "3"
  } else {
    Side[[c]] <- NA
  }
}

sideA <- subset(Side, select = which(Side %in% c(1, 2)))
sideD <- subset(Side, select = which(Side %in% c(2, 3)))

# Eliminar las columnas con solo un valor
one_value_cols <- !apply(cutTag, 2, function(x) length(x[!is.na(x)])) == 1

Decision <- cutTag[, one_value_cols, drop = FALSE]
colnames(Decision) <- colnames(cutTag)[one_value_cols]
rownames(Decision) <- rownames(cutTag)

if (dim(Decision)[2] == 0) {
  message_error <- "Not variants expected"
  write.table(message_error, paste0(inv, "/SNP/",inv, "_error2.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  stop(message_error)
}


SideDecision <- subset(Side, select = which(colnames(Side) %in% colnames(Decision)))

# Graficamos la matriz

png(filename = paste0(inv, "/SNP/", "SNPs_",inv, ".png"), width = 10000, height = 5000)
cols <- c(colorRampPalette(c("red", "white"))(100), colorRampPalette(c("white", "blue"))(100))
colsABCD <- c("lightblue", "lightyellow", "lightgreen")
layout(matrix(c(1,2), nrow = 2),
       heights = c(1,25))
plot(as.matrix(SideDecision), na.col = "white", key = NULL, col = colsABCD, breaks = c("1","2","3"))
plot(as.matrix(Decision), na.col = "white", col = cols, key = NULL)
dev.off()

DecisionA <- as.data.frame(Decision[, colnames(Decision) %in% colnames(sideA)])
oneA <- !apply(DecisionA, 2, function(x) length(x[!is.na(x)])) == 1
DecisionA <- DecisionA[, oneA, drop = FALSE]
DecisionD <- as.data.frame(Decision[, colnames(Decision) %in% colnames(sideD)])
oneD <- !apply(DecisionD, 2, function(x) length(x[!is.na(x)])) == 1
DecisionD <- DecisionD[, oneD, drop = FALSE]

if (dim(Decision)[1] == 0) {
  message_error <- "Only 1 read on the SNPs"
  write.table(message_error, paste0(inv, "/SNP/",inv, "_error3.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  stop(message_error)
} else if (dim(DecisionA)[1] == 1){
  message_error <- "Only 1 read in A"
  write.table(message_error, paste0(inv, "/SNP/",inv, "_error5.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  stop(message_error)
} else if (dim(DecisionD)[1] == 1){
  message_error <- "Only 1 read in D"
  write.table(message_error, paste0(inv, "/SNP/",inv, "_error4.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  stop(message_error)
}

## cluster para zona A,B,C

if (dim(DecisionA)[2] > 1){
  matrizDistA <- matrix(nrow = dim(DecisionA)[1], ncol = dim(DecisionA)[1])
  rownames(matrizDistA) <- rownames(DecisionA)
  colnames(matrizDistA) <- rownames(DecisionA)
  
  for(i in seq(dim(DecisionA)[1])){
    for(j in seq(dim(DecisionA)[1])){
      rows <- DecisionA[i,,drop = F]
      rows <- rbind(rows, DecisionA[j,,drop = F])
      
      rows <- rows[, colSums(is.na(rows)) == 0, drop = F]
      if(dim(rows)[2] == 0){
        matrizDistA[i,j] <- NA
        next
      }
      
      rows[rows == -1] <- 0
      
      distVect <- sapply(colnames(rows), function(x){
        abs(rows[1,x] - rows[2,x])
      })
      
      matrizDistA[i,j] <- sum(distVect)/length(colnames(rows))
      
    }
  }
  
  rows_with_na <- apply(matrizDistA, 1, function(x) all(is.na(x)))
  matrizDistA <- matrizDistA[!rows_with_na, ]
  cols_with_na <- apply(matrizDistA, 2, function(x) all(is.na(x)))
  matrizDistA <- matrizDistA[,!cols_with_na ]
  
  means <- colMeans(matrizDistA, na.rm = TRUE)
  matrizDistA[is.na(matrizDistA)] <- means[col(matrizDistA)[is.na(matrizDistA)]]
  
  clustA <- hclust(as.dist(matrizDistA))
  
  cutTreeA <- cutree(clustA, h = tree_distance)
  
  if (dim(matrizDistA)[1] > 2){ 
    png(paste0(inv, "/SNP/",inv,"_clusterA.png"), width = 1000, height = 1000)
    plot(clustA)
    dev.off()
  } else {
    message_error <- "Matrix of 2X2 for A"
    write.table(message_error, paste0(inv, "/SNP/",inv, "_error6.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
} else {
  message_error <- "No possible matrix, only one SNP position on A"
  write.table(message_error, paste0(inv, "/SNP/",inv, "_error8.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}
## cluster para zona B,C,D

if (dim(DecisionD)[2] > 1){
  matrizDistD <- matrix(nrow = dim(DecisionD)[1], ncol = dim(DecisionD)[1])
  rownames(matrizDistD) <- rownames(DecisionD)
  colnames(matrizDistD) <- rownames(DecisionD)
  
  for(i in seq(dim(DecisionD)[1])){
    for(j in seq(dim(DecisionD)[1])){
      rows <- DecisionD[i,,drop = F]
      rows <- rbind(rows, DecisionD[j,,drop = F])
      
      rows <- rows[, colSums(is.na(rows)) == 0, drop = F]
      if(dim(rows)[2] == 0){
        matrizDistD[i,j] <- NA
        next
      }
      
      rows[rows == -1] <- 0
      
      distVect <- sapply(colnames(rows), function(x){
        abs(rows[1,x] - rows[2,x])
      })
      
      matrizDistD[i,j] <- sum(distVect)/length(colnames(rows))
      
    }
  }
  
  rows_with_na <- apply(matrizDistD, 1, function(x) all(is.na(x)))
  matrizDistD <- matrizDistD[!rows_with_na, ]
  cols_with_na <- apply(matrizDistD, 2, function(x) all(is.na(x)))
  matrizDistD <- matrizDistD[,!cols_with_na ]
  
  means <- colMeans(matrizDistD, na.rm = TRUE)
  matrizDistD[is.na(matrizDistD)] <- means[col(matrizDistD)[is.na(matrizDistD)]]
  
  clustD <- hclust(as.dist(matrizDistD))
  
  cutTreeD <- cutree(clustD, h = tree_distance)
  
  if (dim(matrizDistD)[1] > 2){ 
    png(paste0(inv, "/SNP/",inv,"_clusterD.png"), width = 1000, height = 1000)
    plot(clustD)
    dev.off()
  } else {
    message_error <- "Matrix of 2X2 for D"
    write.table(message_error, paste0(inv, "/SNP/",inv, "_error7.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
} else {
  message_error <- "No possible matrix, only one SNP position on D"
  write.table(message_error, paste0(inv, "/SNP/",inv, "_error9.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

if (length(unique(cutTreeA)) == 2){
  write.table(cutTreeA, paste0(inv, "/SNP/",inv,"_resolvedA.txt"), sep = "\t", col.names = FALSE, row.names = TRUE, quote = FALSE)
} else {
  write.table(cutTreeA, paste0(inv, "/SNP/",inv,"_not_resolvedA.txt"), sep = "\t", col.names = FALSE, row.names = TRUE, quote = FALSE)
}

if (length(unique(cutTreeD)) == 2){
  write.table(cutTreeD, paste0(inv, "/SNP/",inv,"_resolvedD.txt"), sep = "\t", col.names = FALSE, row.names = TRUE, quote = FALSE)
} else {
  write.table(cutTreeD, paste0(inv, "/SNP/",inv,"_not_resolvedD.txt"), sep = "\t", col.names = FALSE, row.names = TRUE, quote = FALSE)
}