#!/usr/bin/Rscript

#repos <- "https://cran.rstudio.com"
#packages <- c("ggplot2")
#
#for (pkg in packages) {
#  if (!requireNamespace(pkg, quietly = TRUE)) {
#    install.packages(pkg, repos = repos, dependencies = TRUE)
#    library(pkg, character.only = TRUE)
#  } else {
#    library(pkg, character.only = TRUE)
#  }
#}

library(ggplot2)

inputs <- commandArgs(T)

#inputs <- c("C:/Users/Rikip/OneDrive/Escritorio", "HsInv0020", "M")

inv <- inputs[2]
gender <- inputs[3]
minor_al <- as.numeric(inputs[4])
low_con <- as.numeric(inputs[5])
diff_exp <- as.numeric(inputs[6])
min_sv <- as.numeric(inputs[7])
splitting_svs <- as.numeric(inputs[8])

setwd(paste(inputs[1], inv, sep = "/"))

## Read all files

plotable <- read.table("seqInfo", sep = "\t", header = F)
colnames(plotable) <- c("chr", "p1", "p2")
rownames(plotable) <- c(inv)

mapinfo <- read.table(file = "Resultados/MapInfo", sep = "\t", header = F)
colnames(mapinfo) <- c("Read", "Sonda", "Signo")
mapinfo$Signo <- as.character(mapinfo$Signo)

if(file.size("Resultados/sondaErronea") > 0){
  sonErr <- unique(read.table("Resultados/sondaErronea", sep = "\t", header = F)$V1)
} else {
  sonErr <- c()
}

dupReads <- unique(mapinfo$Read)
'%!in%' <- Negate('%in%')
dupReads <- dupReads[dupReads %!in% sonErr]

### Redo tc

tc <- as.data.frame(dupReads)
colnames(tc) <- "Read"
rownames(tc) <- tc$Read

tc$A <- sapply(tc$Read, function(x){
  if(file.exists(paste0("Resultados/Distancias/Reads/", x, "/BLASTCoord/A"))){
    T
  } else {F}
})

tc$B <- sapply(tc$Read, function(x){
  if(file.exists(paste0("Resultados/Distancias/Reads/", x, "/BLASTCoord/B"))){
    T
  } else {F}
})

tc$C <- sapply(tc$Read, function(x){
  if(file.exists(paste0("Resultados/Distancias/Reads/", x, "/BLASTCoord/C"))){
    T
  } else {F}
})

tc$D <- sapply(tc$Read, function(x){
  if(file.exists(paste0("Resultados/Distancias/Reads/", x, "/BLASTCoord/D"))){
    T
  } else {F}
})

dupReads <- dupReads[sapply(dupReads, function(x){
  if(sum(tc[tc$Read == x, c("A", "B", "C", "D")]) > 1){
    T
  } else {F}
})]

tc <- tc[dupReads,]

tablaSignos <- tc

InvDist <- read.table("../../../Infor/ExpectedDistProvesInverted.txt",
                      sep = "\t", header = T)

rownames(InvDist) <- InvDist$Inv

## Read information about the probes

sondas <- new.env()
for(sonda in c("A", "B", "C", "D")){
  sondas[[sonda]] <- read.table(file = paste0("CoordenadasRelativas/",sonda), header = F, sep = "\t")
}

DistSondas <- new.env()
relDir <- "Resultados/Distancias/Reads"

for (read in dupReads){
  sonPres <- dir(paste(relDir, read, "BLASTCoord" ,sep = "/"))
  if(length(sonPres) > 1){
    for(sp in sonPres){
      DistSondas[[read]][[sp]] <- read.table(file = paste(relDir, read, "BLASTCoord", sp, sep = "/"), 
                                             sep = "\t", header = F, stringsAsFactors = F)
    }
  }
}

seqInfo <- read.table(file = "seqInfo", sep = "\t", header = F)

dupReads <- ls(DistSondas)

tablaSignos <- tablaSignos[dupReads,]
tc <- tc[dupReads,]

mapinfo <- mapinfo[mapinfo$Read %in% dupReads,]

## Filtro los BC

nd <- c()

for(read in dupReads){
  sMap <- ls(DistSondas[[read]])
  
  chA <- "A" %in% sMap
  chD <- "D" %in% sMap
  
  if(sum(c(chA, chD)) > 0){
    nd <- c(nd, read)
  }
}

DistSondas <- mget(nd, envir = DistSondas)

dupReads <- ls(DistSondas)

tc <- tc[dupReads,]
mapinfo <- mapinfo[mapinfo$Read %in% dupReads,]

## Loop for each probe

for(sonda in c("A", "B", "C", "D")){
  tablaSignos[,paste0("signo", sonda)] <- NA
  
  for(read in dupReads){
    if(tablaSignos[read,sonda]){
      ori <- DistSondas[[read]][[sonda]]$V4
      if(ori == "plus"){
        tablaSignos[read,paste0("signo", sonda)] <- "+"
      }
      if(ori == "minus"){
        tablaSignos[read,paste0("signo", sonda)] <- "-"
      }
    }
    
  }
}


matrixPosibSignos <- as.data.frame(
  matrix(c("Std", "+", "+", "+", "+",
           "Inv", "+", "-", "-", "+",
           "Inv", "-", "+", "+", "-",
           "Std", "-", "-", "-", "-"),
         nrow = 4, byrow = T))

colnames(matrixPosibSignos) <- c("Orient", "A", "B", "C", "D")

tablaSignos$Profile <- sapply(tablaSignos$Read, function(x){
  cualSondas <- colnames(tablaSignos)[2:5][as.logical(tablaSignos[x,c("A", "B", "C", "D")])]
  combSig <- c()
  for(sig in cualSondas){
    combSig <- c(combSig, tablaSignos[x, paste0("signo", sig)])
  }
  
  subPosib <- matrixPosibSignos[,cualSondas]
  
  whichRow <- sapply(seq(4), function(i){
    identical(as.character(subPosib[i,]), combSig)
  })
  
  if(sum(whichRow) == 1){
    matrixPosibSignos$Orient[whichRow]
  } else if(sum(whichRow) == 0){
    "Err"
  } else {
    "AD"
  }
  
})


## Representing the plot of each probe and read

refPlot <- as.data.frame(t(unlist(c(1, as.numeric(seqInfo[3]-seqInfo[2]), 
                                    sondas[["A"]][1], sondas[["A"]][2],
                                    sondas[["B"]][1], sondas[["B"]][2],
                                    sondas[["C"]][1], sondas[["C"]][2],
                                    sondas[["D"]][1], sondas[["D"]][2]))))

colnames(refPlot) <- c("Init", "End", "Ai", "Ae", "Bi", "Be", "Ci", "Ce", "Di", "De")

dataPlot <- as.data.frame(dupReads, stringsAsFactors = F)
rownames(dataPlot) <- dataPlot$dupReads
dataPlot$ReadLength <- NA

for(sonda in c("A", "B", "C", "D")){
  dataPlot[,paste0(sonda, "i")] <- NA
  dataPlot[,paste0(sonda, "e")] <- NA
}

dataPlot$AB <- NA
dataPlot$AC <- NA
dataPlot$AD <- NA
dataPlot$BD <- NA
dataPlot$CD <- NA
dataPlot$ABr <- NA
dataPlot$ACr <- NA
dataPlot$ADr <- NA
dataPlot$BDr <- NA
dataPlot$CDr <- NA

for(read in dupReads){
  sonPres <- ls(DistSondas[[read]])
  dataPlot[read, "ReadLength"] <- DistSondas[[read]][[sonPres[1]]]$V3
  
  for(sonda in c("A", "B", "C", "D")){
    if(sonda %in% sonPres){
      dataPlot[read, paste0(sonda, "i")] <- DistSondas[[read]][[sonda]]$V1
      dataPlot[read, paste0(sonda, "e")] <- DistSondas[[read]][[sonda]]$V2
    }
  }
  
  if(!is.na(dataPlot[read,"Ai"])){
    if(!is.na(dataPlot[read,"Bi"])){
      dataPlot[read, "AB"] <- ((dataPlot[read, "Be"] + dataPlot[read, "Bi"])/2) -
        ((dataPlot[read, "Ae"] + dataPlot[read, "Ai"])/2)
      dataPlot[read, "ABr"] <- abs(dataPlot[read, "AB"])/
        (((refPlot$Be + refPlot$Bi)/2)-((refPlot$Ae + refPlot$Ai)/2))
    }
    if(!is.na(dataPlot[read,"Ci"])){
      dataPlot[read, "AC"] <- ((dataPlot[read, "Ce"] + dataPlot[read, "Ci"])/2) -
        ((dataPlot[read, "Ae"] + dataPlot[read, "Ai"])/2)
      dataPlot[read, "ACr"] <- abs(dataPlot[read, "AC"])/
        (((refPlot$Ce + refPlot$Ci)/2)-((refPlot$Ae + refPlot$Ai)/2))
    }
    if(!is.na(dataPlot[read,"Di"])){
      dataPlot[read, "AD"] <- ((dataPlot[read, "De"] + dataPlot[read, "Di"])/2) -
        ((dataPlot[read, "Ae"] + dataPlot[read, "Ai"])/2)
      dataPlot[read, "ADr"] <- abs(dataPlot[read, "AD"])/
        (((refPlot$De + refPlot$Di)/2)-((refPlot$Ae + refPlot$Ai)/2))
    }
  }
  if(!is.na(dataPlot[read,"Di"])){
    if(!is.na(dataPlot[read,"Bi"])){
      dataPlot[read, "BD"] <- ((dataPlot[read, "De"] + dataPlot[read, "Di"])/2) -
        ((dataPlot[read, "Be"] + dataPlot[read, "Bi"])/2)
      dataPlot[read, "BDr"] <- abs(dataPlot[read, "BD"])/
        (((refPlot$De + refPlot$Di)/2)-((refPlot$Be + refPlot$Bi)/2))
    }
    if(!is.na(dataPlot[read,"Ci"])){
      dataPlot[read, "CD"] <- ((dataPlot[read, "De"] + dataPlot[read, "Di"])/2) -
        ((dataPlot[read, "Ce"] + dataPlot[read, "Ci"])/2)
      dataPlot[read, "CDr"] <- abs(dataPlot[read, "CD"])/
        (((refPlot$De + refPlot$Di)/2)-((refPlot$Ce + refPlot$Ci)/2))
    }
  }
}

repGraf <- data.frame(matrix(nrow = 0, ncol = 13))

repP1 <- plotable$p1 - seqInfo$V2
repP2 <- (plotable$p2 - plotable$p1) + repP1

for(read in dupReads){
  forward <- sum(dataPlot[read, c("AB", "AC", "AD", "BD", "CD")], na.rm = T) > 0
  rRl <- dataPlot[read, "ReadLength"]
  
  if(!is.na(dataPlot[read, "Ai"])){
    rAi <- refPlot$Ai
    rAe <- refPlot$Ae
    
    if(!is.na(dataPlot[read, "Bi"])){
      rBi <- rAi + abs(dataPlot[read, "AB"])
      rBe <- rAe + abs(dataPlot[read, "AB"])
    } else {
      rBi <- NA
      rBe <- NA
    }
    if(!is.na(dataPlot[read, "Ci"])){
      rCi <- rAi + abs(dataPlot[read, "AC"])
      rCe <- rAe + abs(dataPlot[read, "AC"])
    } else {
      rCi <- NA
      rCe <- NA
    }
    if(!is.na(dataPlot[read, "Di"])){
      rDi <- rAi + abs(dataPlot[read, "AD"])
      rDe <- rAe + abs(dataPlot[read, "AD"])
    } else {
      rDi <- NA
      rDe <- NA
    }
    
    if(forward){
      rRi <- refPlot$Ai - dataPlot[read, "Ai"]
      if(rRi < 0){
        rRi <- 0
      }
      rRe <- refPlot$Ae + (rRl - dataPlot[read, "Ae"])
      if(rRe > refPlot$End){
        rRe <- refPlot$End
      }
    } else {
      rRi <- refPlot$Ai - (rRl - dataPlot[read, "Ai"])
      if(rRi < 0){
        rRi <- 0
      }
      rRe <- refPlot$Ae + dataPlot[read, "Ae"]
      if(rRe > refPlot$End){
        rRe <-  refPlot$End
      }
    }
    
  } else {
    rAi <- NA
    rAe <- NA
    
    rDi <- refPlot$Di
    rDe <- refPlot$De
    
    if(!is.na(dataPlot[read, "Bi"])){
      rBi <- rDi - abs(dataPlot[read, "BD"])
      rBe <- rDe - abs(dataPlot[read, "BD"])
    } else {
      rBi <- NA
      rBe <- NA
    }
    if(!is.na(dataPlot[read, "Ci"])){
      rCi <- rDi - abs(dataPlot[read, "CD"])
      rCe <- rDe - abs(dataPlot[read, "CD"])
    } else {
      rCi <- NA
      rCe <- NA
    }
    
    if(forward){
      rRi <- refPlot$Di - dataPlot[read, "Di"]
      if(rRi < 0){
        rRi <-  0
      }
      rRe <- refPlot$De + (rRl - dataPlot[read, "De"])
      if(rRe > refPlot$End){
        rRe <-  refPlot$End
      }
    } else {
      rRi <- refPlot$Di - (rRl - dataPlot[read, "Di"])
      if(rRi < 0){
        rRi <-  0
      }
      rRe <- refPlot$De + dataPlot[read, "De"]
      if(rRe > refPlot$End){
        rRe <-  refPlot$End
      }
    }
  }
  
  Or <- NA
  repGraf <- rbind(repGraf, c(read, rAi, rAe, rBi, rBe, rCi, rCe, rDi, rDe, rRi, rRe, rRl, Or))
  
}

repGraf <- rbind(repGraf, c("Reference", refPlot$Ai, refPlot$Ae, refPlot$Bi, refPlot$Be, refPlot$Ci, refPlot$Ce,
                            refPlot$Di, refPlot$De, NA, NA, NA, NA))
colnames(repGraf) <- c("Read", "Ai", "Ae", "Bi", "Be", "Ci", "Ce", "Di", "De", "Ri", "Re", "Rl", "Or")
rownames(repGraf) <- repGraf$Read

for(cn in colnames(repGraf)){
  repGraf[,cn] <- as.numeric(repGraf[,cn])
}

Recons <- as.data.frame(rownames(dataPlot))
colnames(Recons) <- "Read"
rownames(Recons) <- Recons$Read

Recons$Dist <- sapply(rownames(Recons), function(x){
  ss <- repGraf[x,]
  
  if(!is.na(ss$Ai)){
    if(!is.na(ss$Bi)){
      if(!is.na(ss$Ci)){
        if(ss$Ci > ss$Bi){
          "Std"
        } else {
          "Inv"
        }
      } else {
        "Std"
      }
    } else {
      "Inv"
    }
  } else {
    if(!is.na(ss$Ci)){
      if(!is.na(ss$Bi)){
        if(ss$Ci > ss$Bi){
          "Std"
        } else {
          "Inv"
        }
      } else {
        "Std"
      }
    } else {
      "Inv"
    }
  }
  
})

Recons$SignGT <- sapply(Recons$Read, function(x){
  tablaSignos[x, "Profile"]
})

Recons$ProbGeno <- NA

ELGENOTIPO <- NA
ProbGeno <- NA

gt_freq_geno <- table(Recons$SignGT)
gt_prop_geno <- prop.table(gt_freq_geno)

determine_genotype <- function(gt_prop_geno, Recons, plotable, gender, minor_al, min_low_con) {
  if (length(Recons$SignGT) == 0) {
    return(list(ELGENOTIPO = "ND", ProbGeno = NA))
  }
  
  if (length(Recons$SignGT) == 1) {
    if (plotable$chr == 24 && gender == "W") {
      return(list(ELGENOTIPO = "ND", ProbGeno = NA))
    } else {
      return(list(ELGENOTIPO = "NER", ProbGeno = NA))
    }
  }
  
  if (length(Recons$SignGT) <= 4) {
    if (length(names(gt_prop_geno[gt_prop_geno > minor_al])) > 2) {
      return(list(ELGENOTIPO = "Error", ProbGeno = NA))
    } else {
      return(handle_genotype_under_4(gt_prop_geno, Recons, plotable, gender, minor_al, min_low_con))
    }
  }
  
  if (length(Recons$SignGT) >= 5) {
    return(handle_genotype_over_5(gt_prop_geno, Recons, plotable, gender, minor_al))
  }
  
  return(list(ELGENOTIPO = "Failing_code", ProbGeno = NA))
}

handle_genotype_under_4 <- function(gt_prop_geno, Recons, plotable, gender, minor_al, min_low_con) {
  if (plotable$chr == 23) {
    if (gender == "W") {
      return(Main_condition(gt_prop_geno, minor_al, min_low_con))
    } else {
      return(Det_males(gt_prop_geno, Recons, minor_al))
    }
  }
  
  if (plotable$chr == 24) {
    if (gender == "W") {
      return(list(ELGENOTIPO = "ND", ProbGeno = NA))
    } else {
      return(Det_males(gt_prop_geno, Recons, minor_al))
    }
  }
  
  return(Main_condition(gt_prop_geno, minor_al, min_low_con))
}

handle_genotype_over_5 <- function(gt_prop_geno, Recons, plotable, gender, minor_al) {
  if (plotable$chr == 24) {
    if (gender == "W") {
      return(list(ELGENOTIPO = "ND", ProbGeno = NA))
    } else {
      return(Det_males(gt_prop_geno, Recons, minor_al))
    }
  }
  
  if (plotable$chr == 23 && gender == "M") {
    return(Det_males(gt_prop_geno, Recons, minor_al))
  }
  
  return(Main_condition(gt_prop_geno, minor_al, min_low_con))
}

Det_males <- function(gt_prop_geno, Recons, minor_al) {
  if (length(names(gt_prop_geno[gt_prop_geno >= (1 - minor_al)])) == 1) {
    geno <- names(gt_prop_geno[gt_prop_geno > (1 - minor_al)])
    prob <- round(1 - (0.5 ^ sum(Recons$SignGT == geno)), digits = 3)
    return(list(ELGENOTIPO = geno, ProbGeno = prob))
  } else {
    return(list(ELGENOTIPO = "Error", ProbGeno = NA))
  }
}

Main_condition <- function(gt_prop_geno, minor_al, min_low_con) {
  if ("Std" %in% names(gt_prop_geno) && gt_prop_geno["Std"] > minor_al) {
    if ("Inv" %in% names(gt_prop_geno) && gt_prop_geno["Inv"] > minor_al) {
      return(list(ELGENOTIPO = "Std/Inv", ProbGeno = 1))
    } else if ("Inv" %in% names(gt_prop_geno) && gt_prop_geno["Inv"] >= min_low_con && gt_prop_geno["Inv"] <= minor_al) {
      return(list(ELGENOTIPO = "Low_confi", ProbGeno = NA))
    } else if ("AD" %in% names(gt_prop_geno) && gt_prop_geno["AD"] > minor_al) {
      return(list(ELGENOTIPO = "Std/AD", ProbGeno = 1))
    } else if ("AD" %in% names(gt_prop_geno) && gt_prop_geno["AD"] >= min_low_con && gt_prop_geno["AD"] <= minor_al) {
      return(list(ELGENOTIPO = "Low_confi", ProbGeno = NA))
    } else if ("Error" %in% names(gt_prop_geno) && gt_prop_geno["Error"] > minor_al) {
      return(list(ELGENOTIPO = "Std/Err", ProbGeno = 1))
    } else {
      if (length(Recons$SignGT) <= 4) {
        return(list(ELGENOTIPO = "NER", ProbGeno = NA, SNP = "Y"))
      } else {
        if (length(names(gt_prop_geno[gt_prop_geno > (1 - minor_al)])) == 1) {
          geno <- names(gt_prop_geno[gt_prop_geno > (1 - minor_al)])
          prob <- round(1 - (0.5 ^ sum(Recons$SignGT == geno)), digits = 3)
          return(list(ELGENOTIPO = geno, ProbGeno = prob))
        } else {
          return(list(ELGENOTIPO = "Failing_code", ProbGeno = NA))
        }
      }
    }
  } else if ("Inv" %in% names(gt_prop_geno) && gt_prop_geno["Inv"] > minor_al) {
    if ("Std" %in% names(gt_prop_geno) && gt_prop_geno["Std"] >= min_low_con && gt_prop_geno["Std"] <= minor_al) {
        return(list(ELGENOTIPO = "Low_confi", ProbGeno = NA))
    } else if ("AD" %in% names(gt_prop_geno) && gt_prop_geno["AD"] > minor_al) {
        return(list(ELGENOTIPO = "Inv/AD", ProbGeno = 1))
    } else if ("AD" %in% names(gt_prop_geno) && gt_prop_geno["AD"] >= min_low_con && gt_prop_geno["AD"] <= minor_al) {
      return(list(ELGENOTIPO = "Low_confi", ProbGeno = NA))  
    } else if ("Error" %in% names(gt_prop_geno) && gt_prop_geno["Error"] > minor_al) {
      return(list(ELGENOTIPO = "Std/Err", ProbGeno = 1))
    } else {
      if (length(Recons$SignGT) <= 4) {
        return(list(ELGENOTIPO = "NER", ProbGeno = NA, SNP = "Y"))
      } else {
        if (length(names(gt_prop_geno[gt_prop_geno > (1 - minor_al)])) == 1) {
          geno <- names(gt_prop_geno[gt_prop_geno > (1 - minor_al)])
          prob <- round(1 - (0.5 ^ sum(Recons$SignGT == geno)), digits = 3)
          return(list(ELGENOTIPO = geno, ProbGeno = prob))
        } else {
          return(list(ELGENOTIPO = "Failing_code", ProbGeno = NA))
        }
      }
    }
  } else if ("AD" %in% names(gt_prop_geno) && gt_prop_geno["AD"] > minor_al){
      if ("Std" %in% names(gt_prop_geno) && gt_prop_geno["Std"] >= min_low_con && gt_prop_geno["Std"] <= minor_al) {
        return(list(ELGENOTIPO = "Low_confi", ProbGeno = NA))
      } else if ("Inv" %in% names(gt_prop_geno) && gt_prop_geno["Inv"] >= min_low_con && gt_prop_geno["Inv"] <= minor_al) {
          return(list(ELGENOTIPO = "Low_confi", ProbGeno = NA))
      } else if("Error" %in% names(gt_prop_geno) && gt_prop_geno["Error"] > minor_al) {
        return(list(ELGENOTIPO = "Std/Err", ProbGeno = 1))
      } else {
        return(list(ELGENOTIPO = "NER", ProbGeno = NA, SNP = "Y"))
      }
  } else {
    if (length(names(gt_prop_geno[gt_prop_geno > (1 - minor_al)])) == 1) {
      geno <- names(gt_prop_geno[gt_prop_geno > (1 - minor_al)])
      prob <- round(1 - (0.5 ^ sum(Recons$SignGT == geno)), digits = 3)
      return(list(ELGENOTIPO = geno, ProbGeno = prob))
    } else {
      return(list(ELGENOTIPO = "Failing_code", ProbGeno = NA))
    }
  }
}

result <- determine_genotype(gt_prop_geno, Recons, plotable, gender, minor_al, min_low_con)

ELGENOTIPO <- result$ELGENOTIPO
ProbGeno <- result$ProbGeno
  
FGeno <- as.data.frame(matrix(c("FinalGenotype", ELGENOTIPO, NA, ProbGeno), nrow = 1))
colnames(FGeno) <- colnames(Recons)
rownames(FGeno) <- "FinalGenotype"
Recons <- rbind(Recons, FGeno)

rG <- repGraf
rG[is.na(rG)] <- -40000

hGraph <- 300 + 100*dim(repGraf)[1]
iterGraf <- c("Reference", dupReads)

colReads <- sapply(rownames(repGraf)[1:dim(repGraf)[1]-1], function(x){
  if(Recons[x,"SignGT"] == "Std"){
    "#555555"
  } else if(Recons[x,"SignGT"] == "Inv"){
    "#bbbbbb"
  } else if(Recons[x,"SignGT"] == "Err"){
    "#FF7B6C"
  } else if(Recons[x,"SignGT"] == "AD"){
    "#6CD0FF"
  }
})

colReads <- c(colReads, "#222222")
repGraf$Ri <- sapply(repGraf$Ri, function(x){
  if(!is.na(x)){
    if(x < repP1){
      repP1
    } else {
      x
    }
  } else {
    NA
  }
})

repGraf$Re <- sapply(repGraf$Re, function(x){
  if(!is.na(x)){
    if(x > repP2){
      repP2
    } else {
      x
    }
  } else {
    NA
  }
})

png(filename = paste0(inv, ".png"), height = hGraph, width = 3000)
ggplot(repGraf) +
  geom_rect(aes(xmin = repP1, xmax = repP2, ymin = dim(repGraf)[1]-0.10, ymax = dim(repGraf)[1]+0.10), fill = "#222222", color = NA, linewidth = 1) +
  geom_rect(aes(xmin = Ri, xmax = Re, ymin = seq(length(iterGraf))-0.05, ymax = seq(length(iterGraf))+0.05), fill = colReads, color = NA, linewidth = 1) +
  geom_rect(aes(xmin = Ai, xmax = Ae, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#90af43", color = NA, linewidth = 6) +
  geom_rect(aes(xmin = Bi, xmax = Be, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#932e1f", color = NA, linewidth = 6) +
  geom_rect(aes(xmin = Ci, xmax = Ce, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#e4a21a", color = NA, linewidth = 6) +
  geom_rect(aes(xmin = Di, xmax = De, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#424daa", color = NA, linewidth = 6) +
  theme_classic() + 
  xlim(repP1, repP2 + 1) + 
  ylim(0, dim(repGraf)[1] + 1.5) +
  annotate("text",
           x = (repP1 + repP2) / 2,
           y = dim(repGraf)[1] + 1.5,
           label = paste0(inv, " is ", FGeno$Dist, " with accuracy ", FGeno$ProbGeno),
           size = 30) +
  theme(axis.line = element_blank(), axis.text = element_blank()) 

dev.off()

Recons[is.na(Recons)] <- ""

write.table(Recons, paste0(inv,"_Genotype.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)


##### Checking if BP distinction needed

if(repGraf["Reference", "Ci"] - repGraf["Reference", "Be"] > 150000){
  dAB <- repGraf["Reference", "Be"] - repGraf["Reference", "Ai"]
  dCD <- repGraf["Reference", "De"] - repGraf["Reference", "Ci"]
  
  pasegura <- max(dAB, dCD)
  
  BP1 <- repGraf
  BP2 <- repGraf
  
  ### Separate by BPs
  
  BP1 <- BP1[!is.na(BP1$Ae),]
  BP1 <- BP1[BP1$Ae + pasegura + 100000 > BP1$Re ,]
  BP1 <- BP1[!is.na(BP1$Ri),]
  BP1$Orient <- Recons[rownames(BP1), "SignGT"]
  BP1$Colorinchis <- sapply(BP1$Orient, function(x){
    if(x == "Std"){
      "#555555"
    } else if(x == "Inv"){
      "#bbbbbb"
    } else if(x == "Err"){
      "#FF7B6C"
    } else if(x == "AD"){
      "#6CD0FF"
    }
  })
  
  BP2 <- BP2[!is.na(BP2$Di),]
  BP2 <- BP2[BP2$Di - pasegura - 100000 < BP2$Re ,]
  BP2 <- BP2[!is.na(BP2$Re),]
  BP2$Orient <- Recons[rownames(BP2), "SignGT"]
  BP2$Colorinchis <- sapply(BP2$Orient, function(x){
    if(x == "Std"){
      "#555555"
    } else if(x == "Inv"){
      "#bbbbbb"
    } else if(x == "Err"){
      "#FF7B6C"
    } else if(x == "AD"){
      "#6CD0FF"
    }
  })
  
  BP1_1 <- refPlot$Init
  BP1_2 <- refPlot$Be + 100000
  
  BP2_1 <- refPlot$Ci - 100000
  BP2_2 <- refPlot$End
  
  hBP1 = hGraph <- 300 + 100*dim(BP1)[1]
  hBP2 = hGraph <- 300 + 100*dim(BP2)[1]
  
  refBP <- c(NA, as.numeric(refPlot[1,3:10]), rep(NA, 6))
  
  BP1 <- rbind(BP1, refBP)
  BP2 <- rbind(BP2, refBP)
  
  if(dim(BP1)[1] > 1){
    png(filename = paste0(inv, "_BP1.png"), height = hBP1, width = 3000)
    g1 <- ggplot(BP1) +
      geom_rect(aes(xmin = BP1_1, xmax = BP1_2, ymin = dim(BP1)[1]-0.10, ymax = dim(BP1)[1]+0.10), fill = "#222222", color = NA, linewidth = 1) +
      geom_rect(aes(xmin = Ri, xmax = Re, ymin = seq(dim(BP1)[1])-0.05, ymax = seq(dim(BP1)[1])+0.05), linewidth = 1, fill = BP1$Colorinchis, color = NA) +
      geom_rect(aes(xmin = Ai, xmax = Ae, ymin = seq(dim(BP1)[1])-0.3, ymax = seq(dim(BP1)[1])+0.3), linewidth = 6, fill = "#90af43", color = NA) +
      geom_rect(aes(xmin = Bi, xmax = Be, ymin = seq(dim(BP1)[1])-0.3, ymax = seq(dim(BP1)[1])+0.3), linewidth = 6, fill = "#932e1f", color = NA) +
      geom_rect(aes(xmin = Ci, xmax = Ce, ymin = seq(dim(BP1)[1])-0.3, ymax = seq(dim(BP1)[1])+0.3), linewidth = 6, fill = "#e4a21a", color = NA) +
      geom_rect(aes(xmin = Di, xmax = De, ymin = seq(dim(BP1)[1])-0.3, ymax = seq(dim(BP1)[1])+0.3), linewidth = 6, fill = "#424daa", color = NA) +
      theme_classic() + 
      xlim(BP1_1,BP1_2+1) + 
      ylim(0,dim(BP1)[1]+1.5) +
      annotate("text",
           x = (repP1 + repP2) / 2,
           y = dim(repGraf)[1] + 1.5,
           label = paste0(inv, " is ", FGeno$Dist, " with accuracy ", FGeno$ProbGeno),
           size = 30) +
  theme(axis.line = element_blank(), axis.text = element_blank())
    
    print(g1)
    
    dev.off()
  }
  
  if(dim(BP2)[1] > 1){
    png(filename = paste0(inv, "_BP2.png"), height = hBP2, width = 3000)
    g2 <- ggplot(BP2) +
      geom_rect(aes(xmin = BP2_1, xmax = BP2_2, ymin = dim(BP2)[1]-0.10, ymax = dim(BP2)[1]+0.10), fill = "#222222", color = NA, linewidth = 1) +
      geom_rect(aes(xmin = Ri, xmax = Re, ymin = seq(dim(BP2)[1])-0.05, ymax = seq(dim(BP2)[1])+0.05), linewidth = 1, fill = BP2$Colorinchis, color = NA) +
      geom_rect(aes(xmin = Ai, xmax = Ae, ymin = seq(dim(BP2)[1])-0.3, ymax = seq(dim(BP2)[1])+0.3), linewidth = 6, fill = "#90af43", color = NA) +
      geom_rect(aes(xmin = Bi, xmax = Be, ymin = seq(dim(BP2)[1])-0.3, ymax = seq(dim(BP2)[1])+0.3), linewidth = 6, fill = "#932e1f", color = NA) +
      geom_rect(aes(xmin = Ci, xmax = Ce, ymin = seq(dim(BP2)[1])-0.3, ymax = seq(dim(BP2)[1])+0.3), linewidth = 6, fill = "#e4a21a", color = NA) +
      geom_rect(aes(xmin = Di, xmax = De, ymin = seq(dim(BP2)[1])-0.3, ymax = seq(dim(BP2)[1])+0.3), linewidth = 6, fill = "#424daa", color = NA) +
      theme_classic() + 
      xlim(BP2_1,BP2_2+1) + 
      ylim(0,dim(BP2)[1]+1.5) +
      annotate("text",
           x = (repP1 + repP2) / 2,
           y = dim(repGraf)[1] + 1.5,
           label = paste0(inv, " is ", FGeno$Dist, " with accuracy ", FGeno$ProbGeno),
           size = 30) +
  theme(axis.line = element_blank(), axis.text = element_blank())
    
    print(g2)
    
    dev.off()
  }
  
}


####### Add info about DelIns

CalcDist <- as.data.frame(row.names(repGraf))
colnames(CalcDist) <- "Read"
rownames(CalcDist) <- CalcDist$Read

CalcDist$PosA <- sapply(CalcDist$Read, function(x){
  if(!is.na(repGraf[x, "Ae"])){
    repGraf[x, "Ai"] + (repGraf[x, "Ae"] - repGraf[x, "Ai"])/2
    
  } else {
    NA
  }
  
})

CalcDist$PosB <- sapply(CalcDist$Read, function(x){
  if(!is.na(repGraf[x, "Be"])){
    repGraf[x, "Bi"] + (repGraf[x, "Be"] - repGraf[x, "Bi"])/2
    
  } else {
    NA
  }
  
})

CalcDist$PosC <- sapply(CalcDist$Read, function(x){
  if(!is.na(repGraf[x, "Ce"])){
    repGraf[x, "Ci"] + (repGraf[x, "Ce"] - repGraf[x, "Ci"])/2
    
  } else {
    NA
  }
  
})

CalcDist$PosD <- sapply(CalcDist$Read, function(x){
  if(!is.na(repGraf[x, "De"])){
    repGraf[x, "Di"] + (repGraf[x, "De"] - repGraf[x, "Di"])/2
    
  } else {
    NA
  }
  
})

inv_ance_path <- "../../../Infor/Inv_ancestral.txt"

if (file.exists(inv_ance_path)) {
  Inv_ance <- read.table("../../../Infor/Inv_ancestral.txt")
  if(inv %in% Inv_ance$V1){
    temp <- CalcDist$PosB
    CalcDist$PosB <- CalcDist$PosC
    CalcDist$PosC <- temp
  }
}

CalcDist$GT <- sapply(CalcDist$Read, function(x){
  if(x == "Reference"){
    "Std"
  } else {
    Recons[x, "SignGT"]
  }
  
})

CalcDist <- CalcDist[CalcDist$GT %in% c("Std", "Inv", "AD"),]

CalcDist$BP1 <- sapply(CalcDist$Read, function(x){
  if (CalcDist[x, "GT"] == "Std"){
    if (!is.na(CalcDist[x, "PosA"]) && !is.na(CalcDist[x, "PosB"])){
      CalcDist[x, "PosB"] - CalcDist[x, "PosA"]
    } else if (!is.na(CalcDist[x, "PosA"]) && !is.na(CalcDist[x, "PosC"])){
      CalcDist[x, "PosC"] - CalcDist[x, "PosA"]
    } else {
      NA
    }
  } else if (CalcDist[x, "GT"] == "Inv"){
    if (!is.na(CalcDist[x, "PosC"]) && !is.na(CalcDist[x, "PosA"])){
      CalcDist[x, "PosC"] - CalcDist[x, "PosA"]
    } else if (!is.na(CalcDist[x, "PosB"]) && !is.na(CalcDist[x, "PosA"])){
      CalcDist[x, "PosB"] - CalcDist[x, "PosA"]
    } else {
      NA
    }
  } else {
    NA
  } 
})


CalcDist$BP2 <- sapply(CalcDist$Read, function(x){
  if (CalcDist[x, "GT"] == "Std"){
    if (!is.na(CalcDist[x, "PosC"]) && !is.na(CalcDist[x, "PosD"])){
      CalcDist[x, "PosD"] - CalcDist[x, "PosC"]
    } else if (!is.na(CalcDist[x, "PosB"]) && !is.na(CalcDist[x, "PosD"])){
      CalcDist[x, "PosD"] - CalcDist[x, "PosB"]
    } else {
      NA
    }
  } else if (CalcDist[x, "GT"] == "Inv"){
    if (!is.na(CalcDist[x, "PosD"]) && !is.na(CalcDist[x, "PosB"])){
      CalcDist[x, "PosD"] - CalcDist[x, "PosB"]
    } else if (!is.na(CalcDist[x, "PosC"]) && !is.na(CalcDist[x, "PosD"])){
      CalcDist[x, "PosD"] - CalcDist[x, "PosC"]
    } else {
      NA
    }
  } else {
    NA
  }
})

CalcDist$BC <- sapply(CalcDist$Read, function(x){
  abs(CalcDist[x, "PosC"] - CalcDist[x, "PosB"]) 
})

CalcDist$AD <- sapply(CalcDist$Read, function(x){
  abs(CalcDist[x, "PosD"] - CalcDist[x, "PosA"]) 
})

InvDist <- InvDist[inv, ]
RefDist <- CalcDist["Reference", ]

StdDist <- InvDist[inv, ]
colnames(StdDist) <- c("Inv","BP1","BP2","BC","AD")
StdDist <- StdDist[,!names(StdDist) %in% c("BC","AD")]
StdDist$BC <- RefDist$BC
StdDist$AD <- RefDist$AD

if(RefDist$PosB > RefDist$PosC){
  RefDist <- StdDist
} 

CalcDist <- CalcDist[CalcDist$Read != "Reference", ]

if (dim(InvDist)[1] == 1 && dim(repGraf)[1] > 0){
  
  ### DelIns inside IRs
  CalcDist$DifBP1_Std <- sapply(CalcDist$Read, function (x){
    GT <- CalcDist[x, "GT"]
    PosA <- CalcDist[x, "PosA"]
    PosB <- CalcDist[x, "PosB"]
    PosC <- CalcDist[x, "PosC"]
    
    if (GT == "Std" && !is.na(PosA) && !is.na(PosB)){
      CalcDist[x, "BP1"] - RefDist[ ,"BP1"]
    } else if (GT == "Std" && !is.na(PosA) && !is.na(PosC)){
      CalcDist[x, "BP1"] - (RefDist[ ,"BP1"] + RefDist[ ,"BC"])
    } else {
      NA
    }
  })
  
  CalcDist$DifBP1_Inv <- sapply(CalcDist$Read, function (x){
    GT <- CalcDist[x, "GT"]
    PosA <- CalcDist[x, "PosA"]
    PosB <- CalcDist[x, "PosB"]
    PosC <- CalcDist[x, "PosC"]
    
    if (GT == "Inv" && !is.na(PosA) && !is.na(PosC)){
      CalcDist[x, "BP1"] - InvDist[ , "DistAC"]
    } else if (GT == "Inv" && !is.na(PosA) && !is.na(PosB)){
      CalcDist[x, "BP1"] - InvDist[ , "DistAB"]
    } else {
      NA
    }
  })
  
  
  CalcDist$DifBP2_Std <- sapply(CalcDist$Read, function (x){
    GT <- CalcDist[x, "GT"]
    PosD <- CalcDist[x, "PosD"]
    PosB <- CalcDist[x, "PosB"]
    PosC <- CalcDist[x, "PosC"]
    
    if (GT == "Std" && !is.na(PosC) && !is.na(PosD)){
      CalcDist[x, "BP2"] - RefDist[ ,"BP2"]
    } else if (GT == "Std" && !is.na(PosB) && !is.na(PosD)){
      CalcDist[x, "BP2"] - (RefDist[ ,"BP2"] + RefDist[ ,"BC"])
    } else {
      NA
    }
  })
  
  CalcDist$DifBP2_Inv <- sapply(CalcDist$Read, function (x){
    GT <- CalcDist[x, "GT"]
    PosD <- CalcDist[x, "PosD"]
    PosB <- CalcDist[x, "PosB"]
    PosC <- CalcDist[x, "PosC"]
    
    if (GT == "Inv" && !is.na(PosB) && !is.na(PosD)){
      CalcDist[x, "BP2"] - InvDist[ , "DistBD"]
    } else if (GT == "Inv" && !is.na(PosC) && !is.na(PosD)){
      CalcDist[x, "BP2"] - InvDist[ , "DistCD"]
    } else {
      NA
    }
  })
  
  CalcDist$DifBC <- sapply(CalcDist$Read, function (x){
    CalcDist[x, "BC"] - RefDist[ ,"BC"]
  })
  
  
  CalcDist$DifAD <- sapply(CalcDist$Read, function (x){
    GT <- CalcDist[x, "GT"]
    PosD <- CalcDist[x, "PosD"]
    PosA <- CalcDist[x, "PosA"]
    
    if (GT == "AD"){
      CalcDist[x, "AD"] - RefDist[ ,"AD"]
    } else {
      NA
    }
  })
  
  CalcDist$SVBP1_Std <- NA
  for(x in CalcDist$Read){
    var <- abs(CalcDist[x,"BP1"])*diff_exp
    res <- (abs(CalcDist[x, "DifBP1_Std"]) > var) && (abs(CalcDist[x, "DifBP1_Std"]) > min_sv)
    CalcDist[x,"SVBP1_Std"] <- res
  }
  
  CalcDist$SVBP1_Inv <- NA
  for(x in CalcDist$Read){
    var <- abs(CalcDist[x,"BP1"])*diff_exp
    res <- (abs(CalcDist[x, "DifBP1_Inv"]) > var) && (abs(CalcDist[x, "DifBP1_Inv"]) > min_sv)
    CalcDist[x,"SVBP1_Inv"] <- res
  }
  
  
  CalcDist$SVBP2_Std <- NA
  for(x in CalcDist$Read){
    var <- abs(CalcDist[x,"BP2"])*diff_exp
    res <- (abs(CalcDist[x, "DifBP2_Std"]) > var) && (abs(CalcDist[x, "DifBP2_Std"]) > min_sv)
    CalcDist[x,"SVBP2_Std"] <- res
  }
  
  CalcDist$SVBP2_Inv <- NA
  for(x in CalcDist$Read){
    var <- abs(CalcDist[x,"BP2"])*diff_exp
    res <- (abs(CalcDist[x, "DifBP2_Inv"]) > var) && (abs(CalcDist[x, "DifBP2_Inv"]) > min_sv)
    CalcDist[x,"SVBP2_Inv"] <- res
  }
  
  CalcDist$SVBC <- NA
  for(x in CalcDist$Read){
    var <- abs(CalcDist[x,"BC"])*diff_exp
    res <- (abs(CalcDist[x, "DifBC"]) > var) && (abs(CalcDist[x, "DifBC"]) > min_sv)
    CalcDist[x,"SVBC"] <- res
  }
  
  CalcDist$SVAD <- NA
  for (x in CalcDist$Read){
    var <- abs(CalcDist[x, "AD"])*diff_exp
    res <- (abs(CalcDist[x, "DifAD"]) > var) && (abs(CalcDist[x, "DifAD"]) > min_sv)
    CalcDist[x,"SVAD"] <- res
  }
  
  ### Calculating the proportions of the inversion orientation
  
  mtxOr <- as.data.frame(unique(CalcDist$GT))
  colnames(mtxOr) <- "Orientation"
  
  if(dim(mtxOr)[1] > 1){
    mtxOr$Quant <- sapply(mtxOr$Orientation, function(x){
      sum(CalcDist$GT == x)
    })
    mtxOr <- mtxOr[order(mtxOr$Quant),]
    
    gt_freq <- table(CalcDist$GT)
    gt_prop <- prop.table(gt_freq)
    minor_gt <- names(gt_prop[gt_prop < minor_al])
    CalcDist <- CalcDist[!CalcDist$GT %in% minor_gt, ]  
  }
  
  Props <- as.data.frame(c("SVBP1_Std", "SVBP2_Std", "SVBP1_Inv", "SVBP2_Inv", "SVBC", "SVAD"))
  colnames(Props) <- "Region"
  rownames(Props) <- Props$Region
  
  ### Adding proportions of trues
  
  pr <- c(sum(na.omit(CalcDist$SVBP1_Std))/(length(na.omit(CalcDist$SVBP1_Std))+1e-6),
          sum(na.omit(CalcDist$SVBP2_Std))/(length(na.omit(CalcDist$SVBP2_Std))+1e-6),
          sum(na.omit(CalcDist$SVBP1_Inv))/(length(na.omit(CalcDist$SVBP1_Inv))+1e-6),
          sum(na.omit(CalcDist$SVBP2_Inv))/(length(na.omit(CalcDist$SVBP2_Inv))+1e-6),
          sum(na.omit(CalcDist$SVBC))/(length(na.omit(CalcDist$SVBC))+1e-6),
          sum(na.omit(CalcDist$SVAD))/(length(na.omit(CalcDist$SVAD))+1e-6))
  
  Props$Ratio <- pr
  
  CalcDist$SVBP1_Std <- sapply(CalcDist$SVBP1_Std, function(x){
    if (is.na(x)){
      F
    } else {
      x
    }
  })
  
  CalcDist$SVBP2_Std <- sapply(CalcDist$SVBP2_Std, function(x){
    if (is.na(x)){
      F
    } else {
      x
    }
  })
  
  
  CalcDist$SVBP1_Inv <- sapply(CalcDist$SVBP1_Inv, function(x){
    if (is.na(x)){
      F
    } else {
      x
    }
  })
  
  CalcDist$SVBP2_Inv <- sapply(CalcDist$SVBP2_Inv, function(x){
    if (is.na(x)){
      F
    } else {
      x
    }
  })
  
  CalcDist$SVBC <- sapply(CalcDist$SVBC, function(x){
    if (is.na(x)){
      F
    } else {
      x
    }
  })
  
  CalcDist$SVAD <- sapply(CalcDist$SVAD, function(x){
    if (is.na(x)){
      F
    } else {
      x
    }
  })
  
  write.table(CalcDist, file = "Main_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  Props$SV <- sapply(Props$Ratio, function(x){
    if(x < 0.25){
      "HomoNoSV"
    } else if(x >= 0.25 && x <= 0.75){
      "HetSV"
    } else if(x > 0.75){
      "HomoSV"
    }
  })
  
  checkGroupsHomoSV <- function(dfsv, bp){
    redo <- F
    namedif <- gsub("SV", "Dif", bp)
    processDif <- dfsv[,c("Read",namedif)]
    processDif <- processDif[order(processDif[[namedif]]),]
    rownames(processDif) <- as.character(seq(dim(processDif)[1]))
    
    jumps <- c(0)
    for(i in seq(length(processDif[[namedif]])-1)){
      jumps <- c(jumps, processDif[i+1, namedif]-processDif[i,namedif])
    }
    
    processDif$Jumps <- jumps
    
    jumps <- sort(jumps)
    maxJump <- max(jumps)
    jumps <- jumps[1:length(jumps)-1]
    
    sepSVs <- F
    
    if (length(jumps) > 1) {
      if (maxJump > mean(jumps) + 7 * sd(jumps) && maxJump > 1000) {
        sepSVs <- TRUE
      }
    }
    
    
    if(sepSVs){
      cutIdx <- which.max(processDif$Jumps)
      
      firstSV <- processDif[1:cutIdx-1,]
      secondSV <- processDif[cutIdx:dim(processDif)[1],]
      
      DiffMtx <- as.data.frame(matrix(c(1,2),nrow = 2))
      colnames(DiffMtx) <- "Table"
      
      DiffMtx$NumComp <- NA
      DiffMtx[1,2] <- dim(firstSV)[1]
      DiffMtx[2,2] <- dim(secondSV)[1]
      
      propMtx <- DiffMtx[1,2]/sum(DiffMtx[,2])
      
      if(propMtx < minor_al){
        dfsv <- secondSV
        redo <- T
      } else if(propMtx > (1 - minor_al)){
        dfsv <- firstSV
        redo <- T
      }
      
    }
    
    checkedList <- list(df = dfsv,
                        redo = redo,
                        pD = processDif,
                        sepSVs = sepSVs,
                        nd = namedif)
    
    
    return(checkedList)
  }
  
  ResFinSV <- new.env()
  
  for(bp in c("SVBP1_Std", "SVBP2_Std", "SVBP1_Inv", "SVBP2_Inv", "SVBC", "SVAD")){
    svdet <- Props[bp,"SV"]
    
    if(dim(CalcDist[CalcDist[[bp]], ])[1] < 2){
      next
    }
    
    if(svdet == "HetSV"){
      dfsv <- CalcDist[CalcDist[[bp]], ]
      namedif <- gsub("SV", "Dif", bp)
      
      ### Double checking errors
      ## First by GT
      
      mtxSVOr <- as.data.frame(unique(dfsv$GT))
      colnames(mtxSVOr) <- "Orientation"
      if(dim(mtxSVOr)[1] > 1){
        mtxSVOr$Quant <- sapply(mtxSVOr$Orientation, function(x){
          sum(dfsv$GT == x)
        })
        mtxSVOr <- mtxSVOr[order(mtxSVOr$Quant),]
        
        if(mtxSVOr[1,2]/sum(mtxSVOr[,2]) < minor_al){
          dfsv <- dfsv[dfsv$GT == mtxSVOr[2,1],]
        }
        
      }
      
      # ## Then by 2*sd()
      
      dfAux <- dfsv[,c("Read", namedif)]
      minDet <- mean(dfAux[,namedif])-2*sd(dfAux[,namedif])
      maxDet <- mean(dfAux[,namedif])+2*sd(dfAux[,namedif])
      
      rmrd <- sapply(dfAux[,namedif], function(x){
        if(x < minDet || x > maxDet){
          F
        } else {
          T
        }
        
      })
      
      dfsv <- dfsv[rmrd,]
      
      
      mean <- round(mean(dfsv[[namedif]]),0)
      range <- paste0(round(min(abs(dfsv[[namedif]])),0),"-",
                      round(max(abs(dfsv[[namedif]])),0))
      type <- "Ins"
      if(mean < 0){
        type <- "Del"
      }
      dfsv$Orientation <- dfsv$GT
      dfsv$VarDet <- paste(type, as.character(abs(mean)), sep = "_")
      dfsv$Range <- range
      ResFinSV[[bp]] <- dfsv[,c("Read", "Orientation", "VarDet", "Range")]
      
    } else if(svdet == "HomoSV"){
      dfsv <- CalcDist[CalcDist[[bp]], ]
      
      redo <- T
      
      while(redo){
        res <- checkGroupsHomoSV(dfsv,bp)
        dfsv <- res$df
        processDif <- res$pD
        redo <- res$redo
        sepSVs <- res$sepSVs
        namedif <- res$nd
      }
      
      if(sepSVs){
        cutIdx <- which.max(processDif$Jumps)
        
        firstSV <- processDif[1:cutIdx-1,]
        secondSV <- processDif[cutIdx:dim(processDif)[1],]
        
        mean1 <- round(mean(firstSV[[namedif]]),0)
        mean2 <- round(mean(secondSV[[namedif]]),0)
        
        if(abs(mean1) >= abs(mean2)){
          rest <- abs(mean1)-abs(mean2)
          div <- rest/(mean1)
          if (div <= splitting_svs){
            mean1 <- mean1
            mean2 <- mean1
          } else {
            mean1 <- mean1
            mean2 <- mean2
          }
        } else if (abs(mean2) >= abs(mean1)){
          rest <- abs(mean2)-abs(mean1)
          div <- rest/(mean2)
          if (div <= splitting_svs){
            mean1 <- mean2
            mean2 <- mean2
          } else {
            mean1 <- mean1
            mean2 <- mean2
          }
        }
        
        type1 <- "Ins"
        type2 <- "Ins"
        
        if(mean1 < 0){
          type1 <- "Del"
        }
        if(mean2 < 0){
          type2 <- "Del"
        }
        
        firstSV$Orientation <- sapply(firstSV$Read, function(x){
          CalcDist[x,"GT"]
        })
        secondSV$Orientation <- sapply(secondSV$Read, function(x){
          CalcDist[x,"GT"]
        })
        
        firstSV$VarDet <- paste(type1, abs(mean1), sep = "_")
        secondSV$VarDet <- paste(type2, abs(mean2), sep = "_")
        
        firstSV$Range <- paste0(round(min(abs(firstSV[[namedif]])),0),"-",
                                round(max(abs(firstSV[[namedif]])),0))
        secondSV$Range <- paste0(round(min(abs(secondSV[[namedif]])),0),"-",
                                 round(max(abs(secondSV[[namedif]])),0))
        
        firstSV <- rbind(firstSV, secondSV)
        rownames(firstSV) <- firstSV$Read
        
        ResFinSV[[bp]] <- firstSV[,c("Read", "Orientation", "VarDet", "Range")]
        
      } else {
        mean <- round(mean(processDif[[namedif]]),0)
        type <- "Ins"
        
        range <- paste0(round(min(abs(processDif[[namedif]])),0),"-",
                        round(max(abs(processDif[[namedif]])),0))
        
        if(mean < 0){
          type <- "Del"
        }
        processDif$Orientation <- sapply(processDif$Read, function(x){
          CalcDist[x,"GT"]
        })
        processDif$VarDet <- paste(type, as.character(abs(mean)), sep = "_")
        processDif$Range <- range
        
        ResFinSV[[bp]] <- processDif[,c("Read", "Orientation","VarDet", "Range")]
        
      }
      
    }
    
  }
  
  if(length(ls(ResFinSV)) != 0){
    for(bp in ls(ResFinSV)){
      write.table(ResFinSV[[bp]], paste0("DetSV_",inv,"_",bp,".txt"),
                  row.names = F, col.names = T, quote = F, sep = "\t")
    }
    
    system(paste0("cat DetSV_* > DetSV_", inv, "_AllBPs.txt"))
    
  }
  
  
  
  
}
