#########################
#
# FUNCTIONS
#
#########################
##-----------------------------------------------------##
# onlyNAs_simplifier
##-----------------------------------------------------##
onlyNAs_simplifier <- function(x, separator = ";", NAvalue = "NA"){
  
  cell <- x
  parsed_cell <- unlist(strsplit(cell, separator))
  parsed_cell <- paste0(unique(parsed_cell), collapse = separator)
  
  if(parsed_cell == NAvalue){
    cell <- NA
  } 
  
  return(cell)
}

##-----------------------------------------------------##
# geneBiotype_creator
##-----------------------------------------------------##
geneBiotype_creator <- function(feelnc_file){
  
  type <- switch(EXPR = feelnc_file["type"], intergenic = "linc", genic = "lncg", "ERROR")
  subtype <- switch(EXPR = feelnc_file["subtype"], containing = "Cont", convergent = "Conv", divergent = "Divg", nested = "Nest", overlapping = "Ovlp", same_strand = "SS", "ERROR") 
  location <- switch(EXPR = feelnc_file["location"], downstream = "dw", exonic = "ex", intronic = "in", upstream = "up", "ERROR")
  direction <- switch(EXPR = feelnc_file["direction"], antisense = "AS", sense = "SS", "ERROR")
  
  if(type == "linc"){
    res <- paste0(type, subtype, ifelse(subtype == "SS", location, ""))
  } else if(type == "lncg") {
    res <- paste0(type, direction, location, subtype)
  }
  
  return(res)
}

##-----------------------------------------------------##
# FEELncClass_simplificator
##-----------------------------------------------------##
FEELncClass_simplificator <- function(x){
  
  x <- unlist(strsplit(x, ";"))
  x[grepl("lncgSSex", x)] <- "lncgSSex"
  x[grepl("lncgSSin", x)] <- "lncgSSin"
  x[grepl("lncgASex", x)] <- "lncgASex"
  x[grepl("lncgASin", x)] <- "lncgASin"
  x[grepl("lincDivg", x)] <- "lincDivg"
  x[grepl("lincConv", x)] <- "lincConv"
  x[grepl("lincSS", x)] <- "lincSS"
  x <- paste0(unique(x), collapse = ";")
  
  return(x)
}

##-----------------------------------------------------##
# nbCategories_calculator
##-----------------------------------------------------##
nbCategories_calculator <- function(x, pattern_to_remove){
  
  a <- strsplit(x, ";")
  b <- lapply(a, unique)
  b <- lapply(b, function(x)list_partRemover(x, pattern_to_remove = pattern_to_remove))
  d <- unlist(lapply(b, function(x){length(x) == 1}))
  nbCategories <- ifelse(d, "1", "n")
  
  return(nbCategories)
}

##-----------------------------------------------------##
# list_partRemover
##-----------------------------------------------------##
list_partRemover <- function(x, pattern_to_remove){
  
  isthere_pattern_to_remove <- any(grepl(pattern_to_remove, x))
  
  if(isthere_pattern_to_remove){
    
    x <- x[- which(grepl(pattern_to_remove, x))]
  }
  
  return(x)
}

#########################
#
# DATA
#
#########################
##-----------------------------------------------------##
# Classified genes
##-----------------------------------------------------##
feelnc_results_lncRNArelativeTo_mRNA <- read.table("./FEELnc/1_Ensembl_20180731_all_24881g_38118t_gg5_v87v93_Ensemblv87eqv93_mRNAonly_FEELnc_classes_lncRNA_relativeTo_mRNA.txt",
                                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
feelnc_results_lncRNArelativeTo_mRNA <- feelnc_results_lncRNArelativeTo_mRNA[feelnc_results_lncRNArelativeTo_mRNA$isBest == 1, ] # selection of the is_Best = 1 only

##-----------------------------------------------------##
# Unclassified genes
##-----------------------------------------------------##
feelnc_unclassified_lnc_relativeTo_pcg <- scan("./FEELnc/1_Ensembl_20180731_all_24881g_38118t_gg5_v87v93_Ensemblv87eqv93_mRNAonly_FEELnc_log_lncRNA_relativeTo_mRNA.log",
                                               what = "character")
pos_firsttpID <- grep("GALT", feelnc_unclassified_lnc_relativeTo_pcg)[1] # theorical position of the first transcript id in the vector : 51. To check if re-used
feelnc_unclassified_lnc_relativeTo_pcg <- feelnc_unclassified_lnc_relativeTo_pcg[- (1:(pos_firsttpID - 1))]

#########################
#
# TRANSCRIPT LEVEL
#
#########################
##-----------------------------------------------------##
# Empty transcript-level dataframe creation
##-----------------------------------------------------##
annotation_transcript <- data.frame(gnId = # <gtf_transcripts$gene_id>,
                                    gnName = # <gtf_transcripts$gene_name>,
                                    tpId = # <gtf_transcripts$transcript_id>,
                                    feelLncPcgClassName = NA, 
                                    feelLncPcgGnId = NA, feelLncPcgGnName = NA, feelLncPcgGnDist = NA,
                                    stringsAsFactors = FALSE)

##-----------------------------------------------------##
# Completion of the transcript-level dataframe
##-----------------------------------------------------##
feelnc_results_lncRNArelativeTo_mRNA$biotype <- apply(feelnc_results_lncRNArelativeTo_mRNA, 1, geneBiotype_creator)

annotation_transcript$feelLncPcgClassName <- feelnc_results_lncRNArelativeTo_mRNA$biotype[match(annotation_transcript$tpId, feelnc_results_lncRNArelativeTo_mRNA$lncRNA_transcript)]
annotation_transcript$feelLncPcgClassName[annotation_transcript$tpId %in% feelnc_unclassified_lnc_relativeTo_pcg] <- "unclassified"

annotation_transcript$feelLncPcgGnId <- feelnc_results_lncRNArelativeTo_mRNA$partnerRNA_gene[match(annotation_transcript$tpId, feelnc_results_lncRNArelativeTo_mRNA$lncRNA_transcript)]
annotation_transcript$feelLncPcgGnName <- annotation_transcript$gnName[match(annotation_transcript$feelLncPcgGnId, annotation_transcript$gnId)]
annotation_transcript$feelLncPcgGnDist <- feelnc_results_lncRNArelativeTo_mRNA$distance[match(annotation_transcript$tpId, feelnc_results_lncRNArelativeTo_mRNA$lncRNA_transcript)]

#########################
#
# GENE LEVEL
#
#########################
##-----------------------------------------------------##
# Gene-level dataframe creation
##-----------------------------------------------------##
annotation_gene <- data.frame(gnId = # <gtf_genes$gene_id>,
                              nbTp = NA,
                              feelLncPcgClassName = NA, 
                              feelLncPcgClassType = NA,
                              feelLncPcgGnId = NA, feelLncPcgGnName = NA, feelLncPcgGnDist = NA,
                              feelLncPcgClassNameByTp = NA, 
                              feelLncPcgGnIdByTp = NA, feelLncPcgGnNameByTp = NA, feelLncPcgGnDistByTp = NA,
                              stringsAsFactors = FALSE)

for (i in 1:nrow(annotation_gene)){
  
  if(i %% 5000 == 0){cat(i, "/", nrow(annotation_gene), "\n")}
  gene_id <- annotation_gene$gnId[i]
  lineNumber <- which(annotation_transcript$gnId %in% gene_id)
  
  annotation_gene$feelLncPcgClassNameByTp[i] <- paste0(annotation_transcript$feelLncPcgClassName[lineNumber], collapse = ";")
  annotation_gene$feelLncPcgGnIdByTp[i] <- paste0(annotation_transcript$feelLncPcgGnId[lineNumber], collapse = ";")
  annotation_gene$feelLncPcgGnNameByTp[i] <- paste0(annotation_transcript$feelLncPcgGnName[lineNumber], collapse = ";")
  annotation_gene$feelLncPcgGnDistByTp[i] <- paste0(annotation_transcript$feelLncPcgGnDist[lineNumber], collapse = ";")
  
  a <- annotation_transcript$tpId[lineNumber]
  annotation_gene$nbTp[i] <- length(a)
}

##-----------------------------------------------------##
# Order of classification of the subtypes
##-----------------------------------------------------##
reclassificationOrder <- data.frame(classification = a <- c("lncgSSexNest", "lncgSSexCont", "lncgSSexOvlp", 
                                                            "lncgSSinNest", "lncgSSinCont", "lncgSSinOvlp",
                                                            "lncgASexNest", "lncgASexCont", "lncgASexOvlp", 
                                                            "lncgASinNest", "lncgASinCont", "lncgASinOvlp",
                                                            "lincDivg", "lincSSup", "lincSSdw", "lincConv", 
                                                            "unclassified", "NA"),
                                    index = seq_len(length.out = length(a)),
                                    stringsAsFactors = FALSE)

reclassificationOrder$index[reclassificationOrder$classification == "lincSSdw"] <- reclassificationOrder$index[reclassificationOrder$classification == "lincSSup"]

##-----------------------------------------------------##
# Classification of each gene
##-----------------------------------------------------##
annotation_gene[, "feelLncPcgClassNameByTp"] <- apply(annotation_gene[, "feelLncPcgClassNameByTp", drop = FALSE], 1, onlyNAs_simplifier)
annotation_gene[, "feelLncPcgGnIdByTp"] <- apply(annotation_gene[, "feelLncPcgGnIdByTp", drop = FALSE], 1, onlyNAs_simplifier)
annotation_gene[, "feelLncPcgGnNameByTp"] <- apply(annotation_gene[, "feelLncPcgGnNameByTp", drop = FALSE], 1, onlyNAs_simplifier)
annotation_gene[, "feelLncPcgGnDistByTp"] <- apply(annotation_gene[, "feelLncPcgGnDistByTp", drop = FALSE], 1, onlyNAs_simplifier)

for (i in 1:nrow(annotation_gene)){
  
  if(i %% 5000 == 0){cat(i, "/", nrow(annotation_gene), "\n")}
  
  tmp_lnc <- data.frame(
    class = unlist(strsplit(annotation_gene$feelLncPcgClassNameByTp[i], ";")),
    gnId = unlist(strsplit(annotation_gene$feelLncPcgGnIdByTp[i], ";")),
    gnNm = unlist(strsplit(annotation_gene$feelLncPcgGnNameByTp[i], ";")),
    gnDist = unlist(strsplit(annotation_gene$feelLncPcgGnDistByTp[i], ";")),
    order = NA,
    stringsAsFactors = FALSE)
  
  tmp_lnc$order <- reclassificationOrder$index[match(tmp_lnc$class, reclassificationOrder$classification)]
  tmp_lnc <- tmp_lnc[order(tmp_lnc$order, tmp_lnc$gnDist), ]
  
  annotation_gene$feelLncPcgClassName[i] <- tmp_lnc$class[1]
  annotation_gene$feelLncPcgGnId[i] <- tmp_lnc$gnId[1]
  annotation_gene$feelLncPcgGnName[i] <- tmp_lnc$gnNm[1]
  annotation_gene$feelLncPcgGnDist[i] <- tmp_lnc$gnDist[1]
}

##-----------------------------------------------------##
# Typology of the FEELnc classification
##-----------------------------------------------------##
nbTp <- ifelse(annotation_gene$nbTp == 1, "1", "n")

nbClass_lncRelativeToPcg <- nbCategories_calculator(annotation_gene$feelLncPcgClassNameByTp, pattern_to_remove = "unclassified")
nbGenes_lncRelativeToPcg <- nbCategories_calculator(annotation_gene$feelLncPcgGnIdByTp, pattern_to_remove = "NA")

feelLncPcgClassType <- paste0(nbTp, ".", nbClass_lncRelativeToPcg, ".", nbGenes_lncRelativeToPcg)

annotation_gene$feelLncPcgClassType <- feelLncPcgClassType
annotation_gene$feelLncPcgClassType[is.na(annotation_gene$feelLncPcgGnIdByTp)] <- "unclassified"

##-----------------------------------------------------##
# Modification of the nn1, nnn and n1n types
##-----------------------------------------------------##
pos_nn1 <- grepl("n.n.1", annotation_gene$feelLncPcgClassType)
simplifiedClasses_nn1 <- apply(annotation_gene[pos_nn1, "feelLncPcgClassNameByTp", drop = FALSE], 1, FEELncClass_simplificator)
simplifiedClasses_nn1 <- sapply(strsplit(simplifiedClasses_nn1, ";"), function(x){x <- x[1]})
annotation_gene[pos_nn1, "feelLncPcgClassName"] <- simplifiedClasses_nn1

pos_nnn <- grepl("n.n.n", annotation_gene$feelLncPcgClassType)
pos_n1n <- grepl("n.1.n", annotation_gene$feelLncPcgClassType)

annotation_gene$feelLncPcgClassName[pos_nnn] <- paste0(annotation_gene$feelLncPcgClassName[pos_nnn], "_n.n.n")
annotation_gene$feelLncPcgClassName[pos_n1n] <- paste0(annotation_gene$feelLncPcgClassName[pos_n1n], "_n.1.n")
