# Author: Juraj Michalik
# Date: 03/04/2021

# Script containing various functions to import Mixcr VDJ in format similar to 10X and work with them
# -> import VDJ from mixcr and formatting them to format compatible with 10X
# -> combining VDJ from mixcr and 10X
require(tidyverse)
require(data.table)

##############################
# Mixcr GEX VDJ              #  
##############################

# Transforms and adds GEX (Mixcr) VDJ data to 10X as a preparation for insertion of all data
# Mixcr data are not used to filter VDJ doublets, so it is necessary to get only the maximum number of them
# MINUS those that are already found by 10X

# extracts gene_name from pattern
# - str - pattern of mixcr genes ordered by decreasing score
extract_gene_name_from_mixcr <- function(str){
  if(is.na(str) | str == ''){
    return('')
  }
  first_val <- unlist(strsplit(str, ','))[1] # the first one is the one with best score we want
  first_val <- unlist(strsplit(first_val, '\\('))[1] # remove the score
  return(first_val)
}


# rebuilds GEX VDJ table into format employed by 10X VDJ for given chain type
# - GEX.VDJ.dat - output file from Mixcr for specific chain
# - chain - type of chain in question
# - extend.barcode - adds '-1' at the end of barcode
# - vdj.data - 10X VDJ data, this is only necessary for 
# determining the version of cellranger that was used (<=4 or 5>=)
build.GEX.VDJ.table <- function(GEX.VDJ.dat, 
                                chain, 
                                extend.barcode = TRUE,
                                vdj.data = NULL){
  
  null.value <- 'None' # Cellranger 4
  if(any('' %in% vdj.data$d_gene)){
    null.value <- '' # Cellranger 5
  }
  
  true.value <- 'True' # Cellranger 4
  if(any('true' %in% vdj.data$productive)){
    true.value <- 'true' # Cellranger 5
  }  
  
  barcode <- GEX.VDJ.dat$barcode
  if(extend.barcode){
    barcode <- paste0(barcode, '-1')
  }
  v_gene <- unlist(lapply(GEX.VDJ.dat$allVHitsWithScore, extract_gene_name_from_mixcr))
  d_gene <- unlist(lapply(GEX.VDJ.dat$allDHitsWithScore, extract_gene_name_from_mixcr))
  j_gene <- unlist(lapply(GEX.VDJ.dat$allJHitsWithScore, extract_gene_name_from_mixcr))
  c_gene <- unlist(lapply(GEX.VDJ.dat$allCHitsWithScore, extract_gene_name_from_mixcr))
  
  cdr3 <- GEX.VDJ.dat$aaSeqCDR3
  cdr3_nt <- GEX.VDJ.dat$nSeqCDR3
  umis <- GEX.VDJ.dat$cloneCount
  GEX.VDJ.new <- cbind(data.frame(barcode = barcode), true.value, null.value, true.value, 0, 
                       chain, v_gene, d_gene, j_gene, c_gene, true.value, true.value, cdr3, 
                       cdr3_nt, 0, umis, null.value, null.value)
  colnames(GEX.VDJ.new) <- c('barcode', 'is_cell', 'contig_id', 'high_confidence',
                             'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 
                             'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id', 'raw_consensus_id')
  GEX.VDJ.new <- GEX.VDJ.new %>% mutate_at(c('v_gene','d_gene','j_gene','c_gene'), 
                                           ~replace(., .=='', null.value))
  return(GEX.VDJ.new)
}


# processes the results from Mixcr before joining them to 10X dataset
# - GEX.VDJ.dat.A - list of preprocessed (see above function) results from Mixcr for TRA
# - GEX.VDJ.dat.B - list of preprocessed (see above function) results from Mixcr for TRB
# - VDJ.dat - 10X VDJ to filter results
# - returns combined data
# - mixcr data are not used for doublet filtering, as they seem to e a bit unreliable for that
process.VDJ.GEX <- function(GEX.VDJ.dat.A, 
                            GEX.VDJ.dat.B,
                            VDJ.dat = NULL){
  
  # first keep 2 TCR A and 1 TCR B for each cell with most counts
  # for this we perform ordering and take the last (two) sequences
  
  GEX.VDJ.dat.A <- GEX.VDJ.dat.A %>% arrange(barcode, umis) %>%
    group_by(barcode) %>% dplyr::filter(row_number() %in% c(n(), n()-1)) %>%
    dplyr::arrange(barcode, umis)
  
  GEX.VDJ.dat.B <- GEX.VDJ.dat.B %>% arrange(barcode, umis) %>%
    group_by(barcode) %>% dplyr::filter(row_number() == n())
  
  # if there are 10X VDJ then do the following
  if(!is.null(VDJ.dat)){
    # pre-process data so they contain only productive values (non-productive values 
    # -> are used for filtering only, which is not case here)
    # -> before that, we remove also non/existent values
    int_vdj <- VDJ.dat %>% dplyr::filter(cdr3 != 'None' & (cdr3 != '')) %>%
      group_by(barcode, cdr3) %>% dplyr::filter(row_number() == 1)
    
    VDJ.dat.A <- int_vdj %>% dplyr::filter(chain == 'TRA')
    VDJ.dat.B <- int_vdj %>% dplyr::filter(chain == 'TRB')
    
    count.VDJ.A <- table(VDJ.dat.A$barcode) 
    count.VDJ.B <- table(VDJ.dat.B$barcode)
    
    VDJ.dat.A <- VDJ.dat.A %>% 
      dplyr::filter(grepl('[T|t]rue', productive)) %>%
      dplyr::mutate(barcode.cdr3 = paste0(barcode, cdr3))
    
    VDJ.dat.B <- VDJ.dat.B %>% 
      dplyr::filter(grepl('[T|t]rue', productive)) %>%
      dplyr::mutate(barcode.cdr3 = paste0(barcode, cdr3))
    
    count.VDJ.A.prod <- table(VDJ.dat.A$barcode) 
    count.VDJ.B.prod <- table(VDJ.dat.B$barcode)    
    
    # remove any entries that were previously found in enriched VDJ analysis 
    # (have same barcode and CDR3 for the same chain), this is to prevent duplicated results
    # also, we want to filter out only if there is productive variant so it can fill the spot
    GEX.VDJ.dat.A <- GEX.VDJ.dat.A %>% dplyr::mutate(barcode.cdr3 = paste0(barcode, cdr3))%>%
      dplyr::filter(!(barcode.cdr3 %in% VDJ.dat.A$barcode.cdr3)) %>% dplyr::select(-barcode.cdr3)
    GEX.VDJ.dat.B <- GEX.VDJ.dat.B %>% dplyr::mutate(barcode.cdr3 = paste0(barcode, cdr3))%>%
      dplyr::filter(!(barcode.cdr3 %in% VDJ.dat.B$barcode.cdr3)) %>% dplyr::select(-barcode.cdr3)
    
    
  }
  # now bind the results
  GEX.VDJ.dat <- rbind(GEX.VDJ.dat.A, GEX.VDJ.dat.B)
  GEX.VDJ.dat <- as.data.frame(GEX.VDJ.dat)
  GEX.VDJ.dat <- add.vdj.source(GEX.VDJ.dat, 'mixcr')
  
  if(!is.null(VDJ.dat)){
    GEX.VDJ.dat <- rbind(VDJ.dat, GEX.VDJ.dat)
  }
  return(GEX.VDJ.dat)
}


# builds up a vdj data table that contains 10X and Mixcr VDJs
# - path.GEX.VDJ.A, path.GEX.VDJ.B - path to Mixcr VDJ files for TCRA and TCRB respectively
# - path.10X.data - path to 10X VDJ data
build.VDJ.table <- function(path.GEX.VDJ.A, path.GEX.VDJ.B, path.vdj.data = NULL){
  GEX.VDJ.A <- read.csv(path.GEX.VDJ.A, sep = '\t')
  GEX.VDJ.B <- read.csv(path.GEX.VDJ.B, sep = '\t')
  vdj.data <- NULL
  if(!is.null(path.vdj.data)){
    vdj.data <- load.10X.VDJ(path.vdj.data)
  }
  GEX.VDJ.dat.A <- build.GEX.VDJ.table(GEX.VDJ.A, 'TRA', extend.barcode = TRUE, vdj.dat = vdj.data)
  GEX.VDJ.dat.B <- build.GEX.VDJ.table(GEX.VDJ.B, 'TRB', extend.barcode = TRUE, vdj.dat = vdj.data)
  vdj.data <- process.VDJ.GEX(GEX.VDJ.dat.A, GEX.VDJ.dat.B, vdj.data)
  return(vdj.data)
}

# attempts to correct a (nt) mismatches in a data frame
# for this, v, (d), j, c genes must be identical
# mixcr only, as 10X does corrections itself
# - vdj.data - a data.frame where cdr3 will be corrected
# - max.mismatch - maximum difference for correction
# - min.cdr3.count - minimum number of VDJC + cdr3 to be considered
#   as reference for correction
# - max.corr.count
correct.mixcr.mismatches <- function(vdj.data,
                                     max.mismatch = 1,
                                     min.cdr3.count = 3,
                                     max.corr.count = 1){
  if(max.corr.count > min.cdr3.count){
    warning('min.cdr3.count lower than max.corr.count. 
           Setting max.corr.count to min.cdr3.count - 1.')
    max.corr.count <- min.cdr3.count - 1
  }
  
  # we need to establish the data
  cdr3.classes <- vdj.data %>% dplyr::filter(!(cdr3_nt %in% c('None',''))) %>%
    group_by(chain, cdr3, cdr3_nt) %>% 
    summarise(counts = n(), vdj_source = dplyr::first(vdj_source)) %>%
    ungroup %>% group_by(chain, cdr3) %>%
    mutate(count.cdr3 = sum(counts))
  
  # build table of references and table of corrections
  cdr3.reference <- cdr3.classes %>% dplyr::filter(counts >= min.cdr3.count)
  cdr3.corrected <- cdr3.classes %>% dplyr::filter(counts <= max.corr.count & vdj_source == 'mixcr')
  
  # pre-compute distance
  strdistmat <- stringdist::stringdistmatrix(cdr3.reference$cdr3_nt, cdr3.corrected$cdr3_nt)
  
  # now compare both tables and fix corrections if needed
  for(i in 1:ncol(strdistmat)){
    for(j in 1:nrow(strdistmat)){
      if(strdistmat[j,i] <= max.mismatch &
         cdr3.corrected$cdr3[i] != cdr3.reference$cdr3[j]){
        print(cdr3.reference[j,])
        print(cdr3.corrected[i,])
        print(vdj.data[
          vdj.data$cdr3 %in% cdr3.corrected$cdr3[i],] )
        vdj.data$cdr3[vdj.data$vdj_source == 'mixcr' &
                        vdj.data$cdr3 %in% cdr3.corrected$cdr3[i]] <- cdr3.reference$cdr3[j]
      }
    }
  }
  return(vdj.data)
}
