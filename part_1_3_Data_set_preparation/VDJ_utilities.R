# Author: Juraj Michalik
# Date: 03/04/2021

# Script containing various functions to work with 10X VDJ
# -> mostly importing data to Seurat
# -> computing clones
require(tidyverse)
require(data.table)
require(pbapply)

##############################
# 10X VDJ                    #  
##############################

# performs treatment of 10X VDJs

# adds column source to vdj loaded from 10X so we know from where VDJ comes
# this column does not exist in 10X output so it needs to be added prior to any treatment
# - vdj.data - data to add vdj_source to 10X data
# - vdj.source - notation of source to add
add.vdj.source <- function(vdj.data, vdj.source = '10X'){
  vdj_source <- rep(vdj.source, nrow(vdj.data))
  vdj.data <- cbind(vdj.data, vdj_source)
  return(vdj.data)
}

# filters any NA from x; if there is a single NA, then returns it
# - x - value to filter NAs from
na.simplify <- function(x){
  if(any(!is.na(x))){
    return(na.omit(x))
  }else{
    return(NA)
  }
}


# loads VDJ file and adds vdj_source column to it
# - vdj.csv - VDJ file to load
load.10X.VDJ <- function(vdj.csv){
  vdj.data <- read.csv(vdj.csv)
  vdj.data <- add.vdj.source(vdj.data)
  return(vdj.data)
}


# finds VDJ doublets according to pre-set parameters
# - seurat.obj - seurat object where we want to trace duplicates
# - input.vdj  - a data frame from 10X (see above) containing VDJ information
# - lim.A.all, lim.B.all, lim.A.prod, lim.B.prod are maximum inclusive limits for unique TCRA, TCRB,
#   productive TCRA and productive TCRB. Uniqueness is defined by barcode, chain and CDR3
# - mixcr data are not considered for VDJ filtering as we find them a bit unreliable, 
#   filtering is done only on VDJ -> we exclude them at the beginning
trace.VDJ.duplicates <- function(seurat.obj, 
                                 input.vdj,
                                 lim.A.all = 2,
                                 lim.B.all = 2,
                                 lim.A.prod = 2,
                                 lim.B.prod = 1){
  
  # preprocess 10X table, notable re-order it and remove cdr3 duplicates or mixcr values
  input.vdj <- input.vdj %>% dplyr::arrange(chain, cdr3, productive, reads, v_gene, j_gene, d_gene, c_gene) %>%
    dplyr::filter(cdr3 != 'None' & cdr3 != '' & vdj_source == '10X') %>% 
    dplyr::filter(!duplicated(.[c('chain', 'barcode', 'cdr3')], fromLast = TRUE))
  
  barcode.counts <- input.vdj %>% group_by(chain, barcode) %>%
    summarise(counts = n())
  barcode.counts.prod <- input.vdj %>% group_by(chain, barcode) %>%
    dplyr::filter(grepl('[T|t]rue', productive)) %>% summarise(counts = n())
  
  # get all cells that have TCR sequences over given limits
  multiple.A <- barcode.counts %>% dplyr::filter(chain == 'TRA', counts > lim.A.all) %>% 
    pull(barcode)
  multiple.B <- barcode.counts %>% dplyr::filter(chain == 'TRB', counts > lim.B.all) %>% 
    pull(barcode)
  multiple.A.prod <- barcode.counts.prod %>% dplyr::filter(chain == 'TRA', counts > lim.A.prod) %>% 
    pull(barcode)
  multiple.B.prod <- barcode.counts.prod %>% dplyr::filter(chain == 'TRB', counts > lim.B.prod) %>% 
    pull(barcode)
  
  multiple.chain <- unique(c(multiple.A, multiple.B, multiple.A.prod, multiple.B.prod))
  multiple.chain <- intersect(multiple.chain, colnames(seurat.obj))

  return(multiple.chain)
}


# imports VDJ data from 10X into metadata of seurat Objects
# - seurat.obj - seurat.obj meta data to import
# - input.vdj - a data frame from 10X (see above) containing VDJ information
# - filter.vdj - whether to pre-filter vdj meta.data from doublets or not
#   productive TCRA and productive TCRB. Uniqueness is defined by barcode, chain and CDR3.
#   Considered only if filter.vdj is TRUE
# - generate.clonotypes - whether to create unique identifier for each clonotype -
#   - this takes a long time for big data-sets
# CAUTION: only A alpha chains and one B chain are considered! If not pre-filtered for duplicates,
#  will take the productives with most reads (first)
create.VDJ.metadata <- function(seurat.obj, 
                                input.vdj, 
                                filter.vdj = FALSE,
                                generate.clonotypes = TRUE){
  input.vdj <- input.vdj %>% dplyr::arrange(chain, cdr3, productive, reads, v_gene, j_gene, d_gene, c_gene) %>%
    dplyr::filter(cdr3 != 'None' & cdr3 != '') %>% 
    dplyr::filter(!duplicated(.[c('chain', 'barcode', 'cdr3')], fromLast = TRUE)) %>%
    dplyr::filter(barcode %in% colnames(seurat.obj))
  
  if(filter.vdj){
    barcode.counts <- input.vdj %>% dplyr::filter(vdj_source == '10X') %>%
      group_by(chain, barcode) %>% summarise(counts = n())
    barcode.counts.prod <- input.vdj %>% group_by(chain, barcode) %>%
      dplyr::filter(grepl('[T|t]rue', productive)) %>% summarise(counts = n())
    
    # get all cells that have TCR sequences over given limits
    multiple.A <- barcode.counts %>% dplyr::filter(chain == 'TRA', counts > 2) %>% 
      pull(barcode)
    multiple.B <- barcode.counts %>% dplyr::filter(chain == 'TRB', counts > 2) %>% 
      pull(barcode)
    multiple.A.prod <- barcode.counts.prod %>% dplyr::filter(chain == 'TRA', counts > 2) %>% 
      pull(barcode)
    multiple.B.prod <- barcode.counts.prod %>% dplyr::filter(chain == 'TRB', counts > 1) %>% 
      pull(barcode)
    
    multiple.chain <- unique(c(multiple.A, multiple.B, multiple.A.prod, multiple.B.prod))
    multiple.chain <- intersect(multiple.chain, colnames(seurat.obj))
    input.vdj <- input.vdj %>% dplyr::filter(!(barcode %in% multiple.chain))
  }
  
  # remove all non-productive VDJs
  input.vdj <- input.vdj %>% dplyr::filter(grepl('[T|t]rue', productive)) %>% 
    dplyr::arrange(chain, vdj_source, cdr3, desc(reads), v_gene, j_gene, d_gene, c_gene)
  
  # reformat data
  # first we'll treat columns we can't/don't wan't to just aggregate
  # Multi results will be aggregated
  vdj.data <- input.vdj %>% dplyr::select(barcode, chain, v_gene, d_gene, j_gene, cdr3, cdr3_nt)
  vdj.summary <- vdj.data %>% group_by(barcode, chain) %>% summarise(count = n())
  vdj.summary.A <- vdj.summary %>% dplyr::filter(chain == 'TRA') %>% dplyr::select(-chain) %>%
    dplyr::rename(count_TRA = count)
  vdj.summary.B <- vdj.summary %>% dplyr::filter(chain == 'TRB') %>% dplyr::select(-chain) %>% 
    dplyr::rename(count_TRB = count)
  vdj.summary.Multi <- vdj.summary %>% dplyr::filter(chain == 'Multi') %>% dplyr::select(-chain) %>%
    dplyr::rename(count_Multi = count)
  vdj.data.A1 <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'TRA') %>%  dplyr::filter(row_number()==1) %>%
    dplyr::select(-d_gene, -chain) %>%
    dplyr::rename(cdr3_A1 = cdr3, cdr3_A1_nt = cdr3_nt, 
                  v_gene_A1 = v_gene, j_gene_A1 = j_gene)
  vdj.data.A2 <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'TRA') %>% dplyr::filter(row_number()==2) %>%
    dplyr::select(-d_gene, -chain) %>% 
    dplyr::rename(cdr3_A2 = cdr3, cdr3_A2_nt = cdr3_nt, 
                  v_gene_A2= v_gene, j_gene_A2 = j_gene)
  vdj.data.B <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'TRB') %>% dplyr::filter(row_number()==1) %>%
    dplyr::select(-chain) %>% 
    dplyr::rename(cdr3_B = cdr3, cdr3_B_nt = cdr3_nt, v_gene_B= v_gene,
                  d_gene_B = d_gene, j_gene_B = j_gene)
  vdj.data.Multi <- vdj.data %>% group_by(barcode) %>%
    dplyr::filter(chain == 'Multi') %>% 
    summarise(v_gene = paste(v_gene, collapse = ','), d_gene = paste(d_gene, collapse = ','), 
              j_gene = paste(j_gene, collapse = ','), cdr3 = paste(cdr3, collapse = ','),
              cdr3_nt = paste(cdr3_nt, collapse = ',')) %>%
    dplyr::rename(cdr3_Multi = cdr3, cdr3_Multi_nt = cdr3_nt, v_gene_Multi = v_gene,
                  d_gene_Multi = d_gene, j_gene_Multi = j_gene)
  
  
  
  # descriptors for each cell
  VDJ.cell.descr <- input.vdj %>% dplyr::select(barcode, is_cell) %>% 
    group_by(barcode) %>% dplyr::filter(row_number()==1)
  
  # descriptors for TCR that are not that instrumental
  # filter to keep only first two TRA and single TRB (arranged as previously)
  VDJ.TCR.descr <- input.vdj %>% 
    group_by(barcode, chain) %>%
    dplyr::filter((row_number() <= 2 & chain == 'TRA') | (row_number() <= 1 & chain == 'TRB') | 
           (row_number() <= 1 & chain == 'Multi')) %>% ungroup %>%
    dplyr::select(barcode, chain, contig_id, high_confidence, productive, length, reads, umis, vdj_source) %>%
    group_by(barcode) %>% 
    summarise(chain = paste(chain, collapse = ','), contig_id = paste(contig_id, collapse = ','), 
              high_confidence = paste(high_confidence, collapse = ','), productive = paste(productive, collapse = ','),
              length = paste(length, collapse = ','), reads = paste(reads, collapse = ','),
              umis = paste(umis, collapse = ','), vdj_source = paste(vdj_source, collapse = ','))
  
  VDJ.meta.data <- VDJ.cell.descr %>% left_join(VDJ.TCR.descr, by = "barcode") %>% 
    ungroup %>% dplyr::arrange(barcode) %>% left_join(vdj.summary.A, by = "barcode") %>%
    left_join(vdj.summary.B, by = "barcode") %>%
    left_join(vdj.summary.Multi, by = "barcode") %>%
    left_join(vdj.data.A1, by = "barcode") %>% left_join(vdj.data.A2, by = "barcode") %>%
    left_join(vdj.data.B, by = "barcode") %>% left_join(vdj.data.Multi, by = "barcode")
  
  VDJ.meta.data <- as.data.frame(VDJ.meta.data)
  rownames(VDJ.meta.data) <- VDJ.meta.data$barcode
  VDJ.meta.data <- dplyr::select(VDJ.meta.data, is_cell, chain, contig_id, high_confidence, productive, length, 
                                 reads, umis, vdj_source, count_TRA, count_TRB, count_Multi, cdr3_A1, cdr3_A1_nt, 
                                 cdr3_A2, cdr3_A2_nt, cdr3_B, cdr3_B_nt, cdr3_Multi, cdr3_Multi_nt, v_gene_A1, j_gene_A1,
                                 v_gene_A2, j_gene_A2, v_gene_B, d_gene_B, j_gene_B, v_gene_Multi, d_gene_Multi, 
                                 j_gene_Multi)
  

  seurat.obj <- AddMetaData(seurat.obj, VDJ.meta.data)
  seurat.obj@meta.data <- droplevels(seurat.obj@meta.data)
  
  if(generate.clonotypes){
    clono.numbering <- get.sequence.clonotype(input.vdj)
    seurat.obj <- get.clonotypes(seurat.obj, clono.numbering)
  }
  return(seurat.obj)
}


# gets the clonotype tags for all sequences with given type
# - input.vdj - list for which the tags are created
get.sequence.clonotype <- function(input.vdj){
  clono.numbering <- input.vdj %>% 
    dplyr::select('chain', 'cdr3') %>% group_by(chain, cdr3) %>% 
    dplyr::filter(chain == 'TRA' | chain == 'TRB') %>%
    dplyr::filter(row_number()==1) %>% ungroup %>% group_by(chain) %>%
    mutate(clonotype = ifelse(chain == 'TRA', 'A', 'B')) %>%
    mutate(clonotype = paste0('clonotype-', clonotype, row_number()))
  
  return(clono.numbering)
}


# adds clonotypes to seurat object
# - seurat.obj - seurat.object to add clonotypes to
# - clono.numbering - clonotypes with precomputed numbering
get.clonotypes <- function(seurat.obj, clono.numbering){
  meta.data.table <- seurat.obj@meta.data
  # get unique cdr3 sequences along with their tag
  clono.numbering.A <- clono.numbering %>% dplyr::filter(chain == 'TRA') %>%
    ungroup %>% dplyr::select(cdr3, clonotype) 
  clono.numbering.B <- clono.numbering %>% dplyr::filter(chain == 'TRB') %>%
    ungroup %>% dplyr::select(cdr3, clonotype)
  

  # associate cells and clonotype tags
  meta.data.A1 <- meta.data.table %>% dplyr::select(cdr3_A1) %>%
    left_join(clono.numbering.A, by = c('cdr3_A1' = 'cdr3'))
  meta.data.A2 <- meta.data.table %>% dplyr::select(cdr3_A2) %>%
    left_join(clono.numbering.A, by = c('cdr3_A2' = 'cdr3')) 
  meta.data.B <- meta.data.table %>% dplyr::select(cdr3_B) %>%
    left_join(clono.numbering.B, by = c('cdr3_B' = 'cdr3'))

  # combine all tags, create clonotypes, filter them and put them into a list
  meta.data.tags <- cbind(meta.data.A1$clonotype, meta.data.A2$clonotype,
                          meta.data.B$clonotype)
  meta.data.tags <- as.data.frame(meta.data.tags)
  meta.data.tags <- meta.data.tags  %>% `colnames<-`(c('cdrA1', 'cdrA2', 'cdrB')) %>%
    mutate(clonotype.1 = ifelse(!is.na(cdrA1) & !is.na(cdrB), 
                                paste(cdrA1, cdrB, sep = '_'), NA)) %>% 
    mutate(clonotype.2 = ifelse(!is.na(cdrA2) & !is.na(cdrB), 
                                paste(cdrA2, cdrB, sep = '_'), NA))
  
  # create list with 0, 1 or 2 clones per cell (depending on TCR present)
  clonotypes <- mapply(c, as.list(meta.data.tags$clonotype.1),
                       as.list(meta.data.tags$clonotype.2),
                       SIMPLIFY = F)
  clonotypes <- lapply(clonotypes, na.simplify) # removing NAs
  
  
  # add clonotypes to metadata
  cellnames <- rownames(seurat.obj@meta.data) # transformation to data.table removes names, so we need to save then re-add them
  seurat.obj@meta.data <- data.table(seurat.obj@meta.data)
  seurat.obj@meta.data <- cbind(seurat.obj@meta.data, clonotypes) # a hack to insert list to data frame
  seurat.obj@meta.data <- data.frame(seurat.obj@meta.data)
  rownames(seurat.obj@meta.data) <- cellnames 
  
  return(seurat.obj)
}


# creates clones between various cells:
# -> unifies clones that have the same sequences if there are no conflicts
# -> B and at least one A must be same
# -> if there is a clone with a single A and clone with same A and B + another A, they fall to the same clone
# -> in case of multiple possibilities the most frequent clone wins
# - seurat.obj - object for which clones are calculated
create.clones <- function(seurat.obj){
  
  # pick up all seurat CDR3, then reorder A1 and A2 alphabetically
  cdr3.table <- seurat.obj@meta.data %>% 
    dplyr::select(cdr3_A1, cdr3_A2, cdr3_B)
  for(i in 1:nrow(cdr3.table)){
    cdr3.table[i,1:2] <- cdr3.table[i,order(as.character(cdr3.table[i,1:2]))]
  }
  cdr3.table[is.na(cdr3.table)] <- ''
  cdr3.table$tmp.clone <- ifelse(cdr3.table[,1] == '' | cdr3.table[,3] == '',
    'Incomplete', paste0(cdr3.table[,1], '-', cdr3.table[,2], ':', cdr3.table[,3]))
  
  # we group elements that have CDR3A x 2 first, then we
  # order those by number of occurrences
  # optimization: filter everything with single occurrence of CDR3B, as that can't be a clone
  cdr3.table.opt <- cdr3.table %>% group_by(cdr3_B) %>%
    dplyr::filter(n() > 1) %>% ungroup
  double.A.clones <- cdr3.table.opt %>% 
    dplyr::filter(cdr3_A1 != '' & cdr3_A2 != '' & cdr3_B != '') %>%
    group_by(tmp.clone) %>% 
    summarise(cdr3_A1 = dplyr::first(cdr3_A1), cdr3_A2 = dplyr::first(cdr3_A2),  
              cdr3_B = dplyr::first(cdr3_B), occurences = n()) %>%
    ungroup %>% dplyr::arrange(desc(occurences))
  
  # now do the same for cells with single CDR3A
  single.A.clones <- cdr3.table.opt %>% 
    dplyr::filter(cdr3_A1 != '' & cdr3_A2 == '' & cdr3_B != '') %>%
    group_by(tmp.clone) %>% 
    summarise(cdr3_A1 = dplyr::first(cdr3_A1), cdr3_A2 = dplyr::first(cdr3_A2),  
              cdr3_B = dplyr::first(cdr3_B), occurences = n()) %>%
    ungroup %>% dplyr::arrange(desc(occurences))
  
  # we build a dictionary storing which clononotype is associated with which clone
  # first we add sequences that have 2 x CDR3A as those won't change
  double.A.clones.multi <- double.A.clones %>% dplyr::filter(occurences > 1)
  clonodict = as.list(double.A.clones.multi$tmp.clone)
  names(clonodict) = double.A.clones.multi$tmp.clone

  # now look through list of cells having single CDR3A and look which have 
  # CDR3A and CDR3B in common: in that case add them to clonodict
  if(nrow(single.A.clones) > 0){
    for(i in 1:nrow(single.A.clones)){
      singular.clones <- c() # singular cells get a special treatment
      if(nrow(double.A.clones) > 0){
        for(j in 1:nrow(double.A.clones)){
          if((single.A.clones$cdr3_A1[i] == double.A.clones$cdr3_A1[j] |
              single.A.clones$cdr3_A1[i] == double.A.clones$cdr3_A2[j]) &
             single.A.clones$cdr3_B[i] == double.A.clones$cdr3_B[j]){
            if(double.A.clones$occurences[j] > 1){ # if multiple CDR3A x 2 clones, treat immediately
              clonodict[single.A.clones$tmp.clone[i]] <- double.A.clones$tmp.clone[j]
              break
            }else{ # if singular CDR3A x 2 clone, add it to list for now
              singular.clones <- c(singular.clones, double.A.clones$tmp.clone[j])
            }
          }
        }
      }
      # Basically:
      # if there is only a single singular (appearing once) clone (clone with 2 X CDR3A) that matches with 
      # current clone with single CDR3A then it is most likely that one and we need just to look whether it was
      # already added into list of clones
      # if there are multiple singular clones, then we can't decide and add them just as single CDR3A clone, but with 
      # checking whether they were added before, and if all were added before, we add it to first one because (random)
      if(length(singular.clones) == 1){
        # if we have 2 x CDR3A clone with single occurrence and another without second CDR3A,
        # we have to check whether it was not added previously by something else
        # 1.) no: add it with current clonotype with single CDR3A
        # 2.) yes: check if it is compatible with associated single CDR3A clonotype and only if yes, then add CDR3A with given clone
        if(!(singular.clones[1] %in% names(clonodict))){
          clonodict[singular.clones[1]] <- singular.clones[1]
          clonodict[single.A.clones$tmp.clone[i]] <- singular.clones[1]
        }else if(single.A.clones$cdr3_A1[i] %in% clonodict[singular.clones[1]] &
                 single.A.clones$cdr3_B[i] %in% clonodict[singular.clones[1]]){
          clonodict[singular.clones[1]] <- singular.clones[1]
        }
      }else if(length(singular.clones) > 1){
      # think about this part again ---->
      # similar to above; it is enough if only one is not annotated or is compatible
        added <- 0
        # check how many clonotypes are not added
        for(k in 1:length(singular.clones)){
          if(!(clonodict[singular.clones[k]] %in% names(clonodict))){
            added <- added + 1
          }
        }
        # if at least one added then create group where single clonotype will be considered
        if(added > 0){ # add immediately if at least one clonotype was added previously
          clonodict[single.A.clones$tmp.clone[i]] <- single.A.clones$tmp.clone[i]
          for(k in 1:length(singular.clones)){
            if(!(clonodict[singular.clones[k]] %in% names(clonodict))){
              clonodict[singular.clones[k]] <- single.A.clones$tmp.clone[i]
            }
          }
        # if there is not a single CDR3 x2  clone that was not added yet, check whether there is
        # some clone that is compatible with given CDR3
        }else if(added == 0){
          for(k in 1:length(singular.clones)){
            if(single.A.clones$cdr3_A1[i] %in% clonodict[singular.clones[k]] &
               single.A.clones$cdr3_B[i] %in% clonodict[singular.clones[k]]){
              clonodict[single.A.clones$tmp.clone[i]] <- singular.clones[k]
              break
            }
          }
        }
      }
    }
  }
  
  # finally, add the remaining multiple times occurring single 
  # CDR3A clonotypes not added yet
  clonodict.names <- names(clonodict)
  single.A.clones <- single.A.clones %>% 
    dplyr::filter(occurences > 1 & !(tmp.clone %in% clonodict.names))
  clonodict <- c(clonodict, as.list(single.A.clones$tmp.clone))
  names(clonodict) <- c(clonodict.names, single.A.clones$tmp.clone)
  
  # add to data
  seurat.obj@meta.data$clonotype.repeated <- 
    ifelse(cdr3.table$tmp.clone %in% names(clonodict),
           as.character(clonodict[cdr3.table$tmp.clone]),
           'Other')
  
  return(seurat.obj)
}

##############################
# VDJ MISC                   #  
##############################

# function that artificially copies VDJ to another cells if that sequence is present in another clone and there is second sequence existing for it
# - vdj - list of VDJ where missing sequences will be added
# - CDR3_1 - missing CDR3A
# - CDR3_2 - present CDR3B
# - missing.chain - missing CDR3 is A or B, default 'TRA'

# !!!WARNING!!!: this should NOT be used unless in really specific circumstances!
# - Use at your own risk.
insert.missing <- function(vdj,
                           CDR3_1,
                           CDR3_2,
                           missing.chain = 'TRA'){
  present.chain = 'TRB'
  if(missing.chain == 'TRB'){
    present.chain = 'TRA' # other chain is reference
  }
  # consider productives only
  vdj.missing <- vdj$barcode[vdj$chain == missing.chain & 
                             vdj$cdr3 == CDR3_1 &
                             grepl('[T|t]rue', vdj$productive)]
  vdj.present <- vdj$barcode[vdj$chain == present.chain & 
                             vdj$cdr3 == CDR3_2 &
                             grepl('[T|t]rue', vdj$productive)]
  
  vdj.missing <- vdj.present[!(vdj.present %in% vdj.missing)]
  
  vdj.fill <- vdj %>% 
    dplyr::filter(chain == missing.chain & cdr3 == CDR3_1) %>%
    dplyr::filter(row_number() == 1) %>%
    mutate(count = length(vdj.missing)) %>%
    uncount(count) %>% mutate(barcode = vdj.missing, reads = 0, umis = 1) %>%
    mutate(contig_id = paste0(barcode, '_contig_1'))
  
  vdj <- rbind(vdj, vdj.fill)
  return(vdj)
}

# searches and plots cells having specific CDR3 if found
# - seurat.obj - seurat object to find VDJ in
# - cdr3 - either single sequence or vector of all sequences to find in single cell
# - chain - vector identifying all sequences in CDR3 in same order. Values either CDR3A or CDR3B. 
# - - maximum 2 TRA and 1 TRB. 
find.cdr3s <- function(seurat.obj,
                       cdr3,
                       chain){
  cell.list <- colnames(seurat.obj)
  
  if(is.null(cdr3) || is.null(chain)){
    stop('CDR3 or chain are not provided.')
  }
  
  if('TRB' %in% chain){
    CDR3B.cells <- colnames(seurat.obj)[seurat.obj$cdr3_B %in% cdr3[chain == 'TRB']]
    cell.list <- intersect(CDR3B.cells, cell.list)
  }
  if('TRA' %in% chain){
    if(length(chain[chain == 'TRA'] == 1)){
      CDR3A.cells <- colnames(seurat.obj)[seurat.obj$cdr3_A1 %in% cdr3[chain == 'TRA']|
                                          seurat.obj$cdr3_A2 %in% cdr3[chain == 'TRA']]
      cell.list <- intersect(CDR3A.cells, cell.list)
    }else if(length(chain[chain == 'TRA'] == 2)){
      CDR3A.cells <- colnames(seurat.obj)[(seurat.obj$cdr3_A1 %in% cdr3[chain == 'TRA'][1] &
                                          seurat.obj$cdr3_A2 %in% cdr3[chain == 'TRA'][2]) |
                                          (seurat.obj$cdr3_A1 %in% cdr3[chain == 'TRA'][2] &
                                          seurat.obj$cdr3_A2 %in% cdr3[chain == 'TRA'][1])]
      cell.list <- intersect(CDR3A.cells, cell.list)
    }
  }

  if (length(cell.list) == 0){
    message("No cell with given sequence found. Plot won't be traced.")
  }else{
    cdr3.seq <- paste(c(sort(cdr3[chain == 'TRA']), cdr3[chain == 'TRB']), collapse = '-')
    seurat.obj$sel.clone <- ifelse(colnames(seurat.obj) %in% cell.list, cdr3.seq, 'Other')
    print(DimPlot(seurat.obj, group.by = 'sel.clone', 
                  cols = c('firebrick', adjustcolor('gray82', alpha.f = c(0.3))), raster= F))
  }  
}