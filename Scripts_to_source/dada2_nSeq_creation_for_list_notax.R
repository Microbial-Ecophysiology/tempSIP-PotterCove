## this script is meant to be sourced during the dada2 pipeline to provide functions to create nSeq files

# function to count unique sequences in dadaFR files
 # used in functions below
getN <- function(x) sum(getUniques(x))


# function to create nSeq for each lib using merging by rownames to make sure nothing is mixed up
 # output includes columns for read numbers after demultiplexing, primer clipping, filtering, denoising and merging
 # and all separate numbers for fr and rf reads used to sum them up

create_nSeq <- function(list_input, lib_no_input, nSeq_file_input) { # lib_no_input needs to have '_' in front, makes it possible to give also NULL if only 1 library is analysed
  nSeq <-
    list(
  # 1. Demux read no.
  (
    read.table(nSeq_file_input, h = T, stringsAsFactors = F, sep = "\t") %>% 
      tidyr::separate(SID, into = c("SID_n"), sep = "_")
  ),
  # 2. clipped and filtered read no.
  (
    list_input[["filt_FR.out"]] %>%                            # file with numbers for clipped and filtered reads of fr reads
      base::as.data.frame() %>%                                             # reformat from some kind of list to data.frame
      tibble::rownames_to_column() %>%                                        # new column with former rownames
      tidyr::separate(rowname, into = c("SID_n"), sep = "_") %>%             # create from columns new column with sample IDs without ending
      dplyr::rename("clipped.fr" = "reads.in", "filtered.fr" = "reads.out")  # make new unique colnames to be able to distinguish with rf file
  ),
  (
    list_input[["filt_RF.out"]] %>% 
      base::as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      tidyr::separate(rowname, into = c("SID_n"), sep = "_") %>% 
      dplyr::rename("clipped.rf" = "reads.in", "filtered.rf" = "reads.out")
  ),
  # 3. denoised read no.
  (
    base::sapply(list_input[["dadaFR_R1"]], getN) %>% 
      base::as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_fwd_fr" = ".")
  ),
  (
    base::sapply(list_input[["dadaRF_R1"]], getN) %>% 
      base::as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_fwd_rf" = ".")
  ),
  (
    base::sapply(list_input[["dadaFR_R2"]], getN) %>% 
      base::as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_rev_fr" = ".")
  ),
  (
    base::sapply(list_input[["dadaRF_R2"]], getN) %>% 
      base::as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_rev_rf" = ".")
  ),
  # 4. merged
  (
    base::sapply(list_input[["mergers_FR"]], getN) %>% 
      base::as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "merged_fr" = ".")
  ),
  (
    base::sapply(list_input[["mergers_RF"]], getN) %>% 
      base::as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "merged_rf" = ".")
  )) %>% 
  # 5. join all tables
  purrr::reduce(full_join, by = "SID_n") %>% 
  # 6. sum fr and rf reads 
  dplyr::mutate(Clipped = clipped.fr + clipped.rf,                   # sum fr and rf read numbers for clipped
         Filtered = filtered.fr + filtered.rf,                # ...for filtered            
         Denoised_fwd = denoised_fwd_fr + denoised_fwd_rf,    # ...for denoised fwd (=R1)
         Denoised_rev = denoised_rev_fr + denoised_rev_rf,    # ...for denoised rev (=R2)
         Merged = merged_fr + merged_rf)                      # ...for merged
  return(nSeq)
}


# this function merges the nSeq files created above
 # 1. uses if_else function to check if samples were sequenced on multiple lanes or/and libraries.
 #    if that was the case, no. reads counted for each step have to be summed across lanes and libraries
 #
 # 2. adds read numbers after Chimera removal, removal of too long or short sequences 
 # and number of reads after classification and with different bootstrap values
 # it uses a named list with all the different nSeq objects produced before as input

amend_nSeq <- function(lib_list_input) {
  # combine nSeq files produced individually for each lane/library into one data frame
  nSeq_comb <- bind_rows(lib_list_input, .id = "lib")    
  
  # 1. check if samples were sequenced in multiple libraries or lanes:
  # check if there are duplicates in the sample ID column, this gives a TRUE-FALSE vector. if there are any TRUE in the vector, the sum will by >0
  
  if (sum(duplicated(nSeq_comb[["SID_n"]])) == 0) {
    
    # a) if no samples were done in multiple lanes continue with script not summing anything 
    nSeq2 <- nSeq_comb %>% 
      list(
        (
          # read no. after chimera removal
          base::rowSums(seqtab.nochim) %>% 
            base::as.data.frame() %>%
            tibble::rownames_to_column() %>%
            dplyr::rename("SID_n" = "rowname", "Nochim" = ".")
        ),
        (
          # read no. after discarding too long or too short sequences
          base::rowSums(seqtab.nochim2) %>% 
            base::as.data.frame() %>% 
            tibble::rownames_to_column() %>%
            dplyr::rename("SID_n" = "rowname", "Opt.length" = ".")
        )) %>% 
      # join tables
      purrr::reduce(full_join, by = "SID_n")
    return(nSeq2)
    
  } else {
    
    # b) if any sample was sequenced on multiple lanes, sum numbers for each step counted before for each sample across the different lanes
    nSeq2 <- nSeq_comb %>%
      
      # group per sample
      group_by(SID_n) %>%
      
      dplyr::summarise(
        # change lib column so seq lane IDs are pasted together if sample was sequenced in multiple lanes
        lib = str_c(lib, separate = ',', collapse = ""), .groups = "keep",
        
        # sum all numeric columns per sample
        across(where(is.numeric), ~ sum(.x, na.rm = T))
        ) %>% 
      
      # add read numbers for different steps
      list(
        (
          # read no. after chimera removal
          base::rowSums(seqtab.nochim) %>% 
            base::as.data.frame() %>%
            tibble::rownames_to_column() %>%
            dplyr::rename("SID_n" = "rowname", "Nochim" = ".")
        ),
        (
          # read no. after discarding too long or too short sequences
          base::rowSums(seqtab.nochim2) %>% 
            base::as.data.frame() %>% 
            tibble::rownames_to_column() %>%
            dplyr::rename("SID_n" = "rowname", "Opt.length" = ".")
        )) %>% 
      # join tables
      purrr::reduce(full_join, by = "SID_n")
    
    # also return nSeq for individual lanes
    nSeq_indiv_seqlanes <<- nSeq_comb
    
    # print message
    writeLines("=======================================\nNOTICE!\nSome or all samples were sequenced on multiple sequencing lanes.\nThe read numbers for these samples were summed up from the different lanes.\nAn additional data frame called 'nSeq_indiv_seqlanes' is returned, which gives the individual read numbers for each lane until the merging step.\n=======================================")
    
    # return actual nSeq table
    return(nSeq2)
  }
  
}


# this function calculates percent of reads lost in each step
 # it uses the nSeq table from above as input
perc_nSeq <- function(nSeq_input) {
  nSeq3 <- nSeq_input %>% 
    dplyr::select(lib, SID_n, Demux, Clipped, Filtered, Denoised_fwd, Denoised_rev,   # select only columns with combined reads
                  Merged, Nochim, Opt.length) %>% 
    dplyr::mutate(across(where(is.numeric), ~ (. / Demux) * 100)) %>%                 # calculate percent of reads retained in each step compared to reads after demultiplexing
    dplyr::mutate(across(where(is.numeric), ~ round(.x, 2)))                          # round to 2 digits after decimal
  return(nSeq3) 
}

