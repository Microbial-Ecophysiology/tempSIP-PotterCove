# function for modifying asv tables with only taxa above certain threshold
## now also checks if summed up "others" for that rank level is above threshold - work in progress

## corrected script, previous script did not sum all ASVs per taxon, but checked if any ASV is above threshold

## also includes function to re-order and re-name vector for sorting taxa easily, see at end of script

## info input
# table: input table, long format
# abundance threshold: taxa above which threshold in at least one sample should be kept separate, in %
# Abundance: name of column which has abundance numbers, not in %
 # by default named 'Abundance'
# Sample: name of column which has unique name for every sequenced sample
 # by default named 'Sample'


# function to change taxon names if below set threshold
sort_abundant_taxa <- function(table, abundance_threshold, Abundance = "Abundance", Sample = "Sample") {

  # 1. check if rank names are correct ####
  if (identical(tail(colnames(phylo_melt), n=6), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))) {
    
  # 2. function to check if taxon is above threshold ####
    list_taxa_above_threshold <- function(input_table, rank, threshold, abundance_column, sample_column) {
      df1 <- input_table %>% 
        group_by({{rank}}, {{sample_column}}) %>%                                          ## group table by rank (Genus, Family, etc.) and by sample
        summarise(sum_tax_abundance = sum({{abundance_column}}), .groups = "drop") %>%     ## sum abundance of selected rank together, but separately per sample
        ## check across all samples if taxon is above threshold in at least one
        group_by({{rank}}) %>% 
        summarise(max_abundance = max(sum_tax_abundance), .groups = "keep") %>% 
        filter(max_abundance > threshold/100) %>% 
        arrange(desc(max_abundance))
      out_vector <- pull(df1, {{rank}})
      return(out_vector)
    }
    
    
  # 3. rename NAs in taxonomy to unclassified specifying highest assigned rank, no ASV ID ####
    pf0 <- table %>% 
      mutate(Phylum_uncl = if_else(is.na(Phylum) == FALSE, Phylum, paste0(Kingdom, "_p_unclassified")),
             
             Class_uncl = if_else(is.na(Class) == FALSE, Class, 
                                  if_else(is.na(Phylum) == FALSE, paste0(Phylum, "_c_unclassified"),
                                          paste0(Kingdom, "_p_unclassified"))),
             
             Order_uncl = if_else(is.na(Order) == FALSE, Order,
                                  if_else(is.na(Class) == FALSE, paste0(Class, "_o_unclassified"),
                                          if_else(is.na(Phylum) == FALSE, paste0(Phylum, "_c_unclassified"),
                                                  paste0(Kingdom, "_p_unclassified")))),
             
             Family_uncl = if_else(is.na(Family) == FALSE, Family,
                                   if_else(is.na(Order) == FALSE, paste0(Order, "_f_unclassified"),
                                           if_else(is.na(Class) == FALSE, paste0(Class, "_o_unclassified"),
                                                   if_else(is.na(Phylum) == FALSE, paste0(Phylum, "_c_unclassified"),
                                                           paste0(Kingdom, "_p_unclassified"))))),
             
             Genus_uncl = if_else(is.na(Genus) == FALSE, Genus,
                                  if_else(is.na(Family) == FALSE, paste0(Family, "_g_unclassified"),
                                          if_else(is.na(Order) == FALSE, paste0(Order, "_f_unclassified"),
                                                  if_else(is.na(Class) == FALSE, paste0(Class, "_o_unclassified"),
                                                          if_else(is.na(Phylum) == FALSE, paste0(Phylum, "_c_unclassified"),
                                                                  paste0(Kingdom, "_p_unclassified")))))),
             
             ASV_uncl = if_else(is.na(Genus) == FALSE, paste0(Genus, "_", OTU),
                                if_else(is.na(Family) == FALSE, paste0(Family, "_g_unclassified_", OTU),
                                        if_else(is.na(Order) == FALSE, paste0(Order, "_f_unclassified_", OTU),
                                                if_else(is.na(Class) == FALSE, paste0(Class, "_o_unclassified_", OTU),
                                                        if_else(is.na(Phylum) == FALSE, paste0(Phylum, "_c_unclassified_", OTU),
                                                                paste0(Kingdom, "_p_unclassified_", OTU)))))))
    
    
  # 4. check with taxa are above threshold for each rank ####
    ## 4.0 ASV level ####
    # only find ASVs, don't change table
    ### just ASV ID
    asv_thr <- list_taxa_above_threshold(input_table = pf0,
                                         rank = OTU,
                                         threshold = abundance_threshold,
                                         abundance_column = Abundance,
                                         sample_column = Sample)
    
    ### named ASV by highest classified rank
    asv_thr_named <- list_taxa_above_threshold(input_table = pf0,
                                               rank = ASV_uncl,
                                               threshold = abundance_threshold,
                                               abundance_column = Abundance,
                                               sample_column = Sample)
    
    ## 4.1 genus level ####
    # check which genera are above threshold and save in vector
    genus_thr_g <- list_taxa_above_threshold(input_table = pf0,
                                             rank = Genus_uncl,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use on table to create new column, rename genus to "other_" family if not above threshold
    pf_all <- pf0 %>% 
      mutate(Genus_mod_f = if_else(Genus_uncl %in% genus_thr_g, 
                                   # if genus is in list, keep name as genus
                                   Genus_uncl, 
                                   # if genus not in list, rename as "other_" family
                                   paste0("other_f_", Family_uncl, "_<", abundance_threshold, "%")))
    
    ## family level
    # check again which taxa groups are above threshold
    genus_thr_f <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Genus_mod_f,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by order
    pf_all <- pf_all %>% 
      mutate(Genus_mod_o = if_else(Genus_mod_f %in% genus_thr_f, 
                                   # if is in list, keep name
                                   Genus_mod_f, 
                                   # if not in list, rename as "other_" order
                                   paste0("other_o_", Order_uncl, "_<", abundance_threshold, "%")))
    
    ## order level
    # check again which taxa groups are above threshold
    genus_thr_o <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Genus_mod_o,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by class
    pf_all <- pf_all %>% 
      mutate(Genus_mod_c = if_else(Genus_mod_o %in% genus_thr_o, 
                                   # if is in list, keep name 
                                   Genus_mod_o, 
                                   # if not in list, rename as "other_" class
                                   paste0("other_c_", Class_uncl, "_<", abundance_threshold, "%")))
    
    ## class level
    # check again which taxa groups are above threshold
    genus_thr_c <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Genus_mod_c,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Genus_mod_p = if_else(Genus_mod_c %in% genus_thr_c, 
                                   # if is in list, keep name
                                   Genus_mod_c, 
                                   # if not in list, rename as "other_" phylum
                                   paste0("other_p_", Phylum_uncl, "_<", abundance_threshold, "%")))
    
    ## phylum level
    # check again which taxa groups are above threshold
    genus_thr_p <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Genus_mod_p,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Genus_abt = if_else(Genus_mod_p %in% genus_thr_p, 
                                 # if is in list, keep name
                                 Genus_mod_p, 
                                 # if not in list, rename as "other_" kingdom
                                 paste0("other_", Kingdom, "_<", abundance_threshold, "%")))
    
    
    ## 4.2 family level ####
    # check which families are above threshold and save in vector
    family_thr_f <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Family_uncl,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use on table to create new column, renaming by order
    pf_all <- pf_all %>% 
      mutate(Family_mod_o = if_else(Family_uncl %in% family_thr_f, 
                                   # if is in list, keep name
                                   Family_uncl, 
                                   # if not in list, rename as "other_" order
                                   paste0("other_o_", Order_uncl, "_<", abundance_threshold, "%")))
    
    ## order level
    # check again which taxa groups are above threshold
    family_thr_o <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Family_mod_o,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by class
    pf_all <- pf_all %>% 
      mutate(Family_mod_c = if_else(Family_mod_o %in% family_thr_o, 
                                   # if is in list, keep name 
                                   Family_mod_o, 
                                   # if not in list, rename as "other_" class
                                   paste0("other_c_", Class_uncl, "_<", abundance_threshold, "%")))
    
    ## class level
    # check again which taxa groups are above threshold
    family_thr_c <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Family_mod_c,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Family_mod_p = if_else(Family_mod_c %in% family_thr_c, 
                                   # if is in list, keep name
                                   Family_mod_c, 
                                   # if not in list, rename as "other_" phylum
                                   paste0("other_p_", Phylum_uncl, "_<", abundance_threshold, "%")))
    
    ## phylum level
    # check again which taxa groups are above threshold
    family_thr_p <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Family_mod_p,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Family_abt = if_else(Family_mod_p %in% family_thr_p, 
                                 # if is in list, keep name
                                 Family_mod_p, 
                                 # if not in list, rename as "other_" kingdom
                                 paste0("other_", Kingdom, "_<", abundance_threshold, "%")))
    
    
    ## 4.3 order level ####
    # check which orders are above threshold and save in vector
    order_thr_o <- list_taxa_above_threshold(input_table = pf_all,
                                              rank = Order_uncl,
                                              threshold = abundance_threshold,
                                              abundance_column = Abundance,
                                              sample_column = Sample)
    
    # use on table to create new column, renaming by class
    pf_all <- pf_all %>% 
      mutate(Order_mod_c = if_else(Order_uncl %in% order_thr_o, 
                                    # if is in list, keep name 
                                    Order_uncl, 
                                    # if not in list, rename as "other_" class
                                    paste0("other_c_", Class_uncl, "_<", abundance_threshold, "%")))
    
    ## class level
    # check again which taxa groups are above threshold
    order_thr_c <- list_taxa_above_threshold(input_table = pf_all,
                                              rank = Order_mod_c,
                                              threshold = abundance_threshold,
                                              abundance_column = Abundance,
                                              sample_column = Sample)
    
    # use again on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Order_mod_p = if_else(Order_mod_c %in% order_thr_c, 
                                    # if is in list, keep name
                                    Order_mod_c, 
                                    # if not in list, rename as "other_" phylum
                                    paste0("other_p_", Phylum_uncl, "_<", abundance_threshold, "%")))
    
    ## phylum level
    # check again which taxa groups are above threshold
    order_thr_p <- list_taxa_above_threshold(input_table = pf_all,
                                              rank = Order_mod_p,
                                              threshold = abundance_threshold,
                                              abundance_column = Abundance,
                                              sample_column = Sample)
    
    # use again on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Order_abt = if_else(Order_mod_p %in% order_thr_p, 
                                  # if is in list, keep name
                                  Order_mod_p, 
                                  # if not in list, rename as "other_" kingdom
                                  paste0("other_", Kingdom, "_<", abundance_threshold, "%")))
    
    
    ## 4.4 class level ####
    # check which classes are above threshold and save in vector
    class_thr_c <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Class_uncl,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use on table to create new column, renaming by phylum
    pf_all <- pf_all %>% 
      mutate(Class_mod_p = if_else(Class_uncl %in% class_thr_c, 
                                   # if is in list, keep name
                                   Class_uncl, 
                                   # if not in list, rename as "other_" phylum
                                   paste0("other_p_", Phylum_uncl, "_<", abundance_threshold, "%")))
    
    ## phylum level
    # check again which taxa groups are above threshold
    class_thr_p <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Class_mod_p,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use again on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Class_abt = if_else(Class_mod_p %in% class_thr_p, 
                                 # if is in list, keep name
                                 Class_mod_p, 
                                 # if not in list, rename as "other_" kingdom
                                 paste0("other_", Kingdom, "_<", abundance_threshold, "%")))
    
    
    ## 4.5 phylum level ####
    # check which phyla are above threshold and save in vector
    phylum_thr_p <- list_taxa_above_threshold(input_table = pf_all,
                                             rank = Phylum_uncl,
                                             threshold = abundance_threshold,
                                             abundance_column = Abundance,
                                             sample_column = Sample)
    
    # use on table to create new column, renaming now by phylum
    pf_all <- pf_all %>% 
      mutate(Phylum_abt = if_else(Phylum_uncl %in% phylum_thr_p, 
                                 # if is in list, keep name
                                 Phylum_uncl, 
                                 # if not in list, rename as "other_" kingdom
                                 paste0("other_", Kingdom, "_<", abundance_threshold, "%")))
    
    pf_all <- pf_all %>% 
      arrange(Kingdom, Phylum_abt, Class_abt, Order_abt, Family_abt, Genus_abt) %>% 
      
      ## remove not needed columns
      select(-Genus_mod_f, -Genus_mod_o, -Genus_mod_c, -Genus_mod_p,
             -Family_mod_o, -Family_mod_c, -Family_mod_p,
             -Order_mod_c, -Order_mod_p,
             -Class_mod_p)
    

  # 5. Output #### 
  ## create list which contains individual vectors for each taxon rank
    output <- vector("list", 9)
    output[[1]] <- abundance_threshold
    output[[2]] <- asv_thr
    output[[3]] <- asv_thr_named
    output[[4]] <- genus_thr_g
    output[[5]] <- family_thr_f
    output[[6]] <- order_thr_o
    output[[7]] <- class_thr_c
    output[[8]] <- phylum_thr_p
    # actual long output table
    output[[9]] <- pf_all
    
    # give the elements of the list names which contains information of the threshold used
    names(output) <- c("used_abundance_threshold",
                       "ASVs_abv_thr",
                       "ASVs_named_abv_thr",
                       "genera_abv_thr",
                       "families_abv_thr",
                       "orders_abv_thr",
                       "classes_abv_thr",
                       "phyla_abv_thr",
                       "ASV_table_taxa_abv_thr")
    
    # output list
    return(output)
  }
  else {print("Error: rank names have to be changed to 'Kingdom', 'Phylum', 'Class', 'Order', 'Family',  'Genus'")}
  
}


# reordering phyla and class entries using function from https://github.com/mrdwab/SOfun ####
#' Reorders the Contents of a Vector
#' 
#' Shuffle the order of a vector around using natural language statements.
#' 
#' This can be a useful function for reordering the columns of a
#' \code{data.frame} or \code{data.table} in a convenient manner. In such
#' cases, the \code{invec} would be \code{names(your_data_frame)}. When using
#' \code{data.table}s, remember to use \code{setcolorder} to avoid copying.
#' 
#' The \code{movecommand} argument is specified in the form of \code{"a, b
#' before f"}. The positions to move are: \itemize{ \item \strong{first}: move
#' the specified items to the first postion. \item \strong{last}: move the
#' specified items to the last position. \item \strong{before}: move the
#' specified items before the value mentioned. \item \strong{after}: move the
#' specified items after the value mentioned. } Multiples are allowed:
#' \itemize{ \item Specify multiple values to be moved by separating them with
#' a comma. \item Chain multiple move commands by separating them with a
#' semicolon. }
#' 
#' @param invec The input vector
#' @param movecommand The command that describes how you want to shuffle the
#' vector. See \emph{Details}.
#' @return A vector.
#' @author Ananda Mahto
#' @references \url{http://stackoverflow.com/a/18420673/1270695}
#' @examples
#' 
#' myvec <- letters[1:10]
#' myvec
#' moveMe(myvec, "a last; b, e, g before d; c first; h after j")
#' 
#' x <- names(mtcars)
#' x
#' moveMe(x, "hp first; cyl after drat; vs, am, gear before mpg; wt last")
#' 
#' @export moveMe
moveMe <- function(invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], ",|\\s+"), 
                        function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp)-1
      } else if (A == "after") {
        after <- match(ba, temp)
      }    
    } else if (A == "first") {
      after <- 0
    } else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}
NULL

