## functions mysubset and mydbRDA_data for running dbRDAs with different subsets of large datasets

mysubset <- function(metadata, rel.abd.table, variable, levels_to_compare, id_column) 
{
  if (identical(sort(metadata[[id_column]]), sort(rownames(rel.abd.table)))) {
    sub_mdata <- subset(metadata, metadata[[variable]] %in% levels_to_compare)
    sub_rel.abd.table <- rel.abd.table[rownames(rel.abd.table) %in% sub_mdata[[id_column]],]
    output <- list("rel.abd.table" = sub_rel.abd.table, "mdata" = sub_mdata)
    return(output)
  }
  else {print("Error: rownames of rel.abd.table is not identical to id_column")}
}

mydbRDA_data <- function(metadata, rel.abd.table, var1, var2, id_column) {
  set.seed(96)
  dbRDA <- dbrda(rel.abd.table ~ metadata[[var1]] + metadata[[var2]], 
                 data = metadata, distance = "bray")
  sum_dbrda <- summary(dbRDA)
  dbRDA_stat <- anova.cca(dbRDA)
  message("\nstatistic whole model:")
  print(dbRDA_stat)
  dbRDA_stat_margin <- anova.cca(dbRDA, by = "margin")
  message("\nstatistic individual compounds:\nvar1 = ", var1,
          "\nvar2 = ", var2)
  print(dbRDA_stat_margin)
  message("\nvariation explained per dbRDA axis:")
  print(sum_dbrda$concont)
  dbrda_data <- as.data.frame(scores(dbRDA, display="sites"))
  dbrda_data[[id_column]] <- rownames(dbrda_data)
  dbrda_data_join <- dplyr::left_join(dbrda_data, metadata)
  data <- list("dbRDA_plot" = dbrda_data_join, "dbRDA" = dbRDA, "dbRDA_anova" = dbRDA_stat,
               "dbRDA_anova_margin" = dbRDA_stat_margin)
  return(data)
}
