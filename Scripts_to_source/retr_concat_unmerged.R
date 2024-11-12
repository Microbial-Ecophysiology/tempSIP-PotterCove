# function to retrieve and concatenate unmerged sequences
retr_concat_unmerged <- function (mergers, long_seq, unmerged) {
mergers_final <- mergers
for(i in 1:length(mergers_final)) {
  if(i %in% unique(long_seq$sample_index)) {
    tmp_index <- long_seq$row_index[long_seq$sample_index == i]
    if(length(tmp_index) > 0) {
      mergers_final[[i]]$sequence[tmp_index] <- paste0(
        unmerged[[1]][paste0(long_seq$seqID[long_seq$sample_index == i], "/1")],
        "NNNNNNNNNN", 
        rc(unmerged[[2]][paste0(long_seq$seqID[long_seq$sample_index == i], "/2")])
      )
      mergers_final[[i]]$nmatch[tmp_index] <- 0
      mergers_final[[i]]$nmismatch[tmp_index] <- 0
      mergers_final[[i]]$nindel[tmp_index] <- 0
      mergers_final[[i]]$prefer[tmp_index] <- NA
      mergers_final[[i]]$accept[tmp_index] <- TRUE
      mergers_final[[i]] <- mergers_final[[i]][mergers_final[[i]]$accept, ]
    }
  } else {
    mergers_final[[i]] <- mergers_final[[i]][mergers_final[[i]]$accept, ]
  }
}

return(mergers_final)

}