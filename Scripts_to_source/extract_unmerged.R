# function to extract unmerged reads
extract_unmerged <- function(dadaF, dadaR, mergers) {
  
  # R1 read
  unmergedFs <- lapply(
    1:length(mergers),
    function(x) {
      tmp <- dadaF[[x]]$sequence[mergers[[x]]$forward[!mergers[[x]]$accept]]
      if (length(tmp) != 0) {
      names(tmp) <- paste0(x, "_", which(!mergers[[x]]$accept), "/1")
      return(tmp)
      }
    }
  )
  
  # R2 read
  unmergedRs <- lapply(
    1:length(mergers),
    function(x) {
      tmp <- dadaR[[x]]$sequence[mergers[[x]]$reverse[!mergers[[x]]$accept]]
      if (length(tmp) != 0) {
      names(tmp) <- paste0(x, "_", which(!mergers[[x]]$accept), "/2")
      return(tmp)
      }
    }
  )
  
  # create one fasta with unmerged from all samples
  fastaFs <- unlist(unmergedFs)
  fastaRs <- unlist(unmergedRs)
  
  # output is a list with paired R1 and R2
  return(list(fastaFs, fastaRs))
  
}
