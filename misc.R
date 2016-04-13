# DATAFRAME DESCRIBING TRANSITIONS ---------------------------------------------
to.trans <- function(trans) {
  # This is the internal mstate function to.trans2 which converts a transition
  # matrix into a nice data frame
  #
  # Args:
  #   trans: Square transition matrix whose r, s entry is i if the ith 
  #          transition type is r, s
  # Returns:
  #   Data from with columns trasno (transition number), from (starting state),
  #   to (ending state), fromname (name of starting state), 
  #   toname (name of ending state) and transname (of the form fromname --> toname)
  dm <- dim(trans)
  if (dm[1] != dm[2]) stop("transition matrix should be square")
  S <- dm[1]
  mx <- max(trans, na.rm=TRUE)
  res <- matrix(NA, mx, 3)
  res[, 1] <- 1:mx
  transvec <- as.vector(trans)
  for (i in 1:mx) {
    idx <- which(transvec==i)
    res[i, 2:3] <- c(idx %% S, idx %/% S + 1)
  }
  res <- data.frame(res)
  names(res) <- c("transno", "from", "to")
  statesfrom <- dimnames(trans)[[1]]
  if (is.null(statesfrom)) statesfrom <- 1:S
  statesto <- dimnames(trans)[[2]]
  if (is.null(statesto)) statesto <- 1:S
  res$fromname <- statesfrom[res$from]
  res$toname <- statesto[res$to]
  res$transname <- paste(res$fromname, res$toname, sep=" -> ")
  return(res)
}