# TRANSITION PROBABILITY PLOT --------------------------------------------------
tpPlot <- function(probs, treat.names=NULL, state.names=NULL){
  # Plots probability of being in each state over time
  #
  # Args:
  #   probs: list of dataframes with left-most column as time and remaining
  #          columns as probability of being in each possible state
  #   treat.names: names of treatment associates with each dataframe in probs
  #   state.names: name of each state in the multi-state model
  #
  # Returns:
  #   ggplot figure
  tp <- do.call("rbind", probs)
  ns <- ncol(probs[[1]]) - 1
  nr <- do.call("c", lapply(probs, function(x) nrow(x)))
  if (length(treat.names) == 0){
    treat.names <- paste0("Treatment ", seq(1, length(probs)))
  } 
  tp$treat <- rep(treat.names, times = nr)
  if (length(state.names) == 0){
    state.names = paste0("State ", seq(1, ns))
  }
  colnames(tp) <- c("time", state.names, "Treatment")
 
  tp <- melt(tp, id = c("time", "Treatment"), value.name = "pstate", 
                    variable.name = "State")
  p <-ggplot(tp, aes(x = time, y = pstate, col = Treatment)) +
    geom_line() + facet_wrap(~State) + xlab("Years") + 
    ylab("Probability of being in state") +
    theme(legend.position="bottom")
  return(p)
}

# CUMULATIVE HAZARDS PLOT ------------------------------------------------------
cumhazPlot <- function(cumhaz, model.names, trans){
  x <- do.call("rbind", cumhaz)
  nmods <- length(cumhaz)
  n <- sapply(cumhaz, nrow)
  x$Model <- rep(model.names, times = n)
  trans.dat <- to.trans(trans)
  x$trans <- factor(x$trans, levels = seq(1, max(x$trans)),
                   labels = trans.dat$transname)
  ggplot(x, aes(x = time, y = Haz, col = Model, linetype = Model)) + 
    geom_line() + facet_wrap(~trans, scale = "free_y") + xlab("Years") +
    ylab("Cumulative hazard") +
    scale_linetype_manual(values = c(rep("solid", 1), rep("dashed", nmods -1))) +
    theme(legend.position="bottom", strip.text = element_text(size = 12))
}