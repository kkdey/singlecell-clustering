
######### tSNE on classtpx model probabilities   ################

library(tsne)
model <- get(load("../output/classtpx_sorted_MCC_15.rda"))
omega <- model$omega
out <- tsne::tsne(omega, k=2)
