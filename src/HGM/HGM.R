# R-script for the Halibut Growth Model

require(Riscam)
require(ggplot2)
A <- read.rep("HGM.rep")

A$tx <- array(as.vector(A$data), dim=c(A$nbin, A$nage, A$narea, A$nsex))