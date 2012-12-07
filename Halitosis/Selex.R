# |---------------------------------------------------------------------------|
# | Selex.R
# |---------------------------------------------------------------------------|
# | Compute the age-specific capture probability based on size selectivity.
# | 



# |---------------------------------------------------------------------------|
# | Calculate vulnerability-at-age given P(l|a) and P(l)
# |---------------------------------------------------------------------------|
# | Args: la -> vector of mean length-at-age
# |       sa -> vector of std in length-at-age
# |       pl -> probabilty of capturing an individual of a given length l.
# |       xl -> vector of mid points for length frequency distributions.
# | Pseudocode:
# |  1. Construct an age-length key given la and sa
# |  2. Compute va = P(l)*P(l|a)
.calcPage <- function(la, sa, pl, xl)
{
	A   <- length(la)
	L   <- length(xl)
	hbw <- 0.5*(xl[2] - xl[1])
	ALK <- matrix(0, nrow=L, ncol=A)
	
	ALK <-  sapply(xl+hbw,pnorm,mean=la,sd=sa) - sapply(xl-hbw,pnorm,mean=la,sd=sa)
	
	va  <-  as.vector(pl %*% t(ALK))
	return(va)
}

.calcALK <- function(la, sa, xl)
{
	A   <- length(la)
	L   <- length(xl)
	hbw <- 0.5*(xl[2] - xl[1])
	ALK <- matrix(0, nrow=L, ncol=A)
	
	ALK <-  sapply(xl+hbw,pnorm,mean=la,sd=sa) - sapply(xl-hbw,pnorm,mean=la,sd=sa)
	
	return(ALK)
}


