# -------------------------------------------------------------------------- ##
# Halitosis.R
# Written by Steve Martell,  IPHC
# Date: Jan 8, 2012
# DATE: November 5,  2012.  Re-organization of the code
# CODE ORGANIZATION
# 	Dependencies.
#   Read in external results from growth parameter estimatation,  selectivities.
#   PSEUDOCODE:
# 				1) Construct a LIST object for each regulatory area with parameters. (RegArea)
# 				2) For each RegArea
# 					- append life-table array (length, weight, fecundity)
# 
# 
# Feb 7, 2011.  Added multiple growth groups to allow for cumulative effects
# of size selective fishing.  Need to show how mean weigth-at-age decreases
# with increasing fishing mortality rates.
# 
# Halibut prices (10-20) (20-40) (40+)
# $6.75  $7.30  $7.50  In Homer Alaska.
#
# SCENARIOS FOR HALIBUT PAPER
# 1) Status quo. 83.3cm minimum size limit
# 2) No size limits
# 3) 70cm minimum size limit
# 4) slot limit 81.3-150 cm
# 5) slot limit 70-150 cm
# -------------------------------------------------------------------------- ##

# |---------------------------------------------------------------------------|
# | Libraries
# |---------------------------------------------------------------------------| 
   require(Hmisc)
   require(ggplot2)
   require(reshape2)
   require(Riscam)
   source("GvonBGrowth.R")  #gvonb function
   source("Selex.R")		#calcPage function


# |---------------------------------------------------------------------------|
# | Data and other constants
# |---------------------------------------------------------------------------|
   regArea <- "2B"
   A	<- 35					# maximum age.
   G	<- 11					# number of growth groups
   S	<- 2					# number of sexes
   dim	<- c(A, G, S)			# array dimensions
   age	<- 1:A					# vector of ages
   pg	<- dnorm(seq(-3, 3, length=G), 0, 1); 
   pg  <- pg/sum(pg) 			# proportion assigned to each growth-type group.
   fe  <- seq(0, 1.6, b=0.005) 	#sequence of fishing mortality rates.

# |---------------------------------------------------------------------------|
# | Growth parameter estimates from vonBH
# |---------------------------------------------------------------------------|
   GM  <- read.rep("../src/VONB/vonbh.rep")
   RA  <- c("2A","2B","2C","3A","3B","4A","4B","4C","4D")
   
# |---------------------------------------------------------------------------|
# | Commercial selectivities from Stewart 2012.                               
# |---------------------------------------------------------------------------|
   bin   <- seq(60, 130, by=10)
   CSelL <- matrix(
   (data=c(0,  0.0252252,  0.250685,  0.617268,  1,  1.36809,  1.74292,  2.12022, 
   		0,  0.0151914,  0.144236,  0.513552,  1,  1.48663,  1.97208,  2.45702)
   ), nrow=8, ncol=2)
   
   CSelL <- t(t(CSelL)/apply(CSelL,2,max))
   slim	<- 81.28
   ulim	<- 1500
   cvlm	<- 0.1
   PI		<- data.frame(lhat=lhat, ghat=ghat, slim=slim, ulim=ulim, cvlm=cvlm)
# |---------------------------------------------------------------------------|


# |---------------------------------------------------------------------------|
# | Population parameters 
# |---------------------------------------------------------------------------|
   bo		<- 100.0			# unfished female spawning biomass
   h		<- 0.75				# steepness
   dm		<- 0.16				# discard mortality rate
   THETA <- data.frame(bo=bo, h=h, dm=dm)
   

# |---------------------------------------------------------------------------|
# | Assesemble Regulatory area LIST OBJECT                          
# |---------------------------------------------------------------------------|
# |
.getRegPars <- function(regArea)
{
	# Sex specific parameters (female, male).
	m		<- c(0.15, 0.1439)			# natural mortality rate
	#0.15 0.135474
	a50		<- rep(10.91, 2)			# age at 50% maturity
	k50		<- rep(1.406, 2)			# std at 50% maturity
	a		<- rep(6.821e-6, 2)			# length-weight allometry (Clark 1992)
	b		<- rep(3.24, 2)				# length-weight allometry (CLark 1992)
	
	
	# Get growth parameters by regArea
	ii <- match(regArea, RA)
	linf <- GM$linf[, ii]
	k    <- GM$vbk[, ii]
	to   <- GM$to[, ii]
	pp   <- GM$p[, ii]
	cv   <- GM$cv[, ii]
         

	PHI	<- data.frame(m=m, a50=a50, k50=k50, a=a, b=b, linf=linf, k=k, to=to, pp=pp, CVlinf=cv)
	rownames(PHI)=c("Female", "Male")

	T1 <- c(RA=RA[ii], THETA, PHI, PI)	
	return(T1)
}


# |---------------------------------------------------------------------------|
# | Calculate survivorship, growth and fecundity (unfished)         
# |---------------------------------------------------------------------------|
# |
.calcLifeTable <- function(RegArea)
{
	# Append to RegArea list object the following arrays:
	# survivorship     (lx)
	# length-at-age    (la)
	# weight-at-age    (wa)
	# fecundity-at-age (fa)
	# 
	
	with(RegArea, {
		lx	   <- array(0, dim)
		la	   <- array(0, dim)
		sd_la  <- array(0, dim)
		wa	   <- array(0, dim)
		fa	   <- array(0, dim)
		ma     <- plogis(age, a50, k50)
		for(i in 1:S)
		{
			# Survivorship
			lx[,,i]  <- exp(-m[i])^(age-1)
			lx[A,,i] <- lx[A,,i]/(1-exp(-m[i]))
			
			# Length-at-age
			dev        <- linf[i]*CVlinf[i]
			linf.g     <- seq(linf[i]-dev, linf[i]+dev, length=G)
			la[,,i]    <- sapply(linf.g,gvonb,t=age,vbk=k[i],to=to[i],p=pp[i])
			sd_la[,,i] <- sqrt(1/G*(CVlinf[i]*la[,,i])^2)
			wa[,,i]    <- a[i]*la[,,i]^b[i]
			fa[,,i]    <- ma*wa[,,i]
		}
		
		
		
		RegArea$lx    = lx
		RegArea$la    = la
		RegArea$sd_la = sd_la
		RegArea$wa    = wa
		RegArea$fa    = fa
		return(RegArea)
	})
	
}

# |---------------------------------------------------------------------------|
# | Calculate size-based selectivities and joint capture probability
# |---------------------------------------------------------------------------|
# | This function calls .calcPage(la, sa, pl, xl) in Selex.R
.calcSelectivities <- function(RegArea)
{
	# TODO Fix selectivities such that va is roughly the same with G=1 or G=11 groups.
	# Calculate capture probabilities
	with(RegArea, {
		# Length-interval midpoints for integration
		xl  <- seq(5,200,by=2.5)
		
		
		# Length-based selectivity (length-based -> age-based)
		sc	<- array(0, dim)
		sr	<- array(0, dim)
		sd	<- array(0, dim)
		va	<- array(0, dim)
		std	<- cvlm*slim+1.e-30	
		for(i in 1:S)
		{
			pl       <- approx(bin, CSelL[,i], xl, yright=1, yleft=0)$y
			sc[,,i]  <- .calcPage(la[,,i],sd_la[,,i],pl,xl)
			#sc[,,i]  <- approx(bin, CSelL[,i], la[,,i], yleft=0, yright=1)$y
			
			
			sr[,,i]  <-  plogis(la[,,i],location=slim, scale=std) 
					    - plogis(la[,,i],location=ulim, scale=std)
			sd[,,i]  <- 1-sr[,,i]
			va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm)
		}
		
		RegArea$sc <- sc	# Length-based commercial selectivity.
		RegArea$sr <- sr	# Age-specific retention probability.
		RegArea$sd <- sd	# Age-specific discard probability.
		RegArea$va <- va	# Joint capture probability.
		
		return(RegArea)
	})
}

# |---------------------------------------------------------------------------|
# | Calculate stock recruitment relationship.                       
# |---------------------------------------------------------------------------|
# |
.calcSRR <- function(RegArea)
{
	with(RegArea, {
		# Unfished SPR  (phi.E)
		phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg))
		
		# Unfished recruitment (ro)
		ro		<- bo/phi.E
		
		# Beverton-Holt model
		kap <- 4*h/(1-h)

		# Ricker Model
		# kap <- (5*h)^(5/4)
		
		RegArea$phi.E <- phi.E
		RegArea$ro    <- ro
		RegArea$kap   <- kap
		
		return(RegArea)
	})
}

# |---------------------------------------------------------------------------|
# | Age-structure equilibrium model asem                            
# |---------------------------------------------------------------------------|
# | fe is the equilibrium fishing mortality rate.
.asem <- function(fe=0, RegArea)
{
	# | Psuedocode:
	# | 1. Calculate age-specific total mortality, retention, and discard rates
	# | 2. Calculate survivorship with fe>0.
	# | 3. Calculate equilibrium recruitment (re) and biomass (be)
	# | 4. Calculate yield per recruit,  spawning biomass per recruit,  yield, discards.
	with(RegArea, {
		# 1. Age-specific total mortality, survival, retention, and discard rate.
		za	<- array(0, dim)
		sa	<- array(0, dim)
		qa	<- array(0, dim)
		da	<- array(0, dim)
		for(i in 1:S)
		{
			za[,,i]  <- m[i] + fe*va[,,i]
			sa[,,i]  <- exp(-za[,,i])
			qa[,,i]  <- (sc[,,i]*sr[,,i]) * (1-sa[,,i])/za[,,i]
			da[,,i]  <- (sc[,,i]*sd[,,i]) * (1-sa[,,i])/za[,,i]
		}
		
		# 2. Survivorship under fished conditions lz(A, G, S)
		lz	<- array(1, dim)
		for(i in 1:S)
		{
			for(j in 2:A)
			{
				lz[j,,i] <- lz[j-1,,i]*exp(-za[j-1,,i])
			}
			lz[A,,i] <- lz[A,,i]/(1-exp(-za[A,,i]))
		}
		
		# 3. Calculate equilibrium recruitment and biomass
		phi.e	<- sum( t(lz[,,1]*fa[,,1])*pg )
		
		# Beverton-Holt model
		t1      <- phi.E/phi.e
		t2      <- (kap-t1)
		re      <- max(0, ro*t2/(kap-1))
		
		# Ricker model
		#t1		<- log(phi.E/(kap*phi.e))
		#t2		<- (log(kap)*phi.e)
		#re		<- max(0, -(t1*ro*phi.E)/t2)
		
		be		<- re * phi.e
		
		# 4. Calculate yield per recruit,  spawning biomass per recruit,  yield, discards.
		ye		<- 0
		de		<- 0
		ypr		<- 0
		bpr     <- 0
		spr		<- phi.e/phi.E
		for(i in 1:S)
		{
			ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			de	<- de + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			bpr <- bpr + sum( t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
		}
		
		RegArea$lz  <- lz
		RegArea$re  <- re
		RegArea$be  <- be
		RegArea$ye  <- ye
		RegArea$de  <- de
		RegArea$bpr <- bpr
		RegArea$ypr <- ypr
		RegArea$spr <- spr
		
		return(RegArea)
	})
}

# |---------------------------------------------------------------------------|
# | Calculate equilibrium values for a given fe vector              
# |---------------------------------------------------------------------------|
# | This function calls asem several times and constructs a data 
# | frame with columns: fe ye be re de ypr spr
.calcEquilibrium <- function(RegArea)
{
	with(RegArea, {
		# Proto-type function to get equilibrium vector.
		fn <- function(fe)
		{
			tmp <- .asem(fe, RegArea)
			out <- c(fe=fe, ye=tmp$ye, be=tmp$be, de=tmp$de, 
				re=tmp$re, spr=tmp$spr, ypr=tmp$ypr, bpr=tmp$bpr)
			return(out)
		}
		
		xx <- sapply(fe, fn)
		RegArea$equil <- as.data.frame(t(xx))
		
		return(RegArea)
	})
}

# |---------------------------------------------------------------------------|
# | Construct a dataframe from all RegArea objects for ggplot use   
# |---------------------------------------------------------------------------|
# |
.makeDataFrame <- function(RegArea)
{
	
	n  <- length(RegArea)
	df <- data.frame()
	for(i in 1:n)
	{
		# | Find Fspr = 0.3
		xx      <- RegArea[[i]]$equil$spr
		yy      <- RegArea[[i]]$equil$fe
		fspr.30 <- approx(xx,yy,0.3)$y
		
		# | Find F0.1
		xx      <- RegArea[[i]]$equil$bpr
		yy      <- RegArea[[i]]$equil$fe
		bpr0    <- 0.1*xx[1]
		f0.1    <- approx(xx, yy, bpr0, yright=max(yy))$y
		
		
		tmp <- data.frame(Area  = RegArea[[i]]$RA, 
						Fspr.30 = fspr.30,
						F0.1    = f0.1,  
						RegArea[[i]]$equil)
		df  <- rbind(df, tmp)
	}
	return(df)	
}

# |---------------------------------------------------------------------------|
# | Plot size at age data
# |---------------------------------------------------------------------------|
# |
plotSizeAtAge <- function(RegArea)
{
	# | Construct a dataframe with Area, Sex, G as id.vars and length-at-age
	
}

# |---------------------------------------------------------------------------|
# | Main function calls.                                                      |
# |---------------------------------------------------------------------------|
# | The order of operations is important here:
# | 1. .getRegPars
# | 2. .calcLifeTable
# | 3. .calcSelectivities
# | 4. .calcSRR (Stock recruitment relationship)
# | 5. .calcEquilibrium (data frame of values versus fe)
# | 6. .makeDataFrame (object for ggplot)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium)
	DF      <- .makeDataFrame(RegArea)
# |
# |---------------------------------------------------------------------------|
# |---------------------------------------------------------------------------|



# |-----------------------------------------------------------------|
# | Graphics                                                        |
# |-----------------------------------------------------------------|
# | p.ypr = yield per recruit
# | p.spr = spawning biomass depletion
# | p.ye  = equilibrium yield

# |-----------------------------------------------------------------|
# | Yield per recruit                                               |
# |-----------------------------------------------------------------|
p.ypr <- ggplot(DF, aes(x=fe, y=ypr)) + geom_line() 
p.ypr <- p.ypr + labs(x="Fishing mortality rate", y="Yield per recruit (lb)")
p.ypr <- p.ypr + facet_wrap(~Area)

# |-----------------------------------------------------------------|
# | Spawning biomass depletion                                      |
# |-----------------------------------------------------------------|
p.spr <- ggplot(DF, aes(x=fe, y=spr)) + geom_line() 
p.spr <- p.spr + labs(x="Fishing mortality rate", y="Spawning biomass depletion")
p.spr <- p.spr + geom_segment(aes(x = Fspr.30, y = 0.30, xend = Fspr.30, yend = 0), col="red")
p.spr <- p.spr + facet_wrap(~Area)

# |-----------------------------------------------------------------|
# | Equilibrium Yield                                               |
# |-----------------------------------------------------------------|
p.ye <- ggplot(DF, aes(x=fe, y=ye)) + geom_line() 
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + facet_wrap(~Area)




