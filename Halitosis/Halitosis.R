# -------------------------------------------------------------------------- ##
# Halitosis.R
# Written by Steve Martell,  IPHC
# Date: Jan 8, 2012
# DATE: November 5,  2012.  Re-organization of the code
# CODE ORGANIZATIONn
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
# 2) Maximum size limit of 140 cm
# 3) No size limits
# 4) Selectivity shift to smaller fish (-10cm)
# 5) Bycatch fishing mortality rates from 0.02 to 0.153
# 6) Size-specific natural mortality rates.
# -------------------------------------------------------------------------- ##
  setwd("/Users/stevenmartell1/Documents/IPHC/SizeLimitz/Halitosis/")
# |---------------------------------------------------------------------------|
# | Libraries
# |---------------------------------------------------------------------------| 
   require(Hmisc)
   require(ggplot2)
   require(reshape2)
   require(Riscam)
   require(grid)            #arrow function
   source("GvonBGrowth.R")  #gvonb function
   source("Selex.R")		#calcPage function
   source("IPHC_theme.R")


# |---------------------------------------------------------------------------|
# | Data and other constants
# |---------------------------------------------------------------------------|
   regArea <- "2B"
   A	<- 35					# maximum age.
   G	<- 11					# number of growth groups
   S	<- 2					# number of sexes
   dim	<- c(A, G, S)			# array dimensions
   age	<- 1:A					# vector of ages
   pg	<- dnorm(seq(-1.96, 1.96, length=G), 0, 1); 
   pg  <- pg/sum(pg) 			# proportion assigned to each growth-type group.
   fe  <- seq(0, 1.00, b=0.0051) 	#sequence of fishing mortality rates.

# |---------------------------------------------------------------------------|
# | Growth parameter estimates from vonBH
# |---------------------------------------------------------------------------|
   GM  <- read.rep("../src/VONB/vonbh.rep")
   RA  <- c("2A","2B","2C","3A","3B","4A","4B","4C","4D")
   HR  <- c(rep(0.215, 4), rep(0.161, 5))
   
# |---------------------------------------------------------------------------|
# | Commercial selectivities from Stewart 2012.                               
# |---------------------------------------------------------------------------|
   bin   <- seq(60, 130, by=10)
   CSelL <- matrix(
   (data=c(0,  0.0252252,  0.250685,  0.617268,  1,  1.36809,  1.74292,  2.12022, 
   		0,  0.0151914,  0.144236,  0.513552,  1,  1.48663,  1.97208,  2.45702)
   ), nrow=8, ncol=2)
   

   # THis is from the survey selectivity
   SSelL <- matrix(
	data=c(0,  0.285875,  0.579811,  0.822709,  1,  1.13862,  1.29717,  1.46621, 
		0,  0.246216,  0.454654,  0.68445,  1,  1.29683,  1.58379,  1.86742)
	, nrow=8, ncol=2)
	
   #CSelL <- t(t(CSelL)/apply(CSelL,2,max))
   #CSelL <- t(t(SSelL)/apply(SSelL,2,max))
   CSelL <- CSelL/max(CSelL)
   CSelL <- SSelL/max(SSelL)
   slim	<- 81.28
   ulim	<- 1500
   cvlm	<- 0.1
   
# |---------------------------------------------------------------------------|


# |---------------------------------------------------------------------------|
# | Population parameters 
# |---------------------------------------------------------------------------|
   bo		<- 100.0			# unfished female spawning biomass
   h		<- 0.95				# steepness
   dm		<- 0.16				# discard mortality rate
   cm		<- 0				# Size-dependent natural mortality rate (-0.5, 0.5)
     

# |---------------------------------------------------------------------------|
# | Assesemble Regulatory area LIST OBJECT                          
# |---------------------------------------------------------------------------|
# |
.getRegPars <- function(regArea, THETA=THETA, PI=PI)
{
	# Sex specific parameters (female, male).
	m		<- c(0.15, 0.1439)			# natural mortality rate
	#0.15 0.135474
	a50		<- rep(10.91, 2)			# age at 50% maturity
	k50		<- rep(1.406, 2)			# std at 50% maturity
	a		<- rep(6.821e-6, 2)			# length-weight allometry (Clark 1992)
	b		<- rep(3.24, 2)				# length-weight allometry (CLark 1992)
	
	
	# Get growth parameters by regArea
	j = 4
	i   = 9 
	ii <- match(regArea, RA)
	linf <- GM$linf[, ii]
	k    <- GM$vbk[, ii]
	to   <- GM$to[, ii]
	pp   <- GM$p[, ii]
	cv   <- GM$cv[, ii]
         
	THETA <- data.frame(bo=bo, h=h, dm=dm, cm=cm)

	PHI	<- data.frame(m=m, a50=a50, k50=k50, a=a, b=b, linf=linf, k=k, to=to, pp=pp, CVlinf=cv)
	rownames(PHI) <- c("Female", "Male")
	
	PI  <- data.frame(slim=slim, ulim=ulim, cvlm=cvlm)

	T1 <- c(RA=RA[ii], THETA, PHI, PI)	
	return(T1)
}


# |---------------------------------------------------------------------------|
# | Calculate survivorship, growth and fecundity (unfished)         
# |---------------------------------------------------------------------------|
# | For the length-at-age calculation, use the mean growth curve,  then 
# | base deviates from the mean length-at-age,  not linf as in the Pine paper.
.calcLifeTable <- function(RegArea)
{
	# Append to RegArea list object the following arrays:
	# survivorship     (lx)
	# length-at-age    (la)
	# weight-at-age    (wa)
	# fecundity-at-age (fa)
	# 
	
	with(RegArea, {
		lx	   <- array(1, dim)
		la	   <- array(0, dim)
		sd_la  <- array(0, dim)
		wa	   <- array(0, dim)
		fa	   <- array(0, dim)
		M      <- array(0, dim)
		ma     <- plogis(age, a50, k50)
		for(i in 1:S)
		{
			# Survivorship
			#lx[,,i]  <- exp(-m[i])^(age-1)
			#lx[A,,i] <- lx[A,,i]/(1-exp(-m[i]))
			
			# Length-at-age
			mu         <- gvonb(age, linf[i], k[i], to[i], pp[i])
			sigma      <- CVlinf[i]*mu
			dev        <- seq(-1.96, 1.96, length=G)
			if(G==1) dev <- 0
			
			
			la[,,i]    <- sapply(dev,fn<-function(dev){la=mu+dev*sigma})
			sd_la[,,i] <- sqrt(1/G*(CVlinf[i]*mu)^2)
			wa[,,i]    <- a[i]*la[,,i]^b[i]
			fa[,,i]    <- ma*wa[,,i]
			
			# Size dependent natural mortality rate 
			# M_l = M (l_a/l_r)^c
			l_r     <- 100
			delta   <- (la[,,i]/l_r)^cm / mean((la[,,i]/l_r)^cm)
			M[,,i]  <- m[i] * delta
			
			# Survivorship
			for(j in 2:A)
			{
				lx[j,,i] <- lx[j-1,,i]*exp(-M[j-1,,i])
			}
			lx[A,,i] <- lx[A,,i]/(1-exp(-M[A,,i]))
		}
		
		
		
		RegArea$lx    = lx
		RegArea$la    = la
		RegArea$sd_la = sd_la
		RegArea$wa    = wa
		RegArea$fa    = fa
		RegArea$M     = M
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
		vd	<- array(0, dim)  #discard fishery
		std	<- cvlm*slim+1.e-30	
		for(i in 1:S)
		{
			pl       <- approx(bin, CSelL[,i], xl, yright=1, yleft=0)$y
			sc[,,i]  <- .calcPage(la[,,i],sd_la[,,i],pl,xl)
			#sc[,,i]  <- approx(bin, CSelL[,i], la[,,i], yleft=0, yright=1)$y
			
			# retention proability
			pr       <-  plogis(xl, slim, 0.1) - plogis(xl, ulim, 0.1)
			sr[,,i]  <- .calcPage(la[,,i],sd_la[,,i],pr,xl)
			
			#sr[,,i]  <-  plogis(la[,,i],location=slim, scale=std) - plogis(la[,,i],location=ulim, scale=std)
			sd[,,i]  <- 1-sr[,,i]
			va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm)
			
			# discard fishery selecitvity
			pd       <- plogis(xl, 66,  0.1) - plogis(xl, 81.3, 0.1)
			vd[,,i]  <- .calcPage(la[,,i],sd_la[,,i], pd, xl)
		}
		
		
		RegArea$sc <- sc	# Length-based commercial selectivity.
		RegArea$sr <- sr	# Age-specific retention probability.
		RegArea$sd <- sd	# Age-specific discard probability.
		RegArea$va <- va	# Joint capture probability.
		RegArea$vd <- vd	# Discard probability in trawl fishery.
		
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
.asem <- function(fe=0, RegArea, ct=0)
{
	# | Psuedocode:
	# | 1. Calculate age-specific total mortality, retention, and discard rates
	# | 2. Calculate survivorship with fe>0.
	# | 3. Calculate equilibrium recruitment (re) and biomass (be)
	# | 4. Calculate yield per recruit,  spawning biomass per recruit,  yield, discards.
	# | 5. Calculate average weight-at-age.
	with(RegArea, {
		# 1. Age-specific total mortality, survival, retention, and discard rate.
		za	<- array(0, dim)
		sa	<- array(0, dim)
		qa	<- array(0, dim)
		da	<- array(0, dim)
		
		bycatch <- ct
		bapprox <- bo * 0.15/(0.15+fe)
		fd      <- bycatch/bapprox
		#if(ct>0)
		#cat("fe = ", fe, " fd = ", fd, "\n")
		
		for(i in 1:S)
		{
			za[,,i]  <- M[,,i] + fe*va[,,i] + fd*vd[,,i]
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
		dpr     <- 0
		bpr     <- 0
		spr		<- phi.e/phi.E
		for(i in 1:S)
		{
			ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			de	<- de + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			bpr <- bpr + sum( t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			dpr <- dpr + sum( fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
		}
		
		# 5. Calculate average weight-at-age
		wbar <- matrix(0, nrow=S, ncol=A)
		for(i in 1:S)
		{
			tmp      <- lz[,,i]/rowSums(lz[,,i])
			wbar[i,] <- rowSums(wa[,,i]*tmp)
		}
		
		RegArea$lz   <- lz
		RegArea$re   <- re
		RegArea$be   <- be
		RegArea$ye   <- ye
		RegArea$de   <- de
		RegArea$bpr  <- bpr
		RegArea$ypr  <- ypr
		RegArea$spr  <- spr
		RegArea$dpr  <- dpr
		RegArea$wbar <- wbar
		
		return(RegArea)
	})
}

# |---------------------------------------------------------------------------|
# | Calculate equilibrium values for a given fe vector              
# |---------------------------------------------------------------------------|
# | This function calls asem several times and constructs a data 
# | frame with columns: fe ye be re de ypr spr
.calcEquilibrium <- function(RegArea, Scenario=NULL, bycatch=0)
{
	with(RegArea, {
		# Proto-type function to get equilibrium vector.
		fn <- function(fe)
		{
			tmp <- .asem(fe, RegArea, bycatch)
			out <- c(fe=fe, ye=tmp$ye, be=tmp$be, de=tmp$de, 
				re=tmp$re, spr=tmp$spr, ypr=tmp$ypr, 
				bpr=tmp$bpr, dpr=tmp$dpr)
				
			# average weight arrays
			wbar_f <- c(wbar=tmp$wbar[1,])
			wbar_m <- c(wbar=tmp$wbar[2,])
			
			out <- c(out, wbar_f=wbar_f, wbar_m=wbar_m, Scenario=Scenario)
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
		xx      <- as.double(RegArea[[i]]$equil$spr)
		yy      <- as.double(RegArea[[i]]$equil$fe)
		fspr.30 <- approx(xx,yy,0.3,yleft=0, yright=max(yy))$y
		
		# | Find F0.1
		xx      <- as.double(RegArea[[i]]$equil$bpr)
		yy      <- as.double(RegArea[[i]]$equil$fe)
		bpr0    <- 0.1*xx[1]
		f0.1    <- approx(xx, yy, bpr0, yright=max(yy))$y
		
		# | Find Fmsy
		yy      <- as.double(RegArea[[i]]$equil$ye)
		ii      <- which.max(yy)
		fmsy    <- as.double(RegArea[[i]]$equil$fe)[ii]
		msy     <- as.double(RegArea[[i]]$equil$ye)[ii]
		
		u1 <- round(1-exp(-f0.1), 3)
		u2 <- round(1-exp(-fspr.30), 3)
		u3 <- round(1-exp(-fmsy), 3)
		if(i==1)
		cat( "\t", "Fspr30", "\t", "Umsy\n")
		cat( "\t", u2, "\t", u3, "\n")
		
		
		tmp <- data.frame(Area  = RegArea[[i]]$RA, 
						Fspr.30 = fspr.30,
						F0.1    = f0.1,  
						Fmsy    = fmsy, 
						msy     = msy, 
						RegArea[[i]]$equil)
		df  <- rbind(df, tmp)
	}
	return(df)	
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
	DF      <- NULL
	
	h       <- 0.75
	slim    <- 81.3
	ulim    <- 15000
	cm      <- 0
	dm      <- 0.16
	bin     <- seq(60, 130, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=1)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF      <- .makeDataFrame(RegArea)
	
	
	slim    <- 81.3
	ulim    <- 140
	#cm      <- 0.5
	#dm      <- 0
	#bin     <- seq(50, 120, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=2)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF2      <- .makeDataFrame(RegArea)
	
	slim    <- 0
	ulim    <- 1500
	#cm      <- -0.5
	#dm      <- 1.0
	#bin     <- seq(70, 140, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=3)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF3      <- .makeDataFrame(RegArea)
	
	h       <- 0.75
	slim    <- 81.3
	ulim    <- 15000
	cm      <- 0
	dm      <- 0.16
	bin     <- seq(50, 120, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=4)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF4      <- .makeDataFrame(RegArea)
	
	#CSelL <- t(t(SSelL)/apply(SSelL,2,max))
	#slim    <- 60
	#ulim    <- 140
	#cm      <- -0.5
	#dm      <- 1.0
	#bin     <- seq(70, 140, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=5, bycatch=2)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF5      <- .makeDataFrame(RegArea)
	
	#CSelL <- t(t(SSelL)/apply(SSelL,2,max))
	#slim    <- 60
	#ulim    <- 140
	cm      <- -0.5
	#dm      <- 1.0
	#bin     <- seq(70, 140, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=6, bycatch=0)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF6      <- .makeDataFrame(RegArea)
	
	#CSelL <- t(t(SSelL)/apply(SSelL,2,max))
	#slim    <- 60
	#ulim    <- 140
	cm      <-  0.5
	#dm      <- 1.0
	#bin     <- seq(70, 140, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=7, bycatch=0)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF7      <- .makeDataFrame(RegArea)
	
	
	#CSelL <- t(t(SSelL)/apply(SSelL,2,max))
	#slim    <- 60
	#ulim    <- 140
	cm      <- 0
	h       <- 0.85
	#dm      <- 1.0
	#bin     <- seq(70, 140, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=8, bycatch=0)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF8      <- .makeDataFrame(RegArea)
	
	#CSelL <- t(t(SSelL)/apply(SSelL,2,max))
	#slim    <- 60
	#ulim    <- 140
	cm      <- 0
	h       <- 0.65
	#dm      <- 1.0
	#bin     <- seq(70, 140, by=10)
	RegArea <- lapply(RA, .getRegPars)
	RegArea <- lapply(RegArea, .calcLifeTable)
	RegArea <- lapply(RegArea, .calcSelectivities)
	RegArea <- lapply(RegArea, .calcSRR)
	RegArea <- lapply(RegArea, .calcEquilibrium, Scenario=9, bycatch=0)
	#RegArea <- lapply(RegArea, .asem, fe=0)
	DF9      <- .makeDataFrame(RegArea)
	
	
	DF      <- rbind(DF, DF2, DF3, DF4, DF5, DF6, DF7, DF8, DF9)
	
# |
# |---------------------------------------------------------------------------|
# |---------------------------------------------------------------------------|

# |---------------------------------------------------------------------------|
# | FMSY,  MSY,  BMSY table
# |---------------------------------------------------------------------------|
# |
ss <- unique(DF$Scenario)
Fmsy <- NULL
MSY  <- NULL
Bmsy <- NULL
for(i in ss)
{
	tmp = subset(DF, Scenario==i)
	ff <- NULL
	mm <- NULL
	bb <- NULL
	for(j in RA)
	{
		tmp2 = subset(tmp, Area==j)
		ff   = c(ff, unique(tmp2$Fmsy))
		mm   = c(mm, tmp2$ye[which.max(tmp2$ye)])
		bb   = c(bb, tmp2$be[which.max(tmp2$ye)])
	}
	Fmsy <- rbind(Fmsy, ff)
	MSY  <- rbind(MSY, mm)
	Bmsy <- rbind(Bmsy, bb)
}


# |---------------------------------------------------------------------------|
# | Graphics                                                        	      |
# |---------------------------------------------------------------------------|
# | p.ypr = yield per recruit
# | p.spr = spawning biomass depletion
# | p.ye  = equilibrium yield
   RELSIZE <- 1.5
   graphics.off()
   quartz("Size at age", width=7.92, height=6.72)

# |---------------------------------------------------------------------------|
# | Plot size at age data
# |---------------------------------------------------------------------------|
# |
plotSizeAtAge <- function(RegArea)
{
	# | Construct a dataframe with Area, Sex, G as id.vars and length-at-age.
	# | 
	n  <- length(RegArea)
	df <- data.frame()
	for(i in 1:n)
	{
		la.f <- t(RegArea[[i]]$la[,,1])
		la.m <- t(RegArea[[i]]$la[,,2])

		mdf.f <- cbind(area=RegArea[[i]]$RA, sex="Female", melt(la.f))
		mdf.m <- cbind(area=RegArea[[i]]$RA, sex="Male", melt(la.m))

		df <- rbind(df, rbind(mdf.f, mdf.m))
	}
	colnames(df) <- c("area","sex", "group", "age", "length")
	p.la <- ggplot(df, aes(x=age, y=length))
	p.la <- p.la + geom_line(aes(linetype=factor(group),col=sex),size=.4, alpha=0.85) 
	p.la <- p.la + facet_wrap(~area)
	print(p.la)
	return(df)
}

   
# |---------------------------------------------------------------------------|
# | Yield per recruit                                               
# |---------------------------------------------------------------------------|
p.ypr <- ggplot(subset(DF, Scenario==c(1,2,3))) + geom_line(aes(x=fe, y=ypr, col=factor(Scenario))) 
p.ypr <- p.ypr + labs(x="Fishing mortality rate", y="Yield per recruit (lb)")
p.ypr <- p.ypr + scale_colour_discrete(name="Scenario", labels=c("Status quo", "81.3–140 cm SL","60–140 cm SL"))
#p.ypr <- p.ypr + scale_colour_discrete(name="Scenario", labels=c("Constant M", "M increases with size","M decreases with size"))
#p.ypr <- p.ypr + scale_colour_discrete(name="Scenario", labels=c("Discard mortality = 0.16", "Discard mortality = 0.00","Discard mortality = 1.00"))
p.ypr <- p.ypr + geom_segment(aes(x=F0.1, y=1, xend=F0.1, yend=0), color="darkgrey"
         , arrow=arrow(length=unit(.2, "cm")), size=0.1)
p.ypr <- p.ypr + geom_vline(xintercept=-log(1-c(0.161, 0.215)), size=0.2 )
p.ypr <- p.ypr + theme(axis.title = element_text(size = rel(RELSIZE)))
p.ypr <- p.ypr + facet_wrap(~Area)


# |---------------------------------------------------------------------------|
# | Discards per recruit                                               
# |---------------------------------------------------------------------------|
p.dpr <- ggplot(subset(DF, Scenario==c(1,2,3))) + geom_line(aes(x=fe, y=dpr, col=factor(Scenario))) 
p.dpr <- p.dpr + labs(x="Fishing mortality rate", y="Discard per recruit (lb)")
p.dpr <- p.dpr + scale_colour_discrete(name="Scenario", labels=c("Status quo", "81.3–140 cm SL","60–140 cm SL"))
p.dpr <- p.dpr + theme(axis.title = element_text(size = rel(RELSIZE)))
p.dpr <- p.dpr + facet_wrap(~Area)


# |---------------------------------------------------------------------------|
# | Percent wastage versus fishing mortality
# |---------------------------------------------------------------------------|
jj    <- seq(1,length(fe),by=15)
sDF   <- subset(DF, Scenario==c(1, 2, 4))
pDF   <- subset(sDF, fe%in%fe[jj])
p.de  <- ggplot(sDF) + geom_line(aes(x=fe, y=de/ye*100, shape=factor(Scenario), linetype=factor(Scenario)), size=0.25)
p.de  <- p.de + geom_point(data=pDF, aes(x=fe, y=de/ye*100, shape=factor(Scenario)))
p.de  <- p.de + labs(x="Fishing mortality rate", y="Percent wastage (dead discards/landed catch)", shape="Scenario", linetype="Scenario")
p.de  <- p.de + facet_wrap(~Area) + theme_bw(12)
ggsave(p.de, file="../FIGS/fig:PercentWastage.pdf")
ggsave(p.de, file="../FIGS/fig:PercentWastage.png")


# |---------------------------------------------------------------------------|
# | Yield per recruit & current harvest policy.                                              
# |---------------------------------------------------------------------------|
f.ypr <- ggplot(subset(DF, Scenario==c(1))) + geom_line(aes(x=fe, y=ypr, col=Area), size=1.5)
f.ypr <- f.ypr + geom_vline(xintercept=-log(1-c(0.161, 0.215)), size=1 )
f.ypr <- f.ypr + labs(x="Fishing mortality rate", y="Yield per recruit (lb)")
f.ypr <- f.ypr + theme(axis.title = element_text(size = rel(RELSIZE)))



# |---------------------------------------------------------------------------|
# | Spawning biomass depletion                                      
# |---------------------------------------------------------------------------|
p.spr <- ggplot(subset(DF,Scenario==c(1, 2, 3)) ) + geom_line(aes(x=fe, y=spr, col=factor(Scenario))) 
p.spr <- p.spr + labs(x="Fishing mortality rate", y="Spawning biomass depletion")
p.spr <- p.spr + scale_colour_discrete(name="Scenario", labels=c("Status quo", "81.3–140 cm SL","60–140 cm SL"))
#p.spr <- p.spr + scale_colour_discrete(name="Scenario", labels=c("Constant M", "M increases with size","M decreases with size"))
#p.spr <- p.spr + scale_colour_discrete(name="Scenario", labels=c("Discard mortality = 0.16", "Discard mortality = 0.00","Discard mortality = 1.00"))
p.spr <- p.spr + geom_segment(aes(x=Fspr.30, y=0.30, xend=Fspr.30, yend=0), col="red"
         , arrow=arrow(length=unit(.2, "cm")), size=0.1)
p.spr <- p.spr + geom_vline(xintercept=-log(1-c(0.161, 0.215)), size=0.2 )
p.spr <- p.spr + theme(axis.title = element_text(size = rel(RELSIZE)))
p.spr <- p.spr + facet_wrap(~Area)



# |---------------------------------------------------------------------------|
# | Equilibrium Yield                                               
# |---------------------------------------------------------------------------|
p.ye <- ggplot(subset(subset(DF, Scenario==1), Area=="2B")) + geom_line(aes(x=fe, y=ye), size=1.5)
p.ye <- p.ye + labs(x="Fishing Mortality Rate", y="Equilibrium Yield") + ylim(c(0, 8))
p.ye2B <- p.ye
SIZE <- 1.5
jj   <- seq(1,length(fe),by=10)

DF   <- subset(DF,Area==c("2B","3A","4A")) #Comment out for paper.
sDF  <- subset(DF, Scenario==1)
pDF  <- subset(sDF, fe%in%fe[jj])
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, shape=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(data=pDF, aes(x=fe, y=ye, shape=Area), size=2)
p.ye <- p.ye + scale_shape_manual(values=c(1:9))
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, shape=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye1 <- p.ye + theme_bw(12)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(data=pDF, aes(x=fe, y=ye, col=Area), size=2)
p.ye <- p.ye + scale_shape_manual(values=c(1:9))
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye1a<-p.ye
ggsave(p.ye1, file="../FIGS/fig:YeBase.pdf")
ggsave(p.ye1, file="../FIGS/fig:YeBase.png")


sDF  <- subset(DF, Scenario==2)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye2 <- p.ye

sDF  <- subset(DF, Scenario==3)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye3 <- p.ye

sDF  <- subset(DF, Scenario==4)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye4 <- p.ye 

sDF  <- subset(DF, Scenario==5)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye5 <- p.ye 

sDF  <- subset(DF, Scenario==6)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye6 <- p.ye 

sDF  <- subset(DF, Scenario==7)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye7 <- p.ye 

sDF  <- subset(DF, Scenario==8)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye8 <- p.ye

sDF  <- subset(DF, Scenario==9)
p.ye <- ggplot(sDF ) + geom_line(aes(x=fe, y=ye, col=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(aes(x=fe, y=ye, shape=Area, col=Area), size=2)
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, col=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.ye9 <- p.ye


pDF  <- subset(DF, fe%in%fe[jj])
p.ye <- ggplot(DF) + geom_line(aes(x=fe, y=ye, shape=Area), size=SIZE) + ylim(c(0, 8))
p.ye <- p.ye + geom_point(data=pDF, aes(x=fe, y=ye, shape=Area), size=1.5)
p.ye <- p.ye + scale_shape_manual(values=c(1:9))
p.ye <- p.ye + labs(x="Fishing mortality rate", y="Relative yield")
p.ye <- p.ye + geom_segment(aes(x=Fmsy, y=msy, xend=Fmsy, yend=0, shape=Area), arrow=arrow(length=unit(.2, "cm")), size=0.25, linetype=1)
p.yeAll <- p.ye + facet_wrap(~Scenario) + theme_bw(12)
ggsave(p.yeAll, file="../FIGS/fig:YeALL.pdf")
ggsave(p.yeAll, file="../FIGS/fig:YeALL.png")


#p.ye <- p.ye + geom_vline(xintercept=-log(1-c(0.161, 0.215)), size=0.2,  col=c(1, 2) )
#p.ye <- p.ye #+ facet_wrap(~Scenario) 

#p.ye <- p.ye + scale_colour_discrete(name="Scenario", labels=c("h=0.95","h=0.75"))
#p.ye <- p.ye + scale_colour_discrete(name="Scenario", labels=c("Commercial Selectivity", "Survey Selectivity"))
#p.ye <- p.ye + scale_colour_discrete(name="Scenario", labels=c("Status quo", "81.3–140 cm SL","60–140 cm SL"))
#p.ye <- p.ye + scale_colour_discrete(name="Scenario", labels=c("Discard mortality = 0.16", "Discard mortality = 0.00","Discard mortality = 1.00"))
#p.ye <- p.ye + theme(axis.text.x = element_text(size = 2))

# |---------------------------------------------------------------------------|
# | Mean weight-at-age vs F  for females                                             
# |---------------------------------------------------------------------------|
   iage        <- c(6, 12, 14, 20)
   Scen        <- 1
   nm          <- colnames(DF)
   wbar_f_cols <- nm %in% grep("^wbar_f",nm,value=TRUE)
   sDF         <- cbind(Area=DF$Area, fe=DF$fe,sex="Female",S=DF$Scenario, DF[,wbar_f_cols])
   sDF         <- subset(sDF, S==Scen)
   colnames(sDF) <- c("Area", "fe", "sex", "Scenario", paste(1:A) )
   icol          <- c(1, 2, 3, match(iage, names(sDF)))
   mDF           <- melt(sDF[icol],id.vars=c("Area","fe","sex"))

   p.wbar_f <- ggplot(mDF,aes(x=fe,y=value,group=variable, linetype=variable)) 
   p.wbar_f <- p.wbar_f + geom_line(size=0.5) +xlim(c(0, 0.75))
   p.wbar_f <- p.wbar_f + labs(x="Fishing mortality rate", y="Mean weight-at-age (lb)")
   p.wbar_f <- p.wbar_f + labs(linetype="Age")
   p.wbar_f <- p.wbar_f + facet_wrap(~Area) + theme_bw(12)
   ggsave(p.wbar_f, file="../FIGS/fig:wbar_female.pdf")
   ggsave(p.wbar_f, file="../FIGS/fig:wbar_female.png")


# |---------------------------------------------------------------------------|
# | Mean weight-at-age vs F  for males                                             
# |---------------------------------------------------------------------------|
   nm          <- colnames(DF)
   wbar_m_cols <- nm %in% grep("^wbar_m",nm,value=TRUE)
   sDF         <- cbind(Area=DF$Area, fe=DF$fe, sex="Male",S=DF$Scenario, DF[,wbar_m_cols])
   sDF         <- subset(sDF, S==Scen)
   colnames(sDF) <- c("Area", "fe", "sex", paste(1:A) )
   icol          <- c(1, 2, 3, match(iage, names(sDF)))
   mDF           <- melt(sDF[icol],id.vars=c("Area","fe","sex"))

   p.wbar_m <- ggplot(mDF,aes(x=fe,y=value,group=variable, col=variable))+geom_line()
   p.wbar_m <- p.wbar_m + labs(x="Fishing mortality rate", y="Mean weight-at-age (lb)")
   p.wbar_m <- p.wbar_m + theme(axis.title = element_text(size = rel(RELSIZE)))
   p.wbar_m <- p.wbar_m + labs(color="Age") 
   p.wbar_m <- p.wbar_m + facet_wrap(~Area) 

