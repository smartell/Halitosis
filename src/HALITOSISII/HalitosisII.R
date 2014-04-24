








  setwd("/Users/stevenmartell1/Documents/IPHC/SizeLimitz/SRC/HALITOSISII/")
# |---------------------------------------------------------------------------|
# | Libraries
# |---------------------------------------------------------------------------| 
   require(Hmisc)
   require(ggplot2)
   require(reshape2)
   require(Riscam)
   require(grid)            #arrow function

   source("Selex.R")

# |---------------------------------------------------------------------------|
# | Data and other constants
# |---------------------------------------------------------------------------|
# | Stock -> is a list object with all parameters and output.
   regArea <- "2B"
   A	<- 35					# maximum age.
   G	<- 11					# number of growth groups
   S	<- 2					# number of sexes
   dim	<- c(A, G, S)			# array dimensions
   age	<- 1:A					# vector of ages
   pg	<- dnorm(seq(-1.96, 1.96, length=G), 0, 1); 
   pg  <- pg/sum(pg) 			# proportion assigned to each growth-type group.
   fe  <- seq(0, 2.0, length=51)#sequence of fishing mortality rates.
   Stock <- list(A=A,G=G,S=S,dim=dim,age=age,pg=pg)


# |---------------------------------------------------------------------------|
# | Population parameters 
# |---------------------------------------------------------------------------|
   bo		<- 300.0			# unfished female spawning biomass
   h		<- 0.75				# steepness
   m        <- c(0.15,0.14)     # natural mortality rate
   linf     <- c(165,107)		# asymptotic length
   vonk     <- c(0.075,0.068)	# vonk
   to       <- c(-4.62,-7.45)   # time at zero length
   p        <- c(1.49,0.98)     # vonb Power parameter.
   cv       <- c(0.1,0.1)		# CV in length-at-age.
   a50      <- rep(10.91,2)     # age at 50% maturity.
   k50      <- rep(1.406,2)     # std at 50% maturity.
	a		<- rep(6.821e-6, 2) # length-weight allometry (Clark 1992)
	b		<- rep(3.24, 2)		# length-weight allometry (CLark 1992)

   # dm		<- 0.16				# discard mortality rate
   cm		<- 0				# Size-dependent natural mortality rate (-0.5, 0.5)
   Stock <- c(Stock,list(bo=bo,h=h,m=m,linf=linf,vonk=vonk,to=to,p=p,cv=cv,a50=a50,k50=k50))
   Stock <- c(Stock,list(a=a,b=b,cm=cm))

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
   CSelL <- t(t(SSelL)/apply(SSelL,2,max))
   #CSelL <- CSelL/max(CSelL)
   #CSelL <- SSelL/max(SSelL)
   slim	<- 81.28
   ulim	<- 1500
   cvlm	<- 0.1

   sizelimit <- 2.54 * c(seq(20,65,by=5))
   discardmortality <- c(0.15,0.30,0.60,0.90)
   bycatch  <- c(0,10,20)
# |---------------------------------------------------------------------------|


# |---------------------------------------------------------------------------|
# | Age Schedule Information.
# |---------------------------------------------------------------------------|
.calcLifeTable <- function( Stock ) 
{

	gvonb <- function(t, linf, vbk, to, p=1)
	{
		l <- linf*(1.0-exp(-vbk*(t - to)))^p
		
		return(l)
	}

	with(Stock, {
		lx	   <- array(1, dim)
		la	   <- array(0, dim)
		sd_la  <- array(0, dim)
		wa	   <- array(0, dim)
		fa	   <- array(0, dim)
		M      <- array(0, dim)
		ma     <- plogis(age, a50, k50)

		for(i in 1:S)
		{
			# Lenght-at-age
			mu    <- gvonb(age,linf[i],vonk[i],to[i],p[i])
			sigma <- cv[i] * mu
			dev   <- seq(-1.96, 1.96, length=G)
			if(G==1) dev <- 0

			la[,,i]    <- sapply(dev,fn<-function(dev){la=mu+dev*sigma})
			sd_la[,,i] <- sqrt(1/G*(cv[i]*mu)^2)
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

		Stock$lx	   = lx	
		Stock$la	   = la	
		Stock$sd_la    = sd_la
		Stock$wa	   = wa	
		Stock$fa	   = fa	
		Stock$M        = M
		Stock$ma       = ma

		return(Stock)
	})

}


# |---------------------------------------------------------------------------|
# | Calculate size-based selectivities and joint capture probability
# |---------------------------------------------------------------------------|
# | This function calls .calcPage(la, sa, pl, xl) in Selex.R
.calcSelectivities <- function(Stock,slim=0,ulim=1500,cvlm=0.1,dm=0)
{
	# TODO Fix selectivities such that va is roughly the same with G=1 or G=11 groups.
	# Calculate capture probabilities
	with(Stock, {
		# Length-interval midpoints for integration
		xl  <- seq(5,200,by=2.5)
		
		# Length-based selectivity (length-based -> age-based)
		sc	<- array(0, dim)  #size capture
		sr	<- array(0, dim)  #size retention
		sd	<- array(0, dim)  #size discard
		va	<- array(0, dim)  #retained
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
		
		
		Stock$sc <- sc	# Length-based commercial selectivity.
		Stock$sr <- sr	# Age-specific retention probability.
		Stock$sd <- sd	# Age-specific discard probability.
		Stock$va <- va	# Joint capture probability.
		Stock$vd <- vd	# Discard probability in trawl fishery.
		Stock$dm <- dm  # Discard mortality rate
		
		return(Stock)
	})
}


# |---------------------------------------------------------------------------|
# | Calculate stock recruitment relationship.                       
# |---------------------------------------------------------------------------|
# |
.calcSRR <- function(Stock)
{
	with(Stock, {
		# Unfished SPR  (phi.E)
		phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg))
		
		# Unfished recruitment (ro)
		ro		<- bo/phi.E
		
		# Beverton-Holt model
		kap <- 4*h/(1-h)

		# Ricker Model
		# kap <- (5*h)^(5/4)
		
		Stock$phi.E <- phi.E
		Stock$ro    <- ro
		Stock$kap   <- kap
		
		return(Stock)
	})
}

# |---------------------------------------------------------------------------|
# | Age-structure equilibrium model asem                            
# |---------------------------------------------------------------------------|
# | fe is the equilibrium fishing mortality rate.
.asem <- function(fe=0, Stock, ct=0)
{
	# | Psuedocode:
	# | 1. Calculate age-specific total mortality, retention, and discard rates
	# | 2. Calculate survivorship with fe>0.
	# | 3. Calculate equilibrium recruitment (re) and biomass (be)
	# | 4. Calculate yield per recruit,  spawning biomass per recruit,  yield, discards.
	# | 5. Calculate average weight-at-age.
	with(Stock, {
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
		ne      <- 0
		de		<- 0
		ypr		<- 0
		dpr     <- 0
		bpr     <- 0
		spr		<- phi.e/phi.E
		for(i in 1:S)
		{

			ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ne  <- ne + sum( re * fe * t(lz[,,i]*qa[,,i])*pg )
			de	<- de + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			bpr <- bpr + sum( t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			dpr <- dpr + sum( fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
		}
		cbar <- ye / ne
		
		# 5. Calculate average weight-at-age
		wbar <- matrix(0, nrow=S, ncol=A)
		for(i in 1:S)
		{
			tmp      <- lz[,,i]/rowSums(lz[,,i])
			wbar[i,] <- rowSums(wa[,,i]*tmp)
		}
		
		# 6. Calculate average weight of the landed catch.

		Stock$lz   <- lz
		Stock$re   <- re
		Stock$be   <- be
		Stock$ye   <- ye
		Stock$de   <- de
		Stock$bpr  <- bpr
		Stock$ypr  <- ypr
		Stock$spr  <- spr
		Stock$dpr  <- dpr
		Stock$wbar <- wbar
		Stock$cbar <- cbar
		
		return(Stock)
	})
}

# |---------------------------------------------------------------------------|
# | Calculate equilibrium values for a given fe vector              
# |---------------------------------------------------------------------------|
# | This function calls asem several times and constructs a data 
# | frame with columns: fe ye be re de ypr spr
.calcEquilibrium <- function(Stock, Scenario=NULL, bycatch=0)
{
	with(Stock, {
		# Proto-type function to get equilibrium vector.
		fn <- function(fe)
		{
			tmp <- .asem(fe, Stock, bycatch)
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
		Stock$equil <- as.data.frame(t(xx))
		
		return(Stock)
	})
}

# |---------------------------------------------------------------------------|
# | Run Model and Construct Data Frame for fe, slim, dm, bycatch combinations.
# |---------------------------------------------------------------------------|
# | This function creates a large data frame object with id.vars for fe, slim
# | dm and bycatch levels, and colums for each performance measure (ye, spr)
# | The order of operations is important here:
# | 2. .calcLifeTable
# | 3. .calcSelectivities
# | 4. .calcSRR (Stock recruitment relationship)
# | 5. .calcEquilibrium (data frame of values versus fe)

.runModel <- function(df)
{
	
	with(as.list(df), {
		print(dm)
		M1 <- .calcLifeTable(Stock)
		M1 <- .calcSelectivities(M1,slim=slim,ulim=1500,cvlm=0.1,dm=dm)
		M1 <- .calcSRR(M1)
		M1 <- .asem(fe,M1,bycatch)

		out <- c(fe=fe,slim=slim,dm=dm,bycatch=bycatch,
		         Ye=M1$ye,De=M1$de,Be=M1$be,Re=M1$re,
		         SPR=M1$spr,YPR=M1$ypr,DPR=M1$dpr,BPR=M1$bpr,
		         Cbar=M1$cbar )
		return(out)
	})
}


# |---------------------------------------------------------------------------|
# | Main function calls.                                                      |
# |---------------------------------------------------------------------------|
# | df is a data frame with the output of the equilibrium model.
# | 
	df <- expand.grid(fe=fe,slim=sizelimit,dm=discardmortality,bycatch=bycatch)

	if(!exists("S1"))
	S1 <- as.data.frame(t(apply(df,1,.runModel)))

	
# |---------------------------------------------------------------------------|
# | GRAPHICS.                                                      |
# |---------------------------------------------------------------------------|
# | p1 -> equilibrium yield
# | p2 -> equilibrium wastage
# | p3 -> Efficiency defined as yield /	(yield + wastage)
# | p4 -> spr
# | p5 -> ypr
# | p6 -> Average weight in the catch.

	p  <- ggplot(S1,aes(x=fe,y=slim,z=Ye))
	p  <- p + stat_contour(geom="polygon",aes(fill = ..level..))
	p  <- p + facet_grid(dm~bycatch)
	p  <- p + labs(x="Fishing intensity",y="Minumum Size Limit (cm)",color="Yield\n(lbs)")
	p1 <- p

	p  <- ggplot(S1,aes(x=fe,y=slim,z=De))
	p  <- p + stat_contour(geom="polygon",aes(fill = ..level..))
	p  <- p + facet_grid(dm~bycatch)
	p  <- p + labs(x="Fishing intensity",y="Minumum Size Limit (cm)",color="Wastage\n(lbs)")
	p2 <- p

	p  <- ggplot(S1,aes(x=fe,y=slim,z=(100*Ye/(De+Ye))))
	p  <- p + stat_contour(aes(colour = ..level..))#,breaks=c(10,25,75,90))
	p  <- p + facet_grid(dm~bycatch)
	p  <- p + labs(x="Fishing intensity",y="Minumum Size Limit (cm)",color="Keepers\n Per\n 100 Caught")
	p3 <- p

	p  <- ggplot(S1,aes(x=fe,y=slim,z=SPR))
	p  <- p + stat_contour(aes(colour = ..level..),breaks=c(0.1,0.2,0.3,0.4))
	p  <- p + stat_contour(aes(colour=..level..), colour = "red", breaks=c(0.35))
	p  <- p + facet_grid(dm~bycatch)
	p  <- p + labs(x="Fishing intensity",y="Minumum Size Limit (cm)",col="Spawning\nPotential\nRatio")
	p4 <- p

	p  <- ggplot(S1,aes(x=fe,y=slim,z=YPR))
	p  <- p + stat_contour(aes(colour = ..level..))
	# p  <- p + stat_contour(aes(colour=..level..), colour = "red", breaks=c(0.35))
	p  <- p + facet_grid(dm~bycatch)
	p  <- p + labs(x="Fishing intensity",y="Minumum Size Limit (cm)",col="Yield\nPer\nRecruit")
	p5 <- p

	p  <- ggplot(S1,aes(x=fe,y=slim,z=Cbar))
	p  <- p + stat_contour(aes(colour = ..level..))
	# p  <- p + stat_contour(aes(colour=..level..), colour = "red", breaks=c(0.35))
	p  <- p + facet_grid(dm~bycatch)
	p  <- p + labs(x="Fishing intensity",y="Minumum Size Limit (cm)",col="Average\nWeight\n(lbs)")
	p6 <- p






