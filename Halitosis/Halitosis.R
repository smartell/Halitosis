# -------------------------------------------------------------------------- ##
# Halitosis.R
# Written by Steve Martell,  IPHC
# Date: Jan 8, 2012
# 
# Feb 7, 2011.  Added multiple growth groups to allow for cumulative effects
# of size selective fishing.  Need to show how mean weigth-at-age decreases
# with increasing fishing mortality rates.
# 
# SCENARIOS FOR HALIBUT PAPER
# 1) Status quo. 83.3cm minimum size limit
# 2) No size limits
# 3) 60cm minimum size limit
# 4) slot limit 81.3-150 cm
# 5) slot limit 60-150 cm
# -------------------------------------------------------------------------- ##

# -------------------------------------------------------------------------- ##
# Libraries
require(Hmisc)
require(ggplot2)
require(reshape2)
# -------------------------------------------------------------------------- ##

# Data and other constants
A	<- 35	# maximum age.
G	<- 21	# number of growth groups
S	<- 2	# number of sexes
dim	<- c(A, G, S)
age	<- 1:A	# vector of ages
pg	<- dnorm(seq(-3, 3, length=G), 0, 1); pg <- pg/sum(pg)

# Population parameters 
bo		<- 100.0			# unfished female spawning biomass
h		<- 0.75				# steepness
dm		<- 0.16				# discard mortality rate
THETA <- data.frame(bo=bo, h=h, dm=dm)

# Selectivity parameters (cm)
lhat	<- 97.132
lhat	<- 107.132
ghat	<- 1/0.1667
slim	<- 81.28
ulim	<- 150
cvlm	<- 0.1




# Sex specific parameters (female, male).
m		<- c(0.15, 0.18)			# natural mortality rate
a50		<- rep(10.91, 2)			# age at 50% maturity
k50		<- rep(1.406, 2)			# std at 50% maturity
a		<- rep(6.821e-6, 2)			# length-weight allometry (Clark 1992)
b		<- rep(3.24, 2)				# length-weight allometry (CLark 1992)
linf	<- c(145, 110)				# Range female 145-190, male 110-155 (cm)
k		<- c(0.1, 0.12)				# eyeballed growth pars from Clark & Hare 2002.
CVlinf  <- c(0.1, 0.1)				# CV in the asymptotic length

PHI	<- data.frame(m=m, a50=a50, k50=k50, a=a, b=b, linf=linf, k=k, CVlinf=CVlinf)
rownames(PHI)=c("Female", "Male")


# Halibut prices (10-20) (20-40) (40+)
# $6.75  $7.30  $7.50  In Homer Alaska.
fe		<- seq(0, 0.5, b=0.01)

lifeh	<-
function(fe = 0)
{
	# Life history and age-schedule information
	lx	<- array(0, dim)
	la	<- array(0, dim)
	wa	<- array(0, dim)
	fa	<- array(0, dim)
	pa	<- array(0, dim)	# price per pound 
	ma	<- plogis(age, a50, k50)
	for(i in 1:S)
	{
		lx[,,i]  <- exp(-m[i])^(age-1)
		lx[A,,i] <- lx[A,,i]/(1-exp(-m[i]))
		
		# growth
		'vonb'  <- function(linf,k) len <- linf*(1-exp(-k*age))
		dev     <- linf[i]*CVlinf[i]
		linf.g  <- seq(linf[i]-dev, linf[i]+dev, length=G)
		la[,,i] <- sapply(linf.g, vonb,k=k[i])
		wa[,,i] <- a[i]*la[,,i]^b[i]
		
		# maturity (this assumes maturity at a fixed age)
		fa[,,i] <- ma*wa[,,i]
	}
	
	# Length-based selectivity (length-based -> age-based)
	sc	<- array(0, dim)
	sr	<- array(0, dim)
	sd	<- array(0, dim)
	va	<- array(0, dim)
	std	<- cvlm*slim+1.e-30
	for(i in 1:S)
	{
		sc[,,i]  <- plogis(la[,,i],location=lhat, scale=ghat)
		sr[,,i]  <- plogis(la[,,i],location=slim, scale=std) - plogis(la[,,i],location=ulim, scale=std)
		sd[,,i]  <- 1-sr[,,i]
		va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm)
	}
	
	# Age-specific total mortality, survival, retention, and discard rate.
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
	# za	<- m+fe*va
	# sa	<- exp(-za)
	# qa	<- (sc*sr)*(1.-sa)/za	# fraction retained
	# da	<- (sc*sd)*(1.-sa)/za	# fraction discarded
	
	# Survivorship under fished conditions lz(A, G, S)
	lz	<- array(1, dim)
	for(i in 1:S)
	{
		for(j in 2:A)
		{
			lz[j,,i] <- lz[j-1,,i]*exp(-za[j-1,,i])
		}
		lz[A,,i] <- lz[A,,i]/(1-exp(-za[A,,i]))
	}
	
	
	# price premiums based on fish weight
	pa[wa<10]  <- 0.00
	pa[wa>=10] <- 6.75
	pa[wa>=20] <- 7.30
	pa[wa>=40] <- 7.50
	
	return(list(la=la, wa=wa, fa=fa, lz=lz, pa=pa*wa, va=va))
}

plot.AgeSchedule <- function()
{
	U   <- lifeh(0.2)
	mU  <- melt(U)
	colnames(mU) <- c("Age","G","Sex","value","type")
	mU$pg  <- pg[mU$G]
	mU$Sex[mU$Sex==1] = "Female"
	mU$Sex[mU$Sex==2] = "Male"
	
	p   <- ggplot(mU)
	p   <- p + geom_point(aes(x=Age, y=value, col=Sex, alpha=pg))
	p   <- p + facet_wrap(~type, scale="free")
	
	print(p)
	dev.copy2pdf(file="../FIGS/AgeSchedule.pdf")
	
}





tsasm	<- 
function(fe=0, slim=0, ulim=1000, dm=0.16)
{
	# A two sex age structured model. #
	
	
	# Age-schedule information
	# Survivorship array lx(A, G, S)
	lx	<- array(0, dim)
	la	<- array(0, dim)
	wa	<- array(0, dim)
	fa	<- array(0, dim)
	pa	<- array(0, dim)	# price per pound 
	ma	<- plogis(age, a50, k50)
	for(i in 1:S)
	{
		lx[,,i]  <- exp(-m[i])^(age-1)
		lx[A,,i] <- lx[A,,i]/(1-exp(-m[i]))
		
		# growth
		'vonb'  <- function(linf,k) len <- linf*(1-exp(-k*age))
		dev     <- linf[i]*CVlinf[i]
		linf.g  <- seq(linf[i]-dev, linf[i]+dev, length=G)
		la[,,i] <- sapply(linf.g, vonb,k=k[i])
		wa[,,i] <- a[i]*la[,,i]^b[i]
		
		# maturity (this assumes maturity at a fixed age)
		fa[,,i] <- ma*wa[,,i]
	}
	
	# price premiums based on fish weight
	pa[wa<10]  <- 3.00
	pa[wa>=10] <- 6.75
	pa[wa>=20] <- 7.30
	pa[wa>=40] <- 7.50
	
	# lx	<- sapply(age, function(age) exp(-m)^(age-1) )
	# lx[, A] <- lx[, A]/(1-exp(-m))
	# la	<- sapply(age, function(age) linf*(1-exp(-k*age)))
	# wa	<- a*la^b
	# ma	<- plogis(age, a50, k50)
	# fa	<- t(t(wa)*ma)
	
	# Derivation of S-R parameters (Ricker model)
	# Re = Ro*[(log(kap)-log(phi.E/phi.e))/log(kap)]
	# Assume 50:50 sex ratio at birth
	kap		<- (5*h)^(5/4)
	phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg))
	#phi.E	<- as.double(lx[1, ] %*% fa[1, ])
	ro		<- bo/phi.E
	
	# Length-based selectivity (length-based -> age-based)
	sc	<- array(0, dim)
	sr	<- array(0, dim)
	sd	<- array(0, dim)
	va	<- array(0, dim)
	std	<- cvlm*slim+1.e-30
	for(i in 1:S)
	{
		sc[,,i]  <- plogis(la[,,i],location=lhat, scale=ghat)
		sr[,,i]  <- plogis(la[,,i],location=slim, scale=std) - plogis(la[,,i],location=ulim, scale=std)
		sd[,,i]  <- 1-sr[,,i]
		va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm)
	}
	# sc	<- t(apply(la,1,plogis,location=lhat,scale=ghat))
	# std	<- cvlm*slim+1.e-30
	# sr	<- t(apply(la,1,plogis,location=slim,scale=std))
	# sd	<- 1-sr
	# va	<- sc*(sr+sd*dm)		# age-specific probability of dying due to F
	
	# Age-specific total mortality, survival, retention, and discard rate.
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
	# za	<- m+fe*va
	# sa	<- exp(-za)
	# qa	<- (sc*sr)*(1.-sa)/za	# fraction retained
	# da	<- (sc*sd)*(1.-sa)/za	# fraction discarded
	
	# Survivorship under fished conditions lz(A, G, S)
	lz	<- array(1, dim)
	for(i in 1:S)
	{
		for(j in 2:A)
		{
			lz[j,,i] <- lz[j-1,,i]*exp(-za[j-1,,i])
		}
		lz[A,,i] <- lz[A,,i]/(1-exp(-za[A,,i]))
	}
	# lz	<- matrix(1, 2, A)
	# for(i in 2:A)
	# {
	# 	lz[, i] <- lz[, i-1]*exp(-za[, i-1])
	# 	if(i==A)
	# 	{
	# 		lz[, A] <- lz[, A]/(1-exp(-za[, i-1]))
	# 	}
	# }
	
	# Incidence functions
	phi.e	<- sum( t(lz[,,1]*fa[,,1])*pg )
	#phi.e	<- as.double(lz[1, ] %*% fa[1, ])
	
	# Equilibrium calculations
	t1		<- log(phi.E/(kap*phi.e))
	t2		<- (log(kap)*phi.e)
	re		<- max(0, -(t1*ro*phi.E)/t2)
	be		<- re * phi.e
	ye		<- 0
	de		<- 0
	ypr		<- 0
	wbar	<- rep(0, S)
	wdot	<- matrix(0,nrow=S, ncol=A)
	for(i in 1:S)
	{
		ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
		de	<- de + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
		ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
		
		# Average weigth of a 10-year old fish (female & male)
		tmp		<- t(lz[,,i]*wa[,,i])*pg
		tmpn	<- t(lz[,,i])*pg
		wbar[i] <- weighted.mean(wa[10,,i], tmpn[,10])
		
		# Average weight-at-age for each sex.
		tmp     <- t(tmpn)
		P       <- tmp/rowSums(tmp)
		wdot[i,] <- rowSums(wa[,,i]*P)
	}
	spr		<- phi.e/phi.E
	
	
	# Need to calculate landed value per recruit versus slim and fe
	# using the price and size categories from above.
	# How many millions of dollars in Halibut are discarded each year?
	landed.value	<- 0
	discard.value	<- 0
	for(i in 1:S)
	{
		t1 <- sum(re * fe * t(lz[,,i]*wa[,,i]*qa[,,i]*pa[,,i])*pg)
		landed.value  <- landed.value + t1
		t2 <- sum(re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i]*pa[,,i])*pg)
		discard.value <- discard.value + t2
	}
	
	# b<- seq(0, 2*bo, length=100)
	# r<- kap*b*exp(-log(kap)*b/bo)/phi.E
	# plot(b, r)
	# points(bo, ro, pch=20, col=2)
	# points(be, re, pch=20, col=3)
	return(c(fe=fe, re=re, be=be, ye=ye, 
		de=de, spr=spr, ypr=ypr, 
		dep=be/bo, wbar.f10=wbar[1], wbar.m10=wbar[2], 
		landed.value=landed.value, 
		discard.value=discard.value ))
		# List objects don't work here.
		#, wdot.f=list(wdot[1,]), wdot.m=list(wdot[2,]) ))
}

.equil	<-
function(arg="ye", dm=0.17)
{
	fn	<- function(fe=0, slim=0, dm=dm)
	{
		tmp <- tsasm(fe, slim, dm)
		idx <- which(names(tmp)==arg)
		return(as.double(tmp[idx]))
	}
	V	<- Vectorize(fn, c("fe", "slim"))
	fe	<- seq(0, 0.45, length=20)
	sl	<- seq(60, 100, length=20)
	Z	<- outer(fe, sl, V, dm)
	obj	<- list(x=fe, y=sl, Z=Z)
	class(obj) <- "isopleth"
	return(obj)
}

plot.isopleth <- 
function(obj, ...)
{
	# Plots the contour plot for the contour class
	x <- obj$x
	y <- obj$y
	z <- obj$Z
	
	#contour(x, y, z, ...)
	mx <- melt(z)
	names(mx) <- c("x","y","z")
	g=ggplot(mx,aes(x=x,y=y,z=z))+stat_contour()
	print(g)
}

# SCENARIOS
lambda = 0.16
S1 = data.frame(Scenario="Status quo", t(sapply(fe, tsasm, slim=81.3, dm=lambda)))
S2 = data.frame(Scenario="No size limit", t(sapply(fe, tsasm, slim=0.00, dm=lambda)))
S3 = data.frame(Scenario="Mininum SL = 70cm", t(sapply(fe, tsasm, slim=70.0, dm=lambda)))
S4 = data.frame(Scenario="Slot size (81.3-150)", t(sapply(fe, tsasm, slim=81.3, dm=lambda, ulim=150)))
S5 = data.frame(Scenario="Slot size (70.0-150)", t(sapply(fe, tsasm, slim=70.0, dm=lambda, ulim=150)))

DF = rbind(S1, S2, S3, S4, S5)

p<-ggplot(DF,aes(x=fe,y=de,col=Scenario)) +geom_line()
print(p+opts(title="Wastage"))
dev.copy2pdf(file="../FIGS/fig:wastage.pdf", width=10, height=7.5)

p<-ggplot(DF,aes(x=fe,y=spr,col=Scenario)) +geom_line()
print(p+opts(title="Relative spawning biomass per recruit"))

p<-ggplot(DF,aes(x=fe,y=ypr,col=Scenario)) +geom_line()
print(p+opts(title="Yield per recruit"))

p<-ggplot(DF,aes(x=fe,y=ye,col=Scenario)) +geom_line()
print(p+opts(title="Equilibrium yield"))
dev.copy2pdf(file="../FIGS/fig:yield.pdf", width=10, height=7.5)

p<-ggplot(DF,aes(x=fe,y=landed.value,col=Scenario)) +geom_line()
print(p+opts(title="Landed Value"))

p<-ggplot(DF,aes(x=fe,y=discard.value,col=Scenario)) +geom_line()
print(p+opts(title="Value of dead discards"))


# A key question is for every pound of bycatch what is the corresponding
# yield loss to the directed fishery. This is computed by IPHC as the 
# yield loss ration = (Wt. of future yield loss)/(Wt. of bycatch).  The wt. of 
# the bycatch is straight forward. The yield loss is the difference between
# the yield obtained with discard mortality =0 and discard mortality =0.17

oldfun <-function()
{
	SPR <- .equil("spr", dm=dm)
	SPR0<- .equil("spr", dm=0)
	YE  <- .equil("ye", dm=dm)
	YE0 <- .equil("ye", dm=0)
	BE  <- .equil("be", dm=dm)
	BE0 <- .equil("be", dm=0)
	DE	<- .equil("de", dm=dm)
	W.F <- .equil("wbar.f", dm=dm)
	W.M <- .equil("wbar.m", dm=dm)
	LV	<- .equil("landed.value", dm=dm)
	DV	<- .equil("discard.value", dm=dm)

	# REPORT SECTION
	par(mfcol=c(1, 1), las=1)
	isolvl <- c(0.35, seq(0, 1, by=0.1))
	isolwd <- c(2, rep(1, 11))
	xl     <- "Fishing mortality"
	yl     <- "Size limit (cm)"
	plot(SPR,xlab=xl,ylab=yl,levels=isolvl,lwd=isolwd,main="Spawn potential ratio")
	plot(YE ,xlab=xl,ylab=yl,main="Equilibrium yield")
	plot(DE ,xlab=xl,ylab=yl,main="Discarded yield")
	X = DE
	X$Z = (YE0$Z-YE$Z)/(DE$Z)
	plot(X, xlab=xl,ylab=yl,main="Yield loss ratio")

	#The following in the spawning biomass per recruit lost per 
	#unit of discard. This should be the spawning biomass,  not SPR
	SE = DE
	SE$Z = (BE0$Z-BE$Z)/(DE$Z)
	plot(SE, add=TRUE, col="blue", levels=seq(0, 10, by=.25))

	E=DE
	E$Z = YE$Z/(YE$Z+DE$Z)
	plot(E, add=TRUE, col="red")

	par(mfcol=c(2, 2))
	plot(W.F, xlab=xl, ylab=yl, main="Mean weight of age-10 females")
	plot(W.M, xlab=xl, ylab=yl, main="Mean weight of age-10 males")
	plot(LV, xlab=xl, ylab=yl, main="Landed Value ($$)")
	X = LV
	X$Z = DV$Z/LV$Z
	plot(X, xlab=xl, ylab=yl, main="Discard Value/Landed Value")
}