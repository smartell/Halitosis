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
# 3) 70cm minimum size limit
# 4) slot limit 81.3-150 cm
# 5) slot limit 70-150 cm
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
# Commercial selectivities from Hare 2012.
bin   <- seq(50, 120, by=10)
CSelL <- matrix(data=c(1.63180e-09,  3.25740e-09,  6.03022e-02,  3.00891e-01,  6.30344e-01,  9.13893e-01, 1.00000e+00,  1.00000e+00, 2.17159e-09,  3.96878e-03,  5.67109e-02,  2.81436e-01,  5.85461e-01,  8.35614e-01, 1.00000e+00,  1.00000e+00), nrow=8, ncol=2)

pg	<- dnorm(seq(-3, 3, length=G), 0, 1); pg <- pg/sum(pg)

# Population parameters 
bo		<- 100.0			# unfished female spawning biomass
h		<- 0.85				# steepness
dm		<- 0.16				# discard mortality rate
THETA <- data.frame(bo=bo, h=h, dm=dm)

# Selectivity parameters (cm)
lhat	<- 97.132
ghat	<- 1/0.1667
slim	<- 81.28
ulim	<- 1500
cvlm	<- 0.1
PI		<- data.frame(lhat=lhat, ghat=ghat, slim=slim, ulim=ulim, cvlm=cvlm, bin=bin)



# Sex specific parameters (female, male).
m		<- c(0.15, 0.135)			# natural mortality rate
#0.15 0.135474
a50		<- rep(10.91, 2)			# age at 50% maturity
k50		<- rep(1.406, 2)			# std at 50% maturity
a		<- rep(6.821e-6, 2)			# length-weight allometry (Clark 1992)
b		<- rep(3.24, 2)				# length-weight allometry (CLark 1992)
#linf	<- c(145, 110)				# Range female 145-190, male 110-155 (cm)
linf    <- c(151.568,  99.3607)		# From 2011 Length_age data
#k		<- c(0.1, 0.12)				# eyeballed growth pars from Clark & Hare 2002.
k		<- c(0.0820581, 0.135409)   # Fron 20111 Length_age data
CVlinf  <- c(0.1, 0.1)				# CV in the asymptotic length

PHI	<- data.frame(m=m, a50=a50, k50=k50, a=a, b=b, linf=linf, k=k, CVlinf=CVlinf)
rownames(PHI)=c("Female", "Male")

T1 <- c(THETA, PHI, PI)

# Halibut prices (10-20) (20-40) (40+)
# $6.75  $7.30  $7.50  In Homer Alaska.
fe		<- seq(0, 0.6, b=0.01)

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
		# SM Adding approx function,  using the length-based selectivity coefficients
		# CSelL from Hare's assessment model.
		#approx(bin,CSelL.F,xout=la,yleft=0,yright=1)
		sc[,,i]  <- approx(bin, CSelL[,i], la[,,i], yleft=0, yright=1)
		
		#sc[,,i]  <- plogis(la[,,i],location=lhat, scale=ghat)
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
#function(fe=0, slim=0, ulim=1000, dm=0.16)
function(fe=0, theta)
{
	# A two sex age structured model. #
	# SM Changing arguments (fe,  THETA)
	# Where theta is a list of the global parameters used in the model.
	
	with(as.list(theta), {
		
	
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
	pa[wa<5]   <- 0.00
	pa[wa>=5.] <- 3.00
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
		#sc[,,i]  <- plogis(la[,,i],location=lhat, scale=ghat)
		sc[,,i]  <- approx(bin, CSelL[,i], la[,,i], yleft=0, yright=1)$y
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
		de=de, spr=spr, ypr=ypr, eye=max(0, landed.value-8*fe), 
		dep=be/bo, wbar.f10=wbar[1], wbar.m10=wbar[2], 
		landed.value=landed.value, 
		discard.value=discard.value ))
		# List objects don't work here.
		#, wdot.f=list(wdot[1,]), wdot.m=list(wdot[2,]) ))
	})
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
.scenarios <- function(T1, Model="M1")
{
	T2 = T1; T2$slim=0.00
	T3 = T1; T3$slim=81.3; T3$h=1.0
	T4 = T1; T4$ulim=131.5
	T5 = T1; T5$ulim=131.5; T5$slim=70.0

	S1 = data.frame(Model=Model, Scenario="Status quo", t(sapply(fe, tsasm, theta=T1)))
	S2 = data.frame(Model=Model, Scenario="No size limit", t(sapply(fe, tsasm, theta=T2)))
	S3 = data.frame(Model=Model, Scenario="No SR relationship", t(sapply(fe, tsasm, theta=T3)))
	S4 = data.frame(Model=Model, Scenario="Slot size (81.3-131.5)", t(sapply(fe, tsasm, theta=T4)))
	S5 = data.frame(Model=Model, Scenario="Slot size (70.0-131.5)", t(sapply(fe, tsasm, theta=T5)))

	DF = rbind(S1, S2, S3, S4, S5)
	return(DF)
}

M1 <- T1
M2 <- T1; M2$bin=T1$bin - 10
M3 <- T1; M3$bin=T1$bin + 10
M4 <- T1; M4$linf=0.9*T1$linf
M5 <- T1; M5$linf=0.9*T1$linf; M5$bin=T1$bin-10
M6 <- T1; M5$linf=0.9*T1$linf; M6$bin=T1$bin+10
DF <- rbind(.scenarios(M1, "Base"), .scenarios(M2, "(-) selectivity"), .scenarios(M3, "(+) selectivity"), 
	.scenarios(M4, "(-) Growth"), .scenarios(M5, "(-) growth & (-) selectivity"), .scenarios(M6, "(-) growth & (+) selectivity"))

# Max yield from status quo
imax    <- which.max(subset(subset(DF,Model=="Base"),Scenario=="Status quo")$ye)
ye_max  <- DF$ye[imax]
ypr_max <- DF$ypr[imax]
spr_max <- DF$spr[imax]
eye_max <- DF$eye[imax]
de_max  <- DF$de[imax]
lv_max  <- DF$landed.value[imax]

p1 <- ggplot(DF,aes(x=fe,y=ypr/ypr_max,col=Scenario)) +geom_line()
p1 <- p1 + labs(x = "Fishing mortality", y = "Relative yield per recruit") 
p1 <- p1 + facet_wrap(~Model)

p2 <- ggplot(DF,aes(x=fe,y=spr/spr_max,col=Scenario)) +geom_line()
p2 <- p2 + labs(x = "Fishing mortality", y = "Relative spawning potential ratio") 
p2 <- p2 + facet_wrap(~Model)

p3 <- ggplot(DF,aes(x=fe,y=ye/ye_max, col=Scenario)) +geom_line()
p3 <- p3 + labs(x = "Fishing mortality", y = "Relative yield") 
p3 <- p3 + facet_wrap(~Model)

p4 <-ggplot(DF,aes(x=fe,y=eye/eye_max,col=Scenario)) +geom_line()
p4 <- p4 + labs(x = "Fishing mortality", y = "Relative net landed value ($)") 
p4 <- p4 + facet_wrap(~Model)


p5 <- ggplot(DF,aes(x=fe,y=de/de_max,col=Scenario)) + geom_line() 
p5 <- p5 + labs(x = "Fishing mortality", y = "Wastage (lb)") 
p5 <- p5 + facet_wrap(~Model)

p6 <- ggplot(DF,aes(x=fe,y=landed.value/lv_max,col=Scenario)) +geom_line()
p6 <- p6 + labs(x = "Fishing mortality", y = "Landed Value ($)") 
p6 <- p6 + facet_wrap(~Model)







p<-ggplot(DF,aes(x=fe,y=discard.value,col=Scenario)) +geom_line()+ facet_wrap(~Model)
print(p+opts(title="Value of dead discards"))

draw <- function()
{
	print(p1); ggsave(file="../FIGS/fig:ypr.pdf", width=10, height=7.5)
	print(p2); ggsave(file="../FIGS/fig:spr.pdf", width=10, height=7.5)
	print(p3); ggsave(file="../FIGS/fig:yield.pdf", width=10, height=7.5)
	print(p4); ggsave(file="../FIGS/fig:eye.pdf", width=10, height=7.5)
	print(p5); ggsave(file="../FIGS/fig:waste.pdf", width=10, height=7.5)
}

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