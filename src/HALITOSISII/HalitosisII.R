

####################################
### Size-selective fishing model ###
### a.k.a "HALITOSIS II"         ###
### written by Steve Martell     ###
####################################
.RES = 300 #png file resolution

# |---------------------------------------------------------------------------|
# | Dependencies
# |---------------------------------------------------------------------------| 
if(!require("Hmisc"))   install.packages("Hmisc")  ## latex
if(!require("ggplot2"))   install.packages("ggplot2")
if(!require("reshape2"))   install.packages("reshape2")
if(!require("grid"))   install.packages("grid") ## arrow function
if(!require("plyr"))   install.packages("plyr")
if(!require("extrafont"))   install.packages("extrafont")

source("src/HALITOSISII/Selex.R")

mytheme <- theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(),
        panel.border = element_blank(),
        axis.text=element_text(family="serif",size=12),
        axis.title=element_text(family="serif",size=12),
        legend.text=element_text(size=12,family="serif"),
        legend.title=element_text(size=12,family="serif"),
        panel.background = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.position="top")
# |---------------------------------------------------------------------------|
# | Data and other constants
# |---------------------------------------------------------------------------|
# | Stock -> is a list object with all parameters and output.
regArea <- "2B"
A	<- 30					# maximum age.
G	<- 11					# number of growth groups
S	<- 2					# number of sexes
dim	<- c(A, G, S)			# array dimensions
age	<- 1:A					# vector of ages
pg	<- dnorm(seq(-1.96, 1.96, length=G), 0, 1); #probability of being a group or the
pg  <- pg/sum(pg) 			# proportion assigned to each growth-type group.

Stock <- list(A=A,G=G,S=S,dim=dim,age=age,pg=pg)


# |---------------------------------------------------------------------------|
# | Population parameters 
# |---------------------------------------------------------------------------|
bo		<- 476.891			# unfished female spawning biomass from Stewart 2012
#    r0   <- 1
#    phiE <-
#    b0   <- r0/phiE
h		<- 0.75				# steepness
m        <- c(0.201826,0.169674) # natural mortality rate
linf     <- c(313,172)		# asymptotic length
vonk     <- c(0.042,0.066)	# vonk
to       <- c(-2.6e-8,-0.92)   # time at zero length
p        <- c(1,1)  # vonb Power parameter.
cv       <- c(0.14,0.14)		# CV in length-at-age

a50      <- c(12.8,10.6)     # age at 50% maturity. ## 11.6 yrs (Stewart 2014) 
k50      <- c(0.72,0.53)     # the scale or slope parameter 50% maturity. 
a				<- rep(3.139e-6, 2) # length-weight allometry (Clark 1992)
b				<- rep(3.24, 2)		# length-weight allometry (CLark 1992)

# dm		<- 0.16				# discard mortality rate
cm		<- 0				# Size-dependent natural mortality rate (-0.5, 0.5)
Stock <- c(Stock,list(bo=bo,h=h,m=m,linf=linf,vonk=vonk,to=to,p=p,cv=cv,a50=a50,k50=k50))
Stock <- c(Stock,list(a=a,b=b,cm=cm))


# |---------------------------------------------------------------------------|
# | Commercial selectivities from Stewart 2012 (old from Hare's model).
# | This has been updated in the modern stock assessments.
# |---------------------------------------------------------------------------|
bin   <- seq(60, 130, by=10) #60 to 130 cm by 10 cm intervals
CSelL <- matrix(
  (data=c(0,  0.0252252,  0.250685,  0.617268,  1,  1.36809,  1.74292,  2.12022, 
          0,  0.0151914,  0.144236,  0.513552,  1,  1.48663,  1.97208,  2.45702)
  ), nrow=8, ncol=2) #in order from 60 to 130 cm by 10 cm intervals, females first then males


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

fe  <- seq(0, 1.00, by=0.01)#sequence of fishing mortality rates.
sizelimit <- 81.28 # 2.54 * c(seq(20,65,by=2))
discardmortality <- 0.16 #c(0.10, 0.15, 0.20, 0.25)
bycatch  <- 0 #c(0,10,20)

bycatchSel <- c(0, 0, 0.379083, 0.923116, 1, 0.748264, rep(0.650509,length=29))
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
  
  with(Stock, {  #Stock is a function
    lx	   <- array(1, dim)
    la	   <- array(0, dim)
    sd_la  <- array(0, dim)
    wa	   <- array(0, dim)
    fa	   <- array(0, dim)
    pa     <- array(0, dim)
    M      <- array(0, dim)
    ma     <- plogis(age, a50, k50) #maturity is a function of age, not length
    
    for(i in 1:S) # S = sexes
    {
      # Lenght-at-age
      mu    <- gvonb(age,linf[i],vonk[i],to[i],p[i])
      sigma <- cv[i] * mu #standard deviation
      dev   <- seq(-1.96, 1.96, length=G) #deviance, standard normal deviation
      if(G==1) dev <- 0
      
      la[,,i]    <- sapply(dev,fn<-function(dev){la=mu+dev*sigma}) #la = mean length-at-age
      sd_la[,,i] <- sqrt(1/G*(cv[i]*mu)^2)
      wa[,,i]    <- a[i]*la[,,i]^b[i] #weight at age
      fa[,,i]    <- ma*wa[,,i] #fecundity at age
      
      # Size dependent natural mortality rate 
      # M_l = M (l_a/l_r)^c 
      l_r     <- 100 # reference length
      delta   <- (la[,,i]/l_r)^cm / mean((la[,,i]/l_r)^cm)
      M[,,i]  <- m[i] * delta
      
      # Survivorship
      for(j in 2:A)
      {
        lx[j,,i] <- lx[j-1,,i]*exp(-M[j-1,,i]) #unfished
      }
      lx[A,,i] <- lx[A,,i]/(1-exp(-M[A,,i])) #plus groups
    }
    
    # price premiums based on fish weight
    pa[wa<10]  <- 0.00
    pa[wa>=10] <- 6.75
    pa[wa>=20] <- 7.30
    pa[wa>=40] <- 7.50
    
    Stock$lx	   = lx	#unfished surv.
    Stock$la	   = la	
    Stock$sd_la    = sd_la
    Stock$wa	   = wa
    Stock$pa       = pa
    Stock$fa	   = fa	
    Stock$M        = M #natural mortality (independent of age and size)
    Stock$ma       = ma #maturity at age
    
    return(Stock)
  })
  
}


# |---------------------------------------------------------------------------|
# | Calculate size-based selectivities and joint capture probability
# |---------------------------------------------------------------------------|
# | This function calls .calcPage(la, sa, pl, xl) in Selex.R
.calcSelectivities <- function(Stock,slim=0,ulim=1500,cvlm=0.1,dm=0) #slim=lowersizelimit
  #ulim=upper size limit, cvlm= variability within saa within 
  #gtg plus measurement error or measuring indv fish, dm=discard mort
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
    vd	<- array(0, dim)  #bycatch fishery
    std	<- cvlm*slim+1.e-30	
    for(i in 1:S)
    {
      #probt of catching a fish of a given length:
      pl       <- approx(bin, CSelL[,i], xl, yright=1, yleft=0)$y  # selectivity curve, length-based 
      sc[,,i]  <- .calcPage(la[,,i],sd_la[,,i],pl,xl) #sc = size of capture # selex.R
      #sc[,,i]  <- approx(bin, CSelL[,i], la[,,i], yleft=0, yright=1)$y
      
      # retention proability 
      pr       <-  plogis(xl, slim, 0.1) - plogis(xl, ulim, 0.1) 
      sr[,,i]  <- .calcPage(la[,,i],sd_la[,,i],pr,xl) #calculates the proportion of indiv at age that are vulnerable
      
      #probabity of discarding it..
      #sr[,,i]  <-  plogis(la[,,i],location=slim, scale=std) - plogis(la[,,i],location=ulim, scale=std)
      sd[,,i]  <- 1-sr[,,i]
      va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm) 
      
      # bycatch fishery selecitvity... made up
      pd       <- plogis(xl, 40,  5) - (1-0.65)*plogis(xl, 100.3, 5.0) #subtracting these makes it dome-shaped.
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
    # Unfished eggs per recruit or spawning stock per recruit (phi.E)
    phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg)) #survivorship*fecundity*proportion of recruits that end up in each growth type group
    #"1" - indicates females, only interested in female spawning stock biomass
    # t() is transpose
    # Unfished recruitment (ro)
    ro		<- bo/phi.E
    
    # Beverton-Holt model
    kap <- 4*h/(1-h) #reck, phil goodyears compensation ratio, h=steepness (fraction of unfished recruit when you reduce spawning stock biomass by 20% of its unfished state)
    
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
    za	<- array(0, dim) #Age-specific total mortality
    sa	<- array(0, dim) # age-spec survival
    qa	<- array(0, dim) #age specific yield in weight ??
    da	<- array(0, dim) #age-specific discard
    ta  <- array(0, dim) # yield per recruit in trawl discard fishery.
    
    bycatch <- ct
    # bapprox <- bo * 0.15/(0.15+fe)
    # fd      <- bycatch/bapprox
    fd      <- 0;
    #if(ct>0)
    #cat("fe = ", fe, " fd = ", fd, "\n")
    
    ### 
    for(iter in 1:25)
    {
      
      
      for(i in 1:S)
      {
        za[,,i]  <- M[,,i] + fe*va[,,i] + fd*vd[,,i]
        sa[,,i]  <- exp(-za[,,i])
        qa[,,i]  <- (sc[,,i]*sr[,,i]) * (1-sa[,,i])/za[,,i]
        da[,,i]  <- (sc[,,i]*sd[,,i]) * (1-sa[,,i])/za[,,i]
        ta[,,i]  <- (vd[,,i])         * (1-sa[,,i])/za[,,i]
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
      if(re ==0 ) break();
      
      # 3b. Calculate bycatch per recruit to determine F in bycatch fisheries.
      de      <- 0          
      for(i in 1:S)
      {
        de	<- de + re * sum( t(lz[,,i]*wa[,,i]*ta[,,i])*pg )  # ta = yield per recruit in trawl discard fishery.
      }
      if(bycatch < de)
        fd <- -log(1.0-bycatch/de)
      else
        fd <- -log(0.01) 
      # cat("fd = ",fd,"\t de = ",de,"re = ",re," \n")
    }
    # 4. Calculate yield per recruit, spawning biomass per recruit,  yield, discards.
    ye		<- 0 #eqb yield in dir fishery
    ne    <- 0   #eqb numbers
    de		<- 0 #eqb discards
    we    <- 0   #eqb wastage
    yev   <- 0   #eqb value of yield in dir fishery
    dev   <- 0   # Value of wastage.
    byv   <- 0 	 # Value of discarded bycatch.
    ypr		<- 0 #yield per recruit in dir fishery
    # 		lpr   <- 	#length per recruit
    dpr   <- 0   #discarded fish in dir fishery (live and dead)
    wpr   <- 0   #wastage per recruit in dir fishery
    bpr   <- 0   #biomass per recruit available to the dir fishery (Ebio)
    spr		<- phi.e/phi.E #spawning potential ratio (spr)
    for(i in 1:S)
    {
      
      ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
      ne  <- ne + sum( re * fe * t(lz[,,i]*qa[,,i])*pg )
      we	<- we + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
      de	<- de + sum( re * fe * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
      bpr <- bpr + sum( t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
      ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
      # 			lpr <- lpr + sum( fe * t(lz[,,i]*la[,,i]*qa[,,i])*pg )
      dpr <- dpr + sum( fe * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
      wpr <- wpr + sum( fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
      n18 <- re * fe * t(lz[18,,i]*qa[,,i])*pg
      
      #
      # landed value of retained fish.
      #
      yev <- yev + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i]*pa[,,i])*pg )
      # 
      # value of discarded fish.
      # 
      dev <- dev + sum( re * fe * t(lz[,,i]*wa[,,i]*da[,,i]*pa[,,i])*pg )
      #
      # value of bycaught fish
      #
      byv <- byv + sum( re * fd * t(lz[,,i]*wa[,,i]*ta[,,i]*pa[,,i])*pg )
    }
    cbar <- ye / ne  #average weight of catch
    
    # 5. Calculate average weight-at-age
    wbar <- matrix(0, nrow=S, ncol=A)
    for(i in 1:S)
    {
      tmp      <- lz[,,i]/rowSums(lz[,,i]) #average weighted by survivorship, transforming surv to a probability
      wbar[i,] <- rowSums(wa[,,i]*tmp)
    }
    
    
    # 6. Calculate average weight of the landed catch.
    
    Stock$lz   <- lz
    Stock$re   <- re
    Stock$be   <- be
    Stock$ye   <- ye
    Stock$de   <- de
    Stock$we   <- we
    Stock$fd   <- fd
    Stock$bpr  <- bpr
    Stock$ypr  <- ypr
    Stock$spr  <- spr
    Stock$dpr  <- dpr
    Stock$wpr  <- wpr
    Stock$wbar <- wbar
    Stock$cbar <- cbar
    Stock$yev  <- yev 
    Stock$dev  <- dev
    Stock$byv  <- byv
    
    
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
      out <- c(fe=fe, ye=tmp$ye, be=tmp$be, de=tmp$de, we=tmp$we,
               re=tmp$re, spr=tmp$spr, ypr=tmp$ypr, fd = tmp$fd,
               bpr=tmp$bpr, dpr=tmp$dpr,wpr=tmp$wpr)
      
      # average weight arrays
      wbar_f1 <- c(wbar=tmp$wbar[1,])
      wbar_f2 <- c(wbar=tmp$wbar[1,])
      wbar_f3 <- c(wbar=tmp$wbar[1,])
      wbar_f4 <- c(wbar=tmp$wbar[1,])
      wbar_f5 <- c(wbar=tmp$wbar[1,])
      wbar_f6 <- c(wbar=tmp$wbar[1,])
      wbar_f7 <- c(wbar=tmp$wbar[1,])
      wbar_f8 <- c(wbar=tmp$wbar[1,])
      wbar_f9 <- c(wbar=tmp$wbar[1,])
      wbar_f10 <- c(wbar=tmp$wbar[1,])
      wbar_f11 <- c(wbar=tmp$wbar[1,])
      wbar_f12 <- c(wbar=tmp$wbar[1,])
      wbar_f13 <- c(wbar=tmp$wbar[1,])
      wbar_f14 <- c(wbar=tmp$wbar[1,])
      wbar_f15 <- c(wbar=tmp$wbar[1,])
      wbar_f16 <- c(wbar=tmp$wbar[1,])
      wbar_f17 <- c(wbar=tmp$wbar[1,])
      wbar_f18 <- c(wbar=tmp$wbar[1,])
      wbar_f19 <- c(wbar=tmp$wbar[1,])
      wbar_f20 <- c(wbar=tmp$wbar[1,])
      
      wbar_m1 <- c(wbar=tmp$wbar[2,])
      wbar_m2 <- c(wbar=tmp$wbar[2,])
      wbar_m3 <- c(wbar=tmp$wbar[2,])
      wbar_m4 <- c(wbar=tmp$wbar[2,])
      wbar_m5 <- c(wbar=tmp$wbar[2,])
      wbar_m6 <- c(wbar=tmp$wbar[2,])
      wbar_m7 <- c(wbar=tmp$wbar[2,])
      wbar_m8 <- c(wbar=tmp$wbar[2,])
      wbar_m9 <- c(wbar=tmp$wbar[2,])
      wbar_m10 <- c(wbar=tmp$wbar[2,])
      wbar_m11 <- c(wbar=tmp$wbar[2,])
      wbar_m12 <- c(wbar=tmp$wbar[2,])
      wbar_m13 <- c(wbar=tmp$wbar[2,])
      wbar_m14 <- c(wbar=tmp$wbar[2,])
      wbar_m15 <- c(wbar=tmp$wbar[2,])
      wbar_m16 <- c(wbar=tmp$wbar[2,])
      wbar_m17 <- c(wbar=tmp$wbar[2,])
      wbar_m18 <- c(wbar=tmp$wbar[2,])
      wbar_m19 <- c(wbar=tmp$wbar[2,])
      wbar_m20 <- c(wbar=tmp$wbar[2,])
      
      
      
      out <- c(out, wbar_f1=wbar_f1,wbar_f2=wbar_f2,wbar_f3=wbar_f3,wbar_f4=wbar_f4, wbar_f5=wbar_f5,
               wbar_f6=wbar_f6,wbar_f7=wbar_f7,wbar_f8=wbar_f8,wbar_f9=wbar_f9,wbar_f10=wbar_f10,
               wbar_f11=wbar_f11,wbar_f12=wbar_f12,wbar_f13=wbar_f13,wbar_f14=wbar_f14,wbar_f15=wbar_f15,
               wbar_f16=wbar_f16,wbar_f17=wbar_f17,wbar_f18=wbar_f18,wbar_f19=wbar_f19,wbar_f20=wbar_f20, 
               
               wbar_m1=wbar_m1,wbar_m2=wbar_m2,wbar_m3=wbar_m3,wbar_m4=wbar_m4, wbar_m5=wbar_m5,
               wbar_m6=wbar_m6,wbar_m7=wbar_m7,wbar_m8=wbar_m8,wbar_m9=wbar_m9,wbar_m10=wbar_m10,
               wbar_m11=wbar_m11,wbar_m12=wbar_m12,wbar_m13=wbar_m13,wbar_m14=wbar_m14,wbar_m15=wbar_m15,
               wbar_m16=wbar_m16,wbar_m17=wbar_m17,wbar_m18=wbar_m18,wbar_m19=wbar_m19,wbar_m20=wbar_m20, 
               
               Scenario=Scenario)
      
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
    M1 <- .calcLifeTable(Stock) # list of 
    M1 <- .calcSelectivities(M1,slim=slim,ulim=1500,cvlm=0.1,dm=dm)
    M1 <- .calcSRR(M1) #spawnerperrecruit relationship
    M1 <- .asem(fe,M1,bycatch) #age-structures eqB model
    
    out <- c(fe=fe,slim=slim,dm=dm,bycatch=bycatch,
             Ye=M1$ye,De=M1$de,We=M1$we,Be=M1$be,Re=M1$re,
             SPR=M1$spr,YPR=M1$ypr,DPR=M1$dpr,WPR=M1$wpr,BPR=M1$bpr,
             Cbar=M1$cbar,EE=M1$ye/(M1$ye+M1$de),Fd = M1$fd,
             YEv=M1$yev,DEv=M1$dev,BYv=M1$byv,
             
             wbar_f1=M1$wbar[1,1],wbar_m1=M1$wbar[2,1],
             wbar_f2=M1$wbar[1,2],wbar_m2=M1$wbar[2,2],
             wbar_f3=M1$wbar[1,3],wbar_m3=M1$wbar[2,3],
             wbar_f4=M1$wbar[1,4],wbar_m4=M1$wbar[2,4],
             wbar_f5=M1$wbar[1,5],wbar_m5=M1$wbar[2,5],
             wbar_f6=M1$wbar[1,6],wbar_m6=M1$wbar[2,6],
             wbar_f7=M1$wbar[1,7],wbar_m7=M1$wbar[2,7],
             wbar_f8=M1$wbar[1,8],wbar_m8=M1$wbar[2,8],
             wbar_f9=M1$wbar[1,9],wbar_m9=M1$wbar[2,9],
             wbar_f10=M1$wbar[1,10],wbar_m10=M1$wbar[2,10],
             wbar_f11=M1$wbar[1,11],wbar_m11=M1$wbar[2,11],
             wbar_f12=M1$wbar[1,12],wbar_m12=M1$wbar[2,12],
             wbar_f13=M1$wbar[1,13],wbar_m13=M1$wbar[2,13],
             wbar_f14=M1$wbar[1,14],wbar_m14=M1$wbar[2,14],
             wbar_f15=M1$wbar[1,15],wbar_m15=M1$wbar[2,15],
             wbar_f16=M1$wbar[1,16],wbar_m16=M1$wbar[2,16],
             wbar_f17=M1$wbar[1,17],wbar_m17=M1$wbar[2,17],
             wbar_f18=M1$wbar[1,18],wbar_m18=M1$wbar[2,18],
             wbar_f19=M1$wbar[1,19],wbar_m19=M1$wbar[2,19],
             wbar_f20=M1$wbar[1,20],wbar_m20=M1$wbar[2,20]) #females=sex1, males=sex2
    
    return(out)
  })
}

	# |---------------------------------------------------------------------------|
	# | Main function calls.                                                      |
	# |---------------------------------------------------------------------------|
	# | df is a data frame with the output of the equilibrium model, 
	# /  includes all permutations of the variables below
	# | 
	df <- expand.grid(fe=fe,slim=sizelimit,dm=discardmortality,bycatch=bycatch)
	
	#	if(!exists("S1")) ## so it doesn't rerun all models just to update graphs
	S1 <- as.data.frame(t(apply(df,1,.runModel)))
	
	# require("parallel")
	# cl <- makeCluster(getOption("cl.cores", 4))
	# # parApply(cl = NULL, X, MARGIN, FUN, ...)
	# S2 <- as.data.frame(t(parApply(cl,df,1,.runModel)))
	sl.fspr <- ddply(S1,.(slim,dm,bycatch),plyr::summarize,
	                 fspr = fe[max(which(SPR>=0.3))])
	
	
	sl.fmsy <- ddply(S1,.(slim,dm,bycatch),plyr::summarize,
	                 fmsy=fe[which.max(Ye)],
	                 Fd = round(Fd[which.max(Ye)],2),
	                 bmsy=round(Be[which.max(Ye)],1),
	                 spr =round(SPR[which.max(Ye)],3),
	                 msy=round(max(Ye),1),
	                 cbar=round(Cbar[which.max(Ye)],1) )
	
	t2<-subset(subset(sl.fmsy,slim%in%c(66.04,81.28)),bycatch%in%c(10,20))
	require("Hmisc")
	#	latex(t2,rowname=NULL)
	
	#Economic losses
	econ <- ddply(S1,.(slim,dm,bycatch),plyr::summarize,
	              fe =   fe[which.max(YEv)],
	              sl = slim[which.max(YEv)],
	              YEv = YEv[which.max(YEv)],
	              DEv = DEv[which.max(YEv)],
	              BYv = BYv[which.max(YEv)])
	
	sdf <- subset(S1,slim==81.28)
	t3<-ddply(sdf,.(bycatch,dm),plyr::summarize,
	          msy  =round(Ye[which.max(Ye)],1),
	          fmsy =round(fe[which.max(Ye)],2),
	          Value=round(YEv[which.max(Ye)],1),
	          DEv  =round(DEv[which.max(Ye)],1),
	          BYv  =round(BYv[which.max(Ye)],1))
	#	latex(t3,rowname=NULL)
	
	
	# |---------------------------------------------------------------------------|
	# | GRAPHICS - CHANGES TO SIZE-AT_AGE                                         |
	# |---------------------------------------------------------------------------|
	
	colnames(S1)
	cols <- c("fe", "wbar_f1",  "wbar_m1",  "wbar_f2",  "wbar_m2",
	          "wbar_f3",  "wbar_m3",  "wbar_f4",  "wbar_m4",  "wbar_f5",  "wbar_m5", 
	          "wbar_f6",  "wbar_m6",  "wbar_f7",  "wbar_m7",  "wbar_f8",  "wbar_m8", 
	          "wbar_f9",  "wbar_m9",  "wbar_f10", "wbar_m10", "wbar_f11", "wbar_m11",
	          "wbar_f12", "wbar_m12", "wbar_f13", "wbar_m13", "wbar_f14", "wbar_m14",
	          "wbar_f15", "wbar_m15", "wbar_f16", "wbar_m16", "wbar_f17", "wbar_m17",
	          "wbar_f18", "wbar_m18", "wbar_f19", "wbar_m19", "wbar_f20", "wbar_m20")
	dff <- S1[,cols]
	colnames(dff)

	sexone <- c(rep("F",101),rep("M",101))
	sextwo <- rep(sexone,20)
	
	melted <- melt(dff, id="fe" )
	melted$sex <- sextwo
	
	ageone <- c(rep(1,202),
	            rep(2,202),
	            rep(3,202),
	            rep(4,202),
	            rep(5,202),
	            rep(6,202),
	            rep(7,202),
	            rep(8,202),
	            rep(9,202),
	            rep(10,202),
	            rep(11,202),
	            rep(12,202),
	            rep(13,202),
	            rep(14,202),
	            rep(15,202),
	            rep(16,202),
	            rep(17,202),
	            rep(18,202),
	            rep(19,202),
	            rep(20,202))
	melted$age <- ageone
	
	melted$waa <- melted$value
	
	cols <- c("fe","sex","age","waa")
	dff <- melted[,cols]
	
	head(dff,101)
	# ff <- c("0.00","0.10","0.20","0.30","0.40","0.50",
	#         "0.60","0.70","0.80","0.90","1.00")
	
	ff <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
	
	# melt(dfftry,id_vars="fe")
	
	
# Summaries of observed survey weight-at-age (2014 to 1980s - 1980s data was pretty spotty, so I combined it all.)
	dat <- read.csv("src/HALITOSISII/WAA14vs80.csv")
	dat <- subset(dat, reg=="Full")
	
	cols <- c("age","weight","sex","year")
	dat <- dat[,cols]
	dat$fe <- rep(0.5,69)
	
	levels(dat$sex) <- c("Females","Males")
	colnames(dat) <- c("age","waa","sex","year","fe")
	dat14 <- subset(dat, year=="2014")
	dat80 <- subset(dat, year=="1980s")
	datall <- rbind(dat14,dat80)
	datall <- na.omit(datall)
	levels(dfftry$sex) <- c("Females","Males")

	# BASE PLOT

	tiff("src/HALITOSISII/figs/cw_weightraj_nc.tiff", height = 4, family="serif", pointsize=12,width = 7, units = 'in', res=300) #, compression="lzw"
	
	p <- ggplot(dfftry, aes(x=age,y=waa,size=as.factor(fe))) +
	  geom_line(colour="black") + 
	  facet_wrap(~sex, scales="free") +
	  scale_size_manual(values=c(0,0.2,0.3,0.4,0.6,0.8,0.9,1,1.2,1.4,1.8)) + #, guide=FALSE)
	  
	  # Adding observational data
	  # guides(size = guide_legend(override.aes = list(colour="black"))) +
	  # geom_point(data=datall,aes(x=age,y=waa, colour= year),size=2) +
	  # facet_wrap(~sex, scales="free") +
	  # geom_smooth(data=datall,aes(x=age,y=waa, colour= year),
	  #                    se=FALSE)
	  
	  mytheme +
	  labs(x = "\n Age", y = "Weight (kg) \n", size = 'Fishing \n mortality') +
	  guides(fill=guide_legend(nrow=2,byrow=TRUE))
	
	print(p)
	
	dev.off()


	# Trajectories with observations 

	tiff("src/HALITOSISII/figs/cw_weightrajwithobs_bw.tiff", height = 4, family="serif", pointsize=12,width = 7, units = 'in', res=300, compression="lzw")
	
	p <- p + geom_line(colour="darkgrey") +
	  geom_point(data=datall,aes(x=age,y=waa, shape= year),colour="black",size=2) +
	  # geom_smooth(data=datall,aes(x=age,y=waa, shape= year),colour="black",
	  #                      se=FALSE,size=0.8) +
	  guides(size = guide_legend(override.aes = list(colour="darkgrey"))) +
	  scale_linetype_discrete("Year") +
	  scale_shape_discrete("Year")
	
	print(p)
	
	dev.off()
	
	# Trajectories with observations + best guess of historical F
	
	dffreal <- subset(dfftry, fe=="0.4")
	
	tiff("src/HALITOSISII/figs/cw_obs_histF_bw.tiff", height = 4, family="serif", pointsize=12,width = 7, units = 'in', res=300, compression="lzw")
	
	p <- p + geom_line(data=dffreal, aes(x=age,y=waa),#,linetype=as.factor(fe)),
	                   linetype=2, colour = "black", size = 1.25) + 
	  guides(linetype = guide_legend(override.aes = list(linetype = 2))) +
	  labs(linetype = 'Historical F') 
	
	print(p)
	
	dev.off()
	
	# Base R plot - changes in weight-at-age - never used this one in the manuscript.
	
	cols <- c("fe",  "wbar_f5",  "wbar_m5", 
	          "wbar_f10", "wbar_m10", 
	          "wbar_f15", "wbar_m15", 
	          "wbar_f20", "wbar_m20")
	dff <- S1[,cols]
	
	par(mfrow=c(1,2))
	with(dff, plot(fe, wbar_f20, 'l', col=2, lwd=2, ylim=c(0,max(wbar_f20)),# ylim=c(0,150),
	               main="Females",
	               ylab="Weight (lb)",
	               xlab="Fishing intensity", yaxt="n"))
	axis(2, at=seq(0,150,20))
	with(dff, lines(fe, wbar_f15, 'l', col=3, lwd=2))
	with(dff, lines(fe, wbar_f10, 'l', col=4, lwd=2))
	with(dff, lines(fe, wbar_f5, 'l', col=5, lwd=2))
	abline(v=0.34, lty=2)
	abline(h=43.4, lty=2, col=2)
	abline(h=23.7, lty=2, col=3)
	abline(h=11.0, lty=2, col=4)
	abline(h=4.5, lty=2, col=5)
	grid()
	legend("topright", legend=c("age-20","age-15","age-10","age-5"),
	       col=c(2,3,4,5), lwd=c(2,2,2,2))
	
	with(dff, plot(fe, wbar_m20, 'l', col=2, lwd=2,ylim=c(0,max(wbar_m20)),#ylim=c(0,100),
	               main="Males",
	               ylab="Weight (lb)",
	               xlab="Fishing intensity",yaxt="n"))
	axis(2, at=seq(0,100,20))
	with(dff, lines(fe, wbar_m15, 'l', col=3, lwd=2))
	with(dff, lines(fe, wbar_m10, 'l', col=4, lwd=2))
	with(dff, lines(fe, wbar_m5, 'l', col=5, lwd=2))
	grid()
	abline(v=0.34, lty=2)
	abline(h=15.5, lty=2, col=2)
	abline(h=11.2, lty=2, col=3)
	abline(h=7.1, lty=2, col=4)
	abline(h=4.2, lty=2, col=5)
	legend("topright", legend=c("age-20","age-15","age-10","age-5"),
	       col=c(2,3,4,5), lwd=c(2,2,2,2))
	
	
# |---------------------------------------------------------------------------|
# | MORE GRAPHICS. - for LWS 2015 paper                                       |
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






