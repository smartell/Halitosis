# vonBH.R

# DATA SECTION
require(Riscam)
require(ggplot2)
A 						<- read.admb("vonBH")
A$rp 					<- read.table("vonBH.mcmc", header=TRUE)
colnames(A$post.samp)	<- A$fit$names[1:(A$fit$nopar-1)]
A$post.samp             <- as.data.frame(A$post.samp)

sx          <- c("Female","Male")
sa          <- c("2A","2B","2C","3A","3B","4A","4B","4C","4D")
A$rp$area   <- sa[A$rp$area]
A$rp$sex    <- sx[A$rp$sex]
DF          <- data.frame(
				area=sa[A$area],
				sex=sx[A$sex], 
				age=A$age, 
				fl=A$fl, 
				fl_hat=A$fl_hat, 
				epsilon=A$epsilon)




# GRAPHICAL SECTION
# Index:
# 		p1 <- length-at-age plots for male and female by area
# 		p2 <- length-age by sex overlayed areas
#       p3 <- residual distribution as box plots by sex
#       p4 <- qqplot on residuals by area and sex
# 		p5 <- qqplot on age residuals 
#		p6 <- weight-at-age fits by area and sex
#       p7 <- Marginal posteriors for F0.1 for male and female by area panels
#       p8 <- Marginal posterior plots for F0.1 for female by area.
graphics.off()
quartz("Size at age", width=10, height=7.0)
p  <- ggplot(DF, aes(age, fl, shape=factor(sex))) + geom_point(size=1, alpha=0.3) 
p  <- p + geom_line(aes(age, fl_hat, col=factor(sex)), size=0.5) 
p  <- p + facet_wrap(~area)
p  <- p + labs(shape="Observed", col="Predicted") + xlab("Age (years)") + ylab("Fork length (cm)")
p1 <- p + theme_bw(12)
 ggsave(p1, file="../../FIGS/fig:lengthAgeFit.pdf")
 ggsave(p1, file="../../FIGS/fig:lengthAgeFit.png")

p  <- ggplot(DF, aes(age, fl_hat, shape=area)) +geom_point(size=1.5) + scale_shape_manual(values=c(1:9))
p  <- p + geom_line(aes(age, fl_hat), size=0.25) + facet_wrap(~sex, ncol=1)
p  <- p + labs(shape="Area", linetype="Area") + xlab("Age (years)") + ylab("Fork length (cm)")
p2 <- p + theme_bw(12)
ggsave(p2, file="../../FIGS/fig:lengthAgeFitbySex.pdf")
ggsave(p2, file="../../FIGS/fig:lengthAgeFitbySex.png")

p  <- ggplot(DF, aes(factor(age), epsilon, col=sex))+geom_boxplot(outlier.size = 0.2)
p  <- p + facet_wrap(~area) + labs(col="Sex")+xlab("Age (years)")+ylab("Residual fork length (cm)")
p3 <- p

# qq plots
p  <- ggplot(DF, aes(sample=epsilon, col=sex)) + stat_qq()
p  <- p + facet_wrap(~area)
p4 <- p


p  <- ggplot(DF, aes(sample=epsilon, col=sex)) + stat_qq()
p  <- p + facet_wrap(~age)
p5 <- p



# Allometric relationship.
# W = a L ^b; a 6.82e-6 b=3.24
a  <- 6.83e-6
b  <- 3.24
p  <- ggplot(DF, aes(age, a*fl^b, col=area)) +geom_point(size=1)
p  <- p + geom_point(aes(age, a*fl_hat^b, col=area), size=2) + facet_wrap(~sex)
p  <- p + labs(col="Regulatory\nArea") + xlab("Age (years)") + ylab("Weight-at-age (lb)")
p6 <- p


# Marginal posteriors for F0.1
p  <- ggplot(A$rp, aes(x=F0.1, fill=sex)) +geom_density(adjust=2, size=0.0, alpha=0.5)
p  <- p +labs(fill="sex") + xlim(c(0, 1.0)) + facet_wrap(~area) 
p7 <- p


p  <- ggplot(subset(A$rp, sex=="Female"), aes(x=F0.1, fill=area)) +geom_density(size=0.1, alpha=0.25)
p  <- p  + facet_wrap(~sex, scales="free") 
p8 <- p


p  <- ggplot(subset(A$rp, select=c("sex", "area", "L1")), aes(x=L1, fill=area)) +geom_density(adjust=1.5, size=0.1, alpha=0.25)
p  <- p  + facet_wrap(~sex, scales="free") 
pL1 <- p

p  <- ggplot(subset(A$rp, select=c("sex", "area", "L2")), aes(x=L2, fill=area)) +geom_density(adjust=1.5, size=0.1, alpha=0.25)
p  <- p  + facet_wrap(~sex, scales="free") 
pL2 <- p

p  <- ggplot(subset(A$rp, select=c("sex", "area", "rho")), aes(x=-log(rho), fill=area)) +geom_density(adjust=1.5, size=0.1, alpha=0.25)
p  <- p +labs(x="Growth coefficient (k)") + facet_wrap(~sex, scales="free") 
pk <- p

# Marginal posteriors for hyper parameters L1,  L2,  rho,  b
require(reshape2)
mdf <- melt(A$rp,id.vars=c("sex","area"))
pM1 <- ggplot(subset(mdf,variable=="L1"),aes(x=value)) + geom_density(aes(col=area)) +facet_wrap(~sex)
pM2 <- ggplot(subset(mdf,variable=="L2"),aes(x=value)) + geom_density(aes(col=area)) +facet_wrap(~sex)
pM3 <- ggplot(subset(mdf,variable=="rho"),aes(x=value)) + geom_density(aes(col=area)) +facet_wrap(~sex)
pM4 <- ggplot(subset(mdf,variable=="b"),aes(x=value)) + geom_density(aes(col=area)) +facet_wrap(~sex)
pM5 <- ggplot(subset(mdf,variable=="cv"),aes(x=value)) + geom_density(aes(col=area)) +facet_wrap(~sex)


# Table for growth parameter estimates
# A.post.samp has the posterior samples for l1,  l2,  rho,  cv for each of the stocks.
narea <- length(unique(A$area))


S <- dget("SurveySelectivities.rda")
colnames(S)<-c("Year","Sex",paste(1:30))
#X <- S[, -1:-2]/apply(S[,-1:-2],1,max)
#S[, -1:-2] <- X
mdf <- melt(S, id.vars=c("Year","Sex"))
colnames(mdf)<-c("Year","Sex","Age","value")

age_breaks = seq(5, 30, by=5)

ggplot(mdf)+geom_line(aes(x=Age,y=value,group=Year,col=factor(Year)), size=1.25) +labs(x="Age",y="Selectivity",col="Year") + facet_wrap(~Sex)+ theme_bw(18)+ scale_x_discrete(breaks=age_breaks)

ggplot(subset(mdf, Year==2011))+geom_line(aes(x=Age,y=value,group=Sex, col=Sex), size=1.25) +labs(x="Age",y="Selectivity",col="Sex") + theme_bw(18) + theme(axis.text.x=element_text(size=rel(0.5)))+opts(legend.position="top")
