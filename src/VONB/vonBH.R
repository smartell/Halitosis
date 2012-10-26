# vonBH.R

# DATA SECTION
require(Riscam)
require(ggplot2)
A 						<- read.admb("vonBH")
A$rp 					<- read.table("vonBH.mcmc", header=TRUE)
colnames(A$post.samp)	<- A$fit$names[1:A$fit$nopar]
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

quartz("Size at age", width=10, height=6)

p  <- ggplot(DF, aes(age, fl, col=factor(sex))) + geom_point(size=1, alpha=0.3) 
p  <- p + geom_line(aes(age, fl_hat, col=factor(sex)), size=1) 
p  <- p + facet_wrap(~area)
p  <- p + labs(col="Sex") + xlab("Age (years)") + ylab("Fork length (cm)")
p1 <- p

p  <- ggplot(DF, aes(age, fl, col=area)) +geom_point(size=1)
p  <- p + geom_line(aes(age, fl_hat, col=area), size=1) + facet_wrap(~sex)
p  <- p + labs(col="Regulatory\nArea") + xlab("Age (years)") + ylab("Fork length (cm)")
p2 <- p

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
p  <- ggplot(A$rp, aes(x=F0.1, fill=sex)) +geom_density(adjust=2, size=0.0)
p  <- p +labs(fill="Sex") + xlim(c(0, 1.0)) + facet_wrap(~area) 
p7 <- p


p  <- ggplot(subset(A$rp, sex=="Female"), aes(x=F0.1, fill=area)) +geom_density(adjust=1.5, size=0.1, alpha=0.15)
p  <- p  + facet_wrap(~sex, scales="free") 
p8 <- p

