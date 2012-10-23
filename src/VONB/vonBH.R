# vonBH.R
require(Riscam)
require(ggplot2)
A <- read.admb("vonBH")
A$rp <- read.table("vonBH.mcmc", header=TRUE)

sx          <- c("Female","Male")
sa          <- c("2A","2B","2C","3A","3B","4A","4B","4C","4D")
A$rp$area   <- sa[A$rp$area]
A$rp$sex    <- sx[A$rp$sex]
DF <- data.frame(area=sa[A$area],sex=sx[A$sex], age=A$age, fl=A$fl, fl_hat=A$fl_hat, epsilon=A$epsilon)

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

p  <- ggplot(DF, aes(age, epsilon, col=sex)) +geom_point(size=1)+geom_smooth()
p  <- p + facet_wrap(~area) + labs(col="Sex")+xlab("Age (years)")+ylab("Residual fork length (cm)")
p3 <- p

p1 
p2 
p3

# qq plots
p  <- ggplot(DF, aes(sample=epsilon, col=sex)) + stat_qq()
p  <- p + facet_wrap(~area)
p4 <- p
p4

p  <- ggplot(DF, aes(sample=epsilon, col=sex)) + stat_qq()
p  <- p + facet_wrap(~age)
p5 <- p
p5


# Allometric relationship.
# W = a L ^b; a 6.82e-6 b=3.24
a  <- 6.83e-6
b  <- 3.24
p  <- ggplot(DF, aes(age, a*fl^b, col=area)) +geom_point(size=1)
p  <- p + geom_point(aes(age, a*fl_hat^b, col=area), size=2) + facet_wrap(~sex)
p  <- p + labs(col="Regulatory\nArea") + xlab("Age (years)") + ylab("Weight-at-age (lb)")
p6 <- p


# Marginal posteriors for F0.1

p  <- ggplot(A$rp, aes(x=F0.1, fill=sx[sex])) +geom_density(adjust=2, size=0.0)
p  <- p +labs(fill="Sex") + xlim(c(0, 1.0)) + facet_wrap(~area) 
p7 <- p
p7

p  <- ggplot(A$rp, aes(x=F0.1, fill=area)) +geom_density(adjust=2, size=0.1, alpha=0.5)
p  <- p  + facet_wrap(~sex, scales="free") 
p
