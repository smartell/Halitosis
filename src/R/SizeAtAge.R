`# R-script for looking at the size at age data.

# Thanks to Aaron Ranta for providing the data.


setwd("/Users/Jane/Documents/GitHub/Halitosis/src/R/")
require(ggplot2)
require(reshape2)
hal.data       <- read.table(file="HalLenAge1998-2011bySex.txt", header=TRUE, na.strings="NULL")
hal.data$byear <- hal.data$year - hal.data$bbage


# Halibut length-age data for a given year.
write.LengthAgeData <- function(iyr=2011, fyr=NULL)
{
	#iyr = 2009
	sx          <- c("F","M")
	sa          <- c("2A","2B","2C","3A","3B","4A","4B","4C","4D")
	df          <- subset(hal.data, year==iyr)
	if(!is.null(fyr))
	{
		df          <- subset(hal.data, year>=fyr)
		df          <- subset(df, year<=iyr)
	}
	df          <- subset(df, sex!="U")

	df$statarea <- match(df$RegArea,sa)
	p           <- ggplot(df,aes(bbage,frklen,col=sex)) + geom_point(size=1) + facet_wrap(~RegArea)

	dfile       <- paste("LengthAge", fyr, iyr, ".dat", sep="")
	write(dim(df)[1], dfile)
	data        <- cbind(area = df$statarea, age = df$bbage, fl = df$frklen,  sex = match(df$sex, sx))
	narea		<- length(unique(df$statarea))
	write(narea, dfile, append=TRUE)
	nsex		<- length(unique(data[, 4]))
	write(nsex, dfile, append=TRUE)
	write("#area\t age\t fl\t sex", file=dfile, append=TRUE)
	write.table(data, file=dfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")	
}

write.LengtAgeKey <- function(iyr=2011)
{
	# Construct a table of the numbers at length (row) for a given age (col)
	# for a given year (iyr) from the setline survey. Note that for each
	# RegArea there will be two length-age tables (female and male),  so 
	# total of 18 tables will be produced.
	
	df    <- subset(hal.data, year==iyr)
	df    <- subset(df,  sex!="U")
	df    <- subset(df,  frklen!=0)
	
	xbin  <- seq(40, 215, by=5)
	df$FL <- xbin[findInterval(df$frklen, xbin)]
	
	abin  <- seq(4, 40)
	df$AGE<- abin[findInterval(df$bbage, abin)]
	
	dfm   <- melt(df, id=c("FL","AGE","sex","RegArea"))
	tx = acast(dfm,FL~AGE~RegArea~sex,length, fill=0)
	
	fn <- paste("../HGM/LenAgeKey", iyr, ".dat", sep="")
	write("#Dimensions, fl age area sex", file=fn)
	write(dim(tx), file=fn, append=TRUE)
	write("#Length intervals", file=fn, append=TRUE)
	write(dimnames(tx)[[1]],file=fn, append=TRUE)
	write("#Age vector", file=fn, append=TRUE)
	write(dimnames(tx)[[2]],file=fn, append=TRUE)
	
	ia <- dimnames(tx)[[3]]
	is <- dimnames(tx)[[4]]
	for(i in 1:2)
	for(j in 1:9)
	{
		id <- paste("# area",ia[j]," sex", is[i] )
		write(id, file=fn, append=TRUE)
		write.table(tx[, , j, i], file=fn, append=TRUE, col.names=FALSE, row.names=FALSE)
		write("#", file=fn, append=TRUE)
	}
		
}

sadf     <- subset(hal.data, bbage==5, select=c(year, frklen, sex, RegArea))
sadf     <- subset(sadf, sex!="U")
sadf     <- subset(sadf, RegArea!="4C")
sadf     <- subset(sadf, RegArea!="4E")

pp <- ggplot(sadf)
pp <- pp + geom_boxplot(aes(factor(year), frklen, col=sex), size=0.3, outlier.shape=".") 
pp <- pp + labs(x="Year", y="Fork length (age-15)")
pp <- pp + theme(axis.text.x = element_text(size = rel(0.8)))
pp <- pp + theme(axis.title = element_text(size = rel(1.8)))
pp <- pp + facet_wrap(~RegArea)
pp


# Change in size at age data to compare with equilibrium model.  
# For each reg area, show boxplot of female size-at-age over time 
# for ages 5, 10,  15, 20
saa  <- subset(hal.data, bbage==c(6, 10, 14), select=c(year, bbage, frklen, sex, RegArea))
saaf <- subset(saa, sex=="F")
saaf <- subset(saaf, RegArea==c("2B", "2C", "3A", "3B", "4A", "4B"))

a  <- 6.83e-6
b  <- 3.24

p  <- ggplot(saaf)
#p  <- p + geom_boxplot(aes(factor(year), a*frklen^b, fill=factor(bbage)), size=0.25, outlier.shape="")
p  <- p + geom_point(aes(x=year, a*frklen^b, shape=factor(bbage), col=factor(bbage)))
p  <- p + geom_smooth(aes(x=year, a*frklen^b, linetype=factor(bbage), fill=factor(bbage), shape=factor(bbage)))
p  <- p + labs(x="Year", y="Net weight (lb)", fill="Age", linetype="Age", shape="Age", col="Age")
#p  <- p + scale_x_discrete(breaks=pretty(saaf$year))
p  <- p + theme(axis.text.x = element_text(size = rel(0.8)))
p  <- p + facet_wrap(~RegArea) + theme_bw(12)
ggsave(p, file="../../FIGS/fig:SAA_age6_14.pdf")
ggsave(p, file="../../FIGS/fig:SAA_age6_14.png")


wraper<- function()
{
	#plots of data

df2b     <- subset(hal.data,RegArea=="3A",select=c(year,byear,frklen,bbage,srfage,sex))
df2b     <- subset(df2b, bbage>=5)
df2b     <- subset(df2b, bbage<=20)
df2b	 <- subset(df2b, sex!="U")
p <- ggplot(df2b) + geom_point(aes(x=bbage,y=frklen, col=sex), size=1, position="jitter") 
p <- p + geom_point(aes(x=srfage, y=frklen, col=sex), size=1, pch=3)
p <- p + facet_wrap(~byear)
p

p<-ggplot(df2b,aes(x=bbage, y=srfage)) + geom_point()
p<-p+geom_abline(slope=1)
p


# plot mean length-at-age over time
# FEMALES
p <- ggplot(subset(subset(hal.data, sex=="F", na.rm=TRUE), bbage==10))
#p <- p + geom_point(aes(year, frklen, col=bbage), size=1)
p <- p + geom_boxplot(aes(factor(year), frklen))
p <- p + facet_wrap(~RegArea)
p



p <- ggplot(df2b) #+ geom_violin(aes(x=factor(year), y=length, col=byear))
#p <- p + geom_point(aes(x=factor(year), y=frklen, col=sex), shape=1, size=0.5)
p <- p + geom_boxplot(aes(factor(year), frklen, col=sex), outlier.shape=".", size=0.2)
#p <- p + geom_smooth(aes(year, frklen, col=sex))
p <- p + facet_wrap(~bbage)
p

#plot weight-at-age data over time.
a  <- 6.83e-6
b  <- 3.24
df2b$weight <- a*df2b$frklen^b
p <- ggplot(df2b) + geom_boxplot(aes(factor(year), weight, col=sex), outlier.shape=".")
p <- p  facet_wrap(~bbage)
p

# Construct SRF vs BB age for estimating BB age from SRF.
p<-ggplot(hal.data)+geom_point(aes(bbage,srfage,col=year),size=0.5,position="jitter") +geom_abline(slope=1,col=2)
p<-p+stat_smooth(aes(bbage, srfage), se=TRUE)
p

dt    <- cbind(X$bbage, X$srfage)
write(dim(dt)[1], file="BBSRFage.dat")
write("# BB SRF", file="BBSRFage.dat", append=TRUE)
write.table(dt, file="BBSRFage.dat", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)



}