# R-script for looking at the size at age data.
# Thanks to Aaron Ranta for providing the data.

#hal.data          <- read.csv(file="LenFrequency98-11.CSV", header=TRUE, skip=7, na.strings="NULL")
hal.data       <- read.table(file="HalLenAge1998-2011bySex.txt", header=TRUE, na.strings="NULL")
hal.data$byear    <- hal.data$year - hal.data$bbage


df2b     <- subset(hal.data,RegArea=="2B",select=c(year,byear,frklen,bbage,srfage,sex))
df2b     <- subset(df2b, bbage>=6)
df2b     <- subset(df2b, bbage<=24)
p <- ggplot(df2b) + geom_point(aes(x=bbage,y=frklen, col=sex), size=1) 
p <- p + geom_point(aes(x=srfage, y=frklen, col=sex), size=1, pch=3)
p <- p + facet_wrap(~year)
p

p<-ggplot(df2b,aes(x=bbage, y=srfage)) + geom_point()
p<-p+geom_abline(slope=1)
p


# plot mean length-at-age over time
# FEMALES
p <- ggplot(subset(hal.data, sex=="F", na.rm=TRUE))
p <- p + geom_point(aes(year, frklen, col=bbage), size=1)
p <- p + facet_wrap(~RegArea)


p <- ggplot(df2b) #+ geom_violin(aes(x=factor(year), y=length, col=byear))
p <- p + geom_point(aes(x=year, y=frklen, col=sex), shape=1, size=0.5)
#p <- p + geom_boxplot(aes(factor(year), length), point=".")
p <- p + geom_smooth(aes(year, frklen, col=sex))
p <- p + facet_wrap(~bbage, scales="free_y")


# Construct SRF vs BB age for estimating BB age from SRF.
p<-ggplot(hal.data)+geom_point(aes(bbage,srfage),size=0.5,position="jitter") +geom_abline(slope=1,col=2)
p<-p+stat_smooth(aes(bbage, srfage), se=TRUE)
p

dt    <- cbind(X$bbage, X$srfage)
write(dim(dt)[1], file="BBSRFage.dat")
write("# BB SRF", file="BBSRFage.dat", append=TRUE)
write.table(dt, file="BBSRFage.dat", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)


# Halibut length-age data for a given year.
iyr = 2009
sx          <- c("F","M")
sa          <- c("2A","2B","2C","3A","3B","4A","4B","4C","4D")
df          <- subset(hal.data, year==iyr)
df          <- subset(df, sex!="U")

df$statarea <- match(df$RegArea,sa)
p           <- ggplot(df,aes(bbage,frklen,col=sex)) + geom_point(size=1) + facet_wrap(~RegArea)

dfile       <- paste("LengthAge", iyr, ".dat", sep="")
write(dim(df)[1], dfile)
data        <- cbind(area = df$statarea, age = df$bbage, fl = df$frklen,  sex = match(df$sex, sx))
narea		<- length(unique(df$statarea))
write(narea, dfile, append=TRUE)
nsex		<- length(unique(data[, 4]))
write(nsex, dfile, append=TRUE)
write("#area\t age\t fl\t sex", file=dfile, append=TRUE)
write.table(data, file=dfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")