library(ggplot2)
setwd("Q:/ClusterFiles/Trap2/")
hoss <- read.csv("tools/HOSS4hr4KComp.csv", stringsAsFactors=FALSE, header=TRUE)
names(hoss) <- c("WY", "Period", "Energy", "Peak", "Min")
hoss$Source <- "HOSS"
rods <- read.csv("tools/RODS4hr4KComp.csv", stringsAsFactors=FALSE, header=TRUE)
names(rods) <- c("WY", "Period", "Energy", "Peak", "Min")
rods$Source <- "Historic"
comb <- rbind(hoss, rods)

trap <- read.fwf("outputs/run 4K-4hr-50pct.OUT", widths=c(4, 8, rep(7, 11), 5), skip=3)
names(trap) <- c("PER", "TM_E", "EON", "EOFF", "TM_W", "WON", "WOFF", "TM_I", "IDON", "IDOFF", "TM_FD", "FDON", "FDOFF", "IWY")
trap$OUT <- rep(seq(1,4), length(trap[,1])/4)
trap.avg <- aggregate(cbind(TM_FD, FDON, FDOFF) ~ IWY + PER, data=trap, FUN=mean)
names(trap.avg) <- c("WY", "Period", "Energy", "Peak", "Min")
trap.avg$Source <- "TRAP"
comb <- rbind(comb, trap.avg)
comb$Period <- as.factor(comb$Period)
levels(comb$Period) <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr1", "Apr2", "May", "Jun", "Jul", "Aug1", "Aug2", "Sep")

comp.plot <- ggplot(data=comb) + geom_point(aes(x=Energy, y=Peak, colour=Source, shape="Peak"))
comp.plot <- comp.plot + geom_point(aes(x=Energy, y=Min, colour=Source, shape="Min"))
comp.plot <- comp.plot + facet_wrap(~ Period)
comp.plot <- comp.plot + scale_shape_discrete(name="")
pdf("outputs/FedComp.pdf", width=14, height=8.5, version="1.4")
comp.plot
dev.off()


