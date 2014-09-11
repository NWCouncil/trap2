dam.plot <- function(plant.df, data.col.name, ...) {
  args <- list(...)
  library(grid)
  plot.df <- data.frame(names=c("H HORS", "KERR", "THOM F", "NOXON", "CAB G", "ALBENI", "BOX C", "BOUND", "LIBBY", 
                                "COULEE", "CH JOE", "WELLS", "CHELAN", "R RECH", "ROCK I", "WANAP", "PRIEST", "BRNLEE", 
                                "OXBOW", "HELL C", "DWRSHK", "LR.GRN", "L GOOS", "LR MON", "ICE H", "MCNARY", "J DAY", 
                                "RND B", "PELTON", "DALLES", "BONN", "SWFT 1", "SWFT 2", "YALE", "MERWIN"),
                       xcoords=c(0.9, 0.9, 0.875, 0.825, 0.75, 0.675, 0.6, 0.575, 0.675, 0.475, 0.4, 0.3, 0.175, 0.3, 
                                 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.725, 0.575, 0.5, 0.425, 0.35, 0.3, 0.225, 0.15, 0.15,
                                 0.15, 0.075, 0.2, 0.2, 0.15, 0.075),
                       ycoords=c(0.85, 0.775, 0.7, 0.75, 0.75, 0.75, 0.8, 0.85, 0.9, 0.8, 0.8, 0.75, 0.725, 0.7, 0.65,
                                 0.6, 0.55, 0.2, 0.25, 0.3, 0.45, 0.4, 0.4, 0.4, 0.4, 0.25, 0.25, 0.075, 0.125, 0.25, 0.25,
                                 0.5, 0.45, 0.4, 0.4),
                       branchtop=c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 
                                   FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE,
                                   FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                   TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
                       hjust=c(0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1),
                       rot=c(0, 0, 45, 45, 45, 45, 45, 45, 45, 45, 45, 0, 45, 0, 0, 0, 0, 0, 0, 0, 45, 45, 45, 45, 45, 45, 45, 0, 0, 45, 45, 0, 0, 45, 45))

  plot.df <- merge(plot.df, plant.df, all.x=TRUE, sort=FALSE)

  vp <- viewport()

  direct.line <- function(vp, plot.df, coords, line.weight) {
    grid.lines(plot.df$xcoords[coords], 
             plot.df$ycoords[coords], vp=vp, gp=gpar(lex=line.weight[coords[1]]))
  }
  along.y <- function(vp, plot.df, coords, line.weight) {
    grid.lines(c(rep(plot.df$xcoords[coords[1]], 2), plot.df$xcoords[coords[2]]), 
               c(plot.df$ycoords[coords[1]], rep(plot.df$ycoords[coords[2]],2)), vp=vp,
               gp=gpar(lex=line.weight[coords[1]]))
  }
  along.x <- function(vp, plot.df, coords, line.weight) {
    grid.lines(c(plot.df$xcoords[coords[1]], rep(plot.df$xcoords[coords[2]], 2)), 
               c(rep(plot.df$ycoords[coords[1]], 2), plot.df$ycoords[coords[2]]), vp=vp,
               gp=gpar(lex=line.weight[coords[1]]))
  }
  s.line <- function(vp, plot.df, coords, line.weight, along.x=TRUE) {
    half.x.diff <- (plot.df$xcoords[coords[1]] - plot.df$xcoords[coords[2]])/2
    half.y.diff <- (plot.df$ycoords[coords[1]] - plot.df$ycoords[coords[2]])/2
    if (along.x) {
      grid.lines(c(plot.df$xcoords[coords[1]], rep(plot.df$xcoords[coords[1]] - half.x.diff, 2), plot.df$xcoords[coords[2]]),
                 rep(c(plot.df$ycoords[coords[1]], plot.df$ycoords[coords[2]]),each=2), vp=vp,
                 gp=gpar(lex=line.weight[coords[1]]))
    } else {
      grid.lines(rep(c(plot.df$xcoords[coords[1]], plot.df$xcoords[coords[2]]),each=2),
                 c(plot.df$ycoords[coords[1]], rep(plot.df$ycoords[coords[2]] - half.y.diff, 2), plot.df$ycoords[coords[2]]), vp=vp,
                 gp=gpar(lex=line.weight[coords[1]]))
    }
  }
  off.left <- function(vp, plot.df, coord, line.weight) {
    grid.lines(c(plot.df$xcoords[coord], 0.01),
               rep(plot.df$ycoords[coord], 2), vp=vp, gp=gpar(lex=line.weight[coord]))
  }

  line.weight <- 1 + 7 * plot.df[, data.col.name]/max(plot.df[, data.col.name], na.rm=TRUE)
  if (hasArg(no.lex)) {
    if (args$no.lex == TRUE) {
      line.weight <- rep(1, length(plot.df[, data.col.name]))
    } 
  }

  direct.line(vp, plot.df, 1:2, line.weight)
  along.y(vp, plot.df, 2:3, line.weight)
  s.line(vp, plot.df, 3:4, line.weight)
  direct.line(vp, plot.df, 4:5, line.weight)
  direct.line(vp, plot.df, 5:6, line.weight)
  along.x(vp, plot.df, 6:7, line.weight)
  along.y(vp, plot.df, 7:8, line.weight)
  s.line(vp, plot.df, c(8,10), line.weight)
  along.x(vp, plot.df, 9:10, line.weight)
  direct.line(vp, plot.df, 10:11, line.weight)
  along.x(vp, plot.df, 11:12, line.weight)
  direct.line(vp, plot.df, c(12, 14), line.weight)
  along.x(vp, plot.df, 13:14, line.weight)
  direct.line(vp, plot.df, 14:15, line.weight)
  direct.line(vp, plot.df, 15:16, line.weight)
  direct.line(vp, plot.df, 16:17, line.weight)
  along.y(vp, plot.df, c(17, 26), line.weight)
  direct.line(vp, plot.df, 18:19, line.weight)
  direct.line(vp, plot.df, 19:20, line.weight)
  along.y(vp, plot.df, c(20, 22), line.weight)
  s.line(vp, plot.df, 21:22, line.weight)
  direct.line(vp, plot.df, 22:23, line.weight)
  direct.line(vp, plot.df, 23:24, line.weight)
  direct.line(vp, plot.df, 24:25, line.weight)
  along.x(vp, plot.df, 25:26, line.weight)
  direct.line(vp, plot.df, 26:27, line.weight)
  direct.line(vp, plot.df, 28:29, line.weight)
  direct.line(vp, plot.df, 29:30, line.weight)
  direct.line(vp, plot.df, c(27, 30), line.weight)
  direct.line(vp, plot.df, 30:31, line.weight)
  off.left(vp, plot.df, 31, line.weight)
  direct.line(vp, plot.df, 32:33, line.weight)
  along.y(vp, plot.df, 33:34, line.weight)
  direct.line(vp, plot.df, 34:35, line.weight)
  off.left(vp, plot.df, 35, line.weight)

  grid.rect(plot.df$xcoords, plot.df$ycoords, width=.04, height=.025, gp=gpar(fill="white"), vp=vp)
  y.adder <- as.numeric(plot.df$rot == 45) * .02
  y.adder <- y.adder * sapply(plot.df$hjust, function(x){if(x==0){1} else {-1}})
  x.adder <- as.numeric(plot.df$rot == 0) * .03
  x.adder <- x.adder * sapply(plot.df$hjust, function(x){if(x==0){1} else {-1}})
  grid.text(plot.df$names, plot.df$xcoords + x.adder, plot.df$ycoords + y.adder, 
            hjust=plot.df$hjust, rot=plot.df$rot, vp=vp, gp=gpar(cex=.8))
  grid.text(sprintf("%.1f", plot.df[, data.col.name]), plot.df$xcoords,
            plot.df$ycoords, vp=vp, gp=gpar(cex=.8))
  if (hasArg(title.text)) {
    grid.text(args$title.text, 0.05, 0.95, vp=vp, gp=gpar(cex=1.2), just="left")
  }

}

setwd("Q:/BK/projects/2014/Flex/HydroFlex")

waterfile <- file("LP.OUT", "rb")
varnamefile <- file("LPVARS", "rt")

varnamesraw <- readLines(varnamefile, n=1)
plantnamesraw <- readLines(varnamefile, n=1)

varnames <- unlist(strsplit(varnamesraw, "[[:space:]]+"))
varnames <- varnames[nchar(varnames) > 0]

plantnames <- substring(plantnamesraw, 2)
plantnames <- regmatches(plantnames, gregexpr(".{6}", plantnames))[[1]]
plantnames <- gsub("^ *| *$", "", plantnames)
plantnames <- plantnames[plantnames!=""]

#eof <- FALSE
#water.matrix <- matrix(ncol=length(varnames), nrow=0)
#while(eof == FALSE) {
  line <- readBin(waterfile, numeric(), size=8, n=length(varnames))
#  eof <- length(line) == 0
#  if(!eof) {
#    water.matrix <- rbind(water.matrix, line)
#  }
#}

plant.df <- data.frame(names=plantnames, 
                       nums=substr(varnames[grep("TN", varnames)], 2, 3),
                       peak.turb=line[grep("TN", varnames)]/1000, 
                       off.turb=line[grep("TF", varnames)]/1000,
                       peak.spill=line[grep("SN", varnames)]/1000,
                       off.spill=line[grep("SF", varnames)]/1000)
storage.df <- data.frame(nums=substr(varnames[grep("S0", varnames)], 2, 3),
                         S0=line[grep("S0", varnames)])
storage.df <- merge(storage.df, data.frame(nums=substr(varnames[grep("S1", varnames)], 2, 3),
                                           S1=line[grep("S1", varnames)]), all.x=TRUE)
storage.df <- merge(storage.df, data.frame(nums=substr(varnames[grep("S2", varnames)], 2, 3),
                                           S2=line[grep("S2", varnames)]), all.x=TRUE)
storage.df$off.storage <- (storage.df$S1 - storage.df$S0)/16/1000
storage.df$peak.storage <- (storage.df$S2 - storage.df$S1)/8/1000
plant.df <- merge(plant.df, storage.df[, c("nums", "peak.storage", "off.storage")], all.x=TRUE)

pdf("Per1WY1929.pdf", width=11, height=8.5, version="1.4")
dam.plot(plant.df, "peak.turb", title.text="Peak Turbine Flows - Period 1 - WY 1929")
grid.newpage()
dam.plot(plant.df, "off.turb", title.text="Off Peak Turbine Flows - Period 1 - WY 1929")
grid.newpage()
dam.plot(plant.df, "peak.spill", title.text="Peak Spill Flows - Period 1 - WY 1929")
grid.newpage()
dam.plot(plant.df, "off.spill", title.text="Off Peak Spill Flows - Period 1 - WY 1929")
grid.newpage()
dam.plot(plant.df, "peak.storage", title.text="Peak Hourly Storage - Period 1 - WY 1929", no.lex=TRUE)
grid.newpage()
dam.plot(plant.df, "off.storage", title.text="Off Peak Hourly Storage - Period 1 - WY 1929", no.lex=TRUE)
dev.off()

