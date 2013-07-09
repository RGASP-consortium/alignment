## Make figure comparing transcript performance using all vs only primary alignments

source(file.path(Sys.getenv("SPALE_HOME"), "scripts", "R", "spale_func.R"))

sim.dir <- get.spale.data.fn("sim_logs")
stats.dir <- get.spale.data.fn("output", "cufflinks", "eval", "ce.stranded")

read.general.stats <- function(cufflinks.ver, aln.set, data.set, method.info, ...) {

  ## Read data
  stats.fn <- file.path(stats.dir, paste("ce", data.set, cufflinks.ver, aln.set, ..., "stats.txt", sep="."))
  print(stats.fn)  
  stats <- read.delim(stats.fn, quote="", as.is=TRUE)

  ## Set method names
  x <- strsplit(stats$submission, "/", fixed=TRUE)
  team <- sapply(x, function(x) x[length(x)-1])
  i <- substring(stats$submission, nchar(stats$submission))
  i <- match(paste(method.info$team, method.info$number), paste(team, i))
  stats <- stats[i, ]
  stats$method <- method.info$name

  return(stats)
}

compute.tx.sens <- function(x) {
  n.pred.true <- tapply(x$ref.tx.corr, x$method, sum)
  n.pred <- tapply(x$ref.tx, x$method, sum)
  return(n.pred.true / n.pred)
}

compute.tx.prec <- function(x) {
  n.pred.true <- tapply(x$pred.tx.corr, x$method, sum)
  n.pred <- tapply(x$pred.tx, x$method, sum)
  return(n.pred.true / n.pred)
}

compute.exon.sens <- function(x) {
  n.pred.true <- tapply(x$ref.exons.corr, x$method, sum)
  n.pred <- tapply(x$ref.exons, x$method, sum)
  return(n.pred.true / n.pred)
}

compute.exon.prec <- function(x) {
  n.pred.true <- tapply(x$pred.exons.corr, x$method, sum)
  n.pred <- tapply(x$pred.exons, x$method, sum)
  return(n.pred.true / n.pred)
}

## Read method info
method.info <- read.method.info(include.exploratory=FALSE)

## Read dataset info
datasets <- read.dataset.info()
datasets <- datasets[ datasets$id != "mouse", ] # RGASP3-specific; we did not run Cufflinks for mouse data  

# Open figure file for output
pdf("cufflinks_abacus-plot_pri_vs_all.pdf", height=11)

## Process each dataset
for(i in 1:length(datasets)) {

  prec.stats.all <- read.general.stats("v2", "all", datasets$id[i], method.info)
  prec.stats.pri <- read.general.stats("v2", "primary", datasets$id[i], method.info)

  if(substring(datasets$id[i], 1, 3) == "sim") {
    sens.stats.all <- prec.stats.all
    sens.stats.pri <- prec.stats.pri
  } else {
    sens.stats.all <- read.general.stats("v2", "all", datasets$id[i], method.info, "coding")
    sens.stats.pri <- read.general.stats("v2", "primary", datasets$id[i], method.info, "coding")
  }
  
  ## Transcript stats
  
  sens.all <- compute.tx.sens(sens.stats.all)
  sens.pri <- compute.tx.sens(sens.stats.pri)
  prec.all <- compute.tx.prec(prec.stats.all)
  prec.pri <- compute.tx.prec(prec.stats.pri)
 
  j <- length(prec.all):1
  par(mar=c(5, 10, 4, 2))
  plot(NULL, type="n", xlim=c(0, max(prec.all, prec.pri, sens.all, sens.pri, na.rm=TRUE)),
       ylim=range(j), xlab="", ylab="", yaxt="n", main=paste(datasets$name[i], "tx"))
  axis(2, j, rownames(prec.all), las=1, tick=FALSE)
  abline(h=j, col="gray", lwd=2)
  points(sens.all, j, pch=21, col="#FF0000", bg="#FF0000", cex=3)
  points(sens.pri, j, pch=1, cex=3)
  points(prec.all, j, pch=23, col="#FF0000", bg="#FF0000", cex=3)
  points(prec.pri, j, pch=5, cex=3)

  ## Exon stats

  sens.all <- compute.exon.sens(sens.stats.all)
  sens.pri <- compute.exon.sens(sens.stats.pri)
  prec.all <- compute.exon.prec(prec.stats.all)
  prec.pri <- compute.exon.prec(prec.stats.pri)
 
  j <- length(prec.all):1
  par(mar=c(5, 10, 4, 2))
  plot(NULL, type="n", xlim=c(0, max(prec.all, prec.pri, sens.all, sens.pri, na.rm=TRUE)),
       ylim=range(j), xlab="", ylab="", yaxt="n", main=paste(datasets$name[i], "exon"))
  axis(2, j, rownames(prec.all), las=1, tick=FALSE)
  abline(h=j, col="gray", lwd=2)
  points(sens.all, j, pch=21, col="#FF0000", bg="#FF0000", cex=3)
  points(sens.pri, j, pch=1, cex=3)
  points(prec.all, j, pch=23, col="#FF0000", bg="#FF0000", cex=3)
  points(prec.pri, j, pch=5, cex=3)

  ## Report exon spec and sens for primary alignments
  cat("[", datasets$name[i], "] Max exon sens =", max(sens.pri, na.rm=TRUE), "\n", sep=" ")
  cat("[", datasets$name[i], "] Max exon prec =", max(prec.pri, na.rm=TRUE), "\n", sep=" ")
  
}
 
## Close PDFs
dev.off()

