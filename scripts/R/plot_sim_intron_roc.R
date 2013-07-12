## Plot ROC-like curve for splice junction discovery
## The counts of true and false junctions at each threshold are also written to a text file

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))
data.dir <- get.rgasp.ali.data.fn("output", "gap_stats", "primary", "iac")
glyph.cex <- 1.5

read.iac.stratified <- function(fn, set="all") {

  x <- read.delim(fn, skip=1)

  if(set == "known") {
    x <- x[x$class1 == "Known splice", ]
  } else if(set == "novel" ) {
    x <- x[x$class1 != "Known splice", ]
  } else {
    stopifnot(set == "all")
  }

  f <- factor(x$class2 == "Known splice")
  x <- x[, grepl("n\\d+", colnames(x), perl=TRUE)]
  counts <- sapply(split(x, f, drop=FALSE), colSums)
  counts <- apply(counts, 2, function(x) rev(cumsum(rev(x))))

  return( counts )
}


get.iac.method.info <- function(data.dir, dataset) {
  files <- dir(path=data.dir, pattern=glob2rx(paste("iac.*.",dataset,"_*.txt",sep="")))
  method.info <- get.method.info.for.files(files, order.methods=TRUE)
  method.info <- method.info[ method.info$style < 5, ]
  return(method.info)
}


make.roc.like.plot <- function(data.dir, dataset, xlim=rep(NA, 2), ylim=rep(NA, 2), legend=FALSE, intron.set="all") {
  ## Get filenames and method info
  method.info <- get.iac.method.info(data.dir, dataset)
 
  ## Read data
  data <- lapply(file.path(data.dir, method.info$file), read.iac.stratified, set=intron.set)
 
  ## Determine scale
  xlim.full <- range(sapply(data, function(x) x[,1]))
  ylim.full <- range(sapply(data, function(x) x[,2]))
  if(is.na(xlim[1])) xlim[1] <- xlim.full[1]
  if(is.na(xlim[2])) xlim[2] <- xlim.full[2]
  if(is.na(ylim[1])) ylim[1] <- ylim.full[1]
  if(is.na(ylim[2])) ylim[2] <- ylim.full[2]
  
  ## Plot data
  plot(0, xlim=xlim, ylim=ylim, type="n", xlab="False junctions", ylab="True junctions", main=paste(dataset, intron.set))
  truth <- method.info$name == "Truth"
  for(i in which(!truth)) {
    lines(data[[i]][,1], data[[i]][,2], col=method.info$color[i], lty=method.info$lty[i], lwd=1)
    points(data[[i]][1,1], data[[i]][1,2], cex=glyph.cex, pch=method.info$pch[i], bg=method.info$color[i])
  }
  abline(h=data[[which(truth)]][1,2], col="darkgray")
  
  ## Plot legend
  if(legend) {
    plot.new()
    legend("topleft", legend=method.info$name, col=method.info$color, lty=method.info$lty, lwd=1, seg.len=3, xpd=NA)
    legend("topright", legend=method.info$name, pt.bg=method.info$color, pch=method.info$pch, pt.cex=glyph.cex, xpd=NA)
  }
}


write.true.false.table <- function(data.dir, dataset, fn) {
  ## Get filenames and method info
  method.info <- get.iac.method.info(data.dir, dataset)
 
  ## Read data
  data <- lapply(file.path(data.dir, method.info$file), read.iac.stratified)

  ## Extraft true and false counts
  false.junc <- t(sapply(data, function(x) x[,1]))
  true.junc <- t(sapply(data, function(x) x[,2]))
  rownames(false.junc) <- method.info$name
  rownames(true.junc) <- method.info$name

  ## Write to file
  out <- file(fn, "w")
  cat("True junctions\n", file=out)
  write.table(true.junc, file=out, col.names=NA, sep="\t", quote=FALSE)
  cat("False junctions\n", file=out)
  write.table(false.junc, file=out, col.names=NA, sep="\t", quote=FALSE)
  close(out)
}


## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Process each dataset
for(dataset in datasets$id) {

  ## Make plots
  out.fn <- paste("intron_roc_", dataset, ".pdf", sep="")
  pdf(out.fn, width=4.5, height=6.75, pointsize=6)
  par(mfrow=c(3,2), lwd=0.6, cex=1, cex.main=1, font.main=1, cex.axis=5/6)
  make.roc.like.plot(data.dir, dataset)
  make.roc.like.plot(data.dir, dataset, xlim=c(0,30000), ylim=c(95000, NA))
  make.roc.like.plot(data.dir, dataset, intron.set="known")
  make.roc.like.plot(data.dir, dataset, xlim=c(0, 30000), ylim=c(68000, NA), intron.set="known")
  make.roc.like.plot(data.dir, dataset, intron.set="novel")
  make.roc.like.plot(data.dir, dataset, xlim=c(0,30000), ylim=c(26000, NA), intron.set="novel", legend=TRUE)
  dev.off()

  ## Write table of true and false junctions
  out.fn <- paste("intron_discovery_", dataset, ".txt", sep="")
  write.true.false.table(data.dir, dataset, out.fn)
}
