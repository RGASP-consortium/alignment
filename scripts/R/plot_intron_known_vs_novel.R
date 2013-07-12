## Plot ROC-like curves for annotated vs. novel junctions

options(width=150)

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))
data.dir <- get.rgasp.ali.data.fn("output", "gap_stats", "primary", "iac")
glyph.cex <- 1.5

read.iac.stratified <- function(fn) {
  x <- read.delim(fn, skip=1)

  f <- factor(x$class1 == "Known splice")
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


make.roc.like.plot <- function(data.dir, dataset, xlim=rep(NA, 2), ylim=rep(NA, 2), legend=FALSE) {
  ## Get filenames and method info
  method.info <- get.iac.method.info(data.dir, dataset)
 
  ## Read data
  data <- lapply(file.path(data.dir, method.info$file), read.iac.stratified)
 
  ## Determine scale
  xlim.full <- range(sapply(data, function(x) x[,1]))
  ylim.full <- range(sapply(data, function(x) x[,2]))
  if(is.na(xlim[1])) xlim[1] <- xlim.full[1]
  if(is.na(xlim[2])) xlim[2] <- xlim.full[2]
  if(is.na(ylim[1])) ylim[1] <- ylim.full[1]
  if(is.na(ylim[2])) ylim[2] <- ylim.full[2]
  
  ## Plot data
  plot(0, xlim=xlim, ylim=ylim, type="n", xlab="Novel junctions", ylab="Known junctions", main=dataset)
  for(i in 1:length(data)) {
    lines(data[[i]][,1], data[[i]][,2], col=method.info$color[i], lty=method.info$lty[i], lwd=1)
    points(data[[i]][1,1], data[[i]][1,2], cex=glyph.cex, pch=method.info$pch[i], bg=method.info$color[i])
  }
  
  ## Plot legend
  if(legend) {
    plot.new()
    legend("topleft", legend=method.info$name, col=method.info$color, lty=method.info$lty, lwd=1, seg.len=3, xpd=NA)
    legend("topright", legend=method.info$name, pt.bg=method.info$color, pch=method.info$pch, pt.cex=glyph.cex, xpd=NA)
  }
}


## Read dataset info
datasets <- read.dataset.info()

## Open figure file
pdf("intron_known_vs_novel.pdf", width=4.9, height=7.35, pointsize=6)
par(mfrow=c(3,2), lwd=0.6, cex=1, cex.main=1, font.main=1, cex.axis=5/6)

## Process each dataset
for(dataset in datasets$id) {
  make.roc.like.plot(data.dir, dataset, legend=(dataset==datasets$id[nrow(datasets)]))
}
dev.off()
