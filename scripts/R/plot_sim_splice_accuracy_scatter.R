## Scatter plots of splice accuracy for simulated data

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))

data.dir <- get.rgasp.ali.data.fn("output", "sim_accuracy")
glyph.cex <- 1.8
xlim <- c(0, 4); ylim <- c(60, 100)

read.splice.stats.file <- function(method.id, in.dir, metric, mapping.cat) {

  ## Construct filename
  in.file <- paste("sa_intron_", metric, "_", mapping.cat, ".", method.id, ".txt", sep="")
  in.file <- file.path(in.dir, in.file)

  ## Read data
  x <- read.delim(in.file, row.names=1)
  stopifnot(rownames(x) == 1:nrow(x))
  stopifnot(colnames(x) == paste("X", 1:ncol(x)-1, sep=""))

  return(rowSums(x))
}

read.splice.stats <- function(method.id, in.dir, metric) {
  read.splice.stats.file(method.id, in.dir, metric, "uni") +
  read.splice.stats.file(method.id, in.dir, metric, "amb")
}

calc.subset.recall <- function(i, tp, fn) {
  i <- (0+i):(76-i)
  colSums(tp[i,]) / ( colSums(tp[i,]) + colSums(fn[i,]) )
}

calc.fdr <- function(i, fp, tp) {
  i <- (0+i):(76-i)
  colSums(fp[i,]) / ( colSums(fp[i,]) + colSums(tp[i,]) )
}


read.fdr.recall <- function(dataset, data.dir, method.info) {
  ## Read data
  method.ids <- paste(method.info$team, ".", dataset, "_", method.info$number, sep="")
  names(method.ids) <- method.info$name
  tp <- sapply(method.ids, read.splice.stats, data.dir, "tp")
  fp <- sapply(method.ids, read.splice.stats, data.dir, "fp")
  fn <- sapply(method.ids, read.splice.stats, data.dir, "fn")
  
  ## Consistency check
  n.actual <- colSums(tp) + colSums(fn)
  stopifnot(n.actual[1] == n.actual)

  ## Compute recall and FDR, and return
  recall <- 100 * sapply(1:30, calc.subset.recall, tp, fn)
  fdr <- 100 * sapply(1:30, calc.fdr, fp, tp)
  return(list(recall=recall, fdr=fdr))
}

plot.fdr.recall <- function(fdr, recall, method.info, ...)
  plot(fdr, recall, pch=method.info$pch, lty=method.info$lty, bg=method.info$color, las=1,
       xlab="Splice FDR (%)", ylab="Splice recall (%)", cex=glyph.cex, ...)

## Read method info
method.info <- read.method.info(include.exploratory=FALSE)
method.info <- method.info[ method.info$name != "Truth", ] # We have no stats for Truth

## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Process each dataset
for(dataset in datasets$id) {

  stats <- read.fdr.recall(dataset, data.dir, method.info)

  out.fn <- paste("splice_accuracy_scatter_", dataset, ".pdf", sep="")
  pdf(out.fn, width=6, height=9, pointsize=6)
  par(mfrow=c(3,2), lwd=0.6, cex=1, cex.main=1, font.main=1, cex.axis=5/6)

  for(min.pos in c(1,20)) {
    title.text <- paste(dataset, " min.pos=", min.pos, sep="")
    plot.fdr.recall(stats$fdr[, min.pos], stats$recall[, min.pos], ylim=c(0, 100),
                    method.info, main=title.text)
    rect(xlim[1], ylim[1], xlim[2], ylim[2], border="grey")
    plot.fdr.recall(stats$fdr[, min.pos], stats$recall[, min.pos], xlim=xlim, ylim=ylim,
                    method.info, main=title.text)
  }   

  plot.new()
  legend(x="topleft", legend=method.info$name, pt.bg=method.info$col, pch=method.info$pch, pt.cex=glyph.cex, y.intersp=1.6, xpd=NA)

  dev.off()
}
