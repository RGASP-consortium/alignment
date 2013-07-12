## Plot splice recall for simulated data stratified by true junction coverage
library(plotrix) # For function axis.break

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))

data.dir <- get.rgasp.ali.data.fn("output", "sim_accuracy")


read.splice.stats.file <- function(method.id, in.dir, metric, mapping.cat) {

  ## Construct filename
  in.file <- paste("sa_intron_", metric, "_", mapping.cat, ".", method.id, ".txt", sep="")
  in.file <- file.path(in.dir, in.file)

  ## Read data
  x <- read.delim(in.file, row.names=1)
  stopifnot(rownames(x) == 1:nrow(x))
  stopifnot(colnames(x) == paste("X", 1:ncol(x)-1, sep=""))

  ## Return result, making sure it is a matrix
  return(as.matrix(x))
}

read.splice.stats <- function(method.id, in.dir, metric) {
  read.splice.stats.file(method.id, in.dir, metric, "uni") +
  read.splice.stats.file(method.id, in.dir, metric, "amb")
}

calc.subset.recall <- function(min.splice.pos, tp, fn) {
  i <- (0+min.splice.pos):(76-min.splice.pos)
  tp <- apply(tp[i, -1, ], 2:3, sum)
  fn <- apply(fn[i, -1, ], 2:3, sum)
  apply(tp, 2, cumsum) / apply(tp + fn, 2, cumsum)
}


read.recall <- function(dataset, data.dir, method.info) {
  ## Read data
  method.ids <- paste(method.info$team, ".", dataset, "_", method.info$number, sep="")
  names(method.ids) <- method.info$name
  tp <- sapply(method.ids, read.splice.stats, data.dir, "tp", simplify="array")
  fn <- sapply(method.ids, read.splice.stats, data.dir, "fn", simplify="array")
  
  ## Consistency check
  n.actual <- apply(tp+fn, 3, sum)
  stopifnot(n.actual[1] == n.actual)

  ## Compute and return recall
  recall <- 100 * sapply(1:30, calc.subset.recall, tp, fn, simplify="array")

  return(recall)
}


## Read method info
method.info <- read.method.info(include.exploratory=FALSE)
method.info <- method.info[ method.info$name != "Truth", ] # We have no stats for Truth

## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Open plot file
pdf("sim_splice_recall_by_coverage.pdf", width=6.69, height=6.69, pointsize=6)
par(mfrow=c(2,2), lwd=0.6, cex=1, cex.main=1.16, font.main=1, cex.axis=5/6)

## Process each dataset
for(dataset in datasets$id) {

  recall <- read.recall(dataset, data.dir, method.info)
  stopifnot(nrow(recall) == 20)
  
  ## Make plot
  xlim <- c(1, 20)
  ylim <- c(0, 100)
  
  plot(0, xlim=xlim, ylim=ylim, type="n", xlab="True read coverage", ylab="Recall (cumulative)", main=dataset, xaxt="n")
  axis(1, c(1,5,10,15,20), c(1,5,10,15,"20+"))
  axis.break(1, 19.5)
  for(i in 1:ncol(recall)) {
    lines(1:20, recall[, i, 1], col=method.info$color[i], lty=method.info$lty[i], lwd=1)
  }

}
  
par(xpd=NA)
plot.new()
legend("topleft", legend=method.info$name, col=method.info$color, lty=method.info$lty, lwd=1)
par(xpd=FALSE)

dev.off()
