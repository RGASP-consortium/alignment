## Plot barcharts showing accuracy of read placement

library(RColorBrewer)

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))

in.dir <- get.rgasp.ali.data.fn("output", "sim_accuracy")

bar.cols <- c(perfect = "forestgreen",
              some.right = "lightgreen",
              elsewhere = "red3",
              overlap = "indianred1")

bar.desc <- c("Perfectly mapped", "Part correctly mapped", "Mapped, no base correct",
              "No base correctly mapped but intersecting correct location")

read.stats <- function(fn) as.matrix(read.delim(fn, row.names=1))


make.read.placement.plot <- function(stats, ...) {

  ## Get total reads sequenced
  total.reads <- stats["total.reads", ]
  stopifnot(total.reads == total.reads[1]) # Check that it is the same for all methods
  total.reads <- total.reads[1]

  ## Compute read placement percentages
  x <- stats[c("perfect", "some.right", "overlap"), ]
  x <- rbind(x, elsewhere = stats[ "aln", ] - colSums(x))
  x <- 100 * x / total.reads

  ## Plot
  x <- x[, ncol(x):1]
  i <- c("perfect", "some.right")
  barplot(x[i, ], horiz=T, las=1, xlim=c(0,100), cex.names=0.7, col=bar.cols[i], ...)
  i <- c("elsewhere", "overlap")
  barplot(x[i, ], horiz=T, las=1, xlim=c(0,18), cex.names=0.7, col=bar.cols[i], ...)
}


## Read method info
method.info <- read.method.info(include.exploratory=FALSE)
method.info <- method.info[ method.info$name != "Truth", ] # We have no stats for Truth

## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Open plot file
pdf("sim_read_placement.pdf", width=14, height=8)
par(mfrow=c(2,4),  mar=c(5, 8, 4, 2)+.1)

## Process each dataset
for(dataset in datasets$id) {

  ## Read data
  in.files <- file.path(in.dir, paste("sa_general.", method.info$team, ".", dataset, "_", method.info$number, ".txt", sep=""))
  names(in.files) <- method.info$name
  stats <- sapply(in.files, read.stats, simplify="array")

  ## Make plots
  make.read.placement.plot(stats["uni_contig",,] + stats["amb_contig",,], main=dataset, xlab="Percent of simulated unspliced reads")
  make.read.placement.plot(stats["uni_spliced",,] + stats["amb_spliced",,], main=dataset, xlab="Percent of simulated spliced reads")

}

## Plot legend
plot(0, type="n", ann=FALSE, axes=FALSE)
legend(x="topleft", fill=bar.cols, xpd=NA, legend=bar.desc)

## Close plot file
dev.off()

