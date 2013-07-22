## Plot positional distribution of gaps and mismatches in read sequences

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))
data.dir <- get.rgasp.ali.data.fn("output", "gap_mm_pos", "primary")

count.types <- c(Mismatches = "mm", Insertions = "ins", Deletions = "del", Introns = "intron")


## Function to compute coefficient of variation
cv <- function(x) {
   pop.sd <- sqrt( var(x) * (length(x)-1) / length(x) )
   pop.sd / mean(x)
}


## Function that finds the index ranges for which the input vector (x) exceeds a specified cutoff (max.x)
## Returns a list of ranges (coordinate pairs)
get.saturated.ranges <- function(x, max.x) {
  i <- x > max.x # Find which values in x are above cutoff
  if(all(!i)) return(list()) # If none, return emtpty list
  i.rle <- rle(i) # Convert to run-length encoding (RLE), where each element represents a segment of x that is entirely above or below cutoff
  i.rle$values <- 1:length(i.rle$values) # Give each segment a distinct value
  tapply(which(i), inverse.rle(i.rle)[i], range) # Compute range for each segement above the cutoff
}


## Function that plots positional distributions
plot.pos.curves <- function(counts, dataset, method.info, qual.scores, max.y=5, abs.val=FALSE) {

  par(mfrow=c(length(counts)+1, 5), mar=c(0.5,0.5,0.5,0.5))
  plot.new()

  legend("topleft", dataset, bty="n", cex=0.8)
  for(i in 1:4) {
    plot.new()
    legend("topleft", names(count.types)[i], bty="n")
  } 

  for (i in 1:length(counts)) {
    plot.new()
    legend("topleft", names(counts)[[i]], bty="n")
    for(j in count.types) {
      y <- counts[[i]][, j]
      stopifnot(length(y) == 76)
      if(j %in% c("del", "intron")) {
        stopifnot(y[76] == 0)
        y <- y[1:75]
      }
      if(abs.val) {
        y <- y/1000
        plot(y, type="l", ann=F, axes=F, ylim=trunc(range(y)))
        axis(2, trunc(range(y)), las=1)      
      } else {
        y <- 100 * y / sum(y)
        y[is.na(y)] <- 0
        sat.ranges.list <- get.saturated.ranges(y, max.y)
        y[y > max.y] <- max.y
        plot(y, type="l", ann=F, axes=F, ylim=c(0, max.y))
        for(sat.range in sat.ranges.list) 
          lines(sat.range, c(max.y, max.y), col="red")
        if(j == count.types[1]) axis(2, c(0, max.y), las=1)
      }
    }   
  }

}


## Read dataset info
datasets <- read.dataset.info()

## Open plot files
pdf("indel_mm_pos_distrib_percent.pdf", width=7, height=12)
pdf("indel_mm_pos_distrib_absolute.pdf", width=7, height=12)
pdf("indel_pos_cv.pdf")

## Read and plot results for each data set
for(i in 1:nrow(datasets)) {

  ## Read data
  files <- dir(path=data.dir, pattern=glob2rx(paste("imp.*", datasets$id[i], "_*.txt",sep="")))
  methods <- get.method.info.for.files(files, order.methods=TRUE)
  methods <- methods[ methods$style < 5, ]  # Exclude exploratory methods
  files.full.path <- file.path(data.dir, methods$file)
  x <- lapply(files.full.path, read.delim)
  names(x) <- methods$name

  ## Compute coefficient of variations for indel positional distributions
  indel.cv <- sort(sapply(x, function(x) cv(x$del + x$ins)))
  
  ## Make plots
  ## Plot positional distributions (percentages)
  dev.set()
  plot.pos.curves(x, datasets$name[i], method.info, qual.scores)
  ## Plot positional distributions (absolute counts)
  dev.set()
  plot.pos.curves(x, datasets$name[i], method.info, qual.scores, abs.val=TRUE)
  ## Make bar plot showing indel CVs
  dev.set()
  par(mar=c(5, 8, 4, 2))
  barplot(indel.cv, main=datasets$name[i], horiz=TRUE, las=1, xlab="Coefficient of variation")
}

## Close plot files
dev.off()
dev.off()
dev.off()

## Done
sessionInfo()
