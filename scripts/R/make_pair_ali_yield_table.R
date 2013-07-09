## Get table of alignment yield for read pairs
## Input: files produced by validate_bam.pl

source(file.path(Sys.getenv("SPALE_HOME"), "scripts", "R", "spale_func.R"))
source(file.path(Sys.getenv("SPALE_HOME"), "scripts", "R", "spale_valbam_func.R"))

in.dir <- get.spale.data.fn("output", "validate_bam", "std")
out.fn <- "pair_ali_counts.txt"

read.counts <- function(data.set, paths, method.info) {
  ## Get counts
  x <- read.pair.ali.counts(data.set, paths)
  ## Match to methods
  team <- basename(dirname(colnames(x)))
  method.i <- substring(colnames(x), nchar(colnames(x)))
  j <- match(paste(method.info$team, method.info$number), paste(team, method.i))
  x <- x[, j]
  colnames(x) <- method.info$name  
  ## Return result
  return(x)
}


datasets <- read.dataset.info()
method.info <- read.method.info(include.exploratory=FALSE)

paths <- file.path(in.dir, unique(method.info$team))

out.file <- file(out.fn, "w")

for(i in 1:nrow(datasets)) {
  x <- read.counts(datasets$id[i], paths, method.info)
  x <- t(x) / 100
  x <- cbind(x, total = rowSums(x), reads = rowSums(x[, 1:3]) + rowSums(x[, 4:5]/2))
  cat("\n", datasets$name[i], "\n", file=out.file, sep="")
  write.table(x, file=out.file, sep="\t", quote=FALSE, col.names=NA)
}

close(out.file)
