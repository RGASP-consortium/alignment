## Get tables of multi-intron precision and recall from sim accuracy stats

source(file.path(Sys.getenv("SPALE_HOME"), "scripts", "R", "spale_func.R"))

data.dir <- get.spale.data.fn("output", "sim_accuracy")
col.uni <- c("uni.mi.read.actual", "uni.mi.read.reported", "uni.mi.read.correct",
             "uni.mi.pair.actual", "uni.mi.pair.reported", "uni.mi.pair.correct")
col.amb <- c("amb.mi.read.actual", "amb.mi.read.reported", "amb.mi.read.correct",
             "amb.mi.pair.actual", "amb.mi.pair.reported", "amb.mi.pair.correct")
col.names <-  c("actual.read", "reported.read", "correct.read",
                "actual.pair", "reported.pair", "correct.pair")

## Function to read multi-intron stats
read.nr.introns <- function(fn, n=11) {
  ## Read data file
  x <- read.delim(fn)

  ## Check that intron count column has values 1 ... nrow(x)
  stopifnot(all(x[, 1] == 1:nrow(x)-1))

  ## Combine results for unique and ambiguous mappings
  x <- x[, col.uni] + x[, col.amb]
  colnames(x) <- col.names

  ## Limit to n rows
  if(nrow(x) > n) x[n, ] <- colSums(x[n:nrow(x), ])
  x <- x[1:n, ]
  x[ is.na(x) ] <- 0

  ## Set row names
  rownames(x) <- 1:n-1

  ## Return result
  return(as.matrix(x))
}


## Function to compute the reverse cumulative sum for each column of a matrix
rcs <- function(x) apply(x, 2, function(x) rev(cumsum(rev(x))))


## Function to construct and write a result table
write.result.table <- function(nr.introns, out.file) {
  
  ## Consistency check
  x <- nr.introns[, "actual", ]
  stopifnot(x == x[,1])
  
  ## Compute precision and recall
  prec <- rcs(nr.introns[-1, "correct", ]) / rcs(nr.introns[-1, "reported", ])
  recall <- rcs(nr.introns[-1, "correct", ]) / rcs(nr.introns[-1, "actual", ])
  
  ## Write table
  out.file <- file(out.file, "w")
  cat("Precision\n", file=out.file)
  write.table(t(prec), file=out.file, col.names=NA, quote=FALSE, sep="\t")
  cat("Recall\n", file=out.file)
  write.table(t(recall), file=out.file, col.names=NA, quote=FALSE, sep="\t")
  cat("Simulated reads", nr.introns[-1, "actual", 1], file=out.file, sep="\t") 
  cat("\n", file=out.file)
  close(out.file)
}


## Read method info
method.info <- read.method.info(include.exploratory=FALSE)
method.info <- method.info[ method.info$name != "Truth", ] # We have no stats for Truth

## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Process each dataset
for(dataset in datasets$id) {
  ## Read intron counts
  in.files <- file.path(data.dir, paste("sa_intron_multi.", method.info$team, ".", dataset, "_", method.info$number, ".txt", sep=""))
  names(in.files) <- method.info$name
  nr.introns <- sapply(in.files, read.nr.introns, simplify="array")
  
  ## Separate read and pair data
  counts.read <- nr.introns[, c("actual.read", "reported.read", "correct.read"), ]
  counts.pair <- nr.introns[, c("actual.pair", "reported.pair", "correct.pair"), ]
  colnames(counts.read) <- colnames(counts.pair) <- c("actual", "reported", "correct")
  
  ## Make result tables
  write.result.table(counts.read, paste("multi-intron_read_accuracy_", dataset, ".txt", sep=""))
  write.result.table(counts.pair, paste("multi-intron_pair_accuracy_", dataset, ".txt", sep=""))
}
