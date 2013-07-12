## Make tables with fundamental simulation accuracy stats

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))

in.dir <- get.spale.data.fn("output", "sim_accuracy")


read.stats <- function(fn) {
  x <- read.delim(fn, row.names=1)
  x <- as.matrix(x)
  return(x)
}


write.accuracy.table <- function(stats.uni, stats.amb, out.fn) {

  stats.all <- stats.uni + stats.amb
  rm(stats.amb)

  ## Get total reads sequenced
  total.reads <- stats.all["total.reads", ]
  stopifnot(total.reads == total.reads[1]) # Check that it is the same for all methods
  total.reads <- total.reads[1]
  
  ## Get total bases sequenced
  total.bases <- stats.all["total.bases",]
  stopifnot(total.bases == total.bases[1]) # Check that it is the same for all methods
  
  ## Compute fractions, and combine into a result table
  i <- c("aln", "perfect", "some.right", "overlap")
  j <- c("bases.right", "bases.wrong")
  x <- rbind(stats.uni[i, ] / total.reads,
             stats.uni[j, ] / total.bases,
             stats.all[i, ] / total.reads,
             stats.all[j, ] / total.bases)
  rownames(x) <- paste(rep(c("uni", "all"), each=6), rownames(x)) 

  ## Write table to file
  write.table(t(x), file=out.fn, sep="\t", quote=FALSE, col.names=NA)

} 


write.ambig.table <- function(stats.uni, stats.amb, out.fn) {

  ## Get total reads sequenced
  stats.all <- stats.uni + stats.amb
  total.reads <- stats.all["total.reads", ]
  stopifnot(total.reads == total.reads[1]) # Check that it is the same for all methods
  total.reads <- total.reads[1]
  
  ## Compute fractions, and combine into a result table
  i <- c("aln", "perfect")
  j <- c("bases.right")
  x <- rbind("uni.mapped" = stats.uni["aln", ] / total.reads,
             "amb.mapped" = stats.amb["aln", ] / total.reads,
             "uni.perfect" = stats.uni["perfect", ] / stats.uni["aln", ],
             "amb.perfect" = stats.amb["perfect", ] / stats.amb["aln", ],
             "uni.bases.right" = stats.uni["bases.right", ] / (stats.uni["bases.right", ] + stats.uni["bases.wrong", ]),
             "amb.bases.right" = stats.amb["bases.right", ] / (stats.amb["bases.right", ] + stats.amb["bases.wrong", ]))

  ## Write table to file
  write.table(t(x), file=out.fn, sep="\t", quote=FALSE, col.names=NA)
} 


## Read method info
method.info <- read.method.info(include.exploratory=FALSE)
method.info <- method.info[ method.info$name != "Truth", ] # We have no stats for Truth

## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Process each dataset
for(dataset in datasets$id) {

  ## Read data
  in.files <- file.path(in.dir, paste("sa_general.", method.info$team, ".", dataset, "_", method.info$number, ".txt", sep=""))
  names(in.files) <- method.info$name
  stats <- sapply(in.files, read.stats, simplify="array")
  
  ## Write result tables
  write.accuracy.table(stats["uni_contig",,] + stats["uni_spliced",,],
                     stats["amb_contig",,] + stats["amb_spliced",,],
                     paste("accuracy_overall_summary_", dataset, ".txt", sep=""))
  write.accuracy.table(stats["uni_contig",,], stats["amb_contig",,],
                     paste("accuracy_unspliced_summary_", dataset, ".txt", sep=""))
  write.accuracy.table(stats["uni_spliced",,], stats["amb_spliced",,],
                     paste("accuracy_spliced_summary_", dataset, ".txt", sep=""))
  write.ambig.table(stats["uni_contig",,] + stats["uni_spliced",,],
                    stats["amb_contig",,] + stats["amb_spliced",,],
                    paste("ambig_overall_summary_", dataset, ".txt", sep=""))
}
