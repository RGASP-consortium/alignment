## Tabulate frequencies at which true and false introns from alignments are incorporated into transcripts by Cufflinks

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))

data.dir <- get.spale.data.fn("output", "cufflinks", "eval", "acg")

read.acg.freq <- function(fn) {
  counts <- read.delim(fn, row.names=1)
  true.freq <- counts["true.pred", ] / (counts["true.skipped", ] + counts["true.pred", ])
  false.freq <- counts["false.pred", ] / (counts["false.skipped", ] + counts["false.pred", ])
  sums <- rowSums(counts)
  true.data <- c(sums["true.pred"], sums["true.skipped"],
                 sums["true.pred"] / (sums["true.pred"] + sums["true.skipped"]),
                 true.freq)
  false.data <- c(sums["false.pred"], sums["false.skipped"],
                 sums["false.pred"] / (sums["false.pred"] + sums["false.skipped"]),
                 false.freq)
  x <- rbind(True=as.numeric(true.data), False=as.numeric(false.data))
  colnames(x) <- c("pred", "skipped", "freq.all", colnames(counts))
  return(x)
}

datasets <- read.dataset.info(id.pattern="^sim")

for(dataset in datasets$id) {

  ## Get filenames and method info
  files <- dir(data.dir, pattern=glob2rx(paste("acg.*.", dataset, "_?.txt", sep="")))
  method.info <- get.method.info.for.files(files, order.methods=TRUE)
  method.info <- method.info[method.info$style < 5, ] # Exclude exploratory results
  rm(files)

  ## Read data
  x <- sapply(file.path(data.dir, method.info$file), read.acg.freq, simplify="array")
  dimnames(x)[[3]] <- method.info$name

  ## Write to file
  out.file <- paste("aln_vs_cuff_gaps_", dataset, ".txt", sep="")
  for(i in 1:nrow(method.info)) {
    y <- data.frame(method=c(method.info$name[i], ""), intron.type=rownames(x[, , i]), x[, , i])
    write.table(y, out.file, sep="\t", col.names=(i==1), row.names=FALSE, append=(i>1), quote=FALSE)
  }
}

