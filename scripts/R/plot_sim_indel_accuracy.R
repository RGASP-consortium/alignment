## Plot indel accuracy stats for simulated data

library(gplots)
library(RColorBrewer)

source(file.path(Sys.getenv("RGASP_ALI_HOME"), "scripts", "R", "general_functions.R"))

data.dir <- get.rgasp.ali.data.fn("output", "sim_accuracy")
hm.cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(255)
hm.cols <- rev(hm.cols)
hm.breaks <- seq(0, 100, length=256)


rcs <- function(x) rev(cumsum(rev(x)))


calc.accuracy <- function(correct, true, reported) {
  recall <- correct / true
  prec <- correct / reported
  recall.cumul <- rcs(correct) / rcs(true)
  prec.cumul <- rcs(correct) / rcs(reported)
  f <- 2 * recall.cumul * prec.cumul / (recall.cumul + prec.cumul)
  cbind(recall, prec, recall.cumul, prec.cumul, f)
}


read.indel.stats <- function(fn, n=8) {

  ## Columns to extract deletion data from
  del.correct.col <- c("uni.del.correct", "amb.del.correct")
  del.false.col <- c("uni.del.false", "amb.del.false")
  del.true.col <- c(del.correct.col, "uni.del.missed", "amb.del.missed")
  del.reported.col <- c(del.correct.col, del.false.col)

  ## Columns to extract insertion data from
  ins.correct.col <- c("uni.ins.correct", "amb.ins.correct")
  ins.false.col <- c("uni.ins.false", "amb.ins.false")
  ins.true.col <- c(ins.correct.col, "uni.ins.missed", "amb.ins.missed")
  ins.reported.col <- c(ins.correct.col, ins.false.col)

  ## Read data
  x <- read.delim(fn)
  stopifnot(all(x[, 1] == 1:nrow(x)-1))

  ## Discard counts for indel size 0
  x <- x[-1, ]

  ## Set number of rows
  if(nrow(x) > n) x[n, ] <- colSums(x[n:nrow(x), ])
  x <- x[1:n, ]
  x[ is.na(x) ] <- 0
  
  ## Compute precision, recall and F-score
  del.correct <- rowSums(x[, del.correct.col])
  del.true <- rowSums(x[, del.true.col])
  del.reported <- rowSums(x[, del.reported.col])
  ins.correct <- rowSums(x[, ins.correct.col])
  ins.true <- rowSums(x[, ins.true.col])
  ins.reported <- rowSums(x[, ins.reported.col])
  del.acc <- calc.accuracy(del.correct, del.true, del.reported)
  ins.acc <- calc.accuracy(ins.correct, ins.true, ins.reported)
  indel.acc <- calc.accuracy(del.correct + ins.correct, del.true + ins.true, del.reported + ins.reported)
  
  ## Assemble into one matrix
  colnames(del.acc) <- paste("del", colnames(del.acc), sep=".")
  colnames(ins.acc) <- paste("ins", colnames(ins.acc), sep=".")
  colnames(indel.acc) <- paste("indel", colnames(indel.acc), sep=".")
  x <- cbind(del.acc, ins.acc, indel.acc)
  rownames(x) <- 1:n

  ## Return result, making sure it is a matrix
  return(as.matrix(x))
}


## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Process each dataset
for(dataset in datasets$id) {
  
  ## Get filenames and method info
  files <- dir(path=data.dir, pattern=glob2rx(paste("sa_indel.*.",dataset,"_?.txt",sep="")))
  method.info <- get.method.info.for.files(files, order.methods=TRUE)
  method.info <- method.info[method.info$style < 5, ] # Exclude exploratory results
  rm(files)
  
  ## Read data
  x <- sapply(file.path(data.dir, method.info$file), read.indel.stats, simplify="array")
  dimnames(x)[[3]] <- method.info$name
  
  save(x, file=paste("indel_accuracy_", dataset, ".RData", sep=""))
  
  ## Plot
  pdf(paste("indel_accuracy_", dataset, ".pdf", sep=""))
  par(oma=c(0,0,0,3))
  for(i in c("del.recall", "del.prec", "ins.recall", "ins.prec"))
    heatmap.2(t(100*x[, i, ]), col=hm.cols, breaks=hm.breaks,
              Colv=NULL, Rowv=NULL, dendrogram="none", density.info="none", trace="none",
              main=paste(dataset, i), na.color="red")
  dev.off()
  
  ## Print F-scores
  cat(dataset, "del.f\n")
  print(round(100*t(x[,"del.f",])))
  cat(dataset, "ins.f\n")
  print(round(100*t(x[,"ins.f",])))
  
}


quit()

## Stats can be extracted from the .RData files gerated by this script as follows:

load("indel_accuracy_sim1.RData")
x[1, "indel.prec.cumul", ]
apply(x[,"del.recall",], 2, min)
sort(x[5, "ins.recall.cumul", ])

load("indel_accuracy_sim2.RData")
(x[1, "indel.prec.cumul", ])
apply(x[,"del.recall",], 2, min)
sort(x[5, "ins.recall.cumul", ])
