## Plot cufflinks sim accuracy as abacus plot, and make corresponding tables
## This script uses the output from cufflinks_eval.pl

library(RColorBrewer)
abacus.col <- brewer.pal(4, "Set1")[-1]

source(file.path(Sys.getenv("SPALE_HOME"), "scripts", "R", "spale_func.R"))

sim.dir <- get.spale.data.fn("sim_logs")
stats.dir <- get.spale.data.fn("output", "cufflinks", "eval")

read.tx.stats <- function(strand.mode, cufflinks.ver, aln.set, data.set, method.info) {

  ## Read data
  stats.fn <- file.path(stats.dir, paste("ce", strand.mode, sep="."), paste("ce", data.set, cufflinks.ver, aln.set, "tx_match.txt", sep="."))
  stats <- read.delim(stats.fn, row.names=1)
  stats <- as.matrix(stats)

  ## Sort by transcript ID
  stats <- stats[ order(as.numeric(rownames(stats))), ]

  ## Set method names
  x <- strsplit(colnames(stats), ".", fixed=TRUE)
  team <- sapply(x, function(x) x[length(x)-1])
  i <- substring(colnames(stats), nchar(colnames(stats)))
  i <- match(paste(method.info$team, method.info$number), paste(team, i))
  method.info <- method.info[ !is.na(i), ]
  stats <- stats[, i[!is.na(i)] ]
  colnames(stats) <- ifelse(is.na(method.info$name), colnames(stats), method.info$name)

  return(stats)
}


read.tx.variants <- function(strand.mode, cufflinks.ver, aln.set, data.set) {
  data.fn <- file.path(stats.dir, paste("ce", strand.mode, sep="."), paste("ce", data.set, cufflinks.ver, aln.set, "tx_groups.txt", sep="."))
  x <- read.delim(data.fn, quote="", as.is=TRUE, header=FALSE)
  stopifnot(ncol(x) == 2)
  colnames(x) <- c("tx", "variant")
  x <- x[order(x$variant), ]
  rownames(x) <- NULL
  return(x)
}


read.tx.known <- function(strand.mode, cufflinks.ver, aln.set, data.set) {
  data.fn <- file.path(stats.dir, paste("ce", strand.mode, sep="."), paste("ce", data.set, cufflinks.ver, aln.set, "tx_known.txt", sep="."))
  x <- read.delim(data.fn)
  stopifnot(ncol(x) == 2)
  stopifnot(x$id == 1:nrow(x))
  return(x$known == 1)
}


read.general.stats <- function(strand.mode, cufflinks.ver, aln.set, data.set, method.info) {

  ## Read data
  stats.fn <- file.path(stats.dir, paste("ce", strand.mode, sep="."), paste("ce", data.set, cufflinks.ver, aln.set, "stats.txt", sep="."))
  stats <- read.delim(stats.fn, quote="", as.is=TRUE)

  ## Set method names
  x <- strsplit(stats$submission, "/", fixed=TRUE)
  team <- sapply(x, function(x) x[length(x)-1])
  i <- substring(stats$submission, nchar(stats$submission))
  i <- match(paste(method.info$team, method.info$number), paste(team, i))
  method.info <- method.info[ !is.na(i), ]
  stats <- stats[i[!is.na(i)], ]
  rownames(stats) <- method.info$name

  return(stats)
}


## Compute sensitivity for groups of transcripts
calc.sens.stratified <- function(stats, groups) {
  ## Split into bins
  stats.split <- split.data.frame(stats, groups)
  ## Compute sensitivity for each bin
  sens <- sapply(stats.split, function(x) colSums(x) / nrow(x))
  ## Return result
  return(sens)
}


compute.prec <- function(general.stats, tx.stats, tx.known) {

  ## Get numbers for computing overall precision
  n.pred <- general.stats$pred.tx
  n.pred.true <- general.stats$pred.tx.corr
  
  ## Get numbers for computing precision specifically for known and novel transcripts
  ## We get the number of correctly predicted reference transcripts in known and novel categories
  ## However what we really need is the number of predictions matching known and novel transcripts
  ## For stranded comparisons this should be the same, but not necessarily for unstranded comparisons
  ## So we check that the totals match, and if they don't, we report NA
  tx.stats <- tx.stats == 3
  i <- n.pred.true == colSums(tx.stats)
  n.pred.true.known <- ifelse(i, colSums(tx.stats[tx.known, ]), NA)
  n.pred.true.novel <- ifelse(i, colSums(tx.stats[!tx.known, ]), NA)
  prec <- cbind(all = n.pred.true / n.pred,
                known = n.pred.true.known / (n.pred - n.pred.true.novel),
                novel = n.pred.true.novel / (n.pred - n.pred.true.known))
  
  return(prec)
}


write.tx.accuracy.table <- function(stats, groups, prec, out.file, ...) {

  ## Make logical matrix for sensitvity calculcations
  ## TRUE means that the transcript was correctly predicted
  stats <- stats == 3

  ## Compute stratified sensitivity
  sens.strat <- calc.sens.stratified(stats, groups)
  
  ## Put overall sensitivity and precision in the same matrix
  sens <- colSums(stats) / nrow(stats)
  x <- cbind(sens.strat, sens, prec)

  ## Output
  write.table(x, file=out.file, sep="\t", col.names=NA, ...)
}  


write.exon.accuracy.table <- function(x, out.file, ...) {

  ref.exons.corr.novel <- x$ref.exons.corr - x$ref.exons.corr.known
  ref.exons.novel <- x$ref.exons - x$ref.exons.known
  pred.exons.corr.novel <- x$pred.exons.corr - x$pred.exons.corr.known
  pred.exons.false <- x$pred.exons - x$pred.exons.corr

  result <- cbind(sens.all = x$ref.exons.corr / x$ref.exons,
                  sens.known = x$ref.exons.corr.known / x$ref.exons.known,
                  sens.novel = ref.exons.corr.novel / ref.exons.novel,
                  prec.all = x$pred.exons.corr / x$pred.exons,
                  prec.known = x$pred.exons.corr.known / (x$pred.exons.corr.known + pred.exons.false),
                  prec.novel = pred.exons.corr.novel / (pred.exons.corr.novel + pred.exons.false))
  rownames(result) <- rownames(x)
  
  ## Output
  write.table(result, file=out.file, sep="\t", col.names=NA, ...)
}  


make.abacus.plot <- function(stats, fpkm, groups, prec, plot.title) {

  ## Make logical matrix for sensitvity calculcations
  ## TRUE means that the transcript was correctly predicted
  stats <- stats == 3
  
  ## Compute stratified sensitivity
  sens <- calc.sens.stratified(stats, groups)

  ## Create a string with FPKM range for each group
  group.string <- paste(tapply(round(fpkm,2), groups, range), collapse=" ")  
  
  ## Plot
  i <- nrow(sens):1
  par(mar=c(5, 10, 4, 2))
  plot(NULL, type="n", xlim=c(0, max(sens, prec, na.rm=TRUE)), ylim=range(i), xlab=group.string, ylab="", yaxt="n", main=plot.title)
  axis(2, i, rownames(sens), las=1, tick=FALSE)
  abline(h=i, col="gray", lwd=2)
  for(j in 1:ncol(sens)) {
    points(sens[, j], i, pch=16, col=abacus.col[j], cex=3)
  }
  if(!is.null(prec)) points(prec, i, pch=5, cex=3)

  ## Report accuracy for true alignments
  cat("[ ", plot.title, " ] Truth sens = ", sens["Truth", ], "\n")
  cat("[ ", plot.title, " ] Truth prec = ", prec["Truth"], "\n")
 
}  


make.barplot <- function(stats, plot.title) {

  x <- rbind(n.corr = colSums(stats == 3),
             n.part = colSums(stats == 2),
             n.wrong = colSums(stats == 1),
             n.missed = colSums(stats == 0))
  x <- x[, ncol(x):1]
  
  barplot(x, horiz=T, las=1, cex.names=0.5, col=brewer.pal(4, "Paired"), main=plot.title)
 
}  


## Read method info
method.info <- read.method.info(include.exploratory=FALSE)

## Read dataset info
datasets <- read.dataset.info(id.pattern="^sim")

## Process each dataset
for(dataset in datasets$id) {

  expr.fn <- file.path(sim.dir, dataset, paste("expr_", dataset, ".txt", sep=""))
  all.expr <- read.delim(expr.fn)
  all.expr <- all.expr[all.expr$exons > 1 & all.expr$fpkm > 0, ]
  
  for(strand.mode in c("stranded", "unstranded")) {
  
    for (cufflinks.ver in c("v2")) {  # Version 1 is not of interest anymore
      
      for(aln.set in c("all", "primary")) {    
        
        cat("**", dataset, strand.mode, cufflinks.ver, aln.set, "**\n")
          
        ## Read data
        tx.stats <- read.tx.stats(strand.mode, cufflinks.ver, aln.set, dataset, method.info)
        tx.variants <- read.tx.variants(strand.mode, cufflinks.ver, aln.set, dataset)
        tx.known <- read.tx.known(strand.mode, cufflinks.ver, aln.set, dataset)
        general.stats <- read.general.stats(strand.mode, cufflinks.ver, aln.set, dataset, method.info)

        ## Compute precision
        prec <- compute.prec(general.stats, tx.stats, tx.known)

        ## Compute expression for each splice form
        ## For each splice form, the may be multiple variants with different starts and ends
        stopifnot(all.expr$id == tx.variants$variant)  # Check that data structures are in the same order    
        tx.fpkm <- tapply(all.expr$fpkm, tx.variants$tx, sum)

        ## Check row order
        stopifnot(rownames(tx.stats) == 1:nrow(tx.stats))
        stopifnot(names(tx.fpkm) == rownames(tx.stats))

        ## Divide transcripts into groups of low, mid and high expression
        expr.class <- cut(tx.fpkm, quantile(tx.fpkm, probs=c(0, 1/3, 2/3, 1)), include.lowest=TRUE)
        
        ## Make plots
        pdf(paste("cufflinks_abacus-plot_",strand.mode,"_",cufflinks.ver,"_", aln.set,"_",dataset,".pdf", sep=""), height=11)
        make.abacus.plot(tx.stats, tx.fpkm, expr.class, prec[,"all"], dataset)
        make.abacus.plot(tx.stats[ tx.known, ], tx.fpkm[ tx.known], expr.class[ tx.known], prec[,"known"], paste(dataset, "known"))
        make.abacus.plot(tx.stats[!tx.known, ], tx.fpkm[!tx.known], expr.class[!tx.known], prec[,"novel"], paste(dataset, "novel"))
        make.barplot(tx.stats[ tx.known, ], paste(dataset, "known"))
        make.barplot(tx.stats[!tx.known, ], paste(dataset, "novel"))
        dev.off()

        ## Write accuracy table
        out.table <- paste("cufflinks_accuracy_",strand.mode,"_",cufflinks.ver,"_", aln.set,"_",dataset,".txt", sep="")
        out.table <- file(out.table, "w")
        cat("[all transcripts]\n", file=out.table)
        write.tx.accuracy.table(tx.stats, expr.class, prec[,"all"], out.table)
        cat("[known transcripts]\n", file=out.table)
        write.tx.accuracy.table(tx.stats[tx.known, ], expr.class[tx.known], prec[,"known"], out.table)
        cat("[novel transcripts]\n", file=out.table)
        write.tx.accuracy.table(tx.stats[!tx.known, ], expr.class[!tx.known], prec[,"novel"], out.table)
        cat("[exons]\n", file=out.table)
        write.exon.accuracy.table(general.stats, out.table)
        close(out.table)

      }
    }
  }
}
