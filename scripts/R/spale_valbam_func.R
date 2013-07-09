## Functions for handling output from validate_bam.pl

read.tabular.matrix <- function(fn, data.column) {
  x <- read.delim(fn)
  n <- max(x[,1]) + 1
  counts <- 0:(n-1)
  stopifnot(data.column > 2)
  stopifnot(ncol(x) >= data.column)
  stopifnot(nrow(x) == n*n)
  stopifnot(x[,1] == rep(counts, each=n))
  stopifnot(x[,2] == rep(counts, n))
  matrix(x[, data.column], nrow=n, byrow=T, dimnames=list(counts, counts))
}


read.count.table <- function(pattern, path=".") {
  files <- dir(path=path, pattern=glob2rx(pattern))
  if(length(files) == 0) return(NULL)
  all.list <- lapply(paste(path, files, sep="/"), read.delim, header=FALSE, row.names=1, quote="")
  names(all.list) <- sapply(strsplit(files, ".", fixed=TRUE), function(x) x[1])
  metrics <- unique(c(unlist(sapply(all.list, row.names))))
  all.mat <- matrix(NA, nrow=length(all.list), ncol=length(metrics), dimnames=list(names(all.list), metrics))
  for(i in 1:length(all.list)) all.mat[i, rownames(all.list[[i]])] <- all.list[[i]][,1]
  return(all.mat)
}


read.count.array <- function(pattern, path=".", tabular=FALSE, data.column=3) {
  files <- dir(path=path, pattern=glob2rx(pattern), full.names=TRUE)
  if(length(files) == 0) return(NULL)
  if(tabular)
    all.list <- lapply(files, read.tabular.matrix, data.column=data.column)
  else
    all.list <- lapply(files, read.delim, row.names=1, quote="", check.names=FALSE, as.is=TRUE)
  ids <- list(sapply(strsplit(files, ".", fixed=TRUE), function(x) x[1]),
              unique(c(unlist(sapply(all.list, rownames)))),
              unique(c(unlist(sapply(all.list, colnames)))))
  all.mat <- array(NA, dim=sapply(ids, length), dimnames=ids)
  for(i in 1:length(all.list)) all.mat[i, rownames(all.list[[i]]), colnames(all.list[[i]])] <- as.matrix(all.list[[i]])
  return(all.mat)
}


read.bool.count.table <- function(pattern, path=".") {
  files <- unique(unlist(sapply(pattern, function(x) dir(path=path, pattern=glob2rx(x)))))
  if(length(files) == 0) return(NULL)
  all.list <- lapply(paste(path, files, sep="/"), read.delim, quote="", as.is=TRUE)
  names(all.list) <- sapply(strsplit(files, ".", fixed=TRUE), function(x) x[1])
  ids <- unique(c(unlist(sapply(all.list, function(x) x$id))))
  all.mat <- matrix(NA, nrow=length(all.list), ncol=length(ids), dimnames=list(names(all.list), ids))
  for(i in 1:length(all.list)) all.mat[i, all.list[[i]]$id] <- all.list[[i]]$count
  invalid.ids <- unique(c(unlist(sapply(all.list, function(x) x$id[!x$valid]))))
  return(list(counts=all.mat, invalid.ids=invalid.ids))
}


read.pair.ali.counts <- function(data.set, paths) {
  x <- read.count.array(paste(data.set, "_*.aln_counts.txt", sep=""), paths, tabular=TRUE, data.column=3)
  total <- sum(x[1, , ])
  stopifnot(all(apply(x, 1, sum) == total))
  y <- matrix(NA, nrow=dim(x)[1], ncol=5)
  rownames(y) <- dimnames(x)[[1]]
  colnames(y) <- c("Both mates uniquely mapped", "Both mates multi-mapped", "One unique & one multi", "One unique & one unaligned", "One multi & one unaligned")
  for(j in 1:dim(x)[1]) {
    y[j, 1] <- x[j, 2, 2]                             # UU
    y[j, 2] <- sum(x[j, -(1:2), -(1:2)])              # MM
    y[j, 3] <- sum(x[j, 2, -(1:2)], x[j, -(1:2), 2])  # UM
    y[j, 4] <- sum(x[j, 1, 2], x[j, 2, 1])            # UN
    y[j, 5] <- sum(x[j, 1, -(1:2)], x[j, -(1:2), 1])  # MN
    stopifnot(sum(y[j, ]) + x[j, 1, 1] == total)      # Paranoia
  }
  y <- 100 * t(y / total)
  y <- y[, ncol(y):1]
  return(y)
}
