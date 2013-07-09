get.spale.data.fn <- function(...) {
  base.dir <- Sys.getenv("SPALE_DATA")
  if(base.dir == "") stop("Environment variable SPALE_DATA undefined")
  file.path(base.dir, ...)
}


read.dataset.info <- function(fn = get.spale.data.fn("reads", "datasets.txt"),
                              id.pattern = NULL) {
  x <- read.delim(fn, quote="", as.is=TRUE)
  if(!is.null(id.pattern)) x <- x[ grepl(id.pattern, x), ]
  return(x)
}


## deprecate this:
read.total.frag.count <- function(...) {
  warn("deprecated!")
  x <- read.dataset.info(...)
  counts <- x$fragments
  names(counts) <- x$id
  return(counts)
}


## deprecate this:
read.dataset.ids <- function(...) {
  warn("deprecated!")
  x <- read.dataset.info(...)
  ids <- x$id
  names(ids) <- x$name
  return(ids)
}


read.method.info <- function(method.fn = get.spale.data.fn("aligners", "methods.txt"),
                             style.fn = get.spale.data.fn("aligners", "styles.txt"),
                             include.exploratory=TRUE )
{

  ## Read tables
  method.info <- read.delim(method.fn, quote="", as.is=TRUE)
  styles <- read.delim(style.fn)

  ## Exclude exploratory?
  if(!include.exploratory)
    method.info <- method.info[ method.info$style < 5, ]
  
  ## Add an order column
  method.info$order <- 1:nrow(method.info)
  
  ## Set style
  j <- match(method.info$style, styles$style)
  stopifnot( all(!is.na(j)) )
  method.info$lty <- styles$lty[j]
  method.info$pch <- styles$pch[j]

  return(method.info)
}


match.methods.to.filenames <- function(method.info, files) {
  ## Parse filenames
  parts <- strsplit(files,".",fixed=T)
  team <- sapply(parts, function(x) x[length(x)-2])
  i <- as.integer( sapply(parts, function(x) substring(x[length(x)-1], nchar(x[length(x)-1]))) )
  
  ## Match filenames to method info
  j <- match(paste(method.info$team, method.info$number), paste(team, i))

  return(j)
}


get.method.info.for.files <- function(files, order.methods=FALSE, ...) { 

  method.info <- read.method.info(...)

  ## Parse filenames
  parts <- strsplit(files,".",fixed=T)
  team <- sapply(parts, function(x) x[length(x)-2])
  i <- as.integer( sapply(parts, function(x) substring(x[length(x)-1], nchar(x[length(x)-1]))) )

  ## Match filenames to method info
  j <- match(paste(team, i), paste(method.info$team, method.info$number))
  stopifnot( all(!is.na(j)) )
  method.info <- method.info[j,]
  method.info$file <- files
  row.names(method.info) <- files
  if(order.methods) method.info <- method.info[order(method.info$order), ]
  
  return(method.info)
}


draw.ecdf.curve <- function(x, y, xmin=1, xmax=max(x), sample.x=1000, log.x=FALSE, ...) {
  y <- cumsum(y)
  if(!is.null(sample.x)) {
    fn <- stepfun(x, c(0, y))
    if(log.x) xmax <- log(xmax)
    x <- exp(seq(1, xmax, length.out=sample.x))
    y <- fn(x)
  }
  lines(x, y, ...)
}


plot.gap.size.ecdf <- function(sizes, methods, count.column, leg.pos="topleft", scale.factors=NULL,
                               xmin=1, xmax=NULL, ymax=NULL, log.x=TRUE, sample.x=1000, main="", xlab="Gap size (bp)",
                               ylab="Cumulative frequency") {

  if(is.null(xmax)) xmax <- max(sapply(sizes, function(x) { max(x$size,0) }))
  totals <- sapply(sizes, function(x) { sum(x[x$size <= xmax , count.column]) })

  if(is.null(scale.factors)) {
    plot(NULL, xlim=c(xmin, xmax), ylim=c(0, 100), log=if(log.x) "x" else "", las=1, xlab=xlab, ylab=ylab, main=main)
    for(i in 1:length(sizes)) {
      a <- sizes[[i]]
      count.sum <- sum(a[, count.column])
      if(count.sum > 0)
        draw.ecdf.curve(a$size, 100 * a[, count.column] / count.sum, xmin, xmax,
                        col=methods$color[i], lty=methods$lty[i], lwd=1, sample.x=sample.x, log.x=log.x)
    }
  } else {
    if(length(scale.factors) == 1) scale.factors <- rep(scale.factors, length(sizes))
    if(is.null(ymax)) ymax <- max(totals / scale.factors)
    plot(NULL, xlim=c(xmin, xmax), ylim=c(0, ymax), log=if(log.x) "x" else "", las=1, xlab=xlab, ylab=ylab, main=main)
    for(i in 1:length(sizes)) {
      a <- sizes[[i]]
      draw.ecdf.curve(a$size, a[,count.column] / scale.factors[i], xmin, xmax,
                      col=methods$color[i], lty=methods$lty[i], lwd=1, sample.x=sample.x, log.x=log.x)
    }
  }

  if(leg.pos != FALSE) legend(x=leg.pos, legend=methods$name, col=methods$color, lty=methods$lty, lwd=1)
}



