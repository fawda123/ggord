#' Ordination plots with ggplot2
#'
#' Create an ordination biplot using ggplot2 including options for selecting axes, group color aesthetics, and selection of variables to plot.
#'
#' @param ord_in input ordination object
#' @param grp_in vector of grouping objects for the biplot, must have the same number of observations as the original matrix used for the ordination
#' @param obs matrix or data frame of axis scores for each observation
#' @param vecs matrix or data frame of axis scores for each variable
#' @param axes chr string indicating which axes to plot
#' @param arrow numeric indicating length of the arrow heads on the vectors
#' @param ext numeric indicating scalar distance of the labels from the arrow ends
#' @param size numeric indicating size of the observatoin points
#' @param txt numeric indicating size of the text labels for the vectors
#' @param xlims two numeric values indicating x-axis limits
#' @param ylims two numeric values indicating y-axis limits
#' @param var_sub chr string indcating which labels to show.  Regular expression matching is used.
#' @param ... arguments passed to or from other methods
#'
#' @export
#'
#' @import ggplot2 grid
#'
#' @return A \code{\link[ggplot2]{ggplot}} object that can be further modified
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#'
#' # principal components analysis with the iris data set
#' # prcomp
#' ord <- prcomp(iris[, 1:4])
#'
#' p <- ggord(ord, iris$Species)
#' p
#'
#' p + scale_colour_manual('Species', values = c('purple', 'orange', 'blue'))
#' p + theme_classic()
#' p + theme(legend.position = 'top')
#' p + scale_x_continuous(limits = c(-2, 2))
#'
#' # principal components analysis with the iris dataset
#' # princomp
#' ord <- princomp(iris[, 1:4])
#'
#' ggord(ord, iris$Species)
#'
#' # principal components analysis with the iris dataset
#' # PCA
#' library(FactoMineR)
#'
#' ord <- PCA(iris[, 1:4], graph = FALSE)
#'
#' ggord(ord, iris$Species)
#'
#' # principal components analysis with the iris dataset
#' # dudi.pca
#' library(ade4)
#'
#' ord <- dudi.pca(iris[, 1:4], scannf = FALSE, nf = 4)
#'
#' ggord(ord, iris$Species)
#'
#' # multiple correspondence analysis with the tea dataset
#' # MCA
#' data(tea)
#' tea <- tea[, c('Tea', 'sugar', 'price', 'age_Q', 'sex')]
#'
#' ord <- MCA(tea[, -1], graph = FALSE)
#'
#' ggord(ord, tea$Tea)
#'
#' # multiple correspondence analysis with the tea dataset
#' # mca
#' library(MASS)
#'
#' ord <- mca(tea[, -1])
#'
#' ggord(ord, tea$Tea)
#'
#' # multiple correspondence analysis with the tea dataset
#' # acm
#' ord <- dudi.acm(tea[, -1], scannf = FALSE)
#'
#' ggord(ord, tea$Tea)
#'
#' # nonmetric multidimensional scaling with the iris dataset
#' # metaMDS
#' library(vegan)
#' ord <- metaMDS(iris[, 1:4])
#'
#' ggord(ord, iris$Species)
#'
#' # linear discriminant analysis
#' # example from lda in MASS package
#' ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)
#'
#' ggord(ord, iris$Species)
#'
#' # correspondence analysis
#' # dudi.coa
#' ord <- dudi.coa(iris[, 1:4], scannf = FALSE, nf = 4)
#'
#' ggord(ord, iris$Species)
#'
#' # correspondence analysis
#' # ca
#' library(ca)
#' ord <- ca(iris[, 1:4])
#'
#' ggord(ord, iris$Species)
#'
ggord <- function(...) UseMethod('ggord')

#' @rdname ggord
#'
#' @export
#'
#' @method ggord default
ggord.default <- function(obs, vecs, axes = c('1', '2'),
                      arrow = 0.4, ext = 1.2, size = 4, txt = 6, xlims = NULL,
                      ylims = NULL, var_sub = NULL, ...){

  # tweaks to vecs for plotting
  # create vecs label  from vecs for labels
  names(vecs) <- c('one', 'two')
  vecs_lab <- ext * vecs
  vecs_lab$labs <- row.names(vecs_lab)
  vecs$lab <- row.names(vecs)

  # remove vectors for easier viz
  if(!is.null(var_sub)){

    var_sub <- paste(var_sub, collapse = '|')
    vecs <- vecs[grepl(var_sub, vecs$lab), ]
    vecs_lab <- vecs_lab[grepl(var_sub, vecs_lab$lab), ]

  }

  ## plots

  # individual points
  nms <- names(obs)[1:2]
  names(obs)[1:2] <- c('one', 'two')
  p <- ggplot(obs, aes_string(x = 'one', y = 'two')) +
    geom_point(size = size) +
    scale_x_continuous(name = nms[1], limits = xlims) +
    scale_y_continuous(name = nms[2], limits = ylims) +
    theme_bw()

  if(!is.null(obs$Groups))
    p <- p + geom_point(aes_string(colour = 'Groups'), size = size)

  # add vectors
  if(!is.null(arrow))
    p <- p + geom_segment(
      data = vecs,
      aes_string(x = 0, y = 0, xend = 'one', yend = 'two'),
      arrow = arrow(length = unit(arrow, "cm"))
    )

  # add labels
  if(!is.null(txt))
    p <- p + geom_text(data = vecs_lab, aes_string(x = 'one', y = 'two', label = 'labs'),
                       size = txt)

  return(p)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord PCA
ggord.PCA <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- ord_in$eig[as.numeric(axes), 'percentage of variance']
  axes <- paste0('Dim.', axes)
  obs <- data.frame(ord_in$ind$coord)
  obs <- obs[, names(obs) %in% axes]
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$var$coord)
  vecs <- vecs[, names(vecs) %in% axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord MCA
ggord.MCA <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- ord_in$eig[as.numeric(axes), 'percentage of variance']
  axes <- paste0('Dim.', axes)
  obs <- data.frame(ord_in$ind$coord)
  obs <- obs[, names(obs) %in% axes]
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$var$coord)
  vecs <- vecs[, names(vecs) %in% axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord mca
ggord.mca <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- ord_in$eig[as.numeric(axes), 'percentage of variance']
  obs <- data.frame(ord_in$rs[, as.numeric(axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$cs[, as.numeric(axes)])
  exp_var <- 100 * ord_in$d^2/sum(ord_in$d^2)
  exp_var <- exp_var[as.numeric(axes)]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord acm
ggord.acm <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- 100 * ord_in$eig^2 / sum(ord_in$eig^2)
  exp_var <- exp_var[as.numeric(axes)]
  obs <- data.frame(ord_in$li[, paste0('Axis', axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$co[, paste0('Comp', axes)])
  axes <- paste0('Axis', axes)
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord prcomp
ggord.prcomp <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  axes <- paste0('PC', axes)
  obs <- data.frame(ord_in$x[, axes])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$rotation[, axes])
  exp_var <- 100 * summary(ord_in)$importance[2, ][axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord princomp
ggord.princomp <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  axes <- paste0('Comp.', axes)
  obs <- data.frame(ord_in$scores[, axes])
  obs$Groups <- grp_in
  vecs <- loadings(ord_in)
  dims <- dim(vecs)
  rownms <- row.names(loadings(ord_in))
  colnms <- colnames(loadings(ord_in))
  vecs <- matrix(vecs, nrow = dims[1], ncol = dims[2])
  vecs <- data.frame(vecs, row.names = rownms)
  names(vecs) <- colnms
  vecs <- vecs[, axes]
  exp_var <- 100 * (ord_in$sdev^2)/sum(ord_in$sdev^2)
  exp_var <- exp_var[axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord metaMDS
ggord.metaMDS <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  axes <- paste0('MDS', axes)
  obs <- data.frame(ord_in$points[, axes])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$species[, axes])

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord lda
ggord.lda <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- 100 * ord_in$svd^2 / sum(ord_in$svd^2)
  exp_var <- exp_var[as.numeric(axes)]
  axes <- paste0('LD', axes)
  obs <- data.frame(predict(ord_in)$x[, axes])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$scaling[, axes])
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord pca
ggord.pca <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- 100 * ord_in$eig^2 / sum(ord_in$eig^2)
  exp_var <- exp_var[as.numeric(axes)]
  obs <- data.frame(ord_in$li[, paste0('Axis', axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$co[, paste0('Comp', axes)])
  axes <- paste0('Axis', axes)
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord coa
ggord.coa <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- 100 * ord_in$eig^2 / sum(ord_in$eig^2)
  exp_var <- exp_var[as.numeric(axes)]
  obs <- data.frame(ord_in$li[, paste0('Axis', axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$co[, paste0('Comp', axes)])
  axes <- paste0('Axis', axes)
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord ca
ggord.ca <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  exp_var <- 100 * ord_in$sv^2 / sum(ord_in$sv^2)
  exp_var <- exp_var[as.numeric(axes)]
  axes <- paste0('Dim', axes)
  obs <- data.frame(ord_in$rowcoord[, axes])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$colcoord[, axes])
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}
