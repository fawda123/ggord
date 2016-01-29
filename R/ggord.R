#' Ordination plots with ggplot2
#'
#' Create an ordination biplot using ggplot2 including options for selecting axes, group color aesthetics, and selection of variables to plot.
#'
#' @param ord_in input ordination object
#' @param grp_in vector of grouping objects for the biplot, must have the same number of observations as the original matrix used for the ordination
#' @param obs matrix or data frame of axis scores for each observation
#' @param vecs matrix or data frame of axis scores for each variable
#' @param axes chr string indicating which axes to plot
#' @param addpts optional matrix or data.frame of additional points if constrained ordination is used (e.g., species locations in cca, rda)
#' @param obslab logical if the row names for the observations in \code{obs} are plotted rather than points
#' @param ellipse logical if confidence ellipses are shown for each group, method from the ggbiplot package
#' @param ellipse_pro numeric indicating confidence value for the ellipses
#' @param arrow numeric indicating length of the arrow heads on the vectors, use \code{NULL} to suppress arrows
#' @param ext numeric indicating scalar distance of the labels from the arrow ends
#' @param vec_ext numeric indicating a scalar extension for the ordination vectors
#' @param vec_lab list of optional labels for vectors, defaults to names from input data.  The input list must be named using the existing variables in the input data.  Each element of the list will have the desired name change.
#' @param size numeric indicating size of the observation points
#' @param addsize numeric indicating size of the species points if addpts is not \code{NULL}
#' @param addcol numeric indicating color of the species points if addpts is not \code{NULL}
#' @param addpch numeric indicating point type of the species points if addpts is not \code{NULL}
#' @param txt numeric indicating size of the text labels for the vectors, use \code{NULL} to suppress labels
#' @param alpha numeric transparency of points and ellipses from 0 to 1
#' @param xlims two numeric values indicating x-axis limits
#' @param ylims two numeric values indicating y-axis limits
#' @param var_sub chr string indcating which labels to show.  Regular expression matching is used.
#' @param coord_fix logical indicating fixed, equal scaling for axes
#' @param parse logical indicating if optional vector labels are expressions to parse with \code{\link[ggplot2]{geom_text}}
#' @param ... arguments passed to or from other methods
#'
#' @export
#'
#' @import ggplot2 plyr
#'
#' @return A \code{\link[ggplot2]{ggplot}} object that can be further modified
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#'
#' library(ggplot2)
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
#'
#' # change the vector labels with vec_lab
#' new_lab <- list(Sepal.Length = 'SL', Sepal.Width = 'SW', Petal.Width = 'PW',
#'  Petal.Length = 'PL')
#' p <- ggord(ord, iris$Species, vec_lab = new_lab)
#' p
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
#' data(tea, package = 'FactoMineR')
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
#' # rda triplot
#' data(varespec)
#' data(varechem)
#' ord <- rda(varespec, varechem)
#'
#' ggord(ord)
ggord <- function(...) UseMethod('ggord')

#' @rdname ggord
#'
#' @export
#'
#' @method ggord default
ggord.default <- function(obs, vecs, axes = c('1', '2'), addpts = NULL, obslab = FALSE,
                      ellipse = TRUE, ellipse_pro = 0.95, arrow = 0.4, ext = 1.2, vec_ext = 1,
                      vec_lab = NULL, size = 4, addsize = size/2, addcol = 'blue', addpch = 19,
                      txt = 4, alpha = 1, xlims = NULL, ylims = NULL, var_sub = NULL,
                      coord_fix = TRUE, parse = FALSE, ...){

  # extend vectors by scale
  vecs <- vecs * vec_ext

  # tweaks to vecs for plotting
  # create vecs label  from vecs for labels
  names(vecs) <- c('one', 'two')
  vecs_lab <- ext * vecs
  if(is.null(vec_lab)) vecs_lab$labs <- as.character(row.names(vecs_lab))
  else{
    vecs_lab$labs <- vec_lab[row.names(vecs_lab)]
  }
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
  obs$lab <- row.names(obs)
  p <- ggplot(obs, aes_string(x = 'one', y = 'two')) +
    scale_x_continuous(name = nms[1], limits = xlims) +
    scale_y_continuous(name = nms[2], limits = ylims) +
    theme_bw()

  # observations as points or text, colour if groups provided
  if(obslab){
    if(!is.null(obs$Groups))
      p <- p + geom_text(aes_string(colour = 'Groups', label = 'lab')) # size = size, alpha = alpha)
    else
      p <- p + geom_text(label = row.names(obs), size = size, alpha = alpha)
  } else {
    if(!is.null(obs$Groups))
      p <- p + geom_point(aes_string(colour = 'Groups'), size = size, alpha = alpha)
    else
      p <- p + geom_point(size = size, alpha = alpha)
  }

  # add species scores if addpts not null, for triplot
  if(!is.null(addpts)){

    nms <- names(addpts)[1:2]
    names(addpts)[1:2] <- c('one', 'two')
    p <- p +
      geom_point(data = addpts, aes_string(x = 'one', y = 'two'),
        size = addsize, col = addcol, alpha = alpha, shape = addpch)

  }

  # fixed coordiantes if TRUE
  if(coord_fix)
    p <- p + coord_fixed()

  # concentration ellipse if there are groups, from ggbiplot
  if(!is.null(obs$Groups) & ellipse) {

    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))

    ell <- ddply(obs, 'Groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$one, x$two))
      mu <- c(mean(x$one), mean(x$two))
      ed <- sqrt(qchisq(ellipse_pro, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
    })
    names(ell)[2:3] <- c('one', 'two')

    p <- p + geom_path(data = ell, aes_string(color = 'Groups', group = 'Groups'), alpha = alpha)

  }

  # add vectors
  if(!is.null(arrow))
    p <- p + geom_segment(
      data = vecs,
      aes_string(x = 0, y = 0, xend = 'one', yend = 'two'),
      arrow = grid::arrow(length = grid::unit(arrow, "cm"))
    )

  # add labels
  if(!is.null(txt))
    p <- p + geom_text(data = vecs_lab, aes_string(x = 'one', y = 'two'),
      label = unlist(lapply(vecs_lab$labs, function(x) as.character(as.expression(x)))),
      size = txt, parse = parse)

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
  obs <- data.frame(ord_in$rs[, as.numeric(axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$cs[, as.numeric(axes)])
  exp_var <-  100 * ord_in$d/(ord_in$p - 1)
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
  exp_var <- 100 * ord_in$eig / sum(ord_in$eig)
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
  exp_var <- 100 * ord_in$eig / sum(ord_in$eig)
  exp_var <- exp_var[as.numeric(axes)]
  ax_nms <- paste0('Axis', axes)
  obs <- data.frame(ord_in$li[, ax_nms])
  obs$Groups <- grp_in

  # vector locations, uses scaling from factoextra
  vc_nms <- paste0('Comp', axes)
  vecs <- data.frame(ord_in$co[, vc_nms])
  r <- min((max(obs[, ax_nms[1]]) - min(obs[, ax_nms[1]])/(max(vecs[, vc_nms[1]]) -
    min(vecs[, vc_nms[1]]))), (max(obs[, ax_nms[2]]) - min(obs[, ax_nms[2]])/(max(vecs[,
    vc_nms[2]]) - min(vecs[, vc_nms[2]]))))
  vecs[, vc_nms] <- vecs[, vc_nms] * r * 0.7

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
  exp_var <- 100 * ord_in$eig / sum(ord_in$eig)
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

  # axis contribution
  exp_var <- with(ord_in, 100 * sv^2/sum(sv^2))
  exp_var <- exp_var[as.numeric(axes)]

  # data to plot
  sv <- ord_in$sv[as.numeric(axes)]
  colmass <- ord_in$colmass[as.numeric(axes)]
  axes <- paste0('Dim', axes)
  obs <- ord_in$rowcoord[, axes]

  # scale the observations and vectors correctly
  obs <- as.data.frame(t(apply(obs, 1, "*", sv)))
  obs$Groups <- grp_in
  vecs <- ord_in$colcoord[, axes]
  vecs <- as.data.frame(t(apply(vecs, 1, "*", sv)))

  # make a nice label for the axes
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes

  ggord.default(obs, vecs, axes, ...)

}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord rda
ggord.rda <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  axes <- paste0('RDA', axes)
  obs <- data.frame(ord_in$CCA$wa[, axes])
  obs$Groups <- grp_in
  addpts <- data.frame(ord_in$CCA$v[, axes])

  # vectors for constraining matrix
  constr <- data.frame(ord_in$CCA$biplot[, axes])

  ggord.default(obs, vecs = constr, axes, addpts = addpts, ...)

}
