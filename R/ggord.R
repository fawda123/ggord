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
#' # pca with the iris data set
#' ord <- prcomp(iris[, 1:4])
#'
#' ggord(ord, iris$Species)
#'
#' # mca with farms data set
#' library(FactoMineR)
#' library(MASS)
#' ord <- MCA(farms, graph = FALSE)
#'
#' ggord(ord)
ggord <- function(...) UseMethod('ggord')

#' @rdname ggord
#'
#' @export
#'
#' @method ggord default
ggord.default <- function(obs, vecs, axes = c('1', '2'),
                      arrow = 0.4, ext = 1.5, size = 4, txt = 6, xlims = NULL,
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
  p <- ggplot(obs, aes_string(x = axes[1], y = axes[2])) +
    geom_point(size = size) +
    scale_x_continuous(limits = xlims) +
    scale_y_continuous(limits = ylims) +
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
#' @method ggord MCA
ggord.MCA <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){

  # data to plot
  axes <- paste0('Dim.', axes)
  obs <- data.frame(ord_in$ind$coord)
  obs <- obs[, names(obs) %in% axes]
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$var$coord)
  vecs <- vecs[, names(vecs) %in% axes]

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

  ggord.default(obs, vecs, axes, ...)

}

