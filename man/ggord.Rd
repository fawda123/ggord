% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/ggord.R
\name{ggord}
\alias{ggord}
\alias{ggord.MCA}
\alias{ggord.PCA}
\alias{ggord.default}
\alias{ggord.metaMDS}
\alias{ggord.prcomp}
\alias{ggord.princomp}
\title{Ordination plots with ggplot2}
\usage{
ggord(...)

\method{ggord}{default}(obs, vecs, axes = c("1", "2"), arrow = 0.4,
  ext = 1.2, size = 4, txt = 6, xlims = NULL, ylims = NULL,
  var_sub = NULL, ...)

\method{ggord}{PCA}(ord_in, grp_in = NULL, axes = c("1", "2"), ...)

\method{ggord}{MCA}(ord_in, grp_in = NULL, axes = c("1", "2"), ...)

\method{ggord}{prcomp}(ord_in, grp_in = NULL, axes = c("1", "2"), ...)

\method{ggord}{princomp}(ord_in, grp_in = NULL, axes = c("1", "2"), ...)

\method{ggord}{metaMDS}(ord_in, grp_in = NULL, axes = c("1", "2"), ...)
}
\arguments{
\item{...}{arguments passed to or from other methods}

\item{obs}{matrix or data frame of axis scores for each observation}

\item{vecs}{matrix or data frame of axis scores for each variable}

\item{axes}{chr string indicating which axes to plot}

\item{arrow}{numeric indicating length of the arrow heads on the vectors}

\item{ext}{numeric indicating scalar distance of the labels from the arrow ends}

\item{size}{numeric indicating size of the observatoin points}

\item{txt}{numeric indicating size of the text labels for the vectors}

\item{xlims}{two numeric values indicating x-axis limits}

\item{ylims}{two numeric values indicating y-axis limits}

\item{var_sub}{chr string indcating which labels to show.  Regular expression matching is used.}

\item{ord_in}{input ordination object}

\item{grp_in}{vector of grouping objects for the biplot, must have the same number of observations as the original matrix used for the ordination}
}
\value{
A \code{\link[ggplot2]{ggplot}} object that can be further modified
}
\description{
Create an ordination biplot using ggplot2 including options for selecting axes, group color aesthetics, and selection of variables to plot.
}
\examples{
# principal components analysis with the iris data set
# prcomp
ord <- prcomp(iris[, 1:4])

p <- ggord(ord, iris$Species)
p

p + scale_colour_manual('Species', values = c('purple', 'orange', 'blue'))
p + theme_classic()
p + theme(legend.position = 'top')
p + scale_x_continuous(limits = c(-2, 2))

# principal components analysis with the iris dataset
# princomp
ord <- princomp(iris[, 1:4])

ggord(ord, iris$Species)

# principal components analysis with the iris dataset
# PCA
library(FactoMineR)

ord <- PCA(iris[, 1:4], graph = FALSE)

ggord(ord, iris$Species)

# multiple correspondence analysis with farms data set
library(MASS)
ord <- MCA(farms, graph = FALSE)

ggord(ord)

# nonmetric multidimensional scaling with the iris dataset
# metaMDS
library(vegan)
ord <- metaMDS(iris[, 1:4])

ggord(ord, iris$Species)
}
\seealso{
\code{\link[ggplot2]{ggplot}}
}
