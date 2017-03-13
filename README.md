
## ggord

#### *Marcus W. Beck, mbafs2012@gmail.com*

[![Travis-CI Build Status](https://travis-ci.org/fawda123/ggord.png?branch=master)](https://travis-ci.org/fawda123/ggord)

A simple package for creating ordination plots with ggplot2 (aka reinventing the wheel, see [this](https://github.com/vqv/ggbiplot) and [this](https://github.com/kassambara/factoextra)).  Install the package as follows:


```r
install.packages('devtools')
library(devtools)
install_github('fawda123/ggord')
library(ggord)
```

The following shows some examples of creating biplots using the methods available with ggord.  These methods were developed independently from the [ggbiplot](https://github.com/vqv/ggbiplot) and [factoextra](https://github.com/kassambara/factoextra) packages, though the biplots are practically identical.  Most methods are for results from principal components analysis, although methods are available for nonmetric multidimensional scaling, multiple correspondence analysis, correspondence analysis, and linear discriminant analysis.  Available methods are as follows:

```
##  [1] ggord.acm      ggord.ca       ggord.cca      ggord.coa     
##  [5] ggord.default  ggord.dpcoa    ggord.lda      ggord.mca     
##  [9] ggord.MCA      ggord.metaMDS  ggord.pca      ggord.PCA     
## [13] ggord.prcomp   ggord.princomp ggord.rda     
## see '?methods' for accessing help and source code
```

```r
# principal components analysis with the iris data set
# prcomp
ord <- prcomp(iris[, 1:4])

p <- ggord(ord, iris$Species)
p
```

![](README_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
p <- ggord(ord, iris$Species, cols = c('purple', 'orange', 'blue'))
p
```

![](README_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
library(ggplot2)
p + scale_shape_manual('Groups', values = c(1, 2, 3))
```

![](README_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```r
p + theme_classic()
```

![](README_files/figure-html/unnamed-chunk-3-4.png)<!-- -->

```r
p + theme(legend.position = 'top')
```

![](README_files/figure-html/unnamed-chunk-3-5.png)<!-- -->

```r
# change the vector labels with vec_lab
new_lab <- list(Sepal.Length = 'SL', Sepal.Width = 'SW', Petal.Width = 'PW',
  Petal.Length = 'PL')
p <- ggord(ord, iris$Species, vec_lab = new_lab)
p
```

![](README_files/figure-html/unnamed-chunk-3-6.png)<!-- -->

```r
# observations as labels from row names
p <- ggord(ord, iris$Species, obslab = TRUE)
p
```

![](README_files/figure-html/unnamed-chunk-3-7.png)<!-- -->

```r
# principal components analysis with the iris dataset
# princomp
ord <- princomp(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-8.png)<!-- -->

```r
# principal components analysis with the iris dataset
# PCA
library(FactoMineR)

ord <- PCA(iris[, 1:4], graph = FALSE)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-9.png)<!-- -->

```r
# principal components analysis with the iris dataset
# dudi.pca
library(ade4)

ord <- dudi.pca(iris[, 1:4], scannf = FALSE, nf = 4)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-10.png)<!-- -->

```r
# multiple correspondence analysis with the tea dataset
# MCA
data(tea, package = 'FactoMineR')
tea <- tea[, c('Tea', 'sugar', 'price', 'age_Q', 'sex')]

ord <- MCA(tea[, -1], graph = FALSE)

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-11.png)<!-- -->

```r
# multiple correspondence analysis with the tea dataset
# mca
library(MASS)

ord <- mca(tea[, -1])

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-12.png)<!-- -->

```r
# multiple correspondence analysis with the tea dataset
# acm
ord <- dudi.acm(tea[, -1], scannf = FALSE)

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-13.png)<!-- -->

```r
# nonmetric multidimensional scaling with the iris dataset
# metaMDS
library(vegan)
ord <- metaMDS(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-14.png)<!-- -->

```r
# linear discriminant analysis
# example from lda in MASS package
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-15.png)<!-- -->

```r
# correspondence analysis
# dudi.coa
ord <- dudi.coa(iris[, 1:4], scannf = FALSE, nf = 4)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-16.png)<!-- -->

```r
# correspondence analysis
# ca
library(ca)
ord <- ca(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-17.png)<!-- -->

```r
# double principle coordinate analysis (DPCoA)
# dpcoa
library(ade4)
data(ecomor)
grp <- rep(c("Bu", "Ca", "Ch", "Pr"), each = 4)    # sample groups
dtaxo <- dist.taxo(ecomor$taxo)                    # taxonomic distance between species
ord <- dpcoa(data.frame(t(ecomor$habitat)), dtaxo, scan = FALSE, nf = 2)
 
ggord(ord, grp_in = grp, ellipse = FALSE, arrow = 0.2, txt = 3)
```

![](README_files/figure-html/unnamed-chunk-3-18.png)<!-- -->

```r
######
# triplots

# redundancy analysis
# rda from vegan
data(varespec)
data(varechem)
ord <- rda(varespec, varechem)

ggord(ord)
```

![](README_files/figure-html/unnamed-chunk-3-19.png)<!-- -->

```r
# canonical correspondence analysis
# cca from vegan
ord <- cca(varespec, varechem)

ggord(ord)
```

![](README_files/figure-html/unnamed-chunk-3-20.png)<!-- -->

```r
# species points as text
# suppress site points
ggord(ord, ptslab = TRUE, size = NA, addsize = 5)
```

![](README_files/figure-html/unnamed-chunk-3-21.png)<!-- -->

