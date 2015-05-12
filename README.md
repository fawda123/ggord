
## ggord

#### *Marcus W. Beck, mbafs2012@gmail.com*

[![Travis-CI Build Status](https://travis-ci.org/fawda123/ggord.png?branch=master)](https://travis-ci.org/fawda123/ggord)

A simple package for creating ordination plots with ggplot2.  Install the package as follows:


```r
install.packages(devtools)
library(devtools)
install_github('fawda123/ggord')
library(ggord)
```

Examples of ordination plots for pcord (principal components), MCA (multiple correspondence analysis), and metaMDS (nonmetric multidimensional scaling) are shown below.  Additional methods not shown are also available for princomp and PCA.



```r
# principal components analysis with the iris data set
# prcomp
ord <- prcomp(iris[, 1:4])

p <- ggord(ord, iris$Species)
p
```

![](README_files/figure-html/unnamed-chunk-3-1.png) 

```r
p + scale_colour_manual('Species', values = c('purple', 'orange', 'blue'))
```

![](README_files/figure-html/unnamed-chunk-3-2.png) 

```r
p + theme_classic()
```

![](README_files/figure-html/unnamed-chunk-3-3.png) 

```r
p + theme(legend.position = 'top')
```

![](README_files/figure-html/unnamed-chunk-3-4.png) 

```r
p + scale_x_continuous(limits = c(-2, 2))
```

```
## Warning: Removed 75 rows containing missing values (geom_point).
```

```
## Warning: Removed 75 rows containing missing values (geom_point).
```

![](README_files/figure-html/unnamed-chunk-3-5.png) 

```r
# principal components analysis with the iris dataset
# princomp
ord <- princomp(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-6.png) 

```r
# principal components analysis with the iris dataset
# PCA
library(FactoMineR)

ord <- PCA(iris[, 1:4], graph = FALSE)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-7.png) 

```r
# principal components analysis with the iris dataset
# dudi.pca
library(ade4)

ord <- dudi.pca(iris[, 1:4], scannf = FALSE, nf = 4)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-8.png) 

```r
# multiple correspondence analysis with the tea dataset
# MCA
data(tea)
tea <- tea[, c('Tea', 'sugar', 'price', 'age_Q', 'sex')]
ord <- MCA(tea[, -1], graph = FALSE)

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-9.png) 

```r
# multiple correspondence analysis with the tea dataset
# mca
library(MASS)
ord <- mca(tea[, -1])

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-10.png) 

```r
# nonmetric multidimensional scaling with the iris dataset
# metaMDS
library(vegan)
ord <- metaMDS(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-11.png) 

```r
# linear discriminant analysis
# example from lda in MASS package
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-12.png) 

