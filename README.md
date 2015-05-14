
## ggord

#### *Marcus W. Beck, mbafs2012@gmail.com*

[![Travis-CI Build Status](https://travis-ci.org/fawda123/ggord.png?branch=master)](https://travis-ci.org/fawda123/ggord)

A simple package for creating ordination plots with ggplot2 (aka reinventing the wheel, see [this](https://github.com/vqv/ggbiplot) and [this](https://github.com/kassambara/factoextra)).  Install the package as follows:


```r
install.packages(devtools)
library(devtools)
install_github('fawda123/ggord')
library(ggord)
```

The following shows some examples of creating biplots using the methods available with ggord.  These methods were developed independently from the [ggbiplot](https://github.com/vqv/ggbiplot) and [factoextra](https://github.com/kassambara/factoextra), though the biplots are practically identical.  I made liberal use of the ellipses feature from ggbiplot, so credit is given where credit is due.  Most methods are for results from principal components analysis, although methods are available for nonmetric multidimensional scaling, multiple correspondence analysis, correspondence analysis, and linear discriminant analysis.  Available methods are as follows:

```
##  [1] ggord.acm      ggord.ca       ggord.coa      ggord.default 
##  [5] ggord.lda      ggord.mca      ggord.MCA      ggord.metaMDS 
##  [9] ggord.pca      ggord.PCA      ggord.prcomp   ggord.princomp
```

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
```

```
## Warning: package 'FactoMineR' was built under R version 3.1.3
```

```r
ord <- PCA(iris[, 1:4], graph = FALSE)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-7.png) 

```r
# principal components analysis with the iris dataset
# dudi.pca
library(ade4)
```

```
## Warning: package 'ade4' was built under R version 3.1.3
```

```r
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
# multiple correspondence analysis with the tea dataset
# acm
ord <- dudi.acm(tea[, -1], scannf = FALSE)

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-11.png) 

```r
# nonmetric multidimensional scaling with the iris dataset
# metaMDS
library(vegan)
```

```
## Warning: package 'vegan' was built under R version 3.1.3
```

```r
ord <- metaMDS(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-12.png) 

```r
# linear discriminant analysis
# example from lda in MASS package
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-13.png) 

```r
# correspondence analysis
# dudi.coa
ord <- dudi.coa(iris[, 1:4], scannf = FALSE, nf = 4)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-14.png) 

```r
# correspondence analysis
# ca
library(ca)
```

```
## Warning: package 'ca' was built under R version 3.1.3
```

```r
ord <- ca(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-15.png) 

