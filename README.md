
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
ord <- prcomp(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-1.png) 

```r
# multiple correspondence analysis with farms data set
library(FactoMineR)
library(MASS)
ord <- MCA(farms, graph = FALSE)

ggord(ord)
```

![](README_files/figure-html/unnamed-chunk-3-2.png) 

```r
# nonmetric multidimensional scaling with the iris dataset
# metaMDS
library(vegan)
ord <- metaMDS(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-3.png) 

