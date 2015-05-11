
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

Examples of ordination plots for pcord (principal components) and MCA (multiple correspondence analysis) are shown below.



```r
# pca with the iris data set
ord <- prcomp(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-1.png) 

```r
# mca with the farms dataset
library(FactoMineR)
library(MASS)

mod <- MCA(farms, graph = FALSE)
ggord(mod)
```

![](README_files/figure-html/unnamed-chunk-3-2.png) 

