
## ggord

#### *Marcus W. Beck, mbafs2012@gmail.com*

[![Travis-CI Build Status](https://travis-ci.org/fawda123/ggord.png?branch=master)](https://travis-ci.org/fawda123/ggord)[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/fawda123/ggord?branch=master&svg=true)](https://ci.appveyor.com/project/fawda123/ggord)[![DOI](https://zenodo.org/badge/35334615.svg)](https://zenodo.org/badge/latestdoi/35334615)

A simple package for creating ordination plots with ggplot2 (aka reinventing the wheel, see [this](https://github.com/vqv/ggbiplot) and [this](https://github.com/kassambara/factoextra)).  

### Installation

Install the package as follows:


```r
install.packages('devtools')
library(devtools)
install_github('fawda123/ggord')
library(ggord)
```

### Citation

Please cite the current release as follows:

Marcus W. Beck (2017). ggord: Ordination Plots with ggplot2. R package version 1.0.0. [https://zenodo.org/badge/latestdoi/35334615](https://zenodo.org/badge/latestdoi/35334615)

### Usage

The following shows some examples of creating biplots using the methods available with ggord.  These methods were developed independently from the [ggbiplot](https://github.com/vqv/ggbiplot) and [factoextra](https://github.com/kassambara/factoextra) packages, though the biplots are practically identical.  Most methods are for results from principal components analysis, although methods are available for nonmetric multidimensional scaling, multiple correspondence analysis, correspondence analysis, and linear discriminant analysis.  Available methods are as follows:

```
##  [1] ggord.acm      ggord.ca       ggord.capscale ggord.cca     
##  [5] ggord.coa      ggord.default  ggord.dpcoa    ggord.lda     
##  [9] ggord.mca      ggord.MCA      ggord.metaMDS  ggord.pca     
## [13] ggord.PCA      ggord.ppca     ggord.prcomp   ggord.princomp
## [17] ggord.rda     
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
# transparent ellipses
p <- ggord(ord, iris$Species, poly = FALSE)
p
```

![](README_files/figure-html/unnamed-chunk-3-6.png)<!-- -->

```r
# change linetype for transparent ellipses
p <- ggord(ord, iris$Species, poly = FALSE, polylntyp = iris$Species)
p
```

![](README_files/figure-html/unnamed-chunk-3-7.png)<!-- -->

```r
# convex hulls 
p <- ggord(ord, iris$Species, ellipse = FALSE, hull = TRUE)
p
```

![](README_files/figure-html/unnamed-chunk-3-8.png)<!-- -->

```r
# change the vector labels with vec_lab
new_lab <- list(Sepal.Length = 'SL', Sepal.Width = 'SW', Petal.Width = 'PW',
  Petal.Length = 'PL')
p <- ggord(ord, iris$Species, vec_lab = new_lab)
p
```

![](README_files/figure-html/unnamed-chunk-3-9.png)<!-- -->

```r
# observations as labels from row names
p <- ggord(ord, iris$Species, obslab = TRUE)
p
```

![](README_files/figure-html/unnamed-chunk-3-10.png)<!-- -->

```r
# map a variable to point sizes
p <- ggord(ord, grp_in = iris$Species, size = iris$Sepal.Length, sizelab = 'Sepal\nlength')
p
```

![](README_files/figure-html/unnamed-chunk-3-11.png)<!-- -->

```r
# change vector scaling, arrow length, line color, size, and type
p <- ggord(ord, grp_in = iris$Species, arrow = 1, vec_ext = 3, veccol = 'red', veclsz = 1, vectyp = 'dotted')
p
```

![](README_files/figure-html/unnamed-chunk-3-12.png)<!-- -->

```r
# change color of text labels on vectors, use ggrepel to prevent text overlap
p <- ggord(ord, grp_in = iris$Species, labcol = 'purple', repel = TRUE)
p
```

![](README_files/figure-html/unnamed-chunk-3-13.png)<!-- -->

```r
# faceted by group
p <- ggord(ord, iris$Species, facet = TRUE, nfac = 1)
p
```

![](README_files/figure-html/unnamed-chunk-3-14.png)<!-- -->

```r
# principal components analysis with the iris dataset
# princomp
ord <- princomp(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-15.png)<!-- -->

```r
# principal components analysis with the iris dataset
# PCA
library(FactoMineR)

ord <- PCA(iris[, 1:4], graph = FALSE)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-16.png)<!-- -->

```r
# principal components analysis with the iris dataset
# dudi.pca
library(ade4)

ord <- dudi.pca(iris[, 1:4], scannf = FALSE, nf = 4)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-17.png)<!-- -->

```r
# multiple correspondence analysis with the tea dataset
# MCA
data(tea, package = 'FactoMineR')
tea <- tea[, c('Tea', 'sugar', 'price', 'age_Q', 'sex')]

ord <- MCA(tea[, -1], graph = FALSE)

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-18.png)<!-- -->

```r
# multiple correspondence analysis with the tea dataset
# mca
library(MASS)

ord <- mca(tea[, -1])

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-19.png)<!-- -->

```r
# multiple correspondence analysis with the tea dataset
# acm
ord <- dudi.acm(tea[, -1], scannf = FALSE)

ggord(ord, tea$Tea)
```

![](README_files/figure-html/unnamed-chunk-3-20.png)<!-- -->

```r
# nonmetric multidimensional scaling with the iris dataset
# metaMDS
library(vegan)
ord <- metaMDS(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-21.png)<!-- -->

```r
# linear discriminant analysis
# example from lda in MASS package
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-22.png)<!-- -->

```r
# correspondence analysis
# dudi.coa
ord <- dudi.coa(iris[, 1:4], scannf = FALSE, nf = 4)

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-23.png)<!-- -->

```r
# correspondence analysis
# ca
library(ca)
ord <- ca(iris[, 1:4])

ggord(ord, iris$Species)
```

![](README_files/figure-html/unnamed-chunk-3-24.png)<!-- -->

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

![](README_files/figure-html/unnamed-chunk-3-25.png)<!-- -->

```r
# phylogenetic PCA
# ppca
library(adephylo)
library(phylobase)
library(ape)

data(lizards)

# example from help file, adephylo::ppca
# original example from JOMBART ET AL 2010

# build a tree and phylo4d object
liz.tre <- read.tree(tex=lizards$hprA)
liz.4d <- phylobase::phylo4d(liz.tre, lizards$traits)

# remove duplicated populations
liz.4d <- phylobase::prune(liz.4d, c(7,14))


# correct labels
lab <- c("Pa", "Ph", "Ll", "Lmca", "Lmcy", "Phha", "Pha",
         "Pb", "Pm", "Ae", "Tt", "Ts", "Lviv", "La", "Ls", "Lvir")
tipLabels(liz.4d) <- lab

# remove size effect
dat <- tdata(liz.4d, type="tip")
dat <- log(dat)
newdat <- data.frame(lapply(dat, function(v) residuals(lm(v~dat$mean.L))))
rownames(newdat) <- rownames(dat)
tdata(liz.4d, type="tip") <- newdat[,-1] # replace data in the phylo4d object

# create ppca
liz.ppca <- ppca(liz.4d,scale=FALSE,scannf=FALSE,nfposi=1,nfnega=1, method="Abouheif")

# plot
ggord(liz.ppca)
```

![](README_files/figure-html/unnamed-chunk-3-26.png)<!-- -->

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

![](README_files/figure-html/unnamed-chunk-3-27.png)<!-- -->

```r
# distance-based redundancy analysis, from vegan
ord <- capscale(varespec ~ N + P + K + Condition(Al), varechem, dist = "bray")

ggord(ord)
```

![](README_files/figure-html/unnamed-chunk-3-28.png)<!-- -->

```r
# canonical correspondence analysis
# cca from vegan
ord <- cca(varespec, varechem)

ggord(ord)
```

![](README_files/figure-html/unnamed-chunk-3-29.png)<!-- -->

```r
# species points as text
# suppress site points
ggord(ord, ptslab = TRUE, size = NA, addsize = 5, parse = TRUE)
```

![](README_files/figure-html/unnamed-chunk-3-30.png)<!-- -->

