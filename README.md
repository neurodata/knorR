[![Build
Status](https://travis-ci.org/flashxio/knorR.svg?branch=master)](https://travis-ci.org/flashxio/knorR)

# knor

R bindings for *k-means* NUMA optimized routines. This package is supported for **Linux**, **Mac OSX** and **Windows**

## Best Performance configuration

For the best performance on **Linux** make sure the `numa` system package is installed via

```
apt-get install -y libnuma-dbg libnuma-dev libnuma1
```

#### `R` Dependencies

- We require a recent version of `Rcpp` (`install.packages("Rcpp")`)
- We recommend the `testthat` package if you want to run unit-tests (`install.packages("testthat")`)

### Stable builds

Install from CRAN directly.

```
install.packages("knor")
```

### Bleeding edge install

Install directly from Github. This has dependency on the following system packages:
- git
- autoconf

```
git clone --recursive https://github.com/flashxio/knorR.git
cd knorR
./install.sh
```

**NOTE:** The command may require administrator privileges (i.e., `sudo`)

# Docker

A Docker images with all dependencies installed can be obtained by:

```
docker pull flashxio/knorr-base
```

**NOTE**: The knor `R` package must still be installed on this image via:
`install.packages("knor")`

If you prefer to build the image yourself, you can use this
[Dockerfile](https://github.com/flashxio/knor/tree/dev/R/Dockerfile)

# Examples:

## Work with data already in-memory
```
iris.mat <- as.matrix(iris[,1:4])
k <- length(unique(iris[, dim(iris)[2]])) # Number of unique classes
kms <- Kmeans(iris.mat, k)
```
## Work with data from disk

To work with data from disk simply use binary row-major data. Please see [this link](TODO) for a detailed description.

```
fn <- "/path/to/file.bin" # Use real file
k <- 2 # The number of clusters
nrow <- 50 # The number of rows
ncol <- 5 # The number of columns
kms <-Kmeans(fn, nrow, ncol, k, init="kmeanspp", nthread=2)
```

## Help
Check the R docs

```
?knor::Kmeans
```
