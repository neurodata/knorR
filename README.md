[![Build Status](https://travis-ci.org/neurodata/knorR.svg?branch=master)](https://travis-ci.org/neurodata/knorR)[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ggplot2)](https://cran.r-project.org/package=knor)[![Downloads](http://cranlogs.r-pkg.org/badges/knor?color=brightgreen)](http://www.r-pkg.org/pkg/knor)

# knor (K-means NUMA optimized routines)

- [Repo contents](#repo-contents)
- [Best Performance configuration](#best-performance-configuration)
- [Dependencies](#r-dependencies)
- [Installation](#stable-builds)
- [Docker](#docker)
- [Examples](#examples)
- [Test data](#test-data)
- [Reproduction and Verification](#reproduction-and-verification)

## Repo contents

- [**R**](https://github.com/flashxio/knorR/tree/master/R): `R` building blocks for user interface code. Internally called by user interface.
- [**data**](https://github.com/flashxio/knorR/tree/master/data): Data files for testing.
- [**inst**](https://github.com/flashxio/knorR/tree/master/inst): Citation files
- [**man**](https://github.com/flashxio/knorR/tree/master/man): Package documentation
- [**src**](https://github.com/flashxio/knorR/tree/master/src): `R` bindings interface and C++ submodule to base repo.
- [**tests**](https://github.com/flashxio/knorR/tree/master/tests): `R` unit tests written using the `testthat` package.

R bindings for *k-means* NUMA optimized routines. This package is supported for **Linux**, **Mac OSX** and **Windows**.

**NOTE**: This is a package from C++ source that will compile using your
`gcc` compiler.

## Tested on
- Mac OSX: 10.11 (El Capitan), 10.12 (Sierra), 10.13 (High Sierra)
- Linux: Ubuntu 14.04, 16.04, CentOS 6, Fedora 25, Fedora 26
- Windows: 8.1, 10

## Hardware requirements
- Any machine with >= 2 GB RAM

## License

Our software is licensed under the [Apache version 2.0 license](https://github.com/flashxio/knor/blob/master/LICENSE).

## Best Performance configuration

For the best performance on **Linux** make sure the `numa` system package is installed via

```
apt-get install -y build-essential libnuma-dbg libnuma-dev libnuma1
```

#### `R` Dependencies

- We require a recent version of `Rcpp` (`install.packages("Rcpp")`)
- We recommend the `testthat` package if you want to run unit-tests (`install.packages("testthat")`)

### Stable builds

Install from CRAN directly. Installation time is normally **~2min**.

```
install.packages("knor")
```

### Bleeding edge install

Install directly from Github. This has dependency on the following system packages:

- `git`
- `autoconf`

```
git clone --recursive https://github.com/flashxio/knorR.git
cd knorR
./install.sh
```

**Mac:** Install via `brew install autoconf`

**Ubuntu:** Install via `apt-get install autoconf`

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

To work with data from disk simply use binary row-major data. Please see [this link](https://github.com/flashxio/knor#data-format) for a detailed description.

```
fn <- "/path/to/file.bin" # Use real file
k <- 2 # The number of clusters
nrow <- 50 # The number of rows
ncol <- 5 # The number of columns
kms <-Kmeans(fn, nrow, ncol, k, init="kmeanspp", nthread=2)
```

# Test data

We provide [test data](https://github.com/flashxio/knorR/tree/master/data) that is included as part of the package and can be accessed directly via [this link](https://github.com/flashxio/knorR/tree/master/data) or through the `R` interpreter after the package is `require`d in `R` as `knor::test_data`.


# Reproduction and Verification

```
require(knor)
kms <- Kmeans(knor::test_data, knor::test_centroids)
```

**Expected output**:

Runtime for this action should be nearly **instantaneous** on any machine:

```
> kms
$nrow
[1] 50

$ncol
[1] 5

$iters
[1] 5

$k
[1] 8

$centers
         [,1]     [,2]     [,3]     [,4]     [,5]
[1,] 2.881889 4.079735 4.243061 1.953790 2.690649
[2,] 2.494522 2.334093 2.204031 4.161763 2.444349
[3,] 3.630086 2.398294 3.793616 2.404824 4.490043
[4,] 3.909759 3.991190 2.947161 3.762090 1.950588
[5,] 4.574327 3.645658 3.975175 4.505870 3.595890
[6,] 3.190091 4.267428 1.643788 3.229366 3.700539
[7,] 2.110254 3.147714 2.153235 1.581510 3.102312
[8,] 2.186852 2.027695 3.938736 1.410910 2.383727

$cluster
 [1] 3 2 3 3 6 8 8 3 3 2 3 4 7 7 5 4 2 1 2 1 2 7 7 5 1 1 8 7 5 2 6 2 4 6 6 8 2 5
[39] 7 4 6 5 6 4 7 4 5 4 2 5

$size
[1] 4 9 6 7 7 6 7 4
```

## Help
Check the R docs provided:

```
?knor::Kmeans
```
