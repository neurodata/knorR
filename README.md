[![Build
Status](https://travis-ci.org/flashxio/knorR.svg?branch=master)](https://travis-ci.org/flashxio/knorR)

# knor

R bindings for *k-means* NUMA optimized routines. This package runs on **Linux**
and **Mac OSX** only!

## Best Performance configuration

For the best performance make sure the `numa` system package is installed via

```
apt-get install -y libnuma-dbg libnuma-dev libnuma1
```

# Installation

The following packages are required to support and install R.

```
libssl-dev libxml2-dev libcurl4-openssl-dev r-base-core
```

A one-liner for Ubuntu 14.04 - 16.04 is:

```
apt-get install -y libssl-dev libxml2-dev libcurl4-openssl-dev libnuma-dbg libnuma-dev libnuma1 r-base-core
```

### Bleeding edge install

Install directly from github

```
# Install dependencies if you don't have them
install.packages("devtools", dependencies=TRUE)

# Load and install the package
require(devtools)
install_github("flashxio/knorR")
```

### Stable builds (NOT yet available)

Install from CRAN directly.

```
install.packages("knor")
```

# Docker

A Docker images with all dependencies installed can be obtained by:

```
docker pull flashxio/knorr-base
```

**NOTE**: The knor package must still be installed on this image via:
`install.packages("knor")`

If you prefer to build the image yourself, you can use this
[Dockerfile](https://github.com/flashxio/knor/tree/dev/R/Dockerfile)

### Help
Check the R docs

```
?knor::kmeans
```
