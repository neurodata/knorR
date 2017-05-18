# knorR

R bindings for *k-means* NUMA optimized routines. This package runs on **Linux** only!

# Installation

The following packages must be installed as follows on the system in order to compile knorR.

```
libssl-dev libxml2-dev libcurl4-openssl-dev libnuma-dbg libnuma-dev libnuma1 libboost-all-dev r-base-core
```

For Ubuntu 16.04 this can be done as follows:

```
apt install -y libssl-dev libxml2-dev libcurl4-openssl-dev libnuma-dbg libnuma-dev libnuma1 libboost-all-dev r-base-core
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
install.packages("knorR")
```

### Help
Check the R docs

```
??knorR::kmeans
```
