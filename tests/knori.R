require(knorR)

fn <- "~/Research/knor/test-data/matrix_r50_c5_rrw.bin"
k <- 8
nrow <- 50
ncol <- 5
nthread <- 2


# Data in memory, compute centroids
test.data.in.mem <- function() {
    cat("Data ==> memory, centroids ==> compute\n\n")
    data <- replicate(ncol, rnorm(nrow))
    print(kmeans(data, k, nrow, ncol, nthread=nthread))
}

# Data on disk, compute centroids
test.data.ex.mem <- function() {
    cat("Data ==> disk, centroids ==> compute\n\n")
    print(kmeans(fn, k, nrow, ncol, nthread=nthread))
}

# Data on disk, centroids in memory
test.centroids.in.mem <- function() {
    cat("Data ==> disk, centroids ==> memory\n\n")
    centroids <- replicate(ncol, rnorm(k))
    print(kmeans(fn, centroids, nrow, nthread=nthread))
}

# Data in memory, centroids in memory
test.data.centroids.in.mem <- function() {
    cat("Data ==> memory, centroids ==> memory\n\n")
    data <- replicate(ncol, rnorm(nrow))
    centroids <- replicate(ncol, rnorm(k))

    print(kmeans(data, centroids, nthread=nthread))
}

# Data on disk, centroids on disk
# TODO

# Data in memory, centroids on disk
# TODO


# Main
test.data.in.mem()
test.data.in.mem()
test.data.ex.mem()
test.centroids.in.mem()
test.data.centroids.in.mem()
