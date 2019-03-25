require(clusternor)

cat("Running Kmeans\n")
example("Kmeans", package="clusternor")
print(kms)
cat("\n\n************************************************************\n\n")

cat("Running KmeansPP\n")
example("KmeansPP", package="clusternor")
print(km)
cat("\n\n************************************************************\n\n")

cat("Running MiniBatchKmeans\n")
example("MiniBatchKmeans", package="clusternor")
print(kms)
cat("\n\n************************************************************\n\n")


cat("Running Skmeans\n")
example("Skmeans", package="clusternor")
print(km)
cat("\n\n************************************************************\n\n")

cat("Running Kmedoids\n")
example("Kmedoids", package="clusternor")
print(km)
cat("\n\n************************************************************\n\n")

cat("Running FuzzyCMeans\n")
example("FuzzyCMeans", package="clusternor")
print(fcm)
cat("\n\n************************************************************\n\n")

cat("Running Hmeans\n")
example("Hmeans", package="clusternor")
print(kms)
cat("\n\n************************************************************\n\n")

cat("Running Xmeans\n")
example("Xmeans", package="clusternor")
print(xms)
cat("\n\n************************************************************\n\n")
