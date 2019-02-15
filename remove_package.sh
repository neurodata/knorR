#!/bin/bash

for loc in `Rscript -e "cat(.libPaths())"`; do
    path=$loc/clusternor
    echo "Removing $path"
    rm -rvf $path
done
