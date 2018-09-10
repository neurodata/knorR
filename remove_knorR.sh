#!/bin/bash

for loc in `Rscript -e "cat(.libPaths())"`; do
    path=$loc/knorR
    echo "Removing $path"
    rm -rvf $path
done
