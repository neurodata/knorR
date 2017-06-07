#!/bin/bash

#echo 'roxygen2::roxygenize("./")' | R --no-save
#R -e 'library(roxygen2);roxygenise()'
R -e 'library(devtools);document()'
