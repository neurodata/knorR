#!/usr/bin/env bash

PACK="knor*.tar.gz"
./cleanup.sh
autoconf
./configure
rm -f $PACK
R CMD build .
R CMD INSTALL $PACK
