#!/usr/bin/env bash

PACK="knorR*.tar.gz"
./cleanup.sh
autoconf
./configure
rm -f $PACK
R CMD build .
R CMD INSTALL $PACK
