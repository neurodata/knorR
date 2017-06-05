#!/usr/bin/env bash

PACK="knorR*.tar.gz"
rm -rf config.* configure
autoconf
rm -f $PACK
R CMD build .
R CMD INSTALL $PACK
