#!/usr/bin/env bash

PACK="knorR*.tar.gz"
rm -rf config.* configure
autoconf
configure
rm -f $PACK
R CMD build .
R CMD INSTALL $PACK
