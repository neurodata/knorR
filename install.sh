#!/usr/bin/env bash

PACK="knorR*.tar.gz"
rm -f configure
aclocal; autoconf
rm -f $PACK
R CMD build .
R CMD INSTALL $PACK
