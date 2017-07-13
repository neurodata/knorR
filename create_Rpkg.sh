#!/bin/sh

set -e

CDir=../knor
RDir=.

cd $RDir
cp -rvf $CDir/libkcommon/*.*pp $RDir/src/libkcommon/.
cp -rvf $CDir/libman/*.*pp $RDir/src/libman/.
cp -rvf $CDir/libauto/*.*pp $RDir/src/libauto/.
cp -rvf $CDir/binding/*.*pp $RDir/src/binding/.
aclocal; autoconf
#./install.sh
