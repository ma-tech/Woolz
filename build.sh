#!/bin/sh
# This script will configure and build Woolz. Uncomment the appropriate
# configure command line for the build you want. The easiest way to use
# this script is probably to copy it to mybuild.sh and the edit that script.

set -x
# In most cases a simple autoreconf should be sufficient
autoreconf
# If you hit problems with missing files or libtool use the following
# autoreconf
# autoreconf -i --force

#export MA=$HOME
#export MA=$HOME/MouseAtlas/Build/
export MA=/opt/MouseAtlas

# Build; with debug and external file formats.
# This requires jpeg, nifti and tiff libraries.
#./configure --prefix=$MA --enable-debug --enable-extff --with-tiff=$MA --with-nifti=$MA --with-jpeg=$MA

# The default build; optimised, using openmp and external file formats.
# This requires jpeg, nifti and tiff libraries.
./configure --prefix=$MA --enable-optimise --enable-openmp --enable-extff --with-nifti=$MA --with-jpeg=$MA --with-tiff=$MA --enable-test

# Build the core Woolz code unoptimised with debug support.
# This requires only standard system libraries.
#./configure --prefix=$MA --enable-debug

# Build the core and test Woolz code unoptimised with debug support.
# This requires only standard system libraries.
#./configure --prefix=$MA --enable-debug --enable-test

# Build the core Woolz code optimised
# This requires only standard system libraries.
#./configure --prefix=$MA --enable-optimise --enable-openmp

# Build the core Woolz code optimised with test code.
# This requires only standard system libraries.
#./configure --prefix=$MA --enable-optimise --enable-openmp --enable-test

