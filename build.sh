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

# The default build; optimised, using openmp and external file formats.
# This requires jpeg, nifti and tiff libraries.
./configure --prefix=/opt/MouseAtlas --enable-optimise --enable-openmp --enable-extff --with-nifti=/opt/MouseAtlas

# Build optimised, using openmp, external file formats and the Java binding.
# This requires javacc; the jpeg, nifti and tiff libraries.
#./configure --prefix=/opt/MouseAtlas --enable-optimise --enable-openmp --enable-extff --with-nifti=/opt/MouseAtlas --enable-java --enable-shared --disable-static --with-jdk=/opt/java --with-javacc=/opt/JavaCC/bin/javacc --with-pic

# Build the Woolz code unoptimised with debug support.
# This requires only standard system libraries.
#./configure --prefix=/opt/MouseAtlas --enable-debug

# Build the Woolz code unoptimised with debug support, external file formats,
# the Java binding and tests.
# This requires javacc; the jpeg, nifti and tiff libraries.
#./configure --prefix=/opt/MouseAtlas --enable-debug --enable-extff --enable-test --with-nifti=/opt/MouseAtlas --enable-java --enable-shared --disable-static --with-jdk=/opt/java --with-javacc=/opt/JavaCC/bin/javacc --with-pic
