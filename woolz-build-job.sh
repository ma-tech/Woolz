#!/bin/sh
#$ -N woolz-build-job
#$ -cwd
#$ -l h_rt=05:00:00
#$ -l h_vmem=4G
# This is a build script for Woolz which is intended to be submitted as a
# SGE job on a a cluster, however since it's just a shell script it may
# well be a possible route to downloading/building/installing Woolz on
# other Unix/Linux systems. For other builds look at build.sh too.

# If the cluster uses modules uncomment and add modules as appropriate
#. /etc/profile.d/modules.sh
# An assembler such as nasm is required to build the external jpeg library.
# module add igmm/apps/nasm/2.12.02

set -x

# Set the environment variable WOOLZ to the install/build root
# sources will be pulled from github and built in $WOOLZ/src
# Keep Woolz out of your system directories as it has it's own
# copy of some libraries such as jpeg and tiff.
export WOOLZ='/exports/igmm/software/pkg/el7/apps/Woolz/1.7.7'
export PATH="$WOOLZ"/bin:$PATH
export LD_LIBRARY_PATH="$WOOLZ"/lib:$LD_LIBRARY_PATH

cd $WOOLZ
mkdir -p src lib bin share include
ln -s lib lib64

cd "$WOOLZ"/src
git clone https://github.com/ma-tech/Woolz.git
git clone https://github.com/ma-tech/External.git
git clone https://github.com/ma-tech/ctypesgen.git
git clone https://github.com/ma-tech/PyWoolz.git

##External - libraries required to support non-Woolz file formats.
##  Jpeg
cd "$WOOLZ"/src/External/Jpeg
rm -rf libjpeg-turbo-1.5.1
tar -zxf libjpeg-turbo-1.5.1.tar.gz
cd libjpeg-turbo-1.5.1
autoreconf -fi
./configure --prefix $WOOLZ \
            --disable-shared \
            --enable-static \
            --with-jpeg7 \
            --with-pic
make
make install
##  Tiff
cd "$WOOLZ"/src/External/Tiff
rm -rf tiff-4.0.8
tar -zxf tiff-4.0.8.tar.gz
cd tiff-4.0.8
./configure --prefix=$WOOLZ \
            --disable-shared \
            --enable-static \
            --with-pic \
            --with-jpeg-include-dir=$WOOLZ/include \
            --with-jpeg-lib-dir==$WOOLZ/lib
make
make install
##  NIfTI
cd "$WOOLZ"/src/External/NIfTI
tar -zxf nifticlib-2.0.0.tar.gz
cd nifticlib-2.0.0
cmake -DCMAKE_INSTALL_PREFIX:PATH=$WOOLZ -DBUILD_SHARED_LIBS:BOOL=OFF
make VERBOSE=1
make install
cd $WOOLZ/include/
ln -s nifti/* .

## Woolz - the Woolz libraries and commandline applications.
cd "$WOOLZ"/src/Woolz
git checkout -b release-1.7.7 remotes/origin/release-1.7.7
autoreconf -f -i
./configure --prefix=$WOOLZ --enable-optimise --enable-openmp \
    --enable-lto --enable-extff --with-nifti=$WOOLZ \
    --with-jpeg=$WOOLZ --with-tiff=$WOOLZ --enable-test \
    --with-pic --enable-static --disable-shared
make
make install

## ctypesgen - PyWoolz needs a modified version of ctypesgen.
export PYTHONPATH=$WOOLZ/lib/python2.7/site-packages/:$PYTHONPATH
cd "$WOOLZ"/src/ctypesgen
./setup.py install --prefix=$WOOLZ

## PyWoolz
cd "$WOOLZ"/src/PyWoolz
make WLZDIR=$WOOLZ
cp libPyWlz.so "$WOOLZ"/lib
cp Wlz.py "$WOOLZ"/lib/python2.7/site-packages/


