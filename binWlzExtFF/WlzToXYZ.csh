#! /bin/tcsh -f

# set the data files
set reconstruction = $1
set data = $2
set dstdir = "$reconstruction:r"_XYZ
mkdir $dstdir

# get the bounds
set bounds = `WlzBoundingBox $reconstruction`

# convert the reconstruction to gifs
WlzExtFFConvert -fwlz -Fpgm -o$dstdir/plane.pgm $reconstruction
cd $dstdir
foreach i (*.pgm)
convert $i $i:r.gif
rm $i
end
cd ..

