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

# now generate the XYZ files - taking care of the numbering
# from WlzExtFFConvert
set plane = $bounds[3]
set xOff = $bounds[1]
set yOff = $bounds[2]
@ xOff *= -1
@ yOff *= -1

while( $plane < $bounds[6] )

  if( $plane < 10 ) then
    WlzExtFFWlzToXYZ -o$xOff,$yOff -p$plane $data > $dstdir/plane0000000$plane.xyz
  else if( $plane < 100 ) then
    WlzExtFFWlzToXYZ -o$xOff,$yOff -p$plane $data > $dstdir/plane000000$plane.xyz
  else if( $plane < 1000  )then
    WlzExtFFWlzToXYZ -o$xOff,$yOff -p$plane $data > $dstdir/plane00000$plane.xyz
  else if( $plane < 10000 ) then
    WlzExtFFWlzToXYZ -o$xOff,$yOff -p$plane $data > $dstdir/plane0000$plane.xyz
  endif

  echo done plane = $plane
  @ plane += 1
end

