#! /bin/tcsh -f
##
# \file         binWlzExtFF/WlzToXYZ.csh
# \author       Richard Baldock
# \date         November 2000
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2012],
# The University Court of the University of Edinburgh,
# Old College, Edinburgh, UK.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief	Woolz to x, y, z files.
# \ingroup	BinWlzExtFF
##


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

