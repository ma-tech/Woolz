#! /bin/tcsh -f
##
# \file         binWlzApp/WlzObjCompareSpecial_01.csh
# \author       Richard Baldock
# \date         July 2000
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
# \brief	Calculates domain comparison values.
# \ingroup	BinWlzApp
##


set filelist = $1

echo " "
echo Calculating domain comparison values for $filelist
echo Path: `pwd`
echo User: `whoami`
echo Date: `date`
echo " "
echo Output from script: $0

set files = `cat $1`
set i = 1

echo " "
echo Domain 1, Domain 2, volume 1, volume 2, mass 1, mass 2, cm-cm, cm-surf, surf-cm, surf-surf

while ( $i <  $#files )
    set cosmidfile = $files[$i]
    @ i += 1
    set terrfile = $files[$i]
    @ i += 1
    echo -n $cosmidfile, $terrfile, " "
    cat $cosmidfile $terrfile | WlzObjCompareSpecial_01
end
