#!/bin/sh
# \file         binWlzApp/WlzDomainsByPlane.sh
# \author       Bill Hill
# \date         November 2012
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
# \brief        Computes a table of domain occupancy by plane.
# \ingroup      BinWlz

P='WlzDomainsByPlane'
D=$1
O=/tmp/$P-$$
if [ \( -z "$D" \) -o \( ! -d "$D" \)  -o \( "$D" = "-h" \) ]
then
  echo "Usage: $P <dir>"
  echo "Given a directory containing only Woolz domains (all one domain per"
  echo "file and all files ending with .wlz) this script outputs a two"
  echo "column file with the first column having the plane coordinate and"
  echo "the second a comma seperated list of the domain files which intersect"
  echo "the plane. If no planes intersect a plane the string \"empty\" is"
  echo "used."
  exit 1
else
  echo "" >$O
  for W in $D/*.wlz
  do
    echo -n "$W ">>$O
    WlzDomainOccupancy -c -r $W >>$O
  done
 # Find the first and last plane in the file
  awk '
  BEGIN \
  {
    first = 1;
  }
  NF == 2 \
  {
    dom = $1;
    npl = split($2, pln, ",");
    if(first != 0)
    {
      pl1 = pln[1];
      lpl = pln[npl];
      first = 0;
    }
    else
    {
      if(pl1 > pln[1])
      {
        pl1 = pln[1];
      }
      if(lpl < pln[npl])
      {
        lpl = pln[npl];
      }
    }
    for(i = 1; i <= npl; ++i)
    {
      if("z" == domains[pln[i]] "z")
      {
        domains[pln[i]] = dom;
      }
      else
      {
        domains[pln[i]] = domains[pln[i]] "," dom;
      }
    }
  }
  END \
  {
    for(p = pl1; p <= lpl; ++p)
    {
      if("z" == domains[p] "z")
      {
        print p, "empty"
      }
      else
      {
        print p, domains[p]
      }
    }
  }
  {
  }' $O
  rm -f $O
fi
