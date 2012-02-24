#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstGeomArcLength2D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzTstGeomArcLength2D.c
* \author       Bill Hill
* \date         September 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	Test for WlzGeomArcLength2D().
* \ingroup	BinWlzTst
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

int		main(int argc, char *argv[])
{
  double	d;
  WlzDVertex2	a,
  		b,
		c;

  c.vtX = 200.0;
  c.vtY = 200.0;
  a.vtX = 100.0;
  a.vtY = 200.0;
  b.vtX = 300.0;
  b.vtY = 200.0;
  d = WlzGeomArcLength2D(a, b, c); 
  (void )printf("% 4.4lf\n", d);
  exit(0);
}
