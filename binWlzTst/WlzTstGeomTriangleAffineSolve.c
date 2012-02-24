#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstGeomTriangleAffineSolve_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstGeomTriangleAffineSolve.c
* \author       Bill Hill
* \date         June 2007
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
* \brief	Simple test for WlzGeomTriangleSnArea2().
* \ingroup	BinWlzTst
*/

#include <stdio.h>
#include <Wlz.h>

int		main(int argc, char *argv[])
{
  int		squashed;
  double	dd,
  		thr;
  double	xTr[3],
  		yTr[3];
  WlzDVertex2	sVx[3],
  		dVx[3];


  thr = 1.0e-6;
  sVx[0].vtX = 260; sVx[0].vtY = 60;
  sVx[1].vtX = 240; sVx[1].vtY = 60;
  sVx[2].vtX = 250; sVx[2].vtY = 40;
  dVx[0].vtX = 19;  dVx[0].vtY = 175;
  dVx[1].vtX = 17;  dVx[1].vtY = 190;
  dVx[2].vtX = 5;   dVx[2].vtY = 178;
  (void )printf("Results should be\n");
  (void )printf("  dd = 400\n");
  (void )printf("  squashed = 0\n");
  (void )printf("  xTr = 0.1, 0.65, -46\n");
  (void )printf("  yTr = -0.75, 0.225, 356.5\n");
  (void )printf("Results are\n");
  dd = WlzGeomTriangleSnArea2(sVx[0], sVx[1], sVx[2]);
  (void )printf("  dd = %g\n", dd);
  squashed = WlzGeomTriangleAffineSolve(xTr, yTr, dd, sVx, dVx, thr);
  (void )printf("  squashed = %d\n", squashed);
  (void )printf("  xTr = %g, %g, %g\n", xTr[0], xTr[1], xTr[2]);
  (void )printf("  yTr = %g, %g, %g\n", yTr[0], yTr[1], yTr[2]);
  exit(0);
}

