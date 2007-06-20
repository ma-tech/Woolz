#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTstGeomTetraAffineSolve_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTstGeomTetraAffineSolve.c
* \author       Bill Hill
* \date         May 2007
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2007 Medical research Council, UK.
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
* \brief	Test program for libWlz/WlzGeomTetraAffineSolve().
* \ingroup	WlzTst
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

double		MyRandom(void)
{
  double	scale = 3.0,
  		rnd;

  rnd = exp(exp(scale * AlgRandUniform())) - exp(exp(scale * 0.5));
  return(rnd);
}

int		main(int argc, char *argv[])
{
  int		idN,
		idR,
		squashed,
  		repeats = 10;
  double	ss,
  		delta = 0.000001,
		scale;
  WlzDVertex3	eVx;
  WlzDVertex3	sVx[4],
		tVx[4],
  		dVx[4];
  WlzAffineTransform *tr = NULL;

  (void )AlgRandSeed(0);
  tr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, NULL);
  for(idR = 0; idR < repeats; ++idR)
  {
   for(idN = 0; idN < 4; ++idN)
   {
     scale = MyRandom();
     sVx[idN].vtX = scale * MyRandom();
     sVx[idN].vtY = scale * MyRandom();
     sVx[idN].vtZ = scale * MyRandom();
     scale = MyRandom();
     dVx[idN].vtX = scale * MyRandom();
     dVx[idN].vtY = scale * MyRandom();
     dVx[idN].vtZ = scale * MyRandom();
    }
    ss = 0.0;
    squashed = WlzGeomTetraAffineSolve(*(tr->mat), sVx, dVx, delta);
    for(idN = 0; idN < 4; ++idN)
    {
      tVx[idN] = WlzAffineTransformVertexD3(tr, sVx[idN], NULL);
      WLZ_VTX_3_SUB(eVx, tVx[idN], dVx[idN]);
      ss += WLZ_VTX_3_SQRLEN(eVx);
    }
    if(ss > DBL_EPSILON)
    {
      ss = sqrt(ss);
    }
    printf("%d %012e | %012e %012e %012e | %012e %012e %012e | "
           "%012e %012e %012e|\n",
           squashed, ss,
	   sVx[0].vtX, sVx[0].vtY, sVx[0].vtZ, 
	   dVx[0].vtX, dVx[0].vtY, dVx[0].vtZ,
	   tVx[0].vtX - dVx[0].vtX,
	   tVx[0].vtY - dVx[0].vtY,
	   tVx[0].vtZ - dVx[0].vtZ);
    for(idN = 1; idN < 4; ++idN)
    {
      printf("               | %012e %012e %012e | %012e %012e %012e | "
	     "%012e %012e %012e|\n",
	     sVx[idN].vtX, sVx[idN].vtY, sVx[idN].vtZ, 
	     dVx[idN].vtX, dVx[idN].vtY, dVx[idN].vtZ,
	     tVx[idN].vtX - dVx[idN].vtX,
	     tVx[idN].vtY - dVx[idN].vtY,
	     tVx[idN].vtZ - dVx[idN].vtZ);
    }
    for(idN = 0; idN < 4; ++idN)
    {
      printf("%012e %012e %012e %012e\n",
	     tr->mat[idN][0], tr->mat[idN][1],
	     tr->mat[idN][2], tr->mat[idN][3]);
    }
  }
  (void )WlzFreeAffineTransform(tr);
}

