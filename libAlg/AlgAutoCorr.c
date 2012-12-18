#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgAutoCorr_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgAutoCorr.c
* \author       Bill Hill
* \date         November 2007
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
* \brief	Frequency domain auto correlation functions.
* \ingroup	AlgCorr
*/

#include <Alg.h>
#include <math.h>
#include <float.h>

AlgError	AlgAutoCorrelate2D(double **data, int nX, int nY)
{
  int		tI0,
  		tI1,
		idX,
                idY,
                nX2,
                nY2;
  double        tD1,
                tD2,
                tD3,
                tD4;
  double        *tDP1,
  		*tDP2;
  AlgError      errNum = ALG_ERR_NONE;
  const int     minN = 8,
		maxN = 1048576;

  if((data == NULL) ||
     (nX < minN) || (nX > maxN) || (nY < minN) || (nY > maxN))
  {
    errNum = ALG_ERR_FUNC;
  }
  else
  {
    (void )AlgBitNextPowerOfTwo((unsigned int *)&tI0, nX);
    (void )AlgBitNextPowerOfTwo((unsigned int *)&tI1, nY);
    if((tI0 != nX) || (tI1 != nY))
    {
      errNum = ALG_ERR_FUNC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    AlgFourReal2D(data, 1, nX, nY);
    nX2 = nX / 2;
    nY2 = nY / 2;
    for(idY = 0; idY < nY; ++idY)
    {
      tDP1 = *(data + idY) + 1;
      for(idX = 1; idX < nX2; ++idX)
      {
        tD1 = *tDP1;
        tD2 = *(tDP1 + nX2);
        *tDP1 = (tD1 * tD1) + (tD2 * tD2);
        *(tDP1 + nX2) = 0.0;
        ++tDP1;
      }
    }
    for(idX = 0; idX < nX; idX += nX2)
    {
      for(idY = 1; idY < nY2; ++idY)
      {
        tDP1 = *(data + idY) + idX;
        tDP2 = *(data + nY2 + idY) + idX;
        tD1 = *tDP1;
        tD2 = *tDP2;
        tD3 = *(*(data + idY) + idX);
        tD4 = -*(*(data + nY2 + idY) + idX);
        *tDP1 = tD1 * tD3 - tD2 * tD4;
        *tDP2 = tD1 * tD4 + tD2 * tD3;
      }
    }
    **data *= **data;
    **(data + nY2) *= **(data + nY2);
    *(*data + nX2) *= *(*data + nX2);
    *(*(data + nY2) + nX2) *= *(*(data + nY2) + nX2);
    AlgFourRealInv2D(data, 1, nX, nY);
  }
  return(errNum);
}
