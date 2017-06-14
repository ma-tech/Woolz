#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixRSTDiag_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixRSTDiag.c
* \author       Bill Hill
* \date         May 2001
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
* \brief        Reduces a real symmetric matrix to symmetric
*		tridiagonal form by orthogonal similarity transformation
*		and construction of the right operator of the reduction.
* \ingroup      AlgMatrix
*/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

/*!
* \return       Error code.
* \ingroup      AlgMatrix
* \brief        An implementation of Householder's alorithm which
*		reduces a real aSz x aSz symmetric matrix to symmetric
*		tridiagonal form by aSz - 2 orthogonal similarity
*		transformations.
*		This is a translation of the functions/subroutines
*		'tred2' in the fortran EISPACK library and the
*		Numerical Recipies book, both of these are descended
*		from the algol procedure tred2 in Wilkinson J.H and
*		Reinsch C. Linear Algebra, Volume II Handbook for
*		Automatic Computation. Springer-Verlag, 1971.
*		See Numerical Recipies in C: The Art of Scientific
*		Computing. Cambridge University Press, 1992.
* \param        aMat 			Given real symmetric matrix.
* \param	dVec 			Given array for the return of the
* 					diagonal elements.
* \param	oVec 			Given array for the return of the
* 					off diagonal elements with the first
*					element set to 0.0.
*/
AlgError	AlgMatrixRSTDiag(AlgMatrix aMat, double *dVec, double *oVec)
{
  int 		id0,
  		id1,
  		id2,
		id3;
  double	scale,
  		hh,
		h,
		g,
		f;
  AlgError	errCode = ALG_ERR_NONE;

  if((aMat.core == NULL) || (aMat.core->type != ALG_MATRIX_RECT) ||
     (aMat.core->nR < 2) || (aMat.core->nR != aMat.core->nC) ||
     (dVec == NULL) || (oVec == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    int	   aSz;;
    double **aAry;

    aSz = aMat.rect->nR;
    aAry = aMat.rect->array;
    for(id0 = aSz - 1; id0 > 0; --id0)
    {
      id1 = id0 - 1;
      h = scale = 0.0;
      if(id1 > 0)
      {
	/* Do scaling for the computation of the L2 norm */
	for(id2 = 0; id2 <= id1; ++id2)
	{
	  scale += fabs(aAry[id0][id2]);
	}
	/* Bypass the transformation if column is already in reduced form */
	if(scale < DBL_EPSILON)
	{
	  oVec[id0] = aAry[id0][id1];
	}
	else
	{
	  for(id2 = 0; id2 <= id1; ++id2)
	  {
	    aAry[id0][id2]  /= scale;
	    h += aAry[id0][id2] * aAry[id0][id2];
	  }
	  f = aAry[id0][id1];
	  g = (f > 0)? -sqrt(h): sqrt(h);
	  oVec[id0] = scale * g;
	  h -= f * g;
	  aAry[id0][id1] = f - g;
	  f = 0.0;
	  for(id3 = 0; id3 <= id1; ++id3)
	  {
	    aAry[id3][id0] = aAry[id0][id3] / h;
	    g = 0.0;
	    for(id2 = 0; id2 <= id3; ++id2)
	    {
	      g += aAry[id3][id2] * aAry[id0][id2];
	    }
	    for(id2 = id3 + 1; id2 <= id1; ++id2)
	    {
	      g += aAry[id2][id3] * aAry[id0][id2];
	    }
	    oVec[id3] = g / h;
	    f += oVec[id3] * aAry[id0][id3];
	  }
	  hh = f / (h + h);
	  for(id3 = 0; id3 <= id1; ++id3)
	  {
	    f = aAry[id0][id3];
	    oVec[id3] = g = oVec[id3] - (hh * f);
	    for(id2 = 0; id2 <= id3; ++id2)
	    {
	      aAry[id3][id2] -= (f * oVec[id2]) + (g * aAry[id0][id2]);
	    }
	  }
	}
      }
      else
      {
	oVec[id0] = aAry[id0][id1];
      }
      dVec[id0] = h;
    }
    dVec[0] = 0.0;
    oVec[0] = 0.0;
    /* Contents of this loop can be omitted if eigenvectors not
       wanted except for statement dVec[id0] = aMat[id0][id0]; */
    for(id0 = 0; id0 < aSz; ++id0)
    {
      id1 = id0 - 1;
      if(fabs(dVec[id0]) > DBL_EPSILON)
      {
	for(id3 = 0; id3 <= id1; ++id3)
	{
	  g = 0.0;
	  for(id2 = 0; id2 <= id1; ++id2)
	  {
	    g += aAry[id0][id2] * aAry[id2][id3];
	  }
	  for(id2 = 0; id2 <= id1; ++id2)
	  {
	    aAry[id2][id3] -= g * aAry[id2][id0];
	  }
	}
      }
      dVec[id0] = aAry[id0][id0];
      /* Reset the row and column id3 of aMat to the values of the identity
       * matrix ready for the next itteration. */
      aAry[id0][id0] = 1.0;
      for(id3 = 0; id3 <= id1; ++id3)
      {
	aAry[id3][id0] = aAry[id0][id3] = 0.0;
      }
    }
  }
  return(errCode);
}
