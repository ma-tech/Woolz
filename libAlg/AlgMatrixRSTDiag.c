#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgMatrixRSTDiag.c
* \author       Bill Hill
* \date         May 2001
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Reduces a real symmetric matrix to symmetric
*		tridiagonal form by orthogonal similarity transformation
*		and construction of the right operator of the reduction.
* \bug          None known.
* \note
* Maintenance log with most recent changes at top of list.
*/

/*!
* \ingroup      AlgMatrix
* @{
*/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

/*!
* \return       	                  Error code.
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
* \param        aM 			Given real symetric matrix.
* \param        aSz 		    	Size of the array.
* \param	dM 			Given array for the return of the
* 					diagonal elements.
* \param	oM 			Given array for the return of the
* 					off diagonal elements.
*/
AlgError	AlgMatrixRSTDiag(double **aM, int aSz, double *dM, double *oM)
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

  if((aM == NULL) || (*aM == NULL) || (dM == NULL) || (oM == NULL) ||
     (aSz < 2))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    for(id0 = aSz - 1; id0 > 0; --id0)
    {
      id1 = id0 - 1;
      h = scale = 0.0;
      if(id1 > 0)
      {
	/* Do scaling for the computation of the L2 norm */
	for(id2 = 0; id2 <= id1; ++id2)
	{
	  scale += fabs(aM[id0][id2]);
	}
	/* Bypass the transformation if column is already in reduced form */
	if(scale < DBL_EPSILON)
	{
	  oM[id0] = aM[id0][id1];
	}
	else
	{
	  for(id2 = 0; id2 <= id1; ++id2)
	  {
	    aM[id0][id2]  /= scale;
	    h += aM[id0][id2] * aM[id0][id2];
	  }
	  f = aM[id0][id1];
	  g = (f > 0)? -sqrt(h): sqrt(h);
	  oM[id0] = scale * g;
	  h -= f * g;
	  aM[id0][id1] = f - g;
	  f = 0.0;
	  for(id3 = 0; id3 <= id1; ++id3)
	  {
	    aM[id3][id0] = aM[id0][id3] / h;
	    g = 0.0;
	    for(id2 = 0; id2 <= id3; ++id2)
	    {
	      g += aM[id3][id2] * aM[id0][id2];
	    }
	    for(id2 = id3 + 1; id2 <= id1; ++id2)
	    {
	      g += aM[id2][id3] * aM[id0][id2];
	    }
	    oM[id3] = g / h;
	    f += oM[id3] * aM[id0][id3];
	  }
	  hh = f / (h + h);
	  for(id3 = 0; id3 <= id1; ++id3)
	  {
	    f = aM[id0][id3];
	    oM[id3] = g = oM[id3] - (hh * f);
	    for(id2 = 0; id2 <= id3; ++id2)
	    {
	      aM[id3][id2] -= (f * oM[id2]) + (g * aM[id0][id2]);
	    }
	  }
	}
      }
      else
      {
	oM[id0] = aM[id0][id1];
      }
      dM[id0] = h;
    }
    dM[0] = 0.0;
    oM[0] = 0.0;
    /* Contents of this loop can be omitted if eigenvectors not
       wanted except for statement dM[id0] = aM[id0][id0]; */
    for(id0 = 0; id0 < aSz; ++id0)
    {
      id1 = id0 - 1;
      if(dM[id0])
      {
	for(id3 = 0; id3 <= id1; ++id3)
	{
	  g = 0.0;
	  for(id2 = 0; id2 <= id1; ++id2)
	  {
	    g += aM[id0][id2] * aM[id2][id3];
	  }
	  for(id2 = 0; id2 <= id1; ++id2)
	  {
	    aM[id2][id3] -= g * aM[id2][id0];
	  }
	}
      }
      dM[id0] = aM[id0][id0];
      /* Reset the row and column id3 of aM to the values of the identity#
       * matrix ready for the next itteration. */
      aM[id0][id0] = 1.0;
      for(id3 = 0; id3 <= id1; ++id3)
      {
	aM[id3][id0] = aM[id0][id3] = 0.0;
      }
    }
  }
  return(errCode);
}

/*!
* @}
*/
