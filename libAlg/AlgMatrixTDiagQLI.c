#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgMatrixTDiagQLI.c
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
* \brief        Determines the eigenvalues and eigenvectors of a
*		real symmetric tridiagonal matrix using implicit shifts.
* \ingroup      AlgMatrix
* \bug          None known.
* \todo
*/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

/*!
* \return       Error code.
* \ingroup      AlgMatrix
* \brief        Determines the eigenvalues and eigenvectors of a
*		real symmetric tridiagonal matrix using the QL
*		algorithm with implicit shifts.
*		This is a based on the function 'tqli' in Numerical
*		Recipies in C: The Art of Scientific Computiung.
*		Cambridge University Press, 1992.

* \param        dM 			Given diagonal elements of the
*					tridiagonal matrix. On return
*					it contains the eigenvalues.
* \param	oM 			Given off diagonal elements of the
*					tridiagonal matrix, with an
*					arbitrary first element. On return
*					it's contents are destroyed.
* \param        aSz			Size of the array.
* \param	zM 			Given array for the return of the
* 					eigenvectors, may be NULL if the
*					eigenvectors are not required.
*					If required the i'th eigenvector
*					is returned in the i'th column
*					of zM.
*/
AlgError	AlgMatrixTDiagQLI(double *dM, double *oM, int aSz, double **zM)
{
  int 		id0,
  		id1,
  		id2,
		id3,
		itr = 0;
  double 	b,
		c,
		dd,
		f,
		g,
		p,
  		r,
  		s;
  AlgError	errCode = ALG_ERR_NONE;
  const int	maxItr = 100;


  if((dM == NULL) || (oM == NULL) || (aSz < 2) || (zM && (*zM == NULL)))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    for(id0 = 1; id0 < aSz; ++id0)
    {
      oM[id0 - 1] = oM[id0];
    }
    oM[aSz - 1] = 0.0;
    for(id1 = 0; id1 < aSz; ++id1)
    {
      do
      {
	/* Find single small subdiagional element to split the matrix */
	id3 = id1;
	id2 = aSz - 2;
	while(id3 <= id2)
	{
	  dd = fabs(dM[id3]) + fabs(dM[id3 + 1]);
	  if((fabs(oM[id3]) + dd) == dd)
	  {
	    break;
	  }
	  ++id3;
	}
	if(id3 != id1)
	{
	  if(++itr > maxItr)
	  {
	    errCode = ALG_ERR_CONVERGENCE;
	  }
	  else
	  {
	    g = (dM[id1 + 1] - dM[id1]) / (2.0 * oM[id1]);
	    r = sqrt((g * g) + 1.0);
	    s = fabs(r);
	    g = dM[id3] - dM[id1] + (oM[id1] / (g + ((g < 0)? -s: s)));
	    s = c = 1.0;
	    p = 0.0;
	    /* Plane rotation followed by Givens rotations to restore
	     * tridiagonal form. */
	    for(id0 = id3 - 1; id0 >= id1; --id0)
	    {
	      f = s * oM[id0];
	      b = c * oM[id0];
	      if(fabs(f) >= fabs(g))
	      {
		c = g / f;
		r = sqrt((c * c) + 1.0);
		oM[id0 + 1] = f * r;
		c *= (s = 1.0 / r);
	      }
	      else
	      {
		s = f / g;
		r = sqrt((s * s) + 1.0);
		oM[id0 + 1] = g * r;
		s *= (c = 1.0 / r);
	      }
	      g = dM[id0 + 1] - p;
	      r = ((dM[id0] - g) * s) + (2.0 * c * b);
	      p = s * r;
	      dM[id0 + 1] = g + p;
	      g = (c * r) - b;
	      if(zM)
	      {
		for(id2 = 0 ; id2 < aSz; ++id2)
		{
		  f = zM[id2][id0 + 1];
		  zM[id2][id0 + 1] = s * zM[id2][id0] + (c * f);
		  zM[id2][id0] = c * zM[id2][id0] - (s * f);
		}
	      }
	    }
	    dM[id1] -= p;
	    oM[id1] = g;
	    oM[id3] = 0.0;
	  }
	}
      } while((errCode == ALG_ERR_NONE) && (id3 != id1));
    }
  }
  return(errCode);
}
