#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixRSEigen_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixRSEigen.c
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
* \brief        Functions to find the eigenvalues and eigenvectors of a
*		real symmetric matrix.
* \ingroup      AlgMatrix
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

static void			AlgMatrixRSEigenSort(
				  double **vM,
				  double *xM,
				  int n,
				  int reqEV);

/*!
* \return       Error code.
* \ingroup      AlgMatrix
* \brief        Determines the eigenvalues and eigenvectors of a
*		real symmetric matrix by calling AlgMatrixRSTDiag()
*		to create a tridiagonal symmetric matrix and then
*		AlgMatrixTDiagQLI() to compute its eigenvalues and
*		eigenvectors. The eigenvectors and eigenvalues 
*		are returned in descending eigenvalue order.
*		For efficiency, the eigenvectors should only be
*		computed if required.
* \param        aM 			Given real symmetric matrix
*					which contains the eigenvectors
*					in it's columns on return if
*					required.
* \param	vM 			Given vector for the return of the
* 					eigenvalues.
* \param	reqEV			Non zero if the eigenvectors are
*					required.
*/
AlgError	AlgMatrixRSEigen(AlgMatrix aM, double *vM, int reqEV)
{
  double	*oM = NULL;
  AlgError	errCode = ALG_ERR_NONE;


  if((aM.core == NULL) || (aM.core->type != ALG_MATRIX_RECT) ||
     (aM.core->nR <= 0) || (aM.core->nR != aM.core->nC) || (vM == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    if((oM = (double *)AlcMalloc(sizeof(double) * aM.core->nR)) == NULL)
    {
      errCode = ALG_ERR_MALLOC;
    }
    if(errCode == ALG_ERR_NONE)
    {
      if((errCode = AlgMatrixRSTDiag(aM, vM, oM)) == ALG_ERR_NONE)
      {
	AlgMatrix rM;

	rM.core = (reqEV == 0)? NULL: aM.core;
	errCode = AlgMatrixTDiagQLI(vM, oM, aM.core->nR, rM);
      }
    }
    if(errCode == ALG_ERR_NONE)
    {
      AlgMatrixRSEigenSort(aM.rect->array, vM, aM.core->nR, reqEV);
    }
    if(oM)
    {
      AlcFree(oM);
    }
  }
  return(errCode);
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief	Sorts the eigenvectors and eigenvalues into descending
* 		eigenvalue order. Because AlgMatrixRSTDiag() runs in
*		O(N^3) an O(N^2) insertion sort algorithm is used
*		without it significantly impacting on the execution
*		time, but it shouldn't be used elsewhere!
*		This function is based on 'eigsrt' in Numerical Recipies
*		in C: The Art of ScientificComputing. Cambridge University
*		Press, 1992.
*		in 
* \param	vM			The eigenvectors.
* \param	xM			The eigenvalues.
* \param	n			The number of eigenvectors or
* 					eigenvalues.
* \param	reqEV			Non zero if the eigenvectors are
*					required.
*/
static void	AlgMatrixRSEigenSort(double **vM, double *xM, int n, int reqEV)
{
  int		id0,
  		id1,
		id2,
		n1;
  double 	keyVal;

  n1 = n - 1;
  for(id0 = 0; id0 < n1; ++id0)
  {
    keyVal = xM[id2 = id0];
    for(id1 = id0 + 1; id1 < n; ++id1)
    {
      if(xM[id1] >= keyVal)
      {
	keyVal = xM[id2 = id1];
      }
    }
    if(id2 != id0)
    {
      xM[id2] = xM[id0];
      xM[id0] = keyVal;
      if(reqEV)
      {
	for(id1 = 0; id1 < n; id1++)
	{
	  keyVal = vM[id1][id0];
	  vM[id1][id0] = vM[id1][id2];
	  vM[id1][id2] = keyVal;
	}
      }
    }
  }
}
