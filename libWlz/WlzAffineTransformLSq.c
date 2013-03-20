#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzAffineTransformLSq_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzAffineTransformLSq.c
* \author       Bill Hill, Richard Baldock
* \date         March 1999
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
* \brief	Functions for computing Woolz affine transforms that
* 		give the best fit, in a least squares sense, when
* 		used to transform one set of vertices to another.
* \ingroup	WlzTransform
*/
#include <stdio.h>
#include <float.h>
#include <Wlz.h>


static WlzErrorNum 		WlzAffineTransformLSqLinSysSolve(
				  AlgMatrix aM,
				  double *bV,
				  double tol);

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		source vertices onto the target vertices.
*		The weights are optional but if givn must correspond
*		to the vertices.
*		This function calls either WlzAffineTransformLSq2D()
*		or WlzAffineTransformLSq3D() depending on the given
*		vertex type, see these functions of greater detail.
* \param	vType			Type of vertices.
* \param	nVT			Number of target vertices.
* \param	vT			Target vertices.
* \param	nVS			Number source vertices, which
*					must be the same as the number
*					of target vertices.
* \param	vS			Target vertices.
* \param	nVW			number of vertex weights, which
*					must either be zero or the same
*					as the number target vertices.
* \param	vW			Vertex pair weights which may
*					be NULL if the number of weights
*					is zero.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSq(WlzVertexType vType,
				          int nVT, WlzVertexP vT,
					  int nVS, WlzVertexP vS,
					  int nVW, double *vW,
					  WlzTransformType trType,
					  WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nVT != nVS) || (nVT <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vT.v == NULL) || (vS.v == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(vType)
    {
      case WLZ_VERTEX_D2:
	tr = WlzAffineTransformLSq2D(nVT, vT.d2, nVS, vS.d2, nVW, vW,
				     trType, &errNum);
	break;
      case WLZ_VERTEX_D3:
	tr = WlzAffineTransformLSq3D(nVT, vT.d3, nVS, vS.d3, nVW, vW,
				     trType, &errNum);
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the affine transform which gives the
*		best (least squares) fit when used to transform the
*		2D source vertices onto the target vertices.
*		This function calls the appropriate 2D least
*		squares affine transform function:
*		<ul>
*		<li> WlzAffineTransformLSqTrans2D() - if the requested
*		     transform type is WLZ_TRANSFORM_2D_TRANS or
*		     only a single vertex pair is given.
*		<li> WlzAffineTransformLSqReg2D() - if the requested
*		     transform type is WLZ_TRANSFORM_2D_REG,
*		     WLZ_TRANSFORM_2D_NOSHEAR or the number of vertex
*		     pairs is \f$<\f$ 3.
*		<li> WlzAffineTransformLSqGen2D() - all other cases.
*		</ul>
* \param	nVT			Number of target vertices.
* \param	vT			First vector of vertices.
* \param	nVTS			Number of vertices in second
*					vector (MUST be same as first).
* \param	vS			Second vector of vertices.
* \param	nVW			number of vertex weights, which
*					must either be zero or the same
*					as the number target vertices.
* \param	vW			Vertex pair weights which may
*					be NULL if the number of weights
*					is zero.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSq2D(int nVT, WlzDVertex2 *vT,
					    int nVTS, WlzDVertex2 *vS,
					    int nVW, double *vW,
					    WlzTransformType trType,
					    WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(nVW == 0)
  {
    vW = NULL;
  }
  if((nVT != nVTS) || (nVT <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vT == NULL) || (vS == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(trType)
    {
      case WLZ_TRANSFORM_2D_TRANS:
      case WLZ_TRANSFORM_2D_REG:
      case WLZ_TRANSFORM_2D_NOSHEAR:
      case WLZ_TRANSFORM_2D_AFFINE:
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((trType == WLZ_TRANSFORM_2D_TRANS) || (nVT == 1))
    {
      tr = WlzAffineTransformLSqTrans2D(vT, vS, vW, nVT, &errNum);
    }
    else if((trType == WLZ_TRANSFORM_2D_REG))
    {
      tr = WlzAffineTransformLSqReg2D(vT, vS, vW, nVT, &errNum);
    }
    else if((trType == WLZ_TRANSFORM_2D_NOSHEAR) || (nVT < 3))
    {
      WlzAffineTransform	*tr1, *tr2;
      WlzDVertex2		*vS1;
      int			i;

      /* kludge for now - scale then rigid body */
      if((tr1 = WlzAffineTransformLSqScale2D(vT, vS, vW, nVT,
                                             &errNum)) != NULL ){
	/* apply to the source vertices */
	if((vS1 = (WlzDVertex2 *) AlcMalloc(sizeof(WlzDVertex2) *
	                                    nVT)) != NULL ){
	  for(i=0; (i < nVT) && (errNum == WLZ_ERR_NONE); i++){
	    vS1[i] = WlzAffineTransformVertexD2(tr1, vS[i], &errNum);
	  }
	  /* now rigid body from re-scaled vertices */
	  if( errNum == WLZ_ERR_NONE ){
	    if((tr2 = WlzAffineTransformLSqReg2D(vT, vS1, vW, nVT,
	                                         &errNum)) != NULL){
	      tr = WlzAffineTransformProduct(tr1, tr2, &errNum);
	      WlzFreeAffineTransform(tr2);
	    }
	  }
	  AlcFree(vS1);
	}
	else {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	WlzFreeAffineTransform(tr1);
      }
    }
    else
    {
      tr = WlzAffineTransformLSqGen2D(vT, vS, vW, nVT, &errNum);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the affine transform which gives the
*		best (least squares) fit when used to transform the
*		3D source vertices onto the target vertices.
*		This function calls the appropriate 3D least
*		squares affine transform function:
*		<ul>
*		<li> WlzAffineTransformLSqTrans3D() - if the requested
*		     transform type is WLZ_TRANSFORM_3D_TRANS or
*		     only a single vertex pair is given.
*		<li> WlzAffineTransformLSqReg3D() - if the requested
*		     transform type is WLZ_TRANSFORM_3D_REG, or the
*		     number of vertex pairs is \f$<\f$ 6.
*		<li> WlzAffineTransformLSqGen3D() - all other cases.
*		</ul>
* \param	nVT			Number of target vertices.
* \param	vT			First vector of vertices.
* \param	nVS			Number of vertices in second
*					vector (MUST be same as first).
* \param	vS			Second vector of vertices.
* \param	nVW			number of vertex weights, which
*					must either be zero or the same
*					as the number target vertices.
* \param	vW			Vertex pair weights which may
*					be NULL if the number of weights
*					is zero.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSq3D(int nVT, WlzDVertex3 *vT,
					    int nVS, WlzDVertex3 *vS,
					    int nVW, double *vW,
					    WlzTransformType trType,
					    WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(nVW == 0)
  {
    vW = NULL;
  }
  if((nVT != nVS) || (nVT <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vT == NULL) || (vS == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(trType)
    {
      case WLZ_TRANSFORM_3D_TRANS:
      case WLZ_TRANSFORM_3D_REG:
      case WLZ_TRANSFORM_3D_AFFINE:
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((trType == WLZ_TRANSFORM_3D_TRANS) || (nVT == 1))
    {
      tr = WlzAffineTransformLSqTrans3D(vT, vS, vW, nVT, &errNum);
    }
    else if((trType == WLZ_TRANSFORM_3D_REG) || (nVT < 6))
    {
      tr = WlzAffineTransformLSqReg3D(vT, vS, vW, nVT, &errNum);
    }
    else
    {
      tr = WlzAffineTransformLSqGen3D(vT, vS, vW, nVT, &errNum);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the 2D translation transform which gives the
*		best (least squares) fit when used to transform the
*		source vertices onto the target vertices.
* \param	vT			Target vertices.
* \param	vS			Source vertices.
* \param	vW			Vertex pair weights.
* \param	nV			Number of vertex pairs.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqTrans2D(WlzDVertex2 *vT,
		    		WlzDVertex2 *vS, double *vW, int nV,
				WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex2	sum;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(nV <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vT == NULL) || (vS == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if(vW)
    {
      sum.vtX = (vT->vtX - vS->vtX) * *vW;
      sum.vtY = (vT->vtY - vS->vtY) * *vW;
      for(idx = 1; idx < nV; ++idx)
      {
	++vT;
	++vS;
	++vW;
	sum.vtX += (vT->vtX - vS->vtX) * *vW;
	sum.vtY += (vT->vtY - vS->vtY) * *vW;
      }
    }
    else
    {
      sum.vtX = vT->vtX - vS->vtX;
      sum.vtY = vT->vtY - vS->vtY;
      for(idx = 1; idx < nV; ++idx)
      {
	++vT;
	++vS;
	sum.vtX += vT->vtX - vS->vtX;
	sum.vtY += vT->vtY - vS->vtY;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tr = WlzAffineTransformFromTranslation(WLZ_TRANSFORM_2D_AFFINE,
				           sum.vtX / nV, sum.vtY / nV, 0.0,
				           &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the 3D translation transform which gives the
*		best (least squares) fit when used to transform the
*		source vertices onto the target vertices.
* \param	vT			Target vertices.
* \param	vS			Source vertices.
* \param	vW			Vertex pair weights.
* \param	nV			Number of vertex pairs.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform  *WlzAffineTransformLSqTrans3D(WlzDVertex3 *vT,
		    		WlzDVertex3 *vS, double *vW, int nV,
				WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex3	sum;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  sum.vtX = sum.vtY = sum.vtZ = 0.0;
  if(nV <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vT == NULL) || (vS == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if(vW)
    {
      sum.vtX = (vT->vtX - vS->vtX) * *vW;
      sum.vtY = (vT->vtY - vS->vtY) * *vW;
      sum.vtZ = (vT->vtZ - vS->vtZ) * *vW;
      for(idx = 1; idx < nV; ++idx)
      {
	++vT;
	++vS;
	++vW;
	sum.vtX += (vT->vtX - vS->vtX) * *vW;
	sum.vtY += (vT->vtY - vS->vtY) * *vW;
	sum.vtZ += (vT->vtZ - vS->vtZ) * *vW;
      }
    }
    else
    {
      sum.vtX = vT->vtX - vS->vtX;
      sum.vtY = vT->vtY - vS->vtY;
      sum.vtZ = vT->vtZ - vS->vtZ;
      for(idx = 1; idx < nV; ++idx)
      {
	++vT;
	++vS;
	sum.vtX += vT->vtX - vS->vtX;
	sum.vtY += vT->vtY - vS->vtY;
	sum.vtZ += vT->vtZ - vS->vtZ;
      }
    }
  }
  tr = WlzAffineTransformFromTranslation(WLZ_TRANSFORM_3D_AFFINE,
				   sum.vtX / nV, sum.vtY / nV, sum.vtZ / nV,
					 &errNum);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the general 2D affine transform which
*		gives the best (least squares) fit when used to
*		transform the target vertices onto the source
*		vertices.
*
*		Given \f$N\f$ source, target vertex pairs
*		\f$\mathbf{x_i}\f$, \f$\mathbf{\acute{x_i}}\f$ respectively,
*		each with a weight \f$w_i\f$, the least
*		squares general affine transform \f$T\f$, with
*		\f[
		  \mathbf{\acute{x_i}} = T \mathbf{x_i}
		\f]
*		can be found by minimizing \f$\chi^2\f$ in:
*		\f[
		  \chi^2 = \sum_{i=0}^{N-1}
		           {(w_i(\mathbf{x_i} - \mathbf{\acute{x_i}})^2)}
	        \f]
*		By taking partial derivatives of \f$\chi^2\f$ w.r.t.
*		the transform array elements and setting these to zero
*		the least squares affine transform is found through solving:
*		\f[
		  A \mathbf{x} = \mathbf{b}
		\f]
*		with
*		\f[
		A = \left( \begin{array}{ccc}
		    \sum_{i=0}^{N-1}{w_i^2 x_i^2} &
		      \sum_{i=0}^{N-1}{w_i^2 x_i y_i} &
			\sum_{i=0}^{N-1}{w_i^2 x_i} \\
		    \sum_{i=0}^{N-1}{w_i^2 x_i y_i} &
		      \sum_{i=0}^{N-1}{w_i^2 y_i^2} &
			\sum_{i=0}^{N-1}{w_i^2 y_i} \\
		    \sum_{i=0}^{N-1}{w_i^2 x_i} &
		      \sum_{i=0}^{N-1}{w_i^2 y_i} &
		        \sum_{i=0}^{N-1}{w_i^2}
		    \end{array} \right)
		\f]
*		For \f$\mathbf{x} = (t00, t01, t02)^T\f$:
*		\f[
		  \mathbf{b} = \left( \begin{array}{c}
				 \sum_{i=0}^{N-1}{w_i^2 x_i \acute{x_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 y_i \acute{x_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 \acute{x_i}}
		               \end{array} \right)
		\f]
*		and for \f$\mathbf{y} = (t10, t11, t12)^T\f$:
*		\f[
		  \mathbf{b} = \left( \begin{array}{c}
				 \sum_{i=0}^{N-1}{w_i^2 x_i \acute{y_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 y_i \acute{y_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 \acute{y_i}}
		               \end{array} \right)
		\f]
* \param	vT			Target vertices.
* \param	vS			Source vertices.
* \param	vW			Vertex pair weights in range
*					[0.0-1.0], may be NULL in which case
*					all weights have value 1.0.
* \param	nV			Number of vertex pairs.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqGen2D(WlzDVertex2 *vT,
				WlzDVertex2 *vS, double *vW, int nV,
				WlzErrorNum *dstErr)
{
  int		idx;
  double	wSq;
  double	**aA,
		**trA;
  double	bV[4],
  		sums[12];
  AlgMatrix	aM,
  		trM;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  aM.core = NULL;
  trM.core = NULL;
  if(nV <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vT == NULL) || (vS == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    /* Initialise the array */
    for(idx = 0; idx < 12; ++idx)
    {
      sums[idx] = 0;
    }
    /* Accumulate values */
    if(vW)
    {
      for(idx = 0; idx < nV; ++idx)
      {
	wSq = vW[idx] * vW[idx];
	sums[0]  += vS[idx].vtX * vS[idx].vtX * wSq;
	sums[1]  += vS[idx].vtX * vS[idx].vtY * wSq;
	sums[2]  += vS[idx].vtX * wSq;
	sums[3]  += vS[idx].vtY * vS[idx].vtY * wSq;
	sums[4]  += vS[idx].vtY * wSq;
	sums[5]  += wSq;
	sums[6]  += vS[idx].vtX * vT[idx].vtX * wSq;
	sums[7]  += vS[idx].vtY * vT[idx].vtX * wSq;
	sums[8]  += vT[idx].vtX * wSq;
	sums[9]  += vS[idx].vtX * vT[idx].vtY * wSq;
	sums[10] += vS[idx].vtY * vT[idx].vtY * wSq;
	sums[11] += vT[idx].vtY * wSq;
      }
    }
    else
    {
      sums[5]  = (double )nV;
      for(idx = 0; idx < nV; ++idx)
      {
	sums[0]  += vS[idx].vtX * vS[idx].vtX;
	sums[1]  += vS[idx].vtX * vS[idx].vtY;
	sums[2]  += vS[idx].vtX;
	sums[3]  += vS[idx].vtY * vS[idx].vtY;
	sums[4]  += vS[idx].vtY;
	sums[6]  += vS[idx].vtX * vT[idx].vtX;
	sums[7]  += vS[idx].vtY * vT[idx].vtX;
	sums[8]  += vT[idx].vtX;
	sums[9]  += vS[idx].vtX * vT[idx].vtY;
	sums[10] += vS[idx].vtY * vT[idx].vtY;
	sums[11] += vT[idx].vtY;
      }
    }
    /* Allocate workspace */
    if(((aM.rect = AlgMatrixRectNew(3, 3, NULL)) == NULL) ||
       ((trM.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    trA = trM.rect->array;
    /* Determine the least square transformation matrix values */
    /* x parameters first */
    aA[0][0] = sums[0];  aA[0][1] = sums[1];  aA[0][2] = sums[2];
    aA[1][0]  = sums[1];  aA[1][1] = sums[3];  aA[1][2] = sums[4];
    aA[2][0]  = sums[2];  aA[2][1] = sums[4];  aA[2][2] = sums[5];
    bV[0] = sums[6];     bV[1] = sums[7];     bV[2] = sums[8];
    errNum = WlzAffineTransformLSqLinSysSolve(aM, bV, 0.000001);
    trA[0][0] = bV[0]; trA[0][1] = bV[1]; trA[0][2] = bV[2];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Now y parameters */
    aA[0][0] = sums[0];  aA[0][1] = sums[1];  aA[0][2] = sums[2];
    aA[1][0] = sums[1];  aA[1][1] = sums[3];  aA[1][2] = sums[4];
    aA[2][0] = sums[2];  aA[2][1] = sums[4];  aA[2][2] = sums[5];
    bV[0] = sums[9];     bV[1] = sums[10];    bV[2] = sums[11];
    errNum = WlzAffineTransformLSqLinSysSolve(aM, bV, 0.000001);
    trA[1][0] = bV[0]; trA[1][1] = bV[1]; trA[1][2] = bV[2];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make the transform */
    trA[2][0] = 0; trA[2][1] = 0; trA[2][2] = 1;
    tr = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, trA,
				      &errNum);
  }
  AlgMatrixFree(aM);
  AlgMatrixFree(trM);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the general 3D affine transform which
*		gives the best (least squares) fit when used to
*		transform the target vertices onto the source
*		vertices.
*
*		Given \f$N\f$ source, target vertex pairs
*		\f$\mathbf{x_i}\f$, \f$\mathbf{\acute{x_i}}\f$ respectively,
*		each with a weight \f$w_i\f$, the least
*		squares general affine transform \f$T\f$, with
*		\f[
		  \mathbf{\acute{x_i}} = T \mathbf{x_i}
		\f]
*		can be found by minimizing \f$\chi^2\f$ in:
*		\f[
		  \chi^2 = \sum_{i=0}^{N-1}
		           {(w_i(\mathbf{x_i} - \mathbf{\acute{x_i}})^2)}
	        \f]
*		By taking partial derivatives of \f$\chi^2\f$ w.r.t.
*		the transform array elements and setting these to zero
*		the least squares affine transform is found through solving:
*		\f[
		  A \mathbf{x} = \mathbf{b}
		\f]
*		with
*		\f[
		A = \left( \begin{array}{cccc}
		    \sum_{i=0}^{N-1}{w_i^2 x_i^2} &
		      \sum_{i=0}^{N-1}{w_i^2 x_i y_i} &
			\sum_{i=0}^{N-1}{w_i^2 x_i z_i} &
			  \sum_{i=0}^{N-1}{w_i^2 x_i} \\
		    \sum_{i=0}^{N-1}{w_i^2 x_i y_i} &
		      \sum_{i=0}^{N-1}{w_i^2 y_i^2} &
			\sum_{i=0}^{N-1}{w_i^2 y_i z_i} &
			  \sum_{i=0}^{N-1}{w_i^2 y_i} \\
		    \sum_{i=0}^{N-1}{w_i^2 x_i z_i} &
		      \sum_{i=0}^{N-1}{w_i^2 y_i z_i} &
			\sum_{i=0}^{N-1}{w_i^2 z_i^2} &
			  \sum_{i=0}^{N-1}{w_i^2 z_i} \\
		    \sum_{i=0}^{N-1}{w_i^2 x_i} &
		      \sum_{i=0}^{N-1}{w_i^2 y_i} &
			\sum_{i=0}^{N-1}{w_i^2 z_i} &
			  \sum_{i=0}^{N-1}{w_i^2}
		    \end{array} \right)
		\f]
*		For \f$\mathbf{x} = (t00, t01, t02, t03)^T\f$:
*		\f[
		  \mathbf{b} = \left( \begin{array}{c}
				 \sum_{i=0}^{N-1}{w_i^2 x_i \acute{x_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 y_i \acute{x_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 z_i \acute{x_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 \acute{x_i}}
		               \end{array} \right)
		\f]
*		for \f$\mathbf{y} = (t10, t11, t12, t13)^T\f$:
*		\f[
		  \mathbf{b} = \left( \begin{array}{c}
				 \sum_{i=0}^{N-1}{w_i^2 x_i \acute{y_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 y_i \acute{y_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 z_i \acute{y_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 \acute{y_i}}
		               \end{array} \right)
		\f]
*		and for \f$\mathbf{z} = (t20, t21, t22, t23)^T\f$:
*		\f[
		  \mathbf{b} = \left( \begin{array}{c}
				 \sum_{i=0}^{N-1}{w_i^2 x_i \acute{z_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 y_i \acute{z_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 z_i \acute{z_i}} \\
				 \sum_{i=0}^{N-1}{w_i^2 \acute{z_i}}
		               \end{array} \right)
		\f]
* \param	vT			Target vertices.
* \param	vS			Source vertices.
* \param	vW			Vertex pair weights in range
*					[0.0-1.0], may be NULL in which case
*					all weights have value 1.0.
* \param	nV			Number of vertex pairs.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqGen3D(WlzDVertex3 *vT,
				WlzDVertex3 *vS, double *vW, int nV,
				WlzErrorNum *dstErr)
{
  int		idx;
  double	wSq;
  double	**aA,
		**trA;
  double	bV[4],
  		sums[22];
  AlgMatrix	aM,
  		trM;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = 0.000000001;

  aM.core = NULL;
  trM.core = NULL;
  if(nV <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vT == NULL) || (vS == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    /* Initialise the array */
    for(idx = 0; idx < 22; ++idx)
    {
      sums[idx] = 0;
    }
    /* Accumulate values */
    if(vW)
    {
      for(idx = 0; idx < nV; ++idx)
      {
	wSq = vW[idx] * vW[idx];
	sums[0]  += vS[idx].vtX * vS[idx].vtX * wSq;
	sums[1]  += vS[idx].vtX * vS[idx].vtY * wSq;
	sums[2]  += vS[idx].vtX * vS[idx].vtZ * wSq;
	sums[3]  += vS[idx].vtX * wSq;
	sums[4]  += vS[idx].vtY * vS[idx].vtY * wSq;
	sums[5]  += vS[idx].vtY * vS[idx].vtZ * wSq;
	sums[6]  += vS[idx].vtY * wSq;
	sums[7]  += vS[idx].vtZ * vS[idx].vtZ * wSq;
	sums[8]  += vS[idx].vtZ * wSq;
	sums[9]  += wSq;
	sums[10] += vS[idx].vtX * vT[idx].vtX * wSq;
	sums[11] += vS[idx].vtY * vT[idx].vtX * wSq;
	sums[12] += vS[idx].vtZ * vT[idx].vtX * wSq;
	sums[13] += vT[idx].vtX * wSq;
	sums[14]  += vS[idx].vtX * vT[idx].vtY * wSq;
	sums[15]  += vS[idx].vtY * vT[idx].vtY * wSq;
	sums[16]  += vS[idx].vtZ * vT[idx].vtY * wSq;
	sums[17]  += vT[idx].vtY * wSq;
	sums[18]  += vS[idx].vtX * vT[idx].vtZ * wSq;
	sums[19]  += vS[idx].vtY * vT[idx].vtZ * wSq;
	sums[20]  += vS[idx].vtZ * vT[idx].vtZ * wSq;
	sums[21]  += vT[idx].vtZ * wSq;

      }
    }
    else
    {
      sums[9]  = (double )nV;
      for(idx = 0; idx < nV; ++idx)
      {
	sums[0]  += vS[idx].vtX * vS[idx].vtX;
	sums[1]  += vS[idx].vtX * vS[idx].vtY;
	sums[2]  += vS[idx].vtX * vS[idx].vtZ;
	sums[3]  += vS[idx].vtX;
	sums[4]  += vS[idx].vtY * vS[idx].vtY;
	sums[5]  += vS[idx].vtY * vS[idx].vtZ;
	sums[6]  += vS[idx].vtY;
	sums[7]  += vS[idx].vtZ * vS[idx].vtZ;
	sums[8]  += vS[idx].vtZ;
	sums[10] += vS[idx].vtX * vT[idx].vtX;
	sums[11] += vS[idx].vtY * vT[idx].vtX;
	sums[12] += vS[idx].vtZ * vT[idx].vtX;
	sums[13] += vT[idx].vtX;
	sums[14] += vS[idx].vtX * vT[idx].vtY;
	sums[15] += vS[idx].vtY * vT[idx].vtY;
	sums[16] += vS[idx].vtZ * vT[idx].vtY;
	sums[17] += vT[idx].vtY;
	sums[18] += vS[idx].vtX * vT[idx].vtZ;
	sums[19] += vS[idx].vtY * vT[idx].vtZ;
	sums[20] += vS[idx].vtZ * vT[idx].vtZ;
	sums[21] += vT[idx].vtZ;
      }
    }
    /* Allocate workspace */
    if(((aM.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL) ||
       ((trM.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    trA = trM.rect->array;
    /* Determine the least square transformation matrix values */
    /* x parameters first */
    aA[0][0] = sums[0];  aA[0][1] = sums[1];
    aA[0][2] = sums[2];  aA[0][3] = sums[3];
    aA[1][0] = sums[1];  aA[1][1] = sums[4];
    aA[1][2] = sums[5];  aA[1][3] = sums[6];
    aA[2][0] = sums[2];  aA[2][1] = sums[5];
    aA[2][2] = sums[7];  aA[2][3] = sums[8];
    aA[3][0] = sums[3];  aA[3][1] = sums[6];
    aA[3][2] = sums[8];  aA[3][3] = sums[9];
    bV[0] =    sums[10]; bV[1] = sums[11];
    bV[2] =    sums[12]; bV[3] = sums[13];
    errNum = WlzAffineTransformLSqLinSysSolve(aM, bV, eps);
    trA[0][0] = bV[0]; trA[0][1] = bV[1];
    trA[0][2] = bV[2]; trA[0][3] = bV[3];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Now y parameters */
    aA[0][0] = sums[0];  aA[0][1] = sums[1];
    aA[0][2] = sums[2];  aA[0][3] = sums[3];
    aA[1][0] = sums[1];  aA[1][1] = sums[4];
    aA[1][2] = sums[5];  aA[1][3] = sums[6];
    aA[2][0] = sums[2];  aA[2][1] = sums[5];
    aA[2][2] = sums[7];  aA[2][3] = sums[8];
    aA[3][0] = sums[3];  aA[3][1] = sums[6];
    aA[3][2] = sums[8];  aA[3][3] = sums[9];
    bV[0] =    sums[14]; bV[1] = sums[15];
    bV[2] =    sums[16]; bV[3] = sums[17];
    errNum = WlzAffineTransformLSqLinSysSolve(aM, bV, eps);
    trA[1][0] = bV[0]; trA[1][1] = bV[1];
    trA[1][2] = bV[2]; trA[1][3] = bV[3];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Finaly z parameters */
    aA[0][0] = sums[0];  aA[0][1] = sums[1];
    aA[0][2] = sums[2];  aA[0][3] = sums[3];
    aA[1][0] = sums[1];  aA[1][1] = sums[4];
    aA[1][2] = sums[5];  aA[1][3] = sums[6];
    aA[2][0] = sums[2];  aA[2][1] = sums[5];
    aA[2][2] = sums[7];  aA[2][3] = sums[8];
    aA[3][0] = sums[3];  aA[3][1] = sums[6];
    aA[3][2] = sums[8];  aA[3][3] = sums[9];
    bV[0] =    sums[18]; bV[1] = sums[19];
    bV[2] =    sums[20]; bV[3] = sums[21];
    errNum = WlzAffineTransformLSqLinSysSolve(aM, bV, eps);
    trA[2][0] = bV[0]; trA[2][1] = bV[1];
    trA[2][2] = bV[2]; trA[2][3] = bV[3];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make the transform */
    trA[3][0] = 0; trA[3][1] = 0; trA[3][2] = 0; trA[3][3] = 1;
    tr = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_3D_AFFINE, trA,
				      &errNum);
  }
  AlgMatrixFree(aM);
  AlgMatrixFree(trM);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the 2D registration transform which
*		gives the best (weighted least squares) fit when used to
*		transform the first set of vertices onto the second
*		using the given match weights.
*		The transform is constrained to rotation and
*		translation only.
*		The algorithm has been addapted from the algorithm:
*		Arun K.S., Huang T.T. and Blostein S.D. "Least-Squares
*		Fitting of Two 3-D Point Sets" PAMI 9(5), 698-700, 1987,
*		by transforming the problem to a 2D space, including
*		vertex weights and computing the translation component
*		of the transform using
	        \f[
	          \mathbf{T} = \frac{\sum_{i=0}^{N-1}{w_i^2 \mathbf{x}_i}}
				    {\sum_{i=0}^{N-1}{w_i^2}}
	        \f]
	        \f[
	          \mathbf{x}_i = {\mathbf{p}_i}' - \mathbf{R}\mathbf{p}_i
	        \f]
*	        where \f$w_i\f$ are the weights, \f${\mathbf{p}_i}'\f$ are
*		the target vertices, \f$\mathbf{p}_i\f$ are the source
*		vertices and \f$\mathbf{R}\f$ is the rotation matrix.
* \param	vT			Target vertices.
* \param	vS			Source vertices.
* \param	vW			Vertex pair weights, range [0-1],
*					may be NULL in which case all
*					weights are 1.
* \param	nVtx			Number of vertices in each array.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqReg2D(WlzDVertex2 *vT,
		    		WlzDVertex2 *vS, double *vW, int nVtx,
				WlzErrorNum *dstErr)
{
  int		tI0,
  		idN;
  double	tD0,
  		meanSqD;
  WlzDVertex2	p0,
  		p1,
		cen0,
  		cen1,
		rel0,
		rel1;
  double	wV[2];
  double	**hA,
		**vA,
  		**trA = NULL;
  AlgMatrix	hM,
  		vM;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0E-06;

  wV[0] = wV[1] = 0.0;
  hM.core = vM.core = NULL;
  if(((hM.rect = AlgMatrixRectNew(2, 2, NULL)) == NULL) ||
     ((vM.rect = AlgMatrixRectNew(2, 2, NULL)) == NULL) ||
     (AlcDouble2Malloc(&trA, 4, 4) !=  ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    hA = hM.rect->array;
    vA = vM.rect->array;
    /* Compute weighted centroids (cen0 and cen1) and a mean of squares of
     * distance * between the weighted vertices. */
    if(vW)
    {
      WLZ_VTX_2_SCALE(p0, *vT, *vW);
      WLZ_VTX_2_SCALE(p1, *vS, *vW);
    }
    else
    {
      p0 = *vT;
      p1 = *vS;
    }
    cen0 = p0;
    cen1 = p1;
    WLZ_VTX_2_SUB(rel0, p0, p1);
    meanSqD = WLZ_VTX_2_DOT(rel0, rel0);
    for(idN = 1; idN < nVtx; ++idN)
    {
      if(vW)
      {
	WLZ_VTX_2_SCALE(p0, *(vT + idN), *(vW + idN));
	WLZ_VTX_2_SCALE(p1, *(vS + idN), *(vW + idN));
      }
      else
      {
        p0 = *(vT + idN);
	p1 = *(vS + idN);
      }
      WLZ_VTX_2_ADD(cen0, cen0, p0);
      WLZ_VTX_2_ADD(cen1, cen1, p1);
      WLZ_VTX_2_SUB(rel0, p0, p1);
      meanSqD += WLZ_VTX_2_DOT(rel0, rel0);
    }
    tD0 = 1.0 / nVtx;
    meanSqD *= tD0;
    if(meanSqD < tol)
    {
      /* The vertices coincide, make an identity transform. */
      tr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum);
    }
    else
    {
      WLZ_VTX_2_SCALE(cen0, cen0, tD0);
      WLZ_VTX_2_SCALE(cen1, cen1, tD0);
      /* Compute the 2x2 matrix hM, which is the sum of tensor products of
       * the weighted vertices relative to their centroids. */
      for(idN = 0; idN < nVtx; ++idN)
      {
	if(vW)
	{
	  WLZ_VTX_2_SCALE(p0, *(vT + idN), *(vW + idN));
	  WLZ_VTX_2_SCALE(p1, *(vS + idN), *(vW + idN));
	}
	else
	{
	  p0 = *(vT + idN);
	  p1 = *(vS + idN);
	}
	WLZ_VTX_2_SUB(rel0, p0, cen0);
	WLZ_VTX_2_SUB(rel1, p1, cen1);
	hA[0][0] += rel0.vtX * rel1.vtX;
	hA[0][1] += rel0.vtX * rel1.vtY;
	hA[1][0] += rel0.vtY * rel1.vtX;
	hA[1][1] += rel0.vtY * rel1.vtY;
      }
      /* Compute the SVD of the 2x2 matrix, hM = hM.wV.vM. */
      errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(hM, wV, vM));
      if(errNum == WLZ_ERR_NONE)
      {
	/* Compute 2x2 rotation matrix trA = vM'.hM, where vM' is the
	 * transpose of vM. */
	trA[0][0] = (vA[0][0] * hA[0][0]) + (vA[0][1] * hA[0][1]);
	trA[0][1] = (vA[1][0] * hA[0][0]) + (vA[1][1] * hA[0][1]);
	trA[1][0] = (vA[0][0] * hA[1][0]) + (vA[0][1] * hA[1][1]);
	trA[1][1] = (vA[1][0] * hA[1][0]) + (vA[1][1] * hA[1][1]);
	/* Test for degeneracy using the determinant of the rotation matrix. */
	tD0 = (trA[0][0] * trA[1][1]) - (trA[0][1] * trA[1][0]);
	if(tD0 < 0.0)
	{
	  /* Are source vertices (rel0) colinear? They are iff one of the 2
	   * singular values of hM in wV is zero. If the source vertices
	   * are not coplanar, the the solution of the SVD is correct. */
	  tI0 = (fabs(*(wV + 0)) >= DBL_EPSILON) |
		((fabs(*(wV + 1)) >= DBL_EPSILON) << 1);
	  if(tI0)
	  {
	    /* Source vertices are colinear and there exists an infinity of
	     * solutions. But select the identity rotation matrix */
	    trA[0][0] = trA[1][1] = trA[2][2] = 1.0;
	    trA[0][1] = trA[0][2] =
	    trA[1][0] = trA[1][2] =
	    trA[2][0] = trA[2][1] = 0.0;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* Fill in other matrix elements. */
	trA[2][0] = trA[2][1] = 0.0;
	trA[2][2] = 1.0;
	/* Compute the translation by applying the rotation to the source
	 * tie points before computing the translation \f$\mathbf{T}\f$
	 * using:
	 * \f[
	   \mathbf{T} = \frac{\sum_{i=0}^{N-1}{w_i^2 \mathbf{x}_i}}
	                     {\sum_{i=0}^{N-1}{w_i^2}}
	   \f]
	 * \f[
	   \mathbf{x}_i = {\mathbf{p}_i}' - \mathbf{R}\mathbf{p}_i
	   \f]
	 * where \f$w_i\f$ are the weights, \f${\mathbf{p}_i}'\f$ are the
	 * target vertices, \f$\mathbf{p}_i\f$ are the source vertices
	 and \f$\mathbf{R}\f$ is the rotation matrix.
	 */
        trA[0][2] = cen0.vtX - (cen1.vtX * trA[0][0] + cen1.vtY * trA[0][1]);
        trA[1][2] = cen0.vtY - (cen1.vtX * trA[1][0] + cen1.vtY * trA[1][1]);
	/* Build affine transform. */
	tr = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, trA,
	    				  &errNum);
      }
    }
  }
  AlgMatrixFree(hM);
  AlgMatrixFree(vM);
  (void )AlcDouble2Free(trA);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}


/*! 
* \return       Affine transform or NULL on error.
* \ingroup      WlzTransform
* \brief        Computes the 2D transform to rescale the source vertices.
* 		The assumption is that the source vertices have a different
* 		"spread" to the target and this re-scaling transform can be
* 		used in conjunction with the rigid-body (registration)
* 		tansform to determine a re-scaled shape-preserving transform
* 		i.e. no-shear. It is called by WlzAffineTransformLSq2D when
* 		transform type WLZ_TRANSFORM_2D_NOSHEAR is requested. The
* 		algorithm compares the mean distance from the centroid of
* 		each set of vertices.
* \param    	vT			Target vertices
* \param    	vS			Source vertices
* \param    	vW			Vertex weights
* \param    	nVtx			Number of vertices
* \param    	dstErr			Destination error pointer, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqScale2D(WlzDVertex2 *vT,
		    		WlzDVertex2 *vS, double *vW, int nVtx,
				WlzErrorNum *dstErr)
{
  int		idN;
  double	tD0;
  WlzDVertex2	p0,
  		p1,
		cen0,
                cen1;
  double	dist0, dist1;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0E-06;

  /* Compute weighted centroids (cen0 and cen1) */
  if(vW)
  {
    WLZ_VTX_2_SCALE(p0, *vT, *vW);
    WLZ_VTX_2_SCALE(p1, *vS, *vW);
  }
  else
  {
    p0 = *vT;
    p1 = *vS;
  }
  cen0 = p0;
  cen1 = p1;
  for(idN = 1; idN < nVtx; ++idN)
  {
    if(vW)
    {
      WLZ_VTX_2_SCALE(p0, *(vT + idN), *(vW + idN));
      WLZ_VTX_2_SCALE(p1, *(vS + idN), *(vW + idN));
    }
    else
    {
      p0 = *(vT + idN);
      p1 = *(vS + idN);
    }
    WLZ_VTX_2_ADD(cen0, cen0, p0);
    WLZ_VTX_2_ADD(cen1, cen1, p1);
  }
  tD0 = 1.0 / nVtx;
  WLZ_VTX_2_SCALE(cen0, cen0, tD0);
  WLZ_VTX_2_SCALE(cen1, cen1, tD0);
  /* Compute the sum of weighted distances from each centroid 
   * dist0 and dist1. */
  dist0 = 0.0;
  dist1 = 0.0;
  for(idN=0; idN < nVtx; idN++)
  {
    if(vW)
    {
      WLZ_VTX_2_SCALE(p0, *(vT + idN), *(vW + idN));
      WLZ_VTX_2_SCALE(p1, *(vS + idN), *(vW + idN));
    }
    else
    {
      p0 = *(vT + idN);
      p1 = *(vS + idN);
    }
    WLZ_VTX_2_SUB(p0, p0, cen0);
    WLZ_VTX_2_SUB(p1, p1, cen1);
    if((fabs(p0.vtX) > DBL_EPSILON) || (fabs(p0.vtY) > DBL_EPSILON))
    {
      dist0 += WLZ_VTX_2_LENGTH(p0);
    }
    if((fabs(p1.vtX) > DBL_EPSILON) || (fabs(p1.vtY) > DBL_EPSILON))
    {
      dist1 += WLZ_VTX_2_LENGTH(p1);
    }
  }
  /* Now make the transform */
  if(dist1 > tol)
  {
    tr = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
				       0.0, 0.0, 0.0,
				       dist0/dist1,
				       0.0, 0.0, 0.0, 0.0, 0.0,
				       0, &errNum);
  }
  else
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the Woolz 3D registration transform which
*		gives the best (weighted least squares) fit when used to
*		transform the first set of vertices onto the second
*		using the given match weights.
*		The transform is constrained to rotation and
*		translation only.
*		The algorithm has been addapted from the algorithm:
*		Arun K.S., Huang T.T. and Blostein S.D. "Least-Squares
*		Fitting of Two 3-D Point Sets" PAMI 9(5), 698-700, 1987,
*		by including vertex weights and computing the translation
*		component of the transform using
	        \f[
	          \mathbf{T} = \frac{\sum_{i=0}^{N-1}{w_i^2 \mathbf{x}_i}}
				    {\sum_{i=0}^{N-1}{w_i^2}}
	        \f]
	        \f[
	          \mathbf{x}_i = {\mathbf{p}_i}' - \mathbf{R}\mathbf{p}_i
	        \f]
*	        where \f$w_i\f$ are the weights, \f${\mathbf{p}_i}'\f$ are
*		the target vertices, \f$\mathbf{p}_i\f$ are the source
*		vertices and \f$\mathbf{R}\f$ is the rotation matrix.
* \param	vT			Target vertices.
* \param	vS			First array of vertices.
* \param	vW			Vertex pair weights, range [0-1],
*					may be NULL in which case all
*					weights are 1..
* \param	nV			Number of vertices in each array.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform  *WlzAffineTransformLSqReg3D(WlzDVertex3 *vT,
		    		WlzDVertex3 *vS, double *vW, int nV,
				WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
  		idK,
		idR;
  double	tD0,
  		meanSqD;
  WlzDVertex3	p0,
  		p1,
		cen0,
  		cen1,
		rel0,
		rel1;
  double	*wV = NULL;
  double	**hA,
		**vA,
  		**trM = NULL;
  AlgMatrix	hM,
  		vM;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0E-06;

  hM.core = vM.core = NULL;
  if(((hM.rect = AlgMatrixRectNew(3, 3, NULL)) == NULL) ||
     ((vM.rect = AlgMatrixRectNew(3, 3, NULL)) == NULL) ||
     ((wV = (double *)AlcCalloc(sizeof(double), 3)) == NULL) ||
     (AlcDouble2Malloc(&trM, 4, 4) !=  ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    hA = hM.rect->array;
    vA = vM.rect->array;
    /* Compute weighted centroids (cen0 and cen1) and a mean of squares of
     * distance * between the weighted vertices. */
    if(vW)
    {
      WLZ_VTX_3_SCALE(p0, *vT, *vW);
      WLZ_VTX_3_SCALE(p1, *vS, *vW);
    }
    else
    {
      p0 = *vT;
      p1 = *vS;
    }
    cen0 = p0;
    cen1 = p1;
    WLZ_VTX_3_SUB(rel0, p0, p1);
    meanSqD = WLZ_VTX_3_DOT(rel0, rel0);
    for(idN = 1; idN < nV; ++idN)
    {
      if(vW)
      {
	WLZ_VTX_3_SCALE(p0, *(vT + idN), *(vW + idN));
	WLZ_VTX_3_SCALE(p1, *(vS + idN), *(vW + idN));
      }
      else
      {
        p0 = *(vT + idN);
	p1 = *(vS + idN);
      }
      WLZ_VTX_3_ADD(cen0, cen0, p0);
      WLZ_VTX_3_ADD(cen1, cen1, p1);
      WLZ_VTX_3_SUB(rel0, p0, p1);
      meanSqD += WLZ_VTX_3_DOT(rel0, rel0);
    }
    tD0 = 1.0 / nV;
    meanSqD *= tD0;
    if(meanSqD < tol)
    {
      /* The vertices coincide, make an identity transform. */
      tr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, &errNum);
    }
    else
    {
      WLZ_VTX_3_SCALE(cen0, cen0, tD0);
      WLZ_VTX_3_SCALE(cen1, cen1, tD0);
      /* Compute the 3x3 matrix hM, which is the sum of tensor products of
       * the weighted vertices relative to their centroids. */
      for(idN = 0; idN < nV; ++idN)
      {
	if(vW)
	{
	  WLZ_VTX_3_SCALE(p0, *(vT + idN), *(vW + idN));
	  WLZ_VTX_3_SCALE(p1, *(vS + idN), *(vW + idN));
	}
	else
	{
	  p0 = *(vT + idN);
	  p1 = *(vS + idN);
	}
	WLZ_VTX_3_SUB(rel0, p0, cen0);
	WLZ_VTX_3_SUB(rel1, p1, cen1);
	hA[0][0] += rel0.vtX * rel1.vtX;
	hA[0][1] += rel0.vtX * rel1.vtY;
	hA[0][2] += rel0.vtX * rel1.vtZ;
	hA[1][0] += rel0.vtY * rel1.vtX;
	hA[1][1] += rel0.vtY * rel1.vtY;
	hA[1][2] += rel0.vtY * rel1.vtZ;
	hA[2][0] += rel0.vtZ * rel1.vtX;
	hA[2][1] += rel0.vtZ * rel1.vtY;
	hA[2][2] += rel0.vtZ * rel1.vtZ;
      }
      /* Compute the SVD of the 3x3 matrix, hM = hM.wV.vM. */
      errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(hM, wV, vM));
      if(errNum == WLZ_ERR_NONE)
      {
	/* Compute 3x3 rotation matrix trM = vM'.hM, where vM' is the
	 * transpose of vM. */
	trM[0][0] = (vA[0][0] * hA[0][0]) + (vA[0][1] * hA[0][1]) +
	            (vA[0][2] * hA[0][2]);
	trM[0][1] = (vA[1][0] * hA[0][0]) + (vA[1][1] * hA[0][1]) +
	            (vA[1][2] * hA[0][2]);
	trM[0][2] = (vA[2][0] * hA[0][0]) + (vA[2][1] * hA[0][1]) +
	            (vA[2][2] * hA[0][2]);
	trM[1][0] = (vA[0][0] * hA[1][0]) + (vA[0][1] * hA[1][1]) +
	            (vA[0][2] * hA[1][2]);
	trM[1][1] = (vA[1][0] * hA[1][0]) + (vA[1][1] * hA[1][1]) +
	            (vA[1][2] * hA[1][2]);
	trM[1][2] = (vA[2][0] * hA[1][0]) + (vA[2][1] * hA[1][1]) +
	            (vA[2][2] * hA[1][2]);
	trM[2][0] = (vA[0][0] * hA[2][0]) + (vA[0][1] * hA[2][1]) +
	            (vA[0][2] * hA[2][2]);
	trM[2][1] = (vA[1][0] * hA[2][0]) + (vA[1][1] * hA[2][1]) +
	            (vA[1][2] * hA[2][2]);
	trM[2][2] = (vA[2][0] * hA[2][0]) + (vA[2][1] * hA[2][1]) +
	            (vA[2][2] * hA[2][2]);
	/* Test for degeneracy using the determinant of the rotation matrix. */
	tD0 = (trM[0][0] * trM[1][1] * trM[2][2]) -
	      (trM[0][0] * trM[1][2] * trM[2][1]) +
	      (trM[0][1] * trM[1][2] * trM[2][0]) -
	      (trM[0][1] * trM[1][0] * trM[2][2]) +
	      (trM[0][2] * trM[1][0] * trM[2][1]) -
	      (trM[0][2] * trM[1][1] * trM[2][0]);
	if(tD0 < 0.0)
	{
	  /* Are source vertices (rel0) coplanar? They are iff one of the 3
	   * singular values of hM in wV is zero. If the source vertices
	   * are not coplanar, the the solution of the SVD is correct. */
	  tI0 = (fabs(*(wV + 0)) >= DBL_EPSILON) |
	        ((fabs(*(wV + 1)) >= DBL_EPSILON) << 1) |
	        ((fabs(*(wV + 2)) >= DBL_EPSILON) << 2);
	  switch(tI0)
	  {
	    case 0:
	      /* Source vertices are not coplanar or colinear. The SVD gives
	       * the correct rotation. */
	      break;
	    case 1:
	    case 2:
	    case 4:
	      /* Source vertices are coplanar, but not colinear. There is a
	       * unique reflection as well as a unique rotation. The SVD
	       * may give either BUT in this case it has found the reflection
	       * so need to recompute for the rotation.
	       * If the singular values of \f$\mathbf{W}\f$ are
	       * \f$w_0 > w_1 > w_2 = 0\f$, then
	       * \f[
	         \mathbf{H} = w_0 u_0 v_0^t + w_1 u_1 v_1^t + 0.u_2 v_2^t
		 \f]
	       * where \f$u_i\f$ and \f$v_i\f$ are the columns of the
	       * matricies \f$\mathbf{H}\f$ and \f$\mathbf{V}\f$
	       * respectively (hM is used for both \f$\mathbf{H}\f$
	       * and \f$\mathbf{U}\f$).
	       * So to get the rotation rather than the reflection we
	       * recalculate the transform matrix with \f$v_2 = -v_2\f$.
	       */
	      for(idR = 0; idR < 3; ++idR)
	      {
		for(idK = 0; idK < 3; ++idK)
		{
		  trM[idR][idK] = vA[idR][0] * hA[idK][0] + 
		    vA[idR][1] * hA[idK][1] -
		    vA[idR][2] * hA[idK][2];
		}
	      }
	      break;
	    case 3:
	    case 5:
	    case 6:
	      /* Source vertices are colinear and there exists an infinity of
	       * solutions. But select the identity rotation matrix */
	      trM[0][0] = trM[1][1] = trM[2][2] = 1.0;
	      trM[0][1] = trM[0][2] =
	      trM[1][0] = trM[1][2] =
	      trM[2][0] = trM[2][1] = 0.0;
	      break;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* Fill in other matrix elements. */
	trM[3][0] = trM[3][1] = trM[3][2] = 0.0;
	trM[3][3] = 1.0;
	/* Compute the translation by applying the rotation to the source
	 * tie points before computing the translation \f$\mathbf{T}\f$
	 * using:
	 * \f[
	   \mathbf{T} = \frac{\sum_{i=0}^{N-1}{w_i^2 \mathbf{x}_i}}
	                     {\sum_{i=0}^{N-1}{w_i^2}}
	   \f]
	 * \f[
	   \mathbf{x}_i = {\mathbf{p}_i}' - \mathbf{R}\mathbf{p}_i
	   \f]
	 * where \f$w_i\f$ are the weights, \f${\mathbf{p}_i}'\f$ are the
	 * target vertices, \f$\mathbf{p}_i\f$ are the source vertices
	 and \f$\mathbf{R}\f$ is the rotation matrix.
	 */
        trM[0][3] = cen0.vtX -
	           (cen1.vtX * trM[0][0] + cen1.vtY * trM[0][1] +
		    cen1.vtZ * trM[0][2]);
        trM[1][3] = cen0.vtY -
	           (cen1.vtX * trM[1][0] + cen1.vtY * trM[1][1] +
		    cen1.vtZ * trM[1][2]);
        trM[2][3] = cen0.vtZ -
	           (cen1.vtX * trM[2][0] + cen1.vtY * trM[2][1] +
		    cen1.vtZ * trM[2][2]);
	/* Build affine transform. */
	tr = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_3D_AFFINE, trM,
	    &errNum);
      }
    }
  }
  (void )AlcFree(wV);
  (void )AlcDouble2Free(trM);
  AlgMatrixFree(hM);
  AlgMatrixFree(vM);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return	Computed affine transform, may be NULL on error.
* \brief	Computes the Woolz 2D registration transform which
*		gives the best (least squares) fit when used to
*		transform the source vertices onto the target
*		vertices. This is an old function which should
*		not be used in new code, use WlzAffineTransformLSqReg2D()
*		instead.
* \param	vT			Target vertices.
* \param	vS			Source vertices.
* \param	nV			Number of vertex pairs.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqRegWlz2D(WlzDVertex2 *vT,
		    		WlzDVertex2 *vS, int nV,
				WlzErrorNum *dstErr)
{
  int		idx;
  double	s;
  double	**aA = NULL,
		**trA = NULL;
  double	*bV = NULL;
  WlzAffineTransform *tr = NULL;
  double	aV[12];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Initialise the array */
  for(idx=0; idx < 12; idx++)
  {
    aV[idx] = 0;
  }
  /* accumulate values */
  for(idx = 0; idx < nV; ++idx, ++vS, ++vT)
  {
    aV[0]  += 1;
    aV[1]  += vS->vtX;
    aV[2]  += vS->vtY;
    aV[3]  += vS->vtX * vS->vtX;
    /* aV[4]  += vS->vtX * vS->vtY;*/
    aV[5]  += vS->vtY * vS->vtY;
    aV[6]  += vT->vtX;
    aV[7]  += vS->vtX * vT->vtX;
    aV[8]  += vS->vtY * vT->vtX;
    aV[9]  += vT->vtY;
    aV[10] += vS->vtX * vT->vtY;
    aV[11] += vS->vtY * vT->vtY;
  }
  /* Allocate workspace */
  if((AlcDouble2Malloc(&aA, 4, 4) != ALC_ER_NONE) ||
      (AlcDouble1Malloc(&bV, 4) != ALC_ER_NONE) ||
      (AlcDouble2Malloc(&trA, 4, 4) != ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Determine the least square transformation matrix values */
    aA[0][0] = aV[0];          aA[0][1] = aV[1];
    aA[0][2] = aV[2];          aA[0][3] = 0;
    aA[1][0] = aV[1];          aA[1][1] = aV[3] + aV[5];
    aA[1][2] = 0;              aA[1][3] = aV[2];
    aA[2][0] = aV[2];          aA[2][1] = 0;
    aA[2][2] = aV[3] + aV[5];  aA[2][3] = -aV[1];
    aA[3][0] = 0;              aA[3][1] = aV[2];
    aA[3][2] = -aV[1];         aA[3][3] = aV[0];
    bV[0]    = aV[6];          bV[1]    = aV[7] + aV[11];
    bV[2]    = aV[8] - aV[10]; bV[3]    = aV[9];
    errNum = WlzErrorFromAlg(AlgMatrixLUSolveRaw4(aA, bV, 1));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check for scale constraint - this is a kludge */
    s = sqrt(bV[1]*bV[1] + bV[2]*bV[2]);
    bV[1] /= s;
    bV[2] /= s;
    bV[0] = (aV[6] - aV[1] * bV[1] - aV[2] * bV[2]) / aV[0];
    bV[3] = (aV[9] + aV[1] * bV[2] - aV[2] * bV[1]) / aV[0];
    /* Make the transformation */
    trA[0][0] = bV[1];  trA[0][1] = bV[2]; trA[0][2] = bV[0];
    trA[1][0] = -bV[2]; trA[1][1] = bV[1]; trA[1][2] = bV[3];
    trA[2][0] = 0;      trA[2][1] = 0;     trA[2][2] = 1;
    tr = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, trA,
					 &errNum);
  }
  (void )AlcDouble2Free(aA);
  (void )AlcFree((void *)bV);
  (void )AlcDouble2Free(trA);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \ingroup	WlzTransform
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of vertices and normals onto the second
*		set. The vertex and normal weighting factors must
*		be in the range [0-1].
*		See WlzAffineTransformLSqDQ3D() from which this
*		function has been derived.
*		This algorithm may be less stable than the SVD algorithm
*		particularly when the data are co--linear, use
*		WlzAffineTransformLSq2D() instead.
* \param	nV			Number of vertices.
* \param	vW			Vertex weights (Walker's beta),
*					may be NULL which implies that
*					all the weights are 1.0.
* \param	vT			Vertices of the target set.
* \param	vS			Vertices of the source set.
* \param	nN			Number of normals.
* \param	nW			Normal weights (Walker's alpha).,
*					may be NULL which implies that
*					all the weights are 1.0
* \param	nT			Target normals, may be NULL.
* \param	nS			Source normals, may be NULL.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqDQ2D(int nV, double *vW,
					WlzDVertex2 *vT, WlzDVertex2 *vS,
					int nN, double *nW,
					WlzDVertex2 *nT, WlzDVertex2 *nS,
					WlzErrorNum *dstErr)
{
  int		id0;
  double	wt,
		sumVW,
		t0D,
		t11D,
		t12D,
		t21D,
		t22D,
		t33D,
		t34D,
		t44D;
  double	rV[4];
  double	**aA,
  		**t0A,
  		**t1A,
		**c1A,
		**c3A;
  AlgMatrix	aM,	 		/* Walker's A */
  		t0M,	 		/* Working matrix */
  		t1M,	 		/* Working matrix */
  		t2M,	 		/* Working matrix */
		c1M,	 		/* Walker's C_1 */
		c3M;	 		/* Walker's C_3 */
  WlzDVertex3	t0V,
  		t1V;
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  aM.core = NULL;
  c1M.core = NULL;
  c3M.core = NULL;
  t0M.core = NULL;
  t1M.core = NULL;
  t2M.core = NULL;
  /* Allocate matricies required to compute c1M, c2M and c3M. */
  if(((aM.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL) ||
     ((c1M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL) ||
     ((c3M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL) ||
     ((t0M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL) ||
     ((t1M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL) ||
     ((t2M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  /* Compute c1M, c2M and c3M (see Walker's paper). */
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    c1A = c1M.rect->array;
    c3A = c3M.rect->array;
    t0A = t0M.rect->array;
    t1A = t1M.rect->array;
    if(vW)
    {
      sumVW = 0.0;
    }
    else
    {
      sumVW = 1.0;
    }
    for(id0 = 0; id0 < nV; ++id0)
    {
      /* Update sumVW */
      if(vW)
      {
        sumVW += vW[id0];
      }
      /* Update c1M: Compute \beta_i{Q(v_i^t)}^TW(v_i_s) */
      wt = (vW)? vW[id0]: 1.0;
      t0V.vtX = (vT + id0)->vtX; t0V.vtY = (vT + id0)->vtY; t0V.vtZ = 0.0;
      t1V.vtX = (vS + id0)->vtX; t1V.vtY = (vS + id0)->vtY; t1V.vtZ = 0.0;
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      c1A[0][0] += wt * ( t11D - t22D);
      c1A[1][1] += wt * (-t11D + t22D);
      c1A[2][2] += wt * (-t11D - t22D);
      c1A[3][3] += wt * ( t11D + t22D);
      t0D = wt * ( t12D + t21D); c1A[0][1] += t0D; c1A[1][0] += t0D;
      t0D = wt * (-t12D + t21D); c1A[2][3] += t0D; c1A[3][2] += t0D;
      /* Update c3M: Compute \beta_i(W(v_i_s) - Q(v_i^t)) */
      t0D = wt * (-t1V.vtY - t0V.vtY); c3A[0][2] += t0D; c3A[2][0] -= t0D; 
      t0D = wt * ( t1V.vtX - t0V.vtX); c3A[0][3] += t0D; c3A[3][0] -= t0D;
      t0D = wt * ( t1V.vtX + t0V.vtX); c3A[1][2] += t0D; c3A[2][1] -= t0D;
      t0D = wt * ( t1V.vtY - t0V.vtY); c3A[1][3] += t0D; c3A[3][1] -= t0D;
    }
    for(id0 = 0; id0 < nN; ++id0)
    {
      /* Update c1M: Compute \alpha_i{Q(n_i^t)}^TW(n_i_s) */
      wt = (nW)? nW[id0]: 1.0;
      t0V.vtX = (nT + id0)->vtX; t0V.vtY = (nT + id0)->vtY; t0V.vtZ = 0.0;
      t1V.vtX = (nS + id0)->vtX; t1V.vtY = (nS + id0)->vtY; t1V.vtZ = 0.0;
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      c1A[0][0] += wt * ( t11D - t22D);
      c1A[1][1] += wt * (-t11D + t22D);
      c1A[2][2] += wt * (-t11D - t22D);
      c1A[3][3] += wt * ( t11D + t22D);
      t0D = wt * ( t12D + t21D); c1A[0][1] += t0D; c1A[1][0] += t0D;
      t0D = wt * (-t12D + t21D); c1A[2][3] += t0D; c1A[3][2] += t0D;
    }
    AlgMatrixScale(c1M, c1M, -2.0);
    AlgMatrixScale(c3M, c3M, 2.0);
    if(sumVW < DBL_EPSILON)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Now compute A from C1, C3 and \sum_i^n{\beta_i} */
  if(errNum == WLZ_ERR_NONE)
  {
    AlgMatrixTranspose(t0M, c3M);			/* T0 = C3^T */
    AlgMatrixMul(t1M, t0M, c3M);	  		/* T1 = C3^T C3 */
    t0D = 1.0 / (4.0 * sumVW);
    AlgMatrixScale(t1M, t1M, t0D);
    AlgMatrixSub(aM, t1M, c1M);
    /* Find the eigenvector of A which has the greatest eigenvalue.
     * This is returned in the first column of aM. */
    errNum = WlzErrorFromAlg(AlgMatrixRSEigen(aM, rV, 1));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rV[0] = aA[0][0]; rV[1] = aA[1][0];
    rV[2] = aA[2][0]; rV[3] = aA[3][0];
    /* Compute the transform's rotation elements from the eigen vector. */
    t11D = rV[0] * rV[0];
    t22D = rV[1] * rV[1];
    t33D = rV[2] * rV[2];
    t44D = rV[3] * rV[3];
    t12D = 2.0 * rV[0] * rV[1];
    t34D = 2.0 * rV[2] * rV[3];
    aA[0][0] = t44D + t11D - t22D - t33D;
    aA[0][1] = t12D - t34D;
    aA[1][0] = t12D + t34D;
    aA[1][1] = t44D - t11D + t22D - t33D;
    /* Compute the translation elements from the eigen vector r, the sum of the
     * vertex weights \sum_i{\beta_i} and matrix C3. */
    /* Set t0M[0] to r (see Walker's paper). */
    t0A[0][0] = rV[0]; t0A[1][0] = rV[1]; t0A[2][0] = rV[2]; t0A[3][0] = rV[3];
    /* Compute s in t1M[0], but don't scale with -1/(2\sum_i{\beta_i} yet (see
     * Walker's paper). */
    t0M.rect->nC = 1; AlgMatrixMul(t1M, c3M, t0M); t0M.rect->nC = 4;
    /* Set t0M to be W^T(r) (see Walker's paper). */
    t0A[0][0] =  rV[3]; t0A[0][1] =  rV[2];
    t0A[0][2] =  rV[1]; t0A[0][3] = -rV[0];
    t0A[1][0] = -rV[2]; t0A[1][1] =  rV[3];
    t0A[1][2] =  rV[0]; t0A[1][3] = -rV[1];
    t0A[2][0] =  rV[1]; t0A[2][1] = -rV[0];
    t0A[2][2] =  rV[3]; t0A[2][3] = -rV[2];
    t0A[3][0] =  rV[0]; t0A[3][1] =  rV[1];
    t0A[3][2] =  rV[2]; t0A[3][3] =  rV[3];
    /* Compute the product W^T(r) s (see Walker's paper). */
    t1M.rect->nC = 1; AlgMatrixMul(t1M, t0M, t1M); t1M.rect->nC = 4;
    /* Extract the translation components and scale by -1/(2\sum_i{\beta_i}
     * (see Walker's paper). */
    t0D = -1.0 / (2.0 * sumVW);
    aA[0][2] = t0D * t1A[0][0];
    aA[1][2] = t0D * t1A[1][0];
    /* Set the transform's perspective and scale elements. */
    aA[2][0] = aA[2][1] = 0.0;
    aA[2][2] = 1.0;
    /* Create the affine transform. */
    trans = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, aA, &errNum);
  }
  /* Free the matricies. */
  AlgMatrixFree(aM);
  AlgMatrixFree(c1M);
  AlgMatrixFree(c3M);
  AlgMatrixFree(t0M);
  AlgMatrixFree(t1M);
  AlgMatrixFree(t2M);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trans);
}

/*!
* \ingroup	WlzTransform
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		source vertices and normals onto the target vertices.
*		The vertex and normal weighting factors must
*		be in the range [0-1].
*		This algorithm may be less stable than the SVD algorithm
*		particularly when the data are co--linear, use
*		WlzAffineTransformLSq3D() instead.
*
*		This function is based on the dual quaternion
*		algorithm: M.W. Walker  and Shao L. Estimating 3-D
*		Location Parameters Using Dual Number Quaternions,
*		CVGIP 54(3), 1991.
*		Given a set on \f$k\f$ source points
*		\f$\bar{\mathbf{p}}^0_i\f$, \f$l\f$ source unit
*		normals \f$\bar{\mathbf{n}}^0_i\f$;
*		the corresponding target points \f$\mathbf{p}^0_i\f$
*		and unit normals \f$\mathbf{n}^0_i\f$, with weights
*		\f$\alpha_i\f$ and \f$\beta_i\f$ the algorithm performs
*		the following steps:
*		<ol>
*		  <li> Compute \f$\mathbf{C}_1\f$, \f$\mathbf{C}_2\f$,
*	               \f$\mathbf{C}_3\f$:
*	               \f[
  \mathbf{C}_1 = -2 \sum_{i=1}^k{
    \alpha_i \mathbf{Q}(\bar{\mathbf{n}}_i)^T\mathbf{W}(\mathbf{n}^0_i)} -
                    2 \sum_{i=1}^l{
    \beta_i \mathbf{Q}(\bar{\mathbf{p}}_i)^T\mathbf{W}(\mathbf{p}^0_i)}
		       \f]
*	               \f[
  \mathbf{C}_2 = \left(\sum_{i=1}^l{\beta_i} \right)\mathbf{I}
		       \f]
*	               \f[
  \mathbf{C}_3 = 2 \sum_{i=1}^l{\beta_i
    (\mathbf{W}(\mathbf{p}^0_i) - \mathbf{Q}(\bar{\mathbf{p}}_i)) }
		       \f]
*		  <li> Compute the \f$4 \times 4\f$ symetric matrix
*		       \f$\mathbf{A}\f$:
*		       \f[
  \mathbf{A} = \frac{1}{2}
    (\mathbf{C}_3^T(\mathbf{C}_2 + \mathbf{C}_2^T)^{-1}\mathbf{C}_3 -
     \mathbf{C}_1 - \mathbf{C}_1^T)
		       \f]
*                     This can be simplified to:
*		       \f[
		 \mathbf{A} = {\frac{1}{4\sum_{i=1}^{n}{\beta_i}}}
			      {\mathbf{C}_3^T\mathbf{C}_3} -
			      \mathbf{C}_1
		       \f]
*    		      because \f$\mathbf{C}_1\f$ is real and symetric,
*		      \f$\mathbf{C}_2\f$ is scalar and \f$\mathbf{C}_3\f$
*		      is anti-symetric.
*		  <li> Compute the eigenvector of \f$r\f$ corresponding to
*		       the largest positive eigenvalue of matrix
*		       \f$\mathbf{A}\f$ and derive \f$\mathbf{s}\f$ from
*		       \f$\mathbf{r}\f$ using
*		       \f[
  \mathbf{t} = \mathbf{W}(\mathbf{r})^T \mathbf{s}
		       \f]
*		   <li> Compute the translation vector \f$\mathbf{t}\f$ and
*			rotation matrix \f$\mathbf{R}\f$ using:
*		       \f[
  \mathbf{t} = \mathbf{W}(\mathbf{r})^T \mathbf{s}
                       \f]
*		       and
*		       \f[
  \mathbf{R} = (r^2_4 - \mathbf{r}^T\mathbf{r}) \mathbf{I} +
               2 \mathbf{r} \mathbf{r}^T +
	       2 r_4 \mathbf{K}(\mathbf{r})
                       \f]
*		</ol>
* \param	nV			Number of vertices.
* \param	vW			Vertex weights (Walker's beta),
*					may be NULL which implies that
*					all the weights are 1.0.
* \param	vT			Vertices of the target set.
* \param	vS			Vertices of the source set.
* \param	nN			Number of normals.
* \param	nW			Normal weights (Walker's alpha).,
*					may be NULL which implies that
*					all the weights are 1.0
* \param	nT			Target normals, may be NULL.
* \param	nS			Source normals, may be NULL.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqDQ3D(int nV, double *vW,
					WlzDVertex3 *vT, WlzDVertex3 *vS,
					int nN, double *nW,
					WlzDVertex3 *nT, WlzDVertex3 *nS,
					WlzErrorNum *dstErr)
{
  int		id0;
  double	wt,
		sumVW,
		t0D,
		t11D,
		t12D,
		t13D,
		t14D,
		t21D,
		t22D,
		t23D,
		t24D,
		t31D,
		t32D,
		t33D,
		t34D,
		t44D;
  double	rV[4];
  double	**aA,
  		**t0A,
  		**t1A,
		**c1A,
		**c3A;
  WlzDVertex3	t0V,
  		t1V;
  AlgMatrix	aM, 			/* Walker's A */
  		t0M,			/* Working matrix */
  		t1M,			/* Working matrix */
  		t2M,			/* Working matrix */
		c1M,			/* Walker's C_1 */
		c3M;			/* Walker's C_3 */
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  aM.core = t0M.core = t1M.core = t2M.core = c1M.core = c3M.core = NULL;
  /* Allocate matricies required to compute c1M, c2M and c3M. */
  if(((aM.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL)||
     ((c1M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL)||
     ((c3M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL)||
     ((t0M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL)||
     ((t1M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL)||
     ((t2M.rect = AlgMatrixRectNew(4, 4, NULL)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  /* Compute c1M, c2M and c3M (see Walker's paper). */
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    t0A = t0M.rect->array;
    t1A = t1M.rect->array;
    c1A = c1M.rect->array;
    c3A = c3M.rect->array;
    sumVW = 0.0;
    for(id0 = 0; id0 < nV; ++id0)
    {
      /* Update sumVW */
      sumVW += (vW)? vW[id0]: 1.0;
      /* Update c1M: Compute \beta_i{Q(v_i^t)}^T W(v_i_s) */
      wt = (vW)? vW[id0]: 1.0;
      t0V = vT[id0];
      t1V = vS[id0];
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t13D = t0V.vtX * t1V.vtZ;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      t23D = t0V.vtY * t1V.vtZ;
      t31D = t0V.vtZ * t1V.vtX;
      t32D = t0V.vtZ * t1V.vtY;
      t33D = t0V.vtZ * t1V.vtZ;
      c1A[0][0] += wt * ( t11D - t22D - t33D);
      c1A[1][1] += wt * (-t11D + t22D - t33D);
      c1A[2][2] += wt * (-t11D - t22D + t33D);
      c1A[3][3] += wt * ( t11D + t22D + t33D);
      t0D = wt * ( t12D + t21D); c1A[0][1] += t0D; c1A[1][0] += t0D;
      t0D = wt * ( t13D + t31D); c1A[0][2] += t0D; c1A[2][0] += t0D;
      t0D = wt * (-t23D + t32D); c1A[0][3] += t0D; c1A[3][0] += t0D;
      t0D = wt * ( t23D + t32D); c1A[1][2] += t0D; c1A[2][1] += t0D;
      t0D = wt * ( t13D - t31D); c1A[1][3] += t0D; c1A[3][1] += t0D;
      t0D = wt * (-t12D + t21D); c1A[2][3] += t0D; c1A[3][2] += t0D;
      /* Update c3M: Compute \beta_i(W(v_i_s) - Q(v_i^t)) */
      t0D = wt * ( t1V.vtZ + t0V.vtZ); c3A[0][1] += t0D; c3A[1][0] -= t0D; 
      t0D = wt * (-t1V.vtY - t0V.vtY); c3A[0][2] += t0D; c3A[2][0] -= t0D; 
      t0D = wt * ( t1V.vtX - t0V.vtX); c3A[0][3] += t0D; c3A[3][0] -= t0D;
      t0D = wt * ( t1V.vtX + t0V.vtX); c3A[1][2] += t0D; c3A[2][1] -= t0D;
      t0D = wt * ( t1V.vtY - t0V.vtY); c3A[1][3] += t0D; c3A[3][1] -= t0D;
      t0D = wt * ( t1V.vtZ - t0V.vtZ); c3A[2][3] += t0D; c3A[3][2] -= t0D;
    }
    for(id0 = 0; id0 < nN; ++id0)
    {
      /* Update c1M: Compute \alpha_i{Q(n_i^t)}^T W(n_i_s) */
      wt = (nW)? nW[id0]: 1.0;
      t0V = nT[id0];
      t1V = nS[id0];
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t13D = t0V.vtX * t1V.vtZ;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      t23D = t0V.vtY * t1V.vtZ;
      t31D = t0V.vtZ * t1V.vtX;
      t32D = t0V.vtZ * t1V.vtY;
      t33D = t0V.vtZ * t1V.vtZ;
      c1A[0][0] += wt * ( t11D - t22D - t33D);
      c1A[1][1] += wt * (-t11D + t22D - t33D);
      c1A[2][2] += wt * (-t11D - t22D + t33D);
      c1A[3][3] += wt * ( t11D + t22D + t33D);
      t0D = wt * ( t12D + t21D); c1A[0][1] += t0D; c1A[1][0] += t0D;
      t0D = wt * ( t13D + t31D); c1A[0][2] += t0D; c1A[2][0] += t0D;
      t0D = wt * (-t23D + t32D); c1A[0][3] += t0D; c1A[3][0] += t0D;
      t0D = wt * ( t23D + t32D); c1A[1][2] += t0D; c1A[2][1] += t0D;
      t0D = wt * ( t13D - t31D); c1A[1][3] += t0D; c1A[3][1] += t0D;
      t0D = wt * (-t12D + t21D); c1A[2][3] += t0D; c1A[3][2] += t0D;
    }
    AlgMatrixScale(c1M, c1M, -2.0);
    AlgMatrixScale(c3M, c3M, 2.0);
    if(sumVW < DBL_EPSILON)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Now compute A from C1, C3 and \sum_i^n{\beta_i} */
  if(errNum == WLZ_ERR_NONE)
  {
    AlgMatrixTranspose(t0M, c3M);		/* T0 = C3^T */
    AlgMatrixMul(t1M, t0M, c3M);	 	/* T1 = C3^T C3 */
    t0D = 1.0 / (4.0 * sumVW);
    AlgMatrixScale(t1M, t1M, t0D);		
    AlgMatrixSub(aM, t1M, c1M);
    /* Find the eigenvector of A which has the greatest eigenvalue.
     * This is returned in the first column of aM. */
    errNum = WlzErrorFromAlg(AlgMatrixRSEigen(aM, rV, 1));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rV[0] = aA[0][0]; rV[1] = aA[1][0];
    rV[2] = aA[2][0]; rV[3] = aA[3][0];
    /* Compute the transform's rotation elements from the eigen vector. */
    t11D = rV[0] * rV[0];
    t22D = rV[1] * rV[1];
    t33D = rV[2] * rV[2];
    t44D = rV[3] * rV[3];
    t12D = 2.0 * rV[0] * rV[1];
    t13D = 2.0 * rV[0] * rV[2];
    t14D = 2.0 * rV[0] * rV[3];
    t23D = 2.0 * rV[1] * rV[2];
    t24D = 2.0 * rV[1] * rV[3];
    t34D = 2.0 * rV[2] * rV[3];
    aA[0][0] = t44D + t11D - t22D - t33D;
    aA[0][1] = t12D - t34D;
    aA[0][2] = t13D + t24D;
    aA[1][0] = t12D + t34D;
    aA[1][1] = t44D - t11D + t22D - t33D;
    aA[1][2] = t23D - t14D;
    aA[2][0] = t13D - t24D;
    aA[2][1] = t23D + t14D;
    aA[2][2] = t44D - t11D - t22D + t33D;
    /* Compute the translation elements from the eigen vector r, the sum of the
     * vertex weights \sum_i{\beta_i} and matrix C3. */
    /* Set t0M[0] to r (see Walker's paper). */
    t0A[0][0] = rV[0]; t0A[1][0] = rV[1]; t0A[2][0] = rV[2]; t0A[3][0] = rV[3];
    /* Compute s in t1M[0], but don't scale with -1/(2\sum_i{\beta_i} yet (see
     * Walker's paper). */
    t0M.rect->nC = 1; AlgMatrixMul(t1M, c3M, t0M); t0M.rect->nC = 4;
    /* Set t0M to be W^T(r) (see Walker's paper). */
    t0A[0][0] =  rV[3]; t0A[0][1] =  rV[2];
    t0A[0][2] =  rV[1]; t0A[0][3] = -rV[0];
    t0A[1][0] = -rV[2]; t0A[1][1] =  rV[3];
    t0A[1][2] =  rV[0]; t0A[1][3] = -rV[1];
    t0A[2][0] =  rV[1]; t0A[2][1] = -rV[0];
    t0A[2][2] =  rV[3]; t0A[2][3] = -rV[2];
    t0A[3][0] =  rV[0]; t0A[3][1] =  rV[1];
    t0A[3][2] =  rV[2]; t0A[3][3] =  rV[3];
    /* Compute the product W^T(r) s (see Walker's paper). */
    t1M.rect->nC = 1; AlgMatrixMul(t1M, t0M, t1M); t1M.rect->nC = 4;
    /* Extract the translation components and scale by -1/(2\sum_i{\beta_i}
     * (see Walker's paper). */
    t0D = -1.0 / (2.0 * sumVW);
    aA[0][3] = t0D * t1A[0][0];
    aA[1][3] = t0D * t1A[1][0];
    aA[2][3] = t0D * t1A[2][0];
    /* Set the transform's perspective and scale elements. */
    aA[3][0] = aA[3][1] = aA[3][2] = 0.0;
    aA[3][3] = 1.0;
    /* Create the affine transform. */
    trans = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_3D_AFFINE, aA, &errNum);
  }
  /* Free the matricies. */
  AlgMatrixFree(aM);
  AlgMatrixFree(c1M);
  AlgMatrixFree(c3M);
  AlgMatrixFree(t0M);
  AlgMatrixFree(t1M);
  AlgMatrixFree(t2M);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trans);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Solves the given linear system \f$\mathbf{Ax} = \mathbf{b}\f$
*		first using an SVD solver. If the design matrix is
*		ill-conditioned then the estimate is improved using an
*		itterative conjugate gradient solver.
* \param	aM			Square matrix A.
* \param	bV			Column matrix b, overwritten
*					by matrix x on return.
* \param	tol			Tolerance for singular values,
					1.0e-06 should be suitable as
					a default value.
*/
static WlzErrorNum WlzAffineTransformLSqLinSysSolve(AlgMatrix aM,
					double *bV, double tol)
{
  int		iCon;
  double	*bCV = NULL;
  AlgMatrix	aCM,
  		wCM;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	itr = 1000;

  wCM.core = NULL;
  if(((aCM.rect = AlgMatrixRectNew(aM.core->nR, aM.core->nC, NULL)) == NULL) ||
     ((wCM.rect = AlgMatrixRectNew(4, aM.core->nC, NULL)) == NULL) ||
     ((bCV = (double *)AlcMalloc(aM.core->nR * sizeof(double))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    AlgVectorCopy(bCV, bV, aM.core->nR);
    AlgMatrixCopy(aCM, aM);
    errNum = WlzErrorFromAlg(AlgMatrixSVSolve(aM, bV, tol, &iCon));
  }
  if((errNum == WLZ_ERR_NONE) && (iCon > 0))
  {
    errNum = WlzErrorFromAlg(AlgMatrixCGSolve(aCM, bV, bCV, wCM,
    					      NULL, NULL, tol, itr,
					      NULL, NULL));
  }
  AlcFree(bCV);
  AlgMatrixFree(aCM);
  AlgMatrixFree(wCM);
  return(errNum);
}
