#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzAffineTransformLSq.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions for computing Woolz affine transforms that
*		give the best fit, in a least squares sense, when
*		used to transform one set of vertices to another.
* \ingroup      WlzTransform
* \todo         -
* \bug          None known.
* \note
* Maintenance log with most recent changes at top of list.
*/
#include <stdio.h>
#include <float.h>
#include <Wlz.h>

static WlzAffineTransform 	*WlzAffineTransformLSqWgt3D(
				  int nVtx,
				  double *wgt,
				  WlzDVertex3 *vtxVec0,
				  WlzDVertex3 *vtxVec1,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzAffineTransformLSqWgt2D(
				  int nVtx,
				  double *wgt,
				  WlzDVertex2 *vtxVec0,
				  WlzDVertex2 *vtxVec1,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzAffineTransformLSqTrans3D(
				  WlzDVertex3 *pos0,
				  WlzDVertex3 *pos1,
				  int nVtx,
				  WlzErrorNum *dstErr);
static WlzAffineTransform	*WlzAffineTransformLSqReg3D(
				  WlzDVertex3 *vtxVec0,
				  WlzDVertex3 *vtxVec1,
				  int nVtx,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzAffineTransformLSqGen2D(
				  WlzDVertex2 *vtxVec0,
				  WlzDVertex2 *vtxVec1,
				  int nVtx,
				  WlzErrorNum *dstErr);
static WlzAffineTransform	*WlzAffineTransformLSqReg2D(
				  WlzDVertex2 *vtxVec0,
				  WlzDVertex2 *vtxVec1,
				  int nVtx,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
static WlzAffineTransform	*WlzAffineTransformLSqTrans2D(
				  WlzDVertex2 *vtxVec0,
				  WlzDVertex2 *vtxVec1,
				  int nVtx,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzAffineTransformLSqWgtGen2D(
				  WlzDVertex2 *vtxVec0,
				  WlzDVertex2 *vtxVec1,
				  double *wgt,
				  int nVtx,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzAffineTransformLSqDQ3D(
				  int nV,
				  double *vW,
				  WlzDVertex3 *v0,
				  WlzDVertex3 *v1,
				  int nN,
				  double *nW,
				  WlzDVertex3 *n0,
				  WlzDVertex3 *n1,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzAffineTransformLSqDQ2D(
				  int nV,
				  double *vW,
				  WlzDVertex2 *v1,
				  WlzDVertex2 *v0,
				  int nN,
				  double *nW,
				  WlzDVertex2 *n1,
				  WlzDVertex2 *n0,
				  WlzErrorNum *dstErr);

/*!
* \ingroup	WlzTransform
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of vertices and normals onto the second
*		set. The vertex and normal weighting factors must
*		be in the range [0-1].
* \param	vtxType			Type of vertices.
* \param	nV			Number of vertices.
* \param	vW			Vertex weights, may be NULL which
*					implies that all the weights are 1.0.
* \param	v0			Vertices of the first set.
* \param	v1			Vertices of the second set.
* \param	nW			Normal weights, may be NULL which
*					implies that all the weights are 1.0
* \param	nN			Number of normals, may be zero.
* \param	n0			Normals of the first set, may be NULL.
* \param	n1			Normals of the second set, may be NULL.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSq2(WlzVertexType vtxType,
					   int nV, double *vW,
					   WlzVertexP v0, WlzVertexP v1,
					   int nN, double *nW,
					   WlzVertexP n0, WlzVertexP n1,
					   WlzTransformType trType,
					   WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nV <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((v0.v == NULL) || (v1.v == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((vtxType != WLZ_VERTEX_D2) && (vtxType != WLZ_VERTEX_D3))
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  else
  {
    switch(trType)
    {
      case WLZ_TRANSFORM_2D_REG:
	tr = WlzAffineTransformLSqDQ2D(nV, vW, v0.d2, v1.d2,
				       nN, nW, n0.d2, n1.d2, &errNum);
        break;
      case WLZ_TRANSFORM_3D_REG: 
	tr = WlzAffineTransformLSqDQ3D(nV, vW, v0.d3, v1.d3,
				       nN, nW, n0.d3, n1.d3, &errNum);
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
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of vertices onto the second.
* \param	vtxType			Type of vertices.
* \param	nVtx0			Number of vertices in first vector.
* \param	vtxVec0			First vector of vertices.
* \param	nVtx1			Number of vertices in second
*					vector (MUST be same as first).
* \param	vtxVec1			Second vector of vertices.
* \param	trType			 Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSq(WlzVertexType vtxType,
				          int nVtx0, WlzVertexP vtxVec0,
					  int nVtx1, WlzVertexP vtxVec1,
					  WlzTransformType trType,
					  WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nVtx0 != nVtx1) || (nVtx0 <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vtxVec0.v == NULL) || (vtxVec1.v == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(vtxType)
    {
      case WLZ_VERTEX_D2:
	tr = WlzAffineTransformLSq2D(nVtx0, vtxVec0.d2,
				     nVtx1, vtxVec1.d2, 
				     trType, &errNum);
	break;
      case WLZ_VERTEX_D3:
	tr = WlzAffineTransformLSq3D(nVtx0, vtxVec0.d3, 
				     nVtx1, vtxVec1.d3,
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
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of vertices onto the second, while using
*		the given vertex weights.
* \param	vtxType			Type of vertices.
* \param	nVtx			Number of vertex pairs.
* \param	wgt			Weights for the matched pairs
*					of vertices.
* \param	vtxVec0			First vector of vertices.
* \param	vtxVec1			Second vector of vertices.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSqWgt(WlzVertexType vtxType,
				          int nVtx, double *wgt,
					  WlzVertexP vtxVec0,
					  WlzVertexP vtxVec1,
					  WlzTransformType trType,
					  WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(nVtx <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vtxVec0.v == NULL) || (vtxVec1.v == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(vtxType)
    {
      case WLZ_VERTEX_D2:
	tr = WlzAffineTransformLSqWgt2D(nVtx, wgt, vtxVec0.d2, vtxVec1.d2, 
					trType, &errNum);
	break;
      case WLZ_VERTEX_D3:
	tr = WlzAffineTransformLSqWgt3D(nVtx, wgt, vtxVec0.d3, vtxVec1.d3, 
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
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of 2D double vertices onto the second,
*		while using the given vertex weights.
* \param	nVtx			Number of vertex pairs.
* \param	wgt			Weights for the matched pairs
*					of vertices.
* \param	vtxVec0			First vector of vertices.
* \param	vtxVec1			Second vector of vertices.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqWgt2D(int nVtx, double *wgt,
					WlzDVertex2 *vtxVec0,
					WlzDVertex2 *vtxVec1,
					WlzTransformType trType,
					WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum errNum = WLZ_ERR_UNIMPLEMENTED;

  switch(trType)
  {
    case WLZ_TRANSFORM_2D_AFFINE:
      tr = WlzAffineTransformLSqWgtGen2D(vtxVec0, vtxVec1, wgt, nVtx, &errNum);
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
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
*		first set of 3D double vertices onto the second,
*		while using the given vertex weights.
* \todo		UNIMPLEMENTED
* \param	nVtx			Number of vertex pairs.
* \param	wgt			Weights for the matched pairs
*					of vertices.
* \param	vtxVec0			First vector of vertices.
* \param	vtxVec1			Second vector of vertices.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqWgt3D(int nVtx, double *wgt,
					WlzDVertex3 *vtxVec0,
					WlzDVertex3 *vtxVec1,
					WlzTransformType trType,
					WlzErrorNum *dstErr)
{
  WlzAffineTransform *tr = NULL;
  WlzErrorNum errNum = WLZ_ERR_UNIMPLEMENTED;

  /* TODO UNIMPLEMENTED. */
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
*		first set of 3D vertices onto the second.
* \param	nVtx0			Number of vertices in first
*					vector.
* \param	vtxVec0			First vector of vertices.
* \param	nVtx1			Number of vertices in second
*					vector (MUST be same as first).
* \param	vtxVec1			Second vector of vertices.
* \param	trType			 Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSq3D(int nVtx0, WlzDVertex3 *vtxVec0,
					    int nVtx1, WlzDVertex3 *vtxVec1,
					    WlzTransformType trType,
					    WlzErrorNum *dstErr)
{
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nVtx0 != nVtx1) || (nVtx0 <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vtxVec0 == NULL) || (vtxVec1 == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(trType)
    {
      case WLZ_TRANSFORM_3D_TRANS:
	trans = WlzAffineTransformLSqTrans3D(vtxVec0, vtxVec1, nVtx0, &errNum);
	break;
      case WLZ_TRANSFORM_3D_REG:
	if(nVtx0 == 1)
	{
	  trans = WlzAffineTransformLSqTrans3D(vtxVec0, vtxVec1, nVtx0,
	  				       &errNum);
	}
	else
	{
	  trans = WlzAffineTransformLSqReg3D(vtxVec0, vtxVec1, nVtx0, &errNum);
	}
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
  return(trans);
}

/*!
* \ingroup	WlzTransform
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of 2D vertices onto the second.
* \param	nVtx0			Number of vertices in first
*					vector.
* \param	vtxVec0			First vector of vertices.
* \param	nVtx1			Number of vertices in second
*					vector (MUST be same as first).
* \param	vtxVec1			Second vector of vertices.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzAffineTransform *WlzAffineTransformLSq2D(int nVtx0, WlzDVertex2 *vtxVec0,
					    int nVtx1, WlzDVertex2 *vtxVec1,
					    WlzTransformType trType,
					    WlzErrorNum *dstErr)
{
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nVtx0 != nVtx1) || (nVtx0 <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA; 
  }
  else if((vtxVec0 == NULL) || (vtxVec1 == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(trType)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
      case WLZ_TRANSFORM_2D_REG:
      case WLZ_TRANSFORM_2D_NOSHEAR:
      case WLZ_TRANSFORM_2D_TRANS:
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((trType == WLZ_TRANSFORM_2D_TRANS) || (nVtx0 == 1))
    {
      trans = WlzAffineTransformLSqTrans2D(vtxVec0, vtxVec1, nVtx0, &errNum);
    }
    else if((trType == WLZ_TRANSFORM_2D_REG) ||
	    (trType == WLZ_TRANSFORM_2D_NOSHEAR) || (nVtx0 < 3) )
    {
      trans = WlzAffineTransformLSqReg2D(vtxVec0, vtxVec1, nVtx0, trType,
      					 &errNum);
    }
    else
    {
      trans = WlzAffineTransformLSqGen2D(vtxVec0, vtxVec1, nVtx0, &errNum);
    }
  }
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
* \brief	Computes the Woolz 3D registration transform which
*		gives the best (least squares) fit when used to
*		transform the first set of vertices onto the second.
*		The transform is constrained to rotation and
*		translation only.
*		The algorithm is based on Arun K.S., Huang T.T.
*		and Blostein S.D. "Least-Squares Fitting of Two 3-D
*		Point Sets" PAMI 9(5), 698-700, 1987.
* \param	pos0			First vector of vertices.
* \param	pos1			Second vector of vertices.
* \param	nVtx			Number of vertices in vectors.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqReg3D(WlzDVertex3 *pos0,
		    		WlzDVertex3 *pos1, int nVtx,
				WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
  		idK,
		idR;
  double	tD0,
  		meanSqD;
  WlzDVertex3	cen0,
  		cen1,
		rel0,
		rel1;
  double	*wMx = NULL;
  double	**hMx = NULL,
  		**vMx = NULL,
		**trMx = NULL;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0E-06;

  if(((wMx = (double *)AlcCalloc(sizeof(double), 3)) == NULL) ||
     (AlcDouble2Calloc(&hMx, 3, 3) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&vMx, 3, 3) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&trMx, 4, 4) !=  ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute centroids (cen0 and cen1) and a mean of squares of distance
     * between the vertices. */
    meanSqD = 0.0;
    cen0 = *pos0;
    cen1 = *pos1;
    for(idN = 1; idN < nVtx; ++idN)
    {
      WLZ_VTX_3_ADD(cen0, cen0, *(pos0 + idN));
      WLZ_VTX_3_ADD(cen1, cen1, *(pos1 + idN));
      WLZ_VTX_3_SUB(rel0, *(pos0 + idN), *(pos1 + idN));
      meanSqD += WLZ_VTX_3_DOT(rel0, rel0);
    }
    tD0 = 1.0 / nVtx;
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
      /* Compute the 3x3 matrix hMx, which is the sum of tensor products of
       * the vertices relative to their centroids. */
      for(idN = 0; idN < nVtx; ++idN)
      {
	WLZ_VTX_3_SUB(rel0, *(pos0 + idN), cen0);
	WLZ_VTX_3_SUB(rel1, *(pos1 + idN), cen1);
	hMx[0][0] += rel0.vtX * rel1.vtX;
	hMx[0][1] += rel0.vtX * rel1.vtY;
	hMx[0][2] += rel0.vtX * rel1.vtZ;
	hMx[1][0] += rel0.vtY * rel1.vtX;
	hMx[1][1] += rel0.vtY * rel1.vtY;
	hMx[1][2] += rel0.vtY * rel1.vtZ;
	hMx[2][0] += rel0.vtZ * rel1.vtX;
	hMx[2][1] += rel0.vtZ * rel1.vtY;
	hMx[2][2] += rel0.vtZ * rel1.vtZ;
      }
      /* Compute the SVD of the 3x3 matrix, hMx = hMx.mMx.vMx. */
      errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(hMx, 3, 3, wMx, vMx));
      if(errNum == WLZ_ERR_NONE)
      {
	/* Compute 3x3 rotation matrix trMx = vMx.hMx', where hMx' is the
	 * transpose of hMx. */
	for(idR = 0; idR < 3; ++idR)
	{
	  for(idK = 0; idK < 3; ++idK)
	  {
	    trMx[idR][idK] = vMx[idR][0] * hMx[idK][0] + 
	      vMx[idR][1] * hMx[idK][1] + 
	      vMx[idR][2] * hMx[idK][2];
	  }
	}
	/* Test for degeneracy using the determinant of the rotation matrix. */
	tD0 = (trMx[0][0] * trMx[1][1] * trMx[2][2]) -
	  (trMx[0][0] * trMx[1][2] * trMx[2][1]) +
	  (trMx[0][1] * trMx[1][2] * trMx[2][0]) -
	  (trMx[0][1] * trMx[1][0] * trMx[2][2]) +
	  (trMx[0][2] * trMx[1][0] * trMx[2][1]) -
	  (trMx[0][2] * trMx[1][1] * trMx[2][0]);
	if(tD0 < 0.0)
	{
	  /* Are source vertices (rel0) coplanar? They are iff one of the 3
	   * singular values of hMx in wMx is zero. If the source vertices
	   * are not coplanar, the the solution of the SVD is correct. */
	  tI0 = ((fabs(*(wMx + 2)) >= DBL_EPSILON) << 2) |
	    ((fabs(*(wMx + 1)) >= DBL_EPSILON) << 1) |
	    (fabs(*(wMx + 2)) >= DBL_EPSILON);
	  switch(tI0)
	  {
	    case 0:
	      /* Source vertices are not coplanar or colinear. The SVD gives the
	       * correct solution. */
	      break;
	    case 1:
	    case 2:
	    case 4:
	      /* Source vertices are coplanar, but not colinear. There is a
	       * unique reflection as well as a unique rotation. The SVD
	       * may give either BUT in this case it has found the reflection
	       * so need to recompute for the rotation.
	       * If the singular values of wMx are w0 > w1 > w2 = 0, then
	       *
	       *              t        t        t
	       *   hMx = w u v  + w u v  + 0.u v 
	       *          0 0 0    1 1 1      2 2
	       *
	       * where ui and vi are the columns of the matricies hMx and vMx
	       * respectively (hMx is used for both H and U).
	       * So to get the rotation rather than the reflection we recalculate
	       * the transform matrix with v2 = -v2.
	       */
	      for(idR = 0; idR < 3; ++idR)
	      {
		for(idK = 0; idK < 3; ++idK)
		{
		  trMx[idR][idK] = vMx[idR][0] * hMx[idK][0] + 
		    vMx[idR][1] * hMx[idK][1] -
		    vMx[idR][2] * hMx[idK][2];
		}
	      }
	      break;
	    case 3:
	    case 5:
	    case 6:
	      /* Source vertices are colinear and there exists an infinity of
	       * solutions! */
	      errNum = WLZ_ERR_ALG_SINGULAR;
	      break;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* Fill in other matrix elements. */
	trMx[3][0] = trMx[3][1] = trMx[3][2] = 0.0;
	trMx[3][3] = 1.0;
	/* Compute translation by applying the rotation transform to the first
	 * centroid and subtracting this from the second centroid. */
	trMx[0][3] = cen1.vtX -
	  ((trMx[0][0] * cen0.vtX) + (trMx[0][1] * cen0.vtY) +
	   (trMx[0][2] * cen0.vtZ));
	trMx[1][3] = cen1.vtY -
	  ((trMx[1][0] * cen0.vtX) + (trMx[1][1] * cen0.vtY) +
	   (trMx[1][2] * cen0.vtZ));
	trMx[2][3] = cen1.vtZ -
	  ((trMx[2][0] * cen0.vtX) + (trMx[2][1] * cen0.vtY) +
	   (trMx[2][2] * cen0.vtZ));
	/* Build affine transform. */
	tr = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_3D_AFFINE, trMx,
	    &errNum);
      }
      /* Clear up on error. */
      if(wMx)
      {
	(void )AlcFree(wMx);
      }
      if(hMx)
      {
	(void )AlcDouble2Free(hMx);
      }
      if(vMx)
      {
	(void )AlcDouble2Free(vMx);
      }
      if(trMx)
      {
	(void )AlcDouble2Free(trMx);
      }
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
* \return				Computed affine transform, may
*					be NULL on error.
* \brief	Computes the Woolz 3D registration transform which
*		gives the best (least squares) fit when used to
*		transform the first set of vertices onto the second.
*		The transform is constrained to translation only.
* \param	pos0			First vector of vertices.
* \param	pos1			Second vector of vertices.
* \param	nVtx			Number of vertices in vectors.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqTrans3D(WlzDVertex3 *pos0,
		    		WlzDVertex3 *pos1, int nVtx,
				WlzErrorNum *dstErr)
{
  int		idN;
  double	tD0;
  WlzDVertex3	cen0,
  		cen1;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Compute centroids, cen0 and cen1. */
  cen0 = *pos0;
  cen1 = *pos1;
  for(idN = 1; idN < nVtx; ++idN)
  {
    WLZ_VTX_3_ADD(cen0, cen0, *(pos0 + idN));
    WLZ_VTX_3_ADD(cen1, cen1, *(pos1 + idN));
  }
  tD0 = 1.0 / nVtx;
  WLZ_VTX_3_SCALE(cen0, cen0, tD0);
  WLZ_VTX_3_SCALE(cen1, cen1, tD0);
  /* Build the affine transform which brings the first centroid to the
   * position of the second. */
  tr = WlzAffineTransformFromTranslation(WLZ_TRANSFORM_3D_AFFINE,
					 cen1.vtX - cen0.vtX,
					 cen1.vtY - cen0.vtY,
					 cen1.vtZ - cen0.vtZ,
					 &errNum);
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
* \brief	Computes the Woolz general 2D affine transform which
*		gives the best (least squares) fit when used to
*		transform the first set of vertices onto the second.
* \param	vtxVec0			First vector of vertices.
* \param	vtxVec1			Second vector of vertices.
* \param	nVtx			Number of vertices in vectors.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqGen2D(WlzDVertex2 *vtxVec0,
				WlzDVertex2 *vtxVec1, int nVtx,
				WlzErrorNum *dstErr)
{
  int		idx;
  double 	*b = NULL;
  double	**a = NULL,
		**trMat = NULL;
  double	A[12];
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Initialise the array */
  for(idx = 0; idx < 12; ++idx)
  {
    A[idx] = 0;
  }
  /* Accumulate values */
  for(idx = 0; idx < nVtx; ++idx, ++vtxVec0, ++vtxVec1)
  {
    A[0] += 1;
    A[1] += vtxVec0->vtX;
    A[2] += vtxVec0->vtY;
    A[3] += vtxVec0->vtX * vtxVec0->vtX;
    A[4] += vtxVec0->vtX * vtxVec0->vtY;
    A[5] += vtxVec0->vtY * vtxVec0->vtY;
    A[6] += vtxVec1->vtX;
    A[7] += vtxVec0->vtX * vtxVec1->vtX;
    A[8] += vtxVec0->vtY * vtxVec1->vtX;
    A[9] += vtxVec1->vtY;
    A[10] += vtxVec0->vtX * vtxVec1->vtY;
    A[11] += vtxVec0->vtY * vtxVec1->vtY;
  }
  /* Allocate workspace */
  if((AlcDouble2Malloc(&a, 3, 3) != ALC_ER_NONE) ||
      (AlcDouble1Malloc(&b, 3) != ALC_ER_NONE) ||
      (AlcDouble2Malloc(&trMat, 4, 4) != ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Determine the least square transformation matrix values */
    /* x parameters first */
    a[0][0] = A[0];  a[0][1] = A[1];  a[0][2] = A[2];
    a[1][0] = A[1];  a[1][1] = A[3];  a[1][2] = A[4];
    a[2][0] = A[2];  a[2][1] = A[4];  a[2][2] = A[5];
    b[0] = A[6];     b[1] = A[7];     b[2] = A[8];
    errNum = WlzErrorFromAlg(AlgMatrixLUSolve(a, 3, b, 1));
    trMat[0][0] = b[1]; trMat[0][1] = b[2]; trMat[0][2] = b[0];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* now y parameters */
    a[0][0] = A[0];  a[0][1] = A[1];  a[0][2] = A[2];
    a[1][0] = A[1];  a[1][1] = A[3];  a[1][2] = A[4];
    a[2][0] = A[2];  a[2][1] = A[4];  a[2][2] = A[5];
    b[0] = A[9];     b[1] = A[10];    b[2] = A[11];
    errNum = WlzErrorFromAlg(AlgMatrixLUSolve(a, 3, b, 1));
    trMat[1][0] = b[1]; trMat[1][1] = b[2]; trMat[1][2] = b[0];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make the transform */
    trMat[2][0] = 0; trMat[2][1] = 0; trMat[2][2] = 1;
    trans = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, trMat,
					 &errNum);
  }
  if(a)
  {
    (void )AlcDouble2Free(a);
  }
  if(b)
  {
    (void )AlcFree((void *)b);
  }
  if(trMat)
  {
    (void )AlcDouble2Free(trMat);
  }
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
* \brief	Computes the Woolz 2D registration transform which
*		gives the best (least squares) fit when used to
*		transform the first set of vertices onto the second.
*		The transform is constrained to rotation and
*		translation only.
* \param	vtxVec0			First vector of vertices.
* \param	vtxVec1			Second vector of vertices.
* \param	nVtx			Number of vertices in vectors.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqReg2D(WlzDVertex2 *vtxVec0,
		    		WlzDVertex2 *vtxVec1, int nVtx,
				WlzTransformType trType,
				WlzErrorNum *dstErr)
{
  int		idx;
  double	**a,
		 **trMat;
  double	*b;
  WlzAffineTransform *trans;
  double	A[12];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Initialise the array */
  for(idx=0; idx < 12; idx++)
  {
    A[idx] = 0;
  }
  /* accumulate values */
  for(idx = 0; idx < nVtx; ++idx, ++vtxVec0, ++vtxVec1)
  {
    A[0] += 1;
    A[1] += vtxVec0->vtX;
    A[2] += vtxVec0->vtY;
    A[3] += vtxVec0->vtX * vtxVec0->vtX;
    /* A[4] += vtxVec0->vtX * vtxVec0->vtY;*/
    A[5] += vtxVec0->vtY * vtxVec0->vtY;
    A[6] += vtxVec1->vtX;
    A[7] += vtxVec0->vtX * vtxVec1->vtX;
    A[8] += vtxVec0->vtY * vtxVec1->vtX;
    A[9] += vtxVec1->vtY;
    A[10] += vtxVec0->vtX * vtxVec1->vtY;
    A[11] += vtxVec0->vtY * vtxVec1->vtY;
  }
  /* Allocate workspace */
  if((AlcDouble2Malloc(&a, 4, 4) != ALC_ER_NONE) ||
      (AlcDouble1Malloc(&b, 4) != ALC_ER_NONE) ||
      (AlcDouble2Malloc(&trMat, 4, 4) != ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Determine the least square transformation matrix values */
    a[0][0]=A[0];a[0][1]=A[1];       a[0][2]= A[2];       a[0][3]= 0;
    a[1][0]=A[1];a[1][1]=A[3] + A[5];a[1][2]= 0;          a[1][3]= A[2];
    a[2][0]=A[2];a[2][1]=0;          a[2][2]= A[3] + A[5];a[2][3]= -A[1];
    a[3][0]=0;   a[3][1]=A[2];       a[3][2]= -A[1];      a[3][3]= A[0];
    b[0]=A[6];   b[1]=A[7] + A[11];  b[2]=A[8] - A[10];   b[3]= A[9];
    errNum = WlzErrorFromAlg(AlgMatrixLUSolve(a, 4, b, 1));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check for scale constraint - this is a kludge */
    if(trType == WLZ_TRANSFORM_2D_REG)
    {
      double	s;
      
      s = sqrt(b[1]*b[1] + b[2]*b[2]);
      b[1] /= s;
      b[2] /= s;
      b[0] = (A[6] - A[1]*b[1] - A[2]*b[2]) / A[0];
      b[3] = (A[9] + A[1]*b[2] - A[2]*b[1]) / A[0];
    }
    /* Make the transformation */
    trMat[0][0] = b[1];  trMat[0][1] = b[2]; trMat[0][2] = b[0];
    trMat[1][0] = -b[2]; trMat[1][1] = b[1]; trMat[1][2] = b[3];
    trMat[2][0] = 0; trMat[2][1] = 0; trMat[2][2] = 1;
    trans = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, trMat,
					 &errNum);
  }
  if(a)
  {
    (void )AlcDouble2Free(a);
  }
  if(b)
  {
    (void )AlcFree((void *)b);
  }
  if(trMat)
  {
    (void )AlcDouble2Free(trMat);
  }
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
* \brief	Computes the Woolz general 2D affine transform which
*		gives the best (least squares) fit when used to
*		transform the first set of vertices onto the second,
*		while using the given vertex weights.
* \param	vtxVec0			First vector of vertices.
* \param	vtxVec1			Second vector of vertices.
* \param	nVtx			Number of vertices in vectors.
* \param	wgt			Weights for the matched pairs
*					of vertices.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqWgtGen2D(WlzDVertex2 *vtxVec0,
				WlzDVertex2 *vtxVec1, double *wgt, int nVtx,
				WlzErrorNum *dstErr)
{
  int		idx;
  double 	*b = NULL;
  double	**a = NULL,
		**trMat = NULL;
  double	A[12];
  WlzDVertex2	wV0,
  		wV1;
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Initialise the array */
  for(idx = 0; idx < 12; ++idx)
  {
    A[idx] = 0;
  }
  /* Accumulate values */
  for(idx = 0; idx < nVtx; ++idx)
  {
    WLZ_VTX_2_SCALE(wV0, *vtxVec0, wgt[idx]);
    WLZ_VTX_2_SCALE(wV1, *vtxVec1, wgt[idx]);
    A[0] += 1;
    A[1] += wV0.vtX;
    A[2] += wV0.vtY;
    A[3] += wV0.vtX * wV0.vtX;
    A[4] += wV0.vtX * wV0.vtY;
    A[5] += wV0.vtY * wV0.vtY;
    A[6] += wV1.vtX;
    A[7] += wV0.vtX * wV1.vtX;
    A[8] += wV0.vtY * wV1.vtX;
    A[9] += wV1.vtY;
    A[10] += wV0.vtX * wV1.vtY;
    A[11] += wV0.vtY * wV1.vtY;
    ++vtxVec0;
    ++vtxVec1;
  }
  /* Allocate workspace */
  if((AlcDouble2Malloc(&a, 3, 3) != ALC_ER_NONE) ||
      (AlcDouble1Malloc(&b, 3) != ALC_ER_NONE) ||
      (AlcDouble2Malloc(&trMat, 4, 4) != ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Determine the least square transformation matrix values */
    /* x parameters first */
    a[0][0] = A[0];  a[0][1] = A[1];  a[0][2] = A[2];
    a[1][0] = A[1];  a[1][1] = A[3];  a[1][2] = A[4];
    a[2][0] = A[2];  a[2][1] = A[4];  a[2][2] = A[5];
    b[0] = A[6];     b[1] = A[7];     b[2] = A[8];
    errNum = WlzErrorFromAlg(AlgMatrixLUSolve(a, 3, b, 1));
    trMat[0][0] = b[1]; trMat[0][1] = b[2]; trMat[0][2] = b[0];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* now y parameters */
    a[0][0] = A[0];  a[0][1] = A[1];  a[0][2] = A[2];
    a[1][0] = A[1];  a[1][1] = A[3];  a[1][2] = A[4];
    a[2][0] = A[2];  a[2][1] = A[4];  a[2][2] = A[5];
    b[0] = A[9];     b[1] = A[10];    b[2] = A[11];
    errNum = WlzErrorFromAlg(AlgMatrixLUSolve(a, 3, b, 1));
    trMat[1][0] = b[1]; trMat[1][1] = b[2]; trMat[1][2] = b[0];
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make the transform */
    trMat[2][0] = 0; trMat[2][1] = 0; trMat[2][2] = 1;
    trans = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, trMat,
					 &errNum);
  }
  if(a)
  {
    (void )AlcDouble2Free(a);
  }
  if(b)
  {
    (void )AlcFree((void *)b);
  }
  if(trMat)
  {
    (void )AlcDouble2Free(trMat);
  }
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
* \brief	Computes the Woolz 2D translation transform which
*		gives the best (least squares) fit when used to
*		transform the first set of vertices onto the second.
*		The transform is constrained to rotation and
*		translation only.
* \param	vtxVec0			First vector of vertices.
* \param	vtxVec1			Second vector of vertices.
* \param	nVtx			Number of vertices in vectors.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqTrans2D(WlzDVertex2 *vtxVec0,
		    		WlzDVertex2 *vtxVec1, int nVtx,
				WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex2	sum;
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  sum.vtX = 0.0;
  sum.vtY = 0.0;
  for(idx = 0; idx < nVtx; ++idx)
  {
    sum.vtX += vtxVec1->vtX - vtxVec0->vtX;
    sum.vtY += vtxVec1->vtY - vtxVec0->vtY;
    ++vtxVec0;
    ++vtxVec1;
  }
  trans = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE, 
					sum.vtX / nVtx, sum.vtY / nVtx, 0.0,
					1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
					&errNum);
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
*		first set of vertices and normals onto the second
*		set. The vertex and normal weighting factors must
*		be in the range [0-1].
*		This function is based on the dual quaternion
*		algorithm: M.W. Walker  and Shao L. Estimating 3-D
*		Location Parameters Using Dual Number Quaternions,
*		CVGIP 54(3), 1991.
*		Equation 47 from this paper can be simplified to
*		\f[
		 \mathbf{A} = {\frac{1}{4\sum{i=1}{n}{\beta_i}}}
			      {\mathbf{C}_3^T\mathbf{C}_3} -
			      \mathbf{C}_1
		\f]
*    		because \f$\mathbf{C}_1\f$ is real and symetric,
*		\f$\mathbf{C}_2\f$ is scalar and \f$\mathbf{C}_3\f$
*		is anti-symetric.
* \param	vtxType			Type of vertices.
* \param	nV			Number of vertices.
* \param	vW			Vertex weights (Walker's beta),
*					may be NULL which implies that
*					all the weights are 1.0.
* \param	v1			Vertices of the target set.
* \param	v0			Vertices of the source set.
* \param	nN			Numbr of normals.
* \param	nW			Normal weights (Walker's alpha).,
*					may be NULL which implies that
*					all the weights are 1.0
* \param	n1			Normals of the target set, may
*					be NULL.
* \param	n0			Normals of the source set, may
*					be NULL.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqDQ3D(int nV, double *vW,
					WlzDVertex3 *v1, WlzDVertex3 *v0,
					int nN, double *nW,
					WlzDVertex3 *n1, WlzDVertex3 *n0,
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
  double	rM[4];
  double	**aM = NULL,		/* Walker's A */
  		**t0M = NULL,		/* Working matrix */
  		**t1M = NULL,		/* Working matrix */
  		**t2M = NULL,		/* Working matrix */
		**c1M = NULL,		/* Walker's C_1 */
		**c3M = NULL;		/* Walker's C_3 */
  WlzDVertex3	t0V,
  		t1V;
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  /* Allocate matricies required to compute c1M, c2M and c3M. */
  if((AlcDouble2Calloc(&aM, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&c1M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&c3M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&t0M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&t1M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&t2M, 4, 4) != ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  /* Compute c1M, c2M and c3M (see Walker's paper). */
  if(errNum == WLZ_ERR_NONE)
  {
    sumVW = 0.0;
    for(id0 = 0; id0 < nV; ++id0)
    {
      /* Update sumVW */
      sumVW += (vW)? vW[id0]: 1.0;
      /* Update c1M: Compute \beta_i{Q(v_i^t)}^TW(v_i_s) */
      wt = (vW)? vW[id0]: 1.0;
      t0V = v0[id0];
      t1V = v1[id0];
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t13D = t0V.vtX * t1V.vtZ;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      t23D = t0V.vtY * t1V.vtZ;
      t31D = t0V.vtZ * t1V.vtX;
      t32D = t0V.vtZ * t1V.vtY;
      t33D = t0V.vtZ * t1V.vtZ;
      c1M[0][0] += wt * ( t11D - t22D - t33D);
      c1M[1][1] += wt * (-t11D + t22D - t33D);
      c1M[2][2] += wt * (-t11D - t22D + t33D);
      c1M[3][3] += wt * ( t11D + t22D + t33D);
      t0D = wt * ( t12D + t21D); c1M[0][1] += t0D; c1M[1][0] += t0D;
      t0D = wt * ( t13D + t31D); c1M[0][2] += t0D; c1M[2][0] += t0D;
      t0D = wt * (-t23D + t32D); c1M[0][3] += t0D; c1M[3][0] += t0D;
      t0D = wt * ( t23D + t32D); c1M[1][2] += t0D; c1M[2][1] += t0D;
      t0D = wt * ( t13D - t31D); c1M[1][3] += t0D; c1M[3][1] += t0D;
      t0D = wt * (-t12D + t21D); c1M[2][3] += t0D; c1M[3][2] += t0D;
      /* Update c3M: Compute \beta_i(W(v_i_s) - Q(v_i^t)) */
      t0D = wt * ( t1V.vtZ + t0V.vtZ); c3M[0][1] += t0D; c3M[1][0] -= t0D; 
      t0D = wt * (-t1V.vtY - t0V.vtY); c3M[0][2] += t0D; c3M[2][0] -= t0D; 
      t0D = wt * ( t1V.vtX - t0V.vtX); c3M[0][3] += t0D; c3M[3][0] -= t0D;
      t0D = wt * ( t1V.vtX + t0V.vtX); c3M[1][2] += t0D; c3M[2][1] -= t0D;
      t0D = wt * ( t1V.vtY - t0V.vtY); c3M[1][3] += t0D; c3M[3][1] -= t0D;
      t0D = wt * ( t1V.vtZ - t0V.vtZ); c3M[2][3] += t0D; c3M[3][2] -= t0D;
    }
    for(id0 = 0; id0 < nN; ++id0)
    {
      /* Update c1M: Compute \alpha_i{Q(n_i^t)}^TW(n_i_s) */
      wt = (nW)? nW[id0]: 1.0;
      t0V = n0[id0];
      t1V = n1[id0];
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t13D = t0V.vtX * t1V.vtZ;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      t23D = t0V.vtY * t1V.vtZ;
      t31D = t0V.vtZ * t1V.vtX;
      t32D = t0V.vtZ * t1V.vtY;
      t33D = t0V.vtZ * t1V.vtZ;
      c1M[0][0] += wt * ( t11D - t22D - t33D);
      c1M[1][1] += wt * (-t11D + t22D - t33D);
      c1M[2][2] += wt * (-t11D - t22D + t33D);
      c1M[3][3] += wt * ( t11D + t22D + t33D);
      t0D = wt * ( t12D + t21D); c1M[0][1] += t0D; c1M[1][0] += t0D;
      t0D = wt * ( t13D + t31D); c1M[0][2] += t0D; c1M[2][0] += t0D;
      t0D = wt * (-t23D + t32D); c1M[0][3] += t0D; c1M[3][0] += t0D;
      t0D = wt * ( t23D + t32D); c1M[1][2] += t0D; c1M[2][1] += t0D;
      t0D = wt * ( t13D - t31D); c1M[1][3] += t0D; c1M[3][1] += t0D;
      t0D = wt * (-t12D + t21D); c1M[2][3] += t0D; c1M[3][2] += t0D;
    }
    AlgMatrixScale(c1M, c1M, -2.0, 4, 4);
    AlgMatrixScale(c3M, c3M, 2.0, 4, 4);
    if(sumVW < DBL_EPSILON)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Now compute A from C1, C3 and \sum_i^n{\beta_i} */
  if(errNum == WLZ_ERR_NONE)
  {
    AlgMatrixTranspose(t0M, c3M, 4, 4);		/* T0 = C3^T */
    AlgMatrixMul(t1M, t0M, c3M, 4, 4, 4); 	/* T1 = C3^T C3 */
    t0D = 1.0 / (4.0 * sumVW);
    AlgMatrixScale(t1M, t1M, t0D, 4, 4);		
    AlgMatrixSub(aM, t1M, c1M, 4, 4);
    /* Find the eigenvector of A which has the greatest eigenvalue.
     * This is returned in the first column of aM. */
    AlgMatrixRSEigen(aM, 4, rM, 1);
    rM[0] = aM[0][0]; rM[1] = aM[1][0];
    rM[2] = aM[2][0]; rM[3] = aM[3][0];
    /* Compute the transform's rotation elements from the eigen vector. */
    t11D = rM[0] * rM[0];
    t22D = rM[1] * rM[1];
    t33D = rM[2] * rM[2];
    t44D = rM[3] * rM[3];
    t12D = 2.0 * rM[0] * rM[1];
    t13D = 2.0 * rM[0] * rM[2];
    t14D = 2.0 * rM[0] * rM[3];
    t23D = 2.0 * rM[1] * rM[2];
    t24D = 2.0 * rM[1] * rM[3];
    t34D = 2.0 * rM[2] * rM[3];
    aM[0][0] = t44D + t11D - t22D - t33D;
    aM[0][1] = t12D - t34D;
    aM[0][2] = t13D + t24D;
    aM[1][0] = t12D + t34D;
    aM[1][1] = t44D - t11D + t22D - t33D;
    aM[1][2] = t23D - t14D;
    aM[2][0] = t13D - t24D;
    aM[2][1] = t23D + t14D;
    aM[2][2] = t44D - t11D - t22D + t33D;
    /* Compute the translation elements from the eigen vector r, the sum of the
     * vertex weights \sum_i{\beta_i} and matrix C3. */
    /* Set t0M[0] to r (see Walker's paper). */
    t0M[0][0] = rM[0]; t0M[1][0] = rM[1]; t0M[2][0] = rM[2]; t0M[3][0] = rM[3];
    /* Compute s in t1M[0], but don't scale with -1/(2\sum_i{\beta_i} yet (see
     * Walker's paper). */
    AlgMatrixMul(t1M, c3M, t0M, 4, 4, 1);
    /* Set t0M to be W^T(r) (see Walker's paper). */
    t0M[0][0] =  rM[3]; t0M[0][1] =  rM[2];
    t0M[0][2] =  rM[1]; t0M[0][3] = -rM[0];
    t0M[1][0] = -rM[2]; t0M[1][1] =  rM[3];
    t0M[1][2] =  rM[0]; t0M[1][3] = -rM[1];
    t0M[2][0] =  rM[1]; t0M[2][1] = -rM[0];
    t0M[2][2] =  rM[3]; t0M[2][3] = -rM[2];
    t0M[3][0] =  rM[0]; t0M[3][1] =  rM[1];
    t0M[3][2] =  rM[2]; t0M[3][3] =  rM[3];
    /* Compute the product W^T(r) s (see Walker's paper). */
    AlgMatrixMul(t1M, t0M, t1M, 4, 4, 1);
    /* Extract the translation components and scale by -1/(2\sum_i{\beta_i}
     * (see Walker's paper). */
    t0D = -1.0 / (2.0 * sumVW);
    aM[0][3] = t0D * t1M[0][0];
    aM[1][3] = t0D * t1M[1][0];
    aM[2][3] = t0D * t1M[2][0];
    /* Set the transform's perspective and scale elements. */
    aM[3][0] = aM[3][1] = aM[3][2] = 0.0;
    aM[3][3] = 1.0;
    /* Create the affine transform. */
    trans = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_3D_AFFINE, aM, &errNum);
  }
  /* Free the matricies. */
  (void )AlcDouble2Free(aM);
  (void )AlcDouble2Free(c1M);
  (void )AlcDouble2Free(c3M);
  (void )AlcDouble2Free(t0M);
  (void )AlcDouble2Free(t1M);
  (void )AlcDouble2Free(t2M);
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
*		first set of vertices and normals onto the second
*		set. The vertex and normal weighting factors must
*		be in the range [0-1].
*		See WlzAffineTransformLSqDQ3D() from which this
*		function has been derived.
* \param	vtxType			Type of vertices.
* \param	nV			Number of vertices.
* \param	vW			Vertex weights (Walker's beta),
*					may be NULL which implies that
*					all the weights are 1.0.
* \param	v1			Vertices of the target set.
* \param	v0			Vertices of the source set.
* \param	nN			Numbr of normals.
* \param	nW			Normal weights (Walker's alpha).,
*					may be NULL which implies that
*					all the weights are 1.0
* \param	n1			Normals of the target set, may
*					be NULL.
* \param	n0			Normals of the source set, may
*					be NULL.
* \param	trType			Required transform type.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
static WlzAffineTransform *WlzAffineTransformLSqDQ2D(int nV, double *vW,
					WlzDVertex2 *v1, WlzDVertex2 *v0,
					int nN, double *nW,
					WlzDVertex2 *n1, WlzDVertex2 *n0,
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
  double	rM[4];
  double	**aM = NULL,		/* Walker's A */
  		**t0M = NULL,		/* Working matrix */
  		**t1M = NULL,		/* Working matrix */
  		**t2M = NULL,		/* Working matrix */
		**c1M = NULL,		/* Walker's C_1 */
		**c3M = NULL;		/* Walker's C_3 */
  WlzDVertex3	t0V,
  		t1V;
  WlzAffineTransform *trans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  /* Allocate matricies required to compute c1M, c2M and c3M. */
  if((AlcDouble2Calloc(&aM, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&c1M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&c3M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&t0M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&t1M, 4, 4) != ALC_ER_NONE) ||
     (AlcDouble2Calloc(&t2M, 4, 4) != ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  /* Compute c1M, c2M and c3M (see Walker's paper). */
  if(errNum == WLZ_ERR_NONE)
  {
    sumVW = 0.0;
    for(id0 = 0; id0 < nV; ++id0)
    {
      /* Update sumVW */
      sumVW += (vW)? vW[id0]: 1.0;
      /* Update c1M: Compute \beta_i{Q(v_i^t)}^TW(v_i_s) */
      wt = (vW)? vW[id0]: 1.0;
      t0V.vtX = (v0 + id0)->vtX; t0V.vtY = (v0 + id0)->vtY; t0V.vtZ = 0.0;
      t1V.vtX = (v1 + id0)->vtX; t1V.vtY = (v1 + id0)->vtY; t1V.vtZ = 0.0;
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      c1M[0][0] += wt * ( t11D - t22D);
      c1M[1][1] += wt * (-t11D + t22D);
      c1M[2][2] += wt * (-t11D - t22D);
      c1M[3][3] += wt * ( t11D + t22D);
      t0D = wt * ( t12D + t21D); c1M[0][1] += t0D; c1M[1][0] += t0D;
      t0D = wt * (-t12D + t21D); c1M[2][3] += t0D; c1M[3][2] += t0D;
      /* Update c3M: Compute \beta_i(W(v_i_s) - Q(v_i^t)) */
      t0D = wt * (-t1V.vtY - t0V.vtY); c3M[0][2] += t0D; c3M[2][0] -= t0D; 
      t0D = wt * ( t1V.vtX - t0V.vtX); c3M[0][3] += t0D; c3M[3][0] -= t0D;
      t0D = wt * ( t1V.vtX + t0V.vtX); c3M[1][2] += t0D; c3M[2][1] -= t0D;
      t0D = wt * ( t1V.vtY - t0V.vtY); c3M[1][3] += t0D; c3M[3][1] -= t0D;
    }
    for(id0 = 0; id0 < nN; ++id0)
    {
      /* Update c1M: Compute \alpha_i{Q(n_i^t)}^TW(n_i_s) */
      wt = (nW)? nW[id0]: 1.0;
      t0V.vtX = (n0 + id0)->vtX; t0V.vtY = (n0 + id0)->vtY; t0V.vtZ = 0.0;
      t1V.vtX = (n1 + id0)->vtX; t1V.vtY = (n1 + id0)->vtY; t1V.vtZ = 0.0;
      t11D = t0V.vtX * t1V.vtX;
      t12D = t0V.vtX * t1V.vtY;
      t21D = t0V.vtY * t1V.vtX;
      t22D = t0V.vtY * t1V.vtY;
      c1M[0][0] += wt * ( t11D - t22D);
      c1M[1][1] += wt * (-t11D + t22D);
      c1M[2][2] += wt * (-t11D - t22D);
      c1M[3][3] += wt * ( t11D + t22D);
      t0D = wt * ( t12D + t21D); c1M[0][1] += t0D; c1M[1][0] += t0D;
      t0D = wt * (-t12D + t21D); c1M[2][3] += t0D; c1M[3][2] += t0D;
    }
    AlgMatrixScale(c1M, c1M, -2.0, 4, 4);
    AlgMatrixScale(c3M, c3M, 2.0, 4, 4);
    if(sumVW < DBL_EPSILON)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Now compute A from C1, C3 and \sum_i^n{\beta_i} */
  if(errNum == WLZ_ERR_NONE)
  {
    AlgMatrixTranspose(t0M, c3M, 4, 4);		/* T0 = C3^T */
    AlgMatrixMul(t1M, t0M, c3M, 4, 4, 4); 	/* T1 = C3^T C3 */
    t0D = 1.0 / (4.0 * sumVW);
    AlgMatrixScale(t1M, t1M, t0D, 4, 4);		
    AlgMatrixSub(aM, t1M, c1M, 4, 4);
    /* Find the eigenvector of A which has the greatest eigenvalue.
     * This is returned in the first column of aM. */
    AlgMatrixRSEigen(aM, 4, rM, 1);
    rM[0] = aM[0][0]; rM[1] = aM[1][0];
    rM[2] = aM[2][0]; rM[3] = aM[3][0];
    /* Compute the transform's rotation elements from the eigen vector. */
    t11D = rM[0] * rM[0];
    t22D = rM[1] * rM[1];
    t33D = rM[2] * rM[2];
    t44D = rM[3] * rM[3];
    t12D = 2.0 * rM[0] * rM[1];
    t34D = 2.0 * rM[2] * rM[3];
    aM[0][0] = t44D + t11D - t22D - t33D;
    aM[0][1] = t12D - t34D;
    aM[1][0] = t12D + t34D;
    aM[1][1] = t44D - t11D + t22D - t33D;
    /* Compute the translation elements from the eigen vector r, the sum of the
     * vertex weights \sum_i{\beta_i} and matrix C3. */
    /* Set t0M[0] to r (see Walker's paper). */
    t0M[0][0] = rM[0]; t0M[1][0] = rM[1]; t0M[2][0] = rM[2]; t0M[3][0] = rM[3];
    /* Compute s in t1M[0], but don't scale with -1/(2\sum_i{\beta_i} yet (see
     * Walker's paper). */
    AlgMatrixMul(t1M, c3M, t0M, 4, 4, 1);
    /* Set t0M to be W^T(r) (see Walker's paper). */
    t0M[0][0] =  rM[3]; t0M[0][1] =  rM[2];
    t0M[0][2] =  rM[1]; t0M[0][3] = -rM[0];
    t0M[1][0] = -rM[2]; t0M[1][1] =  rM[3];
    t0M[1][2] =  rM[0]; t0M[1][3] = -rM[1];
    t0M[2][0] =  rM[1]; t0M[2][1] = -rM[0];
    t0M[2][2] =  rM[3]; t0M[2][3] = -rM[2];
    t0M[3][0] =  rM[0]; t0M[3][1] =  rM[1];
    t0M[3][2] =  rM[2]; t0M[3][3] =  rM[3];
    /* Compute the product W^T(r) s (see Walker's paper). */
    AlgMatrixMul(t1M, t0M, t1M, 4, 4, 1);
    /* Extract the translation components and scale by -1/(2\sum_i{\beta_i}
     * (see Walker's paper). */
    t0D = -1.0 / (2.0 * sumVW);
    aM[0][2] = t0D * t1M[0][0];
    aM[1][2] = t0D * t1M[1][0];
    /* Set the transform's perspective and scale elements. */
    aM[2][0] = aM[2][1] = 0.0;
    aM[2][2] = 1.0;
    /* Create the affine transform. */
    trans = WlzAffineTransformFromMatrix(WLZ_TRANSFORM_2D_AFFINE, aM, &errNum);
  }
  /* Free the matricies. */
  (void )AlcDouble2Free(aM);
  (void )AlcDouble2Free(c1M);
  (void )AlcDouble2Free(c3M);
  (void )AlcDouble2Free(t0M);
  (void )AlcDouble2Free(t1M);
  (void )AlcDouble2Free(t2M);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trans);
}
