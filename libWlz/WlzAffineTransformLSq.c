#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzAffineTransformLSq.c
* Date:         March 1999
* Author:       Richard Baldock, Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for computing Woolz affine transforms that
*		give the best fit, in a least squares sense, when
*		used to transform one set of verticies to another.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 13-12-00 bill Change members of WlzVertex and WlzVertexP.
* 05-12-00 bill In WlzAffineTransformLSqReg3D add code for degenerate
*		solutions.
* 29-11-00 bill Rename WlzAffineTransformLSq to WlzAffineTransformLSq2D
*		add WlzAffineTransformLSq3D, WlzAffineTransformLSqReg3D,
*		WlzAffineTransformLSqTrans3D and a new WlzAffineTransformLSq.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

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

/************************************************************************
* Function:	WlzAffineTransformLSq
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of verticies onto the second.
* Global refs:	-
* Parameters:	WlzVertexType vtxType:	Type of verticies.
*		int nVtx0:		Number of verticies in first
*					vector.
*		WlzVertexP vtxVec0:	First vector of verticies.
*		int nVtx1:		Number of verticies in second
*					vector (MUST be same as first).
*		WlzVertexP vtxVec1:	Second vector of verticies.
*		WlzTransformType trType: Required transform type.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzAffineTransformLSq3D
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of 3D verticies onto the second.
* Global refs:	-
* Parameters:	int nVtx0:		Number of verticies in first
*					vector.
*		WlzDVertex3 *vtxVec0:	First vector of verticies.
*		int nVtx1:		Number of verticies in second
*					vector (MUST be same as first).
*		WlzDVertex3 *vtxVec1:	Second vector of verticies.
*		WlzTransformType trType: Required transform type.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzAffineTransformLSq2D
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz affine transform which gives the
*		best (least squares) fit when used to transform the
*		first set of 2D verticies onto the second.
* Global refs:	-
* Parameters:	int nVtx0:		Number of verticies in first
*					vector.
*		WlzDVertex2 *vtxVec0:	First vector of verticies.
*		int nVtx1:		Number of verticies in second
*					vector (MUST be same as first).
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.
*		WlzTransformType trType: Required transform type.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzAffineTransformLSqReg3D
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz 3D registration transform which
*		gives the best (least squares) fit when used to
*		transform the first set of verticies onto the second.
*		The transform is constrained to rotation and
*		translation only.
*		The algorithm is based on Arun K.S., Huang T.T.
*		and Blostein S.D. "Least-Squares Fitting of Two 3-D
*		Point Sets" PAMI 9(5), 698-700, 1987.
* Global refs:	-
* Parameters:	WlzDVertex3 *pos0:	First vector of verticies.
*		WlzDVertex3 *pos1:	Second vector of verticies.
*		int nVtx:		Number of verticies in vectors.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
static WlzAffineTransform *WlzAffineTransformLSqReg3D(WlzDVertex3 *pos0,
		    		WlzDVertex3 *pos1, int nVtx,
				WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
  		idK,
		idR;
  double	tD0;
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

  if(((wMx = (double *)AlcCalloc(sizeof(double), 3)) == NULL) ||
     (AlcDouble2Calloc(&hMx, 3, 3) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&vMx, 3, 3) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&trMx, 4, 4) !=  ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
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
    /* Compute the 3x3 matrix hMx, which is the sum of tensor products of the
     * verticies relative to their centroids. */
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
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute 3x3 rotation matrix trMx = vMx.hMx', where hMx' is the transpose
     * of hMx. */
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
      /* Are source verticies (rel0) coplanar? They are iff one of the 3
       * singular values of hMx in wMx is zero. If the source verticies
       * are not coplanar, the the solution of the SVD is correct. */
      tI0 = ((fabs(*(wMx + 2)) >= DBL_EPSILON) << 2) |
            ((fabs(*(wMx + 1)) >= DBL_EPSILON) << 1) |
            (fabs(*(wMx + 2)) >= DBL_EPSILON);
      switch(tI0)
      {
	case 0:
	  /* Source verticies are not coplanar or colinear. The SVD gives the
	   * correct solution. */
	  break;
        case 1:
	case 2:
	case 4:
	  /* Source verticies are coplanar, but not colinear. There is a
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
	  /* Source verticies are colinear and there exists an infinity of
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/************************************************************************
* Function:	WlzAffineTransformLSqTrans3D
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz 3D registration transform which
*		gives the best (least squares) fit when used to
*		transform the first set of verticies onto the second.
*		The transform is constrained to translation only.
* Global refs:	-
* Parameters:	WlzDVertex3 *pos0:	First vector of verticies.
*		WlzDVertex3 *pos1:	Second vector of verticies.
*		int nVtx:		Number of verticies in vectors.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzAffineTransformLSqGen2D
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz general 2D affine transform which
*		gives the best (least squares) fit when used to
*		transform the first set of verticies onto the second.
* Global refs:	-
* Parameters:	WlzDVertex2 *vtxVec0:	First vector of verticies.
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.
*		int nVtx:		Number of verticies in vectors.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzAffineTransformLSqReg2D
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz 2D registration transform which
*		gives the best (least squares) fit when used to
*		transform the first set of verticies onto the second.
*		The transform is constrained to rotation and
*		translation only.
* Global refs:	-
* Parameters:	WlzDVertex2 *vtxVec0:	First vector of verticies.
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.
*		int nVtx:		Number of verticies in vectors.
*		WlzTransformType trType: Required transform type.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzAffineTransformLSqTrans2D
* Returns:	WlzAffineTransform:	Computed affine transform, may
*					be NULL on error.
* Purpose:	Computes the Woolz 2D translation transform which
*		gives the best (least squares) fit when used to
*		transform the first set of verticies onto the second.
*		The transform is constrained to rotation and
*		translation only.
* Global refs:	-
* Parameters:	WlzDVertex2 *vtxVec0:	First vector of verticies.
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.
*		int nVtx:		Number of verticies in vectors.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
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
