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
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static WlzAffineTransform *WlzAffineTransformLSqGen2D(WlzDVertex2 *vtxVec0,
				WlzDVertex2 *vtxVec1, int nVtx,
				WlzErrorNum *dstErr),
		*WlzAffineTransformLSqReg2D(WlzDVertex2 *vtxVec0,
				WlzDVertex2 *vtxVec1, int nVtx,
				WlzTransformType trType,
				WlzErrorNum *dstErr),
		*WlzAffineTransformLSqTrans2D(WlzDVertex2 *vtxVec0,
				WlzDVertex2 *vtxVec1, int nVtx,
				WlzErrorNum *dstErr);

/************************************************************************
* Function:	WlzAffineTransformLSq					*
* Returns:	WlzAffineTransform:	Computed affine transform, may	*
*					be NULL on error.		*
* Purpose:	Computes the Woolz affine transform which gives the	*
*		best (least squares) fit when used to transform the	*
*		first set of verticies onto the second.			*
* Global refs:	-							*
* Parameters:	int nVtx0:		Number of verticies in first	*
*					vector.				*
*		WlzDVertex2 *vtxVec0:	First vector of verticies.	*
*		int nVtx1:		Number of verticies in second	*
*					vector (MUST be same as first).	*
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.	*
*		WlzTransformType trType: Required transform type.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
************************************************************************/
WlzAffineTransform *WlzAffineTransformLSq(int nVtx0, WlzDVertex2 *vtxVec0,
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
  return( trans );
}

/************************************************************************
* Function:	WlzAffineTransformLSqGen2D				*
* Returns:	WlzAffineTransform:	Computed affine transform, may	*
*					be NULL on error.		*
* Purpose:	Computes the Woolz general 2D affine transform which	*
*		gives the best (least squares) fit when used to 	*
*		transform the first set of verticies onto the second.	*
* Global refs:	-							*
* Parameters:	WlzDVertex2 *vtxVec0:	First vector of verticies.	*
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.	*
*		int nVtx:		Number of verticies in vectors.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
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
    trans = WlzAffineTransformFromMatrix4X4(WLZ_TRANSFORM_2D_AFFINE, trMat,
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
* Function:	WlzAffineTransformLSqReg2D				*
* Returns:	WlzAffineTransform:	Computed affine transform, may	*
*					be NULL on error.		*
* Purpose:	Computes the Woolz 2D registration transform which	*
*		gives the best (least squares) fit when used to 	*
*		transform the first set of verticies onto the second.	*
*		The transform is constrained to rotation and 		*
*		translation only.					*
* Global refs:	-							*
* Parameters:	WlzDVertex2 *vtxVec0:	First vector of verticies.	*
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.	*
*		int nVtx:		Number of verticies in vectors.	*
*		WlzTransformType trType: Required transform type.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
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
    trans = WlzAffineTransformFromMatrix4X4(WLZ_TRANSFORM_2D_AFFINE, trMat,
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
* Function:	WlzAffineTransformLSqTrans2D				*
* Returns:	WlzAffineTransform:	Computed affine transform, may	*
*					be NULL on error.		*
* Purpose:	Computes the Woolz 2D translation transform which	*
*		gives the best (least squares) fit when used to 	*
*		transform the first set of verticies onto the second.	*
*		The transform is constrained to rotation and 		*
*		translation only.					*
* Global refs:	-							*
* Parameters:	WlzDVertex2 *vtxVec0:	First vector of verticies.	*
*		WlzDVertex2 *vtxVec1:	Second vector of verticies.	*
*		int nVtx:		Number of verticies in vectors.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
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
  trans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE, 
  				     sum.vtX / nVtx, sum.vtY / nVtx, 0.0,
				     1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, &errNum);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trans);
}
