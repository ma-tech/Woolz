#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzAffineTransform.c
* Date:         March 1999
* Author:       Richard Baldock, Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for computing Woolz affine transforms and
*		applying them to Woolz objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 27-09-00 bill Make modifications for 3D affine transforms.
*		Add WlzAffineTransformDimension(),
*               WlzAffineTransformTranslationSet(),
*               WlzAffineTransformFromTranslation(),
*               WlzAffineTransformRotationSet(),
*               WlzAffineTransformFromScale(),
*               WlzAffineTransformScaleSet() and
*               WlzAffineTransformFromRotation().
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static int			WlzAffineTransformIsTranslate2(
				  WlzAffineTransform *trans,
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
static int			WlzAffineTransformIsTranslate3(
				  WlzAffineTransform *trans,
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzAffineTransformIntTranslate(
				  WlzObject *,
				  WlzAffineTransform *,
				  WlzErrorNum *);
static WlzPolygonDomain 	*WlzAffineTransformPoly2(
				  WlzPolygonDomain *,
				  WlzAffineTransform *,
				  WlzErrorNum *);
static WlzBoundList 		*WlzAffineTransformBoundList(WlzBoundList *,
				  WlzAffineTransform *,
				  WlzErrorNum *);
static WlzErrorNum 		WlzAffineTransformValues2(WlzObject *,
				  WlzObject *,
				  WlzAffineTransform *,
				  WlzInterpolationType);
static WlzErrorNum 		WlzAffineTransformPrimSet2(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim prim);
static WlzErrorNum 		WlzAffineTransformPrimSet3(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim prim);
static void			WlzAffineTransformPrimGet2(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim *prim);

/************************************************************************
* Function:	WlzAffineTransformDimension
* Returns:	int:			2 or 3 for a 2D or 3D affine
*					transform, 0 on error.
* Purpose:	Computes the dimension of the given affine transform.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given affine transform.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
int		WlzAffineTransformDimension(WlzAffineTransform *tr,
					    WlzErrorNum *dstErr)
{
  int		dim = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(tr->type)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
      case WLZ_TRANSFORM_2D_REG:
      case WLZ_TRANSFORM_2D_TRANS:
      case WLZ_TRANSFORM_2D_NOSHEAR:
	dim = 2;
        break;
      case WLZ_TRANSFORM_3D_AFFINE:
      case WLZ_TRANSFORM_3D_REG:
      case WLZ_TRANSFORM_3D_TRANS:
      case WLZ_TRANSFORM_3D_NOSHEAR:
	dim = 3;
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
  return(dim);
}

/************************************************************************
* Function:	WlzAffineTransformIsTranslate
* Returns:	int:			Non-zero if translation.
* Purpose:	Tests whether the given affine transform is a simple
*		integer translation.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given affine transform.
*		WlzObject *obj:		Optional object, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
int		WlzAffineTransformIsTranslate(WlzAffineTransform *tr,
					      WlzObject *obj,
					      WlzErrorNum *dstErr)
{
  int		transFlg = 1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(WlzAffineTransformDimension(tr, NULL))
    {
      case 2:
	transFlg = WlzAffineTransformIsTranslate2(tr, obj, &errNum);
        break;
      case 3:
	transFlg = WlzAffineTransformIsTranslate3(tr, obj, &errNum);
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
  return(transFlg);
}

/************************************************************************
* Function:	WlzAffineTransformIsTranslate2D
* Returns:	int:			Non-zero if translation.
* Purpose:	Tests wether the given 2D affine transform is a simple
*		integer translation.
*		Because this is a static function the parameters are
*		not checked.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D affine transform.
*		WlzObject *obj:		Optional object, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static int	WlzAffineTransformIsTranslate2(WlzAffineTransform *tr,
						WlzObject *obj,
						WlzErrorNum *dstErr)
{
  int		idx,
  		transFlg = 1;
  double	trX,
  		trY;
  WlzDVertex2	*tV;
  WlzIBox2	box;
  WlzDVertex2	tstV[4];
  const double	transDelta = 0.1,
  		nonTransDelta = 1.0E-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* First check for integer translation matrix elements. */
  trX = tr->mat[0][2];
  trY = tr->mat[1][2];
  if((fabs(trX - WLZ_NINT(trX)) > transDelta) ||
     (fabs(trY - WLZ_NINT(trY)) > transDelta))
  {
    transFlg = 0;
  }
  else
  {
    if(obj)
    {
      /* Given an object: Check for integer translation within the bounding
       * box of the given object. */
      box = WlzBoundingBox2D(obj, &errNum);
      tstV[0].vtX = box.xMin; tstV[0].vtY = box.yMin;
      tstV[1].vtX = box.xMin; tstV[1].vtY = box.yMax;
      tstV[2].vtX = box.xMax; tstV[2].vtY = box.yMin;
      tstV[3].vtX = box.xMax; tstV[3].vtY = box.yMax;
      idx = 0;
      do
      {
	tV = tstV + idx;
	trX = (tV->vtX * tr->mat[0][0]) + (tV->vtY * tr->mat[0][1]);
	trY = (tV->vtX * tr->mat[1][0]) + (tV->vtY * tr->mat[1][1]);
	transFlg = (fabs(trX - tV->vtX) <= transDelta) &&
	           (fabs(trY - tV->vtY) <= transDelta);
      }
      while(transFlg && (++idx < 4));
    }
    else
    {
      /* Not given an object: Check the rest of the matrix elemets. */
      transFlg = (fabs(tr->mat[0][0] - 1.0) <= nonTransDelta) &&
		 (fabs(tr->mat[0][1]) <= nonTransDelta) &&
		 (fabs(tr->mat[1][0]) <= nonTransDelta) &&
		 (fabs(tr->mat[1][1] - 1.0) <= nonTransDelta) &&
		 (fabs(tr->mat[2][0]) <= nonTransDelta) &&
		 (fabs(tr->mat[2][1]) <= nonTransDelta) &&
		 (fabs(tr->mat[2][2] - 1.0) <= nonTransDelta);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(transFlg);
}

/************************************************************************
* Function:	WlzAffineTransformIsTranslate3
* Returns:	int:			Non-zero if translation.
* Purpose:	Tests wether the given 3D affine transform is a simple
*		integer translation.
*		Because this is a static function the parameters are
*		not checked.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 3D affine transform.
*		WlzObject *obj:		Optional object, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static int	WlzAffineTransformIsTranslate3(WlzAffineTransform *tr,
					       WlzObject *obj,
					       WlzErrorNum *dstErr)
{
  int		idx,
  		transFlg = 1;
  double	trX,
  		trY,
		trZ;
  WlzDVertex3	*tV;
  WlzIBox3	box;
  WlzDVertex3	tstV[8];
  const double	transDelta = 0.1,
  		nonTransDelta = 1.0E-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* First check for integer translation matrix elements. */
  trX = tr->mat[0][3];
  trY = tr->mat[1][3];
  trZ = tr->mat[2][3];
  if((fabs(trX - WLZ_NINT(trX)) > transDelta) ||
     (fabs(trY - WLZ_NINT(trY)) > transDelta) ||
     (fabs(trZ - WLZ_NINT(trZ)) > transDelta))
  {
    transFlg = 0;
  }
  else
  {
    if(obj)
    {
      /* Given an object: Check for integer translation within the bounding
       * box of the given object. */
      box = WlzBoundingBox3D(obj, &errNum);
      tstV[0].vtX = box.xMin; tstV[0].vtY = box.yMin; tstV[0].vtZ = box.zMin;
      tstV[1].vtX = box.xMax; tstV[1].vtY = box.yMin; tstV[1].vtZ = box.zMin;
      tstV[2].vtX = box.xMin; tstV[2].vtY = box.yMax; tstV[2].vtZ = box.zMin;
      tstV[3].vtX = box.xMax; tstV[3].vtY = box.yMax; tstV[3].vtZ = box.zMin;
      tstV[4].vtX = box.xMin; tstV[4].vtY = box.yMin; tstV[4].vtZ = box.zMax;
      tstV[5].vtX = box.xMax; tstV[5].vtY = box.yMin; tstV[5].vtZ = box.zMax;
      tstV[6].vtX = box.xMin; tstV[6].vtY = box.yMax; tstV[6].vtZ = box.zMax;
      tstV[7].vtX = box.xMin; tstV[7].vtY = box.yMin; tstV[7].vtZ = box.zMax;
      idx = 0;
      do
      {
	tV = tstV + idx;
	trX = (tV->vtX * tr->mat[0][0]) + (tV->vtY * tr->mat[0][1]) +
	      (tV->vtZ * tr->mat[0][2]);
	trY = (tV->vtX * tr->mat[1][0]) + (tV->vtY * tr->mat[1][1]) +
	      (tV->vtZ * tr->mat[1][2]);
	trZ = (tV->vtX * tr->mat[2][0]) + (tV->vtY * tr->mat[2][1]) +
	      (tV->vtZ * tr->mat[2][2]);
	transFlg = (fabs(trX - tV->vtX) <= transDelta) &&
	           (fabs(trY - tV->vtY) <= transDelta) &&
		   (fabs(trZ - tV->vtZ) <= transDelta);
      }
      while(transFlg && (++idx < 8));
    }
    else
    {
      /* Not given an object: Check the rest of the matrix elemets. */
      transFlg = (fabs(tr->mat[0][0] - 1.0) <= nonTransDelta) &&
		 (fabs(tr->mat[0][1]) <= nonTransDelta) &&
		 (fabs(tr->mat[0][2]) <= nonTransDelta) &&
		 (fabs(tr->mat[1][0]) <= nonTransDelta) &&
		 (fabs(tr->mat[1][1] - 1.0) <= nonTransDelta) &&
		 (fabs(tr->mat[1][3]) <= nonTransDelta) &&
		 (fabs(tr->mat[2][0]) <= nonTransDelta) &&
		 (fabs(tr->mat[2][1]) <= nonTransDelta) &&
		 (fabs(tr->mat[2][2] - 1.0) <= nonTransDelta) &&
		 (fabs(tr->mat[3][0]) <= nonTransDelta) &&
		 (fabs(tr->mat[3][1]) <= nonTransDelta) &&
		 (fabs(tr->mat[3][2]) <= nonTransDelta) &&
		 (fabs(tr->mat[3][3] - 1.0) <= nonTransDelta);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(transFlg);
}

/************************************************************************
* Function:	WlzAffineTransformIntTranslate
* Returns:	WlzObject *:		Translated object or NULL on
*					error.
* Purpose:	Translates the given 2D domain object with an integral
*		translation.
*		Because this is a static function the parameters are
*		not checked.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given 2D domain object.
*		WlzAffineTransform *tr: Given affine transform.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
static WlzObject *WlzAffineTransformIntTranslate(WlzObject *srcObj,
					         WlzAffineTransform *tr,
					      	 WlzErrorNum *dstErr)
{
  int		trX,
  		trY,
		trZ;
  WlzObject	*newObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(WlzAffineTransformDimension(tr, NULL))
  {
    case 2:
      trX = WLZ_NINT(tr->mat[0][2]);
      trY = WLZ_NINT(tr->mat[1][2]);
      trZ = 0;
      newObj = WlzShiftObject(srcObj, trX, trY, trZ, &errNum);
      break;
    case 3:
      trX = WLZ_NINT(tr->mat[0][3]);
      trY = WLZ_NINT(tr->mat[1][3]);
      trZ = WLZ_NINT(tr->mat[2][3]);
      newObj = WlzShiftObject(srcObj, trX, trY, trZ, &errNum);
      break;
    default:
      errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newObj);
}

/************************************************************************
* Function:     WlzAffineTransformPoly2
* Returns:      WlzPolygonDomain *:     Transformed polygon domain or
*                                       NULL on error.
* Purpose:      Transforms the given polygon domain.
*               Because this is a static function the parameters (other
*               than polygon type) are not checked.
* Global refs:  -
* Parameters:   WlzPolygonDomain *srcPoly: Given polygon domain.
*               WlzAffineTransform *trans: Given affine transform.
*               WlzErrorNum *dstErr:    Destination pointer for error
*                                       number.
************************************************************************/
static WlzPolygonDomain *WlzAffineTransformPoly2(WlzPolygonDomain *srcPoly,
                                                 WlzAffineTransform *trans,
                                                 WlzErrorNum *dstErr)
{
  int           count;
  double        cx,
                cy,
                dx,
                dy,
                sx,
                sy,
                tx,
                ty;
  WlzIVertex2   *srcVtxI,
                *dstVtxI;
  WlzFVertex2   *srcVtxF,
                *dstVtxF;
  WlzDVertex2   *srcVtxD,
                *dstVtxD;
  WlzPolygonDomain *dstPoly = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((srcPoly->type != WLZ_POLYGON_INT) &&
     (srcPoly->type != WLZ_POLYGON_FLOAT) &&
     (srcPoly->type != WLZ_POLYGON_DOUBLE))
  {
    errNum = WLZ_ERR_POLYGON_TYPE;
  }
  else if(WlzAffineTransformDimension(trans, NULL) != 2)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if((dstPoly = WlzMakePolyDmn(srcPoly->type, NULL, 0,
                                    srcPoly->nvertices, 1, &errNum)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    cx = trans->mat[0][0];
    cy = trans->mat[1][1];
    sx = trans->mat[0][1];
    sy = trans->mat[1][0];
    tx = trans->mat[0][2];
    ty = trans->mat[1][2];
    dstPoly->nvertices = srcPoly->nvertices;
    switch(srcPoly->type)
    {
      case WLZ_POLYGON_INT:
        srcVtxI = srcPoly->vtx;
        dstVtxI = dstPoly->vtx;
        count = srcPoly->nvertices;
        while(count-- > 0)
        {
          dx = (cx * srcVtxI->vtX) + (sx * srcVtxI->vtY) + tx;
          dy = (sy * srcVtxI->vtX) + (cy * srcVtxI->vtY) + ty;
          dstVtxI->vtX = WLZ_NINT(dx);
          dstVtxI->vtY = WLZ_NINT(dy);
          ++srcVtxI;
          ++dstVtxI;
        }
        break;
      case WLZ_POLYGON_FLOAT:
        srcVtxF = (WlzFVertex2 *)(srcPoly->vtx);
        dstVtxF = (WlzFVertex2 *)(dstPoly->vtx);
        count = srcPoly->nvertices;
        while(count-- > 0)
        {
          dstVtxF->vtX = (cx * srcVtxF->vtX) + (sx * srcVtxF->vtY) + tx;
          dstVtxF->vtY = (sy * srcVtxF->vtX) + (cy * srcVtxF->vtY) + ty;
          ++srcVtxF;
          ++dstVtxF;
        }
        break;
      case WLZ_POLYGON_DOUBLE:
        srcVtxD = (WlzDVertex2 *)(srcPoly->vtx);
        dstVtxD = (WlzDVertex2 *)(dstPoly->vtx);
        count = srcPoly->nvertices;
        while(count-- > 0)
        {
          dstVtxF->vtX = (cx * srcVtxF->vtX) + (sx * srcVtxF->vtY) + tx;
          dstVtxF->vtY = (sy * srcVtxF->vtX) + (cy * srcVtxF->vtY) + ty;
          ++srcVtxD;
          ++dstVtxD;
        }
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstPoly);
}

/************************************************************************
* Function:	WlzAffineTransformBoundList
* Returns:	WlzBoundList *:		Transformed boundary list or
*					NULL on error.
* Purpose:	Transforms the given boundary list.
*		Because this is a static function the parameters are
*		not checked.
* Global refs:	-
* Parameters:	WlzBoundList *srcBound: Given boundary list.
*		WlzAffineTransform *trans: Given affine transform.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
static WlzBoundList *WlzAffineTransformBoundList(WlzBoundList *srcBound,
					        WlzAffineTransform *trans,
					      	WlzErrorNum *dstErr)
{
  WlzDomain	dumDom;
  WlzBoundList	*dstBound = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstBound = (WlzBoundList *)AlcCalloc(sizeof(WlzBoundList), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else if(WlzAffineTransformDimension(trans, NULL) != 2)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    dstBound->type = srcBound->type;
    dstBound->wrap = srcBound->wrap;
    /* transform the polygon */
    if((dstBound->poly = WlzAffineTransformPoly2(srcBound->poly, trans,
    						  &errNum)) != NULL)
    {
      /* transform next */
      if(srcBound->next)
      {
	if((dumDom.b = WlzAffineTransformBoundList(srcBound->next, trans,
					           &errNum)) != NULL)
	{
	  (void )WlzAssignDomain(dumDom, &errNum);
	  dstBound->next = dumDom.b;
	}
      }
      /* transform down */
      if(srcBound->down && (errNum == WLZ_ERR_NONE))
      {
	if((dumDom.b = WlzAffineTransformBoundList(srcBound->down, trans,
						   &errNum)) != NULL)
	{
	  (void )WlzAssignDomain(dumDom, &errNum);
	  dstBound->down = dumDom.b;
	  dstBound->down->up = dstBound;
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreeBoundList(dstBound);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBound);
}


/************************************************************************
* Function:	WlzAffineTransformValues2
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Creates a new value table, fills in the values and
*		adds it to the given new object.
*		Because this is a static function the parameters are
*		not checked.
* Global refs:	-
* Parameters:	WlzObject *newObj:	Partialy transformed object
*					with a valid domain.
*		WlzObject *srcObj:	2D domain object which is being
*					transformed.
*		WlzAffineTransform *trans: Given affine transform.
*		WlzInterpolationType interp: Level of interpolation to
*					use.
************************************************************************/
static WlzErrorNum WlzAffineTransformValues2(WlzObject *newObj,
					     WlzObject *srcObj,
					     WlzAffineTransform *trans,
					     WlzInterpolationType interp)
{
  int		count;
  double	tD0,
  		tD1,
		dx,
		dy,
		tx,
  		ty,
		sx,
		cx,
		sy,
		cy,
		lxyy,
		lyyy;
  WlzIVertex2	posI;
  WlzGreyType	newGreyType;
  WlzPixelV	bkdV;
  WlzValues	newValues;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzAffineTransform *invTrans = NULL;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  newValues.core = NULL;
  if(WlzAffineTransformDimension(trans, NULL) != 2)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if((invTrans = WlzAffineTransformInverse(trans, &errNum)) != NULL)
  {
    bkdV = WlzGetBackground(srcObj, &errNum);
    newGreyType = WlzGreyTableTypeToGreyType(srcObj->values.v->type,
					     &errNum);
    if(bkdV.type != newGreyType)
    {
      errNum = WlzValueConvertPixel(&bkdV, bkdV, newGreyType);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newValues.v = WlzNewValueTb(newObj, srcObj->values.v->type,
				bkdV, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newObj->values = WlzAssignValues(newValues, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Set up back transformation parameters */
    cx = invTrans->mat[0][0];
    cy = invTrans->mat[1][1];
    sx = invTrans->mat[0][1];
    sy = invTrans->mat[1][0];
    tx = invTrans->mat[0][2];
    ty = invTrans->mat[1][2];
    /* Fill in grey values */
    errNum = WlzInitGreyScan(newObj, &iWSp, &gWSp);
    if(errNum == WLZ_ERR_NONE)
    {
      gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);	
    }
    while((errNum == WLZ_ERR_NONE) &&
	  (WlzNextGreyInterval(&iWSp) == 0))
    {
      posI.vtX = iWSp.lftpos;
      posI.vtY = iWSp.linpos;
      lxyy = (sx * posI.vtY) + tx;
      lyyy = (cy * posI.vtY) + ty;
      count = iWSp.rgtpos - iWSp.lftpos + 1;
      switch(interp)
      {
	case WLZ_INTERPOLATION_NEAREST:
	  while(count-- > 0)
	  {
	    dx = lxyy + (cx * posI.vtX);
	    dy = lyyy + (sy * posI.vtX);
	    WlzGreyValueGet(gVWSp, 0, dy, dx);
	    switch(gWSp.pixeltype)
	    {
	      case WLZ_GREY_INT:
	        *(gWSp.u_grintptr.inp)++ = (*(gVWSp->gVal)).inv;
		break;
	      case WLZ_GREY_SHORT:
	        *(gWSp.u_grintptr.shp)++ = (*(gVWSp->gVal)).shv;
		break;
	      case WLZ_GREY_UBYTE:
	        *(gWSp.u_grintptr.ubp)++ = (*(gVWSp->gVal)).ubv;
		break;
	      case WLZ_GREY_FLOAT:
	        *(gWSp.u_grintptr.flp)++ = (*(gVWSp->gVal)).flv;
		break;
	      case WLZ_GREY_DOUBLE:
	        *(gWSp.u_grintptr.dbp)++ = (*(gVWSp->gVal)).dbv;
		break;
	    }
	    ++(posI.vtX);
	  }
	  break;
	case WLZ_INTERPOLATION_LINEAR:
	  while(count-- > 0)
	  {
	    dx = lxyy + (cx * posI.vtX);
	    dy = lyyy + (sy * posI.vtX);
	    WlzGreyValueGetCon(gVWSp, 0, dy, dx);
	    tD0 = dx - WLZ_NINT(dx-0.5);
	    tD1 = dy - WLZ_NINT(dy-0.5);
	    switch(gWSp.pixeltype)
	    {
	      case WLZ_GREY_INT:
		tD0 = (((gVWSp->gVal[0]).inv * (1.0 - tD0) * (1.0 - tD1)) +
		       ((gVWSp->gVal[1]).inv * tD0 * (1.0 - tD1)) +
		       ((gVWSp->gVal[2]).inv * (1.0 - tD0) * tD1) +
		       ((gVWSp->gVal[3]).inv * tD0 * tD1));
		*(gWSp.u_grintptr.inp)++ = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_SHORT:
		tD0 = (((gVWSp->gVal[0]).shv * (1.0 - tD0) * (1.0 - tD1)) +
		       ((gVWSp->gVal[1]).shv * tD0 * (1.0 - tD1)) +
		       ((gVWSp->gVal[2]).shv * (1.0 - tD0) * tD1) +
		       ((gVWSp->gVal[3]).shv * tD0 * tD1));
		*(gWSp.u_grintptr.shp)++ = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_UBYTE:
		tD0 = (((gVWSp->gVal[0]).ubv * (1.0 - tD0) * (1.0 - tD1)) +
		       ((gVWSp->gVal[1]).ubv * tD0 * (1.0 - tD1)) +
		       ((gVWSp->gVal[2]).ubv * (1.0 - tD0) * tD1) +
		       ((gVWSp->gVal[3]).ubv * tD0 * tD1));
		*(gWSp.u_grintptr.ubp)++ = WLZ_CLAMP(tD0, 0, 255);
		break;
	      case WLZ_GREY_FLOAT:
		tD0 = (((gVWSp->gVal[0]).flv * (1.0 - tD0) * (1.0 - tD1)) +
		       ((gVWSp->gVal[1]).flv * tD0 * (1.0 - tD1)) +
		       ((gVWSp->gVal[2]).flv * (1.0 - tD0) * tD1) +
		       ((gVWSp->gVal[3]).flv * tD0 * tD1));
		*(gWSp.u_grintptr.flp)++ = tD0;
		break;
	      case WLZ_GREY_DOUBLE:
		tD0 = (((gVWSp->gVal[0]).dbv * (1.0 - tD0) * (1.0 - tD1)) +
		       ((gVWSp->gVal[1]).dbv * tD0 * (1.0 - tD1)) +
		       ((gVWSp->gVal[2]).dbv * (1.0 - tD0) * tD1) +
		       ((gVWSp->gVal[3]).dbv * tD0 * tD1));
		*(gWSp.u_grintptr.dbp)++ = tD0;
		break;
	    }
	    ++(posI.vtX);
	  }
	  break;
	default:
	  errNum = WLZ_ERR_INTERPOLATION_TYPE;
	  break;
      }
    }
    if(errNum == WLZ_ERR_EOO)	        /* Reset error from end of intervals */ 
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  WlzGreyValueFreeWSp(gVWSp);
  if(invTrans)
  {
    (void )WlzFreeAffineTransform(invTrans);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformPrimSet
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Sets the given transform's matrix from an
*		affine transform primitives data structure.
*		A composite transform is built from the primitives
*		with the order of composition being scale (applied first),
*		shear, rotation and then translation (applied last),
*		ie:
*		  A = T.R.Sh.Sc, x' = A.x
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given affine transform.
*		WlzAffineTransformPrim: Given primitives.
************************************************************************/
WlzErrorNum	WlzAffineTransformPrimSet(WlzAffineTransform *tr,
					  WlzAffineTransformPrim prim)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(WlzAffineTransformDimension(tr, NULL))
    {
      case 2:
	errNum = WlzAffineTransformPrimSet2(tr, prim);
        break;
      case 3:
	errNum = WlzAffineTransformPrimSet3(tr, prim);
        break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
        break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformPrimSet2
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Sets the given transform's matrix from an
*		affine transform primitives data structure.
*		A composite transform is built from the primitives
*		with the order of composition being scale (applied first),
*		shear, rotation and then translation (applied last),
*		ie:
*		  A = T.R.Sh.Sc, x' = A.x
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D affine transform.
*		WlzAffineTransformPrim: Given primitives.
************************************************************************/
static WlzErrorNum WlzAffineTransformPrimSet2(WlzAffineTransform *tr,
					      WlzAffineTransformPrim prim)
{
  int		idx0;
  double	tS,
  		tCos1,
		tCos2,
  		tSin1,
		tSin2,
		tCosSin2;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tCos1 = cos(prim.theta);
  tSin1 = sin(prim.theta);
  tCos2 = cos(prim.psi);
  tSin2 = sin(prim.psi);
  tCosSin2 = tCos2 * tSin2 * prim.alpha;
  tSin2 *= tSin2 * prim.alpha;
  tCos2 *= tCos2 * prim.alpha;
  tS = prim.scale;
  tr->mat[0][0] = tS * ((tCos1 * (1 - tCosSin2)) + (tSin1 * tSin2));
  tr->mat[0][1] = tS * ((-tSin1 * (1 + tCosSin2)) + (tCos1 * tCos2));
  tr->mat[1][0] = tS * ((tSin1 * (1 - tCosSin2)) - (tCos1 * tSin2));
  tr->mat[1][1] = tS * ((tCos1 * (1 + tCosSin2)) + (tSin1 * tCos2));
  tr->mat[0][2] = prim.tx;
  tr->mat[1][2] = prim.ty;
  for(idx0 = 0; idx0 < 4; ++idx0)
  {
    tr->mat[2][idx0] = 0.0;
    tr->mat[3][idx0] = 0.0;
    tr->mat[idx0][3] = 0.0;
  }
  tr->mat[2][2] = 1.0;
  tr->mat[3][3] = 1.0;
  if(prim.invert)
  {
    for(idx0 = 0; idx0 < 4; ++idx0)
    {
      tr->mat[0][idx0] *= -1.0;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformPrimSet3
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Sets the given transform's matrix from an
*		affine transform primitives data structure.
*		A composite transform is built from the primitives
*		with the order of composition being scale (applied first),
*		shear, rotation and then translation (applied last),
*		ie:
*		  A = T.R.Sh.Sc, x' = A.x
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D affine transform.
*		WlzAffineTransformPrim: Given primitives.
************************************************************************/
static WlzErrorNum WlzAffineTransformPrimSet3(WlzAffineTransform *tr,
					      WlzAffineTransformPrim prim)
{
  double	cx,
		sx,
		cy,
		sy,
		cz,
  		sz;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((fabs(prim.alpha) > DBL_EPSILON) || (fabs(prim.psi) > DBL_EPSILON) ||
     (fabs(prim.xsi) > DBL_EPSILON))
  {
   /* TODO: Shear not yet implemented for 3D! So if a shear is specified
    * return an error. */
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    sy = sin(prim.phi); cy = cos(prim.phi);
    sz = sin(prim.theta); cz = cos(prim.theta);
    tr->mat[0][0] = prim.scale * cy * cz;
    tr->mat[0][1] = -(prim.scale * sz);
    tr->mat[0][2] = prim.scale * sy * cz;
    tr->mat[0][3] = prim.scale * prim.tx;
    tr->mat[1][0] = prim.scale * cy * sz;
    tr->mat[1][1] = prim.scale * cz;
    tr->mat[1][2] = prim.scale * sy * sz;
    tr->mat[1][3] = prim.scale * prim.ty;
    tr->mat[2][0] = -(prim.scale * sy);
    tr->mat[2][1] = 0.0;
    tr->mat[2][2] = prim.scale * cy;
    tr->mat[2][3] = prim.scale * prim.tz;
    tr->mat[3][0] = tr->mat[3][1] = tr->mat[3][2] = 0.0;
    tr->mat[3][3] = 1.0;
    if(prim.invert)
    {
      tr->mat[0][0] = -(tr->mat[0][0]);
      tr->mat[0][1] = -(tr->mat[0][1]);
      tr->mat[0][2] = -(tr->mat[0][2]);
      tr->mat[0][3] = -(tr->mat[0][3]);
    }
  }
  return(errNum);
}


/************************************************************************
* Function:	WlzAffineTransformTranslationSet
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Sets the given transform's matrix from the given
*		translations.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D or 3D affine transform.
*		double tx: 		Translation along the x-axis.
*		double ty: 		Translation along the y-axis.
*		double tz: 		Translation along the z-axis,
*					ignored for 2D transforms.
************************************************************************/
WlzErrorNum	WlzAffineTransformTranslationSet(WlzAffineTransform *tr,
					   double tx, double ty, double tz)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(WlzAffineTransformDimension(tr, NULL))
    {
      case 2:
	tr->mat[0][0] = tr->mat[1][1] = tr->mat[2][2] = 1.0;
	tr->mat[0][1] = tr->mat[1][0] = tr->mat[2][0] = tr->mat[2][1] = 0.0;
	tr->mat[0][2] = tx;
	tr->mat[1][2] = ty;
        break;
      case 3:
	tr->mat[0][0] = tr->mat[1][1] = tr->mat[2][2] = tr->mat[3][3] = 1.0;
	tr->mat[0][1] = tr->mat[0][2] = tr->mat[1][0] = 0.0;
	tr->mat[1][2] = tr->mat[2][0] = tr->mat[2][1] = 0.0;
	tr->mat[3][0] = tr->mat[3][1] = tr->mat[3][2] = 0.0;
	tr->mat[0][3] = tx;
	tr->mat[1][3] = ty;
	tr->mat[2][3] = tz;
        break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformFromTranslation
* Returns:	WlzAffineTransform *:	New affine transform,
*					NULL on error.
* Purpose:	Constructs a new affine transform from the given 
*		translations.
* Global refs:	-
* Parameters:	WlzTransformType type:	Required transform type.
* 		double tx: 		Translation along the x-axis.
*		double ty: 		Translation along the y-axis.
*		double tz: 		Translation along the z-axis,
*					ignored for 2D transforms.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromTranslation(WlzTransformType type,
					double tx, double ty, double tz,
				        WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((newTr = WlzMakeAffineTransform(type, &errNum)) != NULL)
  {
    if((errNum = WlzAffineTransformTranslationSet(newTr,
    					     tx, ty, tz)) != WLZ_ERR_NONE)
    {
      WlzFreeAffineTransform(newTr);
      newTr = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newTr);
}

/************************************************************************
* Function:	WlzAffineTransformScaleSet
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Sets the given transform's matrix from the given
*		scales.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D or 3D affine transform.
*		double sx: 		Scale along the x-axis.
*		double sy: 		Scale along the y-axis.
*		double sz: 		Scale along the z-axis,
*					ignored for 2D transforms.
************************************************************************/
WlzErrorNum	WlzAffineTransformScaleSet(WlzAffineTransform *tr,
					   double sx, double sy, double sz)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(WlzAffineTransformDimension(tr, NULL))
    {
      case 2:
	tr->mat[0][0] = sx;
	tr->mat[1][1] = sy;
	tr->mat[0][1] = tr->mat[0][2] = 0.0;
	tr->mat[1][0] = tr->mat[1][2] = 0.0;
	tr->mat[2][0] = tr->mat[2][1] = 0.0;
	tr->mat[2][2] = 1.0;
        break;
      case 3:
	tr->mat[0][0] = sx;
	tr->mat[1][1] = sy;
	tr->mat[2][2] = sz;
	tr->mat[0][1] = tr->mat[0][2] = tr->mat[0][3] = 0.0;
	tr->mat[1][0] = tr->mat[1][2] = tr->mat[1][3] = 0.0;
	tr->mat[2][0] = tr->mat[2][1] = tr->mat[2][3] = 0.0;
	tr->mat[2][0] = tr->mat[2][1] = tr->mat[2][2] = 0.0;

        break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  return(errNum);
}


/************************************************************************
* Function:	WlzAffineTransformFromScale
* Returns:	WlzAffineTransform *:	New affine transform,
*					NULL on error.
* Purpose:	Constructs a new affine transform from the given 
*		scales.
* Global refs:	-
* Parameters:	WlzTransformType type:	Required transform type.
* 		double sx: 		Scale along the x-axis.
*		double sy: 		Scale along the y-axis.
*		double sz: 		Scale along the z-axis,
*					ignored for 2D transforms.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromScale(WlzTransformType type,
					double sx, double sy, double sz,
				        WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((newTr = WlzMakeAffineTransform(type, &errNum)) != NULL)
  {
    if((errNum = WlzAffineTransformScaleSet(newTr,
					    sx, sy, sz)) != WLZ_ERR_NONE)
    {
      WlzFreeAffineTransform(newTr);
      newTr = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newTr);
}

/************************************************************************
* Function:	WlzAffineTransformRotationSet
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Sets the given transform's matrix from the given
*		rotations. Although the 3 rotations contain redundant
*		information this may be a useful method for setting
*		rotation transforms. The order of composition is
*		R = Rz.Ry.Rx, x' = R.x.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D or 3D affine transform.
*		double rx: 		Rotation about the x-axis,
*					ignored for 2D transforms.
*		double ry: 		Rotation about the y-axis,
*					ignored for 2D transforms.
*		double rz: 		Rotation about the z-axis.
************************************************************************/
WlzErrorNum	WlzAffineTransformRotationSet(WlzAffineTransform *tr,
					    double rx, double ry, double rz)
{
  double	cx,
		sx,
		cy,
		sy,
		cz,
  		sz;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(WlzAffineTransformDimension(tr, NULL))
    {
      case 2:
	sz = sin(rz); cz = cos(rz);
	tr->mat[0][0] = cz; tr->mat[0][1] = -sz; tr->mat[0][2] = 0.0;
	tr->mat[1][0] = sz; tr->mat[1][1] = cz; tr->mat[1][2] = 0.0;
	tr->mat[2][0] = tr->mat[2][1] = 0.0; tr->mat[2][2] = 1.0;
        break;
      case 3:
	sx = sin(rx); cx = cos(rx);
	sy = sin(ry); cy = cos(ry);
	sz = sin(rz); cz = cos(rz);
	tr->mat[0][0] = cy * cz;
	tr->mat[0][1] = (sx * sy * cz) - (cx * sz);
	tr->mat[0][2] = (cx * sy * cz) + (sx * sz);
	tr->mat[1][0] = cy * sz;
	tr->mat[1][1] = (sx * sy * sz) + (cx * cz);
	tr->mat[1][2] = (cx * sy * sz) - (sx * cz);
	tr->mat[2][0] = -sy;
	tr->mat[2][1] = sx * cy;
	tr->mat[2][2] = cx * cy;
	tr->mat[0][3] = tr->mat[1][3] = tr->mat[2][3] = 0.0;
	tr->mat[3][0] = tr->mat[3][1] = tr->mat[3][2] = 0.0;
	tr->mat[3][3] = 1.0;
        break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformFromRotation
* Returns:	WlzAffineTransform *:	New affine transform,
*					NULL on error.
* Purpose:	Constructs a new affine transform from the given 
*		rotations.
* Global refs:	-
* Parameters:	WlzTransformType type:	Required transform type.
*		double rx: 		Rotation about the x-axis,
*					ignored for 2D transforms.
*		double ry: 		Rotation about the y-axis,
*					ignored for 2D transforms.
*		double rz: 		Rotation about the z-axis.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromRotation(WlzTransformType type,
					double rx, double ry, double rz,
				        WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((newTr = WlzMakeAffineTransform(type, &errNum)) != NULL)
  {
    if((errNum = WlzAffineTransformRotationSet(newTr,
    					     rx, ry, rz)) != WLZ_ERR_NONE)
    {
      WlzFreeAffineTransform(newTr);
      newTr = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newTr);
}

/************************************************************************
* Function:	WlzAffineTransformPrimGet
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Gets the given 2D transform's primitives from it's
*		matrix.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D affine transform.
*		WlzAffineTransformPrim *prim: Primitives data
*					structure to be set.
************************************************************************/
WlzErrorNum	WlzAffineTransformPrimGet(WlzAffineTransform *tr,
					   WlzAffineTransformPrim *prim)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (prim == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(WlzAffineTransformDimension(tr, NULL))
    {
      case 2:
        WlzAffineTransformPrimGet2(tr, prim);
	break;
      case 3:
	/* TODO don't know how to get the primitives from a 3D affine
	 * transform. */
        errNum = WLZ_ERR_TRANSFORM_TYPE;
        break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformPrimGet2
* Returns:	void
* Purpose:	Gets the given 2D transform's primitives from it's
*		matrix.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given 2D affine transform.
*		WlzAffineTransformPrim *prim: Primitives data
*					structure to be set.
************************************************************************/
static void	WlzAffineTransformPrimGet2(WlzAffineTransform *tr,
					   WlzAffineTransformPrim *prim)
{
  double  	s2,
  		tD0,
		tD1,
		tD2;

  /* Test for inversion */
  s2 = (tr->mat[0][0] * tr->mat[1][1]) - (tr->mat[0][1] * tr->mat[1][0]);
  if(s2 < 0.0)
  {
    prim->invert = 1;
    tr->mat[0][0] *= -1.0;
    tr->mat[0][1] *= -1.0;
    tr->mat[0][2] *= -1.0;
    s2 *= -1.0;
  }
  else
  {
    prim->invert = 0;
  }
  tD0 = tr->mat[0][0] + tr->mat[1][1];
  tD1 = tr->mat[0][1] - tr->mat[1][0];
  /* Scale and shear strength */ 
  if(fabs(s2) > DBL_EPSILON)
  {
    prim->scale = sqrt(s2);
    tD2 = (((tD0 * tD0) + (tD1 * tD1)) / s2)  - 4.0;
    prim->alpha = (tD2 > DBL_EPSILON)? sqrt(tD2): 0.0;
  }
  else
  {
    prim->scale =  0.0;
    prim->alpha = 0.0;
  }
  /* Rotation */
  prim->theta = atan2((prim->alpha * tD0 / 2.0) - tD1,
		      (prim->alpha * tD1 / 2.0) + tD0);
  /* Shear angle */
  prim->phi = 0.0;
  prim->xsi = 0.0;
  if(fabs(prim->alpha) > DBL_EPSILON)
  {
    tD2 = tan(prim->theta);
    prim->psi = -atan2(tr->mat[0][0] - tr->mat[1][1] +
		       ((tr->mat[0][1] + tr->mat[1][0]) * tD2),
		       (tr->mat[0][1] + (tr->mat[1][1] * tD2)) * 2.0);
  }
  else
  {
    prim->psi = 0.0;
  }
  /* Translation */
  prim->tx = tr->mat[0][2];
  prim->ty = tr->mat[1][2];
  prim->tz = 0.0;
  /* Restore matrix if inversion */
  if(prim->invert)
  {
    tr->mat[0][0] *= -1.0;
    tr->mat[0][1] *= -1.0;
    tr->mat[0][2] *= -1.0;
  }
}

/************************************************************************
* Function:	WlzAffineTransformMatrixSet
* Returns:	WlzErrorNum:		Error number.
* Purpose:	Sets the given transform from the given matrix.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Given affine transform.
*		double **matrix:	4x4 transform matrix values to
*					be copied.
************************************************************************/
WlzErrorNum	WlzAffineTransformMatrixSet(WlzAffineTransform *trans,
					    double **matrix)
{
  int		idx0,
  		idx1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformMatrixSet FE 0x%lx 0x%lx\n",
	   (unsigned long )trans, (unsigned long )matrix));
  if(trans == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(matrix == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch( trans->type )
    {
      case WLZ_TRANSFORM_2D_AFFINE:
	for(idx0 = 0; idx0 < 3; ++idx0)
	{
	  for(idx1 = 0; idx1 < 3; ++idx1)
	  {
	    trans->mat[idx0][idx1] = matrix[idx0][idx1];
	  }
	}
	break;
      case WLZ_TRANSFORM_3D_AFFINE:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformMatrixSet FX %d\n",
	   (int )errNum));
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformPrimValSet
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Sets a 2D affine transform from the given primitives.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr:	Given 2D affine transform.
*		double trX:		Column (x) translation.
*		double trY:		Line (y) translation.
*		double trZ:		Plane (z) translation.
*		double trScale:		Scale transformation.
*		double trTheta:		Rotation about z-axis.
*		double trPhi:		Rotation about y-axis.
*		double trAlpha:		Shear strength.
*		double trPsi:		Shear angle in x-y plane.
*		double trXsi:		3D shear angle.
*		int trInvert:		Reflection about y-axis if
*					non-zero.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzErrorNum	WlzAffineTransformPrimValSet(WlzAffineTransform *tr,
					     double trX,
					     double trY,
					     double trZ,
					     double trScale,
					     double trTheta,
					     double trPhi,
					     double trAlpha,
					     double trPsi,
					     double trXsi,
					     int trInvert)
{
  WlzAffineTransformPrim prim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformPrimValSet2D FE "
	  "0x%lx %g %g %g %g %g %g %g %g %g %d\n",
	   (unsigned long )tr, trX, trY, trZ, trScale, trTheta, trPhi,
	   trAlpha, trPsi, trXsi, trInvert));
  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    prim.tx = trX;
    prim.ty = trY;
    prim.tz = trZ;
    prim.scale = trScale;
    prim.theta = trTheta;
    prim.phi = trPhi;
    prim.alpha = trAlpha;
    prim.psi = trPsi;
    prim.xsi = trXsi;
    prim.invert = trInvert;
    errNum = WlzAffineTransformPrimSet(tr, prim);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformPrimValSet2D FX %d\n",
	   (int )errNum));
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformFromMatrix
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Makes a new affine transform of the given type and
*		then sets it's matrix.
* Global refs:	-
* Parameters:	WlzTransformType type:	Required transform type.
*		double **matrix:	Given matrix.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromMatrix(WlzTransformType type,
						 double **matrix,
						 WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromMatrix FE %d 0x%lx 0x%lx\n",
	   (int )type, (unsigned long )matrix,
	   (unsigned long )dstErr));
  if(matrix == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if((newTr = WlzMakeAffineTransform(type, &errNum)) != NULL)
    {
      if((errNum = WlzAffineTransformMatrixSet(newTr, matrix)) != WLZ_ERR_NONE)
      {
	WlzFreeAffineTransform(newTr);
	newTr = NULL;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromMatrix FX 0x%lx\n",
	   (unsigned long )newTr));
  return(newTr);
}

/************************************************************************
* Function:	WlzAffineTransformFromPrimVal
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Makes a new affine transform from the given primitive
*		transform properties.
* Global refs:	-
* Parameters:	WlzTransformType type:	Required transform type.
*		double trX:		Column (x) translation.
*		double trY:		Line (y) translation.
*		double trZ:		Plane (z) translation.
*		double trScale:		Scale transformation.
*		double trTheta:		Rotation about z-axis.
*		double trPhi:		Rotation about y-axis.
*		double trAlpha:		Shear strength.
*		double trPsi:		Shear angle in x-y plane.
*		double trXsi:		3D shear angle.
*		int trInvert:		Reflection about y-axis if
*					non-zero.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromPrimVal(WlzTransformType type,
				    	double trX, double trY, double trZ,
				    	double trScale, double trTheta,
					double trPhi, double trAlpha,
					double trPsi, double trXsi,
					int trInvert, WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromPrimVal FE "
	  "%d %g %g %g %g %g %g %g %g %g %d 0x%lx\n",
	   (int )type, trX, trY, trZ, trScale, trTheta, trPhi,
	   trAlpha, trPsi, trXsi, trInvert, (unsigned long )dstErr));
  if((newTr = WlzMakeAffineTransform(type, &errNum)) != NULL)
  {
    errNum = WlzAffineTransformPrimValSet(newTr, trX, trY, trZ,
					 trScale, trTheta, trPhi,
					 trAlpha, trPsi, trXsi,
					 trInvert);
    if(errNum != WLZ_ERR_NONE)
    {
      WlzFreeAffineTransform(newTr);
      newTr = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromPrimVal FX 0x%lx\n",
	   (unsigned long )newTr));
  return(newTr);
}
/************************************************************************
* Function:	WlzAffineTransformFromSpin
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Makes a new 2D affine transform from the given spin
*		angle and centre of rotation.
* Global refs:	-
* Parameters:	double spX:		Spin centre column (x).
*		double spY:		Spin centre line (y).
*		double spTheta:		Spin rotation about centre.
*					number.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromSpin(double spX, double spY,
					       double spTheta,
					       WlzErrorNum *dstErr)
{
  double	sinTheta,
  		cosTheta,
		trX,
		trY;
  WlzAffineTransform *newTrans;

  WLZ_DBG((WLZ_DBG_LVL_1),
          ("WlzAffineTransformFromSpin %g %g %g 0x%lx\n",
	   spX, spY, spTheta, dstErr));
  sinTheta = sin(spTheta);
  cosTheta = cos(spTheta);
  trX = spX - (spX * cosTheta) + (spY * sinTheta);
  trY = spY - (spX * sinTheta) - (spY * cosTheta);
  newTrans = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					   trX, trY, 0.0,
					   1.0, spTheta, 0.0, 0.0,
					   0.0, 0.0, 0, dstErr);
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromSpin FX 0x%lx\n",
	   (unsigned long )newTrans));
  return(newTrans);
}

/************************************************************************
* Function:	WlzAffineTransformFromSpinSqueeze
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Makes a new 2D affine transform from the given spin
*		angle, centre of rotation and scale factors.
* Global refs:	-
* Parameters:	double spX:		Spin centre column (x).
*		double spY:		Spin centre line (y).
*		double spTheta:		Spin rotation about centre.
*					number.
*		double sqX:		Squeeze (x) factor.
*		double sqY:		Squeeze (y) factor.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromSpinSqueeze(double spX, double spY,
					       double spTheta,
					       double sqX, double sqY,
					       WlzErrorNum *dstErr)
{
  double	sinTheta,
  		cosTheta;
  double	**matrix;
  WlzAffineTransform *newTr;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
          ("WlzAffineTransformFromSpinSqueeze %g %g %g %g %g 0x%lx\n",
	   spX, spY, spTheta, sqX, sqY, dstErr));
  if((newTr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE,
  				     &errNum)) == NULL)
  {
    errNum = WLZ_ERR_UNSPECIFIED;
  }
  else {
    matrix = newTr->mat;
    sinTheta = sin(spTheta);
    cosTheta = cos(spTheta);
    matrix[0][0] =  sqX * cosTheta;
    matrix[0][1] = -sqX * sinTheta;
    matrix[1][0] =  sqY * sinTheta;
    matrix[1][1] =  sqY * cosTheta;
    matrix[0][2] =  (spX * (1 - matrix[0][0])) - (spY * matrix[0][1]);
    matrix[1][2] = (-spX * matrix[1][0]) + (spY * (1.0 - matrix[1][1]));
    matrix[2][0] = 0.0;
    matrix[2][1] = 0.0;
    matrix[2][2] = 1.0;
  }

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromSpinSqueeze FX 0x%lx\n",
	   (unsigned long )newTr));
  if( dstErr ){
    *dstErr = errNum;
  }
  return(newTr);
}

/************************************************************************
* Function:	WlzAffineTransformCopy
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Copies the given affine transform.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Given affine transform.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformCopy(WlzAffineTransform *tr,
					   WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformCopy FE 0x%lx 0x%lx\n",
	   (unsigned long )tr, (unsigned long )dstErr));
  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    newTr = WlzAffineTransformFromMatrix(tr->type, tr->mat, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformCopy FX 0x%lx\n",
	   (unsigned long )newTr));
  return(newTr);
}

/************************************************************************
* Function:	WlzAffineTransformProduct
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Computes the product of the two given affine
*		transforms.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr0: First affine transform.
*		WlzAffineTransform *tr1: Second affine transform.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform *WlzAffineTransformProduct(WlzAffineTransform *tr0,
					      WlzAffineTransform *tr1,
					      WlzErrorNum *dstErr)
{
  int		idx0,
  		idx1,
		idx2;
  double	tD0;
  WlzAffineTransform *prodTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformProduct FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )tr0, (unsigned long )tr1,
	   (unsigned long )dstErr));
  if((tr0 == NULL) || (tr1 == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(tr0->type != tr1->type)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    prodTr = WlzMakeAffineTransform(tr0->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(tr0->type)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
	for(idx0=0; idx0 < 3; idx0++)
	{
	  for(idx1=0; idx1 < 3; idx1++)
	  {
	    tD0 = 0.0;
	    for(idx2=0; idx2 <3; idx2++)
	    {
	      tD0 += tr1->mat[idx0][idx2] * tr0->mat[idx2][idx1];
	    }
	    prodTr->mat[idx0][idx1] = tD0;
	  }
	}
	break;
      case WLZ_TRANSFORM_3D_AFFINE:
	for(idx0=0; idx0 < 4; idx0++)
	{
	  for(idx1=0; idx1 < 4; idx1++)
	  {
	    tD0 = 0.0;
	    for(idx2=0; idx2 <4; idx2++)
	    {
	      tD0 += tr1->mat[idx0][idx2] * tr0->mat[idx2][idx1];
	    }
	    prodTr->mat[idx0][idx1] = tD0;
	  }
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
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformProduct FX 0x%lx\n",
	   (unsigned long )prodTr));
  return(prodTr);
}

/************************************************************************
* Function:	WlzAffineTransformInverse
* Returns:	WlzAffineTransform *:	New affine transform, or NULL
*					on error.
* Purpose:	Computes the inverse of the given affine transform.
* Global refs:	-
* Parameters:	WlzAffineTransform *tr: Given affine transform.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number.
************************************************************************/
WlzAffineTransform	*WlzAffineTransformInverse(WlzAffineTransform *tr,
					           WlzErrorNum *dstErr)
{
  int		dim;
  WlzAffineTransform *invTr = NULL;
  double	*matrix[4];
  double	buffer[16];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformInverse FE 0x%lx 0x%lx\n",
	   (unsigned long )tr, (unsigned long )dstErr));
  if(tr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    invTr = WlzAffineTransformCopy(tr, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dim = WlzAffineTransformDimension(tr, NULL);
    if(AlgMatrixLUInvert(invTr->mat, dim + 1) != ALG_ERR_NONE)
    {
      errNum = WLZ_ERR_TRANSFORM_DATA;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformInverse FX 0x%lx\n",
	   (unsigned long )invTr));
  return(invTr);
}

/************************************************************************
* Function:     WlzAffineTransformIsIdentity
* Returns:      int:                    Non-zero if the given transform
*                                       is an identity transform.
* Purpose:      Checks whether the given transform is an identity
*               transform.
* Global refs:  -
* Parameters:   WlzAffineTransform *trans: Given affine transform.
*               WlzErrorNum *dstErr:    Destination pointer for error
*                                       number.
************************************************************************/
int		WlzAffineTransformIsIdentity(WlzAffineTransform *trans,
					     WlzErrorNum *dstErr)
{
  int           dim,
  		isIdentity = 0;
  double        tD0;
  double        **mat;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
 
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformIsIdentity FE 0x%lx 0x%lx\n",
	   (unsigned long )trans, (unsigned long )dstErr));
  if(trans == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    mat = trans->mat;
    dim = WlzAffineTransformDimension(trans, NULL);
    switch(dim)
    {
      case 2:
	if((fabs(mat[0][0] - 1.0) <= DBL_EPSILON) &&
	   (fabs(mat[0][1]) <= DBL_EPSILON) &&
	   (fabs(mat[0][2]) <= DBL_EPSILON) &&
	   (fabs(mat[1][0]) <= DBL_EPSILON) &&
	   (fabs(mat[1][1] - 1.0) <= DBL_EPSILON) &&
	   (fabs(mat[1][2]) <= DBL_EPSILON) &&
	   (fabs(mat[2][0]) <= DBL_EPSILON) &&
	   (fabs(mat[2][1]) <= DBL_EPSILON) &&
	   (fabs(mat[2][2] - 1.0) <= DBL_EPSILON)) 
	{
	  isIdentity = 1;
	}
        break;
      case 3:
	if((fabs(mat[0][0] - 1.0) <= DBL_EPSILON) &&
	   (fabs(mat[0][1]) <= DBL_EPSILON) &&
	   (fabs(mat[0][2]) <= DBL_EPSILON) &&
	   (fabs(mat[0][3]) <= DBL_EPSILON) &&
	   (fabs(mat[1][0]) <= DBL_EPSILON) &&
	   (fabs(mat[1][1] - 1.0) <= DBL_EPSILON) &&
	   (fabs(mat[1][2]) <= DBL_EPSILON) &&
	   (fabs(mat[1][3]) <= DBL_EPSILON) &&
	   (fabs(mat[2][0]) <= DBL_EPSILON) &&
	   (fabs(mat[2][1]) <= DBL_EPSILON) &&
	   (fabs(mat[2][2] - 1.0) <= DBL_EPSILON) &&
	   (fabs(mat[2][3]) <= DBL_EPSILON) &&
	   (fabs(mat[3][0]) <= DBL_EPSILON) &&
	   (fabs(mat[3][1]) <= DBL_EPSILON) &&
	   (fabs(mat[3][2]) <= DBL_EPSILON) &&
	   (fabs(mat[3][3] - 1.0) <= DBL_EPSILON))
	{
	  isIdentity = 1;
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
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformIsIdentity FX %d\n",
	   isIdentity));
  return(isIdentity);
}

/************************************************************************
* Function:	WlzAffineTransformObj
* Returns:	WlzObject *:		Transformed object, NULL on
*					error.
* Purpose:	Applies the given affine transform to the given Woolz
*		object.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Object to be transformed.
*		WlzAffineTransform *trans: Affine transform to apply.
*		WlzInterpolationType interp: Level of interpolation to
*					use.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
WlzObject	*WlzAffineTransformObj(WlzObject *srcObj,
				       WlzAffineTransform *trans,
				       WlzInterpolationType interp,
				       WlzErrorNum *dstErr)
{
  int		planeCount,
  		planeIdx;
  WlzDomain	srcDom,
  		dstDom;
  WlzValues	srcValues,
  		dumValues,
  		dstValues;
  WlzObject	*tObj0,
  		*tObj1,
		*dstObj = NULL;
  WlzAffineTransformPrim prim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformObj FE 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )srcObj, (unsigned long )trans,
	   (int )interp, (unsigned long )dstErr));
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    dstDom.core = NULL;
    dumValues.core = NULL;
    dstValues.core = NULL;
    srcValues.core = NULL;
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	if((dstObj = WlzMakeEmpty(&errNum)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
        break;
      case WLZ_2D_POLYGON:
      case WLZ_BOUNDLIST:
      case WLZ_TRANS_OBJ:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(srcObj->type)
	  {
	    case WLZ_2D_POLYGON:
	      dstDom.poly = WlzAffineTransformPoly2(srcObj->domain.poly,
						     trans, &errNum);
	      break;
	    case WLZ_BOUNDLIST:
	      dstDom.b = WlzAffineTransformBoundList(srcObj->domain.b,
						     trans, &errNum);
	      break;
	    case WLZ_TRANS_OBJ:
	      dstDom.t = WlzAffineTransformProduct(srcObj->domain.t,
	      					   trans, &errNum);
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if((dstObj = WlzMakeMain(srcObj->type, dstDom,
				   srcValues, NULL, NULL, &errNum)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    (void )WlzFreeDomain(dstDom);
	  }
	}
        break;
      case WLZ_2D_DOMAINOBJ:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	} 
	else if(WlzAffineTransformIsTranslate(trans, srcObj, NULL))
	{
	  dstObj = WlzAffineTransformIntTranslate(srcObj, trans, &errNum);
	}
	else
	{
	  tObj0 = NULL;
	  tObj1 = NULL;
	  tObj0 = WlzObjToBoundary(srcObj, 1, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tObj1 = WlzAffineTransformObj(tObj0, trans, interp, &errNum);
	    WlzFreeObj(tObj0);
	    tObj0 = NULL;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj = WlzBoundToObj(tObj1->domain.b, WLZ_EVEN_ODD_FILL,
				   &errNum);
	    WlzFreeObj(tObj1);
	    tObj1 = NULL;
	  }
	  if((errNum == WLZ_ERR_NONE) &&
	     (srcObj->values.core) )
	  {
	    errNum = WlzAffineTransformValues2(dstObj, srcObj, trans,
						interp);
	  }
	  if(tObj0)
	  {
	    WlzFreeObj(tObj0);
	  }
	  if(tObj1)
	  {
	    WlzFreeObj(tObj1);
	  }
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	srcDom = srcObj->domain;
	srcValues = srcObj->values;
	if(srcDom.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	} 
	else if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if(WlzAffineTransformIsTranslate(trans, srcObj, NULL))
	{
	  dstObj = WlzAffineTransformIntTranslate(srcObj, trans, &errNum);
	}
	else
	{
	  /* TODO. This is a kludge. Currently a full 3D transform has not been
	   * implemented and in this code the transform is applied to each
	   * plane in turn. If the transform has a z translation set then
	   * the resultant object is shifted in z - BUT NOT SCALED!!!! 
	   * In principle the z-transformation associated with 2D transform 
	   * is meaningless and should have been flagged as an error. Clearly
	   * it does not contribute to the tranformation matrix. */
	  dstDom.p = WlzMakePlaneDomain(srcDom.p->type,
	  				srcDom.p->plane1, srcDom.p->lastpl,
					srcDom.p->line1, srcDom.p->lastln,
					srcDom.p->kol1, srcDom.p->lastkl,
					&errNum);
	  			 /* Need to fix the line column bounds later */
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzAffineTransformPrimGet(trans, &prim);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstDom.p->voxel_size[0] = srcDom.p->voxel_size[0];
	    dstDom.p->voxel_size[1] = srcDom.p->voxel_size[1];
	    dstDom.p->voxel_size[2] = srcDom.p->voxel_size[2];
	    if(srcValues.core)
	    {
	      dstValues.vox = WlzMakeVoxelValueTb(srcObj->values.vox->type,
						  srcDom.p->plane1,
						  srcDom.p->lastpl,
						  WlzGetBackground(srcObj,
								   NULL),
						  NULL, &errNum);
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    planeIdx = 0;
	    planeCount = srcDom.p->lastpl - srcDom.p->plane1 + 1;
	    while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	    {
	      tObj0 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				  *(srcDom.p->domains + planeIdx),
				  dumValues, NULL, NULL, &errNum);
	      if(tObj0->domain.core &&
		 (tObj0->domain.core->type != WLZ_EMPTY_OBJ))
	      {
		if(srcValues.core &&
		   (srcValues.core->type != WLZ_EMPTY_OBJ))
		{
		  tObj0->values = WlzAssignValues(*(srcValues.vox->values +
						    planeIdx), &errNum);
		}
		else
		{
		  tObj0->values.core = NULL;
		}
		tObj1 = WlzAffineTransformObj(tObj0, trans, interp, &errNum);
	      }
	      else
	      {
		tObj1 = NULL;
	      }
	      if( tObj0 ){
		(void) WlzFreeObj(tObj0);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(tObj1)
		{
		  *(dstDom.p->domains + planeIdx) = tObj1->domain;
		  if( dstValues.vox ){
		    *(dstValues.vox->values + planeIdx) = tObj1->values;
		  }
		  tObj1->domain.core = NULL;
		  tObj1->values.core = NULL;
		  (void )WlzFreeObj(tObj1);
		  tObj1 = NULL;
		}
		else
		{
		  (dstDom.p->domains + planeIdx)->core = NULL;
		  if( dstValues.vox ){
		    (dstValues.vox->values + planeIdx)->core = NULL;
		  }
		}
	      }
	      ++planeIdx;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzStandardPlaneDomain(dstDom.p, dstValues.vox);
	    if( dstDom.p ){
	      dstDom.p->plane1 += prim.tz;
	      dstDom.p->lastpl += prim.tz;
	    }
	    if( dstValues.vox ){
	      dstValues.vox->plane1 += prim.tz;
	      dstValues.vox->lastpl += prim.tz;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues,
				 NULL, NULL, &errNum);
	  }
	}
        break;
      default:
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformObj FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/************************************************************************
* Function:	WlzAffineTransformVertexI2
* Returns:	WlzIVertex2:		Transformed vertex.
* Purpose:	Transforms the given WlzIVertex2.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Affine transform to apply.
*		WlzIVertex2 srcVtx:	Vertex to be transformed.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
WlzIVertex2	WlzAffineTransformVertexI2(WlzAffineTransform *trans,
					   WlzIVertex2 srcVtx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex2	dVtx;
  WlzIVertex2	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dVtx.vtX = srcVtx.vtX;
  dVtx.vtY = srcVtx.vtY;
  dVtx = WlzAffineTransformVertexD2(trans, dVtx, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dstVtx.vtX = WLZ_NINT(dVtx.vtX);
    dstVtx.vtY = WLZ_NINT(dVtx.vtY);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVtx);
}

/************************************************************************
* Function:	WlzAffineTransformVertexI3
* Returns:	WlzIVertex3:		Transformed vertex.
* Purpose:	Transforms the given WlzIVertex3.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Affine transform to apply.
*		WlzIVertex3 srcVtx:	Vertex to be transformed.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
WlzIVertex3	WlzAffineTransformVertexI3(WlzAffineTransform *trans,
					   WlzIVertex3 srcVtx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex3	dVtx;
  WlzIVertex3	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dVtx.vtX = srcVtx.vtX;
  dVtx.vtY = srcVtx.vtY;
  dVtx.vtZ = srcVtx.vtZ;
  dVtx = WlzAffineTransformVertexD3(trans, dVtx, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dstVtx.vtX = WLZ_NINT(dVtx.vtX);
    dstVtx.vtY = WLZ_NINT(dVtx.vtY);
    dstVtx.vtZ = WLZ_NINT(dVtx.vtZ);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVtx);
}

/************************************************************************
* Function:	WlzAffineTransformVertexF2
* Returns:	WlzFVertex2:		Transformed vertex.
* Purpose:	Transforms the given WlzFVertex2.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Affine transform to apply.
*		WlzFVertex2 srcVtx:	Vertex to be transformed.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
WlzFVertex2	WlzAffineTransformVertexF2(WlzAffineTransform *trans,
					   WlzFVertex2 srcVtx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex2	dVtx;
  WlzFVertex2	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dVtx.vtX = srcVtx.vtX;
  dVtx.vtY = srcVtx.vtY;
  dVtx = WlzAffineTransformVertexD2(trans, dVtx, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dstVtx.vtX = dVtx.vtX;
    dstVtx.vtY = dVtx.vtY;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVtx);
}

/************************************************************************
* Function:	WlzAffineTransformVertexF3
* Returns:	WlzFVertex3:		Transformed vertex.
* Purpose:	Transforms the given WlzFVertex3.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Affine transform to apply.
*		WlzFVertex3 srcVtx:	Vertex to be transformed.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
WlzFVertex3	WlzAffineTransformVertexF3(WlzAffineTransform *trans,
					   WlzFVertex3 srcVtx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex3	dVtx;
  WlzFVertex3	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dVtx.vtX = srcVtx.vtX;
  dVtx.vtY = srcVtx.vtY;
  dVtx.vtZ = srcVtx.vtZ;
  dVtx = WlzAffineTransformVertexD3(trans, dVtx, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dstVtx.vtX = dVtx.vtX;
    dstVtx.vtY = dVtx.vtY;
    dstVtx.vtZ = dVtx.vtZ;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVtx);
}

/************************************************************************
* Function:	WlzAffineTransformVertexD2
* Returns:	WlzDVertex2:		Transformed vertex.
* Purpose:	Transforms the given WlzDVertex2.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Affine transform to apply.
*		WlzDVertex2 srcVtx:	Vertex to be transformed.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
WlzDVertex2	WlzAffineTransformVertexD2(WlzAffineTransform *trans,
					   WlzDVertex2 srcVtx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex2	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(WlzAffineTransformDimension(trans, NULL) != 2)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    dstVtx.vtX = (trans->mat[0][0] * srcVtx.vtX) +
		 (trans->mat[0][1] * srcVtx.vtY) + trans->mat[0][2];
    dstVtx.vtY = (trans->mat[1][0] * srcVtx.vtX) +
		 (trans->mat[1][1] * srcVtx.vtY) + trans->mat[1][2];
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVtx);
}

/************************************************************************
* Function:	WlzAffineTransformVertexD3
* Returns:	WlzDVertex3:		Transformed vertex.
* Purpose:	Transforms the given WlzDVertex3.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Affine transform to apply.
*		WlzDVertex3 srcVtx:	Vertex to be transformed.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL.
************************************************************************/
WlzDVertex3	WlzAffineTransformVertexD3(WlzAffineTransform *trans,
					   WlzDVertex3 srcVtx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex3	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(WlzAffineTransformDimension(trans, NULL))
    {
      case 2:
        dstVtx.vtX = (trans->mat[0][0] * srcVtx.vtX) +
		     (trans->mat[0][1] * srcVtx.vtY) + trans->mat[0][2];
        dstVtx.vtY = (trans->mat[1][0] * srcVtx.vtX) +
		     (trans->mat[1][1] * srcVtx.vtY) + trans->mat[1][2];
        break;
      case 3:
	dstVtx.vtX = (trans->mat[0][0] * srcVtx.vtX) +
		     (trans->mat[0][1] * srcVtx.vtY) +
		     (trans->mat[0][2] * srcVtx.vtZ) + trans->mat[0][3];
	dstVtx.vtY = (trans->mat[1][0] * srcVtx.vtX) +
		     (trans->mat[1][1] * srcVtx.vtY) +
		     (trans->mat[1][2] * srcVtx.vtZ) + trans->mat[1][3];
	dstVtx.vtZ = (trans->mat[2][0] * srcVtx.vtX) +
		     (trans->mat[2][1] * srcVtx.vtY) +
		     (trans->mat[2][2] * srcVtx.vtZ) + trans->mat[2][3];
        break;
      default:
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVtx);
}
