#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzAffineTransform.c
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
* \brief        Functions for computing Woolz affine transforms and
*		applying them to Woolz objects.
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <limits.h>
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
static WlzPlaneDomain 		*WlzAffineTransformPDom(
				  WlzObject *srcObj,
				  WlzAffineTransform *trans,
				  WlzErrorNum *dstErr);
static WlzPolygonDomain 	*WlzAffineTransformPoly2(
				  WlzPolygonDomain *,
				  WlzAffineTransform *,
				  WlzErrorNum *);
static WlzBoundList 		*WlzAffineTransformBoundList(WlzBoundList *,
				  WlzAffineTransform *,
				  WlzErrorNum *);
static WlzErrorNum 		WlzAffineTransformValues2(WlzObject *newObj,
				  WlzObject *srcObj,
				  WlzAffineTransform *tr,
				  WlzInterpolationType interp);
static WlzErrorNum 		WlzAffineTransformValues3(WlzObject *newObj,
				  WlzObject *srcObj,
				  WlzAffineTransform *tr,
				  WlzInterpolationType interp);
static WlzErrorNum 		WlzAffineTransformPrimSet2(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim prim);
static WlzErrorNum 		WlzAffineTransformPrimSet3(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim prim);
static void			WlzAffineTransformPrimGet2(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim *prim);

/*!
* \ingroup	Wlz
* \defgroup	WlzAffineTransform
* @{
*/

/*!
* \return				2 or 3 for a 2D or 3D affine
*					transform, 0 on error.
* \brief	Computes the dimension of the given affine transform.
* \param 	tr			Given affine transform.
* \param	dstErr			Destination error pointer, may
*					be null.
*/
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

/*!
* \return				Non-zero if translation.
* \brief	Tests whether the given affine transform is a simple
*		integer translation.
* \param	tr			Given affine transform.
* \param	obj			Optional object, may be NULL.
* \param	dstErr			Destination error pointer, may
*					be null.
*/
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

/*!
* \return				Non-zero if translation.
* \brief	Tests wether the given 2D affine transform is a simple
*		integer translation.
*		Because this is a static function the parameters are
*		not checked.
* \param	tr			Given 2D affine transform.
* \param	obj			Optional object, may be NULL.
* \param	dstErr			Destination error pointer, may
*					be null.
*/
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

/*!
* \return				Non-zero if translation.
* \brief	Tests wether the given 3D affine transform is a simple
*		integer translation.
*		Because this is a static function the parameters are
*		not checked.
* \param	tr			Given 3D affine transform.
* \param	obj			Optional object, may be NULL.
* \param	dstErr			Destination error pointer, may
*					be null.
*/
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

/*!
* \return				Translated object or NULL on
*					error.
* \brief	Translates the given 2D domain object with an integral
*		translation.
*		Because this is a static function the parameters are
*		not checked.
* \param	srcObj			Given 2D domain object.
* \param	tr			Given affine transform.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return			     	Transformed polygon domain or
*                                    	NULL on error.
* \brief      	Transforms the given polygon domain.
*             	Because this is a static function the parameters (other
*             	than polygon type) are not checked.
* \param	srcPoly		 	Given polygon domain.
* \param        trans		 	Given affine transform.
* \param        dstErr		 	Destination pointer for error number.
*/
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

/*!
* \return				Transformed boundary list or
*					NULL on error.
* \brief	Transforms the given boundary list.
*		Because this is a static function the parameters are
*		not checked.
* \param	srcBound		Given boundary list.
* \param	trans			Given affine transform.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				Transformed contour or
*					NULL on error.
* \brief	Transforms the given contour.
* \param	srcCtr		 	Given contour.
* \param	tr			Given affine transform.
* \param	newModFlg		Make a new model if non-zero,
*					otherwise transform the given
*					model in place.
* \param	dstErr			Destination pointer for error
*					number.
*/
WlzContour	*WlzAffineTransformContour(WlzContour *srcCtr,
					   WlzAffineTransform *tr,
					   int newModFlg,
					   WlzErrorNum *dstErr)
{
  WlzContour 	*dstCtr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstCtr = WlzMakeContour(&errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dstCtr->model = WlzAssignGMModel(
    		    WlzAffineTransformGMModel(srcCtr->model, tr, newModFlg,
		    			      &errNum), NULL);
  }
  if((errNum != WLZ_ERR_NONE) && dstCtr)
  {
    (void )WlzFreeContour(dstCtr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstCtr);
}

/*!
* \return				Transformed model or
*					NULL on error.
* \brief	Transforms the given geometric model.
* \param	srcM		 	Given geometric model.
* \param	tr			Given affine transform.
* \param	newModFlg		Make a new model if non-zero,
*					otherwise transform the given
*					model in place.
* \param	dstErr			Destination pointer for error
*					number.
*/
WlzGMModel	*WlzAffineTransformGMModel(WlzGMModel *srcM,
					   WlzAffineTransform *tr,
					   int newModFlg,
					   WlzErrorNum *dstErr)
{
  int		idx,
  		cnt;
  AlcVector	*vec;
  WlzGMElemP	elmP;
  WlzGMModel 	*dstM = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcM == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  dstM = (newModFlg)? WlzGMModelCopy(srcM, &errNum):
  		      WlzAssignGMModel(srcM, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Transform vertex geometries. */
    idx = 0;
    vec = dstM->res.vertexG.vec;
    cnt = dstM->res.vertexG.numIdx;
    while((idx < cnt) && (errNum == WLZ_ERR_NONE))
    {
      elmP.core = (WlzGMCore *)AlcVectorItemGet(vec, idx);
      if(elmP.core && (elmP.core->idx >= 0))
      {
        switch(dstM->type)
	{
	  case WLZ_GMMOD_2I:
	    elmP.vertexG2I->vtx = WlzAffineTransformVertexI2(tr,
	    				elmP.vertexG2I->vtx, &errNum);
	    break;
	  case WLZ_GMMOD_2D:
	    elmP.vertexG2D->vtx = WlzAffineTransformVertexD2(tr,
	    				elmP.vertexG2D->vtx, &errNum);
	    break;
	  case WLZ_GMMOD_3I:
	    elmP.vertexG3I->vtx = WlzAffineTransformVertexI3(tr,
	    				elmP.vertexG3I->vtx, &errNum);
	    break;
	  case WLZ_GMMOD_3D:
	    elmP.vertexG3D->vtx = WlzAffineTransformVertexD3(tr,
	    				elmP.vertexG3D->vtx, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
      }
      ++idx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Transform shell geometries. */
    idx = 0;
    vec = dstM->res.shellG.vec;
    cnt = dstM->res.shellG.numIdx;
    while((idx < cnt) && (errNum == WLZ_ERR_NONE))
    {
      elmP.core = (WlzGMCore *)AlcVectorItemGet(vec, idx);
      if(elmP.core && (elmP.core->idx >= 0))
      {
        switch(dstM->type)
	{
	  case WLZ_GMMOD_2I:
	    elmP.shellG2I->bBox = WlzAffineTransformBBoxI2(tr,
	    				elmP.shellG2I->bBox, &errNum);
	    break;
	  case WLZ_GMMOD_2D:
	    elmP.shellG2D->bBox = WlzAffineTransformBBoxD2(tr,
	    				elmP.shellG2D->bBox, &errNum);
	    break;
	  case WLZ_GMMOD_3I:
	    elmP.shellG3I->bBox = WlzAffineTransformBBoxI3(tr,
	    				elmP.shellG3I->bBox, &errNum);
	    break;
	  case WLZ_GMMOD_3D:
	    elmP.shellG3D->bBox = WlzAffineTransformBBoxD3(tr,
	    				elmP.shellG3D->bBox, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
      }
      ++idx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMModelRehashVHT(dstM, 0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstM);
}

/*!
* \return				Error number.
* \brief	Creates a new value table, fills in the values and
*		adds it to the given new object.
*		Because this is a static function the parameters are
*		not checked.
* \param	newObj			Partialy transformed object
*					with a valid domain.
* \param	srcObj			2D domain object which is being
*					transformed.
* \param	trans			Given affine transform.
* \param	interp			Level of interpolation to
*					use.
*/
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

/*!
* \return				Error number.
* \brief	Creates new value, fills in the values and adds it
*		to the given new object.
*		not checked.
* \param	newObj			Partialy transformed object
*					with a valid domain.
* \param	srcObj			3D domain object which is being
*					transformed.
* \param	trans			Given affine transform.
* \param	interp			Level of interpolation to
*					use.
*/
static WlzErrorNum WlzAffineTransformValues3(WlzObject *newObj,
					     WlzObject *srcObj,
					     WlzAffineTransform *trans,
					     WlzInterpolationType interp)
{
  int		tI0,
  		count;
  double	tD0, x, y, z;
  WlzIVertex3	idx,
  		sPos,
		dPos;
  WlzDVertex3	tDV0,
  		tDV1;
  WlzIBox3	bBox;
  WlzPixelV	bkdV;
  WlzValues	tVal,
  		dstValues,
  		emptyValues;
  WlzObject 	*tObj0;
  WlzGreyWSpace	gWSp;
  WlzIntervalWSpace iWSp;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzAffineTransform *invTrans = NULL;
  WlzGreyType	gType;
  WlzDomain	dom2D;
  double	tMat[3][3];
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  dstValues.core = NULL;
  emptyValues.core = NULL;
  /* Make a new voxel value table. */
  bkdV = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    gType = WlzGreyTypeFromObj(srcObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox = WlzBoundingBox3D(newObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					bBox.zMin, bBox.zMax,
					bkdV, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newObj->values = WlzAssignValues(dstValues, NULL);
  }
  /* Invert the transform. */
  if(errNum == WLZ_ERR_NONE)
  {
    invTrans = WlzAffineTransformInverse(trans, &errNum);
  }
  /* For each plane in the new object make a new value table and
   * then fill it in. */
  if(errNum == WLZ_ERR_NONE)
  {
    gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idx.vtZ = 0;
    dPos.vtZ = bBox.zMin;
    tMat[0][0] = invTrans->mat[0][0];
    tMat[1][0] = invTrans->mat[1][0];
    tMat[2][0] = invTrans->mat[2][0];
    while((errNum == WLZ_ERR_NONE) && (dPos.vtZ <= bBox.zMax))
    {
      if(((dom2D = *(newObj->domain.p->domains + idx.vtZ)).core != NULL) &&
         (dom2D.core->type != WLZ_EMPTY_DOMAIN))
      {
	tMat[0][2] = invTrans->mat[0][3] + (invTrans->mat[0][2] * dPos.vtZ);
	tMat[1][2] = invTrans->mat[1][3] + (invTrans->mat[1][2] * dPos.vtZ);
	tMat[2][2] = invTrans->mat[2][3] + (invTrans->mat[2][2] * dPos.vtZ);
	/* Make a 2D domain object for the plane. */
	tObj0 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	    *(newObj->domain.p->domains + idx.vtZ),
	    emptyValues, NULL, NULL, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  tVal.v = WlzNewValueTb(tObj0,
	  		         WlzGreyTableType(WLZ_GREY_TAB_RAGR, gType,
				 		  NULL),
				 bkdV, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  tObj0->values = WlzAssignValues(tVal, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzInitGreyScan(tObj0, &iWSp, &gWSp);
	}
	/* Fill in the values of the new 2D object. */
	while((errNum == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
	{
	  dPos.vtX = iWSp.lftpos;
	  dPos.vtY = iWSp.linpos;
	  tMat[0][1] = tMat[0][2] + (invTrans->mat[0][1] * dPos.vtY);
	  tMat[1][1] = tMat[1][2] + (invTrans->mat[1][1] * dPos.vtY);
	  tMat[2][1] = tMat[2][2] + (invTrans->mat[2][1] * dPos.vtY);
	  count = iWSp.rgtpos - iWSp.lftpos + 1;
	  switch(interp)
	  {
	    case WLZ_INTERPOLATION_NEAREST:
	      while(count-- > 0)
	      {
		sPos.vtX = tMat[0][1] + (tMat[0][0] * dPos.vtX);
		sPos.vtY = tMat[1][1] + (tMat[1][0] * dPos.vtX);
		sPos.vtZ = tMat[2][1] + (tMat[2][0] * dPos.vtX);
		WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
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
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		++(dPos.vtX);
	      }
	      break;
	    case WLZ_INTERPOLATION_LINEAR:
	      while(count-- > 0)
	      {
		x = tMat[0][1] + (tMat[0][0] * dPos.vtX);
		y = tMat[1][1] + (tMat[1][0] * dPos.vtX);
		z = tMat[2][1] + (tMat[2][0] * dPos.vtX);
		WlzGreyValueGetCon(gVWSp, z, y, x);
		tDV0.vtX = x - WLZ_NINT(x - 0.5);
		tDV0.vtY = y - WLZ_NINT(y - 0.5);
		tDV0.vtZ = z - WLZ_NINT(z - 0.5);
		tDV1.vtX = 1.0 - tDV0.vtX;
		tDV1.vtY = 1.0 - tDV0.vtY;
		tDV1.vtZ = 1.0 - tDV0.vtZ;
		switch(gWSp.pixeltype)
		{
		  case WLZ_GREY_INT:
		    tD0 = ((gVWSp->gVal[0]).inv *
			tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).inv *
		       tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).inv *
		       tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).inv *
		       tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).inv *
		       tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).inv *
		       tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).inv *
		       tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).inv *
		       tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    tD0 = WLZ_CLAMP(tD0, INT_MIN, INT_MAX);
		    tI0 = WLZ_NINT(tD0);
		    *(gWSp.u_grintptr.inp)++ = tI0;
		    break;
		  case WLZ_GREY_SHORT:
		    tD0 = ((gVWSp->gVal[0]).shv *
			tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).shv *
		       tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).shv *
		       tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).shv *
		       tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).shv *
		       tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).shv *
		       tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).shv *
		       tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).shv *
		       tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    tD0 = WLZ_CLAMP(tD0, SHRT_MIN, SHRT_MAX);
		    tI0 = WLZ_NINT(tD0);
		    *(gWSp.u_grintptr.shp)++ = tI0;
		    break;
		  case WLZ_GREY_UBYTE:
		    tD0 = ((gVWSp->gVal[0]).ubv *
			tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).ubv *
		       tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).ubv *
		       tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).ubv *
		       tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).ubv *
		       tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).ubv *
		       tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).ubv *
		       tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).ubv *
		       tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    tD0 = WLZ_CLAMP(tD0, 0, 255);
		    tI0 = WLZ_NINT(tD0);
		    *(gWSp.u_grintptr.ubp)++ = tI0;
		    break;
		  case WLZ_GREY_FLOAT:
		    tD0 = ((gVWSp->gVal[0]).flv *
			tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).flv *
		       tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).flv *
		       tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).flv *
		       tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).flv *
		       tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).flv *
		       tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).flv *
		       tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).flv *
		       tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    *(gWSp.u_grintptr.flp)++ = tD0;
		    break;
		  case WLZ_GREY_DOUBLE:
		    tD0 = ((gVWSp->gVal[0]).dbv *
			tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).dbv *
		       tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).dbv *
		       tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).dbv *
		       tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).dbv *
		       tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).dbv *
		       tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).dbv *
		       tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).dbv *
		       tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    tD0 = WLZ_CLAMP(tD0, FLT_MIN, FLT_MAX);
		    *(gWSp.u_grintptr.dbp)++ = tD0;
		    break;
		  default:
		    errNum = WLZ_ERR_GREY_TYPE;
		    break;
		}
		++(dPos.vtX);
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_INTERPOLATION_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *(newObj->values.vox->values + idx.vtZ) =
	    WlzAssignValues(tObj0->values, NULL);
	}
	if(tObj0)
	{
	  (void )WlzFreeObj(tObj0);
	}
      }
      ++(idx.vtZ);
      ++(dPos.vtZ);
    }
  }
  if(invTrans)
  {
    (void )WlzFreeAffineTransform(invTrans);
  }
  return(errNum);
}

/*!
* \return				Transformed plane domain,
*					NULL on error.
* \brief	Creates a new plane domain which is the transformed
*		plane domain of the given source object.
*		The algorithm used by this function isn't very
*		efficient. It would be better to forward transform
*		the source objects domain, but that's not easy!
*  \param	srcObj			Given source object.
* \param	trans			Given affine transform.
* \param	dstErr			Destination ptr for Woolz error
*					code, may be NULL.
*/
static WlzPlaneDomain *WlzAffineTransformPDom(WlzObject *srcObj,
					      WlzAffineTransform *trans,
					      WlzErrorNum *dstErr)
{
  int		insideFlg;
  WlzIVertex2	pMskOrg,
  		pMskSz;
  WlzIVertex3	idx,
  		sPos,
  		dPos;
  WlzIBox3	bBox;
  UBYTE		*lMsk;
  UBYTE		**pMsk = NULL;
  WlzObject	*tObj;
  WlzPlaneDomain *dstPDom = NULL;
  WlzAffineTransform *invTrans = NULL;
  double	tMat[3][3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Forward trasform the bounding box. */
  bBox = WlzBoundingBox3D(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    bBox = WlzAffineTransformBBoxI3(trans, bBox, &errNum);
  }
  /* Create a plane domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstPDom = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
    				 bBox.zMin, bBox.zMax,
    				 bBox.yMin, bBox.yMax,
    				 bBox.xMin, bBox.xMax, &errNum);
  }
  /* Create a bit mask that's the maximum size of any plane. */
  if(errNum == WLZ_ERR_NONE)
  {
    pMskOrg.vtX = bBox.xMin;
    pMskOrg.vtY = bBox.yMin;
    pMskSz.vtX = bBox.xMax - bBox.xMin + 1;
    pMskSz.vtY = bBox.yMax - bBox.yMin + 1;
    if(AlcBit2Malloc(&pMsk, pMskSz.vtY, pMskSz.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Invert the transform. */
  if(errNum == WLZ_ERR_NONE)
  {
    invTrans = WlzAffineTransformInverse(trans, &errNum);
  }
  /* For each plane. */
  if(errNum == WLZ_ERR_NONE)
  {
    idx.vtZ = 0;
    dPos.vtZ = bBox.zMin;
    tMat[0][0] = invTrans->mat[0][0];
    tMat[1][0] = invTrans->mat[1][0];
    tMat[2][0] = invTrans->mat[2][0];
    while((errNum == WLZ_ERR_NONE) && (dPos.vtZ <= bBox.zMax))
    {
      /* Set the plane bit mask. */
      idx.vtY = 0;
      dPos.vtY = bBox.yMin;
      tMat[0][2] = invTrans->mat[0][3] + (invTrans->mat[0][2] * dPos.vtZ);
      tMat[1][2] = invTrans->mat[1][3] + (invTrans->mat[1][2] * dPos.vtZ);
      tMat[2][2] = invTrans->mat[2][3] + (invTrans->mat[2][2] * dPos.vtZ);
      while(dPos.vtY <= bBox.yMax)
      {
	idx.vtX = 0;
	dPos.vtX = bBox.xMin;
	lMsk = *(pMsk + idx.vtY);
	tMat[0][1] = tMat[0][2] + (invTrans->mat[0][1] * dPos.vtY);
	tMat[1][1] = tMat[1][2] + (invTrans->mat[1][1] * dPos.vtY);
	tMat[2][1] = tMat[2][2] + (invTrans->mat[2][1] * dPos.vtY);
	while(dPos.vtX <= bBox.xMax)
	{
	  sPos.vtX = tMat[0][1] + (tMat[0][0] * dPos.vtX);
          sPos.vtY = tMat[1][1] + (tMat[1][0] * dPos.vtX);
          sPos.vtZ = tMat[2][1] + (tMat[2][0] * dPos.vtX);
	  insideFlg = WlzInsideDomain(srcObj, sPos.vtZ, sPos.vtY, sPos.vtX, 
	  			      NULL) != 0;
	  *(lMsk + (idx.vtX / 8)) |= insideFlg << (idx.vtX % 8);
	  ++(idx.vtX);
	  ++(dPos.vtX);
	}
        ++(idx.vtY);
	++(dPos.vtY);
      }
      /* Create a domain from the bit mask. */
      tObj = WlzFromArray2D((void **)pMsk, pMskSz, pMskOrg, WLZ_GREY_BIT,
      			    WLZ_GREY_BIT, 0.0, 1.0, 0, 0, &errNum);
      if(tObj)
      {
        *(dstPDom->domains + idx.vtZ) = WlzAssignDomain(tObj->domain, NULL);
        (void )WlzFreeObj(tObj);
	tObj = NULL;
      }
      ++(idx.vtZ);
      ++(dPos.vtZ);
    }
  }
  if(invTrans)
  {
    (void )WlzFreeAffineTransform(invTrans);
  }
  if(pMsk)
  {
    Alc2Free((void **)pMsk);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(dstPDom, NULL);
  }
  else
  {
    if(dstPDom)
    {
      (void )WlzFreePlaneDomain(dstPDom);
    }
    dstPDom = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstPDom);
}

/*!
* \return				Error number.
* \brief	Sets the given transform's matrix from an
*		affine transform primitives data structure.
*		A composite transform is built from the primitives
*		with the order of composition being scale (applied first),
*		shear, rotation and then translation (applied last),
*		ie:
*		  A = T.R.Sh.Sc, x' = A.x
* \param	tr			 Given affine transform.
* \param	prim			 Given primitives.
*/
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

/*!
* \return				Error number.
* \brief	Sets the given transform's matrix from an
*		affine transform primitives data structure.
*		A composite transform is built from the primitives
*		with the order of composition being scale (applied first),
*		shear, rotation and then translation (applied last),
*		ie:
*		  A = T.R.Sh.Sc, x' = A.x
* \param	tr			Given 2D affine transform.
* \param	prim			Given primitives.
*/
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

/*!
* \return				Error number.
* \brief	Sets the given transform's matrix from an
*		affine transform primitives data structure.
*		A composite transform is built from the primitives
*		with the order of composition being scale (applied first),
*		shear, rotation and then translation (applied last),
*		ie:
*		  A = T.R.Sh.Sc, x' = A.x
* \prim		tr			Given 2D affine transform.
* \param	prim			Given primitives.
*/
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


/*!
* \return				Error number.
* \brief	Sets the given transform's matrix from the given
*		translations.
* \param	tr			Given 2D or 3D affine transform.
* \param	tx	 		Translation along the x-axis.
* \param	ty	 		Translation along the y-axis.
* \param	tz	 		Translation along the z-axis,
*					ignored for 2D transforms.
*/
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

/*!
* \return				New affine transform,
*					NULL on error.
* \brief	Constructs a new affine transform from the given 
*		translations.
* \param	type			Required transform type.
* \param	tx			Translation along the x-axis.
* \param	ty			Translation along the y-axis.
* \param	tz			Translation along the z-axis,
*					ignored for 2D transforms.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				Error number.
* \brief	Sets the given transform's matrix from the given
*		scales.
* \param	tr			Given 2D or 3D affine transform.
* \param	sx			Scale along the x-axis.
* \param	sy			Scale along the y-axis.
* \param	sz			Scale along the z-axis,
*					ignored for 2D transforms.
*/
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


/*!
* \return				New affine transform,
*					NULL on error.
* \brief	Constructs a new affine transform from the given 
*		scales.
* \param	type			Required transform type.
* \param	sx	 		Scale along the x-axis.
* \param	sy	 		Scale along the y-axis.
* \param	sz	 		Scale along the z-axis,
*					ignored for 2D transforms.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				Error number.
* \brief	Sets the given transform's matrix from the given
*		rotations. Although the 3 rotations contain redundant
*		information this may be a useful method for setting
*		rotation transforms. The order of composition is
*		R = Rz.Ry.Rx, x' = R.x.
* \param	tr			Given 2D or 3D affine transform.
* \param	rx			Rotation about the x-axis,
*					ignored for 2D transforms.
* \param	ry			Rotation about the y-axis,
*					ignored for 2D transforms.
* \param	rz			Rotation about the z-axis.
*/
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

/*!
* \return				New affine transform,
*					NULL on error.
* \brief	Constructs a new affine transform from the given 
*		rotations.
* \param	type			Required transform type.
* \param	rx	 		Rotation about the x-axis,
*					ignored for 2D transforms.
* \param	ry	 		Rotation about the y-axis,
*					ignored for 2D transforms.
* \param	rz	 		Rotation about the z-axis.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				Error number.
* \brief	Gets the given 2D transform's primitives from it's
*		matrix.
* \param	tr			Given 2D affine transform.
* \param	prim			Primitives data
*					structure to be set.
*/
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

/*!
* \return	<void>
* \brief	Gets the given 2D transform's primitives from it's
*		matrix.
* \param	tr			Given 2D affine transform.
* \param	prim			Primitives data
*					structure to be set.
*/
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

/*!
* \return				Error number.
* \brief	Sets the given transform from the given matrix.
* \param	trans			Given affine transform.
* \param	matrix			4x4 transform matrix values to
*					be copied.
*/
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
	for(idx0 = 0; idx0 < 4; ++idx0)
	{
	  for(idx1 = 0; idx1 < 4; ++idx1)
	  {
	    trans->mat[idx0][idx1] = matrix[idx0][idx1];
	  }
	}
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

/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Sets a 2D affine transform from the given primitives.
* \param	tr			Given 2D affine transform.
* \param	trX			Column (x) translation.
* \param	trY			Line (y) translation.
* \param	trZ			Plane (z) translation.
* \param	trScale			Scale transformation.
* \param	trTheta			Rotation about z-axis.
* \param	trPhi			Rotation about y-axis.
* \param	trAlpha			Shear strength.
* \param	trPsi			Shear angle in x-y plane.
* \param	trXsi			3D shear angle.
* \param	trInvert		Reflection about y-axis if
*					non-zero.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Makes a new affine transform of the given type and
*		then sets it's matrix.
* \param	type			Required transform type.
* \param	matrix			Given matrix.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Makes a new affine transform from the given primitive
*		transform properties.
* \param	type			Required transform type.
* \param	trX			Column (x) translation.
* \param	trY			Line (y) translation.
* \param	trZ			Plane (z) translation.
* \param	trScale			Scale transformation.
* \param	trTheta			Rotation about z-axis.
* \param	trPhi			Rotation about y-axis.
* \param	trAlpha			Shear strength.
* \param	trPsi			Shear angle in x-y plane.
* \param	trXsi			3D shear angle.
* \param	trInvert		Reflection about y-axis if
*					non-zero.
* \param	dstErr			Destination pointer for error
*					number.
*/
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
/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Makes a new 2D affine transform from the given spin
*		angle and centre of rotation.
* \param	spX			Spin centre column (x).
* \param	spY			Spin centre line (y).
* \param	spTheta			Spin rotation about centre.
*					number.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Makes a new 2D affine transform from the given spin
*		angle, centre of rotation and scale factors.
* \param	spX			Spin centre column (x).
* \param	spY			Spin centre line (y).
* \param	spTheta			Spin rotation about centre.
*					number.
* \param	sqX			Squeeze (x) factor.
* \param	sqY			Squeeze (y) factor.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Copies the given affine transform.
* \param	trans			Given affine transform.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Computes the product of the two given affine
*		transforms.
* \param	tr0			First affine transform.
* \param	tr1			Second affine transform.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return				New affine transform, or NULL
*					on error.
* \brief	Computes the inverse of the given affine transform.
* \param	tr			Given affine transform.
* \param	dstErr			Destination pointer for error
*					number.
*/
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

/*!
* \return		                Non-zero if the given transform
*                                       is an identity transform.
* \brief      	Checks whether the given transform is an identity
*               transform.
* \param	trans			Given affine transform.
* \param        dstErr			Destination pointer for error
*                                       number.
*/
int		WlzAffineTransformIsIdentity(WlzAffineTransform *trans,
					     WlzErrorNum *dstErr)
{
  int           dim,
  		isIdentity = 0;
  double        tD0;
  double        **mat;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const double	delta = 1.0e-06;
 
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
	if((fabs(mat[0][0] - 1.0) <= delta) &&
	   (fabs(mat[0][1]) <= delta) &&
	   (fabs(mat[0][2]) <= delta) &&
	   (fabs(mat[1][0]) <= delta) &&
	   (fabs(mat[1][1] - 1.0) <= delta) &&
	   (fabs(mat[1][2]) <= delta) &&
	   (fabs(mat[2][0]) <= delta) &&
	   (fabs(mat[2][1]) <= delta) &&
	   (fabs(mat[2][2] - 1.0) <= delta)) 
	{
	  isIdentity = 1;
	}
        break;
      case 3:
	if((fabs(mat[0][0] - 1.0) <= delta) &&
	   (fabs(mat[0][1]) <= delta) &&
	   (fabs(mat[0][2]) <= delta) &&
	   (fabs(mat[0][3]) <= delta) &&
	   (fabs(mat[1][0]) <= delta) &&
	   (fabs(mat[1][1] - 1.0) <= delta) &&
	   (fabs(mat[1][2]) <= delta) &&
	   (fabs(mat[1][3]) <= delta) &&
	   (fabs(mat[2][0]) <= delta) &&
	   (fabs(mat[2][1]) <= delta) &&
	   (fabs(mat[2][2] - 1.0) <= delta) &&
	   (fabs(mat[2][3]) <= delta) &&
	   (fabs(mat[3][0]) <= delta) &&
	   (fabs(mat[3][1]) <= delta) &&
	   (fabs(mat[3][2]) <= delta) &&
	   (fabs(mat[3][3] - 1.0) <= delta))
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

/*!
* \return				Transformed object, NULL on
*					error.
* \brief	Applies the given affine transform to the given Woolz
*		object.
* \param	srcObj			Object to be transformed.
* \param	trans			Affine transform to apply.
* \param	interp			Level of interpolation to
*					use.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
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
      case WLZ_CONTOUR:
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
	    case WLZ_CONTOUR:
	      dstDom.ctr = WlzAffineTransformContour(srcObj->domain.ctr,
	      					     trans, 1, &errNum);
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
	else if(WlzAffineTransformIsTranslate(trans, srcObj, NULL))
	{
	  dstObj = WlzAffineTransformIntTranslate(srcObj, trans, &errNum);
	}
	else if(srcDom.core->type == WLZ_PLANEDOMAIN_DOMAIN)
	{
	  dstDom.p = WlzAffineTransformPDom(srcObj, trans, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues,
				 NULL, NULL, &errNum);
	  }
	  if((errNum == WLZ_ERR_NONE) && (srcValues.core != NULL))
	  {
	    errNum = WlzAffineTransformValues3(dstObj, srcObj, trans,
	    				       interp);
	  }
	}
	else
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
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

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzIVertex2.
* \param	trans			Affine transform to apply.
* \param	srcVtx			Vertex to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
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

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzIVertex3.
* \param	trans			Affine transform to apply.
* \param	srcVtx			Vertex to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
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

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzFVertex2.
* \param	trans			Affine transform to apply.
* \param	srcVtx			Vertex to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
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

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzFVertex3.
* \param	trans			Affine transform to apply.
* \param	srcVtx			Vertex to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
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

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzDVertex2.
* \param	trans			Affine transform to apply.
* \param	srcVtx			Vertex to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
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

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzDVertex3.
* \param	trans			Affine transform to apply.
* \param	srcVtx			Vertex to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
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

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzDVertex2 which is a normal.
* \param	trans			Affine transform to apply.
* \param	srcNrm			Normal to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzDVertex2	WlzAffineTransformNormalD2(WlzAffineTransform *trans,
					   WlzDVertex2 srcNrm,
					   WlzErrorNum *dstErr)
{
  WlzDVertex2	dstNrm;
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
    dstNrm.vtX = (trans->mat[0][0] * srcNrm.vtX) +
		 (trans->mat[0][1] * srcNrm.vtY);
    dstNrm.vtY = (trans->mat[1][0] * srcNrm.vtX) +
		 (trans->mat[1][1] * srcNrm.vtY);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstNrm);
}

/*!
* \return				Transformed vertex.
* \brief	Transforms the given WlzDVertex3 which is a normal.
* \param	trans			Affine transform to apply.
* \param	srcNrm			Normal to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzDVertex3	WlzAffineTransformNormalD3(WlzAffineTransform *trans,
					   WlzDVertex3 srcNrm,
					   WlzErrorNum *dstErr)
{
  WlzDVertex3	dstNrm;
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
        dstNrm.vtX = (trans->mat[0][0] * srcNrm.vtX) +
		     (trans->mat[0][1] * srcNrm.vtY);
        dstNrm.vtY = (trans->mat[1][0] * srcNrm.vtX) +
		     (trans->mat[1][1] * srcNrm.vtY);
        break;
      case 3:
	dstNrm.vtX = (trans->mat[0][0] * srcNrm.vtX) +
		     (trans->mat[0][1] * srcNrm.vtY) +
		     (trans->mat[0][2] * srcNrm.vtZ);
	dstNrm.vtY = (trans->mat[1][0] * srcNrm.vtX) +
		     (trans->mat[1][1] * srcNrm.vtY) +
		     (trans->mat[1][2] * srcNrm.vtZ);
	dstNrm.vtZ = (trans->mat[2][0] * srcNrm.vtX) +
		     (trans->mat[2][1] * srcNrm.vtY) +
		     (trans->mat[2][2] * srcNrm.vtZ);
        break;
      default:
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstNrm);
}

/*!
* \return				Transformed bounding box.
* \brief	Transforms the given WlzIBox2.
* \param	tr			Affine transform to apply.
* \param	srcBox			Bounding box to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzIBox2	WlzAffineTransformBBoxI2(WlzAffineTransform *tr,
				         WlzIBox2 srcBox,
					 WlzErrorNum *dstErr)
{
  int		idx;
  WlzIVertex2	tV[4];
  WlzIBox2	dstBox;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  tV[0].vtX = srcBox.xMin; tV[0].vtY = srcBox.yMin;
  tV[1].vtX = srcBox.xMax; tV[1].vtY = srcBox.yMin;
  tV[2].vtX = srcBox.xMin; tV[2].vtY = srcBox.yMax;
  tV[3].vtX = srcBox.xMax; tV[3].vtY = srcBox.yMax;
  idx = 0;
  while((errNum == WLZ_ERR_NONE) && (idx < 4))
  {
    tV[idx] = WlzAffineTransformVertexI2(tr, tV[idx], &errNum);
    ++idx;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstBox.xMin = dstBox.xMax = tV[0].vtX;
    dstBox.yMin = dstBox.yMax = tV[0].vtY;
    for(idx = 1; idx < 4; ++idx)
    {
      if(tV[idx].vtX < dstBox.xMin)
      {
	dstBox.xMin = tV[idx].vtX;
      }
      else if(tV[idx].vtX > dstBox.xMax)
      {
        dstBox.xMax = tV[idx].vtX;
      }
      if(tV[idx].vtY < dstBox.yMin)
      {
	dstBox.yMin = tV[idx].vtY;
      }
      else if(tV[idx].vtY > dstBox.yMax)
      {
        dstBox.yMax = tV[idx].vtY;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBox);
}

/*!
* \return				Transformed bounding box.
* \brief	Transforms the given WlzDBox2.
* \param	tr			Affine transform to apply.
* \param	srcBox			Bounding box to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzDBox2	WlzAffineTransformBBoxD2(WlzAffineTransform *tr,
				         WlzDBox2 srcBox,
					 WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex2	tV[4];
  WlzDBox2	dstBox;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  tV[0].vtX = srcBox.xMin; tV[0].vtY = srcBox.yMin;
  tV[1].vtX = srcBox.xMax; tV[1].vtY = srcBox.yMin;
  tV[2].vtX = srcBox.xMin; tV[2].vtY = srcBox.yMax;
  tV[3].vtX = srcBox.xMax; tV[3].vtY = srcBox.yMax;
  idx = 0;
  while((errNum == WLZ_ERR_NONE) && (idx < 4))
  {
    tV[idx] = WlzAffineTransformVertexD2(tr, tV[idx], &errNum);
    ++idx;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstBox.xMin = dstBox.xMax = tV[0].vtX;
    dstBox.yMin = dstBox.yMax = tV[0].vtY;
    for(idx = 1; idx < 4; ++idx)
    {
      if(tV[idx].vtX < dstBox.xMin)
      {
	dstBox.xMin = tV[idx].vtX;
      }
      else if(tV[idx].vtX > dstBox.xMax)
      {
        dstBox.xMax = tV[idx].vtX;
      }
      if(tV[idx].vtY < dstBox.yMin)
      {
	dstBox.yMin = tV[idx].vtY;
      }
      else if(tV[idx].vtY > dstBox.yMax)
      {
        dstBox.yMax = tV[idx].vtY;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBox);
}

/*!
* \return				Transformed bounding box.
* \brief	Transforms the given WlzIBox3.
* \param	tr			Affine transform to apply.
* \param	srcBox			Bounding box to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzIBox3	WlzAffineTransformBBoxI3(WlzAffineTransform *tr,
				         WlzIBox3 srcBox,
					 WlzErrorNum *dstErr)
{
  int		idx;
  WlzIVertex3	tV[8];
  WlzIBox3	dstBox;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  tV[0].vtX = srcBox.xMin; tV[0].vtY = srcBox.yMin; tV[0].vtZ = srcBox.zMin;
  tV[1].vtX = srcBox.xMax; tV[1].vtY = srcBox.yMin; tV[1].vtZ = srcBox.zMin;
  tV[2].vtX = srcBox.xMin; tV[2].vtY = srcBox.yMax; tV[2].vtZ = srcBox.zMin;
  tV[3].vtX = srcBox.xMax; tV[3].vtY = srcBox.yMax; tV[3].vtZ = srcBox.zMin;
  tV[4].vtX = srcBox.xMin; tV[4].vtY = srcBox.yMin; tV[4].vtZ = srcBox.zMax;
  tV[5].vtX = srcBox.xMax; tV[5].vtY = srcBox.yMin; tV[5].vtZ = srcBox.zMax;
  tV[6].vtX = srcBox.xMin; tV[6].vtY = srcBox.yMax; tV[6].vtZ = srcBox.zMax;
  tV[7].vtX = srcBox.xMax; tV[7].vtY = srcBox.yMax; tV[7].vtZ = srcBox.zMax;
  idx = 0;
  while((errNum == WLZ_ERR_NONE) && (idx < 8))
  {
    tV[idx] = WlzAffineTransformVertexI3(tr, tV[idx], &errNum);
    ++idx;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstBox.xMin = dstBox.xMax = tV[0].vtX;
    dstBox.yMin = dstBox.yMax = tV[0].vtY;
    dstBox.zMin = dstBox.zMax = tV[0].vtZ;
    for(idx = 1; idx < 8; ++idx)
    {
      if(tV[idx].vtX < dstBox.xMin)
      {
	dstBox.xMin = tV[idx].vtX;
      }
      else if(tV[idx].vtX > dstBox.xMax)
      {
        dstBox.xMax = tV[idx].vtX;
      }
      if(tV[idx].vtY < dstBox.yMin)
      {
	dstBox.yMin = tV[idx].vtY;
      }
      else if(tV[idx].vtY > dstBox.yMax)
      {
        dstBox.yMax = tV[idx].vtY;
      }
      if(tV[idx].vtZ < dstBox.zMin)
      {
	dstBox.zMin = tV[idx].vtZ;
      }
      else if(tV[idx].vtZ > dstBox.zMax)
      {
        dstBox.zMax = tV[idx].vtZ;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBox);
}

/*!
* \return				Transformed bounding box.
* \brief	Transforms the given WlzDBox3.
* \param	tr			Affine transform to apply.
* \param	srcBox			Bounding box to be transformed.
* \param	dstErr			Destination pointer for error
*					number, may be NULL.
*/
WlzDBox3	WlzAffineTransformBBoxD3(WlzAffineTransform *tr,
				         WlzDBox3 srcBox,
					 WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex3	tV[8];
  WlzDBox3	dstBox;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  tV[0].vtX = srcBox.xMin; tV[0].vtY = srcBox.yMin; tV[0].vtZ = srcBox.zMin;
  tV[1].vtX = srcBox.xMax; tV[1].vtY = srcBox.yMin; tV[1].vtZ = srcBox.zMin;
  tV[2].vtX = srcBox.xMin; tV[2].vtY = srcBox.yMax; tV[2].vtZ = srcBox.zMin;
  tV[3].vtX = srcBox.xMax; tV[3].vtY = srcBox.yMax; tV[3].vtZ = srcBox.zMin;
  tV[4].vtX = srcBox.xMin; tV[4].vtY = srcBox.yMin; tV[4].vtZ = srcBox.zMax;
  tV[5].vtX = srcBox.xMax; tV[5].vtY = srcBox.yMin; tV[5].vtZ = srcBox.zMax;
  tV[6].vtX = srcBox.xMin; tV[6].vtY = srcBox.yMax; tV[6].vtZ = srcBox.zMax;
  tV[7].vtX = srcBox.xMax; tV[7].vtY = srcBox.yMax; tV[7].vtZ = srcBox.zMax;
  idx = 0;
  while((errNum == WLZ_ERR_NONE) && (idx < 8))
  {
    tV[idx] = WlzAffineTransformVertexD3(tr, tV[idx], &errNum);
    ++idx;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstBox.xMin = dstBox.xMax = tV[0].vtX;
    dstBox.yMin = dstBox.yMax = tV[0].vtY;
    dstBox.zMin = dstBox.zMax = tV[0].vtZ;
    for(idx = 1; idx < 8; ++idx)
    {
      if(tV[idx].vtX < dstBox.xMin)
      {
	dstBox.xMin = tV[idx].vtX;
      }
      else if(tV[idx].vtX > dstBox.xMax)
      {
        dstBox.xMax = tV[idx].vtX;
      }
      if(tV[idx].vtY < dstBox.yMin)
      {
	dstBox.yMin = tV[idx].vtY;
      }
      else if(tV[idx].vtY > dstBox.yMax)
      {
        dstBox.yMax = tV[idx].vtY;
      }
      if(tV[idx].vtZ < dstBox.zMin)
      {
	dstBox.zMin = tV[idx].vtZ;
      }
      else if(tV[idx].vtZ > dstBox.zMax)
      {
        dstBox.zMax = tV[idx].vtZ;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBox);
}

/*!
* @}
*/

