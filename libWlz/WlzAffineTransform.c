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
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static WlzObject *WlzAffineTransformIntTranslate(WlzObject *,
						 WlzAffineTransform *,
						 WlzErrorNum *);
static WlzPolygonDomain *WlzAffineTransformPoly(WlzPolygonDomain *,
					        WlzAffineTransform *,
						WlzErrorNum *);
static WlzBoundList *WlzAffineTransformBoundList(WlzBoundList *,
						 WlzAffineTransform *,
					    WlzErrorNum *);
static WlzErrorNum WlzAffineTransformValues2D(WlzObject *,
					      WlzObject *,
					      WlzAffineTransform *,
					      WlzInterpolationType);

/************************************************************************
* Function:	WlzAffineTransformIsTranslate				*
* Returns:	int:			Non-zero if translation.	*
* Purpose:	Tests wether the given affine transform is a simple	*
*		translation.						*
*		Because this is a static function the parameters are	*
*		not checked.						*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
*		WlzObject *obj:		Given 2D domain object.		*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
int		WlzAffineTransformIsTranslate(WlzAffineTransform *trans,
					      WlzObject *obj,
					      WlzErrorNum *dstErr)
{
  int		translate = 1;
  double	trX,
  		trY;
  WlzIBox2	box;
  const double	translateDelta = 0.1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(trans->type != WLZ_TRANSFORM_2D_AFFINE)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    if(obj && (obj->type == WLZ_2D_DOMAINOBJ) && obj->domain.core)
    {
      box.xMin = obj->domain.i->kol1;
      box.yMin = obj->domain.i->line1;
      box.xMax = obj->domain.i->lastkl;
      box.yMax = obj->domain.i->lastln;
    }
    else
    {
      /* Should be something better here */
      box.xMin = 0;
      box.yMin = 0;
      box.xMax = 1;
      box.yMax = 1;
    }
    /* Check translation */
    trX = trans->mat[0][2];
    trY = trans->mat[1][2];
    if((fabs(trX - WLZ_NINT(trX)) > translateDelta) ||
       (fabs(trY - WLZ_NINT(trY)) > translateDelta))
    {
      translate = 0;
    }
    else
    {
      /* Check rotation, scale, shear and invert */
      trX = box.xMin * (trans->mat[0][0]) + box.yMin * (trans->mat[0][1]);
      trY = box.xMin * (trans->mat[1][0]) + box.yMin * (trans->mat[1][1]);
      if((fabs(trX - box.xMin) > translateDelta) ||
	 (fabs(trY - box.yMin) > translateDelta))
      {
	translate = 0;
      }
      else
      {
	trX = box.xMin * (trans->mat[0][0]) + box.yMax * (trans->mat[0][1]);
	trY = box.xMin * (trans->mat[1][0]) + box.yMax * (trans->mat[1][1]);
	if((fabs(trX - box.xMin) > translateDelta) ||
	   (fabs(trY - box.yMax) > translateDelta))
	{
	  translate = 0;
	}
	else
	{
	  trX = box.xMax * (trans->mat[0][0]) + box.yMin * (trans->mat[0][1]);
	  trY = box.xMax * (trans->mat[1][0]) + box.yMin * (trans->mat[1][1]);
	  if((fabs(trX - box.xMax) > translateDelta) ||
	     (fabs(trY - box.yMin) > translateDelta))
	  {
	    translate = 0;
	  }
	  else
	  {
	    trX = box.xMax * (trans->mat[0][0]) + box.yMax * (trans->mat[0][1]);
	    trY = box.xMax * (trans->mat[1][0]) + box.yMax * (trans->mat[1][1]);
	    if((fabs(trX - box.xMax) > translateDelta) ||
	       (fabs(trY - box.yMax) > translateDelta))

	    {
	      translate = 0;
	    }
	  }
	}
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(translate);
}

/************************************************************************
* Function:	WlzAffineTransformIntTranslate				*
* Returns:	WlzObject *:		Translated object or NULL on	*
*					error.				*
* Purpose:	Translates the given 2D domain object with an integral	*
*		translation.						*
*		Because this is a static function the parameters are	*
*		not checked.						*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
*		WlzObject *obj:		   Given 2D domain object.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
static WlzObject *WlzAffineTransformIntTranslate(WlzObject *srcObj,
					         WlzAffineTransform *trans,
					      	 WlzErrorNum *dstErr)
{
  int		trX,
  		trY;
  WlzObject	*newObj = NULL;
  WlzDomain	newDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  trX = WLZ_NINT(trans->mat[0][2]);
  trY = WLZ_NINT(trans->mat[1][2]);
  newObj = WlzNewGrey(srcObj, &errNum);
  if( errNum == WLZ_ERR_NONE ){
    newDom.i = WlzNewIDomain(srcObj->domain.i->type, srcObj->domain.i,
    			     &errNum);
  }

  if( errNum == WLZ_ERR_NONE ){

    /* Translate the interval domain */
    newDom.i->kol1 += trX;
    newDom.i->lastkl += trX;
    newDom.i->line1 += trY;
    newDom.i->lastln += trY;
    /* Replace domain with the new translated domain */
    (void )WlzFreeDomain(newObj->domain);
    newObj->domain = WlzAssignDomain(newDom, &errNum);
    /* Translate the valuetable */
    if(newObj->values.core)
    {
      switch(WlzGreyTableTypeToTableType(newObj->values.core->type,
					 &errNum))
      {
	case WLZ_GREY_TAB_RAGR:
	  newObj->values.v->kol1   += trX;
	  newObj->values.v->line1  += trY;
	  newObj->values.v->lastln += trY;
	  break;
	case WLZ_GREY_TAB_RECT:
	  newObj->values.r->kol1   += trX;
	  newObj->values.r->line1  += trY;
	  newObj->values.r->lastln += trY;
	  break;
	case WLZ_GREY_TAB_INTL:
	  newObj->values.i->kol1   += trX;
	  newObj->values.i->line1  += trY;
	  newObj->values.i->lastln += trY;
	  break;
	default:
	  errNum = WLZ_ERR_VALUES_TYPE;
	  break;
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(newObj)
    {
      (void )WlzFreeObj(newObj);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newObj);
}

/************************************************************************
* Function:	WlzAffineTransformPoly					*
* Returns:	WlzPolygonDomain *:	Transformed polygon domain or	*
*					NULL on error.			*
* Purpose:	Transforms the given polygon domain.			*
*		Because this is a static function the parameters (other	*
*		than polygon type) are not checked.			*
* Global refs:	-							*
* Parameters:	WlzPolygonDomain *srcPoly: Given polygon domain.	*
*		WlzAffineTransform *trans: Given affine transform.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
static WlzPolygonDomain *WlzAffineTransformPoly(WlzPolygonDomain *srcPoly,
					        WlzAffineTransform *trans,
					      	WlzErrorNum *dstErr)
{
  int		count;
  double	cx,
  		cy,
		dx,
		dy,
		sx,
		sy,
		tx,
		ty;
  WlzIVertex2	*srcVtxI,
  		*dstVtxI;
  WlzFVertex2	*srcVtxF,
  		*dstVtxF;
  WlzDVertex2	*srcVtxD,
  		*dstVtxD;
  WlzPolygonDomain *dstPoly = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((srcPoly->type != WLZ_POLYGON_INT) &&
     (srcPoly->type != WLZ_POLYGON_FLOAT) &&
     (srcPoly->type != WLZ_POLYGON_DOUBLE))
  {
    errNum = WLZ_ERR_POLYGON_TYPE;
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
* Function:	WlzAffineTransformBoundList				*
* Returns:	WlzBoundList *:		Transformed boundary list or	*
*					NULL on error.			*
* Purpose:	Transforms the given boundary list.			*
*		Because this is a static function the parameters are	*
*		not checked.						*
* Global refs:	-							*
* Parameters:	WlzBoundList *srcBound: Given boundary list.		*
*		WlzAffineTransform *trans: Given affine transform.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
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
  else
  {
    dstBound->type = srcBound->type;
    dstBound->wrap = srcBound->wrap;
    /* transform the polygon */
    if((dstBound->poly = WlzAffineTransformPoly(srcBound->poly, trans,
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
* Function:	WlzAffineTransformValues2D				*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Creates a new value table, fills in the values and	*
*		adds it to the given new object.			*
*		Because this is a static function the parameters are	*
*		not checked.						*
* Global refs:	-							*
* Parameters:	WlzObject *newObj:	Partialy transformed object	*
*					with a valid domain.		*
*		WlzObject *srcObj:	2D domain object which is being	*
*					transformed.			*
*		WlzAffineTransform *trans: Given affine transform.	*
*		WlzInterpolationType interp: Level of interpolation to	*
*					use.				*
************************************************************************/
static WlzErrorNum WlzAffineTransformValues2D(WlzObject *newObj,
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
  if((invTrans = WlzAffineTransformInverse(trans, &errNum)) != NULL)
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
* Function:	WlzAffineTransformMatrixUpdate				*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Updates the given transform's matrix from it's		*
*		primitives.						*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
************************************************************************/
WlzErrorNum	WlzAffineTransformMatrixUpdate(WlzAffineTransform *trans)
{
  int		idx0;
  double	tS,
  		tCos1,
		tCos2,
  		tSin1,
		tSin2,
		tCosSin2;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(trans->type)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
	tCos1 = cos(trans->theta);
	tSin1 = sin(trans->theta);
	tCos2 = cos(trans->psi);
	tSin2 = sin(trans->psi);
	tCosSin2 = tCos2 * tSin2 * trans->alpha;
	tSin2 *= tSin2 * trans->alpha;
	tCos2 *= tCos2 * trans->alpha;
	tS = trans->scale;
	trans->mat[0][0] = tS * ((tCos1 * (1 - tCosSin2)) + (tSin1 * tSin2));
	trans->mat[0][1] = tS * ((-tSin1 * (1 + tCosSin2)) + (tCos1 * tCos2));
	trans->mat[1][0] = tS * ((tSin1 * (1 - tCosSin2)) - (tCos1 * tSin2));
	trans->mat[1][1] = tS * ((tCos1 * (1 + tCosSin2)) + (tSin1 * tCos2));
	trans->mat[0][2] = trans->tx;
	trans->mat[1][2] = trans->ty;
	for(idx0 = 0; idx0 < 4; ++idx0)
	{
	  trans->mat[2][idx0] = 0.0;
	  trans->mat[3][idx0] = 0.0;
	  trans->mat[idx0][3] = 0.0;
	}
	trans->mat[2][2] = 1.0;
	trans->mat[3][3] = 1.0;
	if(trans->invert)
	{
	  for(idx0 = 0; idx0 < 4; ++idx0)
	  {
	    trans->mat[0][idx0] *= -1.0;
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
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformPrimUpdate				*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Updates the given transform's primitives from it's	*
*		matrix.							*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
************************************************************************/
WlzErrorNum	WlzAffineTransformPrimUpdate(WlzAffineTransform *trans)
{
  double  	s2,
  		tD0,
		tD1,
		tD2;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(trans->type)
    { 
      case WLZ_TRANSFORM_2D_AFFINE:
	/* Test for inversion */
	s2 = (trans->mat[0][0] * trans->mat[1][1]) - 
	     (trans->mat[0][1] * trans->mat[1][0]);
	if(s2 < 0.0)
	{
	  trans->invert = 1;
	  trans->mat[0][0] *= -1.0;
	  trans->mat[0][1] *= -1.0;
	  trans->mat[0][2] *= -1.0;
	  s2 *= -1.0;
	}
	else
	{
	  trans->invert = 0;
	}
	tD0 = trans->mat[0][0] + trans->mat[1][1];
	tD1 = trans->mat[0][1] - trans->mat[1][0];
	/* Scale and shear strength */ 
	if(fabs(s2) > DBL_EPSILON)
	{
	  trans->scale = sqrt(s2);
	  tD2 = (((tD0 * tD0) + (tD1 * tD1)) / s2)  - 4.0;
	  trans->alpha = (tD2 > DBL_EPSILON)? sqrt(tD2): 0.0;
	}
	else
	{
	  trans->scale =  0.0;
	  trans->alpha = 0.0;
	}
	/* Rotation */
	trans->theta = atan2((trans->alpha * tD0 / 2.0) - tD1,
			     (trans->alpha * tD1 / 2.0) + tD0);
	/* Shear angle */
	if(fabs(trans->alpha) > DBL_EPSILON)
	{
	  tD2 = tan(trans->theta);
	  trans->psi = -atan2(trans->mat[0][0] - trans->mat[1][1] +
			      ((trans->mat[0][1] + trans->mat[1][0]) * tD2),
			      (trans->mat[0][1] +
			       (trans->mat[1][1] * tD2)) * 2.0);
	}
	else
	{
	  trans->psi = 0.0;
	}
	/* Translation */
	trans->tx = trans->mat[0][2];
	trans->ty = trans->mat[1][2];
	/* Restore matrix if inversion */
	if(trans->invert)
	{
	  trans->mat[0][0] *= -1.0;
	  trans->mat[0][1] *= -1.0;
	  trans->mat[0][2] *= -1.0;
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
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformMatrixSet				*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Sets the given transform from the given matrix.		*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
*		WlzIVertex2 arraySizeMat: Matrix size (4 x 4).		*
*		double ** arrayMat:	Transform matrix values to be	*
*					copied.				*
************************************************************************/
WlzErrorNum	WlzAffineTransformMatrixSet(WlzAffineTransform *trans,
					    WlzIVertex2 arraySizeMat,
					    double **arrayMat)
{
  return(WlzAffineTransformMatrixSet4X4(trans, arrayMat));
}

/************************************************************************
* Function:	WlzAffineTransformMatrixSet4X4				*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Sets the given transform from the given matrix.		*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
*		double **matrix:	Transform matrix values to be	*
*					copied.				*
************************************************************************/
WlzErrorNum	WlzAffineTransformMatrixSet4X4(WlzAffineTransform *trans,
					       double **matrix)
{
  int		idx0,
  		idx1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformMatrixSet4X4 FE 0x%lx 0x%lx\n",
	   (unsigned long )trans, (unsigned long )matrix));
  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
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
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzAffineTransformPrimUpdate(trans);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformMatrixSet4X4 FX %d\n",
	   (int )errNum));
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformPrimSet				*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Makes a new affine transform from the given primitive	*
*		transform properties.					*
* Global refs:	-							*
* Parameters:	WlzTransformType type:	Required transform type.	*
*		double trX:		Column (x) translation.		*
*		double trY:		Line (y) translation.		*
*		double trZ:		Plane (z) translation.		*
*		double trScale:		Scale transformation.		*
*		double trTheta:		Rotation about z-axis.		*
*		double trPhi:		Rotation about y-axis.		*
*		double trAlpha:		Shear strength.			*
*		double trPsi:		Shear angle in x-y plane.	*
*		double trXsi:		3D shear angle.			*
*		int trInvert:		Reflection about y-axis if 	*
*					non-zero.			*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzErrorNum	WlzAffineTransformPrimSet(WlzAffineTransform *trans,
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
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformPrimSet FE "
	  "0x%lx %g %g %g %g %g %g %g %g %g %d\n",
	   (unsigned long )trans, trX, trY, trZ, trScale, trTheta, trPhi,
	   trAlpha, trPsi, trXsi, trInvert));
  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    trans->tx = trX;
    trans->ty = trY;
    trans->tz = trZ;
    trans->scale = trScale;
    trans->theta = trTheta;
    trans->phi = trPhi;
    trans->alpha = trAlpha;
    trans->psi = trPsi;
    trans->xsi = trXsi;
    trans->invert = trInvert;
    errNum = WlzAffineTransformMatrixUpdate(trans);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformPrimSet FX %d\n",
	   (int )errNum));
  return(errNum);
}

/************************************************************************
* Function:	WlzAffineTransformFromMatrix				*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Makes a new affine transform of the given type and	*
*		then sets it's matrix.					*
* Global refs:	-							*
* Parameters:	WlzTransformType type:	Required transform type.	*
*		WlzIVertex2 arraySizeMat: Size of array (4 x 4).	*
*		double **arrayMat:	Given matrix.			*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromMatrix(WlzTransformType type,
						 WlzIVertex2 arraySizeMat,
						 double **arrayMat,
						 WlzErrorNum *dstErr)
{
  return(WlzAffineTransformFromMatrix4X4(type, arrayMat, dstErr));
}

/************************************************************************
* Function:	WlzAffineTransformFromMatrix4X4				*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Makes a new affine transform of the given type and	*
*		then sets it's matrix.					*
* Global refs:	-							*
* Parameters:	WlzTransformType type:	Required transform type.	*
*		double **matrix:	Given matrix.			*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromMatrix4X4(WlzTransformType type,
						    double **matrix,
						    WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTrans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromMatrix4X4 FE %d 0x%lx 0x%lx\n",
	   (int )type, (unsigned long )matrix,
	   (unsigned long )dstErr));
  if(matrix == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if((newTrans = WlzMakeAffineTransform(type, &errNum)) == NULL)
    {
      errNum = WLZ_ERR_UNSPECIFIED;
    }
    else if((errNum = WlzAffineTransformMatrixSet4X4(newTrans,
    						     matrix)) != WLZ_ERR_NONE)
    {
      WlzFreeAffineTransform(newTrans);
      newTrans = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromMatrix4X4 FX 0x%lx\n",
	   (unsigned long )newTrans));
  return(newTrans);
}

/************************************************************************
* Function:	WlzAffineTransformFromPrim				*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Makes a new affine transform from the given primitive	*
*		transform properties.					*
* Global refs:	-							*
* Parameters:	WlzTransformType type:	Required transform type.	*
*		double trX:		Column (x) translation.		*
*		double trY:		Line (y) translation.		*
*		double trZ:		Plane (z) translation.		*
*		double trScale:		Scale transformation.		*
*		double trTheta:		Rotation about z-axis.		*
*		double trPhi:		Rotation about y-axis.		*
*		double trAlpha:		Shear strength.			*
*		double trPsi:		Shear angle in x-y plane.	*
*		double trXsi:		3D shear angle.			*
*		int trInvert:		Reflection about y-axis if 	*
*					non-zero.			*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromPrim(WlzTransformType type,
						double trX,
						double trY,
						double trZ,
						double trScale,
						double trTheta,
						double trPhi,
						double trAlpha,
						double trPsi,
						double trXsi,
						int trInvert,
						WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTrans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromPrim FE "
	  "%d %g %g %g %g %g %g %g %g %g %d 0x%lx\n",
	   (int )type, trX, trY, trZ, trScale, trTheta, trPhi,
	   trAlpha, trPsi, trXsi, trInvert, (unsigned long )dstErr));
  if((newTrans = WlzMakeAffineTransform(type, &errNum)) == NULL)
  {
/*    errNum = WlzGetErrno();*/
    errNum = WLZ_ERR_UNSPECIFIED;
  }
  else if((errNum = WlzAffineTransformPrimSet(newTrans,
					      trX, trY, trZ,
					      trScale, trTheta, trPhi,
					      trAlpha, trPsi, trXsi,
					      trInvert)) != WLZ_ERR_NONE)
  {
    WlzFreeAffineTransform(newTrans);
    newTrans = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromPrim FX 0x%lx\n",
	   (unsigned long )newTrans));
  return(newTrans);
}
/************************************************************************
* Function:	WlzAffineTransformFromSpin				*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Makes a new 2D affine transform from the given spin	*
*		angle and centre of rotation.				*
* Global refs:	-							*
* Parameters:	double spX:		Spin centre column (x).		*
*		double spY:		Spin centre line (y).		*
*		double spTheta:		Spin rotation about centre.	*
*					number.				*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
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
  newTrans = WlzAffineTransformFromPrim(WLZ_TRANSFORM_2D_AFFINE, trX, trY, 0.0,
					1.0, spTheta, 0.0, 0.0, 0.0, 0.0, 0,
					dstErr);
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromSpin FX 0x%lx\n",
	   (unsigned long )newTrans));
  return(newTrans);
}

/************************************************************************
* Function:	WlzAffineTransformFromSpinSqueeze			*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Makes a new 2D affine transform from the given spin	*
*		angle, centre of rotation and scale factors.		*
* Global refs:	-							*
* Parameters:	double spX:		Spin centre column (x).		*
*		double spY:		Spin centre line (y).		*
*		double spTheta:		Spin rotation about centre.	*
*					number.				*
*		double sqX:		Squeeze (x) factor.		*
*		double sqY:		Squeeze (y) factor.		*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzAffineTransform *WlzAffineTransformFromSpinSqueeze(double spX, double spY,
					       double spTheta,
					       double sqX, double sqY,
					       WlzErrorNum *dstErr)
{
  double	sinTheta,
  		cosTheta;
  double	**matrix;
  WlzAffineTransform *newTrans;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
          ("WlzAffineTransformFromSpinSqueeze %g %g %g %g %g 0x%lx\n",
	   spX, spY, spTheta, sqX, sqY, dstErr));
  if((newTrans = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE,
					&errNum)) == NULL)
  {
    errNum = WLZ_ERR_UNSPECIFIED;
  }
  else {
    matrix = newTrans->mat;
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
    errNum = WlzAffineTransformPrimUpdate(newTrans);
  }

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformFromSpinSqueeze FX 0x%lx\n",
	   (unsigned long )newTrans));
  if( dstErr ){
    *dstErr = errNum;
  }
  return(newTrans);
}

/************************************************************************
* Function:	WlzAffineTransformCopy					*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Copies the given affine transform.			*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzAffineTransform *WlzAffineTransformCopy(WlzAffineTransform *trans,
					   WlzErrorNum *dstErr)
{
  WlzAffineTransform *newTrans = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformCopy FE 0x%lx 0x%lx\n",
	   (unsigned long )trans, (unsigned long )dstErr));
  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    newTrans = WlzAffineTransformFromMatrix4X4(trans->type, trans->mat,
    					       &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformCopy FX 0x%lx\n",
	   (unsigned long )newTrans));
  return(newTrans);
}

/************************************************************************
* Function:	WlzAffineTransformProduct				*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Computes the product of the two given affine		*
*		transforms.						*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans0: First affine transform.	*
*		WlzAffineTransform *trans1: Second affine transform.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzAffineTransform *WlzAffineTransformProduct(WlzAffineTransform *trans0,
					      WlzAffineTransform *trans1,
					      WlzErrorNum *dstErr)
{
  int		idx0,
  		idx1,
		idx2;
  double	*matrix[4];
  double	buffer[16];
  WlzAffineTransform *product = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformProduct FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )trans0, (unsigned long )trans1,
	   (unsigned long )dstErr));
  if((trans0 == NULL) || (trans1 == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(trans0->type != trans1->type)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    /* Set up the matrix pointers and initialise the array */
    for(idx0=0; idx0 < 4; idx0++)
    {
      matrix[idx0] = &buffer[idx0*4];
    }
    for(idx0=0; idx0 < 16; idx0++)
    {
      buffer[idx0] = 0.0;
    }
    switch(trans0->type)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
	for(idx0=0; idx0 < 3; idx0++)
	{
	  for(idx1=0; idx1 < 3; idx1++)
	  {
	    for(idx2=0; idx2 <3; idx2++)
	    {
	      matrix[idx0][idx1] += trans1->mat[idx0][idx2] *
	      			    trans0->mat[idx2][idx1];
	    }
	  }
	}
	break;
      case WLZ_TRANSFORM_3D_AFFINE:
	for(idx0=0; idx0 < 4; idx0++)
	{
	  for(idx1=0; idx1 < 4; idx1++)
	  {
	    for(idx2=0; idx2 <4; idx2++)
	    {
	      matrix[idx0][idx1] += trans1->mat[idx0][idx2] *
	      			    trans0->mat[idx2][idx1];
	    }
	  }
	}
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    product = WlzAffineTransformFromMatrix4X4(trans0->type, matrix, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformProduct FX 0x%lx\n",
	   (unsigned long )product));
  return(product);
}

/************************************************************************
* Function:	WlzAffineTransformInverse				*
* Returns:	WlzAffineTransform *:	New affine transform, or NULL	*
*					on error.			*
* Purpose:	Computes the inverse of the given affine transform.	*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Given affine transform.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
WlzAffineTransform	*WlzAffineTransformInverse(WlzAffineTransform *trans,
					           WlzErrorNum *dstErr)
{
  WlzAffineTransform *invTrans = NULL;
  double	*matrix[4];
  double	buffer[16];
  int		idx0,
  		idx1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformInverse FE 0x%lx 0x%lx\n",
	   (unsigned long )trans, (unsigned long )dstErr));
  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    /* Set up the matrix pointers and initialise the array */
    *(matrix + 0) = buffer + 0;
    *(matrix + 1) = buffer + 4;
    *(matrix + 2) = buffer + 8;
    *(matrix + 3) = buffer + 12;
    /* Copy the matrix */
    for(idx0 = 0; idx0 < 4; ++idx0)
    {
      for(idx1 = 0; idx1 < 4; ++idx1)
      {
	matrix[idx0][idx1] = trans->mat[idx0][idx1];
      }
    }
    /* Find the inverse */
    switch(trans->type)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
	if(AlgMatrixLUInvert(matrix, 3) != ALG_ERR_NONE)
	{
	  errNum = WLZ_ERR_TRANSFORM_DATA;
	}
	break;
      case WLZ_TRANSFORM_3D_AFFINE:
	if(AlgMatrixLUInvert(matrix, 4) != ALG_ERR_NONE)
	{
	  errNum = WLZ_ERR_TRANSFORM_DATA;
	}
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      invTrans = WlzAffineTransformFromMatrix4X4(trans->type, matrix,
					         &errNum);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzAffineTransformInverse FX 0x%lx\n",
	   (unsigned long )invTrans));
  return(invTrans);
}

/************************************************************************
* Function:     WlzAffineTransformIsIdentity                            *
* Returns:      int:                    Non-zero if the given transform *
*                                       is an identity transform.       *
* Purpose:      Checks whether the given transform is an identity       *
*               transform.                                              *
* Global refs:  -                                                       *
* Parameters:   WlzAffineTransform *trans: Given affine transform.      *
*               WlzErrorNum *dstErr:    Destination pointer for error   *
*                                       number.                         *
************************************************************************/
int		WlzAffineTransformIsIdentity(WlzAffineTransform *trans,
					     WlzErrorNum *dstErr)
{
  int           isIdentity = 0;
  double        tD0;
  double        **mat;
  const double  zeroM = -(DBL_EPSILON),
		zeroP = DBL_EPSILON,
		oneM = 1.0 - DBL_EPSILON,
		oneP = 1.0 + DBL_EPSILON;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
 
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformIsIdentity FE 0x%lx 0x%lx\n",
	   (unsigned long )trans, (unsigned long )dstErr));
  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((trans->type != WLZ_TRANSFORM_2D_AFFINE) &&
	  (trans->type != WLZ_TRANSFORM_2D_AFFINE))
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if(((tD0 = (mat = trans->mat)[0][2]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[1][2]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[0][0]) >= oneM) && (tD0 <= oneP) &&
	  ((tD0 = mat[0][1]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[0][3]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[1][0]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[1][1]) >= oneM) && (tD0 <= oneP) &&
	  ((tD0 = mat[1][3]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[2][0]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[2][1]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[2][2]) >= oneM) && (tD0 <= oneP) &&
	  ((tD0 = mat[2][3]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[3][0]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[3][1]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[3][2]) >= zeroM) && (tD0 <= zeroP) &&
	  ((tD0 = mat[3][3]) >= oneM) && (tD0 <= oneP))
  {
    isIdentity = 1;
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
* Function:	WlzAffineTransformObj					*
* Returns:	WlzObject *:		Transformed object, NULL on	*
*					error.				*
* Purpose:	Applies the given affine transform to the given Woolz	*
*		object.							*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Object to be transformed.	*
*		WlzAffineTransform *trans: Affine transform to apply. 	*
*		WlzInterpolationType interp: Level of interpolation to	*
*					use.				*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
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
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzAffineTransformObj FE 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )srcObj, (unsigned long )trans,
	   (int )interp, (unsigned long )dstErr));
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
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
	      dstDom.poly = WlzAffineTransformPoly(srcObj->domain.poly,
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
	else
	{
	  if(WlzAffineTransformIsTranslate(trans, srcObj, NULL))
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
	      dstObj = WlzBoundToObj(tObj1->domain.b, WLZ_SIMPLE_FILL,
				     &errNum);
	      WlzFreeObj(tObj1);
	      tObj1 = NULL;
	    }
	    if((errNum == WLZ_ERR_NONE) &&
	       (srcObj->values.core) )
	    {
	      errNum = WlzAffineTransformValues2D(dstObj, srcObj, trans,
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
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	/* this is a kludge. Currently a full 3D transform has not been
	   implemented and in this code the transform is applied to each
	   plane in turn. If the transform has a z translation set then
	   the resultant object is shifted in z - BUT NOT SCALED!!!! 
	   In principle the z-transformation associated with 2D transform 
	   is meaningless and should have been flagged as an error. Clearly
	   it does not contribute to the tranformation matrix. */
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
	else
	{
	  dstDom.p = WlzMakePlaneDomain(srcDom.p->type,
	  				srcDom.p->plane1, srcDom.p->lastpl,
					srcDom.p->line1, srcDom.p->lastln,
					srcDom.p->kol1, srcDom.p->lastkl,
					&errNum);
	  			 /* Need to fix the line column bounds later */
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
	    dstDom.p->plane1 += trans->tz;
	    dstDom.p->lastpl += trans->tz;
	  }
	  if( dstValues.vox ){
	    dstValues.vox->plane1 += trans->tz;
	    dstValues.vox->lastpl += trans->tz;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues,
			       NULL, NULL, &errNum);
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
* Function:	WlzAffineTransformVertexI				*
* Returns:	WlzIVertex2:		Transformed vertex.		*
* Purpose:	Transforms the given WlzIVertex2.			*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Affine transform to apply. 	*
*		WlzIVertex2 srcVtx:	Vertex to be transformed.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
************************************************************************/
WlzIVertex2	WlzAffineTransformVertexI(WlzAffineTransform *trans,
					  WlzIVertex2 srcVtx,
					  WlzErrorNum *dstErr)
{
  double	tD0;
  WlzIVertex2	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(trans->type != WLZ_TRANSFORM_2D_AFFINE)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    tD0 = (trans->mat[0][0] * srcVtx.vtX) +
    	  (trans->mat[0][1] * srcVtx.vtY) + trans->mat[0][2];
    dstVtx.vtX = WLZ_NINT(tD0);
    tD0= (trans->mat[1][0] * srcVtx.vtX) +
         (trans->mat[1][1] * srcVtx.vtY) + trans->mat[1][2];
    dstVtx.vtY = WLZ_NINT(tD0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVtx);
}

/************************************************************************
* Function:	WlzAffineTransformVertexF				*
* Returns:	WlzFVertex2:		Transformed vertex.		*
* Purpose:	Transforms the given WlzFVertex2.			*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Affine transform to apply. 	*
*		WlzFVertex2 srcVtx:	Vertex to be transformed.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
************************************************************************/
WlzFVertex2	WlzAffineTransformVertexF(WlzAffineTransform *trans,
					  WlzFVertex2 srcVtx,
					  WlzErrorNum *dstErr)
{
  WlzFVertex2	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(trans->type != WLZ_TRANSFORM_2D_AFFINE)
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
* Function:	WlzAffineTransformVertexD				*
* Returns:	WlzDVertex2:		Transformed vertex.		*
* Purpose:	Transforms the given WlzDVertex2.			*
* Global refs:	-							*
* Parameters:	WlzAffineTransform *trans: Affine transform to apply. 	*
*		WlzDVertex2 srcVtx:	Vertex to be transformed.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number, may be NULL.		*
************************************************************************/
WlzDVertex2	WlzAffineTransformVertexD(WlzAffineTransform *trans,
					  WlzDVertex2 srcVtx,
					  WlzErrorNum *dstErr)
{
  WlzDVertex2	dstVtx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(trans->type != WLZ_TRANSFORM_2D_AFFINE)
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
