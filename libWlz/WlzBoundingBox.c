#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzBoundingBox.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for computing the bounding box of woolz
*		objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 15-08-00 bill	Add WlzBoundingBoxContour(). Remove obsolete types:
*		WLZ_VECTOR_(INT)|(FLOAT) and WLZ_VECTOR_(INT)|FLOAT).
************************************************************************/
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <Wlz.h>

static WlzIBox3			WlzBoundingBoxBound3D(
				  WlzBoundList *,
				  WlzErrorNum *);
static WlzIBox3			WlzBoundingBoxPoly3D(
				  WlzPolygonDomain *,
				  WlzErrorNum *);
static WlzIBox3			WlzBoundingBoxTransObj3D(
				  WlzObject *,
				  WlzAffineTransform *,
				  WlzErrorNum *);
static WlzIBox3 		WlzBoundingBoxContour(
				  WlzContour *ctr,
				  WlzErrorNum *dstErr);

/************************************************************************
* Function:	WlzBoundingBox2D					*
* Returns:	WlzIBox2 *:		2D bounding box.		*
* Purpose:	Computes the 2D bounding box of the given object.	*
* Global refs:	-							*
* Parameters:	WlzObject *inObj:	The given object.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
WlzIBox2	WlzBoundingBox2D(WlzObject *inObj, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIBox2	bBox2D;
  WlzIBox3	bBox3D;

  bBox3D = WlzBoundingBox3D(inObj, &errNum);
  bBox2D.xMin = bBox3D.xMin;
  bBox2D.yMin = bBox3D.yMin;
  bBox2D.xMax = bBox3D.xMax;
  bBox2D.yMax = bBox3D.yMax;
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox2D);
}

/************************************************************************
* Function:	WlzBoundingBox3D					*
* Returns:	WlzIBox3 *:		3D bounding box.		*
* Purpose:	Computes the 3D bounding box of the given object.	*
* Global refs:	-							*
* Parameters:	WlzObject *inObj:	The given object.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
WlzIBox3	WlzBoundingBox3D(WlzObject *inObj, WlzErrorNum *dstErr)
{
  WlzIBox3	bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D.xMin = 0;
  bBox3D.xMax = 0;
  bBox3D.yMin = 0;
  bBox3D.yMax = 0;
  bBox3D.zMin = 0;
  bBox3D.zMax = 0;
  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(inObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(inObj->domain.core->type)
	  {
	    case WLZ_INTERVALDOMAIN_INTVL: /* FALLTHROUGH */
            case WLZ_INTERVALDOMAIN_RECT:
	      bBox3D.xMin = inObj->domain.i->kol1;
	      bBox3D.yMin = inObj->domain.i->line1;
	      bBox3D.xMax = inObj->domain.i->lastkl;
	      bBox3D.yMax = inObj->domain.i->lastln;
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	}
        break;
      case WLZ_3D_DOMAINOBJ:
        if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  if(inObj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	  {
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	  }
	  else
	  {
	    bBox3D.xMin = inObj->domain.p->kol1;
	    bBox3D.yMin = inObj->domain.p->line1;
	    bBox3D.zMin = inObj->domain.p->plane1;
	    bBox3D.xMax = inObj->domain.p->lastkl;
	    bBox3D.yMax = inObj->domain.p->lastln;
	    bBox3D.zMax = inObj->domain.p->lastpl;
	  }
	}
	break;
      case WLZ_TRANS_OBJ:
        if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(inObj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else if(inObj->domain.core->type != WLZ_TRANSFORM_2D_AFFINE)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  bBox3D = WlzBoundingBoxTransObj3D(inObj->values.obj,
	  				    inObj->domain.t, &errNum);
	}
	break;
      case WLZ_2D_POLYGON:
        if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
          bBox3D = WlzBoundingBoxPoly3D(inObj->domain.poly, &errNum);
	}
        break;
      case WLZ_BOUNDLIST:
        if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
          bBox3D = WlzBoundingBoxBound3D(inObj->domain.b, &errNum);
	}
	break;
      case WLZ_CONTOUR:
	if(inObj->domain.ctr == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  bBox3D = WlzBoundingBoxContour(inObj->domain.ctr, &errNum);
	}
        break;
      case WLZ_EMPTY_OBJ:
      case WLZ_AFFINE_TRANS:
      case WLZ_HISTOGRAM:
      case WLZ_PROPERTY_OBJ:
      case WLZ_CONV_HULL:
      case WLZ_3D_WARP_TRANS:
      case WLZ_3D_POLYGON:
      case WLZ_RECTANGLE:
      case WLZ_CONVOLVE_INT:
      case WLZ_CONVOLVE_FLOAT:
      case WLZ_WARP_TRANS:
      case WLZ_FMATCHOBJ:
      case WLZ_TEXT:
      case WLZ_COMPOUND_ARR_1:
      case WLZ_COMPOUND_ARR_2:
      case WLZ_COMPOUND_LIST_1:
      case WLZ_COMPOUND_LIST_2:
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/************************************************************************
* Function:	WlzBoundingBoxTransObj3D				*
* Returns:	WlzIBox3 *:		3D bounding box.		*
* Purpose:	Computes the 3D bounding box of the transformed object.	*
* Global refs:	-							*
* Parameters:	WlzObject *inObj:	Input object.			*
*		WlzAffineTransform *trans: Given transform.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzIBox3	WlzBoundingBoxTransObj3D(WlzObject *inObj,
					 WlzAffineTransform *trans,
				         WlzErrorNum *dstErr)
{
  int		tI0,
  		idx;
  WlzIVertex2	bufVx[4];
  WlzIBox3	tBox3D,
  		bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D.xMin = 0;
  bBox3D.xMax = 0;
  bBox3D.yMin = 0;
  bBox3D.yMax = 0;
  bBox3D.zMin = 0;
  bBox3D.zMax = 0;
  tBox3D = WlzBoundingBox3D(inObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    bufVx[0].vtX = tBox3D.xMin;
    bufVx[0].vtY = tBox3D.yMin;
    bufVx[1].vtX = tBox3D.xMax;
    bufVx[1].vtY = tBox3D.yMin;
    bufVx[2].vtX = tBox3D.xMax;
    bufVx[2].vtY = tBox3D.yMax;
    bufVx[3].vtX = tBox3D.xMin;
    bufVx[3].vtY = tBox3D.yMax;
    idx = 0;
    while((idx < 4) && (errNum == WLZ_ERR_NONE))
    {
      bufVx[idx] = WlzAffineTransformVertexI2(trans, bufVx[idx], &errNum);
      ++idx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox3D.xMax = bBox3D.xMin = bufVx[0].vtX;
    bBox3D.yMax = bBox3D.yMin = bufVx[0].vtY;
    bBox3D.zMin = tBox3D.zMin;
    bBox3D.zMax = tBox3D.zMax;
    idx = 1;
    while(idx < 4)
    {
      if((tI0 = bufVx[idx].vtX) < bBox3D.xMin)
      {
	bBox3D.xMin = tI0;
      }
      else if(tI0 > bBox3D.xMax)
      {
	bBox3D.xMax = tI0;
      }
      if((tI0 = bufVx[idx].vtY) < bBox3D.yMin)
      {
	bBox3D.yMin = tI0;
      }
      else if(tI0 > bBox3D.yMax)
      {
	bBox3D.yMax = tI0;
      }
      ++idx;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/************************************************************************
* Function:	WlzBoundingBoxPoly3D					*
* Returns:	WlzIBox3 *:		3D bounding box.		*
* Purpose:	Computes the 3D bounding box of the polygon domain.	*
* Global refs:	-							*
* Parameters:	WlzPolygonDomain *poly:	The given polygon domain.	*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzIBox3	WlzBoundingBoxPoly3D(WlzPolygonDomain *poly,
				     WlzErrorNum *dstErr)
{
  int		tI0,
  		idx;
  WlzIVertex2	*tIVxP;
  WlzFVertex2	*tFVxP;
  WlzDVertex2	*tDVxP;
  WlzIBox3	bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D.xMin = 0;
  bBox3D.xMax = 0;
  bBox3D.yMin = 0;
  bBox3D.yMax = 0;
  bBox3D.zMin = 0;
  bBox3D.zMax = 0;
  switch(poly->type)
  {
    case WLZ_POLYGON_INT:
      if((idx = poly->nvertices) > 0)
      {
	tIVxP = poly->vtx;
	bBox3D.xMax = bBox3D.xMin = tIVxP->vtX;
	bBox3D.yMax = bBox3D.yMin = tIVxP->vtY;
	while(--idx > 0)
	{
	  if((tI0 = tIVxP->vtX) < bBox3D.xMin)
	  {
	    bBox3D.xMin = tI0;
	  }
	  else if(tI0 > bBox3D.xMax)
	  {
	    bBox3D.xMax = tI0;
	  }
	  if((tI0 = tIVxP->vtY) < bBox3D.yMin)
	  {
	    bBox3D.yMin = tI0;
	  }
	  else if(tI0 > bBox3D.yMax)
	  {
	    bBox3D.yMax = tI0;
	  }
	  ++tIVxP;
	}
      }
      break;
    case WLZ_POLYGON_FLOAT:
      if((idx = poly->nvertices) > 0)
      {
	tFVxP = (WlzFVertex2 *)(poly->vtx);
	bBox3D.xMax = bBox3D.xMin = WLZ_NINT(tIVxP->vtX);
	bBox3D.yMax = bBox3D.yMin = WLZ_NINT(tIVxP->vtY);
	while(--idx > 0)
	{
	  if((tI0 = WLZ_NINT(tFVxP->vtX)) < bBox3D.xMin)
	  {
	    bBox3D.xMin = tI0;
	  }
	  else if(tI0 > bBox3D.xMax)
	  {
	    bBox3D.xMax = tI0;
	  }
	  if((tI0 = WLZ_NINT(tFVxP->vtY)) < bBox3D.yMin)
	  {
	    bBox3D.yMin = tI0;
	  }
	  else if(tI0 > bBox3D.yMax)
	  {
	    bBox3D.yMax = tI0;
	  }
	  ++tFVxP;
	}
      }
      break;
    case WLZ_POLYGON_DOUBLE:
      if((idx = poly->nvertices) > 0)
      {
	tDVxP = (WlzDVertex2 *)(poly->vtx);
	bBox3D.xMax = bBox3D.xMin = WLZ_NINT(tIVxP->vtX);
	bBox3D.yMax = bBox3D.yMin = WLZ_NINT(tIVxP->vtY);
	while(--idx > 0)
	{
	  if((tI0 = WLZ_NINT(tDVxP->vtX)) < bBox3D.xMin)
	  {
	    bBox3D.xMin = tI0;
	  }
	  else if(tI0 > bBox3D.xMax)
	  {
	    bBox3D.xMax = tI0;
	  }
	  if((tI0 = WLZ_NINT(tDVxP->vtY)) < bBox3D.yMin)
	  {
	    bBox3D.yMin = tI0;
	  }
	  else if(tI0 > bBox3D.yMax)
	  {
	    bBox3D.yMax = tI0;
	  }
	  ++tDVxP;
	}
      }
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/************************************************************************
* Function:	WlzBoundingBoxBound3D					*
* Returns:	WlzIBox3 *:		3D bounding box.		*
* Purpose:	Computes the 3D bounding box of the boundary list.	*
* Global refs:	-							*
* Parameters:	WlzBoundList *bound:	The given boundary list.	*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzIBox3	WlzBoundingBoxBound3D(WlzBoundList *bound,
				      WlzErrorNum *dstErr)
{
  WlzIBox3	tBox3D,
  		bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D = WlzBoundingBoxPoly3D(bound->poly, &errNum);
  if(bound->next && (errNum == WLZ_ERR_NONE))
  {
    tBox3D = WlzBoundingBoxBound3D(bound->next, &errNum);
    if(tBox3D.xMin < bBox3D.xMin)
    {
      bBox3D.xMin = tBox3D.xMin;
    }
    if(tBox3D.yMin < bBox3D.yMin)
    {
      bBox3D.yMin = tBox3D.yMin;
    }
    if(tBox3D.zMin < bBox3D.zMin)
    {
      bBox3D.zMin = tBox3D.zMin;
    }
    if(tBox3D.xMax < bBox3D.xMax)
    {
      bBox3D.xMax = tBox3D.xMax;
    }
    if(tBox3D.yMax < bBox3D.yMax)
    {
      bBox3D.yMax = tBox3D.yMax;
    }
    if(tBox3D.zMax < bBox3D.zMax)
    {
      bBox3D.zMax = tBox3D.zMax;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/************************************************************************
* Function:	WlzBoundingBoxContour					*
* Returns:	WlzIBox3 *:		3D bounding box.		*
* Purpose:	Computes the 3D bounding box of the contour.		*
* Global refs:	-							*
* Parameters:	WlzContour *ctr:	The given contour.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzIBox3 WlzBoundingBoxContour(WlzContour *ctr, WlzErrorNum *dstErr)
{
  WlzGMShell	*shell0,
  		*shell1;
  WlzDBox3	sBox3D,
  		bBox3D;
  WlzIBox3	bBox3I;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr->model == NULL)
  {
    errNum =  WLZ_ERR_DOMAIN_NULL;
  }
  else if((shell0 = ctr->model->child) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    /* Find the bounding box of all the GM's shell's bounding boxes. */
    shell1 = shell0->next;
    errNum = WlzGMShellGetGBB3D(shell0, &bBox3D);
    while((errNum == WLZ_ERR_NONE) && (shell1->idx != shell0->idx))
    {
      if((errNum = WlzGMShellGetGBB3D(shell0, &sBox3D)) == WLZ_ERR_NONE)
      {
	if(sBox3D.xMin < bBox3D.xMin)
	{
	  bBox3D.xMin = sBox3D.xMin;
	}
	if(sBox3D.yMin < bBox3D.yMin)
	{
	  bBox3D.yMin = sBox3D.yMin;
	}
	if(sBox3D.zMin < bBox3D.zMin)
	{
	  bBox3D.zMin = sBox3D.zMin;
	}
	if(sBox3D.xMax < bBox3D.xMax)
	{
	  bBox3D.xMax = sBox3D.xMax;
	}
	if(sBox3D.yMax < bBox3D.yMax)
	{
	  bBox3D.yMax = sBox3D.yMax;
	}
	if(sBox3D.zMax < bBox3D.zMax)
	{
	  bBox3D.zMax = sBox3D.zMax;
	}
        shell1 = shell1->next;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox3I.xMin = (int )floor(bBox3D.xMin);
    bBox3I.yMin = (int )floor(bBox3D.yMin);
    bBox3I.zMin = (int )floor(bBox3D.zMin);
    bBox3I.xMax = (int )ceil(bBox3D.xMax);
    bBox3I.yMax = (int )ceil(bBox3D.yMax);
    bBox3I.zMax = (int )ceil(bBox3D.zMax);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3I);
}
