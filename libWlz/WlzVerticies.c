#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        WlzVerticies.c
* Date:         November 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for extracting verticies from objects
*		represented by verticies, eg WLZ_2D_POLYGON.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 13-12-00 bill Change members of WlzVertex and WlzVertexP.
************************************************************************/
#include <Wlz.h>

static WlzVertexP 		WlzVerticiesFromPoly2(
				  WlzPolygonDomain *poly,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzVerticiesFromBound(
				  WlzBoundList *bound,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzVerticiesFromCtr(
				  WlzContour *ctr,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static int			WlzVerticiesCntBound(
				  WlzBoundList *bound);
static WlzErrorNum 		WlzVerticiesCpBound(
				  WlzVertexP vData,
				  WlzVertexType vType,
				  int *off,
				  WlzBoundList *bound);
static WlzVertexP 		WlzVerticiesAlcPoly(
				  WlzObjectType polyType,
				  int cnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);

/************************************************************************
* Function:	WlzVerticiesFromObj
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer which it fills with the verticies
*		from the given object. The object must be one of the
*		types that is represented by verticies, eg
*		WLZ_2D_POLYGON.
* Global refs:	-
* Parameters:	WlzObject *obj:		Given polygon domain object.
*		int *dstCnt:		Destination ptr for the number
*					of verticies. Can NOT be NULL.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies. Can NOT be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
WlzVertexP	WlzVerticiesFromObj(WlzObject *obj, int *dstCnt,
				    WlzVertexType *dstType, 
				    WlzErrorNum *dstErr)
{
  WlzVertexP	vData;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_POLYGON:
	vData = WlzVerticiesFromPoly2(obj->domain.poly,
				      dstCnt, dstType, &errNum);
	break;
      case WLZ_BOUNDLIST:
	vData = WlzVerticiesFromBound(obj->domain.b, dstCnt,
				      dstType, &errNum);
	break;
      case WLZ_CONTOUR:
	vData = WlzVerticiesFromCtr(obj->domain.ctr, dstCnt,
				    dstType, &errNum);
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/************************************************************************
* Function:	WlzVerticiesFromPoly2
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer which it fills with the verticies
*		from a 3D polygon domain.
* Global refs:	-
* Parameters:	WlzPolygonDomain *poly:	Given polygon domain.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromPoly2(WlzPolygonDomain *poly, int *dstCnt,
					WlzVertexType *dstType,
					WlzErrorNum *dstErr)
{
  int		cnt;
  WlzVertexType	type;
  WlzVertexP	vData;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  if(poly && ((cnt = poly->nvertices) > 0))
  {
    vData = WlzVerticiesAlcPoly(poly->type, cnt, &type, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(poly->type)
      {
	case WLZ_POLYGON_INT:
	  type = WLZ_VERTEX_I2;
	  WlzValueCopyIVertexToIVertex(vData.i2,
				       (WlzIVertex2 *)(poly->vtx), cnt);
	  break;
	case WLZ_POLYGON_FLOAT:
	  type = WLZ_VERTEX_F2;
	  WlzValueCopyFVertexToFVertex(vData.f2,
				       (WlzFVertex2 *)(poly->vtx), cnt);
	  break;
	case WLZ_POLYGON_DOUBLE:
	  type = WLZ_VERTEX_D2;
	  WlzValueCopyDVertexToDVertex(vData.d2,
				       (WlzDVertex2 *)(poly->vtx), cnt);
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstType = type;
    *dstCnt = cnt;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/************************************************************************
* Function:	WlzVerticiesFromBound
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer which it fills with the verticies
*		from a boundary domain.
* Global refs:	-
* Parameters:	WlzBoundList *bound:	Given boundary domain.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromBound(WlzBoundList *bound, int *dstCnt,
					WlzVertexType *dstType,
					WlzErrorNum *dstErr)
{
  int		off,
  		cnt;
  WlzVertexType	type;
  WlzVertexP    vData;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  if((cnt = WlzVerticiesCntBound(bound)) > 0)
  {
    vData = WlzVerticiesAlcPoly(bound->poly->type, cnt, &type, &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && vData.v)
  {
    off = 0;
    errNum = WlzVerticiesCpBound(vData, type, &off, bound);
  }
  *dstCnt = cnt;
  *dstType = type;
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/************************************************************************
* Function:	WlzVerticiesFromCtr
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer which it fills with the verticies
*		from a convex hull.
* Global refs:	-
* Parameters:	WlzContour *ctr:	Given contour.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromCtr(WlzContour *ctr, int *dstCnt,
				       WlzVertexType *dstType,
				       WlzErrorNum *dstErr)
{
  int		idx,
  		cnt,
		vIdx;
  WlzVertexType	type;
  WlzVertexP    vData;
  WlzGMVertex	*cV;
  AlcVector	*vec;
  WlzGMModel	*model;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  if(ctr && ((model = ctr->model) != NULL))
  {
    vIdx = 0;
    cnt = model->res.vertex.numElm;
    vec = model->res.vertex.vec;
    switch(model->type)
    {
      case WLZ_GMMOD_2I:
	type = WLZ_VERTEX_I2;
	if((vData.v = AlcMalloc(sizeof(WlzIVertex2) * cnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idx = 0; idx < cnt; ++idx)
	  {
	    cV = (WlzGMVertex *)AlcVectorItemGet(vec, vIdx++);
	    if(cV->idx >= 0)
	    {
	      *(vData.i2 + idx) = cV->geo.vg2I->vtx;
	    }
	  }
	}
        break;
      case WLZ_GMMOD_2D:
	type = WLZ_VERTEX_D2;
	if((vData.v = AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idx = 0; idx < cnt; ++idx)
	  {
	    cV = (WlzGMVertex *)AlcVectorItemGet(vec, vIdx++);
	    if(cV->idx >= 0)
	    {
	      *(vData.d2 + idx) = cV->geo.vg2D->vtx;
	    }
	  }
	}
        break;
      case WLZ_GMMOD_3I:
	type = WLZ_VERTEX_I3;
	if((vData.v = AlcMalloc(sizeof(WlzIVertex3) * cnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idx = 0; idx < cnt; ++idx)
	  {
	    cV = (WlzGMVertex *)AlcVectorItemGet(vec, vIdx++);
	    if(cV->idx >= 0)
	    {
	      *(vData.i3 + idx) = cV->geo.vg3I->vtx;
	    }
	  }
	}
        break;
      case WLZ_GMMOD_3D:
	type = WLZ_VERTEX_D3;
	if((vData.v = AlcMalloc(sizeof(WlzDVertex3) * cnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idx = 0; idx < cnt; ++idx)
	  {
	    cV = (WlzGMVertex *)AlcVectorItemGet(vec, vIdx++);
	    if(cV->idx >= 0)
	    {
	      *(vData.d3 + idx) = cV->geo.vg3D->vtx;
	    }
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  *dstCnt = cnt;
  *dstType = type;
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/************************************************************************
* Function:	WlzVerticiesCntBound
* Returns:	int:			Number of verticies in boundary.
* Purpose:	Counts the number of verticies in the polygon domains
*		of the boundary. This is a recursive function.
* Global refs:	-
* Parameters:	WlzBoundList *bound:	Given boundary domain.
************************************************************************/
static int	WlzVerticiesCntBound(WlzBoundList *bound)
{
  int		cnt;

  cnt = bound->poly->nvertices + WlzVerticiesCntBound(bound->next) +
  	WlzVerticiesCntBound(bound->down);
  return(cnt);
}

/************************************************************************
* Function:	WlzVerticiesCpBound
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Copies verticies from the boundaries polygon domain
*		to the buffer.
* Global refs:	-
* Parameters:	WlzVertexP vData:	Given buffer.
*		WlzVertexType vType:	Type of verticies.
*		int *off:		Ptr to offset into buffer.
*		WlzBoundList *bound:	Given boundary domain.
************************************************************************/
static WlzErrorNum WlzVerticiesCpBound(WlzVertexP vData,
				       WlzVertexType vType,
				       int *off, WlzBoundList *bound)
{
  int		cnt;
  WlzPolygonDomain *poly;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((poly = bound->poly) != NULL) && ((cnt = poly->nvertices) > 0))
  {
    switch(poly->type)
    {
      case WLZ_POLYGON_INT:
	if(vType != WLZ_VERTEX_I2)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  WlzValueCopyIVertexToIVertex(vData.i2 + *off,
				       (WlzIVertex2 *)(poly->vtx), cnt);
	  *off += cnt;
	}
	break;
      case WLZ_POLYGON_FLOAT:
	if(vType != WLZ_VERTEX_F2)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  WlzValueCopyFVertexToFVertex(vData.f2 + *off,
				       (WlzFVertex2 *)(poly->vtx), cnt);
	  *off += cnt;
	}
	break;
      case WLZ_POLYGON_DOUBLE:
	if(vType != WLZ_VERTEX_D2)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  WlzValueCopyDVertexToDVertex(vData.d2 + *off,
				       (WlzDVertex2 *)(poly->vtx), cnt);
	  *off += cnt;
	}
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && bound->next)
  {
    errNum = WlzVerticiesCpBound(vData, vType, off, bound->next);
  }
  if((errNum == WLZ_ERR_NONE) && bound->down)
  {
    errNum = WlzVerticiesCpBound(vData, vType, off, bound->down);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzVerticiesAlcPoly
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer for copting the verticies of a
*		polygon domain.
* Global refs:	-
* Parameters:	WlzObjectType polyType:	Type of polygon domain.
*		int cnt:		Number of verticies to allocate
*					room for.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesAlcPoly(WlzObjectType polyType, int cnt,
				      WlzVertexType *dstType,
				      WlzErrorNum *dstErr)
{
  int		vSize;
  WlzVertexType	type;
  WlzVertexP	vData;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  switch(polyType)
  {
    case WLZ_POLYGON_INT:
      type = WLZ_VERTEX_I2;
      vSize = sizeof(WlzIVertex2);
      break;
    case WLZ_POLYGON_FLOAT:
      type = WLZ_VERTEX_F2;
      vSize = sizeof(WlzFVertex2);
      break;
    case WLZ_POLYGON_DOUBLE:
      type = WLZ_VERTEX_D2;
      vSize = sizeof(WlzDVertex2);
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((vData.v = AlcMalloc(vSize * cnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      *dstType = type;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}
