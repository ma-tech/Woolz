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
* 22-12-00 bill Add normals.
* 13-12-00 bill Change members of WlzVertex and WlzVertexP.
************************************************************************/
#include <float.h>
#include <Wlz.h>

static WlzVertexP 		WlzVerticiesFromPoly2(
				  WlzPolygonDomain *poly,
				  WlzVertexP *dstNr,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzVerticiesFromBound(
				  WlzBoundList *bound,
				  WlzVertexP *dstNr,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzVerticiesFromCtr(
				  WlzContour *ctr,
				  WlzVertexP *dstNr,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzVerticiesFromGM2(
				  WlzGMModel *model,
				  WlzVertexP *dstNr,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzVerticiesFromGM3(
				  WlzGMModel *model,
				  WlzVertexP *dstNr,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static int			WlzVerticiesCntBound(
				  WlzBoundList *bound);
static WlzErrorNum 		WlzVerticiesCpBound(
				  WlzVertexP vData,
				  WlzDVertex2 *vNorm,
				  WlzVertexType vType,
				  int *off,
				  WlzBoundList *bound);
static WlzVertexP 		WlzVerticiesAlcPoly(
				  WlzObjectType polyType,
				  int cnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static void			WlzVerticiesNorm2(
				  WlzDVertex2 *nrm,
				  WlzVertexP vtx,
				  int cnt,
				  WlzObjectType type);
static WlzDVertex2 		WlzVerticiesNormPair2(
				  WlzDVertex2 *vtx);

/************************************************************************
* Function:	WlzVerticiesFromObj
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer which it fills with the verticies
*		from the given object. The object must be one of the
*		types that is represented by verticies, eg
*		WLZ_2D_POLYGON.
* Global refs:	-
* Parameters:	WlzObject *obj:		Given polygon domain object.
*		WlzVertexP *dstNr:	Destination ptr for normals.
*					The normals will always be
*					either WlzDVertex2 or
*					WlzDVertex3. May be NULL.
*		int *dstCnt:		Destination ptr for the number
*					of verticies. Can NOT be NULL.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies. Can NOT be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
WlzVertexP	WlzVerticiesFromObj(WlzObject *obj, WlzVertexP *dstNr,
				    int *dstCnt, WlzVertexType *dstType, 
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
	vData = WlzVerticiesFromPoly2(obj->domain.poly, dstNr,
				      dstCnt, dstType, &errNum);
	break;
      case WLZ_BOUNDLIST:
	vData = WlzVerticiesFromBound(obj->domain.b, dstNr,
				      dstCnt, dstType, &errNum);
	break;
      case WLZ_CONTOUR:
	vData = WlzVerticiesFromCtr(obj->domain.ctr, dstNr,
				    dstCnt, dstType, &errNum);
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
*		WlzVertexP *dstNr:	Destination ptr for normals,
*					may be NULL.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromPoly2(WlzPolygonDomain *poly,
				        WlzVertexP *dstNr, int *dstCnt,
					WlzVertexType *dstType,
					WlzErrorNum *dstErr)
{
  int		cnt;
  WlzVertexType	type;
  WlzVertexP	vData;
  WlzDVertex2	*vNorm = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  if(poly && ((cnt = poly->nvertices) > 0))
  {
    vData = WlzVerticiesAlcPoly(poly->type, cnt, &type, &errNum);
    if((errNum == WLZ_ERR_NONE) && dstNr)
    {
      if((vNorm = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
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
    if(dstNr)
    {
      WlzVerticiesNorm2(vNorm, vData, cnt, poly->type);
      (*dstNr).d2 = vNorm;
    }
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
*		WlzVertexP *dstNr:	Destination ptr for normals,
*					may be NULL.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromBound(WlzBoundList *bound,
				        WlzVertexP *dstNr, int *dstCnt,
					WlzVertexType *dstType,
					WlzErrorNum *dstErr)
{
  int		off,
  		cnt;
  WlzVertexType	type;
  WlzVertexP    vData;
  WlzDVertex2	*vNorm = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  if((cnt = WlzVerticiesCntBound(bound)) > 0)
  {
    vData = WlzVerticiesAlcPoly(bound->poly->type, cnt, &type, &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && dstNr)
  {
    if((vNorm = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && vData.v)
  {
    off = 0;
    errNum = WlzVerticiesCpBound(vData, vNorm, type, &off, bound);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstCnt = cnt;
    *dstType = type;
    if(dstNr)
    {
      (*dstNr).d2 = vNorm;
    }
  }
  else
  {
    if(vData.v)
    {
      AlcFree(vData.v);
      vData.v = NULL;
    }
    if(vNorm)
    {
      AlcFree(vNorm);
    }
  }
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
*		WlzVertexP *dstNr:	Destination ptr for normals,
*					may be NULL.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromCtr(WlzContour *ctr,
				      WlzVertexP *dstNr, int *dstCnt,
				      WlzVertexType *dstType,
				      WlzErrorNum *dstErr)
{
  WlzVertexP    vData;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  if(ctr && (ctr->model != NULL))
  {
    switch(ctr->model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	vData = WlzVerticiesFromGM2(ctr->model, dstNr, dstCnt, dstType,
				    &errNum);
        break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	vData = WlzVerticiesFromGM3(ctr->model, dstNr, dstCnt, dstType,
				    &errNum);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
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
* Function:	WlzVerticiesFromGM2
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer which it fills with the verticies
*		from a 2D GM.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Given model.
*		WlzVertexP *dstNr:	Destination ptr for normals,
*					may be NULL.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromGM2(WlzGMModel *model,
				      WlzVertexP *dstNr, int *dstCnt,
				      WlzVertexType *dstType,
				      WlzErrorNum *dstErr)
{
  int		idx,
  		cnt,
		nIdx,
		vIdx;
  WlzVertexType	type;
  WlzVertexP    vData;
  WlzGMVertex	*cV;
  WlzGMVertexT	*cVT;
  AlcVector	*vec;
  WlzDVertex2	*vNorm = NULL;
  WlzGMVertex	*nV[2];
  WlzVertexP	tVP[3];
  WlzDVertex2	nrmV[2],
  		segV[3];
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vData.v = NULL;
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
      break;
    case WLZ_GMMOD_2D:
      type = WLZ_VERTEX_D2;
      if((vData.v = AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      break;
  }
  if((errNum == WLZ_ERR_NONE) && dstNr)
  {
    if((vNorm = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idx = 0; idx < cnt; ++idx)
    {
      cV = (WlzGMVertex *)AlcVectorItemGet(vec, vIdx++);
      if(cV->idx >= 0)
      {
	if(model->type == WLZ_GMMOD_2I)
	{
	  *(vData.i2 + idx) = cV->geo.vg2I->vtx;
	}
	else /* model->type == WLZ_GMMOD_2D */
	{
	  *(vData.d2 + idx) = cV->geo.vg2D->vtx;
	}
	if(vNorm)
	{
	  cVT = cV->diskT->vertexT;
	  if((cVT == cVT->next) || (cVT->prev != cVT->next))
	  {
	    /* Vertex is either an isolated vertex or used by more than
	     * two edges. Normal undefined. */
	    (vNorm + idx)->vtX = 0.0;
	    (vNorm + idx)->vtY = 0.0;
	  }
	  else
	  {
	    /* Vertex is used by two edges. Find the other two verticies
	     * that are used by these two edges. */
	    nV[0] = cVT->prev->diskT->vertex;
	    nV[1] = cVT->next->diskT->vertex;
	    if(model->type == WLZ_GMMOD_2I)
	    { 
	      tVP[0].i2 = &(nV[0]->geo.vg2I->vtx);
	      tVP[1].i2 = &(cV->geo.vg2I->vtx);
	      tVP[2].i2 = &(nV[1]->geo.vg2I->vtx);
	      segV[0].vtX = tVP[0].i2->vtX;
	      segV[0].vtY = tVP[0].i2->vtY;
	      segV[1].vtX = tVP[1].i2->vtX;
	      segV[1].vtY = tVP[1].i2->vtY;
	      segV[2].vtX = tVP[2].i2->vtX;
	      segV[2].vtY = tVP[2].i2->vtY;
	    }
	    else /* model->type == WLZ_GMMOD_2D */
	    {
	      tVP[0].d2 = &(nV[0]->geo.vg2D->vtx);
	      tVP[1].d2 = &(cV->geo.vg2D->vtx);
	      tVP[2].d2 = &(nV[1]->geo.vg2D->vtx);
	      segV[0] = *(tVP[0].d2);
	      segV[1] = *(tVP[1].d2);
	      segV[2] = *(tVP[2].d2);
	    }
	    nrmV[0] = WlzVerticiesNormPair2(segV);
	    nrmV[1] = WlzVerticiesNormPair2(segV + 1);
	    (vNorm + idx)->vtX = (nrmV[0].vtX + nrmV[1].vtX) / 2.0;
	    (vNorm + idx)->vtY = (nrmV[0].vtY + nrmV[1].vtY) / 2.0;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstCnt = cnt;
    *dstType = type;
    if(dstNr)
    {
      (*dstNr).d2 = vNorm;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/************************************************************************
* Function:	WlzVerticiesFromGM3
* Returns:	WlzVertexP:		Allocated verticies.
* Purpose:	Allocates a buffer which it fills with the verticies
*		from a 3D GM.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Given model.
*		WlzVertexP *dstNr:	Destination ptr for normals,
*					may be NULL.
*		int *dstCnt:		Destination ptr for the number
*					of verticies.
*		WlzVertexType *dstType:	Destination ptr for the type
*					of verticies.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzVertexP WlzVerticiesFromGM3(WlzGMModel *model,
				      WlzVertexP *dstNr, int *dstCnt,
				      WlzVertexType *dstType,
				      WlzErrorNum *dstErr)
{
  int		vIdx,
		vecIdx,
		sIdx,
		sCnt,
		sMax = 0,
  		vCnt,
		manifold;
  double	tD0;
  AlcVector	*vec;
  WlzVertexType	type;
  WlzVertexP    vData;
  WlzGMVertex	*cV;
  WlzGMVertex	**sVBuf = NULL;
  WlzGMVertexT	*vT0,
  		*vT1;
  WlzGMEdgeT	*eT1;
  WlzDVertex3	*vNorm = NULL;
  WlzDVertex3	sVG[3];
  WlzDVertex3	cNrm,
  		sNrm;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  vData.v = NULL;
  vecIdx = 0;
  vCnt = model->res.vertex.numElm;
  vec = model->res.vertex.vec;
  if(model->type == WLZ_GMMOD_3I)
  {
    type = WLZ_VERTEX_I3;
    if((vData.v = AlcMalloc(sizeof(WlzIVertex3) * vCnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  else /* model->type == WLZ_GMMOD_3D */
  {
    type = WLZ_VERTEX_D3;
    if((vData.v = AlcMalloc(sizeof(WlzDVertex3) * vCnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstNr)
  {
    if((vNorm = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * vCnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(vIdx = 0; vIdx < vCnt; ++vIdx)
    {
      manifold = 1;
      cV = (WlzGMVertex *)AlcVectorItemGet(vec, vecIdx++);
      if(cV->idx >= 0)
      {
	if(model->type == WLZ_GMMOD_3I)
	{
	  *(vData.i3 + vIdx) = cV->geo.vg3I->vtx;
	}
	else /* model->type == WLZ_GMMOD_3D */
	{
	  *(vData.d3 + vIdx) = cV->geo.vg3D->vtx;
	}
	if(dstNr)
	{
	  *(vNorm + vIdx) = WlzGMVertexNormal3D(model, cV, &sMax, &sVBuf,
	  					&errNum);
	}
      }
    }
  }
  if(sVBuf)
  {
    AlcFree(sVBuf);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstCnt = vCnt;
    *dstType = type;
    if(dstNr)
    {
      (*dstNr).d3 = vNorm;
    }
  }
  else
  {
    if(vNorm)
    {
      AlcFree(vNorm);
    }
  }
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
*		WlzDVertex2 *vNorm:	Given buffer for normals, may
*					be NULL.
*		WlzVertexType vType:	Type of verticies.
*		int *off:		Ptr to offset into buffer.
*		WlzBoundList *bound:	Given boundary domain.
************************************************************************/
static WlzErrorNum WlzVerticiesCpBound(WlzVertexP vData, WlzDVertex2 *vNorm,
				       WlzVertexType vType,
				       int *off, WlzBoundList *bound)
{
  int		cnt;
  WlzPolygonDomain *poly;
  WlzVertexP	vPtr;
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
	  vPtr.i2 = vData.i2 + *off;
	  WlzValueCopyIVertexToIVertex(vPtr.i2,
				       (WlzIVertex2 *)(poly->vtx), cnt);
	}
	break;
      case WLZ_POLYGON_FLOAT:
	if(vType != WLZ_VERTEX_F2)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  vPtr.f2 = vData.f2 + *off;
	  WlzValueCopyFVertexToFVertex(vPtr.f2,
				       (WlzFVertex2 *)(poly->vtx), cnt);
	}
	break;
      case WLZ_POLYGON_DOUBLE:
	if(vType != WLZ_VERTEX_D2)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  vPtr.d2 = vData.d2 + *off;
	  WlzValueCopyDVertexToDVertex(vPtr.d2,
				       (WlzDVertex2 *)(poly->vtx), cnt);
	}
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(vNorm)
      {
	WlzVerticiesNorm2(vNorm + *off, vPtr, cnt, poly->type);
      }
      *off += cnt;
    }
  }
  if((errNum == WLZ_ERR_NONE) && bound->next)
  {
    errNum = WlzVerticiesCpBound(vData, vNorm, vType, off, bound->next);
  }
  if((errNum == WLZ_ERR_NONE) && bound->down)
  {
    errNum = WlzVerticiesCpBound(vData, vNorm, vType, off, bound->down);
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

/************************************************************************
* Function:	WlzVerticiesNorm2
* Returns:	void
* Purpose:	Computes the normals of the given verticies which are
*		assumed to lie in a 2D polygon.
*		The normals all have +ve x components. 
*		There is no need for any type checking in this function
*		because it is static and all types have been checked.
* Global refs:	-
* Parameters:	WlzDVertex2 *nrm:	Buffer for the normals.
*		WlzVertexP vtx:		The given verticies.
*		int cnt:		The number of verticies (and
*					normals).
*		WlzObjectType pType:	Polygon type.
************************************************************************/
static void	WlzVerticiesNorm2(WlzDVertex2 *nrm, WlzVertexP vtx, int cnt,
				  WlzObjectType pType)
{
  int		idx,
  		idx1;
  WlzVertexType vType;
  WlzDVertex2	nrmV[2],
  		segV[2];

  switch(cnt)
  {
    case 1:
      /* Normal doesn't have a meaning, set it to (0,0). */
      nrm->vtX = nrm->vtY = 0.0;
      break;
    case 2:
      /* Find normal to the line segment. */
      switch(pType)
      {
	case WLZ_POLYGON_INT:
	  segV[0].vtX = (vtx.i2 + 0)->vtX;
	  segV[0].vtY = (vtx.i2 + 0)->vtY;
	  segV[1].vtX = (vtx.i2 + 1)->vtX;
	  segV[1].vtY = (vtx.i2 + 1)->vtY;
	  break;
	case WLZ_POLYGON_FLOAT:
	  segV[0].vtX = (vtx.f2 + 0)->vtX;
	  segV[0].vtY = (vtx.f2 + 0)->vtY;
	  segV[1].vtX = (vtx.f2 + 1)->vtX;
	  segV[1].vtY = (vtx.f2 + 1)->vtY;
	  break;
	case WLZ_POLYGON_DOUBLE:
	  segV[0] = *(vtx.d2 + 0);
	  segV[1] = *(vtx.d2 + 1);
	  break;
      }
      *nrm = WlzVerticiesNormPair2(segV);
      break;
    default:
      /* There are more than two verticies.
       * Normals are computed for a vertex by averaging the normals
       * of the line segments previous to and after the vertex. The
       * resulting normal's length is 1.0. */
      idx = 0;
      idx1 = 1;
      switch(pType)
      {
	case WLZ_POLYGON_INT:
	  segV[0].vtX = (vtx.i2 + cnt - 1)->vtX; 
	  segV[0].vtY = (vtx.i2 + cnt - 1)->vtY; 
	  segV[1].vtX = (vtx.i2 + 0)->vtX;
	  segV[1].vtY = (vtx.i2 + 0)->vtY;
	  break;
	case WLZ_POLYGON_FLOAT:
	  segV[0].vtX = (vtx.f2 + cnt - 1)->vtX; 
	  segV[0].vtY = (vtx.f2 + cnt - 1)->vtY; 
	  segV[1].vtX = (vtx.f2 + 0)->vtX;
	  segV[1].vtY = (vtx.f2 + 0)->vtY;
	  break;
	case WLZ_POLYGON_DOUBLE:
	  segV[0] = *(vtx.d2 + cnt - 1);
	  segV[1] = *(vtx.d2 + 0);
	  break;
      }
      nrmV[1] = WlzVerticiesNormPair2(segV);
      while(idx1 < cnt)
      {
        segV[0] = segV[1];
	nrmV[0] = nrmV[1];
	switch(pType)
	{
	  case WLZ_POLYGON_INT:
	    segV[1].vtX = (vtx.i2 + idx1)->vtX;
	    segV[1].vtY = (vtx.i2 + idx1)->vtY;
	    break;
	  case WLZ_POLYGON_FLOAT:
	    segV[1].vtX = (vtx.f2 + idx1)->vtX;
	    segV[1].vtY = (vtx.f2 + idx1)->vtY;
	    break;
	  case WLZ_POLYGON_DOUBLE:
	    segV[1] = *(vtx.d2 + idx1);
	    break;
	}
	nrmV[1] = WlzVerticiesNormPair2(segV);
	(nrm + idx)->vtX = (nrmV[0].vtX +  nrmV[1].vtX) / 2.0;
	(nrm + idx)->vtY = (nrmV[0].vtY +  nrmV[1].vtY) / 2.0;
        idx = idx1++;
      }
      segV[0] = segV[1];
      nrmV[0] = nrmV[1];
      switch(pType)
      {
	case WLZ_POLYGON_INT:
	  segV[1].vtX = (vtx.i2 + 0)->vtX;
	  segV[1].vtY = (vtx.i2 + 0)->vtY;
	  break;
	case WLZ_POLYGON_FLOAT:
	  segV[1].vtX = (vtx.f2 + 0)->vtX;
	  segV[1].vtY = (vtx.f2 + 0)->vtY;
	  break;
	case WLZ_POLYGON_DOUBLE:
	  segV[1] = *(vtx.d2 + 0);
	  break;
      }
      nrmV[1] = WlzVerticiesNormPair2(segV);
      (nrm + idx)->vtX = (nrmV[0].vtX +  nrmV[1].vtX) / 2.0;
      (nrm + idx)->vtY = (nrmV[0].vtY +  nrmV[1].vtY) / 2.0;
      break;
  }
}

/************************************************************************
* Function:	WlzVerticiesNormPair2
* Returns:	void
* Purpose:	Computes the normal (n) to a segment (g) between the
*		given pair of verticies. There are clearly two solutions
*		to the problem of finding a normal to a line segment,
*		but this function always finds the normal vector with a
*		+ve x component.
*		If the two verticies are coincident then the normal
*		vector is set to {0, 0}.
*		With two non-coincident verticies the normal vector is
*		computed using the relationships g.n = 0 and |n|^2 = 1.
*		Giving nx = 1/sqrt(1 + (gx/gy)^2), ny = -nx gx/gy.
*		There is no need for any type checking in this function
*		because it is static and all types have been checked.
* Global refs:	-
* Parameters:	WlzDVertex2 *vtx:	The given pair of verticies.
************************************************************************/
static WlzDVertex2 WlzVerticiesNormPair2(WlzDVertex2 *vtx)
{
  WlzDVertex2	tV0,
  		tV1,
		nrm;

  WLZ_VTX_2_SUB(tV0, *(vtx + 1), *(vtx + 0));
  tV1.vtX = tV0.vtX * tV0.vtX;
  tV1.vtY = tV0.vtY * tV0.vtY; 
  if(tV1.vtY < DBL_EPSILON)
  {
    nrm.vtX = 0.0;
    if(tV1.vtX < DBL_EPSILON)
    {
      nrm.vtY = 0.0;
    }
    else
    {
      nrm.vtY = 1.0;
    }
  }
  else if(tV1.vtX < DBL_EPSILON)
  {
    nrm.vtX = 1.0;
    nrm.vtY = 0.0;
  }
  else
  {
    nrm.vtX = 1.0 / sqrt(1.0 + (tV1.vtX / tV1.vtY));
    nrm.vtY = -((tV0.vtX) * nrm.vtX) / tV0.vtY;
  }
  return(nrm);
}
