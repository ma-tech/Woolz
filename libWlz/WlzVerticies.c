#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzVerticies.c
* \author       Bill Hill
* \date         November 2000
* \version      $Id$
* \note
*               Copyright:
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par  Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions for extracting verticies from objects represented
*		by verticies, eg polylines, boundlists and contours.
* \ingroup      WlzFeatures
* \todo
* \bug
*/

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
static WlzVertexP 		WlzVerticiesFromGM2(
				  WlzGMModel *model,
				  WlzVertexP *dstNr,
				  int **dstVId,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
static WlzVertexP 		WlzVerticiesFromGM3(
				  WlzGMModel *model,
				  WlzVertexP *dstNr,
				  int **dstVId,
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
				  WlzDVertex2 v0,
				  WlzDVertex2 v1);
static WlzDVertex2 		WlzVerticiesNormTriple2(
				  WlzDVertex2 v0,
				  WlzDVertex2 v1,
				  WlzDVertex2 v2);

/*!
* \ingroup      WlzFeatures
* \return				Allocated verticies.
* \brief	Allocates a buffer which it fills with the verticies
*		from the given object. The object must be one of the
*		types that is represented by verticies, eg
*		polylines, boundlists and contours.
* \param	obj			Given polygon domain object.
* \param	dstNr			Destination ptr for normals.
*					The normals will always be
*					either WlzDVertex2 or
*					WlzDVertex3. May be NULL.
* \param	dstCnt			Destination ptr for the number
*					of verticies. Can NOT be NULL.
* \param	dstType			Destination ptr for the type
*					of verticies. Can NOT be NULL.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
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
	vData = WlzVerticiesFromCtr(obj->domain.ctr, dstNr, NULL,
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

/*!
* \ingroup      WlzFeatures
* \return				Allocated verticies.
* \brief	Allocates a buffer which it fills with the verticies
*		from a 3D polygon domain.
* \param	poly			Given polygon domain.
* \param	dstNr			Destination ptr for normals,
*					may be NULL.
* \param	dstCnt			Destination ptr for the number
*					of verticies.
* \param	dstType			Destination ptr for the type
*					of verticies.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
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

/*!
* \ingroup      WlzFeatures
* \return				Allocated verticies.
* \brief	Allocates a buffer which it fills with the verticies
*		from a boundary domain.
* \param	bound			Given boundary domain.
* \param	dstNr			Destination ptr for normals,
*					may be NULL.
* \param	dstCnt			Destination ptr for the number
*					of verticies.
* \param	dstType			Destination ptr for the type
*					of verticies.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
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

/*!
* \ingroup	WlzFeatures
* \return				Allocated verticies.
* \brief	Allocates a buffer which it fills with the verticies
*		from a contour.
* \param	ctr			Given contour.
* \param	dstNr			Destination ptr for normals,
*					may be NULL.
* \param	dstVId			Destination ptr for the GM
*					vertex indicies, may be NULL.
* \param	dstCnt			Destination ptr for the number
*					of verticies.
* \param	dstType			Destination ptr for the type
*					of verticies.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzVertexP 	WlzVerticiesFromCtr(WlzContour *ctr,
				    WlzVertexP *dstNr, int **dstVId,
				    int *dstCnt, WlzVertexType *dstType,
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
	vData = WlzVerticiesFromGM2(ctr->model, dstNr, dstVId, dstCnt, dstType,
				    &errNum);
        break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	vData = WlzVerticiesFromGM3(ctr->model, dstNr, dstVId, dstCnt, dstType,
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

/*!
* \ingroup	WlzFeatures
* \return				Allocated verticies.
* \brief	Allocates a buffer which it fills with the verticies
*		from a 2D GM.
* \param	model			Given model.
* \param	dstNr			Destination ptr for normals,
*					may be NULL.
* \param	dstVId			Destination ptr for the GM
*					vertex indicies, may be NULL.
* \param	dstCnt			Destination ptr for the number
*					of verticies.
* \param	dstType			Destination ptr for the type
*					of verticies.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzVertexP WlzVerticiesFromGM2(WlzGMModel *model,
				      WlzVertexP *dstNr, int **dstVId,
				      int *dstCnt, WlzVertexType *dstType,
				      WlzErrorNum *dstErr)
{
  int		idx,
  		cnt,
		nIdx,
		vIdx;
  int		*vId = NULL;
  WlzVertexType	type;
  WlzVertexP    vData;
  WlzGMVertex	*cV;
  WlzGMEdgeT	*cET;
  WlzGMVertexT	*cVT,
  		*nVT,
		*pVT;
  AlcVector	*vec;
  WlzDVertex2	*vNorm = NULL;
  WlzGMVertex	*nV[2];
  WlzVertexP	tVP[3];
  WlzDVertex2	segV[3];
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
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstNr)
    {
      if((vNorm = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * cnt)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if((errNum == WLZ_ERR_NONE) && dstVId)
    {
      if((vId = (int *)AlcMalloc(sizeof(int) * cnt)) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idx = 0; idx < cnt; ++idx)
    {
      cV = (WlzGMVertex *)AlcVectorItemGet(vec, vIdx++);
      if(cV->idx >= 0)
      {
	*(vId + idx) = cV->idx;
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
	  cET = cVT->parent;
	  if((cET == NULL) || (cVT->prev != cVT->next))
	  {
	    /* Vertex is either an isolated vertex or used by more than
	     * two edges. Normal undefined. */
	    (vNorm + idx)->vtX = 0.0;
	    (vNorm + idx)->vtY = 0.0;
	  }
	  else if((nV[0] = cET->prev->vertexT->diskT->vertex) ==
	          (nV[1] = cET->next->vertexT->diskT->vertex))
	  {
	    /* Vertex is on the end of a contour segment. */
	    if(model->type == WLZ_GMMOD_2I)
	    { 
	      tVP[0].i2 = &(nV[0]->geo.vg2I->vtx);
	      tVP[1].i2 = &(cV->geo.vg2I->vtx);
	      segV[0].vtX = tVP[0].i2->vtX;
	      segV[0].vtY = tVP[0].i2->vtY;
	      segV[1].vtX = tVP[1].i2->vtX;
	      segV[1].vtY = tVP[1].i2->vtY;
	    }
	    else /* model->type == WLZ_GMMOD_2D */
	    {
	      tVP[0].d2 = &(nV[0]->geo.vg2D->vtx);
	      tVP[1].d2 = &(cV->geo.vg2D->vtx);
	      segV[0] = *(tVP[0].d2);
	      segV[1] = *(tVP[1].d2);
	    }
	    *(vNorm + idx) = WlzVerticiesNormPair2(segV[0], segV[1]);
	  }
	  else
	  {
	    /* Vertex is used by two edges. Find the other two verticies
	     * that are used by these two edges. */
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
	    *(vNorm + idx) = WlzVerticiesNormTriple2(segV[0], segV[1],
	    					     segV[2]);
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
    if(dstVId)
    {
      *dstVId = vId;
    }
  }
  else
  {
    AlcFree(vNorm);
    AlcFree(vId);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/*!
* \ingroup WlzFeatures
* \return				Allocated verticies.
* \brief	Allocates a buffer which it fills with the verticies
*		from a 3D GM.
* \param	model			Given model.
* \param	dstNr			Destination ptr for normals,
*					may be NULL.
* \param	dstVId			Destination ptr for the GM
*					vertex indicies, may be NULL.
* \param	dstCnt			Destination ptr for the number
*					of verticies.
* \param	dstType			Destination ptr for the type
*					of verticies.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzVertexP WlzVerticiesFromGM3(WlzGMModel *model,
				      WlzVertexP *dstNr, int **dstVId,
				      int *dstCnt, WlzVertexType *dstType,
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
  int		*vId = NULL;
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
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstNr)
    {
      if((vNorm = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) *
      					   vCnt)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if((errNum == WLZ_ERR_NONE) && dstVId)
    {
      if((vId = (int *)AlcMalloc(sizeof(int) * vCnt)) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
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
	*(vId + vIdx) = cV->idx;
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
    if(dstVId)
    {
      (*dstVId) = vId;
    }
  }
  else
  {
    AlcFree(vNorm);
    AlcFree(vId);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vData);
}

/*!
* \ingroup 	WlzFeatures
* \return				Number of verticies in boundary.
* \brief	Counts the number of verticies in the polygon domains
*		of the boundary. This is a recursive function.
*  \param	bound			Given boundary domain.
*/
static int	WlzVerticiesCntBound(WlzBoundList *bound)
{
  int		cnt;

  cnt = bound->poly->nvertices + WlzVerticiesCntBound(bound->next) +
  	WlzVerticiesCntBound(bound->down);
  return(cnt);
}

/*!
* \ingroup 	WlzFeatures
* \return				Woolz error code.
* \brief	Copies verticies from the boundaries polygon domain
*		to the buffer.
* \param	vData			Given buffer.
* \param	vNorm			Given buffer for normals, may
*					be NULL.
* \param	vType			Type of verticies.
* \param	off			Ptr to offset into buffer.
* \param	bound			Given boundary domain.
*/
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

/*!
* \ingroup 	WlzFeatures
* \return				Allocated verticies.
* \brief	Allocates a buffer for copting the verticies of a
*		polygon domain.
* \param	polyType		Type of polygon domain.
* \param	cnt			Number of verticies to allocate
*					room for.
* \param	dstType			Destination ptr for the type
*					of verticies.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
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

/*!
* \ingroup	WlzFeatures
* \return	<void>
* \brief	Computes the normals of the given verticies which are
*		assumed to lie in a 2D polygon.
*		The normals all have +ve x components. 
*		There is no need for any type checking in this function
*		because it is static and all types have been checked.
* \param	nrm			Buffer for the normals.
* \param	vtx			The given verticies.
*		cnt			The number of verticies (and
*					normals).
* \param	pType			Polygon type.
*/
static void	WlzVerticiesNorm2(WlzDVertex2 *nrm, WlzVertexP vtx, int cnt,
				  WlzObjectType pType)
{
  int		idx,
  		idx1;
  WlzVertexType vType;
  WlzDVertex2	segV[3];

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
      *nrm = WlzVerticiesNormPair2(segV[0], segV[1]);
      break;
    default:
      /* There are more than two verticies.
       * Normals are computed using WlzVerticiesNormTriple2(). */
      idx = 0;
      idx1 = 2;
      switch(pType)
      {
	case WLZ_POLYGON_INT:
	  segV[0].vtX = (vtx.i2 + cnt - 1)->vtX; 
	  segV[0].vtY = (vtx.i2 + cnt - 1)->vtY; 
	  segV[1].vtX = (vtx.i2 + 0)->vtX;
	  segV[1].vtY = (vtx.i2 + 0)->vtY;
	  segV[2].vtX = (vtx.i2 + 1)->vtX;
	  segV[2].vtY = (vtx.i2 + 1)->vtY;
	  break;
	case WLZ_POLYGON_FLOAT:
	  segV[0].vtX = (vtx.f2 + cnt - 1)->vtX; 
	  segV[0].vtY = (vtx.f2 + cnt - 1)->vtY; 
	  segV[1].vtX = (vtx.f2 + 0)->vtX;
	  segV[1].vtY = (vtx.f2 + 0)->vtY;
	  segV[2].vtX = (vtx.f2 + 1)->vtX;
	  segV[2].vtY = (vtx.f2 + 1)->vtY;
	  break;
	case WLZ_POLYGON_DOUBLE:
	  segV[0] = *(vtx.d2 + cnt - 1);
	  segV[1] = *(vtx.d2 + 0);
	  segV[2] = *(vtx.d2 + 1);
	  break;
      }
      *(nrm + 0) = WlzVerticiesNormTriple2(segV[0], segV[1], segV[2]);
      while(idx1 < cnt)
      {
        segV[0] = segV[1];
        segV[1] = segV[2];
	switch(pType)
	{
	  case WLZ_POLYGON_INT:
	    segV[2].vtX = (vtx.i2 + idx1)->vtX;
	    segV[2].vtY = (vtx.i2 + idx1)->vtY;
	    break;
	  case WLZ_POLYGON_FLOAT:
	    segV[2].vtX = (vtx.f2 + idx1)->vtX;
	    segV[2].vtY = (vtx.f2 + idx1)->vtY;
	    break;
	  case WLZ_POLYGON_DOUBLE:
	    segV[2] = *(vtx.d2 + idx1);
	    break;
	}
	*(nrm + ++idx) = WlzVerticiesNormTriple2(segV[0], segV[1], segV[2]);
        ++idx1;
      }
      segV[0] = segV[1];
      segV[1] = segV[2];
      switch(pType)
      {
	case WLZ_POLYGON_INT:
	  segV[2].vtX = (vtx.i2 + 0)->vtX;
	  segV[2].vtY = (vtx.i2 + 0)->vtY;
	  break;
	case WLZ_POLYGON_FLOAT:
	  segV[2].vtX = (vtx.f2 + 0)->vtX;
	  segV[2].vtY = (vtx.f2 + 0)->vtY;
	  break;
	case WLZ_POLYGON_DOUBLE:
	  segV[2] = *(vtx.d2 + 0);
	  break;
      }
      *(nrm + ++idx) = WlzVerticiesNormTriple2(segV[0], segV[1], segV[2]);
      break;
  }
}

/*!
* \ingroup	WlzFeatures
* \return	<void>
* \brief	Computes the normal (n) to a segment (g) between the
*		given pair of verticies. There are clearly two solutions
*		to the problem of finding a normal to a line segment,
*		but this function always finds the normal vector with a
*		+ve x component.
*		If the two verticies are coincident then the normal
*		vector is set to {0, 0}.
*		With two non-coincident verticies the normal vector is
*		computed using the relationships g.n = 0 and \|n\|^2 = 1.
*		Giving:
*		\f$ nx = \frac{1}{\sqrt(1+(gx/gy)^2)} \f$ ,
*		\f$ ny = -nx\frac{gx}{gy} \f$
*		There is no need for any type checking in this function
*		because it is static and all types have been checked.
* \param	v0			First of the given pair.
* \param	v1			Second of the given pair.
*/
static WlzDVertex2 WlzVerticiesNormPair2(WlzDVertex2 v0, WlzDVertex2 v1)
{
  WlzDVertex2	tV0,
  		tV1,
		nrm;

  WLZ_VTX_2_SUB(tV0, v1, v0);
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

/*!
* \return				Normal vector.
* \brief	Computes the normal (n) at a vertex. This is chosen to
*		be the mean of normals of the two line segments which the
*  		vertex is common to. This normal also bisects the two
*		angles of the line segments at the vertex.
*		Considering the three given verticies (A,B and C) to
*		form a triangle ABC, a line which bisects the angles
*		at B from some point D on the line segment CA to B.
*		The point D is given by:
*		\f$ p_D = p_A + (p_C - p_A)
                                \frac{\|p_A - p_B\|}{\|p_B - p_C\|} \f$
*		So the normal n is given by:
*		\f$ n = \frac{(p_B - p_D)}{\|p_B - p_D\|} \f$
*		unless \f$ \|p_B - p_D\| < \epsilon \f$ in which case
*		WlzVerticiesNormPair2() is used to compute the normal
*		from A and C.
* \param	v0			First of the given triple (A).
* \param	v1			Second of the given triple (B).
* \param	v2			Third of the given triple (C).
*/
static WlzDVertex2 WlzVerticiesNormTriple2(WlzDVertex2 v0, WlzDVertex2 v1,
					   WlzDVertex2 v2)
{
  double	tD0,
  		tD1;
  WlzDVertex2	tV0,
  		tV1,
		nrm;

  WLZ_VTX_2_SUB(tV0, v0, v1);
  tD0 = WLZ_VTX_2_SQRLEN(tV0);
  WLZ_VTX_2_SUB(tV0, v1, v2);
  tD1 = WLZ_VTX_2_SQRLEN(tV0);
  if(tD0 < DBL_EPSILON)
  {
    nrm = WlzVerticiesNormPair2(v1, v2);
  }
  else if(tD1 < DBL_EPSILON)
  {
    nrm = WlzVerticiesNormPair2(v0, v1);
  }
  else
  {
    tD0 /= tD0 + tD1;
    WLZ_VTX_2_SUB(tV0, v2, v0);
    WLZ_VTX_2_SCALE(tV1, tV0, tD0);
    WLZ_VTX_2_ADD(tV0, tV1, v0);
    WLZ_VTX_2_SUB(tV1, v1, tV0);
    tD0 = WLZ_VTX_2_LENGTH(tV1);
    if(tD0 < DBL_EPSILON)
    {
      nrm = WlzVerticiesNormPair2(v0, v2);
    }
    else
    {
      tD1 = 1.0 / tD0;
      WLZ_VTX_2_SCALE(nrm, tV1, tD1);
    }
  }
  return(nrm);
}
