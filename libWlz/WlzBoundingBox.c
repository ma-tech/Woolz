#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBoundingBox_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzBoundingBox.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	Functions for computing the axis aligned bounding box
* 		of objects.
* \ingroup	WlzFeatures
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <Wlz.h>

static WlzDBox3			WlzBoundingBoxBound3D(
				  WlzBoundList *,
				  WlzErrorNum *);
static WlzDBox3			WlzBoundingBoxPoly3D(
				  WlzPolygonDomain *,
				  WlzErrorNum *);
static WlzDBox3			WlzBoundingBoxTransObj3D(
				  WlzObject *,
				  WlzAffineTransform *,
				  WlzErrorNum *);
static WlzDBox3 		WlzBoundingBoxContour3D(
				  WlzContour *ctr,
				  WlzErrorNum *dstErr);
static WlzDBox3 		WlzBoundingBoxCompound3D(
				  WlzCompoundArray *cObj,
				  WlzErrorNum *dstErr);

/*!
* \return	2D integer bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the 2D integer axis aligned bounding box of the
*		given object.
* \param	inObj			The given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIBox2	WlzBoundingBox2I(WlzObject *inObj, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIBox2	bBox2I;
  WlzIBox3	bBox3I;

  bBox3I = WlzBoundingBox3I(inObj, &errNum);
  bBox2I.xMin = bBox3I.xMin;
  bBox2I.yMin = bBox3I.yMin;
  bBox2I.xMax = bBox3I.xMax;
  bBox2I.yMax = bBox3I.yMax;
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox2I);
}

/*!
* \return	2D double precision bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the 2D double precision axis aligned bounding box
*		of the given object.
* \param	inObj			The given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDBox2	WlzBoundingBox2D(WlzObject *inObj, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDBox2	bBox2D;
  WlzDBox3	bBox3D;

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

/*!
* \return	3D integer bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the integer 3D axis aligned bounding box of the
*		given object.
* \param	inObj			The given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIBox3	WlzBoundingBox3I(WlzObject *inObj, WlzErrorNum *dstErr)
{
  WlzIBox3	bBox3I;
  WlzDBox3	bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3I.xMin = 0;
  bBox3I.xMax = 0;
  bBox3I.yMin = 0;
  bBox3I.yMax = 0;
  bBox3I.zMin = 0;
  bBox3I.zMax = 0;
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
	      bBox3I.xMin = inObj->domain.i->kol1;
	      bBox3I.yMin = inObj->domain.i->line1;
	      bBox3I.xMax = inObj->domain.i->lastkl;
	      bBox3I.yMax = inObj->domain.i->lastln;
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
	    bBox3I.xMin = inObj->domain.p->kol1;
	    bBox3I.yMin = inObj->domain.p->line1;
	    bBox3I.zMin = inObj->domain.p->plane1;
	    bBox3I.xMax = inObj->domain.p->lastkl;
	    bBox3I.yMax = inObj->domain.p->lastln;
	    bBox3I.zMax = inObj->domain.p->lastpl;
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
	else if(inObj->domain.core->type != (WlzObjectType)
	                                    WLZ_TRANSFORM_2D_AFFINE)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  bBox3D = WlzBoundingBoxTransObj3D(inObj->values.obj,
	  				    inObj->domain.t, &errNum);
	  bBox3I = WlzBoundingBox3DTo3I(bBox3D);
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
	  bBox3I = WlzBoundingBox3DTo3I(bBox3D);
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
	  bBox3I = WlzBoundingBox3DTo3I(bBox3D);
	}
	break;
      case WLZ_CONTOUR:
	if(inObj->domain.ctr == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  bBox3D = WlzBoundingBoxContour3D(inObj->domain.ctr, &errNum);
	  bBox3I = WlzBoundingBox3DTo3I(bBox3D);
	}
        break;
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
        bBox3D = WlzBoundingBoxCompound3D((WlzCompoundArray *)inObj, &errNum);
	bBox3I = WlzBoundingBox3DTo3I(bBox3D);
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
  return(bBox3I);
}

/*!
* \return	3D double precision bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the double precision 3D axis aligned bounding box
*		of the given object.
* \param	inObj			The given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDBox3	WlzBoundingBox3D(WlzObject *inObj, WlzErrorNum *dstErr)
{
  WlzDBox3	bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D.xMin = 0.0;
  bBox3D.xMax = 0.0;
  bBox3D.yMin = 0.0;
  bBox3D.yMax = 0.0;
  bBox3D.zMin = 0.0;
  bBox3D.zMax = 0.0;
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
	else if(inObj->domain.core->type != (WlzObjectType )
	                                    WLZ_TRANSFORM_2D_AFFINE)
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
	  bBox3D = WlzBoundingBoxContour3D(inObj->domain.ctr, &errNum);
	}
        break;
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
        bBox3D = WlzBoundingBoxCompound3D((WlzCompoundArray *)inObj, &errNum);
	break;
      case WLZ_CMESH_2D:
	if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(inObj->domain.core->type != WLZ_CMESH_2D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  WlzCMeshUpdateBBox2D(inObj->domain.cm2);
          bBox3D.xMin = inObj->domain.cm2->bBox.xMin;
          bBox3D.yMin = inObj->domain.cm2->bBox.yMin;
          bBox3D.zMin = 0.0;
          bBox3D.xMax = inObj->domain.cm2->bBox.xMax;
          bBox3D.yMax = inObj->domain.cm2->bBox.yMax;
          bBox3D.zMax = 0.0;
	}
	break;
      case WLZ_CMESH_2D5:
	if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(inObj->domain.core->type != WLZ_CMESH_2D5)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  WlzCMeshUpdateBBox2D5(inObj->domain.cm2d5);
          bBox3D.xMin = inObj->domain.cm2d5->bBox.xMin;
          bBox3D.yMin = inObj->domain.cm2d5->bBox.yMin;
          bBox3D.zMin = inObj->domain.cm2d5->bBox.zMin;
          bBox3D.xMax = inObj->domain.cm2d5->bBox.xMax;
          bBox3D.yMax = inObj->domain.cm2d5->bBox.yMax;
          bBox3D.zMax = inObj->domain.cm2d5->bBox.zMax;
	}
	break;
      case WLZ_CMESH_3D:
	if(inObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(inObj->domain.core->type != WLZ_CMESH_3D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  WlzCMeshUpdateBBox3D(inObj->domain.cm3);
          bBox3D.xMin = inObj->domain.cm3->bBox.xMin;
          bBox3D.yMin = inObj->domain.cm3->bBox.yMin;
          bBox3D.zMin = inObj->domain.cm3->bBox.zMin;
          bBox3D.xMax = inObj->domain.cm3->bBox.xMax;
          bBox3D.yMax = inObj->domain.cm3->bBox.yMax;
          bBox3D.zMax = inObj->domain.cm3->bBox.zMax;
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

/*!
* \return	Double precision 3D bounding box.
* \brief	Computes the bounding box of the given 3D double precision
*		vertices.
* \param	nVtx			Number of vertices.
* \param	vtx			Vertices.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDBox3	WlzBoundingBoxVtx3D(int nVtx, WlzDVertex3 *vtx,
				    WlzErrorNum *dstErr)
{
  int		cnt;
  WlzDBox3	bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nVtx <= 0) || (vtx == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    cnt = nVtx;
    bBox3D.xMin = bBox3D.xMax = vtx->vtX;
    bBox3D.yMin = bBox3D.yMax = vtx->vtY;
    bBox3D.zMin = bBox3D.zMax = vtx->vtZ;
    while(--cnt > 0)
    {
      ++vtx;
      if(vtx->vtX < bBox3D.xMin)
      {
	bBox3D.xMin = vtx->vtX;
      }
      else if(vtx->vtX > bBox3D.xMax)
      {
        bBox3D.xMax = vtx->vtX;
      }
      if(vtx->vtY < bBox3D.yMin)
      {
	bBox3D.yMin = vtx->vtY;
      }
      else if(vtx->vtY > bBox3D.yMax)
      {
        bBox3D.yMax = vtx->vtY;
      }
      if(vtx->vtZ < bBox3D.zMin)
      {
	bBox3D.zMin = vtx->vtZ;
      }
      else if(vtx->vtZ > bBox3D.zMax)
      {
        bBox3D.zMax = vtx->vtZ;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/*!
* \return	Integer 3D bounding box.
* \brief	Computes the bounding box of the given 3D integer vertices.
* \param	nVtx			Number of vertices.
* \param	vtx			Vertices.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIBox3	WlzBoundingBoxVtx3I(int nVtx, WlzIVertex3 *vtx,
				    WlzErrorNum *dstErr)
{
  int		cnt;
  WlzIBox3	bBox3I;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nVtx <= 0) || (vtx == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    cnt = nVtx;
    bBox3I.xMin = bBox3I.xMax = vtx->vtX;
    bBox3I.yMin = bBox3I.yMax = vtx->vtY;
    bBox3I.zMin = bBox3I.zMax = vtx->vtZ;
    while(--cnt > 0)
    {
      ++vtx;
      if(vtx->vtX < bBox3I.xMin)
      {
	bBox3I.xMin = vtx->vtX;
      }
      else if(vtx->vtX > bBox3I.xMax)
      {
        bBox3I.xMax = vtx->vtX;
      }
      if(vtx->vtY < bBox3I.yMin)
      {
	bBox3I.yMin = vtx->vtY;
      }
      else if(vtx->vtY > bBox3I.yMax)
      {
        bBox3I.yMax = vtx->vtY;
      }
      if(vtx->vtZ < bBox3I.zMin)
      {
	bBox3I.zMin = vtx->vtZ;
      }
      else if(vtx->vtZ > bBox3I.zMax)
      {
        bBox3I.zMax = vtx->vtZ;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3I);
}

/*!
* \return	3D bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the 3D axis aligned bounding box of the
*		transformed object.
* \param	inObj			Input object.
* \param	trans			Given transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzDBox3	WlzBoundingBoxTransObj3D(WlzObject *inObj,
					 WlzAffineTransform *trans,
				         WlzErrorNum *dstErr)
{
  WlzDBox3	tBox3D,
  		bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D.xMin = 0.0;
  bBox3D.xMax = 0.0;
  bBox3D.yMin = 0.0;
  bBox3D.yMax = 0.0;
  bBox3D.zMin = 0.0;
  bBox3D.zMax = 0.0;
  tBox3D = WlzBoundingBox3D(inObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    bBox3D = WlzAffineTransformBBoxD3(trans, tBox3D, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/*!
* \return	3D bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 3D axis aligned bounding box of the
* 		compound array object.
* \param	WlzCompoundArray *cObjdstErr
*/
static WlzDBox3 WlzBoundingBoxCompound3D(WlzCompoundArray *cObj,
					 WlzErrorNum *dstErr)
{
  WlzDBox3	tBox3D,
  		bBox3D;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(cObj->n < 1)
  {
    errNum = WLZ_ERR_OBJECT_DATA;
  }
  else
  {
    int		idx = 0;

    while((idx < cObj->n) &&
          ((cObj->o[idx] == NULL) || (cObj->o[idx]->type == WLZ_EMPTY_OBJ)))
    {
      ++idx;
    }
    if(idx >= cObj->n)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else
    {
      bBox3D = WlzBoundingBox3D(cObj->o[idx], &errNum);
      while((errNum == WLZ_ERR_NONE) && (++idx < cObj->n))
      {
	if((cObj->o[idx] != NULL) && (cObj->o[idx]->type != WLZ_EMPTY_OBJ))
	{
	  tBox3D = WlzBoundingBox3D(cObj->o[idx], &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    bBox3D = WlzBoundingBoxUnion3D(bBox3D, tBox3D);
	  }
	}
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/*!
* \return	3D bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the 3D axis aligned bounding box of the polygon
* 		domain.
* \param	poly			The given polygon domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzDBox3	WlzBoundingBoxPoly3D(WlzPolygonDomain *poly,
				     WlzErrorNum *dstErr)
{
  int		tI0,
  		idx;
  WlzIVertex2	*tIVxP;
  WlzFVertex2	*tFVxP;
  WlzDVertex2	*tDVxP;
  WlzDBox3	bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D.xMin = 0.0;
  bBox3D.xMax = 0.0;
  bBox3D.yMin = 0.0;
  bBox3D.yMax = 0.0;
  bBox3D.zMin = 0.0;
  bBox3D.zMax = 0.0;
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
	bBox3D.xMax = bBox3D.xMin = WLZ_NINT(tFVxP->vtX);
	bBox3D.yMax = bBox3D.yMin = WLZ_NINT(tFVxP->vtY);
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
	bBox3D.xMax = bBox3D.xMin = WLZ_NINT(tDVxP->vtX);
	bBox3D.yMax = bBox3D.yMin = WLZ_NINT(tDVxP->vtY);
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

/*!
* \return	3D bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the 3D axis aligned bounding box of the
*		boundary list.
* \param	bound			The given boundary list.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzDBox3	WlzBoundingBoxBound3D(WlzBoundList *bound,
				      WlzErrorNum *dstErr)
{
  WlzDBox3	bBox3D,
  		tBox3D;
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
    if(tBox3D.xMax > bBox3D.xMax)
    {
      bBox3D.xMax = tBox3D.xMax;
    }
    if(tBox3D.yMax > bBox3D.yMax)
    {
      bBox3D.yMax = tBox3D.yMax;
    }
    if(tBox3D.zMax > bBox3D.zMax)
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

/*!
* \return	3D bounding box.
* \ingroup      WlzFeatures
* \brief	Computes the 3D axis aligned bounding box of the contour.
* \param	ctr			The given boundary list.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzDBox3 WlzBoundingBoxContour3D(WlzContour *ctr, WlzErrorNum *dstErr)
{
  WlzDBox3	bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox3D = WlzBoundingBoxGModel3D(ctr->model, &errNum);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/*!
* \return	3D bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 3D axis aligned bounding box of the 3D geometric
* 		model. It is an error if the model is not a 3D model.
* \param	model			Given model.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDBox3 	WlzBoundingBoxGModel3D(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMShell	*shell0,
  		*shell1;
  WlzDBox3	sBox3D,
  		bBox3D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum =  WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(model->type)
    {
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((shell0 = model->child) == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find the bounding box of all the GM's shell's bounding boxes. */
    shell1 = shell0->next;
    errNum = WlzGMShellGetGBB3D(shell0, &bBox3D);
    while((errNum == WLZ_ERR_NONE) && (shell1->idx != shell0->idx))
    {
      if((errNum = WlzGMShellGetGBB3D(shell1, &sBox3D)) == WLZ_ERR_NONE)
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
	if(sBox3D.xMax > bBox3D.xMax)
	{
	  bBox3D.xMax = sBox3D.xMax;
	}
	if(sBox3D.yMax > bBox3D.yMax)
	{
	  bBox3D.yMax = sBox3D.yMax;
	}
	if(sBox3D.zMax > bBox3D.zMax)
	{
	  bBox3D.zMax = sBox3D.zMax;
	}
        shell1 = shell1->next;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox3D);
}

/*!
* \return	2D bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 2D axis aligned bounding box of the 2D geometric
* 		model. It is an error if the model is not a 2D model.
* \param	model			Given model.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDBox2 	WlzBoundingBoxGModel2D(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMShell	*shell0,
  		*shell1;
  WlzDBox2	sBox2D,
  		bBox2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum =  WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N:
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((shell0 = model->child) == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find the bounding box of all the GM's shell's bounding boxes. */
    shell1 = shell0->next;
    errNum = WlzGMShellGetGBB2D(shell0, &bBox2D);
    while((errNum == WLZ_ERR_NONE) && (shell1->idx != shell0->idx))
    {
      if((errNum = WlzGMShellGetGBB2D(shell1, &sBox2D)) == WLZ_ERR_NONE)
      {
	if(sBox2D.xMin < bBox2D.xMin)
	{
	  bBox2D.xMin = sBox2D.xMin;
	}
	if(sBox2D.yMin < bBox2D.yMin)
	{
	  bBox2D.yMin = sBox2D.yMin;
	}
	if(sBox2D.xMax > bBox2D.xMax)
	{
	  bBox2D.xMax = sBox2D.xMax;
	}
	if(sBox2D.yMax > bBox2D.yMax)
	{
	  bBox2D.yMax = sBox2D.yMax;
	}
        shell1 = shell1->next;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bBox2D);
}

/*!
* \return	Integer bounding box.
* \ingroup	WlzFeatures
* \brief	Converts a 2D double precision bounding box to an integer
*		bounding box. There is no checking done for underflows or
*		overflows.
* \param	bBox2D			Double precision bounding box.
*/
WlzIBox2	WlzBoundingBox2DTo2I(WlzDBox2 bBox2D)
{
  WlzIBox2	bBox2I;

  bBox2I.xMin = (int )floor(bBox2D.xMin);
  bBox2I.yMin = (int )floor(bBox2D.yMin);
  bBox2I.xMax = (int )ceil(bBox2D.xMax);
  bBox2I.yMax = (int )ceil(bBox2D.yMax);
  return(bBox2I);
}

/*!
* \return	Double precision bounding box.
* \ingroup	WlzFeatures
* \brief	Converts a 2D integer bounding box to an double precision
*		bounding box.
* \param	bBox2I			Double precision bounding box.
*/
WlzDBox2	WlzBoundingBox2ITo2D(WlzIBox2 bBox2I)
{
  WlzDBox2	bBox2D;

  bBox2D.xMin = bBox2I.xMin;
  bBox2D.yMin = bBox2I.yMin;
  bBox2D.xMax = bBox2I.xMax;
  bBox2D.yMax = bBox2I.yMax;
  return(bBox2D);
}

/*!
* \return	Integer bounding box.
* \ingroup	WlzFeatures
* \brief	Converts a 3D double precision bounding box to an integer
*		bounding box. There is no checking done for underflows or
*		overflows.
* \param	bBox3D			Double precision bounding box.
*/
WlzIBox3	WlzBoundingBox3DTo3I(WlzDBox3 bBox3D)
{
  WlzIBox3	bBox3I;

  bBox3I.xMin = (int )floor(bBox3D.xMin);
  bBox3I.yMin = (int )floor(bBox3D.yMin);
  bBox3I.zMin = (int )floor(bBox3D.zMin);
  bBox3I.xMax = (int )ceil(bBox3D.xMax);
  bBox3I.yMax = (int )ceil(bBox3D.yMax);
  bBox3I.zMax = (int )ceil(bBox3D.zMax);
  return(bBox3I);
}

/*!
* \return	Double precision bounding box.
* \ingroup	WlzFeatures
* \brief	Converts a 3D integer bounding box to an double precision
*		bounding box.
* \param	bBox3I			Integer bounding box.
*/
WlzDBox3	WlzBoundingBox3ITo3D(WlzIBox3 bBox3I)
{
  WlzDBox3	bBox3D;

  bBox3D.xMin = bBox3I.xMin;
  bBox3D.yMin = bBox3I.yMin;
  bBox3D.zMin = bBox3I.zMin;
  bBox3D.xMax = bBox3I.xMax;
  bBox3D.yMax = bBox3I.yMax;
  bBox3D.zMax = bBox3I.zMax;
  return(bBox3D);
}

/*!
* \return	Single precision bounding box.
* \ingroup	WlzFeatures
* \brief	Converts a 3D double bounding box to a single precision
*		bounding box.
* \param	bBox3D			Double precision bounding box.
*/
WlzFBox3	WlzBoundingBox3DTo3F(WlzDBox3 bBox3D)
{
  WlzFBox3	bBox3F;

  bBox3F.xMin = (float )(bBox3D.xMin);
  bBox3F.yMin = (float )(bBox3D.yMin);
  bBox3F.zMin = (float )(bBox3D.zMin);
  bBox3F.xMax = (float )(bBox3D.xMax);
  bBox3F.yMax = (float )(bBox3D.yMax);
  bBox3F.zMax = (float )(bBox3D.zMax);
  return(bBox3F);
}

/*!
* \return	Double precision bounding box.
* \ingroup	WlzFeatures
* \brief	Converts a 3D single bounding box to a double precision
*		bounding box.
* \param	bBox3F			Single precision bounding box.
*/
WlzDBox3	WlzBoundingBox3FTo3D(WlzFBox3 bBox3F)
{
  WlzDBox3	bBox3D;

  bBox3D.xMin = bBox3F.xMin;
  bBox3D.yMin = bBox3F.yMin;
  bBox3D.zMin = bBox3F.zMin;
  bBox3D.xMax = bBox3F.xMax;
  bBox3D.yMax = bBox3F.yMax;
  bBox3D.zMax = bBox3F.zMax;
  return(bBox3D);
}

/*!
* \return	Axis alligned 2D integer bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 2D integer bounding box which
*		encloses the given pair of bounding boxes.
* \param	box0			First bounding box.
* \param	box1			Second bounding box.
*/
WlzIBox2	WlzBoundingBoxUnion2I(WlzIBox2 box0, WlzIBox2 box1)
{
  if(box1.xMin < box0.xMin)
  {
    box0.xMin = box1.xMin;
  }
  if(box1.xMax > box0.xMax)
  {
    box0.xMax = box1.xMax;
  }
  if(box1.yMin < box0.yMin)
  {
    box0.yMin = box1.yMin;
  }
  if(box1.yMax > box0.yMax)
  {
    box0.yMax = box1.yMax;
  }
  return(box0);
}

/*!
* \return	Axis alligned 2D single precision bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 2D single precision bounding box which
*		encloses the given pair of bounding boxes.
* \param	box0			First bounding box.
* \param	box1			Second bounding box.
*/
WlzFBox2	WlzBoundingBoxUnion2F(WlzFBox2 box0, WlzFBox2 box1)
{
  if(box1.xMin < box0.xMin)
  {
    box0.xMin = box1.xMin;
  }
  if(box1.xMax > box0.xMax)
  {
    box0.xMax = box1.xMax;
  }
  if(box1.yMin < box0.yMin)
  {
    box0.yMin = box1.yMin;
  }
  if(box1.yMax > box0.yMax)
  {
    box0.yMax = box1.yMax;
  }
  return(box0);
}

/*!
* \return	Axis alligned 2D double precision bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 2D double precision bounding box which
*		encloses the given pair of bounding boxes.
* \param	box0			First bounding box.
* \param	box1			Second bounding box.
*/
WlzDBox2	WlzBoundingBoxUnion2D(WlzDBox2 box0, WlzDBox2 box1)
{
  if(box1.xMin < box0.xMin)
  {
    box0.xMin = box1.xMin;
  }
  if(box1.xMax > box0.xMax)
  {
    box0.xMax = box1.xMax;
  }
  if(box1.yMin < box0.yMin)
  {
    box0.yMin = box1.yMin;
  }
  if(box1.yMax > box0.yMax)
  {
    box0.yMax = box1.yMax;
  }
  return(box0);
}

/*!
* \return	Axis alligned 3D integer bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 3D integer bounding box which
*		encloses the given pair of bounding boxes.
* \param	box0			First bounding box.
* \param	box1			Second bounding box.
*/
WlzIBox3	WlzBoundingBoxUnion3I(WlzIBox3 box0, WlzIBox3 box1)
{
  if(box1.xMin < box0.xMin)
  {
    box0.xMin = box1.xMin;
  }
  if(box1.xMax > box0.xMax)
  {
    box0.xMax = box1.xMax;
  }
  if(box1.yMin < box0.yMin)
  {
    box0.yMin = box1.yMin;
  }
  if(box1.yMax > box0.yMax)
  {
    box0.yMax = box1.yMax;
  }
  if(box1.zMin < box0.zMin)
  {
    box0.zMin = box1.zMin;
  }
  if(box1.zMax > box0.zMax)
  {
    box0.zMax = box1.zMax;
  }
  return(box0);
}

/*!
* \return	Axis alligned 3D float precision bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 3D single precision bounding box which
*		encloses the given pair of bounding boxes.
* \param	box0			First bounding box.
* \param	box1			Second bounding box.
*/
WlzFBox3	WlzBoundingBoxUnion3F(WlzFBox3 box0, WlzFBox3 box1)
{
  if(box1.xMin < box0.xMin)
  {
    box0.xMin = box1.xMin;
  }
  if(box1.xMax > box0.xMax)
  {
    box0.xMax = box1.xMax;
  }
  if(box1.yMin < box0.yMin)
  {
    box0.yMin = box1.yMin;
  }
  if(box1.yMax > box0.yMax)
  {
    box0.yMax = box1.yMax;
  }
  if(box1.zMin < box0.zMin)
  {
    box0.zMin = box1.zMin;
  }
  if(box1.zMax > box0.zMax)
  {
    box0.zMax = box1.zMax;
  }
  return(box0);
}

/*!
* \return	Axis alligned 3D double precision bounding box.
* \ingroup	WlzFeatures
* \brief	Computes the 3D double precision bounding box which
*		encloses the given pair of bounding boxes.
* \param	box0			First bounding box.
* \param	box1			Second bounding box.
*/
WlzDBox3	WlzBoundingBoxUnion3D(WlzDBox3 box0, WlzDBox3 box1)
{
  if(box1.xMin < box0.xMin)
  {
    box0.xMin = box1.xMin;
  }
  if(box1.xMax > box0.xMax)
  {
    box0.xMax = box1.xMax;
  }
  if(box1.yMin < box0.yMin)
  {
    box0.yMin = box1.yMin;
  }
  if(box1.yMax > box0.yMax)
  {
    box0.yMax = box1.yMax;
  }
  if(box1.zMin < box0.zMin)
  {
    box0.zMin = box1.zMin;
  }
  if(box1.zMax > box0.zMax)
  {
    box0.zMax = box1.zMax;
  }
  return(box0);
}
