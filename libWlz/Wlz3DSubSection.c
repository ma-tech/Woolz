#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _Wlz3DSubSection_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/Wlz3DSubSection.c
* \author       Richard Baldock, Bill Hill
* \date         May 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \brief        Return a sub-region of a 3D section via a 3D section transform.
* \ingroup      WlzSectionTransform
*/

#include <limits.h>
#include <float.h>

#include <Wlz.h>

#ifndef WLZ_FAST_CODE
#define WLZ_FAST_CODE
#endif

static WlzObject	*WlzGetSubSectionFrom3DTiledValueObj(
  WlzObject		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct	*viewStr,
  WlzErrorNum		*dstErr);
static WlzObject *WlzGetSubSectionFrom3DDomObj(
  WlzObject 		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct 	*viewStr,
  WlzInterpolationType	interp,
  WlzObject		**maskRtn,
  WlzErrorNum 		*dstErr);


/*!
* \return	New sub-section object.
* \ingroup	WlzSectionTransform
* \brief	Computes a section through the given 3D object.
* \param	obj			Given 3D object.
* \param	subDomain		Given 2D domain within which to
* 					restrict the section. If NULL
* 					returned section will be have a
* 					rectangular domain which is the maximum
* 					for the given object and view
* 					transform when the given 3D object
* 					is a spatial domain object.
* \param	view			Given view transform.
* \param	interp			Interpolation, should be either
* 					WLZ_INTERPOLATION_NEAREST or
* 					WLZ_INTERPOLATION_LINEAR.
* \param	maskRtn			Destination pointer for returned
* 					domain mask.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzGetSubSectionFromObject(
  WlzObject 		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct 	*view,
  WlzInterpolationType	interp,
  WlzObject		**maskRtn,
  WlzErrorNum 		*dstErr)
{
  WlzObject	*newObj = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

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
      case WLZ_3D_DOMAINOBJ:
	if(obj->values.core &&
	   (maskRtn == NULL) &&
	   (interp == WLZ_INTERPOLATION_NEAREST) &&
	   WlzGreyTableIsTiled(obj->values.core->type) &&
	   (obj->values.t->tileWidth == 16) &&
	   (obj->values.t->original_table.core == NULL))
	{
	  newObj = WlzGetSubSectionFrom3DTiledValueObj(obj, subDomain, view,
	                                               &errNum);
	}
	else
	{
	  newObj = WlzGetSubSectionFrom3DDomObj(obj, subDomain, view,
						interp, maskRtn, &errNum);
	}
        break;
      default:
	newObj = WlzGetSectionFromObject(obj, view, interp, &errNum);
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newObj);
}

/*!
* \return	New sub-section object.
* \ingroup	WlzSectionTransform
* \brief	Computes a section through the given 3D domain object
* 		where that object has tiled values and only the section
* 		values are required. Assumes a tile width of 16 voxels
* 		and that there is a valid background value in that part
* 		of each tile which falls outside of the object's domain.
* \param	obj			Given 3D object which is known to
* 					be non-NULL, have non-NULL domain
* 					and values with the values being
* 					tiled values.
* \param	subDomain		Given 2D domain within which to
* 					restrict the section. If NULL
* 					returned section will be have a
* 					rectangular domain which is the maximum
* 					for the given object and view
* 					transform.
* \param	view			Given view transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject	*WlzGetSubSectionFrom3DTiledValueObj(
  WlzObject		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct	*view,
  WlzErrorNum		*dstErr)
{
  WlzObject		*dstObj = NULL;
  WlzDomain		dstDom;
  WlzValues		dstVal;
  WlzGreyType		gType;
  WlzPixelV		bgd;
  WlzTiledValues        *tv;
  WlzIBox2 		subBox;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dstVal.core = NULL;
  if(view == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(view->type != WLZ_3D_VIEW_STRUCT)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(!(view->initialised))
  {
    errNum = WLZ_ERR_OBJECT_DATA;
  }
  else
  {
    tv = obj->values.t;
  }
  /* Make rectangular domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    subBox.xMin = WLZ_NINT(view->minvals.vtX);
    subBox.xMax = WLZ_NINT(view->maxvals.vtX);
    subBox.yMin = WLZ_NINT(view->minvals.vtY);
    subBox.yMax = WLZ_NINT(view->maxvals.vtY);
    if(subDomain)
    {
      switch(subDomain->type)
      {
	case WLZ_2D_DOMAINOBJ:
	  dstObj = WlzClipObjToBox2D(subDomain, subBox, &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    else
    {
      if((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					   subBox.yMin, subBox.yMax,
					   subBox.xMin, subBox.xMax,
					   &errNum)) != NULL)
      {
	dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, dstVal,
	                     NULL, NULL, &errNum);
      }
    }
  }
  /* Make and add a rectangular value table */
  if(errNum == WLZ_ERR_NONE)
  {
    subBox.xMin = dstObj->domain.i->kol1;
    subBox.xMax = dstObj->domain.i->lastkl;
    subBox.yMin = dstObj->domain.i->line1;
    subBox.yMax = dstObj->domain.i->lastln;
    bgd = WlzGetBackground(obj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gType = WlzGreyTypeFromObj(obj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make sure that the grey type is one that's handled here. This
     * simplifies/removes testing latter. */
    switch(gType)
    {
      case WLZ_GREY_INT:    /* FALLTHROUGH */
      case WLZ_GREY_SHORT:  /* FALLTHROUGH */
      case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
      case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE: /* FALLTHROUGH */
      case WLZ_GREY_RGBA:
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((dstVal.v = WlzNewValueTb(dstObj,
                                 WlzGreyTableType(WLZ_GREY_TAB_RECT,
				                  gType, NULL),
				 bgd, &errNum)))
    {
      dstObj->values = WlzAssignValues(dstVal, &errNum);
    }
  }
  /* Set section values from tiles. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		ln;
    size_t	width;
    WlzGreyP	dstGP;
    WlzIVertex3	ntw;

    dstGP = dstObj->values.r->values;
    width = subBox.xMax - subBox.xMin + 1;
/*
 * Throughout this function WLZ_FAST_CODE assumes that the tiles are
 * 16 x 16 x 16 voxels.
 */
#ifdef WLZ_FAST_CODE
    ntw.vtX = tv->nIdx[0] << 4;
    ntw.vtY = tv->nIdx[1] << 4;
    ntw.vtZ = tv->nIdx[2] << 4;
#else
    ntw.vtX = tv->nIdx[0] * tv->tileWidth;
    ntw.vtY = tv->nIdx[1] * tv->tileWidth;
    ntw.vtZ = tv->nIdx[2] * tv->tileWidth;
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
    for(ln = subBox.yMin; ln <= subBox.yMax; ++ln)
    {
      int         y,
                  kl,
                  lnRel;
      WlzDVertex3 vty;
      WlzGreyP	  lnGP;

      lnRel = ln - subBox.yMin;
      switch(gType)
      {
	case WLZ_GREY_INT:
	  lnGP.inp = dstGP.inp + (width * lnRel);
	  break;
	case WLZ_GREY_SHORT:
	  lnGP.shp = dstGP.shp + (width * lnRel);
	  break;
	case WLZ_GREY_UBYTE:
	  lnGP.ubp = dstGP.ubp + (width * lnRel);
	  break;
	case WLZ_GREY_FLOAT:
	  lnGP.flp = dstGP.flp + (width * lnRel);
	  break;
	case WLZ_GREY_DOUBLE:
	  lnGP.dbp = dstGP.dbp + (width * lnRel);
	  break;
	case WLZ_GREY_RGBA:
	  lnGP.rgbp = dstGP.rgbp + (width * lnRel);
	  break;
	default:
	  lnGP.v = NULL; 		    /* Just to keep compilers happy. */
	  break;
      }
      y = ln - WLZ_NINT(view->minvals.vtY);
      vty.vtX = view->yp_to_x[y];
      vty.vtY = view->yp_to_y[y];
      vty.vtZ = view->yp_to_z[y];
      WlzValueSetGrey(lnGP, 0, bgd.v, gType, width);
      for(kl = subBox.xMin; kl <= subBox.xMax; ++kl)
      {
	int   		x,
			pz,
			klRel;
	double	        tz;

	klRel = kl - subBox.xMin;
	x = kl - WLZ_NINT(view->minvals.vtX);
	tz = view->xp_to_z[x] + vty.vtZ;
	pz = WLZ_NINT(tz) - tv->plane1;
#ifdef WLZ_FAST_CODE
	if((unsigned int)(pz) < (unsigned int)(ntw.vtZ))
#else
	if((pz >= 0) && (pz < ntw.vtZ))
#endif
	{
	  int		iz,
	                py;
	  double	ty;

#ifdef WLZ_FAST_CODE
          iz = pz >> 4;
#else
	  iz = pz / tv->tileWidth;
#endif
	  ty = view->xp_to_y[x] + vty.vtY;
	  py = WLZ_NINT(ty) - tv->line1;
#ifdef WLZ_FAST_CODE
	  if((unsigned int)(py) < (unsigned int)(ntw.vtY))
#else
	  if((py >= 0) && (py < ntw.vtY))
#endif
	  {
	    int 	iy,
	                px;
	    double	tx;

#ifdef WLZ_FAST_CODE
	    iy = py >> 4;
#else
	    iy = py / tv->tileWidth;
#endif
	    tx = view->xp_to_x[x] + vty.vtX;
	    px = WLZ_NINT(tx) - tv->kol1;
#ifdef WLZ_FAST_CODE
	  if((unsigned int)(px) < (unsigned int)(ntw.vtX))
#else
	    if((px >= 0) && (px < ntw.vtX))
#endif
	    {
	      size_t	ix,
	      		idx;

#ifdef WLZ_FAST_CODE
              ix = px >> 4;
#else
	      ix = px / tv->tileWidth;
#endif
	      idx = ((iz * tv->nIdx[1] + iy) * tv->nIdx[0]) + ix;
	      idx = *(tv->indices + idx);
#ifdef WLZ_FAST_CODE
              if((unsigned int)(idx) < (unsigned int)(tv->numTiles))
#else
	      if((idx >= 0) && (idx < tv->numTiles))
#endif
	      {
		size_t   off;

#ifdef WLZ_FAST_CODE
		off = (idx  << 12) +
		      ((((pz & 15) << 4) + (py & 15)) << 4) + (px & 15);
#else
		WlzIVertex3 tOff;

		tOff.vtX = px % tv->tileWidth;
		tOff.vtY = py % tv->tileWidth;
		tOff.vtZ = pz % tv->tileWidth;
		off = (idx * tv->tileSz) +
		      ((tOff.vtZ * tv->tileWidth + tOff.vtY) *
		       tv->tileWidth) + tOff.vtX;
#endif
		switch(gType)
		{
	          case WLZ_GREY_INT:
		    lnGP.inp[klRel] = tv->tiles.inp[off];
		    break;
	          case WLZ_GREY_SHORT:
		    lnGP.shp[klRel] = tv->tiles.shp[off];
		    break;
	          case WLZ_GREY_UBYTE:
		    lnGP.ubp[klRel] = tv->tiles.ubp[off];
		    break;
	          case WLZ_GREY_FLOAT:
		    lnGP.flp[klRel] = tv->tiles.flp[off];
		    break;
	          case WLZ_GREY_DOUBLE:
		    lnGP.dbp[klRel] = tv->tiles.dbp[off];
		    break;
	          case WLZ_GREY_RGBA:
		    lnGP.rgbp[klRel] = tv->tiles.rgbp[off];
		    break;
		  default:
		    break;
		}
	      }
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
  return(dstObj);
}

#define WLZ_GETSUBSEC_POS(P,V,X,Y) \
  (P).vtX = (V)->xp_to_x[(X)] + (Y).vtX; \
  (P).vtY = (V)->xp_to_y[(X)] + (Y).vtY; \
  (P).vtZ = (V)->xp_to_z[(X)] + (Y).vtZ; \

#define WLZ_GETSUBSEC_VAL(G,V,K,Y) \
{ \
  int		x; \
  WlzDVertex3   p; \
 \
  x = k - WLZ_NINT((V)->minvals.vtX); \
  WLZ_GETSUBSEC_POS(p,(V),x,(Y)) \
  WlzGreyValueGet((G), WLZ_NINT(p.vtZ), WLZ_NINT(p.vtY), WLZ_NINT(p.vtX)); \
}

#define WLZ_GETSUBSEC_CONVAL(G,F0,F1,V,K,Y) \
{ \
  int		x; \
  WlzDVertex3   p; \
 \
  x = k - WLZ_NINT((V)->minvals.vtX); \
  WLZ_GETSUBSEC_POS(p,(V),x,(Y)) \
  WlzGreyValueGetCon((G), p.vtZ, p.vtY, p.vtX); \
  (F0).vtX = p.vtX - WLZ_NINT(p.vtX - 0.5); \
  (F0).vtY = p.vtY - WLZ_NINT(p.vtY - 0.5); \
  (F0).vtZ = p.vtZ - WLZ_NINT(p.vtZ - 0.5); \
  (F1).vtX = 1.0 - (F0).vtX; \
  (F1).vtY = 1.0 - (F0).vtY; \
  (F1).vtZ = 1.0 - (F0).vtZ; \
}

/*!
* \return	New sub-section object.
* \ingroup	WlzSectionTransform
* \brief	Computes a section through the given 3D domain object.
* \param	obj			Given 3D object.
* \param	subDomain		Given 2D domain within which to
* 					restrict the section. If NULL
* 					returned section will be have a
* 					rectangular domain which is the maximum
* 					for the given object and view
* 					transform.
* \param	view			Given view transform.
* \param	interp			Interpolation, should be either
* 					WLZ_INTERPOLATION_NEAREST or
* 					WLZ_INTERPOLATION_LINEAR.
* \param	maskRtn			Destination pointer for returned
* 					domain mask.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzGetSubSectionFrom3DDomObj(
  WlzObject		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct	*viewStr,
  WlzInterpolationType	interp,
  WlzObject		**maskRtn,
  WlzErrorNum		*dstErr)
{
  WlzObject		*newObj = NULL,
			*mask = NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzGreyValueWSpace	*gVWSp = NULL;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  int			maskFlg = 0,
  			greyFlg = 0;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  newObj = NULL;
  domain.core = NULL;
  values.core = NULL;
  /* Check parameters */
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(viewStr == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(viewStr->type != WLZ_3D_VIEW_STRUCT)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(!(viewStr->initialised))
  {
    errNum = WLZ_ERR_OBJECT_DATA;
  }
  else
  {
    switch(interp)
    {
      case WLZ_INTERPOLATION_NEAREST:
      case WLZ_INTERPOLATION_LINEAR:
	break;
      default:
	errNum = WLZ_ERR_INTERPOLATION_TYPE;
	break;
    }
  }

  if(errNum == WLZ_ERR_NONE)
  {
    WlzIBox2 	subBox;


    subBox.xMin = WLZ_NINT(viewStr->minvals.vtX);
    subBox.xMax = WLZ_NINT(viewStr->maxvals.vtX);
    subBox.yMin = WLZ_NINT(viewStr->minvals.vtY);
    subBox.yMax = WLZ_NINT(viewStr->maxvals.vtY);
    /* Sort out return object requirements */
    greyFlg = (obj->values.core != NULL);
    maskFlg = (obj->values.core == NULL) || (maskRtn != NULL);
    /* Create a new return object - domain only. */
    if(subDomain)
    {
      switch(subDomain->type)
      {
	case WLZ_2D_DOMAINOBJ:
	  newObj = WlzClipObjToBox2D(subDomain, subBox, &errNum);
	  break;
	case WLZ_2D_POLYGON:
	case WLZ_BOUNDLIST:
	case WLZ_RECTANGLE:
	  errNum = WLZ_ERR_UNIMPLEMENTED;
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    else
    {
      if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					   subBox.yMin, subBox.yMax,
					   subBox.xMin, subBox.xMax,
					   &errNum)) != NULL)
      {
	newObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
	                     NULL, NULL, &errNum);
      }
    }
  }
  /* Add a value table */
  if((errNum == WLZ_ERR_NONE) && greyFlg)
  {
    WlzPixelV	pixval;

    pixval = WlzGetBackground(obj, &errNum);
    if((values.v = WlzNewValueTb(newObj,
	    WlzGreyTableType(WLZ_GREY_TAB_RECT,
	                     WlzGreyTypeFromObj(obj, NULL), NULL),
	    pixval, &errNum)))
    {
      newObj->values = WlzAssignValues(values, &errNum);
    }
  }
  /* Check if mask required */
  if((errNum == WLZ_ERR_NONE) && maskFlg)
  {
    WlzPixelV	pixval;

    pixval.type = WLZ_GREY_UBYTE;
    pixval.v.ubv = (WlzUByte )0;
    if((values.v = WlzNewValueTb(newObj,
	    WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_UBYTE, NULL),
	    pixval, &errNum)))
    {
      mask = WlzMakeMain(WLZ_2D_DOMAINOBJ, newObj->domain, values, NULL, NULL,
	                 &errNum);
    }

    /* newObj is redundant if mask only required */
    if(!greyFlg)
    {
      WlzFreeObj(newObj);
      newObj = NULL;
    }
  }
  /* Scan object setting values */
  if((errNum == WLZ_ERR_NONE) && greyFlg)
  {
    if((gVWSp = WlzGreyValueMakeWSp(obj, &errNum)))
    {
      errNum = WlzInitGreyScan(newObj, &iwsp, &gwsp);
      while((errNum == WLZ_ERR_NONE) && 
	    ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
      {
        int 	k,
		yp;
	WlzDVertex3 vty;

	yp = iwsp.linpos - WLZ_NINT(viewStr->minvals.vtY);
	vty.vtX = viewStr->yp_to_x[yp];
	vty.vtY = viewStr->yp_to_y[yp];
	vty.vtZ = viewStr->yp_to_z[yp];
	switch(interp)
	{
	  case WLZ_INTERPOLATION_NEAREST:
	    switch(gVWSp->gType){
	      case WLZ_GREY_INT:
	        for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		{
                  WLZ_GETSUBSEC_VAL(gVWSp, viewStr, k, vty)
	          *(gwsp.u_grintptr.inp)++ = gVWSp->gVal[0].inv;
		}
		break;
	      case WLZ_GREY_SHORT:
	        for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++)
		{
                  WLZ_GETSUBSEC_VAL(gVWSp, viewStr, k, vty)
	          *(gwsp.u_grintptr.shp)++ = gVWSp->gVal[0].shv;
		}
		break;
	      case WLZ_GREY_UBYTE:
	        for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		{
                  WLZ_GETSUBSEC_VAL(gVWSp, viewStr, k, vty)
	          *(gwsp.u_grintptr.ubp)++ = gVWSp->gVal[0].ubv;
		}
		break;
	      case WLZ_GREY_FLOAT:
	        for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		{
                  WLZ_GETSUBSEC_VAL(gVWSp, viewStr, k, vty)
	          *(gwsp.u_grintptr.flp)++ = gVWSp->gVal[0].flv;
		}
		break;
	      case WLZ_GREY_DOUBLE:
	        for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		{
                  WLZ_GETSUBSEC_VAL(gVWSp, viewStr, k, vty)
	          *(gwsp.u_grintptr.dbp)++ = gVWSp->gVal[0].dbv;
		}
		break;
	      case WLZ_GREY_RGBA:
	        for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		{
                  WLZ_GETSUBSEC_VAL(gVWSp, viewStr, k, vty)
	          *(gwsp.u_grintptr.rgbp)++ = gVWSp->gVal[0].rgbv;
		}
		break;
	      default:
		break;
	    }
	    break;
	  case WLZ_INTERPOLATION_LINEAR:
	    {
	      double		tD0;
	      WlzDVertex3	tDV0,
				tDV1;

	      switch(gVWSp->gType){
		case WLZ_GREY_INT:
		  for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		  {
		    WLZ_GETSUBSEC_CONVAL(gVWSp, tDV0, tDV1, viewStr, k, vty)
		    tD0 =
		      ((gVWSp->gVal[0]).inv * tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).inv * tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).inv * tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).inv * tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).inv * tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).inv * tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).inv * tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).inv * tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    tD0 = WLZ_CLAMP(tD0, INT_MIN, INT_MAX);
		    *(gwsp.u_grintptr.inp)++ = WLZ_NINT(tD0);
		  }
		  break;
		case WLZ_GREY_SHORT:
		  for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		  {
		    WLZ_GETSUBSEC_CONVAL(gVWSp, tDV0, tDV1, viewStr, k, vty)
		    tD0 =
		      ((gVWSp->gVal[0]).shv * tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).shv * tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).shv * tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).shv * tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).shv * tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).shv * tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).shv * tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).shv * tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    tD0 = WLZ_CLAMP(tD0, SHRT_MIN, SHRT_MAX);
		    *(gwsp.u_grintptr.shp)++ = WLZ_NINT(tD0);
		  }
		  break;
		case WLZ_GREY_UBYTE:
		  for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		  {
		    WLZ_GETSUBSEC_CONVAL(gVWSp, tDV0, tDV1, viewStr, k, vty)
		    tD0 =
		      ((gVWSp->gVal[0]).ubv * tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).ubv * tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).ubv * tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).ubv * tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).ubv * tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).ubv * tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).ubv * tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).ubv * tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    tD0 = WLZ_CLAMP(tD0, 0, 255);
		    *(gwsp.u_grintptr.ubp)++ = WLZ_NINT(tD0);
		  }
		  break;
		case WLZ_GREY_FLOAT:
		  for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		  {
		    WLZ_GETSUBSEC_CONVAL(gVWSp, tDV0, tDV1, viewStr, k, vty)
		    tD0 =
		      ((gVWSp->gVal[0]).flv * tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).flv * tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).flv * tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).flv * tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).flv * tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).flv * tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).flv * tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).flv * tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    *(gwsp.u_grintptr.flp)++ = WLZ_CLAMP(tD0, FLT_MIN, FLT_MAX);
		  }
		  break;
		case WLZ_GREY_DOUBLE:
		  for(k = iwsp.lftpos; k <= iwsp.rgtpos; ++k)
		  {
		    WLZ_GETSUBSEC_CONVAL(gVWSp, tDV0, tDV1, viewStr, k, vty)
		    tD0 =
		      ((gVWSp->gVal[0]).dbv * tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[1]).dbv * tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[2]).dbv * tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[3]).dbv * tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		      ((gVWSp->gVal[4]).dbv * tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[5]).dbv * tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[6]).dbv * tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		      ((gVWSp->gVal[7]).dbv * tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		    *(gwsp.u_grintptr.dbp)++ = WLZ_NINT(tD0);
		  }
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	}
      }
      if(errNum == WLZ_ERR_EOO)	   /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
      WlzGreyValueFreeWSp(gVWSp);
    }
  }

  /* Check if mask required */
  if((errNum == WLZ_ERR_NONE) && maskFlg)
  {
    errNum = WlzInitGreyScan(mask, &iwsp, &gwsp);
    while((errNum == WLZ_ERR_NONE) && 
	  ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
    {
      int 	k,
      		yp;
      WlzDVertex3 vty;

      yp = iwsp.linpos - WLZ_NINT(viewStr->minvals.vtY);
      vty.vtX = viewStr->yp_to_x[yp];
      vty.vtY = viewStr->yp_to_y[yp];
      vty.vtZ = viewStr->yp_to_z[yp];
      for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++)
      {
	int	xp;
        WlzDVertex3 vtx;

	xp = k - WLZ_NINT(viewStr->minvals.vtX);
        WLZ_GETSUBSEC_POS(vtx, viewStr, xp, vty)
	*(gwsp.u_grintptr.ubp)++ = WlzInsideDomain(obj,
	    WLZ_NINT(vtx.vtZ), WLZ_NINT(vtx.vtY), WLZ_NINT(vtx.vtX), NULL);
      }
    }
    if(errNum == WLZ_ERR_EOO)	   /* Reset EOO error from end of intervals */ 
    {
      errNum = WLZ_ERR_NONE;
    }
    /* Threshold to determine the mask */
    if(errNum == WLZ_ERR_NONE)
    {
      WlzPixelV	pixval;
      WlzObject *tmpObj;

      pixval.type = WLZ_GREY_INT;
      pixval.v.inv = 1;
      mask = WlzAssignObject(mask, NULL);
      tmpObj = WlzThreshold(mask, pixval, WLZ_THRESH_HIGH, &errNum);
      WlzFreeObj(mask);
      values.core = NULL;
      mask = WlzMakeMain(tmpObj->type, tmpObj->domain, values,
			 NULL, NULL, &errNum);
      WlzFreeObj(tmpObj);

      /* Set returns */
      if(!greyFlg)
      {
	/* New object if mask is returned twice */
	if(maskRtn)
	{
	  *maskRtn = WlzMakeMain(mask->type, mask->domain, mask->values,
			 NULL, NULL, &errNum);
	}
	newObj = mask;
      }
      else if(maskRtn)
      {
	*maskRtn = mask;
      }
      else
      {
	/* Should not get here */
	WlzFreeObj(mask);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newObj);
}
