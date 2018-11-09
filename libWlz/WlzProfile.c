#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzProfile_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzProfile.c
* \author       Bill Hill
* \date         September 2018
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2018],
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
* \brief	Functions for extracting a profile (with or without)
* 		grey values along a ray through an object.
* \ingroup	WlzTransform
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

typedef struct _WlzProfileWalkWSp
{
  int		index;
  int		inside;
  WlzIVertex3	start;
  WlzIVertex3	end;
  WlzDomain 	domain;
  WlzValues 	values;
  WlzGreyValueWSpace *gVWSp;
} WlzProfileWalkWSp;

typedef WlzErrorNum (*WlzProfileWalkFn)(
  WlzObject *, WlzIVertex3, WlzProfileWalkWSp *);


static WlzErrorNum 		WlzProfileSetIDom(
				  WlzObject *obj,
				  WlzIVertex3 pos,
				  WlzProfileWalkWSp *wSp);
static WlzErrorNum 		WlzProfileWalk(
				  WlzObject *obj,
				  WlzProfileWalkFn fn,
				  WlzProfileWalkWSp *wSp);
static WlzObject 		*WlzProfileLnSD(
  				  WlzObject *gObj,
				  WlzIVertex3 sPos,
				  WlzIVertex3 ePos,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzProfileWalkValues(
				  WlzObject *gObj,
				  WlzObject *rObj,
				  WlzGreyValueWSpace *gVWSp,
				  WlzIVertex3 sPos,
				  WlzIVertex3 ePos,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzProfileSetPixels(
				  WlzObject *obj,
				  WlzIVertex3 pos,
				  WlzProfileWalkWSp *pWSp);

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzTransform
* \brief	Creates a new 2D Woolz object with a single line which
* 		is the profile from the given start position to the
* 		given end position.
* 		Empty objects are allowed and in this case a new empty object
* 		is returned. In other cases (where there is no error) the
* 		returned object will always be a 2D spatial domain object
* 		(WLZ_2D_DOMAINOBJ) which has an interval domain and where
* 		the given object has values a rectangular value table.
* 		These will have a single line with the pixel at 0,0.
* 		When the line segment between the given points does not
* 		intersect the given object a spatial domain object is
* 		still returned but it will have no intervals and it
* 		may be necessary to call WlzIsEmpty().
* \param	gObj		Given Woolz object.
* \param	sPos		Start position.
* \param	ePos		End position.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject *WlzProfileLine(
  WlzObject *gObj,
  WlzIVertex3 sPos,
  WlzIVertex3 ePos,
  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_EMPTY_OBJ:
        rObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	if(gObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
          rObj = WlzProfileLnSD(gObj, sPos, ePos, &errNum);
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New Woolz interval domain or NULL on error.
* \ingroup	WlzTransform
* \brief	Creates a new 2D Woolz interval domain with a single line
* 		which is the profile from the given start position to the
* 		given end position for the given spatial domain object.
* \param	gObj		Given Woolz object.
* \param	sPos		Start position.
* \param	ePos		End position.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzIntervalDomain *WlzProfileLineIDom(
  WlzObject *gObj,
  WlzIVertex3 sPos,
  WlzIVertex3 ePos,
  WlzErrorNum *dstErr)
{
  WlzIntervalDomain *iDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check given object. */
  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
	if(gObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		len;
    WlzIVertex3	dPos;
    WlzProfileWalkWSp pWSp = {0};
    WlzInterval *rItv = NULL;

    WLZ_VTX_3_SUB(dPos, ePos, sPos);
    len = (int )ceil(WLZ_VTX_3_LENGTH(dPos));
    /* Create return object. */
    iDom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
                                 0, 0, 0, len - 1, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if((rItv = (WlzInterval *)
		 AlcCalloc((len / 2) + 1, sizeof(WlzInterval))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
	(void )WlzFreeIntervalDomain(iDom);
	iDom = NULL;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzIntervalLine *itvLn;

      iDom->freeptr = AlcFreeStackPush(iDom->freeptr, (void *)rItv, NULL);
      pWSp.start = sPos;
      pWSp.end = ePos;
      pWSp.domain.i = iDom;
      itvLn = pWSp.domain.i->intvlines + 0;
      itvLn->intvs = rItv;
      itvLn->nintvs = 0;
      errNum = WlzProfileWalk(gObj, WlzProfileSetIDom, &pWSp);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    iDom = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(iDom);
}

/*!
* \return	New 2D Woolz object or NULL on error.
* \ingroup	WlzTransform
* \brief	Creates a new 2D Woolz object with a single line which
* 		is the profile from the given start position to the
* 		given end position. This function assumes that it is
* 		given either a 2 or 3D spatial domain object and that
* 		this object has a non-NULL domain. See WlzProfileLine().
* \param	gObj		Given Woolz object.
* \param	sPos		Start position.
* \param	ePos		End position.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static WlzObject *WlzProfileLnSD(
  WlzObject *gObj,
  WlzIVertex3 sPos,
  WlzIVertex3 ePos,
  WlzErrorNum *dstErr)
{
  WlzDomain	rDom = {0};
  WlzValues	rVal = {0};
  WlzObject	*rObj = NULL;
  WlzGreyP	rPix = {0};
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rDom.i = WlzProfileLineIDom(gObj, sPos, ePos, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    if(gObj->values.core)
    {
      gVWSp = WlzGreyValueMakeWSp(gObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	int	len;
	size_t	gSz;

	len = rDom.i->lastkl - rDom.i->kol1 +1;
	if((gSz = WlzGreySize(gVWSp->gType)) == 0)
	{
	  errNum = WLZ_ERR_GREY_TYPE;
	}
	else if((rPix.v = AlcMalloc(len * gSz)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	WlzPixelV bgd;
	WlzObjectType tType;

	bgd.type = gVWSp->gType;
	bgd.v = gVWSp->gBkd;
	tType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RECT, gVWSp->gType,
	                              NULL);
	rVal.r = WlzMakeRectValueTb(tType,
	                            rDom.i->line1, rDom.i->lastln, rDom.i->kol1,
				    (rDom.i->lastkl - rDom.i->kol1) + 1,
				    bgd, rPix.inp, &errNum);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(gObj->type, rDom, rVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzProfileWalkValues(gObj, rObj, gVWSp, sPos, ePos, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else
    {
      if(rVal.core)
      {
	(void )WlzFreeDomain(rDom);
	(void )WlzFreeValues(rVal);
      }
      else
      {
        AlcFree(rPix.v);
      }
    }
  }
  if(dstErr)    
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Sets the values in a profiles previously defined valuetable.
* \param	gObj		Given object.
* \param	rObj		Profile object with values to set.
* \param	gVWSp		Grey workspace initialised for the given object
* 				through which the profile runs.
* \param	sPos		Start position.
* \param	ePos		End position.
* \param	dstErr		Destination error pointer, my be NULL.
*/
static WlzErrorNum	WlzProfileWalkValues(
  WlzObject *gObj,
  WlzObject *rObj,
  WlzGreyValueWSpace *gVWSp,
  WlzIVertex3 sPos,
  WlzIVertex3 ePos,
  WlzErrorNum *dstErr)
{
  WlzProfileWalkWSp pWSp = {0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  pWSp.start = sPos;
  pWSp.end = ePos;  
  pWSp.gVWSp = gVWSp;
  pWSp.values = rObj->values;
  errNum = WlzProfileWalk(gObj, WlzProfileSetPixels, &pWSp);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Sets the interval for the profile for the given point
* 		with respect to the given object.
* \param	obj		Given object.
* \param	pos		Given point with respect to the given
* 				object.
* \param	wSp		profile walking work space.
*/
static WlzErrorNum WlzProfileSetIDom(
  WlzObject *obj,
  WlzIVertex3 pos,
  WlzProfileWalkWSp *pWSp)
{
  int		msk;
  WlzInterval	*itv;
  WlzIntervalLine *itvLn;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  msk = (pWSp->inside << 1) |
        (WlzInsideDomain(obj, pos.vtZ, pos.vtY, pos.vtX, &errNum) != 0);
  if(errNum == WLZ_ERR_NONE)
  {
    itvLn = pWSp->domain.i->intvlines + 0;
    switch(msk)
    {
      case 1:  /* Was out, now in. Start new interval. */
	itv = itvLn->intvs + itvLn->nintvs;
	itv->ileft = itv->iright = pWSp->index;
	pWSp->inside = 1;
	++(itvLn->nintvs);
	break;
      case 3:  /* Was in,  still in. Just update interval. */
	itv = itvLn->intvs + itvLn->nintvs - 1;
	itv->iright = pWSp->index;
	break;
      default: /* Either was in, now out or was out, still out. */
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Sets the profile value for walk function.
* \param	obj		Given object from which to set the value.
* \param	pos		Position of value in given object.
* \param	pWSp		Profile workspace.
*/
static WlzErrorNum WlzProfileSetPixels(
  WlzObject *obj,
  WlzIVertex3 pos,
  WlzProfileWalkWSp *pWSp)
{
  WlzGreyP	gP;
  WlzGreyV	gV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WlzGreyValueGet(pWSp->gVWSp, pos.vtZ, pos.vtY, pos.vtX);
  gP = pWSp->values.r->values;
  gV = pWSp->gVWSp->gVal[0];
  switch(pWSp->gVWSp->gType)
  {
    case WLZ_GREY_LONG:
      *(gP.lnp + pWSp->index) = gV.lnv;
      break;
    case WLZ_GREY_INT:
      *(gP.inp + pWSp->index) = gV.inv;
      break;
    case WLZ_GREY_SHORT:
      *(gP.shp + pWSp->index) = gV.shv;
      break;
    case WLZ_GREY_UBYTE:
      *(gP.ubp + pWSp->index) = gV.ubv;
      break;
    case WLZ_GREY_FLOAT:
      *(gP.flp + pWSp->index) = gV.flv;
      break;
    case WLZ_GREY_DOUBLE:
      *(gP.dbp + pWSp->index) = gV.dbv;
      break;
    case WLZ_GREY_RGBA:
      *(gP.rgbp + pWSp->index) = gV.rgbv;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Walks along a line segment calling the given function.
* \param	obj		Object passed to function.
* \param	fn		Function to call.
* \param	pWSp		Workspace in which the start and end positions,
* 				a valid interval domain and sufficent
* 				intervals have been set.
*/
static WlzErrorNum WlzProfileWalk(
  WlzObject *obj,
  WlzProfileWalkFn fn,
  WlzProfileWalkWSp *pWSp)
{
  int		i,
  		m;
  WlzIVertex3	d,
  		e,
		p,
  		s;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_VTX_3_SUB(d, pWSp->end, pWSp->start);
  WLZ_VTX_3_SIGN(s, d);
  WLZ_VTX_3_ABS(d, d);
  m = ALG_MAX3(d.vtX, d.vtY, d.vtZ);
  WLZ_VTX_3_SET(e, m / 2, m / 2, m / 2);
  i = m;
  pWSp->index = 0;
  pWSp->inside = 0;
  p = pWSp->start;
  while(1)
  {
    WlzIVertex3 t0,
    		t1;

    errNum = (*fn)(obj, p, pWSp);
    if((errNum != WLZ_ERR_NONE) || (--i < 0))
    {
      break;
    }
    ++(pWSp->index);
    WLZ_VTX_3_SUB(e, e, d);
    WLZ_VTX_3_SET(t0, (e.vtX < 0), (e.vtY < 0), (e.vtZ < 0));
    WLZ_VTX_3_SCALE_ADD(e, t0, m, e);
    WLZ_VTX_3_HAD(t1, t0, s);
    WLZ_VTX_3_ADD(p, p, t1);
  }
  return(errNum);
}
