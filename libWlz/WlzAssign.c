#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzAssign_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzAssign.c
* \author       Bill Hill, Richard Baldock, Christophe Dubreuil
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
* \brief	Woolz objects domains and values maintain a linkcount,
* 		which records it's usage by other objects, domains or
* 		values. To increment a linkcount the appropriate assignment
* 		function should be used.
* \ingroup	WlzAllocation
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

/*!
* \return	Given object with incremented linkcount or NULL on error.
* \ingroup	WlzAllocation
* \brief	Assign an object (increment it's linkcount) by first
* 		checking for NULL, then the value of linkcount, before
*		incrementing the linkcount. If used concientiously,
*		assignment should avoid memory errors.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzAssignObject(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj)
  {
    if(obj->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnObj = obj;
      obj->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return rtnObj;
}

/*!
* \return	Given domain with incremented linkcount or NULL on error.
* \ingroup	WlzAllocation
* \brief	Assign a domain by incrementing it's linkcount.
* \param	domain			Given domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDomain 	WlzAssignDomain(WlzDomain domain, WlzErrorNum *dstErr)
{
  WlzDomain	rtnDomain;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  rtnDomain.core = NULL;
  if(domain.core)
  {
    if(domain.core->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnDomain = domain;
      domain.core->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return rtnDomain;
}

/*!
* \return	Given values with incremented linkcount or NULL on error.
* \ingroup	WlzAllocation
* \brief	Assign a values by incrementing it's linkcount.
* \param	values			Given values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzValues 	WlzAssignValues(WlzValues values, WlzErrorNum	*dstErr)
{
  WlzValues	rtnValues;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rtnValues.core = NULL;
  if(values.core)
  {
    if(values.core->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnValues = values;
      values.core->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return rtnValues;
}

/*!
* \return	Given property with incremented linkcount or the core
*		property set to NULL on error.
* \ingroup	WlzAllocation
* \brief	Assign a property by incrementing it's linkcount.
* \param	property		Given property.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzProperty WlzAssignProperty(WlzProperty property, WlzErrorNum *dstErr)
{
  WlzProperty	rtnProp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rtnProp.core = NULL;
  if(property.core)
  {
    if(property.core->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnProp = property;
      property.core->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return rtnProp;
}

/*!
* \return	Property list with incremented link count or NULL on error.
* \brief	Assigned a Woolz property list, incrementing the link count
		or the number of times the property list is used.
* \param	pList			Given property list.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPropertyList *WlzAssignPropertyList(WlzPropertyList *pList,
				       WlzErrorNum *dstErr)
{
  WlzPropertyList *rtnPList = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(pList)
  {
    if(pList->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnPList = pList;
      pList->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnPList);
}

/*!
* \return	Given affine transform with incremented linkcount or NULL on
*		error.
* \ingroup      WlzAllocation
* \brief	Assign an affine transform by incrementing it's linkcount.
* \param	trans			Given affine transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzAffineTransform *WlzAssignAffineTransform(
  WlzAffineTransform *trans,
  WlzErrorNum	*dstErr)
{
  WlzAffineTransform	*rtnTrans=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  if(trans)
  {
    if(trans->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnTrans = trans;
      ++(trans->linkcount);
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return rtnTrans;
}

/*!
* \return	Given boundary list with incremented linkcount or NULL on
*		error.
* \ingroup	WlzAllocation
* \brief	Assign a boundary list by incrementing it's linkcount.
* \param	blist			Given boundary list.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBoundList 	*WlzAssignBoundList(WlzBoundList *blist, WlzErrorNum *dstErr)
{
  WlzBoundList	*rtnBlist = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(blist)
  {
    if(blist->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnBlist = blist;
      ++(blist->linkcount);
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return rtnBlist;
}

/*!
* \return	Given polygon domain with incremented linkcount or NULL on
*		error.
* \ingroup	WlzAllocation
* \brief	Assign a polygon domain by incrementing it's linkcount.
* \param	poly			Given boundary list.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPolygonDomain *WlzAssignPolygonDomain(
  WlzPolygonDomain 	*poly,
  WlzErrorNum	*dstErr)
{
  WlzPolygonDomain	*rtnPoly=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if(poly)
  {
    if(poly->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnPoly = poly;
      ++(poly->linkcount);
    }
  }
#ifdef WLZ_NO_NULL
  else
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return rtnPoly;
}

/*!
* \return	Given geometric model with incremented linkcount or NULL on
* 		error.
* \ingroup      WlzAllocation
* \brief	Assign a geometric model by incrementing it's linkcount.
* \param	model			Given geometric model.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGMModel 	*WlzAssignGMModel(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMModel	*rtnModel = NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if(model)
  {
    if(model->linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      rtnModel = model;
      ++(model->linkcount);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnModel);
}

/*!
* \return	Non-zero if object can be free'd.
* \ingroup      WlzAllocation
* \brief	Unlink an object, domain or values by decrementing
*		and testing it's linkcount.
* \param	linkcount		Given linkcount pointer.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzUnlink(int *linkcount, WlzErrorNum *dstErr)
{
  int		canFree = 0;
  WlzErrorNum	errNum = WLZ_ERR_PARAM_NULL;

  if(linkcount)
  {
    if(*linkcount < 0)
    {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else
    {
      errNum = WLZ_ERR_NONE;
      if(--*linkcount <= 0)
      {
	*linkcount = -1;
	canFree = 1;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(canFree);
}
