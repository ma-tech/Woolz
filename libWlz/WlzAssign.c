#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzAssign.c
* Date:         March 1999
* Author:       Christophe Dubreuil
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Assignemnt of Woolz objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

/************************************************************************
*   Function   : WlzAssignObject					*
*   Date       : Fri Oct 18 12:09:13 1996				*
*************************************************************************
*   Synopsis    : Assign an object pointer (increment linkcount) by     *
*		first checking for NULL, then value of linkcount then	*
*		incrementing the linkcount before return. If used 	*
*		concientiously this should avoid memory free conflicts	*
*   Returns    : WlzObject *: given object pointer			*
*   Parameters : WlzObject *obj: object to be assigned			*
*   Global refs: None							*
************************************************************************/

WlzObject *WlzAssignObject(
  WlzObject *obj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if (obj != (WlzObject *) NULL) {
    if (obj->linkcount < 0) {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else {
      rtnObj = obj;
      obj->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

/************************************************************************
*   Function   : WlzAssignDomain					*
*   Date       : Fri Oct 18 12:07:56 1996				*
*************************************************************************
*   Synopsis   : assign a domain pointer (increment linkcount)		*
*   Returns    : WlzDomain: the given domain				*
*   Parameters : WlzDomain domain: domain to be assigned		*
*   Global refs: None							*
************************************************************************/

WlzDomain WlzAssignDomain(
  WlzDomain domain,
  WlzErrorNum	*dstErr)
{
  WlzDomain	rtnDomain;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  rtnDomain.core = NULL;
  if (domain.core != NULL) {
    if (domain.core->linkcount < 0) {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else {
      rtnDomain = domain;
      domain.core->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
#endif /* WLZ_NO_NULL */

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnDomain;
}

/************************************************************************
*   Function   : WlzAssignValues					*
*   Date       : Fri Oct 18 12:08:17 1996				*
*************************************************************************
*   Synopsis   : assign a valutable pointer (increment linkcount)	*
*   Returns    : WlzValues: the given valuetable			*
*   Parameters : WlzValues values: valuetable to be assigned		*
*   Global refs: None							*
************************************************************************/

WlzValues WlzAssignValues(
  WlzValues values,
  WlzErrorNum	*dstErr)
{
  WlzValues	rtnValues;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  rtnValues.core = NULL;
  if (values.core != NULL) {
    if (values.core->linkcount < 0) {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else {
      rtnValues = values;
      values.core->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else {
    errNum = WLZ_ERR_VALUES_NULL;
  }
#endif /* WLZ_NO_NULL */

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnValues;
}

/************************************************************************
*   Function   : WlzAssignProperty					*
*   Date       : Fri Oct 18 12:08:43 1996				*
*************************************************************************
*   Synopsis   : assign a property pointer (increment linkcount)	*
*   Returns    : WlzSimpleProperty *: the given property		*
*   Parameters : WlzSimpleProperty *property: the property to be	*
*		 assigned						*
*   Global refs: None							*
************************************************************************/

WlzSimpleProperty *WlzAssignProperty(
  WlzSimpleProperty *property,
  WlzErrorNum	*dstErr)
{
  WlzSimpleProperty	*rtnProp=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  if (property != (WlzSimpleProperty *) NULL) {
    if (property->linkcount < 0) {
      errNum = WLZ_ERR_LINKCOUNT_DATA;
    }
    else {
      rtnProp = property;
      property->linkcount++;
    }
  }
#ifdef WLZ_NO_NULL
  else {
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
#endif /* WLZ_NO_NULL */

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnProp;
}

/************************************************************************
*   Function   : WlzAssignAffineTransform				*
*   Date       : Thu Jul 10 09:14:55 BST 1997				*
*************************************************************************
*   Synopsis   : assign an affine transform pointer (ie increment	*
*		 linkcount)						*
*   Returns    : WlzAffineTransform *: the given affine transform	*
*   Parameters : WlzAffineTransform *trans: the affine transform to be	*
*		 assigned						*
*   Global refs: None							*
************************************************************************/

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
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnTrans;
}


/************************************************************************
*   Function   : WlzAssignBoundList					*
*   Date       : Thu Jul 10 09:14:55 BST 1997				*
*************************************************************************
*   Synopsis   : assign an boundary list pointer (ie increment		*
*		 linkcount)						*
*   Returns    : WlzBoundList *: the given boundary list		*
*   Parameters : WlzBoundList *blist: the boundary list to be assigned	*
*   Global refs: None							*
************************************************************************/

WlzBoundList *WlzAssignBoundList(
  WlzBoundList 	*blist,
  WlzErrorNum	*dstErr)
{
  WlzBoundList	*rtnBlist=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

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
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnBlist;
}

/************************************************************************
*   Function   : WlzAssignPolygonDomain					*
*   Date       : Thu Jul 10 09:14:55 BST 1997				*
*************************************************************************
*   Synopsis   : assign an polygon domain pointer (ie increment		*
*		 linkcount)						*
*   Returns    : WlzPolygonDomain *: the given polygon domain		*
*   Parameters : WlzPolygonDomain *poly: the polgon to be assigned	*
*   Global refs: None							*
************************************************************************/

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
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
#endif /* WLZ_NO_NULL */

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnPoly;
}


/************************************************************************
* Function:	WlzUnlink						*
* Returns:	int:			Non-zero if object can be	*
*					free'd.				*
* Purpose:	Given a pointer to an object's linkcount: Decrements	*
*		and tests the linkcount to check if the object can be 	*
*		free'd.							*
* Global refs:	-							*
* Parameters:	int *linkcount:		Given linkcount pointer.	*
*		WlzErrorNum *dstErr:	Destination pointer for error	*
*					number.				*
************************************************************************/
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
