#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzVolume.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes the volume of a Woolz domain object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

 
int		WlzVolume(WlzObject *obj, WlzErrorNum *wlzErr)
{
  WlzObject		*tmpobj;
  WlzPlaneDomain	*pldom;
  WlzDomain		domain;
  WlzValues		values;
  int			vol = -1,
  			p;
  WlzErrorNum		errNum = WLZ_ERR_NONE; 

  /* check the object */
  if( obj == NULL )
  {
    if(*wlzErr)
    {
      *wlzErr = WLZ_ERR_OBJECT_NULL;
    }
    return( -1 );
  }

  switch( obj->type ){

  case WLZ_3D_DOMAINOBJ:
    if( obj->domain.core == NULL ){
      if(*wlzErr)
      {
        *wlzErr = WLZ_ERR_DOMAIN_NULL;
      }
      return -1;
    }

    if( obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN )
    {
      if(*wlzErr)
      {
        *wlzErr = WLZ_ERR_DOMAIN_TYPE;
      }
      return -1;
    }
    break;

  case WLZ_EMPTY_OBJ:
    return 0;

  default:
    if(*wlzErr)
    {
      *wlzErr = WLZ_ERR_OBJECT_TYPE;
    }
    return( -1 );

  }

  /* scan through planes calculating the area */
  vol = 0;
  values.core = NULL;
  pldom = obj->domain.p;
  for(p=0; (p <= (pldom->lastpl - pldom->plane1)) && (errNum == WLZ_ERR_NONE);
      p++){
    domain = pldom->domains[p];
    if( domain.core ){
      tmpobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL,
      			   &errNum);
      if(tmpobj) {
	vol += WlzArea(tmpobj , &errNum);
	WlzFreeObj( tmpobj );
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    vol = -1;
  }
  if(wlzErr)
  {
    *wlzErr = errNum;
  }

  return( vol );
}
