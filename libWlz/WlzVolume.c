#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzVolume.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Sep 24 17:36:44 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzDomainOps
* \brief        Compute volume of a Woolz domain object.
*               
* \todo         Should extend to include voxel size/ affine transform
 if requested.
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
* 03-02-2k bill Fix null pointer reference bug in WlzVolume().
*/
#include <Wlz.h>


/* function:     WlzVolume    */
/*! 
* \ingroup      WlzDomainOps
* \brief        Calculate the volume of the input 3D domain object.
*
* \return       Volume of input 3D object, -1 on error.
* \param    obj	Input object pointer.
* \param    wlzErr	Error return.
* \par      Source:
*                WlzVolume.c
*/
int WlzVolume(
  WlzObject 	*obj,
  WlzErrorNum 	*wlzErr)
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
    if(*wlzErr)
    {
      *wlzErr = WLZ_ERR_NONE;
    }
    return 0;

  default:
    if(wlzErr)
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
