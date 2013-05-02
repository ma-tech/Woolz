#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzVolume_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzVolume.c
* \author       Richard Baldock
* \date         September 2003
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
* \brief	Computes the volume of domain objects.
* \ingroup	WlzDomainOps
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
WlzLong WlzVolume(
  WlzObject 	*obj,
  WlzErrorNum 	*wlzErr)
{
  WlzObject		*tmpobj;
  WlzPlaneDomain	*pldom;
  WlzDomain		domain;
  WlzValues		values;
  WlzLong		vol = -1,
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

  case WLZ_2D_DOMAINOBJ:
    return(WlzArea(obj, wlzErr));

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
