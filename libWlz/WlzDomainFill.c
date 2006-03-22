#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDomainFill_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzDomainFill.c
* \author       Richard Baldock
* \date         September 2003
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
* \brief	Functions to fill holes in domain objects.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>


/* function:     WlzDomainFill    */
/*! 
* \ingroup      WlzDomainOps
* \brief        fill holes in a woolz domain object domain. The returned
 object will have a NULL valuetable.
*
* \return       Domain object with holes filled.
* \param    obj	Input domain object.
* \param    dstErr	Error return
* \par      Source:
*                WlzDomainFill.c
*/
WlzObject *WlzDomainFill(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL, *obj1, *obj2;
  WlzBoundList	*bndList;
  WlzDomain	domain;
  WlzValues	values;
  WlzPlaneDomain	*pdom, *rtnpdom;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		p;

  /* check the object pointer and type */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      /* check domain */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( obj->domain.core->type ){
	case WLZ_INTERVALDOMAIN_INTVL:
	  if( obj1 = WlzObjToBoundary(obj, 1, &errNum) ){
	    bndList = obj1->domain.b;
	    while( bndList != NULL ){
	      if( bndList->down ){
		WlzFreeBoundList(bndList->down);
		bndList->down = NULL;
	      }
	      bndList = bndList->next;
	    }
	    rtnObj = WlzBoundToObj(obj1->domain.b, WLZ_SIMPLE_FILL, &errNum);
	    WlzFreeObj(obj1);
	  }
	  break;

	case WLZ_INTERVALDOMAIN_RECT:
	  values.core = NULL;
	  return WlzMakeMain(obj->type, obj->domain, values,
			     NULL, NULL, dstErr);

	case WLZ_EMPTY_DOMAIN:
	  return WlzMakeEmpty(dstErr);

	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
	}
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      switch( obj->domain.core->type ){
      case WLZ_PLANEDOMAIN_DOMAIN:
	pdom = obj->domain.p;
	if( !(rtnpdom = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					   pdom->plane1, pdom->lastpl,
					   pdom->line1, pdom->lastln,
					   pdom->kol1, pdom->lastkl,
					   &errNum)) ){
	  break;
	}
	rtnpdom->voxel_size[0] = pdom->voxel_size[0];
	rtnpdom->voxel_size[1] = pdom->voxel_size[1];
	rtnpdom->voxel_size[2] = pdom->voxel_size[2];
	values.core = NULL;
	for(p=pdom->plane1; p <= pdom->lastpl; p++){
	  if( domain.i = (pdom->domains)[p - pdom->plane1].i ){
	    obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			       NULL, NULL, NULL);
	    if( obj2 = WlzDomainFill(obj1, &errNum) ){
	      rtnpdom->domains[p - pdom->plane1] =
		WlzAssignDomain(obj2->domain, NULL);
	      WlzFreeObj(obj2);
	    }
	    WlzFreeObj(obj1);
	  }
	  else {
	    rtnpdom->domains[p - pdom->plane1].core = NULL;
	  }
	}
	domain.p = rtnpdom;
	return WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			   NULL, NULL, dstErr);

      case WLZ_EMPTY_DOMAIN:
	return WlzMakeEmpty(dstErr);

      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      }
      break;

    case WLZ_TRANS_OBJ:
      rtnObj = WlzDomainFill(obj->values.obj, &errNum);
      if( errNum == WLZ_ERR_NONE ){
	values.obj = rtnObj;
	return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( dstErr ){
    *dstErr=errNum;
  }
  return rtnObj;
}
