#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzDomainFill.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Sep 24 08:35:58 2003
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
* \brief        Functions to fill holes in woolz domain objects.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
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
