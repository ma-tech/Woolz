#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzBoundToObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Boundary to interval domain conversion.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

WlzObject *WlzBoundaryToObj(
  WlzObject		*boundary,
  WlzPolyFillMode	fillMode,
  WlzErrorNum		*dstErr)
{
  WlzObject	*rtnObj=NULL, *tmpObj;
  WlzDomain	domain, *domains, *bnddmns;
  WlzValues	values;
  int		p;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object */
  if( boundary == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( boundary->type ){
    case WLZ_3D_DOMAINOBJ:
      /* check plane domain */
      if( boundary->domain.p == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( boundary->domain.p->type ){
	case WLZ_PLANEDOMAIN_BOUNDLIST:
	  if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					    boundary->domain.p->plane1,
					    boundary->domain.p->lastpl,
					    boundary->domain.p->line1,
					    boundary->domain.p->lastln,
					    boundary->domain.p->kol1,
					    boundary->domain.p->lastkl,
					    &errNum) ){
	    domain.p->voxel_size[0] = boundary->domain.p->voxel_size[0];
	    domain.p->voxel_size[1] = boundary->domain.p->voxel_size[1];
	    domain.p->voxel_size[2] = boundary->domain.p->voxel_size[2];
	    domains = domain.p->domains;
	    bnddmns = boundary->domain.p->domains;
	    for(p=domain.p->plane1; p <= domain.p->lastpl;
		p++, domains++, bnddmns++){
	      if( (*bnddmns).poly ){
		if( tmpObj = WlzBoundToObj((*bnddmns).b, fillMode,
					  &errNum) ){
		  *domains = WlzAssignDomain(tmpObj->domain, NULL);
		  WlzFreeObj(tmpObj);
		}
		else {
		  WlzFreePlaneDomain(domain.p);
		  domain.p = NULL;
		  break;
		}
	      }
	    }
	    if( domain.p ){
	      values.core = NULL;
	      if( (rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
					NULL, NULL, &errNum)) == NULL ){
		WlzFreePlaneDomain(domain.p);
	      }
	    }
	  }
	  break;
	  
	case WLZ_EMPTY_DOMAIN:
	  return WlzMakeEmpty(dstErr);

	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
	}
      }
      break;

    case WLZ_BOUNDLIST:
      return WlzBoundToObj(boundary->domain.b, fillMode, dstErr);

    case WLZ_TRANS_OBJ:
      if( values.obj = WlzBoundaryToObj(boundary->values.obj, fillMode,
				       &errNum) ){
	return WlzMakeMain(WLZ_TRANS_OBJ, boundary->domain, values,
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
    *dstErr = errNum;
  }
  return rtnObj;
}
WlzObject *WlzBoundToObj(
  WlzBoundList		*bound,
  WlzPolyFillMode	fillMode,
  WlzErrorNum		*dstNum)
{
  WlzObject *obj, *selfobj, *nextobj, *downobj;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check input */
  if( bound == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* find object corresponding to the current boundary */
  if( errNum == WLZ_ERR_NONE ){
    if( bound->poly != NULL ){
      selfobj = WlzPolyToObj(bound->poly, fillMode, &errNum);
    }
    else {
      selfobj = NULL;
    }
  }
  if( (errNum == WLZ_ERR_NONE) && (selfobj != NULL)
     && (bound->type != WLZ_BOUNDLIST_PIECE)
     && (fillMode != WLZ_VERTEX_FILL) ){
    obj = WlzErosion(selfobj, WLZ_4_CONNECTED, &errNum);
    WlzFreeObj(selfobj);
    selfobj = obj;
  }

  /* add object corresponding to next */
  if( (errNum == WLZ_ERR_NONE) && bound->next ){
    nextobj = WlzBoundToObj(bound->next, fillMode, &errNum);
    if( selfobj == NULL ){
      selfobj = nextobj;
    }
    else {
      if( (errNum == WLZ_ERR_NONE) && (nextobj != NULL) ){
	obj = WlzUnion2(selfobj, nextobj, &errNum);
	WlzFreeObj(nextobj);
	WlzFreeObj(selfobj);
	selfobj = obj;
      }
    }
  }

  /* remove object corresponding to down */
  if( (errNum == WLZ_ERR_NONE) && bound->down ){
    downobj = WlzBoundToObj(bound->down, fillMode, &errNum);
    if( downobj != NULL ){
      if( selfobj != NULL ){
	if( fillMode != WLZ_VERTEX_FILL ){
	  obj = WlzDiffDomain(selfobj, downobj, &errNum);
	}
	else {
	  obj = WlzUnion2(selfobj, downobj, &errNum);
	}
	WlzFreeObj(selfobj);
	selfobj = obj;
      }
      WlzFreeObj(downobj);
    }
  }

  /* return object */
  if( dstNum ){
    *dstNum = errNum;
  }
  return selfobj;
}
