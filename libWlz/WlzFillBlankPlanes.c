#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzFillBlankPlanes.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Fills blank planes of a 3D object. Originaly for
*		MAPaint to allow painting of intermediate planes.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>

#include <Wlz.h>

WlzErrorNum WlzFillBlankPlanes(
  WlzObject	*obj,
  int		min_domain)
{
  WlzObject		*fill_obj, *obj1, *obj2;
  WlzPlaneDomain	*planedmn;
  WlzVoxelValues	*voxtab;
  WlzDomain		domain, *domains;
  WlzValues		values, *valuess;
  WlzObjectType		valTbType;
  WlzGreyType		gType;
  int			p;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    return WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.p == NULL ){
	return WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( obj->domain.p->type ){

	case WLZ_PLANEDOMAIN_DOMAIN:
	  planedmn = obj->domain.p;
	  domains = planedmn->domains;
	  valuess = NULL;
	  if( obj->values.vox ){
	    valuess = obj->values.vox->values;
	  }
	  break;

	default:
	  return WLZ_ERR_DOMAIN_TYPE;
	}
      }
      break;

    case WLZ_EMPTY_OBJ:
      return errNum;

    default:
      return WLZ_ERR_OBJECT_TYPE;
    }
  }

  /* make a domain to fill each plane, two options:
     min_domain = 0  :- use bounding box
     min_domain != 0 :- use union of non-blank domains
     */
  fill_obj = NULL;
  values.core = NULL;
  if( min_domain ){
    for(p=planedmn->plane1; p <= planedmn->lastpl; p++){
      if( domains[p - planedmn->plane1].core == NULL ){
	continue;
      }

      if( fill_obj == NULL ){
	fill_obj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			       domains[p - planedmn->plane1],
			       values, NULL, NULL, NULL);
      }
      else {
	obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			   domains[p - planedmn->plane1],
			   values, NULL, NULL, NULL);
	obj2 = WlzUnion2(fill_obj, obj1, NULL);
	WlzFreeObj(fill_obj);
	WlzFreeObj(obj1);
	fill_obj = obj2;
      }

      if( valuess ){
	gType = WlzGreyTableTypeToGreyType(
	  valuess[p - planedmn->plane1].core->type, NULL);
	valTbType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, gType, NULL);
      }
    }
  }
  else {
    domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				     planedmn->line1, planedmn->lastln,
				     planedmn->kol1, planedmn->lastkl,
				     NULL);
    fill_obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			   NULL, NULL, NULL);
  }
				   

  /* loop through planes filling blanks */
  for(p=planedmn->plane1; p <= planedmn->lastpl; p++){
    if( domains[p - planedmn->plane1].core != NULL ){
      continue;
    }

    if( fill_obj == NULL ){
      continue;
    }

    domain.i = WlzNewIDomain(fill_obj->domain.i->type,
			     fill_obj->domain.i, NULL);
    domains[p - planedmn->plane1] = WlzAssignDomain(domain, NULL);

    if( valuess != NULL ){
      WlzPixelV background = WlzGetBackground(obj, NULL);

      values.core = NULL;
      obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			 NULL, NULL, NULL);
      
      values.v = WlzNewValueTb(obj1, valTbType, background, NULL);
      valuess[p - planedmn->plane1] = WlzAssignValues(values, NULL);
      WlzFreeObj( obj1 );
    }
  }

  return errNum;
}
