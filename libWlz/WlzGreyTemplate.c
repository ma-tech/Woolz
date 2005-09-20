#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzGreyTemplate.c
* \author       richard Baldock
* \date         March 1999.
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
* \brief	Attach a grey table to a template object.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

static WlzObject *WlzGreyTemplate3d(WlzObject	*obj,
				WlzObject	*tmpl,
				WlzPixelV	tmplVal,
				WlzErrorNum	*dstErr);


/* function:     WlzGreyTemplate    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        
*
* \return       New object with the same domain <tt>tmpl</tt> but
 values in the intersection with <tt>obj</tt> set to those of the object.
 Returns NULL on error.
* \param    obj	Input object to which the template is applied
* \param    tmpl	Template object
* \param    tmplVal	Template value for regions in the template
 not in the original object
* \param    dstErr	Error return.
* \par      Source:
*                WlzGreyTemplate.c
*/
WlzObject *WlzGreyTemplate(
  WlzObject	*obj,
  WlzObject	*tmpl,
  WlzPixelV	tmplVal,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzObject	*obj1, *obj2;
  WlzValues	values;
  WlzPixelV	bckgrnd;
  WlzObjectType	type;
  WlzGreyType	gtype=WLZ_GREY_UBYTE;
  WlzIntervalWSpace	iwsp1, iwsp2;
  WlzGreyWSpace		gwsp1, gwsp2;
  int			size;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check obj */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      bckgrnd = WlzGetBackground(obj, &errNum);
      gtype = WlzGreyTableTypeToGreyType(obj->values.core->type, NULL);
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzGreyTemplate3d(obj, tmpl, tmplVal, dstErr);

    case WLZ_TRANS_OBJ:
      if( values.obj = WlzGreyTemplate(obj->values.obj, tmpl,
				       tmplVal, &errNum) ){
	return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      bckgrnd.type = WLZ_GREY_UBYTE;
      bckgrnd.v.ubv = 0;
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* check the template */
  if( errNum == WLZ_ERR_NONE ){
    if( tmpl == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      values.core = NULL;
      switch( tmpl->type ){
      case WLZ_2D_DOMAINOBJ:
	rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, tmpl->domain, values,
			      NULL, NULL, &errNum);
	break;

      case WLZ_TRANS_OBJ:
	rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, tmpl->values.obj->domain,
			      values, NULL, NULL, &errNum);
	break;

      case WLZ_EMPTY_OBJ:
	return WlzMakeEmpty(dstErr);

      case WLZ_2D_POLYGON:
	rtnObj = WlzPolyToObj(tmpl->domain.poly, WLZ_SIMPLE_FILL, &errNum);
	break;

      case WLZ_BOUNDLIST:
	rtnObj = WlzBoundToObj(tmpl->domain.b, WLZ_SIMPLE_FILL, &errNum);
	break;

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* attach a value table to the template and set to the template value,
     note the background is set to the input object or zero if empty */
  if( errNum == WLZ_ERR_NONE ){
    type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, gtype, NULL);
    if( values.v = WlzNewValueTb(rtnObj, type, bckgrnd, &errNum) ){
      rtnObj->values = WlzAssignValues(values, NULL);
      errNum = WlzGreySetValue(rtnObj, tmplVal);
    }
  }

  /* copy input obj values within the intersection */
  if( errNum == WLZ_ERR_NONE ){
    if((obj->type != WLZ_EMPTY_OBJ) ){
      if( (obj1 = WlzIntersect2(obj, rtnObj, &errNum)) ){
	obj1->values = WlzAssignValues(rtnObj->values, NULL);
	obj2 = WlzMakeMain(obj1->type, obj1->domain, obj->values,
			   NULL, NULL, NULL);

	errNum = WlzInitGreyScan(obj1, &iwsp1, &gwsp1);
	errNum = WlzInitGreyScan(obj2, &iwsp2, &gwsp2);
	switch( gwsp1.pixeltype ){
	case WLZ_GREY_INT:
	  size = sizeof(int);
	  break;
	case WLZ_GREY_SHORT:
	  size = sizeof(short);
	  break;
	case WLZ_GREY_UBYTE:
	  size = sizeof(UBYTE);
	  break;
	case WLZ_GREY_FLOAT:
	  size = sizeof(float);
	  break;
	case WLZ_GREY_DOUBLE:
	  size = sizeof(double);
	  break;
	case WLZ_GREY_RGBA:
	  size = sizeof(UINT);
	  break;
	}

	while( (errNum = WlzNextGreyInterval(&iwsp1)) == WLZ_ERR_NONE ){
	  (void) WlzNextGreyInterval(&iwsp2);
	  memcpy((void *) gwsp1.u_grintptr.inp,
		 (const void *) gwsp2.u_grintptr.inp,
		 size * iwsp1.colrmn);
	}
	if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	}
	WlzFreeObj(obj2);
	WlzFreeObj(obj1);
      }
      else {
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

static WlzObject *WlzGreyTemplate3d(
  WlzObject	*obj,
  WlzObject	*tmpl,
  WlzPixelV	tmplVal,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzObject	*tmpObj, *obj1, *obj2;
  WlzDomain	domain, *domains;
  WlzValues	values, *valuess;
  WlzPlaneDomain	*pdom;
  int		p;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object - it is non-NULL and 3D but the
     domain needs checking */
  if( obj->domain.p == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else {
    switch( obj->domain.p->type ){
    case WLZ_2D_DOMAINOBJ:
      /* check there is a valuetable */
      if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      break;

    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }

  /* check the template and create the return object */
  if( errNum == WLZ_ERR_NONE ){
    if( tmpl == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      values.core = NULL;
      switch( tmpl->type ){
      case WLZ_2D_DOMAINOBJ:
	pdom = obj->domain.p;
	if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					  pdom->plane1, pdom->lastpl,
					  pdom->line1, pdom->lastpl,
					  pdom->kol1, pdom->lastkl,
					  &errNum) ){
	  domain.p->voxel_size[0] = pdom->voxel_size[0];
	  domain.p->voxel_size[1] = pdom->voxel_size[1];
	  domain.p->voxel_size[2] = pdom->voxel_size[2];
	  for(p=pdom->plane1; p <= pdom->lastpl; p++){
	    domain.p->domains[p - pdom->plane1] = WlzAssignDomain(tmpl->domain,
								  NULL);
	  }
	  rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			       NULL, NULL, &errNum);
	}
	break;

      case WLZ_2D_POLYGON:
	pdom = obj->domain.p;
	if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					  pdom->plane1, pdom->lastpl,
					  pdom->line1, pdom->lastpl,
					  pdom->kol1, pdom->lastkl,
					  &errNum) ){
	  domain.p->voxel_size[0] = pdom->voxel_size[0];
	  domain.p->voxel_size[1] = pdom->voxel_size[1];
	  domain.p->voxel_size[2] = pdom->voxel_size[2];
	  obj1 = WlzPolyToObj(tmpl->domain.poly, WLZ_SIMPLE_FILL, &errNum);
	  for(p=pdom->plane1; p <= pdom->lastpl; p++){
	    domain.p->domains[p - pdom->plane1] = WlzAssignDomain(obj1->domain,
								  NULL);
	  }
	  WlzFreeObj(obj1);
	  rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			       NULL, NULL, &errNum);
	}
	break;

      case WLZ_BOUNDLIST:
	pdom = obj->domain.p;
	if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					  pdom->plane1, pdom->lastpl,
					  pdom->line1, pdom->lastpl,
					  pdom->kol1, pdom->lastkl,
					  &errNum) ){
	  domain.p->voxel_size[0] = pdom->voxel_size[0];
	  domain.p->voxel_size[1] = pdom->voxel_size[1];
	  domain.p->voxel_size[2] = pdom->voxel_size[2];
	  obj1 = WlzBoundToObj(tmpl->domain.b, WLZ_SIMPLE_FILL, &errNum);
	  for(p=pdom->plane1; p <= pdom->lastpl; p++){
	    domain.p->domains[p - pdom->plane1] = WlzAssignDomain(obj1->domain,
								  NULL);
	  }
	  WlzFreeObj(obj1);
	  rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			       NULL, NULL, &errNum);
	}
	break;

      case WLZ_3D_DOMAINOBJ:
	if( tmpl->domain.p ){
	  switch( tmpl->domain.p->type ){
	  case WLZ_2D_DOMAINOBJ:
	    domain.p = tmpl->domain.p;
	    break;

	  case WLZ_PLANEDOMAIN_POLYGON:
	  case WLZ_PLANEDOMAIN_CONV_HULL:
	    pdom = tmpl->domain.p;
	    if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					      pdom->plane1, pdom->lastpl,
					      pdom->line1, pdom->lastpl,
					      pdom->kol1, pdom->lastkl,
					      &errNum) ){
	      domain.p->voxel_size[0] = pdom->voxel_size[0];
	      domain.p->voxel_size[1] = pdom->voxel_size[1];
	      domain.p->voxel_size[2] = pdom->voxel_size[2];
	      for(p=pdom->plane1; p <= pdom->lastpl; p++){
		if( pdom->domains[p-pdom->plane1].core ){
		  obj1 = WlzPolyToObj(pdom->domains[p-pdom->plane1].poly,
				      WLZ_SIMPLE_FILL, &errNum);
		  domain.p->domains[p - pdom->plane1] =
		    WlzAssignDomain(obj1->domain, NULL);
		}
		WlzFreeObj(obj1);
	      }
	      values.core = NULL;
	      rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
				   NULL, NULL, &errNum);
	    }
	    break;

	  case WLZ_PLANEDOMAIN_BOUNDLIST:
	    pdom = tmpl->domain.p;
	    if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					      pdom->plane1, pdom->lastpl,
					      pdom->line1, pdom->lastpl,
					      pdom->kol1, pdom->lastkl,
					      &errNum) ){
	      domain.p->voxel_size[0] = pdom->voxel_size[0];
	      domain.p->voxel_size[1] = pdom->voxel_size[1];
	      domain.p->voxel_size[2] = pdom->voxel_size[2];
	      for(p=pdom->plane1; p <= pdom->lastpl; p++){
		if( pdom->domains[p-pdom->plane1].core ){
		  obj1 = WlzBoundToObj(pdom->domains[p-pdom->plane1].b,
				      WLZ_SIMPLE_FILL, &errNum);
		  domain.p->domains[p - pdom->plane1] =
		    WlzAssignDomain(obj1->domain, NULL);
		}
		WlzFreeObj(obj1);
	      }
	      values.core = NULL;
	      rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
				   NULL, NULL, &errNum);
	    }
	    break;

	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	  }
	  if( errNum == WLZ_ERR_NONE ){
	    rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
				 NULL, NULL, &errNum);
	  }
	}
	else {
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	break;

      case WLZ_EMPTY_OBJ:
	return WlzMakeEmpty(dstErr);

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* now we have a 3D obj and 3D template so run through the template
     and map values as required, note we must check that all the valuetables
     have the same type ie switch to obj type if necessary */
  if( errNum == WLZ_ERR_NONE ){
    WlzDomain	*objDoms;
    WlzValues	*objVals;
    WlzGreyType	gtype=WLZ_GREY_UBYTE;

    /* attach a voxel table with empty values list */
    values.vox = WlzMakeVoxelValueTb(obj->values.vox->type,
				     rtnObj->domain.p->plane1,
				     rtnObj->domain.p->lastpl,
				     obj->values.vox->bckgrnd,
				     NULL, NULL);
    rtnObj->values = WlzAssignValues(values, NULL);

    /* set some local variables */
    pdom = rtnObj->domain.p;
    domains = rtnObj->domain.p->domains;
    valuess = rtnObj->values.vox->values;
    objDoms = obj->domain.p->domains;
    objVals = obj->values.vox->values;

    /* calculate the new valuetables */
    for(p=pdom->plane1; p <= pdom->lastpl; p++, domains++, valuess++){
      if(((*domains).core)){
	if((p >= obj->domain.p->plane1) &&
	   (p <= obj->domain.p->lastpl) &&
	   (objDoms[p - obj->domain.p->plane1].core) ){
	  tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				objDoms[p - obj->domain.p->plane1],
				objVals[p - obj->domain.p->plane1],
				NULL, NULL, NULL);
	  gtype = WlzGreyTableTypeToGreyType(tmpObj->values.core->type, NULL);
	}
	else {
	  tmpObj = WlzMakeEmpty(NULL);
	}
	tmpObj = WlzAssignObject(tmpObj, NULL);
	values.core = NULL;
	obj1 = WlzAssignObject(
	  WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
		      NULL, NULL, NULL), NULL);
	if( obj2 = WlzGreyTemplate(tmpObj, obj1, tmplVal, &errNum) ){
	  *valuess = WlzAssignValues(obj2->values, NULL);
	  WlzFreeObj(obj2);
	}
	WlzFreeObj(obj1);
	WlzFreeObj(tmpObj);
      }
    }

    /* now check all valuetables have the same grey type */
    domains = rtnObj->domain.p->domains;
    valuess = rtnObj->values.vox->values;
    for(p=pdom->plane1; p <= pdom->lastpl; p++, domains++, valuess++){
      if((*domains).core &&
	 (WlzGreyTableTypeToGreyType((*valuess).core->type, NULL) != gtype) ){
	obj1 = WlzAssignObject(
	  WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *valuess,
		      NULL, NULL, NULL), NULL);
	if( obj2 = WlzConvertPix(obj1, gtype, &errNum) ){
	  /* substitute the valuetable in the voxel table array */
	  WlzFreeValues(*valuess);
	  *valuess = WlzAssignValues(obj2->values, NULL);
	  WlzFreeObj(obj2);
	}
	WlzFreeObj(obj1);
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}


