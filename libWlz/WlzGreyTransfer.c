#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyTransfer_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreyTransfer.c
* \author       Richard Baldock
* \date         September 2002
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
* \brief	Transfers grey values from a source to a destination
* 		object. The source object grey values are set in the
* 		intersection domain between source and destination.
* 		Destination domain and the destination values outside of the
* 		intersection are unchanged.
* \ingroup	WlzValuesUtils
*/

#include <stdio.h>
#include <string.h>

#include <Wlz.h>

static WlzObject *WlzGreyTransfer3d(WlzObject	*obj,
				    WlzObject	*tmpl,
				    WlzErrorNum	*dstErr);



/* function:     WlzGreyTransfer    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Transfer grey values from the source object to the
 destination object. Currently it is assumed that the objects are
 of the same type (2D/3D) and have the same grey-value type.
*
* \return       Woolz object with transferred grey values
* \param    obj	destination object
* \param    srcObj	source object
* \param    dstErr	error return
* \par      Source:
*                WlzGreyTransfer.c
*/
WlzObject *WlzGreyTransfer(
  WlzObject	*obj,
  WlzObject	*srcObj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzObject	*obj1, *obj2;
  WlzValues	values;
  WlzIntervalWSpace	iwsp1, iwsp2;
  WlzGreyWSpace		gwsp1, gwsp2;
  int			size;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check destination obj */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      rtnObj = WlzCopyObject(obj, &errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzGreyTransfer3d(obj, srcObj, dstErr);

    case WLZ_TRANS_OBJ:
      if((values.obj = WlzGreyTransfer(obj->values.obj, srcObj,
				       &errNum)) != NULL){
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

  /* check the source object */
  if( errNum == WLZ_ERR_NONE ){
    if( srcObj == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      switch( srcObj->type ){
      case WLZ_2D_DOMAINOBJ:
	break;

      case WLZ_TRANS_OBJ:
	srcObj = srcObj->values.obj;
	break;

      case WLZ_EMPTY_OBJ:
	if( dstErr ){
	  *dstErr = errNum;
	}
	return rtnObj;

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* copy source obj values within the intersection */
  if( errNum == WLZ_ERR_NONE ){
    if((srcObj->type != WLZ_EMPTY_OBJ) ){
      if( (obj1 = WlzIntersect2(srcObj, rtnObj, &errNum)) ){
	obj1->values = WlzAssignValues(rtnObj->values, NULL);
	obj2 = WlzMakeMain(obj1->type, obj1->domain, srcObj->values,
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
	  size = sizeof(WlzUByte);
	  break;
	case WLZ_GREY_FLOAT:
	  size = sizeof(float);
	  break;
	case WLZ_GREY_DOUBLE:
	  size = sizeof(double);
	  break;
	case WLZ_GREY_RGBA:
	  size = sizeof(WlzUInt);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
	}

	while((errNum == WLZ_ERR_NONE) &&
	      (errNum = WlzNextGreyInterval(&iwsp1)) == WLZ_ERR_NONE){
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

/* function:     WlzGreyTransfer3d    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        static function to implement WlzGreyTransfer for 3D objects.
*
* \return       woolz object
* \param    obj	destination object
* \param    srcObj	source object
* \param    dstErr	error return
* \par      Source:
*                WlzGreyTransfer.c
*/
static WlzObject *WlzGreyTransfer3d(
  WlzObject	*obj,
  WlzObject	*srcObj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzObject	*obj1, *obj2, *tmpObj;
  WlzDomain	*domains;
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

  /* check the source object */
  if( errNum == WLZ_ERR_NONE ){
    if( srcObj == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      switch( srcObj->type ){

      case WLZ_3D_DOMAINOBJ:
	if( srcObj->domain.p ){
	  switch( srcObj->domain.p->type ){
	  case WLZ_2D_DOMAINOBJ:
	    break;

	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	  }
	}
	else {
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	break;

      case WLZ_EMPTY_OBJ:
	return WlzCopyObject(obj, dstErr);

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* now we have a 3D obj and 3D srcObject so run through the source
     and map values as required */
  if( errNum == WLZ_ERR_NONE ){
    WlzDomain	*objDoms;
    WlzValues	*objVals;

    /* attach a voxel table with empty values list */
    values.vox = WlzMakeVoxelValueTb(obj->values.vox->type,
				     obj->domain.p->plane1,
				     obj->domain.p->lastpl,
				     obj->values.vox->bckgrnd,
				     NULL, NULL);
    rtnObj = WlzMakeMain(obj->type, obj->domain, values, NULL, NULL, &errNum);

    /* set some local variables */
    pdom = rtnObj->domain.p;
    domains = rtnObj->domain.p->domains;
    valuess = rtnObj->values.vox->values;
    objDoms = obj->domain.p->domains;
    objVals = obj->values.vox->values;

    /* calculate the new valuetables */
    for(p=pdom->plane1; p <= pdom->lastpl;
	p++, domains++, valuess++, objDoms++, objVals++){
      if(((*domains).core)){
	obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *objDoms, *objVals,
			   NULL, NULL, &errNum);
	obj1 = WlzAssignObject(obj1, &errNum);
			   
	if((p >= srcObj->domain.p->plane1) &&
	   (p <= srcObj->domain.p->lastpl) &&
	   (srcObj->domain.p->domains[p-srcObj->domain.p->plane1].core)){
	  obj2 = 
	    WlzMakeMain(WLZ_2D_DOMAINOBJ,
			srcObj->domain.p->domains[p-srcObj->domain.p->plane1],
			srcObj->values.vox->values[p-srcObj->domain.p->plane1],
			NULL, NULL, &errNum);
	}
	else {
	  obj2 = WlzMakeEmpty(NULL);
	}
	obj2 = WlzAssignObject(obj2, &errNum);

	tmpObj = WlzGreyTransfer(obj1, obj2, &errNum);
	*valuess = WlzAssignValues(tmpObj->values, &errNum);
	WlzFreeObj(obj1);
	WlzFreeObj(obj2);
	WlzFreeObj(tmpObj);
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}


