#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzGreyTransfer.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Jan 16 10:02:22 2002
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzValuesUtils
* \brief        Procedures to transfer grey values from a source to
 a destination object. The source object grey values are set in the
 intersection domain between source and destination. Destination
 domain and the destination values outside of the intersection are
 unchanged.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
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
  WlzPixelV	bckgrnd;
  WlzObjectType	type;
  WlzGreyType	gtype=WLZ_GREY_UBYTE;
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
      bckgrnd = WlzGetBackground(obj, &errNum);
      gtype = WlzGreyTableTypeToGreyType(obj->values.core->type, NULL);
      rtnObj = WlzCopyObject(obj, &errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzGreyTransfer3d(obj, srcObj, dstErr);

    case WLZ_TRANS_OBJ:
      if( values.obj = WlzGreyTransfer(obj->values.obj, srcObj,
				       &errNum) ){
	return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      bckgrnd.type = WLZ_GREY_UBYTE;
      bckgrnd.v.ubv = 0;
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
    WlzGreyType	gtype=WLZ_GREY_UBYTE;

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


