#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzConvertPix.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for converting between the various pixel
*		types.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************
* 6/7/1		richard	put in WlzConverVtx
************************************************************************/
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

static WlzObject *WlzConvertPix3d(WlzObject	*obj,
				  WlzGreyType	newpixtype,
				  WlzErrorNum	*dstErr);
static WlzObject *WlzConvertVtx3d(WlzObject	*obj,
				  WlzVertexType	newVtxType,
				  WlzErrorNum	*dstErr);

/************************************************************************
*   Function   : WlzConvertPix						*
*   Date       : Sat Nov  2 13:30:12 1996				*
*************************************************************************
*   Synopsis   :Convert the pixel type of an image object		*
*   Returns    :WlzObject *: New object with converted valuetable. The	*
*		domain is identical to the input. Returns NULL on error	*
*		The error can be obtained with WlzGetError or by setting*
*		an error handler. Possible errors:			*
*		WLZ_ERR_OBJECT_NULL, WLZ_ERR_OBJECT_TYPE, 		*
*		WLZ_ERR_DOMAIN_NULL, WLZ_ERR_VALUES_NULL,		*
*		WLZ_ERR_GREY_TYPE				*
*   Parameters :WlzObject	*obj: the object for conversion		*
*		WlzGreyType	newpixtype: the required grey-value type*
*   Global refs:None							*
************************************************************************/

WlzObject *WlzConvertPix(
  WlzObject	*obj,
  WlzGreyType	newpixtype,
  WlzErrorNum	*dstErr)
{
  WlzGreyType		oldpixtype;
  WlzGreyV		g;
  WlzGreyP		go, gn;
  WlzIntervalWSpace	oldiwsp, newiwsp;
  WlzGreyWSpace		oldgwsp, newgwsp;
  WlzObjectType		newvtbltype;
  WlzObject 		*newobj=NULL;
  WlzPixelV		bg;
  WlzValues		newvalues,
  			values;
  int 			k;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzConvertPix3d(obj, newpixtype, dstErr);

    case WLZ_TRANS_OBJ:
      newobj = WlzConvertPix(obj->values.obj,
			     newpixtype, &errNum);
      if( errNum == WLZ_ERR_NONE ){
	newvalues.obj = newobj;
	return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, newvalues,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeMain(WLZ_EMPTY_OBJ, obj->domain, obj->values,
			 NULL, NULL, dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    }
  }

  /* check domain and valuetable */
  if( (errNum == WLZ_ERR_NONE) && (obj->domain.core == NULL) ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  if( (errNum == WLZ_ERR_NONE) && (obj->values.core == NULL) ){
    errNum = WLZ_ERR_VALUES_NULL;
  }
  if( errNum == WLZ_ERR_NONE ){
    oldpixtype = WlzGreyTableTypeToGreyType(obj->values.core->type, NULL);
  }
	
  /*
   * Set type of new value table so as to preserve
   * rectangular/single interval/multiple interval
   * type.
   */
  if( errNum == WLZ_ERR_NONE ){
    newvtbltype = WlzGreyTableTypeToTableType(obj->values.core->type,
					      &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    newvtbltype = WlzGreyTableType(newvtbltype, newpixtype, &errNum);
  }

  /* get the background  - note background now carries its own type */
  if( errNum == WLZ_ERR_NONE ){
    bg = WlzGetBackground(obj, &errNum);
    WlzValueConvertPixel(&bg, bg, newpixtype);
  }

  /*
   * Make the new object with new value table type and value table
   * allocated (but blank).  Share original idom.
   */
  if( errNum == WLZ_ERR_NONE ){
    values.v = WlzNewValueTb(obj, newvtbltype, bg, &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    newobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj->domain, values,
			 obj->plist, obj->assoc, &errNum);
  }

  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(obj, &oldiwsp, &oldgwsp);
  }
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(newobj, &newiwsp, &newgwsp);
  }

  while( ((errNum = WlzNextGreyInterval(&oldiwsp)) == WLZ_ERR_NONE)
	 && ((errNum = WlzNextGreyInterval(&newiwsp)) == WLZ_ERR_NONE) ){
    go = oldgwsp.u_grintptr;	
    gn = newgwsp.u_grintptr;

    for (k=oldiwsp.lftpos; k<=oldiwsp.rgtpos; k++){

      /* First copy into a temporary grey value - use int for
       int, short & UBYTE, use double for float and double */		   
      switch (oldpixtype) {
      case WLZ_GREY_INT:
	g.inv = *go.inp++;
	switch (newpixtype) {
	case WLZ_GREY_INT:
	  *gn.inp++ = g.inv; break;
	case WLZ_GREY_SHORT:
	  *gn.shp++ = WLZ_CLAMP(g.inv, SHRT_MIN, SHRT_MAX); break;
	case WLZ_GREY_UBYTE:
	  *gn.ubp++ = (UBYTE) (WLZ_CLAMP(g.inv, 0, UCHAR_MAX)); break;
	case WLZ_GREY_FLOAT:
	  *gn.flp++ = g.inv; break;
	case WLZ_GREY_DOUBLE:
	  *gn.dbp++ = g.inv; break;
	}
	break;
      case WLZ_GREY_SHORT:
	g.shv = *go.shp++;
	switch (newpixtype) {
	case WLZ_GREY_INT:
	  *gn.inp++ = g.shv; break;
	case WLZ_GREY_SHORT:
	  *gn.shp++ = g.shv; break;
	case WLZ_GREY_UBYTE:
	  *gn.ubp++ = (UBYTE) (WLZ_CLAMP(g.shv, 0, UCHAR_MAX)); break;
	case WLZ_GREY_FLOAT:
	  *gn.flp++ = g.shv; break;
	case WLZ_GREY_DOUBLE:
	  *gn.dbp++ = g.shv; break;
	}
	break;
      case WLZ_GREY_UBYTE:
	g.ubv = *go.ubp++;
	switch (newpixtype) {
	case WLZ_GREY_INT:
	  *gn.inp++ = g.ubv; break;
	case WLZ_GREY_SHORT:
	  *gn.shp++ = g.ubv; break;
	case WLZ_GREY_UBYTE:
	  *gn.ubp++ = g.ubv; break;
	case WLZ_GREY_FLOAT:
	  *gn.flp++ = g.ubv; break;
	case WLZ_GREY_DOUBLE:
	  *gn.dbp++ = g.ubv; break;
	}
	break;
      case WLZ_GREY_FLOAT:
	g.flv = *go.flp++;
	switch (newpixtype) {
	case WLZ_GREY_INT:
	  *gn.inp++ = WLZ_CLAMP(g.flv, INT_MIN, INT_MAX); break;
	case WLZ_GREY_SHORT:
	  *gn.shp++ = WLZ_CLAMP(g.flv, SHRT_MIN, SHRT_MAX); break;
	case WLZ_GREY_UBYTE:
	  *gn.ubp++ = (UBYTE) (WLZ_CLAMP(g.flv, 0, UCHAR_MAX)); break;
	case WLZ_GREY_FLOAT:
	  *gn.flp++ = g.flv; break;
	case WLZ_GREY_DOUBLE:
	  *gn.dbp++ = g.flv; break;
	}
	break;
      case WLZ_GREY_DOUBLE:
	g.dbv = *go.dbp++;
	switch (newpixtype) {
	case WLZ_GREY_INT:
	  *gn.inp++ = WLZ_CLAMP(g.dbv, INT_MIN, INT_MAX); break;
	case WLZ_GREY_SHORT:
	  *gn.shp++ = WLZ_CLAMP(g.dbv, SHRT_MIN, SHRT_MAX); break;
	case WLZ_GREY_UBYTE:
	  *gn.ubp++ = (UBYTE) (WLZ_CLAMP(g.dbv, 0, UCHAR_MAX)); break;
	case WLZ_GREY_FLOAT:
	  *gn.flp++ = WLZ_CLAMP(g.dbv, FLT_MIN, FLT_MAX); break;
	case WLZ_GREY_DOUBLE:
	  *gn.dbp++ = g.dbv; break;
	}
	break;
      }
    } /* for */
  } /* while */
  if(errNum == WLZ_ERR_EOO)	        /* Reset error from end of intervals */ 
  {
    errNum = WLZ_ERR_NONE;
  }
  if( dstErr ){
    *dstErr = errNum;
  }
  return newobj;
}	
		
/************************************************************************
*   Function   : WlzConvertPix3d					*
*   Date       : Tue Nov 12 14:23:07 1996				*
*************************************************************************
*   Synopsis   :Convert the pixel type of a 3D image object. This is	*
*		static to this module and should only be accessed via	*
*		WlzConvertPix.						*
*   Returns    :WlzObject *: the converted object or NULL in error	*
*		Possible errors: WLZ_ERR_DOMAIN_NULL, WLZ_ERR_VALUES_NULL,		*
*		WLZ_ERR_DOMAIN_TYPE, WLZ_ERR_VALUES_TYPE.		*
*   Parameters :WlzObject	*obj: object for conversion		*
*		WlzGreyType	newpixtype: the required pixel type	*
*   Global refs:None							*
************************************************************************/

static WlzObject *WlzConvertPix3d(
  WlzObject	*obj,
  WlzGreyType	newpixtype,
  WlzErrorNum	*dstErr)
{
  WlzObject		*tmp_obj1, *tmp_obj2;
  WlzPlaneDomain 	*pdom, *new_pdom;
  WlzVoxelValues	*voxtab, *new_voxtab;
  WlzDomain		domain;
  WlzValues		values;
  WlzPixelV		bg;
  int	       		p;
  WlzErrorNum		errNum=WLZ_ERR_NONE;


  /* the object, object type have been checked by WlzConvertPix but
     must still check domain, domain type and the voxel table */
  if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
    
  if( (errNum == WLZ_ERR_NONE) && (obj->values.core == NULL) ){
    errNum = WLZ_ERR_VALUES_NULL;
  }
    
  if( (errNum == WLZ_ERR_NONE) && 
     (obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN) ){
    errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
  }
    
  if( (errNum == WLZ_ERR_NONE) && 
     (obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY) ){
    errNum = WLZ_ERR_VOXELVALUES_TYPE;
  }

  /* set up new 3D domains */
  if( errNum == WLZ_ERR_NONE ){
    pdom = obj->domain.p;
    new_pdom = WlzMakePlaneDomain(pdom->type,
				  pdom->plane1, pdom->lastpl,
				  pdom->line1, pdom->lastln,
				  pdom->kol1, pdom->lastkl, &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    for(p=0; p < 3; p++){
      new_pdom->voxel_size[p] = pdom->voxel_size[p];
    }

    voxtab = obj->values.vox;
    /* get the background  - note background now carries its own type */
    bg = voxtab->bckgrnd;
    WlzValueConvertPixel(&bg, bg, newpixtype);

    new_voxtab = WlzMakeVoxelValueTb(voxtab->type,
				     voxtab->plane1, voxtab->lastpl,
				     bg, NULL, &errNum);
  }

  /* set up the interval and value tables */
  if( errNum == WLZ_ERR_NONE ){
    for(p=0; (p < (pdom->lastpl - pdom->plane1 + 1)) &&
	  (errNum == WLZ_ERR_NONE); p++){

      if( (pdom->domains[p]).core == NULL ){
	(new_pdom->domains[p]).i = NULL;
	(new_voxtab->values[p]).v = NULL;
	continue;
      }

      if( (voxtab->values[p]).core == NULL ){
	new_pdom->domains[p] = WlzAssignDomain(pdom->domains[p], NULL);
	(new_voxtab->values[p]).v = NULL;
	continue;
      }
	    
      tmp_obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, pdom->domains[p],
			     voxtab->values[p], NULL, NULL, &errNum);
      if( errNum != WLZ_ERR_NONE ){break;}

      tmp_obj2 = WlzConvertPix(tmp_obj1, newpixtype, &errNum);
      if( errNum != WLZ_ERR_NONE ){
	WlzFreeObj(tmp_obj1);
	break;
      }

      new_pdom->domains[p] = WlzAssignDomain(tmp_obj2->domain, NULL);
      new_voxtab->values[p] = WlzAssignValues(tmp_obj2->values, NULL);
      WlzFreeObj(tmp_obj2);
      WlzFreeObj(tmp_obj1);
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    domain.p = new_pdom;
    values.vox = new_voxtab;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values, NULL, NULL, dstErr);
}


WlzPolygonDomain *WlzConvertPolyType(
  WlzPolygonDomain	*pdom,
  WlzObjectType		type,
  WlzErrorNum		*dstErr)
{
  WlzPolygonDomain	*rtnDom=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzIVertex2		*iVtxs;
  WlzFVertex2		*fVtxs;
  WlzDVertex2		*dVtxs;

  if( pdom == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( pdom->nvertices <= 0 ){
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else {
    if( rtnDom = WlzMakePolyDmn(type, NULL, pdom->nvertices,
				pdom->nvertices, 1, &errNum) ){
      switch( type ){
      case WLZ_POLYGON_INT:
	iVtxs = rtnDom->vtx;
	switch( pdom->type ){
	case WLZ_POLYGON_INT:
	  WlzValueCopyIVertexToIVertex(iVtxs, pdom->vtx, pdom->nvertices);
	  break;
	case WLZ_POLYGON_FLOAT:
	  WlzValueCopyFVertexToIVertex(iVtxs, (WlzFVertex2 *) pdom->vtx,
				       pdom->nvertices);
	  break;
	case WLZ_POLYGON_DOUBLE:
	  WlzValueCopyDVertexToIVertex(iVtxs, (WlzDVertex2 *) pdom->vtx,
				       pdom->nvertices);
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  WlzFreePolyDmn(rtnDom);
	  rtnDom = NULL;
	  break;
	}
	break;

      case WLZ_POLYGON_FLOAT:
	fVtxs = (WlzFVertex2 *) rtnDom->vtx;
	switch( pdom->type ){
	case WLZ_POLYGON_INT:
	  WlzValueCopyIVertexToFVertex(fVtxs, pdom->vtx, pdom->nvertices);
	  break;
	case WLZ_POLYGON_FLOAT:
	  WlzValueCopyFVertexToFVertex(fVtxs, (WlzFVertex2 *) pdom->vtx,
				       pdom->nvertices);
	  break;
	case WLZ_POLYGON_DOUBLE:
	  WlzValueCopyDVertexToFVertex(fVtxs, (WlzDVertex2 *) pdom->vtx,
				       pdom->nvertices);
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  WlzFreePolyDmn(rtnDom);
	  rtnDom = NULL;
	  break;
	}
	break;

      case WLZ_POLYGON_DOUBLE:
	dVtxs = (WlzDVertex2 *) rtnDom->vtx;
	switch( pdom->type ){
	case WLZ_POLYGON_INT:
	  WlzValueCopyIVertexToDVertex(dVtxs, pdom->vtx, pdom->nvertices);
	  break;
	case WLZ_POLYGON_FLOAT:
	  WlzValueCopyFVertexToDVertex(dVtxs, (WlzFVertex2 *) pdom->vtx,
				       pdom->nvertices);
	  break;
	case WLZ_POLYGON_DOUBLE:
	  WlzValueCopyDVertexToDVertex(dVtxs, (WlzDVertex2 *) pdom->vtx,
				       pdom->nvertices);
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  WlzFreePolyDmn(rtnDom);
	  rtnDom = NULL;
	  break;
	}
	break;

      default:
	errNum = WLZ_ERR_PARAM_TYPE;
	break;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnDom;
}

WlzBoundList *WlzConvertBoundType(
  WlzBoundList		*bound,
  WlzObjectType		type,
  WlzErrorNum		*dstErr)
{
  WlzBoundList	*rtnBound=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( bound == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( rtnBound = WlzMakeBoundList(bound->type, bound->wrap,
				       NULL, &errNum) ){
    if( (errNum == WLZ_ERR_NONE) && bound->next ){
      rtnBound->next =
	WlzAssignBoundList(WlzConvertBoundType(bound->next, type,
					       &errNum), NULL);
    }

    if( (errNum == WLZ_ERR_NONE) && bound->down ){
      rtnBound->down =
	WlzAssignBoundList(WlzConvertBoundType(bound->down, type,
					       &errNum), NULL);
    }

    if( (errNum == WLZ_ERR_NONE) && bound->poly ){
      rtnBound->poly =
	WlzAssignPolygonDomain(WlzConvertPolyType(bound->poly,
						  type, &errNum), NULL);
    }
    if( errNum != WLZ_ERR_NONE ){
      WlzFreeBoundList(rtnBound);
      rtnBound = NULL;
    }
      
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnBound;
}

WlzObject *WlzConvertVtx(
  WlzObject	*obj,
  WlzVertexType	newVtxType,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzDomain	domain;
  WlzValues	values;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE){
    switch( obj->type ){

    case WLZ_2D_POLYGON:
      switch( newVtxType ){
      case WLZ_VERTEX_I2:
	domain.poly = WlzConvertPolyType(obj->domain.poly, WLZ_POLYGON_INT,
					 &errNum);
	break;

      case WLZ_VERTEX_F2:
	domain.poly = WlzConvertPolyType(obj->domain.poly, WLZ_POLYGON_FLOAT,
					 &errNum);
	break;

      case WLZ_VERTEX_D2:
	domain.poly = WlzConvertPolyType(obj->domain.poly, WLZ_POLYGON_DOUBLE,
					 &errNum);
	break;

      default:
	domain.poly = NULL;
	errNum = WLZ_ERR_PARAM_TYPE;
      }
      if( errNum == WLZ_ERR_NONE ){
	values.core = NULL;
	rtnObj = WlzMakeMain(obj->type, domain, values, NULL, NULL, NULL);
      }
      break;

    case WLZ_BOUNDLIST:
      switch( newVtxType ){
      case WLZ_VERTEX_I2:
	domain.b = WlzConvertBoundType(obj->domain.b, WLZ_POLYGON_INT,
				       &errNum);
	break;

      case WLZ_VERTEX_F2:
	domain.b = WlzConvertBoundType(obj->domain.b, WLZ_POLYGON_FLOAT,
				       &errNum);
	break;

      case WLZ_VERTEX_D2:
	domain.b = WlzConvertBoundType(obj->domain.b, WLZ_POLYGON_DOUBLE,
				       &errNum);
	break;

      default:
	domain.b = NULL;
	errNum = WLZ_ERR_PARAM_TYPE;
      }
      if( errNum == WLZ_ERR_NONE ){
	values.core = NULL;
	rtnObj = WlzMakeMain(obj->type, domain, values, NULL, NULL, NULL);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      rtnObj = WlzConvertVtx3d(obj, newVtxType, &errNum);
      break;

    case WLZ_TRANS_OBJ:
      values.obj = WlzConvertVtx(obj->values.obj,
				 newVtxType, &errNum);
      if( errNum == WLZ_ERR_NONE ){
	rtnObj = WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, values,
			     NULL, NULL, &errNum);
      }
      break;

    case WLZ_EMPTY_OBJ:
      rtnObj = WlzMakeMain(WLZ_EMPTY_OBJ, obj->domain, obj->values,
			   NULL, NULL, &errNum);

    case WLZ_3D_POLYGON:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;

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

static WlzObject *WlzConvertVtx3d(
  WlzObject	*obj,
  WlzVertexType	newVtxType,
  WlzErrorNum	*dstErr)
{
  WlzObject		*rtnObj=NULL;
  WlzObject		*tmp_obj1, *tmp_obj2;
  WlzPlaneDomain 	*pdom, *new_pdom;
  WlzDomain		domain;
  WlzValues		values;
  int	       		p;
  WlzErrorNum		errNum=WLZ_ERR_NONE;


  /* the object, object type have been checked by WlzConvertVtx but
     must still check domain, domain type */
  if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else {
    switch( obj->domain.p->type ){
    case WLZ_PLANEDOMAIN_POLYGON:
    case WLZ_PLANEDOMAIN_BOUNDLIST:
      break;

    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }
    
  /* set up new 3D domains */
  if( errNum == WLZ_ERR_NONE ){
    pdom = obj->domain.p;
    new_pdom = WlzMakePlaneDomain(pdom->type,
				  pdom->plane1, pdom->lastpl,
				  pdom->line1, pdom->lastln,
				  pdom->kol1, pdom->lastkl, &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    for(p=0; p < 3; p++){
      new_pdom->voxel_size[p] = pdom->voxel_size[p];
    }
  }

  /* set up the interval and value tables */
  if( errNum == WLZ_ERR_NONE ){
    values.core = NULL;
    for(p=0; (p < (pdom->lastpl - pdom->plane1 + 1)) &&
	  (errNum == WLZ_ERR_NONE); p++){

      if( (pdom->domains[p]).core == NULL ){
	(new_pdom->domains[p]).i = NULL;
	continue;
      }

      if( pdom->type == WLZ_PLANEDOMAIN_POLYGON ){
	tmp_obj1 = WlzMakeMain(WLZ_2D_POLYGON, pdom->domains[p],
			       values, NULL, NULL, &errNum);
      }
      else {
	tmp_obj1 = WlzMakeMain(WLZ_BOUNDLIST, pdom->domains[p],
			       values, NULL, NULL, &errNum);
      }

      if( errNum != WLZ_ERR_NONE ){break;}

      tmp_obj2 = WlzConvertVtx(tmp_obj1, newVtxType, &errNum);
      if( errNum != WLZ_ERR_NONE ){
	WlzFreeObj(tmp_obj1);
	break;
      }

      new_pdom->domains[p] = WlzAssignDomain(tmp_obj2->domain, NULL);
      WlzFreeObj(tmp_obj2);
      WlzFreeObj(tmp_obj1);
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    domain.p = new_pdom;
    values.core = NULL;
    rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values, NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

