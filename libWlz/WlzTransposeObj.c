#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTransposeObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzTransposeObj.c
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
* \brief	Transposes objects (ie interchanges row and column
* 		coordinates).
* \ingroup	WlzTransform
*/

#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

static WlzObject *WlzTransposeRectObj(WlzObject *obj,
				      WlzErrorNum *dstErr);
static WlzBoundList *WlzTransposeBound(WlzBoundList	*blist,
				       WlzErrorNum *dstErr);
static WlzPolygonDomain *WlzTransposePolygon(WlzPolygonDomain *poly,
					     WlzErrorNum *dstErr);
static WlzObject *WlzTranspose3DObj(WlzObject *obj,
				      WlzErrorNum *dstErr);

/*! 
* \return       Transposed object, NULL on error.
* \ingroup      WlzTransform
* \brief        Transpose a woolz object. Currently implemented for
*		2D domain, boundary and polyline objects. Transpose
*		means pixels/vertices at (k,l) are moved to (l,k).
*		Domain, boundary and polyline objects can all be transposed.
*		This is used in the WlzSepTrans procedure.
* \par      Source:
*                WlzTransposeObj.c
* \param    	obj			Input object
* \param    	dstErr			Error return.
*/
WlzObject *WlzTransposeObj(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject		*bobj = NULL, *nobj = NULL;
  WlzGreyValueWSpace	*gVWSp = NULL;
  WlzDomain		domain;
  WlzValues		nvalues,
  			values;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  int			i;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  domain.core = NULL;
  values.core = NULL;
  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    /* switch on object type, transpose the domain now, values later */
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }

      switch( obj->domain.core->type ){

      case WLZ_INTERVALDOMAIN_INTVL:
	if((bobj = WlzObjToBoundary(obj, 1, &errNum)) != NULL){
	  if((domain.b = WlzTransposeBound(bobj->domain.b, &errNum)) != NULL){
	    nobj = WlzBoundToObj(domain.b, WLZ_SIMPLE_FILL, &errNum);
	    WlzFreeObj(bobj);
	    WlzFreeBoundList(domain.b);
	  }
	  else {
	    WlzFreeObj(bobj);
	  }
	}
	break;

      case WLZ_INTERVALDOMAIN_RECT:
	return WlzTransposeRectObj(obj, dstErr);

      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzTranspose3DObj(obj, dstErr);

    case WLZ_TRANS_OBJ:
      if((nobj = WlzTransposeObj(obj->values.obj, &errNum)) != NULL){
	nvalues.obj = nobj;
	return WlzMakeMain(obj->type, obj->domain, nvalues,
			   NULL, obj, dstErr);
      }
      break;

    case WLZ_2D_POLYGON:
      if((domain.poly = WlzTransposePolygon(obj->domain.poly,
                                            &errNum)) != NULL){
	values.core = NULL;
	return WlzMakeMain(WLZ_2D_POLYGON, domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_BOUNDLIST:
      if((domain.b = WlzTransposeBound(obj->domain.b, &errNum)) != NULL){
	values.core = NULL;
	return WlzMakeMain(WLZ_BOUNDLIST, domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_RECTANGLE:
      if( obj->domain.r == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      switch( obj->domain.r->type ){

      case WLZ_RECTANGLE_DOMAIN_INT:
	domain.r = (WlzIRect *) AlcMalloc(sizeof(WlzIRect));
	domain.r->type = WLZ_RECTANGLE_DOMAIN_INT;
	domain.r->linkcount = 0;
	domain.r->freeptr = NULL;
	for(i=0; i < 4; i++){
	  domain.r->irk[i] = obj->domain.r->irl[i];
	  domain.r->irl[i] = obj->domain.r->irk[i];
	}
	domain.r->rangle = (float )(WLZ_M_PI_2 - obj->domain.r->rangle);
	return WlzMakeMain(WLZ_RECTANGLE, domain, values,
			   NULL, NULL, dstErr);

      case WLZ_RECTANGLE_DOMAIN_FLOAT:
	domain.fr = (WlzFRect *) AlcMalloc(sizeof(WlzFRect));
	domain.fr->type = WLZ_RECTANGLE_DOMAIN_FLOAT;
	domain.fr->linkcount = 0;
	domain.fr->freeptr = NULL;
	for(i=0; i < 4; i++){
	  domain.fr->frk[i] = obj->domain.fr->frl[i];
	  domain.fr->frl[i] = obj->domain.fr->frk[i];
	}
	domain.fr->rangle = (float )(WLZ_M_PI_2 - obj->domain.fr->rangle);
	return WlzMakeMain(WLZ_RECTANGLE, domain, values,
			   NULL, NULL, dstErr);

      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* now attach a grey-table */
  if((errNum == WLZ_ERR_NONE) && obj->values.core ){
    WlzObjectType tt;
    WlzGreyType	gt;

    gt = WlzGreyTableTypeToGreyType(obj->values.core->type, NULL);
    tt = WlzGreyTableTypeToTableType(obj->values.core->type, NULL);
    if(tt == WLZ_GREY_TAB_TILED)
    {
      tt = WLZ_GREY_TAB_RAGR;
    }
    if((values.v = WlzNewValueTb(nobj, WlzGreyTableType(tt, gt, NULL),
				 WlzGetBackground(obj,
				                  NULL), &errNum)) != NULL){
      nobj->values = WlzAssignValues(values, NULL);

      /* transfer values */
      errNum = WlzInitGreyScan(nobj, &iwsp, &gwsp);
      if(errNum == WLZ_ERR_NONE){
	gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
	if(errNum == WLZ_ERR_NONE){
	  while( WlzNextGreyInterval(&iwsp) == WLZ_ERR_NONE ){
	    for(i=iwsp.lftpos; i<=iwsp.rgtpos; i++){
	      WlzGreyValueGet(gVWSp, 0, i, iwsp.linpos);
	      switch(gVWSp->gType) {
	      case WLZ_GREY_INT:
		*gwsp.u_grintptr.inp++ = (*(gVWSp->gVal)).inv;
		break;
	      case WLZ_GREY_SHORT:
		*gwsp.u_grintptr.shp++ = (*(gVWSp->gVal)).shv;
		break;
	      case WLZ_GREY_UBYTE:
		*gwsp.u_grintptr.ubp++ = (*(gVWSp->gVal)).ubv;
		break;
	      case WLZ_GREY_FLOAT:
		*gwsp.u_grintptr.flp++ = (*(gVWSp->gVal)).flv;
		break;
	      case WLZ_GREY_DOUBLE:
		*gwsp.u_grintptr.dbp++ = (*(gVWSp->gVal)).dbv;
		break;
	      case WLZ_GREY_RGBA:
		*gwsp.u_grintptr.rgbp++ = (*(gVWSp->gVal)).rgbv;
		break;
	      default:
		break;
	      }
	    }
	  }
	}
	(void )WlzEndGreyScan(&iwsp, &gwsp);
      }
      if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
      {
	errNum = WLZ_ERR_NONE;
      }
      WlzGreyValueFreeWSp(gVWSp);
    }
  }
  if(errNum != WLZ_ERR_NONE){
    (void )WlzFreeObj(nobj);
    nobj = NULL;
  }
  /* return new object */
  if( dstErr ){
    *dstErr = errNum;
  }
  return(nobj);
}

/* function:     WlzTransposeBound    */
/*! 
* \ingroup      WlzTransform
* \brief        Transpose a boundlist structure, only accessed via
 WlzTransposeObj().
*
* \return       Transposed boundary list
* \param    blist	Input boundary list domain
* \param    dstErr	Error return.
* \par      Source:
*                WlzTransposeObj.c
*/
static WlzBoundList *WlzTransposeBound(
  WlzBoundList	*blist,
  WlzErrorNum	*dstErr)
{
  WlzBoundList		*newblist=NULL;
  WlzPolygonDomain	*poly;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* minimal checking since this is only called from WlzTransposeObj */
  if( blist == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if( blist->poly == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else {
    /* make a new polygon with transposed vertices */
    if((poly = WlzTransposePolygon(blist->poly, &errNum)) != NULL){
      /* make new boundlist structure */
      if((newblist = WlzMakeBoundList(blist->type,
				      blist->wrap, poly, &errNum)) == NULL){
	WlzFreePolyDmn( poly );
      }
    }

    /* now rotate next and down */
    if( (errNum == WLZ_ERR_NONE) && blist->next ){
      newblist->next = WlzAssignBoundList(
	WlzTransposeBound(blist->next, &errNum), NULL);
    }
    if( errNum != WLZ_ERR_NONE ){
      WlzFreeBoundList(newblist);
      newblist = NULL;
    }

    if( (errNum == WLZ_ERR_NONE) && blist->down ){
      newblist->down = WlzAssignBoundList(
	WlzTransposeBound(blist->down, &errNum), NULL);
    }
    if( errNum != WLZ_ERR_NONE ){
      WlzFreeBoundList(newblist);
      newblist = NULL;
    }
  }

  /* return newlist */
  if( dstErr ){
    *dstErr = errNum;
  }
  return newblist;
}

/* function:     WlzTransposePolygon    */
/*! 
* \ingroup      WlzTransform
* \brief        Transpose a polyline, only access from WlzTransposeObj().
*
* \return       Transposed polygon domain.
* \param    poly	Input polygon domain
* \param    dstErr	Error return.
* \par      Source:
*                WlzTransposeObj.c
*/
static WlzPolygonDomain *WlzTransposePolygon(
  WlzPolygonDomain *poly,
  WlzErrorNum	*dstErr)
{
  WlzPolygonDomain	*newpoly=NULL;
  int 			i;
  int			inx;
  double		dbx;
  WlzIVertex2 		*ivtx;
  WlzFVertex2		*fvtx;
  WlzDVertex2		*dvtx;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the domain */
  if( poly == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  /* first make a new polygon domain, allocating space */
  else {
    newpoly = WlzMakePolygonDomain(poly->type, poly->nvertices, poly->vtx,
				poly->nvertices, 1, &errNum);
  }

  /* now transpose it */
  if( errNum == WLZ_ERR_NONE ){
    switch( poly->type ){

    case WLZ_POLYGON_INT:
      ivtx = newpoly->vtx;
      for(i=0; i < newpoly->nvertices; i++, ivtx++){
	inx = ivtx->vtX;
	ivtx->vtX = ivtx->vtY;
	ivtx->vtY = inx;
      }
      break;

    case WLZ_POLYGON_FLOAT:
      fvtx = (WlzFVertex2 *) newpoly->vtx;
      for(i=0; i < newpoly->nvertices; i++, fvtx++){
	dbx = fvtx->vtX;
	fvtx->vtX = fvtx->vtY;
	fvtx->vtY = (float )dbx;
      }
      break;

    case WLZ_POLYGON_DOUBLE:
      dvtx = (WlzDVertex2 *) newpoly->vtx;
      for(i=0; i < newpoly->nvertices; i++, dvtx++){
	dbx = dvtx->vtX;
	dvtx->vtX = dvtx->vtY;
	dvtx->vtY = dbx;
      }
      break;
    default:
      errNum = WLZ_ERR_POLYGON_TYPE;
      break;
    }
  }

  /* return newpoly */
  if( dstErr ){
    *dstErr = errNum;
  }
  return(newpoly);
}


/* function:     WlzTransposeRectObj    */
/*! 
* \ingroup      WlzTransform
* \brief        Transpose a rectangular object. Static function anly accessed
 via WlzTransposeObj().
*
* \return       Transposed object.
* \param    obj	Inout rectangular object.
* \param    dstErr	Error return.
* \par      Source:
*                WlzTransposeObj.c
*/
static WlzObject *WlzTransposeRectObj(
  WlzObject *obj,
  WlzErrorNum *dstErr)
{
  WlzObject	*nobj=NULL;
  WlzDomain	domain;
  WlzValues	values;
  int 		i, width, height, size = 0;
  WlzGreyP 	newvals;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* minimal checking since this is only called from WlzTransposeObj,
     make the new domain */
  domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				   obj->domain.i->kol1,
				   obj->domain.i->lastkl,
				   obj->domain.i->line1,
				   obj->domain.i->lastln,
				   &errNum);

  /* now the valuetable - convert to rectangular - note rectangular domain
     does not imply rectangular valuetable */
  if( errNum == WLZ_ERR_NONE ){
    values.core = NULL;
    if( obj->values.core ){
      WlzIntervalWSpace	iwsp;
      WlzGreyWSpace	gwsp;
      WlzObjectType	newtype = WLZ_NULL;

      width = obj->domain.i->lastln - obj->domain.i->line1 + 1;
      height = obj->domain.i->lastkl - obj->domain.i->kol1 + 1;

      switch( WlzGreyTableTypeToGreyType(obj->values.core->type, NULL) ){
      case WLZ_GREY_INT:
	size = sizeof(int);
	newtype = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_INT, NULL);
	break;
      case WLZ_GREY_SHORT:
	size = sizeof(short);
	newtype = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_SHORT, NULL);
	break;
      case WLZ_GREY_UBYTE:
	size = sizeof(WlzUByte);
	newtype = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_UBYTE, NULL);
	break;
      case WLZ_GREY_FLOAT:
	size = sizeof(float);
	newtype = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_FLOAT, NULL);
	break;
      case WLZ_GREY_DOUBLE:
	size = sizeof(double);
	newtype = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_DOUBLE, NULL);
	break;
      case WLZ_GREY_RGBA:
	size = sizeof(WlzUInt);
	newtype = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_RGBA, NULL);
	break;
      default:
        break;
      }

      /* allocate space */
      if( (newvals.inp = (int *) AlcMalloc(size*width*height)) == NULL ){
	WlzFreeIntervalDomain(domain.i);
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else {
	/* fill in values */
	width *= size;
	errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
	if(errNum == WLZ_ERR_NONE){
	  while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	    int	off1, off2;
	    off1 = (iwsp.linpos - iwsp.linbot) * size;
	    off2 = 0;
	    for(i=iwsp.lftpos; i<=iwsp.rgtpos; i++, off1 += width,
		  off2 += size){
	      (void) memcpy((void *) (newvals.ubp+off1),
			    (const void *) (gwsp.u_grintptr.ubp+off2), size);
	    }
	  }
	  (void )WlzEndGreyScan(&iwsp, &gwsp);
	}
	if(errNum == WLZ_ERR_EOO)	/* Reset error from end of intervals */
	{
	  errNum = WLZ_ERR_NONE;
	}

	/* create the valuetable */
	if( (values.r = WlzMakeRectValueTb(newtype,
					   obj->domain.i->kol1,
					   obj->domain.i->lastkl,
					   obj->domain.i->line1,
					   width / size,
					   WlzGetBackground(obj, NULL),
					   newvals.inp,
					   &errNum)) == NULL ){
	  if(errNum == WLZ_ERR_NONE) {
	    values.r->freeptr = AlcFreeStackPush(
	    	values.r->freeptr, (void *)(newvals.inp), NULL);
	  }
	  WlzFreeIntervalDomain(domain.i);
	  return NULL;
	}
      }
    }
  }

  /* make the new object */
  if( errNum == WLZ_ERR_NONE ){
    nobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
		       NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return nobj;
}

static WlzObject *WlzTranspose3DObj(
  WlzObject *obj,
  WlzErrorNum *dstErr)
{
  WlzObject	*rtnObj = NULL;
  WlzObject	*obj1 = NULL, *obj2 = NULL;
  WlzDomain	domain;
  WlzValues	values;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		p, indx;

  domain.core = NULL;
  values.core = NULL;
  /* only check for the plane domain - other chacks are done */
  if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core && WlzGreyTableIsTiled(obj->values.core->type)) {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else {
    if((domain.p = WlzMakePlaneDomain(obj->domain.p->type,
				      obj->domain.p->plane1,
				      obj->domain.p->lastpl,
				      obj->domain.p->kol1,
				      obj->domain.p->lastkl,
				      obj->domain.p->line1,
				      obj->domain.p->lastln,
				      &errNum)) != NULL){
      domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
      domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
      domain.p->voxel_size[2] = obj->domain.p->voxel_size[2];
      if( obj->values.core ){
	values.vox = WlzMakeVoxelValueTb(obj->values.vox->type,
					 obj->values.vox->plane1,
					 obj->values.vox->lastpl,
					 obj->values.vox->bckgrnd,
					 NULL, &errNum);
      }
      else {
	values.vox = NULL;
      }

      if( errNum == WLZ_ERR_NONE ){
	rtnObj = WlzMakeMain(obj->type, domain, values,
			     NULL, NULL, &errNum);
      }
    
      if( errNum == WLZ_ERR_NONE ){
	for(p=obj->domain.p->plane1, indx=0;
	    p <= obj->domain.p->lastpl; p++, indx++){
	  if( obj->domain.p->domains[indx].core ){
	    domain = obj->domain.p->domains[indx];
	    if( obj->values.core ){
	      values = obj->values.vox->values[indx];
	    }
	    else {
	      values.core = NULL;
	    }
	    obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			       NULL, NULL, NULL);
	    obj2 = WlzTransposeObj(obj1, &errNum);
	    rtnObj->domain.p->domains[indx] =
	      WlzAssignDomain(obj2->domain, &errNum);
	    if( values.core ){
	      rtnObj->values.vox->values[indx] =
		WlzAssignValues(obj2->values, &errNum);
	    }
	    WlzFreeObj(obj1);
	    WlzFreeObj(obj2);
	  }
	  else {
	    rtnObj->domain.p->domains[indx].core = NULL;
	    if( obj->values.core ){
	      rtnObj->values.vox->values[indx].core = NULL;
	    }
	  }
	}
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
