#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConvexHull_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzConvexHull.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Functions for computing the convex hull of objects.
* \ingroup	WlzConvexHull
*/

#include <stdlib.h>
#include <math.h>
#include <Wlz.h>

static WlzConvHullValues *WlzMakeConvexHullValues(WlzObject *cvh,
						      WlzObject *obj,
						      WlzErrorNum *dstErr);
static WlzConvHullValues *WlzMakeConvexHullValues3d(WlzObject *cvh,
							WlzObject *obj,
							WlzErrorNum *dstErr);


static WlzObject *WlzObjToConvexPolygon3d(WlzObject	*obj,
					  WlzErrorNum	*dstErr);


/*!
* \return	Convex hull object, NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes the convex hull of the given object.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzObjToConvexHull(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*cvh=NULL;
  WlzValues	values;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* the convex hull is a polygon domain with values which are
     a set of chords with pre-calculated parameters which can be used
     by other procedures */
  if((cvh = WlzObjToConvexPolygon(obj, &errNum)) != NULL){
    if((values.c = WlzMakeConvexHullValues(cvh, obj, &errNum)) != NULL){
      /* assign values and reset object type which is now
	 WLZ_CONV_HULL rather than WLZ_2D_POLYGON */
      cvh->values = WlzAssignValues(values, NULL);
      if( cvh->type == WLZ_2D_POLYGON ){
	cvh->type = WLZ_CONV_HULL;
      }
      if( cvh->type == WLZ_3D_DOMAINOBJ ){
	cvh->domain.p->type = WLZ_PLANEDOMAIN_CONV_HULL;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return cvh;
}

/*!
* \return	Convex hull object.
* \ingroup      WlzConvexHull
* \brief	Construct the minimal convex polygonal cover from interval
* 		domain.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzObjToConvexPolygon(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject 		*cvh=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzIVertex2 		*wtx, *w1tx, *w2tx;
  WlzIntervalWSpace 	iwsp;
  int 			nhalfway;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj ){
    switch( obj->type )
    {
      
    case WLZ_2D_DOMAINOBJ:
      /* check the domain */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check the planedomain type */
      if( obj->domain.p ){
	switch( obj->domain.p->type ){
	case WLZ_PLANEDOMAIN_DOMAIN:
	  return WlzObjToConvexPolygon3d(obj, dstErr);

	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_TRANS_OBJ:
      if((obj->values.core) && 
	 (values.obj = WlzObjToConvexPolygon(obj->values.obj, &errNum))){
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
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* now the real algorithm, 2D only, definitely a domain at this point.
     Make a polygon object with the maximum number of vertices for the
     convex polygon = (2 * num_lines + 1) the extra is to allow the polygon
     to be closed i.e. first = last */
  if( errNum == WLZ_ERR_NONE ){
    if((domain.poly = WlzMakePolygonDomain(WLZ_POLYGON_INT, 0, NULL,
				      3+2*(obj->domain.i->lastln -
					   obj->domain.i->line1),
				      1, &errNum))){
      values.core = NULL;
      cvh = WlzMakeMain(WLZ_2D_POLYGON, domain, values, NULL, NULL, &errNum);
    }
  }
  
  if( errNum == WLZ_ERR_NONE ){
    wtx = cvh->domain.poly->vtx;
    /*
     * proceed down right hand side of object
     */
    if( (errNum = WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_ILIC))
       == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
	/*
	 * set up first chord
	 */
	if (iwsp.linpos == obj->domain.i->line1) {
	  if (iwsp.nwlpos == 1) {
	    wtx->vtX = iwsp.lftpos;
	    wtx->vtY = iwsp.linpos;
	    wtx++;
	  }
	  if (iwsp.intrmn == 0) {
	    wtx->vtX = iwsp.rgtpos;
	    wtx->vtY = iwsp.linpos;
	    wtx++;
	    cvh->domain.poly->nvertices = 2;
	    w1tx = wtx-1;
	    w2tx = wtx-2;
	  }
	} else {
	  /*
	     * add extra chords, checking concavity condition
	     */
	  if (iwsp.intrmn == 0) {
	    wtx->vtX = iwsp.rgtpos;
	    wtx->vtY = iwsp.linpos;
	    cvh->domain.poly->nvertices++;
	    /*
	     * Concavity condition (may propagate backwards).
	     * Also deals satisfactorily with the case that first
	     * line consists of a single interval, itself a single point.
	     */
	    while ((cvh->domain.poly->nvertices >= 3) &&
		   (wtx->vtY-w2tx->vtY)*(w1tx->vtX-w2tx->vtX) <=
		   (w1tx->vtY-w2tx->vtY)*(wtx->vtX-w2tx->vtX)) {
	      w1tx->vtX = wtx->vtX;
	      w1tx->vtY = wtx->vtY;
	      wtx--;
	      w1tx--;
	      w2tx--;
	      cvh->domain.poly->nvertices--;
	    }
	    wtx++;
	    w1tx++;
	    w2tx++;
	  }
	}
      }
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    /*
     * now proceed up left hand side of object
     */
    if( (errNum = WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_DLDC))
       == WLZ_ERR_NONE ){
      nhalfway = cvh->domain.poly->nvertices + 2;
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
	if (iwsp.intrmn == 0) {
	  wtx->vtX = iwsp.lftpos;
	  wtx->vtY = iwsp.linpos;
	  cvh->domain.poly->nvertices++;
	  /*
	   * Concavity condition (may propagate backwards).
	   * Also deals satisfactorily with the case that last
	   * line consists of a single interval, itself a single point.
	   */
	  while ((cvh->domain.poly->nvertices >= nhalfway) &&
		 (wtx->vtY-w2tx->vtY)*(w1tx->vtX-w2tx->vtX) <=
		 (w1tx->vtY-w2tx->vtY)*(wtx->vtX-w2tx->vtX)) {
	    w1tx->vtX = wtx->vtX;
	    w1tx->vtY = wtx->vtY;
	    wtx--;
	    w1tx--;
	    w2tx--;
	    cvh->domain.poly->nvertices--;
	  }
	  wtx++;
	  w1tx++;
	  w2tx++;
	}
      }
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return cvh;
}


/*!
* \return
* \brief	
* \param	obj
* \param	dstErr
*/
static WlzObject *WlzObjToConvexPolygon3d(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*polygon=NULL, *obj1, *obj2;
  WlzDomain	domain, *domains, *new_domains;
  WlzValues	values;
  int		p;
  WlzErrorNum	errNum=WLZ_ERR_NONE;


  /* the object and domain have been checked therefore can create the
     new straight away and fill each plane appropriately */
  if((domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_POLYGON,
				    obj->domain.p->plane1,
				    obj->domain.p->lastpl,
				    obj->domain.p->line1,
				    obj->domain.p->lastln,
				    obj->domain.p->kol1,
				    obj->domain.p->lastkl,
				    &errNum)) != NULL){
    domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
    domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
    domain.p->voxel_size[2] = obj->domain.p->voxel_size[2];
    values.core = NULL;
    polygon = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values, NULL, NULL,
			  &errNum);
  }

  if( errNum == WLZ_ERR_NONE ){
    domains = obj->domain.p->domains;
    new_domains = domain.p->domains;
    values.core = NULL;
    for(p=obj->domain.p->plane1; p <= obj->domain.p->lastpl;
	p++, domains++, new_domains++){
      if( (*domains).core ){
	obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
			   NULL, NULL, NULL);
	if((obj2 = WlzObjToConvexPolygon(obj1, &errNum)) != NULL){
	  *new_domains = WlzAssignDomain(obj2->domain, NULL);
	  WlzFreeObj(obj2);
	}
	else {
	  WlzFreeObj(polygon);
	  polygon = NULL;
	  break;
	}
	WlzFreeObj(obj1);
      }
      else {
	(*new_domains).core = NULL;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return polygon;
}


/*!
* \return	New convex hull values.
* \ingroup	WlzConvexHull
* \brief	Fill in parameters of the convex hull into the values table.
*		Compute line equation parameters of chords plus 8*length.
* \param	cvh			Given convex hull.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzConvHullValues *WlzMakeConvexHullValues(
  WlzObject *cvh,
  WlzObject *obj,
  WlzErrorNum *dstErr)
{
  WlzConvHullValues	*cdom=NULL;
  WlzPolygonDomain	*cvhpdom;
  WlzChord 		*chord;
  WlzIVertex2 		*vtx, *wtx;
  int 		i;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* this only gets called if the call to WlzObjToConvexPolygon succeeds
     therefore no need to check objects except for 3D type */
  if( cvh->type == WLZ_3D_DOMAINOBJ ){
    return WlzMakeConvexHullValues3d(cvh, obj, dstErr);
  }

  cvhpdom = cvh->domain.poly;
  /*
   * allocate space
   */
  if((cdom = (WlzConvHullValues *)
             AlcCalloc(1, sizeof(WlzConvHullValues) +
	               (cvhpdom->nvertices-1) * sizeof(WlzChord))) != NULL){
    cdom->ch = (WlzChord *)(cdom + 1);

    cdom->type = WLZ_CONVHULL_VALUES;
    cdom->linkcount = 0;
    cdom->freeptr = NULL;
    cdom->original_table.core = NULL;
    cdom->mdlin = (obj->domain.i->line1 + obj->domain.i->lastln) / 2;
    cdom->mdkol = (obj->domain.i->kol1 + obj->domain.i->lastkl) / 2;
    cdom->nchords = cvhpdom->nvertices - 1;
    cdom->nsigchords = 0;

    chord = cdom->ch;
    vtx = cvhpdom->vtx;
    wtx = vtx + 1;
    for (i=0; i< cdom->nchords; i++) {
      chord->sig = 0;
      chord->acon = wtx->vtY - vtx->vtY;
      chord->bcon = wtx->vtX - vtx->vtX;
      chord->ccon = (wtx->vtX - cdom->mdkol)*chord->acon -
	(wtx->vtY - cdom->mdlin)*chord->bcon;
      chord->cl = 8.0 * sqrt((double)(chord->acon*chord->acon +
				      chord->bcon*chord->bcon));
      chord++;
      vtx++;
      wtx++;
    }
  }
  else {
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return cdom;
}

/*!
* \return	New convex hull values.
* \ingroup	WlzConvexHull
* \brief	Fill in parameters of the convex hull into the values table.
* \param	cvh			Given convex hull.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzConvHullValues *WlzMakeConvexHullValues3d(
  WlzObject *cvh,
  WlzObject *obj,
  WlzErrorNum *dstErr)
{
  WlzValues		rtnvalues, *valuess, values;
  WlzObject		*obj1, *obj2;
  WlzDomain		*domains1, *domains2;
  WlzPixelV		bckgrnd;
  int			p;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* this is only called after a successful call to WlzObjToConvexPolygon
     therefore the given convex hull object and object have been checked
     and match */
  rtnvalues.c = NULL;
  bckgrnd.type = WLZ_GREY_UBYTE;
  bckgrnd.v.ubv = 0;

  /* make a voxeltable and calculate the convex hull values for each plane */
  if((rtnvalues.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_CONV_HULL,
					  cvh->domain.p->plane1,
					  cvh->domain.p->lastpl,
					  bckgrnd, NULL, &errNum)) != NULL){
    domains1 = cvh->domain.p->domains;
    domains2 = obj->domain.p->domains;
    valuess = rtnvalues.vox->values;
    for(p=cvh->domain.p->plane1; p <= cvh->domain.p->lastpl;
	p++, domains1++, domains2++, valuess++){
      if( (*domains1).core != NULL ){
	values.core = NULL;
	obj1 = WlzMakeMain(WLZ_2D_POLYGON, *domains1, values,
			   NULL, NULL, NULL);
	obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains2, values,
			   NULL, NULL, NULL);
	if((values.c = WlzMakeConvexHullValues(obj1, obj2, &errNum)) != NULL){
	  *valuess = WlzAssignValues(values, NULL);
	}
	else {
	  WlzFreeObj(obj2);
	  WlzFreeObj(obj1);
	  WlzFreeVoxelValueTb(rtnvalues.vox);
	  rtnvalues.vox = NULL;
	  break;
	}
	WlzFreeObj(obj2);
	WlzFreeObj(obj1);
      }
    }
  }


  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnvalues.c;
}

