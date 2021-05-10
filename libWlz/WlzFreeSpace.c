#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFreeSpace_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzFreeSpace.c
* \author       Bill Hill, Richard Baldock
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
* \brief	Functions for freeing objects and their components.
* \ingroup	WlzAllocation
*/

#include <stdlib.h>
#include <Wlz.h>


/*! 
* \ingroup      WlzAllocation
* \brief        Free space allocated to a woolz object.
*
* \return       Error number, values: WLZ_ERR_NONE, WLZ_ERR_MEM_FREE
* \param    obj	Object to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreeObj(WlzObject *obj)
{
  int 			i;
  WlzCompoundArray	*ca = (WlzCompoundArray *) obj;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzFreeObj FE %p\n",
	   obj));

  /* check the object pointer and linkcount */
  if (obj == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(obj->linkcount), &errNum) ){    /* Check linkcount */

    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	      ("WlzFreeObj %p WLZ_2D_DOMAINOBJ "
	       "%p %d %p %d %p\n",
	       obj, (obj->domain.i),
	       (obj->domain.i?obj->domain.i->linkcount: 0),
	       (obj->values.core),
	       ((obj->values.core)? obj->values.core->linkcount: 0),
	       (obj->plist)));
      errNum = WlzFreeDomain(obj->domain);
      if((errNum == WLZ_ERR_NONE) && (obj->values.core != NULL)) {
	if(WlzGreyTableIsTiled(obj->values.core->type) == WLZ_GREY_TAB_TILED) {
	  errNum = WlzFreeTiledValues(obj->values.t);
	}
	else {
	  errNum = WlzFreeValueTb(obj->values.v);
	}
      }
      if(errNum == WLZ_ERR_NONE) {
        errNum = WlzFreePropertyList(obj->plist);
      }
      if(errNum == WLZ_ERR_NONE) {
        errNum = WlzFreeObj(obj->assoc);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_3D_DOMAINOBJ %p, "
	       "%d %p %d %p\n",
	       obj, obj->domain.i,
	       (obj->domain.p? obj->domain.p->linkcount: 0),
	       obj->values.core,
	       (obj->values.core? obj->values.core->linkcount: 0),
	       obj->plist));
      errNum = WlzFreePlaneDomain(obj->domain.p);
      if((errNum == WLZ_ERR_NONE) && (obj->values.core != NULL)){
	if(WlzGreyTableIsTiled(obj->values.core->type) == WLZ_GREY_TAB_TILED) {
	  errNum = WlzFreeTiledValues(obj->values.t);
	}
	else {
	  errNum = WlzFreeVoxelValueTb(obj->values.vox);
	}
      }
      if(errNum == WLZ_ERR_NONE) {
        errNum = WlzFreePropertyList(obj->plist);
      }
      if(errNum == WLZ_ERR_NONE) {
        errNum = WlzFreeObj(obj->assoc);
      }
      break;

    case WLZ_TRANS_OBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_TRANS_OBJ %p, "
	       "%d %p %d %p\n",
	       obj, obj->domain.t,
	       ((obj->domain.t)?(obj->domain.t)->linkcount: 0),
	       obj->values.obj,
	       ((obj->values.obj)?(obj->values.obj)->linkcount: 0),
	       obj->plist));
      if( WlzFreeAffineTransform(obj->domain.t) ||
	  WlzFreeObj(obj->values.obj) ||
	  WlzFreePropertyList(obj->plist) ||
	  WlzFreeObj(obj->assoc) ){
	errNum = WLZ_ERR_MEM_FREE;
      }
      break;

    case WLZ_SPLINE:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_SPLINE %p\n",
	       obj, obj->domain.bs));
      errNum = WlzFreeBSpline(obj->domain.bs);
      break;

    case WLZ_2D_POLYGON:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_2D_POLYGON %p\n",
	       obj, obj->domain.poly));
      errNum = WlzFreePolyDmn(obj->domain.poly);
      break;

    case WLZ_BOUNDLIST:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_BOUNDLIST %p\n",
	       obj, obj->domain.b));
      errNum = WlzFreeBoundList(obj->domain.b);
      break;

    case WLZ_CONTOUR:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_CONTOUR %p\n",
	       obj, obj->domain.ctr));
      errNum = WlzFreeContour(obj->domain.ctr);
      break;

    case WLZ_CONV_HULL:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_CONV_HULL %p %p\n",
	       obj, obj->domain.core, obj->values.core));
      if(obj->domain.core) {
        switch(obj->domain.core->type) {
	  case WLZ_CONVHULL_DOMAIN_2D:
	    errNum = WlzFreeConvexHullDomain2(obj->domain.cvh2);
	    break;
	  case WLZ_CONVHULL_DOMAIN_3D:
	    errNum = WlzFreeConvexHullDomain3(obj->domain.cvh3);
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
      }
      break;

    case WLZ_CMESH_2D:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_CMESH_2D %p, "
	       "%d %p %d %p\n",
	       obj, obj->domain.cm2,
	       ((obj->domain.cm2)? (obj->domain.cm2)->linkcount: 0),
	       obj->values.x,
	       ((obj->values.x)? (obj->values.x)->linkcount: 0),
	       obj->plist));
      errNum = WlzCMeshFree2D(obj->domain.cm2);
      if((errNum == WLZ_ERR_NONE) && (obj->values.core != NULL))
      {
	errNum = WlzFreeIndexedValues(obj->values.x);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->plist != NULL))
      {
	errNum = WlzFreePropertyList(obj->plist);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->assoc != NULL))
      {
	errNum = WlzFreeObj(obj->assoc);
      }
      break;
    case WLZ_CMESH_2D5:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_CMESH_2D5 %p, "
	       "%d %p %d %p\n",
	       obj, obj->domain.cm2d5,
	       ((obj->domain.cm2d5)? (obj->domain.cm2d5)->linkcount: 0),
	       obj->values.x,
	       ((obj->values.x)? (obj->values.x)->linkcount: 0),
	       obj->plist));
      errNum = WlzCMeshFree2D5(obj->domain.cm2d5);
      if((errNum == WLZ_ERR_NONE) && (obj->values.core != NULL))
      {
	errNum = WlzFreeIndexedValues(obj->values.x);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->plist != NULL))
      {
	errNum = WlzFreePropertyList(obj->plist);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->assoc != NULL))
      {
	errNum = WlzFreeObj(obj->assoc);
      }
      break;

    case WLZ_CMESH_3D:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_CMESH_3D %p, "
	       "%d %p %d %p\n",
	       obj, obj->domain.cm3,
	       ((obj->domain.cm3)?(obj->domain.cm3)->linkcount: 0),
	       obj->values.x,
	       ((obj->values.x)?(obj->values.x)->linkcount: 0),
	       obj->plist));
      errNum = WlzCMeshFree3D(obj->domain.cm3);
      if((errNum == WLZ_ERR_NONE) && (obj->values.core != NULL))
      {
	errNum = WlzFreeIndexedValues(obj->values.x);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->plist != NULL))
      {
	errNum = WlzFreePropertyList(obj->plist);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->assoc != NULL))
      {
	errNum = WlzFreeObj(obj->assoc);
      }
      break;

    case WLZ_POINTS:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_POINTS %p\n",
	       obj, obj->domain.pts));
      errNum = WlzFreeDomain(obj->domain);
      if((errNum == WLZ_ERR_NONE) && (obj->values.core != NULL))
      {
        errNum = WlzFreePointValues(obj->values.pts);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->plist != NULL))
      {
        errNum = WlzFreePropertyList(obj->plist);
      }
      if((errNum == WLZ_ERR_NONE) && (obj->assoc != NULL))
      {
	errNum = WlzFreeObj(obj->assoc);
      }
      break;

    case WLZ_HISTOGRAM:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_CONV_HULL %p\n",
	       obj, obj->domain.hist));
      errNum = WlzFreeDomain(obj->domain);
      break;

    case WLZ_RECTANGLE:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_RECTANGLE %p\n",
	       obj, obj->domain.r));
      errNum = WlzFreeDomain(obj->domain);
      break;

    case WLZ_AFFINE_TRANS:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_AFFINE_TRANS\n",
	       obj));
      errNum = WlzFreeAffineTransform(obj->domain.t);
      break;

    case WLZ_LUT:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
              ("WlzFreeObj %p WLZ_LUT\n",
	      obj));
      errNum = WlzFreeDomain(obj->domain);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzFreeLUTValues(obj->values);
      }
      break;

    case WLZ_COMPOUND_ARR_1:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_COMPOUND_ARR_1\n",
	       ca));
      for (i=0; i<ca->n; i++){
	if( WlzFreeObj(ca->o[i]) != WLZ_ERR_NONE ){
	  errNum = WLZ_ERR_MEM_FREE;
	}
      }
      break;

    case WLZ_COMPOUND_ARR_2:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_COMPOUND_ARR_2\n",
	       ca));
      for (i=0; i<ca->n; i++){
	if( WlzFreeObj(ca->o[i]) != WLZ_ERR_NONE ){
	  errNum = WLZ_ERR_MEM_FREE;
	}
      }
      break;

    case WLZ_PROPERTY_OBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_PROPERTY_OBJ\n",
	       obj));
      errNum = WlzFreePropertyList(obj->plist);
      break;

    case WLZ_EMPTY_OBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p WLZ_EMPTY_OBJ\n",
	       obj));
      errNum = WlzFreePropertyList(obj->plist);
      break;

    default:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj %p %d\n",
	       obj, (int )(obj->type)));
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    }    /* End of switch */

    AlcFree((void *) obj);
  }

  return( errNum );
}

/*! 
* \ingroup      WlzAllocation
* \brief        Free space allocated to a woolz compound array object.
*
* \return       Error number, values: WLZ_ERR_NONE, WLZ_ERR_MEM_FREE
* \param    	obj			Compound array object to be freed.
*/
WlzErrorNum 	WlzFreeCompoundArray(WlzCompoundArray *obj)
{
  return(WlzFreeObj((WlzObject *)obj));
}

/* function:     WlzFreeIntervalDomain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free an interval domain - convenience link to
 WlzFreeDomain()
*
* \return       Error number, values: from WlzFreeDomain().
* \param    idom	interval domain pointer to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreeIntervalDomain(WlzIntervalDomain *idom)
{
  WlzDomain	domain;

  domain.i = idom;

  return( WlzFreeDomain(domain) );
}

/* function:     WlzFreeHistogramDomain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free a histogram domain.
*
* \return       Error number, values: from WlzFreeDomain().
* \param    hist	Histogram domain to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreeHistogramDomain(WlzHistogramDomain *hist)
{
  WlzDomain	domain;

  domain.hist  = hist;

  return( WlzFreeDomain(domain) );
}

/*! 
* \return       Woolz error code.
* \ingroup      WlzAllocation
* \brief        Free a domain structure of any type. All domain
*		structures must have a type and linkcount. Most
*		also have a freeptr by which, if set, all the space
*		allocated can be freed, however there are some special
*		cases.
* \param    domain			Domain to be freed.
*/
WlzErrorNum WlzFreeDomain(WlzDomain domain)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  AlcErrno	alcErrNum = ALC_ER_NONE;

  if(domain.core != NULL)
  {
    if(WlzUnlink(&(domain.core->linkcount), &errNum))
    {
      switch(domain.core->type)
      {
	case WLZ_CMESH_2D:
	  errNum = WlzCMeshFree2D(domain.cm2);
	  break;
	case WLZ_CMESH_3D:
	  errNum = WlzCMeshFree3D(domain.cm3);
	  break;
	case WLZ_3D_VIEW_STRUCT:
	  errNum = WlzFree3DViewStruct(domain.vs3d);
	  break;
	default:
	  /* Most domains are are freed in the same way. */
	  if(domain.core->freeptr != NULL)
	  {
	    alcErrNum = AlcFreeStackFree(domain.core->freeptr);
	  }
	  AlcFree((void *)domain.core);
	  if(alcErrNum != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_FREE;
	  }
	  break;
      }
    }
  }
  return(errNum);
}

/* function:     WlzFreePlaneDomain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free a planedomain
*
* \return       Error number, values: WLZ_ERR_NONE and errors from
 WlzFreeAffineTransform(), WlzFreeDomain(), WlzFreeBoundList(),
 WlzFreePolyDmn().
* \param    planedm	Pointer to planedomain structure to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreePlaneDomain(WlzPlaneDomain *planedm)
{
  WlzDomain	*domains;
  int 		nplanes;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (planedm == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(planedm->linkcount), &errNum) ){

    nplanes = planedm->lastpl - planedm->plane1 + 1;
    domains = planedm->domains;

    switch( planedm->type ){

    case WLZ_PLANEDOMAIN_DOMAIN:
      while( nplanes-- ){
	errNum |= WlzFreeDomain( *domains );
	domains++;
      }
      break;

    case WLZ_PLANEDOMAIN_POLYGON:
      while( nplanes-- ){
	errNum |= WlzFreePolyDmn((*domains).poly);
	domains++;
      }
      break;

    case WLZ_PLANEDOMAIN_BOUNDLIST:
      while( nplanes-- ){
	errNum |= WlzFreeBoundList((*domains).b);
	domains++;
      }
      break;

    case WLZ_PLANEDOMAIN_HISTOGRAM:
      while( nplanes-- ){
	errNum |= WlzFreeDomain(*domains);
	domains++;
      }
      break;

    case WLZ_PLANEDOMAIN_AFFINE:
      while( nplanes-- ){
	errNum |= WlzFreeAffineTransform((*domains).t);
	domains++;
      }
      break;

    default:
      return( WLZ_ERR_PLANEDOMAIN_TYPE );

    }
    AlcFree((void *) planedm->domains);
    AlcFree((void *) planedm);
  }

  return( errNum );
}

/* function:     WlzFreeValueTb    */
/*! 
* \ingroup      WlzAllocation
* \brief        Convenience routine to free a ragged rect valuetable.
*
* \return       Error number, values: WLZ_ERR_NONE and from WlzFreeValues().
* \param    vdmn	Value domain to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreeValueTb(WlzRagRValues *vdmn)
{       
  WlzValues	values;

  /* check the object pointer and linkcount */
  if (vdmn == NULL){
    return( WLZ_ERR_NONE );
  }

  values.v = vdmn;

  return( WlzFreeValues( values ) );
}

/* function:     WlzFreeValues    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free a values structure, currently only WlzRagRValues
and WlzRectValues DO NOT call this function with any
other values structure types!
*
* \return       Error number, values: WLZ_ERR_NONE, WLZ_ERR_VALUES_DATA.
* \param    values	Values union to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreeValues(WlzValues values)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (values.v == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(values.v->linkcount), &errNum) ){
			
    /* if there is a freeptr then free it */
    if (values.v->freeptr != NULL){

      /* it is illegal for a table to point to itself */
      if( values.v->original_table.v != NULL ){
	return( WLZ_ERR_VALUES_DATA );
      }
      (void )AlcFreeStackFree(values.v->freeptr);

    }
    
    if( values.v->original_table.v ){
      errNum = WlzFreeValues( values.v->original_table );
    }

    AlcFree((void *) values.v);
  }

  return( errNum );
}

/* function:     WlzFreeVoxelValueTb    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free a voxel value table
*
* \return       Error number, values: WLZ_ERR_NONE,
 WLZ_ERR_VOXELVALUES_TYPE and from WlzFreeValues().
* \param    voxtab	
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreeVoxelValueTb(WlzVoxelValues *voxtab)
{
  WlzValues	*values;
  int 		nplanes;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (voxtab == NULL){
    return( WLZ_ERR_NONE );
  }

  /* check the type */
  if( voxtab->type != WLZ_VOXELVALUETABLE_GREY ){
    return WLZ_ERR_VOXELVALUES_TYPE;
  }

  if( WlzUnlink(&(voxtab->linkcount), &errNum) ){

    nplanes = voxtab->lastpl - voxtab->plane1 + 1;
    values = voxtab->values;
    while( nplanes-- ){
      if((errNum = WlzFreeValues(*values)) != WLZ_ERR_NONE)
      {
        break;
      }
      values++;
    }
    WlzFreeVoxelValueTb(voxtab->original_table.vox);
    AlcFreeStackFree(voxtab->freeptr);
    AlcFree((char *) voxtab);
  }

  return( errNum );
}

/* function:     WlzFreePolyDmn    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free a polygon domain.
*
* \return       Error number, values: WLZ_ERR_NONE and from WlzUnlink().
* \param    poly	Polygon domain to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreePolyDmn(WlzPolygonDomain *poly)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (poly == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(poly->linkcount), &errNum) ){
    AlcFree((void *) poly);
  }

  return errNum;
}

/* function:     WlzFreeBoundList    */
/*! 
* \ingroup      WlzAllocation
* \brief        Recursively free a boundary list.
*
* \return       Error number, values: WLZ_ERR_NONE and from WlzUnlink().
* \param    b	Boundary list structure to be freed (note this will call WlzFreeBoundList recursively).
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum WlzFreeBoundList(WlzBoundList *b)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (b == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(b->linkcount), &errNum) ){
    errNum |= WlzFreePolyDmn(b->poly);
    errNum |= WlzFreeBoundList(b->next);
    errNum |= WlzFreeBoundList(b->down);
    AlcFree((void *) b);
  }

  return( errNum );
}

/* function:     WlzFree3DWarpTrans    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free a 3D warp transform.
*
* \return       Error number, values: WLZ_ERR_NONE and from WlzFreePlaneDomain().
* \param    obj	3D warp transform object to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum	WlzFree3DWarpTrans(Wlz3DWarpTrans *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
     errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    if(obj->intptdoms)
    {
      (void )AlcFree(obj->intptdoms);
    }
    if(obj->pdom)
    {
      errNum = WlzFreePlaneDomain(obj->pdom);
    }
    (void )AlcFree(obj);
  }
  return(errNum);
}

/* function:     WlzFreeContour    */
/*! 
* \ingroup      WlzAllocation
* \brief        Free's a WlzContour data structure.
*
* \return       Error number, values: WLZ_ERR_NONE, WLZ_ERR_DOMAIN_NULL, WLZ_ERR_DOMAIN_TYPE and from WlzUnlink().
* \param    ctr				Contour to be freed.
* \par      Source:
*                WlzFreeSpace.c
*/
WlzErrorNum	WlzFreeContour(WlzContour *ctr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctr->type != WLZ_CONTOUR)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if(WlzUnlink(&(ctr->linkcount), &errNum))
    {
      if(ctr->model && WlzUnlink(&(ctr->model->linkcount), &errNum))
      {
        (void )WlzGMModelFree(ctr->model);
      }
      AlcFree((void *)ctr);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Frees an indexed valuetable.
* \param	ixv			Given indexed valuetable.
*/
WlzErrorNum	WlzFreeIndexedValues(WlzIndexedValues *ixv)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ixv == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(ixv->type != (WlzObjectType )WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    (void )AlcVectorFree(ixv->values);
    if(ixv->rank > 0)
    {
      AlcFree(ixv->dim);
    }
    AlcFree(ixv);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Frees a points valuetable.
* \param	pv			Given  points valuetable.
*/
WlzErrorNum	WlzFreePointValues(WlzPointValues *pv)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(pv == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(pv->type != (WlzObjectType )WLZ_POINT_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    (void )AlcFree(pv->values.v);
    if(pv->rank > 0)
    {
      AlcFree(pv->dim);
    }
    AlcFree(pv);
  }
  return(errNum);
}
