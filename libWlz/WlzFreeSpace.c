#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzFreeSpace.c
* Date:         March 1999
* Author:       Bill Hill, Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions to free objects and their data.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 03-03-2K bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzFreeObj						*
*   Date       : Mon Oct 21 19:54:29 1996				*
*************************************************************************
*   Synopsis   :Free space allocated to a woolz object			*
*   Returns    :WlzErrorNum: error return				*
*   Parameters :WlzObject *obj: object to be freed. Note a NULL object	*
*		is assumed legal.					*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzFreeObj(WlzObject *obj)
{
  WlzCompoundArray	*ca = (WlzCompoundArray *) obj;
  WlzErrorNum		errNum = WLZ_ERR_NONE;
  int 			i;

  WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzFreeObj FE 0x%lx\n",
	   (unsigned long )obj));

  /* check the object pointer and linkcount */
  if (obj == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(obj->linkcount), &errNum) ){    /* Check linkcount */

    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	      ("WlzFreeObj 01 0x%lx WLZ_2D_DOMAINOBJ "
	       "0x%lx %d 0x%lx %d 0x%lx\n",
	       (unsigned long )obj,
	       (unsigned long )(obj->domain.i),
	       (obj->domain.i?obj->domain.i->linkcount: 0),
	       (unsigned long )(obj->values.v),
	       ((obj->values.v)? obj->values.v->linkcount: 0),
	       (unsigned long )(obj->plist)));
      if( WlzFreeDomain(obj->domain) ||
	  WlzFreeValueTb(obj->values.v) ||
	  WlzFreeProperty(obj->plist) ||
	  WlzFreeObj(obj->assoc) ){
	errNum = WLZ_ERR_MEM_FREE;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 02 0x%lx WLZ_3D_DOMAINOBJ 0x%lx, "
	       "%d 0x%lx %d 0x%lx\n",
	       (unsigned long )obj, (unsigned long )(obj->domain.i),
	       (obj->domain.p?obj->domain.p->linkcount: 0),
	       (unsigned long )(obj->values.vox),
	       (obj->values.vox?obj->values.vox->linkcount: 0),
	       (unsigned long )(obj->plist)));
      if( WlzFreePlaneDomain(obj->domain.p) ||
	  WlzFreeVoxelValueTb(obj->values.vox) ||
	  WlzFreeProperty(obj->plist) ||
	  WlzFreeObj(obj->assoc) ){
	errNum = WLZ_ERR_MEM_FREE;
      }
      break;

    case WLZ_TRANS_OBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 03 0x%lx WLZ_TRANS_OBJ 0x%lx, "
	       "%d 0x%lx %d 0x%lx\n",
	       (unsigned long )obj,
	       (unsigned long )(obj->domain.t),
	       ((obj->domain.t)?(obj->domain.t)->linkcount: 0),
	       (unsigned long )(obj->values.obj),
	       ((obj->values.obj)?(obj->values.obj)->linkcount: 0),
	       (unsigned long )(obj->plist)));
      if( WlzFreeAffineTransform(obj->domain.t) ||
	  WlzFreeObj(obj->values.obj) ||
	  WlzFreeProperty(obj->plist) ||
	  WlzFreeObj(obj->assoc) ){
	errNum = WLZ_ERR_MEM_FREE;
      }
      break;

    case WLZ_2D_POLYGON:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 05 0x%lx WLZ_2D_POLYGON 0x%lx\n",
	       (unsigned long )obj, (unsigned long )(obj->domain.poly)));
      errNum = WlzFreePolyDmn(obj->domain.poly);
      break;

    case WLZ_BOUNDLIST:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 06 0x%lx WLZ_BOUNDLIST 0x%lx\n",
	       (unsigned long )obj, (unsigned long )(obj->domain.b)));
      errNum = WlzFreeBoundList(obj->domain.b);
      break;

    case WLZ_CONV_HULL:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 07 0x%lx WLZ_CONV_HULL 0x%lx 0x%lx\n",
	       (unsigned long )obj,
	       (unsigned long )(obj->domain.poly),
	       (unsigned long)(obj->values.c)));
      if( WlzFreePolyDmn(obj->domain.poly) ||
	  WlzFreeConvHull(obj->values.c) ){
	errNum = WLZ_ERR_MEM_FREE;
      }
      break;

    case WLZ_HISTOGRAM:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 08 0x%lx WLZ_CONV_HULL 0x%lx\n",
	       (unsigned long )obj, (unsigned long )(obj->domain.hist)));
      errNum = WlzFreeDomain(obj->domain);
      break;

    case WLZ_RECTANGLE:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 09 0x%lx WLZ_RECTANGLE 0x%lx\n",
	       (unsigned long )obj, (unsigned long )(obj->domain.r)));
      errNum = WlzFreeDomain(obj->domain);
      break;

    case WLZ_POINT_INT:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 10 0x%lx WLZ_POINT_INT\n",
	       (unsigned long )obj));
      break;

    case WLZ_POINT_FLOAT:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 11 0x%lx WLZ_POINT_INT\n",
	       (unsigned long )obj));
      break;

    case WLZ_AFFINE_TRANS:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 12 0x%lx WLZ_AFFINE_TRANS\n",
	       (unsigned long )obj));
      errNum = WlzFreeAffineTransform(obj->domain.t);
      break;

    case WLZ_COMPOUND_ARR_1:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 15 0x%lx WLZ_COMPOUND_ARR_1\n",
	       (unsigned long )ca));
      for (i=0; i<ca->n; i++){
	if( WlzFreeObj(ca->o[i]) != WLZ_ERR_NONE ){
	  errNum = WLZ_ERR_MEM_FREE;
	}
      }
      break;

    case WLZ_COMPOUND_ARR_2:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 16 0x%lx WLZ_COMPOUND_ARR_2\n",
	       (unsigned long )ca));
      for (i=0; i<ca->n; i++){
	if( WlzFreeObj(ca->o[i]) != WLZ_ERR_NONE ){
	  errNum = WLZ_ERR_MEM_FREE;
	}
      }
      break;

    case WLZ_PROPERTY_OBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 17 0x%lx WLZ_PROPERTY_OBJ\n",
	       (unsigned long )obj));
      errNum = WlzFreeProperty(obj->plist);
      break;

    case WLZ_EMPTY_OBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 18 0x%lx WLZ_EMPTY_OBJ\n",
	       (unsigned long )obj));
      break;

    default:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
      	      ("WlzFreeObj 18 0x%lx %d\n",
	       (unsigned long )obj, (int )(obj->type)));
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    }    /* End of switch */

    AlcFree((void *) obj);
  }

  return( errNum );
}

/************************************************************************
*   Function   : WlzFreeIntervalDomain					*
*   Date       : Fri Oct 18 22:39:21 1996				*
*************************************************************************
*   Synopsis   :Free an interval domain - convenience link to		*
*		WlzFreeDomain						*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzIntervalDomain *idom: domain to be freed		*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzFreeIntervalDomain(WlzIntervalDomain *idom)
{
  WlzDomain	domain;

  domain.i = idom;

  return( WlzFreeDomain(domain) );
}

/************************************************************************
*   Function   : WlzFreeHistogramDomain					*
*   Date       : Sun Oct 20 10:04:31 1996				*
*************************************************************************
*   Synopsis   :Free a histogram domain					*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzHistogram *hist: the histogram domain to be freed	*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzFreeHistogramDomain(WlzHistogramDomain *hist)
{
  WlzDomain	domain;

  domain.hist  = hist;

  return( WlzFreeDomain(domain) );
}

/************************************************************************
*   Function   : WlzFreeDomain						*
*   Date       : Fri Oct 18 22:39:34 1996				*
*************************************************************************
*   Synopsis   :Free a domain structure of any type. All domain		*
*		structures must have a type, linkcount and freeptr by	*
*		which, if set, all allocated space can be freed. 	*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzDomain domain: domain to be freed			*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzFreeDomain(WlzDomain domain)
{
  WlzErrorNum errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (domain.core == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(domain.core->linkcount), &errNum) ){

    if (domain.core->freeptr != NULL){
      errNum = AlcFreeStackFree(domain.core->freeptr);
    }

    AlcFree((void *) domain.core);
  }

  return errNum;
}

/************************************************************************
*   Function   : WlzFreePlaneDomain					*
*   Date       : Fri Oct 18 22:40:51 1996				*
*************************************************************************
*   Synopsis   :Free a planedomain					*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzPlaneDomain *planedm: domain to be freed		*
*   Global refs:None.							*
************************************************************************/

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

/************************************************************************
*   Function   : WlzFreeValueTb						*
*   Date       : Fri Oct 18 23:15:13 1996				*
*************************************************************************
*   Synopsis   :Convenience routine to free a ragged rect valuetable.	*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzRagRValues *vdmn: valuetable to be freed - type	*
*		WlzRagRValues only					*
*   Global refs:None							*
************************************************************************/

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

/************************************************************************
*   Function   : WlzFreeValues						*
*   Date       : Thu Nov 14 10:50:58 1996				*
*************************************************************************
*   Synopsis   :Free a values structure, currently only WlzRagRValues	*
*		and WlzRectValues DO NOT call this function with any	*
*		other values structure types!				*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzValues values: values table to be freed		*
*   Global refs:None.							*
************************************************************************/

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

/************************************************************************
*   Function   : WlzFreeVoxelValueTb					*
*   Date       : Sun Oct 20 09:47:26 1996				*
*************************************************************************
*   Synopsis   :Free a voxel value table				*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzVoxelValues *voxtab: voxel table to be freed		*
*   Global refs:None.							*
************************************************************************/

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
      errNum |= WlzFreeValues(*values);
      values++;
    }
    AlcFree((char *) voxtab);
  }

  return( errNum );
}

/************************************************************************
*   Function   : WlzFreeConvHull					*
*   Date       : Sun Oct 20 10:18:55 1996				*
*************************************************************************
*   Synopsis   :Free convex hull values.				*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzConvHullValues *c: the values to be freed		*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzFreeConvHull(WlzConvHullValues *c)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (c == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(c->linkcount), &errNum) ){
    AlcFree((void *) c);
  }

  return errNum;
}

/************************************************************************
*   Function   : WlzFreePolyDmn						*
*   Date       : Sun Oct 20 10:21:50 1996				*
*************************************************************************
*   Synopsis   :Free a polygon domain					*
*   Returns    :WlzErrorNum:						*
*   Parameters :WlzPolygonDomain *poly: the polygon domain to be freed	*
*   Global refs:None.							*
************************************************************************/

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

/************************************************************************
*   Function   : WlzFreeBoundList					*
*   Date       : Thu Nov 14 10:59:12 1996				*
*************************************************************************
*   Synopsis   :Recursively free a boundary list			*
*   Returns    :WlzErrorNum: 						*
*   Parameters :WlzBoundList *b: the boundary list			*
*   Global refs:None.							*
************************************************************************/

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

/************************************************************************
*   Function   : WlzFree3DWarpTrans					*
*   Date       : Thu Oct  2 10:05:19 BST 1997				*
*************************************************************************
*   Synopsis   :Free a 3D warp transform.				*
*   Returns    :WlzErrorNum: 						*
*   Parameters :Wlz3DWarpTrans *obj: The 3D warp transform.		*
*   Global refs:None.							*
************************************************************************/
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
