#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzMakeStructs.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Makes Woolz object types.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 15-08-00 bill	Move WlzMakeContour() from WlzContour.c to here.
*		Add WLZ_CONTOUR to WlzMakemain().
* 03-03-00 bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzMakeIntervalDomain					*
*   Date       : Fri Oct 18 15:39:43 1996				*
*************************************************************************
*   Synopsis   :Allocate space for an interval domain structure. If the	*
*		type is WLZ_INTERVALDOMAIN_INTVL then allocate space for*
*		the interval line array and set the pointer.		*
*   Returns    :WlzIntervalDomain *: pointer to the domain structure
*   Parameters :WlzObjectType type: interval domain type		*
*		int l1: first line of required domain			*
*		int ll: last line - this is used with l1 to allocate	*
*		interval line structures.				*
*		int k1, kl: first and last columns			*
*   Global refs:None.							*
************************************************************************/

WlzIntervalDomain *
WlzMakeIntervalDomain(WlzObjectType type,
		      int l1,
		      int ll,
		      int k1,
		      int kl,
		      WlzErrorNum *dstErr)
{
  WlzIntervalDomain	*idom=NULL;
  int 			nlines, intervallinespace;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check bounding box */
  if( (ll < l1) || (kl < k1) )
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }

  if( errNum == WLZ_ERR_NONE ){
    nlines = ll - l1 + 1;
    intervallinespace = nlines * (int) sizeof(WlzIntervalLine);
  
    switch( type ){

    case WLZ_INTERVALDOMAIN_INTVL:
      if( (idom = (WlzIntervalDomain *)
	   AlcCalloc(sizeof(WlzIntervalDomain) + 
		     intervallinespace, 1)) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else {
        idom->intvlines = (WlzIntervalLine *) (idom + 1);
      }
      break;

    case WLZ_INTERVALDOMAIN_RECT:
      if( (idom = (WlzIntervalDomain *)
	   AlcMalloc(sizeof(WlzIntervalDomain))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    idom->type = type;
    idom->line1 = l1;
    idom->lastln = ll;
    idom->kol1 = k1;
    idom->lastkl = kl;
    idom->freeptr = NULL;
    idom->linkcount = 0;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(idom);
}

/************************************************************************
*   Function   : WlzMakePlaneDomain					*
*   Date       : Fri Oct 18 15:40:04 1996				*
*************************************************************************
*   Synopsis   :Allocate space for a plane domain and the domain	*
*		pointers						*
*   Returns    :WlzPlaneDomain *: pointer to initialised domain		*
*   Parameters :WlzObjectType type: must be a valid planedomain type	*
*		int p1, pl: first and last planes.			*
*		int l1, ll: first and last lines.			*
*		int k1, kl: first and last columns.			*
*   Global refs:None.							*
************************************************************************/

WlzPlaneDomain *
WlzMakePlaneDomain(WlzObjectType type,
		   int p1, int pl,
		   int l1, int ll,
		   int k1, int kl,
		   WlzErrorNum *dstErr)
{
  WlzPlaneDomain	*planedm=NULL;
  int			nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check bounding box */
  if( (pl < p1) || (ll < l1) || (kl < k1) ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  
  if( errNum == WLZ_ERR_NONE ){
    nplanes = pl - p1 + 1 ;
    switch(type) {

    case WLZ_PLANEDOMAIN_DOMAIN:
    case WLZ_PLANEDOMAIN_POLYGON:
    case WLZ_PLANEDOMAIN_BOUNDLIST:
    case WLZ_PLANEDOMAIN_HISTOGRAM:
    case WLZ_PLANEDOMAIN_AFFINE:
    case WLZ_PLANEDOMAIN_WARP:
      if( (planedm = (WlzPlaneDomain *)
	   AlcMalloc(sizeof(WlzPlaneDomain))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }

      if( (planedm->domains = (WlzDomain *) 
	   AlcCalloc(nplanes, sizeof(WlzDomain))) == NULL ){
	AlcFree((void *) planedm);
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    planedm->type = type;
    planedm->plane1 = p1;
    planedm->lastpl = pl;
    planedm->line1 = l1;
    planedm->lastln = ll;
    planedm->kol1 = k1;
    planedm->lastkl = kl;
    planedm->freeptr = NULL;
    planedm->linkcount = 0;
    planedm->voxel_size[0] = 1.0;
    planedm->voxel_size[1] = 1.0;
    planedm->voxel_size[2] = 1.0;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(planedm);
}

/************************************************************************
*   Function   : WlzMakeMain						*
*   Date       : Fri Oct 18 15:40:24 1996				*
*************************************************************************
*   Synopsis   :Make a top-level woolz object assigning domain, values	*
*		and other pointers as required . The type is checked.	*
*		The domain is not checked for NULL although this should	*
*		 only apply to the WLZ_EMPTY_OBJ.			*
*   Returns    :WlzObject *: pointer to the initialised object, NULL on *
*		error							*
*   Parameters :WlzObjectType 	type: Object type, one of:		*
*		WLZ_2D_DOMAINOBJ, WLZ_3D_DOMAINOBJ, WLZ_2D_POLYGON,	*
*		WLZ_BOUNDLIST, WLZ_CONV_HULL, WLZ_HISTOGRAM, 		*
*		WLZ_RECTANGLE, WLZ_AFFINE_TRANS, WLZ_PROPERTY_OBJ,	*
*		WLZ_EMPTY_OBJ.						*
*   Global refs:None.							*
************************************************************************/

WlzObject *
WlzMakeMain(WlzObjectType 	type,
	    WlzDomain 		domain,
	    WlzValues 		values,
	    WlzSimpleProperty 	*plist,
	    WlzObject 		*assoc,
	    WlzErrorNum 	*dstErr)
{
  WlzObject 	*obj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( (obj = (WlzObject *) AlcMalloc(sizeof(WlzObject))) == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  /* check the type */
  if( errNum == WLZ_ERR_NONE ){
    switch( type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_2D_POLYGON:
    case WLZ_BOUNDLIST:
    case WLZ_CONV_HULL:
    case WLZ_HISTOGRAM:
    case WLZ_CONTOUR:
    case WLZ_RECTANGLE:
    case WLZ_AFFINE_TRANS:
    case WLZ_PROPERTY_OBJ:
    case WLZ_EMPTY_OBJ:
    case WLZ_TRANS_OBJ:
      obj->type = type;
      obj->linkcount = 0;
      obj->domain = WlzAssignDomain(domain, &errNum);
      if( errNum == WLZ_ERR_NONE ){
	obj->values = WlzAssignValues(values, &errNum);
      }
      if( errNum == WLZ_ERR_NONE ){
	obj->plist = WlzAssignProperty(plist, &errNum);
      }
      if( errNum == WLZ_ERR_NONE ){
	obj->assoc = WlzAssignObject(assoc, &errNum);
      }
      break;

    case WLZ_WARP_TRANS:
    case WLZ_3D_POLYGON:
    case WLZ_3D_WARP_TRANS:
    default:
      AlcFree((void *) obj);
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* debugging statements */
  if( obj ){
    switch (obj->type) {

    case WLZ_2D_DOMAINOBJ:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	      ("Makestructs - Obj 0x%lx type %d idom 0x%lx dl %d val 0x%lx"
	       " vl %d plist 0x%lx ass 0x%lx\n",
	       (unsigned long) obj, obj->type,
	       (unsigned long) (obj->domain.core), 
	       (obj->domain.core ? obj->domain.core->linkcount: 0),
	       (unsigned long) (obj->values.core), 
	       (obj->values.core ? obj->values.core->linkcount: 0),
	       (unsigned long) (obj->plist),
	       (unsigned long) (obj->assoc)));
      break;

    case WLZ_2D_POLYGON:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	      ("Makestructs - Obj 0x%lx type %d pdom 0x%lx ass 0x%lx\n",
	       (unsigned long) obj, obj->type,
	       (unsigned long) (obj->domain.core),
	       (unsigned long) (obj->assoc)));
      break;

    default:
      break;

    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
*   Function   : WlzMakeValueTb						*
*   Date       : Fri Oct 18 15:40:45 1996				*
*************************************************************************
*   Synopsis   :Allocate and initialise space for a ragged-rectangle	*
*		value table only					*
*   Returns    :WlzRagRValues *: pointer to new structure		*
*   Parameters :WlzObjectType	type: type - defines grey type		*
*		int l1, ll: first and last line				*
*		int kl: first column - the width is set later when grey	*
*		values are attached.					*
*		WlzPixelV backgrnd: background pixel value		*
*		WlzObject *orig: originating object for the values (not	*
*		used).
*   Global refs:None.							*
************************************************************************/

WlzRagRValues *
WlzMakeValueTb(WlzObjectType	type,
	       int 		l1,
	       int 		ll,
	       int		k1,
	       WlzPixelV	backgrnd,
	       WlzObject	*orig,
	       WlzErrorNum	*dstErr)
{
  int		nlines;
  WlzRagRValues	*vtb=NULL;
  int 		valuelinespace;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the table type - the grey type is checked later */
  switch( WlzGreyTableTypeToTableType(type, NULL) ){

  case WLZ_GREY_TAB_RAGR:
    break;

  case WLZ_GREY_TAB_RECT:
  case WLZ_GREY_TAB_INTL:
  default:
    errNum = WLZ_ERR_PARAM_DATA;
    break;
  }

  /* check the line bounds */
  if( errNum == WLZ_ERR_NONE ){
    if( ll < l1 ){
      errNum = WLZ_ERR_LINE_DATA;
    }
    else {
      nlines = ll - l1 + 1;
      valuelinespace = nlines * (int) sizeof(WlzValueLine);
    }
  }
  
  /* allocate space */
  if( errNum == WLZ_ERR_NONE ){
    if( (vtb = (WlzRagRValues *)
	 AlcCalloc(sizeof(WlzRagRValues) + valuelinespace, 1)) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      vtb->type = type;
      vtb->line1 = l1;
      vtb->lastln = ll;
      vtb->kol1 = k1;
      vtb->width = 0;
    }
  }

  /* check the grey type and set the background */
  if( errNum == WLZ_ERR_NONE ){
    switch( WlzGreyTableTypeToGreyType(type, NULL) ){

    case WLZ_GREY_INT:
    case WLZ_GREY_SHORT:
    case WLZ_GREY_UBYTE:
    case WLZ_GREY_FLOAT:
    case WLZ_GREY_DOUBLE:
      WlzValueConvertPixel(&(vtb->bckgrnd), backgrnd,
			   WlzGreyTableTypeToGreyType(type, NULL));
      break;

    default:
      AlcFree((void *) vtb);
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    vtb->freeptr = NULL;
    vtb->linkcount = 0;
    vtb->vtblines = (WlzValueLine *) (vtb + 1);
    vtb->original_table.core = NULL;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(vtb);
}

/************************************************************************
*   Function   : WlzMakeVoxelValueTb					*
*   Date       : Fri Oct 18 15:41:08 1996				*
*************************************************************************
*   Synopsis   :Allocate space for a voxel table			*
*   Returns    :WlzVoxelValues *: a voxel values structure		*
*   Parameters :WlzObjectType	type: must be WLZ_VOXELVALUETABLE_GREY	*
*		int p1, pl: first and last plane values			*
*		WlzPixelV backgrnd: background pixel value		*
*		WlzObject *orig: originating object for the values (not	*
*		used)							*
*   Global refs:None.							*
************************************************************************/

WlzVoxelValues *
WlzMakeVoxelValueTb(WlzObjectType	type,
		    int 		p1,
		    int 		pl,
		    WlzPixelV		backgrnd,
		    WlzObject 		*orig,
		    WlzErrorNum		*dstErr)
{
  int 			nplanes, p;
  WlzVoxelValues 	*voxtab=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check type */
  switch( type ){

  case WLZ_VOXELVALUETABLE_GREY:
    break;

  default:
    errNum = WLZ_ERR_PARAM_DATA;
    break;
  }

  /* check plane bounds */
  if( errNum == WLZ_ERR_NONE ){
    if( pl < p1 ){
      errNum = WLZ_ERR_PLANE_DATA;
    }
    else {
      nplanes = pl - p1 + 1;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    if( (voxtab = (WlzVoxelValues *)
	 AlcCalloc(1, sizeof(WlzVoxelValues)	
		   + nplanes*sizeof(WlzValues *))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      voxtab->type = type;
      voxtab->plane1 = p1;
      voxtab->lastpl = pl;
      voxtab->bckgrnd = backgrnd;
      voxtab->freeptr = NULL;
      voxtab->original_table.core = NULL;
      voxtab->linkcount = 0;
      voxtab->values = (WlzValues *) (voxtab + 1);
      for(p=0; p < nplanes; p++){
	voxtab->values[p].core = NULL;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(voxtab);
}

/************************************************************************
*   Function   : WlzMakeRectValueTb					*
*   Date       : Fri Oct 18 15:41:23 1996				*
*************************************************************************
*   Synopsis   :Make rectangular value table attaching values if set.	*
*		Note the values pointer is just copied and will not be	*
*		freed when the object is freed unless the freeptr is set*
*		to the required value.					*
*   Returns    :WlzRectValues *: pointer to the value table.		*
*   Parameters :WlzObjectType	type: must be 	WLZ_GREY_TAB_RECT with	*
*		a valid grey-type.					*
*		int line1, lastln, kol1, width: defines the bounding	*
*		box and are checked.					*
*		WlzPixelV backgrnd: background pixel value		*
*		int	*values: pointer to array of values cast to	*
*		int * for historical reasons. The actual type is encoded*
*		in the type.						*
*   Global refs:None.							*
************************************************************************/

WlzRectValues *
WlzMakeRectValueTb(WlzObjectType	type, 
		   int 			line1,
		   int 			lastln,
		   int 			kol1,
		   int 			width,
		   WlzPixelV		backgrnd,
		   int 			*values,
		   WlzErrorNum		*dstErr)
{
  WlzRectValues *vtb=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the table type - the grey type is checked later */
  switch( WlzGreyTableTypeToTableType(type, NULL) ){

  case WLZ_GREY_TAB_RECT:
    break;

  case WLZ_GREY_TAB_RAGR:
  case WLZ_GREY_TAB_INTL:
  default:
    errNum = WLZ_ERR_PARAM_DATA;
    break;
  }

  /* check the line bounds & width */
  if( errNum == WLZ_ERR_NONE ){
    if( (lastln < line1) || (width <= 0) ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }

  /* allocate space */
  if( errNum == WLZ_ERR_NONE ){
    if( (vtb = (WlzRectValues *) AlcMalloc(sizeof(WlzRectValues)))
       == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      vtb->type = type;
      vtb->freeptr = NULL;
      vtb->linkcount = 0;
      vtb->line1 = line1;
      vtb->lastln = lastln;
      vtb->kol1 = kol1;
      vtb->width = width;
    }
  }

  /* check grey type and set the background */
  if( errNum == WLZ_ERR_NONE ){
    switch( WlzGreyTableTypeToGreyType(type, NULL) ){

    case WLZ_GREY_INT:
    case WLZ_GREY_SHORT:
    case WLZ_GREY_UBYTE:
    case WLZ_GREY_FLOAT:
    case WLZ_GREY_DOUBLE:
      WlzValueConvertPixel(&(vtb->bckgrnd), backgrnd,
			   WlzGreyTableTypeToGreyType(type, NULL));
      break;

    default:
      AlcFree((void *) vtb);
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  if( vtb ){
    vtb->values.inp = values;
    vtb->original_table.core = NULL;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(vtb);
}

/************************************************************************
*   Function   : WlzMakeInterval					*
*   Date       : Fri Oct 18 15:41:40 1996				*
*************************************************************************
*   Synopsis   :Attach and interval pointer to a an interval domain 	*
*		in the appropriate place and set the number of intervals*
*		on the line. Note this procedure makes no checks on the	*
*		arguments because it is often used in deeply nested 	*
*		loops.							*
*   Returns    :WlzErrorNum: always succeeds!				*
*   Parameters :int line: line containing the intervals - must be legal	*
*		WlzIntervalDomain *idom: domain pointer			*
*		int nints: number of intervals				*
*		WlzInterval *intptr: intervals pointer			*
*   Global refs:None							*
************************************************************************/

WlzErrorNum
WlzMakeInterval(int 			line,
		WlzIntervalDomain	*idom,
		int 			nints,
		WlzInterval 		*intptr)
{
  WlzIntervalLine *ivln;

  /* Note no pointer or value checking because this is a key
     procedure at the core of many loops */
  ivln = idom->intvlines + line - idom->line1;

  if( idom->intvlines ){
    ivln->nintvs = nints;
    ivln->intvs = intptr;
  }

  return( WLZ_ERR_NONE );
}

/************************************************************************
*   Function   : WlzMakeValueLine					*
*   Date       : Fri Oct 18 15:42:18 1996				*
*************************************************************************
*   Synopsis   :Attach the grey values for a valueline within a ragged-	*
*		rectangle value table. Not this procedure does not check*
*		the input arguments because it is often at the core of	*
*		nested loops.						*
*   Returns    :WlzErrorNum: always succeeds				*
*   Parameters :WlzRagRValues 	*vtb: the value table pointer		*
*		int	line: line for the values to be set.		*
*		int k1, kl: column interval for which the greys are set *
*		int *greyptr: grey values pointer cast type int *	*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum 
WlzMakeValueLine(WlzRagRValues 	*vtb,
		 int 		line,
		 int 		k1,
		 int 		kl,
		 int 		*greyptr)
{
  WlzValueLine *vlln;

  /* Note no pointer or value checking because this is a key
     procedure at the core of many loops */
  vlln = vtb->vtblines + line - vtb->line1;
  vlln->vkol1 = k1 - vtb->kol1;
  vlln->vlastkl = kl - vtb->kol1;
  vlln->values.inp = greyptr;

  return( WLZ_ERR_NONE );
}

/************************************************************************
*   Function   : WlzMakePolyDmn						*
*   Date       : Fri Oct 18 15:42:42 1996				*
*************************************************************************
*   Synopsis   :Make a polygon domain, allocating space and copying as	*
*		required:						*
*		vertices != NULL, copy=0 - just plant the pointer	*
*		vertices != NULL, copy=1 - allocate space and copy	*
*		vertices == NULL, copy=0 - nor vertex space allocated -	*
*			probably an error!!				*
*		vertices == NULL, copy=1 - allocate space for maxv 	*
*			vertices.					*
*   Returns    :WlzPolygonDomain *: pointer to initialised structure.	*
*   Parameters :WlzObjectType	type: one of WLZ_POLYGON_INT, 		*
*		WLZ_POLYGON_FLOAT, WLZ_POLYGON_DOUBLE			*
*		WlzIVertex2 	*vertices: vertices to be set, see above*
*		for value options.					*
*		int n: number of vertices if vertices!=NULL		*
*		int maxv: size of array if vertices!=NULL else number	*
*		of vertices for which space it to be allocated.		*
*		int copy: copy flag see above for values.		*
*   Global refs:None.							*
************************************************************************/

WlzPolygonDomain *
WlzMakePolyDmn(WlzObjectType	type,
	       WlzIVertex2 	*vertices,
	       int 		n,
	       int 		maxv,
	       int 		copy,
	       WlzErrorNum		*dstErr)
{
  WlzPolygonDomain	*p=NULL;
  int 			vertexsize;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check type and set vertex size */
  switch( type ){

  case WLZ_POLYGON_INT:
    vertexsize = sizeof(WlzIVertex2);
    break;

  case WLZ_POLYGON_FLOAT:
    vertexsize = sizeof(WlzFVertex2);
    break;

  case WLZ_POLYGON_DOUBLE:
    vertexsize = sizeof(WlzDVertex2);
    break;

  default:
    errNum = WLZ_ERR_PARAM_DATA;
    break;
  }
  
  if (copy == 0){
    vertexsize = 0;
  }

  if( errNum == WLZ_ERR_NONE ){
    if( (p = (WlzPolygonDomain *)
	 AlcCalloc(sizeof(WlzPolygonDomain) + maxv*vertexsize, 1))
       == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      p->type = type;
      p->nvertices = n;
      p->maxvertices = maxv;
      p->linkcount = 0;
      if( copy == 0 ){
	p->vtx = vertices;
      } 
      else {
	p->vtx = (WlzIVertex2 *) (p + 1);
	memcpy((void *) p->vtx, (void *) vertices, n*vertexsize);
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(p);
}

/************************************************************************
*   Function   : WlzMakeBoundList					*
*   Date       : Fri Oct 18 15:42:58 1996				*
*************************************************************************
*   Synopsis   :Allocate space and initialise a boundlist structure.	*
*   Returns    :WlzBoundList *: pointer to the structre.		*
*   Parameters :WlzObjectType type: boundlist type - WLZ_BOUNDLIST_PIECE*
*		or 	WLZ_BOUNDLIST_HOLE				*
*		int wrap: number of vertices by which the polygon is	*
*		"wrapped" ie number of vertices overlapping at the	*
*		beginning and end.					*
*		WlzPolygonDomain *poly: polygon for this boundary struct*
*   Global refs:None.							*
************************************************************************/

WlzBoundList *
WlzMakeBoundList(WlzObjectType		type,
		 int			wrap,
		 WlzPolygonDomain 	*poly,
		 WlzErrorNum		*dstErr)
{
  WlzBoundList *blist=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check type */
  if( type != WLZ_BOUNDLIST_PIECE && type != WLZ_BOUNDLIST_HOLE ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  /* check wrap */
  if( wrap < 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  /* allocate space */
  if( errNum == WLZ_ERR_NONE ){
    if( (blist = ( WlzBoundList *) AlcMalloc(sizeof( WlzBoundList)))
       == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      blist->type = type;
      blist->linkcount = 0;
      blist->freeptr = NULL;
      blist->up = NULL;
      blist->next = NULL;
      blist->down = NULL;
      blist->wrap = wrap;
      if(poly &&
	 ((blist->poly = WlzAssignPolygonDomain(poly, &errNum)) == NULL) ){
	AlcFree((void *) blist);
	blist = NULL;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( blist );
}

/************************************************************************
*   Function   : WlzMakeIVertex						*
*   Date       : Fri Oct 18 15:43:24 1996				*
*************************************************************************
*   Synopsis   :Make an integer vertex array				*
*   Returns    :WlzIVertex2 *: pointer to the array.			*
*   Parameters :int nverts: number of vertices.				*
*   Global refs:None.							*
************************************************************************/

WlzIVertex2 *WlzMakeIVertex(int nverts,
			   WlzErrorNum		*dstErr)
{
  WlzIVertex2	*vtx=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
	
  if( (vtx = (WlzIVertex2 *)
       AlcMalloc(sizeof(WlzIVertex2) * nverts)) == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }
    
  if( dstErr ){
    *dstErr = errNum;
  }
  return(vtx);
}

/************************************************************************
*   Function   : WlzMakeRect						*
*   Date       : Fri Oct 18 15:43:45 1996				*
*************************************************************************
*   Synopsis   :make a top-level rectangular object, setting values	*
*		if non-NULL, uses WlzMakeRectValueTb to assign values	*
*		which by default will not be freed when the object is	*
*		freed. The freeptr needs to be explicitly set by the	*
*		calling procedure.					*
*		This is a convenience procedure calling 		*
*		WlzMakeIntervalDomain then WlzMakeRectValueTb then	*
*		WlzMakeMain.						*
*   Returns    :WlzObject *: pointer to the top-level object		*
*   Parameters :int line1, lastln, kol1, lastkl: bounding box for the	*
*		data.							*
*		WlzGreyType pixeltype: data pixel type.			*
*		int *grey_values: pointer to grey-values array		*
*		WlzPixelV backgrnd: background pixel value.		*
*		WlzSimpleProperty *prop: associated properties		*
*		WlzObject *assoc_obj: associated object			*
*   Global refs:None.							*
************************************************************************/

WlzObject *WlzMakeRect(int 			line1,
		       int 			lastln,
		       int 			kol1,
		       int 			lastkl,
		       WlzGreyType		pixeltype,
		       int 			*grey_values,
		       WlzPixelV		backgrnd,
		       WlzSimpleProperty	*prop,
		       WlzObject		*assoc_obj,
		       WlzErrorNum		*dstErr)
{
  WlzDomain	domain;
  WlzValues	values;
  WlzObject 	*obj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				   line1, lastln,
				   kol1, lastkl, &errNum);

  if((errNum == WLZ_ERR_NONE) && 
     ((values.r = 
       WlzMakeRectValueTb(WlzGreyTableType(WLZ_GREY_TAB_RECT, pixeltype,
					   NULL),
			  line1, lastln, kol1, lastkl - kol1 + 1,
			  backgrnd, grey_values, &errNum)) == NULL) ){
    WlzFreeDomain( domain );
  }

  if((errNum == WLZ_ERR_NONE) &&
     ((obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			 prop, assoc_obj, &errNum)) == NULL) ){
    WlzFreeDomain( domain );
    WlzFreeValues( values );
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
*   Function   : WlzMakeHistogramDomain					*
*   Date       : Tue May 20 11:29:03 BST 1997				*
*************************************************************************
*   Synopsis   :Allocates space for a histogram domain with the space	*
*		the given maximum number of bins of the appropriate	*
*		type.							*
*   Returns    :WlzHistogramDomain *: Pointer to object, NULL on error.	*
*   Parameters :WlzObjectType type: Can be WLZ_HISTOGRAMDOMAIN_INT or	*
*		WLZ_HISTOGRAMDOMAIN_FLOAT only.				*
*		int maxBins: Maximum number of histogram bins.		*
*   Global refs:None.							*
************************************************************************/
WlzHistogramDomain *WlzMakeHistogramDomain(WlzObjectType type, int maxBins,
					   WlzErrorNum		*dstErr)
{
  WlzHistogramDomain *hist = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  AlcErrno	alcErr = ALC_ER_NONE;

  if(((type != WLZ_HISTOGRAMDOMAIN_INT) &&
      (type != WLZ_HISTOGRAMDOMAIN_FLOAT)) || (maxBins < 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((hist = (WlzHistogramDomain *)AlcCalloc(sizeof(WlzHistogramDomain),
  					          1)) == NULL )
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else if(maxBins)
  {
    switch(type)
    {
      case WLZ_HISTOGRAMDOMAIN_INT:
	if((hist->binValues.inp = (int *)AlcCalloc(sizeof(int),
						   maxBins)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  hist->freeptr = AlcFreeStackPush(hist->freeptr,
					   (void *)(hist->binValues.inp),
					   &alcErr);
	  if(alcErr != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	break;
      case WLZ_HISTOGRAMDOMAIN_FLOAT:
	if((hist->binValues.dbp = (double *)AlcCalloc(sizeof(double),
						      maxBins)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  hist->freeptr = AlcFreeStackPush(hist->freeptr, 
					   (void *)(hist->binValues.dbp),
					   &alcErr);
	  if(alcErr != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    hist->type = type;
    hist->maxBins = maxBins;
  }
  else
  {
    if(hist)
    {
      AlcFree(hist);
      hist = NULL;
    }
  }
  if( dstErr ){
    *dstErr = errNum;
  }
  return(hist);
}

/************************************************************************
*   Function   : WlzMakeEmpty						*
*   Date       : Fri Mar  7 14:04:49 GMT 1997				*
*************************************************************************
*   Synopsis   :Convenience procedure to make a top-level empty object.	*
*   Returns    :WlzObject *: NULL on error.				*
*   Parameters :void							*
*   Global refs:None.							*
************************************************************************/

WlzObject *WlzMakeEmpty(WlzErrorNum *dstErr)
{
  WlzDomain	domain;
  WlzValues	values;

  domain.core = NULL;
  values.core = NULL;
  return WlzMakeMain(WLZ_EMPTY_OBJ, domain, values, NULL, NULL, dstErr);
}

/************************************************************************
*   Function   : WlzMakeHistogram					*
*   Date       : Tue May 20 11:38:11 BST 1997				*
*************************************************************************
*   Synopsis   :Convenience procedure to make a top-level object with	*
*		a histogram domain.					*
*   Returns    :WlzObject *: NULL on error.				*
*   Parameters :WlzObjectType type: can be WLZ_HISTOGRAMDOMAIN_INT or	*
*		WLZ_HISTOGRAMDOMAIN_FLOAT only.				*
*		int maxBins: Maximum number of histogram bins.		*
*   Global refs:None.							*
************************************************************************/
WlzObject	*WlzMakeHistogram(WlzObjectType	type, int maxBins,
				  WlzErrorNum		*dstErr)
{
  WlzObject	*obj = NULL;
  WlzDomain	domain;
  WlzValues	values;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if((domain.hist = WlzMakeHistogramDomain(type, maxBins, &errNum)) != NULL)
  {
    values.core = NULL;
    obj = WlzMakeMain(WLZ_HISTOGRAM, domain, values, NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
*   Function   : WlzNewGrey						*
*   Date       : Fri Oct 18 15:44:40 1996				*
*************************************************************************
*   Synopsis   :Make a top-level object with the same domain as iobj 	*
*		(pointer copied and linkcount incremented) and new	*
*		valuetable with values copied from iobj. If iobj has no	*
*		valuetable then the returned object will have no value-	*
*		table. This allows a copy of a 2D grey-value object	*
*   Returns    :WlzObject *: NULL on error				*
*   Parameters :WlzObject *iobj: object to be copied			*
*   Global refs:None.							*
************************************************************************/

WlzObject *WlzNewGrey(WlzObject *iobj,
		      WlzErrorNum		*dstErr)
{
  WlzObject		*jobj=NULL;
  WlzValues 		v;
  WlzGreyP 		g1, g2;
  WlzObjectType		vtype;
  WlzIntervalWSpace	iwsp1,iwsp2;
  WlzGreyWSpace		gwsp1,gwsp2;
  WlzPixelV 		backgrnd;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check input object */
  if( iobj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( iobj->type ){

    case WLZ_2D_DOMAINOBJ:
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( (errNum == WLZ_ERR_NONE) && (iobj->values.v == NULL) ){
    return WlzMakeMain(iobj->type, iobj->domain, iobj->values,
		       iobj->plist, iobj, dstErr);
  }

  /* get the background value */
  if( errNum == WLZ_ERR_NONE ){
    backgrnd = WlzGetBackground( iobj, &errNum );
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( iobj->domain.core->type ){

    case WLZ_INTERVALDOMAIN_INTVL:
      vtype =
	WlzGreyTableType(
	  WLZ_GREY_TAB_RAGR, 
	  WlzGreyTableTypeToGreyType(iobj->values.core->type, NULL), NULL);
      break;

    case WLZ_INTERVALDOMAIN_RECT:
      vtype = 
	WlzGreyTableType(
	  WLZ_GREY_TAB_RECT, 
	  WlzGreyTableTypeToGreyType(iobj->values.core->type, NULL), NULL);
      break;

    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }

  if((errNum == WLZ_ERR_NONE) &&
     (v.v = WlzNewValueTb(iobj, vtype, backgrnd, &errNum)) ){
    jobj = WlzMakeMain(iobj->type, iobj->domain, v, iobj->plist,
		       iobj, &errNum);
  }

  if((errNum == WLZ_ERR_NONE) &&
     ((errNum = WlzInitGreyScan(iobj, &iwsp1, &gwsp1)) == WLZ_ERR_NONE) &&
     ((errNum = WlzInitGreyScan(jobj, &iwsp2, &gwsp2)) == WLZ_ERR_NONE))
  {
    while( (errNum = WlzNextGreyInterval(&iwsp1)) == WLZ_ERR_NONE ){
      size_t	data_size;

      (void) WlzNextGreyInterval(&iwsp2);
      g1 = gwsp1.u_grintptr;
      g2 = gwsp2.u_grintptr;

      switch( gwsp1.pixeltype ){

      case WLZ_GREY_INT:
	data_size = sizeof(int);
	break;

      case WLZ_GREY_SHORT:
	data_size = sizeof(short);
	break;

      case WLZ_GREY_UBYTE:
	data_size = sizeof(UBYTE);
	break;

      case WLZ_GREY_FLOAT:
	data_size = sizeof(float);
	break;

      case WLZ_GREY_DOUBLE:
	data_size = sizeof(double);
	break;

      default:
	WlzFreeObj( jobj );
	return( NULL );

      }
      data_size *= (iwsp1.rgtpos-iwsp1.lftpos+1);
      memcpy((void *) g2.ubp, (void *) g1.ubp, data_size);
    }
    switch( errNum ){
    case WLZ_ERR_NONE:
    case WLZ_ERR_EOO:
      errNum = WLZ_ERR_NONE;
      break;
    default:
      (void) WlzFreeObj(jobj);
      jobj = NULL;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(jobj);
}

/************************************************************************
*   Function   : WlzNewValueTb						*
*   Date       : Fri Oct 18 15:44:55 1996				*
*************************************************************************
*   Synopsis   :Creat a value table of required type with the same size	*
*		and shape as the domain of obj. This must be of type	*
*		WLZ_2D_DOMAINOBJ.					*
*   Returns    :WlzRagRValues *: pointer to the created table.		*
*   Parameters :WlzObject *obj: object pointer whose domain is matched	*
*		WlzObjectType type: type for the valuetable		*
*		WlzPixelV backgrnd: background pixel value		*
*   Global refs:None.							*
************************************************************************/

WlzRagRValues *WlzNewValueTb(WlzObject		*obj,
			     WlzObjectType	type,
			     WlzPixelV		backgrnd,
			     WlzErrorNum	*dstErr)
{
  WlzValues 		v;
  WlzDomain		idom;
  WlzIntervalWSpace	iwsp;
  WlzGreyP		g;
  int 			k1, table_size, bgd_val;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  v.v = NULL;
  if( obj == NULL || obj->domain.core == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) && (obj->type != WLZ_2D_DOMAINOBJ) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else {
    idom = obj->domain;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( WlzGreyTableTypeToTableType(type, &errNum) ){

    case WLZ_GREY_TAB_RAGR:
      switch( WlzGreyTableTypeToGreyType(type, &errNum) ){

      case WLZ_GREY_INT:
	table_size = WlzLineArea(obj, NULL) * sizeof(int);
	bgd_val = (int) backgrnd.v.inv;
	break;

      case WLZ_GREY_SHORT:
	table_size = WlzLineArea(obj, NULL) * sizeof(short);
	bgd_val = (int) backgrnd.v.shv;
	break;

      case WLZ_GREY_UBYTE:
	table_size = WlzLineArea(obj, NULL) * sizeof(UBYTE);
	bgd_val = (int) backgrnd.v.ubv;
	break;

      case WLZ_GREY_FLOAT:
	table_size = WlzLineArea(obj, NULL) * sizeof(float);
	bgd_val = (int) backgrnd.v.flv;
	break;

      case WLZ_GREY_DOUBLE:
	table_size = WlzLineArea(obj, NULL) * sizeof(double);
	bgd_val = (int) backgrnd.v.dbv;
	break;

      default:
	break;
      }

      if( errNum == WLZ_ERR_NONE ){
	if( v.v = WlzMakeValueTb(type, idom.i->line1, idom.i->lastln,
				 idom.i->kol1, backgrnd, obj, &errNum) ){
	  if( (g.inp = (int *) AlcMalloc(table_size)) == NULL ){
	    WlzFreeValueTb(v.v);
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else {
	    memset((void *) g.inp, bgd_val, table_size);
	    v.v->freeptr = AlcFreeStackPush(v.v->freeptr, (void *)g.inp, NULL);
	    v.v->width = idom.i->lastkl - idom.i->kol1 + 1;
	  }
	}
      }

      if( errNum == WLZ_ERR_NONE ){
	errNum = WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_ILIC);
	while( (errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
	  if (iwsp.nwlpos){
	    k1 = iwsp.lftpos;
	  }
	  if (iwsp.intrmn == 0) {
	    WlzMakeValueLine(v.v, iwsp.linpos, k1, iwsp.rgtpos, g.inp);
	    switch( WlzGreyTableTypeToGreyType(type, NULL) ){

	    case WLZ_GREY_INT:
	      g.inp += iwsp.rgtpos - k1 +1;
	      break;

	    case WLZ_GREY_SHORT:
	      g.shp += iwsp.rgtpos - k1 +1;
	      break;

	    case WLZ_GREY_UBYTE:
	      g.ubp += iwsp.rgtpos - k1 +1;
	      break;

	    case WLZ_GREY_FLOAT:
	      g.flp += iwsp.rgtpos - k1 +1;
	      break;

	    case WLZ_GREY_DOUBLE:
	      g.dbp += iwsp.rgtpos - k1 +1;
	      break;

	    }
	  }
	}
	switch( errNum ){
	case WLZ_ERR_NONE:
	case WLZ_ERR_EOO:
	  errNum = WLZ_ERR_NONE;
	  break;
	default:
	  WlzFreeValueTb(v.v);
	  v.v = NULL;
	  break;
	}
      }
      break;

    case WLZ_GREY_TAB_RECT:
      if( (v.r = WlzMakeRectValueTb(type, idom.i->line1, idom.i->lastln,
				    idom.i->kol1,
				    idom.i->lastkl - idom.i->kol1 + 1,
				    backgrnd, NULL, &errNum)) == NULL ){
	break;
      }

      switch( WlzGreyTableTypeToGreyType(type, NULL) ){

      case WLZ_GREY_INT:
	v.r->values.inp  =
	  (int *) AlcMalloc(sizeof(int)*(v.r->lastln - v.r->line1 + 1)
			    * v.r->width);
	break;

      case WLZ_GREY_SHORT:
	v.r->values.shp  =
	  (short *) AlcMalloc(sizeof(short)*(v.r->lastln - v.r->line1 + 1)
			      * v.r->width);
	break;

      case WLZ_GREY_UBYTE:
	v.r->values.ubp  =
	  (UBYTE *) AlcMalloc(sizeof(UBYTE)*(v.r->lastln - v.r->line1 + 1)
			      * v.r->width);
	break;

      case WLZ_GREY_FLOAT:
	v.r->values.flp  = 
	  (float *) AlcMalloc(sizeof(float)*(v.r->lastln - v.r->line1 + 1)
			      * v.r->width);
	break;

      case WLZ_GREY_DOUBLE:
	v.r->values.dbp  =
	  (double *) AlcMalloc(sizeof(double)*(v.r->lastln - v.r->line1 + 1)
			       * v.r->width);
	break;

      }
    
      if( v.r->values.ubp == NULL ){
	WlzFreeValueTb(v.v);
	v.v = NULL;
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else {
	v.r->freeptr = AlcFreeStackPush(v.r->freeptr, (void *)v.r->values.inp,
				        NULL);
      }
      break;

    case WLZ_GREY_TAB_INTL:
      v.i = WlzMakeIntervalValues(type, obj, backgrnd, &errNum);
      break;

    default:
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(v.v);
}

/************************************************************************
*   Function   : WlzNewIDomain						*
*   Date       : Fri Oct 18 15:45:09 1996				*
*************************************************************************
*   Synopsis   :Make a copy of an intervaldomain			*
*   Returns    :WlzIntervalDomain *: NULL on error			*
*   Parameters :WlzObjectType outDomType: Required interval domain	*
*					  type.				*
*		WlzIntervalDomain *inDom: Domain to be copied.		*
*		WlzErrorNum *dstErr:	Destination error pointer, 	*
*					may be NULL.			*
*   Global refs:None.							*
************************************************************************/

WlzIntervalDomain *
WlzNewIDomain(WlzObjectType outDomType,
	      WlzIntervalDomain *inDom,
	      WlzErrorNum	*dstErr)
{
  WlzIntervalDomain	*outDom = NULL;
  WlzInterval		*itvl,
  			*jtvl,
			*ktvl;
  WlzIntervalLine	*ivLn;
  int			line,
  			ivCnt,
			n;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* check the domain */
  if(inDom == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(((inDom->type != WLZ_INTERVALDOMAIN_INTVL) &&
           (inDom->type != WLZ_INTERVALDOMAIN_RECT)) ||
          ((outDomType != WLZ_INTERVALDOMAIN_INTVL) &&
           (outDomType != WLZ_INTERVALDOMAIN_RECT)))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  /* create the new domain */
  if( errNum == WLZ_ERR_NONE )
  {
    outDom = WlzMakeIntervalDomain(outDomType, inDom->line1, inDom->lastln,
				   inDom->kol1, inDom->lastkl, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(outDomType == WLZ_INTERVALDOMAIN_INTVL)
    {
      /* Count the intervals. */
      if(inDom->type == WLZ_INTERVALDOMAIN_INTVL)
      {
        ivLn = inDom->intvlines;
	ivCnt = 0;
	for(line = inDom->line1; line <= inDom->lastln; ++line)
	{
	  ivCnt += ivLn->nintvs;
	  ++ivLn;
	}
      }
      else /* inDom->type == WLZ_INTERVALDOMAIN_RECT */
      {
        ivCnt = inDom->lastln - inDom->line1 + 1;
      }
      /* Allocate space for the intervals. */
      if((itvl = (WlzInterval *)AlcMalloc(ivCnt * sizeof(WlzInterval))) == NULL)
      {
	WlzFreeIntervalDomain(outDom);
	errNum = WLZ_ERR_MEM_ALLOC;
      } 
      else
      {
	jtvl = itvl;
	outDom->freeptr = AlcFreeStackPush(outDom->freeptr, (void *)itvl, NULL);
	/* Set the new intervals. */
	if(inDom->type == WLZ_INTERVALDOMAIN_INTVL)
	{
	  ivLn = inDom->intvlines;
	  for(line = inDom->line1; line <= inDom->lastln; ++line)
	  {
	    ktvl = ivLn->intvs;
	    for(n = 0; n<ivLn->nintvs; ++n)
	    {
	      jtvl->ileft = ktvl->ileft;
	      jtvl->iright = ktvl->iright;
	      ++jtvl;
	      ++ktvl;
	    }
	    (void )WlzMakeInterval(line, outDom, ivLn->nintvs, itvl);
	    itvl = jtvl;
	    ++ivLn;
	  }
	}
	else /* inDom->type == WLZ_INTERVALDOMAIN_RECT */
	{
	  for(line = inDom->line1; line <= inDom->lastln; ++line)
	  { 
	    jtvl->ileft = inDom->kol1;
	    jtvl->iright = inDom->lastkl;
	    (void )WlzMakeInterval(line, outDom, 1, jtvl);
	    ++jtvl;
	  }
	}
      }
    }
  }
  if( dstErr )
  {
    *dstErr = errNum;
  }
  return(outDom);
}

/************************************************************************
* Function:	WlzMakeContour
* Returns:	WlzContour *:		New contour.
* Purpose:	Makes a new contour data structure.
* Global refs:	-
* Parameters:	WlzErrorNum *dstErr:	Destination error pointer, may
					be null.
************************************************************************/
WlzContour	*WlzMakeContour(WlzErrorNum *dstErr)
{
  WlzContour	*ctr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ctr = AlcCalloc(1, sizeof(WlzContour))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    ctr->type = WLZ_CONTOUR;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}
