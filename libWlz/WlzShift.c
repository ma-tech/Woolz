#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzShift.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz functions for shifting (applying integer
*		translations) to Woolz objects, domains and values.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzShiftObject						*
* Returns:	WlzObject *:		Shifted object.			*
* Purpose:	The external object shift interface function.		*
*		Shifts a Woolz object in place, cf WlzAffineTransform()	*
*		which always creates a new object with both a new	*
*		domain and a new value table.				*
*		WlzShiftObject always makes a new domain but keeps	*
*		as much of the given object's value table as possible.	*
* Global refs:	-							*
* Parameters:	WlzObject *inObj:	The given object.		*
*		int xShift:		Column shift.			*
*		int yShift:		Line shift.			*
*		int zShift:		Plane shift (only used for 3D	*
*					objects).			*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
WlzObject	*WlzShiftObject(WlzObject *inObj,
				int xShift, int yShift, int zShift,
				WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*outObj = NULL;

  dom.core = NULL;
  val.core = NULL;
  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(inObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_TRANS_OBJ:    /* FALLTHROUGH */
      case WLZ_AFFINE_TRANS: /* FALLTHROUGH */
      case WLZ_PROPERTY_OBJ: /* FALLTHROUGH */
      case WLZ_2D_POLYGON:   /* FALLTHROUGH */
      case WLZ_BOUNDLIST:
	dom = WlzShiftDomain(inObj->type, inObj->domain,
			     xShift, yShift, zShift, &errNum);
	if(inObj->values.core)
	{
	  val = WlzShiftValues(inObj->type, inObj->values, inObj->domain,
	  		       xShift, yShift, zShift, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzMakeMain(inObj->type, dom, val, NULL, NULL, &errNum);
	}
	break;
      case WLZ_EMPTY_OBJ:
        outObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_HISTOGRAM:
      case WLZ_CONV_HULL:
      case WLZ_3D_WARP_TRANS:
      case WLZ_3D_POLYGON:
      case WLZ_RECTANGLE:
      case WLZ_VECTOR_INT:
      case WLZ_VECTOR_FLOAT:
      case WLZ_POINT_INT:
      case WLZ_POINT_FLOAT:
      case WLZ_CONVOLVE_INT:
      case WLZ_CONVOLVE_FLOAT:
      case WLZ_WARP_TRANS:
      case WLZ_FMATCHOBJ:
      case WLZ_TEXT:
      case WLZ_COMPOUND_ARR_1:
      case WLZ_COMPOUND_ARR_2:
      case WLZ_COMPOUND_LIST_1:
      case WLZ_COMPOUND_LIST_2:
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dom.core)
    {
      (void )WlzFreeDomain(dom);
    }
    if(val.core)
    {
      (void )WlzFreeValues(val);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(outObj);
}

/************************************************************************
* Function:	WlzShiftDomain						*
* Returns:	WlzDomain:		Shifted domain, NULL on error.	*
* Purpose:	Creates a new shifted domain.				*
* Global refs:	-							*
* Parameters:	WlzObjectType inObjType: Type of given domain's parent	*
*					object.				*
*		WlzDomain inDom:	Domain to be shifted.		*
*		int xShift:		Column shift.			*
*		int yShift:		Line shift.			*
*		int zShift:		Plane shift (only used for 3D	*
*					objects).			*
*		WlzErrorNum *dstErr:	Destination error pointer, may	*
*					be NULL.			*
************************************************************************/
WlzDomain	 WlzShiftDomain(WlzObjectType inObjType, WlzDomain inDom,
			        int xShift, int yShift, int zShift,
			        WlzErrorNum *dstErr)
{
  int		idx0;
  WlzDomain	tDom0,
  		tDom1,
  		outDom,
  		nullDom;
  WlzDomain	*domP0,
  		*domP1;
  WlzIVertex2	*tIVP0;
  WlzFVertex2	*tFVP0;
  WlzDVertex2	*tDVP0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  outDom.core = NULL;
  nullDom.core = NULL;
  if(inDom.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(inObjType)
    {
      case WLZ_2D_DOMAINOBJ:
        outDom.i = WlzNewIDomain(inDom.i->type, inDom.i, &errNum);
	if(errNum == WLZ_ERR_NONE )
	{   
	  outDom.i->kol1 += xShift;
	  outDom.i->lastkl += xShift;
	  outDom.i->line1 += yShift;
	  outDom.i->lastln += yShift;
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	outDom.p = WlzMakePlaneDomain(inDom.p->type,
				      inDom.p->plane1 + zShift,
				      inDom.p->lastpl + zShift,
				      inDom.p->line1  + yShift,
				      inDom.p->lastln  + yShift,
				      inDom.p->kol1  + xShift,
				      inDom.p->lastkl  + xShift,
				      &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  outDom.p->voxel_size[0] = inDom.p->voxel_size[0];
	  outDom.p->voxel_size[1] = inDom.p->voxel_size[1];
	  outDom.p->voxel_size[2] = inDom.p->voxel_size[2];
	  idx0 = inDom.p->plane1;
	  domP0 = inDom.p->domains;
	  domP1 = outDom.p->domains;
	  while((idx0 <= inDom.p->lastpl) && (errNum == WLZ_ERR_NONE))
	  {
	    *domP1++ = (domP0->core)?
	    	       WlzAssignDomain(
	    	       WlzShiftDomain(WLZ_2D_DOMAINOBJ, *domP0,
		       		      xShift, yShift, zShift, &errNum), NULL):
		       nullDom;
	    ++domP0;
	    ++idx0;
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    while(--idx0 >= inDom.p->plane1)
	    {
	      if(domP1->core)
	      {
	        (void )WlzFreeDomain(*domP1);
	      }
	      --domP1;
	    }
	  }
	}
	if((errNum != WLZ_ERR_NONE) && outDom.core)
	{
	  (void )WlzFreeDomain(outDom);
	  outDom.core = NULL;
	}
	break;
      case WLZ_TRANS_OBJ: /* FALLTHROUGH */
      case WLZ_AFFINE_TRANS:
        tDom0.t = WlzAffineTransformFromPrim(inDom.t->type, (double )xShift,
					     (double )yShift, (double )zShift,
					     1.0, 0.0, 0.0, 0.0,
					     0.0, 0.0, 0, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  outDom.t = WlzAffineTransformProduct(inDom.t, tDom0.t, &errNum);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    (void )WlzFreeAffineTransform(tDom0.t);
	    outDom.t = NULL;
	  }
	}
	break;
      case WLZ_2D_POLYGON:
        outDom.poly = WlzMakePolyDmn(inDom.poly->type, inDom.poly->vtx,
				     inDom.poly->nvertices,
				     inDom.poly->maxvertices,
				     1, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
          switch(outDom.poly->type)
	  {
	    case WLZ_POLYGON_INT:
	      tIVP0 = outDom.poly->vtx;
	      idx0 = outDom.poly->nvertices;
	      while(idx0-- > 0)
	      {
	        tIVP0->vtX += xShift;
		tIVP0->vtY += yShift;
		++tIVP0;
	      }
	      break;
	    case WLZ_POLYGON_FLOAT:
	      tFVP0 = (WlzFVertex2 *)(outDom.poly->vtx);
	      idx0 = outDom.poly->nvertices;
	      while(idx0-- > 0)
	      {
	        tFVP0->vtX += xShift;
		tFVP0->vtY += yShift;
		++tFVP0;
	      }
	      break;
	    case WLZ_POLYGON_DOUBLE:
	      tDVP0 = (WlzDVertex2 *)(outDom.poly->vtx);
	      idx0 = outDom.poly->nvertices;
	      while(idx0-- > 0)
	      {
	        tDVP0->vtX += xShift;
		tDVP0->vtY += yShift;
		++tDVP0;
	      }
	      break;
	  }
	}
        break;
      case WLZ_BOUNDLIST:
	if(inDom.b->poly)
	{
	  tDom0.poly = inDom.b->poly;
	  tDom1 = WlzShiftDomain(WLZ_2D_POLYGON, tDom0,
				 xShift, yShift, zShift, &errNum);
	}
	else
	{
	  tDom1.poly = NULL;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  outDom.b = WlzMakeBoundList(inDom.b->type, inDom.b->wrap, tDom1.poly,
	  			      &errNum);
	  if((errNum != WLZ_ERR_NONE) && tDom1.poly)
	  {
	    (void )WlzFreePolyDmn(tDom1.poly);
	  }
	}

	if((errNum == WLZ_ERR_NONE) && inDom.b->up)
	{
	  tDom0.b = inDom.b->up;
	  tDom1 = WlzAssignDomain(WlzShiftDomain(WLZ_BOUNDLIST, tDom0,
	  					 xShift, yShift, zShift,
						 &errNum), NULL);
	  outDom.b->up = tDom1.b;
	}
	if((errNum == WLZ_ERR_NONE) && inDom.b->next)
	{
	  tDom0.b = inDom.b->next;
	  tDom1 = WlzAssignDomain(WlzShiftDomain(WLZ_BOUNDLIST, tDom0,
						 xShift, yShift, zShift,
	  			  		 &errNum), NULL);
	  outDom.b->next = tDom1.b;
	}
	if((errNum == WLZ_ERR_NONE) && inDom.b->down)
	{
	  tDom0.b = inDom.b->down;
	  tDom1 = WlzAssignDomain(WlzShiftDomain(WLZ_BOUNDLIST, tDom0,
						 xShift, yShift, zShift,
	  			  		 &errNum), NULL);
	  outDom.b->down = tDom1.b;
	}
	if((errNum != WLZ_ERR_NONE) && outDom.core)
	{
	  (void )WlzFreeBoundList(outDom.b);
	  outDom.core = NULL;
	}
        break;
      case WLZ_HISTOGRAM:
      case WLZ_3D_WARP_TRANS:
      case WLZ_CONV_HULL:
      case WLZ_3D_POLYGON:
      case WLZ_RECTANGLE:
      case WLZ_VECTOR_INT:
      case WLZ_VECTOR_FLOAT:
      case WLZ_POINT_INT:
      case WLZ_POINT_FLOAT:
      case WLZ_CONVOLVE_INT:
      case WLZ_CONVOLVE_FLOAT:
      case WLZ_WARP_TRANS:
      case WLZ_FMATCHOBJ:
      case WLZ_TEXT:
      case WLZ_COMPOUND_ARR_1:
      case WLZ_COMPOUND_ARR_2:
      case WLZ_COMPOUND_LIST_1:
      case WLZ_COMPOUND_LIST_2:
        break;
      case WLZ_PROPERTY_OBJ: /* FALLTHROUGH */
      case WLZ_EMPTY_OBJ:
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(outDom);
}

/************************************************************************
* Function:	WlzShiftValues						*
* Returns:	WlzValues:		Shifted values, NULL on error.	*
* Purpose:	Shifts the given objects values.			*
* Global refs:	-							*
* Parameters:	WlzObjectType inObjType: Type of given domain's parent	*
*					object.				*
*		WlzValues inVal:	Values to be shifted.		*
*		WlzDomain inDom:	Domain over which values are	*
*					defined (parent object's 	*
*					domain).			*
*		int xShift:		Column shift.			*
*		int yShift:		Line shift.			*
*		int zShift:		Plane shift (only used for 3D	*
*					objects).			*
*		WlzErrorNum *dstErr:	Destination error pointer, may	*
*					be NULL.			*
************************************************************************/
WlzValues	 WlzShiftValues(WlzObjectType inObjType, WlzValues inVal,
			       WlzDomain inDom,
			       int xShift, int yShift, int zShift,
			       WlzErrorNum *dstErr)
{
  int		idx0,
		idx1,
		lnCount,
		vIvCount;
  WlzValues	outVal,
  		nullVal;
  WlzObjectType vTabType;
  WlzValues	*valP0,
  		*valP1;
  WlzDomain	*domP0;
  WlzValueIntervalLine *inVIvLn,
  		*outVIvLn;
  WlzValueLine	*inVLn,
  		*outVLn;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  outVal.core = NULL;
  if(inVal.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(inObjType)
    {
      case WLZ_2D_DOMAINOBJ:
 	vTabType = WlzGreyTableTypeToTableType(inVal.core->type, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  switch(vTabType)
	  {
	    case WLZ_GREY_TAB_RAGR:
	      outVal.v = WlzMakeValueTb(inVal.v->type,
	      				inVal.v->line1 + yShift,
					inVal.v->lastln + yShift,
					inVal.v->kol1 + xShift,
					inVal.v->bckgrnd,
					NULL, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		outVal.v->width = inVal.v->width;
	        outVal.v->original_table = WlzAssignValues(inVal, NULL);
		lnCount = inVal.v->lastln - inVal.v->line1 + 1;
		idx0 = lnCount;
		inVLn = inVal.v->vtblines;
		outVLn = outVal.v->vtblines;
		while(idx0-- > 0)
		{
		  *outVLn++ = *inVLn++;
		}
	      }
	      break;
	    case WLZ_GREY_TAB_RECT:
	      outVal.r = WlzMakeRectValueTb(inVal.r->type,
	      				    inVal.r->line1 + yShift,
					    inVal.r->lastln + yShift,
					    inVal.r->kol1 + xShift,
					    inVal.r->width,
					    inVal.r->bckgrnd,
					    NULL,
					    &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        outVal.r->original_table = WlzAssignValues(inVal, NULL);
		outVal.r->values = inVal.r->values;
	      }
	      break;
	    case WLZ_GREY_TAB_INTL:
	      vIvCount = 0;
	      inVIvLn = inVal.i->vil;
	      lnCount = inVal.i->lastln - inVal.i->line1 + 1;
	      idx0 = lnCount;
	      while(idx0-- > 0)
	      {
		vIvCount += inVIvLn->nintvs;
		++inVIvLn;
	      }
	      if((outVal.i = (WlzIntervalValues *)
	      		     AlcCalloc(sizeof(WlzIntervalValues) +
			 	       (sizeof(WlzValueIntervalLine) *
				        lnCount) +
			 	       (sizeof(WlzValueLine) * vIvCount),
				       1)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        outVal.i->type = inVal.i->type;
		outVal.i->linkcount = 0;
		outVal.i->freeptr = NULL;
		outVal.i->original_table = WlzAssignValues(inVal, NULL);
		outVal.i->line1 = inVal.i->line1 + yShift;
		outVal.i->lastln = inVal.i->lastln + yShift;
		outVal.i->kol1 = inVal.i->kol1 + xShift;
		outVal.i->width = inVal.i->width;
		outVIvLn = (WlzValueIntervalLine *)(outVal.i + 1);
		outVLn = (WlzValueLine *)(outVIvLn + lnCount);
		outVal.i->vil = outVIvLn;
		inVIvLn = inVal.i->vil;
		idx0 = lnCount;
		while(idx0-- > 0)
		{
		  idx1 = outVIvLn->nintvs = inVIvLn->nintvs;
		  inVLn = inVIvLn->vtbint;
		  outVIvLn->vtbint = outVLn;
		  while(idx1-- > 0)
		  {
		    *outVLn++ = *inVLn++;
		  }
		  ++inVIvLn;
		  ++outVIvLn;
		}
	      }
	      break;
	    default:
      	      errNum = WLZ_ERR_OBJECT_TYPE;
	      break;
	  }
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	if(inVal.vox->type != WLZ_VOXELVALUETABLE_GREY)
	{
	  errNum = WLZ_ERR_VALUES_TYPE;
	}
	else if(inDom.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(inDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if((outVal.vox = WlzMakeVoxelValueTb(inVal.vox->type,
					       inVal.vox->plane1,
					       inVal.vox->lastpl,
					       inVal.vox->bckgrnd, NULL,
					       &errNum)) != NULL)
	  {
	    idx0 = inVal.vox->plane1;
	    domP0 = inDom.p->domains;
	    valP0 = inVal.vox->values;
	    valP1 = outVal.vox->values;
	    while((idx0 <= inVal.vox->lastpl) && (errNum == WLZ_ERR_NONE))
	    {
	      *valP1++ = (valP0->core)?
			 WlzAssignValues(
			 WlzCopyValues(WLZ_2D_DOMAINOBJ, *valP0, *domP0,
			 	       &errNum), NULL):
			 nullVal;
	      ++valP0;
	      ++domP0;
	      ++idx0;
	    }
	    if(errNum != WLZ_ERR_NONE)
	    {
	      while(--idx0 >= inVal.vox->plane1)
	      {
		if(valP1->core)
		{
		  (void )WlzFreeValues(*valP1);
		}
		--valP1;
	      }
	    }
	  }
	}
	break;
      case WLZ_TRANS_OBJ:
        outVal.obj = inVal.obj;
	break;
      case WLZ_3D_WARP_TRANS:
      case WLZ_CONV_HULL:
      case WLZ_3D_POLYGON:
      case WLZ_RECTANGLE:
      case WLZ_VECTOR_INT:
      case WLZ_VECTOR_FLOAT:
      case WLZ_POINT_INT:
      case WLZ_POINT_FLOAT:
      case WLZ_CONVOLVE_INT:
      case WLZ_CONVOLVE_FLOAT:
      case WLZ_WARP_TRANS:
      case WLZ_FMATCHOBJ:
      case WLZ_TEXT:
      case WLZ_COMPOUND_ARR_1:
      case WLZ_COMPOUND_ARR_2:
      case WLZ_COMPOUND_LIST_1:
      case WLZ_COMPOUND_LIST_2:
      case WLZ_EMPTY_OBJ:
        break;
      case WLZ_AFFINE_TRANS: /* FALLTHROUGH */
      case WLZ_HISTOGRAM:    /* FALLTHROUGH */
      case WLZ_PROPERTY_OBJ: /* FALLTHROUGH */
      case WLZ_2D_POLYGON:   /* FALLTHROUGH */
      case WLZ_BOUNDLIST:
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(outVal);
}
