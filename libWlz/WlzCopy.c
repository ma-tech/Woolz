#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzCopy.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions to make 'deep' copies of Woolz objects.
* \ingroup	WlzAllocation
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <Wlz.h>

WlzDomain	WlzCopyDomain(WlzObjectType, WlzDomain, WlzErrorNum *);
WlzValues	WlzCopyValues(WlzObjectType, WlzValues, WlzDomain,
			       WlzErrorNum *);
WlzSimpleProperty *WlzCopySimpleProperty(WlzSimpleProperty *,
					 	WlzErrorNum *);

/*!
* \return	Copy of given object.
* \ingroup	WlzAllocation
* \brief	Copies a Woolz object together with it's domain, values
*		and properties and then returns the copy.
* \param	inObj			The given object.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzObject	*WlzCopyObject(WlzObject *inObj, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*outObj = NULL;
  WlzSimpleProperty *pLst = NULL;

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
      case WLZ_AFFINE_TRANS:
      case WLZ_HISTOGRAM:
      case WLZ_PROPERTY_OBJ:
      case WLZ_2D_POLYGON:
      case WLZ_BOUNDLIST:
      case WLZ_CONTOUR:
	dom = WlzCopyDomain(inObj->type, inObj->domain, &errNum);
	if(inObj->values.core)
	{
	  val = WlzCopyValues(inObj->type, inObj->values, inObj->domain,
	  		      &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && inObj->plist)
	{
	  pLst = WlzCopySimpleProperty(inObj->plist, &errNum);
        }
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzMakeMain(inObj->type, dom, val, pLst, NULL, &errNum);
	}
	break;
      case WLZ_EMPTY_OBJ:
        outObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_CONV_HULL:
      case WLZ_3D_WARP_TRANS:
      case WLZ_3D_POLYGON:
      case WLZ_RECTANGLE:
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
    if(pLst)
    {
      WlzFreeProperty(pLst);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(outObj);
}

/*!
* \return	Copied domain, NULL on error.
* \ingroup      WlzAllocation
* \brief	Copies the given objects domain.
* \param	inObjType		Type of given domain's parent
*					object.
* \param	inDom			Domain to be copied.
* \param	dstErr			Destination error pointer, may
*					be NULL.
*/
WlzDomain	 WlzCopyDomain(WlzObjectType inObjType, WlzDomain inDom,
			       WlzErrorNum *dstErr)
{
  int		idx0;
  WlzDomain	tDom0,
  		tDom1,
  		outDom,
  		nullDom;
  WlzDomain	*domP0,
  		*domP1;
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
	break;
      case WLZ_3D_DOMAINOBJ:
	outDom.p = WlzMakePlaneDomain(inDom.p->type,
				      inDom.p->plane1, inDom.p->lastpl,
				      inDom.p->line1, inDom.p->lastln,
				      inDom.p->kol1, inDom.p->lastkl,
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
		       WlzCopyDomain(WLZ_2D_DOMAINOBJ, *domP0, &errNum), NULL):
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
        outDom.t = WlzAffineTransformCopy(inDom.t, &errNum);
	break;
      case WLZ_HISTOGRAM:
        outDom.hist = WlzMakeHistogramDomain(inDom.hist->type,
					     inDom.hist->maxBins, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
	  outDom.hist->nBins = inDom.hist->nBins;
	  outDom.hist->origin = inDom.hist->origin;
	  outDom.hist->binSize = inDom.hist->binSize;
	  WlzValueCopyGreyToGrey(outDom.hist->binValues, 0,
			       (outDom.hist->type == WLZ_HISTOGRAMDOMAIN_INT)?
			       (WLZ_GREY_INT): (WLZ_GREY_DOUBLE),
			       inDom.hist->binValues, 0,
			       (inDom.hist->type == WLZ_HISTOGRAMDOMAIN_INT)?
			       (WLZ_GREY_INT): (WLZ_GREY_DOUBLE),
			       inDom.hist->nBins);
	}
	break;
      case WLZ_2D_POLYGON:
        outDom.poly = WlzMakePolyDmn(inDom.poly->type, inDom.poly->vtx,
				     inDom.poly->nvertices,
				     inDom.poly->maxvertices,
				     1, &errNum);
        break;
      case WLZ_BOUNDLIST:
	tDom0.poly = (inDom.b->poly)?
		     WlzMakePolyDmn(inDom.b->poly->type, inDom.b->poly->vtx,
		     		    inDom.b->poly->nvertices,
				    inDom.b->poly->maxvertices,
				    1, &errNum): NULL;
	if(errNum == WLZ_ERR_NONE)
	{
	  outDom.b = WlzMakeBoundList(inDom.b->type, inDom.b->wrap, tDom0.poly,
	  			      &errNum);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    (void )WlzFreePolyDmn(tDom0.poly);
	  }
	}
	if((errNum == WLZ_ERR_NONE) && inDom.b->next)
	{
	  tDom0.b = inDom.b->next;
	  tDom1 = WlzAssignDomain(WlzCopyDomain(WLZ_BOUNDLIST, tDom0,
	  			  	        &errNum), NULL);
	  outDom.b->next = tDom1.b;
	}
	if((errNum == WLZ_ERR_NONE) && inDom.b->down)
	{
	  tDom0.b = inDom.b->down;
	  tDom1 = WlzAssignDomain(WlzCopyDomain(WLZ_BOUNDLIST, tDom0,
	  			  	        &errNum), NULL);
	  outDom.b->down = tDom1.b;
	}
	if((errNum != WLZ_ERR_NONE) && outDom.core)
	{
	  (void )WlzFreeBoundList(outDom.b);
	  outDom.core = NULL;
	}
        break;
      case WLZ_CONTOUR:
	if((outDom.ctr = WlzMakeContour(&errNum)) != NULL)
	{
	  outDom.ctr->model = WlzAssignGMModel(
	  		      WlzGMModelCopy(inDom.ctr->model, &errNum), NULL);
	}
        break;
      case WLZ_3D_WARP_TRANS:
      case WLZ_CONV_HULL:
      case WLZ_3D_POLYGON:
      case WLZ_RECTANGLE:
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

/*!
* \return	Copied values, NULL on error.
* \ingroup      WlzAllocation
* \brief	Copies the given values.
* \param	inObjType		Type of given values parent
*					object.
* \param	inVal			Values to be copied.
* \param	inDom			Domain over which values are
*					defined (parent object's domain).
* \param	dstErr			Destination error pointer, may
*					be NULL.
*/
WlzValues	 WlzCopyValues(WlzObjectType inObjType, WlzValues inVal,
			       WlzDomain inDom, WlzErrorNum *dstErr)
{
  int		idx0;
  WlzObject	*tObj0,
  		*tObj1;
  WlzValues	outVal,
  		nullVal;
  WlzValues	*valP0,
  		*valP1;
  WlzDomain	*domP0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nullVal.core = outVal.core = NULL;
  if(inVal.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(inObjType)
    {
      case WLZ_2D_DOMAINOBJ:
	if( (tObj0 = WlzMakeMain(inObjType, inDom, inVal,
				 NULL, NULL, &errNum)) != NULL ){
	  tObj0 = WlzAssignObject(tObj0, NULL);
	  if( (tObj1 = WlzNewGrey(tObj0, &errNum)) != NULL ){
	    outVal = WlzAssignValues(tObj1->values, NULL);
	    (void )WlzFreeObj(tObj1);
	  }
	  (void )WlzFreeObj(tObj0);
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
        outVal.obj = WlzCopyObject(inVal.obj, &errNum);
	break;
      case WLZ_3D_WARP_TRANS:
      case WLZ_CONV_HULL:
      case WLZ_3D_POLYGON:
      case WLZ_CONTOUR:
      case WLZ_RECTANGLE:
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

/*!
* \return	Copied property list.
* \ingroup	WlzAllocation
* \brief	Copies the given simple property list.
* \param	inPLst			Given property list.
* \param	dstErr			Destination error pointer, may
*					be NULL.
*/
WlzSimpleProperty *WlzCopySimpleProperty(WlzSimpleProperty *inPLst,
					 WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzSimpleProperty *outPLst = NULL;

  if(inPLst == NULL)
  {
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
  else if(inPLst->type != WLZ_PROPERTY_SIMPLE)
  {
    errNum = WLZ_ERR_PROPERTY_TYPE;
  }
  else
  {
    if((outPLst = (WlzSimpleProperty *)
	       AlcCalloc(1, sizeof(WlzSimpleProperty))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else if (inPLst->prop && (inPLst->size > 0))
    {
      if((outPLst->prop = AlcMalloc(inPLst->size)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )memcpy(outPLst->prop, inPLst->prop, inPLst->size);
      outPLst->type = inPLst->type;
      outPLst->size = inPLst->size;
    }
    else
    {
      if(outPLst)
      {
        AlcFree(outPLst);
	outPLst = NULL;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(outPLst);
}
