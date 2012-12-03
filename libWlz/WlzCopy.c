#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCopy_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCopy.c
* \author       Bill Hill
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
* \brief	Functions to make 'deep' copies of objects.
* \ingroup	WlzAllocation
*/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <Wlz.h>

static WlzSimpleProperty  	*WlzCopySimpleProperty(
				  WlzSimpleProperty *PList,
				  WlzErrorNum *dstErr);
static WlzLUTValues 		*WlzCopyLUTValues(
				  WlzLUTValues *inVal,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzCopyObjectGreyValues2D(
				  WlzObject *dObj,
				  WlzObject *sObj);
static WlzErrorNum 		WlzCopyObjectGreyValues3D(
				  WlzObject *dObj,
				  WlzObject *sObj);
static WlzErrorNum 		WlzCopyObjectGreyValuesAny2D(
				  WlzObject *dObj,
				  WlzObject *sObj);
static WlzErrorNum 		WlzCopyObjectGreyValuesAny3D(
				  WlzObject *dObj,
				  WlzObject *sObj);
static WlzErrorNum 		WlzCopyObjectGreyValuesScan2D(
				  WlzObject *dObj,
				  WlzObject *sObj);
static WlzErrorNum 		WlzCopyObjectGreyValuesScan3D(
				  WlzObject *dObj,
				  WlzObject *sObj);
static WlzErrorNum 		WlzCopyObjectGreyValuesGVWSp2D(
				  WlzObject *dObj,
				  WlzGreyValueWSpace *dGVWSp,
				  WlzObject *sObj,
				  WlzGreyValueWSpace *sGVWSp,
				  int pln);

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
  WlzPropertyList *pLst = NULL;

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
      case WLZ_2D_DOMAINOBJ:   /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:   /* FALLTHROUGH */
      case WLZ_3D_VIEW_STRUCT: /* FALLTHROUGH */
      case WLZ_TRANS_OBJ:      /* FALLTHROUGH */
      case WLZ_AFFINE_TRANS:   /* FALLTHROUGH */
      case WLZ_HISTOGRAM:      /* FALLTHROUGH */
      case WLZ_PROPERTY_OBJ:   /* FALLTHROUGH */
      case WLZ_2D_POLYGON:     /* FALLTHROUGH */
      case WLZ_BOUNDLIST:      /* FALLTHROUGH */
      case WLZ_CONTOUR:        /* FALLTHROUGH */
      case WLZ_LUT:
	dom = WlzCopyDomain(inObj->type, inObj->domain, &errNum);
	if(inObj->values.core)
	{
	  val = WlzCopyValues(inObj->type, inObj->values, inObj->domain,
	  		      &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && inObj->plist)
	{
	  pLst = WlzCopyPropertyList(inObj->plist, &errNum);
        }
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzMakeMain(inObj->type, dom, val, pLst, NULL, &errNum);
	}
	break;
      case WLZ_CMESH_2D:	/* FALLTHROUGH */
      case WLZ_CMESH_2D5:	/* FALLTHROUGH */
      case WLZ_CMESH_3D:
	if(inObj->values.core != NULL)
	{
	  if(inObj->values.core->type != WLZ_INDEXED_VALUES)
	  {
	    errNum = WLZ_ERR_VALUES_TYPE;
	  }
	  else
	  {
	    WlzIndexedValues *ixv;

	    ixv = inObj->values.x;
	    val.x = WlzMakeIndexedValues(inObj, ixv->rank, ixv->dim,
	                                 ixv->vType, ixv->attach, &errNum);
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  size_t datSz;
	  WlzCMeshP meshP;
	  AlcVector *newVec = NULL;

	  datSz = WlzIndexedValueSize(val.x, NULL);
	  /* Although the 2D conforming mesh union member is used all that
	   * really matters is that the pointer is passed. */
	  meshP.m2 = inObj->domain.cm2;
	  meshP =  WlzCMeshCopy(meshP, 1, datSz,
				(val.core)? &newVec: NULL,
				(val.core)? (inObj->values.x->values): NULL,
				&errNum);
	  dom.cm2 = meshP.m2;
	  if(newVec)
	  {
	    (void )AlcVectorFree(val.x->values);
	    val.x->values = newVec;
	  }
	}
	if((errNum == WLZ_ERR_NONE) && inObj->plist)
	{
	  pLst = WlzCopyPropertyList(inObj->plist, &errNum);
        }
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzMakeMain(inObj->type, dom, val, pLst, NULL, &errNum);
	}
	break;
      case WLZ_EMPTY_OBJ:
        outObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_CONV_HULL:       /* FALLTHROUGH */
      case WLZ_3D_WARP_TRANS:   /* FALLTHROUGH */
      case WLZ_3D_POLYGON:      /* FALLTHROUGH */
      case WLZ_RECTANGLE:       /* FALLTHROUGH */
      case WLZ_CONVOLVE_INT:    /* FALLTHROUGH */
      case WLZ_CONVOLVE_FLOAT:  /* FALLTHROUGH */
      case WLZ_WARP_TRANS:      /* FALLTHROUGH */
      case WLZ_FMATCHOBJ:       /* FALLTHROUGH */
      case WLZ_TEXT:            /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_1:  /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:  /* FALLTHROUGH */
      case WLZ_COMPOUND_LIST_1: /* FALLTHROUGH */
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
      WlzFreePropertyList(pLst);
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
      case WLZ_3D_VIEW_STRUCT:
	outDom.vs3d = WlzMake3DViewStructCopy(inDom.vs3d, &errNum);
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
        outDom.poly = WlzMakePolygonDomain(inDom.poly->type,
				     inDom.poly->nvertices,
				     inDom.poly->vtx,
				     inDom.poly->maxvertices,
				     1, &errNum);
        break;
      case WLZ_BOUNDLIST:
	tDom0.poly = (inDom.b->poly)?
		     WlzMakePolygonDomain(inDom.b->poly->type,
		     		    inDom.b->poly->nvertices,
		     		    inDom.b->poly->vtx,
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
      case WLZ_CMESH_2D:
        outDom.cm2 = WlzCMeshCopy2D(inDom.cm2, 1, 0, NULL, NULL, &errNum);
	break;
      case WLZ_CMESH_2D5:
        outDom.cm2d5 = WlzCMeshCopy2D5(inDom.cm2d5, 1, 0, NULL, NULL, &errNum);
	break;
      case WLZ_CMESH_3D:
        outDom.cm3 = WlzCMeshCopy3D(inDom.cm3, 0, 0, NULL, NULL, &errNum);
	break;
      case WLZ_LUT:
	outDom.lut = WlzMakeLUTDomain(inDom.lut->bin1, inDom.lut->lastbin,
	                              &errNum);
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
	if((tObj0 = WlzMakeMain(inObjType, inDom, inVal,
				 NULL, NULL, &errNum)) != NULL)
        {
	  tObj0 = WlzAssignObject(tObj0, NULL);
	  if((tObj1 = WlzNewGrey(tObj0, &errNum)) != NULL)
	  {
	    outVal = WlzAssignValues(tObj1->values, NULL);
	    (void )WlzFreeObj(tObj1);
	  }
	  (void )WlzFreeObj(tObj0);
	  if(outVal.core && (outVal.core->linkcount == 1))
	  {
	    /* The linkcount is 1 but returned objects, values, etc should
	     * have a linkcount of 0 unless they are used more than once.
	     */
	    outVal.core->linkcount = 0;
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
        outVal.obj = WlzCopyObject(inVal.obj, &errNum);
	break;
      case WLZ_LUT:
        outVal.lut = WlzCopyLUTValues(inVal.lut, &errNum);
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
* \brief	Copies the given property list.
*
*		A new property list is created with a zero link count
*		and a linked list for the properties. Each property of
*		given list is then copied and assigned to the new list.
*		The order of the properties in the new list is the same
*		as that in the given list, although a property that
*		occurs more than once in the given list will have many
*		copies created.
* \param	gList			Given property list.
* \param	dstErr			Destination error pointer, may
*					be NULL.
*/
WlzPropertyList	*WlzCopyPropertyList(WlzPropertyList *gList,
					WlzErrorNum *dstErr)
{
  WlzPropertyList *nList = NULL;
  AlcDLPItem	*gItem;
  WlzProperty	gProp,
  		nProp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nList = WlzMakePropertyList(NULL)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else if((gItem = gList->list->head) != NULL)
  {
    do
    {
      if(gItem->entry == NULL)
      {
        errNum = WLZ_ERR_PROPERTY_NULL;
      }
      else
      {
	gProp.core = (WlzCoreProperty *)(gItem->entry);
        switch(gProp.core->type)
	{
	  case WLZ_PROPERTY_SIMPLE:
	    nProp.simple = WlzCopySimpleProperty(gProp.simple, &errNum);
	    break;
	  case WLZ_PROPERTY_EMAP:
	    nProp.emap = WlzMakeEMAPProperty(gProp.emap->emapType,
					     gProp.emap->modelUID,
					     gProp.emap->anatomyUID,
					     gProp.emap->targetUID,
					     gProp.emap->targetVersion,
					     gProp.emap->stage,
					     gProp.emap->subStage,
					     gProp.emap->modelName,
					     gProp.emap->version,
					     gProp.emap->fileName,
					     gProp.emap->comment,
					     &errNum);
	    break;
	  case WLZ_PROPERTY_NAME:
	    nProp.name = WlzMakeNameProperty(gProp.name->name, &errNum);
	    break;
	  case WLZ_PROPERTY_GREY:
	    nProp.greyV = WlzMakeGreyProperty(gProp.greyV->name,
	    				       gProp.greyV->value, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_PROPERTY_TYPE;
	    break;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzAssignProperty(nProp, NULL);
	if(AlcDLPListEntryAppend(nList->list, NULL, (void *)(nProp.core),
	 			 WlzFreePropertyListEntry) != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	gItem = gItem->next;
      }
    } while((errNum == WLZ_ERR_NONE) && (gItem != gList->list->head));
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreePropertyList(nList);
    nList = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nList);
}

/*!
* \return	Copied simple property.
* \ingroup	WlzAllocation
* \brief	Copies the given simple property.
* \param	inPLst			Given property list.
* \param	dstErr			Destination error pointer, may
*					be NULL.
*/
static WlzSimpleProperty *WlzCopySimpleProperty(WlzSimpleProperty *inPLst,
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

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the source object to the destination
* 		object, overwriting the values of the destination object.
* \param	dObj			Destination object.
* \param	sObj			Source object.
*/
WlzErrorNum	WlzCopyObjectGreyValues(WlzObject *dObj, WlzObject *sObj)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((dObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(dObj->type != sObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((dObj->domain.core == NULL) || (sObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dObj->values.core == NULL) || (sObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(dObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	errNum = WlzCopyObjectGreyValues2D(dObj, sObj);
        break;
      case WLZ_3D_DOMAINOBJ:
	errNum = WlzCopyObjectGreyValues3D(dObj, sObj);
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the 2D source object to the 2D
* 		destination object, overwriting the values of the
* 		destination object. Both the destination and source objects
* 		are known to be valid 2D domain objects with valid non-NULL
* 		values tables.
* \param	dObj			Destination object.
* \param	sObj			Source object.
*/
static WlzErrorNum WlzCopyObjectGreyValues2D(WlzObject *dObj, WlzObject *sObj)
{
  WlzObjectType dGTType = WLZ_NULL,
  		sGTType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dGTType = WlzGreyTableTypeToTableType(dObj->values.core->type, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    sGTType = WlzGreyTableTypeToTableType(sObj->values.core->type, &errNum);
  }
  if((dGTType == WLZ_GREY_TAB_TILED) || (sGTType == WLZ_GREY_TAB_TILED))
  {
    errNum = WlzCopyObjectGreyValuesAny2D(dObj, sObj);
  }
  else
  {
    errNum = WlzCopyObjectGreyValuesScan2D(dObj, sObj);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the 3D source object to the 3D
* 		destination object, overwriting the values of the
* 		destination object. Both the destination and source objects
* 		are known to be valid 3D domain objects with valid non-NULL
* 		values tables.
* \param	dObj			Destination object.
* \param	sObj			Source object.
*/
static WlzErrorNum WlzCopyObjectGreyValues3D(WlzObject *dObj, WlzObject *sObj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dObj->values.core->type == WLZ_VOXELVALUETABLE_GREY) &&
     (sObj->values.core->type == WLZ_VOXELVALUETABLE_GREY))
  {
    errNum = WlzCopyObjectGreyValuesScan3D(dObj, sObj);
  }
  else
  {
    errNum = WlzCopyObjectGreyValuesAny3D(dObj, sObj);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the 2D source object to the 2D
* 		destination object, overwriting the values of the
* 		destination object. Both the destination and source objects
* 		are known to be valid 2D domain objects with valid values
* 		tables of any type suitable for WlzGreyValueGet().
* \param	dObj			Destination object.
* \param	sObj			Source object.
* 					
*/
static WlzErrorNum WlzCopyObjectGreyValuesAny2D(WlzObject *dObj,
					WlzObject *sObj)
{
  WlzGreyValueWSpace *dGVWSp = NULL,
  		*sGVWSp = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  dGVWSp = WlzGreyValueMakeWSp(dObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    sGVWSp = WlzGreyValueMakeWSp(sObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCopyObjectGreyValuesGVWSp2D(dObj, dGVWSp, sObj, sGVWSp, 0);
  }
  WlzGreyValueFreeWSp(dGVWSp);
  WlzGreyValueFreeWSp(sGVWSp);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the 3D source object to the 3D
* 		destination object, overwriting the values of the
* 		destination object. Both the destination and source objects
* 		are known to be valid 3D domain objects with valid values
* 		tables of any type suitable for WlzGreyValueGet().
* \param	dObj			Destination object.
* \param	sObj			Source object.
*/
static WlzErrorNum WlzCopyObjectGreyValuesAny3D(WlzObject *dObj,
					WlzObject *sObj)
{
  int		idP,
  		plMin,
		plMax;
  WlzPlaneDomain *dPDom,
  		 *sPDom;
  WlzValues	nullVal;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  nullVal.core = NULL;
  dPDom = dObj->domain.p;
  sPDom = sObj->domain.p;
  plMin = ALG_MAX(dPDom->plane1, sPDom->plane1);
  plMax = ALG_MIN(dPDom->lastpl, sPDom->lastpl);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(idP = plMin; idP <= plMax; ++idP)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain	*dDom2D,
		*sDom2D;
      WlzErrorNum errNum2D = WLZ_ERR_NONE;

      if(((dDom2D = dPDom->domains + idP - dPDom->plane1) != NULL) &&
	 ((sDom2D = sPDom->domains + idP - sPDom->plane1) != NULL))
      {
	WlzObject *dObj2D = NULL,
		  *sObj2D = NULL;
	WlzGreyValueWSpace *dGVWSp = NULL,
  		*sGVWSp = NULL;
	
	if(((dObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dDom2D, nullVal,
				 NULL, NULL, &errNum2D)) != NULL) &&
	   ((sObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *sDom2D, nullVal,
				 NULL, NULL, &errNum2D)) != NULL) &&
	   ((dGVWSp = WlzGreyValueMakeWSp(dObj, &errNum2D)) != NULL) &&
	   ((sGVWSp = WlzGreyValueMakeWSp(sObj, &errNum2D)) != NULL))
	{
	  errNum2D = WlzCopyObjectGreyValuesGVWSp2D(dObj2D, dGVWSp,
						    sObj2D, sGVWSp, idP);
	}
	WlzGreyValueFreeWSp(dGVWSp);
	WlzGreyValueFreeWSp(sGVWSp);
	(void )WlzFreeObj(dObj2D);
	(void )WlzFreeObj(sObj2D);
      }
#ifdef _OPENMP
#pragma omp critical
      {
#endif
	if(errNum2D != WLZ_ERR_NONE)
	{
	  errNum = errNum2D;
	}
#ifdef _OPENMP
      }
#endif
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the 2D source object to the 2D
* 		destination object, overwriting the values of the
* 		destination object. Both the destination and source objects
* 		are known to be valid 2D domain objects with values for which
* 		the WlzInitGreyScan() and WlzNextGreyInterval() access
* 		methods work.
* \param	dObj			Destination object.
* \param	sObj			Source object.
*/
static WlzErrorNum WlzCopyObjectGreyValuesScan2D(WlzObject *dObj,
					WlzObject *sObj)
{
  int		lnMin,
		lnMax;
  WlzIntervalWSpace dIWSp,
  		sIWSp;
  WlzGreyWSpace dGWSp,
  		sGWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  lnMin = ALG_MAX(dObj->domain.i->line1, sObj->domain.i->line1);
  lnMax = ALG_MIN(dObj->domain.i->lastln, sObj->domain.i->lastln);
  if(lnMin <= lnMax)
  {
    if((errNum = WlzInitGreyScan(dObj, &dIWSp, &dGWSp) == WLZ_ERR_NONE) &&
       (errNum = WlzInitGreyScan(sObj, &sIWSp, &sGWSp) == WLZ_ERR_NONE))
    {
      do
      {
	errNum = WlzNextGreyInterval(&dIWSp);
      } while((errNum == WLZ_ERR_NONE) && (dIWSp.linpos < lnMin));
      while((errNum == WLZ_ERR_NONE) && (dIWSp.linpos <= lnMax))
      {
	while((errNum == WLZ_ERR_NONE) && (sIWSp.linpos < dIWSp.linpos))
	{
	  errNum = WlzNextGreyInterval(&sIWSp);
	}
	while((errNum == WLZ_ERR_NONE) && (sIWSp.linpos == dIWSp.linpos))
	{
	  int   klMin,
		klMax;

	  klMin = ALG_MAX(dIWSp.lftpos, sIWSp.lftpos);
	  klMax = ALG_MIN(dIWSp.rgtpos, sIWSp.rgtpos);
	  if(klMin <= klMax)
	  {
	    WlzValueCopyGreyToGrey(dGWSp.u_grintptr, dIWSp.lftpos - klMin,
				   dGWSp.pixeltype,
	                           sGWSp.u_grintptr, sIWSp.lftpos - klMin,
				   sGWSp.pixeltype,
				   klMax - klMin + 1);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzNextGreyInterval(&sIWSp);
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzNextGreyInterval(&dIWSp);
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the 3D source object to the 3D
* 		destination object, overwriting the values of the
* 		destination object. Both the destination and source objects
* 		are known to be valid 3D domain objects with valid voxel
* 		valuses tables.
* \param	dObj			Destination object.
* \param	sObj			Source object.
*/
static WlzErrorNum WlzCopyObjectGreyValuesScan3D(WlzObject *dObj,
					WlzObject *sObj)
{
  int		idP,
  		plMin,
		plMax;
  WlzVoxelValues *dVVal,
  		 *sVVal;
  WlzPlaneDomain *dPDom,
  		 *sPDom;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  dPDom = dObj->domain.p;
  sPDom = sObj->domain.p;
  dVVal = dObj->values.vox;
  sVVal = sObj->values.vox;
  plMin = ALG_MAX(dPDom->plane1, sPDom->plane1);
  plMax = ALG_MIN(dPDom->lastpl, sPDom->lastpl);
  for(idP = plMin; idP <= plMax; ++idP)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain	*dDom2D,
		*sDom2D;
      WlzValues *dVal2D,
      		*sVal2D;
      WlzErrorNum errNum2D = WLZ_ERR_NONE;

      if(((dDom2D = dPDom->domains + idP - dPDom->plane1) != NULL) &&
	 ((sDom2D = sPDom->domains + idP - sPDom->plane1) != NULL) &&
	 ((dVal2D = dVVal->values + idP - dPDom->plane1) != NULL) &&
	 ((sVal2D = sVVal->values + idP - sPDom->plane1) != NULL))
      {
	WlzObject *dObj2D = NULL,
		  *sObj2D = NULL;
	
	if(((dObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dDom2D, *dVal2D,
				 NULL, NULL, &errNum2D)) != NULL) &&
	   ((sObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *sDom2D, *sVal2D,
				 NULL, NULL, &errNum2D)) != NULL))
	{
	  errNum2D = WlzCopyObjectGreyValuesScan2D(dObj2D, sObj2D);
	}
	(void )WlzFreeObj(dObj2D);
	(void )WlzFreeObj(sObj2D);
      }
      if(errNum2D != WLZ_ERR_NONE)
      {
        errNum = errNum2D;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the grey values of the 2 or 3D source object to the
* 		2 or 3D destination object, overwriting the values of the
* 		destination object. Both the destination and source objects
* 		are known to be valid 2 or 3D domain objects with values
* 		tables of any type suitable for WlzGreyValueGet(). The
* 		given objects are used to scan through the object's domains
* 		and the grey value workspaces for their values. Both
* 		workspaces must have been initialised for the objects.
* \param	dObj			Destination object.
* \param	dGVWSp			Destination grey value workspace.
* \param	sObj			Source object.
* \param	sGVWSp			Source grey value workspace.
* \param	pln			The current plane, may be zero
* 					for 2D objects.
* 					
*/
static WlzErrorNum WlzCopyObjectGreyValuesGVWSp2D(
				WlzObject *dObj, WlzGreyValueWSpace *dGVWSp,
				WlzObject *sObj, WlzGreyValueWSpace *sGVWSp,
				int pln)
{
  int		lnMin,
		lnMax;
  WlzIntervalWSpace dIWSp,
  		sIWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  lnMin = ALG_MAX(dObj->domain.i->line1, sObj->domain.i->line1);
  lnMax = ALG_MIN(dObj->domain.i->lastln, sObj->domain.i->lastln);
  if(lnMin <= lnMax)
  {
    if((errNum = WlzInitRasterScan(dObj, &dIWSp,
				   WLZ_RASTERDIR_ILIC) == WLZ_ERR_NONE) &&
       (errNum = WlzInitRasterScan(sObj, &sIWSp,
				   WLZ_RASTERDIR_ILIC) == WLZ_ERR_NONE))
    {
      do
      {
	errNum = WlzNextInterval(&dIWSp);
      } while((errNum == WLZ_ERR_NONE) && (dIWSp.linpos < lnMin));
      while((errNum == WLZ_ERR_NONE) && (dIWSp.linpos <= lnMax))
      {
	while((errNum == WLZ_ERR_NONE) && (sIWSp.linpos < dIWSp.linpos))
	{
	  errNum = WlzNextInterval(&sIWSp);
	}
	while((errNum == WLZ_ERR_NONE) && (sIWSp.linpos == dIWSp.linpos))
	{
	  int		t0,
	  		t1;

	  /* Classify the possible intersections. Here s (source), d
	   * (destination) and o (overlap) are used in comment strings
	   * to represent the six possible interval intersection cases. */
	  if(sIWSp.rgtpos < dIWSp.lftpos)                       /* ssss dddd */
	  {
	    t0 = 1;
	  }
	  else if(sIWSp.lftpos > dIWSp.rgtpos)                  /* dddd ssss */
	  {
	    t0 = 0;
	  }
	  else
	  {
	    int		idK;
	    WlzInterval itv;

            t0 = sIWSp.lftpos >= dIWSp.lftpos;
	    t1 = sIWSp.rgtpos <= dIWSp.rgtpos;
	    if((t0 != 0) && (t1 != 0))                     /* dooood || oooo */
	    {
	      itv.ileft = sIWSp.lftpos;
	      itv.iright = sIWSp.rgtpos;
	    }
	    else if((t0 == 0) && (t1 == 0))                        /* soooos */
	    {
	      itv.ileft = dIWSp.lftpos;
	      itv.iright = dIWSp.rgtpos;
	    }
	    else
	    {
	      if(t0 == 0)                                           /* soood */
	      {
	        itv.ileft = dIWSp.lftpos;
		itv.iright = sIWSp.rgtpos;
	      }
	      else                                                  /* dooos */
	      {
	        itv.ileft = sIWSp.lftpos;
		itv.iright = dIWSp.rgtpos;
	      }
	      t0 = !t0;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      for(idK = itv.ileft; idK <= itv.iright; ++idK)
	      {
		WlzGreyValueGet(dGVWSp, pln, dIWSp.linpos, idK);
		WlzGreyValueGet(sGVWSp, pln, dIWSp.linpos, idK);
		WlzValueCopyGreyToGrey(dGVWSp->gPtr[0], 0, dGVWSp->gType,
				       sGVWSp->gPtr[0], 0, sGVWSp->gType, 1);
	      }
	    }
	  }
	  if(t0 == 0)
	  {
	    errNum = WlzNextInterval(&dIWSp);
	  }
	  else
	  {
	    errNum = WlzNextInterval(&sIWSp);
	  }
	}
	while((errNum == WLZ_ERR_NONE) && (dIWSp.linpos < sIWSp.linpos))
	{
	  errNum = WlzNextInterval(&dIWSp);
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Copies the LUT values of the source object.
* \param	obj			Source LUT object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzLUTValues *WlzCopyLUTValues(WlzLUTValues *inVal, WlzErrorNum *dstErr)
{
  WlzLUTValues	*outVal = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(inVal == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(inVal->type != WLZ_LUT)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    outVal = WlzMakeLUTValues(inVal->vType, inVal->maxVal, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValueCopyGreyToGrey(outVal->val, 0, outVal->vType,
			   inVal->val,  0, inVal->vType, inVal->maxVal);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(outVal);
}
