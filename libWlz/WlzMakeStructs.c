#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMakeStructs_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzMakeStructs.c
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
* \brief	Functions for allocating woolz structures.
* \ingroup	WlzAllocation
*/

#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

/* function:     WlzMakeIntervalDomain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Allocate space for an interval domain structure. If the
type is WLZ_INTERVALDOMAIN_INTVL then allocate space for
the interval line array and set the pointer.
*
* \return       Pointer to allocated interval domain.
* \param    type	Required interval domain type.
* \param    l1	First line
* \param    ll	Last line.
* \param    k1	First column
* \param    kl	Last column
* \param    dstErr	error return values: WLZ_ERR_NONE, WLZ_ERR_PARAM_DATA, WLZ_ERR_MEM_ALLOC
* \par      Source:
*                WlzMakeStructs.c
*/
WlzIntervalDomain *WlzMakeIntervalDomain(WlzObjectType type,
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

/* function:     WlzMakePlaneDomain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Allocate space for a plane domain and the domain
pointers
*
* \return       POinter to the allocated planedomain.
* \param    type	PlaneDomain type.
* \param    p1	First plane
* \param    pl	Last plane.
* \param    l1	First line.
* \param    ll	Last line.
* \param    k1	First column.
* \param    kl	Last column.
* \param    dstErr	error return values: WLZ_ERR_NONE, WLZ_ERR_MEM_ALLOC, WLZ_ERR_PARAM_DATA.
* \par      Source:
*                WlzMakeStructs.c
*/

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

/*!
* \return	New object or NULL on error.
* \ingroup	WlzAllocation
* \brief	Creaes a new object with the domain of the given object
*		and a new values.
*
*		This function conveniently wraps up WlzNewValueTb(),
*		WlzNewValuesVox(), WlzMakeMain() and WlzGreySetValue().
* \param	sObj			Given object.
* \param	tType			Grey table type.
* \param	bgdV			Background value.
* \param	setFG			Set forground value if non zero.
* \param	fgdV			Foreground value to be set.
* \param	dstErr			Destinaition error pointer, may be NULL.
*/
WlzObject	*WlzNewObjectValues(WlzObject *sObj, WlzObjectType tType,
				    WlzPixelV bgdV, int setFG, WlzPixelV fgdV,
				    WlzErrorNum *dstErr)
{
  WlzValues 	val;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(sObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(sObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(sObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	if(sObj->type == WLZ_2D_DOMAINOBJ)
	{
	  val.v = WlzNewValueTb(sObj, tType, bgdV, &errNum);
	}
	else
	{
	  val.vox = WlzNewValuesVox(sObj, tType, bgdV, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  rObj = WlzMakeMain(sObj->type, sObj->domain, val, NULL, NULL,
			     &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && setFG)
	{
	  errNum = WlzGreySetValue(rObj, fgdV);
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  (void )WlzFreeObj(rObj);
	  rObj = NULL;
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New voxel value table.
* \ingroup      WlzAllocation
* \brief	From the domain of the given source object a new voxel
*		value table is created with the given grey type and
*		background value.
* \param	sObj			Source object.
* \param	gTType			Given grey table type for the plane
*					value tables.
* \param	bgdV			Background value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzVoxelValues	*WlzNewValuesVox(WlzObject *sObj, WlzObjectType gTType,
				 WlzPixelV bgdV, WlzErrorNum *dstErr)
{
  int		idx0;
  WlzDomain	*domP0;
  WlzValues	*valP0;
  WlzObject	*tObj;
  WlzVoxelValues *vox = NULL;
  WlzValues	tVal,
  		dumVal;
  WlzPlaneDomain *pDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dumVal.core = NULL;
  if(sObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(sObj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(sObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(sObj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    pDom = sObj->domain.p;
    vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
			      pDom->plane1, pDom->lastpl,
			      bgdV, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idx0 = pDom->plane1;
    domP0 = pDom->domains;
    valP0 = vox->values;
    while((idx0 <= pDom->lastpl) && (errNum == WLZ_ERR_NONE))
    {
      if((*domP0).core)
      {
	tObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domP0, dumVal,
			   NULL, NULL, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  tVal.v = WlzNewValueTb(tObj, gTType, bgdV, &errNum);
	  *valP0  = WlzAssignValues(tVal, NULL);
	  (void )WlzFreeObj(tObj);
	}
      }
      ++valP0;
      ++domP0;
      ++idx0;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vox);
}

/* function:     WlzMakeMain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Make a top-level woolz object assigning domain, values
and other pointers as required . The type is checked.
The domain is not checked for NULL although this should
only apply to the WLZ_EMPTY_OBJ.
*
* \return       Pointer to object structure.
* \param    type	Object type.
* \param    domain	Domain to be assigned using WlzAssignDomain()
* \param    values	Values to be attached using WlzAssignValues()
* \param    plist	Property list attached using WlzAssignPropertyList()
* \param    assoc	Associated Object attached using WlzAssignObject().
* \param    dstErr	error return values: WLZ_ERR_NONE,
 WLZ_ERR_MEM_ALLOC, WLZ_ERR_PARAM_DATA snd error values from WlzAssign
procedures.
* \par      Source:
*                WlzMakeStructs.c
*/

WlzObject *
WlzMakeMain(WlzObjectType 	type,
	    WlzDomain 		domain,
	    WlzValues 		values,
	    WlzPropertyList 	*plist,
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
    case WLZ_2D_POLYGON:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_3D_VIEW_STRUCT:
    case WLZ_AFFINE_TRANS:
    case WLZ_BOUNDLIST:
    case WLZ_CMESH_TRANS:
    case WLZ_CONTOUR:
    case WLZ_CMESH_2D:
    case WLZ_CMESH_2D5:
    case WLZ_CMESH_3D:
    case WLZ_CONV_HULL:
    case WLZ_EMPTY_OBJ:
    case WLZ_HISTOGRAM:
    case WLZ_LUT:
    case WLZ_MESH_TRANS:
    case WLZ_POINTS:
    case WLZ_PROPERTY_OBJ:
    case WLZ_RECTANGLE:
    case WLZ_TRANS_OBJ:
      obj->type = type;
      obj->linkcount = 0;
      obj->domain = WlzAssignDomain(domain, &errNum);
      if( errNum == WLZ_ERR_NONE ){
	obj->values = WlzAssignValues(values, &errNum);
      }
      /* property lists now more complicated  */
      if( errNum == WLZ_ERR_NONE ){
	if( plist ){
	  obj->plist = WlzAssignPropertyList(plist, &errNum);
	}
	else {
	  obj->plist = NULL;
	}
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
	      ("Makestructs - Obj %p type %d idom %p dl %d val %p"
	       " vl %d plist %p ass %p\n",
	       obj, obj->type, obj->domain.core, 
	       (obj->domain.core ? obj->domain.core->linkcount: 0),
	       obj->values.core, 
	       (obj->values.core ? obj->values.core->linkcount: 0),
	       obj->plist, obj->assoc));
      break;

    case WLZ_2D_POLYGON:
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	      ("Makestructs - Obj %p type %d pdom %p ass %p\n",
	       obj, obj->type, obj->domain.core, obj->assoc));
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

/* function:     WlzMakeValueTb    */
/*! 
* \ingroup      WlzAllocation
* \brief        Allocate and initialise space for a ragged-rectangle
value table only
*
* \return       Pointer to a ragged rectangle values table.
* \param    type	Values structure type. Must be a correct type as given
by WlzGreyTableType() with table_type = WLZ_GREY_TAB_RAGR.
* \param    l1	First line.
* \param    ll	Last line.
* \param    k1	First column.
* \param    backgrnd	Background pixel value
* \param    orig	Original object holding the grey-value data.
* \param    dstErr	Error return values: WLZ_ERR_NONE,
 WLZ_ERR_PARAM_DATA, WLZ_ERR_LINE_DATA, WLZ_ERR_MEM_ALLOC
* \par      Source:
*                WlzMakeStructs.c
*/
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
	 AlcCalloc(1, sizeof(WlzRagRValues) + valuelinespace)) == NULL ){
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
    case WLZ_GREY_RGBA:
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

/* function:     WlzMakeVoxelValueTb    */
/*! 
* \ingroup      WlzAllocation
* \brief        Allocate space for a voxel table
*
* \return       Pointer to a voxel values table.
* \param    type	Voxel table type.
* \param    p1	Fist plane.
* \param    pl	Last plane.
* \param    backgrnd	Background pixel value.
* \param    orig	Original object holding the voxel table
* \param    dstErr	Error return values: WLZ_ERR_NONE,
 WLZ_ERR_PARAM_DATA, WLZ_ERR_MEM_ALLOC
* \par      Source:
*                WlzMakeStructs.c
*/

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

/* function:     WlzMakeRectValueTb    */
/*! 
* \ingroup      WlzAllocation
* \brief        Make rectangular value table attaching values if set.
Note the values pointer is just copied and will not be
freed when the object is freed unless the freeptr is set
using AlcFreePointerPush().
*
* \return       Pointer to a rectangular value table.
* \param    type	Value table type. Must be a valid table type e.g.
as returned by WlzGreyTableType() with table_type = WLZ_GREY_TAB_RECT.
* \param    line1	First line
* \param    lastln	Last line.
* \param    kol1	First column
* \param    width	Width
* \param    backgrnd	Background pixel value.
* \param    values	Pointer to array of pixel values.
* \param    dstErr	Error return values: WLZ_ERR_NONE,
WLZ_ERR_PARAM_DATA, WLZ_ERR_MEM_ALLOC
* \par      Source:
*                WlzMakeStructs.c
*/

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
    if( (vtb = (WlzRectValues *) AlcMalloc(sizeof(WlzRectValues))) == NULL ){
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
    case WLZ_GREY_RGBA:
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


/* function:     WlzMakeInterval    */
/*! 
* \ingroup      WlzAllocation
* \brief        Attach and interval pointer to a an interval domain
in the appropriate place and set the number of intervals
on the line. Note this procedure makes no checks on the
arguments because it is often used in deeply nested loops.
*
* \return       Woolz error, values: WLZ_ERR_NONE
* \param    line	Line on which the interval array is to be attached.
* \param    idom	Pointer to the interval domain to hold the intervals.
* \param    nints	Number of intervals in the array.
* \param    intptr	Pointer to the set of intervals.
* \par      Source:
*                WlzMakeStructs.c
*/
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

/* function:     WlzMakeValueLine    */
/*! 
* \ingroup      WlzAllocation
* \brief        Attach the grey values for a valueline within a ragged-
rectangle value table. Not this procedure does not check
the input arguments because it is often at the core of
nested loops.
*
* \return       Woolz error, values: WLZ_ERR_NONE.
* \param    vtb	Pointer to a ragged rectangle value table.
* \param    line	Line for the values to be set.
* \param    k1	First column of grey interval.
* \param    kl	Last column of grey interval.
* \param    greyptr	Grey values pointer cast type int *.
* \par      Source:
*                WlzMakeStructs.c
*/
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

/*!
* \return       New domain object without values.
* \ingroup	WlzFeatures
* \brief        Constructs a domain from the union of marker domains with
*               a marker domain at each of the given vertex positions.
* \param        nVtx                    Number of vertices.
* \param        vtx                     Given vertices.
* \param        mType                   Marker type.
* \param        mSz                     Marker size. This is the radius of a
*					sphere marker. The marker size is
*					ignored for point markers.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject 	*WlzMakeMarkers(WlzVertexType vType,
                                int nVtx, WlzVertexP vtx,
                                WlzMarkerType mType, int mSz,
                                WlzErrorNum *dstErr)
{
  int           idx,
                dim = 0;
  WlzObjectType oType;
  WlzObject     *mObj = NULL;
  WlzObject     *tObj[4];
  WlzIVertex3   off;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  off.vtX = off.vtY = off.vtZ = 0;
  tObj[0] = tObj[1] = tObj[2] = tObj[3] = NULL;
  if((nVtx <= 0) || (mSz <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(vtx.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL; 
  } 
  else
  {
    switch(vType)
    { 
      case WLZ_VERTEX_I2:
        dim = 2;
        oType = WLZ_2D_DOMAINOBJ;
        break;
      case WLZ_VERTEX_I3:
        dim = 3;
        oType = WLZ_3D_DOMAINOBJ;
        break;  
      default:  
        errNum = WLZ_ERR_PARAM_DATA;
        break;
    }
  }                             
  if(errNum == WLZ_ERR_NONE)    
  {
    switch(mType)
    {
      case WLZ_MARKER_POINT:
        tObj[0] = WlzMakeSinglePixelObject(oType, 0, 0, 0, &errNum);
        break;
      case WLZ_MARKER_SPHERE:
        tObj[0] = WlzMakeSphereObject(oType, mSz, 0.0, 0.0, 0.0, &errNum);
        break;
      default:
        errNum = WLZ_ERR_PARAM_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj[1] = WlzMakeEmpty(&errNum);
  }
  for(idx = 0; (idx < nVtx) && (errNum == WLZ_ERR_NONE); ++idx)
  {
    if(dim == 2)
    {
      off.vtX = (vtx.i2 + idx)->vtX;
      off.vtY = (vtx.i2 + idx)->vtY;
    }
    else
    {
      off = *(vtx.i3 + idx);
    }
    tObj[2] = WlzShiftObject(tObj[0], off.vtX, off.vtY, off.vtZ, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      tObj[3] = WlzUnion2(tObj[1], tObj[2], &errNum);
    }
    WlzFreeObj(tObj[1]); tObj[1] = NULL;
    WlzFreeObj(tObj[2]); tObj[2] = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      tObj[1] = tObj[3];
      tObj[3] = NULL;
    }
  }
  WlzFreeObj(tObj[0]);
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = tObj[1];
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return       New domain object without values.
* \ingroup	WlzFeatures
* \brief	Creates a new domain object that is a formed from a lattice
* 		of markers covering the given domain.
* \param	gObj			Given spatial domain object.
* \param	mType			Marker type.
* \param	mSz			Marker size.
* \param	mSep			Marker separation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzMarkerLattice(WlzObject *gObj, WlzMarkerType mType,
                                  int  mSz, int mSep, WlzErrorNum *dstErr)
{
  int		nMrk = 0;
  WlzVertexP	mPos;
  WlzObject	*mObj = NULL,
  		*sObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  mPos.v = NULL;
  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((gObj->type != WLZ_2D_DOMAINOBJ) && (gObj->type != WLZ_3D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    WlzValues	val;
    WlzObject	*tObj = NULL;

    val.core = NULL;
    tObj = WlzMakeMain(gObj->type, gObj->domain, val, NULL, NULL, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzTransformType trType;
      WlzAffineTransform *tr = NULL;

      trType = (gObj->type == WLZ_2D_DOMAINOBJ)?
               WLZ_TRANSFORM_2D_AFFINE: WLZ_TRANSFORM_3D_AFFINE;
      tr = WlzAffineTransformFromPrimVal(trType, 0.0, 0.0, 0.0,
					 1.0 / mSep, 0.0, 0.0, 0.0, 0.0, 0.0,
      				         0, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        sObj = WlzAffineTransformObj(tObj, tr, WLZ_INTERPOLATION_NEAREST,
	                             &errNum);
      }
      (void )WlzFreeAffineTransform(tr);
    }
    (void )WlzFreeObj(tObj);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(gObj->type == WLZ_2D_DOMAINOBJ)
    {
      errNum = WlzVerticesFromObj2I(sObj, &nMrk, &(mPos.i2));
    }
    else
    {
      errNum = WlzVerticesFromObj3I(sObj, &nMrk, &(mPos.i3));
    }
  }
  (void )WlzFreeObj(sObj);
  if(errNum == WLZ_ERR_NONE)
  {
    int idx;
    WlzVertexType vType;

    if(gObj->type == WLZ_2D_DOMAINOBJ)
    {
      vType = WLZ_VERTEX_I2;
      for(idx = 0; idx < nMrk; ++idx)
      {
	WlzIVertex2 *p;

        p = mPos.i2 + idx;
	p->vtX *= mSep;
	p->vtY *= mSep;
      }
    }
    else
    {
      vType = WLZ_VERTEX_I3;
      for(idx = 0; idx < nMrk; ++idx)
      {
	WlzIVertex3 *p;

        p = mPos.i3 + idx;
	p->vtX *= mSep;
	p->vtY *= mSep;
	p->vtZ *= mSep;
      }
    }
    mObj = WlzMakeMarkers(vType, nMrk, mPos, mType, mSz, &errNum);
  }
  AlcFree(mPos.v);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(mObj);
    mObj = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	New point domain.
* \ingroup	WlzFeatures
* \brief	Creates a new point domain. A point domain consists of an
* 		array of vertices which are treated as seperate points.
* \param	type		Type of point vertices, which must be one of
* 				WLZ_POINTS_2I, WLZ_POINTS_2D, WLZ_POINTS_3I
* 				or WLZ_POINTS_3D.
* \param	nVtx		Number of points to copy to the new domain.
* \param	vtxP		Points to copy to the new domain. These
* 				must be of the correct type. If NULL then
* 				no points are copied.
* \param	maxVtx		The number of vertices for which space is
* 				allocated.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzPoints 	*WlzMakePoints(WlzObjectType type, int nVtx, WlzVertexP vtxP,
    			       int maxVtx, WlzErrorNum *dstErr)
{
  size_t	pntSz;
  WlzPoints	*pnt = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(maxVtx < 1)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    switch(type)
    {
      case WLZ_POINTS_2I:
	pntSz = sizeof(WlzIVertex2);
	break;
      case WLZ_POINTS_2D:
	pntSz = sizeof(WlzDVertex2);
	break;
      case WLZ_POINTS_3I:
	pntSz = sizeof(WlzIVertex3);
	break;
      case WLZ_POINTS_3D:
	pntSz = sizeof(WlzDVertex3);
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((pnt = (WlzPoints *)
	      AlcCalloc(sizeof(WlzPoints) + pntSz * maxVtx, 1)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    pnt->type = type;
    pnt->maxPoints = maxVtx;
    pnt->points.v = (void *)(pnt + 1);
    if(vtxP.v && (nVtx > 0))
    {
      pnt->nPoints = nVtx;
      (void )memcpy(pnt->points.v, vtxP.v, pntSz * nVtx);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pnt);
}

/*!
* \return	New point value table or NULL on error.
* \ingroup	WlzAllocation
* \brief	Makes a new point value table which covers the points of
* 		the given object.
* \param	obj			Given object, the domain of which is
* 					used to determine the value allocation.
* 					This must have type WLZ_POINTS.
* \param	rank			The rank of the individual values.
* \param	dim			The dimensions of individual indexed
* 					values.
* \param	vType			The type of the data in the individual
* 					values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPointValues	*WlzNewPointValues(WlzObject *obj,
                                   int rank, int *dim, WlzGreyType vType,
				   WlzErrorNum *dstErr)
{
  WlzPointValues *pVal = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_POINTS)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    pVal = WlzMakePointValues(obj->domain.pts->nPoints, rank, dim, vType,
    			      &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pVal);
}

/*!
* \return	New point value table or NULL on error.
* \ingroup	WlzAllocation
* \brief	Makes a new point value table with space for the requested
* 		number of points.
* \param	nP:			Number of points to allocate space
* 					for.
* \param	rank			The rank of the individual values.
* \param	dim			The dimensions of individual indexed
* 					values.
* \param	vType			The type of the data in the individual
* 					values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPointValues	*WlzMakePointValues(size_t nP,
                                    int rank, int *dim, WlzGreyType vType,
				    WlzErrorNum *dstErr)
{
  int		pSz,
  		gSz = 0,
  		nDat = 0;
  WlzPointValues *pv = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check parameters and compute the number of individual values per point. */
  if((nP == 0) || (rank < 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(rank == 0)
  {
    nDat = 1;
  }
  else
  {
    if(dim[0] < 1)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
    else
    {
      int	idx;

      nDat = dim[0];
      for(idx = 1; idx < rank; ++idx)
      {
	if(dim[idx] < 1)
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
	}
	nDat *= dim[idx];
      }
    }
  }
  /* Compute individual value size. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((vType == WLZ_GREY_BIT) ||
       ((gSz = WlzGreySize(vType)) == 0))
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
  }
  /* Allocate the point value table. */
  if(errNum == WLZ_ERR_NONE)
  {
    pSz = gSz * nDat;
    /* No dimension array allocated if rank == 0. */
    if(((pv = (WlzPointValues *)
              AlcCalloc(1, sizeof(WlzPointValues))) == NULL) ||
       ((pv->values.v = AlcCalloc(nP, pSz)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else if((rank > 0) &&
            ((pv->dim = (int *)AlcMalloc(rank * sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx;

    pv->type = WLZ_POINT_VALUES;
    pv->rank = rank;
    pv->vType = vType;
    pv->pSz = pSz;
    for(idx = 0; idx < rank; ++idx)
    {
      pv->dim[idx] = dim[idx];
    }
    pv->maxPoints = nP;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreePointValues(pv);
    pv = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pv);
}

/*! 
* \ingroup      WlzAllocation
* \brief        Make a polygon domain, allocating space and copying as
required:
<ul>
<li>vertices != NULL,  copy=0 - just plant the pointer </li>
<li>vertices != NULL,  copy=1 - allocate space and copy </li>
<li>vertices == NULL,  copy=0 - no vertex space allocated\n
probably an error!!</li>
<li>vertices == NULL,  copy=1 - allocate space for maxv vertices</li>
</ul>

*
* \return       Pointer to the initialised structure.
* \param    type	one of WLZ_POLYGON_INT, WLZ_POLYGON_FLOAT,
 WLZ_POLYGON_DOUBLE
 * \param    n	number of vertices if vertices!=NULL.
 * \param    vertices	vertices to be set, see type for value options.
 * \param    maxv	size of array if vertices!=NULL else number of vertices for which space it to be allocated.
 * \param    copy	copy flag see description for values.
 * \param    dstErr	Error return, values: WLZ_ERR_NONE, WLZ_ERR_PARAM_DATA, WLZ_ERR_MEM_ALLOC.
* \par      Source:
*                WlzMakeStructs.c
*/
WlzPolygonDomain *
WlzMakePolygonDomain(WlzObjectType	type,
	       int 		n,
	       WlzIVertex2 	*vertices,
	       int 		maxv,
	       int 		copy,
	       WlzErrorNum		*dstErr)
{
  WlzPolygonDomain	*p = NULL;
  int 			vertexsize = 0;
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
	if( vertices ){
	  (void )memcpy((void * )p->vtx, (void * )vertices, n * vertexsize);
	}
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(p);
}

/* function:     WlzMakeBoundList    */
/*! 
* \ingroup      WlzAllocation
* \brief        Allocate space and initialise a boundlist structure.
*
* \return       Pointer to a BoundList structure.
* \param    type	BoundList type, one of: WLZ_BOUNDLIST_PIECE or
 WLZ_BOUNDLIST_HOLE.
* \param    wrap	number of vertices by which the polygon is "wrapped"
 ie number of vertices overlapping at the beginning and end.
* \param    poly	polygon for this boundary structure
* \param    dstErr	Error return, values: WLZ_ERR_NONE,
 WLZ_ERR_PARAM_DATA, WLZ_ERR_MEM_ALLOC.
* \par      Source:
*                WlzMakeStructs.c
*/
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
    if( (blist = ( WlzBoundList *) AlcMalloc(sizeof( WlzBoundList))) == NULL ){
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

/* function:     WlzMakeIVertex    */
/*! 
* \ingroup      WlzAllocation
* \brief        Make an integer vertex array.
*
* \return       Pointer to the array of vertex structures.
* \param    nverts	Number of vertices.
* \param    dstErr	Error return, values: WLZ_ERR_NONE, WLZ_ERR_MEM_ALLOC.
* \par      Source:
*                WlzMakeStructs.c
*/
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

/* function:     WlzMakeRect    */
/*! 
* \ingroup      WlzAllocation
* \brief        Make a top-level rectangular object, setting values
if non-NULL, uses WlzMakeRectValueTb() to assign values 
which by default will not be freed when the object is
freed. The freeptr needs to be explicitly set by the
calling procedure. This is a convenience procedure calling
WlzMakeIntervalDomain() then WlzMakeRectValueTb() then WlzMakeMain().
*
* \return       Pointer to the top-level object.
* \param    line1	First line.
* \param    lastln	Last line
* \param    kol1	First column
* \param    lastkl	last column
* \param    pixeltype	Pixel type for the grey values. If WLZ_GREY_ERROR
*                       is given then no values are created.
* \param    grey_values	Pointer to the grey values array.
* \param    backgrnd	Background pixel value.
* \param    plist	Property list to be attached.
* \param    assoc_obj	Associated object.
* \param    dstErr	Error return, values: WLZ_ERR_NONE and valuea
 from WlzMakeRectValueTb() and WlzMakeMain().
* \par      Source:
*                WlzMakeStructs.c
*/
WlzObject *WlzMakeRect(int 			line1,
		       int 			lastln,
		       int 			kol1,
		       int 			lastkl,
		       WlzGreyType		pixeltype,
		       int 			*grey_values,
		       WlzPixelV		backgrnd,
		       WlzPropertyList		*plist,
		       WlzObject		*assoc_obj,
		       WlzErrorNum		*dstErr)
{
  WlzDomain	domain;
  WlzValues	values;
  WlzObject 	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  values.core = NULL;
  domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				   line1, lastln,
				   kol1, lastkl, &errNum);
  if((errNum == WLZ_ERR_NONE) &&
     (pixeltype != WLZ_GREY_ERROR) &&
     ((values.r = 
       WlzMakeRectValueTb(WlzGreyTableType(WLZ_GREY_TAB_RECT, pixeltype,
					   NULL),
			  line1, lastln, kol1, lastkl - kol1 + 1,
			  backgrnd, grey_values, &errNum)) == NULL)){
    WlzFreeDomain(domain);
  }

  if((errNum == WLZ_ERR_NONE) &&
     ((obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			 plist, assoc_obj, &errNum)) == NULL)){
    (void )WlzFreeDomain(domain);
    (void )WlzFreeValues(values);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New 3D domain object with values or NULL on error.
* \ingroup      WlzAllocation
* \brief	Creates a 3D domain object with values, for which the
*		domain is a cuboid and the values are of the given
*		type and initialized to have value zero.
* \param	plane1			First plane.
* \param	lastpl			Last plane.
* \param	line1			First line.
* \param	lastln			Last line.
* \param	kol1			First column.
* \param	lastkl			Last column.
* \param	pixType			Pixel type for the grey values. If
* 					WLZ_GREY_ERROR is given then no values
* 					are created.
* \param	bgdV			Background pixel value.
* \param	plist			Property list to be attached.
* \param	assocObj		Associated object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzMakeCuboid(int plane1, int lastpl,
			       int line1, int lastln,
			       int kol1, int lastkl,
			       WlzGreyType pixType, WlzPixelV bgdV,
			       WlzPropertyList *plist, WlzObject *assocObj,
			       WlzErrorNum *dstErr)
{
  int		pPos;
  size_t	arSz,
  		arElmSz = 0;
  void		*arDat = NULL;
  WlzDomain	dom,
  		dom2D;
  WlzDomain	*dom2DP = NULL;
  WlzObjectType	tbType;
  WlzValues	val,
  		val2D;
  WlzValues	*val2DP = NULL;
  WlzObject	*obj = NULL;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  dom2D.core = NULL;
  val2D.core = NULL;
  if(pixType != WLZ_GREY_ERROR)
  {
    if((arElmSz = WlzGreySize(pixType)) == 0)
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    arSz = (lastln - line1 + 1) * (lastkl - kol1 + 1);
    dom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
			       plane1, lastpl, line1, lastln, kol1, lastkl,
			       &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && (pixType != WLZ_GREY_ERROR))
  {
    val.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
    				  plane1, lastpl, bgdV, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    pPos = plane1				;
    dom2DP = dom.p->domains;
    if(val.core != NULL)
    {
      val2DP = val.vox->values;
      tbType = WlzGreyTableType(WLZ_GREY_TAB_RECT, pixType, &errNum);
    }
    while((errNum == WLZ_ERR_NONE) && (pPos <= lastpl))
    {
      dom2D.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
      				      line1, lastln, kol1, lastkl, &errNum);
      if(val.core != NULL)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  if((arDat = AlcCalloc(arSz, arElmSz)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  val2D.r = WlzMakeRectValueTb(tbType, line1, lastln,
				      kol1, lastkl - kol1 + 1,
				      bgdV, (int *)arDat, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  val2D.r->freeptr = AlcFreeStackPush(val2D.r->freeptr, arDat, &alcErr);
	  if(alcErr != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	++pPos;
        *dom2DP++ = WlzAssignDomain(dom2D, NULL);
	if(val.core != NULL)
	{
	  *val2DP++ = WlzAssignValues(val2D, NULL);
	  dom2D.core = NULL;
	  val2D.core = NULL;
	  arDat = NULL;
	}
      }
      else
      {
        (void )WlzFreeDomain(dom2D);
	if(val2D.core)
	{
	  (void )WlzFreeValues(val2D);
	}
	else
	{
	  AlcFree(arDat);
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.p->voxel_size[0] = 1.0;
    dom.p->voxel_size[1] = 1.0;
    dom.p->voxel_size[2] = 1.0;
    obj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, val, plist, assocObj, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    obj = NULL;
    (void )WlzFreeDomain(dom);
    (void )WlzFreeValues(val);
  }
  return(obj);
}

/* function:     WlzMakeHistogramDomain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Allocates space for a histogram domain with the space
for the given maximum number of bins of the appropriate type.
*
* \return       Pointer to object, NULL on error.
* \param    type	histogram type, one of: WLZ_HISTOGRAMDOMAIN_INT
 or WLZ_HISTOGRAMDOMAIN_FLOAT
* \param    maxBins	Maximum number of histogram bins.
* \param    dstErr	Error return, values: WLZ_ERR_NONE,
 WLZ_ERR_PARAM_DATA, WLZ_ERR_MEM_ALLOC.
* \par      Source:
*                WlzMakeStructs.c
*/
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
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
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

/* function:     WlzMakeEmpty    */
/*! 
* \ingroup      WlzAllocation
* \brief        Convenience procedure to make a top-level empty object.
*
* \return       Pointer to the Object structure, NULL on error.
* \param    dstErr	Error return, values from WlzMakeMain().
* \par      Source:
*                WlzMakeStructs.c
*/
WlzObject *WlzMakeEmpty(WlzErrorNum *dstErr)
{
  WlzDomain	domain;
  WlzValues	values;

  domain.core = NULL;
  values.core = NULL;
  return WlzMakeMain(WLZ_EMPTY_OBJ, domain, values, NULL, NULL, dstErr);
}

/* function:     WlzMakeHistogram    */
/*! 
* \ingroup      WlzAllocation
* \brief        Convenience procedure to make a top-level object with
a histogram domain.
*
* \return       Pointer to a Histogram object, NULL on error.
* \param    type	Type of the histogram domain, one of:
 WLZ_HISTOGRAMDOMAIN_INT or WLZ_HISTOGRAMDOMAIN_FLOAT.
* \param    maxBins	Maximum number of histogram bins.
* \param    dstErr	Error return, values: from WlzMakeHistogramDomain()
 and WlzMakeMain().
* \par      Source:
*                WlzMakeStructs.c
*/
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

/* function:     WlzNewGrey    */
/*! 
* \ingroup      WlzAllocation
* \brief        Make a top-level object with the same domain as iobj
(pointer copied and linkcount incremented) and new
valuetable with values copied from iobj. If iobj has no
valuetable then the returned object will have no value-
table. This allows a copy of a 2D grey-value object.
*
* \return       Pointer to the top-level object, NULL on error.
* \param    iobj	Input object which defines the domain and grey
values for which the new grey table will be defined.
 * \param    dstErr	Error return, values: WLZ_ERR_NONE,
 WLZ_ERR_OBJECT_TYPE and errors from WlzNewValueTb(), WlzMakeMain(),
 WlzInitGreyScan() and WlzNextGreyInterval().
* \par      Source:
*                WlzMakeStructs.c
*/
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

  /* check that the values aren't tiled. */
  if((errNum == WLZ_ERR_NONE) && iobj->values.core &&
     WlzGreyTableIsTiled(iobj->values.core->type)){
    errNum = WLZ_ERR_VALUES_TYPE;
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
	data_size = sizeof(WlzUByte);
	break;

      case WLZ_GREY_FLOAT:
	data_size = sizeof(float);
	break;

      case WLZ_GREY_DOUBLE:
	data_size = sizeof(double);
	break;

      case WLZ_GREY_RGBA:
	data_size = sizeof(WlzUInt);
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

/*! 
* \return       Pointer to new ragged-rectangle value table structure.
* \ingroup      WlzAllocation
* \brief        Create a value table of required type with the same size
*		and shape as the domain of obj. This must be of type
*		WLZ_2D_DOMAINOBJ.
* \param    obj				Given object in which the domain
* 					defines minimum coverage of the
* 					new value table.
* \param    type			Value table type.
* \param    backgrnd			Background pixel value.
* \param    dstErr			Destination error pointer, may be NULL.
*/
WlzRagRValues *WlzNewValueTb(WlzObject		*obj,
			     WlzObjectType	type,
			     WlzPixelV		backgrnd,
			     WlzErrorNum	*dstErr)
{
  WlzValues 		v;
  WlzDomain		idom;
  WlzIntervalWSpace	iwsp;
  WlzGreyP		g;
  int 			k1 = 0, table_size, bgd_val;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  v.v = NULL;
  g.v = NULL;
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
	table_size = WlzLineArea(obj, NULL) * sizeof(WlzUByte);
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

      case WLZ_GREY_RGBA:
	table_size = WlzLineArea(obj, NULL) * sizeof(WlzUInt);
	bgd_val = (int) backgrnd.v.rgbv;
	break;

      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
      }

      if( errNum == WLZ_ERR_NONE ){
	if((v.v = WlzMakeValueTb(type, idom.i->line1, idom.i->lastln,
				 idom.i->kol1, backgrnd, obj,
				 &errNum)) != NULL){
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

	    case WLZ_GREY_RGBA:
	      g.rgbp += iwsp.rgtpos - k1 +1;
	      break;

	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
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
				    backgrnd, NULL,
				    &errNum)) == NULL ){
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
	  (WlzUByte *) AlcMalloc(sizeof(WlzUByte)*(v.r->lastln - v.r->line1 +
	                                             1) * v.r->width);
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

      case WLZ_GREY_RGBA:
	v.r->values.rgbp  =
	  (WlzUInt *) AlcMalloc(sizeof(WlzUInt) *
	                         (v.r->lastln - v.r->line1 + 1) * v.r->width);
	break;

      default:
        errNum = WLZ_ERR_GREY_TYPE;
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

/* function:     WlzNewIDomain    */
/*! 
* \ingroup      WlzAllocation
* \brief        Make a copy of an intervaldomain.
*
* \return       Pointer to the new interval domain.
* \param    outDomType	Interval domain type, one of:
 WLZ_INTERVALDOMAIN_INTVL or WLZ_INTERVALDOMAIN_RECT.
* \param    inDom	Input domain to be copied.
* \param    dstErr	Error return, values: WLZ_ERR_NONE,
 WLZ_ERR_DOMAIN_TYPE, WLZ_ERR_MEM_ALLOC or errors from WlzMakeInterval(),
 WlzMakeIntervalDomain().
* \par      Source:
*                WlzMakeStructs.c
*/
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
	    jtvl->ileft = 0;
	    jtvl->iright = inDom->lastkl - inDom->kol1;
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

/* function:     WlzMakeContour    */
/*! 
* \ingroup      WlzAllocation
* \brief        Makes a new contour data structure.
*
* \return       Pointer to a WLZContour object.
* \param    dstErr	Error return, values: WLZ_ERR_NONE, WLZ_ERR_MEM_ALLOC.
* \par      Source:
*                WlzMakeStructs.c
*/
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

/*!
* \return	New index value table or NULL on error.
* \ingroup	WlzAllocation
* \brief	Makes a new indexed value table.
* \param	obj			Given object, the domain of which is
* 					used to determine the value allocation.
* \param	rank			The rank of the individual values.
* \param	dim			The dimensions of individual indexed
* 					values.
* \param	vType			The type of the data in the individual
* 					values.
* \param	attach			Specifies what the values are to be
* 					attached to.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIndexedValues *WlzMakeIndexedValues(WlzObject *obj,
                                    int rank, int *dim, WlzGreyType vType,
				    WlzValueAttach attach,
				    WlzErrorNum *dstErr)
{
  int		idx;
  size_t	bSz = 0,	                    /* AlcVector block size. */
		bCnt = 0,		           /* AlcVector block count. */
  		gSz = 0, 	        /* Size of data in individual value. */
                vSz = 0,                     /* Size of an individual value. */
		nDat = 0;	   /* Number of data in an individual value. */
  WlzIndexedValues *ixv = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check parameters and compute AlcVector allocation unit sizes. */
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(rank < 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(rank == 0)
  {
    nDat = 1;
  }
  else
  {
    if(dim[0] < 1)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
    else
    {
      nDat = dim[0];
      for(idx = 1; idx < rank; ++idx)
      {
	if(dim[idx] < 1)
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
	}
	nDat *= dim[idx];
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj->domain.core->type)
    {
      case WLZ_CMESH_2D:
	switch(attach)
	{
	  case WLZ_VALUE_ATTACH_NOD:
	    bSz = obj->domain.cm2->res.nod.vec->blkSz;
	    bCnt= obj->domain.cm2->res.nod.vec->blkCnt;
	    break;
	  case WLZ_VALUE_ATTACH_ELM:
	    bSz = obj->domain.cm2->res.elm.vec->blkSz;
	    bCnt= obj->domain.cm2->res.elm.vec->blkCnt;
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      case WLZ_CMESH_2D5:
	switch(attach)
	{
	  case WLZ_VALUE_ATTACH_NOD:
	    bSz = obj->domain.cm2d5->res.nod.vec->blkSz;
	    bCnt= obj->domain.cm2d5->res.nod.vec->blkCnt;
	    break;
	  case WLZ_VALUE_ATTACH_ELM:
	    bSz = obj->domain.cm2d5->res.elm.vec->blkSz;
	    bCnt= obj->domain.cm2d5->res.elm.vec->blkCnt;
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      case WLZ_CMESH_3D:
	switch(attach)
	{
	  case WLZ_VALUE_ATTACH_NOD:
	    bSz = obj->domain.cm3->res.nod.vec->blkSz;
	    bCnt= obj->domain.cm3->res.nod.vec->blkCnt;
	    break;
	  case WLZ_VALUE_ATTACH_ELM:
	    bSz = obj->domain.cm3->res.elm.vec->blkSz;
	    bCnt= obj->domain.cm3->res.elm.vec->blkCnt;
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  /* Compute value size. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((vType == WLZ_GREY_BIT) ||
       ((gSz = WlzGreySize(vType)) == 0))
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
  }
  /* Allocate the indexed value table. */
  if(errNum == WLZ_ERR_NONE)
  {
    vSz = gSz * nDat;
    /* No dimension array allocated if rank == 0. */
    if(((ixv = (WlzIndexedValues *)
              AlcCalloc(1, sizeof(WlzIndexedValues))) == NULL) ||
       ((ixv->values = AlcVectorNew(bSz * bCnt, vSz, bSz, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else if((rank > 0) &&
            ((ixv->dim = (int *)AlcMalloc(rank * sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ixv->type = WLZ_INDEXED_VALUES;
    ixv->rank = rank;
    ixv->vType = vType;
    ixv->attach = attach;
    for(idx = 0; idx < rank; ++idx)
    {
      ixv->dim[idx] = dim[idx];
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeIndexedValues(ixv);
    ixv = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ixv);
}
