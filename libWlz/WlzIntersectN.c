#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzIntersectN_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzIntersectN.c
* \author       Richard Baldock
* \date         August 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Intersection of N woolz domain objects.
* \ingroup	WlzBinaryOps
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

extern WlzObject *WlzIntersect3d(WlzObject	**objs,
				 int		n,
				 int		uvt,
				 WlzErrorNum	 *dstErr);


/* function:     WlzIntersectN    */
/*! 
* \ingroup      WlzBinaryOps
* \brief        Calculate the intersection of a set of objects. If
 uvt=0 calculate domain only, uvt=1 calculate the mmean grey-value at
 each point. Input objects must be all non-NULL and domain objects of
 the same type i.e. either 2D or 3D otherwise an error is returned.
*
* \return       Intersection object with grey-table as required, if the intersection is empty returns WLZ_EMPTY_OBJ, NULL on error.
* \param    n	number of input objects
* \param    objs	input object array
* \param    uvt	grey-table copy flag (1 - copy, 0 - no copy)
* \param    dstErr	error return.
* \par      Source:
*                WlzIntersectN.c
*/
WlzObject *WlzIntersectN(
  int 	n,
  WlzObject **objs,
  int 	uvt,
  WlzErrorNum *dstErr)
{
  WlzObject 		*obj = NULL;
  WlzIntervalDomain 	*idom;
  WlzInterval 		*itvl, *jtvl;
  WlzIntervalWSpace 	*iwsp;
  WlzIntervalWSpace 	*biwsp,*tiwsp,niwsp;
  WlzGreyWSpace 	*gwsp,ngwsp;
  WlzDomain		domain;
  WlzValues		values;
  WlzObjectType		type;
  WlzPixelV		backg;
  WlzGreyP		greyptr;
  WlzGreyV		gv;
  int 			i, j, k, l, inttot, change, lwas, nints;
  int 			line1, lastln;
  int 			kol1, lastkl;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /*
   * check pointers
   */
  /* intersecction of no objects is an empty domain */
  domain.i = NULL;
  values.v = NULL;
  if( n < 1 )
  {
    return WlzMakeEmpty(dstErr);
  }

  /* array pointer == NULL or any object pointer == NULL is an error */
  if(objs == NULL){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    for(i=0; i<n; i++){
      if(objs[i] == NULL){
	errNum = WLZ_ERR_OBJECT_NULL;
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    obj = NULL;
    if(dstErr) {
      *dstErr = errNum;
    }
    return(obj);
  }

  /* check type of first object */
  switch( objs[0]->type ){

  case WLZ_2D_DOMAINOBJ:
    break;

  case WLZ_3D_DOMAINOBJ:
    return WlzIntersect3d(objs, n, uvt, &errNum);

  case WLZ_EMPTY_OBJ:
    return WlzMakeEmpty(dstErr);

  default:
    errNum = WLZ_ERR_OBJECT_TYPE;
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;

  }

  /* check number */
  if (n == 1){
    obj = WlzMakeMain(objs[0]->type, objs[0]->domain, objs[0]->values,
		      NULL, NULL, &errNum);
    if(dstErr) {
      *dstErr = errNum;
    }
    return(obj);
  }

  /* check all objects are non-empty and have the same type
     Note an empty object is not an error */
  for (i=0; i<n; i++){
    if( objs[i]->type != objs[0]->type ){
      if( objs[i]->type == WLZ_EMPTY_OBJ ){
	return WlzMakeEmpty(dstErr);
      }
      else {
	errNum = WLZ_ERR_OBJECT_TYPE;
	if(dstErr) {
	  *dstErr = errNum;
	}
	return NULL;
      }
    }

    /* check for size */
    if( WlzIsEmpty(objs[i], &errNum) ){
      return WlzMakeEmpty(dstErr);
    }
    else {
      if( errNum != WLZ_ERR_NONE ){
	if(dstErr) {
	  *dstErr = errNum;
	}
	return NULL;
      }
    }
  }

  /*
   * Find the line and column bounds of the intersection.
   */
  line1 = objs[0]->domain.i->line1;
  lastln = objs[0]->domain.i->lastln;
  kol1 = objs[0]->domain.i->kol1;
  lastkl = objs[0]->domain.i->lastkl;
  for (i=1; i<n; i++) {
    idom = objs[i]->domain.i;
    if (line1 < idom->line1)
    {
      line1 = idom->line1;
    }
    if (lastln > idom->lastln)
    {
      lastln = idom->lastln;
    }
    if (kol1 < idom->kol1)
    {
      kol1 = idom->kol1;
    }
    if (lastkl > idom->lastkl)
    {
      lastkl = idom->lastkl;
    }
  }
  if( lastkl < kol1 || lastln < line1 ){
    return WlzMakeEmpty(dstErr);
  }

  /*
   * Count the individual intervals so that sufficient space
   * for the intersection may be allocated.
   */
  inttot=0;
  for (i=0; (i < n) && (errNum == WLZ_ERR_NONE); i++) {
    inttot += WlzIntervalCount(objs[i]->domain.i, &errNum);
  }
  if(errNum != WLZ_ERR_NONE) {
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;
  }

  /*
   * Set up domain, value table structures, and object.
   */

  if( (idom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				    line1, lastln, 0, lastkl-kol1,
				    &errNum)) == NULL ){
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;
  }
  if( (itvl = (WlzInterval *)
       AlcMalloc (inttot * sizeof(WlzInterval))) == NULL ){
    WlzFreeIntervalDomain(idom);
    errNum = WLZ_ERR_MEM_ALLOC;
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;
  }

  idom->freeptr = AlcFreeStackPush(idom->freeptr, (void *)itvl, NULL);
  lwas = line1;
  jtvl = itvl;
  nints = 0;
  domain.i = idom;
  values.v = NULL;
  if( (obj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			 domain, values, NULL, NULL, &errNum)) == NULL ){
    WlzFreeIntervalDomain(idom);
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;
  }

  /*
   * allocate space for workspaces
   */
  if( (iwsp = (WlzIntervalWSpace *)
       AlcMalloc(n * sizeof(WlzIntervalWSpace))) == NULL ){
    WlzFreeObj( obj );
    errNum = WLZ_ERR_MEM_ALLOC;
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;
  }
  biwsp = iwsp;
  tiwsp = biwsp + n;

  /*
   * Construct the intersection object's table of intervals.
   * Initialise scanning on each object/workspace combination.
   * Scan synchronously, setting up the intersection of
   * overlapping intervals.  Needs a clear head !!
   */
  for (i=0; (i < n) && (errNum == WLZ_ERR_NONE); i++) {
    errNum = WlzInitRasterScan(objs[i],iwsp,WLZ_RASTERDIR_ILIC);
    if(errNum == WLZ_ERR_NONE) {
      errNum = WlzNextInterval(iwsp++);
    }
  }
  if(errNum != WLZ_ERR_NONE) {
    WlzFreeObj( obj );
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;
  }

  l = lwas;
  for (;;) {
    /*
     * find next line of intersection
     */
    do {
      change = 0;
      for (iwsp=biwsp; iwsp<tiwsp; iwsp++) {
	if (iwsp->linpos > l) {
	  l = iwsp->linpos;
        }
      }
      for (iwsp=biwsp; iwsp<tiwsp; iwsp++) {
	while (iwsp->linpos < l) {
	  if ((errNum = WlzNextInterval(iwsp)) != WLZ_ERR_NONE) {
	    if(errNum == WLZ_ERR_EOO)
	    {
	      errNum = WLZ_ERR_NONE;
	    }
	    goto firstfinished;
	  }
	  if (iwsp->linpos > l) {
	    change = 1;
	  }
	}
      }
    } while (change != 0);
    /*
     * find next interval of intersection
     */
    kol1 = biwsp->lftpos;
    lastkl = biwsp->rgtpos;
    while(errNum == WLZ_ERR_NONE) {
      do {
	change = 0;
	for (iwsp=biwsp; iwsp<tiwsp; iwsp++)
	  if (iwsp->lftpos > lastkl) {
	    kol1 = iwsp->lftpos;
	    lastkl = iwsp->rgtpos;
	  } else if (iwsp->lftpos > kol1)
	    kol1 = iwsp->lftpos;
	for (iwsp=biwsp; iwsp<tiwsp; iwsp++) {
	  while (iwsp->rgtpos < kol1) {
	    if ((errNum = WlzNextInterval(iwsp)) != WLZ_ERR_NONE) {
	      if(errNum == WLZ_ERR_EOO)
	      {
	        errNum = WLZ_ERR_NONE;
	      }
	      goto firstfinished;
	    }
	    if (iwsp->linpos != l) {
	      l = iwsp->linpos;
	      goto jumpline;
	    }
	    if (iwsp->lftpos > kol1)
	      change = 1;
	  }
	}
      } while ((change != 0) && (errNum == WLZ_ERR_NONE));
      if(errNum == WLZ_ERR_NONE)
      {
	for (iwsp=biwsp; iwsp < tiwsp; iwsp++) {
	  if (iwsp->rgtpos <= lastkl) {
	    lastkl = iwsp->rgtpos;
	  }
	}
	if (lastkl >= kol1) {
	  itvl->ileft = kol1 - idom->kol1;
	  itvl->iright = lastkl - idom->kol1;
	  if (l == lwas) {
	    nints++;
	  } else {
	    errNum = WlzMakeInterval(lwas,idom,nints,jtvl);
	    for (j = lwas+1; (j < l) && (errNum == WLZ_ERR_NONE); j++) {
	      errNum = WlzMakeInterval(j,idom,0,NULL);
	    }
	    if(errNum == WLZ_ERR_NONE) {
	      lwas = l;
	      nints = 1;
	      jtvl = itvl;
	    }
	  }
	  itvl++;
	}
	kol1 = lastkl+1;
      }
    }

  jumpline:
    ;
  }

firstfinished:
  if(errNum != WLZ_ERR_NONE) {
    WlzFreeObj(obj);
    if(dstErr) {
      *dstErr = errNum;
    }
    return NULL;
  }

  errNum = WlzMakeInterval(lwas,idom,nints,jtvl);
  for (j = lwas+1; (j <= lastln) && (errNum == WLZ_ERR_NONE); j++) {
    errNum = WlzMakeInterval(j,idom,0,NULL);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(obj);
    AlcFree((void *) biwsp);
    if(dstErr)
    {
      *dstErr = errNum;
    }
    return NULL;
  }

  /*
   * standardise the interval list (remove leading and trailing
   * empty lines, set kol1 so that at least one interval commences
   * at zero, set lastkl correctly)
   */
  (void )WlzStandardIntervalDomain(idom);

  /* check for error or empty */
  nints = WlzIntervalCount(idom, &errNum);
  if( nints < 0 ){
    /* error */
    WlzFreeObj(obj);
    AlcFree((void *) biwsp);
    if(dstErr)
    {
      *dstErr = errNum;
    }
    return NULL;
  }
  else if( nints == 0 ){
    /* empty object */
    WlzFreeObj(obj);
    AlcFree((void *) biwsp);
    return WlzMakeEmpty(dstErr);
  }

  if (uvt != 0) {
    WlzGreyType	grey_type;
    if( (gwsp = (WlzGreyWSpace *)
	 AlcMalloc (n * sizeof (WlzGreyWSpace))) == NULL ){
      WlzFreeObj(obj);
      AlcFree((void *) biwsp);
      errNum = WLZ_ERR_MEM_ALLOC;
      if(dstErr)
      {
        *dstErr = errNum;
      }
      return NULL;
    }

    /* construct an empty "ragged-rectangle" (type 1  or 2) greytable
       choosing the grey-value type from the first object in the list */
    backg = WlzGetBackground(objs[0], &errNum);
    if(errNum == WLZ_ERR_NONE) {
      grey_type = WlzGreyTableTypeToGreyType(objs[0]->values.core->type,
      					     &errNum);
    }
    if(errNum == WLZ_ERR_NONE) {
      type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, grey_type, &errNum);
    }
    if((errNum != WLZ_ERR_NONE) ||
       (values.v = WlzNewValueTb(obj, type, backg, &errNum)) == NULL ){
      WlzFreeObj( obj );
      AlcFree((void *) gwsp);
      AlcFree((void *) biwsp);
      if(dstErr)
      {
        *dstErr = errNum;
      }
      return NULL;
    }
    obj->values = WlzAssignValues(values, &errNum);
    if(errNum != WLZ_ERR_NONE) {
      WlzFreeObj( obj );
      AlcFree((void *) gwsp);
      AlcFree((void *) biwsp);
      if(dstErr)
      { 
	*dstErr = errNum;
      }
      return NULL;
    }
	
    /*
     * fill the grey table with mean of grey values.
     */
    /* initialise the work-spaces and check pixel type */
    errNum = WlzInitGreyScan(obj, &niwsp, &ngwsp);
    iwsp = biwsp;
    for (i=0; (i < n) && (errNum == WLZ_ERR_NONE); i++) {
      errNum = WlzInitGreyScan(objs[i], iwsp, &gwsp[i]);
      if(errNum == WLZ_ERR_NONE) {
	errNum = WlzNextGreyInterval(iwsp++);
      }
      if( gwsp[i].pixeltype != grey_type ){
        errNum = WLZ_ERR_GREY_TYPE;
      }
    }
    if(errNum != WLZ_ERR_NONE) {
      WlzFreeObj( obj );
      AlcFree((void *) gwsp);
      AlcFree((void *) biwsp);
      if(dstErr) {
        *dstErr = errNum;
      }
      return NULL;
    }

    while (WlzNextGreyInterval(&niwsp) == WLZ_ERR_NONE) {
      l = niwsp.linpos;
      greyptr = ngwsp.u_grintptr;
      switch( ngwsp.pixeltype ){

      case WLZ_GREY_INT:
	for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	  gv.inv = 0;
	  for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	    while(iwsp->linrmn >= 0 &&
		  (iwsp->linpos < l ||
		   (iwsp->linpos == l && iwsp->rgtpos < k))){
	      (void )WlzNextGreyInterval(iwsp);
	    }
	    gv.inv += *(gwsp[i].u_grintptr.inp + k - iwsp->lftpos);
	  }
	  *greyptr.inp = gv.inv / n;
	  greyptr.inp++;
	}
	break;

      case WLZ_GREY_SHORT:
	for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	  gv.shv = 0;
	  for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	    while(iwsp->linrmn >= 0 &&
		  (iwsp->linpos < l ||
		   (iwsp->linpos == l && iwsp->rgtpos < k))){
	      (void )WlzNextGreyInterval(iwsp);
	    }
	    gv.shv += *(gwsp[i].u_grintptr.shp + k - iwsp->lftpos);
	  }
	  *greyptr.shp = (short )(gv.shv / n);
	  greyptr.shp++;
	}
	break;

      case WLZ_GREY_UBYTE:
	for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	  gv.ubv = 0;
	  for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	    while(iwsp->linrmn >= 0 &&
		  (iwsp->linpos < l ||
		   (iwsp->linpos == l && iwsp->rgtpos < k))){
	      (void )WlzNextGreyInterval(iwsp);
	    }
	    gv.ubv += *(gwsp[i].u_grintptr.ubp + k - iwsp->lftpos);
	  }
	  *greyptr.ubp = (WlzUByte )(gv.ubv / n);
	  greyptr.ubp++;
	}
	break;

      case WLZ_GREY_FLOAT:
	for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	  gv.flv = 0;
	  for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	    while(iwsp->linrmn >= 0 &&
		  (iwsp->linpos < l ||
		   (iwsp->linpos == l && iwsp->rgtpos < k))){
	      (void )WlzNextGreyInterval(iwsp);
	    }
	    gv.flv += *(gwsp[i].u_grintptr.flp + k - iwsp->lftpos);
	  }
	  *greyptr.flp = gv.flv / n;
	  greyptr.flp++;
	}
	break;

      case WLZ_GREY_DOUBLE:
	for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	  gv.dbv = 0;
	  for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	    while(iwsp->linrmn >= 0 &&
		  (iwsp->linpos < l ||
		   (iwsp->linpos == l && iwsp->rgtpos < k))){
	      (void )WlzNextGreyInterval(iwsp);
	    }
	    gv.dbv += *(gwsp[i].u_grintptr.dbp + k - iwsp->lftpos);
	  }
	  *greyptr.dbp = gv.dbv / n;
	  greyptr.dbp++;
	}
	break;

      case WLZ_GREY_RGBA: /* RGBA to be done again RAB */
	for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	  gv.rgbv = 0;
	  for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	    while(iwsp->linrmn >= 0 &&
		  (iwsp->linpos < l ||
		   (iwsp->linpos == l && iwsp->rgtpos < k))){
	      (void )WlzNextGreyInterval(iwsp);
	    }
	    gv.rgbv += *(gwsp[i].u_grintptr.rgbp + k - iwsp->lftpos);
	  }
	  *greyptr.rgbp = gv.rgbv / n;
	  greyptr.rgbp++;
	}
	break;

      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;

      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
    AlcFree((void *) gwsp);
  }
  AlcFree((void *) biwsp);

  if(dstErr) {
    *dstErr = errNum;
  }
  return obj;
}
