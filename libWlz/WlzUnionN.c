#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzUnionN_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzUnionN.c
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
* \brief	Computes the set union of N objects.
* \ingroup	WlzBinaryOps
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

extern WlzObject *WlzUnion3d(int 	n,
			     WlzObject 	**objs,
			     int 	uvt,
			     WlzErrorNum *dstErr);

/* function:     WlzUnionN    */
/*! 
* \ingroup      WlzBinaryOps
* \brief        Calculate the set union of an array of domain objects.
 Domians only unless uvt non-zero in which case make an average grey
 table. Note background values are used in the averaging process. All
 objects must be domain objects of the same type (2D or 3D) unless
 WLZ_EMPTY_OBJ, NULL input objects are an error.
*
* \return       Union of the array of object.
* \param    n	number of input objects
* \param    objs	input object array
* \param    uvt	grey-table copy flag, copy if non-zero.
* \param    dstErr	error return.
* \par      Source:
*                WlzUnionN.c
*/
WlzObject *WlzUnionN(
  int		n,
  WlzObject 	**objs,
  int 		uvt,
  WlzErrorNum	*dstErr)
{
  WlzObject		*obj=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzIntervalDomain	*idom;
  WlzInterval		*itvl, *jtvl;
  WlzIntervalWSpace	*iwsp;
  WlzIntervalWSpace	*biwsp, *tiwsp, niwsp;
  WlzGreyWSpace		*gwsp, ngwsp;
  WlzObjectType		type;
  int 			i, j, k, l;
  int			inttot, numactive, change, lwas, nints, noverlap;
  WlzPixelV		backg;
  int			line1, lastln;
  int			kol1,lastkl;
  WlzGreyV		gv;
  WlzGreyP		greyptr;
  int			**locbuff;
  int			*locbuffs;
  int			span;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* preliminary stuff - count of non-NULL objects, note WLZ_EMPTY_OBJs
     are ignored but NULL objects are an error */
  if( n < 1 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else {
    for (i=0; i<n ; i++ ){
      if ( objs[i] == NULL ){
	errNum = WLZ_ERR_OBJECT_NULL;
	break;
      }
      if( objs[i]->type == WLZ_EMPTY_OBJ ){
	obj = objs[i];
	for ( j=i; j<n-1 ; j++ ){
	  objs[j] = objs[j+1];
	}
	objs[n-1] = obj;
	n--;
	i--;
      }
    }
  }
	
  /* n has been checked therefore no objects implies all empty */
  if( (errNum == WLZ_ERR_NONE) && (n < 1) ){
    return WlzMakeEmpty(dstErr);
  }

  /* now check they are all of the same type */
  if( errNum == WLZ_ERR_NONE ){
    for(i=1; i < n; i++){
      if( objs[i]->type != objs[0]->type ){
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* now test the type note empty objects have been discarded */
  if( errNum == WLZ_ERR_NONE ){
    switch( objs[0]->type ){

    case WLZ_2D_DOMAINOBJ:
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzUnion3d(n, objs, uvt, dstErr);

    case WLZ_TRANS_OBJ:
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* now discard empty objects */
  if( errNum == WLZ_ERR_NONE ){
    for (i=0; i<n ; i++ ){
      if( WlzIsEmpty(objs[i], NULL) ){
	obj = objs[i];
	for ( j=i; j<n-1 ; j++ ){
	  objs[j] = objs[j+1];
	}
	objs[n-1] = obj;
	n--;
	i--;
      }
    }
  }
  obj = NULL;

  /* recheck number of objects */
  if( (errNum == WLZ_ERR_NONE) && (n < 1) ){
    return WlzMakeEmpty(dstErr);
  }
  if( (errNum == WLZ_ERR_NONE) && (n == 1) ){
    return WlzMakeMain(objs[0]->type, objs[0]->domain,
		       objs[0]->values, NULL, NULL, dstErr);
  }

  /*
   * find the line and column bounds of the union.
   */
  if( errNum == WLZ_ERR_NONE ){
    line1 = objs[0]->domain.i->line1;
    lastln = objs[0]->domain.i->lastln;
    kol1 = objs[0]->domain.i->kol1;
    lastkl = objs[0]->domain.i->lastkl;
    for (i=1; i<n; i++) {
      idom = objs[i]->domain.i;
      if (line1 > idom->line1)
	line1 = idom->line1;
      if (lastln < idom->lastln)
	lastln = idom->lastln;
      if (kol1 > idom->kol1)
	kol1 = idom->kol1;
      if (lastkl < idom->lastkl)
	lastkl = idom->lastkl;
    }
    span = lastkl - kol1 +1 ;
    if( (locbuff = (int **) AlcMalloc((n+1)*sizeof(int *))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  /* space must be allocated for the largest variety of grey-value */
  if( errNum == WLZ_ERR_NONE ){
    if( (locbuffs = (int *) AlcMalloc(sizeof(double)*span*(n+1))) == NULL ){
      AlcFree((void *) locbuff);
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      for(i=0; i <= n; i++){
	locbuff[i] = locbuffs + i*span;
      }
    }
  }
  /*
   * count the individual intervals so that sufficient space
   * for the union may be allocated.
   */
  if( errNum == WLZ_ERR_NONE ){
    inttot=0;
    for(i=0; i < n; i++){
      inttot += WlzIntervalCount(objs[i]->domain.i, &errNum);
    }
  }

  /*
   * set up domain, value table structures, and object.
   */
  if( errNum == WLZ_ERR_NONE ){
    if( (idom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				      line1,lastln,kol1,lastkl,
				      &errNum)) == NULL ){
      AlcFree((void *) locbuffs);
      AlcFree((void *) locbuff);
    }
    else if( (itvl = (WlzInterval *)
	      AlcMalloc(inttot * sizeof(WlzInterval))) == NULL){
      AlcFree((void *) locbuffs);
      AlcFree((void *) locbuff);
      WlzFreeIntervalDomain(idom);
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      idom->freeptr = AlcFreeStackPush(idom->freeptr, (void *)itvl, NULL);
      lwas = line1;
      jtvl = itvl;
      nints = 0;
      domain.i = idom;
      values.v = NULL;
      if( (obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			     NULL, NULL, &errNum)) == NULL ){
	WlzFreeIntervalDomain( idom );
	AlcFree((void *) locbuffs);
	AlcFree((void *) locbuff);
	return(NULL);
      }
    }
  }    

  /*
   * allocate space for workspaces
   */
  if( errNum == WLZ_ERR_NONE ){
    if( (iwsp = (WlzIntervalWSpace *)
        AlcMalloc (n * sizeof (WlzIntervalWSpace))) == NULL ){
      WlzFreeObj( obj );
      AlcFree((void *) locbuffs);
      AlcFree((void *) locbuff);
      errNum = WLZ_ERR_MEM_ALLOC;
      obj = NULL;
    }
    else {
      biwsp = iwsp;
      tiwsp = iwsp + n;
    }
  }

  /*
   * Construct the union object's table of intervals.
   * Initialise scanning on each object/workspace combination.
   * Scan synchronously, setting up the union of adjacent and
   * overlapping intervals.  Needs a clear head !!
   */
  if( errNum == WLZ_ERR_NONE ){
    for (i=0; i<n; i++) {
      WlzInitRasterScan(objs[i], iwsp, WLZ_RASTERDIR_ILIC);
      WlzNextInterval(iwsp++);
    }
    numactive = n;

    /*
     * find next line and left hand end of next interval of union
     */
    while (numactive > 0) {
      /*
       * find first remaining active object
       */
      iwsp = biwsp;
      while( iwsp->linrmn < 0 ){
	iwsp++;
      }
      /*
       * find minimum line number of remaining active intervals
       */
      l = iwsp->linpos;
      kol1 = iwsp->lftpos;
      lastkl = iwsp->rgtpos;
      for (iwsp++; iwsp<tiwsp; iwsp++)
	if (iwsp->linrmn >= 0 && iwsp->linpos < l) {
	  l = iwsp->linpos;
	  kol1 = iwsp->lftpos;
	  lastkl = iwsp->rgtpos;
	}
      /*
       * find left-most interval in this line
       */
      for (iwsp=biwsp; iwsp<tiwsp; iwsp++)
	if (iwsp->linrmn >= 0 && iwsp->linpos == l && iwsp->lftpos < kol1) {
	  kol1 = iwsp->lftpos;
	  lastkl = iwsp->rgtpos;
	}
      /*
       * construct maximal interval with current left end-point
       */
      do {
	change = 0;
	for (iwsp=biwsp; iwsp<tiwsp; iwsp++) {
	  while( iwsp->linrmn >= 0 && iwsp->linpos == l
		&& iwsp->lftpos <= lastkl+1 ){
	    if (iwsp->rgtpos > lastkl) {
	      lastkl = iwsp->rgtpos;
	      change = 1;
	    }
	    if (WlzNextInterval(iwsp) != WLZ_ERR_NONE) {
	      numactive--;
	    }
	  }
	}
      } while (change == 1);
      itvl->ileft = kol1 - idom->kol1;
      itvl->iright = lastkl - idom->kol1;
      if (l == lwas)
	nints++;
      else {
	(void) WlzMakeInterval(lwas,idom,nints,jtvl);
	for (j = lwas+1; j<l; j++) {
	  (void) WlzMakeInterval(j,idom,0,NULL);
	}
	lwas = l;
	nints = 1;
	jtvl = itvl;
      }
      itvl++;
    }
    (void) WlzMakeInterval(lwas,idom,nints,jtvl);
    for (j = lwas+1; j<=lastln; j++) {
      (void) WlzMakeInterval(j,idom,0,NULL);
    }
  }

  /* now deal with the grey-values if required */
  if( (errNum == WLZ_ERR_NONE) && (uvt != 0) ){
    WlzGreyType	grey_type;

    if( (gwsp = (WlzGreyWSpace *)
	 AlcMalloc (n * sizeof (WlzGreyWSpace))) == NULL){
      WlzFreeObj( obj );
      AlcFree((void *) locbuffs);
      AlcFree((void *) locbuff);
      AlcFree((void *) biwsp);
      errNum = WLZ_ERR_MEM_ALLOC;
      obj = NULL;
    }
    
    /* construct an empty "ragged-rectangle" greytable
       with appropriate grey-type */
    if( errNum == WLZ_ERR_NONE ){
      backg = WlzGetBackground(objs[0], NULL);
      grey_type = WlzGreyTableTypeToGreyType(objs[0]->values.core->type,
					     NULL);
      type = WlzGreyTableType(WLZ_GREY_TAB_RAGR, grey_type, NULL);
      if( (values.v = WlzNewValueTb(obj, type, backg, &errNum)) == NULL ){
	WlzFreeObj( obj );
	AlcFree((void *) locbuffs);
	AlcFree((void *) locbuff);
	AlcFree((void *) biwsp);
	obj = NULL;
      }
      else {
	obj->values = WlzAssignValues(values, NULL);
      }
    }

    /* fill the grey table.  Where more than one input objects
       overlap, take mean of grey values. */
    if( errNum == WLZ_ERR_NONE ){
      WlzInitGreyScan(obj, &niwsp, &ngwsp);
      iwsp = biwsp;
      for (i=0; i<n; i++) {
	WlzInitGreyScan(objs[i], iwsp, &gwsp[i]);
	WlzNextGreyInterval(iwsp++);
	if( gwsp[i].pixeltype != grey_type ){
	  AlcFree((void *) gwsp);
	  AlcFree((void *) locbuffs);
	  AlcFree((void *) locbuff);
	  AlcFree((void *) biwsp);
	  WlzFreeObj( obj );
	  obj = NULL;
	  errNum = WLZ_ERR_GREY_TYPE;
	}      
      }
    }

    if( errNum == WLZ_ERR_NONE ){
      while (WlzNextGreyInterval(&niwsp) == WLZ_ERR_NONE) {
	l = niwsp.linpos;
	greyptr = ngwsp.u_grintptr;
	switch( ngwsp.pixeltype ){

	case WLZ_GREY_INT:
	  for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	    noverlap = 0;
	    gv.inv = 0;
	    for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	      while (iwsp->linrmn >= 0
		     && (iwsp->linpos < l
			 || (iwsp->linpos == l && iwsp->rgtpos < k))){
		WlzNextGreyInterval(iwsp);
	      }
	      if (iwsp->linrmn >= 0 && iwsp->linpos == l &&
		  iwsp->lftpos <= k) {
		noverlap++;
		gv.inv += *(gwsp[i].u_grintptr.inp + k - iwsp->lftpos);
	      }
	    }
	    *greyptr.inp = gv.inv / noverlap;
	    greyptr.inp++;
	  }
	  break;

	case WLZ_GREY_SHORT:
	  for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	    noverlap = 0;
	    gv.shv = 0;
	    for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	      while (iwsp->linrmn >= 0
		     && (iwsp->linpos < l
			 || (iwsp->linpos == l && iwsp->rgtpos < k))){
		WlzNextGreyInterval(iwsp);
	      }
	      if (iwsp->linrmn >= 0 && iwsp->linpos == l &&
		  iwsp->lftpos <= k) {
		noverlap++;
		gv.shv += *(gwsp[i].u_grintptr.shp + k - iwsp->lftpos);
	      }
	    }
	    *greyptr.shp = (short )(gv.shv / noverlap);
	    greyptr.shp++;
	  }
	  break;

	case WLZ_GREY_UBYTE:
	  for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	    noverlap = 0;
	    gv.inv = 0;
	    for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	      while (iwsp->linrmn >= 0
		     && (iwsp->linpos < l
			 || (iwsp->linpos == l && iwsp->rgtpos < k))){
		WlzNextGreyInterval(iwsp);
	      }
	      if (iwsp->linrmn >= 0 && iwsp->linpos == l &&
		  iwsp->lftpos <= k) {
		noverlap++;
		gv.inv += *(gwsp[i].u_grintptr.ubp + k - iwsp->lftpos);
	      }
	    }
	    *greyptr.ubp = (WlzUByte )(gv.inv / noverlap);
 	    greyptr.ubp++;
	  }
	  break;

	case WLZ_GREY_FLOAT:
	  for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	    noverlap = 0;
	    gv.flv = 0;
	    for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	      while (iwsp->linrmn >= 0
		     && (iwsp->linpos < l
			 || (iwsp->linpos == l && iwsp->rgtpos < k))){
		WlzNextGreyInterval(iwsp);
	      }
	      if (iwsp->linrmn >= 0 && iwsp->linpos == l &&
		  iwsp->lftpos <= k) {
		noverlap++;
		gv.flv += *(gwsp[i].u_grintptr.flp + k - iwsp->lftpos);
	      }
	    }
	    *greyptr.flp = gv.flv / noverlap;
	    greyptr.flp++;
	  }
	  break;

	case WLZ_GREY_DOUBLE:
	  for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	    noverlap = 0;
	    gv.dbv = 0;
	    for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	      while (iwsp->linrmn >= 0
		     && (iwsp->linpos < l
			 || (iwsp->linpos == l && iwsp->rgtpos < k))){
		WlzNextGreyInterval(iwsp);
	      }
	      if (iwsp->linrmn >= 0 && iwsp->linpos == l &&
		  iwsp->lftpos <= k) {
		noverlap++;
		gv.dbv += *(gwsp[i].u_grintptr.dbp + k - iwsp->lftpos);
	      }
	    }
	    *greyptr.dbp = gv.dbv / noverlap;
	    greyptr.dbp++;
	  }
	  break;
	case WLZ_GREY_RGBA: /* RGBA to be done - do properly RAB */
	  for (k = niwsp.lftpos; k <= niwsp.rgtpos; k++) {
	    noverlap = 0;
	    gv.rgbv = 0;
	    for (iwsp=biwsp,i=0; iwsp<tiwsp; iwsp++,i++) {
	      while (iwsp->linrmn >= 0
		     && (iwsp->linpos < l
			 || (iwsp->linpos == l && iwsp->rgtpos < k))){
		WlzNextGreyInterval(iwsp);
	      }
	      if (iwsp->linrmn >= 0 && iwsp->linpos == l &&
		  iwsp->lftpos <= k) {
		noverlap++;
		gv.rgbv += *(gwsp[i].u_grintptr.rgbp + k - iwsp->lftpos);
	      }
	    }
	    *greyptr.rgbp = gv.rgbv / noverlap;
	    greyptr.rgbp++;
	  }
	  break;

	default:
	  break;

	}
      }
    }
    AlcFree((void *) gwsp);
  }

  if( errNum == WLZ_ERR_NONE ){
    AlcFree( (void *) biwsp);
    AlcFree( (void *)locbuff );
    AlcFree( (void *)locbuffs );
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( obj );
}
