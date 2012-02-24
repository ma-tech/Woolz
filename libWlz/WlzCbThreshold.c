#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCbThreshold_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCbThreshold.c
* \author       Richard Baldock
* \date         August 2003
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
* \brief	Threshold an object using a callback function
* 		to determine if the pixel value satifies the threshold
* 		criteria.
* \ingroup	WlzThreshold
*/

#include <stdlib.h>
#include <Wlz.h>

static WlzObject *WlzCbThreshold3d(WlzObject	*obj,
				   WlzThreshCbFn	threshCb,
				   void		*clientData,
				   WlzErrorNum	*dstErr);

WlzObject *WlzCbThreshold(
  WlzObject	*obj,
  WlzThreshCbFn	threshCb,
  void		*clientData,
  WlzErrorNum	*dstErr)
{
  WlzObject		*nobj=NULL;
  WlzIntervalDomain	*idom;
  WlzGreyP		g;
  WlzThreshCbStr	callData;
  int			colno, nints;
  int			over;
  int			nl1,nll,nk1,nkl;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzInterval		*itvl, *jtvl;
  WlzDomain		domain;
  WlzValues		values;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      /* check object 2D domain and valuetable */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzCbThreshold3d(obj, threshCb, clientData, dstErr);

    case WLZ_TRANS_OBJ:
      if((nobj = WlzCbThreshold(obj->values.obj, threshCb, clientData,
				&errNum)) != NULL){
	values.obj = nobj;
	return WlzMakeMain(obj->type, obj->domain, values,
			   NULL, obj, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }

  /*
   * first pass - find line and column bounds of thresholded
   * object and number of intervals.
   */
  if( errNum == WLZ_ERR_NONE ){
    idom = obj->domain.i;
    nl1 = idom->lastln;
    nll = idom->line1;
    nk1 = idom->lastkl;
    nkl = idom->kol1;
    (void) WlzInitGreyScan(obj, &iwsp, &gwsp);
    callData.pix.type = gwsp.pixeltype;
    nints = 0;
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
      callData.pos.vtY = iwsp.linpos;
      g = gwsp.u_grintptr;
      over = 0;

      for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	callData.pix.p = g;
	callData.pos.vtX = colno;
	if((*threshCb)(obj, clientData, &callData) ){
	  if (over == 0) {
	    over = 1;
	    if (iwsp.linpos < nl1)
	      nl1 = iwsp.linpos;
	    if (iwsp.linpos > nll)
	      nll = iwsp.linpos;
	    if (colno < nk1)
	      nk1 = colno;
	  }
	} else {
	  if (over == 1) {
	    if (colno > nkl)
	      nkl = colno;
	    over = 0;
	    nints++;
	  }
	}

	switch( gwsp.pixeltype ){
	case WLZ_GREY_INT:
	  g.inp++;
	  break;
	case WLZ_GREY_SHORT:
	  g.shp++;
	  break;
	case WLZ_GREY_UBYTE:
	  g.ubp++;
	  break;
	case WLZ_GREY_FLOAT:
	  g.flp++;
	  break;
	case WLZ_GREY_DOUBLE:
	  g.dbp++;
	  break;
	case WLZ_GREY_RGBA:
	  g.rgbp++;
	  break;
	default:
	  break;
	}
      }

      if (over == 1) {
	if (colno > nkl)
	  nkl = colno;
	over = 0;
	nints++;
      }
    }
    nkl--;	/* since we have looked at points beyond interval ends */
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  /* domain structure */
  if( errNum == WLZ_ERR_NONE ){
    if( nints > 0 ){
      if((idom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				       nl1, nll, nk1, nkl, &errNum)) != NULL){
	if( (itvl = (WlzInterval *)
	     AlcMalloc(nints * sizeof(WlzInterval))) == NULL ){
	  errNum = WLZ_ERR_MEM_ALLOC;
	  WlzFreeIntervalDomain(idom);
	}
	else {
	  idom->freeptr = AlcFreeStackPush(idom->freeptr, (void *)itvl, NULL);
	}
      }

      /*
       * second pass - construct intervals
       */
      if( errNum == WLZ_ERR_NONE ){
	errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
	callData.pix.type = gwsp.pixeltype;
	nints = 0;
	jtvl = itvl;
      }

      /* find thresholded endpoints */
      if( errNum == WLZ_ERR_NONE ){
	while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	  if( iwsp.linpos < nl1 || iwsp.linpos > nll ){
	    continue;
	  }
	  callData.pos.vtY = iwsp.linpos;
	  g = gwsp.u_grintptr;
	  over = 0;

	  for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	    callData.pix.p = g;
	    callData.pos.vtX = colno;
	    if((*threshCb)(obj, clientData, &callData) ){
	      if (over == 0) {
		over = 1;
		itvl->ileft = colno - nk1;
	      }
	    } else {
	      if (over == 1) {
		over = 0;
		itvl->iright = colno - nk1 - 1;
		nints++;
		itvl++;
	      }
	    }

	    switch( gwsp.pixeltype ){
	    case WLZ_GREY_INT:
	      g.inp++;
	      break;
	    case WLZ_GREY_SHORT:
	      g.shp++;
	      break;
	    case WLZ_GREY_UBYTE:
	      g.ubp++;
	      break;
	    case WLZ_GREY_FLOAT:
	      g.flp++;
	      break;
	    case WLZ_GREY_DOUBLE:
	      g.dbp++;
	      break;
	    case WLZ_GREY_RGBA:
	      g.rgbp++;
	      break;
	    default:
	      break;
	    }
	  }

	  if (over == 1) {
	    over = 0;
	    itvl->iright = colno - nk1 - 1;
	    nints++;
	    itvl++;
	  }
	  /*
	   * end of line ?
	   */
	  if (iwsp.intrmn == 0) {
	    WlzMakeInterval(iwsp.linpos, idom, nints, jtvl);
	    jtvl = itvl;
	    nints = 0;
	  }
	}
	if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	}
      }
    } else {
      /* no thresholded points - make a dummy domain anyway */
      return WlzMakeEmpty(dstErr);
    }
  }

  /* main object */
  if( errNum == WLZ_ERR_NONE ){
    domain.i = idom;
    nobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, obj->values,
		       obj->plist, obj, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(nobj);
}

static WlzObject *WlzCbThreshold3d(
  WlzObject	*obj,
  WlzThreshCbFn	threshCb,
  void		*clientData,
  WlzErrorNum	*dstErr)
{
  WlzObject		*obj1, *temp;
  WlzPlaneDomain	*pdom, *npdom;
  WlzVoxelValues	*voxtab, *nvoxtab;
  WlzDomain		*domains, *ndomains, domain;
  WlzValues		*values, *nvalues, vals;
  int			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* no need to check the object pointer or type because this procedure
     can only be accessed via WlzThreshold. The domain and valuetable
     must be checked however */
  obj1 = NULL;
  if( obj->domain.p == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->values.vox == NULL ){
    errNum = WLZ_ERR_VALUES_NULL;
  }

  /* check types */
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->domain.p->type ){

    case WLZ_PLANEDOMAIN_DOMAIN:
      break;

    default:
      errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
      break;
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->values.vox->type ){

    case WLZ_VOXELVALUETABLE_GREY:
      break;

    default:
      errNum = WLZ_ERR_VOXELVALUES_TYPE;
      break;
    }
  }

  /* make new planedomain and voxelvaluetable */
  if( errNum == WLZ_ERR_NONE ){
    pdom = obj->domain.p;
    voxtab = obj->values.vox;
    npdom = WlzMakePlaneDomain(pdom->type,
			       pdom->plane1, pdom->lastpl,
			       pdom->line1, pdom->lastln,
			       pdom->kol1, pdom->lastkl, &errNum);
  }
    
  if((errNum == WLZ_ERR_NONE) &&
     ((nvoxtab = WlzMakeVoxelValueTb(voxtab->type, voxtab->plane1,
				     voxtab->lastpl, voxtab->bckgrnd,
				     NULL, &errNum)) == NULL) ){
    WlzFreePlaneDomain(npdom);
  }

  if( errNum == WLZ_ERR_NONE ){
    /* set up variables */
    domains = pdom->domains;
    ndomains = npdom->domains;
    values = voxtab->values;
    nvalues = nvoxtab->values;
    nplanes = pdom->lastpl - pdom->plane1 + 1;

    /* copy voxel_sizes */
    for(i=0; i < 3; i++){
      npdom->voxel_size[i] = pdom->voxel_size[i];
    }
  }

  /* threshold each plane */
  while( (errNum == WLZ_ERR_NONE) && nplanes-- ){
    if(((*domains).core == NULL) || ((*values).core == NULL)){
      (*ndomains).core = NULL;
      (*nvalues).core = NULL;
    }
    else if((temp = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				NULL, NULL, &errNum)) != NULL){
      
      if( (temp->domain.i != NULL) && 
	 (obj1 = WlzCbThreshold(temp, threshCb, clientData, &errNum)) ){
	*ndomains = WlzAssignDomain(obj1->domain, NULL);
	*nvalues = WlzAssignValues(obj1->values, NULL);
	WlzFreeObj(obj1);
      } else {
	(*ndomains).core = NULL;
	(*nvalues).core = NULL;
      }
    }
    else {
      WlzFreePlaneDomain(npdom);
      WlzFreeVoxelValueTb( nvoxtab );
      break;
    }

    domains++;
    ndomains++;
    values++;
    nvalues++;
  }

  /* standardise the plane domain */
  if((errNum == WLZ_ERR_NONE) &&
     ((errNum = WlzStandardPlaneDomain(npdom, nvoxtab)) != WLZ_ERR_NONE) ){
    WlzFreePlaneDomain( npdom );
    WlzFreeVoxelValueTb( nvoxtab );
  }

  /* return a new object */
  if( errNum == WLZ_ERR_NONE ){
    domain.p = npdom;
    vals.vox = nvoxtab;
    if((obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, vals,
			   NULL, obj, &errNum)) != NULL){
      /*	nvoxtab->original = obj1; */
      nvoxtab->original_table = WlzAssignValues(obj->values, NULL);
    }
    else {
      WlzFreePlaneDomain( npdom );
      WlzFreeVoxelValueTb( nvoxtab );
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return obj1;
}
