#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzThreshold.c
* \author       Richard Baldock
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
* \brief	Thresholds a Woolz grey-level object, 2D or 3D.
* \ingroup	WlzThreshold
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <Wlz.h>

static WlzObject *WlzThreshold3d(WlzObject	*obj,
				 WlzPixelV	threshV,
				 WlzThresholdType highlow,
				 WlzErrorNum	*dstErr);

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzThreshold
* \brief	Thresholds a woolz grey-level object, 2D or 3D.
* \param	obj			Object to be thresholded.
* \param	threshV			Threshold pixel value.
* \param	highlow			Mode parameter with possible values:
*					<ul>
*					<li> WLZ_THRESH_HIGH - thresholded
*					object is of values >= threshold value.
*					</li>
*					<li> WLZ_THRESH_LOW - thresholded
*					object is of values < threshold value.
*					</li>
*					</ul>
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzObject *WlzThreshold(WlzObject	*obj,
			WlzPixelV	threshV,
			WlzThresholdType highlow,
			WlzErrorNum	*dstErr)
{
  WlzObject		*nobj=NULL;
  WlzIntervalDomain	*idom;
  WlzGreyP		g;
  int			colno, nints;
  int			over;
  int			nl1,nll,nk1,nkl;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzInterval		*itvl, *jtvl;
  int			thresh_i;
  float			thresh_f;
  double		thresh_d;
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
      return WlzThreshold3d(obj, threshV, highlow, dstErr);

    case WLZ_TRANS_OBJ:
      if( nobj = WlzThreshold(obj->values.obj, threshV, highlow,
			      &errNum) ){
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

  /* check the highlow flag */
  if( errNum == WLZ_ERR_NONE ){
    switch( highlow ){

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;

    case WLZ_THRESH_HIGH:
    case WLZ_THRESH_LOW:
      break;
    }
  }
    
  /* get the threshold value - this does not need to be the same
     as the valuetable pixel type */
  if( errNum == WLZ_ERR_NONE ){
    switch( threshV.type ){
    case WLZ_GREY_INT:
      thresh_i = thresh_f = thresh_d = threshV.v.inv;
      break;
    case WLZ_GREY_SHORT:
      thresh_i = thresh_f = thresh_d = (int) threshV.v.shv;
      break;
    case WLZ_GREY_UBYTE:
      thresh_i = thresh_f = thresh_d = (int) threshV.v.ubv;
      break;
    case WLZ_GREY_FLOAT:
      thresh_i = thresh_f = thresh_d = threshV.v.flv;
      break;
    case WLZ_GREY_DOUBLE:
      thresh_i = thresh_f = thresh_d = threshV.v.dbv;
      break;
    case WLZ_GREY_RGBA:
      thresh_i = WLZ_RGBA_MODULUS(threshV.v.rgbv);
      thresh_f = thresh_d = thresh_i;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
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
    nints = 0;
    if( gwsp.pixeltype == WLZ_GREY_RGBA ){
      thresh_i *= thresh_i;
    }
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
      g = gwsp.u_grintptr;
      over = 0;
      switch( gwsp.pixeltype ){

      case WLZ_GREY_INT:
	for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	  if((highlow == WLZ_THRESH_HIGH && *g.inp >= thresh_i) ||
	     (highlow == WLZ_THRESH_LOW && *g.inp < thresh_i) ){
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
	  g.inp++;
	}
	break;

      case WLZ_GREY_SHORT:
	for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	  if((highlow == WLZ_THRESH_HIGH && *g.shp >= thresh_i) ||
	     (highlow == WLZ_THRESH_LOW && *g.shp < thresh_i) ){
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
	  g.shp++;
	}
	break;

      case WLZ_GREY_UBYTE:
	for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	  if((highlow == WLZ_THRESH_HIGH && ((int) (*g.ubp)) >= thresh_i) ||
	     (highlow == WLZ_THRESH_LOW && ((int) (*g.ubp)) < thresh_i) ){
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
	  g.ubp++;
	}
	break;

      case WLZ_GREY_FLOAT:
	for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	  if((highlow == WLZ_THRESH_HIGH && (*g.flp >= thresh_f)) ||
	     (highlow == WLZ_THRESH_LOW && (*g.flp < thresh_f)) ){
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
	  g.flp++;
	}
	break;

      case WLZ_GREY_DOUBLE:
	for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	  if((highlow == WLZ_THRESH_HIGH && (*g.dbp >= thresh_d)) ||
	     (highlow == WLZ_THRESH_LOW && (*g.dbp < thresh_d)) ){
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
	  g.dbp++;
	}
	break;

      case WLZ_GREY_RGBA: /* what to do - OR or AND ? choose MODULUS */
	for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	  if(((highlow == WLZ_THRESH_HIGH) && 
	      (WLZ_RGBA_MODULUS_2(*g.rgbp) >= thresh_i)) ||
	     ((highlow == WLZ_THRESH_LOW) && 
	      (WLZ_RGBA_MODULUS_2(*g.rgbp) < thresh_i)) ){
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
	  g.rgbp++;
	}
	break;

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
      if( idom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				       nl1, nll, nk1, nkl, &errNum) ){
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
	nints = 0;
	jtvl = itvl;
      }

      /* find thresholded endpoints */
      if( errNum == WLZ_ERR_NONE ){
	while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	  if( iwsp.linpos < nl1 || iwsp.linpos > nll ){
	    continue;
	  }
	  g = gwsp.u_grintptr;
	  over = 0;
	  switch( gwsp.pixeltype ){

	  case WLZ_GREY_INT:
	    for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	      if ((highlow == WLZ_THRESH_HIGH && *g.inp >= thresh_i) ||
		  (highlow == WLZ_THRESH_LOW && *g.inp < thresh_i)) {
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
	      g.inp++;
	    }
	    break;

	  case WLZ_GREY_SHORT:
	    for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	      if ((highlow == WLZ_THRESH_HIGH && *g.shp >= thresh_i) ||
		  (highlow == WLZ_THRESH_LOW && *g.shp < thresh_i)) {
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
	      g.shp++;
	    }
	    break;

	  case WLZ_GREY_UBYTE:
	    for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	      if((highlow == WLZ_THRESH_HIGH &&
		  ((int) (*g.ubp)) >= thresh_i) ||
		 (highlow == WLZ_THRESH_LOW &&
		  ((int) (*g.ubp)) < thresh_i)) {
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
	      g.ubp++;
	    }
	    break;

	  case WLZ_GREY_FLOAT:
	    for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	      if ((highlow == WLZ_THRESH_HIGH && *g.flp >= thresh_f) ||
		  (highlow == WLZ_THRESH_LOW && *g.flp < thresh_f)) {
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
	      g.flp++;
	    }
	    break;

	  case WLZ_GREY_DOUBLE:
	    for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	      if ((highlow == WLZ_THRESH_HIGH && *g.dbp >= thresh_d) ||
		  (highlow == WLZ_THRESH_LOW && *g.dbp < thresh_d)) {
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
	      g.dbp++;
	    }
	    break;

	  case WLZ_GREY_RGBA: /* what to do - OR or AND ? choose MODULUS */
	    for (colno = iwsp.lftpos; colno <= iwsp.rgtpos; colno++) {
	      if(((highlow == WLZ_THRESH_HIGH) && 
		  (WLZ_RGBA_MODULUS_2(*g.rgbp) >= thresh_i)) ||
		 ((highlow == WLZ_THRESH_LOW) && 
		  (WLZ_RGBA_MODULUS_2(*g.rgbp) < thresh_i)) ){
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
	      g.rgbp++;
	    }
	    break;

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

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzThreshold
* \brief	Private function used to threshold 3D objects.
* \param	obj			Object to be thresholded.
* \param	threshV			Threshold pixel value.
* \param	highlow			Mode parameter with possible values:
*					<ul>
*					<li> WLZ_THRESH_HIGH - thresholded
*					object is of values >= threshold value.
*					</li>
*					<li> WLZ_THRESH_LOW - thresholded
*					object is of values < threshold value.
*					</li>
*					</ul>
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
static WlzObject *WlzThreshold3d(WlzObject	*obj,
				 WlzPixelV	threshV,
				 WlzThresholdType highlow,
				 WlzErrorNum	*dstErr)
{
  /*----LOCAL AUTOMATIC VARIABLES-----*/
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

  /* Threshold each plane */
#ifdef _OPENMP
  {
    int		pnIdx;
    WlzErrorNum pvErrNum;
    /* (*(ndomains + pnIdx)).core and (*(nvalues + pnIdx)).core are
     * initialized to NULL by WlzMakePlaneDomain() and WlzMakeVoxelValueTb()
     * when they are created. */
    #pragma omp parallel for default(shared) private(pvErrNum,temp,obj1)
    for(pnIdx = 0; pnIdx < nplanes; ++pnIdx)
    {
      if((errNum == WLZ_ERR_NONE) &&
         (*(domains + pnIdx)).core && (*(values + pnIdx)).core)
      {
	temp = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			   *(domains + pnIdx), *(values + pnIdx),
			   NULL, NULL, &pvErrNum);
	obj1 = NULL;
	if(temp && temp->domain.core)
	{
	  obj1 = WlzThreshold(temp, threshV, highlow, &pvErrNum);
	}
	if(obj1)
	{
	  if( obj1->type == WLZ_2D_DOMAINOBJ ){
	    *(ndomains + pnIdx) = WlzAssignDomain(obj1->domain, NULL);
	    *(nvalues + pnIdx) = WlzAssignValues(obj1->values, NULL);
	  }
	  else {
	    (*(ndomains + pnIdx)).core = NULL;
	    (*(nvalues + pnIdx).core = NULL;
	  WlzFreeObj(obj1);
	}
	#pragma omp critical
	{
	  if((pvErrNum != WLZ_ERR_NONE) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = pvErrNum;
	  }
	}
      }
    }
  }
#else
  while( (errNum == WLZ_ERR_NONE) && nplanes-- ){
    if(((*domains).core == NULL) || ((*values).core == NULL)){
      (*ndomains).core = NULL;
      (*nvalues).core = NULL;
    }
    else if( temp = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				NULL, NULL, &errNum) ){

      if( temp->domain.i != NULL ){
	if( obj1 = WlzThreshold(temp, threshV, highlow, &errNum) ){
	  if( obj1->type == WLZ_2D_DOMAINOBJ ){
	    *ndomains = WlzAssignDomain(obj1->domain, NULL);
	    *nvalues = WlzAssignValues(obj1->values, NULL);
	  }
	  else {
	    (*ndomains).core = NULL;
	    (*nvalues).core = NULL;
	  }
	  WlzFreeObj(obj1);
	}
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
#endif
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
    if( obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, vals,
			   NULL, obj, &errNum) ){
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
