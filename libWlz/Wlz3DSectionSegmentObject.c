#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        Wlz3DSectionSegmentObject.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Segments a sectionthrough a 3D object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 03-03-2K bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
************************************************************************/
#include <Wlz.h>

WlzErrorNum Wlz3DSectionSegmentObject(
  WlzObject		*obj,
  WlzThreeDViewStruct	*viewStr,
  int			*numObjs,
  WlzObject		***rtnObjs)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*tmpObj, *obj0, *obj1, **objs;
  WlzDomain	domain;
  WlzValues	values;
  WlzIntervalDomain	*idom0, *idom1;
  WlzInterval		*intvls0, *intvls1;
  int			nints0, nints1;
  WlzIntervalWSpace	iwsp;
  WlzDVertex3	vtxL, vtxR;
  double	fdelta;
  int		i, j, p, pp;
  int		delta;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_3D_DOMAINOBJ:
      switch( obj->domain.core->type ){
      case WLZ_PLANEDOMAIN_DOMAIN:
	break;

      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      }
      break;

    case WLZ_EMPTY_OBJ:
      if( objs = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * 2) ){
	objs[0] = WlzAssignObject(WlzMakeEmpty(&errNum), NULL);
	objs[1] = WlzAssignObject(WlzMakeEmpty(&errNum), NULL);
	*rtnObjs = objs;
	*numObjs = 2;
      }
      else {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      return errNum;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* check the view structure */
  if( errNum == WLZ_ERR_NONE ){
    if( viewStr == NULL ){
      errNum = WLZ_ERR_VALUES_NULL;
    }
    else {
      if( viewStr->initialised == 0 ){
	errNum = WLZ_ERR_VALUES_DATA;
      }
    }
  }

  /* allocate space and set number of objects */
  if( objs = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * 2) ){
    *rtnObjs = objs;
    *numObjs = 2;
  }
  else {
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  /* create a new 3d objects for each side */
  if( errNum == WLZ_ERR_NONE ){
    values.core = NULL;
    if(domain.p = 
       WlzMakePlaneDomain(obj->domain.p->type,
			  obj->domain.p->plane1, obj->domain.p->lastpl,
			  obj->domain.p->line1, obj->domain.p->lastln,
			  obj->domain.p->kol1, obj->domain.p->lastkl,
			  &errNum) ){
      if( (objs[0] = WlzMakeMain(obj->type, domain, values,
				 NULL, NULL, &errNum)) == NULL ){
	WlzFreePlaneDomain(domain.p);
	AlcFree((void *) objs);
      }
      else {
	objs[0] = WlzAssignObject(objs[0], &errNum);
      }
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    values.core = NULL;
    if(domain.p = 
       WlzMakePlaneDomain(obj->domain.p->type,
			  obj->domain.p->plane1, obj->domain.p->lastpl,
			  obj->domain.p->line1, obj->domain.p->lastln,
			  obj->domain.p->kol1, obj->domain.p->lastkl,
			  &errNum) ){
      if( (objs[1] = WlzMakeMain(obj->type, domain, values,
				 NULL, NULL, &errNum)) == NULL ){
	WlzFreePlaneDomain(domain.p);
	WlzFreeObj(objs[0]);
	AlcFree((void *) objs);
      }
      else {
	objs[1] = WlzAssignObject(objs[1], &errNum);
      }
    }
    else {
      WlzFreeObj(objs[0]);
      AlcFree((void *) objs);
    }
  }


  /* place intervals as required in -ve and +ve object */
  if( errNum == WLZ_ERR_NONE ){
    for(p=obj->domain.p->plane1, pp=0; p <= obj->domain.p->lastpl;
	p++, pp++){
      if((obj->domain.p->domains[pp].core == NULL) ||
	 (obj->domain.p->domains[pp].core->type == WLZ_EMPTY_OBJ) ||
	 (obj->domain.p->domains[pp].core->type == WLZ_EMPTY_DOMAIN)){
	objs[0]->domain.p->domains[pp].core = NULL;
	objs[1]->domain.p->domains[pp].core = NULL;
	continue;
      }

      /* temporary plane object from the original */
      tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj->domain.p->domains[pp],
			   values, NULL, NULL, NULL);

      /* make interval domains and intervals for the two new objects */
      nints0 = WlzIntervalCount(tmpObj->domain.i, NULL);
      idom0 = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				    tmpObj->domain.i->line1,
				    tmpObj->domain.i->lastln,
				    tmpObj->domain.i->kol1,
				    tmpObj->domain.i->lastkl,
				    NULL);
      intvls0 = AlcMalloc(nints0*sizeof(WlzInterval));
      idom0->freeptr = AlcFreeStackPush(idom0->freeptr,
				        (void *)intvls0, NULL);

      idom1 = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				    tmpObj->domain.i->line1,
				    tmpObj->domain.i->lastln,
				    tmpObj->domain.i->kol1,
				    tmpObj->domain.i->lastkl,
				    NULL);
      intvls1 = AlcMalloc(nints0*sizeof(WlzInterval));
      idom1->freeptr = AlcFreeStackPush(idom1->freeptr,
				        (void *)intvls1, NULL);

      /* scan the original object adding/splitting intervals as required */
      nints0 = 0;
      nints1 = 0;
      WlzInitRasterScan(tmpObj, &iwsp, WLZ_RASTERDIR_ILIC);
      while( (errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){

	/* check interval, mark as required */
	vtxL.vtX = iwsp.colpos + 0.05;
	vtxL.vtY = iwsp.linpos + 0.05;
	vtxL.vtZ = p + 0.05;
	Wlz3DSectionTransformVtx(&vtxL, viewStr);
	vtxL.vtZ -= viewStr->dist;

	vtxR.vtX = iwsp.colpos + iwsp.colrmn - 1 + 0.05;
	vtxR.vtY = iwsp.linpos + 0.05;
	vtxR.vtZ = p + 0.05;
	Wlz3DSectionTransformVtx(&vtxR, viewStr);
	vtxR.vtZ -= viewStr->dist;

	if( vtxL.vtZ < 0.0 ){
	  if( vtxR.vtZ < 0.0 ){
	    /* pure negative object only */
	    intvls0[nints0++] = *iwsp.intpos;
	  }
	  else {
	    /* mixed */
	    fdelta = vtxR.vtZ / (vtxR.vtZ - vtxL.vtZ) * (iwsp.colrmn-1);
	    delta = WLZ_NINT(fdelta);
	    intvls0[nints0].ileft = iwsp.intpos->ileft;
	    intvls0[nints0].iright = iwsp.intpos->iright - delta;
	    if( delta ){
	      intvls1[nints1].ileft = intvls0[nints0].iright + 1;
	      intvls1[nints1].iright = iwsp.intpos->iright;
	      nints1++;
	    }
	    nints0++;
	  }
	}
	else {
	  if( vtxR.vtZ < 0.0 ){
	    /* mixed */
	    fdelta = vtxR.vtZ / (vtxR.vtZ - vtxL.vtZ) * (iwsp.colrmn-1);
	    delta = WLZ_NINT(fdelta);
	    intvls1[nints1].ileft = iwsp.intpos->ileft;
	    intvls1[nints1].iright = iwsp.intpos->iright - delta;
	    if( delta ){
	      intvls0[nints0].ileft = intvls1[nints1].iright + 1;
	      intvls0[nints0].iright = iwsp.intpos->iright;
	      nints0++;
	    }
	    nints1++;
	  }
	  else {
	    /* pure positive object only */
	    intvls1[nints1++] = *iwsp.intpos;
	  }
	}

	/* if it is the end of line then add the intervals to
	   their domains */
	if( iwsp.intrmn == 0 ){
	  if( nints0 ){
	    WlzMakeInterval(iwsp.linpos, idom0, nints0, intvls0);
	    intvls0 += nints0;
	    nints0 = 0;
	  }
	  if( nints1 ){
	    WlzMakeInterval(iwsp.linpos, idom1, nints1, intvls1);
	    intvls1 += nints1;
	    nints1 = 0;
	  }
	}

      }
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
      WlzFreeObj(tmpObj);

      /* standardise the domains, check for empty planes */
      WlzStandardIntervalDomain(idom0);
      if( WlzIntervalCount(idom0, NULL) ){
	domain.i = idom0;
	objs[0]->domain.p->domains[pp] =
	  WlzAssignDomain(domain, NULL);
      }
      else {
	WlzFreeIntervalDomain(idom0);
	objs[0]->domain.p->domains[pp].core = NULL;
      }

      WlzStandardIntervalDomain(idom1);
      if( WlzIntervalCount(idom1, NULL) ){
	domain.i = idom1;
	objs[1]->domain.p->domains[pp] =
	  WlzAssignDomain(domain, NULL);
      }
      else {
	WlzFreeIntervalDomain(idom1);
	objs[1]->domain.p->domains[pp].core = NULL;
      }
      
    }
  }

  /* standardise the planedomains */
  WlzStandardPlaneDomain(objs[0]->domain.p, NULL);
  WlzStandardPlaneDomain(objs[1]->domain.p, NULL);

  /* free empty objects */
  if( errNum == WLZ_ERR_NONE ){
    if( WlzIsEmpty(objs[0], NULL) ){
      WlzFreeObj(objs[0]);
      objs[0] = NULL;
    }
    if( WlzIsEmpty(objs[1], NULL) ){
      WlzFreeObj(objs[1]);
      objs[1] = NULL;
    }
  }

  return errNum;
}
