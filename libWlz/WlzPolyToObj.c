#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzFilledPolyToObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Makes a Woolz domain object from a polygon.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 03-03-2K bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>


/* static vertex comparison procedures for qsort */
static int vtx_compare(
  const void	*p1,
  const void	*p2)
{
  WlzIVertex2 *v1 = (WlzIVertex2 *) p1;
  WlzIVertex2 *v2 = (WlzIVertex2 *) p2;

  if( (v1)->vtY > (v2)->vtY ){
    return(1);
  } else if( (v1)->vtY < (v2)->vtY ){
    return(-1);
  } else {
    if( v1->vtX > v2->vtX ){
      return(1);
    } else if( v1->vtX < v2->vtX ){
      return(-1);
    } else {
      return(0);
    }
  }
}

WlzObject *WlzPolygonToObj(
  WlzObject		*polygon,
  WlzPolyFillMode	fillMode,
  WlzErrorNum		*dstErr)
{
  WlzObject	*rtnObj=NULL, *tmpObj;
  WlzDomain	domain, *domains, *polydmns;
  WlzValues	values;
  int		p;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object */
  if( polygon == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( polygon->type ){
    case WLZ_3D_DOMAINOBJ:
      /* check plane domain */
      if( polygon->domain.p == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( polygon->domain.p->type ){
	case WLZ_PLANEDOMAIN_POLYGON:
	  if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					    polygon->domain.p->plane1,
					    polygon->domain.p->lastpl,
					    polygon->domain.p->line1,
					    polygon->domain.p->lastln,
					    polygon->domain.p->kol1,
					    polygon->domain.p->lastkl,
					    &errNum) ){
	    domain.p->voxel_size[0] = polygon->domain.p->voxel_size[0];
	    domain.p->voxel_size[1] = polygon->domain.p->voxel_size[1];
	    domain.p->voxel_size[2] = polygon->domain.p->voxel_size[2];
	    domains = domain.p->domains;
	    polydmns = polygon->domain.p->domains;
	    for(p=domain.p->plane1; p <= domain.p->lastpl;
		p++, domains++, polydmns++){
	      if( (*polydmns).poly ){
		if( tmpObj = WlzPolyToObj((*polydmns).poly, fillMode,
					  &errNum) ){
		  *domains = WlzAssignDomain(tmpObj->domain, NULL);
		  WlzFreeObj(tmpObj);
		}
		else {
		  WlzFreePlaneDomain(domain.p);
		  domain.p = NULL;
		  break;
		}
	      }
	    }
	    if( domain.p ){
	      values.core = NULL;
	      if( (rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
					NULL, NULL, &errNum)) == NULL ){
		WlzFreePlaneDomain(domain.p);
	      }
	    }
	  }
	  break;
	  
	case WLZ_EMPTY_DOMAIN:
	  return WlzMakeEmpty(dstErr);

	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
	}
      }
      break;

    case WLZ_2D_POLYGON:
      return WlzPolyToObj(polygon->domain.poly, fillMode, dstErr);

    case WLZ_TRANS_OBJ:
      if( values.obj = WlzPolygonToObj(polygon->values.obj, fillMode,
				       &errNum) ){
	return WlzMakeMain(WLZ_TRANS_OBJ, polygon->domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

WlzObject *WlzPolyToObj(
  WlzPolygonDomain	*pgdm,
  WlzPolyFillMode	fillMode,
  WlzErrorNum		*dstErr)
{
  /* local variables */
  WlzObject		*new_poly, *obj1, *obj2, *obj3, *obj=NULL;
  WlzObject		**objs;
  WlzDomain		domain;
  WlzValues		values;
  WlzIVertex2		*vtxs;
  WlzInterval		*intptr;
  int			i, j, n, width, height, ignlns;
  int			l1, ll, k1, lk, nints;
  int			num_vtxs, nv;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* find a new polygon eight-connected and closed */
  if( new_poly = WlzPolyTo8Polygon(pgdm, 1, &errNum) ){
    vtxs     = ((WlzPolygonDomain *) new_poly->domain.poly)->vtx;
    num_vtxs = ((WlzPolygonDomain *) new_poly->domain.poly)->nvertices - 1;
  }

  /* find line and column bounds */
  if( errNum == WLZ_ERR_NONE ){
    l1 = ll = vtxs->vtY;
    k1 = lk = vtxs->vtX;
    for(i=1, vtxs++; i < num_vtxs; i++, vtxs++){
      if( vtxs->vtX < k1 ){
	k1 = vtxs->vtX;
      } else if( vtxs->vtX > lk ){
	lk = vtxs->vtX;
      }
      if( vtxs->vtY < l1 ){
	l1 = vtxs->vtY;
      } else if( vtxs->vtY > ll ){
	ll = vtxs->vtY;
      }
    }
    vtxs -= num_vtxs;
    l1--; ll++; k1--; lk++;
    height = ll - l1 + 1;
    width = lk - k1 + 1;

    /* order the vertices first wrt line number then wrt kol number */
    qsort((void *) vtxs, num_vtxs, sizeof(WlzIVertex2), vtx_compare);
  }

  /* build an object with intervals given by non-polynomial points */
  if((errNum == WLZ_ERR_NONE) &&
     (domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				       l1, ll, k1, lk, &errNum)) ){
    values.core = NULL;
    obj1 = WlzAssignObject(WlzMakeMain(WLZ_2D_DOMAINOBJ, domain,
				       values, NULL, NULL, &errNum), NULL);
  }

  /* maximum number of intervals = num_vtxs + height */
  if((errNum == WLZ_ERR_NONE) &&
     (intptr = (WlzInterval *)
       AlcMalloc(sizeof(WlzInterval)*(num_vtxs + height + 1))) ){
    domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)intptr,
    				         NULL);
  }
  else {
    WlzFreeObj(obj1);
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  /* first and last lines - definitely empty */
  if( errNum == WLZ_ERR_NONE ){
    intptr->ileft = 0;
    intptr->iright = lk - k1;
    WlzMakeInterval(l1, domain.i, 1, intptr);	/* first line - empty */
    intptr++;
  }
  if( errNum == WLZ_ERR_NONE ){
    intptr->ileft = 0;
    intptr->iright = lk - k1;
    WlzMakeInterval(ll, domain.i, 1, intptr);	/* last line - empty */
    intptr++;
  }

  /* other lines */
  if( errNum == WLZ_ERR_NONE ){
    for(j=l1+1, nv=0; j < ll; j++){
      n = 0;
      while( (vtxs->vtY == j) && (nv < num_vtxs) ){
	n++;		 /* count vertices on line */
	nv++;
	vtxs++;
      }
      vtxs -= n;

      /* count intervals */
      for(i=1, nints=2; i < n; i++){
	if( ((vtxs+i)->vtX - (vtxs+i-1)->vtX) > 1 ){
	  nints++;
	}
      }
      intptr->ileft = 0;
      intptr++->iright = vtxs++->vtX - 1 - k1;  /* first interval */
      for(i=1; i < n; i++, vtxs++){
	if( (vtxs->vtX - (vtxs-1)->vtX) > 1 ){
	  intptr->ileft = (vtxs-1)->vtX + 1 - k1;
	  intptr++->iright = vtxs->vtX - 1 - k1;
	}
      }
      intptr->ileft = (vtxs-1)->vtX + 1 - k1;
      intptr++->iright = lk - k1;			/* last interval */
      WlzMakeInterval(j, domain.i, nints, intptr-nints);
    }

    /* switch on fill mode */
    switch( fillMode ){

    case WLZ_SIMPLE_FILL:
      /* label to find external object */
      ignlns = (height < width ? height : width) - 1;
      errNum = WlzLabel(obj1, &n, &objs, 1, ignlns, WLZ_4_CONNECTED);
      if( (errNum == WLZ_ERR_NONE) && n ){
	obj2 = *objs;
	AlcFree((void *) objs);
      }
      break;

    case WLZ_EVEN_ODD_FILL:
      WlzFreeObj( obj1 );
      errNum = WLZ_ERR_PARAM_DATA;
      break;

    case WLZ_VERTEX_FILL:
      obj2 = WlzAssignObject(
	WlzMakeMain(obj1->type, obj1->domain, obj1->values,
		    NULL, NULL, &errNum), NULL);
      break;

    }
  }

  /* find residual object */
  if((errNum == WLZ_ERR_NONE) &&
     (domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				       l1, ll, k1, lk, &errNum)) &&
     (obj3 = WlzAssignObject(
       WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL,
		   NULL, &errNum), NULL)) ){
    obj = WlzDiffDomain(obj3, obj2, &errNum);
  }

  /* clean up and return */
  if( errNum == WLZ_ERR_NONE ){
    WlzFreeObj( obj1 );
    WlzFreeObj( obj2 );
    WlzFreeObj( obj3 );
    WlzFreeObj( new_poly );
    if( obj != NULL ){
      switch( obj->type ){

      case WLZ_2D_DOMAINOBJ:
	if( obj->domain.core->type == WLZ_INTERVALDOMAIN_INTVL ){
	  WlzStandardIntervalDomain( obj->domain.i );
	}
	break;

      case WLZ_EMPTY_OBJ:
      default:
	break;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return obj;
}


WlzObject *WlzPolyTo8Polygon(
  WlzPolygonDomain	*pgdm,
  int			wrap,
  WlzErrorNum		*dstErr)
{
  /* local variables */
  WlzObject		*obj=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzPolygonDomain	*npgdm;
  WlzIVertex2		*vtxs, *nvtxs;
  void			*freeptr=NULL;
  int			i, j, n, length, k;
  int			x, y, x0, y0, lx, ly;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the domain and wrap */
  if( pgdm == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if( wrap < 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  /* check type nad copy non-integer vertices */
  if( errNum == WLZ_ERR_NONE ){
    domain.core = NULL;
    values.core = NULL;
    n = pgdm->nvertices;
    switch( pgdm->type ){

    case WLZ_EMPTY_OBJ:
      return WlzMakeMain(WLZ_EMPTY_OBJ, domain, values, NULL, NULL, dstErr);

    case WLZ_POLYGON_INT:
      vtxs = pgdm->vtx;
      break;

    case WLZ_POLYGON_FLOAT:
      if( (freeptr = AlcMalloc(sizeof(WlzIVertex2)*n)) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      vtxs = (WlzIVertex2 *) freeptr;
      WlzValueCopyFVertexToIVertex(vtxs, (WlzFVertex2 *) pgdm->vtx, n);
      break;

    case WLZ_POLYGON_DOUBLE:
      if( (freeptr = AlcMalloc(sizeof(WlzIVertex2)*n)) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      vtxs = (WlzIVertex2 *) freeptr;
      WlzValueCopyDVertexToIVertex(vtxs, (WlzDVertex2 *) pgdm->vtx, n);
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* check eight connected length */
  if( errNum == WLZ_ERR_NONE ){
    x0 = vtxs->vtX;
    y0 = vtxs->vtY;
    for(i=1, vtxs++, length=0; i < n; i++, vtxs++){
      x = vtxs->vtX;
      y = vtxs->vtY;
      lx = abs(x - x0);
      ly = abs(y - y0);
      length += (lx > ly) ? lx : ly;
      x0 = x; y0 = y;
    }
    vtxs -= n;
    x = vtxs->vtX; y = vtxs->vtY;
    lx = abs(x - x0); ly = abs(y - y0);
    length += (lx > ly) ? lx : ly;
    x0 = x; y0 = y;
  }

  /* make an eight-connected polyline */
  if( errNum == WLZ_ERR_NONE ){
    if( (npgdm = WlzMakePolyDmn(WLZ_POLYGON_INT,
				NULL, 0, length+10+wrap,
				1, &errNum)) == NULL ){
      if( freeptr ){
	AlcFree( freeptr );
      }
    }
    else {
      nvtxs = npgdm->vtx;
      nvtxs->vtX = x;
      nvtxs->vtY = y;
      nvtxs++;
      k=1;
      for(j=1, vtxs++; j < n; j++, vtxs++){
	x = vtxs->vtX;
	y = vtxs->vtY;
	lx = x - x0;
	ly = y - y0;
	length = (abs(lx) > abs(ly)) ? abs(lx) : abs(ly);
	for(i=1; i <= length; i++, nvtxs++, k++){
	  nvtxs->vtX = x0 + i * lx / length;
	  nvtxs->vtY = y0 + i * ly / length;
	}
	x0 = x;
	y0 = y;
      }
      vtxs -= n;
    }
  }

  /* close the polyline */
  if( errNum == WLZ_ERR_NONE ){
    if( wrap > 0 ){
      x = vtxs->vtX;
      y = vtxs->vtY;
      lx = x - x0;
      ly = y - y0;
      length = (abs(lx) > abs(ly)) ? abs(lx) : abs(ly);
      for(i=1; i <= length; i++, nvtxs++, k++){
	nvtxs->vtX = x0 + i * lx / length;
	nvtxs->vtY = y0 + i * ly / length;
      }
    }
    for(i=1; i < wrap; i++, nvtxs++, k++){
      *nvtxs = npgdm->vtx[i];
    }

    npgdm->nvertices = k;
    if( freeptr ){
      AlcFree( freeptr );
    }

    domain.poly = npgdm;
    obj = WlzMakeMain(WLZ_2D_POLYGON, domain, values, NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return obj;
}


