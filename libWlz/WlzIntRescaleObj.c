#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzIntRescaleObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Rescale a Woolz object using an integral scale.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 03-03-2K bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

static int check_intvs(
  WlzInterval	*intvs,
  int		nintvs)
{
  int	i;

  /* stopping condition */
  if( nintvs == 0 )
    return( 0 );

  /* check tail of the list */
  nintvs = check_intvs(intvs+1, nintvs-1) + 1;

  /* check this interval */
  if( intvs->iright < intvs->ileft ){
    intvs->iright = intvs->ileft;
  }

  /* return if one interval */
  if( nintvs == 1 )
    return( nintvs );

  /* check next interval */
  if( (intvs[1].ileft - intvs[0].iright) < 2 ){ /* merge condition */
    intvs[0].iright = intvs[1].iright;
    for(i=2; i < nintvs; i++)
      intvs[i-1] = intvs[i];
    return( nintvs - 1 );
  }

  return( nintvs );
}

WlzObject *WlzIntRescaleObj(
  WlzObject	*obj,
  int		scale,
  int		expand,
  WlzErrorNum	*dstErr)
{
  WlzObject		*rtnObj=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzInterval		*intvls;
  int			k1, kl, l1, ll, l, num_intvls;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL )
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      /* check the domain */
      if( obj->domain.i == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.i->type == WLZ_EMPTY_DOMAIN ){
	return WlzMakeEmpty(dstErr);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    case WLZ_2D_POLYGON:
    case WLZ_BOUNDLIST:
    case WLZ_3D_POLYGON:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* check the scale and for no change */
  if( errNum == WLZ_ERR_NONE ){
    if( scale < 1 ){
      errNum = WLZ_ERR_PARAM_DATA;
    }
    else if( scale == 1 )
    {
      return WlzMakeMain(obj->type, obj->domain, obj->values,
			  NULL, NULL, dstErr);
    }
  }

  /* check expand or contract */
  if( errNum == WLZ_ERR_NONE ){
    if( expand )
    {
      k1 = obj->domain.i->kol1   * scale;
      kl = obj->domain.i->lastkl * scale + scale - 1;
      l1 = obj->domain.i->line1  * scale;
      ll = obj->domain.i->lastln * scale + scale - 1;
    }
    else {
      k1 = obj->domain.i->kol1   / scale;
      kl = obj->domain.i->lastkl / scale;
      l1 = obj->domain.i->line1  / scale;
      ll = obj->domain.i->lastln / scale;
    }
  }

  /* create a new object */
  if( errNum == WLZ_ERR_NONE ){
    if( domain.i = WlzMakeIntervalDomain(obj->domain.i->type,
					 l1, ll, k1, kl, &errNum) ){
      values.core = NULL;
      rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, 
			   NULL, NULL, NULL);
    }
  }

  /* fill in the intervals */
  if( errNum == WLZ_ERR_NONE){
    if( domain.i->type == WLZ_INTERVALDOMAIN_INTVL )
    {
      int intvline_offset;
      WlzIntervalLine *intvline;

      num_intvls = WlzIntervalCount(obj->domain.i, NULL);
      num_intvls = expand ? num_intvls * scale : num_intvls;
      intvls = (WlzInterval *)AlcMalloc(sizeof(WlzInterval) * num_intvls);
      domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)intvls,
      					   NULL);

      for(l=l1; l <= ll; l++)
      {
	int		i;

	intvline_offset = (expand?l/scale:l*scale) - obj->domain.i->line1;
	intvline_offset = WLZ_MAX(intvline_offset, 0);
	intvline = obj->domain.i->intvlines + intvline_offset;

	for(i=0; i < intvline->nintvs; i++)
	{
	  intvls[i].ileft  = (intvline->intvs + i)->ileft;
	  intvls[i].iright = (intvline->intvs + i)->iright;
	  if( expand )
	  {
	    intvls[i].ileft  *= scale;
	    intvls[i].iright *= scale;
	    intvls[i].iright += scale - 1;
	  }
	  else
	  {
	    intvls[i].ileft  /= scale;
	    intvls[i].iright /= scale;
	  }
	}

	i = check_intvs(intvls, i);
	WlzMakeInterval(l, domain.i, i, intvls);
	intvls += i;
      }
      (void) WlzStandardIntervalDomain( domain.i );
    }
  }

  /* create the valuetable */
  if( (errNum == WLZ_ERR_NONE) && obj->values.core )
  {
    WlzIntervalWSpace	iwsp;
    WlzGreyWSpace	gwsp;
    WlzPixelV 		backgrnd;
    WlzGreyValueWSpace	*gVWSp = NULL;
    WlzGreyType		gtype;

    backgrnd = WlzGetBackground(obj, NULL);
    if( values.v = WlzNewValueTb(rtnObj, obj->values.v->type,
				 backgrnd, &errNum) ){

      rtnObj->values = WlzAssignValues(values, NULL);

      /* fill in the grey-values */
      WlzInitGreyScan(rtnObj, &iwsp, &gwsp);
      gVWSp = WlzGreyValueMakeWSp(obj, NULL);
      gtype = WlzGreyTableTypeToGreyType(obj->values.v->type, NULL);
      while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE )
      {
	int k;
	int lp = expand ? iwsp.linpos/scale : iwsp.linpos*scale;

	for( k=0; k <= (iwsp.rgtpos - iwsp.lftpos); k++ )
	{
	  int	kp = expand ? k/scale : k*scale;

	  WlzGreyValueGet(gVWSp, 0, (double) lp, (double) kp);

	  switch(gtype)
	  {
	  case WLZ_GREY_INT: 
	    gwsp.u_grintptr.inp[k] = (*(gVWSp->gVal)).inv;
	    break;
	  case WLZ_GREY_SHORT:
	    gwsp.u_grintptr.shp[k] = (*(gVWSp->gVal)).shv;
	    break;
	  case WLZ_GREY_UBYTE:
	    gwsp.u_grintptr.ubp[k] = (*(gVWSp->gVal)).ubv;
	    break;
	  case WLZ_GREY_FLOAT:
	    gwsp.u_grintptr.flp[k] = (*(gVWSp->gVal)).flv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    gwsp.u_grintptr.dbp[k] = (*(gVWSp->gVal)).dbv;
	    break;
	  }
	}
      }
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
      WlzGreyValueFreeWSp(gVWSp);
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
