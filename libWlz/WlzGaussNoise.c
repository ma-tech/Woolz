#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGaussNoise.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for making Gaussian noise filled Woolz
*		objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

WlzErrorNum WlzGaussNoise(
  WlzObject	*obj,
  WlzPixelV	val)
{
  WlzObject		*tmpObj;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzPixelV		tmpVal;
  double		mu, sigma;
  WlzDomain		*domains;
  WlzValues		*values;
  int			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.i == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->values.v == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check planedomain and voxeltable */
      if( obj->domain.p == NULL ){
	errNum =  WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if( obj->values.vox == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else {
	/* set range of each plane if non-empty - indicated by NULL */
	domains = obj->domain.p->domains;
	values = obj->values.vox->values;
	nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
	for(i=0; (errNum == WLZ_ERR_NONE) && (i < nplanes);
	    i++, domains++, values++){

	  if( (*domains).core == NULL || (*values).core == NULL ){
	    continue;
	  }

	  if( tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				    NULL, NULL,	NULL) ){
	    errNum = WlzGaussNoise(tmpObj, val);
	    WlzFreeObj( tmpObj );
	  }
	}
      }
      return errNum;

    case WLZ_TRANS_OBJ:
      return WlzGaussNoise(obj->values.obj, val);

    case WLZ_EMPTY_OBJ:
      return errNum;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
    WlzValueConvertPixel(&tmpVal, val, WLZ_GREY_DOUBLE);
    sigma = tmpVal.v.dbv;

    /* set the seed */
    AlgRandSeed((long )obj);

    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	  mu = (double) *gptr.inp;
	  tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	  *gptr.inp = WLZ_CLAMP(tmpVal.v.dbv, INT_MIN, INT_MAX);
	}
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp.colrmn; i++, gptr.shp++){
	  mu = (double) *gptr.shp;
	  tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	  *gptr.shp = WLZ_CLAMP(tmpVal.v.dbv, SHRT_MIN, SHRT_MAX);
	}
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp.colrmn; i++, gptr.ubp++){
	  mu = (double) *gptr.ubp;
	  tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	  *gptr.ubp = WLZ_CLAMP(tmpVal.v.dbv, 0, 255);
	}
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp.colrmn; i++, gptr.flp++){
	  mu = (double) *gptr.flp;
	  tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	  *gptr.flp = WLZ_CLAMP(tmpVal.v.dbv, FLT_MIN, FLT_MAX);
	}
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	  mu = (double) *gptr.inp;
	  tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	  *gptr.inp = tmpVal.v.dbv;
	}
	break;

      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  return errNum;
}
