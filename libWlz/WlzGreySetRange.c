#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreySetRange.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Set new grey range for a Woolz object either using
*		simple linear interpolation or nearest value.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzGreySetRange					*
*   Date       : Wed Sep 17 12:00:41 1997				*
*************************************************************************
*   Synopsis   :Set new grey-range by simple linear interpolation. It	*
*		assumed that the min and max values enclose the true	*
*		range of grey-values in the object. Failure to check	*
*		this could result in a segmentation fault.		*
*		The transform function is:				*
*		g' = ((gMax - gMin)/(gmax-gmin))*(g - gmin) + gMin	*
*   Returns    :							*
*   Parameters :							*
*   Global refs:							*
************************************************************************/

WlzErrorNum WlzGreySetRange(
  WlzObject	*obj,
  WlzPixelV	min,
  WlzPixelV	max,
  WlzPixelV	Min,
  WlzPixelV	Max)
{
  double		factor;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzObject		*tempobj;
  WlzValues 		*values;
  WlzDomain		*domains;
  int			i, nplanes;
  WlzErrorNum		wlzErrno=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    wlzErrno = WLZ_ERR_OBJECT_NULL;
  }

  if( wlzErrno == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.i == NULL ){
	return WLZ_ERR_DOMAIN_NULL;
      }
      if( obj->values.v == NULL ){
	return WLZ_ERR_VALUES_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check planedomain and voxeltable */
      if( obj->domain.p == NULL ){
	return WLZ_ERR_DOMAIN_NULL;
      }
      if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	return WLZ_ERR_PLANEDOMAIN_TYPE;
      }
      if( obj->values.vox == NULL ){
	return WLZ_ERR_VALUES_NULL;
      }
      if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	return WLZ_ERR_VOXELVALUES_TYPE;
      }

      /* set range of each plane if non-empty - indicated by NULL */
      domains = obj->domain.p->domains;
      values = obj->values.vox->values;
      nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
      for(i=0; i < nplanes; i++, domains++, values++){

	if( (*domains).core == NULL || (*values).core == NULL ){
	  continue;
	}

	tempobj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			      *domains, *values, NULL, NULL,
			      &wlzErrno);
	if((tempobj == NULL) && (wlzErrno == WLZ_ERR_NONE) ){
	  wlzErrno = WLZ_ERR_UNSPECIFIED;
	  break;
	}

	wlzErrno = WlzGreySetRange(tempobj, min, max, Min, Max);
	WlzFreeObj( tempobj );
	if( wlzErrno != WLZ_ERR_NONE ){
	  break;
	}
      }
      
      return wlzErrno;

    case WLZ_TRANS_OBJ:
      return WlzGreySetRange(obj->values.obj, min, max, Min, Max);

    case WLZ_EMPTY_OBJ:
      return wlzErrno;

    default:
      wlzErrno = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( wlzErrno == WLZ_ERR_NONE ){
    WlzValueConvertPixel(&min, min, WLZ_GREY_DOUBLE);
    WlzValueConvertPixel(&max, max, WLZ_GREY_DOUBLE);
    WlzValueConvertPixel(&Min, Min, WLZ_GREY_DOUBLE);
    WlzValueConvertPixel(&Max, Max, WLZ_GREY_DOUBLE);
    if( fabs(max.v.dbv - min.v.dbv) < DBL_EPSILON ){
      return WLZ_ERR_FLOAT_DATA;
    }
    factor = (Max.v.dbv - Min.v.dbv) / (max.v.dbv - min.v.dbv);

    WlzInitGreyScan(obj, &iwsp, &gwsp);
    while( WlzNextGreyInterval(&iwsp) == 0 ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++)
	  *gptr.inp = (int) (factor * (*gptr.inp - min.v.dbv) + Min.v.dbv);
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp.colrmn; i++, gptr.shp++)
	  *gptr.shp = (short) (factor * (*gptr.shp - min.v.dbv) + Min.v.dbv);
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp.colrmn; i++, gptr.ubp++)
	  *gptr.ubp = (UBYTE) (factor * (*gptr.ubp - min.v.dbv) + Min.v.dbv);
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp.colrmn; i++, gptr.flp++)
	  *gptr.flp = (float) (factor * (*gptr.flp - min.v.dbv) + Min.v.dbv);
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp.colrmn; i++, gptr.dbp++)
	  *gptr.dbp = (double) (factor * (*gptr.dbp - min.v.dbv) + Min.v.dbv);
	break;

      default:
	wlzErrno = WLZ_ERR_GREY_TYPE;
	break;
      }
    }
  }

  return wlzErrno;
}
