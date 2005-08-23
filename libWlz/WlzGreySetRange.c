#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzGreySetRange.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri May 23 08:09:14 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzValuesFilters
* \brief        Set new grey range for a Woolz object using
 simple linear interpolation
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>


/* function:     WlzGreySetRange    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        Set new grey-range by simple linear interpolation. It
 assumes that the min and max values enclose the true range of
 grey-values in the object. Failure to check this could result in a
 segmentation fault.

The transform function is:
\f[g' = \frac{(gMax - gMin)}{(gmax - gmin)} (g - gmin) + gMin
\f]
*
* \return       Woolz error number
* \param    obj	Input grey-level object whose values are to be reset.
* \param    min	Initial minimum value
* \param    max	Initial maximum value
* \param    Min	Final minimum value
* \param    Max	Final maximum value
* \par      Source:
*                WlzGreySetRange.c
*/
WlzErrorNum WlzGreySetRange(
  WlzObject	*obj,
  WlzPixelV	min,
  WlzPixelV	max,
  WlzPixelV	Min,
  WlzPixelV	Max)
{
  double		factor, val;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzObject		*tempobj;
  WlzValues 		*values;
  WlzDomain		*domains;
  int			i, j, nplanes;
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
    /* get conversion function - should use LUT
       use 4 LUTS for rgb type since bounded */
    if( WlzGreyTypeFromObj(obj, &wlzErrno) == WLZ_GREY_RGBA ){
      UBYTE	rgbaLut[4][256];
      UINT	rgbamin[4], rgbaMin[4];
      double	rgbaFactor[4], val;
      UINT	red, green, blue, alpha;

      rgbamin[0] = WLZ_RGBA_RED_GET(min.v.rgbv);
      rgbaMin[0] = WLZ_RGBA_RED_GET(Min.v.rgbv);
      if(WLZ_RGBA_RED_GET(max.v.rgbv) > WLZ_RGBA_RED_GET(min.v.rgbv)){
	rgbaFactor[0] = (((double) WLZ_RGBA_RED_GET(Max.v.rgbv) - 
			  WLZ_RGBA_RED_GET(Min.v.rgbv))/
			 (WLZ_RGBA_RED_GET(max.v.rgbv) - 
			  WLZ_RGBA_RED_GET(min.v.rgbv)));
      }
      else {
	rgbaFactor[0] = 0.0;
      }

      rgbamin[1] = WLZ_RGBA_GREEN_GET(min.v.rgbv);
      rgbaMin[1] = WLZ_RGBA_GREEN_GET(Min.v.rgbv);
      if(WLZ_RGBA_GREEN_GET(max.v.rgbv) > WLZ_RGBA_GREEN_GET(min.v.rgbv)){
	rgbaFactor[1] = (((double) WLZ_RGBA_GREEN_GET(Max.v.rgbv) - 
			  WLZ_RGBA_GREEN_GET(Min.v.rgbv))/
			 (WLZ_RGBA_GREEN_GET(max.v.rgbv) - 
			  WLZ_RGBA_GREEN_GET(min.v.rgbv)));
      }
      else {
	rgbaFactor[1] = 0.0;
      }

      rgbamin[2] = WLZ_RGBA_BLUE_GET(min.v.rgbv);
      rgbaMin[2] = WLZ_RGBA_BLUE_GET(Min.v.rgbv);
      if(WLZ_RGBA_BLUE_GET(max.v.rgbv) > WLZ_RGBA_BLUE_GET(min.v.rgbv)){
	rgbaFactor[2] = (((double) WLZ_RGBA_BLUE_GET(Max.v.rgbv) - 
			  WLZ_RGBA_BLUE_GET(Min.v.rgbv))/
			 (WLZ_RGBA_BLUE_GET(max.v.rgbv) - 
			  WLZ_RGBA_BLUE_GET(min.v.rgbv)));
      }
      else {
	rgbaFactor[2] = 0.0;
      }

      rgbamin[3] = WLZ_RGBA_ALPHA_GET(min.v.rgbv);
      rgbaMin[3] = WLZ_RGBA_ALPHA_GET(Min.v.rgbv);
      if(WLZ_RGBA_ALPHA_GET(max.v.rgbv) > WLZ_RGBA_ALPHA_GET(min.v.rgbv)){
	rgbaFactor[3] = (((double) WLZ_RGBA_ALPHA_GET(Max.v.rgbv) - 
			  WLZ_RGBA_ALPHA_GET(Min.v.rgbv))/
			 (WLZ_RGBA_ALPHA_GET(max.v.rgbv) - 
			  WLZ_RGBA_ALPHA_GET(min.v.rgbv)));
      }
      else {
	rgbaFactor[3] = 0.0;
      }

      /* now set up the LUTS */
      for(i=0; i < 4; i++){
	for(j=0; j < 256; j++){
	  val = rgbaFactor[i] * (j - rgbamin[i]) + rgbaMin[i];
	  rgbaLut[i][j] = WLZ_CLAMP(val, 0, 255);
	}
      }

      /* set values - can assume rgba grey-type */
      WlzInitGreyScan(obj, &iwsp, &gwsp);
      while( WlzNextGreyInterval(&iwsp) == WLZ_ERR_NONE ){

	gptr = gwsp.u_grintptr;
	for (i=0; i<iwsp.colrmn; i++, gptr.rgbp++){
	  red = rgbaLut[0][WLZ_RGBA_RED_GET(*gptr.rgbp)];
	  green = rgbaLut[0][WLZ_RGBA_GREEN_GET(*gptr.rgbp)];
	  blue = rgbaLut[0][WLZ_RGBA_BLUE_GET(*gptr.rgbp)];
	  alpha = rgbaLut[0][WLZ_RGBA_ALPHA_GET(*gptr.rgbp)];
	  WLZ_RGBA_RGBA_SET(*gptr.rgbp, red, green, blue, alpha);
	}
      }
    }
    else {
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
	  for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
            val = factor * (*gptr.inp - min.v.dbv) + Min.v.dbv;
	    *gptr.inp = WLZ_NINT(val);
          }
	  break;

	case WLZ_GREY_SHORT:
	  for (i=0; i<iwsp.colrmn; i++, gptr.shp++){
            val = factor * (*gptr.shp - min.v.dbv) + Min.v.dbv;
	    *gptr.shp = WLZ_NINT(val);
          }
	  break;

	case WLZ_GREY_UBYTE:
	  for (i=0; i<iwsp.colrmn; i++, gptr.ubp++){
            val = factor * (*gptr.ubp - min.v.dbv) + Min.v.dbv;
	    *gptr.ubp = WLZ_NINT(val);
          }
	  break;

	case WLZ_GREY_FLOAT:
	  for (i=0; i<iwsp.colrmn; i++, gptr.flp++){
	    *gptr.flp = (float) (factor * (*gptr.flp - min.v.dbv) + Min.v.dbv);
          }
	  break;

	case WLZ_GREY_DOUBLE:
	  for (i=0; i<iwsp.colrmn; i++, gptr.dbp++){
	    *gptr.dbp = (double) (factor * (*gptr.dbp - min.v.dbv) + Min.v.dbv);
          }
	  break;

	default:
	  wlzErrno = WLZ_ERR_GREY_TYPE;
	  break;
	}
      }
    }
  }

  return wlzErrno;
}
