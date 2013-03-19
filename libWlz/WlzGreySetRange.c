#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreySetRange_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreySetRange.c
* \author       Richard Baldock
* \date         September 2003
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
* \brief	Sets the new grey range for an object using
*  		simple linear interpolation.
* \ingroup	WlzValuesFilters
*/

#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <Wlz.h>


/* function:     WlzGreySetRange    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        Set new grey-range by simple linear interpolation. It
 assumes that the min and max values enclose the true range of
 grey-values in the object. Failure to check this could result in a
 segmentation fault.

The transform function is:
\f[g' = \frac{(gMax - gMin)}{(gmax - gmin)} (g - gmin) + gMin + \delta
\f]
Here \f$\delta\f$ is the dither value.
*
* \return       Woolz error number
* \param    obj	Input grey-level object whose values are to be reset.
* \param    min	Initial minimum value
* \param    max	Initial maximum value
* \param    Min	Final minimum value
* \param    Max	Final maximum value
* \param    Dither values if destination range is greater than source range
*           and this flag is non-zero.
* \par      Source:
*                WlzGreySetRange.c
*/
WlzErrorNum WlzGreySetRange(
  WlzObject	*obj,
  WlzPixelV	min,
  WlzPixelV	max,
  WlzPixelV	Min,
  WlzPixelV	Max,
  int		dither)
{
  double		gMin = 0.0,
  			gMax = 0.0,
			sigma = 0.0,
  			factor,
			val;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzObject		*tempobj;
  WlzValues 		*values;
  WlzDomain		*domains;
  int			i, j, nplanes;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
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
			      &errNum);
	if((tempobj == NULL) && (errNum == WLZ_ERR_NONE) ){
	  errNum = WLZ_ERR_UNSPECIFIED;
	  break;
	}

	errNum = WlzGreySetRange(tempobj, min, max, Min, Max, dither);
	WlzFreeObj( tempobj );
	if( errNum != WLZ_ERR_NONE ){
	  break;
	}
      }
      
      return errNum;

    case WLZ_TRANS_OBJ:
      return WlzGreySetRange(obj->values.obj, min, max, Min, Max, dither);

    case WLZ_EMPTY_OBJ:
      return errNum;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    /* get conversion function - should use LUT
       use 4 LUTS for rgb type since bounded */
    if( WlzGreyTypeFromObj(obj, &errNum) == WLZ_GREY_RGBA ){
      WlzUByte	rgbaLut[4][256];
      WlzUInt	rgbamin[4], rgbaMin[4];
      double	rgbaFactor[4], val1;
      WlzUInt	red, green, blue, alpha;

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
	  val1 = rgbaFactor[i] * (j - rgbamin[i]) + rgbaMin[i];
	  rgbaLut[i][j] = (WlzUByte )WLZ_CLAMP(val1, 0, 255);
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
      WlzInitGreyScan(obj, &iwsp, &gwsp);
      factor = (Max.v.dbv - Min.v.dbv) / (max.v.dbv - min.v.dbv);
      if(dither)
      {
	if(fabs(factor) < 1.0 + DBL_EPSILON)
	{
	  dither = 0;
	}
	else
	{
	  sigma = fabs(2.0 / factor);
	  AlgRandSeed(101);
	  switch(gwsp.pixeltype) {
	    case WLZ_GREY_INT:
	      gMin = ALG_MAX(Min.v.dbv, INT_MIN);
	      gMax = ALG_MIN(Max.v.dbv, INT_MAX);
	      break;
	    case WLZ_GREY_SHORT:
	      gMin = ALG_MAX(Min.v.dbv, SHRT_MIN);
	      gMax = ALG_MIN(Max.v.dbv, SHRT_MAX);
	      break;
	    case WLZ_GREY_UBYTE:
	      gMin = ALG_MAX(Min.v.dbv, 0);
	      gMax = ALG_MIN(Max.v.dbv, 255);
	      break;
	    case WLZ_GREY_FLOAT:
	      gMin = ALG_MAX(Min.v.dbv, FLT_MIN);
	      gMax = ALG_MIN(Max.v.dbv, FLT_MAX);
	      break;
	    case WLZ_GREY_DOUBLE:
	      gMin = Min.v.dbv;
	      gMax = Max.v.dbv;
	      break;
	    default:
	      break;
	  }
	}
      }
      while((errNum == WLZ_ERR_NONE) &&
            ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE)){
	gptr = gwsp.u_grintptr;
	switch (gwsp.pixeltype) {
	  case WLZ_GREY_INT:
	    for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	      if(dither){
		val = factor * (*gptr.inp +
				AlgRandZigNormal(0.0, sigma) -
				min.v.dbv) + Min.v.dbv;
		val = WLZ_CLAMP(val, gMin, gMax);
	      } else {
		val = factor * (*gptr.inp - min.v.dbv) + Min.v.dbv;
	      }
	      *gptr.inp = WLZ_NINT(val);
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    for (i=0; i<iwsp.colrmn; i++, gptr.shp++){
	      if(dither){
		val = factor * (*gptr.shp +
				AlgRandZigNormal(0.0, sigma) -
				min.v.dbv) + Min.v.dbv;
		val = WLZ_CLAMP(val, gMin, gMax);
	      } else {
		val = factor * (*gptr.shp - min.v.dbv) + Min.v.dbv;
	      }
	      *gptr.shp = (short )WLZ_NINT(val);
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    for (i=0; i<iwsp.colrmn; i++, gptr.ubp++){
	      if(dither){
		val = factor * (*gptr.ubp +
				AlgRandZigNormal(0.0, sigma) -
				min.v.dbv) + Min.v.dbv;
		val = WLZ_CLAMP(val, gMin, gMax);
	      } else {
		val = factor * (*gptr.ubp - min.v.dbv) + Min.v.dbv;
	      }
	      *gptr.ubp = (WlzUByte )WLZ_NINT(val);
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    for (i=0; i<iwsp.colrmn; i++, gptr.flp++){
	      if(dither){
		val = factor * (*gptr.flp +
				AlgRandZigNormal(0.0, sigma) -
				min.v.dbv) + Min.v.dbv;
		val = WLZ_CLAMP(val, gMin, gMax);
	      } else {
		val = factor * (*gptr.flp - min.v.dbv) + Min.v.dbv;
	      }
	      *gptr.flp = (float )val;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    for (i=0; i<iwsp.colrmn; i++, gptr.dbp++){
	      if(dither){
		val = factor * (*gptr.dbp +
				AlgRandZigNormal(0.0, sigma) -
				min.v.dbv) + Min.v.dbv;
		val = WLZ_CLAMP(val, gMin, gMax);
	      } else {
		val = factor * (*gptr.dbp - min.v.dbv) + Min.v.dbv;
	      }
	      *gptr.dbp = val;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
      if(errNum == WLZ_ERR_EOO){
        errNum = WLZ_ERR_NONE;
      }
    }
  }
  return errNum;
}
