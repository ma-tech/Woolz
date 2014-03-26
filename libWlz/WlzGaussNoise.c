#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGaussNoise_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGaussNoise.c
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
* \brief	Functions for making Gaussian noise filled objects.
* \ingroup	WlzValuesFilters
*/

#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>


/* function:     WlzGaussNoise    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        Add Gaussian noise to each pixel with mean xero and sigma given
 by parameter <tt>val</tt>.
*
* \return       Woolz error.
* \param    obj	Input object
* \param    val	Sigma value for the aditive Guassian noise.
* \par      Source:
*                WlzGaussNoise.c
*/
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
      else if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if(WlzGreyTableIsTiled(obj->values.core->type)) {
        errNum = WLZ_ERR_VALUES_TYPE;
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

	  if((tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				    NULL, NULL,	NULL)) != NULL){
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
    if(errNum == WLZ_ERR_NONE) {
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
	    *gptr.inp = (int )WLZ_CLAMP(tmpVal.v.dbv, INT_MIN, INT_MAX);
	  }
	  break;

	case WLZ_GREY_SHORT:
	  for (i=0; i<iwsp.colrmn; i++, gptr.shp++){
	    mu = (double) *gptr.shp;
	    tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	    *gptr.shp = (short )WLZ_CLAMP(tmpVal.v.dbv, SHRT_MIN, SHRT_MAX);
	  }
	  break;

	case WLZ_GREY_UBYTE:
	  for (i=0; i<iwsp.colrmn; i++, gptr.ubp++){
	    mu = (double) *gptr.ubp;
	    tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	    *gptr.ubp = (WlzUByte )WLZ_CLAMP(tmpVal.v.dbv, 0, 255);
	  }
	  break;

	case WLZ_GREY_FLOAT:
	  for (i=0; i<iwsp.colrmn; i++, gptr.flp++){
	    mu = (double) *gptr.flp;
	    tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	    *gptr.flp = (float )WLZ_CLAMP(tmpVal.v.dbv, FLT_MIN, FLT_MAX);
	  }
	  break;

	case WLZ_GREY_DOUBLE:
	  for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	    mu = (double) *gptr.inp;
	    tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	    *gptr.inp = (int )(tmpVal.v.dbv);
	  }
	  break;

	case WLZ_GREY_RGBA:
	  for (i=0; i<iwsp.colrmn; i++, gptr.rgbp++){
	    mu = (double) *gptr.inp;
	    tmpVal.v.dbv = AlgRandNormal(mu, sigma);
	    *gptr.rgbp = (WlzUInt )WLZ_CLAMP(tmpVal.v.dbv, 0, 0xffffff);
	    *gptr.rgbp |= 0xff000000;
	  }
	  break;

	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
	}
      }
      (void )WlzEndGreyScan(&gwsp);
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  return errNum;
}
