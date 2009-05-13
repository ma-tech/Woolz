#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreyInvertMinMax_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzGreyInvertMinMax.c
* \author       Margaret Stark
* \date         September 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Grey value inversion.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

static WlzErrorNum WlzGreyInvertMinMax3d(WlzObject	*obj,
					 WlzPixelV	min,
					 WlzPixelV	max);

/* function:     WlzGreyInvertMinMax    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Invert the Grey values of an object within a given	
*		range. The supplied min and max values define the grey	
*		transformation by: 		       	       		
*		\f[
		g' = g_{max} + g_{min} - g				
		\f]
*		which means that g=gmax gives g'=gmin and g=gmin gives	
*		g'=gmax. This is a generalisation of the original invert
*		and allows for the "normal" inversion of 255-g for byte	
*		images. The user must ensure the grey-value type of the	
*		given object can hold the result e.g. if the result is	
*		negative. Note it is assumed that min and max have 	
*		a grey type which can be converted to the grey-type	
*		of the given image.
*
* \return       Woolz error.
* \param    obj	Input object.
* \param    min	Minimu grey value for the inversion function.
* \param    max	Maximun value for the inversion function.
* \par      Source:
*                WlzGreyInvertMinMax.c
*/
WlzErrorNum WlzGreyInvertMinMax(
  WlzObject	*obj,
  WlzPixelV	min,
  WlzPixelV	max)
{
 
  WlzGreyP		g;
  WlzIntervalWSpace 	iwsp;
  WlzGreyWSpace 	gwsp;
  int			irange;
  double		drange;
  WlzUInt		redrange, greenrange, bluerange;
  int			i;
 
  /* check for NULL object */
  if( obj == NULL ){
    return WLZ_ERR_OBJECT_NULL;
  }

  /* check the object types */
  switch( obj->type ){
 
  case WLZ_2D_DOMAINOBJ:
    if(obj->domain.core == NULL){ 
      return WLZ_ERR_DOMAIN_NULL;
    }
    if(obj->values.core == NULL){ 
      return WLZ_ERR_VALUES_NULL;
    }
    break;

  case WLZ_3D_DOMAINOBJ:
    if(obj->domain.core == NULL){ 
      return WLZ_ERR_DOMAIN_NULL;
    }
    if(obj->values.core == NULL){ 
      return WLZ_ERR_VALUES_NULL;
    }
    return WlzGreyInvertMinMax3d(obj, min, max);

  case WLZ_TRANS_OBJ:
    return WlzGreyInvertMinMax(obj->values.obj, min, max);

  case WLZ_EMPTY_OBJ:
    return WLZ_ERR_NONE;

  default:
    return WLZ_ERR_OBJECT_TYPE;
  }

  /* calculate the grey invert parameter */
  WlzInitGreyScan(obj,&iwsp,&gwsp);
  WlzValueConvertPixel(&min, min, gwsp.pixeltype);
  WlzValueConvertPixel(&max, max, gwsp.pixeltype);
  switch( gwsp.pixeltype ){
  case WLZ_GREY_INT:
    irange = min.v.inv + max.v.inv;
    break;
  case WLZ_GREY_SHORT:
    irange = min.v.shv + max.v.shv;
    break;
  case WLZ_GREY_UBYTE:
    irange = min.v.ubv + max.v.ubv;
    break;
  case WLZ_GREY_FLOAT:
    drange = min.v.flv + max.v.flv;
    break;
  case WLZ_GREY_DOUBLE:
    drange = min.v.dbv + max.v.dbv;
    break;
  case WLZ_GREY_RGBA:
    redrange = WLZ_RGBA_RED_GET(min.v.rgbv) + WLZ_RGBA_RED_GET(max.v.rgbv);
    greenrange = WLZ_RGBA_GREEN_GET(min.v.rgbv) +
                 WLZ_RGBA_GREEN_GET(max.v.rgbv);
    bluerange = WLZ_RGBA_BLUE_GET(min.v.rgbv) + WLZ_RGBA_BLUE_GET(max.v.rgbv);
    break;
  default:
    return WLZ_ERR_GREY_TYPE;
    break;
  }

  while( WlzNextGreyInterval(&iwsp) == 0 ){
    g = gwsp.u_grintptr;
    switch( gwsp.pixeltype ){
    case WLZ_GREY_INT:
      for (i=0; i<iwsp.colrmn; i++, g.inp++){
	*g.inp = irange - *g.inp;
      }
      break;

    case WLZ_GREY_SHORT:
      for (i=0; i<iwsp.colrmn; i++, g.shp++){
	*g.shp = (short )(irange - *g.shp);
      }
      break;

    case WLZ_GREY_UBYTE: 
      for (i=0; i<iwsp.colrmn; i++, g.ubp++){
	*g.ubp = (WlzUByte )(irange - *g.ubp);
      }
      break;

    case WLZ_GREY_FLOAT: 
      for (i=0; i<iwsp.colrmn; i++, g.flp++){
	*g.flp = (float )(drange - *g.flp);
      }
      break;

    case WLZ_GREY_DOUBLE: 
      for (i=0; i<iwsp.colrmn; i++, g.dbp++){
	*g.dbp = drange - *g.dbp;
      }
      break;

    case WLZ_GREY_RGBA: 
      for (i=0; i<iwsp.colrmn; i++, g.rgbp++){
	WlzUInt red, green, blue;
	red = redrange - WLZ_RGBA_RED_GET(*g.rgbp);
	green = greenrange - WLZ_RGBA_GREEN_GET(*g.rgbp);
	blue = bluerange - WLZ_RGBA_BLUE_GET(*g.rgbp);
	WLZ_RGBA_RED_SET(*g.rgbp, red);
	WLZ_RGBA_GREEN_SET(*g.rgbp, green);
	WLZ_RGBA_BLUE_SET(*g.rgbp, blue);
      }
      break;

    default:
      return( WLZ_ERR_GREY_TYPE );

    }
  }
  return WLZ_ERR_NONE;
}


/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	Private version of WlzGreyInvertMinMax() which omits
* 		the parameter checks.
* \param	obj			Given object.
* \param	min			Minimum grey value.
* \param	max			Maximum grey value.
*/
static WlzErrorNum WlzGreyInvertMinMax3d(WlzObject	*obj,
					 WlzPixelV	min,
					 WlzPixelV	max)
{
					 
		 
  WlzObject		o;
  WlzPlaneDomain	*pdom;
  WlzDomain		*domains;
  WlzValues		*values;
  int			nplanes;
 
  /* no need to check the object pointer or type because this procedure
     can only be accessed via WlzGreyInvertMinMax. The domain and valuetable
     types must be checked however */

  switch( obj->domain.p->type ){

  case WLZ_PLANEDOMAIN_DOMAIN:
    break;

  default:
    return WLZ_ERR_PLANEDOMAIN_TYPE;

  }
  switch( obj->values.vox->type ){

  case WLZ_VOXELVALUETABLE_GREY:
    break;

  default:
    return WLZ_ERR_VOXELVALUES_TYPE;

  }
  /* initialise variables */
  pdom =  obj->domain.p;
  nplanes = pdom->lastpl - pdom->plane1 + 1;
  domains = pdom->domains;
  values =  obj->values.vox->values;

  
  /* set up the temporary object */
  o.type = WLZ_2D_DOMAINOBJ;
  o.linkcount = 0;
  o.plist = NULL;
  o.assoc = NULL;
    
  /* loop through each plane setting values */
  while( nplanes > 0 ){

    o.domain = *domains;
    o.values = *values;
    if( (*domains).core != NULL && (*values).core != NULL ){
      WlzGreyInvertMinMax( &o, min, max );
    }
    nplanes--;
    domains++;
    values++;
  }
    
  return WLZ_ERR_NONE;
}
