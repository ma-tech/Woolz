#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyInvertMinMax.c
* Date:         March 1999
* Author:       Margaret Stark
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Grey value inversion.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

static WlzErrorNum WlzGreyInvertMinMax3d(WlzObject	*obj,
					 WlzPixelV	min,
					 WlzPixelV	max);

/************************************************************************
*   Function   : WlzGreyInvertMinMax					*
*   Date       : Tue Jan  7 15:37:35 1997				*
*************************************************************************
*   Synopsis   : Invert the Grey values of an object within a given	*
*		range. The supplied min and max values define the grey	*
*		transformation by: 		       	       		*
*			g' = gmax + gmin - g				*
*		which means that g=gmax gives g'=gmin and g=gmin gives	*
*		g'=gmax. This is a generalisation of the original invert*
*		and allows for the "normal" inversion of 255-g for byte	*
*		images. The user must ensure the grey-value type of the	*
*		given object can hold the result e.g. if the result is	*
*		negative. Note it is assumed that min and max have 	*
*		a grey type which can be converted to the grey-type	*
*		of the given image.					*
*   Returns    : WlzErrorNum:    WLZ_ERR_NONE on success                *
*                Possible Errors:                                       *
*              :                  WLZ_ERR_INT_DATA,                     *
*              :                  WLZ_ERR_OBJECT_NULL     		*
*              :                  WLZ_ERR_OBJECT_TYPE                   *
*              :                  WLZ_ERR_DOMAIN_NULL                   *
*              :                  WLZ_ERR_VALUES_NULL                   *
*              :                  WLZ_ERR_GREY_TYPE                     *
*   Parameters : Wlzobject	  *obj: Object Pointer	                *
*		 WlzPixelV	  min: minimum grey of transform	*
*				  max: maximum grey of transform	*
*   Global refs: None.							*
************************************************************************/

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
	*g.shp = irange - *g.shp;
      }
      break;

    case WLZ_GREY_UBYTE: 
      for (i=0; i<iwsp.colrmn; i++, g.ubp++){
	*g.ubp = irange - *g.ubp;
      }
      break;

    case WLZ_GREY_FLOAT: 
      for (i=0; i<iwsp.colrmn; i++, g.flp++){
	*g.flp = drange - *g.flp;
      }
      break;

    case WLZ_GREY_DOUBLE: 
      for (i=0; i<iwsp.colrmn; i++, g.dbp++){
	*g.dbp = drange - *g.dbp;
      }
      break;

    default:
      return( WLZ_ERR_GREY_TYPE );

    }
  }
  return WLZ_ERR_NONE;
}


/************************************************************************
*   Function   : WlzGreyInvertMinMax3d					*
*   Date       : Tue Jan 14 12:14:09 1997				*
*************************************************************************
*   Synopsis   :Private procedure for WlzGreyInvertMinMax to implement	*
*	        the function for 3D domain objects. This procedure omits*
*	        most of the object checks and therefore should never be	*
*	        called directly.					*
*   Returns    :WlzErrorNum:  see WlzGreyInvertMinMax		       	*
*   Parameters :see WlzGreyInvertMinMax	       				*
*   Global refs:None							*
************************************************************************/

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
