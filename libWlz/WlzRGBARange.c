#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzRGBARange.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri May 23 13:44:34 2003
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
* \ingroup      WlzValuesUtils
* \brief        Find range of values in a RGBA type image. Currently
 implemented is modulus range. Use WlzGreyRange to get the individual
 colour ranges.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>

WlzErrorNum WlzRGBAModulusRange(
  WlzObject	*obj,
  double	*min,
  double	*max)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  double	val, lmin, lmax;
  int		i, nplanes, initFlg;
  WlzGreyP		g;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzObject		tempobj;
  WlzPlaneDomain	*planedm;
  WlzValues 		*values;
  WlzDomain		*domains;

  /* object checks */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->values.core == NULL){
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if( (min == NULL) || (max == NULL) ){
    errNum = WLZ_ERR_PARAM_NULL;
  }

  /* object type checks */
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      WlzInitGreyScan(obj, &iwsp, &gwsp);
      while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	g = gwsp.u_grintptr;
	for(i=0; i<iwsp.colrmn; i++, g.rgbp++){
	  val = WLZ_RGBA_MODULUS(*g.rgbp);
	  if( !initFlg ){
	    lmin = val;
	    lmax = val;
	    initFlg = 1;
	  }
	  else {
	    if( val < lmin ){
	      lmin = val;
	    }
	    if( val > lmax ){
	      lmax = val;
	    }
	  }
	}
      }
      if( errNum = WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
      *min = lmin;
      *max = lmax;
      break;

    case WLZ_3D_DOMAINOBJ:
      planedm = obj->domain.p;
      switch( planedm->type ){

      case WLZ_PLANEDOMAIN_DOMAIN:
	nplanes = planedm->lastpl - planedm->plane1 + 1;
	domains = planedm->domains;
	if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	  return( WLZ_ERR_VOXELVALUES_TYPE );
	}
	values = obj->values.vox->values;
	tempobj.type = WLZ_2D_DOMAINOBJ;
	tempobj.plist = NULL;
	tempobj.assoc = NULL;
	for(i=0; i < nplanes; i++, domains++, values++){
	  if( (*domains).core == NULL || (*values).core == NULL ){
	    continue;
	  }

	  tempobj.domain = *domains;
	  tempobj.values = *values;
	  errNum = WlzRGBAModulusRange(&tempobj, min, max);
	  if( errNum != WLZ_ERR_NONE ){
	    return( errNum );
	  }
	  if( !initFlg ){
	    lmin = *min;
	    lmax = *max;
	    initFlg = 1;
	    continue;
	  }
	  else {
	    if( *min < lmin ){
	      lmin = *min;
	    }
	    if( *max > lmax ){
	      lmax = *max;
	    }
	  }
	}
	break;

      default:
	errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
	break;

      }
      *min = lmin;
      *max = lmax;
      break;

    case WLZ_TRANS_OBJ:
      return( WlzRGBAModulusRange(obj->values.obj, min, max) );

    case WLZ_EMPTY_OBJ:
      *min = 0.0;
      *max = 0.0;
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  return errNum;
}
