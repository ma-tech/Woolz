#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzGreySetValue.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Sep 26 11:08:55 2003
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
* \brief        Functions to set the grey values of Woolz objects.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdlib.h>
#include <Wlz.h>


/* function:     WlzGreySetValue    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Set the grey value of every pixel/voxel tp
 <tt>val</tt>.
*
* \return       Woolz error.
* \param    obj	Input object.
* \param    val	New grey value.
* \par      Source:
*                WlzGreySetValue.c
*/
WlzErrorNum WlzGreySetValue(
  WlzObject	*obj,
  WlzPixelV	val)
{
  WlzObject		*tmpObj;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;
  WlzPixelV		tmpVal;
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
	    errNum = WlzGreySetValue(tmpObj, val);
	    WlzFreeObj( tmpObj );
	  }
	}
      }
      return errNum;

    case WLZ_TRANS_OBJ:
      return WlzGreySetValue(obj->values.obj, val);

    case WLZ_EMPTY_OBJ:
      return errNum;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
    WlzValueConvertPixel(&tmpVal, val, gwsp.pixeltype);
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++)
	  *gptr.inp = tmpVal.v.inv;
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp.colrmn; i++, gptr.shp++)
	  *gptr.shp = tmpVal.v.shv;
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp.colrmn; i++, gptr.ubp++)
	  *gptr.ubp = tmpVal.v.ubv;
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp.colrmn; i++, gptr.flp++)
	  *gptr.flp = tmpVal.v.flv;
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp.colrmn; i++, gptr.dbp++)
	  *gptr.dbp = tmpVal.v.dbv;
	break;

      case WLZ_GREY_RGBA:
	for (i=0; i<iwsp.colrmn; i++, gptr.rgbp++)
	  *gptr.rgbp = tmpVal.v.rgbv;
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
