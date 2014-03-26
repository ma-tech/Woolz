#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreySetValue_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreySetValue.c
* \author       Richard Baldock, Bill Hill
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
* \brief	Sets the grey values of objects.
* \ingroup	WlzValuesUtils
*/

#include <stdlib.h>
#include <Wlz.h>


/*! 
* \return       Woolz error code.
* \ingroup      WlzValuesUtils
* \brief        Set the grey value of every pixel/voxel to the given value.
* \param    	obj				Input object.
* \param    	val				New grey value.
*/
WlzErrorNum WlzGreySetValue(
  WlzObject	*obj,
  WlzPixelV	val)
{
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
      else if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( WlzGreyTableIsTiled(obj->values.core->type) ){
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
      else if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.core->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else {
	/* set range of each plane if non-empty - indicated by NULL */
	domains = obj->domain.p->domains;
	values = obj->values.vox->values;
	nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i = 0; i < nplanes; ++i)
	{
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzDomain    *doms;
	    WlzValues    *vals;
	    WlzErrorNum  errNum2D = WLZ_ERR_NONE;

	    doms = domains + i;
	    vals = values + i;
	    if(((*doms).core != NULL) && ((*vals).core != NULL))
	    {
	      WlzObject *obj2D = NULL;

	      if((obj2D = WlzAssignObject(
		      WlzMakeMain(WLZ_2D_DOMAINOBJ, *doms, *vals,
			NULL, NULL,
			&errNum2D), NULL)) != NULL)
	      {
		errNum2D = WlzGreySetValue(obj2D, val);
	        (void )WlzFreeObj(obj2D);
	      }
#ifdef _OPENMP
	      {
		if(errNum2D != WLZ_ERR_NONE)
		{
#pragma omp critical
		  {
		    if(errNum == WLZ_ERR_NONE)
		    {
		      errNum = errNum2D;
		    }
		  }
		}
	      }
#endif
	    }
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
  }
  if( errNum == WLZ_ERR_NONE ){
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
    (void )WlzEndGreyScan(&gwsp);
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  return errNum;
}
