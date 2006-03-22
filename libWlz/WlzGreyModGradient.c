#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreyModGradient_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzGreyModGradient.c
* \author       Richard Baldock
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
* \brief	Functions to calculate the modulus of the grey-level
* 		gradient of Woolz objects.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/* function:     WlzGreyModGradient    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Calculate the modulus of the grey-level gradient at each
 point. The gradient images are calculated using WlzGauss2() with width
parameter set to <tt>width</tt>. Will now calculate the modulus for each
object of a compound object if appropriate. 
*
* \return       Object with values set to the gradient modulus at each pixel.
* \param    obj	Input object.
* \param    width	Width parameter for the gaussian gradient operator.
* \param    dstErr	Error return.
* \par      Source:
*                WlzGreyModGradient.c
*/
WlzObject *WlzGreyModGradient(
  WlzObject	*obj,
  double	width,
  WlzErrorNum	*dstErr)
{
  WlzObject		*xobj, *yobj, *returnobj=NULL;
  WlzIntervalWSpace	iwsp1, iwsp2, iwsp3;
  WlzGreyWSpace		gwsp1, gwsp2, gwsp3;
  WlzCompoundArray	*cobj1, *cobj2;
  int			i;
  double		g1, g2, g3;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core ){
	switch( obj->domain.core->type ){
	case WLZ_EMPTY_DOMAIN:
	  returnobj = WlzMakeEmpty(&errNum);
	  break;

	default:
	  if( obj->values.core ){
	    if( obj->values.core->type == WLZ_EMPTY_VALUES ){
	      errNum = WLZ_ERR_VALUES_TYPE;
	    }
	  }
	  else {
	    errNum = WLZ_ERR_VALUES_NULL;
	  }
	  break;
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      /* one last check for rgb values */
      if(WlzGreyTableTypeToGreyType(obj->values.core->type, NULL)
	 == WLZ_GREY_RGBA){
	return WlzRGBAModGradient(obj, width, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      returnobj = WlzMakeEmpty(&errNum);
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      cobj1 = (WlzCompoundArray *) obj;
      if( cobj2 = WlzMakeCompoundArray(cobj1->type, 1, cobj1->n, NULL,
				       cobj1->otype, &errNum) ){
	/* transform each object, ignore type errors */
	for(i=0; i < cobj1->n; i++){
	  cobj2->o[i] =
	    WlzAssignObject(WlzGreyModGradient(cobj1->o[i],
					       width, &errNum), NULL);
	}
	return (WlzObject *) cobj2;
      }
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* if UBYTE grey values then copy to int */
  if( (errNum == WLZ_ERR_NONE) && !returnobj ){
    if(WlzGreyTableTypeToGreyType(obj->values.core->type, NULL) == WLZ_GREY_UBYTE)
    {
      returnobj = WlzConvertPix(obj, WLZ_GREY_INT, NULL);
    }
    else
    {
      returnobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj->domain,
			      WlzCopyValues(obj->type, obj->values,
					    obj->domain, &errNum),
			      NULL, NULL, NULL);
    }

    /* calculate gradient images */
    xobj = WlzGauss2(returnobj, width, width, 1, 0, NULL);
    yobj = WlzGauss2(returnobj, width, width, 0, 1, NULL);

    /* calculate modulus  - lockstep raster scan assumes equal domains */
    errNum = WlzInitGreyScan(returnobj, &iwsp1, &gwsp1);
    errNum = WlzInitGreyScan(xobj, &iwsp2, &gwsp2);
    errNum = WlzInitGreyScan(yobj, &iwsp3, &gwsp3);
    while((errNum == WLZ_ERR_NONE) && 
	  (WlzNextGreyInterval(&iwsp1) == WLZ_ERR_NONE) )
    {
      (void) WlzNextGreyInterval(&iwsp2);
      (void) WlzNextGreyInterval(&iwsp3);

      switch( gwsp1.pixeltype )
      {
      default:
      case WLZ_GREY_INT:
	for(i=0; i < iwsp1.colrmn; i++)
	{
	  g2 = *gwsp2.u_grintptr.inp++;
	  g3 = *gwsp3.u_grintptr.inp++;
	  g1 = sqrt( g2*g2 + g3*g3 );
	  *gwsp1.u_grintptr.inp++ = (int) g1;
	}
	break;

      case WLZ_GREY_SHORT:
	for(i=0; i < iwsp1.colrmn; i++)
	{
	  g2 = *gwsp2.u_grintptr.shp++;
	  g3 = *gwsp3.u_grintptr.shp++;
	  g1 = sqrt( g2*g2 + g3*g3 );
	  *gwsp1.u_grintptr.shp++ = (short) g1;
	}
	break;

      case WLZ_GREY_UBYTE:
	for(i=0; i < iwsp1.colrmn; i++)
	{
	  g2 = *gwsp2.u_grintptr.ubp++;
	  g3 = *gwsp3.u_grintptr.ubp++;
	  g1 = sqrt( g2*g2 + g3*g3 );
	  *gwsp1.u_grintptr.ubp++ = (UBYTE) g1;
	}
	break;

      case WLZ_GREY_FLOAT:
	for(i=0; i < iwsp1.colrmn; i++)
	{
	  g2 = *gwsp2.u_grintptr.flp++;
	  g3 = *gwsp3.u_grintptr.flp++;
	  g1 = sqrt( g2*g2 + g3*g3 );
	  *gwsp1.u_grintptr.flp++ = (float) g1;
	}
	break;

      case WLZ_GREY_DOUBLE:
	for(i=0; i < iwsp1.colrmn; i++)
	{
	  g2 = *gwsp2.u_grintptr.dbp++;
	  g3 = *gwsp3.u_grintptr.dbp++;
	  g1 = sqrt( g2*g2 + g3*g3 );
	  *gwsp1.u_grintptr.dbp++ = (double) g1;
	}
	break;

      case WLZ_GREY_RGBA: /* RGBA to be done - should not get here */
	errNum = WLZ_ERR_GREY_TYPE;
	break;

      }
    }

    /* check for normal return */
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }

    /* clean up */
    WlzFreeObj( xobj );
    WlzFreeObj( yobj );
  }

  /* check error return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return( returnobj );
}

