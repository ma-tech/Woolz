#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRGBAImageArithmetic_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzRGBAImageArithmetic.c
* \author       Richard Baldock
* \date         June 2004
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
* \brief	Performs image arithmetic on RGBA data.
* 		Both input files must be RGBA value type.
* \ingroup	WlzArithmetic
*/

#include <Wlz.h>

/*!
* \return	New object.
* \ingroup	WlzArithmetic
* \brief	Performs image arithmetic on objects with RGBA values.
* 		See WlzImageArithmetic().
* \param	obj0			First object.
* \param	obj1			Second object.
* \param	op			Binary operator.
* \param	overwrite		Allow the destination object
* 					to share values with one of
* 					the given objects if non zero.
* \param	dstErr			Destination error pointer, may
* 					be NULL.
*/
WlzObject *WlzRGBAImageArithmetic(
  WlzObject 		*obj0,
  WlzObject 		*obj1,
  WlzBinaryOperatorType	op,
  int 			overwrite,
  WlzErrorNum 		*dstErr)
{
  WlzObject	*rtnObj=NULL, *objs[4];
  WlzCompoundArray	*cmpnd0, *cmpnd1, *rtnCmpnd;
  int		i;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check objects */
  if( (obj0 == NULL) || (obj1 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;;
  }
  else if(obj0->type != obj1->type){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if( (overwrite < 0) || (overwrite > 2) ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((WlzGreyTypeFromObj(obj0, NULL) != WLZ_GREY_RGBA) ||
	  (WlzGreyTypeFromObj(obj1, NULL) != WLZ_GREY_RGBA) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else {
    /* convert to compound */
    cmpnd0 = WlzRGBAToCompound(obj0, WLZ_RGBA_SPACE_RGB, &errNum);
    cmpnd1 = WlzRGBAToCompound(obj1, WLZ_RGBA_SPACE_RGB, &errNum);

    /* operate on each channel */
    for(i=0; (i < 4) && (errNum == WLZ_ERR_NONE); i++){
      objs[i] = WlzImageArithmetic(cmpnd0->o[i], cmpnd1->o[i],
				   op, 0, &errNum);
    }

    /* create return compound object */
    if( errNum == WLZ_ERR_NONE ){
     rtnCmpnd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 4, &(objs[0]),
				     obj0->type, &errNum);
    }
    else {
      rtnCmpnd = NULL;
    }
    
    /* free temporary objects */
    WlzFreeObj((WlzObject *) cmpnd0);
    WlzFreeObj((WlzObject *) cmpnd1);
  }

  /* convert back to RGBA  */
  if( errNum == WLZ_ERR_NONE ){
    rtnObj = WlzCompoundToRGBA(rtnCmpnd, WLZ_RGBA_SPACE_RGB, &errNum);
  }
  if( rtnCmpnd ){
    WlzFreeObj((WlzObject *) rtnCmpnd);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
