#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRGBAScalarBinaryOp_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzRGBAScalarBinaryOp.c
* \author       Richard Baldock
* \date         January 2008
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
* \brief        Apply a scalar binary operation to an RGBA image.
* \ingroup      WlzArithmetic
*/

#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzAritmetic
* \brief	Performs scalar operation on objects with RGBA values.
* \param	o1			Input object.
* \param	pval			Operand value.
* \param	op			Opertor to be applied.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzRGBAScalarBinaryOp(
  WlzObject		*o1,
  WlzPixelV		pval,
  WlzBinaryOperatorType op,
  WlzErrorNum		*dstErr)
{
  WlzObject	*rtnObj;
  WlzCompoundArray	*cmpnd0;
  int		i;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check objects */
  if( (o1 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;;
  }
  else if((WlzGreyTypeFromObj(o1, NULL) != WLZ_GREY_RGBA) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else {
    /* convert to compound */
    cmpnd0 = WlzRGBAToCompound(o1, WLZ_RGBA_SPACE_RGB, &errNum);

    /* operate on each channel */
    for(i=0; (i < 4) && (errNum == WLZ_ERR_NONE); i++){
      WlzScalarBinaryOp(cmpnd0->o[i], pval, cmpnd0->o[i], op);
    }

  }


  /* convert back to RGBA  */
  if( errNum == WLZ_ERR_NONE ){
    rtnObj = WlzCompoundToRGBA(cmpnd0, WLZ_RGBA_SPACE_RGB, &errNum);
  }
  if( cmpnd0 ){
    WlzFreeObj((WlzObject *) cmpnd0);
  }


  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
