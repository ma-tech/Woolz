#if defined(__GNUC__)
#ident "MRC HGU $Id:"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id:"
#else static char _WlzRGBA_calarBinaryOp_c[] = "MRC HGU $Id:";
#endif
#endif
/*!
* \file         WlzRGBAScalarBinaryOp.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Jan 25 10:55:52 2008
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par Copyright:
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
* \ingroup      WlzArithmetic
* \brief        Apply a scalar binary operation to an RGBA image.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>

/*!
* \return       Woolz error.
* \param    o1	Input object.
* \param    pval	Operand value.
* \param    o3	Oject for the return values. Setting equal to
 <tt>o1</tt> means values will be overwritten.
 * \param    op	Opertor to be applied.
* \par      Source:
*                WlzScalarBinaryOp.c
* \ingroup	WlzAritmetic
* \brief	Performsscalar operation on objects with RGBA values.
* 		See WlzScalarBinaryOp().
* 					be NULL.
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
    rtnObj = WlzCompoundToRGBA(cmpnd0, WLZ_RGBA_SPACE_RGB,
			       0, &errNum);
  }
  if( cmpnd0 ){
    WlzFreeObj((WlzObject *) cmpnd0);
  }


  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
