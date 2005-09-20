#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzMakeAffineTransform.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Allocation and freeing of affine transforms.
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	New affine transform initialized to the identity transform.
* \ingroup	WlzTransform
* \brief	Allocates and initialises space for a 2D or 3D affine
*               transform. Sufficient space is always allocated for a
*               3D transform.
*		The transform should be freed using WlzFreeAffineTransform().
* \param	type			Transform type.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzAffineTransform *WlzMakeAffineTransform(WlzTransformType type,
  				           WlzErrorNum *dstErr)
{
  WlzAffineTransform *trans=NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(type)
  {
    case WLZ_TRANSFORM_2D_AFFINE:
    case WLZ_TRANSFORM_2D_REG:
    case WLZ_TRANSFORM_2D_TRANS:
    case WLZ_TRANSFORM_2D_NOSHEAR:
    case WLZ_TRANSFORM_3D_AFFINE:
    case WLZ_TRANSFORM_3D_REG:
    case WLZ_TRANSFORM_3D_TRANS:
    case WLZ_TRANSFORM_3D_NOSHEAR:
      break;
    default:
      errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((trans = (WlzAffineTransform *)
                 AlcCalloc(1, sizeof(WlzAffineTransform))) == NULL) ||
       (AlcDouble2Calloc(&trans->mat, 4, 4) != ALC_ER_NONE))
    {
      if(trans)
      {
        AlcFree(trans);
	trans = NULL;
      }
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    trans->type = type;
    /* Initialize to the identity transform */
    trans->mat[0][0] = 1.0;
    trans->mat[1][1] = 1.0;
    trans->mat[2][2] = 1.0;
    trans->mat[3][3] = 1.0;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trans);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzTransform
* \brief	Frees an affine transform allocated by
*		WlzMakeAffineTransform().
* \param	trans			Affine transform to free.
*/
WlzErrorNum 	WlzFreeAffineTransform(WlzAffineTransform *trans)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(trans && WlzUnlink(&(trans->linkcount), &errNum))
  {
    /* Free the matrix - assumes allocated by AlcDouble2Alloc
       then the structure */
    if(trans->mat)
    {
      AlcDouble2Free(trans->mat);
    }
    AlcFree(trans);
  }
  return(errNum);
}
