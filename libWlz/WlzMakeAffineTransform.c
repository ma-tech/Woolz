#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMakeAffineTransform.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Allocation and freeing of Woolz affine transforms.
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
