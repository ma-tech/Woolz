#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzMakeAffineTransform.c
* Date:         March 1999
* Author:       Richard Baldock, Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Make and free Woolz affine transforms.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 27-09-00 bill Make modifications for 3D affine transforms (more
* 		checking). Update comments.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzMakeAffineTransform
* Returns:	WlzAffineTransform *:	New affine transform initialized
*					to the identity transform.
* Purpose:	Allocates and initialises space for a 2D or 3D affine
*		transform. Sufficient space is always allocated for a
*		3D transform.
* Global refs:	-
* Parameters:	WlzTransformType type:	Transform type.
*		WlzErrorNum   *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
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
    trans->scale = 1.0;
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

/************************************************************************
* Function:	WlzFreeAffineTransform
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Free's an affine transform allocated by
*		WlzMakeAffineTransform.
* Global refs:	-
* Parameters:	WlzAffineTransform *trans: Affine transform to free.
************************************************************************/
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
