#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzMakeAffineTransform.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Make and free Woolz affine transforms.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzMakeAffineTransform					*
*   Returns    : WlzAffineTransform *:	required structure initialized	*
*					to the identity transform.	*
*   Parameters : WlzTransformType type: 2D or 3D type			*
*   Date       : Wed Oct 16 23:49:50 1996				*
*   Synopsis   : Allocates and initialises space for a transform	*
*		 structure						*
************************************************************************/

WlzAffineTransform *WlzMakeAffineTransform(
  WlzTransformType type,
  WlzErrorNum	*dstErr)
{
  WlzAffineTransform	*trans=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if((trans = (WlzAffineTransform *)
      AlcCalloc(1, sizeof(WlzAffineTransform))) == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  if((errNum == WLZ_ERR_NONE) &&
     (AlcDouble2Calloc(&trans->mat, 4, 4) != ALC_ER_NONE) ){
    AlcFree( (void *) trans );
    trans = NULL;
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  if( errNum == WLZ_ERR_NONE ){
    trans->type = type;

    /* initialize to the identity transform */
    trans->scale = 1.0;
    trans->mat[0][0] = 1.0;
    trans->mat[1][1] = 1.0;
    trans->mat[2][2] = 1.0;
    trans->mat[3][3] = 1.0;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(trans);
}

/************************************************************************
*   Function   : WlzFreeAffineTransform					*
*   Returns    : int:	error return - WLZ_ERR_NONE - success		*
*   Parameters : WlzAffineTransform *trans: transform to be freed	*
*   Date       : Wed Oct 16 23:53:08 1996				*
*   Synopsis   : Free space allocated by WlzMakeAffineTransform		*
************************************************************************/
WlzErrorNum WlzFreeAffineTransform(WlzAffineTransform *trans)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check for NULL pointer */
  if( trans == NULL ){
    return( WLZ_ERR_NONE );
  }

  /* check the linkcount - decrement if > 1
     set to -1 otherwise */
  if( WlzUnlink(&(trans->linkcount), &errNum) ){
    /* free the matrix - assumes allocated by AlcDouble2Alloc
       then the structure */
    AlcDouble2Free(trans->mat);
    AlcFree((void *) trans);
  }

  return errNum;
}
