#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzLaplacian.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Implements a Laplacian edge enhancement filter.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzLaplacian						*
* Returns:	WlzObject *:		Laplacian filtered object.	*
* Purpose:	Applies a Laplacian edge enhancement filter to the	*
*		given woolz object.					*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		int kSize:		Kernel size, must be 3, 5 or 7.	*
*		int newObjFlag:		If zero the convolution is done	*
*					in place, else a new object is	*
*					created.			*
*		int modFlag:		Take the absolute value of the	*
*					convolution if non zero (this	*
*					is always the case for a UBYTE	*
*					object).			*
*		WlzErrorNum *wlzErr:	Destination error pointer, may	*
*					be NULL.			*
************************************************************************/
WlzObject	*WlzLaplacian(WlzObject *srcObj, int kSize,
			      int newObjFlag, int modFlag,
			      WlzErrorNum *wlzErr)
{
  WlzObject	*dstObj = NULL;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	lMask3I[] = 		     /* 3 by 3 with central weight 1 */
  {
    -1, -1, -1,
    -1,  9, -1,
    -1, -1, -1
  },
  		lMask5I[] =		     /* 5 by 5 with central weight 4 */
  {
    -1, -1, -1, -1, -1,
    -1,  0,  0,  0, -1,
    -1,  0, 20,  0, -1,
    -1,  0,  0,  0, -1,
    -1, -1, -1, -1, -1
  },
  		lMask7I[] = 		     /* 7 by 7 with central weight 4 */
  {
     0,  0,  0, -1,  0,  0,  0,
     0, -1,  0,  0,  0, -1,  0,
     0,  0,  0,  1,  0,  0,  0,
    -1,  0,  1,  8,  1,  0, -1,
     0,  0,  0,  1,  0,  0,  0,
     0, -1,  0,  0,  0, -1,  0,
     0,  0,  0, -1,  0,  0,  0
  };
  const double	lMask3D[] = 		     /* 3 by 3 with central weight 1 */
  {
    -1.0, -1.0, -1.0,
    -1.0,  9.0, -1.0,
    -1.0, -1.0, -1.0
  },
  		lMask5D[] =		     /* 5 by 5 with central weight 4 */
  {
    -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0,  0.0,  0.0,  0.0, -1.0,
    -1.0,  0.0, 20.0,  0.0, -1.0,
    -1.0,  0.0,  0.0,  0.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0
  },
  		lMask7D[] = 		     /* 7 by 7 with central weight 4 */
  {
     0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,
     0.0, -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
     0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,
    -1.0,  0.0,  1.0,  8.0,  1.0,  0.0, -1.0,
     0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,
     0.0, -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
     0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0
  };
  WlzConvolution lConv;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzLaplacian FE 0x%lx %d %d %d 0x%lx\n",
	   (unsigned long )srcObj, kSize, newObjFlag, modFlag,
	   (unsigned long)wlzErr));
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->type != WLZ_2D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    gType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    lConv.linkcount = 0;
    lConv.xsize = kSize;
    lConv.ysize = kSize;
    lConv.cv = NULL;
    lConv.divscale = 1;
    lConv.offset = 0;
    lConv.modflag = gType == WLZ_GREY_UBYTE? 1: modFlag;
    switch(gType)
    {
      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_UBYTE:
        lConv.type = WLZ_CONVOLVE_INT;
	switch(kSize)
	{
	  case 3:
    	    lConv.divscale = 1;
	    lConv.cv = (int *)(&(lMask3I[0]));
	    break;
	  case 5:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask5I[0]));
	    break;
	  case 7:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask7I[0]));
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
        break;
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
        lConv.type = WLZ_CONVOLVE_FLOAT;
	switch(kSize)
	{
	  case 3:
    	    lConv.divscale = 1;
	    lConv.cv = (int *)(&(lMask3D[0]));
	    break;
	  case 5:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask5D[0]));
	    break;
	  case 7:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask7D[0]));
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_DATA;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzConvolveObj(srcObj, &lConv, newObjFlag, &errNum);
  }
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzLaplacian FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}
