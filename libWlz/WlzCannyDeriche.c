#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzCannyDeriche.c
* Date:         May 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      A Canny/Deriche edge detection filter.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzCannyDeriche
* Returns:	WlzObject *:		Edge detected object or NULL
*					on error.
* Purpose:	A Canny/Deriche edge detection filter which performs
*		the following steps
*		  1. Recursive filter with the Deriche edge
*		     operator applied in each dimension.
*		  2. Modulus of gradient image generated.
*		  3. Non-maximal suppression of the gradient
*		     image.
*		  4. Hysteresis threshold of NM suppressed gradient
*		     image.
*
* Global refs:	-
* Parameters:	WlzObject **dstGObj:	Destination pointer for the
*					gradient image, may be NULL.
*		WlzObject *srcObj:	Given domain object with grey
*					values.
*		double alpha:		Deriche filter parameter.
*		double mult:		Filter multiplier.
*		WlzPixelV pMinGrdV:	Primary hysteresis threshold
*					value.
*		WlzPixelV sMinGrdV:	Secondary hysteresis threshold
*					value.
************************************************************************/
WlzObject 	*WlzCannyDeriche(WlzObject **dstGObj, WlzObject *srcObj,
				 double alpha, double mult,
				 WlzPixelV pMinGrdV, WlzPixelV sMinGrdV,
			   	 WlzErrorNum *dstErr)
{
  WlzObject     *gMObj = NULL,
		*gZObj = NULL,
		*gYObj = NULL,
		*gXObj = NULL,
		*nMObj = NULL,
		*nMGObj = NULL,
		*hTObj = NULL,
		*dstObj = NULL;
  WlzRsvFilter  *ftr = NULL;
  WlzPixelV     zeroBV;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  zeroBV.type = WLZ_GREY_INT;
  zeroBV.v.inv = 0;
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
  else if((ftr = WlzRsvFilterMakeFilter(WLZ_RSVFILTER_NAME_DERICHE_1,
					alpha, &errNum)) != NULL)
  {
    ftr->c *= mult;
    gMObj = WlzAssignObject(
	    WlzGreyGradient(&gZObj, &gYObj, &gXObj, srcObj, ftr,
			    &errNum), NULL);
    WlzRsvFilterFreeFilter(ftr);
    if(errNum == WLZ_ERR_NONE)
    {
      nMObj = WlzAssignObject(
	      WlzNMSuppress(gMObj, gZObj, gYObj, gXObj, sMinGrdV,
	       		     &errNum), NULL);
    }
    if(gZObj)
    {
      WlzFreeObj(gZObj);
    }
    if(gYObj)
    {
      WlzFreeObj(gYObj);
    }
    if(gXObj)
    {
      WlzFreeObj(gXObj);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nMGObj = WlzAssignObject(
      	       WlzMakeMain(nMObj->type, nMObj->domain, gMObj->values,
	     		   NULL, gMObj, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      hTObj = WlzAssignObject(
      	      WlzHyThreshold(nMGObj, pMinGrdV, sMinGrdV, WLZ_THRESH_HIGH,
	      		     WLZ_8_CONNECTED, &errNum), NULL);
    }
    if(nMGObj)
    {
      WlzFreeObj(nMGObj);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      dstObj = WlzAssignObject(
      	       WlzMakeMain(hTObj->type, hTObj->domain, nMObj->values,
	     		 NULL, nMObj, &errNum), NULL);
    }
    if(hTObj)
    {
      WlzFreeObj(hTObj);
    }
    if(dstGObj)
    {
      *dstGObj = gMObj;
      gMObj = NULL;
    }
  }
  if(gMObj)
  {
    WlzFreeObj(gMObj);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstErr)
    {
      *dstErr = errNum;
    }
    dstObj = NULL;
  }
  return(dstObj);
}

/* #define TEST_WLZNMSUPRESS */
#ifdef TEST_WLZNMSUPRESS
int             main(int argc, char *argv[])
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzPixelV	sMinGrdV,
  		pMinGrdV;
  double        alpha,
  		mult;
  WlzObject     *inObj = NULL,
		*outObj = NULL;

  alpha = 1.0;
  mult = 10.0;
  sMinGrdV.type = WLZ_GREY_INT;
  sMinGrdV.v.inv = 100;
  pMinGrdV.type = WLZ_GREY_INT;
  pMinGrdV.v.inv = 400;
  if(((inObj = WlzAssignObject(WlzReadObj(stdin, &errNum), NULL)) == NULL) ||
     (errNum != WLZ_ERR_NONE))
  {
    (void )fprintf(stderr,
		   "%s: failed to read object from stdin\n",
		   *argv);
  }
  else if(((outObj = WlzCannyDeriche(NULL, inObj, alpha, mult,
  				     pMinGrdV, sMinGrdV,
  			             &errNum)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
  {
    (void )fprintf(stderr,
		   "%s: failed to Canny filter object\n",
		   *argv);
  }
  else if(WlzWriteObj(stdout, outObj) != WLZ_ERR_NONE)
  {
    (void )fprintf(stderr,
		   "%s: failed to write object\n",
		   *argv);
  }
  return(errNum);
}
#endif /* TEST_WLZNMSUPRESS */
