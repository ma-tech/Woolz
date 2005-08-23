#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzGreyCrossing.c
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
* \brief	Computes a domain object with grey values in which the
*		values encode the direction at each pixel of a grey value
*		transition. If the given grey value is zero then this is the
*		zero crossing direction.
* \ingroup	WlzValuesFilters
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	Direction value.
* \ingroup	WlzValuesFilters
* \brief	Tests for a grey value crossing using a 3x3 kernel.
*               This function is called by WlzSeqPar().
*               Given an integer kernel about a given pixel (E):
*		\verbatim
                   A B C
                   D E F
                   G H I
		\endverbatim
*               The sums:
*		\verbatim
                  (A + B + C), (G + H + I); (B + C + F), (D + G + H)
                  (A + D + G), (C + F + I) and (A + B + D), (F + H + I)
		\endverbatim
*               are computed  about the pixel E.
*               The crossing value is then the maximum difference
*               between the sums which have opposite signs.
* \param	spWSpace		Work space data structure
*                                       from WlzSeqPar().
* \param	spData			Data passed on by WlzSeqPar(),
*                                       which is used to pass the
*                                       crossing value.
*/
int 		WlzGreyCrossingSeqParFn(WlzSeqParWSpace *spWSpace,
  			    	        void *spData)
{
  int 		cVal,
		dVal = 0,
  		sum0,
		sum1,
		sum2,
		sum3,
		sum4,
		sum5,
		sum6,
		sum7,
		zVal = 0;
  int		*krn0,
  		*krn1,
		*krn2;
  int		zThresh = 5 * 5;

  if(spWSpace && spData)
  {
    cVal = *(int *)spData * 3;
    krn0 = *(spWSpace->adrptr - 1) - 1;
    krn1 = *(spWSpace->adrptr) - 1;
    krn2 = *(spWSpace->adrptr + 1) - 1;
    /* Compute the sums. */
    sum0 = *krn0 + *(krn0 + 1) +  *(krn0 + 2);
    sum1 = *krn2 + *(krn2 + 1) +  *(krn2 + 2);
    sum2 = *(krn0 + 1) + *(krn0 + 2) + *(krn1 + 2);
    sum3 = *krn1 + *krn2 + *(krn2 + 1);
    sum4 = *krn0 + *krn1 + *krn2;
    sum5 = *(krn0 + 2) + *(krn1 + 2) + *(krn2 + 2);
    sum6 = *krn0 + *(krn0 + 1) + *krn1;
    sum7 = *(krn1 + 2) + *(krn2 + 1) + *(krn2 + 2);
    if(cVal != 0)
    {
      sum0 -= cVal; sum1 -= cVal; sum2 -= cVal; sum3 -= cVal;
      sum4 -= cVal; sum5 -= cVal; sum6 -= cVal; sum7 -= cVal;
    }
    /* Compute maximum crossing value. */
    if((((sum0 > 0) && (sum1 < 0)) || ((sum0 < 0) && (sum1 > 0))) &&
       (sum0 > sum1))
    {
      zVal = sum0 - sum1;
      dVal = 1;
    }
    if((((sum2 > 0) && (sum3 < 0)) || ((sum2 < 0) && (sum3 > 0))) &&
       (sum2 > sum3))
    {
      if(zVal < (sum2 - sum3))
      {
	zVal = sum2 - sum3;
	dVal = 2;
      }
    }
    if((((sum4 > 0) && (sum5 < 0)) || ((sum4 < 0) && (sum5 > 0))) &&
       (sum4 > sum5))
    {
      if(zVal < (sum4 - sum5))
      {
	zVal = sum4 - sum5;
	dVal = 3;
      }
    }
    if((((sum6 > 0) && (sum7 < 0)) || ((sum6 < 0) && (sum7 > 0))) &&
       (sum6 > sum7))
    {
      if(zVal < (sum6 - sum7))
      {
	zVal = sum6 - sum7;
	dVal = 4;
      }
    }
    if((zVal * zVal) < zThresh)
    {
      dVal = 0;
    }
  }
  return(dVal);
}

/*!
* \return	Grey-crossing object or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Uses a 3X3 kernel to examine the given object for
*               transitions (grey value crossings) about the given
*               grey value.
*               Only objects with WLZ_EMPTY_OBJ and WLZ_2D_DOMAINOBJ
*               types are valid. WLZ_2D_DOMAINOBJ objects must have
*               a valid domain and integral (ie WLZ_GREY_INT,
*               WLZ_GREY_SHORT or WLZ_GREY_UBYTE) values.
* \param	inObj			Given object.
* \param	newObjFlag		If zero the convolution is done
*                                       in place, else a new object is
*                                       created.
* \param	cVal			Grey value about which to test
*                                       for transitions.
* \param	dstErr
*/
WlzObject 	*WlzGreyCrossing(WlzObject *inObj, int newObjFlag,
				 int cVal, WlzErrorNum *dstErr)
{
  WlzObject 	*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzPixelV	bkgVal;
  WlzGreyType	gValType;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzGreyCrossing FE 0x%lx %d %d 0x%lx\n",
	   (unsigned long )inObj, newObjFlag, cVal,
	   (unsigned long )dstErr));
  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(inObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(inObj->type)
    {
      case WLZ_EMPTY_OBJ:
	outObj = (newObjFlag)? WlzMakeEmpty(&errNum): inObj;
	break;
      case WLZ_2D_DOMAINOBJ:
	bkgVal = WlzGetBackground(inObj, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzValueConvertPixel(&bkgVal, bkgVal, WLZ_GREY_INT);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  gValType = WlzGreyTableTypeToGreyType(inObj->values.core->type,
	  			                &errNum);
	  switch(gValType)	           /* Check for integral grey values */
	  {
	    case WLZ_GREY_INT:    /* FALLTHROUGH */
	    case WLZ_GREY_SHORT:  /* FALLTHROUGH */
	    case WLZ_GREY_UBYTE:
	      break;
	    case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
	    case WLZ_GREY_DOUBLE: /* FALLTHROUGH */
	    case WLZ_GREY_RGBA: /* FALLTHROUGH */
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzSeqPar(inObj, newObjFlag, 0, WLZ_RASTERDIR_ILIC,
			     3, bkgVal.v.inv, (void *)&cVal,
			     WlzGreyCrossingSeqParFn, &errNum);
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if( dstErr )
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzGreyCrossing FX 0x%lx\n",
	   (unsigned long )outObj));
  return(outObj);
}
