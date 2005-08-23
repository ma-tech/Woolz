#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzRGBAGreyStats.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Jun 18 17:12:20 2004
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzFeatures
* \brief        Calculates simple quick statistics for a
 domain object with RGBA values.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdlib.h>
#include <Wlz.h>

int WlzRGBAGreyStats(
  WlzObject	*srcObj,
  WlzRGBAColorSpace	colSpc,
  WlzGreyType	*dstGType,
  double 	*dstMin,
  double 	*dstMax,
  double 	*dstSum,
  double 	*dstSumSq,
  double 	*dstMean,
  double 	*dstStdDev,
  WlzErrorNum 	*dstErr)
{
  int		area;
  WlzCompoundArray	*cmpnd;
  int		i;
  WlzGreyType	gType;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* all object checks done bt RGBAToCOmpound including
     grey-level type */
  if( cmpnd = WlzRGBAToCompound(srcObj, colSpc, &errNum) ){
    *dstGType = WLZ_GREY_RGBA;
    for(i=0; (i < 4) && (errNum == WLZ_ERR_NONE) ; i++){
      area = WlzGreyStats(cmpnd->o[i], &gType,
			  &(dstMin[i]), &(dstMax[i]),
			  &(dstSum[i]), &(dstSumSq[i]),
			  &(dstMean[i]), &(dstStdDev[i]),
			  &errNum);
    }
    WlzFreeObj((WlzObject *) cmpnd);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return area;
}
