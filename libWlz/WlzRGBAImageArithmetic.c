#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzRGBAImageArithmetic.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Thu Jun 10 16:39:26 2004
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
* \ingroup      WlzValuesFilters
* \brief        Perform image arithmetic on RGBA data. Both input files must be RGBA value type.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>

WlzObject *WlzRGBAImageArithmetic(
  WlzObject 		*obj0,
  WlzObject 		*obj1,
  WlzBinaryOperatorType	op,
  int 			overwrite,
  WlzErrorNum 		*dstErr)
{
  WlzObject	*rtnObj=NULL, *objs[4];
  WlzCompoundArray	*cmpnd0, *cmpnd1, *rtnCmpnd;
  int		i;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check objects */
  if( (obj0 == NULL) || (obj1 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;;
  }
  else if(obj0->type != obj1->type){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if( (overwrite < 0) || (overwrite > 2) ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((WlzGreyTypeFromObj(obj0, NULL) != WLZ_GREY_RGBA) ||
	  (WlzGreyTypeFromObj(obj1, NULL) != WLZ_GREY_RGBA) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else {
    /* convert to compound */
    cmpnd0 = WlzRGBAToCompound(obj0, WLZ_RGBA_SPACE_RGB, &errNum);
    cmpnd1 = WlzRGBAToCompound(obj1, WLZ_RGBA_SPACE_RGB, &errNum);

    /* operate on each channel */
    for(i=0; (i < 4) && (errNum == WLZ_ERR_NONE); i++){
      objs[i] = WlzImageArithmetic(cmpnd0->o[i], cmpnd1->o[i],
				   op, 0, &errNum);
    }

    /* create return compound object */
    if( errNum == WLZ_ERR_NONE ){
     rtnCmpnd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 4, &(objs[0]),
				     obj0->type, &errNum);
    }
    else {
      rtnCmpnd = NULL;
    }
    
    /* free temporary objects */
    WlzFreeObj((WlzObject *) cmpnd0);
    WlzFreeObj((WlzObject *) cmpnd1);
  }

  /* convert back to RGBA  */
  if( errNum == WLZ_ERR_NONE ){
    rtnObj = WlzCompoundToRGBA(rtnCmpnd, WLZ_RGBA_SPACE_RGB,
			       0, &errNum);
  }
  if( rtnCmpnd ){
    WlzFreeObj((WlzObject *) rtnCmpnd);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
