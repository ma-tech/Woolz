#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExplode3D.c
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
* \brief	Explodes a 3D domain object into 2D domain objects.
* \ingroup	WlzSectionTransform
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Error number.
* \ingroup	WlzSectionTransform
* \brief	Explodes the given 3D domain object into 2D domain objects.
* \param	dstExpObjCount		Destination pointer for number
*					of exploded Woolz objects.
* \param	dstExpObjVecP		Destination pointer for vector of
* 					exploded objects.
* \param	srcObj			Given woolz object.
*/
WlzErrorNum	WlzExplode3D(int *dstExpObjCount,
			     WlzObject ***dstExpObjVecP,
			     WlzObject *srcObj)
{
  int		objIdx,
  		objCount = 0;
  WlzDomain	srcDom;
  WlzValues	srcVal;
  WlzObject	**objVec = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzExplode3D FE 0x%lx  0x%lx 0x%lx\n",
	   (unsigned long )dstExpObjCount,
	   (unsigned long )dstExpObjVecP,
	   (unsigned long )srcObj));
  if((dstExpObjVecP == NULL) || (dstExpObjCount == NULL) ||(srcObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(((srcDom = srcObj->domain).core) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
    
  }
  /*
  else if(((srcVal = srcObj->values).core) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  */
  else if((objCount = srcDom.p->lastpl - srcDom.p->plane1 + 1) >= 1)
  {
    if((objVec = (WlzObject **)AlcMalloc((unsigned long )objCount *
                                         sizeof(WlzObject *))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      objIdx = 0;
      while((objIdx < objCount) && (errNum == WLZ_ERR_NONE))
      {
	if(((srcVal = srcObj->values).core) == NULL){
	  WlzValues values;
	  values.core = NULL;
	  *(objVec + objIdx) = WlzMakeMain(WLZ_2D_DOMAINOBJ,
					   *(srcDom.p->domains + objIdx),
					   values,
					   NULL, NULL, &errNum);
	}
	else {
	  *(objVec + objIdx) = WlzMakeMain(WLZ_2D_DOMAINOBJ,
					   *(srcDom.p->domains + objIdx),
					   *(srcVal.vox->values + objIdx),
					   NULL, NULL, &errNum);
	}
        ++objIdx;
      }
      if(errNum != WLZ_ERR_NONE)
      {
	while((objIdx < objCount) && *(objVec + objIdx))
	{
	  WlzFreeObj(*(objVec + objIdx));
	}
        AlcFree(objVec);
	objCount = 0;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstExpObjVecP = objVec;
    *dstExpObjCount = objCount;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzExplode3D FX %d\n",
	   (int )errNum));
  return(errNum);
}
