#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExplode3D_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzExplode3D.c
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
  	  ("WlzExplode3D FE %p  %p %p\n",
	   dstExpObjCount, dstExpObjVecP, srcObj));
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
    if((objVec = (WlzObject **)AlcMalloc((size_t )objCount *
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
