#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzIndexObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzIndexObj.c
* \author       Bill Hill
* \date         June 2011
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \brief	Functions for creating and manipulating objects in which
* 		2 or 3D spatial domains are represented as index grey values
* 		within in a single domain object. Such objects are frequently
* 		refered to as index objects.
* \ingroup	WlzValueUtils
*/
#include <limits.h>
#include <Wlz.h>

/*!
* \return	New spatial domain object with grey values or NULL on error.
* \ingroup	WlzValueUtils
* \brief	Creates a new spatial domain object with grey values;
*  		the grey values being indices of the domains in the given
*  		compound array object.
*
* 		The domain of the new object is the union of the domains
* 		of the objects in the given compound array object. The
* 		values are set so that all values in the domain of the
* 		i'th compound array object are set to i, with i
* 		incremented from 0 through to one less than the number
* 		of objects in the compound array object. If the domains
* 		overlap then higher index objects will overwrite lower
* 		index objects within their intersection.
* 		All objects of the compound array must either be of
* 		type WLZ_EMPTY, WLZ_2D_DOMAINOBJ or WLZ_2D_DOMAINOBJ
* 		and at least one must be non empty. Any NULL pointer
* 		in the compound array will be treated as an empty object.
* 		This function can be considered the inverse of
* 		WlzIndexObjToCompound().
* \param	cObj			Given compound array which must
* 					have all ojects valid.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzIndexObjFromCompound(WlzCompoundArray *cObj,
				         WlzErrorNum *dstErr)
{
  int		nNE = 0;
  WlzObject	*rObj = NULL,
  		*uObj = NULL;
  WlzObject	**uAry = NULL;
  WlzObjectType	gObjType = WLZ_EMPTY_OBJ;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(((cObj->type != WLZ_COMPOUND_ARR_1) &&
           (cObj->type != WLZ_COMPOUND_ARR_2)) ||
          (cObj->n < 1))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((uAry = (WlzObject **)
                  AlcMalloc(sizeof(WlzObject *) * cObj->n)) == NULL)
  {
      errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    int		idx;

    for(idx = 0; (errNum == WLZ_ERR_NONE) && (idx < cObj->n); ++idx)
    {
      if(cObj->o[idx] != NULL)
      {
	switch(cObj->o[idx]->type)
	{
	  case WLZ_EMPTY_OBJ:
	    break;
	  case WLZ_2D_DOMAINOBJ:
	    if(gObjType == WLZ_EMPTY_OBJ)
	    {
	      gObjType = WLZ_2D_DOMAINOBJ;
	    }
	    else if(gObjType != WLZ_2D_DOMAINOBJ)
	    {
	      errNum = WLZ_ERR_OBJECT_TYPE;
	    }
	    uAry[nNE++] = cObj->o[idx];
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    if(gObjType == WLZ_EMPTY_OBJ)
	    {
	      gObjType = WLZ_3D_DOMAINOBJ;
	    }
	    else if(gObjType != WLZ_3D_DOMAINOBJ)
	    {
	      errNum = WLZ_ERR_OBJECT_TYPE;
	    }
	    uAry[nNE++] = cObj->o[idx];
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    uObj = WlzUnionN(nNE, uAry, 0, &errNum);
  }
  AlcFree(uAry);
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues val;
    WlzObjectType gTT;
    WlzPixelV bgV;

    val.core = NULL;
    if(cObj->n <= 255)
    {
      bgV.v.ubv = 0;     
      bgV.type = WLZ_GREY_UBYTE;
    }
    else if(cObj->n <= SHRT_MAX)
    {
      bgV.v.shv = 0;     
      bgV.type = WLZ_GREY_SHORT;
    }
    else
    {
      bgV.v.inv = 0;     
      bgV.type = WLZ_GREY_INT;
    }
    gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, bgV.type, NULL);
    if(gObjType == WLZ_2D_DOMAINOBJ)
    {
      val.v = WlzNewValueTb(uObj, gTT, bgV, &errNum);
      
    }
    else /* gObjType == WLZ_3D_DOMAINOBJ */
    {
      val.vox = WlzNewValuesVox(uObj, gTT, bgV, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rObj = WlzMakeMain(gObjType, uObj->domain, val, NULL, NULL, &errNum);
    }
    else if(val.core != NULL)
    {
      if(gObjType == WLZ_2D_DOMAINOBJ)
      {
        (void )WlzFreeValueTb(val.v);
      }
      else
      {
        (void )WlzFreeVoxelValueTb(val.vox);
      }
    }
  }
  (void )WlzFreeObj(uObj);
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx;
    WlzPixelV	tVal;

    tVal.type = WLZ_GREY_INT;
    for(idx = 0; (errNum == WLZ_ERR_NONE) && (idx < cObj->n); ++idx)
    {
      if((cObj->o[idx] != NULL) && (cObj->o[idx]->type == gObjType))
      {
	WlzObject *tObj;

	tVal.v.inv = idx;
	(void )WlzAssignObject(rObj, NULL);
	tObj = WlzGreyMask(rObj, cObj->o[idx], tVal, &errNum);
	if(tObj != NULL)
	{
	  (void )WlzFreeObj(rObj);
	  rObj = tObj;
	}
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New compound array object or NULL on error.
* \ingroup	WlzValueUtils
* \brief	Creates a new compound array object in which each object
* 		of the array is either an empty object or a spatial domain
* 		without grey values. The spatial domain objects will be of
* 		the same type as the given object and their domains will
* 		correspond those of the grey values of the given object,
* 		where the i'th object of the compound object's array will
* 		have the domain for which all grey values have the value i.
* 		This function can be considered the inverse of
* 		WlzIndexObjFromCompound().
* \param	gObj			Given spatial domain object with
* 					grey values. The grey values must be
* 					of the type WLZ_GREY_UBYTE,
* 					WLZ_GREY_SHORT or WLZ_GREY_INT.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCompoundArray *WlzIndexObjToCompound(WlzObject *gObj, WlzErrorNum *dstErr)
{
  int		idxMin,
  		idxMax;
  WlzPixelV	minV,
  		maxV;
  WlzCompoundArray *cObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((gObj->type != WLZ_2D_DOMAINOBJ) && (gObj->type != WLZ_3D_DOMAINOBJ)) 
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    WlzGreyType	gType;

    gType = WlzGreyTypeFromObj(gObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(gType)
      {
	case WLZ_GREY_UBYTE: /* FALLTHROUGH */
	case WLZ_GREY_SHORT: /* FALLTHROUGH */
	case WLZ_GREY_INT:
	  errNum = WlzGreyRange(gObj,  &minV, &maxV);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValueConvertPixel(&minV, minV, WLZ_GREY_INT);
    WlzValueConvertPixel(&maxV, maxV, WLZ_GREY_INT);
    if(((idxMin = minV.v.inv) < 0) || ((idxMax = maxV.v.inv) < idxMin))
    {
      errNum = WLZ_ERR_GREY_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_2, 1, idxMax + 1, NULL,
                                gObj->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx;
    WlzPixelV	tV;

    tV.type = WLZ_GREY_INT;
    for(idx = 0; (errNum == WLZ_ERR_NONE) && (idx <= idxMax); ++idx)
    {
      if(idx < idxMin)
      {
        cObj->o[idx] = WlzAssignObject(
	               WlzMakeEmpty(&errNum), NULL);
      }
      else
      {
	WlzObject *hObj = NULL,
		  *lObj = NULL;

	tV.v.inv = idx;
	hObj = WlzAssignObject(
	       WlzThreshold(gObj, tV, WLZ_THRESH_HIGH, &errNum), NULL);
	if((errNum == WLZ_ERR_NONE) && (hObj != NULL))
	{
	  tV.v.inv = idx + 1;
	  lObj = WlzAssignObject(
		 WlzThreshold(hObj, tV, WLZ_THRESH_LOW, &errNum), NULL);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if((lObj != NULL) && (lObj->type != WLZ_EMPTY_OBJ))
	  {
	    if(WlzIsEmpty(lObj, NULL))
	    {
	      (void )WlzFreeObj(lObj);
	      lObj = NULL;
	    }
	  }
	  if(lObj == NULL)
	  {
	    lObj = WlzAssignObject(
	           WlzMakeEmpty(&errNum), NULL);
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  WlzValues nullVal;

	  nullVal.core = NULL;
	  cObj->o[idx] = WlzAssignObject(
	                 WlzMakeMain(lObj->type, lObj->domain, nullVal,
			             NULL, NULL, &errNum), NULL);
	}
	(void )WlzFreeObj(lObj);
	(void )WlzFreeObj(hObj);
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj((WlzObject *)cObj);
    cObj = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(cObj);
}
