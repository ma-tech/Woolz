#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzNObjGreyStats_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzNObjGreyStats.c
* \author       Bill Hill
* \date         March 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Calculates statistics for grey values across the
* 		objects of the given compound object.
* \ingroup	WlzFeatures
*/


#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static WlzErrorNum 		WlzNObjGreyStats3D(
				  int nObj,
				  WlzObject **aObj,
				  WlzObject *minObj,
				  WlzObject *maxObj,
				  WlzObject *sumObj,
				  WlzObject *sSqObj);
static WlzErrorNum 		WlzNObjGreyStats2D(
				  int first,
				  int nObj,
				  WlzObject **aObj,
				  WlzObject *minObj,
				  WlzObject *maxObj,
				  WlzObject *sumObj,
				  WlzObject *sSqObj);
static WlzErrorNum		WlzNObjGreyStats2D1(
				  int first,
				  WlzIntervalWSpace *gIWSp,
				  WlzGreyWSpace *gGWSp,
				  WlzObject *minObj,
				  WlzIntervalWSpace *minIWSp,
				  WlzGreyWSpace *minGWSp,
				  WlzObject *maxObj,
				  WlzIntervalWSpace *maxIWSp,
				  WlzGreyWSpace *maxGWSp,
				  WlzObject *sumObj,
				  WlzIntervalWSpace *sumIWSp,
				  WlzGreyWSpace *sumGWSp,
				  WlzObject *sSqObj,
				  WlzIntervalWSpace *sSqIWSp,
				  WlzGreyWSpace *sSqGWSp);

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes a collection of objects from the object (which
* 		should be a compound array object).
* 		The computed objects all have the intersection of the
* 		given object domains as their domain and WLZ_GREY_DOUBLE
* 		values for the minimum, maximum, sum and sum of squares
* 		of the given grey values at each pixel/voxel within
* 		the intersection domain.
* 		Returned objects will have name properties set to: min,
* 		max, sum, ssq, mean and stddev as appropriate.
* \param	gObj			Given object.
* \param	mean			Compute mean not sum if non-zero.
* \param	stddev			Compute standard deviation not sum of
* 					squares if non-zero.
* \param	dstN			Destination pointer for number of
* 					objects, may be NULL.
* \param	dstMinObj		Destination pointer for minimum value
* 					object, may be NULL if not required.
* \param	dstMaxObj		Destination pointer for maximum value
* 					object, may be NULL if not required.
* \param	dstSumObj		Destination pointer for sum of values
* 					object, may be NULL if not required.
* \param	dstSSqObj		Destination pointer for sum of squares
* 					object, may be NULL if not required.
*/
WlzErrorNum	WlzNObjGreyStats(WlzObject *gObj,
				 int mean, int stddev, int *dstN,
				 WlzObject **dstMinObj,
				 WlzObject **dstMaxObj,
				 WlzObject **dstSumObj,
				 WlzObject **dstSSqObj)
{
  int		nDim = 0,
  		nObj = 0;
  WlzObject	**aObj = NULL;
  WlzObject	*isnObj = NULL,
  		*minObj = NULL,
		*maxObj = NULL,
		*sumObj = NULL,
		*sSqObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idN;

    nObj = ((WlzCompoundArray *)gObj)->n;
    aObj = ((WlzCompoundArray *)gObj)->o;
    for(idN = 0; idN < nObj; ++idN)
    {
      WlzObject	*obj;

      obj = aObj[idN];
      if(obj == NULL)
      {
        errNum = WLZ_ERR_OBJECT_NULL;
	break;
      }
      else if((obj->type != WLZ_2D_DOMAINOBJ) &&
              (obj->type != WLZ_3D_DOMAINOBJ))
      {
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
      else if(obj->domain.core == NULL)
      {
        errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      else if(obj->values.core == NULL)
      {
        errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      else if(WlzGreyTableIsTiled(obj->values.core->type))
      {
        errNum = WLZ_ERR_VALUES_TYPE;
	break;
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) &&
     ((dstMinObj != NULL) || (dstMaxObj != NULL) ||
      (dstSumObj != NULL) || (dstSSqObj != NULL)))
  {
    isnObj = WlzIntersectN(nObj, aObj, 0, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(isnObj->type)
      {
	case WLZ_2D_DOMAINOBJ:
	  nDim = 2;
	  break;
	case WLZ_3D_DOMAINOBJ:
	  nDim = 3;
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }

    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzPixelV	bgdV;
      WlzObjectType gTType;

      bgdV.type = WLZ_GREY_DOUBLE;
      bgdV.v.dbv = 0.0;
      gTType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE,
                                     &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	if(dstMinObj != NULL)
	{
	  minObj = WlzNewObjectValues(isnObj, gTType, bgdV, 0, bgdV, &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && (dstMaxObj != NULL))
	{
	  maxObj = WlzNewObjectValues(isnObj, gTType, bgdV, 0, bgdV, &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && (dstSumObj != NULL))
	{
	  sumObj = WlzNewObjectValues(isnObj, gTType, bgdV, 0, bgdV, &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && (dstSSqObj != NULL))
	{
	  sSqObj = WlzNewObjectValues(isnObj, gTType, bgdV, 0, bgdV, &errNum);
	}
      }
    }
    (void )WlzFreeObj(isnObj);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = (nDim == 2)?
	       WlzNObjGreyStats2D(1, nObj, aObj,
	                          minObj, maxObj, sumObj, sSqObj):
	       WlzNObjGreyStats3D(nObj, aObj,
	                          minObj, maxObj, sumObj, sSqObj);
    }
    if((errNum == WLZ_ERR_NONE) && (stddev != 0))
    {
      WlzPixelV tV;

      tV.type = WLZ_GREY_DOUBLE;
      if(nObj <= 1)
      {
	tV.v.dbv = 1.0;
        errNum = WlzGreySetValue(sSqObj, tV);
      }
      else
      {
	WlzObject *tSS = NULL,
		  *tSSN = NULL,
		  *tSSS = NULL,
		  *tSSSN = NULL,
		  *tSqrt = NULL;

	/* stdDev = sqrt((sSq - (sum * sum / n)) / (n - 1)) */
	tV.v.dbv = 1.0 / nObj;
	tSS = WlzScalarFn(sumObj, WLZ_FN_SCALAR_SQR, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  tSSN = WlzScalarMultiply(tSS, tV, &errNum);
	}
	(void )WlzFreeObj(tSS);
	if(errNum == WLZ_ERR_NONE)
	{
	  tSSS = WlzImageArithmetic(sSqObj, tSSN, WLZ_BO_SUBTRACT, 0, &errNum);
	}
	(void )WlzFreeObj(tSSN);
	if(errNum == WLZ_ERR_NONE)
	{
	  tV.v.dbv = 1.0 / (nObj - 1.0);
	  tSSSN = WlzScalarMultiply(tSSS, tV, &errNum);
	}
	(void )WlzFreeObj(tSSS);
	if(errNum == WLZ_ERR_NONE)
	{
	  tSqrt = WlzScalarFn(tSSSN, WLZ_FN_SCALAR_SQRT, &errNum);
	}
	(void )WlzFreeObj(tSSSN);
	if(errNum == WLZ_ERR_NONE)
	{
	  (void )WlzFreeObj(sSqObj);
	  sSqObj = tSqrt;
	  tSqrt = NULL;
	}
	(void )WlzFreeObj(tSqrt);
      }
    }
    if((errNum == WLZ_ERR_NONE) && (mean != 0))
    {
      WlzPixelV tV;
      WlzObject *tObj;

      tV.type = WLZ_GREY_DOUBLE;
      tV.v.dbv = 1.0 / nObj;
      tObj = WlzScalarMultiply(sumObj, tV, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        (void )WlzFreeObj(sumObj);
	sumObj = tObj;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      int	idJ;
      WlzObject	*objs[4] = {NULL};
      char 	*str[4];

      if(minObj)
      {
        objs[0] = minObj;
	str[0] = "min";
      }
      if(maxObj)
      {
        objs[1] = maxObj;
	str[1] = "max";
      }
      if(sumObj)
      {
        objs[2] = sumObj;
	str[2] = (mean)? "mean": "sum";
      }
      if(sSqObj)
      {
        objs[3] = sSqObj;
	str[3] = (stddev)? "stddev": "ssq";
      }
      for(idJ = 0; (idJ < 4) && (errNum == WLZ_ERR_NONE); ++idJ)
      {
        WlzProperty prop;
	WlzObject *obj;
        
	if((obj = objs[idJ]) != NULL)
	{
	  prop.name = WlzMakeNameProperty(str[idJ], &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzPropertyList *pList = NULL;

	    if(((pList = WlzMakePropertyList(NULL)) == NULL) ||
	       (AlcDLPListEntryAppend(pList->list, NULL,
			      (void *)(prop.core),
			      WlzFreePropertyListEntry) != ALC_ER_NONE))
	    {
	      if(prop.core)
	      {
	        (void )WlzFreeProperty(prop);
	      }
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    else
	    {
	      obj->plist = WlzAssignPropertyList(pList, &errNum);
	    }
	  }
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreeObj(minObj);
      (void )WlzFreeObj(maxObj);
      (void )WlzFreeObj(sumObj);
      (void )WlzFreeObj(sSqObj);
    }
    else
    {
      if(dstN)
      {
	*dstN = nObj;
      }
      if(dstMinObj)
      {
	*dstMinObj = minObj;
      }
      if(dstMaxObj)
      {
	*dstMaxObj = maxObj;
      }
      if(dstSumObj)
      {
	*dstSumObj = sumObj;
      }
      if(dstSSqObj)
      {
	*dstSSqObj = sSqObj;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes the minimum, maximum, sum and sum of squares of
* 		all objects in the given compound array object, filling
* 		in the values of these objects. All objects are allocated
* 		before calling this function and all the objects to have
* 		their values filled in share the same domain which is the
* 		intersection of the compound array object domains. All the
* 		object to be filled in have WLZ_GREY_DOUBLE values.
* \param	nObj			Number of given objects.
* \param	aObj			Array of given objects.
* \param	minObj			Minimum value object, may be NULL.
* \param	maxObj			Maximum value object, may be NULL.
* \param	sumObj			Sum of values object, may be NULL.
* \param	sSqObj			Sum of squares object, may be NULL.
*/
static WlzErrorNum WlzNObjGreyStats3D(int nObj, WlzObject **aObj,
				      WlzObject *minObj, WlzObject *maxObj,
				      WlzObject *sumObj, WlzObject *sSqObj)
{
  int		idN,
  		tPlnCnt;
  WlzObject	*tObj = NULL;
  WlzPlaneDomain *tPDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tObj = minObj;
  if(tObj == NULL)
  {
    tObj = maxObj;
  }
  if(tObj == NULL)
  {
    tObj = sumObj;
  }
  if(tObj == NULL)
  {
    tObj = sSqObj;
  }
  tPDom = tObj->domain.p;
  tPlnCnt = tPDom->lastpl - tPDom->plane1 + 1;
  for(idN = 0; idN < nObj; ++idN)
  {
    int		idP,
    		gPlnCnt;
    WlzObject	*gObj;
    WlzPlaneDomain *gPDom;

    gObj = aObj[idN];
    gPDom = gObj->domain.p;
    gPlnCnt = gPDom->lastpl - gPDom->plane1 + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idP = 0; idP < gPlnCnt; ++idP)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	int	idT;
	WlzErrorNum errNum2 = WLZ_ERR_NONE;

	idT = gPDom->plane1 + idP - tPDom->plane1;
	if((idT >= 0) && (idT < tPlnCnt) &&
	   ((gPDom->domains[idP]).core != NULL))
	{
	  WlzObject *gObj2 = NULL,
		    *minObj2 = NULL,
		    *maxObj2 = NULL,
		    *sumObj2 = NULL,
		    *sSqObj2 = NULL;

	  gObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			      gPDom->domains[idP],
			      gObj->values.vox->values[idP],
			      NULL, NULL, &errNum2);
	  if((errNum2 == WLZ_ERR_NONE) && (minObj != NULL))
	  {
	    minObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				  minObj->domain.p->domains[idT], 
				  minObj->values.vox->values[idT], 
				  NULL, NULL, &errNum2);
	  }
	  if((errNum2 == WLZ_ERR_NONE) && (maxObj != NULL))
	  {
	    maxObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				  maxObj->domain.p->domains[idT], 
				  maxObj->values.vox->values[idT], 
				  NULL, NULL, &errNum2);
	  }
	  if((errNum2 == WLZ_ERR_NONE) && (sumObj != NULL))
	  {
	    sumObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				  sumObj->domain.p->domains[idT], 
				  sumObj->values.vox->values[idT], 
				  NULL, NULL, &errNum2);
	  }
	  if((errNum2 == WLZ_ERR_NONE) && (sSqObj != NULL))
	  {
	    sSqObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				  sSqObj->domain.p->domains[idT], 
				  sSqObj->values.vox->values[idT], 
				  NULL, NULL, &errNum2);
	  }
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    int		first;

	    first = (idN == 0);
	    errNum2 =  WlzNObjGreyStats2D(first, 1, &gObj2,
					minObj2, maxObj2,
					sumObj2, sSqObj2);
	  }
	}
	if(errNum2 != WLZ_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp critical (WlzNObjGreyStats3D)
#endif
	  {
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = errNum2;
	    }
	  }
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes the minimum, maximum, sum and sum of squares of
* 		all objects in the given compound array object, filling
* 		in the values of these objects. All objects are allocated
* 		before calling this function and all the objects to have
* 		their values filled in share the same domain which is the
* 		intersection of the compound array object domains. All the
* 		object to be filled in have WLZ_GREY_DOUBLE values.
* \param	first			Non-zero if this is the first
* 					call for this plane of values.
* \param	nObj			Number of given objects.
* \param	aObj			Array of given objects.
* \param	minObj			Minimum value object, may be NULL.
* \param	maxObj			Maximum value object, may be NULL.
* \param	sumObj			Sum of values object, may be NULL.
* \param	sSqObj			Sum of squares object, may be NULL.
*/
static WlzErrorNum WlzNObjGreyStats2D(int first, int nObj, WlzObject **aObj,
				      WlzObject *minObj, WlzObject *maxObj,
				      WlzObject *sumObj, WlzObject *sSqObj)
{
  int		idN;
  WlzIntervalWSpace iIWSp = {0},
  		minIWSp = {0},
  		maxIWSp = {0},
		sumIWSp = {0},
		sSqIWSp = {0};
  WlzGreyWSpace	iGWSp,
  		minGWSp,
  		maxGWSp,
		sumGWSp,
		sSqGWSp;
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  tObj = minObj;
  if(tObj == NULL)
  {
    tObj = maxObj;
  }
  if(tObj == NULL)
  {
    tObj = sumObj;
  }
  if(tObj == NULL)
  {
    tObj = sSqObj;
  }
  for(idN = 0; idN < nObj; ++idN)
  {
    WlzObject	*iObj = NULL;

    iObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, tObj->domain, aObj[idN]->values,
    		       NULL, NULL, &errNum);
    if((errNum == WLZ_ERR_NONE) &&
       (iObj != NULL) && (WlzIsEmpty(iObj, NULL) == 0))
    {
      errNum = WlzInitGreyScan(iObj, &iIWSp, &iGWSp);
      if((errNum == WLZ_ERR_NONE) && (minObj != NULL))
      {
	errNum = WlzInitGreyScan(minObj, &minIWSp, &minGWSp);
      }
      if((errNum == WLZ_ERR_NONE) && (maxObj != NULL))
      {
	errNum = WlzInitGreyScan(maxObj, &maxIWSp, &maxGWSp);
      }
      if((errNum == WLZ_ERR_NONE) && (sumObj != NULL))
      {
	errNum = WlzInitGreyScan(sumObj, &sumIWSp, &sumGWSp);
      }
      if((errNum == WLZ_ERR_NONE) && (sSqObj != NULL))
      {
	errNum = WlzInitGreyScan(sSqObj, &sSqIWSp, &sSqGWSp);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzNObjGreyStats2D1(first, &iIWSp, &iGWSp,
				     minObj, &minIWSp, &minGWSp,
				     maxObj, &maxIWSp, &maxGWSp,
				     sumObj, &sumIWSp, &sumGWSp,
				     sSqObj, &sSqIWSp, &sSqGWSp);
        first = 0;
      }
      (void )WlzEndGreyScan(&iIWSp,   &iGWSp);
      (void )WlzEndGreyScan(&minIWSp, &minGWSp);
      (void )WlzEndGreyScan(&maxIWSp, &maxGWSp);
      (void )WlzEndGreyScan(&sumIWSp, &sumGWSp);
      (void )WlzEndGreyScan(&sSqIWSp, &sSqGWSp);
    }
    (void )WlzFreeObj(iObj);
    if(errNum != WLZ_ERR_NONE)
    {
      break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Sets the value of square of value in the minimum, maximum,
* 		sum and sum of squares objects from a given object using
* 		the given interval and grey workspace pointers. All the
* 		workspace pointers are known to be valid and intialised.
* 		All object domains are known to be the same.
* \param	first			Non-zero if this is the first
* 					object and values are simply to
* 					be set.
* \param 	gIWSp			Sum of squares interval workspace.
* \param	gGWSp			Given object grey value workspace.
* \param	minObj			Minimum value object, may be NULL.
* \param	minIWSp			Minimum value object interval
* 					workspace.
* \param	minGWSp			Minimum value object grey value
* 					workspace.
* \param	maxObj			Maximum value object, may be NULL.
* \param	maxIWSp			Maximum value object interval
* 					workspace.
* \param	maxGWSp			Maximum value object grey value
* 					workspace.
* \param	sumObj			Sum of values object, may be NULL.
* \param	sumIWSp			Sum of values object interval
* 					workspace.
* \param	sumGWSp			Sum of values object grey value
* 					workspace.
* \param	sSqObj			Sum of squares object, may be NULL.
* \param	sSqIWSp			Sum of squares object interval
* 					workspace.
* \param	sSqGWSp			Sum of squares object grey value
* 					workspace.
*/
static WlzErrorNum WlzNObjGreyStats2D1(int first,
				    WlzIntervalWSpace *gIWSp,
				    WlzGreyWSpace *gGWSp,
				    WlzObject *minObj,
				    WlzIntervalWSpace *minIWSp,
				    WlzGreyWSpace *minGWSp,
				    WlzObject *maxObj,
				    WlzIntervalWSpace *maxIWSp,
				    WlzGreyWSpace *maxGWSp,
				    WlzObject *sumObj,
				    WlzIntervalWSpace *sumIWSp,
				    WlzGreyWSpace *sumGWSp,
				    WlzObject *sSqObj,
				    WlzIntervalWSpace *sSqIWSp,
				    WlzGreyWSpace *sSqGWSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) &&
        ((errNum = WlzNextGreyInterval(gIWSp)) == WLZ_ERR_NONE) &&
        ((minObj == NULL) || (WlzNextGreyInterval(minIWSp) == 0)) &&
        ((maxObj == NULL) || (WlzNextGreyInterval(maxIWSp) == 0)) &&
        ((sumObj == NULL) || (WlzNextGreyInterval(sumIWSp) == 0)) &&
        ((sSqObj == NULL) || (WlzNextGreyInterval(sSqIWSp) == 0)))
  {
    int		idK,
    		itvCnt;
    WlzGreyP	gPix;

    gPix = gGWSp->u_grintptr;
    itvCnt = gIWSp->rgtpos - gIWSp->lftpos + 1;
    switch(gGWSp->pixeltype)
    {
      case WLZ_GREY_UBYTE:
	if(first)
	{
	  if(minObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      minGWSp->u_grintptr.dbp[idK] = gPix.ubp[idK];
	    }
	  }
	  if(maxObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      maxGWSp->u_grintptr.dbp[idK] = gPix.ubp[idK];
	    }
	  }
	  if(sumObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sumGWSp->u_grintptr.dbp[idK] = gPix.ubp[idK];
	    }
	  }
	  if(sSqObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sSqGWSp->u_grintptr.dbp[idK] = (double )(gPix.ubp[idK]) *
	                                    (double )(gPix.ubp[idK]);
	    }
	  }
	}
	else
	{
	  if(minObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      if(minGWSp->u_grintptr.dbp[idK] > gPix.ubp[idK])
	      {
	        minGWSp->u_grintptr.dbp[idK] = gPix.ubp[idK];
	      }
	    }
	  }
	  if(maxObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      if(maxGWSp->u_grintptr.dbp[idK] < gPix.ubp[idK])
	      {
	        maxGWSp->u_grintptr.dbp[idK] = gPix.ubp[idK];
	      }
	    }
	  }
	  if(sumObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sumGWSp->u_grintptr.dbp[idK] += gPix.ubp[idK];
	    }
	  }
	  if(sSqObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sSqGWSp->u_grintptr.dbp[idK] += (double )(gPix.ubp[idK]) *
	                                     (double )(gPix.ubp[idK]);
	    }
	  }
	}
	break;
      case WLZ_GREY_DOUBLE:
	if(first)
	{
	  if(minObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      minGWSp->u_grintptr.dbp[idK] = gPix.dbp[idK];
	    }
	  }
	  if(maxObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      maxGWSp->u_grintptr.dbp[idK] = gPix.dbp[idK];
	    }
	  }
	  if(sumObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sumGWSp->u_grintptr.dbp[idK] = gPix.dbp[idK];
	    }
	  }
	  if(sSqObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sSqGWSp->u_grintptr.dbp[idK] = (double )(gPix.dbp[idK]) *
	                                    (double )(gPix.dbp[idK]);
	    }
	  }
	}
	else
	{
	  if(minObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      if(minGWSp->u_grintptr.dbp[idK] > gPix.dbp[idK])
	      {
	        minGWSp->u_grintptr.dbp[idK] = gPix.dbp[idK];
	      }
	    }
	  }
	  if(maxObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      if(maxGWSp->u_grintptr.dbp[idK] < gPix.dbp[idK])
	      {
	        maxGWSp->u_grintptr.dbp[idK] = gPix.dbp[idK];
	      }
	    }
	  }
	  if(sumObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sumGWSp->u_grintptr.dbp[idK] += gPix.dbp[idK];
	    }
	  }
	  if(sSqObj)
	  {
	    for(idK = 0; idK < itvCnt; ++idK)
	    {
	      sSqGWSp->u_grintptr.dbp[idK] += (double )(gPix.dbp[idK]) *
	                                     (double )(gPix.dbp[idK]);
	    }
	  }
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}
