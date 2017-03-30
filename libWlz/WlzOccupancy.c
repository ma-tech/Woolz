#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzOccupancy_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzOccupancy.c
* \author       Bill Hill
* \date         November 2012
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
* \brief	Functions for computing the occupancy of domains.
* \ingroup	WlzFeatures
*/

#include <limits.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes an array of integers which correspond to a
* 		line through the given object's domain perpendicular to the
* 		plane defined by the given view structure. At each point
* 		along the line the occupancy area or number of domains
* 		intersecting the plane is computed.
*		When the given object is a 3D spatial domain object
*		then areas are computed and when it is a compound array
* 		then the occupancy (number of domains) is computed.
* \param	obj			Given object with the domain
* 					or a compound object with many
* 					domains.
* \param	vs			Given view structure which defines the
* 					cutting plane. The view structure
* 					will be initialised within this
* 					function.
* \param	sep			Plane separation distance, must be
* 					greater than ALG_DBL_TOLLERANCE.
* \param	dstFirst		Destination pointer for the first
* 					coordinate of the line in a 1D
* 					coordinate system perpendicular
* 					to the cutting plane. May be NULL.
* \param	dstLast			Destination pointer for the last
* 					coordinate of the line in a 1D
* 					coordinate system perpendicular 
* 					to the cutting plane. May be NULL.
* \param	dstArraySizeOcc		Destination pointer for the size of the
* 					occupancy array. May be NULL.
* \param	dstArrayOcc		Destination pointer for the occupancy
* 					array. The occupancy array should be
* 					freed using AlcFree(). May be NULL.
*/
WlzErrorNum	Wlz3DSectionOcc(WlzObject *obj, WlzThreeDViewStruct *vs,
				double sep,
				double *dstFirst, double *dstLast,
				int *dstArraySizeOcc, int **dstArrayOcc)
{
  int		nOcc = 0;
  double	first = 0.0,
  		last = 0.0,
		eps = ALG_DBL_TOLLERANCE;
  int		*occ = NULL;
  WlzObject	*tObj = NULL;
  WlzCompoundArray *cObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_3D_DOMAINOBJ:
	tObj = obj;
        break;
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
	cObj = (WlzCompoundArray *)obj;
	tObj = (cObj->o)? cObj->o[0]: NULL;
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(tObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    if(tObj->type != WLZ_3D_DOMAINOBJ)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if(tObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else if(vs == NULL)
    {
      errNum = WLZ_ERR_PARAM_NULL;
    }
    else if(vs->type != WLZ_3D_VIEW_STRUCT)
    {
      errNum = WLZ_ERR_PARAM_TYPE;
    }
    else if(sep < eps)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(cObj)
    {
      WlzObject *uObj;

      uObj = WlzAssignObject(
      	     WlzUnionN(cObj->n, cObj->o, 0, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzInit3DViewStruct(vs, uObj);
      }
      (void )WlzFreeObj(uObj);
      vs->ref_obj = NULL;
    }
    else
    {
      errNum = WlzInit3DViewStruct(vs, obj);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    first = vs->minvals.vtZ;
    last = vs->maxvals.vtZ;
    eps *= (fabs(first) + fabs(last));
    nOcc = floor((last - first + sep + eps) / sep);
    if(nOcc < 0)
    {
      errNum = WLZ_ERR_DOUBLE_DATA;
    }
    else if((occ = (int *)AlcCalloc(nOcc, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idP;

    for(idP = 0; idP < nOcc; ++idP)
    {
      WlzObject	*obj2D = NULL;
      const WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;

      vs->dist = first + (idP * sep);
      if(cObj)
      {
        if((first < (vs->dist + eps)) && (vs->dist < (last + eps)))
	{
	  int	idN;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	  for(idN = 0; idN < cObj->n; ++idN)
	  {
	    WlzThreeDViewStruct *vs1;
	    WlzErrorNum errNum2 = WLZ_ERR_NONE;

	    if(errNum == WLZ_ERR_NONE)
	    {
	      vs1 = WlzMake3DViewStructCopy(vs, &errNum);
	      if(errNum2 == WLZ_ERR_NONE)
	      {
		errNum2 = WlzInit3DViewStruct(vs1, cObj->o[idN]);
	      }
	      if(errNum2 == WLZ_ERR_NONE)
	      {
		obj2D = WlzAssignObject(
			WlzGetSubSectionFromObject(cObj->o[idN], NULL,
						   vs1, interp,
						   NULL, &errNum2), NULL);
	      }
	      (void )WlzFree3DViewStruct(vs1);
	      if(errNum2 == WLZ_ERR_NONE)
	      {
		occ[idP] += (WlzIsEmpty(obj2D, &errNum2))? 0: 1;
	      }
	    }
	    if(errNum2 != WLZ_ERR_NONE)
	    {
#ifdef _OPENMP
#pragma omp critical (Wlz3DSectionOcc)
	      {
		if(errNum == WLZ_ERR_NONE)
		{
		  errNum = errNum2;
		}
	      }
#else
	      errNum = errNum2;
#endif
	    }
	  }
	}
      }
      else
      {
	errNum = WlzInit3DViewStruct(vs, obj);
	if(errNum == WLZ_ERR_NONE)
	{
	  obj2D = WlzAssignObject(
		  WlzGetSubSectionFromObject(obj, NULL, vs, interp,
					     NULL, &errNum), NULL);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  occ[idP] = WlzArea(obj2D, &errNum);
	}
	(void )WlzFreeObj(obj2D);
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstFirst)
    {
      *dstFirst = first;
    }
    if(dstLast)
    {
      *dstLast = last;
    }
    if(dstArraySizeOcc)
    {
      *dstArraySizeOcc = nOcc;
    }
    if(dstArrayOcc)
    {
      *dstArrayOcc = occ;
    }
  }
  else
  {
    AlcFree(occ);
    occ = NULL;
  }
  return(errNum);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzFeatures
* \brief	Computes a new Woolz domain object which has the domain
* 		of the given object (union of domains if the object is a
* 		compound array) and a value table with integer values
* 		representing the number of domains present.
* \param	gObj			Given object.
* \param	dstErr			Destination error pointer.
*/
WlzObject	*WlzDomainOccupancy(WlzObject *gObj, WlzErrorNum *dstErr)
{
  WlzObject 	*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    WlzPixelV	gV,
    		zeroV;
    WlzObjectType tType;

    zeroV.v.ubv = 0;
    zeroV.type = gV.type = WLZ_GREY_UBYTE;
    tType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_UBYTE, NULL);
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	/* Easy for a simple spatial domain object, just create a new
	 * object with the same domain and all values in it set to 1. */
	gV.v.inv = 1;
	oObj = WlzNewObjectValues(gObj, tType, zeroV, 1, gV, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzGreySetValue(oObj, zeroV);
	}
        break;
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
	{
	  int		idN;
	  WlzObject	*uObj = NULL;
	  WlzCompoundArray *cObj;

	  cObj = (WlzCompoundArray *)gObj;
	  /* Compute the union of the compound object's domains and set it's
	   * values to 0. */
	  if(cObj->n > SHRT_MAX)
	  {
	    zeroV.v.inv = 0;
	    zeroV.type = gV.type = WLZ_GREY_INT;
	    tType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT,
	                                  NULL);
	  }
	  else if(cObj->n > 255)
	  {
	    zeroV.v.shv = 0;
	    zeroV.type = gV.type = WLZ_GREY_SHORT;
	    tType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_SHORT,
	    				  NULL);
	  }
	  uObj = WlzUnionN(cObj->n, cObj->o, 0, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    oObj = WlzNewObjectValues(uObj, tType, zeroV, 1, gV, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzGreySetValue(oObj, zeroV);
	  }
	  (void )WlzFreeObj(uObj);
	  /* For each object of the compound object, increment all values
	   * of the occupancy object which fall within it's domain. */
	  for(idN = 0; idN < cObj->n; ++idN)
	  {
	    errNum = WlzGreyIncValuesInDomain(oObj, cObj->o[idN]);
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    (void )WlzFreeObj(oObj);
	    oObj = NULL;
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oObj);
}
