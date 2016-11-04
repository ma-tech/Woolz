#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzValuesFromCoords_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzValuesFromCoords.c
* \author       Bill Hill
* \date         April 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Creates a new object with the values set to the
* 		coordinate values.
* \ingroup	WlzValuesFilters
*/

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <Wlz.h>

static WlzCompoundArray		*WlzValuesFromCoords2D(
				  WlzObject *gObj,
				  int radial,
				  WlzGreyType gType,
				  WlzErrorNum *dstErr);
static WlzCompoundArray		*WlzValuesFromCoords3D(
				  WlzObject *gObj,
				  int radial,
				  WlzGreyType gType,
				  WlzErrorNum *dstErr);

/*!
* \return	New compound object with an object per dimension or
* 		NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Given a spatial domain object (with or wthout values)
* 		this function creates a compound array object with the
* 		number of component objects set to the number of the
* 		dimensions of the input object. The domains of the
* 		component objects are those of the given object and
* 		the values are set to the coordinate values within
* 		the object domains. Components are ordered by column,
* 		line, plane or radial distance, angle (in radians), plane.
* \param	gObj			Given spatial domain object.
* \param	radial			Use cylindrical coordinates if
* 					non-zero.
* \param	gType			Required grey type; valid types
* 					are restricted to WLZ_GREY_INT
* 					and WLZ_GREY_DOUBLE any other
* 					grey type is an error.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzValuesFromCoords(
				  WlzObject *gObj,
				  int radial,
				  WlzGreyType gType,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzCompoundArray *cObj = NULL;
  WlzErrorNum  	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gType != WLZ_GREY_INT) && (gType != WLZ_GREY_DOUBLE))
  {
    errNum = WLZ_GREY_ERROR;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	cObj = WlzValuesFromCoords2D(gObj, radial, gType, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	cObj = WlzValuesFromCoords3D(gObj, radial, gType, &errNum);
        break;
      default:
       break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = (WlzObject *)cObj;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New compound object with an object per dimension or
* 		NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Given a 2D spatial domain object this function creates
* 		a compound array object with 2 component objects with
* 		the object domains those of the given object and the
* 		object values the coordinate values.
* \param	gObj			Given spatial domain object.
* \param	radial			Use cylindrical coordinates if
* 					non-zero.
* \param	gType			Required grey type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray		*WlzValuesFromCoords2D(
				  WlzObject *gObj,
				  int radial,
				  WlzGreyType gType,
				  WlzErrorNum *dstErr)
{
  int		idx;
  WlzPixelV	bgdV;
  WlzObjectType tt;
  WlzCompoundArray	*cpd;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(gType == WLZ_GREY_INT)
  {
    bgdV.v.inv = 0;
    bgdV.type = WLZ_GREY_INT;
    tt = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_INT, NULL);
  }
  else /* gType == WLZ_GREY_DOUBLE */
  {
    bgdV.v.dbv = 0.0;
    bgdV.type = WLZ_GREY_DOUBLE;
    tt = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_DOUBLE, NULL);
  }
  cpd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, 2, NULL, gObj->type,
                             &errNum);
  for(idx = 0; (errNum == WLZ_ERR_NONE) && (idx < 2); ++idx)
  {
    WlzValues	val;

    val.v = WlzNewValueTb(gObj, tt, bgdV, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      cpd->o[idx] = WlzAssignObject(
                    WlzMakeMain(WLZ_2D_DOMAINOBJ, gObj->domain, val,
			        NULL, NULL, &errNum), NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzIntervalWSpace iWSp[2];
    WlzGreyWSpace     gWSp[2];

    if(((errNum = WlzInitGreyScan(
                      cpd->o[0], &(iWSp[0]), &(gWSp[0]))) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(
                      cpd->o[1], &(iWSp[1]), &(gWSp[1]))) == WLZ_ERR_NONE))
    {
      while(((errNum = WlzNextGreyInterval(&(iWSp[0]))) == WLZ_ERR_NONE) &&
            ((errNum = WlzNextGreyInterval(&(iWSp[1]))) == WLZ_ERR_NONE))
      {
	int	kl,
		ln;
	WlzGreyP pX,
		 pY;
        
	ln = iWSp[0].linpos;
	pX = gWSp[0].u_grintptr;
	pY = gWSp[1].u_grintptr;
	if(radial)
	{
	  double y,
		 ysq;

	  y = ln;
	  ysq = y * y;
	  for(kl = iWSp[0].lftpos; kl <= iWSp[0].rgtpos; ++kl)
	  {
	    double	a,
	    		r,
			x;

	    x = kl;
	    r = sqrt((x * x) + ysq);
	    a = atan2(y, x);
	    if(gType == WLZ_GREY_INT)
	    {
	      *(pX.inp)++ = WLZ_NINT(r);
	      *(pY.inp)++ = WLZ_NINT(a);
	    }
	    else /* WLZ_GREY_DOUBLE */
	    {
	      *(pX.dbp)++ = r;
	      *(pY.dbp)++ = a;
	    }
	  }
	}
	else
	{
	  if(gType == WLZ_GREY_INT)
	  {
	    for(kl = iWSp[0].lftpos; kl <= iWSp[0].rgtpos; ++kl)
	    {
	      *(pX.inp)++ = kl;
	      *(pY.inp)++ = ln;
	    }
	  }
	  else /* WLZ_GREY_DOUBLE */
	  {
	    *(pX.dbp)++ = kl;
	    *(pY.dbp)++ = ln;
	  }
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
        errNum = WLZ_ERR_NONE;
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj((WlzObject *)cpd);
    cpd = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cpd);
}

/*!
* \return	New compound object with an object per dimension or
* 		NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Given a 3D spatial domain object this function creates
* 		a compound array object with 3 component objects with
* 		the object domains those of the given object and the
* 		object values the coordinate values.
* \param	gObj			Given spatial domain object.
* \param	radial			Use cylindrical coordinates if
* 					non-zero.
* \param	gType			Required grey type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray		*WlzValuesFromCoords3D(
				  WlzObject *gObj,
				  int radial,
				  WlzGreyType gType,
				  WlzErrorNum *dstErr)
{
  WlzDomain 	dom;
  WlzPixelV	bgdV;
  WlzCompoundArray *cpd;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  dom = gObj->domain;
  if(gType == WLZ_GREY_INT)
  {
    bgdV.v.inv = 0;
    bgdV.type = WLZ_GREY_INT;
  }
  else /* gType == WLZ_GREY_DOUBLE */
  {
    bgdV.v.dbv = 0;
    bgdV.type = WLZ_GREY_DOUBLE;
  }
  cpd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1,
                             1, 3, NULL, gObj->type, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx;

    for(idx = 0; (errNum == WLZ_ERR_NONE) && (idx < 3); ++idx)
    {
      WlzValues	val;

      val.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
                                    dom.p->plane1, dom.p->lastpl,
				    bgdV, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	cpd->o[idx] = WlzAssignObject(
		       WlzMakeMain(gObj->type, dom, val, NULL, NULL,
			           &errNum), NULL);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		pln;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(pln = dom.p->plane1; pln <= dom.p->lastpl; ++pln)
    {
      int	idx;

      if(errNum == WLZ_ERR_NONE)
      {
	WlzValues val;
	WlzObject *obj2;
	WlzCompoundArray *cpd2 = NULL;
        WlzErrorNum errNum2 = WLZ_ERR_NONE;

	val.core = NULL;
	idx = pln - dom.p->plane1;
	obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			   gObj->domain.p->domains[idx], val,
			   NULL, NULL, &errNum2);
	if(errNum2 == WLZ_ERR_NONE)
	{
	  cpd2 = WlzValuesFromCoords2D(obj2, radial, gType, &errNum2);
	}
	if(errNum2 == WLZ_ERR_NONE)
	{
	  cpd->o[0]->values.vox->values[idx] = WlzAssignValues(
				    cpd2->o[0]->values, NULL);
	  cpd->o[1]->values.vox->values[idx] = WlzAssignValues(
				    cpd2->o[1]->values, NULL);
	}
	(void )WlzFreeObj((WlzObject *)cpd2);
	if(errNum2 == WLZ_ERR_NONE)
	{
	  WlzObjectType tt;

	  tt = WlzGreyTableType(WLZ_GREY_TAB_RECT, gType, NULL);
	  val.v = WlzNewValueTb(obj2, tt, bgdV, &errNum2);
	}
	if(errNum2 == WLZ_ERR_NONE)
	{
	  WlzPixelV plnV;

	  if(gType == WLZ_GREY_INT)
	  {
	    plnV.type = WLZ_GREY_INT;
	    plnV.v.inv = pln;
	  }
	  else
	  {
	    plnV.type = WLZ_GREY_DOUBLE;
	    plnV.v.dbv = pln;
	  }
	  obj2->values = WlzAssignValues(val, NULL);
	  errNum2 = WlzGreySetValue(obj2, plnV);
	}
	if(errNum2 == WLZ_ERR_NONE)
	{
	  cpd->o[2]->values.vox->values[idx] = WlzAssignValues(
				  obj2->values, NULL);
	}
	(void )WlzFreeObj(obj2);
#ifdef _OPENMP
#pragma omp critical
	{
#endif
	  if((errNum == WLZ_ERR_NONE) && (errNum2 != WLZ_ERR_NONE))
	  {
	    errNum = errNum2;
	  }
#ifdef _OPENMP
	}
#endif
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj((WlzObject *)cpd);
    cpd = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cpd);
}

