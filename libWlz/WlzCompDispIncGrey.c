#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCompDispIncGrey_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCompDispIncGrey.c
* \author       Bill Hill
* \date         November 2009
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
* \brief	Functions to compute vector displacement maps from two
* 		warps of domain objects which originate from a common
* 		object with incrementing integer grey values.
* \ingroup	WlzValueUtils
*/
#include <limits.h>
#include <Wlz.h>

static WlzCompoundArray 	*WlzCompDispIncGrey2D(
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzErrorNum *dstErr);
static WlzCompoundArray 	*WlzCompDispIncGrey3D(
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzErrorNum *dstErr);
static int			*WlzCompDispMakeValAry3D(
				  WlzObject *obj,
				  WlzLong *dstNAry,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzCompDispSetAry(
				  int **ary,
				  WlzObject *obj,
				  int pln,
				  int dim);
static int			WlzCompDispArySortFn3D(
				  const void *p0,
				  const void *p1);
static WlzLong			WlzCompDispFindDsp(
				  WlzLong arraySz,
				  int *array,
				  int vQ,
				  int pad);


/*!
* \return	Compound array with integer displacement values in
* 		x, y, z order.
* \ingroup	WlzValueUtils
* \brief	Computes a compound array object which will have either
* 		two or three members that are the elements of the
* 		computed displacement vector field. This contains the
* 		displacements from the first to the second object.
* 		The two objects must have been derived from the same
* 		domain object with integer incrementing grey values,
* 		see WlzGreySetIncValues().
* \param	obj0			First object.
* \param	obj1			Second object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCompoundArray  *WlzCompDispIncGrey(WlzObject *obj0, WlzObject *obj1,
				      WlzErrorNum *dstErr)
{
  WlzCompoundArray *dsp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((obj0 == NULL) || (obj1 == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj0->type != obj1->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((obj0->domain.core == NULL) || (obj1->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj0->domain.core->type != obj1->domain.core->type)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((obj0->values.core == NULL) || (obj1->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((obj0->values.core->type != obj1->values.core->type) ||
          WlzGreyTableIsTiled(obj0->values.core->type))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    switch(obj0->type)
    {
      case WLZ_2D_DOMAINOBJ:
	dsp = WlzCompDispIncGrey2D(obj0, obj1, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	dsp = WlzCompDispIncGrey3D(obj0, obj1, &errNum);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dsp);
}

/*!
* \return	Compound array with vector displacements.
* \ingroup	WlzValueUtils
* \brief	Computes the 2D displacement field from the first to the
* 		second object, where both are 2D domain objects that are
* 		derived from the same 2D domain object with integer
* 		incremental grey values.
* \param	obj0			First object.
* \param	obj1			Second object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray *WlzCompDispIncGrey2D(WlzObject *obj0, WlzObject *obj1,
				              WlzErrorNum *dstErr)
{
  WlzCompoundArray *dsp = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  errNum = WLZ_ERR_UNIMPLEMENTED;
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dsp);
}

/*!
* \return	Compound array with vector displacements.
* \ingroup	WlzValueUtils
* \brief	Computes the 3D displacement field from the first to the
* 		second object, where both are 3D domain objects that are
* 		derived from the same 3D domain object with integer
* 		incremental grey values.
* \param	obj0			First object.
* \param	obj1			Second object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray *WlzCompDispIncGrey3D(WlzObject *obj0, WlzObject *obj1,
				              WlzErrorNum *dstErr)
{
  int		idN,
  		idO,
		idP,
		iWidth = 0;
  WlzLong 	nAry = 0;
  int		*ary = NULL,
		*valF = NULL;
  int		*val[4];
  WlzObject	*objs[3],
  		*objs2D[4];
  WlzObjectType	gTT;
  WlzPlaneDomain *domP;
  WlzCompoundArray *dsp = NULL;
  WlzIntervalWSpace iWSp[4];
  WlzGreyWSpace gWSp[4];
  WlzPixelV     bgd;
  WlzValues     values;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  /* Create compound object for displacements from the first object. */
  objs[0] = objs[1] = objs[2] = NULL;
  bgd.type = WLZ_GREY_INT;
  bgd.v.inv = 0;
  gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, &errNum);
  for(idN = 0; idN < 3; ++idN)
  {
    values.vox = WlzNewValuesVox(obj0, gTT, bgd, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      objs[idN] = WlzMakeMain(obj0->type, obj0->domain, values, NULL, NULL,
      			      &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
  }
  /* Create sorted (by grey value) array of values and positions for the
   * second object. */
  if(errNum == WLZ_ERR_NONE)
  {
    ary = WlzCompDispMakeValAry3D(obj1, &nAry, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    qsort(ary, nAry, 4 * sizeof(int), WlzCompDispArySortFn3D);
  }
  /* Scan through the first object computing displacements. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzLong	idF = 0;

    domP = obj0->domain.p;
    for(idP = domP->plane1; idP <= domP->lastpl; ++idP)
    {
      idO = idP - domP->plane1;
      objs2D[0] = objs2D[1] = objs2D[2];
      objs2D[3] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
                            *(domP->domains + idO),
                            *(obj0->values.vox->values + idO),
			    NULL, NULL, &errNum);
      for(idN = 0; (errNum == WLZ_ERR_NONE) && (idN < 3); ++idN)
      {
	objs2D[idN] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			     *(objs[idN]->domain.p->domains + idO),
			     *(objs[idN]->values.vox->values + idO),
			     NULL, NULL, &errNum);
      }
      for(idN = 0; (errNum == WLZ_ERR_NONE) && (idN < 4); ++idN)
      {
        errNum = WlzInitGreyScan(objs2D[idN], iWSp + idN, gWSp + idN);
      }
      while(((errNum = WlzNextGreyInterval(iWSp + 0)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(iWSp + 1)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(iWSp + 2)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(iWSp + 3)) == WLZ_ERR_NONE))
      {
	iWidth = iWSp[0].rgtpos - iWSp[0].lftpos + 1;
	val[0] = gWSp[0].u_grintptr.inp;
	val[1] = gWSp[1].u_grintptr.inp;
	val[2] = gWSp[2].u_grintptr.inp;
	val[3] = gWSp[3].u_grintptr.inp;
	for(idN = 0; idN < iWidth; ++idN)
	{
	  /* Find value in array. */
	  idF = WlzCompDispFindDsp(nAry, ary, *(val[3]), 4);
	  if(idF >= 0)
	  {
	    valF = ary + (idF * 4);
	    *(val[0]) = valF[1] - (iWSp[0].lftpos + idN); /* x */
	    *(val[1]) = valF[2] - iWSp[0].linpos;         /* y */
	    *(val[2]) = valF[3] - idP;                    /* z */
	  }
	  else
	  {
	    *(val[0]) = INT_MAX;
	    *(val[1]) = INT_MAX;
	    *(val[2]) = INT_MAX;
	  }
	  ++(val[0]);
	  ++(val[1]);
	  ++(val[2]);
	  ++(val[3]);
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
        errNum = WLZ_ERR_NONE;
      }
      (void )WlzFreeObj(objs2D[0]);
      (void )WlzFreeObj(objs2D[1]);
      (void )WlzFreeObj(objs2D[2]);
      (void )WlzFreeObj(objs2D[3]);
    }
  }
  AlcFree(ary);
  if(errNum == WLZ_ERR_NONE)
  {
    dsp = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 3, objs,
		               WLZ_3D_DOMAINOBJ, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dsp != NULL)
    {
      (void )WlzFreeObj((WlzObject *)dsp);
    }
    else
    {
      for(idN = 0; idN < 3; ++idN)
      {
        (void )WlzFreeObj(objs[idN]);
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dsp);
}

/*!
* \return	Integer array with values and coordinates from 3D object.
* \ingroup	WlzValueUtils
* \brief	Allocates a new array (4 ints per value: 0 = value,
* 		1 = x coordinate, 2 = y coordinate and 3 = z coordinate.
* \param	obj			Given object which must be a valid
* 					3D domain object with integer values.
* \param	dstNAry			Destination pointer for the number of
* 					values, must not be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static int	*WlzCompDispMakeValAry3D(WlzObject *obj, WlzLong *dstNAry,
				         WlzErrorNum *dstErr)
{
  int		idO,
  		idP;
  WlzLong	nAry;
  int		*ary = NULL,
  		*array = NULL;
  WlzObject	*obj2D;
  WlzPlaneDomain *pDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nAry = WlzVolume(obj, &errNum)) <= 0)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((array = AlcMalloc(nAry * 4 * sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ary = array;
    pDom = obj->domain.p;
    for(idP = pDom->plane1; (errNum == WLZ_ERR_NONE) && (idP <= pDom->lastpl);
        ++idP)
    {
      idO = idP - pDom->plane1;
      obj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			  *(obj->domain.p->domains + idO),
			  *(obj->values.vox->values + idO),
			  NULL, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzCompDispSetAry(&ary, obj2D, idP, 3);
	WlzFreeObj(obj2D);
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    AlcFree(ary);
    ary = NULL;
  }
  else
  {
    *dstNAry = nAry;
    if(dstErr != NULL)
    {
      *dstErr = errNum;
    }
  }
  return(array);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValueUtils
* \brief	Sets a portion of a value array using the iven 2D domain
* 		object.
* \param	ary			Current pointer within the array.
* \param	obj			Given object which must be a valid 2D
* 					domain object with integer values.
* \param	pln			The plane coordinate. Not used for
* 					2D.
* \param	dim			The dimension of the object that the
* 					array is being filled from (2 or 3).
*/
static WlzErrorNum WlzCompDispSetAry(int **ary, WlzObject *obj, int pln,
				     int dim)
{
  int		idV,
		iWidth;
  int		*valP;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace gWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzInitGreyScan(obj, &iWSp, &gWSp);
  if(errNum == WLZ_ERR_NONE)
  {
    while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
    {
      valP = gWSp.u_grintptr.inp;
      iWidth = iWSp.rgtpos - iWSp.lftpos + 1;
      for(idV = 0; idV < iWidth; ++idV)
      {
        *(*ary)++ = *valP++;
	*(*ary)++ = iWSp.colpos + idV;
	*(*ary)++ = iWSp.linpos;
	if(dim == 3)
	{
	  *(*ary)++ = pln;
	}
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  return(errNum);
}

/*!
* \return	Comparison value for qsort().
* \ingroup	WlzValueUtils
* \brief	Compares two integer for qsort().
* \param	p0			Pointer for first integer.
* \param	p1			Pointer for second integer.
*/
static int	WlzCompDispArySortFn3D(const void *p0, const void *p1)
{
  int		cmp;
  int		*i0,
  		*i1;
  i0 = (int *)p0;
  i1 = (int *)p1;
  cmp = *i0 - *i1;
  return(cmp);
}

/*!
* \return	Index in terms of the number of values.
* \ingroup	WlzValueUtils
* \brief	Finds the index of the entry in the given array which has
* 		the same value as the given value using a binary search.
* \param	arraySz			Number of values in the array.
* \param	array			The array of values and coordinates.
* \param	vQ			Value to query.
* \param	pad			Number of integers per array value.
*/
static WlzLong	WlzCompDispFindDsp(WlzLong arraySz, int *array,
				   int vQ, int pad)
{
  int		v2;
  WlzLong	l0,
  		l1,
		l2,
		lF = -1;

  l0 = 0;
  l1 = arraySz;
  while(l0 < l1)
  {
    l2 = l0 + ((l1 - l0) / 2);
    v2 = array[l2 * pad];
    if(v2 < vQ)
    {
      l0 = l2 + 1; 
    }
    else
    {
      l1 = l2; 
    }
  }
  if((l0 < arraySz) && (array[l0 * pad] == vQ))
  {
    lF = l0;
  }
  return(lF);
}
