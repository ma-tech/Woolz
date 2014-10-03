#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreySetIncValues_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreySetIncValues.c
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
* \brief	Functions for creating objects with integer grey values
* 		that increment.
* \ingroup	WlzValuesUtils
*/
#include <Wlz.h>

static WlzObject 		*WlzGreyNewIncValues2D(
				  WlzObject *in,
				  int *val,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzGreyNewIncValues3D(
				  WlzObject *in,
				  int *val,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzGreySetIncValuesItr(
				  WlzObject *obj,
				  int *val);

/*!
* \return	New grey value domain object with incrementing values.
* \ingroup	WlzValuesUtils
* \brief	Creates a new domain object with integer values that increment
* 		throughout the object in scan order. Object values start from
* 		1.
* \param	in			Input domain object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzGreyNewIncValues(WlzObject *in, WlzErrorNum *dstErr)
{
  int		val = 1;
  WlzObject	*out = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(in == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(in->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(in->type)
    {
      case WLZ_2D_DOMAINOBJ:
	out = WlzGreyNewIncValues2D(in, &val, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	out = WlzGreyNewIncValues3D(in, &val, &errNum);
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
  return(out);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	Sets the values of a 2 or 3D domain object with int 
* 		values so that they increment throughout the object
* 		in scan order. Object values are set by incrementing
* 		the given value in place. The given object must have
* 		WLZ_GREY_INT values.
* \param	obj			The given object.
* \param	gVal			Pointer to current value, this is
* 					incremented in place. If NULL then
* 					the values will be incremented from
* 					zero.
*/
WlzErrorNum 	WlzGreySetIncValues(WlzObject *obj, int *gVal)
{
  int		val;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    gType = WlzGreyTypeFromObj(obj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(gType != WLZ_GREY_INT)
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    val = (gVal)? *gVal: 0;
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTROUGH */
      case WLZ_3D_DOMAINOBJ:
	errNum = WlzGreySetIncValuesItr(obj, &val);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    if(gVal && (errNum == WLZ_ERR_NONE))
    {
      *gVal = val;
    }
  }
  return(errNum);
}

/*!
* \return	New grey value domain object with incrementing values.
* \ingroup	WlzValuesUtils
* \brief	Creates a new 2D domain object with integer values that
* 		increment throughout the object in scan order. Object
* 		values are set by incrementing the given value in place.
* \param	in			Input domain object.
* \param	val			Pointer to current value, this is
* 					incremented in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzGreyNewIncValues2D(WlzObject *in, int *val,
					WlzErrorNum *dstErr)
{
  WlzObject     *out = NULL;
  WlzObjectType	gTT;
  WlzPixelV	bgd;
  WlzValues	values;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  bgd.type = WLZ_GREY_INT;
  bgd.v.inv = 0;
  values.core = NULL;
  gTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    values.v = WlzNewValueTb(in, gTT, bgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    out = WlzMakeMain(in->type, in->domain, values, in->plist, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetIncValuesItr(out, val);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(out != NULL)
    {
      (void )WlzFreeObj(out);
      out = NULL;
    }
    else if(values.core != NULL)
    {
      (void )WlzFreeValues(values);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(out);
}

/*!
* \return	New grey value domain object with incrementing values.
* \ingroup	WlzValuesUtils
* \brief	Creates a new 3D domain object with integer values that
* 		increment throughout the object in scan order. Object
* 		values are set by incrementing the given value in place.
* \param	in			Input domain object.
* \param	val			Pointer to current value, this is
* 					incremented in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzGreyNewIncValues3D(WlzObject *in, int *val,
					WlzErrorNum *dstErr)
{
  WlzObjectType	gTT;
  WlzObject     *out = NULL;
  WlzPixelV	bgd;
  WlzValues	values;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  bgd.type = WLZ_GREY_INT;
  bgd.v.inv = 0;
  values.core = NULL;
  gTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, &errNum);
  values.vox = WlzNewValuesVox(in, gTT, bgd, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    out = WlzMakeMain(in->type, in->domain, values, in->plist, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetIncValuesItr(out, val);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(out)
    {
      (void )WlzFreeObj(out);
      out = NULL;
    }
    else if(values.core != NULL)
    {
      (void )WlzFreeVoxelValueTb(values.vox);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(out);
}
	
/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	Sets values which increment from the given value throughout
* 		the given objct in WLZ_RASTERDIR_IPILIC scan order.
* \param	obj			Given object.
* \param	val			Pointer to value, with value that is
* 					to be set and on return next value that
* 					would be set.
*/
static WlzErrorNum WlzGreySetIncValuesItr(WlzObject *obj, int *val)
{
  WlzIterateWSpace *itWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  itWSp = WlzIterateInit(obj, WLZ_RASTERDIR_IPILIC, 1, &errNum);
  while(errNum == WLZ_ERR_NONE)
  {
    if((errNum = WlzIterate(itWSp)) == WLZ_ERR_NONE)
    {
      *(itWSp->gP.inp) = *val;
      ++*val;
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  WlzIterateWSpFree(itWSp);
  return(errNum);
}
