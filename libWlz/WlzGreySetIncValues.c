#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreySetIncValues_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzGreySetIncValues.c
* \author       Bill Hill
* \date         November 2009
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2009 Medical research Council, UK.
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

static WlzObject 		*WlzGreySetIncValues2D(
				  WlzObject *in,
				  int *val,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzGreySetIncValues3D(
				  WlzObject *in,
				  int *val,
				  WlzErrorNum *dstErr);

/*!
* \return	New grey value domain object with incrementing values.
* \ingroup	WlzValuesUtils
* \brief	Creates a new domain object with integer values that increment
* 		throughout the object in scan order. Object values start from
* 		1.
* \param	in			Input domain object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzGreySetIncValues(WlzObject *in, WlzErrorNum *dstErr)
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
	out = WlzGreySetIncValues2D(in, &val, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	out = WlzGreySetIncValues3D(in, &val, &errNum);
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
static WlzObject *WlzGreySetIncValues2D(WlzObject *in, int *val,
					WlzErrorNum *dstErr)
{
  int		idV,
  		iWidth;
  int		*valP;
  WlzObject     *out = NULL;
  WlzObjectType	gTT;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace gWsp;
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
    errNum = WlzInitGreyScan(out, &iWSp, &gWsp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
    {
      valP = gWsp.u_grintptr.inp;
      iWidth = iWSp.rgtpos - iWSp.lftpos + 1;
      for(idV = 0; idV < iWidth; ++idV)
      {
        *valP++ = (*val)++;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
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
static WlzObject *WlzGreySetIncValues3D(WlzObject *in, int *val,
					WlzErrorNum *dstErr)
{
  int		idO,
  		idP;
  WlzObject     *in2D,
  		*out2D,
		*out = NULL;
  WlzPixelV	bgd;
  WlzValues	values,
  		nullValues;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  bgd.type = WLZ_GREY_INT;
  bgd.v.inv = 0;
  values.core = nullValues.core =NULL;
  values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
                                   in->domain.p->plane1, in->domain.p->lastpl,
				   bgd, NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    out = WlzMakeMain(in->type, in->domain, values, in->plist, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idP = out->domain.p->plane1; idP <= out->domain.p->lastpl; ++idP)
    {
      in2D = out2D = NULL;
      idO = idP - out->domain.p->plane1;
      in2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *(out->domain.p->domains + idO),
			 nullValues, NULL, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        out2D = WlzGreySetIncValues2D(in2D, val, &errNum);
	(void )WlzFreeObj(in2D);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        *(values.vox->values + idO) = WlzAssignValues(out2D->values, NULL);
	WlzFreeObj(out2D);
      }
      else
      {
        break;
      }
    }
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
