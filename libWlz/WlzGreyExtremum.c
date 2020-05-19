#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzWlzGreyExtremum_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzWlzGreyExtremum.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Finds the location of extremum values within an object.
* \ingroup	WlzFeatures
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static WlzIVertex3		WlzGreyExtremumPos2(
				  WlzObject *gObj,
				  int pln,
				  int isMax,
				  WlzPixelV *dstVal,
				  WlzErrorNum *dstErr);
static WlzIVertex3		WlzGreyExtremumPos3(
				  WlzObject *gObj,
				  int isMax,
				  WlzPixelV *dstVal,
				  WlzErrorNum *dstErr);

/*!
* \return	Position of extremum value.
* \ingroup	WlzFeatures
* \brief	Finds the position of an extremum value within the given
* 		object.
* \param	gObj			Given object which must have a valid
* 					spatial domain and values.
* \param	isMax			Non zero for extremum maximum else
* 					extremum is minimum.
* \param	dstVal			Destination pointer for extremum
* 					value, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIVertex3			WlzGreyExtremumPos(
				  WlzObject *gObj,
				  int isMax,
				  WlzPixelV *dstVal,
				  WlzErrorNum *dstErr)
{
  WlzIVertex3	rPos = {0};
  WlzPixelV	val = {0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        rPos = WlzGreyExtremumPos2(gObj, 0, isMax, &val, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
        rPos = WlzGreyExtremumPos3(gObj, isMax, &val, &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstVal)
  {
    *dstVal = val;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rPos);
}

/*!
* \return	Object with extremum value.
* \ingroup	WlzFeatures
* \brief	Finds the sub-domain of the given object with the extremum
* 		value.
* \param	gObj			Given object which must have a valid
* 					spatial domain and values.
* \param	isMax			Non zero for extremum maximum else
* 					extremum is minimum.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzGreyExtremumObj(
				  WlzObject *gObj,
				  int isMax,
				  WlzErrorNum *dstErr)
{
  WlzPixelV   	minV,
                maxV;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    errNum = WlzGreyRange(gObj, &minV, &maxV);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzThreshold(gObj, (isMax)? maxV: minV, WLZ_THRESH_EQUAL, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Position of extremum value.
* \ingroup	WlzFeatures
* \brief	Finds the position of an extremum value within the given
* 		2 object.
* \param	gObj			Given object which must have a valid
* 					spatial domain and values.
* \param	pln			Plane coordinate, may be zero if
*					the values are 2 and not a 3D tiled
*					object. The plane coordinate will
*					be set to this value on return.
* \param	isMax			Non zero for extremum maximum else
* 					extremum is minimum.
* \param	dstVal			Destination pointer for extremum
* 					value, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzIVertex3		WlzGreyExtremumPos2(
				  WlzObject *gObj,
				  int pln,
				  int isMax,
				  WlzPixelV *dstVal,
				  WlzErrorNum *dstErr)
{
  int		first = 1;
  WlzPixelV	val = {0};
  WlzIVertex3	pos = {0};
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace gWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((errNum = WlzInitGreyScan(gObj, &iWSp, &gWSp)) == WLZ_ERR_NONE)
  {
    if(gWSp.tvb)
    {
      iWSp.plnpos = pln;
    }
    pos.vtZ = pln;
  }
  while((WlzNextGreyInterval(&iWSp) == 0) && (errNum == WLZ_ERR_NONE))
  {
    int		idX;
    WlzGreyP    gPix;

    gPix = gWSp.u_grintptr;                
    if(first)
    {
      first = 0;
      val.type = gWSp.pixeltype;
      pos.vtX = iWSp.lftpos;
      pos.vtY = iWSp.linpos;
      switch(gWSp.pixeltype)
      {
        case WLZ_GREY_INT:                 
	  val.v.inv = *(gPix.inp);
	  break;
        case WLZ_GREY_SHORT:
	  val.v.shv = *(gPix.shp);
	  break;
	case WLZ_GREY_UBYTE:
	  val.v.ubv = *(gPix.ubp);
	  break;
	case WLZ_GREY_FLOAT:
	  val.v.flv = *(gPix.flp);
	  break;
	case WLZ_GREY_DOUBLE:
	  val.v.dbv = *(gPix.dbp);
	  break;
        default:
          errNum = WLZ_ERR_GREY_TYPE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      for(idX = iWSp.lftpos; idX <= iWSp.rgtpos; ++idX)
      {
	switch(gWSp.pixeltype)
	{
	  case WLZ_GREY_INT:                 
	    if((isMax  && (*(gPix.inp) > val.v.inv)) ||
	       (!isMax && (*(gPix.inp) < val.v.inv)))
	    {
	      pos.vtX = idX;
	      pos.vtY = iWSp.linpos;
	      val.v.inv = *(gPix.inp);
	    }
	    ++(gPix.inp);
	    break;
	  case WLZ_GREY_SHORT:
	    if((isMax  && (*(gPix.shp) > val.v.shv)) ||
	       (!isMax && (*(gPix.shp) < val.v.shv)))
	    {
	      pos.vtX = idX;
	      pos.vtY = iWSp.linpos;
	      val.v.shv = *(gPix.shp);
	    }
	    ++(gPix.shp);
	    break;
	  case WLZ_GREY_UBYTE:
	    if((isMax  && (*(gPix.ubp) > val.v.ubv)) ||
	       (!isMax && (*(gPix.ubp) < val.v.ubv)))
	    {
	      pos.vtX = idX;
	      pos.vtY = iWSp.linpos;
	      val.v.ubv = *(gPix.ubp);
	    }
	    ++(gPix.ubp);
	    break;
	  case WLZ_GREY_FLOAT:
	    if((isMax  && (*(gPix.flp) > val.v.flv)) ||
	       (!isMax && (*(gPix.flp) < val.v.flv)))
	    {
	      pos.vtX = idX;
	      pos.vtY = iWSp.linpos;
	      val.v.flv = *(gPix.flp);
	    }
	    ++(gPix.flp);
	    break;
	  case WLZ_GREY_DOUBLE:
	    if((isMax  && (*(gPix.dbp) > val.v.dbv)) ||
	       (!isMax && (*(gPix.dbp) < val.v.dbv)))
	    {
	      pos.vtX = idX;
	      pos.vtY = iWSp.linpos;
	      val.v.dbv = *(gPix.dbp);
	    }
	    ++(gPix.dbp);
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	}
      }
    }
  }
  (void )WlzEndGreyScan(&iWSp, &gWSp);
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  if(dstVal && (errNum == WLZ_ERR_NONE))
  {
    *dstVal = val;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pos);
}

/*!
* \return	Position of extremum value.
* \ingroup	WlzFeatures
* \brief	Finds the position of an extremum value within the given
* 		3D object.
* \param	gObj			Given object which must have a valid
* 					spatial domain and values.
* \param	isMax			Non zero for extremum maximum else
* 					extremum is minimum.
* \param	dstVal			Destination pointer for extremum
* 					value, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzIVertex3		WlzGreyExtremumPos3(
				  WlzObject *gObj,
				  int isMax,
				  WlzPixelV *dstVal,
				  WlzErrorNum *dstErr)
{
  int		tiled = 0;
  WlzPixelV	val = {0};
  WlzIVertex3	pos = {0};
  WlzDomain     *domains;
  WlzValues     *values = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((domains = gObj->domain.p->domains) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    if(WlzGreyTableIsTiled(gObj->values.core->type))
    {
      tiled = 1;
    }
    else
    {
      if(((values = gObj->values.vox->values) == NULL))
      {
        errNum = WLZ_ERR_VALUES_NULL;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		pCnt,
    		pIdx,
		first = 1;

    pCnt = gObj->domain.p->lastpl - gObj->domain.p->plane1 + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(pIdx = 0; pIdx < pCnt; ++pIdx)
    {
      int	pln;
      WlzIVertex3 pos2 = {0};
      WlzPixelV   val2 = {0};
      WlzObject *obj2 = NULL;
      WlzErrorNum errNum2 = WLZ_ERR_NONE;

      pln = gObj->domain.p->plane1 + pIdx;
      if(tiled)
      {
	obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[pIdx],
			    gObj->values, NULL, NULL, &errNum2);
      }
      else
      {
        obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[pIdx], values[pIdx],
	                    NULL, NULL, &errNum2);
      }
      if(errNum2 == WLZ_ERR_NONE)
      {
        pos2 = WlzGreyExtremumPos2(obj2, pln, isMax, &val2, &errNum2);
      }
      WlzFreeObj(obj2);
#ifdef _OPENMP
#pragma omp critical (WlzGreyExtremum)
      {
#endif
        if(errNum2 == WLZ_ERR_NONE)
	{
	  if(first)
	  {
	    first = 0;
	    pos = pos2;
	    val = val2;
	  }
	  else
	  {
	    switch(val2.type)
	    {
	      case WLZ_GREY_INT:
		if((isMax  && (val2.v.inv > val.v.inv)) ||
		   (!isMax && (val2.v.inv < val.v.inv)))
		{
		  pos = pos2;
		  val = val2;
		}
		break;
	      case WLZ_GREY_SHORT:
		if((isMax  && (val2.v.shv > val.v.shv)) ||
		   (!isMax && (val2.v.shv < val.v.shv)))
		{
		  pos = pos2;
		  val = val2;
		}
		break;
	      case WLZ_GREY_UBYTE:
		if((isMax  && (val2.v.ubv > val.v.ubv)) ||
		   (!isMax && (val2.v.ubv < val.v.ubv)))
		{
		  pos = pos2;
		  val = val2;
		}
		break;
	      case WLZ_GREY_FLOAT:
		if((isMax  && (val2.v.flv > val.v.flv)) ||
		   (!isMax && (val2.v.flv < val.v.flv)))
		{
		  pos = pos2;
		  val = val2;
		}
		break;
	      case WLZ_GREY_DOUBLE:
		if((isMax  && (val2.v.dbv > val.v.dbv)) ||
		   (!isMax && (val2.v.dbv < val.v.dbv)))
		{
		  pos = pos2;
		  val = val2;
		}
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
	    }
	  }
	}
	if(errNum2 != WLZ_ERR_NONE)
	{
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = errNum2;
	  }
	}
#ifdef _OPENMP
      }
#endif

    }
  }
  if(dstVal && (errNum == WLZ_ERR_NONE))
  {
    *dstVal = val;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pos);
}

