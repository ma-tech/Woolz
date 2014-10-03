#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFilterNObjValues_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzFilterNObjValues.c
* \author       Bill Hill
* \date         September 2007
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
* \brief	Computes a domain object in which the values are the
*		filtered values of the input objects at the same positions.
* \ingroup	WlzValuesFilters
*/
#include <Wlz.h>

typedef WlzErrorNum (*WlzNObjFilterFn)(WlzGreyType, int, WlzGreyP, int);

static void	 		WlzFilterNObjFreeWBuf(
				  WlzGreyValueWSpace **wBuf,
				  int nObj);
static void			*WlzFilterNObjMakeVBuf(
				  WlzGreyType gType,
				  int nVal);
static WlzErrorNum 		WlzFilterNObjMakeWBuf(
				  WlzGreyValueWSpace ***wBuf,
				  int nObj,
				  WlzObject **objs);
static WlzErrorNum 		WlzFilterValuesMeanFn(
				  WlzGreyType gType,
				  int nVal,
				  WlzGreyP val,
				  int idx);
static WlzErrorNum 		WlzFilterValuesRankFn(
				  WlzGreyType gType,
				  int nVal,
				  WlzGreyP val,
				  int rank);
static WlzErrorNum 		WlzFilterNObjValuesVal(
				  WlzGreyP gP,
				  int idI,
				  WlzGreyType gType,
				  int nObj, WlzObject **objs,
				  WlzGreyValueWSpace **wBuf,
				  WlzNObjFilterFn filterFn,
				  int idF,
				  WlzGreyP vBuf, WlzIVertex3 pos);
static WlzObject 		*WlzFilterNObjValues3DCb(
				  WlzObject *rObj,
				  WlzGreyType gType,
				  int nObj, WlzObject **objs,
				  WlzNObjFilterFn filterFn,
				  int idF, WlzErrorNum *dstErr);
static WlzObject 		*WlzFilterNObjValues2DCb(
				  WlzObject *rObj,
				  WlzGreyType gType,
				  int nObj,
				  WlzObject **objs,
				  WlzNObjFilterFn filterFn,
				  int idF, WlzErrorNum *dstErr);

/*!
* \return	New object in which the values are the filtered values of
*		the input objects or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies the a filter function to all pixels of the
*		input objects, with the filter being applied to all pixels
*		at each coordinate through all objects.
*		The input and reference objects must all be 2D or 3D domain
*		objects of the same type and all the input objects must have
*		valid grey values of the same type.
* \param	rObj			Reference object with domain in which
*					the filtered values are computed.
* \param	nObj			Number of input objects.
* \param	objs			Array on input objects.
* \param	fn			Filter function to apply, with value:
*					0 for rank, 1 for mean.
* \param	rank			If applying a rank filter the required
*					rank, otherwise ignored. Rank value
*					is 0.0 minimum, 0.5 median and
*					1.0 maximum. Intermediate ranks may be
*					given.
* \param	dstErr			Destination error pointer, may be null.
*/
WlzObject 	*WlzFilterNObjValues(WlzObject *rObj,
				    int nObj, WlzObject **objs,
				    int fn, double rank,
				    WlzErrorNum *dstErr)
{
  int		idO,
  		rankI;
  WlzObject	*fObj = NULL;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(rObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((rObj->type != WLZ_2D_DOMAINOBJ) && (rObj->type != WLZ_3D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(rObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(nObj <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idO = 0;
    while((errNum == WLZ_ERR_NONE) && (idO < nObj))
    {
      if(objs[idO] == NULL)
      {
	errNum = WLZ_ERR_OBJECT_NULL;
      }
      else if(objs[idO]->type != rObj->type)
      {
	errNum = WLZ_ERR_OBJECT_TYPE;
      }
      else if(objs[idO]->domain.core == NULL)
      {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if(objs[idO]->values.core == NULL)
      {
        errNum = WLZ_ERR_VALUES_NULL;
      }
      else if(WlzGreyTableIsTiled(objs[idO]->values.core->type))
      {
        errNum = WLZ_ERR_VALUES_TYPE;
      }
      else
      {
	if(idO == 0)
	{
          gType = WlzGreyTypeFromObj(objs[idO], &errNum);
	}
	else
	{
	  if(gType != WlzGreyTypeFromObj(objs[idO], NULL))
	  {
	    errNum = WLZ_ERR_GREY_TYPE;
	  }
	}
      }
      ++idO;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rank = WLZ_CLAMP(rank, 0.0, 1.0);
    rankI = WLZ_NINT(rank * (nObj - 1));
    if(rObj->type == WLZ_2D_DOMAINOBJ)
    {
      switch(fn)
      {
        case 0: /* Rank */
          fObj = WlzFilterNObjValues2DCb(rObj, gType, nObj, objs,
	                                 WlzFilterValuesRankFn, rankI,
					 &errNum);
	  break;
	case 1: /* Mean */
          fObj = WlzFilterNObjValues2DCb(rObj, gType, nObj, objs,
	                                 WlzFilterValuesMeanFn, 0,
					 &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
      }
    }
    else /* rObj->type == WLZ_3D_DOMAINOBJ */
    {
      switch(fn)
      {
        case 0: /* Rank */
          fObj = WlzFilterNObjValues3DCb(rObj, gType, nObj, objs,
	                                 WlzFilterValuesRankFn, rankI,
					 &errNum);
	  break;
	case 1: /* Mean */
          fObj = WlzFilterNObjValues3DCb(rObj, gType, nObj, objs,
	                                 WlzFilterValuesMeanFn, 0,
					 &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(fObj);
}

/*!
* \return	New object in which the values are the filtered values of
*		the input objects or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies the given filter function to all pixels of the
*		input objects, with the filter being applied to all pixels
*		at each coordinate through all objects.
*		The input and reference objects must all be 2D domain objects
*		and all the input objects must have valid values.
* 		Because this is a static function all parameters are assumed
*		to be valid.
* \param	rObj			Reference object with domain in which
*					the filtered values are computed.
* \param	gType			Grey type for all the objects.
* \param	nObj			Number of input objects.
* \param	objs			Array on input objects.
* \param	filterFn		Function to apply.
* \param	idF			Index parameter passed to function.
* \param	dstErr			Destination error pointer, may be null.
*/
static WlzObject *WlzFilterNObjValues2DCb(WlzObject *rObj, WlzGreyType gType,
					  int nObj, WlzObject **objs,
				          WlzNObjFilterFn filterFn,
					  int idF, WlzErrorNum *dstErr)
{
  int		idI;
  WlzObject	*fObj = NULL;
  WlzObjectType	gTblType;
  WlzIVertex3	pos;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace gWSp;
  WlzGreyP	vBuf;
  WlzValues	newVal;
  WlzPixelV	bgdV;
  WlzGreyValueWSpace **wBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  pos.vtZ = 0;
  bgdV.v.inv = 0;
  bgdV.type = WLZ_GREY_INT;
  gTblType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, gType, NULL);
  if((vBuf.v = WlzFilterNObjMakeVBuf(gType, nObj)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzFilterNObjMakeWBuf(&wBuf, nObj, objs);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newVal.v = WlzNewValueTb(rObj, gTblType, bgdV, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    fObj = WlzMakeMain(rObj->type, rObj->domain, newVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetValue(fObj, bgdV);
  }
  if((errNum == WLZ_ERR_NONE) &&
     ((errNum = WlzInitGreyScan(fObj, &iWSp, &gWSp)) == WLZ_ERR_NONE))
  {
    while((WlzNextGreyInterval(&iWSp) == 0) && (errNum == WLZ_ERR_NONE))
    {
      idI = 0;
      pos.vtY = iWSp.linpos;
      pos.vtX = iWSp.lftpos;
      while((errNum == WLZ_ERR_NONE) && (pos.vtX <= iWSp.rgtpos))
      {
	errNum = WlzFilterNObjValuesVal(gWSp.u_grintptr, idI,
	                                gType, nObj, objs, wBuf,
					filterFn, idF,
	                                vBuf, pos);
	++idI;
        ++pos.vtX;
      }
    }
    (void )WlzEndGreyScan(&iWSp, &gWSp);
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  AlcFree(vBuf.v);
  WlzFilterNObjFreeWBuf(wBuf, nObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(fObj);
}

/*!
* \return	New object in which the values are the filtered values of
*		the input objects or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies the given filter function to all pixels of the
*		input objects, with the filter being applied to all pixels
*		at each coordinate through all objects.
*		The input and reference objects must all be 3D domain objects
*		and all the input objects must have valid values.
* 		Because this is a static function all parameters are assumed
*		to be valid.
* \param	rObj			Reference object with domain in which
*					the filtered values are computed.
* \param	gType			Grey type for all the objects.
* \param	nObj			Number of input objects.
* \param	objs			Array on input objects.
* \param	filterFn		Function to apply.
* \param	idF			Index parameter passed to function.
* \param	dstErr			Destination error pointer, may be null.
*/
static WlzObject *WlzFilterNObjValues3DCb(WlzObject *rObj, WlzGreyType gType,
					  int nObj, WlzObject **objs,
				          WlzNObjFilterFn filterFn,
					  int idF, WlzErrorNum *dstErr)
{
  int		idI,
  		idP,
		idQ;
  WlzObject	*obj2D,
  		*fObj = NULL;
  WlzObjectType gTblType;
  WlzIVertex3	pos;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace gWSp;
  WlzGreyP	vBuf;
  WlzDomain	refDom;
  WlzValues	newVal;
  WlzPixelV	bgdV;
  WlzDomain	*domP;
  WlzValues	*valP;
  WlzGreyValueWSpace **wBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.v.inv = 0;
  bgdV.type = WLZ_GREY_INT;
  gTblType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, gType, NULL);
  if((vBuf.v = WlzFilterNObjMakeVBuf(gType, nObj)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzFilterNObjMakeWBuf(&wBuf, nObj, objs);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newVal.vox = WlzNewValuesVox(rObj, gTblType, bgdV, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    fObj = WlzMakeMain(rObj->type, rObj->domain, newVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetValue(fObj, bgdV);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    refDom = fObj->domain;
    idP = refDom.p->plane1;
    domP = refDom.p->domains;
    valP = newVal.vox->values;
    while((errNum == WLZ_ERR_NONE) && (idP <= refDom.p->lastpl))
    {
      idQ = idP - refDom.p->plane1;
      if((*(domP + idQ)).core && (*(valP + idQ)).core)
      {
	pos.vtZ = idP;
        obj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *(domP + idQ), *(valP + idQ),
			    NULL, NULL, &errNum);
	if((errNum == WLZ_ERR_NONE) &&
	   ((errNum = WlzInitGreyScan(obj2D, &iWSp, &gWSp)) == WLZ_ERR_NONE))
	{
	  while((WlzNextGreyInterval(&iWSp) == 0) && (errNum == WLZ_ERR_NONE))
	  {
	    idI = 0;
	    pos.vtY = iWSp.linpos;
	    pos.vtX = iWSp.lftpos;
	    while((errNum == WLZ_ERR_NONE) && (pos.vtX <= iWSp.rgtpos))
	    {
	      errNum = WlzFilterNObjValuesVal(gWSp.u_grintptr, idI,
					      gType, nObj, objs, wBuf,
					      filterFn, idF,
					      vBuf, pos);
	      ++idI;
	      ++pos.vtX;
	    }
	  }
          (void )WlzEndGreyScan(&iWSp, &gWSp);
	  if(errNum == WLZ_ERR_EOO)
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
	(void )WlzFreeObj(obj2D);
      }
      ++idP;
    }
  }
  AlcFree(vBuf.v);
  WlzFilterNObjFreeWBuf(wBuf, nObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(fObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Computes the filtered values for a given coordinate
*		within the given objects.
* \param	gP			Buffer for the filtered value.
* \param	idI			Offset into the buffer for the 
*					filtered value.
* \param	gType			Grey type for all objects and buffer.
* \param	nObj			Number of objects given.
* \param	objs			Given objects.
* \param	wBuf			Buffer of grey value workspaces, one
*					perr object.
* \param	filterFn		Filter function to apply.
* \param	idF			Index parameter passed to function.
* \param	vBuf			Buffer for object values.
* \param	pos			Given coordinate.
*/
static WlzErrorNum WlzFilterNObjValuesVal(WlzGreyP gP, int idI,
				          WlzGreyType gType,
					  int nObj, WlzObject **objs,
					  WlzGreyValueWSpace **wBuf,
				          WlzNObjFilterFn filterFn,
					  int idF,
					  WlzGreyP vBuf, WlzIVertex3 pos)
{
  int		idC,
  		idO;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idC = 0;
  for(idO = 0; idO < nObj; ++idO)
  {
    if(WlzInsideDomain(*(objs + idO), pos.vtZ, pos.vtY, pos.vtX, NULL))
    {
      WlzGreyValueGet(*(wBuf + idO), pos.vtZ, pos.vtY, pos.vtX);
      switch(gType)
      {
	case WLZ_GREY_INT:
	  *(vBuf.inp + idC) = (*(wBuf + idO))->gVal[0].inv;
	  break;
	case WLZ_GREY_SHORT:
	  *(vBuf.shp + idC) = (*(wBuf + idO))->gVal[0].shv;
	  break;
	case WLZ_GREY_UBYTE:
	  *(vBuf.ubp + idC) = (*(wBuf + idO))->gVal[0].ubv;
	  break;
	case WLZ_GREY_FLOAT:
	  *(vBuf.flp + idC) = (*(wBuf + idO))->gVal[0].flv;
	  break;
	case WLZ_GREY_DOUBLE:
	  *(vBuf.dbp + idC) = (*(wBuf + idO))->gVal[0].dbv;
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      ++idC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = (*(filterFn))(gType, idC, vBuf, idF);
    switch(gType)
    {
      case WLZ_GREY_INT:
	*(gP.inp + idI) = *(vBuf.inp + idF);
	break;
      case WLZ_GREY_SHORT:
	*(gP.shp + idI) = *(vBuf.shp + idF);
	break;
      case WLZ_GREY_UBYTE:
	*(gP.ubp + idI) = *(vBuf.ubp + idF);
	break;
      case WLZ_GREY_FLOAT:
	*(gP.flp + idI) = *(vBuf.flp + idF);
	break;
      case WLZ_GREY_DOUBLE:
	*(gP.dbp + idI) = *(vBuf.dbp + idF);
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Ranks the given array of values such that the n'th value
*		has rank n. This is equivalent to sorting but only for the
*		n'th value.
* \param	gType			Grey type.
* \param	nVal			Number of values in the array.
* \param	val			Array of values.
* \param	rank			Required rank, and position in array
* 					for result.
*/
static WlzErrorNum WlzFilterValuesRankFn(WlzGreyType gType,
				         int nVal, WlzGreyP val, int rank)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_INT:
      AlgRankSelectI(val.inp, nVal, rank);
      break;
    case WLZ_GREY_SHORT:
      AlgRankSelectS(val.shp, nVal, rank);
      break;
    case WLZ_GREY_UBYTE:
      AlgRankSelectUB(val.ubp, nVal, rank);
      break;
    case WLZ_GREY_FLOAT:
      AlgRankSelectF(val.flp, nVal, rank);
      break;
    case WLZ_GREY_DOUBLE:
      AlgRankSelectD(val.dbp, nVal, rank);
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Replaces the first value in the array with the mean of the
*		given array values.
* \param	gType			Grey type.
* \param	nVal			Number of values in the array,
*					must be greater than zero.
* \param	val			Array of values.
* \param	idR			Position in array for result.
*/
static WlzErrorNum WlzFilterValuesMeanFn(WlzGreyType gType,
				         int nVal, WlzGreyP val, int idR)
{
  int		idX;
  double	sum = 0.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_INT:
      for(idX = 0; idX < nVal; ++idX)
      {
        sum += val.inp[idX];
      }
      val.inp[idR] = WLZ_NINT(sum / (double )nVal);
      break;

    case WLZ_GREY_SHORT:
      for(idX = 0; idX < nVal; ++idX)
      {
        sum += val.shp[idX];
      }
      val.shp[idR] = (short )WLZ_NINT(sum / (double )nVal);
      break;
    case WLZ_GREY_UBYTE:
      for(idX = 0; idX < nVal; ++idX)
      {
        sum += val.ubp[idX];
      }
      val.ubp[idR] = (WlzUByte )WLZ_NINT(sum / (double )nVal);
      break;
    case WLZ_GREY_FLOAT:
      for(idX = 0; idX < nVal; ++idX)
      {
        sum += val.flp[idX];
      }
      val.flp[idR] = (float )(sum / (double )nVal);
      break;
    case WLZ_GREY_DOUBLE:
      for(idX = 0; idX < nVal; ++idX)
      {
        sum += val.dbp[idX];
      }
      val.dbp[idR] = sum / (double )nVal;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return	Allocated buffer.
* \ingroup	WlzValuesFilters
* \brief	Allocates a buffer of nVal values with type gType.
* \param	gType			Type of values to allocate.
* \param	nVal			Number of values to allocate.
*/
static void	*WlzFilterNObjMakeVBuf(WlzGreyType gType, int nVal)
{
  size_t	sz = 0;
  void		*gP = NULL;

  switch(gType)
  {
    case WLZ_GREY_INT:
      sz = sizeof(int);
      break;
    case WLZ_GREY_SHORT:
      sz = sizeof(short);
      break;
    case WLZ_GREY_UBYTE:
      sz = sizeof(WlzUByte);
      break;
    case WLZ_GREY_FLOAT:
      sz = sizeof(float);
      break;
    case WLZ_GREY_DOUBLE:
      sz = sizeof(double);
      break;
    default:
      break;
  }
  if(sz > 0)
  {
    gP = AlcMalloc(sz * nVal);
  }
  return(gP);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesFilters
* \brief	Allocates a buffer of grey value workspaces, one for each
*		of the given objects.
* \param	wBuf			Destination pointer for grey value
*					workspace buffer.
* \param	nObj			Number of given objects.
* \param	objs			Given objects.
*/
static WlzErrorNum WlzFilterNObjMakeWBuf(WlzGreyValueWSpace ***wBuf,
				         int nObj, WlzObject **objs)
{
  int		idO;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((*wBuf = AlcCalloc(nObj, sizeof(WlzGreyValueWSpace *))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    idO = 0;
    while((errNum == WLZ_ERR_NONE) && (idO < nObj))
    {
      *(*wBuf +idO) = WlzGreyValueMakeWSp(*(objs + idO), &errNum);
      ++idO;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFilterNObjFreeWBuf(*wBuf, nObj);
    *wBuf = NULL;
  }
  return(errNum);
}
  
/*!
* \ingroup	WlzValuesFilters
* \brief	Frees a buffer of grey value workspaces as allocated by
*		WlzFilterNObjMakeWBuf().
* \param	wBuf			Grey value workspace buffer.
* \param	nObj			Number of grey value workspaces.
*/
static void	 WlzFilterNObjFreeWBuf(WlzGreyValueWSpace **wBuf, int nObj)
{
  int		idO;

  if(wBuf)
  {
    for(idO = 0; idO < nObj; ++idO)
    {
      WlzGreyValueFreeWSp(*(wBuf + idO));
    }
    AlcFree(wBuf);
  }
}
