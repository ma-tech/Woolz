#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRank_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzRank.c
* \author       Bill Hill
* \date         March 2002
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
* \brief	Rank filters for woolz objects, these are the
* 		generalization of minimum, maximum and median
* 		value filters.
* \ingroup	WlzValuesFilters
*/

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

static WlzErrorNum 		WlzRankFilterDomObj2D(
				  WlzObject *gObj,
				  int fSz,
				  double rank);
static WlzErrorNum 		WlzRankFilterDomObj3D(
				  WlzObject *gObj,
				  int fSz,
				  double rank);
static void     		WlzRankFilterValLn(WlzObject *gObj,
				  WlzGreyValueWSpace *gVWSp,
				  void **vBuf,
				  WlzUByte **iBuf,
				  void *rBuf,
				  int lnPos,
				  WlzGreyType vType,
				  WlzIVertex2 bufSz,
				  double rank);
static void			WlzRankFilterValPl(
				  WlzObject *gObj,
				  WlzGreyValueWSpace *gVWSp,
				  void ***vBuf,
				  WlzUByte ***iBuf,
				  void *rBuf,
				  int plane,
				  WlzGreyType vType,
				  WlzIVertex3 bufSz,
				  double rank);
static void			WlzRankFilterValues(
				  WlzGreyValueWSpace *gVWSp,
				  int posZ,
				  int posY,
				  int posX,
				  WlzGreyType vType,
				  void *values,
				  int nValues,
				  double rank);

/*!
* \return	Woolz error code.
* \ingroup      WlzValuesFilters
* \brief	Applies a rank filter in place to the given Woolz object.
*		Each value of the given object is replaced by the n'th
*		ranked value of the values in it's immediate neighborhood,
*		where the neighborhood is a simple axis aligned cuboid
*		with the size.
* \param	gObj			Given object.
* \param	fSz			Rank filter size.
* \param	rank			Required rank with values:
*					0.0 minimum, 0.5 median and
*					1.0 maximum.
*/
WlzErrorNum	WlzRankFilter(WlzObject *gObj, int fSz, double rank)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(fSz < 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if(rank < 0.0 + DBL_EPSILON)
    {
      rank = DBL_EPSILON;
    }
    else if(rank > 1.0 - DBL_EPSILON)
    {
      rank = 1.0 - DBL_EPSILON;
    }
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	errNum = WlzRankFilterDomObj2D(gObj, fSz, rank);
	break;
      case WLZ_3D_DOMAINOBJ:
	errNum = WlzRankFilterDomObj3D(gObj, fSz, rank);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzValuesFilters
* \brief	Applies a rank filter in place to the given Woolz object.
*		Each value of the given object is replaced by the n'th
*		ranked value of the values in it's immediate neighborhood,
*		where the neighborhood is a simple axis aligned cube
*		with the given size.
* \param	gObj			Given object.
* \param	fSz			Rank filter size.
* \param	rank			Required rank with values:
*					0.0 minimum, 0.5 median and
*					1.0 maximum.
*/
static WlzErrorNum WlzRankFilterDomObj2D(WlzObject *gObj, int fSz,
				         double rank)
{
  int		bufLft,
		bufRgt,
		inLn,
		fSz2,
		itvWidth,
		iBufWidth, 	      /* No of bytes in interval buffer line */
		bufLn = 0,
		outLn = 0;
  void		*rBuf = NULL;
  void		**vBuf = NULL;
  WlzUByte	**iBuf = NULL;
  WlzIVertex2	bufSz;
  WlzGreyType	vType;
  WlzDomain	gDom;
  WlzGreyWSpace	gWSp;
  WlzIntervalWSpace iWSp;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    fSz2 = fSz / 2;
    gDom = gObj->domain;
    bufSz.vtX = gDom.i->lastkl - gDom.i->kol1 + (2 * fSz2) + 1;
    bufSz.vtY = fSz;
    iBufWidth = (bufSz.vtX + 7) / 8;
    vType = WlzGreyTypeFromObj(gObj, &errNum);
    switch(vType)
    {
      case WLZ_GREY_UBYTE:
	if(((rBuf = AlcMalloc(fSz * fSz * sizeof(WlzUByte))) == NULL) ||
	   (AlcUnchar2Malloc((WlzUByte ***)&vBuf,
	   		     bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_SHORT:
	if(((rBuf = AlcMalloc(fSz * fSz * sizeof(short))) == NULL) ||
	   (AlcShort2Malloc((short ***)&vBuf,
	   		     bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_INT:
	if(((rBuf = AlcMalloc(fSz * fSz * sizeof(int))) == NULL) ||
	   (AlcInt2Malloc((int ***)&vBuf,
	   		     bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_FLOAT:
	if(((rBuf = AlcMalloc(fSz * fSz * sizeof(float))) == NULL) ||
	   (AlcFloat2Malloc((float ***)&vBuf,
	   		     bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_DOUBLE:
	if(((rBuf = AlcMalloc(fSz * fSz * sizeof(double))) == NULL) ||
	   (AlcDouble2Malloc((double ***)&vBuf,
	   		     bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_RGBA:
	if(((rBuf = AlcMalloc(fSz * fSz * sizeof(WlzUInt))) == NULL) ||
	   (AlcInt2Malloc((int ***)&vBuf,
	   		     bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcBit2Calloc(&iBuf, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gVWSp = WlzGreyValueMakeWSp(gObj, &errNum);
  }
  /* Work down through the object. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(gObj, &iWSp, &gWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    inLn = gDom.i->line1 - 1;
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
    {
      if(iWSp.nwlpos)
      {
	while(inLn < iWSp.linpos)
	{
	  outLn = inLn - fSz2 - 1; 	         /* -1 because previous line */
	  if(outLn >= gDom.i->line1)
	  {
	    WlzRankFilterValLn(gObj, gVWSp, vBuf, iBuf, rBuf,
			       outLn, vType, bufSz, rank);
	  }
	  ++inLn;
	  bufLn = (inLn + bufSz.vtY - gDom.i->line1) % bufSz.vtY;
	  WlzValueSetUByte(iBuf[bufLn], 0, iBufWidth);
	}
      }
      bufLft = iWSp.lftpos + fSz2 - gDom.i->kol1;
      bufRgt = iWSp.rgtpos + fSz2 - gDom.i->kol1;
      itvWidth = bufRgt - bufLft + 1;
      /* Set bits in the interval buffer */
      WlzBitLnSetItv(iBuf[bufLn], bufLft, bufRgt, bufSz.vtX);
      /* Copy values to value buffer. */
      switch(vType)
      {
        case WLZ_GREY_UBYTE:
	  (void )memcpy(*(((WlzUByte **)vBuf) + bufLn) + bufLft,
	                gWSp.u_grintptr.ubp, itvWidth);
          break;
        case WLZ_GREY_SHORT:
	  (void )memcpy(*(((short **)vBuf) + bufLn) + bufLft,
	                gWSp.u_grintptr.shp, itvWidth * sizeof(short));
          break;
        case WLZ_GREY_INT:
	  (void )memcpy(*(((int **)vBuf) + bufLn) + bufLft,
	                gWSp.u_grintptr.inp, itvWidth * sizeof(int));
          break;
        case WLZ_GREY_FLOAT:
	  (void )memcpy(*(((float **)vBuf) + bufLn) + bufLft,
	                gWSp.u_grintptr.flp, itvWidth * sizeof(float));
          break;
        case WLZ_GREY_DOUBLE:
	  (void )memcpy(*(((double **)vBuf) + bufLn) + bufLft,
	                gWSp.u_grintptr.dbp, itvWidth * sizeof(double));
          break;
        case WLZ_GREY_RGBA:
	  (void )memcpy(*(((WlzUInt **)vBuf) + bufLn) + bufLft,
	                gWSp.u_grintptr.rgbp, itvWidth * sizeof(WlzUInt));
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
    if(errNum == WLZ_ERR_NONE)
    {
      while(++outLn <= gDom.i->lastln)
      {
	WlzRankFilterValLn(gObj, gVWSp, vBuf, iBuf, rBuf, 
			   outLn, vType, bufSz, rank);
	bufLn = (++inLn + bufSz.vtY - gDom.i->line1) % bufSz.vtY;
        WlzValueSetUByte(iBuf[bufLn], 0, iBufWidth);
      }
    }
  }
  AlcFree(rBuf);
  Alc2Free((void **)iBuf);
  Alc2Free(vBuf);
  WlzGreyValueFreeWSp(gVWSp);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzValuesFilters
* \brief	Applies a rank filter in place to the given Woolz object.
*		Each value of the given object is replaced by the n'th
*		ranked value of the values in it's immediate neighborhood,
*		where the neighborhood is a simple axis aligned cuboid
*		with the given width.
* \param	gObj			Given object.
* \param	fSz			Rank filter size.
* \param	rank			Required rank with values:
*					0.0 minimum, 0.5 median and
*					1.0 maximum.
*/
static WlzErrorNum WlzRankFilterDomObj3D(WlzObject *gObj, int fSz,
				         double rank)
{
  int		bufPl,
		plIdx,
		pnCnt,
		outPl,
		fSz2;
  void		*rBuf = NULL;
  void		***vBuf = NULL;
  WlzUByte	***iBuf = NULL;
  WlzObject	*obj2D = NULL;
  WlzValues	dummyValues;
  WlzDomain	dummyDom,
  		gDom;
  WlzIVertex3	bufSz;
  WlzIVertex2	bufOrg2D,
  		bufSz2D;
  WlzGreyType	vType;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dummyDom.core = NULL;
  dummyValues.core = NULL;
  if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    fSz2 = fSz / 2;
    gDom = gObj->domain;
    bufSz2D.vtX = bufSz.vtX = gDom.p->lastkl - gDom.p->kol1 + (2 * fSz2) + 1;
    bufSz2D.vtY = bufSz.vtY = gDom.p->lastln - gDom.p->line1 + (2 * fSz2) + 1;
    bufSz.vtZ = fSz;
    vType = WlzGreyTypeFromObj(gObj, &errNum);
    switch(vType)
    {
      case WLZ_GREY_UBYTE:
	if(((rBuf = AlcMalloc(fSz * fSz * fSz * sizeof(WlzUByte))) == NULL) ||
	   (AlcUnchar3Malloc((WlzUByte ****)&vBuf,
	   		     bufSz.vtZ, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_SHORT:
	if(((rBuf = AlcMalloc(fSz * fSz * fSz * sizeof(short))) == NULL) ||
	   (AlcShort3Malloc((short ****)&vBuf,
	   		     bufSz.vtZ, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_INT:
	if(((rBuf = AlcMalloc(fSz * fSz * fSz * sizeof(int))) == NULL) ||
	   (AlcInt3Malloc((int ****)&vBuf,
	   		     bufSz.vtZ, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_FLOAT:
	if(((rBuf = AlcMalloc(fSz * fSz * fSz * sizeof(float))) == NULL) ||
	   (AlcFloat3Malloc((float ****)&vBuf,
	   		     bufSz.vtZ, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_DOUBLE:
	if(((rBuf = AlcMalloc(fSz * fSz * fSz * sizeof(double))) == NULL) ||
	   (AlcDouble3Malloc((double ****)&vBuf,
	   		     bufSz.vtZ, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_GREY_RGBA:
	if(((rBuf = AlcMalloc(fSz * fSz * fSz * sizeof(WlzUInt))) == NULL) ||
	   (AlcInt3Malloc((int ****)&vBuf,
	   		     bufSz.vtZ, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcBit3Calloc(&iBuf, bufSz.vtZ, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom, dummyValues,
    			NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gVWSp = WlzGreyValueMakeWSp(gObj, &errNum);
  }
  /* Work down through the object. */
  if(errNum == WLZ_ERR_NONE)
  {
    plIdx = 0;
    outPl = gDom.p->plane1 + plIdx - fSz2;
    pnCnt = gDom.p->lastpl - gDom.p->plane1 + 1;
    while((errNum == WLZ_ERR_NONE) && (pnCnt-- > 0))
    {
      if((errNum == WLZ_ERR_NONE) && (outPl >= gDom.p->plane1))
      {
        WlzRankFilterValPl(gObj, gVWSp, vBuf, iBuf, rBuf,
			   outPl, vType, bufSz, rank);
      }
      bufPl = (gDom.p->plane1 + plIdx) % bufSz.vtZ;
      obj2D->domain = *(gObj->domain.p->domains + plIdx);
      obj2D->values = *(gObj->values.vox->values + plIdx);
      if(errNum == WLZ_ERR_NONE)
      {
        bufOrg2D.vtX = obj2D->domain.i->kol1 - fSz2;
	bufOrg2D.vtY = obj2D->domain.i->line1 - fSz2;
	errNum = WlzToArray2D((void ***)(iBuf + bufPl), obj2D,
			      bufSz2D, bufOrg2D, 0, WLZ_GREY_BIT);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzToArray2D((void ***)(vBuf + bufPl), obj2D,
			      bufSz2D, bufOrg2D, 0, vType);
      }
      ++plIdx;
      ++outPl;
    }
    while((errNum == WLZ_ERR_NONE) && (outPl <= gDom.p->lastpl))
    {
      WlzRankFilterValPl(gObj, gVWSp, vBuf, iBuf, rBuf,
			 outPl, vType, bufSz, rank);
      ++outPl;
    }
  }
  if(obj2D)
  {
    obj2D->domain = dummyDom;
    obj2D->values = dummyValues;
    WlzFreeObj(obj2D);
  }
  AlcFree(rBuf);
  Alc3Free((void ***)iBuf);
  Alc3Free(vBuf);
  WlzGreyValueFreeWSp(gVWSp);
  return(errNum);
}

/*!
* \return	void
* \ingroup      WlzValuesFilters
* \brief	Rank filter the buffer values for a single line.
*		If possible the buffer of values to be ranked just has
*		the values no longer used replaced by new values, rather
*		than all the values in the buffer being re-entered.
* \param	gObj			The object being filtered, needed 
*					for random pixel access.
* \param	gVWSp			Grey value workspace, needed for
*					random pixel access.
* \param	vBuf			Circular buffer of lines of grey
*					values.
* \param	iBuf			Circular buffer of lines of bits
*					used to encode intervals.
* \param	rBuf			A buffer used to rank pixels.
* \param	line			The current line, with both past
*					and future lines valid.
* \param	vType			The grey value type.
* \param	bufSz			The buffer size.
* \param	rank			Required rank with values:
*					0.0 minimum, 0.5 median and
*					1.0 maximum.
*/
static void	WlzRankFilterValLn(WlzObject *gObj, WlzGreyValueWSpace *gVWSp,
				   void **vBuf, WlzUByte **iBuf, void *rBuf,
				   int line, WlzGreyType vType,
				   WlzIVertex2 bufSz, double rank)
{
  int		idK,
		idX,
		idY,
		fSz2,
		sameRCnt,
		rCnt = 0;
  WlzGreyV	gV0,
  		gV1;
  WlzGreyP	gP0,
  		gP1;
  WlzDomain	gDom;
  WlzIVertex2	bufPos0,
		bufPos1,
  		objPos;

  gP0.v = gP1.v = NULL;
  gDom = gObj->domain;
  fSz2 = bufSz.vtY / 2;
  bufPos1.vtY = line - gDom.i->line1 - fSz2 + bufSz.vtY;
  for(idK = fSz2; idK < bufSz.vtX - fSz2; ++idK)
  {
    objPos.vtX = gDom.i->kol1 + idK - fSz2;
    objPos.vtY = line;
    bufPos0.vtX = idK;
    bufPos0.vtY = (line - gDom.i->line1 + bufSz.vtY) % bufSz.vtY;
    if(WLZ_BIT_GET(iBuf[bufPos0.vtY], bufPos0.vtX) != 0)
    {
      bufPos0.vtX = idK - fSz2 - 1;
      sameRCnt = bufPos0.vtX > 0;
      if(sameRCnt)
      {
	idY = 0;
	bufPos1.vtX = bufPos0.vtX + bufSz.vtY;
	while((idY < bufSz.vtY) && sameRCnt)
	{
	  bufPos0.vtY = (bufPos1.vtY + idY) % bufSz.vtY;
	  gP1.ubp = *(iBuf + bufPos0.vtY); 
	  sameRCnt = WLZ_BIT_GET(gP1.ubp, bufPos0.vtX) &&
		     WLZ_BIT_GET(gP1.ubp, bufPos1.vtX);
	  ++idY;
	}
      }
      if(sameRCnt)
      {
	for(idY = 0; idY < bufSz.vtY; ++idY)
	{
	  bufPos0.vtY = (bufPos1.vtY + idY) % bufSz.vtY;
	  switch(vType)
	  {
	    case WLZ_GREY_UBYTE:
	      gP0.ubp = (WlzUByte *)rBuf;
	      gP1.ubp = *((WlzUByte **)vBuf + bufPos0.vtY);
	      gV0.ubv = *(gP1.ubp + bufPos0.vtX);
	      gV1.ubv = *(gP1.ubp + bufPos1.vtX);
	      while(*(gP0.ubp) != gV0.ubv)
	      {
	        ++(gP0.ubp);
	      }
	      *(gP0.ubp) = gV1.ubv;
	      break;
	    case WLZ_GREY_SHORT:
	      gP0.shp = (short *)rBuf;
	      gP1.shp = *((short **)vBuf + bufPos0.vtY);
	      gV0.shv = *(gP1.shp + bufPos0.vtX);
	      gV1.shv = *(gP1.shp + bufPos1.vtX);
	      while(*(gP0.shp) != gV0.shv)
	      {
	        ++(gP0.shp);
	      }
	      *(gP0.shp) = gV1.shv;
	      break;
	    case WLZ_GREY_INT:
	      gP0.inp = (int *)rBuf;
	      gP1.inp = *((int **)vBuf + bufPos0.vtY);
	      gV0.inv = *(gP1.inp + bufPos0.vtX);
	      gV1.inv = *(gP1.inp + bufPos1.vtX);
	      while(*(gP0.inp) != gV0.inv)
	      {
	        ++(gP0.inp);
	      }
	      *(gP0.inp) = gV1.inv;
	      break;
	    case WLZ_GREY_FLOAT:
	      gP0.flp = (float *)rBuf;
	      gP1.flp = *((float **)vBuf + bufPos0.vtY);
	      gV0.flv = *(gP1.flp + bufPos0.vtX);
	      gV1.flv = *(gP1.flp + bufPos1.vtX);
	      while(*(gP0.flp) != gV0.flv)
	      {
	        ++(gP0.flp);
	      }
	      *(gP0.flp) = gV1.flv;
	      break;
	    case WLZ_GREY_DOUBLE:
	      gP0.dbp = (double *)rBuf;
	      gP1.dbp = *((double **)vBuf + bufPos0.vtY);
	      gV0.dbv = *(gP1.dbp + bufPos0.vtX);
	      gV1.dbv = *(gP1.dbp + bufPos1.vtX);
	      while(*(gP0.dbp) != gV0.dbv)
	      {
	        ++(gP0.dbp);
	      }
	      *(gP0.dbp) = gV1.dbv;
	      break;
	    case WLZ_GREY_RGBA:
	      gP0.rgbp = (WlzUInt *)rBuf;
	      gP1.rgbp = *((WlzUInt **)vBuf + bufPos0.vtY);
	      gV0.rgbv = *(gP1.rgbp + bufPos0.vtX);
	      gV1.rgbv = *(gP1.rgbp + bufPos1.vtX);
	      while(*(gP0.rgbp) != gV0.rgbv)
	      {
	        ++(gP0.rgbp);
	      }
	      *(gP0.rgbp) = gV1.rgbv;
	      break;
	    default:
	      break;
	  }
	}
      }
      else
      {
        rCnt = 0;
	for(idY = 0; idY < bufSz.vtY; ++idY)
	{
	  bufPos0.vtY = (bufPos1.vtY + idY) % bufSz.vtY;
	  switch(vType)
	  {
	    case WLZ_GREY_UBYTE:
	      gP0.ubp = *((WlzUByte **)vBuf + bufPos0.vtY);
	      break;
	    case WLZ_GREY_SHORT:
	      gP0.shp = *((short **)vBuf + bufPos0.vtY);
	      break;
	    case WLZ_GREY_INT:
	      gP0.inp = *((int **)vBuf + bufPos0.vtY);
	      break;
	    case WLZ_GREY_FLOAT:
	      gP0.flp = *((float **)vBuf + bufPos0.vtY);
	      break;
	    case WLZ_GREY_DOUBLE:
	      gP0.dbp = *((double **)vBuf + bufPos0.vtY);
	      break;
	    case WLZ_GREY_RGBA:
	      gP0.rgbp = *((WlzUInt **)vBuf + bufPos0.vtY);
	      break;
	    default:
	      break;
	  }
	  bufPos0.vtX = idK - fSz2;
	  switch(vType)
	  {
	    case WLZ_GREY_UBYTE:
	      for(idX = 0; idX < bufSz.vtY; ++idX)
	      {
		if(WLZ_BIT_GET(*(iBuf + bufPos0.vtY), bufPos0.vtX) != 0)
		{
		  *((WlzUByte *)rBuf + rCnt++) = *(gP0.ubp + bufPos0.vtX);
		}
		++(bufPos0.vtX);
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      for(idX = 0; idX < bufSz.vtY; ++idX)
	      {
		if(WLZ_BIT_GET(*(iBuf + bufPos0.vtY), bufPos0.vtX) != 0)
		{
		  *((short *)rBuf + rCnt++) = *(gP0.shp + bufPos0.vtX);
		}
		++(bufPos0.vtX);
	      }
	      break;
	    case WLZ_GREY_INT:
	      for(idX = 0; idX < bufSz.vtY; ++idX)
	      {
		if(WLZ_BIT_GET(*(iBuf + bufPos0.vtY), bufPos0.vtX) != 0)
		{
		  *((int *)rBuf + rCnt++) = *(gP0.inp + bufPos0.vtX);
		}
		++(bufPos0.vtX);
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      for(idX = 0; idX < bufSz.vtY; ++idX)
	      {
		if(WLZ_BIT_GET(*(iBuf + bufPos0.vtY), bufPos0.vtX) != 0)
		{
		  *((float *)rBuf + rCnt++) = *(gP0.flp + bufPos0.vtX);
		}
		++(bufPos0.vtX);
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      for(idX = 0; idX < bufSz.vtY; ++idX)
	      {
		if(WLZ_BIT_GET(*(iBuf + bufPos0.vtY), bufPos0.vtX) != 0)
		{
		  *((double *)rBuf + rCnt++) = *(gP0.dbp + bufPos0.vtX);
		}
		++(bufPos0.vtX);
	      }
	      break;
	  case WLZ_GREY_RGBA: /* RGBA to be done - what to do? RAB */
	    return;
	  default:
	    break;
	  }
	}
      }
    }
    /* Do the rank filtering. */
    if(rCnt > 1)
    {
      WlzRankFilterValues(gVWSp, 0, objPos.vtY, objPos.vtX,
      			  vType, rBuf, rCnt, rank);
    }
  }
}

/*!
* \return	void
* \ingroup      WlzValuesFilters
* \brief	Rank filter the buffer values for a single plane.
*		always popolate the rank buffer afresh, as quick experiments
*		have shown that trying to replace values does not speed up
*		the computation for a 3D buffer.
* \param	gObj			The object being filtered, needed 
*					for random pixel access.
* \param	gVWSp			Grey value workspace, needed for
*					random pixel access.
* \param	vBuf			Circular buffer of planes of grey
*					values.
* \param	iBuf			Circular buffer of planes of bits
*					used to encode intervals.
* \param	rBuf			A buffer used to rank pixels.
* \param	plane			The current plane, with both past
*					and future planes valid.
* \param	vType			The grey value type.
* \param	bufSz			The buffer size.
* \param	rank			Required rank with values:
*					0.0 minimum, 0.5 median and
*					1.0 maximum.
*/
static void	WlzRankFilterValPl(WlzObject *gObj, WlzGreyValueWSpace *gVWSp,
				   void ***vBuf, WlzUByte ***iBuf, void *rBuf,
				   int plane, WlzGreyType vType,
				   WlzIVertex3 bufSz, double rank)
{
  int		idK,
		idL,
		idX,
		idY,
		idZ,
		rCnt,
		fSz2;
  WlzUByte	*iP0;
  WlzGreyP	gP0;
  WlzDomain	gDom;
  WlzIVertex3	bufPos0,
		bufPos1,
  		objPos;

  gP0.v = NULL;
  gDom = gObj->domain;
  fSz2 = bufSz.vtZ / 2;
  bufPos1.vtZ = plane - gDom.p->plane1 - fSz2 + bufSz.vtZ;
  for(idL = fSz2; idL < bufSz.vtY - fSz2; ++idL)
  {
    for(idK = fSz2; idK < bufSz.vtX - fSz2; ++idK)
    {
      rCnt = 0;
      objPos.vtX = gDom.p->kol1 + idK - fSz2;
      objPos.vtY = gDom.p->line1 + idL - fSz2;
      objPos.vtZ = plane;
      bufPos0.vtX = idK;
      bufPos0.vtY = idL;
      bufPos0.vtZ = (plane - gDom.p->plane1 + bufSz.vtZ) % bufSz.vtZ;
      if(WLZ_BIT_GET(*(*(iBuf + bufPos0.vtZ) + bufPos0.vtY), bufPos0.vtX) != 0)
      {
	for(idZ = 0; idZ < bufSz.vtZ; ++idZ)
	{
	  bufPos0.vtZ = (bufPos1.vtZ + idZ) % bufSz.vtZ;
	  for(idY = 0; idY < bufSz.vtZ; ++idY)
	  {
	    bufPos0.vtY = idL + idY - fSz2;
	    switch(vType)
	    {
	      case WLZ_GREY_UBYTE:
		gP0.ubp = *(*((WlzUByte ***)vBuf + bufPos0.vtZ) + bufPos0.vtY);
		break;
	      case WLZ_GREY_SHORT:
		gP0.shp = *(*((short ***)vBuf + bufPos0.vtZ) + bufPos0.vtY);
		break;
	      case WLZ_GREY_INT:
		gP0.inp = *(*((int ***)vBuf + bufPos0.vtZ) + bufPos0.vtY);
		break;
	      case WLZ_GREY_FLOAT:
		gP0.flp = *(*((float ***)vBuf + bufPos0.vtZ) + bufPos0.vtY);
		break;
	      case WLZ_GREY_DOUBLE:
		gP0.dbp = *(*((double ***)vBuf + bufPos0.vtZ) + bufPos0.vtY);
		break;
	      case WLZ_GREY_RGBA:
		gP0.rgbp = *(*((WlzUInt ***)vBuf + bufPos0.vtZ) + bufPos0.vtY);
		break;
	      default:
	        break;
	    }
	    iP0 = *(*(iBuf + bufPos0.vtZ) + bufPos0.vtY);
	    bufPos0.vtX = idK - fSz2;
	    switch(vType)
	    {
	      case WLZ_GREY_UBYTE:
		for(idX = 0; idX < bufSz.vtZ; ++idX)
		{
		  if(WLZ_BIT_GET(iP0, bufPos0.vtX) != 0)
		  {
		    *((WlzUByte *)rBuf + rCnt) = *(gP0.ubp + bufPos0.vtX);
		    ++rCnt;
		  }
		  ++(bufPos0.vtX);
		}
		break;
	      case WLZ_GREY_SHORT:
		for(idX = 0; idX < bufSz.vtZ; ++idX)
		{
		  if(WLZ_BIT_GET(iP0, bufPos0.vtX) != 0)
		  {
		    *((short *)rBuf + rCnt) = *(gP0.shp + bufPos0.vtX);
		    ++rCnt;
		  }
		  ++(bufPos0.vtX);
		}
		break;
	      case WLZ_GREY_INT:
		for(idX = 0; idX < bufSz.vtZ; ++idX)
		{
		  if(WLZ_BIT_GET(iP0, bufPos0.vtX) != 0)
		  {
		    *((int *)rBuf + rCnt) = *(gP0.inp + bufPos0.vtX);
		    ++rCnt;
		  }
		  ++(bufPos0.vtX);
		}
		break;
	      case WLZ_GREY_FLOAT:
		for(idX = 0; idX < bufSz.vtZ; ++idX)
		{
		  if(WLZ_BIT_GET(iP0, bufPos0.vtX) != 0)
		  {
		    *((float *)rBuf + rCnt) = *(gP0.flp + bufPos0.vtX);
		  }
		  ++(bufPos0.vtX);
		}
		break;
	      case WLZ_GREY_DOUBLE:
		for(idX = 0; idX < bufSz.vtZ; ++idX)
		{
		  if(WLZ_BIT_GET(iP0, bufPos0.vtX) != 0)
		  {
		    *((double *)rBuf + rCnt) = *(gP0.dbp + bufPos0.vtX);
		    ++rCnt;
		  }
		  ++(bufPos0.vtX);
		}
		break;
	      case WLZ_GREY_RGBA:
	        /* TODO RGBA to be done - what to do? RAB */
		return;
	      default:
	        break;
	    }
	  }
	}
      }
      /* Do the rank filtering. */
      if(rCnt > 1)
      {
	WlzRankFilterValues(gVWSp, objPos.vtZ, objPos.vtY, objPos.vtX,
	    		    vType, rBuf, rCnt, rank);
      }
    }
  }
}

/*!
* \return	void
* \ingroup      WlzValuesFilters
* \brief	Replace the single pixel/voxel in the object used to
*		initialize the given grey value work space with the
*		rank filtered value.
* \param	gVWSp			Given grey value work space.
* \param	posZ			Plane position.
* \param	posY			Line position.
* \param	posX			Column position.
* \param	vType			Grey value type.
* \param	values			Buffer of values to be ranked.
* \param	nValues			Number of values in buffer.
* \param	rank			Required rank.
*/
static void	WlzRankFilterValues(WlzGreyValueWSpace *gVWSp,
				    int posZ, int posY, int posX,
				    WlzGreyType vType, void *values,
				    int nValues, double rank)
{
  int		rankI;

  rankI = (int )floor(nValues * rank);
  WlzGreyValueGet(gVWSp, posZ, posY, posX);
  switch(vType)
  {
    case WLZ_GREY_UBYTE:
      AlgRankSelectUB((WlzUByte *)values, nValues, rankI);
      *(gVWSp->gPtr[0].ubp) = *((WlzUByte *)values + rankI);
      break;
    case WLZ_GREY_SHORT:
      AlgRankSelectS((short *)values, nValues, rankI);
      *(gVWSp->gPtr[0].shp) = *((short *)values + rankI);
      break;
    case WLZ_GREY_INT:
      AlgRankSelectI((int *)values, nValues, rankI);
      *(gVWSp->gPtr[0].inp) = *((int *)values + rankI);
      break;
    case WLZ_GREY_FLOAT:
      AlgRankSelectF((float *)values, nValues, rankI);
      *(gVWSp->gPtr[0].flp) = *((float *)values + rankI);
      break;
    case WLZ_GREY_DOUBLE:
      AlgRankSelectD((double *)values, nValues, rankI);
      *(gVWSp->gPtr[0].dbp) = *((double *)values + rankI);
      break;
    default:
      break;
  }
}

/* #define WLZ_RANK_TEST */
#ifdef WLZ_RANK_TEST

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0,
		rSz = 3;
  double 	rank = 0.5;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  static char   optList[] = "ho:r:s:";
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'r':
	if((sscanf(optarg, "%lg", &rank) != 1) || (rank < 0.0) || (rank > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
	if((sscanf(optarg, "%d", &rSz) != 1) || (rSz < 0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    errNum = WlzRankFilter(obj, rSz, rank);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed to rank filter object, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outObjFileStr, "-")?
	     fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outObjFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<output file>] [-r#] [-s#] [<input file>]\n"
    	    "Rank filters the grey values of a Woolz domain object.\n"
	    "Options are:\n"
	    "  -r  Required rank. Range [0.0-1.0] with 0.0 minimum, 0.5\n"
	    "      median and 1.0 maximum value. Default 0.5.\n"
	    "  -s  Size of filter region, must be greater than zero.\n"
	    "      Default 3 for 3x3 region.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* WLZ_RANK_TEST */
