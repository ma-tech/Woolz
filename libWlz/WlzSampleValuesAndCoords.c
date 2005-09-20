#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzSampleValuesAndCoords.c
* \author       Bill Hill
* \date         November 2002
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Extracts values and coordinates from a Woolz object
* 		with grey values using a sampoling function.
* \ingroup	WlzFeatures
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <Wlz.h>

static WlzErrorNum 		WlzSampleValuesAndCoords2D(
				  WlzObject *obj,
				  WlzGreyType *dstGType, int *dstNVal,
				  WlzGreyP *dstValP,
				  WlzIVertex2 **dstCoords,
				  WlzSampleFn samFn, int samFac);

/*!
* \return
* \ingroup	WlzFeatures
* \brief	Allocates buffers for both the grey values and the
*		coordinates of the grey values in the given object.
*		On return these buffers contain the sampled object
*		values and the coordinates of the values.
* \param	obj			Given object.
* \param	dstGType		Type of grey value.
* \param	dstNVal			Number of grey values, also the
*					number of coordinates.
* \param	dstValP			Destination pointer for the values.
* \param	dstCoords		Destination pointer for the
* 					coordinates, these are always either
* 					WlzIVertex2 or WlzIVertex3 depending on
* 					the objects dimension.
* \param	samFn			Sampling function.
* \param	samFac			Sampling factor.
*/
WlzErrorNum	WlzSampleValuesAndCoords(WlzObject *obj,
				   WlzGreyType *dstGType, int *dstNVal,
				   WlzGreyP *dstValP, WlzVertexP *dstCoords,
				   WlzSampleFn samFn, int samFac)
{
  WlzIVertex2	*coords2D = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstGType == NULL) || (dstNVal == NULL) ||
     (dstValP == NULL) || (dstCoords == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(samFac < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(obj == NULL)
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
    switch(samFn)
    {
      case WLZ_SAMPLEFN_MIN: /* FALLTHROUGH */
      case WLZ_SAMPLEFN_MAX:
        break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	errNum = WlzSampleValuesAndCoords2D(obj, dstGType, dstNVal,
					    dstValP, &coords2D,
					    samFn, samFac);
	dstCoords->i2 = coords2D;
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return
* \ingroup	WlzFeatures
* \brief	Samples the given 2D object and returns buffers, allocated
*		in this function, with object values and the coordinates of
*		of the values.
*		This function assumes that all it's parameters are valid.
* \param	obj			Given object.
* \param	dstGType		Type of grey value.
* \param	dstNVal			Destination pointer for the number
*					of grey values, also the number
* 					of coordinates.
* \param	dstValP			Destination pointer for the values.
* \param	dstCoords		Destination pointer for the
*					coordinates.
* \param	samFn			Sampling function.
* \param	samFac			Sampling factor.
*/
static WlzErrorNum WlzSampleValuesAndCoords2D(WlzObject *obj,
				        WlzGreyType *dstGType, int *dstNVal,
					WlzGreyP *dstValP,
					WlzIVertex2 **dstCoords,
					WlzSampleFn samFn, int samFac)
{
  int		tI0,
  		gSz,
		iPos,
		bIdx,
		bWidth,
		aSz,
  		nVal = 0;
  WlzGreyType	gType;
  WlzIVertex2	bPos,
  		bStart;
  UBYTE		*dBuf = NULL;
  WlzIVertex2	*cBuf = NULL,
		*cAry = NULL;
  WlzGreyP	gVal,
  		vBuf,
		vAry;

  WlzIntervalWSpace iWsp;
  WlzGreyWSpace gWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vBuf.inp = NULL;
  gType = WlzGreyTypeFromObj(obj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gType)
    {
      case WLZ_GREY_INT:
	gSz = sizeof(int);
        break;
      case WLZ_GREY_SHORT:
	gSz = sizeof(short);
        break;
      case WLZ_GREY_UBYTE:
	gSz = sizeof(UBYTE);
        break;
      case WLZ_GREY_FLOAT:
	gSz = sizeof(float);
        break;
      case WLZ_GREY_DOUBLE:
	gSz = sizeof(double);
        break;
      case WLZ_GREY_RGBA:
	gSz = sizeof(UINT);
        break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  /* Allocate arrays for the sampled values and coordinates. If the size
   * estimated is wrong then either a reallocation will be needed or space
   * will be wasted. */
  if(errNum == WLZ_ERR_NONE)
  {
    aSz = WlzArea(obj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tI0 = (samFac * samFac);
    if(samFac > 1)
    {
      aSz = ((2 * aSz) / (samFac * samFac)) + 1000;
    }
    if(((cAry = (WlzIVertex2 *)
                AlcMalloc(aSz * sizeof(WlzIVertex2))) == NULL) ||
       ((vAry.inp = (int *)AlcMalloc(aSz * gSz)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Allocate working buffers large enough for any line of sampled values or
   * coordinates. */
  if(errNum == WLZ_ERR_NONE)
  {
    bStart.vtX = obj->domain.i->kol1 - (obj->domain.i->kol1 % samFac);
    bStart.vtY = obj->domain.i->line1 - (obj->domain.i->line1 % samFac);
    bWidth = 3 + ((obj->domain.i->lastkl - obj->domain.i->kol1) / samFac);
    if(((cBuf = (WlzIVertex2 *)
    	        AlcMalloc(bWidth * sizeof(WlzIVertex2))) == NULL) ||
       ((vBuf.inp = (int *)AlcMalloc(bWidth * gSz)) == NULL) ||
       ((dBuf = (UBYTE *)AlcCalloc(bWidth, sizeof(UBYTE))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Work through the object's intervals sampling the values. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(obj, &iWsp, &gWsp);
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE))
    {
      gVal = gWsp.u_grintptr;
      bPos.vtX = iWsp.lftpos - bStart.vtX;
      bPos.vtY = iWsp.linpos - bStart.vtY;
      switch(gWsp.pixeltype)
      {
        case WLZ_GREY_INT:
	  switch(samFn)
	  {
	    case WLZ_SAMPLEFN_MIN:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.inp) < *(vBuf.inp + bIdx)))
		{
		  *(vBuf.inp + bIdx) = *(gVal.inp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.inp);
	      }
	      break;
	    case WLZ_SAMPLEFN_MAX:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.inp) > *(vBuf.inp + bIdx)))
		{
		  *(vBuf.inp + bIdx) = *(gVal.inp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.inp);
	      }
	      break;
	  }
	  break;
	case WLZ_GREY_SHORT:
	  switch(samFn)
	  {
	    case WLZ_SAMPLEFN_MIN:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.shp) < *(vBuf.shp + bIdx)))
		{
		  *(vBuf.shp + bIdx) = *(gVal.shp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.shp);
	      }
	      break;
	    case WLZ_SAMPLEFN_MAX:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.shp) > *(vBuf.shp + bIdx)))
		{
		  *(vBuf.shp + bIdx) = *(gVal.shp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.shp);
	      }
	      break;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  switch(samFn)
	  {
	    case WLZ_SAMPLEFN_MIN:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.ubp) < *(vBuf.ubp + bIdx)))
		{
		  *(vBuf.ubp + bIdx) = *(gVal.ubp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.ubp);
	      }
	      break;
	    case WLZ_SAMPLEFN_MAX:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.ubp) > *(vBuf.ubp + bIdx)))
		{
		  *(vBuf.ubp + bIdx) = *(gVal.ubp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.ubp);
	      }
	      break;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  switch(samFn)
	  {
	    case WLZ_SAMPLEFN_MIN:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.flp) < *(vBuf.flp + bIdx)))
		{
		  *(vBuf.flp + bIdx) = *(gVal.flp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.flp);
	      }
	      break;
	    case WLZ_SAMPLEFN_MAX:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.flp) > *(vBuf.flp + bIdx)))
		{
		  *(vBuf.flp + bIdx) = *(gVal.flp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.flp);
	      }
	      break;
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  switch(samFn)
	  {
	    case WLZ_SAMPLEFN_MIN:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.dbp) < *(vBuf.dbp + bIdx)))
		{
		  *(vBuf.dbp + bIdx) = *(gVal.dbp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.dbp);
	      }
	      break;
	    case WLZ_SAMPLEFN_MAX:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.dbp) > *(vBuf.dbp + bIdx)))
		{
		  *(vBuf.dbp + bIdx) = *(gVal.dbp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.dbp);
	      }
	      break;
	  }
	  break;
	case WLZ_GREY_RGBA:
	  switch(samFn)
	  {
	    case WLZ_SAMPLEFN_MIN:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.rgbp) < *(vBuf.rgbp + bIdx)))
		{
		  *(vBuf.rgbp + bIdx) = *(gVal.rgbp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.rgbp);
	      }
	      break;
	    case WLZ_SAMPLEFN_MAX:
	      for(iPos = iWsp.lftpos; iPos <= iWsp.rgtpos; ++iPos)
	      {
		bIdx = bPos.vtX / samFac;
		if((*(dBuf + bIdx) == 0) || (*(gVal.rgbp) > *(vBuf.rgbp + bIdx)))
		{
		  *(vBuf.rgbp + bIdx) = *(gVal.rgbp);
		  (cBuf + bIdx)->vtX = iPos;
		  (cBuf + bIdx)->vtY = iWsp.linpos;
		  *(dBuf + bIdx) = 1;
		}
		++(bPos.vtX);
		++(gVal.rgbp);
	      }
	      break;
	  }
	  break;
      }
      /* Transfer values and coordinates from the buffers to the arrays.
       * It might be necessary to reallocate the arrays. */
      if((iWsp.intrmn == 0) &&
         ((iWsp.linrmn == 0) || (((bPos.vtY + 1) % samFac) == 0)))
      {
	if(aSz <= (nVal + bWidth))
	{
	  aSz *= 2;
	  if(((cAry = (WlzIVertex2 *)
		      AlcRealloc(cAry, aSz * sizeof(WlzIVertex2))) == NULL) ||
	     ((vAry.inp = (int *)AlcRealloc(vAry.inp, aSz * gSz)) == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  bIdx = 0;
	  switch(gType)
	  {
	    case WLZ_GREY_INT:
	      while(bIdx < bWidth)
	      {
		if(*(dBuf + bIdx))
		{
	          *(vAry.inp + nVal) = *(vBuf.inp + bIdx);
	          *(cAry + nVal) = *(cBuf + bIdx);
		  ++nVal;
		}
		++bIdx;
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      while(bIdx < bWidth)
	      {
		if(*(dBuf + bIdx))
		{
	          *(vAry.shp + nVal) = *(vBuf.shp + bIdx);
	          *(cAry + nVal) = *(cBuf + bIdx);
		  ++nVal;
		}
		++bIdx;
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      while(bIdx < bWidth)
	      {
		if(*(dBuf + bIdx))
		{
	          *(vAry.ubp + nVal) = *(vBuf.ubp + bIdx);
	          *(cAry + nVal) = *(cBuf + bIdx);
		  ++nVal;
		}
		++bIdx;
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      while(bIdx < bWidth)
	      {
		if(*(dBuf + bIdx))
		{
	          *(vAry.flp + nVal) = *(vBuf.flp + bIdx);
	          *(cAry + nVal) = *(cBuf + bIdx);
		  ++nVal;
		}
		++bIdx;
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      while(bIdx < bWidth)
	      {
		if(*(dBuf + bIdx))
		{
	          *(vAry.dbp + nVal) = *(vBuf.dbp + bIdx);
	          *(cAry + nVal) = *(cBuf + bIdx);
		  ++nVal;
		}
		++bIdx;
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      while(bIdx < bWidth)
	      {
		if(*(dBuf + bIdx))
		{
	          *(vAry.rgbp + nVal) = *(vBuf.rgbp + bIdx);
	          *(cAry + nVal) = *(cBuf + bIdx);
		  ++nVal;
		}
		++bIdx;
	      }
	      break;
	  }
	  if(iWsp.linrmn)
	  {
	    (void )memset(dBuf, 0, bWidth);
	  }
	}
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstGType = gType;
    *dstNVal = nVal;
    dstValP->inp = vAry.inp;
    *dstCoords = cAry;
  }
  /* Free up the working buffers. */
  AlcFree(cBuf);
  AlcFree(dBuf);
  AlcFree(vBuf.inp);
  return(errNum);
}
