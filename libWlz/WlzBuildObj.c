#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBuildObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzBuildObj.c
* \author       Bill Hill
* \date         July 2013
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
* \brief	Functions to allow spatial domain objects with values
* 		to be built incrementally.
* \ingroup	WlzAllocation
*/

#include <Wlz.h>

static WlzObject 		*WlzBuildObj2(
				  WlzObject *cObj,
				  WlzIVertex2 og,
				  WlzIVertex2 sz,
				  WlzGreyType gType,
				  int bufSz,
				  WlzGreyP bufP,
				  WlzErrorNum *dstErr);

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Wrapper for WlzBuildObj3() with WlzUByte values.
* \param	cObj		Given current object.
* \param	og		Origin of rectangular buffer.
* \param	sz		Buffer size (note 2D).
* \param	bufSz		Number of values in the buffer
* 				(ie sz.vtX * sz.vtY).
* \param	buf		Given buffer of values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3U(WlzObject *cObj,
			       WlzIVertex3 og, WlzIVertex2 sz,
			       int bufSz, WlzUByte *buf,
			       WlzErrorNum *dstErr)
{
  WlzGreyP	bufP;

  bufP.ubp = buf;
  return(WlzBuildObj3(cObj, og, sz, WLZ_GREY_UBYTE, bufSz, bufP, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Wrapper for WlzBuildObj3() with short values.
* \param	cObj		Given current object.
* \param	og		Origin of rectangular buffer.
* \param	sz		Buffer size (note 2D).
* \param	bufSz		Number of values in the buffer
* 				(ie sz.vtX * sz.vtY).
* \param	buf		Given buffer of values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3S(WlzObject *cObj,
			       WlzIVertex3 og, WlzIVertex2 sz,
			       int bufSz, short *buf,
			       WlzErrorNum *dstErr)
{
  WlzGreyP	bufP;

  bufP.shp = buf;
  return(WlzBuildObj3(cObj, og, sz, WLZ_GREY_SHORT, bufSz, bufP, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Wrapper for WlzBuildObj3() with int values.
* \param	cObj		Given current object.
* \param	og		Origin of rectangular buffer.
* \param	sz		Buffer size (note 2D).
* \param	bufSz		Number of values in the buffer
* 				(ie sz.vtX * sz.vtY).
* \param	buf		Given buffer of values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3I(WlzObject *cObj,
			       WlzIVertex3 og, WlzIVertex2 sz,
			       int bufSz, int *buf,
			       WlzErrorNum *dstErr)
{
  WlzGreyP	bufP;

  bufP.inp = buf;
  return(WlzBuildObj3(cObj, og, sz, WLZ_GREY_INT, bufSz, bufP, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Wrapper for WlzBuildObj3() with float values.
* \param	cObj		Given current object.
* \param	og		Origin of rectangular buffer.
* \param	sz		Buffer size (note 2D).
* \param	bufSz		Number of values in the buffer
* 				(ie sz.vtX * sz.vtY).
* \param	buf		Given buffer of values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3F(WlzObject *cObj,
			       WlzIVertex3 og, WlzIVertex2 sz,
			       int bufSz, float *buf,
			       WlzErrorNum *dstErr)
{
  WlzGreyP	bufP;

  bufP.flp = buf;
  return(WlzBuildObj3(cObj, og, sz, WLZ_GREY_FLOAT, bufSz, bufP, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Wrapper for WlzBuildObj3() with double values.
* \param	cObj		Given current object.
* \param	og		Origin of rectangular buffer.
* \param	sz		Buffer size (note 2D).
* \param	bufSz		Number of values in the buffer
* 				(ie sz.vtX * sz.vtY).
* \param	buf		Given buffer of values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3D(WlzObject *cObj,
			       WlzIVertex3 og, WlzIVertex2 sz,
			       int bufSz, double *buf,
			       WlzErrorNum *dstErr)
{
  WlzGreyP	bufP;

  bufP.dbp = buf;
  return(WlzBuildObj3(cObj, og, sz, WLZ_GREY_DOUBLE, bufSz, bufP, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Wrapper for WlzBuildObj3() with RGBA values.
* \param	cObj		Given current object.
* \param	og		Origin of rectangular buffer.
* \param	sz		Buffer size (note 2D).
* \param	bufSz		Number of values in the buffer
* 				(ie sz.vtX * sz.vtY).
* \param	buf		Given buffer of values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3R(WlzObject *cObj,
			       WlzIVertex3 og, WlzIVertex2 sz,
			       int bufSz, WlzUInt *buf,
			       WlzErrorNum *dstErr)
{
  WlzGreyP	bufP;

  bufP.rgbp = buf;
  return(WlzBuildObj3(cObj, og, sz, WLZ_GREY_RGBA, bufSz, bufP, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Wrapper for WlzBuildObj3() with buffer passed as bytes.
* \param	cObj		Given current object.
* \param	og		Origin of rectangular buffer.
* \param	sz		Buffer size (note 2D).
* \param	gType		The grey type.
* \param	bufSz		Number of values in the buffer
* 				(ie sz.vtX * sz.vtY * sizeof(<gType>))).
* \param	buf		Given buffer of values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3B(WlzObject *cObj,
			       WlzIVertex3 og, WlzIVertex2 sz,
			       WlzGreyType gType,
			       int bufSz, char *buf,
			       WlzErrorNum *dstErr)
{
  int		gSz,
  		bufPSz = 0;
  WlzGreyP	bufP;
  WlzObject	*nObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bufP.ubp = (WlzUByte *)buf;
  if((gSz = WlzGreySize(gType)) <= 0)
  {
    errNum = WLZ_ERR_GREY_TYPE;
  }
  else
  {
    bufPSz = bufSz / gSz;
    nObj = WlzBuildObj3(cObj, og, sz, gType, bufPSz, bufP, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nObj);
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Creates a new or adds to an existing 3D spatial domain
*  		object using a rectangular buffer of values to a single
*  		plane of the given current object (which may be NULL or
*  		empty if the object is to be created).
* 		The returned object will share domains and values of
* 		planes other than the given plane with the current object.
* \param	cObj			Given current object.
* \param	og			Origin of rectangular buffer.
* \param	sz			Buffer size (note 2D).
* \param	gType			Grey type which must be consistent
* 					with the current object (if it is
* 					valid) and the buffer of values.
* \param	bufSz			Number of values in the buffer
* 					(ie sz.vtX * sz.vtY).
* \param	bufP			Given buffer of values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzBuildObj3(WlzObject *cObj,
			      WlzIVertex3 og, WlzIVertex2 sz,
			      WlzGreyType gType,
			      int bufSz, WlzGreyP bufP,
			      WlzErrorNum *dstErr)
{
  int		nPlnReq = 1;
  WlzDomain	cDom,
  		nDom;
  WlzValues	cVal,
  		nVal;
  WlzObject	*nObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  cDom.core = NULL;
  nDom.core = NULL;
  cVal.core = NULL;
  nVal.core = NULL;
  if(cObj)
  {
    WlzGreyType cGType = WLZ_GREY_ERROR;;

    switch(cObj->type)
    {
      case WLZ_EMPTY_OBJ:
	cObj = NULL;
        break;
      case WLZ_3D_DOMAINOBJ:
	if(cObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(cObj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else
	{
	  cDom = cObj->domain;
	  cVal = cObj->values;
          cGType = WlzGreyTypeFromObj(cObj, &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && (cGType != gType))
	{
	  errNum = WLZ_ERR_GREY_TYPE;
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Create a new object with new plane domain and voxel values. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzIBox3	nBox;
    WlzPixelV	bgdV;
    float	vxSz[3];

    nBox.xMin = og.vtX;
    nBox.yMin = og.vtY;
    nBox.zMin = og.vtZ;
    nBox.xMax = og.vtX + sz.vtX - 1;
    nBox.yMax = og.vtY + sz.vtY - 1;
    nBox.zMax = og.vtZ;
    if(cObj)
    {
      nPlnReq = (og.vtZ < cDom.p->plane1) ||
                (og.vtZ > cDom.p->lastpl) ||
                ((*(cDom.p->domains + og.vtZ - cDom.p->plane1)).core == NULL);
      nBox.xMin = ALG_MIN(nBox.xMin, cDom.p->kol1);
      nBox.yMin = ALG_MIN(nBox.yMin, cDom.p->line1);
      nBox.zMin = ALG_MIN(nBox.zMin, cDom.p->plane1);
      nBox.xMax = ALG_MAX(nBox.xMax, cDom.p->lastkl);
      nBox.yMax = ALG_MAX(nBox.yMax, cDom.p->lastln);
      nBox.zMax = ALG_MAX(nBox.zMax, cDom.p->lastpl);
      vxSz[0] = cDom.p->voxel_size[0];
      vxSz[1] = cDom.p->voxel_size[1];
      vxSz[2] = cDom.p->voxel_size[2];
      bgdV = WlzGetBackground(cObj, &errNum);
    }
    else
    {
      vxSz[0] = vxSz[1] = vxSz[2] = 1.0f;
      bgdV.type = WLZ_GREY_INT;
      bgdV.v.inv = 0;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nDom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				  nBox.zMin, nBox.zMax,
				  nBox.yMin, nBox.yMax,
				  nBox.xMin, nBox.xMax,
				  &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
				     nBox.zMin, nBox.zMax,
				     bgdV, NULL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nDom.p->voxel_size[0] = vxSz[0];
      nDom.p->voxel_size[1] = vxSz[1];
      nDom.p->voxel_size[2] = vxSz[2];
      nObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, nDom, nVal, NULL, NULL, &errNum);
    }
  }
  /* Set the domain and values on each plane for the new object. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idZ;

    for(idZ = nDom.p->plane1; idZ <= nDom.p->lastpl; ++idZ)
    {
      int	idP;
      WlzDomain	nDom2;
      WlzValues nVal2;

      nDom2.core = NULL;
      nVal2.core = NULL;
      idP = idZ - nDom.p->plane1;
      if(idZ == og.vtZ)
      {
        /* Plane with buffer. */
	WlzIVertex2 og2,
		    sz2;
	WlzObject   *cObj2 = NULL,
		    *nObj2 = NULL;

	og2.vtX = og.vtX;
	og2.vtY = og.vtY;
	sz2.vtX = sz.vtX;
	sz2.vtY = sz.vtY;
        if(nPlnReq == 0)
	{
	  int	idP;

	  idP = idZ - cDom.p->plane1;
	  cObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, 
			      *(cDom.p->domains  + idP),
			      *(cVal.vox->values + idP),
			      NULL, NULL, &errNum);
	}
	nObj2 = WlzBuildObj2(cObj2, og2, sz2, gType, bufSz, bufP, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  nDom2 = WlzAssignDomain(nObj2->domain, NULL);
	  nVal2 = WlzAssignValues(nObj2->values, NULL);
	}
	(void )WlzFreeObj(cObj2);
	(void )WlzFreeObj(nObj2);
      }
      else if((idZ >= cDom.p->plane1) && (idZ <= cDom.p->lastpl))
      {
	/* Not buffer plane, but previously existing plane. */
	int	idQ;

	idQ = idZ - cDom.p->plane1;
	nDom2 = WlzAssignDomain(*(cDom.p->domains  + idQ), NULL);
	nVal2 = WlzAssignValues(*(cVal.vox->values + idQ), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        *(nDom.p->domains  + idP) = nDom2;
        *(nVal.vox->values + idP) = nVal2;
      }
      else
      {
        break;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nObj);
}

/*!
* \return	New Woolz object.
* \ingroup	WlzAllocation
* \brief	Creates a new 2D spatial domain object by adding a
* 		rectangular buffer of values to the given current
* 		object (which may be NULL or empty).
* \param	cObj			Given current object.
* \param	og			Origin of rectangular buffer.
* \param	sz			Buffer size.
* \param	gType			Grey type which must be consistent
* 					with the current object and the
* 					buffer of values.
* \param	bufSz			Number of values in the buffer.
* \param	bufP			Given buffer of values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzBuildObj2(WlzObject *cObj,
			       WlzIVertex2 og, WlzIVertex2 sz,
			       WlzGreyType gType, int bufSz, WlzGreyP bufP,
			       WlzErrorNum *dstErr)
{
  WlzDomain	bDom;
  WlzValues	bVal,
		nVal;
  WlzObject	*bObj = NULL,
  		*nObj = NULL;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bDom.core = NULL;
  bVal.core = NULL;
  nVal.core = NULL;
  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = 0;
  if(cObj)
  {
    WlzGreyType cGType = WLZ_GREY_ERROR;;

    switch(cObj->type)
    {
      case WLZ_EMPTY_OBJ:
	cObj = NULL;
        break;
      case WLZ_2D_DOMAINOBJ:
	if(cObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(cObj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else
	{
          cGType = WlzGreyTypeFromObj(cObj, &errNum);
          bgdV = WlzGetBackground(cObj, &errNum);
	}
	if((errNum == WLZ_ERR_NONE) && (cGType != gType))
	{
	  errNum = WLZ_ERR_GREY_TYPE;
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Create new object with domain and values of given rectangular buffer. */
  if(errNum == WLZ_ERR_NONE)
  {
    bDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				   og.vtY, og.vtY + sz.vtY - 1,
				   og.vtX, og.vtX + sz.vtX - 1,
				   &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObjectType gTT;

    gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RECT, gType, NULL);
    bVal.r = WlzMakeRectValueTb(gTT, bDom.i->line1, bDom.i->lastln,
			        bDom.i->kol1, sz.vtX, bgdV, bufP.inp, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, bDom, bVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(cObj == NULL)
    {
      /* Just copy the buffer object. */
      nObj = WlzCopyObject(bObj, &errNum);
    }
    else
    {
      /* Compute union of current and buffer objects. */
      nObj = (cObj)? WlzUnion2(cObj, bObj, &errNum):
		     WlzMakeMain(WLZ_2D_DOMAINOBJ, bDom, nVal,
				 NULL, NULL, &errNum);
      /* Create new value table. */
      if(errNum == WLZ_ERR_NONE)
      {
	WlzObjectType gTT;

	gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, gType, NULL);
	nVal.v = WlzNewValueTb(nObj, gTT, bgdV, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	nObj->values = WlzAssignValues(nVal, NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	WlzObject 	*tObj;

	/* Copy existing values to new object. */
	tObj = WlzGreyTransfer(nObj, cObj, 0, &errNum);
	(void )WlzFreeObj(nObj);
	nObj = tObj;
	/* Then copy buffer values to new object. */
	if(errNum == WLZ_ERR_NONE)
	{
	  tObj = WlzGreyTransfer(nObj, bObj, 0, &errNum);
	  (void )WlzFreeObj(nObj);
	  nObj = tObj;
	}
      }
    }
  }
  (void )WlzFreeObj(bObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nObj);
}
