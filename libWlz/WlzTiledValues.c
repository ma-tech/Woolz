#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTiledValues_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTiledValues.c
* \author       Bill Hill
* \date         March 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
* \brief	Functions construct and convert tiled objects.
* \ingroup	WlzAllocation
*/

#include <stdlib.h>
#include <limits.h>
#include <Wlz.h>

#define WLZ_USE_MMAP
#ifdef WLZ_USE_MMAP
#include <unistd.h>
#include <sys/mman.h>
#endif /* WLZ_USE_MMAP */

static WlzObject  		*WlzMakeTiledValuesObj2D(
				  WlzObject *gObj,
				  size_t tileSz,
				  int setTiles,
				  WlzGreyType gType,
				  WlzPixelV bgdV,
				  WlzErrorNum *dstErr);
static WlzObject  		*WlzMakeTiledValuesObj3D(
				  WlzObject *gObj,
				  size_t tileSz,
				  int setTiles,
				  WlzGreyType gType,
				  WlzPixelV bgdV,
				  WlzErrorNum *dstErr);

/*!
* \return	New tiled values.
* \ingroup	WlzAllocation
* \brief	Allocates a new tiled value table, but without allocating
* 		any indices or tiles.
* int		dim			Dimension.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzTiledValues	*WlzMakeTiledValues(int dim, WlzErrorNum *dstErr)
{
  WlzTiledValues *tVal = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(dim < 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(((tVal = (WlzTiledValues *)
                   AlcCalloc(sizeof(WlzTiledValues), 1)) == NULL) ||
          ((tVal->nIdx = (int *)AlcCalloc(sizeof(int), dim)) == NULL))
  {
    AlcFree(tVal);
    tVal = NULL;
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    tVal->fd = -1;
    tVal->type = WLZ_VALUETABLE_TILED_INT;
    tVal->dim = dim;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tVal);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Frees a tiled value table.
* \param	tVal			The tiled value table to free.
*/
WlzErrorNum	WlzFreeTiledValues(WlzTiledValues *tVal)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tVal == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(WlzGreyTableIsTiled(tVal->type) == 0)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    if(WlzUnlink(&(tVal->linkcount), &errNum))
    {
      AlcFree(tVal->indices);
      AlcFree(tVal->nIdx);
      if(tVal->tiles.v)
      {
#ifdef WLZ_USE_MMAP
	if(tVal->fd >= 0)
	{
	  size_t        gSz,
	  		tSz;
	  WlzGreyType	gType;

	  gType = WlzGreyTableTypeToGreyType(tVal->type, &errNum);
	  gSz = WlzGreySize(gType);
	  tSz = tVal->numTiles * tVal->tileSz;
	  (void )munmap(tVal->tiles.v, tSz * gSz);
	  (void )close(tVal->fd);
	}
	else
	{
#endif /* WLZ_USE_MMAP */
          AlcFree(tVal->tiles.v);
#ifdef WLZ_USE_MMAP
	}
#endif /* WLZ_USE_MMAP */
      }
      AlcFree(tVal);
    }
  }
  return(errNum);
}

/*!
* \return	New tiled object or NULL on error.
* \ingroup	WlzAllocation
* \brief	Creates a tiled object from the given object.
*
* 		The tile size specifies the number of values that the
* 		tiles will contain and must be an integral power of two.
* 		The actual size in bytes of the tiles will vary with
* 		grey type.
* \param	gObj			Given object which must have an
* 					appropriate type for a tiled
* 					object. The valid types are
* 					currently WLZ_2D_DOMAINOBJ and
* 					WLZ_3D_DOMAINOBJ objects with values.
* \param	tileSz			The required tile size.
* \param	copyValues		Non zero if the grey values should
* 					be copied to the tiled object.
* \param	gType			Required grey type for values table.
* \param	bgdV			required background value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzMakeTiledValuesFromObj(WlzObject *gObj, size_t tileSz,
			        int copyValues, WlzGreyType gType,
				WlzPixelV bgdV, WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	rObj = WlzMakeTiledValuesObj2D(gObj, tileSz, copyValues, gType,
				       bgdV, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	rObj = WlzMakeTiledValuesObj3D(gObj, tileSz, copyValues, gType,
				       bgdV, &errNum);
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
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Allocates tiles for the tiles values but does not set
* 		the values within the tiles.
* \param	tVal			Given tiled values.
*/
WlzErrorNum	WlzMakeTiledValuesTiles(WlzTiledValues *tVal)
{
  size_t	gSz;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gType = WlzGreyTableTypeToGreyType(tVal->type, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    gSz = WlzGreySize(gType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tVal->fd = -1;
    if((tVal->tiles.v = AlcMalloc(tVal->numTiles *
                                  tVal->tileSz * gSz)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  return(errNum);
}

/*!
* \return	New tiled object or NULL on error.
* \ingroup	WlzAllocation
* \brief	Creates a 2D domain object which has the same domain as the
* 		given 2D domain object but a tiled value table.
*
* 		The tile size specifies the number of values that the
* 		tiles will contain and must be an integral power of two.
* 		The actual size in bytes of the tiles will vary with
* 		grey type.
*
* 		If the given object must have valid values but the
* 		tiles these are copied to are not allocated if the
* 		setTiles flag is not set.
* \param	gObj			Given object known to be valid.
* \param	tileSz			The required tile size.
* \param	setTiles		Flag which must be set for the tiles
* 					to be allocated and their values set.
* \param	gType			Required grey type for values table.
* \param	bgdV			required background value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject  *WlzMakeTiledValuesObj2D(WlzObject *gObj, size_t tileSz,
				int setTiles, WlzGreyType gType,
				WlzPixelV bgdV, WlzErrorNum *dstErr)
{
  size_t	width;
  WlzObject	*tObj = NULL;
  WlzTiledValues *tVal = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->type != WLZ_2D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    errNum = WlzValueConvertPixel(&bgdV, bgdV, gType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    unsigned int tmp;
    
    width = AlgBitNextPowerOfTwo(&tmp, tileSz);
    if(tmp == tileSz)
    {
      tmp = width / 2;
      width = 1 << tmp;
      if(width * width != tileSz)
      {
        errNum = WLZ_ERR_PARAM_DATA;
      }
    }
    else
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Create a new tiled values data structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    tVal = WlzMakeTiledValues(2, &errNum);
  }
  /* Set up the fields of the tiled values data structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject *dil = NULL,
	      *idx = NULL,
	      *scl = NULL,
	      *sft = NULL,
    	      *str = NULL;
    WlzValues val;
    WlzDomain dom;
    WlzAffineTransform *tr = NULL;

    dom.core = NULL;
    val.core = NULL;
    tVal->type = WlzGreyTableType(WLZ_GREY_TAB_TILED, gType, NULL);
    tVal->dim = 2;
    tVal->kol1 = gObj->domain.i->kol1;
    tVal->lastkl = gObj->domain.i->lastkl;
    tVal->line1 = gObj->domain.i->line1;
    tVal->lastln = gObj->domain.i->lastln;
    tVal->bckgrnd = bgdV;
    tVal->tileSz = tileSz;
    tVal->tileWidth = width;
    /* Create the tile index. */
    dom = WlzShiftDomain(gObj->type, gObj->domain,
                         -(gObj->domain.i->kol1), -(gObj->domain.i->line1), 0,
                         &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      sft = WlzMakeMain(gObj->type, dom, val, NULL, NULL, &errNum);
      if(sft == NULL)
      {
        (void )WlzFreeDomain(dom);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      str = WlzAssignObject(
	    WlzMakeRectangleObject(width - 1.0, width - 1.0,
				   0.0, 0.0, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      dil = WlzAssignObject(
	    WlzStructDilation(sft, str, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      tr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      double inv;

      inv = 1.0 / (double )width;
      (void )WlzAffineTransformScaleSet(tr, inv, inv, inv);
      scl = WlzAssignObject(
	    WlzAffineTransformObj(dil, tr, WLZ_INTERPOLATION_NEAREST,
				  &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzPixelV bkg;
      WlzObjectType tabType;

      bkg.v.inv = -1;
      bkg.type = WLZ_GREY_INT;
      tabType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
      val.v = WlzNewValueTb(scl, tabType, bkg, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      idx = WlzMakeMain(WLZ_2D_DOMAINOBJ, scl->domain, val,
			NULL, NULL, &errNum);
    }
    (void )WlzFreeObj(sft);
    (void )WlzFreeObj(str);
    (void )WlzFreeObj(dil);
    (void )WlzFreeObj(scl);
    if(idx == NULL)
    {
      (void )WlzFreeValueTb(val.v);
    }
    /* Set tile indices. */
    if(errNum == WLZ_ERR_NONE)
    {
      int tileCnt = 0;

      errNum = WlzGreySetIncValues(idx, &tileCnt);
      tVal->numTiles = tileCnt;
    }
    /* Create index array from the index object. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	**idxArray = NULL;
      WlzIVertex2 sz,
      		  org;

      org.vtX = org.vtY = 0;
      sz.vtX = idx->domain.i->lastkl - idx->domain.i->kol1 + 1;
      sz.vtY = idx->domain.i->lastln - idx->domain.i->line1 + 1;
      errNum = WlzToArray2D((void ***)&idxArray, idx, sz, org, 0,
      			    WLZ_GREY_INT);
      tVal->nIdx[0] = sz.vtX;
      tVal->nIdx[1] = sz.vtY;
      tVal->indices = *idxArray;
      AlcFree(idxArray);
    }
    (void )WlzFreeObj(idx);
  }
  /* Allocate tiles if the given object has values. */
  if((errNum == WLZ_ERR_NONE) && (setTiles != 0))
  {
    errNum = WlzMakeTiledValuesTiles(tVal);
  }
  /* Make the object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues val;

    val.t = tVal;
    tObj = WlzMakeMain(gObj->type, gObj->domain, val, NULL, NULL, &errNum);
  }
  /* Copy the values to the object with the tiles values. */
  if((errNum == WLZ_ERR_NONE) && (setTiles != 0))
  {
    errNum = WlzCopyObjectGreyValues(tObj, gObj);
  }
  /* Clear up on error. */
  if(errNum != WLZ_ERR_NONE)
  {
    if(tObj != NULL)
    {
      WlzFreeObj(tObj);
    }
    else if(tVal != NULL)
    {
      (void )AlcFree(tVal->indices);
      (void )AlcFree(tVal->tiles.v);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	New tiled object or NULL on error.
* \ingroup	WlzAllocation
* \brief	Creates a 3D domain object which has the same domain as the
* 		given 3D domain object but a tiled value table.
*
* 		The tile size specifies the number of values that the
* 		tiles will contain and must be an integral power of two.
* 		The actual size in bytes of the tiles will vary with
* 		grey type.
* \param	gObj			Given object known to be valid.
* \param	tileSz			The required tile size.
* \param	setTiles		Flag which must be set for the tiles
* 					to be allocated and their values set.
* \param	gType			Required grey type for values table.
* \param	bgdV			required background value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject  *WlzMakeTiledValuesObj3D(WlzObject *gObj, size_t tileSz,
				           int setTiles, WlzGreyType gType,
					   WlzPixelV bgdV, WlzErrorNum *dstErr)
{
  size_t	width;
  WlzObject	*tObj = NULL;
  WlzTiledValues *tVal = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzValueConvertPixel(&bgdV, bgdV, gType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    unsigned int tmp;
    
    width = AlgBitNextPowerOfTwo(&tmp, tileSz);
    if(tmp == tileSz)
    {
      tmp = width / 3;
      width = 1 << tmp;
      if(width * width * width != tileSz)
      {
        errNum = WLZ_ERR_PARAM_DATA;
      }
    }
    else
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Create a new tiled values data structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    tVal = WlzMakeTiledValues(3, &errNum);
  }
  /* Set up the fields of the tiled values data structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject *dil = NULL,
	      *idx = NULL,
	      *scl = NULL,
	      *sft = NULL,
    	      *str = NULL;
    WlzValues val;
    WlzDomain dom;
    WlzAffineTransform *tr = NULL;

    dom.core = NULL;
    val.core = NULL;
    tVal->type = WlzGreyTableType(WLZ_GREY_TAB_TILED, gType, NULL);
    tVal->dim = 3;
    tVal->kol1 = gObj->domain.p->kol1;
    tVal->lastkl = gObj->domain.p->lastkl;
    tVal->line1 = gObj->domain.p->line1;
    tVal->lastln = gObj->domain.p->lastln;
    tVal->plane1 = gObj->domain.p->plane1;
    tVal->lastpl = gObj->domain.p->lastpl;
    tVal->bckgrnd = bgdV;
    tVal->tileSz = tileSz;
    tVal->tileWidth = width;
    /* Create the tile index. */
    dom = WlzShiftDomain(gObj->type, gObj->domain,
                         -(gObj->domain.p->kol1),
			 -(gObj->domain.p->line1),
			 -(gObj->domain.p->plane1),
                         &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      sft = WlzMakeMain(gObj->type, dom, val, NULL, NULL, &errNum);
      if(sft == NULL)
      {
        (void )WlzFreeDomain(dom);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      str = WlzAssignObject(
	    WlzMakeCuboidObject(gObj->type,
	                        width - 1.0, width - 1.0, width - 1.0,
				0.0, 0.0, 0.0, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      dil = WlzAssignObject(
	    WlzStructDilation(sft, str, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      tr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      double inv;

      inv = 1.0 / (double )width;
      (void )WlzAffineTransformScaleSet(tr, inv, inv, inv);
      scl = WlzAssignObject(
	    WlzAffineTransformObj(dil, tr, WLZ_INTERPOLATION_NEAREST,
				  &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzPixelV bkg;
      WlzObjectType tabType;

      bkg.v.inv = -1;
      bkg.type = WLZ_GREY_INT;
      tabType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
      val.vox = WlzNewValuesVox(scl, tabType, bkg, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      idx = WlzMakeMain(WLZ_3D_DOMAINOBJ, scl->domain, val,
			NULL, NULL, &errNum);
    }
    (void )WlzFreeObj(sft);
    (void )WlzFreeObj(str);
    (void )WlzFreeObj(dil);
    (void )WlzFreeObj(scl);
    if(idx == NULL)
    {
      (void )WlzFreeVoxelValueTb(val.vox);
    }
    /* Set tile indices. */
    if(errNum == WLZ_ERR_NONE)
    {
      int tileCnt = 0;

      errNum = WlzGreySetIncValues(idx, &tileCnt);
      tVal->numTiles = tileCnt;
    }
    /* Create index array from the index object. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	***idxArray = NULL;
      WlzIVertex3 sz,
      		  org;

      org.vtX = org.vtY = org.vtZ = 0;
      sz.vtX = idx->domain.p->lastkl - idx->domain.p->kol1 + 1;
      sz.vtY = idx->domain.p->lastln - idx->domain.p->line1 + 1;
      sz.vtZ = idx->domain.p->lastpl - idx->domain.p->plane1 + 1;
      errNum = WlzToArray3D((void ****)&idxArray, idx, sz, org, 0,
      			    WLZ_GREY_INT);
      tVal->nIdx[0] = sz.vtX;
      tVal->nIdx[1] = sz.vtY;
      tVal->nIdx[2] = sz.vtZ;
      tVal->indices = **idxArray;
      AlcFree(idxArray);
    }
    (void )WlzFreeObj(idx);
  }
  /* Allocate tiles if the given object has values. */
  if((errNum == WLZ_ERR_NONE) && (setTiles != 0))
  {
    errNum = WlzMakeTiledValuesTiles(tVal);
  }
  /* Make the object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues val;

    val.t = tVal;
    tObj = WlzMakeMain(gObj->type, gObj->domain, val, NULL, NULL, &errNum);
  }
  /* Copy the values to the object with the tiles values. */
  if((errNum == WLZ_ERR_NONE) && (setTiles != 0))
  {
    errNum = WlzCopyObjectGreyValues(tObj, gObj);
  }
  /* Clear up on error. */
  if(errNum != WLZ_ERR_NONE)
  {
    if(tObj != NULL)
    {
      WlzFreeObj(tObj);
    }
    else if(tVal != NULL)
    {
      (void )AlcFree(tVal->indices);
      (void )AlcFree(tVal->tiles.v);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}
