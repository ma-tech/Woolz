#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTiledValues_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzTiledValues.c
* \author       Bill Hill
* \date         March 2010
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
* \brief	Functions construct and convert tiled objects.
* \ingroup	WlzAllocation
*/

#include <stdlib.h>
#include <limits.h>
#include <Wlz.h>

#ifdef HAVE_MMAP
#define WLZ_USE_MMAP
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#endif

/* The use of Hilbert order for tiles need some testing before
real use. */
/* #define WLZ_TILES_USE_HILBERT */

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
* \param	dim			Dimension.
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
      AlcFree(tVal->nIdx);
      if(tVal->original_table.core)
      {
	if(tVal->original_table.core->type != tVal->type)
	{
	  errNum = WLZ_ERR_VALUES_TYPE;
	}
	else
	{
          errNum = WlzFreeTiledValues(tVal->original_table.t);
	}
      }
      else
      {
	AlcFree(tVal->indices);
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
      }
      AlcFree(tVal);
    }
  }
  return(errNum);
}

/*!
* \return	New tiled values or NULL on error.
* \ingroup	WlzAllocation
* \brief	Creates a new tiled values table which shares both
*  		the grey value tiles and indexing with the given
*  		tiled values. This is done by using the original_table
*  		field.
* \param	gVal			The given tiled values.
* \param	bgdV			Background value for the new values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzTiledValues	*WlzNewTiledValues(WlzTiledValues *gVal, WlzPixelV bgdV,
				   WlzErrorNum *dstErr)
{
  WlzTiledValues *rVal = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gVal == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    rVal = WlzMakeTiledValues(gVal->dim, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      int   i;
      WlzValues val;

      rVal->type       = gVal->type;
      rVal->dim        = gVal->dim;
      rVal->kol1       = gVal->kol1;
      rVal->lastkl     = gVal->lastkl;
      rVal->line1      = gVal->line1;
      rVal->lastln     = gVal->lastln;
      rVal->plane1     = gVal->plane1;
      rVal->lastpl     = gVal->lastpl;
      rVal->tileSz     = gVal->tileSz;
      rVal->tileWidth  = gVal->tileWidth;
      rVal->numTiles   = gVal->numTiles;
      rVal->fd         = gVal->fd;
      rVal->tileOffset = gVal->tileOffset;
      rVal->tiles.v    = gVal->tiles.v;
      rVal->indices    = gVal->indices;
      rVal->bckgrnd    = bgdV;
      for(i = 0; i < gVal->dim; ++i)
      {
	rVal->nIdx[i] = gVal->nIdx[i];
      }
      val.t = gVal;
      rVal->original_table = WlzAssignValues(val, NULL);
    }
  }
  if((errNum != WLZ_ERR_NONE) && rVal)
  {
    (void )WlzFreeTiledValues(rVal);
    rVal = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rVal);
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
* \brief	Allocates tiles for the tiles values and sets background value
* 		within the tiles.
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
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValueSetGrey(tVal->tiles, 0, tVal->bckgrnd.v, tVal->bckgrnd.type,
                    tVal->numTiles * tVal->tileSz);
  }
  return(errNum);
}

/*!
* \return	Flags specifying how the value table may be used.
* \ingroup	WlzValuesUtils
* \brief	Determines how the tiled values of a tiles value table
* 		may be used. If memory mapping is in use then the mode
* 		in which the file (containing the memory mapped values)
* 		was opened determines the modes in which the tiled values
* 		can be accessed. Inappropriate access (eg attempting to
* 		change a grey value in a memory mapped tiled values
* 		table that was opened for read only) will generate a
* 		memory fault. The returned value is a bit mask in which
* 		WLZ_IOFLAGS_READ will be set iff grey values can be read
* 		and WLZ_IOFLAGS_WRITE will be set iff the grey values can
* 		be written to (ie modified).
* \param	tv			The given tiled values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzTiledValuesMode(WlzTiledValues *tv, WlzErrorNum *dstErr)
{
  int		flags = WLZ_IOFLAGS_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tv == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(WlzGreyTableIsTiled(tv->type) != WLZ_GREY_TAB_TILED)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
#ifdef WLZ_USE_MMAP
    if((tv->fd < 0) && (tv->tiles.v != NULL))
    {
      flags = WLZ_IOFLAGS_READ | WLZ_IOFLAGS_WRITE;
    }
    else
    {
      int	mod;

      mod = fcntl(tv->fd, F_GETFL);
      if(mod < 0)
      {
	errNum = WLZ_ERR_FILE_OPEN;
      }
      else if((mod & O_ACCMODE) == O_RDONLY)
      {
        flags = WLZ_IOFLAGS_READ;
      }
      else if((mod & O_ACCMODE) == O_WRONLY)
      {
        flags = WLZ_IOFLAGS_WRITE;
      }
      else if((mod & O_ACCMODE) == O_RDWR)
      {
        flags = WLZ_IOFLAGS_READ | WLZ_IOFLAGS_WRITE;
      }
    }
#else
    flags = WLZ_IOFLAGS_READ | WLZ_IOFLAGS_WRITE;
#endif
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(flags);
}


/*!
* \return	New tiled values buffer.
* \ingroup	WlzAllocation
* \brief	Makes a new tiled values buffer for the given tiled values
* 		table.
* \param	tVal			Given tiled values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzTiledValueBuffer *WlzMakeTiledValueBuffer(WlzTiledValues *tVal,
				WlzErrorNum *dstErr)
{
  int		mode = WLZ_IOFLAGS_NONE,
  		lnWidth = 0;
  size_t	gSz = 0;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzTiledValueBuffer *tBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tVal == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(WlzGreyTableIsTiled(tVal->type) != WLZ_GREY_TAB_TILED)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((lnWidth = tVal->lastkl - tVal->kol1 + 1) < 1)
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    mode = WlzTiledValuesMode(tVal, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((gType = WlzGreyTableTypeToGreyType(tVal->type,
						 NULL)) == WLZ_GREY_ERROR) ||
	    ((gSz = WlzGreySize(gType)) <= 0))
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
    else if(((tBuf = (WlzTiledValueBuffer *)
		     AlcCalloc(1, sizeof(WlzTiledValueBuffer))) == NULL) ||
	    ((tBuf->lnbuf.v = AlcMalloc(gSz * lnWidth)) == NULL))
    {
      AlcFree(tBuf);
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      tBuf->gtype = gType;
      tBuf->mode = mode;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tBuf);
}

/*!
* \ingroup	WlzAllocation
* \brief	Frees the given tiled values buffer including it's line
* 		buffer.
* \param	tBuf			Given tiled values, may be NULL.
*/
void		WlzFreeTiledValueBuffer(WlzTiledValueBuffer *tBuf)
{
  if(tBuf)
  {
    AlcFree(tBuf->lnbuf.v);
    AlcFree(tBuf);                                              
  }
}

/*!
* \ingroup	WlzValuesUtils
* \brief	Flushes the values in the tiled values buffer to the tiled
* 		values table.
* \param	tvb			Given tiled values buffer.
* \param	tv			Given tiled values table.
*/
void		WlzTiledValueBufferFlush(WlzTiledValueBuffer *tvb,
				WlzTiledValues *tv)
{
  if(tvb->valid && ((tvb->mode & WLZ_IOFLAGS_WRITE) != 0))
  {
    int		kol,
  		ti,
  		to;

    kol = tvb->kl[0];
    while(kol <= tvb->kl[1])
    {
      int	i,
      		io,
		itc = 0,
		rmn;
      size_t	ii;

      ti = kol / tv->tileWidth;
      to = kol % tv->tileWidth;
      io = tvb->lo + to;
      ii = *(tv->indices + tvb->li + ti);
      if(ii > 0)
      {
	rmn = tvb->kl[1] - kol + 1;
	itc = tv->tileWidth - to;
	if(itc > rmn)
	{
	  itc = rmn;
	}
	switch(tvb->gtype)
	{
	  case WLZ_GREY_LONG:
	    {
	      WlzLong *bp,
		      *tp;

	      tp = tv->tiles.lnp + (ii * tv->tileSz) + io;
	      bp = tvb->lnbuf.lnp + kol;
	      for(i = 0; i < itc; ++i)
	      {
		*tp++ = *bp++;
	      }
	    }
	    break;
	  case WLZ_GREY_INT:
	    {
	      int	 *bp,
		   *tp;

	      tp = tv->tiles.inp + (ii * tv->tileSz) + io;
	      bp = tvb->lnbuf.inp + kol;
	      for(i = 0; i < itc; ++i)
	      {
		*tp++ = *bp++;
	      }
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    {
	      short *bp,
		    *tp;

	      tp = tv->tiles.shp + (ii * tv->tileSz) + io;
	      bp = tvb->lnbuf.shp + kol;
	      for(i = 0; i < itc; ++i)
	      {
		*tp++ = *bp++;
	      }
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    {
	      WlzUByte *bp,
		       *tp;

	      tp = tv->tiles.ubp + (ii * tv->tileSz) + io;
	      bp = tvb->lnbuf.ubp + kol;
	      for(i = 0; i < itc; ++i)
	      {
		*tp++ = *bp++;
	      }
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    {
	      float *bp,
		    *tp;

	      tp = tv->tiles.flp + (ii * tv->tileSz) + io;
	      bp = tvb->lnbuf.flp + kol;
	      for(i = 0; i < itc; ++i)
	      {
		*tp++ = *bp++;
	      }
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    {
	      double *bp,
		     *tp;

	      tp = tv->tiles.dbp + (ii * tv->tileSz) + io;
	      bp = tvb->lnbuf.dbp + kol;
	      for(i = 0; i < itc; ++i)
	      {
		*tp++ = *bp++;
	      }
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    {
	      WlzUInt *bp,
		      *tp;

	      tp = tv->tiles.rgbp + (ii * tv->tileSz) + io;
	      bp = tvb->lnbuf.rgbp + kol;
	      for(i = 0; i < itc; ++i)
	      {
		*tp++ = *bp++;
	      }
	    }
	    break;
	  default:
	    break;
	}
      }
      kol += itc;
    }
  }
}

/*!
* \ingroup	WlzValuesUtils
* \brief	Fills the tiled values buffer using values from the tiled
* 		values table.
* \param	tvb			Given tiled values buffer.
* \param	tv			Given tiled values table.
*/
void		WlzTiledValueBufferFill(WlzTiledValueBuffer *tvb,
				WlzTiledValues *tv)
{
  int		kol;
  int		ti[3],
  		to[3];

  kol = tvb->kl[0];
  ti[1] = tvb->ln / tv->tileWidth;
  to[1] = tvb->ln % tv->tileWidth;
  if(tv->dim == 2)
  {
    tvb->li = ti[1] * tv->nIdx[0];
    tvb->lo = to[1] * tv->tileWidth;
  }
  else /* tv->dim == 3 */
  {
    ti[2] = tvb->pl / tv->tileWidth;
    to[2] = tvb->pl % tv->tileWidth;
    tvb->li = (ti[2] * tv->nIdx[1] + ti[1]) * tv->nIdx[0];
    tvb->lo = ((to[2] * tv->tileWidth) + to[1]) * tv->tileWidth;
  }
  if((tvb->mode & WLZ_IOFLAGS_READ) != 0)
  {
    while(kol <= tvb->kl[1])
    {
      int	i,
		ii,
      		io,
		itc,
		rmn;

      ti[0] = kol / tv->tileWidth;
      to[0] = kol % tv->tileWidth;
      rmn = tvb->kl[1] - kol + 1;
      itc = tv->tileWidth - to[0];
      if(itc > rmn)
      {
	itc = rmn;
      }
      io = tvb->lo + to[0];
      ii = *(tv->indices + tvb->li + ti[0]);
      switch(tvb->gtype)
      {
	case WLZ_GREY_LONG:
	  {
	    WlzLong *bp;

	    bp = tvb->lnbuf.lnp + kol;
	    if(ii < 0)
	    {
	      for(i = 0; i < itc; ++i)
	      {
	        *bp++ = tv->bckgrnd.v.lnv;
	      }
	    }
	    else
	    {
	      WlzLong *tp;

	      tp = tv->tiles.lnp + (ii * tv->tileSz) + io;
	      for(i = 0; i < itc; ++i)
	      {
		*bp++ = *tp++;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_INT:
	  {
	    int *bp;

	    bp = tvb->lnbuf.inp + kol;
	    if(ii < 0)
	    {
	      for(i = 0; i < itc; ++i)
	      {
	        *bp++ = tv->bckgrnd.v.inv;
	      }
	    }
	    else
	    {
	      int *tp;

	      tp = tv->tiles.inp + (ii * tv->tileSz) + io;
	      for(i = 0; i < itc; ++i)
	      {
		*bp++ = *tp++;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_SHORT:
	  {
	    short *bp;

	    bp = tvb->lnbuf.shp + kol;
	    if(ii < 0)
	    {
	      for(i = 0; i < itc; ++i)
	      {
	        *bp++ = tv->bckgrnd.v.shv;
	      }
	    }
	    else
	    {
	      short *tp;

	      tp = tv->tiles.shp + (ii * tv->tileSz) + io;
	      for(i = 0; i < itc; ++i)
	      {
		*bp++ = *tp++;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  {
	    WlzUByte *bp;

	    bp = tvb->lnbuf.ubp + kol;
	    if(ii < 0)
	    {
	      for(i = 0; i < itc; ++i)
	      {
	        *bp++ = tv->bckgrnd.v.ubv;
	      }
	    }
	    else
	    {
	      WlzUByte *tp;

	      tp = tv->tiles.ubp + (ii * tv->tileSz) + io;
	      for(i = 0; i < itc; ++i)
	      {
		*bp++ = *tp++;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  {
	    float *bp;

	    bp = tvb->lnbuf.flp + kol;
	    if(ii < 0)
	    {
	      for(i = 0; i < itc; ++i)
	      {
	        *bp++ = tv->bckgrnd.v.flv;
	      }
	    }
	    else
	    {
	      float *tp;

	      tp = tv->tiles.flp + (ii * tv->tileSz) + io;
	      for(i = 0; i < itc; ++i)
	      {
		*bp++ = *tp++;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  {
	    double *bp;

	    bp = tvb->lnbuf.dbp + kol;
	    if(ii < 0)
	    {
	      for(i = 0; i < itc; ++i)
	      {
	        *bp++ = tv->bckgrnd.v.dbv;
	      }
	    }
	    else
	    {
	      double *tp;

	      tp = tv->tiles.dbp + (ii * tv->tileSz) + io;
	      for(i = 0; i < itc; ++i)
	      {
		*bp++ = *tp++;
	      }
	    }
	  }
	  break;
	case WLZ_GREY_RGBA:
	  {
	    WlzUInt *bp;

	    bp = tvb->lnbuf.rgbp + kol;
	    if(ii < 0)
	    {
	      for(i = 0; i < itc; ++i)
	      {
	        *bp++ = tv->bckgrnd.v.rgbv;
	      }
	    }
	    else
	    {
	      WlzUInt *tp;

	      tp = tv->tiles.rgbp + (ii * tv->tileSz) + io;
	      for(i = 0; i < itc; ++i)
	      {
		*bp++ = *tp++;
	      }
	    }
	  }
	  break;
	default:
	  break;
      }
      kol += itc;
    }
  }
  tvb->valid = 1;
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
    tVal->type = WlzGreyValueTableType(0, WLZ_GREY_TAB_TILED, gType, NULL);
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
      tabType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT,
                                      NULL);
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

#ifdef WLZ_TILES_USE_HILBERT
      errNum = WlzGreySetHilbertRankValues(idx, &tileCnt);
#else /* WLZ_TILES_USE_HILBERT */
      errNum = WlzGreySetIncValues(idx, &tileCnt);
#endif /* WLZ_TILES_USE_HILBERT */
      tVal->numTiles = tileCnt;
    }
    /* Create index array from the index object. */
    if(errNum == WLZ_ERR_NONE)
    {
      unsigned int **idxArray = NULL;
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
  if(errNum == WLZ_ERR_NONE)
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
  if((errNum == WLZ_ERR_NONE) && (setTiles != 0) &&
     (gObj->values.core != NULL))
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
    tVal->type = WlzGreyValueTableType(0, WLZ_GREY_TAB_TILED, gType, NULL);
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
      tabType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT,
                                      NULL);
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

#ifdef WLZ_TILES_USE_HILBERT
      errNum = WlzGreySetHilbertRankValues(idx, &tileCnt);
#else /* WLZ_TILES_USE_HILBERT */
      errNum = WlzGreySetIncValues(idx, &tileCnt);
#endif /* WLZ_TILES_USE_HILBERT */
      tVal->numTiles = tileCnt;
    }
    /* Create index array from the index object. */
    if(errNum == WLZ_ERR_NONE)
    {
      unsigned int ***idxArray = NULL;
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
  if(errNum == WLZ_ERR_NONE)
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
  if((errNum == WLZ_ERR_NONE) && (setTiles != 0) &&
     (gObj->values.core != NULL))
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
