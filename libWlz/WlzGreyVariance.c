#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyVariance_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGreyVariance.c
* \author       Bill Hill
* \date         February 2008
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
* \brief	Applies variance filters of domain objects with values.
* \ingroup	WlzValuesFilters
*/
#include <stdlib.h>
#include <Wlz.h>

/*!
* \struct       _WlzGreyVarianceWSp
* \ingroup      WlzValuesFilters
* \brief        Used to pass size and scale through WlzSeqPar() to
*		WlzGreyVarianceSeqParFn().
*/
typedef struct	_WlzGreyVarianceWSp
{
  int		size;
  WlzGreyType	gType;
  double	scale;
} WlzGreyVarianceWSp;

static int			WlzGreyVarianceSeqParFn(
				  WlzSeqParWSpace *spWSpace,
				  void *spData);

/*!
* \return       Woolz object or NULL on error.
* \ingroup      WlzValuesFilters
* \brief        Uses a small kernel to compute the local variance
*		of the given grey value object.
*               Only objects with WLZ_EMPTY_OBJ and WLZ_2D_DOMAINOBJ
*               types are valid. WLZ_2D_DOMAINOBJ objects must have
*               a valid domain and integral (ie WLZ_GREY_INT,
*               WLZ_GREY_RGBA, WLZ_GREY_SHORT or WLZ_GREY_UBYTE)
*		values.
* \param        inObj                   Given object.
* \param        newObjFlag              If zero the filter is applied in
*					place.
* \param	kernelSz		Kernel size, must be 3, 5 or 7.
* \param	scale			Variance data scale factor.
* \param        dstErr			Destination error pointer, may be NULL.
*/
WlzObject       *WlzGreyVariance(WlzObject *inObj, int newObjFlag,
                                 int kernelSz, double scale,
				 WlzErrorNum *dstErr)
{
  WlzObject     *outObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzPixelV     bkgVal;
  WlzGreyType   gValType;
  WlzGreyVarianceWSp vWSp;

  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(inObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(kernelSz)
    {
      case 3: /* FALLTHROUGH */
      case 5: /* FALLTHROUGH */
      case 7:
        break;
      default:
	errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(inObj->type)
    {
      case WLZ_EMPTY_OBJ:
        outObj = (newObjFlag)? WlzMakeEmpty(&errNum): inObj;
        break;
      case WLZ_2D_DOMAINOBJ:
        bkgVal = WlzGetBackground(inObj, &errNum);
        if(errNum == WLZ_ERR_NONE)
        {
          errNum = WlzValueConvertPixel(&bkgVal, bkgVal, WLZ_GREY_INT);
        }
        if(errNum == WLZ_ERR_NONE)
	{
	  gValType = WlzGreyTableTypeToGreyType(inObj->values.core->type,
	             				&errNum);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  switch(gValType)
	  {
	    case WLZ_GREY_INT:    /* FALLTHROUGH */
	    case WLZ_GREY_SHORT:  /* FALLTHROUGH */
	    case WLZ_GREY_RGBA:   /* FALLTHROUGH */
	    case WLZ_GREY_UBYTE:
	      break;
	    case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
	    case WLZ_GREY_DOUBLE: /* FALLTHROUGH */
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  vWSp.scale = scale;
	  vWSp.size = kernelSz;
	  vWSp.gType = gValType;
	  outObj = WlzSeqPar(inObj, newObjFlag, 0, WLZ_RASTERDIR_ILIC,
			     kernelSz, bkgVal.v.inv, &vWSp,
			     WlzGreyVarianceSeqParFn, &errNum);
	}
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
  return(outObj);
}

/*!
* \return	Variance value.
* \ingroup	WlzValueFilters
* \brief	Computes the variance of the region about a pixel
* 		as a function for WlzSeqPar().
* \param	spWSpace		Work space data structure from
* 					WlzSeqPar().
* \param	spData			Data passed on by WlzSeqPar(),
* 					which is used to pass kernel size
* 					and variance scaling factor.
*/
static int	WlzGreyVarianceSeqParFn(WlzSeqParWSpace *spWSpace,
                                        void *spData)
{
  int		idX,
		idY,
		sz,
		varI;
  double	sum = 0.0,
  		sumSq = 0.0,
  		area,
		varD;
  int		*data;
  WlzGreyVarianceWSp *vWSp;

  vWSp = (WlzGreyVarianceWSp *)spData;
  sz = (vWSp->size - 1) / 2;
  for(idY = -sz; idY <= sz; ++idY)
  {
    data = spWSpace->adrptr[idY] - sz;
    for(idX = -sz; idX <= sz; ++idX)

    {
      sum += *data;
      sumSq += *data * *data;
      ++data;
    }
  }
  area = vWSp->size * vWSp->size;
  varD = vWSp->scale * (sumSq - ((sum * sum) / area)) / (area - 1.0);
  varI = WLZ_NINT(varD);
  return(varI);
}
