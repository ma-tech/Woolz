#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzScalarFn_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzScalarFn.c
* \author       Bill Hill
* \date         August 2006
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
* \brief	Code to apply scalar functions to scalar image values in
*		Woolz objects.
* \ingroup	WlzArithmetic
*/

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Wlz.h>

static WlzObject 		*WlzScalarFn2D(
				  WlzObject *sObj,
				  WlzFnType fn,
			          WlzErrorNum *dstErr);
static WlzObject 		*WlzScalarFn3D(
				  WlzObject *sObj,
				  WlzFnType fn,
			          WlzErrorNum *dstErr);
static WlzGreyType 		WlzScalarFnPromoteGType(
				  WlzFnType fn,
				  WlzGreyType gType,
				  WlzErrorNum *dstErr);
static WlzPixelV 		WlzScalarFnPixel(
				  WlzPixelV sPV,
				  WlzFnType fn,
				  WlzErrorNum *dstErr);
static void			WlzScalarFnItv(
				  WlzGreyP gValP,
				  WlzGreyType gType,
				  int len,
			          WlzFnType fn);
static void			WlzScalarFnItvMod(
				  WlzGreyP gValP,
				  WlzGreyType gType,
				  int len);
static void			WlzScalarFnItvExp(
				  WlzGreyP gValP,
				  WlzGreyType gType,
				  int len);
static void			WlzScalarFnItvLog(
				  WlzGreyP gValP,
				  WlzGreyType gType,
				  int len);
static void			WlzScalarFnItvSqrt(
				  WlzGreyP gValP,
				  WlzGreyType gType,
				  int len);
static void			WlzScalarFnItvInvSqrt(
				  WlzGreyP gValP,
				  WlzGreyType gType,
				  int len);
static void			WlzScalarFnItvSqr(
				  WlzGreyP gValP,
				  WlzGreyType gType,
				  int len);

/*!
* \return	New domain object.
* \ingroup	WlzArithmetic
* \brief	Computes a new object which shares the domain of the given
*		object, but which has grey values that are the result of
*		applying the given function to the grey values of the
*		given object.
* \param	sObj			Given source domainn object with values.
* \param	fn			Scalar function to be applied.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzScalarFn(WlzObject *sObj, WlzFnType fn,
			     WlzErrorNum *dstErr)
{
  WlzObject	*dObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(sObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(sObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(sObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(sObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        dObj = WlzScalarFn2D(sObj, fn, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
        dObj = WlzScalarFn3D(sObj, fn, &errNum);
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
  return(dObj);
}

/*!
* \return	New 2D domain object.
* \ingroup	WlzArithmetic
* \brief	Computes a new 2D object which shares the domain of the given
*		object, but which has grey values that are the result of
*		applying the given function to the grey values of the
*		given object.
* \param	sObj			Given source domain object with values.
* \param	fn			Scalar function to be applied.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzScalarFn2D(WlzObject *sObj, WlzFnType fn,
			        WlzErrorNum *dstErr)
{
  WlzGreyType	sGType,
  		dGType;
  WlzPixelV	sBgd,
  		dBgd;
  WlzObjectType	dVType;
  WlzObject     *dObj = NULL;
  WlzValues	dVal;
  WlzIntervalWSpace sIWSp,
  		dIWSp;
  WlzGreyWSpace	sGWSp,
  		dGWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  sGType = WlzGreyTypeFromObj(sObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dGType = WlzScalarFnPromoteGType(fn, sGType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dVType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, dGType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sBgd = WlzGetBackground(sObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dBgd = WlzScalarFnPixel(sBgd, fn, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dVal.v = WlzNewValueTb(sObj, dVType, dBgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
    		       sObj->domain, dVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(sObj, &sIWSp, &sGWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(dObj, &dIWSp, &dGWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&sIWSp)) == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&dIWSp)) == WLZ_ERR_NONE))
    {
      int	len;
      WlzGreyP	dGP,
      		sGP;

      len = sIWSp.rgtpos - sIWSp.lftpos + 1;
      dGP = dGWSp.u_grintptr;
      sGP = sGWSp.u_grintptr;
      WlzValueCopyGreyToGrey(dGP, 0, dGType, sGP, 0, sGType, len);
      WlzScalarFnItv(dGP, dGType, len, fn);
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dObj);
}

/*!
* \return	New 3D domain object.
* \ingroup	WlzArithmetic
* \brief	Computes a new 3D object which shares the domain of the given
*		object, but which has grey values that are the result of
*		applying the given function to the grey values of the
*		given object.
* \param	sObj			Given source domain object with values.
* \param	fn			Scalar function to be applied.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzScalarFn3D(WlzObject *sObj, WlzFnType fn,
			        WlzErrorNum *dstErr)
{
  WlzPixelV	sBgd,
  		dBgd;
  WlzValues	dVal;
  WlzObject     *dObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(errNum == WLZ_ERR_NONE)
  {
    sBgd = WlzGetBackground(sObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dBgd = WlzScalarFnPixel(sBgd, fn, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
    				   sObj->domain.p->plane1,
				   sObj->domain.p->lastpl,
				   dBgd, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx,
    		plane1,
		lastpl;
    WlzDomain	*doms;
    WlzValues	*dVals,
    		*sVals;
    WlzErrorNum	errNum2 = WLZ_ERR_NONE;

    doms = sObj->domain.p->domains;
    sVals = sObj->values.vox->values;
    dVals = dVal.vox->values;
    plane1 = sObj->domain.p->plane1;
    lastpl = sObj->domain.p->lastpl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idx = plane1; idx <= lastpl; ++idx)
    {
      int	off;
      WlzObject	*sObj2D,
      		*dObj2D = NULL;

      off = idx - sObj->domain.p->plane1;
      sObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, doms[off], sVals[off],
			   NULL, NULL, &errNum2);
      if(errNum2 == WLZ_ERR_NONE)
      {
        dObj2D = WlzScalarFn2D(sObj2D, fn, &errNum2);
      }
      WlzFreeObj(sObj2D);
      if(errNum2 == WLZ_ERR_NONE)
      {
        dVals[off] = WlzAssignValues(dObj2D->values, NULL);
      }
      WlzFreeObj(dObj2D);
      if(errNum2 != WLZ_ERR_NONE)
      {
#ifdef _OPENMP
#pragma omp critical (WlzScalarFn3D)
        {
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = errNum2;
	  }
	}
#else
        errNum = errNum2;
#endif
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, sObj->domain, dVal,
    		       NULL, NULL, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dObj);
}

/*!
* \return	Promoted grey type.
* \ingroup	WlzArithmetic
* \brief	Promotes the given grey value type depending on the function
*		type.
* \param	fn			Scalar function type.
* \param	gType			Grey type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzGreyType WlzScalarFnPromoteGType(WlzFnType fn, WlzGreyType gType,
				WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_LONG:  /* FALLTHROUGH */
    case WLZ_GREY_INT:   /* FALLTHROUGH */
    case WLZ_GREY_SHORT: /* FALLTHROUGH */
    case WLZ_GREY_UBYTE: /* FALLTHROUGH */
    case WLZ_GREY_FLOAT: /* FALLTHROUGH */
    case WLZ_GREY_DOUBLE:
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  switch(fn)
  {
    case WLZ_FN_SCALAR_MOD:
      break;
    case WLZ_FN_SCALAR_EXP:     /* FALLTHROUGH */
    case WLZ_FN_SCALAR_LOG:     /* FALLTHROUGH */
    case WLZ_FN_SCALAR_SQRT:    /* FALLTHROUGH */
    case WLZ_FN_SCALAR_INVSQRT:
    case WLZ_FN_SCALAR_SQR:
      gType = WLZ_GREY_DOUBLE;
      break;
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    gType = WLZ_GREY_ERROR;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(gType);
}

/*!
* \return	Resulting pixel value.
* \ingroup	WlzArithmetic
* \brief	Promotes the grey type of the given pixel value and then
* 		applies the given scalar function to it's value.
* \param	sPV			Single pixel value.
* \param	fn			Scalar function type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzPixelV WlzScalarFnPixel(WlzPixelV sPV, WlzFnType fn,
				  WlzErrorNum *dstErr)
{
  WlzGreyP	dGP;
  WlzPixelV	dPV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dPV.type = WlzScalarFnPromoteGType(fn, sPV.type, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzValueConvertPixel(&dPV, sPV, dPV.type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(dPV.type)
    {
      case WLZ_GREY_LONG:
	dGP.lnp = &(dPV.v.lnv);
	break;
      case WLZ_GREY_INT:
	dGP.inp = &(dPV.v.inv);
	break;
      case WLZ_GREY_SHORT:
	dGP.shp = &(dPV.v.shv);
	break;
      case WLZ_GREY_UBYTE:
	dGP.ubp = &(dPV.v.ubv);
	break;
      case WLZ_GREY_FLOAT:
	dGP.flp = &(dPV.v.flv);
	break;
      case WLZ_GREY_DOUBLE:
	dGP.dbp = &(dPV.v.dbv);
	break;
      case WLZ_GREY_RGBA:
	dGP.rgbp = &(dPV.v.rgbv);
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzScalarFnItv(dGP, dPV.type, 1, fn);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dPV);
}

/*!
* \ingroup	WlzArithmetic
* \brief	Applies the given scalar function to the interval buffer,
*		overwriting the original values.
* \param	gValP			Buffer pointer.
* \param	gType			Grey type.
* \param	len			Interval buffer length.
* \param	fn			Function to apply.
*/
static void	WlzScalarFnItv(WlzGreyP gValP, WlzGreyType gType, int len,
			       WlzFnType fn)
{
  switch(fn)
  {
    case WLZ_FN_SCALAR_MOD:
      WlzScalarFnItvMod(gValP, gType, len);
      break;
    case WLZ_FN_SCALAR_EXP:
      WlzScalarFnItvExp(gValP, gType, len);
      break;
    case WLZ_FN_SCALAR_LOG:
      WlzScalarFnItvLog(gValP, gType, len);
      break;
    case WLZ_FN_SCALAR_SQRT:
      WlzScalarFnItvSqrt(gValP, gType, len);
      break;
    case WLZ_FN_SCALAR_INVSQRT:
      WlzScalarFnItvInvSqrt(gValP, gType, len);
      break;
    case WLZ_FN_SCALAR_SQR:
      WlzScalarFnItvSqr(gValP, gType, len);
      break;
    default:
      break;
  }
}

/*!
* \ingroup	WlzArithmetic
* \brief	Applies a modulus function (eg abs() or fabs()) to data in
*		the buffer.
* \param	gValP			Buffer pointer.
* \param	gType			Grey type.
* \param	len			Interval buffer length.
*/
static void	WlzScalarFnItvMod(WlzGreyP gValP, WlzGreyType gType, int len)
{
  int		idx;

  switch(gType)
  {
    case WLZ_GREY_LONG:
      for(idx = 0; idx < len; ++idx)
      {
        *(gValP.lnp) = labs(*(gValP.lnp));
	(gValP.lnp)++;
      }
      break;
    case WLZ_GREY_INT:
      for(idx = 0; idx < len; ++idx)
      {
        *(gValP.inp) = abs(*(gValP.inp));
	(gValP.inp)++;
      }
      break;
    case WLZ_GREY_SHORT:
      for(idx = 0; idx < len; ++idx)
      {
        *(gValP.shp) = (short )abs(*(gValP.shp));
	(gValP.shp)++;
      }
      break;
    case WLZ_GREY_UBYTE:
      /* Can't be -ve. */
      break;
    case WLZ_GREY_FLOAT:
      for(idx = 0; idx < len; ++idx)
      {
        *(gValP.flp) = fabsf(*(gValP.flp));
        (gValP.flp)++;
      }
      break;
    case WLZ_GREY_DOUBLE:
      for(idx = 0; idx < len; ++idx)
      {
        *(gValP.dbp) = fabs(*(gValP.dbp));
        (gValP.dbp)++;
      }
      break;
    case WLZ_GREY_RGBA:
      /* Can't be -ve. */
      break;
    default:
      break;
  }
}

/*!
* \ingroup	WlzArithmetic
* \brief	Applies an exponential function(ie exp()) to data in
*		the buffer.
* \param	gValP			Buffer pointer.
* \param	gType			Grey type.
* \param	len			Interval buffer length.
*/
static void	WlzScalarFnItvExp(WlzGreyP gValP, WlzGreyType gType, int len)
{
  int		idx;

  switch(gType)
  {
    case WLZ_GREY_DOUBLE:
      for(idx = 0; idx < len; ++idx)
      {
        *(gValP.dbp) = exp(*(gValP.dbp));
	++(gValP.dbp);
      }
    default:
      /* Should only have WLZ_GREY_DOUBLE! */
      break;
  }
}

/*!
* \ingroup	WlzArithmetic
* \brief	Applies a logarithm function(ie log()) to data in
*		the buffer. Values <= to DBL_EPSION result in
*		values of -DBL_MAX.
* \param	gValP			Buffer pointer.
* \param	gType			Grey type.
* \param	len			Interval buffer length.
*/
static void	WlzScalarFnItvLog(WlzGreyP gValP, WlzGreyType gType, int len)
{
  int		idx;
  double	tD0;

  switch(gType)
  {
    case WLZ_GREY_DOUBLE:
      for(idx = 0; idx < len; ++idx)
      {
	tD0 = fabs(*(gValP.dbp));
        *(gValP.dbp)++ = (tD0 > DBL_EPSILON)? log(tD0): -(DBL_MAX);
      }
    default:
      /* Should only have WLZ_GREY_DOUBLE! */
      break;
  }
}

/*!
* \ingroup	WlzArithmetic
* \brief	Applies a square root function(ie sqrt(fabs())) to data in
*		the buffer.
* \param	gValP			Buffer pointer.
* \param	gType			Grey type.
* \param	len			Interval buffer length.
*/
static void	WlzScalarFnItvSqrt(WlzGreyP gValP, WlzGreyType gType, int len)
{
  int		idx;
  double	tD0;

  switch(gType)
  {
    case WLZ_GREY_DOUBLE:
      for(idx = 0; idx < len; ++idx)
      {
	tD0 = fabs(*(gValP.dbp));
        *(gValP.dbp)++ = sqrt(tD0);
      }
    default:
      /* Should only have WLZ_GREY_DOUBLE! */
      break;
  }
}

/*!
* \ingroup	WlzArithmetic
* \brief	Applies a inverse square root function(ie 1/sqrt(fabs())) to
* 		data in the buffer. Values <= to DBL_EPSION result in values
* 		of DBL_MAX.
* \param	gValP			Buffer pointer.
* \param	gType			Grey type.
* \param	len			Interval buffer length.
*/
static void	WlzScalarFnItvInvSqrt(WlzGreyP gValP, WlzGreyType gType, int len)
{
  int		idx;
  double	tD0;

  switch(gType)
  {
    case WLZ_GREY_DOUBLE:
      for(idx = 0; idx < len; ++idx)
      {
	tD0 = sqrt(fabs(*(gValP.dbp)));
        *(gValP.dbp)++ = (tD0 > DBL_EPSILON)? 1.0 / tD0: DBL_MAX;
      }
    default:
      /* Should only have WLZ_GREY_DOUBLE! */
      break;
  }
}

/*!
* \ingroup	WlzArithmetic
* \brief	Applies a square function(ie x^2) to data in the buffer.
* \param	gValP			Buffer pointer.
* \param	gType			Grey type.
* \param	len			Interval buffer length.
*/
static void	WlzScalarFnItvSqr(WlzGreyP gValP, WlzGreyType gType, int len)
{
  int		idx;
  double	tD0;

  switch(gType)
  {
    case WLZ_GREY_DOUBLE:
      for(idx = 0; idx < len; ++idx)
      {
	tD0 = *(gValP.dbp) * *(gValP.dbp);
        *(gValP.dbp)++ = tD0;
      }
    default:
      /* Should only have WLZ_GREY_DOUBLE! */
      break;
  }
}
