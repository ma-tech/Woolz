#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSepFilter_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzSepFilter.c
* \author       Bill Hill
* \date         November 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	A 2/3D implementation of seperable spatial kernel
* 		based convolution functions.
* \ingroup	WlzValuesFilters
*/

#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

#undef _OPENMP

#ifdef _OPENMP
#include <omp.h>
#endif


static WlzObject		*WlzSepFilterX(WlzObject *inObj,
				  int dim,
				  int maxThr,
				  double **iBuf,
				  double **rBuf,
				  int cBufSz,
				  double *cBuf,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzSepFilterY(WlzObject *inObj,
				  int dim,
				  int maxThr,
				  double **iBuf,
				  double **rBuf,
				  int cBufSz,
				  double *cBuf,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzSepFilterZ(WlzObject *inObj,
				  WlzIBox3 bBox,
				  int maxThr,
				  double **iBuf,
				  double **rBuf,
				  int cBufSz,
				  double *cBuf,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr);
/*!
* \return	New Gaussian filtered object with new values or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies a Gaussian filter to the given input spatial domain
* 		object. Only scalar valued image objects are valid for use
* 		with this function.
* \param	inObj			Input 2 or 3D spatial domain object
* 					which must have values.
* \param	sigma			Gaussian sigma values for each of the
* 					Cartesian directions.
* \param	order			Derivative order for each of the
* 					Cartesian directions, range 0-2,
* 					where 0 implies no derivative.
* \param	direc			Required Cartesian directions, with
* 					zero value indicating that no filtering
* 					is to be applied for a direction.
* \param	gType			Required return object grey type.
* 					Passing in WLZ_GREY_ERROR will
* 					request the given input object's grey
* 					type.
* \param	pad			Type of padding.
* \param	padVal			Padding value, only used when
* 					pad == ALG_PAD_VALUE.
* \param	dstErr			Destination error pointer may be NULL.
*/
WlzObject			*WlzGaussFilter(
				  WlzObject *inObj,
				  WlzDVertex3 sigma,
				  WlzIVertex3 order,
				  WlzIVertex3 direc,
				  WlzGreyType gType,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr)
{
  int		dim = 0;
  double	*cBuf[3] = {NULL};
  WlzIVertex3   cBufSz = {0};
  WlzObject	*rnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	maxOrder = 2;
  const double  nSigma = 3.0,
  		sigmaMin = 0.0001;

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
  else if((sigma.vtX < sigmaMin) || (sigma.vtY < sigmaMin) ||
          (order.vtX < 0) || (order.vtX > maxOrder) ||
	  (order.vtY < 0) || (order.vtY > maxOrder))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(inObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        dim = 2;
	break;
      case WLZ_3D_DOMAINOBJ:
        dim = 3;
	if((sigma.vtZ < sigmaMin) ||
	   (order.vtZ < 0) || (order.vtZ > maxOrder))
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (gType == WLZ_GREY_ERROR))
  {
    gType = WlzGreyTypeFromObj(inObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      switch(gType)
      {
        case WLZ_GREY_INT:    /* FALLTHROUGH */
        case WLZ_GREY_SHORT:  /* FALLTHROUGH */
        case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
        case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
        case WLZ_GREY_DOUBLE:
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cBufSz.vtX = (int )floor(nSigma * sigma.vtX);
    cBufSz.vtY = (int )floor(nSigma * sigma.vtY);
    if(errNum == WLZ_ERR_NONE)
    {
      if(dim == 3)
      {
	cBufSz.vtZ = (int )floor(nSigma * sigma.vtZ);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		cSz;

    cSz = (cBufSz.vtX + cBufSz.vtY + cBufSz.vtZ + 3) * 2;
    if((cBuf[0] = (double *)AlcMalloc(sizeof(double) * cSz)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    cBuf[1] = cBuf[0] + 1 + (2 * cBufSz.vtX);
    cBuf[2] = cBuf[1] + 1 + (2 * cBufSz.vtY);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idd;

    for(idd = 0; idd < dim; ++idd)
    {
      int	n,
		odr;
      double	*cb;

      cb = cBuf[idd];
      switch(idd)
      {
        case 0:
	  n = cBufSz.vtX;
	  odr = order.vtX;
	  break;
	case 1:
	  n = cBufSz.vtY;
	  odr = order.vtY;
	  break;
	case 2:
	  n = cBufSz.vtZ;
	  odr = order.vtZ;
	  break;
	default:
	  break;
      }
      /* Set convolution buffers. */
      {
	int	i;
	double	g,
		s,
		t,
		nsn,
		nrm;

	nrm = 0.0;
	nsn = nSigma / n;
	switch(odr)
	{
	  case 0:
	    for(i = 0; i <= n; ++i)
	    {
	      s = n - i;
	      t = s * nsn;
	      g = exp(-0.5 * t * t);
	      nrm += g;
	      cb[i] = g;
	    }
	    nrm = 1.0 / ((2.0 * nrm) - cb[n]);
	    for(i = 0; i <= n; ++i)
	    {
	      cb[i] *= nrm;
	    }
	    for(i = 1; i <= n; ++i)
	    {
	      cb[n + i] = cb[n - i]; 
	    }
	    break;
	  case 1:
	    for(i = 0; i <= n; ++i)
	    {
	      s = n - i;
	      t = s * nsn;
	      g = exp(-0.5 * t * t) * s;
	      nrm += g;
	      cb[i] = g;
	    }
	    nrm = 1.0 / ((2.0 * nrm) - cb[n]);
	    for(i = 0; i <= n; ++i)
	    {
	      cb[i] *= nrm;
	    }
	    for(i = 1; i <= n; ++i)
	    {
	      cb[n + i] = -cb[n - i]; 
	    }
	    break;
	  case 2:
	    for(i = 0; i <= n; ++i)
	    {
	      s = n - i;
	      t = s * nsn;
	      g = exp(-0.5 * t * t) * (s * s - nSigma * nSigma); 
	      nrm += g;
	      cb[i] = g;
	    }
	    nrm = -1.0 / ((2.0 * nrm) - cb[n]);
	    for(i = 0; i <= n; ++i)
	    {
	      cb[i] *= nrm;
	    }
	    for(i = 1; i <= n; ++i)
	    {
	      cb[n + i] = cb[n - i]; 
	    }
	    break;
	  default:
	    break;
	}
#ifdef WLZ_SEPFILTER_DEBUG
	for(i = 0; i < (2 * n) + 1; ++i)
	{
	  (void )fprintf(stderr, "%g\n", cb[i]);
	}
	(void )fprintf(stderr, "\n");
#endif /* WLZ_SEPFILTER_DEBUG */
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rnObj = WlzSepFilter(inObj, cBufSz, cBuf, direc, gType, pad, padVal,
    			 &errNum);
  }
  AlcFree(cBuf[0]);
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(rnObj);
    rnObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rnObj);
}

/*!
* \return	New filtered object with new values or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies a seperable filter to the given object using the given
* 		convolution kernels.
* \param	inObj			Input 2 or 3D spatial domain object
* 					to be filtered which must have scalar
* 					values.
* \param	cBufSz			Convolution kernel sizes (sz), each
* 					kernel buffer is sized (2 * sz) + 1
* 					with the centre indexed sz into the
* 					buffer.
* \param	cBuf			Convolution kernel buffers.
* \param	direc			Set to non-zero in directions for which
* 					the filter is to be applied.
* \param	gType			Required return object grey type.
* 					Passing in WLZ_GREY_ERROR will
* 					request the given input object's grey
* 					type.
* \param	pad			Type of padding.
* \param	padVal			Padding value, only used when
* 					pad == ALG_PAD_VALUE.
* \param	dstErr			Destination error pointer may be NULL.
*/
WlzObject			*WlzSepFilter(WlzObject *inObj,
				  WlzIVertex3 cBufSz,
				  double *cBuf[],
				  WlzIVertex3 direc,
				  WlzGreyType gType,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr)
{
  int		dim = 0,
  		vSz = 0,
  		nThr = 1;
  double	**iBuf = NULL,
  		**rBuf = NULL;
  double	*vBuf = NULL;
  WlzObject	*rnObj = NULL;
  WlzIVertex3	vBufSz = {0};
  WlzIBox3	bBox = {0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp master
    {
      nThr = omp_get_num_threads();
    }
  }
#endif
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
    switch(inObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        dim = 2;
	break;
      case WLZ_3D_DOMAINOBJ:
        dim = 3;
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (gType == WLZ_GREY_ERROR))
  {
    gType = WlzGreyTypeFromObj(inObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      switch(gType)
      {
        case WLZ_GREY_INT:    /* FALLTHROUGH */
        case WLZ_GREY_SHORT:  /* FALLTHROUGH */
        case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
        case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
        case WLZ_GREY_DOUBLE:
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox = WlzBoundingBox3I(inObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      vBufSz.vtX = bBox.xMax - bBox.xMin + 1;
      vBufSz.vtY = bBox.yMax - bBox.yMin + 1;
      if(dim == 3)
      {
	vBufSz.vtZ = bBox.zMax - bBox.zMin + 1;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vSz = ALG_MAX3(vBufSz.vtX, vBufSz.vtY, vBufSz.vtZ);
    if(((iBuf = (double **)AlcMalloc(sizeof(double *) * 2 * nThr)) == NULL) ||
       ((vBuf = (double *)AlcMalloc(sizeof(double) * 2 * nThr * vSz)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int	idt;

      rBuf = iBuf + nThr;
      for(idt = 0; idt < nThr; ++idt)
      {
        iBuf[idt] = vBuf + (idt * vSz);
	rBuf[idt] = vBuf + ((nThr + idt) * vSz);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Convolve the object values. */
    if(direc.vtX)
    {
      rnObj = WlzSepFilterX(inObj, dim, nThr,
			    iBuf, rBuf, cBufSz.vtX, cBuf[0],
			    pad, padVal, &errNum);
    }
    if((errNum == WLZ_ERR_NONE) && direc.vtY)
    {
      WlzObject *tObj;

      tObj = WlzSepFilterY((rnObj)? rnObj: inObj, dim, nThr,
                            iBuf, rBuf, cBufSz.vtY, cBuf[1],
			    pad, padVal, &errNum);
      (void )WlzFreeObj(rnObj);
      rnObj = tObj;
    }
    if((errNum == WLZ_ERR_NONE) && (dim == 3) && direc.vtZ)
    {
      WlzObject *tObj;

      tObj = WlzSepFilterZ((rnObj)? rnObj: inObj, bBox, nThr,
                            iBuf, rBuf, cBufSz.vtZ, cBuf[2],
			    pad, padVal, &errNum);
      (void )WlzFreeObj(rnObj);
      rnObj = tObj;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (rnObj != NULL) && (gType != WLZ_GREY_DOUBLE))
  {
    WlzObject *tObj;

    /* Convert object values to the required grey type. */
    tObj = WlzConvertPix((rnObj)? rnObj: inObj, gType, &errNum);
    (void )WlzFreeObj(rnObj);
    rnObj = tObj;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rnObj);
    rnObj = NULL;
  }
  AlcFree(iBuf);
  AlcFree(vBuf);
  return(rnObj);
}

/*!
* \return	New filtered object with new values or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies a seperable filter along the x axis (columns)
* 		to the given object using the given convolution kernel.
* \param	inObj			Input 2 or 3D spatial domain object
* 					to be filtered which must have scalar
* 					values.
* \param	dim			Object's dimension.
* \param	maxThr			Maximum number of threads to use.
* \param	iBuf			Working buffers large enough for any
* 				        line in the image with one for each
* 				        thread.
* \param	rBuf			Buffers as for iBuf. 			
* \param	cBufSz			Convolution kernel size.
* \param	cBuf			Convolution kernel buffer.
* \param	pad			Type of padding.
* \param	padVal			Padding value.
* \param	dstErr			Destination error pointer may be NULL.
*/
static WlzObject		*WlzSepFilterX(WlzObject *inObj,
				  int dim,
				  int maxThr,
				  double **iBuf,
				  double **rBuf,
				  int cBufSz,
				  double *cBuf,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr)
{
  int		idp,
		poff,
  		nPln;
  WlzObjectType	rGTT;
  WlzDomain	*domains;
  WlzValues	*iVal,
  		*rVal;
  WlzPixelV	zV;
  WlzObject	*rnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  		
  zV.v.dbv = 0.0;
  zV.type = WLZ_GREY_DOUBLE;
  rGTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, NULL);
  if(dim == 3)
  {
    WlzPlaneDomain *pDom;
    WlzVoxelValues *vVal;

    pDom = inObj->domain.p;
    vVal = inObj->values.vox;
    nPln = pDom->lastpl - pDom->plane1 + 1;
    domains = pDom->domains;
    iVal = vVal->values;
    poff = pDom->plane1 - vVal->plane1;
    rnObj = WlzNewObjectValues(inObj, rGTT, zV, 1, zV, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      rVal = rnObj->values.vox->values;
    }
  }
  else
  {
    nPln = 1;
    poff = 0;
    domains = &(inObj->domain);
    iVal = &(inObj->values);
    rnObj = WlzNewObjectValues(inObj, rGTT, zV, 1, zV, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      rVal = &(rnObj->values);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(maxThr)
#endif
    for(idp = 0; idp < nPln; ++idp)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	if(domains[idp].core != NULL)
	{
	  int		thrId = 0;
	  WlzObject 	*iObj = NULL,
			*rObj = NULL;
	  WlzIntervalWSpace iIWSp,
			    rIWSp;
	  WlzGreyWSpace	iGWSp,
			rGWSp;
          WlzErrorNum	errNum2 = WLZ_ERR_NONE;

#ifdef _OPENMP
          thrId = omp_get_thread_num();
#endif
	  iObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[idp], iVal[idp + poff],
			     NULL, NULL, &errNum2);
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    rObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[idp], rVal[idp],
			       NULL, NULL, &errNum2);
	  }
	  if((errNum2 == WLZ_ERR_NONE) &&
	     ((errNum2 = WlzInitGreyScan(iObj, &iIWSp,
	                                 &iGWSp)) == WLZ_ERR_NONE) &&
	     ((errNum2 = WlzInitGreyScan(rObj, &rIWSp,
	                                 &rGWSp)) == WLZ_ERR_NONE))
	  {
	    while(((errNum2 = WlzNextGreyInterval(&iIWSp)) == WLZ_ERR_NONE) &&
		  ((errNum2 = WlzNextGreyInterval(&rIWSp)) == WLZ_ERR_NONE))
	    {
	      size_t	len;
	      WlzGreyP	iBufGP;

	      
	      iBufGP.dbp = iBuf[thrId];
	      len = iIWSp.rgtpos - iIWSp.lftpos + 1;
	      WlzValueCopyGreyToGrey(iBufGP, 0, WLZ_GREY_DOUBLE,
				     iGWSp.u_grintptr, 0, iGWSp.pixeltype,
				     len);
	      AlgConvolveD(len, rBuf[thrId], cBufSz * 2 + 1, cBuf,
                           len, iBuf[thrId], pad, padVal);
	      WlzValueCopyDoubleToDouble(rGWSp.u_grintptr.dbp,
				         rBuf[thrId], len);
	    }
	    (void )WlzEndGreyScan(&iIWSp, &iGWSp);
	    (void )WlzEndGreyScan(&rIWSp, &rGWSp);
	    if(errNum2 == WLZ_ERR_EOO)
	    {
	      errNum2 = WLZ_ERR_NONE;
	    }
	  }
	  (void )WlzFreeObj(iObj);
	  (void )WlzFreeObj(rObj);
	  if(errNum2 != WLZ_ERR_NONE)
	  {
#ifdef _OPENMP
#pragma omp critical
	    {
#endif
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = errNum2;
	      }
#ifdef _OPENMP
	    }
#endif
	  }
	}
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rnObj);
    rnObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rnObj);
}

/*!
* \return	New filtered object with new values or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies a seperable filter along the y axis (lines)
* 		to the given object using the given convolution kernel.
* \param	inObj			Input 2 or 3D spatial domain object
* 					to be filtered which must have scalar
* 					values.
* \param	dim			Object's dimension.
* \param	maxThr			Maximum number of threads to use.
* \param	iBuf			Working buffers large enough for any
* 				        column in the image with one for each
* 				        thread.
* \param	rBuf			Buffers as for iBuf. 			
* \param	cBufSz			Convolution kernel size.
* \param	cBuf			Convolution kernel buffer.
* \param	pad			Type of padding.
* \param	padVal			Padding value.
* \param	dstErr			Destination error pointer may be NULL.
*/
static WlzObject		*WlzSepFilterY(WlzObject *inObj,
				  int dim,
				  int maxThr,
				  double **iBuf,
				  double **rBuf,
				  int cBufSz,
				  double *cBuf,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr)
{
  int		idp,
		poff,
  		nPln;
  WlzObjectType	rGTT;
  WlzDomain	*domains;
  WlzValues	*iVal,
  		*rVal;
  WlzPixelV	zV;
  WlzObject	*rnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  		
  zV.v.dbv = 0.0;
  zV.type = WLZ_GREY_DOUBLE;
  rGTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, NULL);
  if(dim == 3)
  {
    WlzPlaneDomain *pDom;
    WlzVoxelValues *vVal;

    pDom = inObj->domain.p;
    vVal = inObj->values.vox;
    nPln = pDom->lastpl - pDom->plane1 + 1;
    domains = pDom->domains;
    iVal = vVal->values;
    poff = pDom->plane1 - vVal->plane1;
    rnObj = WlzNewObjectValues(inObj, rGTT, zV, 1, zV, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      rVal = rnObj->values.vox->values;
    }
  }
  else
  {
    nPln = 1;
    poff = 0;
    domains = &(inObj->domain);
    iVal = &(inObj->values);
    rnObj = WlzNewObjectValues(inObj, rGTT, zV, 1, zV, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      rVal = &(rnObj->values);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(maxThr)
#endif
    for(idp = 0; idp < nPln; ++idp)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	if(domains[idp].core != NULL)
	{
	  int		thrId = 0;
	  WlzIBox2	bBox;
	  WlzObject	*iObj = NULL,
			*rObj = NULL;
	  WlzGreyValueWSpace *iVWSp = NULL,
			     *rVWSp = NULL;
          WlzErrorNum	errNum2 = WLZ_ERR_NONE;

#ifdef _OPENMP
          thrId = omp_get_thread_num();
#endif
	  iObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[idp], iVal[idp + poff],
			     NULL, NULL, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    rObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[idp], rVal[idp],
			       NULL, NULL, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    bBox = WlzBoundingBox2I(iObj, &errNum);
	  }
	  if((errNum == WLZ_ERR_NONE) &&
	     ((iVWSp = WlzGreyValueMakeWSp(iObj, &errNum)) != NULL) &&
	     ((rVWSp = WlzGreyValueMakeWSp(rObj, &errNum)) != NULL))
	  {
	    int		idx;

	    for(idx = bBox.xMin; idx <= bBox.xMax; ++idx)
	    {
	      int		y0,
			  y1,
			  idy;
	      WlzUByte	in = 0;

	      for(idy = bBox.yMin; idy <= bBox.yMax + 1; ++idy)
	      {
		WlzGreyValueGet(iVWSp, 0, idy, idx);
		if(iVWSp->bkdFlag == 0)
		{
		 double	*ibp;

		 if(in == 0)
		 {
		   y0 = idy;
		   in = 1;
		 }
		 y1 = idy;
		 ibp = iBuf[thrId] + (y1  - y0);
		 switch(iVWSp->gType)
		 {
		   case WLZ_GREY_INT:
		     *ibp = iVWSp->gVal[0].inv;
		     break;
		   case WLZ_GREY_SHORT:
		     *ibp = iVWSp->gVal[0].shv;
		     break;
		   case WLZ_GREY_UBYTE:
		     *ibp = iVWSp->gVal[0].ubv;
		     break;
		   case WLZ_GREY_FLOAT:
		     *ibp = iVWSp->gVal[0].flv;
		     break;
		   case WLZ_GREY_DOUBLE:
		     *ibp = iVWSp->gVal[0].dbv;
		     break;
		   default:
		     break;
		 }
		}
		else if(in)
		{
		  int 	idr,
			len;

		  len = y1 - y0 + 1;
		  AlgConvolveD(len, rBuf[thrId], cBufSz * 2 + 1, cBuf,
		               len, iBuf[thrId], pad, padVal);
		  for(idr = y0; idr <= y1; ++idr)
		  {
		    WlzGreyValueGet(rVWSp, 0, idr, idx);
		    *(rVWSp->gPtr[0].dbp) = rBuf[thrId][idr - y0];
		  }
		}
	      }
	    }
	  }
	  WlzGreyValueFreeWSp(iVWSp);
	  WlzGreyValueFreeWSp(rVWSp);
	  (void )WlzFreeObj(iObj);
	  (void )WlzFreeObj(rObj);
	  if(errNum2 != WLZ_ERR_NONE)
	  {
#ifdef _OPENMP
#pragma omp critical
	    {
#endif
	      if(errNum == WLZ_ERR_NONE)
	      {
	        errNum = errNum2;
	      }
#ifdef _OPENMP
	    }
#endif
	  }
	}
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rnObj);
    rnObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rnObj);
}

/*!
* \return	New filtered object with new values or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Applies a seperable filter along the z axis (planes)
* 		to the given object using the given convolution kernel.
* \param	inObj			Input 3D spatial domain object to be
* 					filtered which must have scalar values.
* \param	bBox			Object's bounding box.
* \param	maxThr			Maximum number of threads to use.
* \param	iBuf			Working buffers large enough for any
* 				        column in the image with one for each
* 				        thread.
* \param	rBuf			Buffers as for iBuf. 			
* \param	cBufSz			Convolution kernel size.
* \param	cBuf			Convolution kernel buffer.
* \param	pad			Type of padding.
* \param	padVal			Padding value.
* \param	dstErr			Destination error pointer may be NULL.
*/
static WlzObject		*WlzSepFilterZ(WlzObject *inObj,
				  WlzIBox3 bBox,
				  int maxThr,
				  double **iBuf,
				  double **rBuf,
				  int cBufSz,
				  double *cBuf,
				  AlgPadType pad,
				  double padVal,
				  WlzErrorNum *dstErr)
{
  WlzObjectType	rGTT;
  WlzPixelV	zV;
  WlzObject	*rnObj = NULL;
  WlzGreyValueWSpace **iVWSp = NULL,
	             **rVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  		
  zV.v.dbv = 0.0;
  zV.type = WLZ_GREY_DOUBLE;
  rGTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, NULL);
  rnObj = WlzNewObjectValues(inObj, rGTT, zV, 1, zV, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    if(((iVWSp = (WlzGreyValueWSpace **)
                 AlcCalloc(maxThr, sizeof(WlzGreyValueWSpace *))) == NULL) ||
       ((rVWSp = (WlzGreyValueWSpace **)
                 AlcCalloc(maxThr, sizeof(WlzGreyValueWSpace *))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idt;

    for(idt = 0; (idt < maxThr) && (errNum == WLZ_ERR_NONE); ++idt)
    {
      iVWSp[idt] = WlzGreyValueMakeWSp(inObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	rVWSp[idt] = WlzGreyValueMakeWSp(rnObj, &errNum);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idy;

#ifdef _OPENMP
#pragma omp parallel for num_threads(maxThr)
#endif
    for(idy = bBox.yMin; idy <= bBox.yMax; ++idy)
    {
      int	idx,
      		thrId = 0;
      WlzGreyValueWSpace *iVWSpT,
			 *rVWSpT;

#ifdef _OPENMP
      thrId = omp_get_thread_num();
#endif
      iVWSpT = iVWSp[thrId];
      rVWSpT = rVWSp[thrId];
      for(idx = bBox.xMin; idx <= bBox.xMax; ++idx)
      {
	int 	z0,
	      	z1,
	      	idz;
	WlzUByte in = 0;

	for(idz = bBox.zMin; idz <= bBox.zMax + 1; ++idz)
	{
	  WlzGreyValueGet(iVWSpT, idz, idy, idx);
	  if(iVWSpT->bkdFlag == 0)
	  {
	   double	*ibp;

	   if(in == 0)
	   {
	     z0 = idz;
	     in = 1;
	   }
	   z1 = idz;
	   ibp = iBuf[thrId] + (z1  - z0);
	   switch(iVWSpT->gType)
	   {
	     case WLZ_GREY_INT:
	       *ibp = iVWSpT->gVal[0].inv;
	       break;
	     case WLZ_GREY_SHORT:
	       *ibp = iVWSpT->gVal[0].shv;
	       break;
	     case WLZ_GREY_UBYTE:
	       *ibp = iVWSpT->gVal[0].ubv;
	       break;
	     case WLZ_GREY_FLOAT:
	       *ibp = iVWSpT->gVal[0].flv;
	       break;
	     case WLZ_GREY_DOUBLE:
	       *ibp = iVWSpT->gVal[0].dbv;
	       break;
	     default:
	       break;
	   }
	  }
	  else if(in)
	  {
	    int	idr,
		len;

	    len = z1 - z0 + 1;
	    AlgConvolveD(len, rBuf[thrId], cBufSz * 2 + 1, cBuf,
			 len, iBuf[thrId], pad, padVal);
	    for(idr = z0; idr <= z1; ++idr)
	    {
	      WlzGreyValueGet(rVWSpT, idr, idy, idx);
	      *(rVWSpT->gPtr[0].dbp) = rBuf[thrId][idr - z0];
	    }
	  }
	}
      }
    }
  }
  if(iVWSp)
  {
    int		idt;

    for(idt = 0; idt < maxThr; ++idt)
    {
      WlzGreyValueFreeWSp(iVWSp[idt]);
    }
    AlcFree(iVWSp);
  }
  if(rVWSp)
  {
    int		idt;

    for(idt = 0; idt < maxThr; ++idt)
    {
      WlzGreyValueFreeWSp(rVWSp[idt]);
    }
    AlcFree(rVWSp);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rnObj);
    rnObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rnObj);
}
