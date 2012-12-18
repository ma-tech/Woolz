#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFourier_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzFourier.c
* \author       Bill Hill
* \date         November 2012
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
* \brief	Functions for computing Fourier transforms and their
* 		inverse.
* \ingroup	WlzFourier
*/

#include <Wlz.h>

static WlzObject 		*WlzFourierTransformObj2D(
				  WlzObject *iObj,
				  int fwd,
				  WlzErrorNum *dstErr);
static WlzCompoundArray 	*WlzFourierTransformObjCpx2D(
				  WlzCompoundArray *iObj,
				  int fwd,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzFourierTransformObj3D(
				  WlzObject *iObj,
				  int fwd,
				  WlzErrorNum *dstErr);
static WlzCompoundArray 	*WlzFourierTransformObjCpx3D(
				  WlzCompoundArray *iObj,
				  int fwd,
				  WlzErrorNum *dstErr);
/*!
* \return	New Woolz domain object or NULL on error.
* \ingroup	WlzFourier
* \brief	Computes either the forward or inverse Fourier transform
* 		of a domain object with real (ie not complex) values.
* 		When computing a transform the object will be padded
* 		to an integer power of two size.
* 		The object's values can have any single valued type
* 		(and therefore RGBA is not acceptable). For forward
* 		transforms the objects frequently have their grey
* 		values modified so that the mean value is zero; this
* 		is not done in this function.
*		The transformed data are scaled, see AlgFour1D(),
*		AlgFourReal1D(), AlgFour2D() AlgFourReal2D(),
*		AlgFour3D() and AlgFourReal3D().
* \param	iObj			Input Woolz object.
* \param	fwd			Non-zero if the transform is a
* 					forward transform, zero if it is to
* 					be an inverse transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzFourierTransformObj(WlzObject *iObj, int fwd,
					WlzErrorNum *dstErr)
{
  WlzObject	*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(iObj->type)
    {
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
        {
	  WlzCompoundArray *cIObj;

	  cIObj = (WlzCompoundArray *)iObj;
	  if(cIObj->n != 2)
	  {
	    errNum = WLZ_ERR_OBJECT_TYPE;
	  }
	  else if((cIObj->o[0] == NULL) || (cIObj->o[1] == NULL))
	  {
	    errNum = WLZ_ERR_OBJECT_NULL;
	  }
	  else if(cIObj->o[0]->type != cIObj->o[1]->type)
	  {
	    errNum = WLZ_ERR_OBJECT_TYPE;
	  }
	  else
	  {
	    switch(cIObj->o[0]->type)
	    {
	      case WLZ_2D_DOMAINOBJ:
	        oObj = (WlzObject *)
		       WlzFourierTransformObjCpx2D(cIObj, fwd, &errNum);
	       break;
	      case WLZ_3D_DOMAINOBJ:
	        oObj = (WlzObject *)
		       WlzFourierTransformObjCpx3D(cIObj, fwd, &errNum);
		break;
	      default:
		errNum = WLZ_ERR_OBJECT_TYPE;
		break;
	    }
	  }
	}
	break;
      case WLZ_2D_DOMAINOBJ:
        oObj = WlzFourierTransformObj2D(iObj, fwd, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
        oObj = WlzFourierTransformObj3D(iObj, fwd, &errNum);
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
  return(oObj);
}

/*!
* \return	New Woolz domain object or NULL on error.
* \ingroup	WlzFourier
* \brief	Computes either the forward or inverse 2D Fourier transform
* 		of a 2D domain object with real (ie not complex) values.
* 		See WlzFourierTransformObj(). The given object is known to
* 		be non-null and have type WLZ_2D_DOMAINOBJ.
* \param	iObj			Input Woolz object.
* \param	fwd			Non-zero for forward transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzFourierTransformObj2D(WlzObject *iObj, int fwd,
					   WlzErrorNum *dstErr)
{
  WlzIVertex2 	org,
  		iSz,
                oSz;
  WlzIBox2  	bBox;
  void		**array = NULL;
  WlzObject	*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(iObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    bBox = WlzBoundingBox2I(iObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    unsigned int t0;

    iSz.vtX = bBox.xMax - bBox.xMin + 1;
    iSz.vtY = bBox.yMax - bBox.yMin + 1;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtX);
    oSz.vtX = t0;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtY);
    oSz.vtY = t0;
    org.vtX = bBox.xMin - (oSz.vtX - iSz.vtX) / 2;
    org.vtY = bBox.yMin - (oSz.vtY - iSz.vtY) / 2;
    errNum = WlzToArray2D(&array, iObj, oSz, org, 0, WLZ_GREY_DOUBLE);
    if(errNum == WLZ_ERR_NONE)
    {
      AlgError	algErr = ALG_ERR_NONE;

      if(fwd)
      {
	{
	  algErr = AlgFourReal2D((double **)array, 1, oSz.vtX, oSz.vtY);
	}
      }
      else
      {
	algErr = AlgFourRealInv2D((double **)array, 1, oSz.vtX, oSz.vtY);
      }
      errNum = WlzErrorFromAlg(algErr);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WLZ_VTX_2_ZERO(org);
    oObj = WlzFromArray2D(array, oSz, org,
			  WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			  0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(oObj)
    {
      (void )WlzFreeObj(oObj);
      oObj = NULL;
    }
    else
    {
      (void )Alc2Free(array);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oObj);
}

/*!
* \return	New Woolz domain object or NULL on error.
* \ingroup	WlzFourier
* \brief	Computes either the forward or inverse 3D Fourier transform
* 		of a 3D domain object with real (ie not complex) values.
* 		See WlzFourierTransformObj(). The given object is known to
* 		be non-null and have type WLZ_3D_DOMAINOBJ.
* \param	iObj			Input Woolz object.
* \param	fwd			Non-zero for forward transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzFourierTransformObj3D(WlzObject *iObj, int fwd,
					   WlzErrorNum *dstErr)
{
  WlzIVertex3 	org,
  		iSz,
                oSz;
  WlzIBox3  	bBox;
  void		***array = NULL;
  WlzObject	*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(iObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    bBox = WlzBoundingBox3I(iObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    unsigned int t0;

    iSz.vtX = bBox.xMax - bBox.xMin + 1;
    iSz.vtY = bBox.yMax - bBox.yMin + 1;
    iSz.vtZ = bBox.zMax - bBox.zMin + 1;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtX);
    oSz.vtX = t0;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtY);
    oSz.vtY = t0;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtZ);
    oSz.vtZ = t0;
    org.vtX = bBox.xMin - (oSz.vtX - iSz.vtX) / 2;
    org.vtY = bBox.yMin - (oSz.vtY - iSz.vtY) / 2;
    org.vtZ = bBox.zMin - (oSz.vtZ - iSz.vtZ) / 2;
    errNum = WlzToArray3D(&array, iObj, oSz, org, 0, WLZ_GREY_DOUBLE);
    if(errNum == WLZ_ERR_NONE)
    {
      AlgError	algErr = ALG_ERR_NONE;

      if(fwd)
      {
	{
	  algErr = AlgFourReal3D((double ***)array, 1,
	                         oSz.vtX, oSz.vtY, oSz.vtZ);
	}
      }
      else
      {
	algErr = AlgFourRealInv3D((double ***)array, 1,
	                          oSz.vtX, oSz.vtY, oSz.vtZ);
      }
      errNum = WlzErrorFromAlg(algErr);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WLZ_VTX_3_ZERO(org);
    oObj = WlzFromArray3D(array, oSz, org,
			  WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			  0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(oObj)
    {
      (void )WlzFreeObj(oObj);
      oObj = NULL;
    }
    else
    {
      (void )Alc3Free(array);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oObj);
}

/*!
* \return	New Woolz domain object or NULL on error.
* \ingroup	WlzFourier
* \brief	Computes either the forward or inverse 2D Fourier transform
* 		of a 2D domain object with real (ie not complex) values.
* 		See WlzFourierTransformObj(). The given object is known to
* 		be non-null and have type WLZ_2D_DOMAINOBJ.
* \param	iObj			Input Woolz object.
* \param	fwd			Non-zero for forward transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray *WlzFourierTransformObjCpx2D(WlzCompoundArray *iObj,
					int fwd, WlzErrorNum *dstErr)
{
  WlzIVertex2 	iSz,
                oSz,
		org;
  WlzIBox2  	bBox[2];
  void		**real = NULL,
  		**imag = NULL;
  WlzCompoundArray *oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((iObj->o[0]->domain.core == NULL) ||
     (iObj->o[1]->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((iObj->o[0]->values.core == NULL) ||
          (iObj->o[1]->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    bBox[0] = WlzBoundingBox2I(iObj->o[0], &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox[1] = WlzBoundingBox2I(iObj->o[1], &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((bBox[0].xMin != bBox[1].xMin) ||
       (bBox[0].yMin != bBox[1].yMin) ||
       (bBox[0].xMax != bBox[1].xMax) ||
       (bBox[0].yMax != bBox[1].yMax))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {

    unsigned int t0;

    iSz.vtX = bBox[0].xMax - bBox[0].xMin + 1;
    iSz.vtY = bBox[0].yMax - bBox[0].yMin + 1;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtX);
    oSz.vtX = t0;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtY);
    oSz.vtY = t0;
    org.vtX = bBox[0].xMin - (oSz.vtX - iSz.vtX) / 2;
    org.vtY = bBox[0].yMin - (oSz.vtY - iSz.vtY) / 2;
    errNum = WlzToArray2D(&real, iObj->o[0], oSz, org, 0, WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzToArray2D(&imag, iObj->o[1], oSz, org, 0, WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    AlgError	algErr = ALG_ERR_NONE;

    if(fwd)
    {
      algErr = AlgFour2D((double **)real, (double **)imag, 1,
			 oSz.vtX, oSz.vtY);
    }
    else
    {
      algErr = AlgFourInv2D((double **)real, (double **)imag, 1,
                            oSz.vtX, oSz.vtY);
    }
    errNum = WlzErrorFromAlg(algErr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject *o[2];

    o[0] = o[1] = NULL;
    WLZ_VTX_2_ZERO(org);
    o[0] = WlzFromArray2D(real, oSz, org,
			  WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			  0.0, 1.0, 0, 1, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      o[1] = WlzFromArray2D(imag, oSz, org,
			    WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			    0.0, 1.0, 0, 1, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      oObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 2, o,
                                  WLZ_2D_DOMAINOBJ, &errNum);
    }
    if(oObj == NULL)
    {
      (void )WlzFreeObj(o[0]);
      (void )WlzFreeObj(o[1]);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(oObj)
    {
      (void )WlzFreeObj((WlzObject *)oObj);
      oObj = NULL;
    }
    else
    {
      (void )Alc2Free(real);
      (void )Alc2Free(imag);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oObj);
}

/*!
* \return	New Woolz domain object or NULL on error.
* \ingroup	WlzFourier
* \brief	Computes either the forward or inverse 3D Fourier transform
* 		of a 3D domain object with complex values represented as a
* 		compound object with the first object having real values
* 		and the second imaginary.
* 		See WlzFourierTransformObj(). The given object is known to
* 		be non-null, have non-null components and the components
* 		are known to have  type WLZ_3D_DOMAINOBJ.
* \param	iObj			Input Woolz object.
* \param	fwd			Non-zero for forward transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray *WlzFourierTransformObjCpx3D(WlzCompoundArray *iObj,
                                 	int fwd, WlzErrorNum *dstErr)
{
  WlzIVertex3 	iSz,
                oSz,
		org;
  WlzIBox3  	bBox[2];
  void		***real = NULL,
  		***imag = NULL;
  WlzCompoundArray *oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((iObj->o[0]->domain.core == NULL) ||
     (iObj->o[1]->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((iObj->o[0]->values.core == NULL) ||
          (iObj->o[1]->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    bBox[0] = WlzBoundingBox3I(iObj->o[0], &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox[1] = WlzBoundingBox3I(iObj->o[1], &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((bBox[0].xMin != bBox[1].xMin) ||
       (bBox[0].yMin != bBox[1].yMin) ||
       (bBox[0].xMax != bBox[1].xMax) ||
       (bBox[0].yMax != bBox[1].yMax))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {

    unsigned int t0;

    iSz.vtX = bBox[0].xMax - bBox[0].xMin + 1;
    iSz.vtY = bBox[0].yMax - bBox[0].yMin + 1;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtX);
    oSz.vtX = t0;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtY);
    oSz.vtY = t0;
    (void )AlgBitNextPowerOfTwo(&t0, iSz.vtZ);
    oSz.vtZ = t0;
    org.vtX = bBox[0].xMin - (oSz.vtX - iSz.vtX) / 2;
    org.vtY = bBox[0].yMin - (oSz.vtY - iSz.vtY) / 2;
    org.vtZ = bBox[0].zMin - (oSz.vtZ - iSz.vtZ) / 2;
    errNum = WlzToArray3D(&real, iObj->o[0], oSz, org, 0, WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzToArray3D(&imag, iObj->o[1], oSz, org, 0, WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    AlgError	algErr = ALG_ERR_NONE;

    if(fwd)
    {
      algErr = AlgFour3D((double ***)real, (double ***)imag, 1,
			 oSz.vtX, oSz.vtY, oSz.vtZ);
    }
    else
    {
      algErr = AlgFourInv3D((double ***)real, (double ***)imag, 1,
                            oSz.vtX, oSz.vtY, oSz.vtZ);
    }
    errNum = WlzErrorFromAlg(algErr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject *o[2];

    o[0] = o[1] = NULL;
    WLZ_VTX_3_ZERO(org);
    o[0] = WlzFromArray3D(real, oSz, org,
			  WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			  0.0, 1.0, 0, 1, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      o[1] = WlzFromArray3D(imag, oSz, org,
			    WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			    0.0, 1.0, 0, 1, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      oObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 2, o,
                                  WLZ_3D_DOMAINOBJ, &errNum);
    }
    if(oObj == NULL)
    {
      (void )WlzFreeObj(o[0]);
      (void )WlzFreeObj(o[1]);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(oObj)
    {
      (void )WlzFreeObj((WlzObject *)oObj);
      oObj = NULL;
    }
    else
    {
      (void )Alc3Free(real);
      (void )Alc3Free(imag);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oObj);
}
