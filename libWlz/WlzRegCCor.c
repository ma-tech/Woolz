#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRegCCor_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzRegCCor.c
* \author       Bill Hill
* \date         January 2005
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
* \brief	Functions to register two objects using frequency
* 		domain cross correlation.
* \ingroup	WlzRegistration
* \todo         -
* \bug          None known.
*/

#include <float.h>
#include <limits.h>
#include <Wlz.h>

/* #define WLZ_REGCCOR_TEST */
/* #define WLZ_REGCCOR_DEBUG */

static WlzAffineTransform 	*WlzRegCCorObjs2D(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  WlzDVertex2 maxTran,
				  double maxRot,
				  int maxItr,
				  int *dstConv,
				  double *dstCCor,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzRegCCorObjs2D1(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzDVertex2 rotCentre,
				  WlzTransformType trType,
				  WlzDVertex2 maxTran,
				  double maxRot,
				  int maxItr,
				  int *dstConv,
				  double *dstCCor,
				  WlzErrorNum *dstErr);
static WlzDVertex2		WlzRegCCorObjs2DTran(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzDVertex2 maxTran,
				  double *dstCCor,
				  WlzErrorNum *dstErr);
static double			WlzRegCCorObjs2DRot(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzIVertex2 rotCentre,
				  WlzAffineTransform *initTr,
				  double maxRot,
				  WlzErrorNum *dstErr);

/*!
* \return	Affine transform which brings the two objects into register.
* \ingroup	WlzRegistration
* \brief	Registers the two given objects using a frequency domain
*               cross correlation.  An affine transform is computed,
*               which when applied to the source object takes it
*               into register with the target object.
* \param	tObj			The target object.
* \param	sObj			The source object to be registered
*					with target object.
* \param	initTr			Initial affine transform to be
*					applied to the source object prior to
*					registration. May be NULL.
* \param	trType			Required transform type.
* \param	maxTran			Maximum translation.
* \param	maxRot			Maximum rotation.
* \param	maxItr			Maximum number of iterations, if
*					\f$\leq\f$ 0 then infinite iterations
*					are allowed. 
* \param	dstConv			Destination ptr for the convergence
* 					flag (non zero on convergence), may be
* 					NULL.
* \param	dstCCor			Destination ptr for the cross
* 					correlation value, may be NULL.
* \param	dstErr			Destination error pointer,
*                                       may be NULL.
*/
WlzAffineTransform *WlzRegCCorObjs(WlzObject *tObj, WlzObject *sObj,
				   WlzAffineTransform *initTr,
				   WlzTransformType trType,
				   WlzDVertex2 maxTran, double maxRot,
				   int maxItr, int *dstConv, double *dstCCor,
				   WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzAffineTransform *regTr = NULL;

  if((tObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((tObj->domain.core == NULL) || (sObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((tObj->values.core == NULL) || (sObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(tObj->type != sObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else
  {
    switch(tObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        regTr = WlzRegCCorObjs2D(tObj, sObj, initTr, trType,
				 maxTran, maxRot, maxItr,
				 dstConv, dstCCor, &errNum);
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
  return(regTr);
}

/*!
* \return	Affine transform which brings the two objects into register.
* \ingroup	WlzRegistration
* \brief	Registers the two given 2D domain objects using a
*               frequency domain cross correlation.  An affine transform
*               is computed, which when applied to the source object
*               takes it into register with the target object.
*               A resolution pyramid is built from the given objects
*               and used to register the objects, progressing from
*               a low resolution towards the full resolution objects.
* \param	tObj			The target object. Must have
*                                       been assigned.
* \param	sObj			The source object to be
*                                       registered with target object.
* \param	initTr			Initial affine transform
*                                       to be applied to the source
*                                       object prior to registration.
*                                       Only translations in x and y
*                                       and rotation about the z axis
*                                       are used. May be NULL which is
*                                       equivalent to an identity transform.
* \param	trType			Required transform type.
* \param	maxTran			Maximum translation.
* \param	maxRot			Maximum rotation.
* \param	maxItr			Maximum number of iterations,
*                                       if \f$leq\f$ 0 then infinite iterations
*                                       are allowed.
* \param	dstConv			Destination ptr for the
*                                       convergence flag (non zero
*                                       on convergence), may be NULL.
* \param	dstCCor			Destination ptr for the cross
*                                       correlation value, may be NULL.
* \param	dstErr			Destination error pointer,
*                                       may be NULL.
*/
static WlzAffineTransform *WlzRegCCorObjs2D(WlzObject *tObj, WlzObject *sObj,
					    WlzAffineTransform *initTr,
					    WlzTransformType trType,
					    WlzDVertex2 maxTran, double maxRot,
					    int maxItr,
					    int *dstConv, double *dstCCor,
					    WlzErrorNum *dstErr)
{
  int		tI1,
		samIdx,
  		nSam,
		conv;
  double	cCor,
  		rot0,
		rot1,
		sMaxRot;
  WlzPixelV	gV[4];
  WlzIVertex2	tIV0,
  		tIV1;
  WlzIVertex3	samFacV;
  WlzDVertex2	tran0,
  		tran1,
		rotCentre,
		samRotCentre,
		sMaxTran;
  WlzIBox2	sBox,
  		tBox;
  int		*samFac = NULL;
  WlzObject	**sTObj = NULL,
  		**sSObj = NULL;
  WlzAffineTransform *samRegTr0 = NULL,
  		*samRegTr1 = NULL,
		*regTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzAffineTransformPrim trPrim;
  WlzPixelV	zeroBgd;
  const int	samFacStep = 4,
  		maxSam = 16,
  		minSamSz = 100;

  zeroBgd.type = WLZ_GREY_INT;
  zeroBgd.v.inv = 0;
  gV[0].type = gV[1].type = gV[2].type = gV[3].type = WLZ_GREY_DOUBLE;
  /* Compute the number of x4 subsampling operations to use. */
  sBox = WlzBoundingBox2I(sObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    tBox = WlzBoundingBox2I(tObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tIV0.vtX = sBox.xMax - sBox.xMin + 1;
    tIV0.vtY = sBox.yMax - sBox.yMin + 1;
    tIV1.vtX = tBox.xMax - tBox.xMin + 1;
    tIV1.vtY = tBox.yMax - tBox.yMin + 1;
    tIV0.vtX = WLZ_MIN(tIV0.vtX, tIV1.vtX);
    tIV0.vtY = WLZ_MIN(tIV0.vtY, tIV1.vtY);
    nSam = 1;
    tI1 = WLZ_MIN(tIV0.vtX, tIV0.vtY);
    while((nSam < maxSam) && (tI1 > minSamSz))
    {
      ++nSam;
      tI1 /= samFacStep;
    }
  }
  /* Allocate space for subsampled objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(((samFac = (int *)AlcMalloc(nSam * sizeof(int))) == NULL) ||
       ((sTObj = (WlzObject **)AlcCalloc(nSam,
    				         sizeof(WlzObject *))) == NULL) ||
       ((sSObj = (WlzObject **)AlcCalloc(nSam,
       					 sizeof(WlzObject *))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Find centre of mass of the target object. */
  if(errNum == WLZ_ERR_NONE)
  {
    rotCentre = WlzCentreOfMass2D(tObj, 0, NULL, &errNum);
  }
  /* Compute subsampled objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    samIdx = 0;
    *samFac = 1;
    samFacV.vtX = samFacV.vtY = samFacStep;
    *(sTObj + 0) = WlzAssignObject(tObj, NULL);
    *(sSObj + 0) = WlzAssignObject(sObj, NULL);
    while((errNum == WLZ_ERR_NONE) && (++samIdx < nSam))
    {
      *(samFac + samIdx) = *(samFac + samIdx - 1) * samFacStep;
      *(sTObj + samIdx) = WlzAssignObject(
      			  WlzSampleObj(*(sTObj + samIdx - 1), samFacV,
      				       WLZ_SAMPLEFN_GAUSS, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
        *(sSObj + samIdx) = WlzAssignObject(
			    WlzSampleObj(*(sSObj + samIdx - 1), samFacV,
			       	         WLZ_SAMPLEFN_GAUSS, &errNum), NULL);
      }
    }
  }
  /* Make sure the background value is zero. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(samIdx = 0; samIdx < nSam; ++samIdx)
    {
      (void )WlzSetBackground(*(sTObj + samIdx), zeroBgd);
      (void )WlzSetBackground(*(sSObj + samIdx), zeroBgd);
    }
  }
  /* Register the subsampled objects starting with the lowest resolution
   * (highest subsampling) and progressing to the unsampled objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(initTr == NULL)
    {
      rot0 = 0.0;
      tran0.vtX = 0.0;
      tran0.vtY = 0.0;
    }
    else
    {
      errNum = WlzAffineTransformPrimGet(initTr, &trPrim);
      rot0 = trPrim.theta;
      tran0.vtX = trPrim.tx;
      tran0.vtY = trPrim.ty;
    }
    conv = 1;
    samIdx = nSam - 1;
    sMaxRot = maxRot;
    sMaxTran.vtX =  maxTran.vtX / *(samFac + nSam - 1);
    sMaxTran.vtY =  maxTran.vtY / *(samFac + nSam - 1);
    while((errNum == WLZ_ERR_NONE) && conv && (samIdx >= 0))
    {
      /* Compute initial transform. */
      rot1 = rot0;
      tran1.vtX = tran0.vtX / *(samFac + samIdx);
      tran1.vtY = tran0.vtY / *(samFac + samIdx);
      samRegTr0 = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
      					        tran1.vtX, tran1.vtY, 0.0, 1.0, 
						rot1, 0.0, 0.0, 0.0, 0.0, 0,
						&errNum);
      /* Compute registration transform. */
      if(errNum == WLZ_ERR_NONE)
      {
        samRotCentre.vtX = rotCentre.vtX / *(samFac + samIdx);
        samRotCentre.vtY = rotCentre.vtY / *(samFac + samIdx);
	samRegTr1 = WlzRegCCorObjs2D1(*(sTObj + samIdx), *(sSObj + samIdx),
				      samRegTr0, samRotCentre,
				      trType, sMaxTran, sMaxRot, maxItr,
				      &conv, &cCor, &errNum);
        
      }
      if(samRegTr0)
      {
        (void )WlzFreeAffineTransform(samRegTr0);
	samRegTr0 = NULL;
      }
      if(samRegTr1)
      {
        samRegTr0 = samRegTr1;
	samRegTr1 = NULL;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzAffineTransformPrimGet(samRegTr0, &trPrim);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	rot0 = trPrim.theta;
	tran0.vtX = trPrim.tx * *(samFac + samIdx);
	tran0.vtY = trPrim.ty * *(samFac + samIdx);
      }
      /* Set registration limits. */
      sMaxRot = WLZ_M_PI / 24.0;
      sMaxTran.vtX = samFacStep * 3;
      sMaxTran.vtY = samFacStep * 3;
      --samIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstConv)
    {
      *dstConv = conv;
    }
    if(dstCCor)
    {
      *dstCCor = cCor;
    }
    regTr = samRegTr0;
  }
  /* Free subsampled objects. */
  if(sTObj)
  {
    for(samIdx = 0; samIdx < nSam; ++samIdx)
    {
      (void )WlzFreeObj(*(sTObj + samIdx));
    }
    AlcFree(sTObj);
  }
  if(sSObj)
  {
    for(samIdx = 0; samIdx < nSam; ++samIdx)
    {
      (void )WlzFreeObj(*(sSObj + samIdx));
    }
    AlcFree(sSObj);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/*!
* \return	Affine transform which brings the two objects into register.
* \ingroup	WlzRegistration
* \brief	Registers the two given 2D domain objects using a
*               frequency domain cross correlation.  An affine transform
*               is computed, which when applied to the source object
*               takes it into register with the target object.
* \param	tObj			The target object. Must have
*                                       been assigned.
* \param	sObj			The source object to be
*                                       registered with target object.
* \param	initTr			Initial affine transform
*                                       to be applied to the source
*                                       object prior to registration.
*                                       Only translations in x and y
*                                       and rotation about the z axis
*                                       are used. May be NULL which is
*                                       equivalent to an identity transform.
* \param	rotCentre		Centre of rotation.
* \param	trType			Required transform type.
* \param	maxTran			Maximum translation.
* \param	maxRot			Maximum rotation.
* \param	maxItr			Maximum number of iterations,
*                                       if \f$\leq\f$ 0 then infinite
					iterations are allowed.
* \param	dstConv			Destination ptr for the
*                                       convergence flag (non zero
*                                       on convergence), may be NULL.
* \param	dstCCor			Destination ptr for the cross
*                                       correlation value, may be NULL.
* \param	dstErr			Destination error pointer,
*                                       may be NULL.
*/
static WlzAffineTransform *WlzRegCCorObjs2D1(WlzObject *tObj, WlzObject *sObj,
					     WlzAffineTransform *initTr,
				  	     WlzDVertex2 rotCentre,
					     WlzTransformType trType,
					     WlzDVertex2 maxTran,
					     double maxRot, int maxItr,
					     int *dstConv, double *dstCCor, 
					     WlzErrorNum *dstErr)
{
  int		itr,
		conv;
  WlzAffineTransform *tTr0 = NULL,
  		*tTr1 = NULL,
		*curTr = NULL,
		*regTr = NULL;
  double	rot,
		lastRot,
  		cCor;
  WlzIVertex2	rotCentreI;
  WlzDVertex2	tran;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tranTol = 1.0,
  		rotTol = 0.01;

  /* Register for translation. */
  tran = WlzRegCCorObjs2DTran(tObj, sObj, initTr, maxTran, &cCor, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    tTr0 = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
    					 tran.vtX, tran.vtY, 0.0,
					 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
					 &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    curTr = WlzAffineTransformProduct(initTr, tTr0, &errNum);
    (void )WlzFreeAffineTransform(tTr0); tTr0 = NULL;
  }
  if(trType == WLZ_TRANSFORM_2D_TRANS)
  {
    conv = 1;
  }
  else /* trType == WLZ_TRANSFORM_2D_REG */
  {
    rotCentreI.vtX = WLZ_NINT(rotCentre.vtX);
    rotCentreI.vtY = WLZ_NINT(rotCentre.vtY);
    /* Iterate until translation is less tahn tollerance value or
     * number of itterations exceeds the maximum. */
    itr = 0;
    rot = 0.0;
    lastRot = 1.0;
    /* The convergence test must only test the translation and not the
     * rotation too, as the rotation has been applied when finding the
     * translation. */
    while((errNum == WLZ_ERR_NONE) &&
	  ((conv = (fabs(tran.vtX) <= tranTol) &&
	           (fabs(tran.vtY) <= tranTol)) == 0) &&
	  ((maxItr < 0) || (itr++ < maxItr)))
    {
      /* Register for rotation. */
      lastRot = rot;
      rot = WlzRegCCorObjs2DRot(tObj, sObj, rotCentreI, curTr, 
				maxRot, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	tTr0 = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					     0.0, 0.0, 0.0,
					     1.0, rot, 0.0, 0.0, 0.0, 0.0, 0,
					     &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tTr1 = WlzAffineTransformProduct(curTr, tTr0, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzFreeAffineTransform(tTr0); tTr0 = NULL;
	(void )WlzFreeAffineTransform(curTr);
	curTr = tTr1; tTr1 = NULL;
      }
      /* Register for translation. */
      if(errNum == WLZ_ERR_NONE)
      {
	tran = WlzRegCCorObjs2DTran(tObj, sObj, curTr, maxTran, &cCor,
				    &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tTr0 = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					     tran.vtX, tran.vtY, 0.0,
					     1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,
					     &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tTr1 = WlzAffineTransformProduct(curTr, tTr0, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzFreeAffineTransform(tTr0); tTr0 = NULL;
	(void )WlzFreeAffineTransform(curTr);
	curTr = tTr1; tTr1 = NULL;
      }
    }
    (void )WlzFreeAffineTransform(tTr0);
    (void )WlzFreeAffineTransform(tTr1);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    regTr = curTr;
    if(dstConv)
    {
      *dstConv = conv;
    }
    if(dstCCor)
    {
      *dstCCor = cCor;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/*!
* \return	Translation.
* \ingroup	WlzRegistration
* \brief	Registers the given 2D domain objects using a
*               frequency domain cross correlation, to find
*               the translation which has the highest cross
*               correlation value.
* \param	tObj			The target object. Must have
*                                       been assigned.
* \param	sObj			The source object to be
*                                       registered with target object.
*                                       Must have been assigned.
* \param	initTr			Initial affine transform
*                                       to be applied to the source
*                                       object prior to registration.
*               WlzDVertex2 maxTran:    Maximum translation.
* \param	maxTran			Maximum translation.
* \param	dstCCor			Destination ptr for the cross
*                                       correlation value, may be NULL.
* \param	dstErr			Destination error pointer,
*                                       may be NULL.
*/
static WlzDVertex2 WlzRegCCorObjs2DTran(WlzObject *tObj, WlzObject *sObj,
					WlzAffineTransform *initTr,
					WlzDVertex2 maxTran, double *dstCCor,
					WlzErrorNum *dstErr)
{
  int		oIdx;
  double	cCor = 0.0;
  double	sSq[2];
  double	**oAr[2];
  WlzIBox2	aBox;
  WlzIBox2	oBox[2];
  WlzIVertex2	aSz,
  		aOrg,
		winOrg,
		winRad,
		tran;
  WlzDVertex2	dstTran;
  WlzObject	*oObj[2],
  		*wObj[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstTran.vtX = 0.0;
  dstTran.vtY = 0.0;
  oAr[0] = oAr[1] = NULL;
  oObj[0] = oObj[1] = NULL;
  wObj[0] = wObj[1] = NULL;
  oObj[0] = WlzAssignObject(tObj, NULL);
  /* Transform source object. */
  if((initTr == NULL) || WlzAffineTransformIsIdentity(initTr, NULL))
  {
    oObj[1] = WlzAssignObject(sObj, NULL);
  }
  else
  {
    oObj[1] = WlzAffineTransformObj(sObj, initTr, WLZ_INTERPOLATION_NEAREST,
				   &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      oBox[oIdx] = WlzBoundingBox2I(oObj[oIdx], &errNum);
      ++oIdx;
    }
  }
  /* Compute windowed objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      winOrg.vtX = (oBox[oIdx].xMax + oBox[oIdx].xMin) / 2;
      winOrg.vtY = (oBox[oIdx].yMax + oBox[oIdx].yMin) / 2;
      winRad.vtX = (oBox[oIdx].xMax - oBox[oIdx].xMin) / 2;
      winRad.vtY = (oBox[oIdx].yMax - oBox[oIdx].yMin) / 2;
      wObj[oIdx] = WlzAssignObject(WlzWindow(oObj[oIdx], WLZ_WINDOWFN_HAMMING,
      					     winOrg, winRad, &errNum), NULL);
      ++oIdx;
    }
  }
  /* Compute the sums of squares for normalizing the cross correlation value if
   * it will be returned through the destination pointer. */
  if(dstCCor)
  {
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      (void )WlzGreyStats(wObj[oIdx], NULL, NULL, NULL, NULL, &(sSq[oIdx]),
      			  NULL, NULL, &errNum);
      ++oIdx;
    }
  }
  /* Create double arrays. */
  if(errNum == WLZ_ERR_NONE)
  {
    aBox.xMin = WLZ_MIN(oBox[0].xMin, oBox[1].xMin) - (int )(maxTran.vtX) + 1;
    aBox.yMin = WLZ_MIN(oBox[0].yMin, oBox[1].yMin) - (int )(maxTran.vtY) + 1;
    aBox.xMax = WLZ_MAX(oBox[0].xMax, oBox[1].xMax) + (int )(maxTran.vtX) + 1;
    aBox.yMax = WLZ_MAX(oBox[0].yMax, oBox[1].yMax) + (int )(maxTran.vtY) + 1;
    aOrg.vtX = aBox.xMin;
    aOrg.vtY = aBox.yMin;
    aSz.vtX = aBox.xMax - aBox.xMin + 1;
    aSz.vtY = aBox.yMax - aBox.yMin + 1;
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtX), aSz.vtX);
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtY), aSz.vtY);
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      errNum = WlzToArray2D((void ***)&(oAr[oIdx]), wObj[oIdx], aSz, aOrg, 0,
			    WLZ_GREY_DOUBLE);
      ++oIdx;
    }
  }
#ifdef WLZ_REGCCOR_DEBUG
  if(errNum == WLZ_ERR_NONE)
  {
    FILE	*fP = NULL;
    WlzObject	*cCObjT = NULL;
    
    cCObjT = WlzFromArray2D((void **)(oAr[0]), aSz, aOrg,
			    WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			    0.0, 1.0, 0, 1, &errNum);
    if(cCObjT)
    {
      if((fP = fopen("wObjT0.wlz", "w")) != NULL)
      {
	(void )WlzWriteObj(fP, cCObjT);
	(void )fclose(fP);
      }
      (void )WlzFreeObj(cCObjT);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    FILE	*fP = NULL;
    WlzObject	*cCObjT = NULL;
    
    cCObjT = WlzFromArray2D((void **)(oAr[1]), aSz, aOrg,
			    WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			    0.0, 1.0, 0, 1, &errNum);
    if(cCObjT)
    {
      if((fP = fopen("wObjT1.wlz", "w")) != NULL)
      {
	(void )WlzWriteObj(fP, cCObjT);
	(void )fclose(fP);
      }
      (void )WlzFreeObj(cCObjT);
    }
  }
#endif /* WLZ_REGCCOR_DEBUG */
  /* Cross correlate. */
  if(errNum == WLZ_ERR_NONE)
  {
    (void )AlgCrossCorrelate2D(oAr[0], oAr[1], aSz.vtX, aSz.vtY);
    AlgCrossCorrPeakXY(&(tran.vtX), &(tran.vtY), &cCor, oAr[0],
		       aSz.vtX, aSz.vtY, maxTran.vtX, maxTran.vtY);
  }
#ifdef WLZ_REGCCOR_DEBUG
  if(errNum == WLZ_ERR_NONE)
  {
    FILE	*fP = NULL;
    WlzObject	*cCObjT = NULL;
    
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      (void )WlzGreyStats(wObj[oIdx], NULL, NULL, NULL, NULL, &(sSq[oIdx]),
			  NULL, NULL, &errNum);
      ++oIdx;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      cCObjT = WlzFromArray2D((void **)(oAr[0]), aSz, aOrg,
			      WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			      0.0,
			      255.0 / (1.0 + (sqrt(sSq[0] * sSq[1]) *
			      		      aSz.vtX * aSz.vtY)),
			      0, 0, &errNum);
    }
    if(cCObjT)
    {
      if((fP = fopen("cCObjT.wlz", "w")) != NULL)
      {
	(void )WlzWriteObj(fP, cCObjT);
	(void )fclose(fP);
      }
      (void )WlzFreeObj(cCObjT);
    }
  }
#endif /* WLZ_REGCCOR_DEBUG */
  for(oIdx = 0; oIdx < 2; ++oIdx)
  {
    (void )WlzFreeObj(oObj[oIdx]);
    (void )WlzFreeObj(wObj[oIdx]);
    AlcDouble2Free(oAr[oIdx]);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstTran.vtX = tran.vtX;
    dstTran.vtY = tran.vtY;
    if(dstCCor)
    {
      cCor = cCor / (1.0 + (sqrt(sSq[0] * sSq[1]) * aSz.vtX * aSz.vtY));
      *dstCCor = cCor;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstTran);
}

/*!
* \return	Angle of roatation.
* \ingroup	WlzRegistration
* \brief	Polar samples then registers the given 2D domain objects
*               using a frequency domain cross correlation, to find
*               the angle of rotation about the given centre of rotation
*               which has the highest cross correlation value.
* \param	tObj			The target object. Must have
*                                       been assigned.
* \param	sObj			The source object to be
*                                       registered with target object.
*                                       Must have been assigned.
* \param	rotCentre		Coordinate of centre of rotation.
* \param	initTr			Initial affine transform
*                                       to be applied to the source
*                                       object prior to registration.
* \param	maxRot			Maximum rotation.
* \param	dstErr			Destination error pointer,
*                                       may be NULL.
*/
static double	WlzRegCCorObjs2DRot(WlzObject *tObj, WlzObject *sObj,
				    WlzIVertex2 rotCentre,
				    WlzAffineTransform *initTr, double maxRot,
				    WlzErrorNum *dstErr)
{
  int		oIdx,
  		angCnt;
  double	angInc,
  		dstRot = 0.0;
  WlzIBox2	aBox;
  WlzIBox2	oBox[2];
  WlzIVertex2	rot,
  		aSz,
  		aOrg,
		winRad,
		winOrg,
		rotPad;
  double	**oAr[2];
  WlzObject	*oObj[2],
  		*pObj[2],
		*wObj[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	rotCnt = 500;
  const double	distInc = 1.0;

  /* Assign the target and transform source objects. */
  oAr[0] = oAr[1] = NULL;
  oAr[0] = oAr[1] = NULL;
  oObj[0] = oObj[1] = NULL;
  pObj[0] = pObj[1] = NULL;
  wObj[0] = wObj[1] = NULL;
  oObj[0] = WlzAssignObject(tObj, NULL);
  if((initTr == NULL) || WlzAffineTransformIsIdentity(initTr, NULL))
  {
    oObj[1] = WlzAssignObject(sObj, NULL);
  }
  else
  {
    oObj[1] = WlzAffineTransformObj(sObj, initTr, WLZ_INTERPOLATION_NEAREST,
				   &errNum);
  }
  /* Compute rectangular to polar transformation. */
  if(errNum == WLZ_ERR_NONE)
  {
    angInc = (2.0 * (maxRot + WLZ_M_PI)) / rotCnt;
    angCnt = (2.0 * WLZ_M_PI) / angInc;
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      pObj[oIdx] = WlzPolarSample(oObj[oIdx], rotCentre, angInc, distInc,
	                          angCnt, 0, &errNum);
      ++oIdx;
    }
  }
  /* Compute windowed objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      oBox[oIdx] = WlzBoundingBox2I(pObj[oIdx], &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	winOrg.vtX = (oBox[oIdx].xMax + oBox[oIdx].xMin) / 2;
	winOrg.vtY = (oBox[oIdx].yMax + oBox[oIdx].yMin) / 2;
	winRad.vtX = (oBox[oIdx].xMax - oBox[oIdx].xMin) / 2;
	winRad.vtY = (oBox[oIdx].yMax - oBox[oIdx].yMin) / 2;
        wObj[oIdx] = WlzAssignObject(WlzWindow(pObj[oIdx],
					       WLZ_WINDOWFN_HAMMING, winOrg,
					       winRad, &errNum), NULL);
      }
      ++oIdx;
    }
  }
#ifdef WLZ_REGCCOR_DEBUG
  if(errNum == WLZ_ERR_NONE)
  {
    FILE	*fP = NULL;
    
    if((fP = fopen("wObjR0.wlz", "w")) != NULL) 
    {
      (void )WlzWriteObj(fP, wObj[0]); 
      (void )fclose(fP);
    }
    if((fP = fopen("wObjR1.wlz", "w")) != NULL)
    {
      (void )WlzWriteObj(fP, wObj[1]);
      (void )fclose(fP);
    }
  }
#endif /* WLZ_REGCCOR_DEBUG */
  /* Create 2D double arrays from the polar sampled objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      oBox[oIdx] = WlzBoundingBox2I(wObj[oIdx], &errNum);
      ++oIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aBox.xMin = WLZ_MIN(oBox[0].xMin, oBox[1].xMin);
    aBox.yMin = WLZ_MIN(oBox[0].yMin, oBox[1].yMin);
    aBox.xMax = WLZ_MAX(oBox[0].xMax, oBox[1].xMax);
    aBox.yMax = WLZ_MAX(oBox[0].yMax, oBox[1].yMax);
    rotPad.vtX = (aBox.xMax - aBox.xMin) / 2;
    rotPad.vtY = 1 + WLZ_NINT(maxRot / angInc);
    aBox.yMin -= rotPad.vtY;
    aBox.yMax += rotPad.vtY;
    aOrg.vtX = aBox.xMin;
    aOrg.vtY = aBox.yMin;
    aSz.vtX = aBox.xMax - aBox.xMin + 1;
    aSz.vtY = aBox.yMax - aBox.yMin + 1;
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtX), aSz.vtX);
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtY), aSz.vtY);
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      errNum = WlzToArray2D((void ***)&(oAr[oIdx]), wObj[oIdx], aSz, aOrg, 0,
			    WLZ_GREY_DOUBLE);
      ++oIdx;
    }
  }
#ifdef WLZ_REGCCOR_DEBUG
  if(errNum == WLZ_ERR_NONE)
  {
    FILE	*fP = NULL;
    WlzObject	*aObj = NULL;

    aObj = WlzFromArray2D((void **)(oAr[0]), aSz, aOrg,
    			  WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE, 0.0, 1.0,
			  0, 0, &errNum);
    if((fP = fopen("wObj0.wlz", "w")) != NULL) 
    {
      (void )WlzWriteObj(fP, aObj); 
      (void )fclose(fP);
    }
    WlzFreeObj(aObj);
    aObj = WlzFromArray2D((void **)(oAr[1]), aSz, aOrg,
    			  WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE, 0.0, 1.0,
			  0, 0, &errNum);
    if((fP = fopen("wObj1.wlz", "w")) != NULL)
    {
      (void )WlzWriteObj(fP, aObj);
      (void )fclose(fP);
    }
    WlzFreeObj(aObj);
  }
#endif /* WLZ_REGCCOR_DEBUG */
  /* Cross correlate. */
  if(errNum == WLZ_ERR_NONE)
  {
    (void )AlgCrossCorrelate2D(oAr[0], oAr[1], aSz.vtX, aSz.vtY);
    AlgCrossCorrPeakXY(&(rot.vtX), &(rot.vtY), NULL, oAr[0],
		       aSz.vtX, aSz.vtY, rotPad.vtX, rotPad.vtY);
    dstRot = rot.vtY * angInc;
    /* dstRot = -(rot.vtY) * angInc; */
  }
#ifdef WLZ_REGCCOR_DEBUG
  if(errNum == WLZ_ERR_NONE)
  {
    double	sSq[2];
    FILE	*fP = NULL;
    WlzObject	*cCObjR = NULL;
    
    oIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (oIdx < 2))
    {
      (void )WlzGreyStats(wObj[oIdx], NULL, NULL, NULL, NULL, &(sSq[oIdx]),
			  NULL, NULL, &errNum);
      ++oIdx;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      cCObjR = WlzFromArray2D((void **)(oAr[0]), aSz, aOrg,
			      WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			      0.0,
			      255.0 / (1.0 + (sqrt(sSq[0] * sSq[1]) *
			      		      aSz.vtX * aSz.vtY)),
			      0, 0, &errNum);
    }
    if(cCObjR)
    {
      if((fP = fopen("cCObjR.wlz", "w")) != NULL)
      {
	(void )WlzWriteObj(fP, cCObjR);
	(void )fclose(fP);
      }
      (void )WlzFreeObj(cCObjR);
    }
  }
#endif /* WLZ_REGCCOR_DEBUG */
  for(oIdx = 0; oIdx < 2; ++oIdx)
  {
    (void )WlzFreeObj(oObj[oIdx]);
    (void )WlzFreeObj(pObj[oIdx]);
    (void )WlzFreeObj(wObj[oIdx]);
    AlcDouble2Free(oAr[oIdx]);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstRot);
}

#ifdef WLZ_REGCCOR_TEST
/* Test main() for WlzRegCCorObjs(). */

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
  		option,
		maxItr = 10,
		verbose = 0,
		conv = 0,
  		ok = 1,
		usage = 0;
  double	cCor,
  		maxRot = 45 * (2.0 * WLZ_M_PI / 360.0);
  FILE		*fP = NULL;
  char		*outObjFileStr;
  char		*inObjFileStr[2];
  WlzDVertex2	maxTran;
  WlzTransformType trType = WLZ_TRANSFORM_2D_REG;
  WlzValues	nullVal;
  WlzDomain	outDom;
  WlzObject	*outObj;
  WlzObject	*inObj[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "htrvo:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  maxTran.vtX = maxTran.vtY = 500.0;
  inObj[0] = inObj[1] = NULL;
  nullVal.core = NULL;
  outDom.core = NULL;
  outObj = NULL;
  outObjFileStr = (char *)outFileStrDef;
  inObjFileStr[0] = inObjFileStr[1] = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'r':
        trType = WLZ_TRANSFORM_2D_REG;
	break;
      case 't':
        trType = WLZ_TRANSFORM_2D_TRANS;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
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
      if((optind + 2) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inObjFileStr[0] = *(argv + optind);
        inObjFileStr[1] = *(argv + optind + 1);
      }
    }
  }
  if(ok)
  {
    idx = 0;
    while(ok && (idx < 2))
    {
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP, &errNum),
	  				 NULL)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read object from file %s\n",
		       *argv, inObjFileStr);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
      }
      ++idx;
    }
  }
  if(ok)
  {
    outDom.t = WlzRegCCorObjs(inObj[0], inObj[1], NULL, trType,
			      maxTran, maxRot, maxItr, &conv,
			      &cCor, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to register objects.\n",
      		     argv[0]);
    }
  }
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_AFFINE_TRANS, outDom, nullVal, NULL, NULL,
    			 &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to make affine transform object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output affine transform object "
		     "to file %s (%s).\n",
		     *argv, outObjFileStr, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    (void )printf("%g\n", cCor);
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(outDom.core)
  {
    (void )WlzFreeAffineTransform(outDom.t);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s",
      *argv,
      " [-o<output object>] [-h] [<input object 0>] [<input object 1>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -v  Be verbose, outputs the cross correlation value to the standard\n"
      "      output, watch out that you dont send the affine transfrom to\n"
      "      the standard output too!\n"
      "  -r  Find registration transform.\n"
      "  -t  Find tranlation only transform.\n"
      "  -o  Output transform object file name.\n"
      "Computes a registration transform which registers the second of the\n"
      "given objects to the first using a cross correlation algorithm.\n");
  }
  return(!ok);
}
#endif /* WLZ_REGCCOR_TEST */
