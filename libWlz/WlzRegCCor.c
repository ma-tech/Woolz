#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzRegCCor.c
* Date:         January 2001
* Author:       Bill Hill
* Copyright:	2001 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions to register two objects using an
*		frequency domain cross correlation.
* TODO renormalise cross correlation values.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <float.h>
#include <Wlz.h>

#define WLZ_REGCCOR_TEST

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
				  double *dstCCor,
				  WlzErrorNum *dstErr);

/************************************************************************
* Function:	WlzRegCCorObjs
* Returns:	WlzAffineTransform *:	Affine transform which brings
*					the two objects into register.
* Purpose:	Registers the two given objects using a frequency domain
*		cross correlation.  An affine transform is computed,
*		which when applied to the source object takes it
*		into register with the target object.
* Global refs:	-
* Parameters:	WlzObject *tObj:	The target object.
*		WlzObject *sObj:	The source object to be
*					registered with target object.
*		WlzAffineTransform *initTr: Initial affine transform
*					to be applied to the source
*					object prior to registration.
*					May be NULL.
*		WlzTransformType trType: Required transform type.
*		WlzDVertex2 maxTran:	Maximum translation.
*		double maxRot:		Maximum rotation.
*		int maxItr:		Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
*		int *dstConv:		Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
*		double *dstCCor:	Destination ptr for the cross
*					correlation value, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzRegCCorObjs2D
* Returns:	WlzAffineTransform *:	Affine transform which brings
*					the two objects into register.
* Purpose:	Registers the two given 2D domain objects using a
*		frequency domain cross correlation.  An affine transform
*		is computed, which when applied to the source object
*		takes it into register with the target object.
*		A resolution pyramid is built from the given objects
*		and used to register the objects, progressing from
*		a low resolution towards the full resolution objects.
* Global refs:	-
* Parameters:	WlzObject *tObj:	The target object. Must have
*					been assigned.
*		WlzObject *sObj:	The source object to be
*					registered with target object.
*		WlzAffineTransform *initTr: Initial affine transform
*					to be applied to the source
*					object prior to registration.
*					Only translations in x and y
*					and rotation about the z axis
*					are used. May be NULL which is
*					equivalent to an identity transform.
*		WlzTransformType trType: Required transform type.
*		WlzDVertex2 maxTran:	Maximum translation.
*		double maxRot:		Maximum rotation.
*		int maxItr:		Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
*		int *dstConv:		Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
*		double *dstCCor:	Destination ptr for the cross
*					correlation value, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzAffineTransform *WlzRegCCorObjs2D(WlzObject *tObj, WlzObject *sObj,
					    WlzAffineTransform *initTr,
					    WlzTransformType trType,
					    WlzDVertex2 maxTran, double maxRot,
					    int maxItr,
					    int *dstConv, double *dstCCor,
					    WlzErrorNum *dstErr)
{
  int		tI0,
  		tI1,
		samIdx,
  		nSam,
		conv;
  double	cCor,
  		rot0,
		rot1,
		sMaxRot;
  WlzGreyType	gType;
  WlzPixelV	gV[4];
  WlzIVertex2	tIV0,
  		tIV1,
		winRad,
		winOrg;
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
  sBox = WlzBoundingBox2D(sObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    tBox = WlzBoundingBox2D(tObj, &errNum);
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
    *(sTObj + 0) = WlzAssignObject(tObj, NULL);
    *(sSObj + 0) = WlzAssignObject(sObj, NULL);
    while((errNum == WLZ_ERR_NONE) && (++samIdx < nSam))
    {
      *(samFac + samIdx) = *(samFac + samIdx - 1) * samFacStep;
      samFacV.vtX = samFacV.vtY = *(samFac + samIdx);
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

/************************************************************************
* Function:	WlzRegCCorObjs2D1
* Returns:	WlzAffineTransform *:	Affine transform which brings
*					the two objects into register.
* Purpose:	Registers the two given 2D domain objects using a
*		frequency domain cross correlation.  An affine transform
*		is computed, which when applied to the source object
*		takes it into register with the target object.
* Global refs:	-
* Parameters:	WlzObject *tObj:	The target object. Must have
*					been assigned.
*		WlzObject *sObj:	The source object to be
*					registered with target object.
*		WlzAffineTransform *initTr: Initial affine transform
*					to be applied to the source
*					object prior to registration.
*					Only translations in x and y
*					and rotation about the z axis
*					are used. May be NULL which is
*					equivalent to an identity transform.
*		WlzDVertex2 rotCentre:	Centre of rotation.
*		WlzTransformType trType: Required transform type.
*		WlzDVertex2 maxTran:	Maximum translation.
*		double maxRot:		Maximum rotation.
*		int maxItr:		Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
*		int *dstConv:		Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
*		double *dstCCor:	Destination ptr for the cross
*					correlation value, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
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
  		cCor;
  WlzIVertex2	rotCentreI;
  WlzDVertex2	tran;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tranTol = 1.0;

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
    while((errNum == WLZ_ERR_NONE) &&
	  ((conv = (fabs(tran.vtX) < tranTol) &&
	           (fabs(tran.vtY) < tranTol)) == 0) &&
	  ((maxItr < 0) || (itr++ < maxItr)))
    {
      /* Register for rotation. */
      if(errNum == WLZ_ERR_NONE)
      {
	rot = WlzRegCCorObjs2DRot(tObj, sObj, rotCentreI, curTr, 
				  maxRot, &cCor, &errNum);
      }
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

/************************************************************************
* Function:	WlzRegCCorObjs2DTran
* Returns:	WlzDVertex2:		Translation.
* Purpose:	Registers the given 2D domain objects using a
*		frequency domain cross correlation, to find
*		the translation which has the highest cross
*		correlation value.
* Global refs:	-
* Parameters:	WlzObject *tObj:	The target object. Must have
*					been assigned.
*		WlzObject *sObj:	The source object to be
*					registered with target object.
*					Must have been assigned.
*		WlzAffineTransform *initTr: Initial affine transform
*					to be applied to the source
*					object prior to registration.
*		WlzDVertex2 maxTran:	Maximum translation.
*		double *dstCCor:	Destination ptr for the cross
*					correlation value, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzDVertex2 WlzRegCCorObjs2DTran(WlzObject *tObj, WlzObject *sObj,
					WlzAffineTransform *initTr,
					WlzDVertex2 maxTran, double *dstCCor,
					WlzErrorNum *dstErr)
{
  double	tSSq,
  		sSSq,
		cCor = 0.0;
  WlzIBox2	aBox,
  		tBox,
  		sBox;
  WlzIVertex2	aSz,
  		aOrg,
		winOrg,
		winRad,
		tranI;
  WlzDVertex2	tran;
  double	**tAr = NULL,
  		**sAr = NULL;
  WlzObject	*sTrObj = NULL,
  		*tWObj = NULL,
		*sWObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tran.vtX = 0.0;
  tran.vtY = 0.0;
  /* Transform source object. */
  if((initTr == NULL) || WlzAffineTransformIsIdentity(initTr, NULL))
  {
    sTrObj = WlzAssignObject(sObj, NULL);
  }
  else
  {
    sTrObj = WlzAffineTransformObj(sObj, initTr, WLZ_INTERPOLATION_NEAREST,
				   &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tBox = WlzBoundingBox2D(tObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sBox = WlzBoundingBox2D(sTrObj, &errNum);
  }
  /* Compute windowed objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    winOrg.vtX = (tBox.xMax + tBox.xMin) / 2;
    winOrg.vtY = (tBox.yMax + tBox.yMin) / 2;
    winRad.vtX = (tBox.xMax - tBox.xMin) / 2;
    winRad.vtY = (tBox.yMax - tBox.yMin) / 2;
    tWObj = WlzAssignObject(WlzWindow(tObj, WLZ_WINDOWFN_HAMMING, winOrg,
    			     	      winRad, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    winOrg.vtX = (sBox.xMax + sBox.xMin) / 2;
    winOrg.vtY = (sBox.yMax + sBox.yMin) / 2;
    winRad.vtX = (sBox.xMax - sBox.xMin) / 2;
    winRad.vtY = (sBox.yMax - sBox.yMin) / 2;
    sWObj = WlzAssignObject(WlzWindow(sTrObj, WLZ_WINDOWFN_HAMMING, winOrg,
    			    winRad, &errNum), NULL);
  }
  /* Create double arrays. */
  if(errNum == WLZ_ERR_NONE)
  {
    aBox.xMin = WLZ_MIN(tBox.xMin, sBox.xMin) - ((int )(maxTran.vtX) + 1);
    aBox.yMin = WLZ_MIN(tBox.yMin, sBox.yMin) - ((int )(maxTran.vtY) + 1);
    aBox.xMax = WLZ_MAX(tBox.xMax, sBox.xMax) + ((int )(maxTran.vtX) + 1);
    aBox.yMax = WLZ_MAX(tBox.yMax, sBox.yMax) + ((int )(maxTran.vtY) + 1);
    aOrg.vtX = aBox.xMin;
    aOrg.vtY = aBox.yMin;
    aSz.vtX = aBox.xMax - aBox.xMin + 1;
    aSz.vtY = aBox.yMax - aBox.yMin + 1;
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtX), aSz.vtX);
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtY), aSz.vtY);
    errNum = WlzToArray2D((void ***)&tAr, tWObj, aSz, aOrg, 0,
    			  WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzToArray2D((void ***)&sAr, sWObj, aSz, aOrg, 0,
    			  WLZ_GREY_DOUBLE);
    
  }
  /* Cross correlate. */
  if(errNum == WLZ_ERR_NONE)
  {
    (void )AlgCrossCorrelate2D(tAr, sAr, aSz.vtX, aSz.vtY);
    AlgCrossCorrPeakXY(&(tranI.vtX), &(tranI.vtY), &cCor, tAr,
		       aSz.vtX, aSz.vtY, maxTran.vtX, maxTran.vtY);
#ifdef WLZ_REGCCOR_DEBUG
    {
      FILE	*fP = NULL;
      WlzObject	*cCObjT = NULL;
      
      cCObjT = WlzFromArray2D((void **)tAr, aSz, aOrg,
      			      WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			      0.0, 1.0, 0, 1, &errNum);
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
  }
  (void )WlzFreeObj(sWObj);
  (void )WlzFreeObj(tWObj);
  (void )WlzFreeObj(sTrObj);
  if(errNum == WLZ_ERR_NONE)
  {
    tran.vtX = tranI.vtX;
    tran.vtY = tranI.vtY;
    if(dstCCor)
    {
      *dstCCor = cCor;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tran);
}

/************************************************************************
* Function:	WlzRegCCorObjs2DRot
* Returns:	double:			Rotation in radians about the
*					given centre of rotation.
* Purpose:	Polar samples then registers the given 2D domain objects
*		using a frequency domain cross correlation, to find
*		the angle of rotation about the given centre of rotation
*		which has the highest cross correlation value.
* Global refs:	-
* Parameters:	WlzObject *tObj:	The target object. Must have
*					been assigned.
*		WlzObject *sObj:	The source object to be
*					registered with target object.
*					Must have been assigned.
*		WlzIVertex2 rotCentre:	Coordinate of centre of rotation.
*		WlzAffineTransform *initTr: Initial affine transform
*					to be applied to the source
*					object prior to registration.
*		double maxRot:		Maximum rotation.
*		double *dstCCor:	Destination ptr for the cross
*					correlation value, may be NULL.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static double	WlzRegCCorObjs2DRot(WlzObject *tObj, WlzObject *sObj,
				    WlzIVertex2 rotCentre,
				    WlzAffineTransform *initTr,
				    double maxRot, double *dstCCor,
				    WlzErrorNum *dstErr)
{
  int		rotY,
  		angCnt,
		rotPad;
  double	angInc,
		tSSq,
		sSSq,
		cCor = 0.0,
  		rot = 0.0;
  WlzIBox2	aBox,
  		tBox,
  		sBox;
  WlzIVertex2	aSz,
  		aOrg,
		winRad,
		winOrg;
  double	**tAr = NULL,
  		**sAr = NULL;
  WlzObject	*sTrObj = NULL,
		*tPObj = NULL,
  		*sPObj = NULL,
		*tWObj = NULL,
		*sWObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minAngCnt = 256;
  const double	distInc = 1.0;

  /* Transform source object. */
  if((initTr == NULL) || WlzAffineTransformIsIdentity(initTr, NULL))
  {
    sTrObj = WlzAssignObject(sObj, NULL);
  }
  else
  {
    sTrObj = WlzAffineTransformObj(sObj, initTr, WLZ_INTERPOLATION_NEAREST,
				   &errNum);
  }
  /* Compute polar transformation. */
  if(errNum == WLZ_ERR_NONE)
  {
    tBox = WlzBoundingBox2D(tObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sBox = WlzBoundingBox2D(sTrObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aBox.xMin = WLZ_MIN(tBox.xMin, sBox.xMin);
    aBox.yMin = WLZ_MIN(tBox.yMin, sBox.yMin);
    aBox.xMax = WLZ_MAX(tBox.xMax, sBox.xMax);
    aBox.yMax = WLZ_MAX(tBox.yMax, sBox.yMax);
    aSz.vtX = aBox.xMax - aBox.xMin + 1;
    aSz.vtY = aBox.yMax - aBox.yMin + 1;
    if((angCnt = WLZ_MIN(aSz.vtX, aSz.vtY)) < minAngCnt)
    {
      angCnt = minAngCnt;
    }
    angInc = (2.0 * WLZ_M_PI) / angCnt;
    tPObj = WlzPolarSample(tObj, rotCentre, angInc, distInc, angCnt, 0,
    			   &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sPObj = WlzPolarSample(sTrObj, rotCentre, angInc, distInc, angCnt, 0,
    			   &errNum);
  }
  /* Compute windowed objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    tBox = WlzBoundingBox2D(tPObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sBox = WlzBoundingBox2D(sPObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    winOrg.vtX = (tBox.xMax + tBox.xMin) / 2;
    winOrg.vtY = (tBox.yMax + tBox.yMin) / 2;
    winRad.vtX = (tBox.xMax - tBox.xMin) / 2;
    winRad.vtY = (tBox.yMax - tBox.yMin) / 2;
    tWObj = WlzAssignObject(WlzWindow(tPObj, WLZ_WINDOWFN_HAMMING, winOrg,
    			     	      winRad, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    winOrg.vtX = (sBox.xMax + sBox.xMin) / 2;
    winOrg.vtY = (sBox.yMax + sBox.yMin) / 2;
    winRad.vtX = (sBox.xMax - sBox.xMin) / 2;
    winRad.vtY = (sBox.yMax - sBox.yMin) / 2;
    sWObj = WlzAssignObject(WlzWindow(sPObj, WLZ_WINDOWFN_HAMMING, winOrg,
    			    winRad, &errNum), NULL);
  }
  /* Create double arrays. */
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef WLZ_REGCCOR_DEBUG
    {
      FILE	*fP = NULL;
      
      if((fP = fopen("sWObj.wlz", "w")) != NULL)
      {
        (void )WlzWriteObj(fP, sWObj);
	(void )fclose(fP);
      }
      if((fP = fopen("tWObj.wlz", "w")) != NULL) 
      {
        (void )WlzWriteObj(fP, tWObj); 
	(void )fclose(fP);
      }
    }
#endif /* WLZ_REGCCOR_DEBUG */
    aBox.xMin = WLZ_MIN(tBox.xMin, sBox.xMin);
    aBox.yMin = WLZ_MIN(tBox.yMin, sBox.yMin);
    aBox.xMax = WLZ_MAX(tBox.xMax, sBox.xMax);
    aBox.yMax = WLZ_MAX(tBox.yMax, sBox.yMax);
    rotPad = 1 + WLZ_NINT(maxRot / angInc);
    aBox.yMin -= rotPad;
    aBox.yMax += rotPad;
    aOrg.vtX = aBox.xMin;
    aOrg.vtY = aBox.yMin;
    aSz.vtX = aBox.xMax - aBox.xMin + 1;
    aSz.vtY = aBox.yMax - aBox.yMin + 1;
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtX), aSz.vtX);
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(aSz.vtY), aSz.vtY);
    errNum = WlzToArray2D((void ***)&tAr, tObj, aSz, aOrg, 0,
    			  WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzToArray2D((void ***)&sAr, sTrObj, aSz, aOrg, 0,
    			  WLZ_GREY_DOUBLE);
  }
  /* Cross correlate. */
  if(errNum == WLZ_ERR_NONE)
  {
    (void )AlgCrossCorrelate2D(tAr, sAr, aSz.vtX, aSz.vtY);
    AlgCrossCorrPeakY(&rotY, &cCor, tAr, aSz.vtY);
    rot = -rotY * angInc;
#ifdef WLZ_REGCCOR_DEBUG
    {
      FILE	*fP = NULL;
      WlzObject	*cCObjR = NULL;
      
      cCObjR = WlzFromArray2D((void **)tAr, aSz, aOrg,
      			      WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
			      0.0, 1.0, 0, 1, &errNum);
      if(cCObjR)
      {
	if((fP = fopen("cCObjR.wlz", "w")) != NULL)
	{
	  (void )WlzWriteObj(fP, cCObjR);
	  (void )fclose(fP);
	}
	(void )WlzFreeObj(cCObjR);
      }
      /* HACK */ exit(0);
    }
#endif /* WLZ_REGCCOR_DEBUG */
  }
  (void )WlzFreeObj(tWObj);
  (void )WlzFreeObj(sWObj);
  (void )WlzFreeObj(tPObj);
  (void )WlzFreeObj(sPObj);
  (void )WlzFreeObj(sTrObj);
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstCCor)
    {
      *dstCCor = cCor;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rot);
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
  static char	optList[] = "htro:";
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
      "  -r  Find registration transform.\n"
      "  -t  Find tranlation only transform.\n"
      "  -o  Output transform object file name.\n"
      "Computes a registration transform which registers the second of the\n"
      "given objects to the first using a cross correlation algorithm.\n");
  }
  return(!ok);
}
#endif /* WLZ_REGCCOR_TEST */
