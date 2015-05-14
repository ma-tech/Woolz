#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructConstruct3D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructConstruct3D.c
* \author       Bill Hill
* \date         April 1999
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
* \brief	Provides functions for the automatic registration of a
*		single pair of serial sections.
* \ingroup	Reconstruct
*/
#include <Reconstruct.h>
#include <string.h>
#include <float.h>


static void			RecConstructForceIntScale(
				  double *scale);
static int			RecConstruct3DPlaneCount(
				  HGUDlpList *secList,
				  int allPlanes,
				  int gvnPlaneMin,
				  int gvnPlaneMax,
				  int *dstPlaneMin,
				  int *dstPlaneMax,
				  RecError *dstErr);
static RecError			RecConstructTargetHist(
				  WlzObject **dstHist,
				  HGUDlpList *secList,
				  int targetIdx,
				  char **eMsg);

/*!
* \return	Reconstruction error code.
* \ingroup	Reconstruct	
* \brief	Constructs a 3D woolz object from the given section list.
*		Each  image is read its section file, transformed and used to
*		build a 3D woolz object. The z component of the scale factor is
*		treated as an integer, ie sections may be either omitted or
*		duplicated.
* \param	dstObj			Destination pointer for new 3D object.
* \param	secList			Given section list.
* \param	confLimit		Confidence limit, any section with
*					a confidence/correlation value less
*					than the limit uses an identity
*					transform.
* \param	gaussSamFlg		Use gaussian subsampling if non-zero.
* \param	interp			Interpolation to use.
* \param	fastSamFlg		Use fast sampling code if non-zero.
* \param	greedyFlg		Be resource greedy if non-zero.
* \param	intScaleFlg		Use integer scale factor if non-zero.
* \param	scale			Scale factor.
* \param	matchHistFlg		Match histograms.
* \param	matchIdx		Section index for target histogram.
* \param	srcReg			Source region in the non-transformed
* 					sections, NULL implies all.
* \param	dstReg			Destination region in the transformed
* 					sections, NULL implies all.
* \param	eMsg			Destination pointer for error
*					messages.
*/
RecError	RecConstruct3DObj(WlzObject **dstObj, HGUDlpList *secList,
			       double confLimit,
			       int gaussSamFlg, WlzInterpolationType interp,
			       int fastSamFlg, int greedyFlg,
			       int intScaleFlg, WlzDVertex3 scale,
			       int matchHistFlg, int matchIdx,
			       WlzDBox3 *srcReg, WlzDBox3 *dstReg,
			       char **eMsg)
{
  int		affineScaleFlg,
  		allSrcFlg,
		firstPlaneFlg = 1,
  		numPlanes,
		plane1 = 0,
		lastPl = 0;
  double	plane1D,
  		lastPlD,
		gaussHWidth = 1.0;
  WlzObject	*tObj0 = NULL,
  		*tObj1 = NULL,
		*obj3D = NULL,
		*dstHistObj = NULL;
  WlzDomain	*doms;
  WlzValues 	*vals;
  WlzAffineTransformPrim prim;
  WlzAffineTransform *scaleTr = NULL,
  		*tTr0 = NULL,
		*tTr1 = NULL;
  RecSection	*sec;
  HGUDlpListItem *secItem;
  WlzDomain	dstDom;
  WlzValues	dstVal;
  WlzIVertex3	samFacI3;
  WlzIBox2      srcClip2I,
  		dstClip2I;
  WlzIBox3	/* dstClip3I, */
  		dstRegI3;
  WlzSampleFn	samFn = WLZ_SAMPLEFN_NONE;
  WlzPixelV	bgdPix;
  RecError	errFlag = REC_ERR_NONE;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecConstruct3DObj FE 0x%lx 0x%lx %g "
	   "%d %d %d %d"
	   "%d {%g %g %g} "
	   "%d %d "
	   "0x%lx 0x%lx ",
	   "0x%lx\n",
	   (unsigned long )dstObj, (unsigned long )secList, confLimit,
	   gaussSamFlg, (int )interp, fastSamFlg, greedyFlg,
	   intScaleFlg, scale.vtX, scale.vtY, scale.vtZ,
	   matchHistFlg, matchIdx,
	   (unsigned long )srcReg, (unsigned long)dstReg,
	   (unsigned long )eMsg));
  dstDom.core = NULL;
  dstVal.core = NULL;
  dstRegI3.xMin = dstRegI3.yMin = dstRegI3.zMin = 0;
  dstRegI3.xMax = dstRegI3.yMax = dstRegI3.zMax = 0;
  if((dstObj == NULL) ||			   /* Check given parameters */
     (secList == NULL) ||
     (HGUDlpListCount(secList) <= 0) ||
     ((interp != WLZ_INTERPOLATION_NEAREST) && 
      (interp != WLZ_INTERPOLATION_LINEAR)) ||
     (scale.vtX < (REC_TRANS_SCALE_MIN - DBL_EPSILON)) ||
     (scale.vtX > (REC_TRANS_SCALE_MAX + DBL_EPSILON)) ||
     (scale.vtY < (REC_TRANS_SCALE_MIN - DBL_EPSILON)) ||
     (scale.vtY > (REC_TRANS_SCALE_MAX + DBL_EPSILON)) ||
     (scale.vtZ < (REC_TRANS_SCALE_MIN - DBL_EPSILON)) ||
     (scale.vtZ > (REC_TRANS_SCALE_MAX + DBL_EPSILON)) ||
     (fabs(scale.vtX - scale.vtY) > DBL_EPSILON) ||
     ((srcReg != NULL) &&
      (((srcReg->xMax - srcReg->xMin) < (DBL_EPSILON)) ||
       ((srcReg->yMax - srcReg->yMin) < (DBL_EPSILON)) ||
       ((srcReg->zMax - srcReg->zMin) < (DBL_EPSILON)))) ||
     ((dstReg != NULL) &&
      (((dstReg->xMax - dstReg->xMin) < (DBL_EPSILON)) ||
       ((dstReg->yMax - dstReg->yMin) < (DBL_EPSILON)) ||
       ((dstReg->zMax - dstReg->zMin) < (DBL_EPSILON)))))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(intScaleFlg)
    {
      RecConstructForceIntScale(&(scale.vtX));
      RecConstructForceIntScale(&(scale.vtY));
      RecConstructForceIntScale(&(scale.vtZ));
    }
    if(srcReg)
    {
      srcClip2I.xMin = WLZ_NINT(srcReg->xMin);
      srcClip2I.xMax = WLZ_NINT(srcReg->xMax);
      srcClip2I.yMin = WLZ_NINT(srcReg->yMin);
      srcClip2I.yMax = WLZ_NINT(srcReg->yMax);
    }
    /*
    if(dstReg)
    {
      dstClip3I.xMin = dstClip2I.xMin = WLZ_NINT(dstReg->xMin);
      dstClip3I.xMax = dstClip2I.xMax = WLZ_NINT(dstReg->xMax);
      dstClip3I.yMin = dstClip2I.yMin = WLZ_NINT(dstReg->yMin);
      dstClip3I.yMax = dstClip2I.yMax = WLZ_NINT(dstReg->yMax);
      dstClip3I.zMin =  WLZ_NINT(dstReg->zMin);
      dstClip3I.zMax =  WLZ_NINT(dstReg->zMax);
    }
    */
    if(confLimit < 0.0)
    {
      confLimit = 0.0;
    }
    else if(confLimit > 1.0)
    {
      confLimit = 1.0;
    }
    if(srcReg || dstReg)
    {
      allSrcFlg = 0;
      if(srcReg && dstReg)
      {
        plane1D = WLZ_MAX(srcReg->zMin, dstReg->zMin);
        lastPlD = WLZ_MIN(srcReg->zMax, dstReg->zMax);
      }
      else if(srcReg)
      {
        plane1D = srcReg->zMin;
        lastPlD = srcReg->zMax;
      }
      else /* if(dstReg) */
      {
        plane1D = dstReg->zMin;
        lastPlD = dstReg->zMax;
      }
      numPlanes = RecConstruct3DPlaneCount(secList, allSrcFlg,
      					   WLZ_NINT(plane1D),
					   WLZ_NINT(lastPlD),
				           &plane1, &lastPl, &errFlag);
    }
    else
    {
      allSrcFlg = 1;
      numPlanes = RecConstruct3DPlaneCount(secList, allSrcFlg, 0, 0,
				           &plane1, &lastPl, &errFlag);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    dstRegI3.zMin = plane1;
    dstRegI3.zMax = lastPl;
    if((numPlanes <= 0) ||
       ((dstDom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
       				       plane1, lastPl,
				       0, 1,
				       0, 1, &wlzErr)) == NULL) ||
       (wlzErr != WLZ_ERR_NONE) ||
       ((dstVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
       					  plane1, lastPl,
					  bgdPix, NULL, &wlzErr)) == NULL) ||
       (wlzErr != WLZ_ERR_NONE) ||
       ((obj3D = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstVal,
       			     NULL, NULL, &wlzErr)) == NULL) ||
       (wlzErr != WLZ_ERR_NONE))
    {
      errFlag = (wlzErr == WLZ_ERR_NONE)? REC_ERR_WLZ: RecErrorFromWlz(wlzErr);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    /* Compute destination histogram */
    if(matchHistFlg)
    {
      errFlag = RecConstructTargetHist(&dstHistObj, secList, matchIdx, eMsg);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(fastSamFlg && intScaleFlg && (scale.vtX < (1.0 - DBL_EPSILON)))
    {
      affineScaleFlg = 0;
      if(gaussSamFlg)
      {
        samFn = WLZ_SAMPLEFN_GAUSS;
      }
      else
      {
        samFn = WLZ_SAMPLEFN_POINT;
      }
      samFacI3.vtX = WLZ_NINT(1.0 / scale.vtX);
      samFacI3.vtY = WLZ_NINT(1.0 / scale.vtY);
      samFacI3.vtZ = 1;
    }
    else
    {
      if(gaussSamFlg)
      {
	if(scale.vtX > 1.0 + DBL_EPSILON)
	{
	  gaussSamFlg = 0;
	}
	else
	{
	  gaussHWidth = 1.0 / scale.vtX; 
	}
      }
      affineScaleFlg = 1;
      scaleTr = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					      0.0, 0.0, 0.0,
					      scale.vtX, 0.0, 0.0,
					      0.0, 0.0, 0.0, 0, &wlzErr);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    secItem = HGUDlpListHead(secList);
    sec = (RecSection *)HGUDlpListEntryGet(secList, secItem);
    doms = dstDom.p->domains;
    vals = dstVal.vox->values;
    dstDom.p->voxel_size[0] = 1.0;
    dstDom.p->voxel_size[1] = 1.0;
    dstDom.p->voxel_size[2] = 1.0;
    if(errFlag == REC_ERR_NONE)
    {
      tTr0 = WlzAssignAffineTransform(
	     WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
					   0.0, 0.0, 0.0,
					   1.0, 0.0, 0.0,
					   0.0, 0.0, 0.0, 0, &wlzErr), NULL);
      errFlag = RecErrorFromWlz(wlzErr);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    /* Loop over all sections in the section list */
    while(secItem && (errFlag == REC_ERR_NONE))
    {
      sec = (RecSection *)HGUDlpListEntryGet(secList, secItem);
      if(allSrcFlg ||
         ((sec->index >= plane1) &&
	  (sec->index <= lastPl)))	               /* Check plane limits */
      {
	/* Compute the transform to be applied to this section */
	if(errFlag == REC_ERR_NONE)
	{
	  if(errFlag == REC_ERR_NONE)
	  {
	    if((confLimit < DBL_EPSILON) ||
	       (sec->correl < confLimit))	  /* Use identity transform? */
	    {
	      if(sec->transform)
	      {
		(void )WlzAffineTransformPrimGet(sec->transform, &prim);
		REC_DBG((REC_DBG_3D|REC_DBG_LVL_2),
			("RecConstruct3DObj 01 0x%lx %d 0x%lx %g %g %g %g\n",
			 (unsigned long )sec,
			 sec->index,
			 (unsigned long )(sec->transform),
			 prim.tx,
			 prim.ty,
			 prim.scale,
			 prim.theta));
	      }
	      else
	      {
		REC_DBG((REC_DBG_3D|REC_DBG_LVL_2),
			("RecConstruct3DObj 01 0x%0\n"));
	      }
	      tTr1 = WlzAssignAffineTransform(
	      	     WlzAffineTransformProduct(sec->transform, tTr0,
		     			       &wlzErr), NULL);
	      errFlag = RecErrorFromWlz(wlzErr);
	      (void )WlzFreeAffineTransform(tTr0);
	      tTr0 = NULL;
	    }
	    else
	    {
	      tTr1 = tTr0;
	      tTr0 = NULL;
	    }
	    if(tTr1)
	    {
	      (void )WlzAffineTransformPrimGet(sec->transform, &prim);
	      REC_DBG((REC_DBG_3D|REC_DBG_LVL_2),
		      ("RecConstruct3DObj 02 0x%lx %g %g %g %g\n",
		       (unsigned long )tTr1,
		       prim.tx, prim.ty, prim.scale, prim.theta));
	    }
	    else
	    {
	      REC_DBG((REC_DBG_3D|REC_DBG_LVL_2),
	      	      ("RecConstruct3DObj 02 0x0\n"));
	    }
	  }
	}
	if(errFlag == REC_ERR_NONE)	        /* Set transform for scaling */
	{
	  if(scaleTr)
	  {
	    tTr0 = WlzAssignAffineTransform(
	    	   WlzAffineTransformProduct(tTr1, scaleTr, &wlzErr), NULL);
	    errFlag = RecErrorFromWlz(wlzErr);
	  }
	  else
	  {
	    tTr0 = WlzAssignAffineTransform(tTr1, NULL);
	  }
	  if(tTr0)
	  {
	    (void )WlzAffineTransformPrimGet(tTr0, &prim);
	    REC_DBG((REC_DBG_3D|REC_DBG_LVL_3),
		    ("RecConstruct3DObj 03 0x%lx %g %g %g %g\n",
		     (unsigned long )tTr0,
		     prim.tx, prim.ty, prim.scale, prim.theta));
	  }
	  else
	  {
	    REC_DBG((REC_DBG_3D|REC_DBG_LVL_2),
		    ("RecConstruct3DObj 03 0x0\n"));
	  }
	}
	/* Pre-process section object (ie prior to transformation) */
	if(errFlag == REC_ERR_NONE)
	{
          if((errFlag = RecFileSecObjRead(sec, eMsg)) == REC_ERR_NONE)
	  {
	    tObj0 = WlzAssignObject(WlzCopyObject(sec->obj, &wlzErr), NULL);
	    errFlag = RecErrorFromWlz(wlzErr);
	    RecFileSecObjFree(sec);
	  }
	}
	if(errFlag == REC_ERR_NONE)
	{
	  /* Match histogram. */
	  if(matchHistFlg)	
	  {
	    errFlag = RecErrorFromWlz(
	    	      WlzHistogramMatchObj(tObj0, dstHistObj, 0, 0,  0.0, 1.0,
		      			   1));
	  }
	}
	if(errFlag == REC_ERR_NONE)
	{
	  /* Clip source region. */
	  if(srcReg)
	  {
	    tObj1 = WlzAssignObject(
	    	   WlzClipObjToBox2D(tObj0, srcClip2I, &wlzErr), NULL);
	    errFlag = RecErrorFromWlz(wlzErr);
	    if(tObj0)
	    {
	      (void )WlzFreeObj(tObj0);
	      tObj0 = tObj1;
	      tObj1 = NULL;
	    }
	  }
	}
	if(errFlag == REC_ERR_NONE)
	{
	  /* Apply gaussian. */
	  if(gaussSamFlg && affineScaleFlg)
	  {
	    tObj1 = WlzAssignObject(
	    	    WlzGauss2(tObj0, gaussHWidth, gaussHWidth, 0, 0,
		              &wlzErr), NULL);
	    errFlag = RecErrorFromWlz(wlzErr);
	    if(tObj0)
	    {
	      (void )WlzFreeObj(tObj0);
	      tObj0 = tObj1;
	      tObj1 = NULL;
	    }
	  }
	}
	/* Apply section transform. */
	if(errFlag == REC_ERR_NONE)
	{
	  tObj1 = WlzAssignObject(
	  	  WlzAffineTransformObj(tObj0, tTr0, interp, &wlzErr), NULL);
	  errFlag = RecErrorFromWlz(wlzErr); 
	  if(tObj0)
	  {
	    (void )WlzFreeObj(tObj0);
	    tObj0 = tObj1;
	    tObj1 = NULL;
	  }
	}
	if(errFlag == REC_ERR_NONE)
	{
	  /* Sample if not using affine scale. */
	  if((affineScaleFlg == 0) && (scale.vtX < 1.0))
	  {
	    tObj1 = WlzAssignObject(
	    	    WlzSampleObj(tObj0, samFacI3, samFn, &wlzErr), NULL);
	    errFlag = RecErrorFromWlz(wlzErr);
	    if(tObj0)
	    {
	      (void )WlzFreeObj(tObj0);
	      tObj0 = tObj1;
	      tObj1 = NULL;
	    }
	  }
	}
	if(errFlag == REC_ERR_NONE)
	{
	  if(tObj0 && (tObj0->type == WLZ_2D_DOMAINOBJ))
	  {
	    if(dstReg)
	    {
	      /* Clip destination. */
	      tObj1 = WlzAssignObject(
		      WlzClipObjToBox2D(tObj0, dstClip2I, &wlzErr), NULL);
	      *doms = WlzAssignDomain(tObj1->domain, NULL);
	      *vals = WlzAssignValues(tObj1->values, NULL);
	      (void )WlzFreeObj(tObj1); 
	      tObj1 = NULL;
	    }
	    else
	    {
	      *doms = WlzAssignDomain(tObj0->domain, NULL);
	      *vals = WlzAssignValues(tObj0->values, NULL);
	    }
	    if(firstPlaneFlg)
	    {
	      dstRegI3.xMin = tObj0->domain.i->kol1;
	      dstRegI3.xMax = tObj0->domain.i->lastkl;
	      dstRegI3.yMin = tObj0->domain.i->line1;
	      dstRegI3.yMax = tObj0->domain.i->lastln;
	    }
	    else
	    {
	      dstRegI3.xMin = WLZ_MIN(dstRegI3.xMin,
				      tObj0->domain.i->kol1);
	      dstRegI3.xMax = WLZ_MAX(dstRegI3.xMax,
				      tObj0->domain.i->lastkl);
	      dstRegI3.yMin = WLZ_MIN(dstRegI3.yMin,
				      tObj0->domain.i->line1);
	      dstRegI3.yMax = WLZ_MAX(dstRegI3.yMax,
				      tObj0->domain.i->lastln);
	    }
	    bgdPix = WlzGetBackground(tObj0, &wlzErr);
	    firstPlaneFlg = 0;
	  }
	  else
	  {
	    (*doms).core = NULL;
	    (*vals).core = NULL;
	  }
	  errFlag = RecErrorFromWlz(wlzErr);
	  if(tObj0)
	  {
	    (void )WlzFreeObj(tObj0);
	    tObj0 = NULL;
	  }
	  if(tObj1)
	  {
	    (void )WlzFreeObj(tObj1);
	    tObj1 = NULL;
	  }
	  ++doms;
	  ++vals;
	  (void )WlzFreeAffineTransform(tTr0);
	  tTr0 = tTr1;
	  tTr1 = NULL;
	}
      }
      secItem = HGUDlpListNext(secList, secItem);
    }
    if(dstHistObj)
    {
      (void )WlzFreeObj(dstHistObj);
      dstHistObj = NULL;
    }
    if(scaleTr)
    {
      (void )WlzFreeAffineTransform(scaleTr);
      scaleTr = NULL;
    }
    if(tTr0)
    {
      (void )WlzFreeAffineTransform(tTr0);
    }
    if(tTr1)
    {
      (void )WlzFreeAffineTransform(tTr1);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    dstDom.p->line1 = dstRegI3.yMin;
    dstDom.p->lastln = dstRegI3.yMax;
    dstDom.p->kol1 = dstRegI3.xMin;
    dstDom.p->lastkl = dstRegI3.xMax;
    dstDom.p->plane1 = dstRegI3.zMin;
    dstDom.p->lastpl = dstRegI3.zMax;
    errFlag = RecErrorFromWlz(WlzSetBackground(obj3D, bgdPix));
  }
  if(errFlag == REC_ERR_NONE)
  {
    *dstObj = obj3D;
  }
  REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecConstruct3DObj FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \ingroup	Reconstruct
* \brief	Forces the scale to be an integer ratio: Eg 0.5 or 2.0.
* \param	scale			Source and destination pointer for
*					the scale.
*/
static void	RecConstructForceIntScale(double *scale)
{
  if(scale)
  {
    if(*scale > (1.0 + DBL_EPSILON))
    {
      *scale = WLZ_NINT(*scale);
    }
    else if((*scale > DBL_EPSILON) && (*scale < (1.0 - DBL_EPSILON)))
    {
      *scale = 1.0 / WLZ_NINT( 1.0 / *scale);
    }
    else
    {
      *scale = 1.0;
    }
  }
}

/*!
* \return	Number of planes within given limits.
* \ingroup	Reconstruct
* \brief	Counts the number of sections for which the image plane lies
* 		within the given plane limits.
* \param	secList			Given section list.
* \param	allPlanes		Count all planes if non zero.
* \param	gvnPlaneMin		Given plane minimum limit.
* \param	gvnPlaneMax		Given plane maximum limit.
* \param	dstPlaneMin		Destination pointer for the minimum
*					plane limit.
* \param	dstPlaneMax		Destination pointer for the maximum
*					plane limit.
* \param	dstErr			Destination pointer for error code.
*/
static int	RecConstruct3DPlaneCount(HGUDlpList *secList, int allPlanes,
			       		 int gvnPlaneMin, int gvnPlaneMax,
					 int *dstPlaneMin, int *dstPlaneMax,
					 RecError *dstErr)
{
  int		lastIdx,
  		planeCount = 0;
  RecSection	*sec;
  HGUDlpListItem *secItem;
  RecError 	errFlag = REC_ERR_NONE;

  if(secList)
  {
    secItem = HGUDlpListHead(secList);
    while((errFlag == REC_ERR_NONE) && secItem && 
          ((sec = (RecSection *)HGUDlpListEntryGet(secList,
	  					   secItem)) != NULL))
    {
      if(allPlanes ||
        ((sec->index >= gvnPlaneMin) && (sec->index <= gvnPlaneMax)))
      {
	if(planeCount == 0)
	{
	  *dstPlaneMin = sec->index;
	  *dstPlaneMax = sec->index;
	  lastIdx = sec->index;
	}
	else if(sec->index <= lastIdx)
	{
	  errFlag = REC_ERR_LIST;
	}
	else
	{
	  if(sec->index > *dstPlaneMax)
	  {
	    *dstPlaneMax = sec->index;
	  }
	  else  if(sec->index < *dstPlaneMin)
	  {
	    *dstPlaneMin = sec->index;
	  }
	}
        ++planeCount;
      }
      secItem = HGUDlpListNext(secList, secItem);
    }
  }
  if(dstErr)
  {
    *dstErr = errFlag;
  }
  return(planeCount);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Creates a histogram which is the mean of the section
* 		histograms. This assumes all objects have been read in.
* \param	dstHist			Destination pointer for mean histogram.
* \param	secList			Given section list.
* \param	targetIdx		Target section index.
* \param	eMsg			Destination pointer for error
*					messages.
*/
static RecError	RecConstructTargetHist(WlzObject **dstHist,
				       HGUDlpList *secList,
				       int targetIdx, char **eMsg)
{
  int		found = 0;
  RecError	errFlag = REC_ERR_NONE;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  HGUDlpListItem *secItem;
  RecSection	*sec;
  WlzObject	*histObj = NULL;
  char		errBuf[256];


  if((secList == NULL) ||
     (HGUDlpListCount(secList) <= 0) ||
     ((secItem = HGUDlpListHead(secList)) == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    /* Find target section */
    while((found == 0) && secItem && (errFlag == REC_ERR_NONE))
    {
      if((sec = (RecSection *)HGUDlpListEntryGet(secList, secItem)) == NULL)
      {
        errFlag = REC_ERR_FUNC;
      }
      else if(sec->index == targetIdx)
      {
        found = 1;
      }
      else
      {
	secItem = HGUDlpListNext(secList, secItem);
      }
    }
    if(found == 0)
    {
      errFlag = REC_ERR_FUNC;
    }
    if(errFlag != REC_ERR_NONE)
    {
      (void )sprintf(errBuf, "index %d not in section list.", targetIdx);
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    /* Read target section image */
    errFlag = RecFileSecObjRead(sec, eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    /* Compute target histogram */
    histObj = WlzHistogramObj(sec->obj, 256, 0.0, 1.0, &wlzErr);
    RecFileSecObjFree(sec);
    errFlag = RecErrorFromWlz(wlzErr);
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(dstHist)
    {
      *dstHist = histObj;
    }
  }
  else
  {
    if(eMsg && (*eMsg == NULL))
    {
      *eMsg = AlcStrDup(errBuf);
    }
  }
  return(errFlag);
}
