#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:	Woolz
* Title:	WlzContour.c
* Date: 	September 1999
* Author:	Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:	Functions for extracting contours from Woolz objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 03-03-2K bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

/* #define WLZ_CONTOUR_DEBUG */
static void 	WlzContourTestOutPSLn2D(FILE *fP,
					WlzDVertex2 org, WlzDVertex2 dst);

/************************************************************************
* WlzContourTriIsn2D: Classification of intersection of line segment
* with a triangle
************************************************************************/
typedef enum
{
  WLZ_CONTOUR_TIC2D_NONE,                                  /*No intersection */
  WLZ_CONTOUR_TIC2D_N1N0,                                 /* Node 1 - node 0 */
  WLZ_CONTOUR_TIC2D_N0N2,                                 /* Node 0 - node 2 */
  WLZ_CONTOUR_TIC2D_N2N1,                                 /* Node 2 - node 1 */
  WLZ_CONTOUR_TIC2D_N1S02,                              /* Node 1 - side 0-2 */
  WLZ_CONTOUR_TIC2D_N0S21,                              /* Node 0 - side 2-1 */
  WLZ_CONTOUR_TIC2D_N2S10,                              /* Node 2 - side 1-0 */
  WLZ_CONTOUR_TIC2D_S10S02,                           /* Side 1-0 - side 0-2 */
  WLZ_CONTOUR_TIC2D_S02S21,                           /* Side 0-2 - side 2-1 */
  WLZ_CONTOUR_TIC2D_S21S10                            /* Side 2-1 - side 1-0 */
} WlzContourTriIsn2D;

static WlzContour	*WlzContourIsoObj2D(
			  WlzObject *srcObj,
			  double isoVal,
			  int minNod,
			  int minEdg,
			  WlzErrorNum *dstErr);
static WlzContour	*WlzContourIsoObj3D(
			  WlzObject *srcObj,
			  double isoVal,
			  int minNod,
			  int minEdg,
			  WlzErrorNum *dstErr);
static WlzContour	*WlzContourGrdObj2D(WlzObject *srcObj,
				double minGrd, double ftrPrm,
				int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContour	*WlzContourGrdObj3D(
			  WlzObject *srcObj,
			  double minGrd,
			  double ftrPrm,
			  int minNod,
			  int minEdg,
			  WlzErrorNum *dstErr);
static WlzDVertex2	WlzContourItpTriSide(
			  double ht0,
			  double ht1,
			  WlzDVertex2 org,
			  WlzDVertex2 dst);
static WlzErrorNum	WlzContourIsoCube2D(
			  WlzContour *ctr,
			  double isoVal,
			  double *ln0,
			  double *ln1,
			  WlzDVertex2 sqOrg,
			  int *vtxSearchIdx);
static WlzErrorNum	WlzContourIsoCube3D(WlzContour *ctr,
			  double isoVal,
			  double *pn0ln0,
			  double *pn0ln1,
			  double *pn1ln0,
			  double *pn1ln1,
			  WlzDVertex3 cbOrg);
static WlzErrorNum	WlzContourIsoTet3D(
			  WlzContour *ctr,
			  double *tVal,
			  WlzDVertex3 *tPos,
			  WlzDVertex3 cbOrg);
static WlzErrorNum	WlzContourGrdLink2D(
			  WlzContour *ctr,
			  UBYTE **grdDBuf,
			  int bufWd,
			  WlzIVertex2 org,
			  int lnOff,
			  int lnIdx[],
			  int klP,
			  int *vtxSearchIdx);
static void		WlzContourTestOutVTK(
			  FILE *fP,
			  WlzContour *ctr);

/************************************************************************
* Function:	WlzMakeContour
* Returns:	WlzContour *:		New contour.
* Purpose:	Makes a new contour data structure.
* Global refs:	-
* Parameters:	WlzContourType ctrType:	Contour type 2D or 3D.
		WlzErrorNum *dstErr:	Destination error pointer, may
					be null.
************************************************************************/
WlzContour	*WlzMakeContour(WlzContourType ctrType, WlzErrorNum *dstErr)
{
  WlzContour	*ctr = NULL;
  WlzGMModelType modType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(ctrType)
  {
    case WLZ_CONTOUR_TYPE_2D:
      modType = WLZ_GMMOD_2D;
      break;
    case WLZ_CONTOUR_TYPE_3D:
      modType = WLZ_GMMOD_3D;
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((ctr = AlcCalloc(1, sizeof(WlzContour))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      ctr->model = WlzGMModelNew(modType, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        AlcFree(ctr);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}

/************************************************************************
* Function:	WlzFreeContour
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Free's a WlzContour data structure.
* Global refs:	-
* Parameters:	WlzContour *ctr:	Given contour to free.
************************************************************************/
WlzErrorNum	WlzFreeContour(WlzContour *ctr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctr->type != WLZ_CONTOUR)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if(WlzUnlink(&(ctr->linkcount), &errNum))
    {
      (void )WlzGMModelFree(ctr->model);
      AlcFree((void *)ctr);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourObj
* Returns:	WlzContour:		Contour, or NULL on error.
* Purpose:	Creates a contour (list of connected edges or surface
*		patches) from a Woolz object's values.
*		The source object should either a 2D or 3D domain
*		object with values. This is the top level contour
*		generation function which calls the appropriate
*		function for the given contour type (dimension) and 
*		generation method.
*		The given contour value is taken to be the iso-value
*		for iso-value contours and the minimum gradient
*		threshold value for maximal gradient contours.
*		The contour width parameter is only used for maximal
*		gradient contours where it is used to generate a
*		recursive Deriche filter (see WlzRsvFilter).
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		WlzContourMethod ctrMtd: Contour generation method.
*		double ctrVal:		Contour value.
*		double ctrWth:		Contour filter width.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzContour	*WlzContourObj(WlzObject *srcObj, WlzContourMethod ctrMtd,
			       double ctrVal, double ctrWth,
			       int minNod, int minEdg,
			       WlzErrorNum *dstErr)
{
  WlzContour	*ctr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	switch(ctrMtd)
	{
	  case WLZ_CONTOUR_MTD_ISO:
	    ctr = WlzContourIsoObj2D(srcObj, ctrVal, minNod, minEdg,
	    				&errNum);
	    break;
	  case WLZ_CONTOUR_MTD_GRD:
	    ctr = WlzContourGrdObj2D(srcObj, ctrVal, ctrWth, minNod,
	    				minEdg, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	switch(ctrMtd)
	{
	  case WLZ_CONTOUR_MTD_ISO:
	    ctr = WlzContourIsoObj3D(srcObj, ctrVal, minNod, minEdg,
	    				&errNum);
	    break;
	  case WLZ_CONTOUR_MTD_GRD:
	    ctr = WlzContourGrdObj3D(srcObj, ctrVal, ctrWth, minNod, 
	    				minEdg, &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
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
  return(ctr);
}

/************************************************************************
* Function:	WlzContourIsoObj2D
* Returns:	WlzContour:		Contour, or NULL on error.
* Purpose:	Creates an iso-value contour (list of edges) from a 2D
*		Woolz object's values.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		double isoVal:		Iso-value.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourIsoObj2D(WlzObject *srcObj, double isoVal,
				      int minNod, int minEdg,
				      WlzErrorNum *dstErr)
{
  int		idX,
  		bufSz,
		bufLft,
		bufRgt,
		bufLnIdx,
  		itvLen,
		itvBufWidth,
		prevFirstVtxOfLine,      /* Previous line first vertex index */
		lastLnVtxSearchIdx;
  WlzDomain	srcDom;
  UBYTE		*itvBuf[2] = {NULL, NULL};
  double	*valBuf[2] = {NULL, NULL};
  WlzContour 	*ctr = NULL;
  WlzDVertex2	sqOrg;
  WlzIntervalWSpace srcIWSp;
  WlzGreyWSpace	srcGWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  prevFirstVtxOfLine = 0;   			       /* Make index invalid */
  /* Create contour. */
  ctr = WlzMakeContour(WLZ_CONTOUR_TYPE_2D, &errNum);
  /* Make buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    srcDom = srcObj->domain;
    bufSz = srcDom.i->lastkl - srcDom.i->kol1 + 1;
    itvBufWidth = (bufSz + 7) / 8;    /* No of bytes in interval buffer line */
    if((AlcBit1Calloc(&(itvBuf[0]), bufSz) != ALC_ER_NONE) ||
       (AlcBit1Calloc(&(itvBuf[1]), bufSz) != ALC_ER_NONE) ||
       (AlcDouble1Malloc(&(valBuf[0]), bufSz) != ALC_ER_NONE) ||
       (AlcDouble1Malloc(&(valBuf[1]), bufSz) != ALC_ER_NONE))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Work down through the object. */
  if((errNum = WlzInitGreyScan(srcObj, &srcIWSp,
  			       &srcGWSp)) == WLZ_ERR_NONE)
  {
    bufLnIdx = 0;
    while((errNum == WLZ_ERR_NONE) &&
    	  ((errNum = WlzNextGreyInterval(&srcIWSp)) == WLZ_ERR_NONE))
    {
      if(srcIWSp.nwlpos > 0)
      {
        bufLnIdx = !bufLnIdx;
	WlzValueSetUByte(itvBuf[bufLnIdx], 0, itvBufWidth);
        if(srcIWSp.nwlpos > 1)
	{
	  WlzValueSetUByte(itvBuf[!bufLnIdx], 0, itvBufWidth);
	}
	lastLnVtxSearchIdx = prevFirstVtxOfLine;
	prevFirstVtxOfLine = ctr->model->res.vertex.nxtIdx;
      }
      itvLen = srcIWSp.rgtpos - srcIWSp.lftpos + 1;
      bufLft = srcIWSp.lftpos - srcDom.i->kol1;
      bufRgt = srcIWSp.rgtpos - srcDom.i->kol1;
      WlzBitLnSetItv(itvBuf[bufLnIdx], bufLft, bufRgt, itvLen);
      switch(srcGWSp.pixeltype)
      {
        case WLZ_GREY_INT:
	  WlzValueCopyIntToDouble(valBuf[bufLnIdx] + bufLft,
	  			  srcGWSp.u_grintptr.inp, itvLen);
	  break;
        case WLZ_GREY_SHORT:
	  WlzValueCopyShortToDouble(valBuf[bufLnIdx] + bufLft, 
	  			    srcGWSp.u_grintptr.shp, itvLen);
	  break;
        case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToDouble(valBuf[bufLnIdx] + bufLft, 
	  			    srcGWSp.u_grintptr.ubp, itvLen);
	  break;
        case WLZ_GREY_FLOAT:
	  WlzValueCopyFloatToDouble(valBuf[bufLnIdx] + bufLft, 
	  			    srcGWSp.u_grintptr.flp, itvLen);
	  break;
        case WLZ_GREY_DOUBLE:
	  WlzValueCopyDoubleToDouble(valBuf[bufLnIdx] + bufLft, 
	  			     srcGWSp.u_grintptr.dbp, itvLen);
	  break;
      }
      idX = bufLft;
      sqOrg.vtY = srcIWSp.linpos - 1.0;
      while((errNum == WLZ_ERR_NONE) && (idX < bufRgt))
      {
	if((WLZ_BIT_GET(itvBuf[!bufLnIdx], idX) != 0) &&
	   (WLZ_BIT_GET(itvBuf[!bufLnIdx], (idX + 1)) != 0))
	{
	  sqOrg.vtX = srcIWSp.lftpos + idX;
	  errNum = WlzContourIsoCube2D(ctr, isoVal,
				       valBuf[!bufLnIdx] + idX,
				       valBuf[bufLnIdx] + idX,
				       sqOrg, &lastLnVtxSearchIdx);
	}
	++idX;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  /* Tidy up on error. */
  if((errNum != WLZ_ERR_NONE) && (ctr != NULL))
  {
    (void )WlzFreeContour(ctr);
    ctr = NULL;
  }
  /* Free buffers. */
  for(idX = 0; idX < 2; ++idX)
  {
    if(itvBuf[idX])
    {
      AlcFree(itvBuf[idX]);
    }
    if(valBuf[idX])
    {
      AlcFree(valBuf[idX]);
    }
  }
  /* Set error code. */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}

/************************************************************************
* Function:	WlzContourIsoObj3D
* Returns:	WlzContour:		Contour , or NULL on error.
* Purpose:	Creates an iso-value contour (list of surface patches)
*		from a 3D Woolz object's values.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		double isoVal:		Iso-value.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourIsoObj3D(WlzObject *srcObj, double isoVal,
				      int minNod, int minEdg,
				      WlzErrorNum *dstErr)
{
  int		klIdx,
  		lnIdx,
		pnIdx,
		klCnt,
  		lnCnt,
		pnCnt,
  		bufIdx0,
  		bufIdx1,
		lastKlIn,
		thisKlIn;
  WlzObject	*obj2D = NULL;
  WlzValues	*srcValues;
  WlzValues	dummyValues;
  WlzDomain	dummyDom,
  		srcDom;
  WlzIVertex2	bufSz,
  		bufOrg,
		bufOff;
  WlzIBox2	bBox2D;
  WlzIBox3	bBox3D;
  WlzDVertex3	cbOrg;
  UBYTE		*tUP0,
  		*tUP1,
		*tUP2,
		*tUP3;
  UBYTE		**itvBuf[2] = {NULL, NULL};
  double	**valBuf[2] = {NULL, NULL};
  WlzContour 	*ctr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dummyDom.core = NULL;
  dummyValues.core = NULL;
  if((srcDom = srcObj->domain).core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((srcValues = srcObj->values.vox->values) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  /* Create contour. */
  ctr = WlzMakeContour(WLZ_CONTOUR_TYPE_3D, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    bBox3D = WlzBoundingBox3D(srcObj, &errNum);
  }
  /* Make buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz.vtX = bBox3D.xMax - bBox3D.xMin + 1;
    bufSz.vtY = bBox3D.yMax - bBox3D.yMin + 1;
    bufOrg.vtX = bBox3D.xMin;
    bufOrg.vtY = bBox3D.yMin;
    if((AlcBit2Calloc(&(itvBuf[0]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcBit2Calloc(&(itvBuf[1]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(valBuf[0]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(valBuf[1]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Sweep down through the object using a pair of plane buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    obj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom, dummyValues,
			NULL, NULL, &errNum);
    if((obj2D == NULL) && (errNum == WLZ_ERR_NONE))
    {
      errNum = WLZ_ERR_UNSPECIFIED;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    pnIdx = 0;
    pnCnt = srcObj->domain.p->lastpl - srcObj->domain.p->plane1 + 1;
    while((errNum == WLZ_ERR_NONE) && (pnCnt-- > 0))
    {
      cbOrg.vtZ = bBox3D.zMin + pnIdx;
      bufIdx0 = (pnIdx + 1) % 2;
      bufIdx1 = pnIdx % 2;
      obj2D->domain = *(srcObj->domain.p->domains + pnIdx);
      obj2D->values = *(srcObj->values.vox->values + pnIdx);
      bBox2D = WlzBoundingBox2D(obj2D, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        bufOff.vtX = bBox2D.xMin - bBox3D.xMin;
        bufOff.vtY = bBox2D.yMin - bBox3D.yMin;
        errNum = WlzToArray2D((void ***)&(itvBuf[bufIdx1]), obj2D,
			      bufSz, bufOff, 0, WLZ_GREY_BIT);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzToArray2D((void ***)&(valBuf[bufIdx1]), obj2D,
			      bufSz, bufOff, 0, WLZ_GREY_DOUBLE);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* Compute the intersection of the iso-value plane with each cube
	 * of values. */
	if(pnIdx > 0)
	{
	  klCnt = bBox2D.xMax - bBox2D.xMin; 			  /* NOT + 1 */
	  lnIdx = bufOff.vtY;
	  lnCnt = bBox2D.yMax - bBox2D.yMin; 			  /* NOT + 1 */
	  while(lnIdx < lnCnt)
	  {
            cbOrg.vtY = bBox3D.yMin + lnIdx;
	    tUP0 = *(itvBuf[bufIdx0] + lnIdx);
	    tUP1 = *(itvBuf[bufIdx0] + lnIdx + 1);
	    tUP2 = *(itvBuf[bufIdx1] + lnIdx);
	    tUP3 = *(itvBuf[bufIdx1] + lnIdx + 1);
	    klIdx = bufOff.vtX;
	    lastKlIn = (WLZ_BIT_GET(tUP0, klIdx) != 0) &&
	    	       (WLZ_BIT_GET(tUP1, klIdx) != 0) &&
		       (WLZ_BIT_GET(tUP2, klIdx) != 0) &&
		       (WLZ_BIT_GET(tUP3, klIdx) != 0);
	    while(klIdx < klCnt)
	    {
	      /* Check if cube is within the 3D object's domain. */
	      thisKlIn = (WLZ_BIT_GET(tUP0, klIdx + 1) != 0) &&
	      		 (WLZ_BIT_GET(tUP1, klIdx + 1) != 0) &&
			 (WLZ_BIT_GET(tUP2, klIdx + 1) != 0) &&
			 (WLZ_BIT_GET(tUP3, klIdx + 1) != 0);
	      if(lastKlIn && thisKlIn)
	      {
                cbOrg.vtX = bBox3D.xMin + klIdx;
		errNum = WlzContourIsoCube3D(ctr, isoVal,
				       *(valBuf[bufIdx0] + lnIdx) + klIdx,
				       *(valBuf[bufIdx0] + lnIdx + 1) + klIdx,
				       *(valBuf[bufIdx1] + lnIdx) + klIdx,
				       *(valBuf[bufIdx1] + lnIdx + 1) + klIdx,
				       cbOrg);
	      }
	      lastKlIn = thisKlIn;
	      ++klIdx;
	    }
	    ++lnIdx;
	  }
	}
	++pnIdx;
      }
    }
    obj2D->domain = dummyDom;
    obj2D->values = dummyValues;
    WlzFreeObj(obj2D);
  }
  /* Tidy up on error. */
  if((errNum != WLZ_ERR_NONE) && (ctr != NULL))
  {
    (void )WlzFreeContour(ctr);
    ctr = NULL;
  }
  /* Free buffers. */
  for(pnIdx = 0; pnIdx < 2; ++pnIdx)
  {
    if(itvBuf[pnIdx])
    {
      Alc2Free((void **)itvBuf[pnIdx]);
    }
    if(valBuf[pnIdx])
    {
      Alc2Free((void **)valBuf[pnIdx]);
    }
  }
  /* Set error code. */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}

/************************************************************************
* Function:	WlzContourGrdObj2D
* Returns:	WlzContour:		Contour , or NULL on error.
* Purpose:	Creates an maximal gradient contour (list of edges)
*		from a 2D Woolz object's values.
*		Direction of gradient is encoded as:
*                 +------+------+
*		  |\   2 | 1   /|
*		  |  \   |   /  |
*		  | 3  \ | /  0 |
*                 +------+------+
*		  | 4  / | \  7 |
*		  |  /   |   \  |
*		  |/   5 | 6   \|
*                 +------+------+
*
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		double minGrd:		Minimum modulus of gradient.
*		double ftrPrm:		Filter width parameter.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourGrdObj2D(WlzObject *srcObj,
				      double minGrd, double ftrPrm,
				      int minNod, int minEdg,
				      WlzErrorNum *dstErr)
{
  int		idX,
		dCode,
		iCnt,
		iLen,
		iBufSz,
		lnInc,
		prevFirstVtxOfLine,      /* Previous line first vertex index */
		lastLnVtxSearchIdx;
  UBYTE		doLn;
  double	tD0,
  		grdM0,
  		grdM1,
		grdM2,
		grdX0,
		grdY0,
		grdMLft,
		grdMRgt,
		sqMinGrd;
  WlzGreyP	grdGP,
  		grdXGP,
  		grdYGP;
  double	*tDP0,
  		*tDP1,
		*tDP2;
  WlzObject	*gXObj = NULL,
  		*gYObj = NULL;
  WlzRsvFilter	*ftr = NULL;
  UBYTE		*itvP[4];
  UBYTE		**grdIBuf = NULL,
  		**grdDBuf = NULL;
  double	**grdMBuf = NULL,	            /* Magnitude of gradient */
  		**grdXBuf = NULL,    	    /* Horizontal gradient component */
  		**grdYBuf = NULL;             /* Vertical gradient component */
  WlzContour 	*ctr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		lnIdx[3],
  		klIdx[3];
  WlzIVertex2	bufSz,
  		org,
  		posAbs,
		posRel;
  WlzIntervalWSpace gXIWSp,
  		gYIWSp;
  WlzGreyWSpace	gXGWSp,
  		gYGWSp;
  const UBYTE   dTable[8] = {3, 2, 0, 1, 4, 5, 7, 6};

  prevFirstVtxOfLine = 0;   			       /* Make index invalid */
  /* Create a new contour. */
  ctr = WlzMakeContour(WLZ_CONTOUR_TYPE_2D, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute the partial derivatives. */
    sqMinGrd = minGrd * minGrd;
    if((ftr = WlzRsvFilterMakeFilter(WLZ_RSVFILTER_NAME_DERICHE_1,
				     ftrPrm, &errNum)) != NULL)
    {
      ftr->c *= 4;
      gXObj = WlzRsvFilterObj(srcObj, ftr, WLZ_RSVFILTER_ACTION_X, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	gYObj = WlzRsvFilterObj(srcObj, ftr, WLZ_RSVFILTER_ACTION_Y, &errNum);
      }
      WlzRsvFilterFreeFilter(ftr);
    }
  }
  /* Make scan line buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz.vtY = 3;			  /* Keep four lines in the buffers. */
    bufSz.vtX = srcObj->domain.i->lastkl - srcObj->domain.i->kol1 + 1;
    iBufSz = (bufSz.vtX + 7) / 8;
    if((AlcBit2Calloc(&grdIBuf, bufSz.vtY + 1, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcUnchar2Calloc(&grdDBuf, bufSz.vtY + 1, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Calloc(&grdMBuf, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Calloc(&grdXBuf, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Calloc(&grdYBuf, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Scan through the partial derivative objects, filling the scan line
   * buffer. After each line compute the maximal gradient intersections. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((errNum = WlzInitGreyScan(gXObj, &gXIWSp, &gXGWSp)) == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(gYObj, &gYIWSp, &gYGWSp);
    }
  }
  org.vtX = srcObj->domain.i->kol1;
  org.vtY = srcObj->domain.i->line1;
  while((errNum == WLZ_ERR_NONE) &&
        ((errNum = WlzNextGreyInterval(&gXIWSp)) == WLZ_ERR_NONE) &&
	((errNum = WlzNextGreyInterval(&gYIWSp)) == WLZ_ERR_NONE))
  {
    posAbs.vtX = gXIWSp.lftpos;
    posRel.vtX = posAbs.vtX - org.vtX;
    /* Set up buffers for new line. */
    if((lnInc = gXIWSp.nwlpos) > 0)
    {
      if(lnInc > bufSz.vtY)
      {
        lnInc = bufSz.vtY;
      }
      posAbs.vtY = gXIWSp.linpos - lnInc;
      while(posAbs.vtY < gXIWSp.linpos)
      {
	++(posAbs.vtY);
	posRel.vtY = posAbs.vtY - org.vtY;
	lnIdx[0] = (1 + posRel.vtY) % 3;	   /* Index of previous line */
        lnIdx[1] = (2 + posRel.vtY) % 3;            /* Index of current line */
        lnIdx[2] = (3 + posRel.vtY) % 3; /* Index of most recent (next) line */
	grdGP.dbp = *(grdMBuf + lnIdx[2]);
	grdXGP.dbp = *(grdXBuf + lnIdx[2]);
	grdYGP.dbp = *(grdYBuf + lnIdx[2]);
	/* Clear buffers. */
	WlzValueSetUByte(*(grdIBuf + lnIdx[2]), 0, iBufSz);
	WlzValueSetUByte(*(grdDBuf + lnIdx[2]), 0, bufSz.vtX);
	WlzValueSetDouble(grdGP.dbp, 0.0, bufSz.vtX);
	WlzValueSetDouble(grdXGP.dbp, 0.0, bufSz.vtX);
	WlzValueSetDouble(grdYGP.dbp, 0.0, bufSz.vtX);
      }
    }
    iLen = gXIWSp.rgtpos - gXIWSp.lftpos + 1;
    /* Add interval to interval mask buffer. */
    WlzBitLnSetItv(*(grdIBuf + lnIdx[2]), posRel.vtX, posRel.vtX + iLen - 1,
    		   bufSz.vtX);
    /* Copy this intervals of gradient data into buffers. */
    WlzValueCopyGreyToGrey(grdXGP, posRel.vtX,
    			   WLZ_GREY_DOUBLE, gXGWSp.u_grintptr, 0,
			   gXGWSp.pixeltype, iLen);
    WlzValueCopyGreyToGrey(grdYGP, posRel.vtX,
    			   WLZ_GREY_DOUBLE, gYGWSp.u_grintptr, 0,
			   gXGWSp.pixeltype, iLen);
    /* Compute the gradient magnitudes in this interval. */
    iCnt = iLen;
    tDP0 = grdXGP.dbp + posRel.vtX;
    tDP1 = grdYGP.dbp + posRel.vtX;
    tDP2 = grdGP.dbp + posRel.vtX;
    while(iCnt-- > 0)
    {
      tD0 = (*tDP0 * *tDP0) + (*tDP1 * *tDP1);
      ++tDP0;
      ++tDP1;
      *tDP2++ = (tD0 > DBL_EPSILON)? sqrt(tD0): 0.0;
    }
    /* Process the lines in the buffer if this is the last interval of
     * the current line. */
    if(gXIWSp.intrmn == 0)
    {
      /* Demand that a pixel has eight neighbours so build mask to test
       * three lines at a time. */
      doLn = 0;
      itvP[0] = *(grdIBuf + lnIdx[0]);
      itvP[1] = *(grdIBuf + lnIdx[1]);
      itvP[2] = *(grdIBuf + lnIdx[2]);
      itvP[3] = *(grdIBuf + 3); /* Used to store mask for all recent lines. */
      for(idX = 0; idX < iBufSz; ++idX)
      {
	doLn |= *(itvP[3])++ = *(itvP[0])++ & *(itvP[1])++ & *(itvP[2])++;
      }
      if(doLn)
      {
	klIdx[0] = 0;
	klIdx[1] = 1;
	klIdx[2] = 2;
        itvP[3] = *(grdIBuf + 3);
	lastLnVtxSearchIdx = prevFirstVtxOfLine;
	prevFirstVtxOfLine = ctr->model->res.vertex.nxtIdx;
	while((errNum == WLZ_ERR_NONE) && (klIdx[2] < bufSz.vtX))
	{
	  if(WLZ_BIT_GET(itvP[3], klIdx[0]) &&
	     WLZ_BIT_GET(itvP[3], klIdx[1]) &&
	     WLZ_BIT_GET(itvP[3], klIdx[2]))
	  {
	    /* Check if gradient magnitude at the centre of the neighbourhood
	     * is above the minimum gradient threshold. */
	    if((grdM0 = *(*(grdMBuf + lnIdx[1]) + klIdx[1])) >= minGrd)
	    {
	      /* Compute and classify the direction of gradient at centre of
	       * neighbourhood. */
	      grdX0 = *(*(grdXBuf + lnIdx[1]) + klIdx[1]);
	      grdY0 = *(*(grdYBuf + lnIdx[1]) + klIdx[1]);
	      dCode = dTable[((grdY0 >= 0.0) << 2) |
	                     ((grdX0 >= 0.0) << 1) |
			     ((grdY0 * grdY0) >= (grdX0 * grdX0))];
	      /* Interpolate gradient magnitudes at the edges of the
	       * neighbourhood. */
	      switch(dCode)
	      {
		case 0: /* \" */
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[0]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[0]);
		  grdMLft = (((grdM0 - grdM1) * grdX0) -
		             ((grdM1 - grdM2) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[2]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[2]);
		  grdMRgt = (((grdM0 - grdM1) * grdX0) -
		             ((grdM1 - grdM2) * grdY0)) / grdM0;
		  break;
		case 1: /* |\ */
		  grdM1 = *(*(grdMBuf + lnIdx[2]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[0]);
		  grdMLft = (((grdM1 - grdM2) * grdX0) -
			     ((grdM0 - grdM1) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[0]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[2]);
		  grdMRgt = (((grdM1 - grdM2) * grdX0) -
			     ((grdM0 - grdM1) * grdY0)) / grdM0;
		  break;
		case 2: /* /| */
		  grdM1 = *(*(grdMBuf + lnIdx[2]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[2]);
		  grdMLft = (((grdM2 - grdM1) * grdX0) -
			     ((grdM0 - grdM1) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[0]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[0]);
		  grdMRgt = (((grdM2 - grdM1) * grdX0) -
			     ((grdM0 - grdM1) * grdY0)) / grdM0;
		  break;
		case 3: /* "/ */
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[2]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[2]);
		  grdMLft = (((grdM1 - grdM0) * grdX0) -
			     ((grdM1 - grdM2) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[0]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[0]);
		  grdMRgt = (((grdM1 - grdM0) * grdX0) -
			     ((grdM1 - grdM2) * grdY0)) / grdM0;
		  break;
		case 4: /* _\ */
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[2]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[2]);
		  grdMLft = (((grdM1 - grdM0) * grdX0) -
			     ((grdM2 - grdM1) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[0]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[0]);
		  grdMRgt = (((grdM1 - grdM0) * grdX0) -
			     ((grdM2 - grdM1) * grdY0)) / grdM0;
		  break;
		case 5: /* \| */
		  grdM1 = *(*(grdMBuf + lnIdx[0]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[2]);
		  grdMLft = (((grdM2 - grdM1) * grdX0) -
			     ((grdM1 - grdM0) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[2]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[0]);
		  grdMRgt = (((grdM2 - grdM1) * grdX0) -
			     ((grdM1 - grdM0) * grdY0)) / grdM0;
		  break;
		case 6: /* |/ */
		  grdM1 = *(*(grdMBuf + lnIdx[0]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[0]);
		  grdMLft = (((grdM1 - grdM2) * grdX0) -
			     ((grdM1 - grdM0) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[2]) + klIdx[1]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[2]);
		  grdMRgt = (((grdM1 - grdM2) * grdX0) -
			     ((grdM1 - grdM0) * grdY0)) / grdM0;
		  break;
		case 7: /* /_ */
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[0]);
		  grdM2 = *(*(grdMBuf + lnIdx[0]) + klIdx[0]);
		  grdMLft = (((grdM0 - grdM1) * grdX0) -
			     ((grdM2 - grdM1) * grdY0)) / grdM0;
		  grdM1 = *(*(grdMBuf + lnIdx[1]) + klIdx[2]);
		  grdM2 = *(*(grdMBuf + lnIdx[2]) + klIdx[2]);
		  grdMRgt = (((grdM0 - grdM1) * grdX0) -
			     ((grdM2 - grdM1) * grdY0)) / grdM0;
		  break;
	      }
	      if((grdMLft > 0.0) && (grdMRgt > 0.0))
	      {
		/* Set element of maximal pixel gradient direction buffer. */
		*(*(grdDBuf + lnIdx[1]) + klIdx[2]) = 0x80 | dCode;
	      }
	    }
	    if(*(*(grdDBuf + lnIdx[1]) + klIdx[1]))
	    {
	      /* Generate and link edge segments. */
	      errNum = WlzContourGrdLink2D(ctr, grdDBuf,
	      				   bufSz.vtX, org, posRel.vtY - 2,
					   lnIdx, klIdx[0],
					   &lastLnVtxSearchIdx);
	    }
	  }
	  klIdx[0] = klIdx[1];
	  klIdx[1] = klIdx[2];
	  ++klIdx[2];
	}
      }
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  /* Tidy up on error. */
  if((errNum != WLZ_ERR_NONE) && (ctr != NULL))
  {
    (void )WlzFreeContour(ctr);
    ctr = NULL;
  }
  /* Free buffers. */
  if(grdIBuf)
  {
    Alc2Free((void **)grdIBuf);
  }
  if(grdDBuf)
  {
    Alc2Free((void **)grdDBuf);
  }
  if(grdMBuf)
  {
    Alc2Free((void **)grdMBuf);
  }
  if(grdXBuf)
  {
    Alc2Free((void **)grdXBuf);
  }
  if(grdYBuf)
  {
    Alc2Free((void **)grdYBuf);
  }
  if(gXObj)
  {
    WlzFreeObj(gXObj);
  }
  if(gYObj)
  {
    WlzFreeObj(gYObj);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}

/************************************************************************
* Function:	WlzContourGrdObj3D
* Returns:	WlzContour:		Contour or NULL on error.
* Purpose:	Creates an maximal gradient contour (list of surface 
*		patches) from a 3D Woolz object's values.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		double minGrd:		Minimum modulus of gradient.
*		double ftrPrm:		Filter width parameter.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourGrdObj3D(WlzObject *srcObj,
				      double minGrd, double ftrPrm,
				      int minNod, int minEdg,
				      WlzErrorNum *dstErr)
{
  WlzContour 	*ctr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* TODO */

  /* Tidy up on error. */
  if((errNum != WLZ_ERR_NONE) && (ctr != NULL))
  {
    (void )WlzFreeContour(ctr);
    ctr = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}

/************************************************************************
* Function:	WlzContourGrdLink2D
* Returns:	WlzErrorNum		Woolz error code.
* Purpose:	Computes edge segment(s) linking the pixel at the
*		centre of a 3x3 neighbourhood to 4 and 8 connected
*		pixels in the neighbourhood.
*		First search for 4 connected neighbours of the
*		central pixel, then search for any 8 connected neighbours
*		which are themselves not 4 connected to any of the
*		4 connected neighbours already found.
*		Indicies of the neighbours of C are:
*		    +    +    +    +
*		    +----+----+----+    lines
*		    + 0  + C  + 2  +    ^
*		    +----+----+----+    |
*		    + 3  + 1  + 4  +    |
*		    +----+----+----+    +---> cols
* Global refs:	-
* Parameters:	WlzContour *ctr: 	Contour being built.
*		UBYTE **grdDBuf:	Buffers containing gradient
*					direction codes for maximal
*					gradient pixels.
*		int bufWd:		Width of the buffers.
*		WlzIVertex2 org:	Origin of the object.
*		int lnOff:		Offset from origin to central
*					pixel line.
*		int lnIdx[]:		Line indicies.
*		int klOff:		Column offset of first column.
*		int *vtxSearchIdx:	Last line vertex search index.
************************************************************************/
static WlzErrorNum WlzContourGrdLink2D(WlzContour *ctr, UBYTE **grdDBuf,
				       int bufWd, WlzIVertex2 org,
				       int lnOff, int lnIdx[], int klOff,
			  	       int *vtxSearchIdx)
{
  int		idN,
		idM,
		conCnt;
  int		con[5],
     		nbr[4];
  WlzDVertex2	nbrOrg,
  		matchDist;
  WlzDVertex2	seg[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	offConKl[5] =
  {
    0, 1, 2, 			     /* 4 connected neighbour column offsets */
    0, 2			     /* 8 connected neighbour column offsets */
  },
  		offConLn[5] =
  {
    1, 0, 1, 			       /* 4 connected neighbour line offsets */
    0, 0			       /* 8 connected neighbour line offsets */
  };

  conCnt = 0;
  matchDist.vtX = 2.01;
  matchDist.vtY = 2.01;
#ifdef WLZ_CONTOUR_DEBUG
  fprintf(stderr, "255 0 0 setrgbcolor\n");
  seg[0].vtX = org.vtX + klOff + 1.0 - 0.1;
  seg[0].vtY = org.vtY + lnOff + 1.0 - 0.1;
  seg[1].vtX = org.vtX + klOff + 1.0 + 0.1;
  seg[1].vtY = org.vtY + lnOff + 1.0 + 0.1;
  WlzContourTestOutPSLn2D(stderr, seg[0], seg[1]);
  seg[0].vtX = org.vtX + klOff + 1.0 - 0.1;
  seg[0].vtY = org.vtY + lnOff + 1.0 + 0.1;
  seg[1].vtX = org.vtX + klOff + 1.0 + 0.1;
  seg[1].vtY = org.vtY + lnOff + 1.0 - 0.1;
  WlzContourTestOutPSLn2D(stderr, seg[0], seg[1]);
  fprintf(stderr, "0 0 0 setrgbcolor\n");
#endif /* WLZ_CONTOUR_DEBUG */
  /* Find 4 connected neighbours. */
  for(idN = 0; idN < 3;  ++idN)
  {
    if(*(*(grdDBuf + lnIdx[offConLn[idN]]) + klOff + offConKl[idN]))
    {
      con[idN] = 1;
      nbr[conCnt++] = idN;
    }
    else
    {
      con[idN] = 0;
    }
  }
  /* Find 8 connected neighbours. */
  if((con[3] = !(con[0] || con[1]) && *(*(grdDBuf + lnIdx[0]) + klOff + 0)))
  {
    nbr[conCnt++] = 3;
  }
  if((con[4] = !(con[1] || con[2]) && *(*(grdDBuf + lnIdx[0]) + klOff + 2)))
  {
    nbr[conCnt++] = 4;
  }
  if(conCnt > 0)
  {
    /* Compute segment end points and add them to the contour workspace. */
    nbrOrg.vtX = org.vtX + klOff;
    nbrOrg.vtY = org.vtY + lnOff;
    seg[0].vtX = org.vtX + klOff + 1.0;
    seg[0].vtY = org.vtY + lnOff + 1.0;
    for(idN = 0; idN < conCnt; ++idN)
    {
      idM = nbr[idN];
      seg[1].vtX = org.vtX + klOff + offConKl[idM];
      seg[1].vtY = org.vtY + lnOff + offConLn[idM];
#ifdef WLZ_CONTOUR_DEBUG
      WlzContourTestOutPSLn2D(stderr, seg[0], seg[1]);
#endif /* WLZ_CONTOUR_DEBUG */
      errNum = WlzGMModelConstructSimplex2D(ctr->model,
      					    matchDist, vtxSearchIdx, seg);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIsoCube2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Computes iso-value contour segments for the given
*		square, creates edge and node elements and inserts
*		them into the contour data structure.
*		From the 4 data values at the corners (1,..4) of
*		the square and a linearly interpolated value at it's
*		centre (0), the square is divided into 4 triangles
*		(0,...3).
*
*                 4-------------------3
*                 | \               / |
*                 |   \     2     /   |
*                 |     \       /     |
*                 |       \   /       |
*                 |   3     0    1 <------ triangle index
*                 |       /   \       |
*                 |     /       \     |
*                 |   /     0     \   |
*                 | /               \ |
*                 1-------------------2 <- square node index
*
*                           0 <----------- triangle node index
*                         /   \        
*                       /       \      
*                     /           \    
*                   /               \  
*                 1-------------------2
*       
* 		Each triangle is classified by generating a above/on/
*		below code for each of it's verticies (above = 2,
*		on 1, below 0) and using these codes to index a
*		look up table.
*		To avoid tests for duplicate edge segments in other
*		parts of the code it is done here.
*		This function is based on Paul Bourke's CONREC.F
*		contouring subroutine, Paul Bourke. CONREC A Contouring
*		Subroutine. BYTE, July 1997.
* Global refs:	-
* Parameters:	WlzContour *ctr: 	Contour being built.
*		double isoVal:		Iso-value to use.
*		double *vLn0:		Ptr to 2 data values at
*					y = yPos and x = xPos,
*					xpos + 1.
*		double *vLn1:		Ptr to 2 data values at
*					y = yPos + 1 and x = xPos,
*					xpos + 1..
*		WlzDVertex2 sqOrg:	The square's origin.
*		int *vtxSearchIdx:	Last line vertex search index.
************************************************************************/
static WlzErrorNum WlzContourIsoCube2D(WlzContour *ctr,
				       double isoVal,
				       double *vLn0, double *vLn1,
				       WlzDVertex2 sqOrg, int *vtxSearchIdx)
{
  int		idT,
		tN1,
  		tN2,
		dupIsn = 0;
  double	tD0;
  WlzContourTriIsn2D iCode;
  WlzDVertex2	tVx0,
  		tVx1,
		matchDist;
  int		lev[3];  /* Rel. triangle node level: above 2, on 1, below 0 */
  double	tV[3],			       /* Rel. triangle node values. */
  		sV[5];	 			 /* Rel. square node values. */
  WlzDVertex2	tIsn[2];		/* Triangle intersection coordinates */
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	tN0 = 0;
  const double	tol = 1.0e-06;
  const WlzContourTriIsn2D iCodeTab[3][3][3] =	  /* Intersection code table */
  {
    {
      {
	WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_S10S02
      },
      {
        WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_N1N0,
	WLZ_CONTOUR_TIC2D_N1S02
      },
      {
        WLZ_CONTOUR_TIC2D_S21S10,
	WLZ_CONTOUR_TIC2D_N0S21,
	WLZ_CONTOUR_TIC2D_S02S21
      }
    },
    {
      {
        WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_N0N2,
	WLZ_CONTOUR_TIC2D_N2S10
      },
      {
        WLZ_CONTOUR_TIC2D_N2N1,
	WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_N2N1
      },
      {
        WLZ_CONTOUR_TIC2D_N2S10,
	WLZ_CONTOUR_TIC2D_N0N2,
	WLZ_CONTOUR_TIC2D_NONE
	}
    },
    {
      {
        WLZ_CONTOUR_TIC2D_S02S21,
	WLZ_CONTOUR_TIC2D_N0S21,
	WLZ_CONTOUR_TIC2D_S21S10
      },
      {
        WLZ_CONTOUR_TIC2D_N1S02,
	WLZ_CONTOUR_TIC2D_N1N0,
	WLZ_CONTOUR_TIC2D_NONE
      },
      {
        WLZ_CONTOUR_TIC2D_S10S02,
	WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_NONE
      }
    }
  };
  const WlzDVertex2 sNOffTab[5] =	   /* Square node horizontal offsets */
  {
    {0.5, 0.5}, {0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0} /* {vtY, vtX} */
  };
  const WlzDVertex2 dupIsnTab[2][4] =    /* Duplicate intersections to avoid */
  {
    {
      {0.0, 0.0}, {0.5, 0.5}, {0.5, 0.5}, {1.0, 0.0}	       /* {vtY, vtX} */
    },
    {
      {0.5, 0.5}, {0.0, 1.0}, {1.0, 1.0}, {0.5, 0.5}	       /* {vtY, vtX} */
    }
  };

  matchDist.vtX = 1.01; /* TODO check match distances */
  matchDist.vtY = 1.01;
  sV[1] = *vLn0 - isoVal; sV[2] = *(vLn0 + 1) - isoVal;
  sV[3] = *(vLn1 + 1) - isoVal; sV[4] = *vLn1 - isoVal;
#ifdef WLZ_CONTOUR_DEBUG
  (void )fprintf(stderr,
		 "%% {%f %f} %f %f %f %f\n",
		 sqOrg.vtX, sqOrg.vtY, sV[1], sV[2], sV[3], sV[4]);
#endif /* WLZ_CONTOUR_DEBUG */
  /* Reject all squares not intersected by the iso-value. */
  if(((sV[1] < 0.0) || (sV[2] < 0.0) || (sV[3] < 0.0) || (sV[4] < 0.0)) &&
     ((sV[1] >= 0.0) || (sV[2] >= 0.0) || (sV[3] >= 0.0) || (sV[4] >= 0.0)))
  {
    /* Interpolate for centre of square */
    sV[0] = (sV[1] + sV[2] + sV[3] + sV[4]) / 4.0;
    tV[0] = sV[0];
    lev[0] = (tV[0] >= DBL_EPSILON) + (tV[0] > -(DBL_EPSILON));
    /* For each triangle: Classify by node levels and compute intersection. */
    idT = 0;
    while((errNum == WLZ_ERR_NONE) && (idT < 4))
    {
      if(idT == 3)
      {
        tN1 = 4;
	tN2 = 1;
      }
      else
      {
        tN1 = idT + 1;
	tN2 = idT + 2;
      }
      tV[1] = sV[tN1];
      tV[2] = sV[tN2];
      lev[1] = (tV[1] >= DBL_EPSILON) + (tV[1] > -(DBL_EPSILON));
      lev[2] = (tV[2] >= DBL_EPSILON) + (tV[2] > -(DBL_EPSILON));
      iCode = iCodeTab[lev[2]][lev[1]][lev[0]];
      if(iCode)
      {
	switch(iCode)
	{
	  case WLZ_CONTOUR_TIC2D_N1N0:
	    tIsn[0] = sNOffTab[tN1];
	    tIsn[1] = sNOffTab[tN0];
	    break;
	  case WLZ_CONTOUR_TIC2D_N0N2:
	    tIsn[0] = sNOffTab[tN0];
	    tIsn[1] = sNOffTab[tN2];
	    break;
	  case WLZ_CONTOUR_TIC2D_N2N1:
	    tIsn[0] = sNOffTab[tN2];
	    tIsn[1] = sNOffTab[tN1];
	    break;
	  case WLZ_CONTOUR_TIC2D_N1S02:
	    tIsn[0] = sNOffTab[tN1];
	    tIsn[1] = WlzContourItpTriSide(tV[0], tV[2],
	    				   sNOffTab[tN0], sNOffTab[tN2]);
	    break;
	  case WLZ_CONTOUR_TIC2D_N0S21:
	    tIsn[0] = sNOffTab[tN0];
	    tIsn[1] = WlzContourItpTriSide(tV[2], tV[1],
	    				   sNOffTab[tN2], sNOffTab[tN1]);
	    break;
	  case WLZ_CONTOUR_TIC2D_N2S10:
	    tIsn[0] = sNOffTab[tN2];
	    tIsn[1] = WlzContourItpTriSide(tV[1], tV[0],
	    				   sNOffTab[tN1], sNOffTab[tN0]);
	    break;
	  case WLZ_CONTOUR_TIC2D_S10S02:
	    tIsn[0] = WlzContourItpTriSide(tV[1], tV[0],
	    				   sNOffTab[tN1], sNOffTab[tN0]);
	    tIsn[1] = WlzContourItpTriSide(tV[0], tV[2],
	    				   sNOffTab[tN0], sNOffTab[tN2]);
	    break;
	  case WLZ_CONTOUR_TIC2D_S02S21:
	    tIsn[0] = WlzContourItpTriSide(tV[0], tV[2],
	    				   sNOffTab[tN0], sNOffTab[tN2]);
	    tIsn[1] = WlzContourItpTriSide(tV[2], tV[1],
	    				   sNOffTab[tN2], sNOffTab[tN1]);
	    break;
	  case WLZ_CONTOUR_TIC2D_S21S10:
	    tIsn[0] = WlzContourItpTriSide(tV[2], tV[1],
	    				   sNOffTab[tN2], sNOffTab[tN1]);
	    tIsn[1] = WlzContourItpTriSide(tV[1], tV[0],
	    				   sNOffTab[tN1], sNOffTab[tN0]);
	    break;
	  default:
	    break;
	}
	if(tIsn[0].vtX > tIsn[1].vtX)
	{
	  tVx0 = tIsn[0];
	  tIsn[0] = tIsn[1];
	  tIsn[1] = tVx0;
	}
	/* Avoid duplicate edge segments along the triangle diagonals! */
	dupIsn = WlzGeomVtxEqual2D(tIsn[0], dupIsnTab[0][idT], tol) &&
		 WlzGeomVtxEqual2D(tIsn[1], dupIsnTab[1][idT], tol);
	if(!dupIsn)
	{
#ifdef WLZ_CONTOUR_DEBUG
	  (void )fprintf(stderr,
		         "%% T%d I%d {%g %g %g} {%g %g} {%g %g}\n",
		         idT, iCode,
		         tV[0], tV[1], tV[2],
		         tIsn[0].vtX, tIsn[0].vtY,
		         tIsn[1].vtX, tIsn[1].vtY);
#endif /* WLZ_CONTOUR_DEBUG */
	  tIsn[0].vtX += sqOrg.vtX;
	  tIsn[0].vtY += sqOrg.vtY;
	  tIsn[1].vtX += sqOrg.vtX;
	  tIsn[1].vtY += sqOrg.vtY;
#ifdef WLZ_CONTOUR_DEBUG
  	  WlzContourTestOutPSLn2D(stderr, tIsn[0], tIsn[1]);
#endif /* WLZ_CONTOUR_DEBUG */
	  errNum = WlzGMModelConstructSimplex2D(ctr->model,
						matchDist, vtxSearchIdx, tIsn);
	}
      }
      ++idT;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIsoCube3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose: 	Checks to see if all above the cube's vertex values are
*		either above or below the iso-value. If they are then
*		there's no intersection between the iso-surface and the
*		cube. If there is a possible intersection then the cube
*		is split into 24 tetrahedra and each tetrahedra is
*		passed to WlzContourIsoTet3D().
*
*		From the 8 data values at the verticies @[1-8] of
*		the cube and interpolated values at it's centre *0,
*		and the centres of the cubes faces (9,...,14) the cube
*		is divided into 24 tetrahedra (0,...,23).
*                                                                      
*                          @8----------------@7                        
*                         /|                /|                        
*                        / |               / |                         
*                       /  |              /  |                         
*                      /   |    +14      /   |                         
*                     /    |        +12 /    |                         
*                    /     |           /     |                         
*                   /      |          /      |                         
*                  /       |         /       |                         
*                 @5----------------@6       |                         
*                 |    +13 |   *0   |    +11 |                         
*                 |        @4-------|--------@3                                 
*                 |       /         |       /                                 
*                 |      /          |      /                                  
*                 |     /           |     /                                  
*                 |    /    +10 +9  |    /                                  
*                 |   /             |   /                                  
*                 |  /              |  /                                  
*                 | /               | /                                  
*                 |/                |/                                  
*                 @1----------------@2                                 
*
*		The tetrahedra are are assigned the following indicies:
*
*		  tetradedron index  	cube vertex indicies
*		   0       		0,  1, 2,  9
*		   1       		0,  2, 3,  9
*		   2       		0,  3, 4,  9
*		   3       		0,  4, 1,  9
*		   4        		0,  1, 2, 10
*		   5        		0,  2, 6, 10
*		   6        		0,  6, 5, 10
*		   7        		0,  5, 1, 10
*		   8        		0,  2, 3, 11
*		   9        		0,  3, 7, 11
*		  10        		0,  7, 6, 11
*		  11       		0,  6, 2, 11
*		  12       		0,  3, 4, 12
*		  13       		0,  4, 8, 12
*		  14       		0,  8, 7, 12
*		  15       		0,  7, 3, 12
*		  16       		0,  4, 1, 13
*		  17       		0,  1, 5, 13
*		  18       		0,  5, 8, 13
*		  19       		0,  8, 4, 13
*		  20       		0,  5, 6, 14
*		  21       		0,  6, 7, 14
*		  22       		0,  7, 8, 14
*		  23       		0,  8, 5, 14
*
*		The values at these vertices are passed onto
*		WlzContourIsoTet3D() using the following indicies:
*		
*                          +3
*                         /|\                                              
*                        / | \                                             
*                       /  |  \                                            
*                      /   |   \                                           
*                     /    |    \                                          
*                    /     |     \                                         
*                   /      *0     \                                        
*                  /    ,/    \,   \                                       
*                 / / '          `\ \                                      
*                 @1-----------------@2
*
*		The tetrahedra verticies are assigned using a look up table
*		in which tetrahedron vertex 0 is always vertex 0 of the cube,
*		tetrahedron vertex 3 is always the cube face vertex
*		(9, 10, 11, 12, 13, 14) and the tetrahedron verticies
*		1 and 2 are cube edge verticies.
* Global refs:	-
* Parameters:	WlzContour *ctr:	Contour being built.
*		double isoVal:		Iso-value to use.
*		double *vPn0Ln0:	Ptr to 2 data values at
*					z = zPos, y = yPos and
*					x = xPos, xpos + 1.
*		double *vPn0Ln1:	Ptr to 2 data values at
*					z = zPos, y = yPos + 1 and
*					x = xPos, xpos + 1.
*		double *vPn1Ln0:	Ptr to 2 data values at
*					z = zPos + 1, y = yPos and
*					x = xPos, xpos + 1.
*		double *vPn1Ln1:	Ptr to 2 data values at
*					z = zPos + 1, y = yPos + 1 and
*					x = xPos, xpos + 1.
*		WlzDVertex3 cbOrg:	The cube's origin.
************************************************************************/
static WlzErrorNum WlzContourIsoCube3D(WlzContour *ctr,
				double isoVal,
				double *vPn0Ln0, double *vPn0Ln1,
				double *vPn1Ln0, double *vPn1Ln1,
				WlzDVertex3 cbOrg)
{
  int		tIdx,
		tI1,
		tI2,
		tI3,
  		intersect;
  double 	cVal[15],	  /* Cube's values relative to the iso-value */
		sVal[8],        /* Temporary 'side' values for interpolation */
  		tVal[4];			       /* Tetrahedron values */
  WlzDVertex3	tPos[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	tVxLUT[24][4] =  /* Tetrahedron to cube vertex look up table */
  {
    {0, 1, 2,  9}, {0, 2, 3,  9}, {0, 3, 4,  9}, {0, 4, 1,  9},
    {0, 1, 2, 10}, {0, 2, 6, 10}, {0, 6, 5, 10}, {0, 5, 1, 10},
    {0, 2, 3, 11}, {0, 3, 7, 11}, {0, 7, 6, 11}, {0, 6, 2, 11},
    {0, 3, 4, 12}, {0, 4, 8, 12}, {0, 8, 7, 12}, {0, 7, 3, 12},
    {0, 4, 1, 13}, {0, 1, 5, 13}, {0, 5, 8, 13}, {0, 8, 4, 13},
    {0, 5, 6, 14}, {0, 6, 7, 14}, {0, 7, 8, 14}, {0, 8, 5, 14}
  };
  const WlzDVertex3 cPos[15] =	 /* Cube positions, order is {vtX, vtY, vtZ} */
  {
    {0.5, 0.5, 0.5},
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0},
    {0.5, 0.5, 0.0},
    {0.5, 0.0, 0.5}, {1.0, 0.5, 0.5}, {0.5, 1.0, 0.5}, {0.0, 0.5, 0.5},
    {0.5, 0.5, 0.1}
  };

  cVal[ 1] = *(vPn0Ln0 + 0) - isoVal;
  cVal[ 2] = *(vPn0Ln0 + 1) - isoVal;
  cVal[ 3] = *(vPn0Ln1 + 1) - isoVal;
  cVal[ 4] = *(vPn0Ln1 + 0) - isoVal;
  cVal[ 5] = *(vPn1Ln0 + 0) - isoVal;
  cVal[ 6] = *(vPn1Ln0 + 1) - isoVal;
  cVal[ 7] = *(vPn1Ln1 + 1) - isoVal;
  cVal[ 8] = *(vPn1Ln1 + 0) - isoVal;
  intersect = !(((cVal[ 1] < 0.0) && (cVal[ 2] < 0.0) &&
		(cVal[ 3] < 0.0) && (cVal[ 4] < 0.0) &&
		(cVal[ 5] < 0.0) && (cVal[ 6] < 0.0) &&
		(cVal[ 7] < 0.0) && (cVal[ 8] < 0.0)) || 
	       ((cVal[ 1] > 0.0) && (cVal[ 2] > 0.0) &&
		(cVal[ 3] > 0.0) && (cVal[ 4] > 0.0) &&
		(cVal[ 5] > 0.0) && (cVal[ 6] > 0.0) &&
		(cVal[ 7] > 0.0) && (cVal[ 8] > 0.0)));
  if(intersect)
  {
    sVal[0] = cVal[ 1] + cVal[ 2];
    sVal[1] = cVal[ 2] + cVal[ 3];
    sVal[2] = cVal[ 3] + cVal[ 4];
    sVal[3] = cVal[ 4] + cVal[ 1];
    sVal[4] = cVal[ 5] + cVal[ 6];
    sVal[5] = cVal[ 6] + cVal[ 7];
    sVal[6] = cVal[ 7] + cVal[ 8];
    sVal[7] = cVal[ 8] + cVal[ 5];
    cVal[ 9] = (sVal[0] + sVal[2]) * 0.25;
    cVal[10] = (sVal[0] + sVal[4]) * 0.25;
    cVal[11] = (sVal[1] + sVal[5]) * 0.25;
    cVal[12] = (sVal[2] + sVal[6]) * 0.25;
    cVal[13] = (sVal[3] + sVal[7]) * 0.25;
    cVal[14] = (sVal[4] + sVal[6]) * 0.25;
    cVal[ 0] = (cVal[ 9] + cVal[14]) * 0.5;
    tVal[0] = cVal[ 0];
    tPos[0] = cPos[ 0];
    for(tIdx = 0; tIdx < 24; ++tIdx)
    {
      tI1 = tVxLUT[tIdx][1];
      tI2 = tVxLUT[tIdx][2];
      tI3 = tVxLUT[tIdx][3];
      tVal[1] = cVal[tI1]; tVal[2] = cVal[tI2]; tVal[3] = cVal[tI3];
      intersect = !(((tVal[0] < 0.0) && (tVal[1] < 0.0) &&
		     (tVal[2] < 0.0) && (tVal[3] < 0.0)) ||
		    ((tVal[0] > 0.0) && (tVal[1] > 0.0) &&
		     (tVal[2] > 0.0) && (tVal[3] > 0.0)));
      if(intersect)
      {
	tPos[1] = cPos[tI1]; tPos[2] = cPos[tI2]; tPos[3] = cPos[tI3];
        errNum = WlzContourIsoTet3D(ctr, tVal, tPos, cbOrg);
      }
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIsoTet3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose: 	Computes the intersection of the given tetrahedron
*		with the isovalue surface. Creates facet, edge and
*		node elements and inserts them into the contour data
*		structure.
* Global refs:	-
* Parameters:	WlzContourWSp *ctr:	Contour being built.
*		double *tVal:		Values wrt the iso-value at the
*					verticies of the tetrahedron.
*		WlzDVertex3 *tPos:	Positions of the tetrahedron
*					verticies wrt the cube's origin.
*		WlzDVertex3 cbOrg:	The cube's origin.
************************************************************************/
static WlzErrorNum WlzContourIsoTet3D(WlzContour *ctr,
				      double *tVal,
				      WlzDVertex3 *tPos,
				      WlzDVertex3 cbOrg)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* TODO */
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIntpTriSide
* Returns:	WlzDVertex2:		Position of intersection with
*					side.
* Purpose:	Calculates the position of the intersection of a
*		line with the zero height plane.
* Global refs:	-
* Parameters:	double ht0:		Height at line origin.
*		double ht1:		Height at line destination.
*		WlzDVertex2 org:	Line origin.
*		WlzDVertex2 dst:	Line destination.
************************************************************************/
static WlzDVertex2 WlzContourItpTriSide(double ht0, double ht1,
				        WlzDVertex2 org, WlzDVertex2 dst)
{
  double	tD0;
  WlzDVertex2	itp;

  tD0 = ht1 - ht0;
  itp.vtX = ((ht1 * org.vtX)  - (ht0 * dst.vtX)) / tD0;
  itp.vtY = ((ht1 * org.vtY)  - (ht0 * dst.vtY)) / tD0;
  return(itp);
}

/************************************************************************
* Function:	WlzContourTestOutPSLn2D
* Returns:	void
* Purpose:	Prints out line segment for testing.
* Global refs:	-
* Parameters:	FILE *fP:		Output file.
		WlzDVertex2 org:	Start of line segment.
*		WlzDVertex2 dst:	End of line segment.
************************************************************************/
static void 	WlzContourTestOutPSLn2D(FILE *fP,
					WlzDVertex2 org, WlzDVertex2 dst)
{
  const double	scaleX = 1.0,
  		scaleY = 1.0,
  		offsetX = 0.0,
		offsetY = 0.0;

  (void )fprintf(fP,
  		 "%f %f moveto "
  		 "%f %f lineto "
		 "stroke\n",
		 (org.vtX * scaleX) + offsetX,
		 (org.vtY * scaleY) + offsetY,
		 (dst.vtX * scaleX) + offsetX,
		 (dst.vtY * scaleY) + offsetY);
}

/************************************************************************
* Function:	WlzContourTestOutPS
* Returns:	WlzErrorNum		Woolz error code.
* Purpose:	Prints out a 2D contour as postscript for testing.
* Global refs:	-
* Parameters:	WlzContour *ctr:	Given contour to print out.
*		FILE *fP:		Output file.
*		int nCol:		Number of colours to cycle
*					through.
*		int headAndFoot:	Print a postscript header and
*					footer if non-zero.
************************************************************************/
WlzErrorNum	WlzContourTestOutPS(WlzContour *ctr, FILE *fP,
				    int nCol, int headAndFoot,
				    WlzDVertex2 offset, WlzDVertex2 scale)
{
  WlzGMModel	*model;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctr->type != WLZ_CONTOUR_TYPE_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(headAndFoot)
    {
      (void )fprintf(fP, "%%!\n"
		     "1 setlinecap\n"
		     "1 setlinejoin\n"
		     "0.01 setlinewidth\n");
    }
    if((model = ctr->model) != NULL)
    {
      errNum = WlzGMModelTestOutPS(model, fP, offset, scale, nCol);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(headAndFoot)
    {
      (void )fprintf(fP, "showpage\n");
    }
  }
  return(errNum);
}

/* #define TEST_WLZCONTOUR */
#ifdef TEST_WLZCONTOUR
/* Test main() for WlzContourObj(). */

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           option,
  		ok = 1,
		usage = 0,
  		colIdx = 0,
		nCol,
		maxCol = 7,
		minNod = 16,
		minEdg = 16;
  double	ctrVal = 100,
  		ctrWth = 1.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  WlzObject     *inObj = NULL;
  WlzContourMethod ctrMtd = WLZ_CONTOUR_MTD_ISO;
  WlzContour	*ctr;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzDVertex2	scale,
  		offset;
  static char	optList[] = "ghic:e:n:o:v:w:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  nCol = maxCol;
  outFileStr = (char *)outFileStrDef;
  inObjFileStr = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'g':
        ctrMtd = WLZ_CONTOUR_MTD_GRD;
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'i':
        ctrMtd = WLZ_CONTOUR_MTD_ISO;
	break;
      case 'c':
        if(sscanf(optarg, "%d", &nCol) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	else
	{
	  if(nCol < 1)
	  {
	    nCol = 1;
	  }
	  else if(nCol > maxCol)
	  {
	    nCol = maxCol;
	  }
	}
	break;
      case 'e':
        if(sscanf(optarg, "%d", &minEdg) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'n':
        if(sscanf(optarg, "%d", &minNod) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'v':
        if(sscanf(optarg, "%lg", &ctrVal) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'w':
        if(sscanf(optarg, "%lg", &ctrWth) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
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
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }

  }
  if(ok)
  {
  if((fP = (strcmp(outFileStr, "-")?
	   fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to open output file %s.\n",
      	      	     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    ctr = WlzContourObj(inObj, ctrMtd, ctrVal, ctrWth,
    			minNod, minEdg, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to compute contour.\n",
      	     	     argv[0]);
    }
  }
  if(ok)
  {
    scale.vtX = 0.5;
    scale.vtY = -0.5;
    offset.vtX = 100.0;
    offset.vtY = 700.0;
    errNum = WlzContourTestOutPS(ctr, fP, nCol, 1, offset, scale);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s Failed to output contour as postscript (%d).\n",
      		     argv[0], (int )errNum);
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
    static char	optList[] = "ghic:e:n:o:v:w:";
      (void )fprintf(stderr,
      "Usage: %s%sExample: %s%s",
      *argv,
      " [-o<output object>] [-h] [-o] [-g] [-i]\n"
      "        [-c#] [-e#] [-n#] [-o#] [-v#] [-w#]\n"
      "        [<input object>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -o  Output object file name.\n"
      "  -g  Compute maximal gradient contours.\n"
      "  -i  Compute iso-value contours.\n"
      "  -c  Cycle contour colours (useful for debuging).\n"
      "  -e  Minimum number of edges for each contour.\n"
      "  -n  Minimum number of nodes for each contour.\n"
      "  -v  Contour iso-value.\n"
      "  -w  Contour (Deriche) gradient operator width.\n"
      "Computes a contour list from the given input object.\n"
      "The input object is read from stdin and output data are written\n"
      "to stdout unless filenames are given.\n",
      *argv,
      " -i -v 0.0 -n 8 in.wlz\n"
      "The input Woolz object is read from in.wlz, and the iso-value\n"
      "(iso-value = 1.0) contour list is written to stdout.\n"
      "Contours with less than 8 nodes are ignored.\n");
  }
  return(!ok);
}

#endif /* TEST_WLZCONTOUR */
