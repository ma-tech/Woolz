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
* 15-12-00 bill Fix offset bug when getting 2D array data in
*		WlzContourIsoObj3D() and WlzContourGrdObj5D().
* 21-11-00 bill Fix bugs in 3D gradient contour generation and removed
*		some unused code.
* 25-08-00 bill	Fix more bugs causing holes in 3D iso-surface.
* 22-08-00 bill	Fix a bug causing holes in 3D iso-surface.
* 15-08-00 bill	Move WlzFreeContour to WlzFreeSpace.c and WlzMakeContour
*		to WlzMakeStructs.c.
* 10-08-00 bill Add WlzContourGrdObj3D() and modify WlzMakeContour().
* 05-06-00 bill Add 3D isosurface generation. Removed unused variables.
* 03-03-00 bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

#define WLZ_CTR_TOLERANCE	(1.0e-06)

/* #define WLZ_CONTOUR_DEBUG */
static void 	WlzContourTestOutPSLn2D(FILE *fP,
					WlzDVertex2 org, WlzDVertex2 dst);

/************************************************************************
* WlzContourTriIsn2D: Classification of intersection of line segment
* with a triangle.
************************************************************************/
typedef enum
{
  WLZ_CONTOUR_TIC2D_NONE,                             /*No intersection */
  WLZ_CONTOUR_TIC2D_V1V0,                             /* Vertex 1 - vertex 0 */
  WLZ_CONTOUR_TIC2D_V0V2,                             /* Vertex 0 - vertex 2 */
  WLZ_CONTOUR_TIC2D_V2V1,                             /* Vertex 2 - vertex 1 */
  WLZ_CONTOUR_TIC2D_V1S02,                            /* Vertex 1 - side 0-2 */
  WLZ_CONTOUR_TIC2D_V0S21,                            /* Vertex 0 - side 2-1 */
  WLZ_CONTOUR_TIC2D_V2S10,                            /* Vertex 2 - side 1-0 */
  WLZ_CONTOUR_TIC2D_S10S02,                           /* Side 1-0 - side 0-2 */
  WLZ_CONTOUR_TIC2D_S02S21,                           /* Side 0-2 - side 2-1 */
  WLZ_CONTOUR_TIC2D_S21S10                            /* Side 2-1 - side 1-0 */
} WlzContourTriIsn2D;

/************************************************************************
* WlzContourTetIsn3D: Classification of intersection of plane with a
* tetradedron.
************************************************************************/
typedef enum
{
  WLZ_CONTOUR_TIC3D_NONE,               /*No intersection */
  WLZ_CONTOUR_TIC3D_V0S12S13,		/* Vertex 0, side 1-2, side 1-3 */
  WLZ_CONTOUR_TIC3D_V0S12S23,		/* Vertex 0, side 1-2, side 2-3 */
  WLZ_CONTOUR_TIC3D_V0S13S23,		/* Vertex 0, side 1-3, side 2-3 */
  WLZ_CONTOUR_TIC3D_V0V1S23,		/* Vertex 0, vertex 1, side 2-3 */
  WLZ_CONTOUR_TIC3D_V0V1V2,		/* Vertex 0, vertex 1, vertex 2 */
  WLZ_CONTOUR_TIC3D_V0V1V3,		/* Vertex 0, vertex 2, side 1-3 */
  WLZ_CONTOUR_TIC3D_V0V2S13,		/* Vertex 0, vertex 2, vertex 3 */
  WLZ_CONTOUR_TIC3D_V0V2V3,		/* Vertex 0, vertex 2, vertex 3 */
  WLZ_CONTOUR_TIC3D_V0V3S12,		/* Vertex 0, side 0-2, side 0-3 */
  WLZ_CONTOUR_TIC3D_V1S02S03,		/* Vertex 1, side 0-2, side 0-3 */
  WLZ_CONTOUR_TIC3D_V1S02S23,		/* Vertex 1, side 0-2, side 2-3 */
  WLZ_CONTOUR_TIC3D_V1S03S23,		/* Vertex 1, side 0-3, side 2-3 */
  WLZ_CONTOUR_TIC3D_V1V2S03,		/* Vertex 1, vertex 2, side 0-3 */
  WLZ_CONTOUR_TIC3D_V1V2V3,		/* Vertex 1, vertex 2, vertex 3 */
  WLZ_CONTOUR_TIC3D_V1V3S02,		/* Vertex 1, vertex 3, side 0-2 */
  WLZ_CONTOUR_TIC3D_V2S01S02,		/* Vertex 2, side 0-1, side 0-2 */
  WLZ_CONTOUR_TIC3D_V2S01S03,		/* Vertex 2, side 0-1, side 0-3 */
  WLZ_CONTOUR_TIC3D_V2S01S13,		/* Vertex 2, side 0-1, side 1-3 */
  WLZ_CONTOUR_TIC3D_V2S03S13,		/* Vertex 2, side 0-3, side 1-3 */
  WLZ_CONTOUR_TIC3D_V2V3S01,		/* Vertex 2, vertex 3, side 0-1 */
  WLZ_CONTOUR_TIC3D_V3S01S02,		/* Vertex 3, side 0-1, side 0-2 */
  WLZ_CONTOUR_TIC3D_V3S01S12,		/* Vertex 3, side 0-1, side 1-2 */
  WLZ_CONTOUR_TIC3D_V3S02S12,		/* Vertex 3, side 0-2, side 1-2 */
  WLZ_CONTOUR_TIC3D_S01S02S03,		/* Side 0-1, side 0-2, side 0-3 */
  WLZ_CONTOUR_TIC3D_S01S02S13S23,  /* Side 0-1, side 0-2, side 1-3, side 2-3 */
  WLZ_CONTOUR_TIC3D_S01S03S12S23,  /* Side 0-1, side 0-3, side 1-2, side 2-3 */
  WLZ_CONTOUR_TIC3D_S01S12S13,		/* Side 0-1, side 1-2, side 1-3 */
  WLZ_CONTOUR_TIC3D_S02S03S12S13,  /* Side 0-2, side 0-3, side 1-2, side 1-3 */
  WLZ_CONTOUR_TIC3D_S02S12S23,		/* Side 0-2, side 1-2, side 2-3 */
  WLZ_CONTOUR_TIC3D_S03S13S23		/* Side 0-3, side 1-3, side 2-3 */
} WlzContourTetIsn3D;

static WlzContour	*WlzContourIsoObj2D(
			  WlzObject *srcObj,
			  double isoVal,
			  WlzErrorNum *dstErr);
static WlzContour	*WlzContourIsoObj3D(
			  WlzObject *srcObj,
			  double isoVal,
			  WlzErrorNum *dstErr);
static WlzContour 	*WlzContourGrdObj2D(
			  WlzObject *srcObj,
			  WlzDVertex2 **dstGrd,
			  double grdLo,
			  double grdHi,
			  double ftrPrm,
			  WlzErrorNum *dstErr);
static WlzContour 	*WlzContourGrdObj3D(
			  WlzObject *srcObj,
			  WlzDVertex3 **dstGrd,
			  double grdLo,
			  double grdHi,
			  double ftrPrm,
			  WlzErrorNum *dstErr);
static WlzDVertex2	WlzContourItpTriSide(
			  double valOrg,
			  double valDst,
			  WlzDVertex2 posOrg,
			  WlzDVertex2 posDst);
static WlzDVertex3 	WlzContourItpTetSide(
			  double valOrg,
			  double valDst,
			  WlzDVertex3 posOrg,
			  WlzDVertex3 posDst);
static WlzErrorNum	WlzContourIsoCube2D(
			  WlzContour *ctr,
			  double isoVal,
			  double *ln0,
			  double *ln1,
			  WlzDVertex2 sqOrg);
static WlzErrorNum 	WlzContourGrdCube3D(
			  WlzContour *ctr,
			  int *dstLnkFlg,
			  UBYTE ***mBuf,
			  double ***zBuf,
			  double ***yBuf,
			  double ***xBuf,
			  int *bufIdx,
			  WlzIVertex3 bufPos,
			  WlzIVertex3 cbOrg);
static WlzErrorNum	WlzContourIsoCube3D24T(WlzContour *ctr,
			  double isoVal,
			  double *pn0ln0,
			  double *pn0ln1,
			  double *pn1ln0,
			  double *pn1ln1,
			  WlzDVertex3 cbOrg);
static WlzErrorNum	WlzContourIsoCube3D6T(WlzContour *ctr,
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
			  int *dstLnkFlg,
			  UBYTE **grdDBuf,
			  WlzIVertex2 org,
			  int lnOff,
			  int lnIdx[],
			  int klP);
static WlzErrorNum	WlzContourGrdLink3D(
			  WlzContour *ctr,
			  int *dstLnkFlg,
			  UBYTE ***mBuf,
			  int *bufIdx,
			  WlzIVertex3 bufPos,
			  WlzIVertex3 cbOrg);

/************************************************************************
* Function:	WlzContourObjGrd
* Returns:	WlzContour:		Contour, or NULL on error.
* Purpose:	Creates a contour (list of connected edges or surface
*		patches) from a Woolz object with values using a
*		maximal gradient algorithm and retains the
*		gradient vectors.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		WlzVertexP *dstGrd:	Destination pointer for the
*					gradients, which are indexed
*					by the contour's GM vertex
*					indicies. The gradients are
*					always either 2D or 3D double
*					verticies. They are only valid
*					for valid vertex indicies.
*		double ctrLo:		Lower maximal gradient
*					threshold value.
*		double ctrHi:		Higher maximal gradient
*					threshold value.
*		double ctrWth:		Contour filter width.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzContour	*WlzContourObjGrd(WlzObject *srcObj,
				  WlzVertexP *dstGrd,
				  double ctrLo, double ctrHi, double ctrWth,
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
	ctr = WlzContourGrdObj2D(srcObj,
				(dstGrd)? &((*dstGrd).d2): NULL,
				ctrLo, ctrHi, ctrWth, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
	ctr = WlzContourGrdObj3D(srcObj,
				 (dstGrd)? &((*dstGrd).d3): NULL,
				 ctrLo, ctrHi, ctrWth, &errNum);
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
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzContour	*WlzContourObj(WlzObject *srcObj, WlzContourMethod ctrMtd,
			       double ctrVal, double ctrWth,
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
	    ctr = WlzContourIsoObj2D(srcObj, ctrVal, &errNum);
	    break;
	  case WLZ_CONTOUR_MTD_GRD:
	    ctr = WlzContourGrdObj2D(srcObj, NULL,
	    			     ctrVal, ctrVal, ctrWth, &errNum);
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
	    ctr = WlzContourIsoObj3D(srcObj, ctrVal, &errNum);
	    break;
	  case WLZ_CONTOUR_MTD_GRD:
	    ctr = WlzContourGrdObj3D(srcObj, NULL,
	    			     ctrVal, ctrVal, ctrWth, &errNum);
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
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourIsoObj2D(WlzObject *srcObj, double isoVal,
				      WlzErrorNum *dstErr)
{
  int		idX,
  		bufSz,
		bufLft,
		bufRgt,
		bufLnIdx,
  		itvLen,
		itvBufWidth;
  WlzDomain	srcDom;
  UBYTE		*itvBuf[2] = {NULL, NULL};
  double	*valBuf[2] = {NULL, NULL};
  WlzContour 	*ctr = NULL;
  WlzDVertex2	sqOrg;
  WlzIntervalWSpace srcIWSp;
  WlzGreyWSpace	srcGWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Create contour. */
  if((ctr = WlzMakeContour(&errNum)) != NULL)
  {
    ctr->model = WlzAssignGMModel(
    		 WlzGMModelNew(WLZ_GMMOD_2D, 0, 0, &errNum), NULL);
  }
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
				       sqOrg);
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
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourIsoObj3D(WlzObject *srcObj, double isoVal,
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
  WlzValues	dummyValues;
  WlzDomain	dummyDom,
  		srcDom;
  WlzIVertex2	bufSz,
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
  else if(srcObj->values.vox->values == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create contour. */
    if((ctr = WlzMakeContour(&errNum)) != NULL)
    {
      ctr->model = WlzAssignGMModel(
      		   WlzGMModelNew(WLZ_GMMOD_3D, 0, 0, &errNum), NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox3D = WlzBoundingBox3D(srcObj, &errNum);
  }
  /* Make buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz.vtX = bBox3D.xMax - bBox3D.xMin + 1;
    bufSz.vtY = bBox3D.yMax - bBox3D.yMin + 1;
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
        bufOff.vtX = bBox2D.xMin;
        bufOff.vtY = bBox2D.yMin;
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
	  lnIdx = bBox2D.yMin - bBox3D.yMin;
	  lnCnt = bBox2D.yMax - bBox2D.yMin; 			  /* NOT + 1 */
	  while(lnIdx < lnCnt)
	  {
            cbOrg.vtY = bBox3D.yMin + lnIdx;
	    tUP0 = *(itvBuf[bufIdx0] + lnIdx);
	    tUP1 = *(itvBuf[bufIdx0] + lnIdx + 1);
	    tUP2 = *(itvBuf[bufIdx1] + lnIdx);
	    tUP3 = *(itvBuf[bufIdx1] + lnIdx + 1);
	    klIdx = bBox2D.xMin - bBox3D.xMin;
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
		errNum = WlzContourIsoCube3D6T(ctr, isoVal,
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
    (void )WlzFreeObj(obj2D);
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
*		WlzDVertex2 **dstGrd:	Destination pointer for gradients,
*					may be NULL.
*		double grdLo:		Lower threshold for modulus of
*					gradient.
*		double grdHi:		Upper threshold for modulus of
*					gradient. TODO Use grdHi!
*		double ftrPrm:		Filter width parameter.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourGrdObj2D(WlzObject *srcObj, WlzDVertex2 **dstGrd,
				      double grdLo, double grdHi,
				      double ftrPrm, WlzErrorNum *dstErr)
{
  int		idX,
		dCode,
		iCnt,
		iLen,
		iBufSz,
		lnInc,
		lnkFlg,
		grdMaxIdx;
  UBYTE		doLn;
  double	tD0,
  		grdM0,
  		grdM1,
		grdM2,
		grdX0,
		grdY0,
		grdMLft,
		grdMRgt;
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
  WlzGMVertex	*mVtx;
  WlzContour 	*ctr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		lnIdx[3],
  		klIdx[3];
  WlzIVertex2	bufSz,
  		org,
  		posAbs,
		posRel;
   WlzDVertex2	pos;
  WlzDVertex2	*grd;
  WlzIntervalWSpace gXIWSp,
  		gYIWSp;
  WlzGreyWSpace	gXGWSp,
  		gYGWSp;
  AlcVector	*grdVec = NULL;
  const UBYTE   dTable[8] = {3, 2, 0, 1, 4, 5, 7, 6};
  const unsigned int grdBlkSz = 1000;       /* Any reasoonable size would do */

  /* Create a new contour. */
  if((ctr = WlzMakeContour(&errNum)) != NULL)
  {
    ctr->model = WlzAssignGMModel(
    		 WlzGMModelNew(WLZ_GMMOD_2D, 0, 0, &errNum), NULL);
  }
  if((errNum == WLZ_ERR_NONE) && dstGrd)
  {
    grdMaxIdx = -1;
    if((grdVec = AlcVectorNew(1, sizeof(WlzDVertex3),
				   grdBlkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute the partial derivatives. */
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
    /* Copy this interval of gradient data into buffers. */
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
	while((errNum == WLZ_ERR_NONE) && (klIdx[2] < bufSz.vtX))
	{
	  if(WLZ_BIT_GET(itvP[3], klIdx[0]) &&
	     WLZ_BIT_GET(itvP[3], klIdx[1]) &&
	     WLZ_BIT_GET(itvP[3], klIdx[2]))
	  {
	    /* Check if gradient magnitude at the centre of the neighbourhood
	     * is above the minimum gradient threshold. */
	    if((grdM0 = *(*(grdMBuf + lnIdx[1]) + klIdx[1])) >= grdLo)
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
	      errNum = WlzContourGrdLink2D(ctr, &lnkFlg,
	      				   grdDBuf, org, posRel.vtY - 2,
					   lnIdx, klIdx[0]);
	      /* If gradients are to be returned set the gradient for this
	       * pixel and if it's in the model, then  add the gradient to
	       * the gradient vector. */
	      if(grdVec && lnkFlg)
	      {
		pos.vtX = org.vtX + klIdx[0] + 1;
		pos.vtY = org.vtY + posRel.vtY - 1;
	        if((mVtx = WlzGMModelMatchVertexG2D(ctr->model, pos)) != NULL)
		{
		  if((grd = (WlzDVertex2 *)
		  	    AlcVectorExtendAndGet(grdVec, mVtx->idx)) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		  else
		  {
		    if(grdMaxIdx < mVtx->idx)
		    {
		      grdMaxIdx = mVtx->idx;
		    }
		    grd->vtX = *(*(grdXBuf + lnIdx[1]) + klIdx[1]);
		    grd->vtY = *(*(grdYBuf + lnIdx[1]) + klIdx[1]);
		  }
		}
	      }
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
    (void )WlzFreeObj(gXObj);
  }
  if(gYObj)
  {
    (void )WlzFreeObj(gYObj);
  }
  if(grdVec)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      *dstGrd = (WlzDVertex2 *)AlcVectorToArray1D(grdVec, 0, grdMaxIdx,
      						  NULL);
    }
    (void )AlcVectorFree(grdVec);
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
*		WlzDVertex3 **dstGrd:	Destination pointer for gradients,
*					may be NULL.
*		double grdLo:		Lower threshold for modulus of
*					gradient.
*		double grdHi:		Upper threshold for modulus of
*					gradient. TODO use grdHi!
*		double ftrPrm:		Filter width parameter.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourGrdObj3D(WlzObject *srcObj, WlzDVertex3 **dstGrd,
				      double grdLo, double grdHi,
				      double ftrPrm, WlzErrorNum *dstErr)
{
  int		cnt,
  		idX,
		idY,
		idZ,
		pnIdx,
  		pnCnt,
		cbInObj,
		grdMaxIdx,
		lnkFlg;
  double	gV,
  		sGV,
		grdLoSq;
  WlzDomain	srcDom;
  UBYTE		*iMP,
	  	*iMC,
	  	*iMN;
  WlzGMVertex	*mVtx;
  WlzRsvFilter	*ftr = NULL;
  WlzObject	*srcObj2D = NULL,
  		*xObj2D = NULL,
  		*yObj2D = NULL,
		*zObj2D = NULL,
  		*zObj = NULL;
  WlzContour 	*ctr = NULL;
  WlzDVertex3	pos;
  WlzDVertex3	*grd;
  WlzIVertex3	bufPos,
  		cbOrg;
  WlzIVertex2	bufSz,
  		bufOff;
  WlzIBox2	bBox2D;
  WlzIBox3	bBox3D;
  int		bufIdx[3],
  		iBufClr[3];
  UBYTE		**iBuf[3],
  		**mBuf[3];
  double	**xBuf[3],
  		**yBuf[3],
		**zBuf[3];
  AlcVector	*grdVec = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const unsigned int grdBlkSz = 10000;      /* Any reasoonable size would do */

  iBufClr[0] = iBufClr[1] = iBufClr[2] = 1;
  iBuf[0] = iBuf[1] = iBuf[2] = NULL;
  mBuf[0] = mBuf[1] = mBuf[2] = NULL;
  xBuf[0] = xBuf[1] = xBuf[2] = NULL;
  yBuf[0] = yBuf[1] = yBuf[2] = NULL;
  zBuf[0] = zBuf[1] = zBuf[2] = NULL;
  if((srcDom = srcObj->domain).core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(srcObj->values.vox->values == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(dstGrd)
  {
    grdMaxIdx = -1;
    if((grdVec = AlcVectorNew(1, sizeof(WlzDVertex3), grdBlkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create contour. */
    if((ctr = WlzMakeContour(&errNum)) != NULL)
    {
      ctr->model = WlzAssignGMModel(
      		   WlzGMModelNew(WLZ_GMMOD_3D, 0, 0, &errNum), NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox3D = WlzBoundingBox3D(srcObj, &errNum);
    bufSz.vtX = bBox3D.xMax - bBox3D.xMin + 1;
    bufSz.vtY = bBox3D.yMax - bBox3D.yMin + 1;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create a recursive filter and compute the partial derivatives
     * through planes. */
    if((ftr = WlzRsvFilterMakeFilter(WLZ_RSVFILTER_NAME_DERICHE_1,
				     ftrPrm, &errNum)) != NULL)
    {
      ftr->c *= 4.0;
      zObj = WlzAssignObject(
             WlzRsvFilterObj(srcObj, ftr, WLZ_RSVFILTER_ACTION_Z,
	     		     &errNum), NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create buffers. */
    if((AlcUnchar2Calloc(&(mBuf[0]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcUnchar2Calloc(&(mBuf[1]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcUnchar2Calloc(&(mBuf[2]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcBit2Calloc(&(iBuf[0]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcBit2Calloc(&(iBuf[1]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcBit2Calloc(&(iBuf[2]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(xBuf[0]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(xBuf[1]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(xBuf[2]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(yBuf[0]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(yBuf[1]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(yBuf[2]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(zBuf[0]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(zBuf[1]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&(zBuf[2]), bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Work down through the object. */
  if(errNum == WLZ_ERR_NONE)
  {
    /* Enforce a minimum gradient. */
    if((grdLoSq = grdLo * grdLo) < 1.0)
    {
      grdLoSq = 1.0;
    }
    pnIdx = 0;
    pnCnt = srcDom.p->lastpl - srcDom.p->plane1 + 1;
    while((errNum == WLZ_ERR_NONE) && (pnCnt-- > 0))
    { 
#ifdef WLZ_CONTOUR_DEBUG
    (void )fprintf(stderr, "WLZ_CONTOUR_DEBUG pnCnt = %d\n", pnCnt);
#endif /* WLZ_CONTOUR_DEBUG */
      bufIdx[0] = (pnIdx + 3 - 2) % 3;
      bufIdx[1] = (pnIdx + 3 - 1) % 3;
      bufIdx[2] = (pnIdx + 3 - 0) % 3;
      srcObj2D = WlzAssignObject(
      	         WlzMakeMain(WLZ_2D_DOMAINOBJ,
      			     *(srcObj->domain.p->domains + pnIdx),
			     *(srcObj->values.vox->values + pnIdx),
			     NULL, NULL, NULL), NULL);
      zObj2D = WlzAssignObject(
      	       WlzMakeMain(WLZ_2D_DOMAINOBJ,
			   *(zObj->domain.p->domains + pnIdx),
			   *(zObj->values.vox->values + pnIdx),
			   NULL, NULL, NULL), NULL);
      /* Clear maximal voxel buffer for this plane. */
      WlzValueSetUByte(*(mBuf[bufIdx[2]]), 0, bufSz.vtX * bufSz.vtY);
      if(srcObj2D->domain.core == NULL)
      {
	/* Clear next plane interval buffer if this is an empty plane. */
	if(iBufClr[bufIdx[2]] == 0)
	{
          WlzValueSetUByte(*(iBuf[bufIdx[2]]), 0, bufSz.vtX * bufSz.vtY);
	  iBufClr[bufIdx[2]] = 1; /* Set plane interval buffer cleared flag. */
	}
      }
      else
      {
	iBufClr[bufIdx[2]] = 0;
	bBox2D = WlzBoundingBox2D(srcObj2D, &errNum);
	/* Compute the partial derivatives through the columns and rows. */
	if(errNum == WLZ_ERR_NONE)
	{
	  xObj2D = WlzAssignObject(
	  	   WlzRsvFilterObj(srcObj2D, ftr, WLZ_RSVFILTER_ACTION_X,
	  			 &errNum), NULL);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  yObj2D = WlzAssignObject(
	  	   WlzRsvFilterObj(srcObj2D, ftr, WLZ_RSVFILTER_ACTION_Y,
	  			 &errNum), NULL);
	}
	/* Update buffers. */
	if(errNum == WLZ_ERR_NONE)
	{
	  bufOff.vtX = bBox2D.xMin;
	  bufOff.vtY = bBox2D.yMin;
	  if((bufOff.vtX < 0) || (bufOff.vtY < 0))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzToArray2D((void ***)&(iBuf[bufIdx[2]]), srcObj2D,
				bufSz, bufOff, 0, WLZ_GREY_BIT);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzToArray2D((void ***)&(xBuf[bufIdx[2]]), xObj2D,
				bufSz, bufOff, 0, WLZ_GREY_DOUBLE);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzToArray2D((void ***)&(yBuf[bufIdx[2]]), yObj2D,
				bufSz, bufOff, 0, WLZ_GREY_DOUBLE);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzToArray2D((void ***)&(zBuf[bufIdx[2]]), zObj2D,
				bufSz, bufOff, 0, WLZ_GREY_DOUBLE);
	}
	/* Check to see if there are 3 consecutive non-empty planes. */
	if((errNum == WLZ_ERR_NONE) &&
	   (iBufClr[bufIdx[0]] == 0) && (iBufClr[bufIdx[1]] == 0))
	{
	  /* Form mask in the previous plane interval mask for the previous,
	   * current and next planes interval masks. */
	  cnt = (bufSz.vtX / 8) * bufSz.vtY;
	  iMP = *(iBuf[bufIdx[0]]);
	  iMC = *(iBuf[bufIdx[1]]);
	  iMN = *(iBuf[bufIdx[2]]);
	  while(cnt-- > 0)
	  {
	    *iMP++ &= *iMC++ & *iMN++;
	  }
	  /* For each cube with all interval mask bits set find the x, y, z
	   * partial derivatives, cumpute the maximal gradient surface
	   * elements and add them to the model. */
	  bufPos.vtZ = bufIdx[0];
          cbOrg.vtZ = bBox3D.zMin + pnIdx - 2;
	  bufPos.vtY = bBox2D.yMin - bBox3D.yMin;
	  cbOrg.vtY = bBox2D.yMin;
	  while((errNum == WLZ_ERR_NONE) && (cbOrg.vtY < (bBox3D.yMax - 2)))
	  {
	    bufPos.vtX = bBox2D.xMin - bBox3D.xMin;
	    cbOrg.vtX = bBox2D.xMin;
	    while((errNum == WLZ_ERR_NONE) && (cbOrg.vtX < (bBox3D.xMax - 2)))
	    {
	      cbInObj = 1;
	      idY = 0;
	      while(cbInObj && (idY < 3))
	      {
	        idX = 0;
	  	iMP = *(iBuf[bufIdx[0]] + bufPos.vtY + idY);
		while(cbInObj && (idX < 3))
		{
		  cbInObj = WLZ_BIT_GET(iMP, bufPos.vtX + idX);
		  ++idX;
		}
		++idY;
	      }
	      if(cbInObj)
	      {
		/* Compute central voxel square of modulus of gradient. */
		gV = *(*(xBuf[bufIdx[1]] + bufPos.vtY + 1) + bufPos.vtX + 1);
		sGV = (gV * gV);
		gV = *(*(yBuf[bufIdx[1]] + bufPos.vtY + 1) + bufPos.vtX + 1);
		sGV += (gV * gV);
		gV = *(*(zBuf[bufIdx[1]] + bufPos.vtY + 1) + bufPos.vtX + 1);
		sGV += (gV * gV);
		/* If central voxel has a modulus of gradient greater than the
		 * threshold value then compute the maximal surface
		 * simplicies and add them to the model. */
		if(sGV >= grdLoSq)
		{
		  errNum = WlzContourGrdCube3D(ctr, &lnkFlg,
		  			       mBuf, zBuf, yBuf, xBuf,
		  			       bufIdx, bufPos, cbOrg);
		/* If gradients are to be returned set the gradient for this
		 * pixel and if it's in the model, then add the gradient to
		 * the gradient vector. */
		  if(grdVec && lnkFlg)
		  {
		    pos.vtX = cbOrg.vtX + 1;
		    pos.vtY = cbOrg.vtY + 1;
		    pos.vtZ = cbOrg.vtZ + 1; 
		    if((mVtx = WlzGMModelMatchVertexG3D(ctr->model,
		    					pos)) != NULL)
		    {
		      if((grd = (WlzDVertex3 *)
				AlcVectorExtendAndGet(grdVec,
						      mVtx->idx)) == NULL)
		      {
			errNum = WLZ_ERR_MEM_ALLOC;
		      }
		      else
		      {
			if(grdMaxIdx < mVtx->idx)
			{
			  grdMaxIdx = mVtx->idx;
			}
			grd->vtX = *(*(xBuf[bufIdx[1]] + bufPos.vtY + 1) +
				     bufPos.vtX + 1);
			grd->vtY = *(*(yBuf[bufIdx[1]] + bufPos.vtY + 1) +
				     bufPos.vtX + 1);
			grd->vtZ = *(*(zBuf[bufIdx[1]] + bufPos.vtY + 1) +
				     bufPos.vtX + 1);
		      }
		    }
		  }
		}
#ifdef WLZ_CONTOUR_DEBUG
		else
		{
		  (void )fprintf(stderr, "%d %d %d 0.0 0.0 0.0\n",
		  		 cbOrg.vtX + 1, cbOrg.vtY + 1, cbOrg.vtZ + 1);
		}
#endif /* WLZ_CONTOUR_DEBUG */
	      }
	      ++(bufPos.vtX);
	      ++(cbOrg.vtX);
	    }
	    ++(bufPos.vtY);
	    ++(cbOrg.vtY);
	  }
	}
	if(xObj2D)
	{
	  (void )WlzFreeObj(xObj2D);
	  xObj2D = NULL;
	}
	if(yObj2D)
	{
	  (void )WlzFreeObj(yObj2D);
	  yObj2D = NULL;
	}
      }
      if(srcObj2D)
      {
	WlzFreeObj(srcObj2D);
	srcObj2D = NULL;
      }
      if(zObj2D)
      {
	WlzFreeObj(zObj2D);
	zObj2D = NULL;
      }
      ++pnIdx;
    }
  }
  /* Tidy up on error. */
  if((errNum != WLZ_ERR_NONE) && (ctr != NULL))
  {
    (void )WlzFreeContour(ctr);
    ctr = NULL;
  }
  /* Free buffer and temporary objects. */
  if(xObj2D)
  {
    (void )WlzFreeObj(xObj2D);
  }
  if(yObj2D)
  {
    (void )WlzFreeObj(yObj2D);
  }
  if(zObj)
  {
    (void )WlzFreeObj(zObj);
  }
  for(idZ = 0; idZ < 3; ++idZ)
  {
    if(iBuf[idZ])
    {
      Alc2Free((void **)(iBuf[idZ]));
    }
    if(mBuf[idZ])
    {
      Alc2Free((void **)(mBuf[idZ]));
    }
    if(xBuf[idZ])
    {
      Alc2Free((void **)(xBuf[idZ]));
    }
    if(yBuf[idZ])
    {
      Alc2Free((void **)(yBuf[idZ]));
    }
    if(zBuf[idZ])
    {
      Alc2Free((void **)(zBuf[idZ]));
    }
  }
  if(grdVec)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      /* Allocate *dstGrd and copy the gradients. */
      *dstGrd = (WlzDVertex3 *)AlcVectorToArray1D(grdVec, 0, grdMaxIdx,
      						  NULL);
    }
    (void )AlcVectorFree(grdVec);
  }
  if(ftr)
  {
    WlzRsvFilterFreeFilter(ftr);
  }
  /* Set error code. */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctr);
}

/************************************************************************
* Function:	WlzContourGrdLink3D
* Returns:	WlzErrorNum		Woolz error code.
* Purpose:	Computes edge simplicies(s) linking the voxel at the
*		centre of a 2x3x3 neighbourhood to it's neighbours.
*
*                z=0 +---+---+---+
*                    | 9 | 5 |10 |
*                    +---+---+---+    z=1 +---+---+                  
*                    | 6 | 2 | 7 |        | 0 | C |                         
*                 y  +---+---+---+        +---+---+---+
*                 ^  |11 | 8 |12 |        | 3 | 1 | 4 |           
*                 |  +---+---+---+        +---+---+---+           
*                 --->x                                   
* Global refs:	-
* Parameters:	WlzContour *ctr: 	Contour being built.
*		int *dstLnkFlg:		Destination pointer for flag
*					set to non-zero value if pixel
*					linked. Must NOT be NULL.
*		UBYTE ***mBuf:		Buffers containing non-zero
*					values for maximal gradient
*					voxels.
*		int *bufIdx:		Z offsets into the buffers.
*		WlzIVertex3 bufPos:	Offset into the buffer for the
*					neighbourhoods origin.
*		WlzIVertex3 pos:	Absolute position of the 
*					neighbourhoods origin.
************************************************************************/
static WlzErrorNum	WlzContourGrdLink3D(WlzContour *ctr, int *dstLnkFlg,
					    UBYTE ***mBuf,
					    int *bufIdx, WlzIVertex3 bufPos,
					    WlzIVertex3 cbOrg)
{
  int		tI0,
  		tI1,
  		idN,
		mCnt,
		spxCnt = 0;
  unsigned int	tU0,
  		mMsk;
  WlzIVertex3	tIV0,
  		tIV1;
  WlzDVertex3	sIsn[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzIVertex3 nbrOffTb[13] =
  {
    /* Offsets into the maximal gradient voxel buffer for the 6, 18 and
     * 26 connected neighbours, used for setting the maximal gradient
     * voxels mask. Vertex members are in vtX, vtY, vtZ order.*/
    {0, 1, 1}, /* index  0,  6 connected neighbour */
    {1, 0, 1}, /* index  1,  6 connected neighbour */
    {1, 1, 0}, /* index  2,  6 connected neighbour */
    {0, 0, 1}, /* index  3, 18 connected neighbour */
    {2, 0, 1}, /* index  4, 18 connected neighbour */
    {1, 2, 0}, /* index  5, 18 connected neighbour */
    {0, 1, 0}, /* index  6, 18 connected neighbour */
    {2, 1, 0}, /* index  7, 18 connected neighbour */
    {1, 0, 0}, /* index  8, 18 connected neighbour */
    {0, 2, 0}, /* index  9, 26 connected neighbour */
    {2, 2, 0}, /* index 10, 26 connected neighbour */
    {0, 0, 0}, /* index 11, 26 connected neighbour */
    {2, 0, 0}  /* index 12, 26 connected neighbour */
  };
  const int sPr[15][2] =
  {
   { 0,  3},
   { 3,  1},
   { 1,  4},
   { 5,  2},
   { 9,  2},
   { 6,  2},
   {11,  2},
   { 8,  2},
   {12,  2},
   { 7,  2},
   {10,  2},
   { 6,  0},
   {11,  3},
   { 8,  1},
   {12,  4}
  };

#ifdef WLZ_CONTOUR_DEBUG
  (void )fprintf(stderr,
  		 "G %d %d %d\n",
  	  	 cbOrg.vtX + 1, cbOrg.vtY + 1, cbOrg.vtZ + 1);
#endif /* WLZ_CONTOUR_DEBUG */

  /* Create a maximal gradient neighbouring voxel mask. */
  idN = 0;
  mCnt = 0;
  mMsk = 0;
  while(idN < 13)
  {
    tIV0 = nbrOffTb[idN];
    tIV0.vtX += bufPos.vtX;
    tIV0.vtY += bufPos.vtY;
    tIV0.vtZ = *(bufIdx + tIV0.vtZ);
    if(*(*(*(mBuf + tIV0.vtZ) + tIV0.vtY) + tIV0.vtX) != 0)
    {
      mMsk |= 1 << idN;
      ++mCnt;
    }
    ++idN;
  }
  if(mCnt > 1)
  {
    sIsn[0].vtX = cbOrg.vtX + 1;
    sIsn[0].vtY = cbOrg.vtY + 1;
    sIsn[0].vtZ = cbOrg.vtZ + 1;
    /* Check for simplicies and add them to the model. */
    idN = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < 15))
    {
      tI0 = sPr[idN][0];
      tI1 = sPr[idN][1];
      tU0 = (1 << tI0) | (1 << tI1);
      if((mMsk & tU0) == tU0)
      {
        tIV0 = nbrOffTb[tI0];
        tIV1 = nbrOffTb[tI1];
        WLZ_VTX_3_ADD(sIsn[1], cbOrg, tIV0);
        WLZ_VTX_3_ADD(sIsn[2], cbOrg, tIV1);
	++spxCnt;;
	errNum = WlzGMModelConstructSimplex3D(ctr->model, sIsn);
      }
      ++idN;
    }
  }
  *dstLnkFlg = spxCnt;
  return(errNum);
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
*		int *dstLnkFlg:		Destination pointer for flag
*					set to non-zero value if pixel
*					linked. Must NOT be NULL.
*		UBYTE **grdDBuf:	Buffers containing gradient
*					direction codes for maximal
*					gradient pixels.
*		WlzIVertex2 org:	Origin of the object.
*		int lnOff:		Offset from origin to central
*					pixel line.
*		int lnIdx[]:		Line indicies.
*		int klOff:		Column offset of first column.
************************************************************************/
static WlzErrorNum WlzContourGrdLink2D(WlzContour *ctr, int *dstLnkFlg, 
				       UBYTE **grdDBuf, WlzIVertex2 org,
				       int lnOff, int lnIdx[], int klOff)
{
  int		idN,
		idM,
		conCnt;
  int		con[5],
     		nbr[4];
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
      errNum = WlzGMModelConstructSimplex2D(ctr->model, seg);
    }
  }
  *dstLnkFlg = conCnt;
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
				       WlzDVertex2 sqOrg)
{
  int		idT,
		tN1,
  		tN2,
		dupIsn = 0;
  WlzContourTriIsn2D iCode;
  WlzDVertex2	tVx0;
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
	WLZ_CONTOUR_TIC2D_V1V0,
	WLZ_CONTOUR_TIC2D_V1S02
      },
      {
        WLZ_CONTOUR_TIC2D_S21S10,
	WLZ_CONTOUR_TIC2D_V0S21,
	WLZ_CONTOUR_TIC2D_S02S21
      }
    },
    {
      {
        WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_V0V2,
	WLZ_CONTOUR_TIC2D_V2S10
      },
      {
        WLZ_CONTOUR_TIC2D_V2V1,
	WLZ_CONTOUR_TIC2D_NONE,
	WLZ_CONTOUR_TIC2D_V2V1
      },
      {
        WLZ_CONTOUR_TIC2D_V2S10,
	WLZ_CONTOUR_TIC2D_V0V2,
	WLZ_CONTOUR_TIC2D_NONE
	}
    },
    {
      {
        WLZ_CONTOUR_TIC2D_S02S21,
	WLZ_CONTOUR_TIC2D_V0S21,
	WLZ_CONTOUR_TIC2D_S21S10
      },
      {
        WLZ_CONTOUR_TIC2D_V1S02,
	WLZ_CONTOUR_TIC2D_V1V0,
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

  sV[1] = *vLn0 - isoVal; sV[2] = *(vLn0 + 1) - isoVal;
  sV[3] = *(vLn1 + 1) - isoVal; sV[4] = *vLn1 - isoVal;
#ifdef WLZ_CONTOUR_DEBUG
  (void )fprintf(stderr,
		 "%% {%f %f} %f %f %f %f\n",
		 sqOrg.vtX, sqOrg.vtY, sV[1], sV[2], sV[3], sV[4]);
#endif /* WLZ_CONTOUR_DEBUG */
  /* Reject all squares not intersected by the iso-value. */
  if(!(((sV[1] < -(WLZ_CTR_TOLERANCE)) && (sV[2] < -(WLZ_CTR_TOLERANCE)) &&
       (sV[3] < -(WLZ_CTR_TOLERANCE)) && (sV[4] < -(WLZ_CTR_TOLERANCE))) ||
      ((sV[1] > WLZ_CTR_TOLERANCE) && (sV[2] > WLZ_CTR_TOLERANCE) &&
       (sV[3] > WLZ_CTR_TOLERANCE) && (sV[4] > WLZ_CTR_TOLERANCE))))
  {
    /* Interpolate for centre of square */
    sV[0] = (sV[1] + sV[2] + sV[3] + sV[4]) / 4.0;
    tV[0] = sV[0];
    lev[0] = (tV[0] > WLZ_CTR_TOLERANCE) + (tV[0] > -(WLZ_CTR_TOLERANCE));
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
      lev[1] = (tV[1] > WLZ_CTR_TOLERANCE) + (tV[1] > -(WLZ_CTR_TOLERANCE));
      lev[2] = (tV[2] > WLZ_CTR_TOLERANCE) + (tV[2] > -(WLZ_CTR_TOLERANCE));
      iCode = iCodeTab[lev[2]][lev[1]][lev[0]];
      if(iCode)
      {
	switch(iCode)
	{
	  case WLZ_CONTOUR_TIC2D_V1V0:
	    tIsn[0] = sNOffTab[tN1];
	    tIsn[1] = sNOffTab[tN0];
	    break;
	  case WLZ_CONTOUR_TIC2D_V0V2:
	    tIsn[0] = sNOffTab[tN0];
	    tIsn[1] = sNOffTab[tN2];
	    break;
	  case WLZ_CONTOUR_TIC2D_V2V1:
	    tIsn[0] = sNOffTab[tN2];
	    tIsn[1] = sNOffTab[tN1];
	    break;
	  case WLZ_CONTOUR_TIC2D_V1S02:
	    tIsn[0] = sNOffTab[tN1];
	    tIsn[1] = WlzContourItpTriSide(tV[0], tV[2],
	    				   sNOffTab[tN0], sNOffTab[tN2]);
	    break;
	  case WLZ_CONTOUR_TIC2D_V0S21:
	    tIsn[0] = sNOffTab[tN0];
	    tIsn[1] = WlzContourItpTriSide(tV[2], tV[1],
	    				   sNOffTab[tN2], sNOffTab[tN1]);
	    break;
	  case WLZ_CONTOUR_TIC2D_V2S10:
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
	  errNum = WlzGMModelConstructSimplex2D(ctr->model, tIsn);
	}
      }
      ++idT;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIsoCube3D24
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
static WlzErrorNum WlzContourIsoCube3D24(WlzContour *ctr,
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

  /* Compute values relative to iso-surface. */
  cVal[ 1] = *(vPn0Ln0 + 0) - isoVal;
  cVal[ 2] = *(vPn0Ln0 + 1) - isoVal;
  cVal[ 3] = *(vPn0Ln1 + 1) - isoVal;
  cVal[ 4] = *(vPn0Ln1 + 0) - isoVal;
  cVal[ 5] = *(vPn1Ln0 + 0) - isoVal;
  cVal[ 6] = *(vPn1Ln0 + 1) - isoVal;
  cVal[ 7] = *(vPn1Ln1 + 1) - isoVal;
  cVal[ 8] = *(vPn1Ln1 + 0) - isoVal;
  /* Test to se if there is an intersection between this cube and the
   * iso-surface. */
  intersect = !(((cVal[ 1] < -(WLZ_CTR_TOLERANCE)) &&
                 (cVal[ 2] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[ 3] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[ 4] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[ 5] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[ 6] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[ 7] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[ 8] < -(WLZ_CTR_TOLERANCE))) || 
	        ((cVal[ 1] > WLZ_CTR_TOLERANCE) &&
		 (cVal[ 2] > WLZ_CTR_TOLERANCE) &&
		 (cVal[ 3] > WLZ_CTR_TOLERANCE) &&
		 (cVal[ 4] > WLZ_CTR_TOLERANCE) &&
		 (cVal[ 5] > WLZ_CTR_TOLERANCE) &&
		 (cVal[ 6] > WLZ_CTR_TOLERANCE) &&
		 (cVal[ 7] > WLZ_CTR_TOLERANCE) &&
		 (cVal[ 8] > WLZ_CTR_TOLERANCE)));
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
    tIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (tIdx < 24))
    {
      tI1 = tVxLUT[tIdx][1];
      tI2 = tVxLUT[tIdx][2];
      tI3 = tVxLUT[tIdx][3];
      tVal[1] = cVal[tI1]; tVal[2] = cVal[tI2]; tVal[3] = cVal[tI3];
      tPos[1] = cPos[tI1]; tPos[2] = cPos[tI2]; tPos[3] = cPos[tI3];
      errNum = WlzContourIsoTet3D(ctr, tVal, tPos, cbOrg);
      ++tIdx;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourGrdCube3D
* Returns:	WlzErrorNum		Woolz error code.
* Purpose:	Computes maximal gradient surface simplicies within a
*		3x3x3 cube.
*		  * Compute central gradient vector's direction (cGV).
*		  * Interpolate gradient at intersection of gradient
*		    vector with the cube's faces.
*		  * If modulus of gradient is greater than the modulus
*		    of the interpolated gradients at the cubes faces.
*		  *   Add position of central voxel to the maximal
*		      gradient voxels buffer.
*		  *   Link maximal gradient voxels forming new
*		      surface elements in the model.
*		The modulus of the central gradient is known to be
*		non-zero.
*		The buffer position wrt z is always modulo 3.
* Global refs:	-
* Parameters:	WlzContour *ctr: 	Contour being built.
*		int *dstLnkFlg:		Destination pointer for flag
*					set to non-zero value if pixel
*					linked. Must NOT be NULL.
*		UBYTE ***mBuf:		Buffer with maximal voxels
*					marked non-zero.
*		double ***zBuf:		Z gradients.
*		double ***yBuf:		Y gradients.
*		double ***xBuf:		X gradients.
*		int *bufIdx:		Z indicies for the buffers.
*		WlzIVertex3 bufPos:	Index into buffers.
*		WlzIVertex3 cbOrg:	Origin of the cube.
************************************************************************/
static WlzErrorNum WlzContourGrdCube3D(WlzContour *ctr, int *dstLnkFlg,
				       UBYTE ***mBuf,
				       double ***zBuf, double ***yBuf,
				       double ***xBuf,
				       int *bufIdx, WlzIVertex3 bufPos,
				       WlzIVertex3 cbOrg)
{
  int		idF,
  		clsIdx;
  double	tD0,
		modCGV,
		tCG;
  WlzIVertex3	aFQ,
		aFQ1,
  		oFQ,
		aCC;
  WlzDVertex3	aCGV,
  		cGV,
		lI,
		pFV,
		tV,
		fS;
  double	gF[2],
  		lIP[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  /* Compute gradient vector at centre of cube as a unit vector and it's
   * modulus. */
  aCC.vtX = bufPos.vtX + 1;
  aCC.vtY = bufPos.vtY + 1;
  aCC.vtZ = *(bufIdx + 1);
  cGV.vtX = *(*(*(xBuf + aCC.vtZ) + aCC.vtY) + aCC.vtX);
  cGV.vtY = *(*(*(yBuf + aCC.vtZ) + aCC.vtY) + aCC.vtX);
  cGV.vtZ = *(*(*(zBuf + aCC.vtZ) + aCC.vtY) + aCC.vtX);
  modCGV = WLZ_VTX_3_LENGTH(cGV);
  WLZ_VTX_3_SCALE(cGV, cGV, 1.0 / modCGV);
  /* Clasify intersection using greatest absolute gradient vector component. */
  aCGV.vtX = fabs(cGV.vtX);
  aCGV.vtY = fabs(cGV.vtY);
  aCGV.vtZ = fabs(cGV.vtZ);
  clsIdx = (aCGV.vtX >= aCGV.vtY) | ((aCGV.vtY >= aCGV.vtZ) << 1) |
	   ((aCGV.vtZ >= aCGV.vtX) << 2);
  /* Switch on the clasification index using symetry to compute the
   * intersection of the central gradient vector with the cubes
   * faces. */
  switch(clsIdx)
  {
    case 1: /* FALLTHROUGH */
    case 3:
    case 7:
      /* Gradient vector passes through +/- x faces. */
      for(idF = 0; idF < 2; ++idF)
      {
	fS.vtX = (idF * 2.0) - 1.0; 	    /* Face: left -1.0 or right +1.0 */
	/* Intersection y, z coordinates. */
	tCG = cGV.vtX * fS.vtX;
	pFV.vtY = cGV.vtY / tCG;
	pFV.vtZ = cGV.vtZ / tCG;
	/* Indicies into the cube array (origin of face quadrant). */
	oFQ.vtX = (tCG >= 0.0) * 2;
	oFQ.vtY = pFV.vtY >= 0.0;
	oFQ.vtZ = pFV.vtZ >= 0.0;
	aFQ.vtX = bufPos.vtX + oFQ.vtX;
	aFQ.vtY = bufPos.vtY + oFQ.vtY;
	aFQ.vtZ = *(bufIdx + oFQ.vtZ);
	aFQ1.vtX = aFQ.vtX;
	aFQ1.vtY = aFQ.vtY + 1;
	aFQ1.vtZ = *(bufIdx + ((oFQ.vtZ + 1) % 3));
	/* Interpolation coefficients lI.vtY and lI.vtZ
	 *
	 *       +ve X face
	 *     O---------------O -
	 *     |               | ^
	 *     |        pFV    | | (1 - lI.vtY)
	 *     |         \     | v
	 *     |          X    | -
	 *     |               | ^
	 *     |               | | lI.vtY
	 *     |           oFQ | |
	 *     |              \| v                    
	 *     O---------------O -                    
	 *     |<-------->|<-->|                    
	 *        \          \             Y      
	 *      (1 - lI.vtZ)  \            ^       
	 *                    lI.vtZ       |             
	 *                                 |                          
	 *                          Z<-----O                          
	 *                                  \     
	 *                                   \                        
	 *                                    X                       
	 *                                                            
	 * Bevasuse pFV.vtZ is in the range [-1.0, +1.0]  and oFQ.vtZ is in
	 * the range [0, 1] then compute the interpolation coeficient by:
	 *   lI.vtZ = (pFV.vtZ + 1.0) - oFQ.vtZ
	 */
	lI.vtY = 1.0 + pFV.vtY - oFQ.vtY;
	lI.vtZ = 1.0 + pFV.vtZ - oFQ.vtZ;
	lIP[0] = (1.0 - lI.vtY) * (1.0 - lI.vtZ);
	lIP[1] = lI.vtY * (1.0 - lI.vtZ);
	lIP[2] = (1.0 - lI.vtY) * lI.vtZ;
	lIP[3] = lI.vtY * lI.vtZ;
	/* Interpolate the modulus of the gradient on the face. */
	tV.vtX = (*(*(*(xBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(xBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[1]) +
		 (*(*(*(xBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[2]) +
		 (*(*(*(xBuf + aFQ1.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[3]);
	tV.vtY = (*(*(*(yBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(yBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[1]) +
		 (*(*(*(yBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[2]) +
		 (*(*(*(yBuf + aFQ1.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[3]);
	tV.vtZ = (*(*(*(zBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(zBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[1]) +
		 (*(*(*(zBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[2]) +
		 (*(*(*(zBuf + aFQ1.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[3]);
	gF[idF] = WLZ_VTX_3_LENGTH(tV);
      }
      break;
    case 2: /* FALLTHROUGH */
    case 6:
      /* Gradient vector passes through +/- y faces. */
      for(idF = 0; idF < 2; ++idF)
      {
	fS.vtY = (idF * 2.0) - 1.0; 	    /* Face: front -1.0 or back +1.0 */
	/* Intersection z, x coordinates. */
	tCG = cGV.vtY * fS.vtY;
	pFV.vtX = cGV.vtX / tCG;
	pFV.vtZ = cGV.vtZ / tCG;
	/* Indicies into the cube array (origin of face quadrant). */
	oFQ.vtX = pFV.vtX >= 0.0;
	oFQ.vtY = (tCG >= 0.0) << 1;
	oFQ.vtZ = pFV.vtZ >= 0.0;
	aFQ.vtX = bufPos.vtX + oFQ.vtX;
	aFQ.vtY = bufPos.vtY + oFQ.vtY;
	aFQ.vtZ = *(bufIdx + oFQ.vtZ);
	aFQ1.vtX = aFQ.vtX + 1;
	aFQ1.vtY = aFQ.vtY;
	aFQ1.vtZ = *(bufIdx + ((oFQ.vtZ + 1) % 3));
	/* Interpolation coefficients lI.vtX and lI.vtY. */
	lI.vtZ = 1.0 + pFV.vtZ - oFQ.vtZ;
	lI.vtX = 1.0 + pFV.vtX - oFQ.vtX;
	lIP[0] = (1.0 - lI.vtZ) * (1.0 - lI.vtX);
	lIP[1] = lI.vtZ * (1.0 - lI.vtX);
	lIP[2] = (1.0 - lI.vtZ) * lI.vtX;
	lIP[3] = lI.vtZ * lI.vtX;
	/* Interpolate the modulus of the gradient on the face. */
	tV.vtX = (*(*(*(xBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(xBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[1]) +
		 (*(*(*(xBuf + aFQ.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[2]) +
		 (*(*(*(xBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[3]);
	tV.vtY = (*(*(*(yBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(yBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[1]) +
		 (*(*(*(yBuf + aFQ.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[2]) +
		 (*(*(*(yBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[3]);
	tV.vtZ = (*(*(*(zBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(zBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[1]) +
		 (*(*(*(zBuf + aFQ.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[2]) +
		 (*(*(*(zBuf + aFQ1.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[3]);
	gF[idF] = WLZ_VTX_3_LENGTH(tV);
      }
      break;
    case 4: /* FALLTHROUGH */
    case 5:
      /* Gradient vector passes through +/- z faces. */
      for(idF = 0; idF < 2; ++idF)
      {
	fS.vtZ = (idF * 2.0) - 1.0; 	    /* Face: bottom -1.0 or top +1.0 */
	/* Intersection x, y coordinates. */
	tCG = cGV.vtZ * fS.vtZ;
	pFV.vtX = cGV.vtX / tCG;
	pFV.vtY = cGV.vtY / tCG;
	/* Indicies into the cube array (origin of face quadrant). */
	oFQ.vtX = pFV.vtX >= 0.0;
	oFQ.vtY = pFV.vtY >= 0.0;
	oFQ.vtZ = (tCG >= 0.0) << 1;
	aFQ.vtX = bufPos.vtX + oFQ.vtX;
	aFQ.vtY = bufPos.vtY + oFQ.vtY;
	aFQ.vtZ = *(bufIdx + oFQ.vtZ);
	aFQ1.vtX = aFQ.vtX + 1;
	aFQ1.vtY = aFQ.vtY + 1;
	aFQ1.vtZ = *(bufIdx + ((oFQ.vtZ + 1) % 3));
	/* Interpolation coefficients lI.vtX and lI.vtY. */
	lI.vtX = 1.0 + pFV.vtX - oFQ.vtX;
	lI.vtY = 1.0 + pFV.vtY - oFQ.vtY;
	lIP[0] = (1.0 - lI.vtX) * (1.0 - lI.vtY);
	lIP[1] = lI.vtX * (1.0 - lI.vtY);
	lIP[2] = (1.0 - lI.vtX) * lI.vtY;
	lIP[3] = (1.0 - lI.vtX) * (1.0 - lI.vtY);
	/* Interpolate the modulus of the gradient on the face. */
	tV.vtX = (*(*(*(xBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(xBuf + aFQ.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[1]) +
		 (*(*(*(xBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[2]) +
		 (*(*(*(xBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ1.vtX) * lIP[3]);
	tV.vtY = (*(*(*(yBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(yBuf + aFQ.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[1]) +
		 (*(*(*(yBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[2]) +
		 (*(*(*(yBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ1.vtX) * lIP[3]);
	tV.vtZ = (*(*(*(zBuf + aFQ.vtZ) + aFQ.vtY) + aFQ.vtX) * lIP[0]) +
		 (*(*(*(zBuf + aFQ.vtZ) + aFQ.vtY) + aFQ1.vtX) * lIP[1]) +
		 (*(*(*(zBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ.vtX) * lIP[2]) +
		 (*(*(*(zBuf + aFQ.vtZ) + aFQ1.vtY) + aFQ1.vtX) * lIP[3]);
	gF[idF] = WLZ_VTX_3_LENGTH(tV);
      }
      break;
  }
  /* Test for maximal gradient along the gradient vector. Is the gradient at
   * the centre of the cube greater than the interpolated gradient at it's
   * faces? */
  if((modCGV > gF[0]) && (modCGV > gF[1]))
  {
    /* Mark this voxel as having a maximal gradient and then link with
     * neighbours to form the model. */
#ifdef WLZ_CONTOUR_DEBUG
    (void )fprintf(stderr, "%d %d %d %g %g %g\n",
    		   cbOrg.vtX + 1, cbOrg.vtY + 1, cbOrg.vtZ + 1,
                   modCGV * cGV.vtX, modCGV * cGV.vtY, modCGV * cGV.vtZ);
#endif /* WLZ_CONTOUR_DEBUG */
    *(*(*(mBuf + aCC.vtZ) + aCC.vtY) + aCC.vtX) = 1;
    errNum = WlzContourGrdLink3D(ctr, dstLnkFlg, mBuf, bufIdx, bufPos, cbOrg);
  }
#ifdef WLZ_CONTOUR_DEBUG
  else
  {
    (void )fprintf(stderr, "%d %d %d 0.0 0.0 0.0\n",
    cbOrg.vtX + 1, cbOrg.vtY + 1, cbOrg.vtZ + 1);
  }
#endif /* WLZ_CONTOUR_DEBUG */
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIsoCube3D6T
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose: 	Checks to see if all above the cube's vertex values are
*		either above or below the iso-value. If they are then
*		there's no intersection between the iso-surface and the
*		cube. If there is a possible intersection then the cube
*		is split into 6 tetrahedra and each tetrahedra is
*		passed to WlzContourIsoTet3D().
*		With the 8 data values at the verticies @[0-7] the
*		cube is split into 6 tetrahedra.
*                                                                      
*                          @7----------------@6                        
*                         /|                /|                        
*                        / |               / |                         
*                       /  |              /  |                         
*                      /   |             /   |                         
*                     /    |            /    |                         
*                    /     |           /     |                         
*                   /      |          /      |                         
*                  /       |         /       |                         
*                 @4----------------@5       |                         
*                 |        |        |        |                         
*                 |        @3-------|--------@2                                 
*                 |       /         |       /                                 
*                 |      /          |      /                                  
*                 |     /           |     /                                  
*                 |    /            |    /                                  
*                 |   /             |   /                                  
*                 |  /              |  /                                  
*                 | /               | /                                  
*                 |/                |/                                  
*                 @0----------------@1                                 
*
*		The tetrahedra are are assigned the following indicies:
*
*		  tetradedron index  	cube vertex indicies
*		   0       		0, 1, 3, 5
*		   1       		1, 2, 3, 5
*		   2       		2, 3, 5, 6
*		   3       		3, 5, 6, 7
*		   4        		3, 4, 5, 7
*		   5        		0, 3, 4, 5
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
static WlzErrorNum WlzContourIsoCube3D6T(WlzContour *ctr,
				double isoVal,
				double *vPn0Ln0, double *vPn0Ln1,
				double *vPn1Ln0, double *vPn1Ln1,
				WlzDVertex3 cbOrg)
{
  int		tI0,
  		tIdx,
		vIdx,
  		intersect;
  double 	cVal[8],	  /* Cube's values relative to the iso-value */
  		tVal[4];			       /* Tetrahedron values */
  WlzDVertex3	tPos[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	tVxLUT[6][4] =  /* Tetrahedron to cube vertex look up table */
  {
    {0, 1, 3, 5}, {1, 2, 3, 5}, {2, 3, 5, 6},
    {3, 5, 6, 7}, {3, 4, 5, 7}, {0, 3, 4, 5}
  };
  const WlzDVertex3 cPos[8] =	 /* Cube positions, order is {vtX, vtY, vtZ} */
  {
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}
  };

  /* Compute values relative to iso-surface. */
  cVal[0] = *(vPn0Ln0 + 0) - isoVal;
  cVal[1] = *(vPn0Ln0 + 1) - isoVal;
  cVal[2] = *(vPn0Ln1 + 1) - isoVal;
  cVal[3] = *(vPn0Ln1 + 0) - isoVal;
  cVal[4] = *(vPn1Ln0 + 0) - isoVal;
  cVal[5] = *(vPn1Ln0 + 1) - isoVal;
  cVal[6] = *(vPn1Ln1 + 1) - isoVal;
  cVal[7] = *(vPn1Ln1 + 0) - isoVal;
  /* Test to se if there is an intersection between this cube and the
   * iso-surface. */
  intersect = !(((cVal[1] < -(WLZ_CTR_TOLERANCE)) &&
                 (cVal[2] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[3] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[4] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[5] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[6] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[7] < -(WLZ_CTR_TOLERANCE)) &&
		 (cVal[8] < -(WLZ_CTR_TOLERANCE))) || 
	        ((cVal[1] > WLZ_CTR_TOLERANCE) &&
		 (cVal[2] > WLZ_CTR_TOLERANCE) &&
	 	 (cVal[3] > WLZ_CTR_TOLERANCE) &&
		 (cVal[4] > WLZ_CTR_TOLERANCE) &&
		 (cVal[5] > WLZ_CTR_TOLERANCE) &&
		 (cVal[6] > WLZ_CTR_TOLERANCE) &&
		 (cVal[7] > WLZ_CTR_TOLERANCE) &&
		 (cVal[8] > WLZ_CTR_TOLERANCE)));
  if(intersect)
  {
    tIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (tIdx < 6))
    {
      for(vIdx = 0; vIdx < 4; ++vIdx)
      {
        tI0 = tVxLUT[tIdx][vIdx];
	tVal[vIdx] = cVal[tI0];
	tPos[vIdx] = cPos[tI0];
      }
      errNum = WlzContourIsoTet3D(ctr, tVal, tPos, cbOrg);
      ++tIdx;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIsoTet3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose: 	Computes the intersection of the given tetrahedron
*		with the isovalue surface. Creates new 3D simplicies
*		for the contour.
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
  int		idx,
  		isnCnt;
  double	tD0,
  		tD1;
  WlzDVertex3	*tVP0;
  WlzDVertex3	tV0;
  WlzContourTetIsn3D iCode;
  int		lev[4];     /* Rel. tetra node level: above 2, on 1, below 0 */
  WlzDVertex3	sIsn[3],
  		tIsn[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzContourTetIsn3D iCodeTab[3][3][3][3] = /* Intersection code table */
  {
    {
      {
	{
	  /* 0000 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 0001 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 0002 */ WLZ_CONTOUR_TIC3D_S01S02S03
	},
	{
	  /* 0010 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 0011 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 0012 */ WLZ_CONTOUR_TIC3D_V1S02S03
	},
	{
	  /* 0020 */ WLZ_CONTOUR_TIC3D_S01S12S13,
	  /* 0021 */ WLZ_CONTOUR_TIC3D_V0S12S13,
	  /* 0022 */ WLZ_CONTOUR_TIC3D_S02S03S12S13
	}
      },
      {
	{
	  /* 0100 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 0101 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 0102 */ WLZ_CONTOUR_TIC3D_V2S01S03
	},
	{
	  /* 0110 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 0111 */ WLZ_CONTOUR_TIC3D_V0V1V2,
	  /* 0112 */ WLZ_CONTOUR_TIC3D_V1V2S03
	},
	{
	  /* 0120 */ WLZ_CONTOUR_TIC3D_V2S01S13,
	  /* 0121 */ WLZ_CONTOUR_TIC3D_V0V2S13,
	  /* 0122 */ WLZ_CONTOUR_TIC3D_V2S03S13
	}
      },
      {
	{
	  /* 0200 */ WLZ_CONTOUR_TIC3D_S02S12S23,
	  /* 0201 */ WLZ_CONTOUR_TIC3D_V0S12S23,
	  /* 0202 */ WLZ_CONTOUR_TIC3D_S01S03S12S23
	},
	{
	  /* 0210 */ WLZ_CONTOUR_TIC3D_V1S02S23,
	  /* 0211 */ WLZ_CONTOUR_TIC3D_V0V1S23,
	  /* 0212 */ WLZ_CONTOUR_TIC3D_V1S03S23
	},
	{
	  /* 0220 */ WLZ_CONTOUR_TIC3D_S01S02S13S23,
	  /* 0221 */ WLZ_CONTOUR_TIC3D_V0S13S23,
	  /* 0222 */ WLZ_CONTOUR_TIC3D_S03S13S23
	}
      }
    },
    {
      {
	{
	  /* 1000 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 1001 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 1002 */ WLZ_CONTOUR_TIC3D_V3S01S02
	},
	{
	  /* 1010 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 1011 */ WLZ_CONTOUR_TIC3D_V0V1V3,
	  /* 1012 */ WLZ_CONTOUR_TIC3D_V1V3S02
	},
	{
	  /* 1020 */ WLZ_CONTOUR_TIC3D_V3S01S12,
	  /* 1021 */ WLZ_CONTOUR_TIC3D_V0V3S12,
	  /* 1022 */ WLZ_CONTOUR_TIC3D_V3S02S12
	}
      },
      {
	{
	  /* 1100 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 1101 */ WLZ_CONTOUR_TIC3D_V0V2V3,
	  /* 1102 */ WLZ_CONTOUR_TIC3D_V2V3S01
	},
	{
	  /* 1110 */ WLZ_CONTOUR_TIC3D_V1V2V3,
	  /* 1111 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 1112 */ WLZ_CONTOUR_TIC3D_V1V2V3
	},
	{
	  /* 1120 */ WLZ_CONTOUR_TIC3D_V2V3S01,
	  /* 1121 */ WLZ_CONTOUR_TIC3D_V0V2V3,
	  /* 1122 */ WLZ_CONTOUR_TIC3D_NONE
	}
      },
      {
	{
	  /* 1200 */ WLZ_CONTOUR_TIC3D_V3S02S12,
	  /* 1201 */ WLZ_CONTOUR_TIC3D_V0V3S12,
	  /* 1202 */ WLZ_CONTOUR_TIC3D_V3S01S12
	},
	{
	  /* 1210 */ WLZ_CONTOUR_TIC3D_V1V3S02,
	  /* 1211 */ WLZ_CONTOUR_TIC3D_V0V1V3,
	  /* 1212 */ WLZ_CONTOUR_TIC3D_NONE
	},
	{
	  /* 1220 */ WLZ_CONTOUR_TIC3D_V3S01S02,
	  /* 1221 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 1222 */ WLZ_CONTOUR_TIC3D_NONE
	}
      }
    },
    {
      {
	{
	  /* 2000 */ WLZ_CONTOUR_TIC3D_S03S13S23,
	  /* 2001 */ WLZ_CONTOUR_TIC3D_V0S13S23,
	  /* 2002 */ WLZ_CONTOUR_TIC3D_S01S02S13S23
	},
	{
	  /* 2010 */ WLZ_CONTOUR_TIC3D_V1S03S23,
	  /* 2011 */ WLZ_CONTOUR_TIC3D_V0V1S23,
	  /* 2012 */ WLZ_CONTOUR_TIC3D_V1S02S23
	},
	{
	  /* 2020 */ WLZ_CONTOUR_TIC3D_S01S03S12S23,
	  /* 2021 */ WLZ_CONTOUR_TIC3D_V0S12S23,
	  /* 2022 */ WLZ_CONTOUR_TIC3D_S02S12S23
	}
      },
      {
	{
	  /* 2100 */ WLZ_CONTOUR_TIC3D_V2S03S13,
	  /* 2101 */ WLZ_CONTOUR_TIC3D_V0V2S13,
	  /* 2102 */ WLZ_CONTOUR_TIC3D_V2S01S13
	},
	{
	  /* 2110 */ WLZ_CONTOUR_TIC3D_V1V2S03,
	  /* 2111 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 2112 */ WLZ_CONTOUR_TIC3D_NONE
	},
	{
	  /* 2120 */ WLZ_CONTOUR_TIC3D_V2S01S03,
	  /* 2121 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 2122 */ WLZ_CONTOUR_TIC3D_NONE
	}
      },
      {
	{
	  /* 2200 */ WLZ_CONTOUR_TIC3D_S02S03S12S13,
	  /* 2201 */ WLZ_CONTOUR_TIC3D_V0S12S13,
	  /* 2202 */ WLZ_CONTOUR_TIC3D_S01S12S13
	},
	{
	  /* 2210 */ WLZ_CONTOUR_TIC3D_V1S02S03,
	  /* 2211 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 2212 */ WLZ_CONTOUR_TIC3D_NONE
	},
	{
	  /* 2220 */ WLZ_CONTOUR_TIC3D_S01S02S03,
	  /* 2221 */ WLZ_CONTOUR_TIC3D_NONE,
	  /* 2222 */ WLZ_CONTOUR_TIC3D_NONE
	}
      }
    }
  };

  for(idx = 0; idx < 4; ++idx)
  {
    lev[idx] = (tVal[idx] >= DBL_EPSILON) + (tVal[idx] > -(DBL_EPSILON));
  }
  iCode = iCodeTab[lev[3]][lev[2]][lev[1]][lev[0]];
  switch(iCode)
  {
    case WLZ_CONTOUR_TIC3D_V0S12S13:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V0S12S23:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V0S13S23:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V0V1S23:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = *(tPos + 1);
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V0V1V2:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = *(tPos + 1);
      *(tIsn + 2) = *(tPos + 2);
      break;
    case WLZ_CONTOUR_TIC3D_V0V1V3:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = *(tPos + 1);
      *(tIsn + 2) = *(tPos + 3);
      break;
    case WLZ_CONTOUR_TIC3D_V0V2S13:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = *(tPos + 2);
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V0V2V3:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = *(tPos + 2);
      *(tIsn + 2) = *(tPos + 3);
      break;
    case WLZ_CONTOUR_TIC3D_V0V3S12:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 0);
      *(tIsn + 1) = *(tPos + 3);
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      break;
    case WLZ_CONTOUR_TIC3D_V1S02S03:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 1);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V1S02S23:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 1);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V1S03S23:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 1);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V1V2S03:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 1);
      *(tIsn + 1) = *(tPos + 2);
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V1V2V3:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 1);
      *(tIsn + 1) = *(tPos + 2);
      *(tIsn + 2) = *(tPos + 3);
      break;
    case WLZ_CONTOUR_TIC3D_V1V3S02:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 1);
      *(tIsn + 1) = *(tPos + 3);
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      break;
    case WLZ_CONTOUR_TIC3D_V2S01S02:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 2);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      break;
    case WLZ_CONTOUR_TIC3D_V2S01S03:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 2);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V2S01S13:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 2);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V2S03S13:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 2);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_V2V3S01:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 2);
      *(tIsn + 1) = *(tPos + 3);
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      break;
    case WLZ_CONTOUR_TIC3D_V3S01S02:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 3);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      break;
    case WLZ_CONTOUR_TIC3D_V3S01S12:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 3);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      break;
    case WLZ_CONTOUR_TIC3D_V3S02S12:
      isnCnt = 3;
      *(tIsn + 0) = *(tPos + 3);
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      break;
    case WLZ_CONTOUR_TIC3D_S01S02S03:
      isnCnt = 3;
      *(tIsn + 0) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_S01S02S13S23:
      isnCnt = 4;
      /* Take care with ordering as later assume ordered verticies when
       * spliting quadrilateral: S01,S02,S23,S13. */
      *(tIsn + 0) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      *(tIsn + 3) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_S01S03S12S23:
      isnCnt = 4;
      /* Take care with ordering as later assume ordered verticies when
       * spliting quadrilateral: S01,S12,S23,S03. */
      *(tIsn + 0) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      *(tIsn + 3) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_S01S12S13:
      isnCnt = 3;
      *(tIsn + 0) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 1),
					 *(tPos + 0), *(tPos + 1));
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_S02S03S12S13:
      isnCnt = 4;
      /* Take care with ordering as later assume ordered verticies when
       * spliting quadrilateral: S03,S13,S12,S02. */
      *(tIsn + 0) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      *(tIsn + 3) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      break;
    case WLZ_CONTOUR_TIC3D_S02S12S23:
      isnCnt = 3;
      *(tIsn + 0) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 2),
					 *(tPos + 0), *(tPos + 2));
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 2),
					 *(tPos + 1), *(tPos + 2));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      break;
    case WLZ_CONTOUR_TIC3D_S03S13S23:
      isnCnt = 3;
      *(tIsn + 0) = WlzContourItpTetSide(*(tVal + 0), *(tVal + 3),
					 *(tPos + 0), *(tPos + 3));
      *(tIsn + 1) = WlzContourItpTetSide(*(tVal + 1), *(tVal + 3),
					 *(tPos + 1), *(tPos + 3));
      *(tIsn + 2) = WlzContourItpTetSide(*(tVal + 2), *(tVal + 3),
					 *(tPos + 2), *(tPos + 3));
      break;
    default:
      isnCnt = 0;
      break;
  }
  if(isnCnt > 0)
  {
    /* Add cube origin */
    idx = isnCnt;
    tVP0 = tIsn;
    while(idx-- > 0)
    {
      tVP0->vtX += cbOrg.vtX;
      tVP0->vtY += cbOrg.vtY;
      tVP0->vtZ += cbOrg.vtZ;
      ++tVP0;
    }
    if(isnCnt == 4)
    {
      /* Split quadrilaterals into triangles along the shortest diagonal.
       * Choose shortest diagonal to try and avoid long thin triangles.
       * Know verticies to be ordered around the quadrilateral. */
      WLZ_VTX_3_SUB(tV0, tIsn[0], tIsn[2]);
      tD0 = WLZ_VTX_3_SQRLEN(tV0);
      WLZ_VTX_3_SUB(tV0, tIsn[1], tIsn[3]);
      tD1 = WLZ_VTX_3_SQRLEN(tV0);
      if(tD0 < tD1)
      {
        sIsn[0] = tIsn[0];
	sIsn[1] = tIsn[2];
	sIsn[2] = tIsn[3];
      }
      else
      {
        sIsn[0] = tIsn[1];
	sIsn[1] = tIsn[2];
	sIsn[2] = tIsn[3];
	tIsn[2] = tIsn[3];
      }
      if((errNum = WlzGMModelConstructSimplex3D(ctr->model,
      						tIsn)) == WLZ_ERR_NONE)
      {
        errNum = WlzGMModelConstructSimplex3D(ctr->model, sIsn);
      }
    }
    else
    {
      errNum = WlzGMModelConstructSimplex3D(ctr->model, tIsn);
    }
#ifdef WLZ_CONTOUR_DEBUG
    (void )fprintf(stderr,
		   "# %d I%d #%g %g %g,%g %g %g,%g %g %g\n",
		   isnCnt, iCode,
		   tIsn[0].vtX, tIsn[0].vtY, tIsn[0].vtZ,
		   tIsn[1].vtX, tIsn[1].vtY, tIsn[1].vtZ,
		   tIsn[2].vtX, tIsn[2].vtY, tIsn[2].vtZ);
    if(isnCnt == 4)
    {
      (void )fprintf(stderr,
		     "# %d I%d #%g %g %g,%g %g %g,%g %g %g\n",
		     isnCnt, iCode,
		     sIsn[0].vtX, sIsn[0].vtY, sIsn[0].vtZ,
		     sIsn[1].vtX, sIsn[1].vtY, sIsn[1].vtZ,
		     sIsn[2].vtX, sIsn[2].vtY, sIsn[2].vtZ);
    }
#endif /* WLZ_CONTOUR_DEBUG */
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourIntpTriSide
* Returns:	WlzDVertex2:		Position of intersection with
*					side.
* Purpose:	Calculates the position of the intersection of a
*		line with the zero height plane in 2D space.
* Global refs:	-
* Parameters:	double valOrg:		Height at line origin.
*		double valDst:		Height at line destination.
*		WlzDVertex2 posOrg:	Line origin.
*		WlzDVertex2 posDst:	Line destination.
************************************************************************/
static WlzDVertex2 WlzContourItpTriSide(double valOrg, double valDst,
				        WlzDVertex2 posOrg, WlzDVertex2 posDst)
{
  double	tD0;
  WlzDVertex2	itp;

  tD0 = valDst - valOrg;
  itp.vtX = ((valDst * posOrg.vtX)  - (valOrg * posDst.vtX)) / tD0;
  itp.vtY = ((valDst * posOrg.vtY)  - (valOrg * posDst.vtY)) / tD0;
  return(itp);
}

/************************************************************************
* Function:	WlzContourIntpTetSide
* Returns:	WlzDVertex3:		Position of intersection with
*					side.
* Purpose:	Calculates the position of the intersection of a
*		line with the zero height plane in 3D space.
* Global refs:	-
* Parameters:	double valOrg:		Height at line origin.
*		double valDst:		Height at line destination.
*		WlzDVertex3 posOrg:	Line origin.
*		WlzDVertex3 posDst:	Line destination.
************************************************************************/
static WlzDVertex3 WlzContourItpTetSide(double valOrg, double valDst,
				        WlzDVertex3 posOrg, WlzDVertex3 posDst)
{
  double	tD0;
  WlzDVertex3	itp;

  tD0 = valDst - valOrg;
  itp.vtX = ((valDst * posOrg.vtX)  - (valOrg * posDst.vtX)) / tD0;
  itp.vtY = ((valDst * posOrg.vtY)  - (valOrg * posDst.vtY)) / tD0;
  itp.vtZ = ((valDst * posOrg.vtZ)  - (valOrg * posDst.vtZ)) / tD0;
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
* Function:	WlzContourTestOutVTK
* Returns:	WlzErrorNum		Woolz error code.
* Purpose:	Prints out a 3D contour as a VTK polydata file
*		for testing.
* Global refs:	-
* Parameters:	WlzContour *ctr:	Given contour to print out.
*		FILE *fP:		Output file.
************************************************************************/
WlzErrorNum	WlzContourTestOutVTK(WlzContour *ctr, FILE *fP)
{
  WlzGMModel	*model;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctr->type != WLZ_CONTOUR)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((model = ctr->model) != NULL)
    {
      errNum = WlzGMModelTestOutVTK(model, fP);
    }
  }
  return(errNum);
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
  else if(ctr->type != WLZ_CONTOUR)
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
  		idN = 0,
		nCol,
		maxCol = 7,
		nmECnt,
		nmEdgeFlg = 0;
  double	ctrVal = 100,
  		ctrWth = 1.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  WlzGMEdge	*tE;
  WlzGMEdgeT	*tET;
  WlzGMEdge	**nmE;
  WlzObject     *inObj = NULL;
  WlzContourMethod ctrMtd = WLZ_CONTOUR_MTD_ISO;
  WlzContour	*ctr = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzDVertex3	tEGV[2];
  WlzDVertex2	scale,
  		offset;
  static char	optList[] = "eghic:o:v:w:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  nCol = maxCol;
  outFileStr = (char *)outFileStrDef;
  inObjFileStr = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'e':
        nmEdgeFlg = 1;
	break;
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
    ctr = WlzContourObj(inObj, ctrMtd, ctrVal, ctrWth, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to compute contour.\n",
      	     	     argv[0]);
    }
  }
  if(ok)
  {
    switch(ctr->model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	scale.vtX = 2.0;
	scale.vtY = -2.0;
	offset.vtX = -240.0;
	offset.vtY = 700.0;
        errNum = WlzContourTestOutPS(ctr, fP, nCol, 1, offset, scale);
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
        errNum = WlzContourTestOutVTK(ctr, fP);
	break;
      default: 
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s Failed to output contour to file (%d).\n",
      		     argv[0], (int )errNum);
    }
  }
  if(ok && nmEdgeFlg)
  {
    nmE = WlzGMModelFindNMEdges(ctr->model, &nmECnt, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s Failed to search for non-manifold edges (%d).\n",
		     argv[0], (int )errNum);
    }
    else
    {
      for(idN = 0; idN < nmECnt; ++idN)
      {
        tE = *(nmE + idN);
	tET = tE->edgeT;
	(void )WlzGMVertexGetG3D(tET->vertexT->diskT->vertex, &tEGV[0]);
	(void )WlzGMVertexGetG3D(tET->opp->vertexT->diskT->vertex, &tEGV[1]);
	(void )fprintf(stderr, "%d {%g, %g, %g} {%g, %g, %g}\n",
		       tE->idx,
		       (tEGV + 0)->vtX, (tEGV + 0)->vtY, (tEGV + 0)->vtZ,
		       (tEGV + 1)->vtX, (tEGV + 1)->vtY, (tEGV + 1)->vtZ);
      }
    }
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(ctr)
  {
    (void )WlzFreeContour(ctr);
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%sExample: %s%s",
      *argv,
      " [-o<output object>] [-h] [-o] [-g] [-i] [-c#] [-o#] [-v#] [-w#]\n"
      "        [<input object>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -o  Output object file name.\n"
      "  -g  Compute maximal gradient contours.\n"
      "  -i  Compute iso-value contours.\n"
      "  -c  Cycle contour colours (useful for debuging).\n"
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
