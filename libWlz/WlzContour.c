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
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

/* #define WLZ_CONTOUR_DEBUG */

typedef enum
{
  WLZ_CONTOURWS_CLEAR_NONE	= (0),
  WLZ_CONTOURWS_CLEAR_DATA	= (1),
  WLZ_CONTOURWS_CLEAR_FLAGS	= (1<<1)
} WlzContourWSpClearMode;

/************************************************************************
* WlzContourWSpace: Workspace for 2D or 3D contour generation.
************************************************************************/
typedef struct _WlzContourWSpace
{
  WlzContourType type;			        /* Type of contour: 2D or 3D */
  WlzContourMethod method;		   /* Iso-value or gradient contour. */
  int		conIdx;				  /* Next connectivity index */
  int		edgIdx;				          /* Next edge index */
  int		nodIdx;				          /* Next node index */
  int		newLn;    /* Non zero until edge created for line, then zero */
  int		minBlkSz; 	     /* Minimum number of elements per block */
  WlzEdge	*lstEdg;		     		        /* Last edge */
  WlzEdge	*edgCurLn;		       /* First edge on current line */
  WlzEdge	*edgLstLn;		          /* First edge on last line */
  AlcBlockStack	*nodStk;				 /* Node block stack */
  AlcBlockStack	*edgStk;				 /* Edge block stack */
  AlcHashTable	*conTable;		 	       /* Connectivity table */
} WlzContourWSpace;

/************************************************************************
* WlzContourTriIsn2D: Classification of intersection of line segment
* with a triangle
************************************************************************/
typedef enum
{
  WLZ_CONTOUR_TIC2D_NONE,	       			   /*No intersection */
  WLZ_CONTOUR_TIC2D_N1N0,	                          /* Node 1 - node 0 */
  WLZ_CONTOUR_TIC2D_N0N2,		                  /* Node 0 - node 2 */
  WLZ_CONTOUR_TIC2D_N2N1,			     	  /* Node 2 - node 1 */
  WLZ_CONTOUR_TIC2D_N1S02,			        /* Node 1 - side 0-2 */
  WLZ_CONTOUR_TIC2D_N0S21,			        /* Node 0 - side 2-1 */
  WLZ_CONTOUR_TIC2D_N2S10,			        /* Node 2 - side 1-0 */
  WLZ_CONTOUR_TIC2D_S10S02,			      /* Side 1-0 - side 0-2 */
  WLZ_CONTOUR_TIC2D_S02S21,			      /* Side 0-2 - side 2-1 */
  WLZ_CONTOUR_TIC2D_S21S10			      /* Side 2-1 - side 1-0 */
} WlzContourTriIsn2D;

static WlzContourWSpace *WlzMakeContourWSpace(WlzContourType type,
				WlzContourMethod method,
				WlzErrorNum *dstErr);
static WlzContourList *WlzContourWSpSegment(WlzContourWSpace *ctrWSp,
				int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContourList *WlzContourIsoObj2D(WlzObject *srcObj,
				double isoVal, int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContourList *WlzContourIsoObj3D(WlzObject *srcObj,
				double isoVal, int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContourList *WlzContourGrdObj2D(WlzObject *srcObj,
				double minGrd, double ftrPrm,
				int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContourList *WlzContourGrdObj3D(WlzObject *srcObj,
				double minGrd, double ftrPrm,
				int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContourList *WlzContourWSpSegment2D(WlzContourWSpace *ctrWSp,
				int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContourList *WlzContourWSpSegment3D(WlzContourWSpace *ctrWSp,
				int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzContour *WlzContourExtractContour2D(WlzContourWSpace *ctrWSp,
				WlzEdge *ctrEdg,
				int minNod, int minEdg,
				WlzErrorNum *dstErr);
static WlzEdge	*WlzContourWSpEdgNextAlc(WlzEdge *edg);
static WlzEdge	*WlzContourWSpEdgPrevAlc(WlzEdge *edg);
static WlzDVertex2 WlzContourItpTriSide(double ht0, double ht1,
				WlzDVertex2 org, WlzDVertex2 dst);
static WlzErrorNum WlzFreeContourWSpace(WlzContourWSpace *ctrWSp);
static WlzErrorNum WlzContourWSpAlcNod2D(WlzContourWSpace  *ctrWSp,
                                WlzEdgeNode2 **dstNod, int reqCnt);
static WlzErrorNum WlzContourWSpAlcEdg(WlzContourWSpace  *ctrWSp,
                                WlzEdge **dstNod, int reqCnt);
static WlzErrorNum WlzContourIsoCube2D(WlzContourWSpace *ctrWSp,
				double isoVal, double *ln0, double *ln1,
				WlzDVertex2 sqOrg);
static WlzErrorNum WlzContourWSpAddEdg2D(WlzContourWSpace *ctrWSp,
				WlzDVertex2 *iPos, WlzDVertex2 sqOrg);
static WlzErrorNum WlzContourGrdLink2D(WlzContourWSpace *ctrWSp,
			        UBYTE **grdDBuf,
				int bufWd, WlzIVertex2 org,
			        int lnOff, int lnIdx[], int klP);
static WlzErrorNum  WlzContourEdgNodCount(WlzEdge *ctrEdg,
				int *dstEdgCnt, int *dstNodCnt);
static WlzErrorNum WlzContourWSpClearFlags(WlzContourWSpace *ctrWSp,
				unsigned mode);
static void	WlzContourExtractContour2DFn(WlzEdge *edgBlk,
				int *edgIdx, WlzEdgeNode2 *nodBlk,
				int *nodIdx, WlzEdge *ctrEdg);
static void	WlzContourExtractContour2DFnNod(WlzEdgeNode2 *nodBlk,
				int *nodIdx, WlzEdgeNode2 *wNod);
static void	WlzContourExtractContour2DFnEdg(WlzEdge *edgBlk,
				int *edgIdx, WlzEdge *wEdg);
static void	WlzContourTestOutPS(FILE *,
				WlzContour *ctr);
static void	WlzContourTestOutPSFn(FILE *fP,
				WlzEdge *edg);
static void	WlzContourMatchEdg2D(WlzContourWSpace *ctrWSp,
				WlzEdge **nxtEdg, WlzDVertex2 *iPos,
				WlzDVertex2 sqOrg);
static void	WlzContourEdgNodCountFn(WlzEdge *edg,
				int *edgCnt, int *nodCnt);
static void 	WlzContourTestOutPSLn2D(FILE *fP, WlzDVertex2 org,
				WlzDVertex2 dst);
static int	WlzContourHashCtrKeyCmpFn(void *key0, void *key1);
static int	WlzContourWSpHashEntryIsNull(AlcHashTable *hTbl,
				AlcHashItem *hItem,
				void *dummy);
static int	WlzContourWSpEdgJoin2D(WlzContourWSpace *ctrWSp, 
				WlzEdge *conEdg0, WlzEdge *conEdg1);
static unsigned	WlzContourHashFn(void *keyP);

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
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ctrType != WLZ_CONTOUR_TYPE_2D) && (ctrType != WLZ_CONTOUR_TYPE_3D))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((ctr = AlcCalloc(1, sizeof(WlzContour))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    ctr->type = WLZ_CONTOUR;
    ctr->ctrType = ctrType;
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
      errNum = WlzFreeFreePtr(ctr->freeptr);
      AlcFree((void *)ctr);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzMakeContourList
* Returns:	WlzContourList *:		New contour list.
* Purpose:	Makes a new contour list data structure.
* Global refs:	-
* Parameters:	WlzErrorNum *dstErr:	Destination error pointer, may
					be null.
************************************************************************/
WlzContourList	*WlzMakeContourList(WlzErrorNum *dstErr)
{
  WlzContourList *ctrLst = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ctrLst = AlcCalloc(1, sizeof(WlzContourList))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    ctrLst->type = WLZ_CONTOUR_LIST;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctrLst);
}

/************************************************************************
* Function:	WlzFreeContourList
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Recursivly free's a WlzContourList data structure.
* Global refs:	-
* Parameters:	WlzContourList *ctrLst:	Given contour to free.
************************************************************************/
WlzErrorNum	WlzFreeContourList(WlzContourList *ctrLst)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctrLst == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctrLst->type != WLZ_CONTOUR_LIST)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if(WlzUnlink(&(ctrLst->linkcount), &errNum))
    {
      if(ctrLst->contour)
      {
        errNum = WlzFreeContour(ctrLst->contour);
      }
      if(ctrLst->next)
      {
        errNum = WlzFreeContourList(ctrLst->next);
      }
      if(ctrLst->next)
      {
        errNum = WlzFreeContourList(ctrLst->prev);
      }
      AlcFree((void *)ctrLst);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourClearFlags
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Clears the flags of all the contour's edges and
*		nodes.
* Global refs:	-
* Parameters:	WlzContour *ctr:	Given contour to free.
************************************************************************/
WlzErrorNum	WlzContourClearFlags(WlzContour *ctr)
{
  int		cnt;
  AlcBlockStack *blk;
  WlzEdge	*edg;
  WlzEdgeNode2	*nod2;
  WlzEdgeNode3	*nod3;
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
    /* Clear all edge flags. */
    blk = ctr->edgeStk;
    while(blk)
    {

      cnt = blk->elmCnt;
      edg = (WlzEdge *)(blk->elements);
      while(cnt-- > 0)
      {
        edg++->flags = WLZ_EDGE_FLAGS_NONE;
      }
      blk = blk->next;
    }
    /* Clear node flags. */
    blk = ctr->nodeStk;
    switch(ctr->ctrType)
    {
      case WLZ_CONTOUR_TYPE_2D:
	while(blk)
	{
	  cnt = blk->elmCnt;
	  nod2 = (WlzEdgeNode2 *)(blk->elements);
	  while(cnt-- > 0)
	  {
	    nod2++->flags = WLZ_EDGE_FLAGS_NONE;
	  }
	  blk = blk->next;
	}
	break;
      case WLZ_CONTOUR_TYPE_3D:
	while(blk)
	{
	  cnt = blk->elmCnt;
	  nod3 = (WlzEdgeNode3 *)(blk->elements);
	  while(cnt-- > 0)
	  {
	    nod3++->flags = WLZ_EDGE_FLAGS_NONE;
	  }
	  blk = blk->next;
	}
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzMakeContourWSpace
* Returns:	WlzContourWSpace *:		New contour.
* Purpose:	Makes a new contour data structure.
* Global refs:	-
* Parameters:	WlzContourType type:	Contour type.
*		WlzContourMethod method: Contour method.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContourWSpace *WlzMakeContourWSpace(WlzContourType type,
				       WlzContourMethod method,
				       WlzErrorNum *dstErr)
{
  WlzContourWSpace *ctrWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minBlkSz = 4096,
  		segTblSz = 4096;

  switch(type)
  {
    case WLZ_CONTOUR_TYPE_2D: /* FALLTHROUGH */
    case WLZ_CONTOUR_TYPE_3D:
      if((ctrWSp = AlcCalloc(1, sizeof(WlzContourWSpace))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	ctrWSp->type = type;
	ctrWSp->method = method;
  	ctrWSp->minBlkSz = minBlkSz;
	ctrWSp->conTable = AlcHashTableNew(segTblSz,
					   WlzContourHashCtrKeyCmpFn,
					   WlzContourHashFn, NULL);
        if(ctrWSp->conTable == NULL)
	{
	  AlcFree(ctrWSp);
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      break;
    default:
      errNum = WLZ_ERR_PARAM_TYPE;
      break;
  }
  return(ctrWSp);
}

/************************************************************************
* Function:	WlzFreeContourWSpace
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Free's a WlzContourWSpace data structure.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Contour workspace to free.
************************************************************************/
static WlzErrorNum WlzFreeContourWSpace(WlzContourWSpace *ctrWSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctrWSp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(ctrWSp->type)
    {
      case WLZ_CONTOUR_TYPE_2D: /* FALLTHROUGH */
      case WLZ_CONTOUR_TYPE_3D:
	if(ctrWSp->conTable)
	{
	  (void )AlcHashTableFree(ctrWSp->conTable);
	}
	AlcFree((void *)ctrWSp);
	break;
      default:
        errNum = WLZ_ERR_PARAM_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourObj
* Returns:	WlzContourList:		Contour list, or NULL on error.
* Purpose:	Creates a contour list from a Woolz object's values.
*		The source object should either a 2D or 3D domain
*		object with values.
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
WlzContourList	*WlzContourObj(WlzObject *srcObj, WlzContourMethod ctrMtd,
			       double ctrVal, double ctrWth,
			       int minNod, int minEdg,
			       WlzErrorNum *dstErr)
{
  WlzContourList *dstCtr = NULL;
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
	    dstCtr = WlzContourIsoObj2D(srcObj, ctrVal, minNod, minEdg,
	    				&errNum);
	    break;
	  case WLZ_CONTOUR_MTD_GRD:
	    dstCtr = WlzContourGrdObj2D(srcObj, ctrVal, ctrWth, minNod,
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
	    dstCtr = WlzContourIsoObj3D(srcObj, ctrVal, minNod, minEdg,
	    				&errNum);
	    break;
	  case WLZ_CONTOUR_MTD_GRD:
	    dstCtr = WlzContourGrdObj3D(srcObj, ctrVal, ctrWth, minNod, 
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
  return(dstCtr);
}

/************************************************************************
* Function:	WlzContourIsoObj2D
* Returns:	WlzContourList:		Contour list, or NULL on error.
* Purpose:	Creates an iso-value contour list from a 2D Woolz
*		object's values.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		double isoVal:		Iso-value.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContourList *WlzContourIsoObj2D(WlzObject *srcObj, double isoVal,
					  int minNod, int minEdg,
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
  WlzContourList *dstCtr = NULL;
  WlzDVertex2	sqOrg;
  WlzIntervalWSpace srcIWSp;
  WlzGreyWSpace	srcGWSp;
  WlzContourWSpace *ctrWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Create contour workspace. */
  ctrWSp = WlzMakeContourWSpace(WLZ_CONTOUR_TYPE_2D, WLZ_CONTOUR_MTD_ISO,
  			        &errNum);
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
	  errNum = WlzContourIsoCube2D(ctrWSp, isoVal,
				       valBuf[!bufLnIdx] + idX,
				       valBuf[bufLnIdx] + idX,
				       sqOrg);
	}
	++idX;
      }
      if(srcIWSp.intrmn == 0)
      {
        ctrWSp->newLn = 1;
	ctrWSp->edgLstLn = ctrWSp->edgCurLn;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstCtr = WlzContourWSpSegment(ctrWSp, minNod, minEdg, &errNum);
  }
  /* Free contour workspace. */
  if(ctrWSp)
  {
    (void )WlzFreeContourWSpace(ctrWSp);
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
  return(dstCtr);
}

/************************************************************************
* Function:	WlzContourIsoObj3D
* Returns:	WlzContourList:		Contour list, or NULL on error.
* Purpose:	Creates an iso-value contour list from a 3D Woolz
*		object's values.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object from which to
*					compute the contours.
*		double isoVal:		Iso-value.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContourList *WlzContourIsoObj3D(WlzObject *srcObj, double isoVal,
					  int minNod, int minEdg,
					  WlzErrorNum *dstErr)
{
  WlzContourList *dstCtr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* TODO */

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstCtr);
}

/************************************************************************
* Function:	WlzContourGrdObj2D
* Returns:	WlzContourList:		Contour list, or NULL on error.
* Purpose:	Creates an maximal gradient contour list from a 2D Woolz
*		object's values.
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
static WlzContourList *WlzContourGrdObj2D(WlzObject *srcObj,
					  double minGrd, double ftrPrm,
					  int minNod, int minEdg,
					  WlzErrorNum *dstErr)
{
  int		idX,
		dCode,
		iCnt,
		iLen,
		iBufSz,
		lnInc;
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
  WlzContourList *dstCtr = NULL;
  WlzContourWSpace *ctrWSp;
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

  /* Create contour workspace. */
  ctrWSp = WlzMakeContourWSpace(WLZ_CONTOUR_TYPE_2D, WLZ_CONTOUR_MTD_GRD,
  				&errNum);
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
	      /* Generate and link edge segments in the contour workspace. */
	      errNum = WlzContourGrdLink2D(ctrWSp, grdDBuf,
	      				   bufSz.vtX, org, posRel.vtY - 2,
					   lnIdx, klIdx[0]);
	    }
	  }
	  klIdx[0] = klIdx[1];
	  klIdx[1] = klIdx[2];
	  ++klIdx[2];
	}
      }
      ctrWSp->newLn = 1;
      ctrWSp->edgLstLn = ctrWSp->edgCurLn;
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstCtr = WlzContourWSpSegment(ctrWSp, minNod, minEdg, &errNum);
  }
  /* Free contour workspace. */
  if(ctrWSp)
  {
    (void )WlzFreeContourWSpace(ctrWSp);
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
  return(dstCtr);
}

/************************************************************************
* Function:	WlzContourGrdObj3D
* Returns:	WlzContourList:		Contour list, or NULL on error.
* Purpose:	Creates an maximal gradient contour list from a 3D Woolz
*		object's values.
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
static WlzContourList *WlzContourGrdObj3D(WlzObject *srcObj,
					  double minGrd, double ftrPrm,
					  int minNod, int minEdg,
					  WlzErrorNum *dstErr)
{
  WlzContourList *dstCtr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* TODO */

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstCtr);
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
* Parameters:	WlzContourWSp *ctrWSp:  Contour generation
*					workspace.
*		UBYTE **grdDBuf:	Buffers containing gradient
*					direction codes for maximal
*					gradient pixels.
*		int bufWd:		Width of the buffers.
*		WlzIVertex2 org:	Origin of the object.
*		int lnOff:		Offset from origin to central
*					pixel line.
*		int lnIdx[]:		Line indicies.
*		int klOff:		Column offset of first column.
************************************************************************/
static WlzErrorNum WlzContourGrdLink2D(WlzContourWSpace *ctrWSp,
				       UBYTE **grdDBuf,
				       int bufWd, WlzIVertex2 org,
				       int lnOff, int lnIdx[], int klOff)
{
  int		idN,
		idM,
		conCnt;
  int		con[5],
     		nbr[4];
  WlzDVertex2	nbrOrg;
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
      errNum = WlzContourWSpAddEdg2D(ctrWSp, seg, nbrOrg);
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
* Parameters:	WlzContourWSp *ctrWSp: Contour generation
					workspace.
*		double isoVal:		Iso-value to use.
*		double *vLn0:		Ptr to 2 data values at
*					y = yPos and x = xPos,
*					xpos + 1.
*		double *vLn1:		Ptr to 2 data values at
*					y = yPos + 1 and x = xPos,
*					xpos + 1..
*		WlzDVertex2 sqOrg:	The square's origin.
************************************************************************/
static WlzErrorNum WlzContourIsoCube2D(WlzContourWSpace *ctrWSp,
				       double isoVal,
				       double *vLn0, double *vLn1,
				       WlzDVertex2 sqOrg)
{
  int		idT,
		tN1,
  		tN2,
		dupIsn = 0;
  double	tD0;
  WlzContourTriIsn2D iCode;
  WlzDVertex2	tVx0,
  		tVx1;
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
  const WlzDVertex2 sNOff[5] =		  /* Square node horizontal offsets. */
  {
    {0.5, 0.5}, {0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0} /* {vtY, vtX} */
  };

  sV[1] = *vLn0 - isoVal; sV[2] = *(vLn0 + 1) - isoVal;
  sV[3] = *(vLn1 + 1) - isoVal; sV[4] = *vLn1 - isoVal;
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
	    tIsn[0] = sNOff[tN1];
	    tIsn[1] = sNOff[tN0];
	    break;
	  case WLZ_CONTOUR_TIC2D_N0N2:
	    tIsn[0] = sNOff[tN0];
	    tIsn[1] = sNOff[tN2];
	    break;
	  case WLZ_CONTOUR_TIC2D_N2N1:
	    tIsn[0] = sNOff[tN2];
	    tIsn[1] = sNOff[tN1];
	    break;
	  case WLZ_CONTOUR_TIC2D_N1S02:
	    tIsn[0] = sNOff[tN1];
	    tIsn[1] = WlzContourItpTriSide(tV[0], tV[2],
	    				   sNOff[tN0], sNOff[tN2]);
	    break;
	  case WLZ_CONTOUR_TIC2D_N0S21:
	    tIsn[0] = sNOff[tN0];
	    tIsn[1] = WlzContourItpTriSide(tV[2], tV[1],
	    				   sNOff[tN2], sNOff[tN1]);
	    break;
	  case WLZ_CONTOUR_TIC2D_N2S10:
	    tIsn[0] = sNOff[tN2];
	    tIsn[1] = WlzContourItpTriSide(tV[1], tV[0],
	    				   sNOff[tN1], sNOff[tN0]);
	    break;
	  case WLZ_CONTOUR_TIC2D_S10S02:
	    tIsn[0] = WlzContourItpTriSide(tV[1], tV[0],
	    				   sNOff[tN1], sNOff[tN0]);
	    tIsn[1] = WlzContourItpTriSide(tV[0], tV[2],
	    				   sNOff[tN0], sNOff[tN2]);
	    break;
	  case WLZ_CONTOUR_TIC2D_S02S21:
	    tIsn[0] = WlzContourItpTriSide(tV[0], tV[2],
	    				   sNOff[tN0], sNOff[tN2]);
	    tIsn[1] = WlzContourItpTriSide(tV[2], tV[1],
	    				   sNOff[tN2], sNOff[tN1]);
	    break;
	  case WLZ_CONTOUR_TIC2D_S21S10:
	    tIsn[0] = WlzContourItpTriSide(tV[2], tV[1],
	    				   sNOff[tN2], sNOff[tN1]);
	    tIsn[1] = WlzContourItpTriSide(tV[1], tV[0],
	    				   sNOff[tN1], sNOff[tN0]);
	    break;
	  default:
	    break;
	}
	tIsn[0].vtX += sqOrg.vtX; tIsn[0].vtY += sqOrg.vtY;
	tIsn[1].vtX += sqOrg.vtX; tIsn[1].vtY += sqOrg.vtY;
	if(tIsn[0].vtX > tIsn[1].vtX)
	{
	  tVx0 = tIsn[0];
	  tIsn[0] = tIsn[1];
	  tIsn[1] = tVx0;
	}
	/* Avoid duplicate edge segments along the triangle diagonals! */
	switch(idT)
	{
	  case 0:
	    tVx0.vtX = sqOrg.vtX;
	    tVx0.vtY = sqOrg.vtY;
	    tVx1.vtX = sqOrg.vtX + 0.5;
	    tVx1.vtY = sqOrg.vtY + 0.5;
	    break;
	  case 1:
	    tVx0.vtX = sqOrg.vtX + 0.5;
	    tVx0.vtY = sqOrg.vtY + 0.5;
	    tVx1.vtX = sqOrg.vtX + 1.0;
	    tVx1.vtY = sqOrg.vtY;
	    break;
	  case 2:
	    tVx0.vtX = sqOrg.vtX + 0.5;
	    tVx0.vtY = sqOrg.vtY + 0.5;
	    tVx1.vtX = sqOrg.vtX + 1.0;
	    tVx1.vtY = sqOrg.vtY + 1.0;
	    break;
	  case 3:
	    tVx0.vtX = sqOrg.vtX;
	    tVx0.vtY = sqOrg.vtY + 1.0;
	    tVx1.vtX = sqOrg.vtX + 0.5;
	    tVx1.vtY = sqOrg.vtY + 0.5;
	    break;
	}
	dupIsn = WlzGeomVtxEqual2D(tIsn[0], tVx0, tol) &&
		 WlzGeomVtxEqual2D(tIsn[1], tVx1, tol);
	if(!dupIsn)
	{
#ifdef WLZ_CONTOUR_DEBUG
	  (void )fprintf(stderr,
		         "%% {%f %f} T%d I%d {%g %g %g} {%g %g} {%g %g}\n",
		         sqOrg.vtX, sqOrg.vtY, idT, iCode,
		         tV[0], tV[1], tV[2],
		         tIsn[0].vtX, tIsn[0].vtY,
		         tIsn[1].vtX, tIsn[1].vtY);
#endif /* WLZ_CONTOUR_DEBUG */
	/* Add new contour element to the contour workspace. */
	  errNum = WlzContourWSpAddEdg2D(ctrWSp, tIsn, sqOrg);
#ifdef WLZ_CONTOUR_DEBUG
  	  WlzContourTestOutPSLn2D(stderr, tIsn[0], tIsn[1]);
#endif /* WLZ_CONTOUR_DEBUG */
	}
      }
      ++idT;
    }
  }
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
* Function:	WlzContourWSpAddEdg2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Adds a line segment to the contour workspace, linking
*		it's node to the workspace's node list and it's edge
*		to the the workspace's edge list.
*
*                                                         / |
*                                                         | |
*                     iPos[0]   newEdg[1]       nxtEdg[1]-| |
*		       |        |                         | /
*		      ### ------------------------------\ ###
*                     # #                                 # #
*		      ### \------------------------------ ###
*                     / |             |                    |
*		      | |-nxtEdg[0]   newEdg[0]           iPos[1]
*		      | |
*		      | /
*
*		Line segment is always such that:
*		  iPos[0].vtX <= iPos[1].vtX.
*
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Contour workspace.
*		WlzDVertex2 *iPos:	Positions of the intersections
*					with the triangle.
************************************************************************/
static WlzErrorNum WlzContourWSpAddEdg2D(WlzContourWSpace *ctrWSp,
				         WlzDVertex2 *iPos,
					 WlzDVertex2 sqOrg)
{
  int		conIdx,
  		done = 0;
  WlzEdge	*tEdg;
  WlzEdgeNode2	*newNod;
  WlzEdge	*newEdg = NULL;
  WlzEdge	*nxtEdg[2],
		*prvEdg[2];
  WlzDVertex2	diff;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	delta = 1.0e-6;

  WlzContourMatchEdg2D(ctrWSp, nxtEdg, iPos, sqOrg);
  if(nxtEdg[0] && nxtEdg[1])
  {
    /* Found edges directed from both the origin and destination nodes:
     * need 0 new nodes and 2 new edges. */
#ifdef WLZ_CONTOUR_DEBUG
   fprintf(stderr,
           "%% WlzContourWSpAddEdg2D "
   	   "e%d {%g %g} {%g %g} e%d {%g %g} {%g %g}\n",
	   (nxtEdg[0])->index,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->node))->pos.vtY,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->opp->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->opp->node))->pos.vtY,
	   (nxtEdg[1])->index,
	   ((WlzEdgeNode2 *)(nxtEdg[1]->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[1]->node))->pos.vtY,
	   ((WlzEdgeNode2 *)(nxtEdg[1]->opp->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[1]->opp->node))->pos.vtY);
#endif /* WLZ_CONTOUR_DEBUG */
    /* Check haven't found this edge segment already. */
    done = nxtEdg[0]->opp == nxtEdg[1];
    if(!done)
    {
      /* Create new edge joining two possibly unconnected contours. */
      errNum = WlzContourWSpAlcEdg(ctrWSp, &newEdg, 2);
      if(errNum == WLZ_ERR_NONE)
      {
	(newEdg + 0)->opp = (newEdg + 1);
	(newEdg + 1)->opp = (newEdg + 0);
	/* Set nodes. */
	(newEdg + 0)->node = nxtEdg[0]->opp->node;
	(newEdg + 1)->node = nxtEdg[1]->opp->node;
	/* Insert edge into existing edge lists. */
	prvEdg[0] = nxtEdg[1]->prev;
	prvEdg[1] = nxtEdg[0]->prev;
	(newEdg + 0)->next = nxtEdg[0];
	nxtEdg[0]->prev = (newEdg + 0);
	(newEdg + 0)->prev = prvEdg[0];
	prvEdg[0]->next = (newEdg + 0);
	(newEdg + 1)->next = nxtEdg[1];
	nxtEdg[1]->prev = (newEdg + 1);
	(newEdg + 1)->prev = prvEdg[1];
	prvEdg[1]->next = (newEdg + 1);
	/* Join two possibly unconnected edges. */
	conIdx = WlzContourWSpEdgJoin2D(ctrWSp, nxtEdg[0], nxtEdg[1]);
	(newEdg + 0)->flags = (newEdg + 1)->flags = conIdx;
      }
    }
  }
  else if(nxtEdg[0])
  {
    /* Only found origin node: need 1 new node and 2 new edges. */
#ifdef WLZ_CONTOUR_DEBUG
   fprintf(stderr,
           "%% WlzContourWSpAddEdg2D "
   	   "e%d {%g %g} {%g %g} NULL {- -} {- -}\n",
	   (nxtEdg[0])->index,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->node))->pos.vtY,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->opp->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[0]->opp->node))->pos.vtY);
#endif /* WLZ_CONTOUR_DEBUG */
    errNum = WlzContourWSpAlcEdg(ctrWSp, &newEdg, 2);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzContourWSpAlcNod2D(ctrWSp, &newNod, 1);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (newEdg + 0)->opp = (newEdg + 1);
      (newEdg + 1)->opp = (newEdg + 0);
      /* Set nodes. */
      (newNod + 0)->pos = iPos[1];
      (newEdg + 0)->node = nxtEdg[0]->opp->node;
      (newEdg + 1)->node.nod2 = newNod + 0;
      /* Insert edge into existing edge lists. */
      prvEdg[1] = nxtEdg[0]->prev;
      (newEdg + 0)->next = nxtEdg[0];
      nxtEdg[0]->prev = (newEdg + 0);
      (newEdg + 0)->prev = (newEdg + 1);
      (newEdg + 1)->next = (newEdg + 0);
      (newEdg + 1)->prev = prvEdg[1];
      prvEdg[1]->next = (newEdg + 1);
      /* Propagate the connectivity index of one of the existing
       * edges connected to the new edges. */
      conIdx = nxtEdg[0]->flags;
      (newEdg + 0)->flags = (newEdg + 1)->flags = conIdx;
    }
  }
  else if(nxtEdg[1])
  {
    /* Only found destination node: need 1 new node and 2 new edges. */
#ifdef WLZ_CONTOUR_DEBUG
   fprintf(stderr,
           "%% WlzContourWSpAddEdg2D "
   	   "NULL {- -} {- -} 0x%x {%g %g} {%g %g}\n",
	   (unsigned )(nxtEdg[1]),
	   ((WlzEdgeNode2 *)(nxtEdg[1]->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[1]->node))->pos.vtY,
	   ((WlzEdgeNode2 *)(nxtEdg[1]->opp->node))->pos.vtX,
	   ((WlzEdgeNode2 *)(nxtEdg[1]->opp->node))->pos.vtY);
#endif /* WLZ_CONTOUR_DEBUG */
    errNum = WlzContourWSpAlcEdg(ctrWSp, &newEdg, 2);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzContourWSpAlcNod2D(ctrWSp, &newNod, 1);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (newEdg + 0)->opp = (newEdg + 1);
      (newEdg + 1)->opp = (newEdg + 0);
      /* Set new node. */
      (newNod + 0)->pos = iPos[0];
      (newEdg + 0)->node.nod2 = newNod + 0;
      (newEdg + 1)->node = nxtEdg[1]->opp->node;
      /* Insert edge into existing edge lists. */
      prvEdg[0] = nxtEdg[1]->prev;
      (newEdg + 0)->next = (newEdg + 1);
      (newEdg + 0)->prev = prvEdg[0];
      prvEdg[0]->next = (newEdg + 0);
      (newEdg + 1)->next = nxtEdg[1];
      nxtEdg[1]->prev = (newEdg + 1);
      (newEdg + 1)->prev = (newEdg + 0);
      /* Propagate the connectivity index of one of the existing
       * edges connected to the new edges. */
      conIdx = nxtEdg[1]->flags;
      (newEdg + 0)->flags = (newEdg + 1)->flags = conIdx;
    }
  }
  else
  {
    /* Didn't find origin or destination nodes: need 2 new nodes and
     * 2 new edges. */
#ifdef WLZ_CONTOUR_DEBUG
   fprintf(stderr,
           "%% WlzContourWSpAddEdg2D "
   	   "NULL {- -} {- -} NULL {- -} {- -}\n");
#endif /* WLZ_CONTOUR_DEBUG */
    errNum = WlzContourWSpAlcEdg(ctrWSp, &newEdg, 2);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzContourWSpAlcNod2D(ctrWSp, &newNod, 2);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (newEdg + 0)->opp = (newEdg + 1);
      (newEdg + 1)->opp = (newEdg + 0);
      /* Set the two new nodes. */
      (newNod + 0)->pos = iPos[0];
      (newNod + 1)->pos = iPos[1];
      (newEdg + 0)->node.nod2 = newNod + 0;
      (newEdg + 1)->node.nod2 = newNod + 1;
      /* Link new edge. */
      (newEdg + 0)->next = (newEdg + 1);
      (newEdg + 0)->prev = (newEdg + 1);
      (newEdg + 1)->next = (newEdg + 0);
      (newEdg + 1)->prev = (newEdg + 0);
      /* Use edge flags to store connectivity index. */
      (newEdg + 0)->flags = (newEdg + 1)->flags = ctrWSp->conIdx;
      /* Have generated a new contour so add a entry to the connectivity
       * hash table: key = connectivity index, value = -1 */
      conIdx = (ctrWSp->conIdx)++;
      alcErr = AlcHashTableEntryInsert(ctrWSp->conTable,
      			      	       (void *)(newEdg + 0), NULL, NULL);
      if(alcErr != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(ctrWSp->newLn && newEdg)
  {
    ctrWSp->edgCurLn = newEdg;
    ctrWSp->newLn = 0;
  }
#ifdef WLZ_CONTOUR_DEBUG
  if(newEdg + 0)
  {
    fprintf(stderr,
	    "%% WlzContourWSpAddEdg2D 0 e%d n%d {%g %g} e%d e%d e%d\n",
	    (newEdg + 0)->index, 
	    (newEdg + 0)->node->index,
	    ((WlzEdgeNode2 *)((newEdg + 0)->node))->pos.vtX,
	    ((WlzEdgeNode2 *)((newEdg + 0)->node))->pos.vtY,
	    ((newEdg + 0)->next->index),
	    ((newEdg + 0)->prev->index),
	    ((newEdg + 0)->opp->index));
    if(newEdg + 1)
    {
      fprintf(stderr,
	      "%% WlzContourWSpAddEdg2D 1 e%d n%d {%g %g} e%d e%d e%d\n",
	      (newEdg + 1)->index, 
	      (newEdg + 1)->node->index,
	      ((WlzEdgeNode2 *)((newEdg + 1)->node))->pos.vtX,
	      ((WlzEdgeNode2 *)((newEdg + 1)->node))->pos.vtY,
	       ((newEdg + 1)->next->index),
	      ((newEdg + 1)->prev->index),
	      ((newEdg + 1)->opp->index));
    }
  }
  else
  {
    fprintf(stderr,
    	    "%% WlzContourWSpAddEdg2D None\n");
  }
#endif /* WLZ_CONTOUR_DEBUG */
  return(errNum);
}

/************************************************************************
* Function:	WlzContourMatchEdg2D
* Returns:	WlzEdge *:		Ptr to matching edge or NULL
*					if no matching edge found.

* Purpose:	Searches for edges directed to the nodes at the
*		given position iPos[0] and iPos[1] the edges opposite
*		to these edges, directed away from the matched edges,
*		would be the next edges for new edges directed to
*		these nodes:
*		  nxtEdg[0] directed from iPos[0] and
*		  nxtEdg[1] directed from iPos[1].
*		The given positions are such that:
*		  iPos[0].vtX <= iPos[1].vtX.
*		Because the data arive in a WLZ_RASTERDIR_ILIC order
*		(at least for the cube origin) and no edge segment
*		has a length greater than 1.0, its easy to do a fast
*		search.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		WlzEdge *nxtEdg:	Return for next edges.
*		WlzDVertex2 *iPos:	Positions of intersections.
*		WlzDVertex2 sqOrg:	Square origin.
************************************************************************/
static void	WlzContourMatchEdg2D(WlzContourWSpace *ctrWSp,
				     WlzEdge **nxtEdg, WlzDVertex2 *iPos,
				     WlzDVertex2 sqOrg)
{
  int		cnt,
		idF,
  		fndMsk,
		schMsk,
		lnIdx,
		stpIdx;
  double	tD0,
  		tD1,
		trX,
		trY,
		trCos,
		trSin,
		trScale;
  WlzDVertex2	tPos,
  		mPos,
		oPos,
		pos0,
		pos1,
  		off;
  WlzEdge	*edg0,
  		*edg1;
  WlzEdge	*fndEdg[2];
  const double	tol = 1.0e-06;

  *(nxtEdg + 0) = NULL;
  *(nxtEdg + 1) = NULL;
  if(ctrWSp->lstEdg)
  {
    fndMsk = 0;
    if(ctrWSp->newLn == 0)
    {
      switch(ctrWSp->method)
      {
        case WLZ_CONTOUR_MTD_ISO:
	  cnt = 14;
	  break;
        case WLZ_CONTOUR_MTD_GRD:
	  cnt = 4;
	  break;
      }
      edg0 = ctrWSp->lstEdg;
      do
      {
	pos0 = (edg0->node.nod2)->pos;
	if(WlzGeomVtxEqual2D(*(iPos + 0), pos0, tol))
	{
	  fndMsk |= 1;
	  *(fndEdg + 0) = edg0;
	}
	if(WlzGeomVtxEqual2D(*(iPos + 1), pos0, tol))
	{
	  fndMsk |= 2;
	  *(fndEdg + 1) = edg0;
	}
      }
      while((--cnt > 0) && (fndMsk != 3) &&
	    ((edg0 = WlzContourWSpEdgPrevAlc(edg0)) != NULL));
    }
    /* If either intersection is not found yet and one of the intersection
     * verticies lies close to the last line then: Update the last line
     * marker and search forward from it for matching nodes. */
    off.vtY = sqOrg.vtY + tol;
    /* Generate a search mask with bits set if vertex is close to
     * box line origin. */
    schMsk = ((((iPos + 0)->vtY < off.vtY) << 0) |
              (((iPos + 1)->vtY < off.vtY) << 1)) & ~fndMsk;
    if((fndMsk != 3) && (schMsk != 0) && ctrWSp->edgLstLn)
    {
      /* Advance the last line edge marker until the next (allocated) edge
       * would be in the box under this box. */
      edg0 = ctrWSp->edgLstLn;
      off.vtX = sqOrg.vtX - tol;
      stpIdx = (ctrWSp->newLn)? ctrWSp->edgIdx: ctrWSp->edgCurLn->index;
      while(edg0 && ((edg0->node.nod2->pos).vtX < off.vtX) &&
	    (edg0->index < stpIdx))
      {
        ctrWSp->edgLstLn = edg0;
	edg0 = WlzContourWSpEdgNextAlc(edg0);
      }
      /* Check all edges from the last line marker until the next (allocated)
       * edge is not in the box under this box. */
      edg0 = ctrWSp->edgLstLn;
      off.vtX = sqOrg.vtX + tol + 2.0;
      stpIdx = (ctrWSp->newLn)? ctrWSp->edgIdx: ctrWSp->edgCurLn->index;
      while(edg0 &&
            ((pos0 = edg0->node.nod2->pos).vtX < off.vtX) &&
	    (edg0->index < stpIdx) &&
            (schMsk != 0))
      {
        if(((schMsk & 1) != 0) &&
	   WlzGeomVtxEqual2D(*(iPos + 0), pos0, tol))
	{
	  fndMsk |= 1;
	  schMsk &= 2;
	  *(fndEdg + 0) = edg0;
	}
        if(((schMsk & 2) != 0) &&
	   WlzGeomVtxEqual2D(*(iPos + 1), pos0, tol))
	{
	  fndMsk |= 2;
	  schMsk &= 1;
	  *(fndEdg + 1) = edg0;
	}
	edg0 = WlzContourWSpEdgNextAlc(edg0);
      }
    }
    if(fndMsk)
    {
      /* Having found matching nodes now find next edges. */
      for(idF = 0; idF < 2; ++idF)
      {
        if(fndMsk & (1 << idF))
	{
	  edg0 = *(fndEdg+ idF);
	  edg1 = edg0->opp->prev;
	  if(edg0 == edg1)
	  {
	    /* Matched node is connected to only one other node, so the
	     * next node is just it's opposite. */
	    *(nxtEdg + idF) = edg0->opp;
#ifdef WLZ_CONTOUR_DEBUG
	    (void *)fprintf(stderr, "%% WlzContourMatchEdg2D %d e%d == e%d\n",
	    		    idF, edg0->index, edg1->index);
#endif /* WLZ_CONTOUR_DEBUG */
	  }
	  else
	  {
#ifdef WLZ_CONTOUR_DEBUG
	    (void *)fprintf(stderr, "%% WlzContourMatchEdg2D %d e%d\n",
	    		    idF, edg0->index);
#endif /* WLZ_CONTOUR_DEBUG */
	    /* Matched node is connected to more than one other node, so need
	     * to search for next edge. */
	    mPos = (*(fndEdg + idF))->node.nod2->pos;
	    oPos = *(iPos + ((idF)? 0: 1));
	    /* Compute the rigid body affine transform which takes the
	     * line segment from the matched node at mPos to the new node at
	     * oPos onto the x-axis and oriented from the origin. 
	     * The transform is of the form:
	     *   x' = scale * (x*cos(theta) - y*sin(theta) + dx)
	     *   y' = scale * (x*sin(theta) + y*cos(theta) + dy) */
	    tD0 = mPos.vtX - oPos.vtX;
	    tD1 = mPos.vtY - oPos.vtY;
	    trScale = (tD0 * tD0) + (tD1 * tD1);
	    if(trScale > tol)
	    {
	      trScale = 1.0 / sqrt(trScale);
	      trCos = (oPos.vtX - mPos.vtX);
	      trSin = (mPos.vtY - oPos.vtY);
	      trX = (mPos.vtX * mPos.vtX) + (mPos.vtY * mPos.vtY) -
		    (mPos.vtX * oPos.vtX) - (mPos.vtY * oPos.vtY);
	      trY = (mPos.vtX * oPos.vtY) - (mPos.vtY * oPos.vtX);
	      /* Search for next edge by walking around the matched
	       * node CCW, comparing the relative polar angles of the line
	       * segments between the nodes until the next CCW node has a 
	       * smaller polar angle than the current node. */
	      tPos = edg0->node.nod2->pos;
	      pos0.vtX = trScale *
	      		 ((tPos.vtX * trCos) - (tPos.vtY * trSin) + trX);
	      pos0.vtY = trScale *
	      		 ((tPos.vtX * trSin) - (tPos.vtY * trCos) + trY);
	      tPos = edg1->node.nod2->pos;
	      pos1.vtX = trScale *
	      		 ((tPos.vtX * trCos) - (tPos.vtY * trSin) + trX);
	      pos1.vtY = trScale *
	      		 ((tPos.vtX * trSin) - (tPos.vtY * trCos) + trY);
	      while(WlzGeomCmpAngle(pos1, pos0) > 0)
	      {
		edg0 = edg1;
		edg1 = edg1->opp->prev;
		pos0 = pos1;
		tPos = edg1->node.nod2->pos;
		pos1.vtX = trScale *
			   ((tPos.vtX * trCos) - (tPos.vtY * trSin) + trX);
	        pos1.vtY = trScale *
			   ((tPos.vtX * trSin) - (tPos.vtY * trCos) + trY);
	      }
	      *(nxtEdg + idF) = edg0->opp;
	    }
	  }
	}
      }
    }
  }
}

/************************************************************************
* Function:	WlzContourWSpEdgJoin2D
* Returns:	int:			Base connectivity index.
* Purpose:	Having found two nodes which join possibly unconnected
*		edges, find the base connectivity of the two given
*		edges. Set the value of the upper of these to that of
*		the lower and return the lower valued index.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		WlzEdge *conEdg0:	First connectivity labeled edge.
*		WlzEdge *conEdg1:	Other connectivity labeled edge.
************************************************************************/
static int	WlzContourWSpEdgJoin2D(WlzContourWSpace *ctrWSp, 
				       WlzEdge *conEdg0, WlzEdge *conEdg1)
{
  int		eIdx;
  WlzEdge	*bConEdg[2];
  AlcHashItem	*conItem[2];

  bConEdg[0] = conEdg0;
  bConEdg[1] = conEdg1;
  if(conEdg0->flags != conEdg1->flags)
  {
    eIdx = 0;
    while(eIdx < 2)
    {
      conEdg1 = bConEdg[eIdx];
      do
      {
	conEdg0 = conEdg1;
	conItem[eIdx] = AlcHashItemGet(ctrWSp->conTable, (void *)conEdg0,
				       NULL);
	conEdg1 = (WlzEdge *)(conItem[eIdx]->entry);
      }
      while(conEdg1 != NULL);
      bConEdg[eIdx] = conEdg0;
      ++eIdx;
    }
    if(bConEdg[0]->flags < bConEdg[1]->flags)
    {
      conItem[1]->entry = (void *)(bConEdg[0]); 
    }
    else if(bConEdg[1]->flags < bConEdg[0]->flags)
    {
      conItem[0]->entry = (void *)(bConEdg[1]);
      bConEdg[0] = bConEdg[1];
    }
  }
  return(bConEdg[0]->flags);
}

/************************************************************************
* Function:	WlzContourWSpAlcNod2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Allocates the requested number of new WlzEdgeNode2
*		nodes.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		WlzEdgeNode2 **dstNod:	Destination pointer for the
*					new nodes.
*		int reqCnt:		Number of new nodes.
************************************************************************/
static WlzErrorNum WlzContourWSpAlcNod2D(WlzContourWSpace  *ctrWSp,
					 WlzEdgeNode2 **dstNod, 
					 int reqCnt)
{
  int		blkCnt;
  WlzEdgeNode2	*nod;
  AlcBlockStack	*blk,
  		*nBlk;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((blk = ctrWSp->nodStk) == NULL) ||
     (blk->elmCnt >= (blk->maxElm - reqCnt)))
  {
    blkCnt = (reqCnt > ctrWSp->minBlkSz)? reqCnt: ctrWSp->minBlkSz;
    if((nBlk = AlcBlockStackNew(blkCnt, sizeof(WlzEdgeNode2),
    			        blk, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      ctrWSp->nodStk = nBlk;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    blk = ctrWSp->nodStk;
    *dstNod = nod = (WlzEdgeNode2 *)(blk->elements) + blk->elmCnt;
    blk->elmCnt += reqCnt;
    while(reqCnt-- > 0)
    {
      nod->index = (ctrWSp->nodIdx)++;
      nod->data = NULL;
      ++nod;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourWSpAlcEdg
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Allocates requested number of  new WlzEdge edges.
*		The allocated edges are returned with ordered indicies
*		(*dstEdg + i)->index < (*dstEdg + j)->index, iff i < j.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		WlzEdge **dstEdg:	Destination pointer for the
*					new edges.
*		int reqCnt:		Number of new edges.
************************************************************************/
static WlzErrorNum WlzContourWSpAlcEdg(WlzContourWSpace  *ctrWSp,
				       WlzEdge **dstEdg, 
				       int reqCnt)
{
  int		blkCnt;
  WlzEdge	*edg;
  AlcBlockStack	*blk,
  		*nBlk;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((blk = ctrWSp->edgStk) == NULL) ||
     (blk->elmCnt >= (blk->maxElm - reqCnt)))
  {
    blkCnt = (reqCnt > ctrWSp->minBlkSz)? reqCnt: ctrWSp->minBlkSz;
    if((nBlk = AlcBlockStackNew(blkCnt, sizeof(WlzEdge),
    			        blk, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      ctrWSp->edgStk = nBlk;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    blk = ctrWSp->edgStk;
    *dstEdg = edg = (WlzEdge *)(blk->elements) + blk->elmCnt;
    blk->elmCnt += reqCnt;
    ctrWSp->lstEdg = edg + reqCnt - 1;
    while(reqCnt-- > 0)
    {
      edg->index = (ctrWSp->edgIdx)++;
      edg->data = (void *)blk;
      ++edg;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourWSpEdgPrevAlc
* Returns:	WlzEdge *:		Previous edge allocated.
* Purpose:	Finds the edge allocated previous to the one given.
* Global refs:	-
* Parameters:	WlzEdge *edg:		Given edge.
************************************************************************/
static WlzEdge	*WlzContourWSpEdgPrevAlc(WlzEdge *edg)
{
  WlzEdge	*nEdg = NULL;
  AlcBlockStack *blk;

  if(edg)
  {
    if(edg > (WlzEdge *)(((AlcBlockStack *)(edg->data))->elements))
    {
      nEdg = edg - 1;
    }
    else
    {
      if((blk = ((AlcBlockStack *)(edg->data))->next) != NULL)
      {
	nEdg = (WlzEdge *)(blk->elements) + blk->elmCnt - 1;
      }
    }
  }
  return(nEdg);
}

/************************************************************************
* Function:	WlzContourWSpEdgNextAlc
* Returns:	WlzEdge *:		Next edge allocated.
* Purpose:	Finds the edge allocated after the one given.
* Global refs:	-
* Parameters:	WlzEdge *edg:		Given edge.
************************************************************************/
static WlzEdge	*WlzContourWSpEdgNextAlc(WlzEdge *edg)
{
  WlzEdge	*nEdg = NULL;
  AlcBlockStack *blk;

  if(edg)
  {
    if(edg < ((WlzEdge *)((blk = (AlcBlockStack *)(edg->data))->elements) +
    				 blk->elmCnt - 1))
    {
      nEdg = edg + 1;
    }
    else
    {
      if((blk = ((AlcBlockStack *)(edg->data))->prev) != NULL)
      {
	nEdg = (WlzEdge *)(blk->elements) + 0;
      }
    }
  }
  return(nEdg);
}

/************************************************************************
* Function:	WlzContourWSpSegment
* Returns:	WlzContourList *:	New contour list.
* Purpose:	Segments the given contour workspace into unconnected
*		contours.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContourList *WlzContourWSpSegment(WlzContourWSpace *ctrWSp,
					    int minNod, int minEdg,
					    WlzErrorNum *dstErr)
{
  WlzContourList *ctrLst = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ctrWSp == NULL) ||
     (ctrWSp->nodStk == NULL) || (ctrWSp->edgStk == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(ctrWSp->type)
    {
      case WLZ_CONTOUR_TYPE_2D:
        ctrLst = WlzContourWSpSegment2D(ctrWSp, minNod,  minEdg, &errNum);
	break;
      case WLZ_CONTOUR_TYPE_3D:
        ctrLst = WlzContourWSpSegment3D(ctrWSp, minNod,  minEdg, &errNum);
	break;
      default:
        errNum = WLZ_ERR_PARAM_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctrLst);
}

/************************************************************************
* Function:	WlzContourWSpSegment2D
* Returns:	WlzContourList *:	New contour list.
* Purpose:	Segments the given 2D contour workspace into
*		unconnected contours.
*		Remove all non-basic contour edges (ie edges which
*		refer to another edge) from the segmented contour
*		hash list.
*		Find number of contours.
*		Compute number of nodes and edges for each contour.
*		Allocate new contours.
*		Copy segmented contours to new contours.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContourList *WlzContourWSpSegment2D(WlzContourWSpace *ctrWSp,
					      int minNod, int minEdg,
				              WlzErrorNum *dstErr)
{
  int		ctrHeadCnt;
  WlzContour	*ctr;
  WlzEdge	*ctrEdg;
  AlcHashItem	*ctrItem;
  AlcHashItem	**ctrHead;
  WlzContourList *ctrLstItem,
  		*ctrLstTop = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Remove all non-basic contour edges from the connectivity hash list. */
  (void )AlcHashTableUnlinkAll(ctrWSp->conTable, WlzContourWSpHashEntryIsNull,
  			       NULL, 1);
  /* Clear the workspace edge and node flags/data. */
  (void *)WlzContourWSpClearFlags(ctrWSp, (WLZ_CONTOURWS_CLEAR_DATA |
					   WLZ_CONTOURWS_CLEAR_FLAGS));
#ifdef WLZ_CONTOUR_DEBUG
  (void )fprintf(stderr, "\n%%Contour Count = %d\n",
  		 AlcHashTableCount(ctrWSp->conTable, NULL));
#endif /* WLZ_CONTOUR_DEBUG */
  /* Copy each segmented contour from the workspace. */
  ctrItem = NULL;
  ctrHead = ctrWSp->conTable->table;
  ctrHeadCnt = ctrWSp->conTable->tableSz;
  while((errNum == WLZ_ERR_NONE) && (ctrHeadCnt > 0))
  {
    ctrItem = *ctrHead;
    while(ctrItem)
    {
      ctrEdg = (WlzEdge *)(ctrItem->key);
      ctr = WlzContourExtractContour2D(ctrWSp, ctrEdg, minNod, minEdg,
      				       &errNum);
      if((errNum == WLZ_ERR_NONE) && (ctr != NULL))
      {
	if((ctrLstItem  = WlzMakeContourList(&errNum)) != NULL)
	{
	  ctrLstItem->contour = ctr;
	  if(ctrLstTop)
	  {
	    ctrLstItem->next = ctrLstTop;
	    ctrLstTop->prev = ctrLstItem;
	  }
	  ctrLstTop = ctrLstItem;
	}
      }
      ctrItem = ctrItem->next;
    }   
    ++ctrHead;
    --ctrHeadCnt;
  } 
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctrLstTop);
}

/************************************************************************
* Function:	WlzContourExtractContour2D
* Returns:	WlzContour *:		New contour.
* Purpose:	Copies the contour which includes the given edge
*		from the workspace.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		WlzEdge *ctrEdg:	An edge of the contour in the
*					workspace.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContour *WlzContourExtractContour2D(WlzContourWSpace *ctrWSp,
				WlzEdge *ctrEdg, int minNod, int minEdg,
				WlzErrorNum *dstErr)
{
  int		eIdx,
  		nIdx,
		nEdg = 0,
  		nNod = 0;
  WlzEdge	*edgBlk = NULL;
  WlzEdgeNode2	*nodBlk = NULL;
  WlzContour	*newCtr = NULL;
  WlzErrorNum	errNum;

  /* Get number of edges and nodes in contour, also sets the node and
   * edge WLZ_EDGE_FLAGS_STATS flag. */
  errNum = WlzContourEdgNodCount(ctrEdg, &nEdg, &nNod);
  if((errNum == WLZ_ERR_NONE) && (nNod >= minNod) && (nEdg >= minEdg))
  {
    /* Make a new contour. */
    newCtr = WlzMakeContour(WLZ_CONTOUR_TYPE_2D, &errNum);
    /* Allocate edge and node space. */
    if(errNum == WLZ_ERR_NONE)
    {
      if(((newCtr->nodeStk = AlcBlockStackNew(nNod, sizeof(WlzEdgeNode2),
				              NULL, NULL)) == NULL) ||
         ((newCtr->edgeStk = AlcBlockStackNew(nEdg, sizeof(WlzEdge),
				              NULL, NULL)) == NULL))
      {
	if(newCtr->nodeStk)
	{
	  AlcBlockStackFree(newCtr->nodeStk);
	}
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    /* Push the edge and node blocks onto the contour's freestack. */
    if(errNum == WLZ_ERR_NONE)
    {
      newCtr->nEdges = nEdg;
      newCtr->edges = (WlzEdge *)(newCtr->edgeStk->elements);
      newCtr->freeptr = WlzPushFreePtr(newCtr->freeptr,
      				       (void *)(newCtr->edgeStk),
				       &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      edgBlk = (WlzEdge *)(newCtr->edgeStk->elements);
      newCtr->freeptr = WlzPushFreePtr(newCtr->freeptr, (void *)edgBlk,
				       &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      newCtr->nNodes = nNod;
      newCtr->nodes.nod2 = (WlzEdgeNode2 *)(newCtr->nodeStk->elements);
      newCtr->freeptr = WlzPushFreePtr(newCtr->freeptr,
      				       (void *)(newCtr->nodeStk),
				       &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nodBlk = (WlzEdgeNode2 *)(newCtr->nodeStk->elements);
      newCtr->freeptr = WlzPushFreePtr(newCtr->freeptr,
      				       (void *)nodBlk, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      eIdx = 0;
      nIdx = 0;
      /* Copy the contour, also sets node and edge WLZ_EDGE_FLAGS_COPY flags. */
      WlzContourExtractContour2DFn(edgBlk, &eIdx, nodBlk, &nIdx, ctrEdg);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(newCtr)
      {
	WlzFreeContour(newCtr);
	newCtr = NULL;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newCtr);
}

/************************************************************************
* Function:	WlzContourExtractContour2DFn
* Returns:	void
* Purpose:	Copies the edges and nodes of a segmented contour
*		from the contour work space to the new contour's,
*		edge and edge node blocks.
*		This is a recursive function which should only be
*		called from WlzContourExtractContour2D().
* Global refs:	-
* Parameters:	WlzEdge *edgBlk:	New edge block.
*		int *edgIdx:		Current index into the new edge
*					block.
*		WlzEdgeNode2 *nodBlk:	New edge node block.
*		int *nodIdx:		Current index into the new node
*					block.
*		WlzEdge *wEdg:		Given workspace contour edge.
************************************************************************/
static void	WlzContourExtractContour2DFn(WlzEdge *edgBlk,
			        int *edgIdx, WlzEdgeNode2 *nodBlk,
				int *nodIdx, WlzEdge *wEdg)
{
  WlzEdgeNode2	*eNod0,
  		*wNod0;
  WlzEdge	*eEdg0,
  		*wEdg0,
  		*wEdg1;

  /* Travel CCW around the loop that includes the given edge. */
  wEdg0 = wEdg;
  wEdg1 = wEdg->next;
  do
  {
    if((wEdg0->flags & WLZ_EDGE_FLAGS_COPY) == WLZ_EDGE_FLAGS_NONE)
    {
      if(wEdg0->data == NULL)
      {
        WlzContourExtractContour2DFnEdg(edgBlk, edgIdx, wEdg0);
      }
      eEdg0 = (WlzEdge *)(wEdg0->data);
      if(wEdg0->node.core->data == NULL)
      {
	WlzContourExtractContour2DFnNod(nodBlk, nodIdx, wEdg0->node.nod2);
      }
      eEdg0->node.nod2 = (WlzEdgeNode2 *)(wEdg0->node.nod2->data);
      if(wEdg0->next->data == NULL)
      {
        WlzContourExtractContour2DFnEdg(edgBlk, edgIdx, wEdg0->next);
      }
      eEdg0->next = (WlzEdge *)(wEdg0->next->data);
      if(wEdg0->prev->data == NULL)
      {
        WlzContourExtractContour2DFnEdg(edgBlk, edgIdx, wEdg0->prev);
      }
      eEdg0->prev = (WlzEdge *)(wEdg0->prev->data);
      if(wEdg0->opp->data == NULL)
      {
        WlzContourExtractContour2DFnEdg(edgBlk, edgIdx, wEdg0->opp);
      }
      eEdg0->opp = (WlzEdge *)(wEdg0->opp->data);
      wEdg0->flags |= WLZ_EDGE_FLAGS_COPY;
    }
    wEdg0 = wEdg1;
    wEdg1 = wEdg1->next;
  } while(wEdg0 != wEdg);
  /* Travel around the same loop recursively calling this function function
   * for each opposite edge which has not been visited yet. */
  wEdg0 = wEdg;
  wEdg1 = wEdg->next;
  do
  {
    if((wEdg0->opp->flags & WLZ_EDGE_FLAGS_COPY) == WLZ_EDGE_FLAGS_NONE)
    {
      WlzContourExtractContour2DFn(edgBlk, edgIdx, nodBlk, nodIdx, wEdg0->opp);
    }
    wEdg0 = wEdg1;
    wEdg1 = wEdg1->next;
  } while(wEdg0 != wEdg);
}

/************************************************************************
* Function:	WlzContourExtractContour2DFnNod
* Returns:	void
* Purpose:	Copies a node into the new contour's node block.
* Global refs:	-
* Parameters:	WlzEdgeNode2 *nodBlk:	New edge node block.
*		int nNod:		Number of nodes in node block.
*		int *nodIdx:		Current index into the new node
*					block.
*		WlzEdgeNode2 *wNod:	Node to copy.
************************************************************************/
static void	WlzContourExtractContour2DFnNod(WlzEdgeNode2 *nodBlk,
				int *nodIdx, WlzEdgeNode2 *wNod)
{
  WlzEdgeNode2 *eNod;

  eNod = nodBlk + *nodIdx;
  wNod->data = (void *)eNod;
  eNod->type = WLZ_EDGENODE_2D;
  eNod->index = *nodIdx;
  eNod->flags = WLZ_EDGE_FLAGS_NONE;
  eNod->data = NULL;
  eNod->pos = wNod->pos;
  ++(*nodIdx);
}

/************************************************************************
* Function:	WlzContourExtractContour2DFnEdg
* Returns:	void
* Purpose:	Copies an edge into the new contour's edge block.
* Global refs:	-
* Parameters:	WlzEdge *edgBlk:	New edge block.
*		int nEdg:		Number of edges in edge block.
*		int *edgIdx:		Current index into the new edge
*					block.
*		WlzEdgeNode2 *wEdg:	Edge to copy.
************************************************************************/
static void	WlzContourExtractContour2DFnEdg(WlzEdge *edgBlk,
				int *edgIdx, WlzEdge *wEdg)
{
  WlzEdge	*eEdg;

  eEdg  = edgBlk + *edgIdx;
  wEdg->data = (void *)eEdg;
  eEdg->index = *edgIdx;
  eEdg->flags = WLZ_EDGE_FLAGS_NONE;
  eEdg->data = NULL;
  eEdg->next = NULL;
  eEdg->prev = NULL;
  eEdg->opp = NULL;
  ++(*edgIdx);
}

/************************************************************************
* Function:	WlzContourEdgNodCount
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Counts the number of edges and nodes that form a
*		connected contour with the given edge.
*		All the node and edge flags must have been cleared
*		before calling this function.
* Global refs:	-
* Parameters:	WlzEdge *ctrEdg:	Given contour edge.
*		int *dstEdgCnt:		Destination pointer for number
*					of edges.
*		int *dstNodCnt:		Destination pointer for number
*					edge nodes.
************************************************************************/
static WlzErrorNum 	WlzContourEdgNodCount(WlzEdge *ctrEdg,
				      	      int *dstEdgCnt, int *dstNodCnt)
{
  int		edgCnt = 0,
  		nodCnt = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctrEdg == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(ctrEdg)
  {
    WlzContourEdgNodCountFn(ctrEdg, &edgCnt, &nodCnt);
    if(dstEdgCnt)
    {
      *dstEdgCnt = edgCnt;
    }
    if(dstNodCnt)
    {
      *dstNodCnt = nodCnt;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourEdgNodCountFn
* Returns:	void
* Purpose:	Counts the number of edges and nodes that form a
*		connected contour with the given edge.
*		This is a recursive function which should only be
*		called from WlzContourEdgNodCount().
* Global refs:	-
* Parameters:	WlzEdge *edg:		Given contour edge.
*		int *edgCnt:		Edge count.
*		int *nodCnt:		Node count.
************************************************************************/
static void	WlzContourEdgNodCountFn(WlzEdge *edg,
					int *edgCnt, int *nodCnt)
{
  WlzEdge 	*edg0,
  		*edg1;

  /* Travel CCW around the loop that includes the given edge. */
  edg0 = edg;
  edg1 = edg->next;
  do
  {
    if((edg0->node.core->flags & WLZ_EDGE_FLAGS_STATS) == WLZ_EDGE_FLAGS_NONE)
    {
      ++*nodCnt;
      edg0->node.core->flags |= WLZ_EDGE_FLAGS_STATS;
    }
    if((edg0->flags & WLZ_EDGE_FLAGS_STATS) == WLZ_EDGE_FLAGS_NONE)
    {
      edg0->flags |= WLZ_EDGE_FLAGS_STATS;
      ++*edgCnt;
    }
    edg0 = edg1;
    edg1 = edg1->next;
  } while(edg0 != edg);
  /* Travel around the same loop recursively calling this function for
   * each opposite edge which has not been visited yet. */
  edg0 = edg;
  edg1 = edg->next;
  do
  {
    if((edg0->opp->flags & WLZ_EDGE_FLAGS_STATS) == WLZ_EDGE_FLAGS_NONE)
    {
      WlzContourEdgNodCountFn(edg0->opp, edgCnt, nodCnt);
    }
    edg0 = edg1;
    edg1 = edg1->next;
  } while(edg0 != edg);
}

/************************************************************************
* Function:	WlzContourWSpClearFlags
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Clears the data and/or flags of all the nodes and edges
*		in the contour workspace.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		unsigned mode:		Clear mode: data and/or flags.
************************************************************************/
static WlzErrorNum WlzContourWSpClearFlags(WlzContourWSpace *ctrWSp,
					   unsigned mode)
{
  int		cnt;
  WlzEdge	*edg;
  WlzEdgeNode	nod;
  WlzEdgeNode2	*nod2;
  WlzEdgeNode3	*nod3;
  AlcBlockStack	*blk;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctrWSp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    /* Clear node flags. */
    blk = ctrWSp->nodStk;
    while(blk)
    {
      cnt = blk->elmCnt;
      nod.core = (WlzEdgeNodeCore *)(blk->elements);
      switch(nod.core->type)
      {
        case WLZ_EDGENODE_2D:
	  nod2 = nod.nod2;
	  while(cnt-- > 0)
	  {
	    if(mode & WLZ_CONTOURWS_CLEAR_DATA)
	    {
	      nod2->data = NULL;
	    }
	    if(mode & WLZ_CONTOURWS_CLEAR_FLAGS)
	    {
	      nod2->flags = WLZ_EDGE_FLAGS_NONE;
	    }
	    ++nod2;
	  }
	  break;
	case WLZ_EDGENODE_3D:
	  nod3 = nod.nod3;
	  while(cnt-- > 0)
	  {
	    if(mode & WLZ_CONTOURWS_CLEAR_DATA)
	    {
	      nod3->data = NULL;
	    }
	    if(mode & WLZ_CONTOURWS_CLEAR_FLAGS)
	    {
	      nod3->flags = WLZ_EDGE_FLAGS_NONE;
	    }
	    ++nod3;
	  }
	  break;
      }
      blk = blk->next;
    }
    /* Clear edge flags. */
    blk = ctrWSp->edgStk;
    while(blk)
    {
      cnt = blk->elmCnt;
      edg  = (WlzEdge *)(blk->elements);
      while(cnt-- > 0)
      {
	if(mode & WLZ_CONTOURWS_CLEAR_DATA)
	{
          edg->data = NULL;
	}
	if(mode & WLZ_CONTOURWS_CLEAR_FLAGS)
	{
          edg->flags = WLZ_EDGE_FLAGS_NONE;
	}
	++edg;
      }
      blk = blk->next;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzContourHashCtrKeyCmpFn
* Returns:	int:			Key comparison <, ==, > 0.
* Purpose:	Compares the two given hash keys.
* Global refs:	-
* Parameters:	void *key0:		First hash key.
*		void *key1:		Second key.
************************************************************************/
static int	WlzContourHashCtrKeyCmpFn(void *key0, void *key1)
{
  return(((WlzEdge *)key0)->flags - ((WlzEdge *)key1)->flags);
}

/************************************************************************
* Function:	WlzContourHashFn
* Returns:	unsigned:		Hash value.
* Purpose:	Generates hash values for the segmented contour hash
*		table. Values are generated using a cheap random
*		number generator.
* Global refs:	-
* Parameters:	void *keyP:		Key from which to generate hash
*					value.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static unsigned WlzContourHashFn(void *keyP)
{
  unsigned      keyU,
                hashVal;
  const unsigned mask0 = 0xaaaaaaaa,
                 mask1 = 0x55555555,
                 mask2 = 0x33333333,
                 prime0 = 16381,
                 prime1 = 10211,
                 prime2 = 13691;

  keyU = (unsigned )(((WlzEdge *)keyP)->flags);
  hashVal = ((keyU & mask0) * prime0) +
            ((keyU & mask1) * prime1) +
            ((keyU & mask2) * prime2);
  return(hashVal);
}

/************************************************************************
* Function:	WlzContourWSpHashEntryIsNull
* Returns:	int:			1 or 0.
* Purpose:	Tests the given item, if it's entry is non NULL then
*		returns 1 else returns 0.
* Global refs:	-
* Parameters:	AlcHashTable *hTbl:	Connectivity hash table.
*		AlcHashItem *hItem:	Given hash item.
*		void *dummy:		Unused data pointer.
************************************************************************/
static int	WlzContourWSpHashEntryIsNull(AlcHashTable *hTbl,
					     AlcHashItem *hItem,
					     void *dummy)
{
  int		tstVal;

  tstVal = (hItem && hItem->entry);
  return(tstVal);
}

/************************************************************************
* Function:	WlzContourWSpSegment3D
* Returns:	WlzContourList *:	New contour list.
* Purpose:	Segments the given 3D contour workspace into 
*		unconnected contours.
* Global refs:	-
* Parameters:	WlzContourWSpace *ctrWSp: Woolz contour workspace.
*		int minNod:		Minimum number of nodes.
*		int minEdg:		Minimum number of edges.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzContourList *WlzContourWSpSegment3D(WlzContourWSpace *ctrWSp,
					      int minNod, int minEdg,
				              WlzErrorNum *dstErr)
{
  WlzContourList *ctrLst = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ctrLst);
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
  		scaleY = -1.0,
  		offsetX = 0.0,
		offsetY = 840.0;

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
* Returns:	void
* Purpose:	Outputs the contour edge segments as postscript.
* Global refs:	-
* Parameters:	FILE *fP:		Output file.
*		WlzContour *ctr:	Given contour.
************************************************************************/
static void	WlzContourTestOutPS(FILE *fP, WlzContour *ctr)
{
  WlzContourTestOutPSFn(fP, ctr->edges);
}

/************************************************************************
* Function:	WlzContourTestOutPSFn
* Returns:	void
* Purpose:	Outputs the contour edge segments as postscript.
* Global refs:	-
* Parameters:	FILE *fP:		Output file.
*		WlzEdge *edg:		Given edge.
************************************************************************/
static void	WlzContourTestOutPSFn(FILE *fP, WlzEdge *edg)
{
  WlzEdge 	*edg0,
  		*edg1;

  /* Travel CCW around the loop that includes the given edge. */
  edg0 = edg;
  edg1 = edg->next;
  do
  {
    if((edg0->node.core->flags & WLZ_EDGE_FLAGS_COPY) ==
       WLZ_EDGE_FLAGS_NONE)
    {
      edg0->node.core->flags |= WLZ_EDGE_FLAGS_COPY;
    }
    if((edg0->flags & WLZ_EDGE_FLAGS_COPY) == WLZ_EDGE_FLAGS_NONE)
    {
      WlzContourTestOutPSLn2D(fP, edg0->node.nod2->pos,
			      edg0->opp->node.nod2->pos);
      edg0->flags |= WLZ_EDGE_FLAGS_COPY;
    }
    edg0 = edg1;
    edg1 = edg1->next;
  } while(edg0 != edg);
  /* Travel around the same loop recursively calling this function function
   * for each opposite edge which has not been visited yet. */
  edg0 = edg;
  edg1 = edg->next;
  do
  {
    if((edg0->opp->flags & WLZ_EDGE_FLAGS_COPY) == WLZ_EDGE_FLAGS_NONE)
    {
      WlzContourTestOutPSFn(fP, edg0->opp); 		       /* Recursive! */
    }
    edg0 = edg1;
    edg1 = edg1->next;
  } while(edg0 != edg);
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
		minNod = 16,
		minEdg = 16;
  double	ctrVal = 16,
  		ctrWth = 1.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  WlzObject     *inObj = NULL;
  WlzContourMethod ctrMtd = WLZ_CONTOUR_MTD_ISO;
  WlzContour	*ctr;
  WlzContourList *ctrLst = NULL,
  		*ctrItem;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  static char	optList[] = "ghic:e:n:o:v:w:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";
  const int	maxCol = 7;			 /* Colours to cycle though. */
  const int	colR[] = {0,   255, 0,   255, 0,   200, 0},
  		colG[] = {0,   0,   255, 200, 0,   0,   255},
		colB[] = {0,   0,   0,   0,   255, 255, 200};

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
      fprintf(stderr, "%s Failed to open output file %s.\n",
      	      argv[0], outFileStr);
    }
  }
  if(ok)
  {
    (void )fprintf(fP, "%%!\n"
    		   "1 setlinecap\n"
		   "1 setlinejoin\n"
		   "0.01 setlinewidth\n");
    ctrItem = ctrLst = WlzContourObj(inObj, ctrMtd, ctrVal, ctrWth,
    				     minNod, minEdg, &errNum);
    while((errNum == WLZ_ERR_NONE) && (ctrItem != NULL))
    {
      /* Set colour. */
      if(nCol > 1)
      {
	(void )fprintf(fP,
		       "%d %d %d setrgbcolor\n",
		       colR[colIdx], colG[colIdx], colB[colIdx]);
	colIdx = ++colIdx % nCol;
      }
      /* Output the contour. */
      WlzContourTestOutPS(fP, ctrItem->contour);
      ctrItem = ctrItem->next;
    }
    (void )fprintf(fP, "showpage\n");
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
