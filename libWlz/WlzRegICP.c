#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzRegICP.c
* Date:         November 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions to register two objects using an
*		implementation of the iterative closest point
*		algorithm.
*		See
*		  * Besel P. J. and McKay N. D. A method for
*   		    registration of 3-D shapes. PAMI, 14(2): 239-256,
*		    1992.
*		  * Zhang Z. Iterative point matching for
*		    registration of free-form curves and surfaces.
*		    International Journal of Computer Vision,
*		    13(2):119-152, 1994.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 22-12-00 bill Add initial affine transform parameter. Add normals.
*		Modify convergence test.
* 13-12-00 bill Change members of WlzVertex and WlzVertexP.
* 06-12-00 bill Remove special case fn WlzRegICPCompTransform2D() and
*		integrate code for both 2D and 3D into
*		WlzRegICPCompTransform().
* 30-11-00 bill Make changes for 3D least squares affine transforms.
************************************************************************/
#include <float.h>
#include <Wlz.h>

typedef struct _WlzRegICPWSp
{
  int		itr;
  double	maxDist;
  int		nT;
  int		nS;
  WlzVertexP	tVx;
  WlzVertexP	tNr;
  WlzVertexP 	sVx;
  WlzVertexP 	sNr;
  WlzVertexType vType;
  AlcKDTTree	*tTree;
  int		*sNN;
  double	*dist;
  double	*cost;
  int		*rank;
  int		nMatch;
  WlzVertexP	tMatchVx;
  WlzVertexP 	sMatchVx;
  double	sumDistCrnt;
  double	sumDistLast;
  WlzAffineTransform *tr;
} WlzRegICPWSp;

static void			WlzRegICPFindNN(
				  WlzRegICPWSp *wSp);
static void			WlzRegICPRankMatch(
				  WlzRegICPWSp *wSp);
static int			WlzRegICPRankCmpFnD(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
static WlzErrorNum		WlzRegICPCompTransform(
				  WlzRegICPWSp *wSp,
				  WlzTransformType trType,
				  int *dstConv);
static WlzErrorNum 		WlzRegICPCheckVerticies(
				  WlzVertexP *vData,
				  int *vCnt,
				  WlzVertexType *vType);
static WlzAffineTransform 	*WlzRegICPVerticies(
				  WlzVertexP tVx,
				  WlzVertexP tNr,
				  int tCnt,
				  WlzVertexP sVx,
				  WlzVertexP sNr,
				  int sCnt,
				  WlzVertexType vType,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv,
				  int *dstItr,
				  int maxItr,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzRegICPBuildTree(
				  WlzRegICPWSp *wSp);
static WlzVertex 		WlzRegICPCompCentroid(
				  WlzVertexType vType,
				  WlzVertexP vx,
				  int nVx);

/************************************************************************
* Function:	WlzRegICPObjs
* Returns:	WlzAffineTransform *:	Affine transform which brings
*					the two objects into register.
* Purpose:	Registers the two given objects using the iterative
* 		closest point algorithm. An affine transform is
*		computed, which when applied to the source object
*		takes it into register with the target object.
* Global refs:	-
* Parameters:	WlzObject *tObj:	The target object.
*		WlzObject *sObj:	The source object to be
*					registered with target object.
*		WlzAffineTransform *initTr: Initial affine transform
*					to be applied to the source
*					object prior to using the ICP
*					algorithm. May be NULL.
*		WlzTransformType trType: Required transform type.
*		int *dstConv:		Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
*		int *dstItr:		Destination ptr for the number
*					of iterations, may be NULL.
*		int maxItr:		Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
WlzAffineTransform *WlzRegICPObjs(WlzObject *tObj, WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv, int *dstItr, int maxItr,
				  WlzErrorNum *dstErr)
{
  int		idx,
		conv,
		itr;
  int		vCnt[2];
  WlzVertexType	vType[2];
  WlzObject	*objs[2];
  WlzVertexP	nData[2],
  		vData[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzAffineTransform *regTr = NULL;

  nData[0].v = nData[1].v = NULL;
  vData[0].v = vData[1].v = NULL;
  objs[0] = tObj; objs[1] = sObj;
  if((tObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      vData[idx] = WlzVerticiesFromObj(objs[idx], nData + idx, vCnt + idx,
      				       vType + idx, &errNum);
      ++idx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check type of verticies and convert to double if required. */
    errNum = WlzRegICPCheckVerticies(vData, vCnt, vType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
     regTr = WlzRegICPVerticies(vData[0], nData[0], vCnt[0],
     				vData[1], nData[1], vCnt[1],
     				vType[0], initTr,
				trType, &conv, &itr, maxItr, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstConv)
    {
      *dstConv = conv;
    }
    if(dstItr)
    {
      *dstItr = itr;
    }
  }
  for(idx = 0; idx < 2; ++idx)
  {
    if(vData[idx].v)
    {
      AlcFree(vData[idx].v);
    }
    if(nData[idx].v)
    {
      AlcFree(nData[idx].v);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/************************************************************************
* Function:	WlzRegICPCheckVerticies
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Checks the two sets of verticies and promotes the
*		vertex type to double (2 or 3D) if not already double
*		verticies.
* Global refs:	-
* Parameters:	WlzVertexP *vData:	The sets of verticies.
*		int *vCnt:		The numbers of verticies.
*		WlzVertexType *vType:	The types of verticies.
************************************************************************/
static WlzErrorNum WlzRegICPCheckVerticies(WlzVertexP *vData, int *vCnt,
					   WlzVertexType *vType)
{
  int		idx,
  		copiedFlg;
  int		dim[2],
  		type[2];
  WlzVertexP	tVP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find vertex set dimensions. */
  idx = 0;
  do
  {
    switch(vType[idx])
    {
      case WLZ_VERTEX_I2: /* FALLTROUGH */
      case WLZ_VERTEX_F2: /* FALLTROUGH */
      case WLZ_VERTEX_D2:
	dim[idx] = 2;
	break;
      case WLZ_VERTEX_I3: /* FALLTROUGH */
      case WLZ_VERTEX_F3: /* FALLTROUGH */
      case WLZ_VERTEX_D3:
	dim[idx] = 3;
	break;
      default:
	errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  } while((++idx < 2) && (errNum == WLZ_ERR_NONE));
  if((errNum == WLZ_ERR_NONE) && (dim[0] != dim[1]))
  {
    /* This is unlikly to work, should it be allowed?  */
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find vertex set primitive types. */
    idx = 0;
    do
    {
      switch(vType[idx])
      {
	case WLZ_VERTEX_I2: /* FALLTROUGH */
	case WLZ_VERTEX_I3:
	  type[idx] = 0;
	case WLZ_VERTEX_F2: /* FALLTROUGH */
	case WLZ_VERTEX_F3:
	  type[idx] = 1;
	case WLZ_VERTEX_D2:
	case WLZ_VERTEX_D3:
	  type[idx] = 2;
	  break;
      }
    } while(++idx < 2);
    idx = 0;
    do
    {
      if((type[0] == 0) || (type[0] == 1) ||
         (type[1] == 0) || (type[1] == 1))
      {
        if(dim[idx] == 2)
	{
	  if((tVP.d2 = (WlzDVertex2 *)
	                AlcMalloc(vCnt[idx] * sizeof(WlzDVertex2))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    copiedFlg = 0;
	    switch(vType[idx])
	    {
	      case WLZ_VERTEX_I2:
		WlzValueCopyIVertexToDVertex(tVP.d2, vData[idx].i2,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	      case WLZ_VERTEX_F2:
		WlzValueCopyFVertexToDVertex(tVP.d2, vData[idx].f2,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	    }
	    if(copiedFlg)
	    {
	      vType[idx] = WLZ_VERTEX_D2;
	      AlcFree(vData[idx].v);
	      vData[idx].d2 = tVP.d2;
	    }
	  }
	}
	else /* dim[idx] == 3 */
	{
	  if((tVP.d3 = (WlzDVertex3 *)
	                AlcMalloc(vCnt[idx] * sizeof(WlzDVertex3))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    copiedFlg = 0;
	    switch(vType[idx])
	    {
	      case WLZ_VERTEX_I3:
		WlzValueCopyIVertexToDVertex3(tVP.d3, vData[idx].i3,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	      case WLZ_VERTEX_F3:
		WlzValueCopyFVertexToDVertex3(tVP.d3, vData[idx].f3,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	    }
	    if(copiedFlg)
	    {
	      vType[idx] = WLZ_VERTEX_D3;
	      AlcFree(vData[idx].v);
	      vData[idx].d3 = tVP.d3;
	    }
	  }
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (++idx < 2));
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzRegICPVerticies
* Returns:	WlzAffineTransform *:	Affine transform which brings
*					the two sets of verticies into
*					register.
* Purpose:	Registers the two given sets of verticies using the
*		iterative closest point algorithm. An affine transform
*		is computed, which when applied to the source verticies
*		takes it into register with the target verticies.
*		The verticies and their normals are known to be either
*		WlzDVertex2 or WlzDVertex3.
* Global refs:	-
* Parameters:	WlzVertexP tVx:		Target verticies.
*		WlzVertexP tNr:		Target normals, may be NULL.
*		int tCnt:		Number of target verticies.
*		WlzVertexP sVx:		Source verticies.
*		WlzVertexP sNr:		Source normals, may be NULL.
*		int sCnt:		Number of source verticies.
*		WlzVertexType vType:	Type of the verticies.
*		WlzAffineTransform *initTr: Initial affine transform
*					to be applied to the source
*					object prior to using the ICP
*					algorithm. May be NULL.
*		WlzTransformType trType: Required transform type.
*		int *dstConv:		Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
*		int *dstItr:		Destination ptr for the number
*					of iterations, may be NULL.
*		int maxItr:		Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
*		WlzErrorNum *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
static WlzAffineTransform *WlzRegICPVerticies(WlzVertexP tVx, WlzVertexP tNr,
					      int tCnt,
					      WlzVertexP sVx, WlzVertexP sNr,
					      int sCnt,
					      WlzVertexType vType,
					      WlzAffineTransform *initTr,
					      WlzTransformType trType,
					      int *dstConv, int *dstItr,
					      int maxItr,
					      WlzErrorNum *dstErr)
{

  int		conv,
		maxCnt;
  WlzAffineTransform *regTr = NULL;
  WlzRegICPWSp	wSp;
  AlcKDTTree	*tTree = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  /* Setup workspace. */
  wSp.itr = 0;
  wSp.maxDist = DBL_MAX;
  wSp.nT = tCnt;
  wSp.vType = vType;
  wSp.tVx = tVx;
  wSp.tNr = tNr;
  wSp.nS = sCnt;
  wSp.sVx = sVx;
  wSp.sNr = sNr;
  wSp.tTree = NULL;
  wSp.sNN = NULL;
  wSp.dist = NULL;
  wSp.cost = NULL;
  wSp.rank = NULL;
  wSp.tMatchVx.v = wSp.sMatchVx.v = NULL;
  wSp.tr = NULL;
  wSp.sumDistCrnt = DBL_MAX;
  maxCnt = WLZ_MAX(tCnt, sCnt);
  if(((wSp.sNN = (int *)AlcMalloc(sizeof(int) * maxCnt)) == NULL) ||
     ((wSp.rank = (int *)AlcMalloc(sizeof(int) * maxCnt)) == NULL) ||
     ((wSp.dist = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL) ||
     ((wSp.cost = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(initTr)
    {
      wSp.tr = WlzAffineTransformCopy(initTr, &errNum);
    }
    else
    {
      if(vType == WLZ_VERTEX_D2)
      {
	wSp.tr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum);
      }
      else /* vType == WLZ_VERTEX_D3 */
      {
        wSp.tr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, &errNum);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(vType == WLZ_VERTEX_D2)
    {
      if(((wSp.tMatchVx.d2 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
      						       maxCnt)) == NULL) ||
         ((wSp.sMatchVx.d2 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
      						       maxCnt)) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    else /* vType == WLZ_VERTEX_D3 */
    {
      if(((wSp.tMatchVx.d3 = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) *
      						       maxCnt)) == NULL) ||
         ((wSp.sMatchVx.d3 = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) *
      						       maxCnt)) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Build k-D tree for the target verticies. */
    errNum = WlzRegICPBuildTree(&wSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute the centroid of the target verticies. */
    /* Iterate to find the registration transform. */
    conv = 0;
    do
    {
      wSp.sumDistLast = wSp.sumDistCrnt;
      /* Find closest points. */
      WlzRegICPFindNN(&wSp);
      /* Rank matched points. */
      WlzRegICPRankMatch(&wSp);
      /* Check for convergence. */
      /* Compute registration transform and check for convergence. */
      errNum = WlzRegICPCompTransform(&wSp, trType, &conv);
    } while((errNum == WLZ_ERR_NONE) && (wSp.itr++ < maxItr) && (conv == 0));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstConv)
    {
      *dstConv = conv;
    }
    if(dstItr)
    {
      *dstItr = wSp.itr;
    }
  }
  if(wSp.sNN)
  {
    AlcFree(wSp.sNN);
  }
  if(wSp.dist)
  {
    AlcFree(wSp.dist);
  }
  if(wSp.cost)
  {
    AlcFree(wSp.cost);
  }
  if(wSp.rank)
  {
    AlcFree(wSp.rank);
  }
  if(wSp.tMatchVx.v == NULL)
  {
    AlcFree(wSp.tMatchVx.v);
  }
  if(wSp.sMatchVx.v == NULL)
  {
    AlcFree(wSp.sMatchVx.v);
  }
  if(wSp.tr)
  {
    if((errNum == WLZ_ERR_NONE) && conv)
    {
      regTr = wSp.tr;
    }
    else
    {
      (void )WlzFreeAffineTransform(wSp.tr);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/************************************************************************
* Function:	WlzRegICPBuildTree
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Allocates and populates a k-D tree from the given
*		verticies. The verticies are either WlzDVertex2 or
*		WlzDVertex3.
* Global refs:	-
* Parameters:	WlzRegICPWSp *wSp:	ICP registration workspace.
************************************************************************/
static WlzErrorNum WlzRegICPBuildTree(WlzRegICPWSp *wSp)
{
  int		idx,
		sIdx,
  		treeDim;
  int		*shuffle = NULL;
  int		datI[3];
  double	datD[3];
  AlcKDTNode	*node;
  WlzVertexP	tVP;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp->vType == WLZ_VERTEX_D2)
  {
    treeDim = 2;
  }
  else /* wSp->vType == WLZ_VERTEX_D3 */
  {
    treeDim = 3;
  }
  /* Create tree. */
  if((wSp->tTree = AlcKDTTreeNew(ALC_POINTTYPE_DBL, treeDim, -1.0, wSp->nT,
  				 NULL)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )AlgShuffleIdx(wSp->nT, wSp->rank, 0);
  }
  /* Populate tree using a shuffle index to get the behaviour of a randomized
   * k-D tree, making sure that the indicies of nodes of the tree are not
   * shuffled too. */
  if(errNum == WLZ_ERR_NONE)
  {
    idx = 0;
    if(wSp->vType == WLZ_VERTEX_D2)
    {
      while((alcErr == ALC_ER_NONE) && (idx < wSp->nT))
      {
	sIdx = *(wSp->rank + idx);
	tVP.d2 = (wSp->tVx.d2 + sIdx);
	datD[0] = tVP.d2->vtX;
	datD[1] = tVP.d2->vtY;
	if((node = AlcKDTInsert(wSp->tTree, datD, NULL, &alcErr)) != NULL)
	{
	  node->idx = sIdx;
	}
        ++idx;
      }
    }
    else /* wSp->vType == WLZ_VERTEX_D3 */
    {
      while((alcErr == ALC_ER_NONE) && (idx < wSp->nT))
      {
	sIdx = *(wSp->rank + idx);
	tVP.d3 = (wSp->tVx.d3 + sIdx);
	datD[0] = tVP.d3->vtX;
	datD[1] = tVP.d3->vtY;
	datD[2] = tVP.d3->vtZ;
	if((node = AlcKDTInsert(wSp->tTree, datD, NULL, &alcErr)) != NULL)
	{
	  node->idx = sIdx;
	}
        ++idx;
      }
    }
    if(alcErr != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzRegICPFindNN
* Returns:	void
* Purpose:	Finds nearest neighbour matches in the target tree for
*		the source verticies.
* Global refs:	-
* Parameters:	WlzRegICPWSp *wSp:	ICP registration workspace.
************************************************************************/
static void	WlzRegICPFindNN(WlzRegICPWSp *wSp)
{
  int		idx;
  WlzDVertex2	tVD2;
  WlzDVertex3	tVD3;
  double	datD[3];
  AlcKDTNode	*node;

  wSp->sumDistCrnt = 0.0;
  for(idx = 0; idx < wSp->nS; ++idx)
  {
    if(wSp->vType == WLZ_VERTEX_D2)
    {
      tVD2 = WlzAffineTransformVertexD2(wSp->tr, *(wSp->sVx.d2 + idx), NULL);
      datD[0] = tVD2.vtX;
      datD[1] = tVD2.vtY;
    }
    else /* wSp->vType == WLZ_VERTEX_D3 */
    {
      tVD3 = WlzAffineTransformVertexD3(wSp->tr, *(wSp->sVx.d3 + idx), NULL);
      datD[0] = tVD3.vtX;
      datD[1] = tVD3.vtY;
      datD[2] = tVD3.vtZ;
    }
    node = AlcKDTGetNN(wSp->tTree, datD, wSp->maxDist, wSp->dist + idx, NULL);
    *(wSp->sNN + idx) = node->idx;
    wSp->sumDistCrnt += *(wSp->dist + idx);
  }
}

/************************************************************************
* Function:	WlzRegICPRankMatch
* Returns:	void
*					quality threshold.
* Purpose:	Ranks matched verticies by cost. The cost of a match
*		is given by
*		  cM = kD.mD + kN.mN
*		and
*		  cM - the cost of the match.
*		  kD - constant for distance cost.
*		  mD - distance between matched verticies / median
*		       of distance between matched verticies.
*		  kN - constant for normal product cost.
*		  mN - 1.0 - scalar product of the normals of the
*		       matched verticies.
*		  
* Global refs:	-
* Parameters:	WlzRegICPWSp *wSp:	ICP registration workspace.
************************************************************************/
static void	WlzRegICPRankMatch(WlzRegICPWSp *wSp)
{
  int		idx,
  		idxL,
		idxH;
  double	mD,
		mN,
  		diff,
  		mDist;
  WlzDVertex3	dist;
  WlzVertex	sV,
  		tV;
  const double	kD = 0.5,
  		kN = 0.5,
		maxCost = 1.0;

  /* Initialize the rank array. */
  for(idx = 0; idx < wSp->nS; ++idx)
  {
    *(wSp->rank + idx) = idx;
  }
  /* Permute the rank array according to distance. */
  (void )AlgHeapSortIdx(wSp->dist, wSp->rank, wSp->nS, 
  			WlzRegICPRankCmpFnD);
  /* Find the median of the distances between the matched verticies. */
  mDist = *(wSp->dist + *(wSp->rank + (wSp->nS / 2)));
  if(mDist < 1.0)
  {
    mDist = 1.0;
  }
  /* Re-initialize the rank array. */
  for(idx = 0; idx < wSp->nS; ++idx)
  {
    *(wSp->rank + idx) = idx;
  }
  /* Compute costs. */
  if(wSp->sNr.v && wSp->tNr.v)
  {
    /* Use normals as they exist. */
    if(wSp->vType == WLZ_VERTEX_D2)
    {
      for(idx = 0; idx < wSp->nS; ++idx)
      {
	mD = *(wSp->dist + idx) / mDist;
	sV.d2 = *(wSp->sNr.d2 + idx);
	tV.d2 = *(wSp->tNr.d2 + *(wSp->sNN + idx));
	mN = WLZ_VTX_2_DOT(sV.d2, tV.d2);
	*(wSp->cost + idx) = (kD * mD) + (kN * mN);
      }
    }
    else /* wSp->vType == WLZ_VERTEX_D3 */
    {
      for(idx = 0; idx < wSp->nS; ++idx)
      {
	mD = *(wSp->dist + idx) / mDist;
	sV.d3 = *(wSp->sNr.d3 + idx);
	tV.d3 = *(wSp->tNr.d3 + *(wSp->sNN + idx));
	mN = WLZ_VTX_2_DOT(sV.d3, tV.d3);
	*(wSp->cost + idx) = (kD * mD) + (kN * mN);
      }
    }
  }
  else
  {
    /* No normals. */
    for(idx = 0; idx < wSp->nS; ++idx)
    {
      mD = *(wSp->dist + idx) / mDist;
      *(wSp->cost + idx) = (kD * mD) + kN;
    }
  }
  /* Permute the rank array according to cost. */
  (void )AlgHeapSortIdx(wSp->cost, wSp->rank, wSp->nS,
  			WlzRegICPRankCmpFnD);
  /* Find the number of matches with a cost less than the maxCost, using
   * binary search for efficiency. */
  idxL = 0;
  idxH = wSp->nS - 1;
  do
  {
    idx = (idxH + idxL) / 2;
    if((diff = (*(wSp->cost + *(wSp->rank + idx)) -
               maxCost)) < -(DBL_EPSILON))
    {
      idxL = idx;
    }
    else
    {
      idxH = idx;
    }
  }
  while(diff && ((idxH - idxL) > 1));
  wSp->nMatch = idxH + 1;
}

/************************************************************************
* Function:	WlzRegICPRankCmpFnD
* Returns:	int:			+ve,0 or -ve.
* Purpose:	Compares indexed double values and returns a signed
*		integer for sorting using AlgHeapSortIdx().
*		The indicies to the values are sorted s.t. smaller
*		values come first.
* Global refs:	-
* Parameters:	void *data:		Used to pass the array of double
*					data.
*		int *idx:		Used to pass the array of
*					indicies.
*		int id0:		First index index.
*		int id1:		Second index index.
************************************************************************/
static int	WlzRegICPRankCmpFnD(void *data, int *idx, int id0, int id1)
{
  double	*valP,
  		cmpD;
  int		cmpI = 0;

  valP = (double *)data;
  cmpD = *(valP + *(idx + id0)) - *(valP + *(idx + id1));
  if(cmpD < 0.0)
  {
    cmpI = -1;
  }
  else if(cmpD > 0.0)
  {
    cmpI = +1;
  }
  return(cmpI);
}

/************************************************************************
* Function:	WlzRegICPCompTransform
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Computes an affine transform from matched verticies
*		and sets a convergence flag.
* Global refs:	-
* Parameters:	WlzRegICPWSp *wSp:	ICP registration workspace.
*		WlzTransformType trType: Required transform type.
*		int dstConv:		Dst ptr for convergence flag,
*					non-zero if converged ie the
*					current affine transform is an
*					identity transform. May NOT be
*					NULL.
************************************************************************/
static WlzErrorNum WlzRegICPCompTransform(WlzRegICPWSp *wSp,
				          WlzTransformType trType,
					  int *dstConv)
{
  int		idx,
  		rIdx,
		conv = 0;
  WlzAffineTransform *curTr = NULL,
  		 *newTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp->vType == WLZ_VERTEX_D2)
  {
    /* Copy verticies to the matched vertex buffers. */
    for(idx = 0; idx < wSp->nMatch; ++idx)
    {
      rIdx = *(wSp->rank + idx);
      *(wSp->tMatchVx.d2 + idx) = *(wSp->tVx.d2 + rIdx);
      *(wSp->sMatchVx.d2 + idx) = WlzAffineTransformVertexD2(wSp->tr,
      						*(wSp->sVx.d2 + rIdx), NULL);
    }
    /* Compute the affine transform. */
    newTr = WlzAffineTransformLSq2D(wSp->nMatch, wSp->sMatchVx.d2,
				    wSp->nMatch, wSp->tMatchVx.d2,
				    trType, &errNum);
  }
  else /* wSp->vType == WLZ_VERTEX_D3 */
  {
    /* Copy verticies to the matched vertex buffers. */
    for(idx = 0; idx < wSp->nMatch; ++idx)
    {
      rIdx = *(wSp->rank + idx);
      *(wSp->tMatchVx.d3 + idx) = *(wSp->tVx.d3 + rIdx);
      *(wSp->sMatchVx.d3 + idx) = WlzAffineTransformVertexD3(wSp->tr,
      				     		*(wSp->sVx.d3 + rIdx), NULL);
    }
    /* Compute the affine transform. */
    newTr = WlzAffineTransformLSq3D(wSp->nMatch, wSp->sMatchVx.d3,
				    wSp->nMatch, wSp->tMatchVx.d3,
				    trType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    conv = WlzAffineTransformIsIdentity(newTr, NULL);
    if(wSp->tr == NULL)
    {
      wSp->tr = newTr;
      newTr = NULL;
    }
    else if(conv == 0)
    {
      curTr = WlzAffineTransformProduct(wSp->tr, newTr, &errNum);
      WlzFreeAffineTransform(wSp->tr);
      wSp->tr = curTr;
      curTr = NULL;
    }
#ifdef WLZ_REGICP_DEBUG
      (void )fprintf(stderr, "WlzRegICP conv = %d\n", conv);
      (void )fprintf(stderr, "WlzRegICP wSp->itr %d\n", wSp->itr);
      (void )fprintf(stderr, "WlzRegICP wSp->tr\n");
      (void )fprintf(stderr, "mat: %-10g %-10g %-10g %-10g\n",
		      wSp->tr->mat[0][0], wSp->tr->mat[0][1],
		      wSp->tr->mat[0][2], wSp->tr->mat[0][3]);
      (void )fprintf(stderr, "     %-10g %-10g %-10g %-10g\n",
		      wSp->tr->mat[1][0], wSp->tr->mat[1][1],
		      wSp->tr->mat[1][2], wSp->tr->mat[1][3]);
      (void )fprintf(stderr, "     %-10g %-10g %-10g %-10g\n",
		      wSp->tr->mat[2][0], wSp->tr->mat[2][1],
		      wSp->tr->mat[2][2], wSp->tr->mat[2][3]);
      (void )fprintf(stderr, "     %-10g %-10g %-10g %-10g\n",
		      wSp->tr->mat[3][0], wSp->tr->mat[3][1],
		      wSp->tr->mat[3][2], wSp->tr->mat[3][3]);
      (void )fprintf(stderr, "\n");
#endif /* WLZ_REGICP_DEBUG */
  }
  *dstConv = conv;
  if(curTr)
  {
    WlzFreeAffineTransform(curTr);
  }
  if(newTr)
  {
    WlzFreeAffineTransform(newTr);
  }
  return(errNum);
}

/* #define WLZ_REGICP_TEST */
#ifdef WLZ_REGICP_TEST

/* Test main() for WlzRegICPObjs(). */

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
  		option,
		tstOut = 0,
  		ok = 1,
		usage = 0;
  FILE		*fP = NULL;
  char		*outObjFileStr;
  char		*inObjFileStr[2];
  WlzValues	nullVal;
  WlzDomain	outDom;
  WlzDomain	ctrDom[2];
  WlzObject	*outObj;
  WlzObject	*inObj[2],
  		*ctrObj[2];
  char		tstOutFileStr[64];
  WlzDVertex2	tstOutOffset2D,
  		tstOutScale2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "ho:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  inObj[0] = inObj[1] = NULL;
  ctrObj[0] = ctrObj[1] = NULL;
  nullVal.core = NULL;
  outDom.core = NULL;
  ctrDom[0].core = ctrDom[1].core = NULL;
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
      case 't':
        tstOut = 1;
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
  if((fP = (strcmp(outObjFileStr, "-")?
	   fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to open output file %s.\n",
      	      	     argv[0], outObjFileStr);
    }
  }
  if(ok)
  {
    idx = 0;
    while(ok && (idx < 2))
    {
      ctrDom[idx].ctr = WlzContourObj(inObj[idx], WLZ_CONTOUR_MTD_GRD,
      				      10.0, 1.0, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )fprintf(stderr, "%s Failed to compute contour.\n",
		       argv[0]);
      }
      ++idx;
    }
  }
  if(ok)
  {
    idx = 0;
    while(ok && (idx < 2))
    {
      ctrObj[idx] = WlzMakeMain(WLZ_CONTOUR, ctrDom[idx], nullVal, NULL, NULL,
      				&errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to make contour object (%s).\n",
		       *argv, errMsg);
      }
      ++idx;
    }
  }
  if(ok && tstOut)
  {
    idx = 0;
    tstOutOffset2D.vtX = 100.0;
    tstOutOffset2D.vtY = 100.0;
    tstOutScale2D.vtX = 1.0;
    tstOutScale2D.vtY = 1.0;
    while(ok && (idx < 2))
    {
      (void )sprintf(tstOutFileStr, "tstOut%d.ps", idx);
      if((fP = fopen(tstOutFileStr, "w")) == NULL)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to open test output file (%s)\n",
		       argv[0], tstOutFileStr);
      }
      if(ok)
      {
	errNum = WlzGMModelTestOutPS(ctrObj[idx]->domain.ctr->model, fP,
				     tstOutOffset2D, tstOutScale2D, 1);
        if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )fprintf(stderr, "%s Failed to write to test file (%s)\n",
	  		 argv[0], tstOutFileStr);
	}
      }
      if(fP)
      {
        fclose(fP);
	fP = NULL;
      }
      ++idx;
    }
  }
  if(ok)
  {
    outDom.t = WlzRegICPObjs(ctrObj[0], ctrObj[1], WLZ_TRANSFORM_2D_REG,
    			     NULL, NULL, 1000, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to register contours.\n",
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
  for(idx = 0; idx < 2; ++idx)
  {
    if(inObj[idx])
    {
      (void )WlzFreeObj(inObj[idx]);
    }
    if(ctrObj[idx])
    {
      (void )WlzFreeObj(ctrObj[idx]);
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
      "Usage: %s%sExample: %s%s",
      *argv,
      " [-o<output object>] [-h] [<input object 0>] [<input object 1>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -o  Output object file name.\n"
      "Computes a registration transform which registers the second of the\n"
      "given objects to the first using an ICP algorithm.\n"
      "list from the given input object.\n");
  }
  return(!ok);
}

#endif /* WLZ_REGICP_TEST */
