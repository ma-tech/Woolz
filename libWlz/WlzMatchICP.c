#pragma ident "MRC HGU $Id$"
/*!
* \file        	WlzMatchICP.c
* \author       Bill Hill
* \date         August 2001
* \version	$Id$
* \note
*		Copyright:
*		2001 Medical Research Council, UK.
*		All rights reserved.
* \par  Address:
*		MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* \brief	Functions to match objects, providing a set of
*		corresponding points, using ICP based registration.
* \ingroup	WlzTransform
* \todo         
* \bug
*/
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>
#include <string.h>

typedef struct _WlzMatchICPShell
{
  int				size;
  WlzGMShell			*shell;
  WlzAffineTransform 		*tr;
} WlzMatchICPShell;

typedef struct _WlzMatchICPShellListElm
{
  struct _WlzMatchICPShell	mShell;
  struct _WlzMatchICPShellListElm *next;
  struct _WlzMatchICPShellListElm *prev;
} WlzMatchICPShellListElm;

typedef struct _WlzMatchICPShellList
{
  struct _WlzMatchICPShellListElm *head;
  struct _WlzMatchICPShellListElm *tail;
} WlzMatchICPShellList;

typedef struct _WlzMatchICPCbData
{
  int		nNS;
  int		nFS;
  WlzErrorNum	errNum;
  struct _WlzMatchICPShellList *list;
} WlzMatchICPCbData;

typedef struct _WlzMatchICPTPPair2D
{
  WlzDVertex2	sVx;
  WlzDVertex2	tVx;
} WlzMatchICPTPPair2D;

static WlzErrorNum		WlzMatchICPRegShellLst(
				  AlcKDTTree *tTree,
				  WlzGMModel *tGM,
				  WlzGMModel *sGM,
				  WlzMatchICPShellList *sLst,
				  WlzAffineTransform *globTr,
				  WlzTransformType trType,
				  WlzVertexType vType,
				  int sgnNrm,
				  int nTV,
				  WlzVertexP tVx,
				  WlzVertexP tNr,
				  int nSV,
				  WlzVertexP sVx,
				  WlzVertexP sNr,
				  int *iBuf,
				  WlzVertexP tVBuf,
				  WlzVertexP sVBuf,
				  double *wBuf,
				  int maxItr,
				  double maxDisp,
				  double maxAng,
				  double maxDeform,
				  int minSpx,
				  WlzAffineTransform *initTr,
				  WlzRegICPUsrWgtFn usrWgtFn,
				  void *usrWgtData);
static void			WlzMatchICPShellListInsertList(
				  WlzMatchICPShellList *list0,
				  WlzMatchICPShellList *list1);
static WlzAffineTransform 	*WlzMatchICPRegModel(
				  AlcKDTTree *tTree,
				  WlzTransformType trType,
				  WlzVertexType vType,
				  int sgnNrm,
				  int nTV,
				  WlzVertexP tVx,
				  WlzVertexP tNr,
				  int nSV,
				  WlzVertexP sVx,
				  WlzVertexP sNr,
				  int *iBuf,
				  WlzVertexP tVBuf,
				  WlzVertexP sVBuf,
				  double *wBuf,
				  int maxItr,
				  WlzAffineTransform *initTr,
				  int *dstConv,
				  WlzRegICPUsrWgtFn usrWgtFn,
				  void *usrWgtData,
				  WlzErrorNum *dstErr);
static WlzAffineTransform 	*WlzMatchICPRegShell(
				  AlcKDTTree *tTree,
				  WlzGMModel *tGM,
				  WlzGMShell *sS,
				  WlzAffineTransform *globTr,
				  WlzTransformType trType,
				  WlzVertexType vType,
				  int sgnNrm,
				  int nTV,
				  WlzVertexP tVx,
				  WlzVertexP tNr,
				  int nSV,
				  WlzVertexP sVx,
				  WlzVertexP sNr,
				  int *iBuf,
				  WlzVertexP tVBuf,
				  WlzVertexP sVBuf,
				  double *wBuf,
				  int maxItr,
				  double maxDisp,
				  double maxAng,
				  double maxDeform,
				  WlzAffineTransform *initTr,
				  int *dstConv,
				  WlzRegICPUsrWgtFn usrWgtFn,
				  void *usrWgtData,
				  WlzErrorNum *dstErr);
static int			WlzMatchICPGetPoints(
				  AlcKDTTree *tTree,
				  WlzMatchICPShell *mS,
				  WlzVertexType vType,
				  WlzVertexP tVx,
				  WlzVertexP tNr,
				  WlzVertexP sNr,
				  WlzVertexP tMatch,
				  WlzVertexP sMatch,
				  int offset,
				  WlzErrorNum *dstErr);
static int			WlzMatchICPGetPoints2D(
				  AlcKDTTree *tTree,
				  WlzMatchICPShell *mS,
				  WlzDVertex2 *tVx,
				  WlzDVertex2 *tNr,
				  WlzDVertex2 *sNr,
				  WlzDVertex2 *tMatch,
				  WlzDVertex2 *sMatch,
				  WlzErrorNum *dstErr);
static int			WlzMatchICPRegCheckConv(
				  WlzGMShell *shell,
				  WlzAffineTransform *globTr,
				  WlzAffineTransform *tr,
				  double maxDisp,
				  double maxAng,
				  double maxDeform);
static int			WlzMatchICPRegCheckConv2D(
				  WlzGMShell *shell,
				  WlzAffineTransform *globTr,
				  WlzAffineTransform *tr,
				  double maxDisp,
				  double maxAng,
				  double maxDeform);
static int			WlzMatchICPRegCheckConv3D(
				  WlzGMShell *shell,
				  WlzAffineTransform *globTr,
				  WlzAffineTransform *tr,
				  double maxDisp,
				  double maxAng,
				  double maxDeform);
static WlzErrorNum		WlzMatchICPBreakShellCon(
				  WlzMatchICPCbData *cbData,
				  WlzGMModel *sGM,
				  WlzGMShell *sS,
				  int *iBuf);
static WlzErrorNum 		WlzMatchICPBreakShellCon2D(
				  WlzMatchICPCbData *cbData,
				  WlzGMModel *sGM,
				  WlzGMShell *sS,
				  int *iBuf);
static WlzErrorNum		WlzMatchICPBreakShellMid(
				  WlzMatchICPCbData *cbData,
				  WlzGMModel *sGM,
				  WlzMatchICPShell *sMS,
				  WlzVertexP sNr,
				  int *iBuf,
				  int minSpx);
static WlzErrorNum		WlzMatchICPBreakShellLoopTs(
				  WlzGMModel *sGM,
				  WlzGMShell *sS,
				  WlzVertexP sNr,
				  int minSpx);
static WlzErrorNum		WlzMatchICPBreakShellLop2D(
				  WlzGMModel *sGM,
				  WlzGMShell *sS,
				  WlzDVertex2 *sNr,
				  int minSpx);
static WlzErrorNum 		WlzMatchICPBreakShellCur2D(
				  WlzMatchICPCbData *cbData,
				  WlzGMModel *sGM,
				  WlzMatchICPShell *sMS,
				  WlzDVertex2 *sNr,
				  int *iBuf,
				  int minSpx);
static WlzErrorNum		WlzMatchICPFilterPts(
				  WlzVertexType vType,
				  WlzVertexP tMatch,
				  WlzVertexP sMatch,
				  int *nMatch,
				  WlzAffineTransform *globTr,
				  int nN,
				  double impThr);
static WlzErrorNum 		WlzMatchICPFilterPtsDisp2D(
				  WlzDVertex2 *tVx,
				  WlzDVertex2 *sVx,
				  int *nVx,
				  WlzAffineTransform *globTr,
				  int nN,
				  double impThr);
static WlzErrorNum		WlzMatchICPFilterPtsRmDup2D(
				  WlzDVertex2 *tVx,
				  WlzDVertex2 *sVx,
				  int *nVxP);
static int			WlzMatchICPTPPairSortFnD(
				  void *p0,
				  void *p1);
static WlzGMVertex 		*WlzMatchICPLoopTMid(
				  WlzGMLoopT *gLT);
static WlzGMVertex 		*WlzMatchICPLoopTMaxCurv2D(
				  WlzGMLoopT *gLT,
				  WlzDVertex2 *gNr,
				  int minEdg,
				  double *dstAngle);
static WlzGMVertex 		*WlzMatchICPLoopTMinCurv2D(
				  WlzGMLoopT *gLT,
				  WlzDVertex2 *gNr,
				  int minEdg,
				  double *dstAngle);
static WlzErrorNum 		WlzMatchICPRemoveVertices(
				  WlzMatchICPCbData *cbData,
				  WlzGMModel *sGM,
				  int *iBuf,
				  int nD);
static WlzErrorNum 		WlzMatchICPMoveListElmToQueue(
				  AlcCPQQueue *queue,
				  WlzMatchICPShellList *list);
static int			WlzMatchICPShellSzCmp(
				  const void *val0,
				  const void *val1);
static int			WlzMatchICPDblSortFnD(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
static double			WlzMatchICPWeightMatches2D(
				  WlzAffineTransform *curTr,
				  AlcKDTTree *tree,
				  WlzDVertex2 *tVx,
				  WlzDVertex2 *sVx,
				  WlzDVertex2 tMVx,
				  WlzDVertex2 sMVx,
				  WlzGMModel *tGM,
				  WlzGMModel *sGM,
				  double maxDisp,
				  int nScatter);
static double			WlzMatchICPWeightMatches3D(
				  WlzAffineTransform *curTr,
				  AlcKDTTree *tree,
				  WlzDVertex3 *tVx,
				  WlzDVertex3 *sVx,
				  WlzDVertex3 tMVx,
				  WlzDVertex3 sMVx,
				  WlzGMModel *tGM,
				  WlzGMModel *sGM,
				  double maxDisp,
				  int nScatter);
static void			WlzMatchICPShellCb(
				  WlzGMModel *model,
				  WlzGMElemP elm,
				  WlzGMCbReason reason,
				  void *data);
static WlzMatchICPShellList *	WlzMatchICPShellListNew(
				  void);
static WlzMatchICPShellListElm *WlzMatchICPShellListElmNew(
				  void);
static void			WlzMatchICPShellListElmFree(
				  WlzMatchICPShellListElm *elm);
static void			WlzMatchICPShellListFree(
				  WlzMatchICPShellList *list);
static void			WlzMatchICPShellListElmInsert(
				  WlzMatchICPShellList *list,
				  WlzMatchICPShellListElm *nElm);
static void			WlzMatchICPShellListElmAppend(
				  WlzMatchICPShellList *list,
				  WlzMatchICPShellListElm *nElm);
static void			WlzMatchICPShellListElmUnlink(
				  WlzMatchICPShellList *list,
				  WlzMatchICPShellListElm *elm);

/*!
* \return				Error code.
* \ingroup	WlzTransform
* \brief	Finds matching points in two objects using an ICP
*		based matching algorithm.
*		On return the match points are ranked by plausibility
*		with the most plausible matches first.
*		If either the target or source object is a contour
*		then it's model is modified in place. If this isn't
*		wanted then use WlzCopyObj() to copy the contours
*		before calling this function.
* \param	tObj			The target object.
* \param	sObj			The source object to be
*					matched with target object.
* \param	initTr			Initial affine transform
*					to be applied to the source
*					object prior to matching. May
*					be NULL.
* \param	dstNMatch		Destination pointer for the number
*					of match points found. Required.
* \param	dstTMatch		Destination pointer for the target
*					match points found. Required.
* \param	dstSMatch		Destination pointer for the source
*					match points found. Required.
* \param	maxItr			Maximum number of iterations.
* \param	minSpx			Minimum number of simplicies in
*					a contour shell for matching.
* \param	brkFlg			Controls the breaking of the source
*					shells. Possible values are:
*					  - 0
*					    Whole model registration, source
*					    shells are never broken.
*					  - 1
*					    Whole source shells are
*					    registered to the target model.
*					    Source shells are never broken.
*					  - 2
*					    Source shells are only ever
*					    broken by connectivity.
*					  - > 2
*					    Source shells are broken by
*					    connectivity and then near
*					    there midpoints brkFlg - 1
*					    times.
* \param	maxDisp			The maximum displacement to the
* 					geometry of a shell, from the global
*					affine transformed position for
*					an acceptable registration.
*					distance, for a registered shell
*					to be used for correspondence
*					computation.
* \param	maxAng			The maximum angle of rotation of a
* 					shell geometry with respect to the
* 					global affine transformed geometry for
* 					an acceptable registration.
* \param	maxDeform		The maximum deformation to the geometry
* 					of a shell, from the global affine
* 					transformed geometry, for an acceptable
* 					registration.
* \param	matchImpNN		Number match points in neighbourhood
*					when removing implausible match
*					points, must be \f$> 2\f$.
* \param	matchImpThr		Implausibility threshold which should
*					be \f$> 0\f$ but the useful range
*					is probably \f$[0.5-2.5]\f$. Higher
*					values allow more implausible matches
*					to be returned.
*/
WlzErrorNum	WlzMatchICPObjs(WlzObject *tObj, WlzObject *sObj,
				WlzAffineTransform *initTr,
				int *dstNMatch, WlzVertexP *dstTMatch,
				WlzVertexP *dstSMatch, int maxItr,
				int minSpx, int brkFlg,
				double maxDisp, double maxAng, 
				double maxDeform,
				int matchImpNN, double matchImpThr)
{
  WlzMatchICPWeightCbData cbData;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	nScatter = 5;

  if((tObj == NULL) | (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if (tObj->type != sObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((tObj->domain.core == NULL) || (sObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstNMatch == NULL) || (dstTMatch == NULL) || (dstSMatch == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(tObj->type)
    {
      case WLZ_CONTOUR:
	if((tObj->domain.ctr->model == NULL) ||
	   (sObj->domain.ctr->model == NULL))
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(tObj->domain.ctr->model->type !=
		sObj->domain.ctr->model->type)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  switch(tObj->domain.ctr->model->type)
	  {
	    case WLZ_GMMOD_2I: /* FALLTHROUGH */
	    case WLZ_GMMOD_2D: /* FALLTHROUGH */
	    case WLZ_GMMOD_2N: /* FALLTHROUGH */
	    case WLZ_GMMOD_3I: /* FALLTHROUGH */
	    case WLZ_GMMOD_3D: /* FALLTHROUGH */
	    case WLZ_GMMOD_3N:
	      /* Set up weighting function callback data. */
	      cbData.tGM = tObj->domain.ctr->model;
	      cbData.sGM = sObj->domain.ctr->model;
	      cbData.maxDisp = maxDisp;
	      cbData.nScatter = nScatter;
	      errNum = WlzMatchICPCtr(tObj->domain.ctr, sObj->domain.ctr,
	      			      initTr, maxItr, minSpx,
				      dstNMatch, dstTMatch, dstSMatch, brkFlg,
				      maxDisp, maxAng, maxDeform,
    				      matchImpNN, matchImpThr,
				      WlzMatchICPWeightMatches, &cbData);
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	}
        break;
      default:
        break;
    }
  }
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	WlzTransform
* \brief	Establishes matching points in two 2D contours using ICP
*		based registration algorithm. The contours have been
*		checked so that their models are neither NULL or
*		of different types.
*		The source contour's model is modified in place.
* \param	tCtr			The target contour.
* \param	sCtr			The source contour to be
*					matched with target contour.
* \param	initTr			Initial affine transform
*					to be applied to the source
*					object prior to matching. May
*					be NULL.
* \param	maxItr			Maximum number of iterations.
* \param	minSpx			Minimum number of simplicies in
*					a contour shell for matching.
* \param	dstNMatch		Destination pointer for the number
*					of match points found.
* \param	dstTMatch		Destination pointer for the target
*					match points found.
* \param	dstSMatch		Destination pointer for the source
*					match points found.
* \param	brkFlg			Controls the breaking of the source
*					shells. Possible values are:
*					  - 0
*					    Whole model registration, source
*					    shells are never broken.
*					  - 1
*					    Whole source shells are
*					    registered to the target model.
*					    Source shells are never broken.
*					  - 2
*					    Source shells are only ever
*					    broken by connectivity.
*					  - > 2
*					    Source shells are broken by
*					    connectivity and then near
*					    there midpoints brkFlg - 1
*					    times.
* \param	maxDisp			The maximum displacement to the
* 					geometry of a shell for an acceptable
* 					registration.
*					distance for a registered shell
*					to be used for correspondence
*					computation.
* \param	maxAng			The maximum angle of rotation of a
* 					shell geometry with respect to the
* 					global affine transformed geometry for
* 					an acceptable registration.
* \param	maxDeform		The maximum deformation to the geometry
* 					of a shell for an acceptable
* 					registration.
* \param	matchImpNN		Number match points in neighbourhood
*					when removing implausible match
*					points, must be \f$> 2\f$.
* \param	matchImpThr		Implausibility threshold which should
*					be \f$> 0\f$ but the useful range
*					is probably \f$[0.5-2.5]\f$. Higher
*					values allow more implausible matches
*					to be returned.
* \param	usrWgtFn		User supplied weighting function.
* \param	usrWgtData		User supplied weighting data.
*/
WlzErrorNum  	WlzMatchICPCtr(WlzContour *tCtr, WlzContour *sCtr,
			       WlzAffineTransform *initTr,
			       int maxItr, int minSpx,
			       int *dstNMatch, WlzVertexP *dstTMatch,
			       WlzVertexP *dstSMatch, int brkFlg,
			       double  maxDisp, double maxAng,
			       double maxDeform,
			       int matchImpNN, double matchImpThr,
			       WlzRegICPUsrWgtFn usrWgtFn, void *usrWgtData)
{
  int		idS,
		n0,
		n1,
  		nTV,
  		nSV,
		nOSS,
		maxVI,
		maxTVI,
		maxSVI,
		sgnNrm = 0,
		vSz,
		brkIdx,
		convFlg,
		nMatch = 0,
  	 	dbgFlg = 0;
  double	maxD,
  		meanD;
  int		*vIBuf = NULL;
  double	*wBuf = NULL;
  WlzGMModel	*tGM,
  		*sGM;
  WlzVertexP	tVx,
		tNr,
  		sVx,
		sNr,
		tVBuf,
		sVBuf;
  WlzVertexType	vType,
  		tVType;
  WlzTransformType trType;
  AlcKDTTree	*tTree = NULL;
  WlzAffineTransform *tTr,
  		*globTr = NULL;
  WlzGMShell 	*cSS,
  		*nSS;
  WlzMatchICPShell *sMS,
  		*sMSBuf = NULL;
  WlzMatchICPShellListElm *lElm0,
  		*lElm1;
  WlzMatchICPShellList *dSList = NULL;
  AlcCPQQueue	*sMSQueue = NULL;
  AlcCPQItem	*qTop;
  WlzMatchICPCbData cbData;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const		trSrcFlg = 0; 	/* This needs to be 0 for valid tie point
  				 * output. */

  tVBuf.v = sVBuf.v = NULL;
  tVx.v = tNr.v = sVx.v = sNr.v = NULL;
  if((tCtr == NULL) || (sCtr == NULL) ||
     ((tGM = tCtr->model) == NULL) || ((sGM = sCtr->model) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(tGM->type != sGM->type)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    switch(tGM->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
        break;
      case WLZ_GMMOD_2N: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
        sgnNrm = 1;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  /* Remove small shells from the models. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMFilterRmSmShells(tGM, minSpx);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMFilterRmSmShells(sGM, minSpx);
  }
  /* Get the vertices and normals from the target and source models. */
  if(errNum == WLZ_ERR_NONE)
  {
    tVx = WlzVerticesFromGM(tCtr->model, &tNr, NULL, &nTV, &tVType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sVx = WlzVerticesFromGM(sCtr->model, &sNr, NULL, &nSV, &vType, &errNum);
    if(tVType != vType)
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
    else if((nTV <= 0) || (nSV <= 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      switch(vType)
      {
        case WLZ_VERTEX_D2:
	  vSz = sizeof(WlzDVertex2);
	  trType = WLZ_TRANSFORM_2D_AFFINE;
	  break;
	case WLZ_VERTEX_D3:
	  vSz = sizeof(WlzDVertex3);
	  trType = WLZ_TRANSFORM_3D_AFFINE;
	  break;
        default:
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  break;
      }
    }
  }
  /* Create some buffers.
   * vIBuf:	A source vertex index buffer to allow registration to be
   *            between the target (using the kd-tree) and the indexed
   *		source vertices.
   * tVBuf:	Temporary buffer for target vertices.
   * sVBuf:	Temporary buffer for source vertices.
   * wBuf:	Temporary for nearest neighbour vertex weights.
   * sMSBuf:	Buffer with source shell pointers, shell sizes and affine
   *            transforms.
   */
  if(errNum == WLZ_ERR_NONE)
  {
    nOSS = sGM->res.shell.numElm;
    maxTVI = tGM->res.vertex.numIdx;
    maxSVI = sGM->res.vertex.numIdx;
    maxVI = WLZ_MAX(maxSVI, maxTVI);
    if(((vIBuf = (int *)AlcMalloc(maxVI * sizeof(int))) == NULL) ||
       ((tVBuf.v = AlcMalloc(maxVI * vSz)) == NULL) ||
       ((sVBuf.v = AlcMalloc(maxVI * vSz)) == NULL) ||
       ((wBuf = (double *)AlcMalloc(maxVI * sizeof(double))) == NULL) ||
       ((sMSBuf = (WlzMatchICPShell *)
       		  AlcCalloc(nOSS, sizeof(WlzMatchICPShell))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Build a kD-tree from the vertices of the the model. */
  if(errNum == WLZ_ERR_NONE)
  {
    tTree = WlzVerticesBuildTree(vType, nTV, tVx, vIBuf, &errNum);
  }
  /* Register the vertices of the source model to those of the target. */
  if(errNum == WLZ_ERR_NONE)
  {
    globTr = WlzAssignAffineTransform(
    	     WlzMatchICPRegModel(tTree, trType, vType, sgnNrm,
	     			 nTV, tVx, tNr, nSV, sVx, sNr,
				 vIBuf, tVBuf, sVBuf, wBuf,
				 maxItr, initTr, &convFlg,
				 NULL, NULL, &errNum), NULL);
    if((errNum == WLZ_ERR_NONE) && (convFlg == 0))
    {
      errNum = WLZ_ERR_ALG_CONVERGENCE;
    }
  }
  /* If only registering whole source model to whole target model
   * then transform the source model using the computed global
   * transform. */
  if((errNum == WLZ_ERR_NONE) && trSrcFlg && (nOSS > 0) && (brkFlg == 0))
  {
    (void )WlzAffineTransformGMModel(sGM, globTr, 0, &errNum);
  }
  /* Register each of the shells of the source model to the target model,
   * deleting any shells which fail to register from the source model and
   * putting all registered shells into match shell entries. */
  if((errNum == WLZ_ERR_NONE) && (nOSS > 0) && (brkFlg > 0))
  {
    nOSS = 0;
    sMS = sMSBuf;
    cSS = sGM->child;
    do
    {
      tTr = WlzAssignAffineTransform(
            WlzMatchICPRegShell(tTree, tGM, cSS, globTr, trType,
	    			vType, sgnNrm,
			        nTV, tVx, tNr, nSV, sVx, sNr,
			        vIBuf, tVBuf, sVBuf, wBuf,
			        maxItr, maxDisp, maxAng, maxDeform,
				globTr, &convFlg,
				usrWgtFn, usrWgtData,
				&errNum), NULL);
      nSS = cSS->next;
      if((errNum == WLZ_ERR_NONE) && (convFlg == 1))
      {
	++nOSS;
        sMS->tr = tTr;
	sMS->shell = cSS;
	sMS->size = WlzGMShellSimplexCnt(cSS);
	++sMS;
      }
      else
      {
        errNum = WlzGMModelDeleteS(sGM, cSS);
      }
      cSS = nSS;
    } while((errNum == WLZ_ERR_NONE) && cSS &&
            sGM->child && (cSS != sGM->child));
  }
  /* If only registering whole source shells to whole target model
   * then transform each of the source shells using the associated
   * transform. */
  if((errNum == WLZ_ERR_NONE) && trSrcFlg && (nOSS > 0) && (brkFlg == 1))
  {
    idS = 0;
    sMS = sMSBuf;
    while((errNum == WLZ_ERR_NONE) && (idS < nOSS))
    {
      (void )WlzAffineTransformGMShell(sMS->shell, sMS->tr);
      ++sMS;
      ++idS;
    }
  }
  /* Create a priority queue. */
  if((errNum == WLZ_ERR_NONE) && (nOSS > 0) && (brkFlg > 1))
  {
    if((sMSQueue = AlcCPQQueueNew(NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Create a list for decomposed shells. */
  if((errNum == WLZ_ERR_NONE) && (nOSS > 0) && (brkFlg > 1))
  {
    if((dSList = WlzMatchICPShellListNew()) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Decompose each of the source shells by removing all high connectivity
   * vertices from the shells and then registering each of the child shells
   * (with a simplex count above the theshold) to the target model. The
   * decomposed shells are entered into a priority queue which is maintained
   * so that the shells with the greatest simplex count are at the top of the
   * queue and have the highest priority. */
  if((errNum == WLZ_ERR_NONE) && (nOSS > 0) && (brkFlg > 1))
  {
    idS = 0;
    sMS = sMSBuf;
    do
    {
      tTr = sMS->tr;
      cSS = sMS->shell;
      cbData.nNS = cbData.nFS = 0;
      cbData.errNum = WLZ_ERR_NONE;
      cbData.list = dSList;
      errNum = WlzGMModelAddResCb(sGM, (WlzGMCbFn )WlzMatchICPShellCb,
				  &cbData);
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzMatchICPBreakShellCon(&cbData, sGM, cSS, vIBuf);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if((lElm0 = WlzMatchICPShellListElmNew()) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  lElm0->mShell.shell = cSS;
	  lElm0->mShell.size = 0;
	  lElm0->mShell.tr = NULL;
	  WlzMatchICPShellListElmInsert(dSList, lElm0);
	  if(dSList->head != NULL)
	  {
	    /* Pass through the list of decomposed shells, registering
	     * those which are above the size threshold to the target
	     * model and removing any small shells or shells that do not
	     * register from both the list and the source model. */
	    errNum = WlzMatchICPRegShellLst(tTree, tGM, sGM, dSList, globTr,
		trType, vType, sgnNrm, nTV, tVx, tNr, nSV, sVx, sNr,
		vIBuf, tVBuf, sVBuf, wBuf, maxItr,
		maxDisp, maxAng, maxDeform, minSpx, tTr,
		usrWgtFn, usrWgtData);
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzMatchICPMoveListElmToQueue(sMSQueue, dSList);
      }
      WlzGMModelRemResCb(sGM, WlzMatchICPShellCb, &cbData);
      ++idS;
      ++sMS;
    } while((errNum == WLZ_ERR_NONE) && (idS < nOSS));
  }
  /* Continue to decompose the source shells which are held in a priority
   * queue until either the break count (this is probably only useful for
   * debuging) or there are no more shells in the queue which have more
   * than three times the threshold number of simplices. */
  if((errNum == WLZ_ERR_NONE) && (nOSS > 0) && (brkFlg > 2))
  {
    brkIdx = 2;
    while((errNum == WLZ_ERR_NONE) &&
	  (brkIdx++ < brkFlg) &&
          ((qTop = AlcCPQItemUnlink(sMSQueue)) != NULL) &&
	  (qTop->entry != NULL) && (qTop->priority > minSpx * 3))
    {
      sMS = (WlzMatchICPShell *)(qTop->entry);
      AlcCPQItemFree(sMSQueue, qTop);
      tTr = sMS->tr;
      cSS = sMS->shell;
      cbData.nNS = cbData.nFS = 0;
      cbData.errNum = WLZ_ERR_NONE;
      cbData.list = dSList;
      errNum = WlzGMModelAddResCb(sGM, (WlzGMCbFn )WlzMatchICPShellCb,
      				  &cbData);
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzMatchICPBreakShellMid(&cbData, sGM, sMS,
					  sNr, vIBuf, minSpx);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if((lElm0 = WlzMatchICPShellListElmNew()) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {

        lElm0->mShell.shell = sMS->shell;
	lElm0->mShell.size = 0;
	lElm0->mShell.tr = NULL;
	AlcFree(sMS);
	WlzMatchICPShellListElmInsert(dSList, lElm0);
	if(dSList->head != NULL)
	{
	  /* Pass through the list of decomposed shells, registering
	   * those which are above the size threshold to the target
	   * model and removing any small shells or shells that do not
	   * register from both the list and the source model. */
	  errNum = WlzMatchICPRegShellLst(tTree, tGM, sGM, dSList, globTr,
	      trType, vType, sgnNrm, nTV, tVx, tNr, nSV, sVx, sNr,
	      vIBuf, tVBuf, sVBuf, wBuf, maxItr,
	      maxDisp, maxAng, maxDeform, minSpx, tTr,
	      usrWgtFn, usrWgtData);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzMatchICPMoveListElmToQueue(sMSQueue, dSList);
      }
      WlzGMModelRemResCb(sGM, WlzMatchICPShellCb, &cbData);
    }
  }
  /* Now a priority queue of decomposed shells each of which has more than
   * the threshold number of simplices and is registered to the target
   * model.  Transform each of the shells in the source model while
   * clearing the queue and building the list of matched points. */
  if((errNum == WLZ_ERR_NONE) && (nOSS > 0) && (brkFlg > 1))
  {
    n0 = 0;
    while((errNum == WLZ_ERR_NONE) &&
          ((qTop = AlcCPQItemUnlink(sMSQueue)) != NULL) &&
	  (qTop->entry != NULL))
    {
      sMS = (WlzMatchICPShell *)(qTop->entry);
      n1 = WlzMatchICPGetPoints(tTree, sMS, vType, tVx,
				tNr, sNr, tVBuf, sVBuf, n0, &errNum);
      if(trSrcFlg)
      {
        (void )WlzAffineTransformGMShell(sMS->shell, sMS->tr);
      }
      n0 += n1;
      AlcFree(sMS);
      AlcCPQItemFree(sMSQueue, qTop);

    }
    nMatch = n0;
  }
  /* Filter out poor matches. */
  if((errNum == WLZ_ERR_NONE) && (dbgFlg == 0) && (brkFlg > 1) && (nMatch > 0))
  {
    errNum = WlzMatchICPFilterPts(vType, tVBuf, sVBuf, &nMatch, globTr,
				  matchImpNN, matchImpThr);
  }
  /* Free the temporary storage that's been allocated. */
  if((errNum == WLZ_ERR_NONE) && (nMatch > 0))
  {
    /* Reallocate the match points. */
    if(((sVBuf.v = AlcRealloc(sVBuf.v, nMatch * vSz)) == NULL) ||
       ((tVBuf.v = AlcRealloc(tVBuf.v, nMatch * vSz)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum != WLZ_ERR_NONE) || (nMatch <= 0))
  {
    AlcFree(tVBuf.v);
    AlcFree(sVBuf.v);
    sVBuf.v = NULL;
    tVBuf.v = NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstNMatch)
    {
      *dstNMatch = nMatch;
    }
    if(dstSMatch)
    {
      *dstSMatch = sVBuf;
    }
    if(dstTMatch)
    {
      *dstTMatch = tVBuf;
    }
  }
  if(dSList)
  {
    lElm0 = dSList->head;
    while(lElm0)
    {
      lElm1 = lElm0->next;
      (void )WlzFreeAffineTransform(lElm0->mShell.tr);
      lElm0 = lElm1;
    }
    WlzMatchICPShellListFree(dSList);
  }
  AlcFree(sMSBuf);
  AlcFree(vIBuf);
  AlcFree(wBuf);
  AlcFree(tVx.v);
  AlcFree(tNr.v);
  AlcFree(sVx.v);
  AlcFree(sNr.v);
  WlzFreeAffineTransform(globTr);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Moves all match shells from the list to the priority
*		queue using the number of simplices in each shell as
*		it's priority.
* \param	queue			Match shell priority queue.
* \param	list			Match shell list.
*/
static WlzErrorNum WlzMatchICPMoveListElmToQueue(AlcCPQQueue *queue,
					      WlzMatchICPShellList *list)
{
  WlzMatchICPShell *qMS;
  WlzMatchICPShellListElm *lElm0,
  		*lElm1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  lElm0 = list->head;
  while((errNum == WLZ_ERR_NONE) && (lElm0 != NULL)) 
  {
    if(((qMS = (WlzMatchICPShell *)
               AlcCalloc(1, sizeof(WlzMatchICPShell))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    *qMS = lElm0->mShell;
    lElm1 = lElm0;
    lElm0 = lElm0->next;
    WlzMatchICPShellListElmFree(lElm1);
    if(AlcCPQEntryInsert(queue, qMS->size, qMS) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  list->head = list->tail = NULL;
  return(errNum);
}

/*!
* \return	Weight value in the range [0.0-1.0].
* \ingroup      WlzTransform
* \brief	Compute an the weighting for the matched pair of
* 		verticies.
*		This function is called from WlzRegICPTreeAndVertices().
* \param	vType			Type of vertices.
* \param	curTr			Current affine transform.
* \param	tree			Given kD-tree populated by the
*					target vertices such that the nodes of
*					the tree have the same indicies as the
*					given target vertices.
* \param	tVx			The target vertices.
* \param	sVx			The source vertices.
* \param	tMVx			The matched target vertex.
* \param	sMVx			The matched source vertex.
* \param	wVx			Distance weight.
* \param	wNr			Normal weight.
* \param	data			Data used to pass the geometric models.
*/
double		WlzMatchICPWeightMatches(WlzVertexType vType,
					 WlzAffineTransform *curTr,
					 AlcKDTTree *tree,
					 WlzVertexP tVx, WlzVertexP sVx,
					 WlzVertex tMVx, WlzVertex sMVx,
					 double wVx, double wNr,
					 void *data)
{
   double	wgt = 1.0;
   WlzMatchICPWeightCbData *cbData;

   if(data)
   {
     cbData = (WlzMatchICPWeightCbData *)data;
     switch(vType)
     {
       case WLZ_VERTEX_D2:
         wgt = WlzMatchICPWeightMatches2D(curTr, tree, tVx.d2, sVx.d2,
					  tMVx.d2, sMVx.d2,
					  cbData->tGM, cbData->sGM,
					  cbData->maxDisp, cbData->nScatter);
	 break;
       case WLZ_VERTEX_D3:
         wgt = WlzMatchICPWeightMatches3D(curTr, tree, tVx.d3, sVx.d3,
	 				  tMVx.d3, sMVx.d3,
					  cbData->tGM, cbData->sGM,
					  cbData->maxDisp, cbData->nScatter);
	 break;
     }
     wgt *= wVx * wNr;
   }
   return(wgt);
}

/*!
* \return	Error code.
* \ingroup	WlzTransform
* \brief	Pass through the given shell list, computing and setting
*		shell sizes within the list, registering shells all
*		above the size threshold to the target model and removing
*		the list elements of any small shells or shells that do not
*		register.
* \param	tTree			Given kD-tree populated by the
*					target vertices such that the
*					nodes of the tree have the same
*					indicies as the given target vertices
*					and normals.
* \param	tGM			Target geometric model.
* \param	sGM			Source geometric model.
* \param	sLst			The given list of shells.
* \param	globTr			Global affine transform.
* \param	trType			The required type of transform,
*					must be either WLZ_TRANSFORM_2D_REG,
*					or WLZ_TRANSFORM_2D_AFFINE.
* \param	vType			Type of vertices.
* \param	sgnNrm			Non zero if normals are consistant
*					in component sign.
* \param	nTV			Number of target vertices.
* \param	tVx			The target vertices.
* \param	tNr			The target normals.
* \param        nSV			Number of source vertices.
* \param        sVx 			The source vertices.
* \param	sNr			The source normals.
*		iBuf			Buffer with at least nSV ints.
* \param	tVBuf			A buffer with room for at least
*					nS vertices.
*		sVBuf			A buffer with room for at least
*					nS vertices.
* \param	wBuf			A buffer with room for at least
*					nS doubles.
* \param	maxItr			Maximum number of iterations.
* \param 	maxDisp			Maximum displacement.
* \param	maxAng			maximum angle (radians).
* \param	maxDeform		Maximum deformation.
* \param	minSpx			Minimum number of simplicies within a
* 					shell.
* \param	gInitTr			Initial affine transform, may be NULL.
* \param	usrWgtFn		User supplied weighting function.
* \param	usrWgtData		User supplied weighting data.
*/
static WlzErrorNum	WlzMatchICPRegShellLst(AlcKDTTree *tTree,
				WlzGMModel *tGM,
				WlzGMModel *sGM,
				WlzMatchICPShellList *sLst,
				WlzAffineTransform *globTr,
				WlzTransformType trType,
				WlzVertexType vType, int sgnNrm,
				int nTV, WlzVertexP tVx, WlzVertexP tNr,
				int nSV, WlzVertexP sVx, WlzVertexP sNr,
				int *iBuf, WlzVertexP tVBuf, WlzVertexP sVBuf,
				double *wBuf, int maxItr, 
				double maxDisp, double maxAng, 
				double maxDeform, int minSpx,
				WlzAffineTransform *gInitTr,
				WlzRegICPUsrWgtFn usrWgtFn, void *usrWgtData)
{
  int		convFlg,
  		remFlg;
  WlzMatchICPShellListElm *lElm0,
  		*lElm1;
  WlzAffineTransform *initTr = NULL,
  		*tTr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  lElm0 = sLst->head;
  /* Copy the given initial affine transform as it may well belong to one of
   * the members of the list. */
  initTr = WlzAffineTransformCopy(gInitTr, &errNum);
  while((errNum == WLZ_ERR_NONE) && (lElm0 != NULL))
  {
    tTr = NULL;
    remFlg = 1;
    lElm0->mShell.size = WlzGMShellSimplexCnt(lElm0->mShell.shell);
    if(lElm0->mShell.size >= minSpx)
    {
      tTr = WlzAssignAffineTransform(
	    WlzMatchICPRegShell(tTree, tGM, lElm0->mShell.shell,
			        globTr, trType, vType, sgnNrm,
			        nTV, tVx, tNr, nSV, sVx, sNr,
			        iBuf, tVBuf, sVBuf, wBuf,
			        maxItr, maxDisp, maxAng, maxDeform,
				initTr, &convFlg,
				usrWgtFn, usrWgtData,
				&errNum),
	    NULL);
      remFlg = !convFlg;
      if(errNum == WLZ_ERR_ALG_CONVERGENCE)
      {
	/* Convergence flag will trap failures to register. */
        errNum = WLZ_ERR_NONE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      lElm1 = lElm0->next;
      if(remFlg)
      {
	/* Either this shell is too small or it has failed to
	 * register with the target model: Remove the element from
	 * the list. */
	WlzMatchICPShellListElmUnlink(sLst, lElm0);
	WlzGMModelDeleteS(sGM, lElm0->mShell.shell);
	WlzMatchICPShellListElmFree(lElm0);
	WlzFreeAffineTransform(tTr);
      }
      else
      {
	/* Shell is larger than threshold size and has been
	 * registered with the target model. */
	lElm0->mShell.tr = tTr;
      }
      lElm0 = lElm1;
    }
  }
  (void )WlzFreeAffineTransform(initTr);
  return(errNum);
}

/*!
* \return				void
* \ingroup	WlzTransform
* \brief	Inserts list1 at the head of list0.
* \param	list0			Pointer to the list of shells.
* \param 	list1			List of shells
*/
static void	WlzMatchICPShellListInsertList(WlzMatchICPShellList *list0,
					       WlzMatchICPShellList *list1)
{
  if(list0 && list1)
  {
    if(list0->head == NULL)
    {
      list0->head = list1->head;
      list0->tail = list1->tail;
    }
    else if(list1->tail)
    {
      list0->head->prev = list1->tail;
      list1->tail->next = list0->head;
      list0->head = list1->head;
    }
    list1->head = list1->tail = NULL;
  }
}

/*!
* \return				Affine transform found, NULL
*					on error.
* \ingroup	WlzTransform
* \brief	Establishes an affine transfrom matching the source
*		model to the target model. The target model's vertices
*		are in an already populated kD-tree in the workspace.
* \param	tTree			Given kD-tree populated by the
*					target vertices such that the
*					nodes of the tree have the same
*					indicies as the given target vertices
*					and normals.
* \param	trType			The required type of transform,
*					must be either WLZ_TRANSFORM_2D_REG,
*					or WLZ_TRANSFORM_2D_AFFINE.
* \param	vType			Type of vertices.
* \param	sgnNrm			Non zero if normal component signs are
*					consistant.
* \param	nTV			Number of target vertices.
* \param	tVx			The target vertices.
* \param	tNr			The target normals.
* \param        nSV			Number of source vertices.
* \param        sVx 			The source vertices.
* \param	sNr			The source normals.
* \param	iBuf			Buffer with at least nSV ints.
* \param	tVBuf			A buffer with room for at least
*					nS vertices.
* \param	sVBuf			A buffer with room for at least
*					nS vertices.
* \param	wBuf			A buffer with room for at least
*					nS doubles.
* \param	maxItr			Maximum number of iterations.
* \param	initTr			Initial affine transform, may be NULL.
* \param	dstConv			Destination pointer for a
*					convergence flag which is set to
*					a non zero value if the registration
*					converges.
* \param	usrWgtFn		User supplied weighting function.
* \param	usrWgtData		User supplied weighting data.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzAffineTransform *WlzMatchICPRegModel(AlcKDTTree *tTree,
				WlzTransformType trType,
				WlzVertexType vType, int sgnNrm,
				int nTV, WlzVertexP tVx, WlzVertexP tNr,
				int nSV, WlzVertexP sVx, WlzVertexP sNr,
				int *iBuf, WlzVertexP tVBuf, WlzVertexP sVBuf,
				double *wBuf, int maxItr, 
				WlzAffineTransform *initTr, int *dstConv,
				WlzRegICPUsrWgtFn usrWgtFn, void *usrWgtData,
				WlzErrorNum *dstErr)
{
  int		idV,
  		conv = 0;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  for(idV = 0; idV < nSV; ++idV)
  {
    *(iBuf + idV) = idV;
  }
  tr = WlzRegICPTreeAndVertices(tTree, trType, vType, sgnNrm,
  			     nTV, tVx, tNr, nSV, iBuf, sVx, sNr,
			     tVBuf, sVBuf, wBuf, maxItr, initTr,
			     &conv, usrWgtFn, usrWgtData, &errNum);
  if(dstConv)
  {
    *dstConv = conv;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \return				Affine transform found, NULL
*					on error.
* \ingroup	WlzTransform
* \brief	Establishes an affine transfrom matching the source
*		shell to the target model. The target model's vertices
*		are in an already populated kD-tree.
* \param	tTree			Given kD-tree populated by the
*					target vertices such that the
*					nodes of the tree have the same
*					indicies as the given target vertices
*					and normals.
* \param	tGM			Target model.
* \param	sS			Source shell.
* \param	globTr			Global affine transform.
* \param	trType			The required type of transform,
*					must be either WLZ_TRANSFORM_2D_REG,
*					or WLZ_TRANSFORM_2D_AFFINE.
* \param	vType			Type of vertices.
* \param	sgnNrm			Non zero if the normals are consistant.
* \param	nTV			Number of target vertices.
* \param	tVx			The target vertices.
* \param	tNr			The target normals.
* \param        nSV			Number of source vertices.
* \param        sVx 			The source vertices.
* \param	sNr			The source normals.
* \param	iBuf			Buffer with at least nSV ints.
* \param	tVBuf			A buffer with room for at least
*					nS vertices.
* \param	sVBuf			A buffer with room for at least
*					nS vertices.
* \param	wBuf			A buffer with room for at least
*					nS doubles.
* \param	maxItr			Maximum number of iterations.
* \param	maxDisp			Maximum displacement.
* \param	maxAng			maximum angle.
* \param	maxDeform		Maximum deformation.
* \param	initTr			Initial affine transform, may be NULL.
* \param	dstConv			Destination pointer for a
*					convergence flag which is set to
*					a non zero value if the registration
*					converges.
* \param	usrWgtFn		User supplied weighting function.
* \param	usrWgtData		User supplied weighting data.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzAffineTransform *WlzMatchICPRegShell(AlcKDTTree *tTree,
				WlzGMModel *tGM, WlzGMShell *sS,
				WlzAffineTransform *globTr,
				WlzTransformType trType,
				WlzVertexType vType, int sgnNrm,
				int nTV, WlzVertexP tVx, WlzVertexP tNr,
				int nSV, WlzVertexP sVx, WlzVertexP sNr,
				int *iBuf,
				WlzVertexP tVBuf, WlzVertexP sVBuf,
				double *wBuf, int maxItr, 
				double maxDisp, double maxAng, 
				double maxDeform,
				WlzAffineTransform *initTr, int *dstConv,
				WlzRegICPUsrWgtFn usrWgtFn, void *usrWgtData,
				WlzErrorNum *dstErr)
{
  int		nSSV,
  		conv = 0;
  WlzGMVertex	*v0;
  WlzGMVertexT	*vT0;
  WlzGMDiskT	*dT0;
  WlzGMEdgeT	*cET,
  		*fET;
  WlzGMLoopT	*cLT,
  		*fLT;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find the indicies of the source shell's vertices. */
  nSSV = 0;
  cLT = fLT = sS->child;
  do
  {
    cET = fET = cLT->edgeT;
    do
    {
      vT0 = cET->vertexT;
      dT0 = vT0->diskT;
      v0 = dT0->vertex;
      if((v0->diskT == dT0) && (dT0->vertexT == vT0))
      {
	*(iBuf + nSSV++) = v0->idx;
      }
      cET = cET->next;
    } while(cET != fET);
    cLT = cLT->next;
  } while(cLT != fLT);
  /* Register the source shell's vertices with the target model using the
   * existing kD-tree. */
  tr = WlzRegICPTreeAndVertices(tTree, trType, vType, sgnNrm,
  			  nTV, tVx, tNr, nSSV, iBuf, sVx, sNr,
			  tVBuf, sVBuf, wBuf,
			  maxItr, initTr, &conv,
			  usrWgtFn, usrWgtData,
			  &errNum);
  if((errNum == WLZ_ERR_NONE) && conv)
  {
    conv = WlzMatchICPRegCheckConv(sS, globTr, tr, maxDisp, maxAng, maxDeform);
  }
  else if((errNum == WLZ_ERR_ALG_CONVERGENCE) && (conv == 0))
  {
    /* Failure to register a shell is not an error if the shell is
     * inappropriate, e.g. normals are opposite. */
    WlzFreeAffineTransform(tr);
    tr = NULL;
    errNum = WLZ_ERR_NONE;
  }
  if(dstConv)
  {
    *dstConv = conv;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tr);
}

/*!
* \return	Non-zero if transform is within limits.
* \ingroup	WlzTransform
* \brief	Check that the transform is within the given limits of
* 		displacement and change in area to the shells axis
*		aligned bounding box.
* \param	sS			Given shell.
* \param	globTr			Global affine transform.
* \param	tr			Given affine transform.
* \param	maxDisp			Maximum displacement.
* \param	maxAng			Maximum angle.
* \param	maxDeform		Maximum deformation.
*/
static int	WlzMatchICPRegCheckConv(WlzGMShell *shell,
				 	WlzAffineTransform *globTr,
				 	WlzAffineTransform *tr,
					double maxDisp, double maxAng,
					double maxDeform)
{
  int 		ok = 0;

  switch(WlzAffineTransformDimension(tr, NULL))
  {
    case 2:
      ok = WlzMatchICPRegCheckConv2D(shell, globTr, tr, maxDisp, maxAng,
      				     maxDeform);
      break;
    case 3:
      ok = WlzMatchICPRegCheckConv3D(shell, globTr, tr, maxDisp, maxAng,
      				     maxDeform);
      break;
  }
  return(ok);
}

/*!
* \return	Non-zero if transform is within limits.
* \ingroup	WlzTransform
* \brief	Check that the transform is within the given limits of
* 		displacement and change in area to the shells axis
*		aligned 2D bounding box.
* \param	sS			Given shell.
* \param	globTr			Global affine transform.
* \param	tr			Given affine transform.
* \param	maxDisp			Maximum displacement.
* \param	maxAng			Maximum angle.
* \param	maxDeform		Minimum area ratio.
*/
static int	WlzMatchICPRegCheckConv2D(WlzGMShell *shell,
					  WlzAffineTransform *globTr,
					  WlzAffineTransform *tr,
					  double maxDisp, double maxAng,
					  double maxDeform)
{
  int		ok = 1;
  double	tD0,
  		tD1,
		param;
  WlzDVertex2	gV,
		pV,
  		s0V,
		s1V,
  		tV;
  WlzDBox2	sBB,
		gBB,
  		tBB;

  (void )WlzGMShellGetGBB2D(shell, &sBB);
  s0V.vtX = sBB.xMin;
  s0V.vtY = sBB.yMin;
  /* Compute the displacement and compare against the given maximum. */
  gV = WlzAffineTransformVertexD2(globTr, s0V, NULL);
  tV = WlzAffineTransformVertexD2(tr, s0V, NULL);
  WLZ_VTX_2_SUB(tV, tV, gV);
  param = WLZ_VTX_2_LENGTH(tV);
  if(param > maxDisp)
  {
    ok = 0;
  }
  if(ok)
  {
    /* Compute the ratio change in area of the axis aligned bounding box and
     * compare against the given maximum. */
    gBB = WlzAffineTransformBBoxD2(globTr, sBB, NULL);
    tBB = WlzAffineTransformBBoxD2(tr, sBB, NULL);
    tD0 = (tBB.xMax - tBB.xMin + 1) * (tBB.yMax - tBB.yMin + 1);
    tD1 = (gBB.xMax - gBB.xMin + 1) * (gBB.yMax - gBB.yMin + 1);
    param = 1.0 - ((tD0 >= tD1)? tD1 / tD0: tD0 / tD1);
    /* The value of param is
     *   small if areas of bounding boxes are similar
     *   large (but <= 1.0) if the areas are very different
     */
    if(param > maxDeform)
    {
      ok = 0;
    }
  }
  if(ok)
  {
    /* Compute the angle using two vectors defined by the minima and maxima of
     * the shell geometry, transform the vectors and compare against the given
     * maximum. */
    s0V.vtX = sBB.xMin;
    s0V.vtY = sBB.yMin;
    s1V.vtX = sBB.xMax;
    s1V.vtY = sBB.yMax;
    pV = WlzAffineTransformVertexD2(globTr, s0V, NULL);
    gV = WlzAffineTransformVertexD2(globTr, s1V, NULL);
    WLZ_VTX_2_SUB(gV, gV, pV);
    pV = WlzAffineTransformVertexD2(tr, s0V, NULL);
    tV = WlzAffineTransformVertexD2(tr, s1V, NULL);
    WLZ_VTX_2_SUB(tV, tV, pV);
    /* Find angle from the scalar product of the two vectors. */
    param = WLZ_VTX_2_DOT(tV, gV);
    tD0 = WLZ_VTX_2_SQRLEN(gV);
    tD0 *= WLZ_VTX_2_SQRLEN(tV);
    if(tD0 < DBL_EPSILON)
    {
      ok = 0;
    }
    else
    {
      param /= sqrt(tD0);
      if(param < cos(maxAng))
      {
	ok = 0;
      }
    }
  }
  return(ok);
}

/*!
* \return	Non-zero if transform is within limits.
* \ingroup	WlzTransform
* \brief	Check that the transform is within the given limits of
* 		displacement and change in area to the shells axis
*		aligned 3D bounding box.
* \param	sS			Given shell.
* \param	globTr			Global affine transform.
* \param	tr			Given affine transform.
* \param	maxDisp			Maximum displacement.
* \param	maxAng			Maximum angle.
* \param	maxDeform		Minimum area ratio.
*/
static int	WlzMatchICPRegCheckConv3D(WlzGMShell *shell,
					  WlzAffineTransform *globTr,
					  WlzAffineTransform *tr,
					  double maxDisp, double maxAng,
					  double maxDeform)
{
  int		ok = 0;

  /* TODO implement this function. */
  return(ok);
}

/*!
* \return				Number of matched points.
* \ingroup	WlzTransform
* \brief	Computes matching points from the broken source shells
*		and their affine transforms.
* \param	tTree			kD-tree populated by the target
*					vertices.
* \param	mS			Match shell.
* \param	vType			Type of vertices.
* \param	tVx			Target vertex geometries.
* \param 	tNr			Target normals.
* \param 	sNr			Source normals.
* \param	tMatch			Allocated space for the target
*					match vertices.
* \param	sMatch			Allocated space for the source
*					match vertices.
* \param	offset			Offset into tMatch and smatch.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static int	WlzMatchICPGetPoints(AlcKDTTree *tTree, WlzMatchICPShell *mS,
				     WlzVertexType vType, WlzVertexP tVx,
				     WlzVertexP tNr, WlzVertexP sNr,
				     WlzVertexP tMatch, WlzVertexP sMatch,
				     int offset, WlzErrorNum *dstErr)
{
  int		nMatch = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(vType)
  {
    case WLZ_VERTEX_D2:
      nMatch = WlzMatchICPGetPoints2D(tTree, mS, tVx.d2, tNr.d2, sNr.d2,
      				      tMatch.d2 + offset, sMatch.d2 + offset,
				      &errNum);
      break;
    case WLZ_VERTEX_D3: /* FALLTHROUGH */
    default:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nMatch);
}

/*!
* \return				Number of matched points.
* \ingroup	WlzTransform
* \brief	Computes matching points from the broken source shell
*		it's associated affine transform and the target kD-tree.
* 		Finds the 2D vertex at the given shell's midpoint, applies
*		the affine transform associated with the shell in the match
*		shell data structure to the vertices geometry then finds the
*		closest point in the target model to the transformed vertex
*		using the target kD-tree.
* \param	tTree			kD-tree populated by the target
*					vertices.
* \param	mS			Match shell.
* \param	tVx			Target vertex geometries.
* \param 	tNr			Target normals.
* \param 	sNr			Source normals.
* \param	tMatch			Allocated space for the target
*					match vertex.
* \param	sMatch			Allocated space for the source
*					match vertex.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static int	WlzMatchICPGetPoints2D(AlcKDTTree *tTree, WlzMatchICPShell *mS,
				       WlzDVertex2 *tVx, 
				       WlzDVertex2 *tNr, WlzDVertex2 *sNr,
				       WlzDVertex2 *tMatch,
				       WlzDVertex2 *sMatch,
				       WlzErrorNum *dstErr)
{
  int		nMatch = 0;
  WlzGMShell	*sS;
  WlzGMVertex	*cV;
  WlzDVertex2	sV,
		tV,
  		tSV;
  double	cDist;
  double	vxD2[2];
  WlzAffineTransform *iTr = NULL;
  AlcKDTNode	*tNode;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minSpx = 5; /* Avoid end regions of shells. */

  if(mS && ((sS = mS->shell) != NULL) && (sS->child == sS->child->next))
  {
    cV = WlzMatchICPLoopTMaxCurv2D(sS->child, sNr, minSpx, NULL);
    /* cV = WlzMatchICPLoopTMid(sS->child); */
    if(cV)
    {
      (void )WlzGMVertexGetG2D(cV, &sV);
      tSV = WlzAffineTransformVertexD2(mS->tr, sV, NULL);
      vxD2[0] = tSV.vtX;
      vxD2[1] = tSV.vtY;
      tNode = AlcKDTGetNN(tTree, vxD2, DBL_MAX, &cDist, NULL);
      if(tNode)
      {
	tV = *(tVx + tNode->idx);
	*(sMatch) = sV;
	*(tMatch) = tV;
	nMatch = 1;
      }
    }
    if(dstErr)
    {
      *dstErr = errNum;
    }
  }
  return(nMatch);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Filter the matched points by ranking the matches by their
*		plausibility and removing those whichtre implausible.
* \param	vType		Vertex type.
* \param	tVx		Given target matches.
* \param	sVx		Given source matches.
* \param	nVx		Pointer to the number of matched point pairs,
				may be modified on return.
* \param	globTr		The global affine transform associated with
*				the source, which maps it onto the target.
* \param	nN		Number of points in neighbourhood, must
*				be \f$> 2\f$.
* \param	impThr		Implausibility threshold which should be
* 				greater than zero.
*/
static WlzErrorNum WlzMatchICPFilterPts(WlzVertexType vType,
					WlzVertexP tVx, WlzVertexP sVx,
					int *nVx,
					WlzAffineTransform *globTr,
					int nN,
					double impThr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(vType)
  {
    case WLZ_VERTEX_D2:
      errNum = WlzMatchICPFilterPtsRmDup2D(tVx.d2, sVx.d2, nVx);
      if(*nVx > 0)
      {
	errNum = WlzMatchICPFilterPtsDisp2D(tVx.d2, sVx.d2, nVx, globTr, nN,
	    			            impThr);
      }
      break;
    case WLZ_VERTEX_D3: /* FALLTHROUGH */
    default:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Filter the matched points to avoid duplicate target
*		points.
*		The match points are sorted by target point position,
*		then scanned through skiping any duplicate target points.
* \param	tVx		Given target matches.
* \param	sVx		Given source matches.
* \param	nVxP		Pointer to the number of matches for both
*				input and return.
*/
static WlzErrorNum WlzMatchICPFilterPtsRmDup2D(WlzDVertex2 *tVx,
					       WlzDVertex2 *sVx, int *nVxP)
{
  int		idN,
  		idM,
		nVx;
  WlzDVertex2	dVx;
  WlzMatchICPTPPair2D *buf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	epsilon = 0.0001;

  nVx = *nVxP;
  if((buf = (WlzMatchICPTPPair2D *)
            AlcMalloc(sizeof(WlzMatchICPTPPair2D) * nVx)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy the match points to the new buffer. */
    for(idN = 0; idN < nVx; ++idN)
    {
      (buf + idN)->sVx = *(sVx + idN);
      (buf + idN)->tVx = *(tVx + idN);
    }
    /* Sort the buffer match points, first by vtY and then by vtX, into
     * increasing order. */
    qsort(buf, nVx, sizeof(WlzMatchICPTPPair2D), 
    	  (int (*)(const void *, const void *))WlzMatchICPTPPairSortFnD);
    /* Now put the buffer match points back skiping over those with
     * near identical target positions. */
    idN = 0;
    idM = 1;
    *sVx = buf->sVx;
    *tVx = buf->tVx;
    while(idM < nVx)
    {
      WLZ_VTX_2_SUB(dVx, *(tVx + idN), (buf + idM)->tVx);
      if((fabs(dVx.vtX) > epsilon) || (fabs(dVx.vtY) > epsilon))
      {
	++idN;
	*(sVx + idN) = (buf + idM)->sVx;
	*(tVx + idN) = (buf + idM)->tVx;
      }
      ++idM;
    }
    nVx = idN;
    *nVxP = nVx;
  }
  AlcFree(buf);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Filter the matched points by ranking the matches by their
*		implausibly and removing those for which the implausibly
*		is above threshold.
*		Affine transformed source points are computed and used to
*		compute the implausibly of the matches. The implausibly
*		of a match is computed by first finding all target points
*		that are near neighbours of the particular target point.
*		The mean and standard deviation of displacements in the
*		neighbourhood of the target point are then computed.
*		The implausibly of a match is computed in 'displacement
*		space'. Given a mean \f$\mu\f$ and a standard deviation
*		\f$\sigma\f$, then \f$\mu\f$ and \f$\sigma\f$ define an
*		ellipse with centre \f$\mu\f$ and half major axes
*		\f$\sigma\f$.
*		The implausibly of a match with displacement \f$d\f$
*		is: \f$\frac{{|d - \mu|}^2}{{|d_e - \mu|}^2}\f$ where
*		\f$d_e\f$ is the distance from \f$\mu\f$ through \f$d\f$
*		to the ellipse.
*		This gives an implausibly value in the range 0.0-\f$\infty\f$,
*		with 0.0 being completely plausible and \f$\infty\f$ being
*		completely implausible. The matches are ranked according
*		to their implausibly with those outside the ellipsoid
*		being removed by reducing the value of *nVx.
* \param	tVx		Given target matches.
* \param	sVx		Given source matches.
* \param	nVxP		Pointer to the number of matches for both
*				input and return.
* \param	globTr		The global affine transform associated with
*				the source, which maps it onto the target.
* \param	nN		Number of points in neighbourhood, must
*				be \f$> 2\f$.
* \param	impThr		Implausibility threshold which should be
*				greater than zero.
*/
static WlzErrorNum WlzMatchICPFilterPtsDisp2D(WlzDVertex2 *tVx,
				              WlzDVertex2 *sVx, int *nVxP,
				              WlzAffineTransform *globTr,
					      int nN, double impThr)
{
  int		tI0,
  		tI1,
  		idN,
  		idS,
  		idT,
		nVx;
  double	tD0,
  		tD1;
  WlzDVertex2	tV0,
		sum,
		sumSq,
  		mean,
		stddev;
  int		*indicies = NULL;
  double	*bufD = NULL,
        	*imp = NULL;
  WlzDVertex2	*bufU = NULL,
  		*bufV = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nVx = *nVxP;
  if(nN > nVx)
  {
    nVx = nN;
  }
  /* Create buffers for the transformed source points/point buffers,
   * implausibly measure, neighbour indicies/rank indicies and neighbour
   * distances. */
  if(((imp = (double *)AlcMalloc(sizeof(double) * nVx)) == NULL) ||
     ((bufU = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * nVx)) == NULL) ||
     ((bufV = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) * nVx)) == NULL) ||
     ((indicies = (int *)AlcMalloc(sizeof(int) * nVx)) == NULL) ||
     ((bufD = (double *)AlcMalloc(sizeof(double) * nVx)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create a set of transformed source points. */
    for(idS = 0; idS < nVx; ++idS)
    {
      *(bufV + idS) = WlzAffineTransformVertexD2(globTr, *(sVx + idS), NULL);
    }
    /* Compute implausibilities for each target vertex indexed by idT. */
    for(idT = 0; idT < nVx; ++idT)
    {
      /* Get indicies of neighbourhood. */
      for(idN = 0; idN < nN; ++idN)
      {
        *(bufD + idN) = DBL_MAX;
      }
      for(idS = 0; idS < nVx; ++idS)
      {
        /* Compute distance to neighbour and insert into the neighbourhood
	 * if it's less than the current maximum distance in it. */
	WLZ_VTX_2_SUB(tV0, *(tVx + idT), *(tVx + idS));
        tD0 = WLZ_VTX_2_LENGTH(tV0);
	if(tD0 < *(bufD + nN - 1))
	{
	  /* Add to point to neighbourhood. */
	  idN = 0;
	  while((*(bufD + idN) < tD0) && (idN < nN))
	  {
	    ++idN;
	  }
	  tI0 = idS;
	  while(idN < nN)
	  {
	    tI1 = *(indicies + idN);
	    tD1 = *(bufD + idN);
	    *(indicies + idN) = tI0;
	    *(bufD + idN) = tD0;
	    tI0 = tI1;
	    tD0 = tD1;
	    ++idN;
	  }
	}
      }
      /* Compute mean and standard deviation of displacements in the
       * neighbourhood while not including the vertex itself. */
      sum.vtX = sum.vtY = sumSq.vtX = sumSq.vtY = 0.0;
      for(idN = 1; idN < nN; ++idN)
      {
	tI0 = *(indicies + idN);
        WLZ_VTX_2_SUB(tV0, *(tVx + tI0), *(bufV + tI0));
	sum.vtX += tV0.vtX;
	sum.vtY += tV0.vtY;
	sumSq.vtX += tV0.vtX * tV0.vtX;
	sumSq.vtY += tV0.vtY * tV0.vtY;
      }
      mean.vtX = sum.vtX / (nN - 1);
      mean.vtY = sum.vtY / (nN - 1);
      stddev.vtX = sqrt((sumSq.vtX - (sum.vtX * sum.vtX / (nN - 1))) /
                        (nN - 2));
      stddev.vtY = sqrt((sumSq.vtY - (sum.vtY * sum.vtY / (nN - 1))) /
                        (nN - 2));
      /* Compute implausibilities. */
      WLZ_VTX_2_SUB(tV0, *(tVx + idT), *(bufV + idT));
      *(imp + idT) = WlzGeomEllipseVxDistSq(mean, stddev, tV0);
    }
    /* Rank matches by plausibility, least plausible last. */
    for(idT = 0; idT < nVx; ++idT)
    {
      *(indicies + idT) = idT;
    }
    (void )AlgHeapSortIdx(imp, indicies, nVx, WlzMatchICPDblSortFnD);
    for(idT = 0; idT < nVx; ++idT)
    {
      *(bufU + idT) = *(sVx + idT);
      *(bufV + idT) = *(tVx + idT);
    }
    for(idT = 0; idT < nVx; ++idT)
    {
      *(sVx + idT) = *(bufU + *(indicies + idT));
      *(tVx + idT) = *(bufV + *(indicies + idT));
      *(bufD + idT) = *(imp + *(indicies + idT));
    }
    /* Reduce the number of matches to filter out any that are implausible. */
    while((nVx > 0) && (*(bufD + nVx - 1) > impThr))
    {
      --nVx;
    }
  }
  *nVxP = nVx;
  AlcFree(bufD);
  AlcFree(bufU);
  AlcFree(bufV);
  AlcFree(imp);
  AlcFree(indicies);
  return(errNum);
}

/*!
* \return				Error code.
* \ingroup	WlzTransform
* \brief	Find all elements within the given shell which have
*		high connectivity and remove them, while at the same
*		time building a list of the shells that are derived from
*		the given shell.
* \param	cbData			Callback data structure in which
*					the list of new shells is being built.
* \param	sGM			Source model.
* \param	sS			Source shell.
* \param	iBuf			Buffer containing the indicies of
*					vertices to be removed from the
*					model.
*/
static WlzErrorNum	WlzMatchICPBreakShellCon(WlzMatchICPCbData *cbData,
				WlzGMModel *sGM, WlzGMShell *sS, int *iBuf)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(sGM->type)
  {
    case WLZ_GMMOD_2I: /* FALLTHROUGH */
    case WLZ_GMMOD_2D: /* FALLTHROUGH */
    case WLZ_GMMOD_2N:
      errNum = WlzMatchICPBreakShellCon2D(cbData, sGM, sS, iBuf);
      break;
    case WLZ_GMMOD_3I: /* FALLTHROUGH */
    case WLZ_GMMOD_3D: /* FALLTHROUGH */
    case WLZ_GMMOD_3N:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return				Error code.
* \ingroup	WlzTransform
* \brief	Find all vertices within the given shell which have
*		high connectivity and remove them, while at the same
*		time building a list of the shells that are derived from
*		the given shell.
* \param	cbData			Callback data structure in which
*					the list of new shells is being built.
* \param	sGM			Source model.
* \param	sS			Source shell.
* \param	iBuf			Buffer containing the indicies of
*					vertices to be removed from the
*					model.
*/
static WlzErrorNum WlzMatchICPBreakShellCon2D(WlzMatchICPCbData *cbData,
				WlzGMModel *sGM, WlzGMShell *sS, int *iBuf)
{
  int		nHCV;
  WlzGMVertex	*v0;
  WlzGMVertexT	*vT0;
  WlzGMDiskT	*dT0;
  WlzGMEdgeT	*cET,
  		*fET;
  WlzGMLoopT	*cLT,
  		*fLT;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find all vertices with high connectivity that are in the given
   * shell. */
  nHCV = 0;
  cLT = fLT = sS->child;
  do
  {
    cET = fET = cLT->edgeT;
    do
    {
      vT0 = cET->vertexT;
      /* Test for high connectivity */
      if((vT0 != vT0->next) && (vT0 != vT0->next->next))
      {
	dT0 = vT0->diskT;
	v0 = dT0->vertex;
	/* Only delete a vertex once. */
	if((v0->diskT == dT0) && (dT0->vertexT == vT0))
	{
	  *(iBuf + nHCV++) = v0->idx;
	}
      }
      cET = cET->next;
    } while(cET != fET);
    cLT = cLT->next;
  } while(cLT != fLT);
  /* Remove all the high connectivity that were found. */
  if(nHCV > 0)
  {
    errNum = WlzMatchICPRemoveVertices(cbData, sGM, iBuf, nHCV);
  }
  return(errNum);
}

/*!
* \return				Error code.
* \ingroup	WlzTransform
* \brief	Break the given shell by removing vertices near the
*		midpoint of the shell. The given shell is known to only
*		have simple connectivity, eg a shell with simple connectivity
*		in a 2D model will no more than two edges using a single
*		vertex.
* \param	cbData			Callback data structure in which
*					the list of new shells is being built.
* \param	sGM			Source model.
* \param	sMS			The match shell.
* \param	sNr			The source shell normals.
* \param	iBuf			Buffer containing the indicies of
*					vertices to be removed from the
*					model.
* \param	minSpx			Minimum number of simplicies to
*					leave in a shell.
*/
static WlzErrorNum	WlzMatchICPBreakShellMid(WlzMatchICPCbData *cbData,
				WlzGMModel *sGM, WlzMatchICPShell *sMS,
				WlzVertexP sNr, int *iBuf, int minSpx)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(sGM->type)
  {
    case WLZ_GMMOD_2I: /* FALLTHROUGH */
    case WLZ_GMMOD_2D: /* FALLTHROUGH */
    case WLZ_GMMOD_2N:
      /* If the shell has more than one loop topology element then it is
       * a closed loop so this needs to be broken into a simple non-cyclic
       * chain of edge segments before it can be broken by removing low
       * curvature vertices. */
      errNum = WlzMatchICPBreakShellCur2D(cbData, sGM, sMS, sNr.d2,
      					  iBuf, minSpx);
      break;
    case WLZ_GMMOD_3I: /* FALLTHROUGH */
    case WLZ_GMMOD_3D: /* FALLTHROUGH */
    case WLZ_GMMOD_3N:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return				Error code.
* \ingroup	WlzTransform
* \brief	Make sure that all simplex topology elements are connected
*		through the next and prev pointers without recourse to the
*		opp pointer for connectivity.
* \param	sGM			Source model.
* \param	sS			Given shell.
* \param	sNr			Source normals.
* \param	minSpx			Minimum number of simplicies to
*					leave in a shell.
*/
static WlzErrorNum	WlzMatchICPBreakShellLoopTs(WlzGMModel *sGM,
					WlzGMShell *sS, WlzVertexP sNr,
					int minSpx)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(sGM->type)
  {
    case WLZ_GMMOD_2I: /* FALLTHROUGH */
    case WLZ_GMMOD_2D: /* FALLTHROUGH */
    case WLZ_GMMOD_2N:
      errNum = WlzMatchICPBreakShellLop2D(sGM, sS, sNr.d2, minSpx);
      break;
    case WLZ_GMMOD_3I: /* FALLTHROUGH */
    case WLZ_GMMOD_3D: /* FALLTHROUGH */
    case WLZ_GMMOD_3N:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  return(errNum);

}

/*!
* \return				Woolz error code.
* \ingroup	WlzTransform
* \brief	Make sure that all simplex topology elements are connected
*		through the next and prev pointers without recourse to the
*		opp pointer for connectivity. In 2D this means that a shell
*		with no high connectivity vertices has two opposite loop
*		topology elements. This is fixed by removing the a vertex
*		that has the lowest curvature and is distant (along the
*		loopT) from the vertex with the highest curvature.
* \param	sGM			Source model.
* \param	sS			Given shell.
* \param	sNr			Source normals.
* \param	minSpx			Minimum number of simplicies to
*					leave in a shell.
*/
static WlzErrorNum	WlzMatchICPBreakShellLop2D(WlzGMModel *sGM,
					WlzGMShell *sS, WlzDVertex2 *sNr,
					int minSpx)
{
  WlzGMVertex	*sV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find the vertex which is at the point of the lowest curvature in a loop
   * topology element of the shell. */
  sV = WlzMatchICPLoopTMinCurv2D(sS->child, sNr, minSpx, NULL);
  errNum = WlzGMModelDeleteV(sGM, sV);
  return(errNum);
}

/*!
* \return				Error code.
* \ingroup	WlzTransform
* \brief	Breaks the shells of the given source model near to their
*		midpoint but avoiding regions of high curvature. This is
*		done by removing vertices at these points from the model.
* \param	cbData			Callback data structure in which
*					the list of new shells is being built.
* \param	sGM			Source model.
* \param	sMS			The match shell with pointers to the
*					shell and it's number of simplicies.
* \param	sNr			Source normals.
* \param	iBuf			Buffer containing the indicies of
*					vertices to be removed from the
*					model.
* \param	minSpx			Minimum number of simplicies to
*					leave in a shell.
*/
static WlzErrorNum WlzMatchICPBreakShellCur2D(WlzMatchICPCbData *cbData,
				WlzGMModel *sGM, WlzMatchICPShell *sMS,
				WlzDVertex2 *sNr, int *iBuf, int minSpx)
{
  int		idD;
  WlzGMVertex	*sV;
  WlzGMLoopT	*sLT;
  WlzGMShell	*sS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idD = 0;
  sS = sMS->shell;
  sLT = sS->child;
  do
  {
    if((sMS->size > (minSpx * 3)) &&
       ((sV = WlzMatchICPLoopTMinCurv2D(sLT, sNr, minSpx, NULL)) != NULL))
    {
      *(iBuf + idD++) = sV->idx;
    }
    sLT = sLT->next;
  } while(sLT != sS->child);
  if(idD > 0)
  {
    errNum = WlzMatchICPRemoveVertices(cbData, sGM, iBuf, idD);
  }
  return(errNum);
}

/*!
* \return
* \ingroup	WlzTransform
* \brief	Find the vertex at the middle of the given loop topology
* 		element, where the given loop topology element is the only
*		one in it's parents shell.
* \param	gLT			Given loop topology element.
*/
static WlzGMVertex *WlzMatchICPLoopTMid(WlzGMLoopT *gLT)
{
  int		cnt;
  WlzGMVertex	*mV = NULL;
  WlzGMEdgeT	*fET,
  		*lET,
		*mET,
		*tET;

  /* Find edge topology element at one end. */
  fET = gLT->edgeT;
  while(((tET = fET->prev) != fET->opp) && (tET != gLT->edgeT))
  {
    fET = tET;
  }
  /* Find edge topology element at other end, counting elements. */
  cnt = 0;
  lET = fET;
  while(((tET = lET->next) != lET->opp) && (tET != fET))
  {
    lET = tET;
    ++cnt;
  }
  /* Go back half the edget topology elements. */
  cnt /= 2;
  mET = fET;
  while(cnt-- > 0)
  {
    mET = mET->next;
  }
  mV = mET->vertexT->diskT->vertex;
  return(mV);
}

/*!
* \return       Vertex at point of maximum curvature which avoids the
*		neighbourhoods of the loopT ends vertices (if they exist).
* \ingroup	WlzTransform
* \brief	Given a loop topology element and an array of normals
*		for the vertices of the loop topology element. Finds
*		the vertex at the point of maximum  curvature.
* 		The actual curvature is computed using the given
*		normals, using the normals of a vertex's previous
*		and next vertices in the loopT: For vetrex \f$mathbf{v_i}\f$
*		having previous and next vertices in the loopT
*		\f$mathbf{v_{i-1}}\f$ and \f$mathbf{v_{i+1}}\f$ respectively.
*		The curvature at \f$mathbf{v_i}\f$ is given by
*		  \f[c = \frac{1}{2}{(1 + \cos\theta)}\f]
*		where
*		  \f[\cos\theta = \frac{\mathbf{n_{i-1}} \cdot
                                        \mathbf{n_{i+1}}}
                                       {|\mathbf{n_{i-1}}|
				        |\mathbf{n_{i+1}}|}\f]
*		Because the normals are all unit normals this simplifies
*		to:
*		  \f[\cos\theta = \mathbf{n_{i-1}} \cdot \mathbf{n_{i+1}}\f]
*		Because it is the vertex at the position of maximum
*		curvature that is wanted the search is for the minimum value
*		of \f$\cos\theta\f$.
* \param	gLT			Given loop topology element.
* \param	gNr			Given vertex normals.
* \param	minEdg			Minimum number of edges to
*					skip at the terminal ends of an
*					loop topology element that is
*					an open edge chain.
* \param	dstAngle		Destination pointer for maximum
*					or minimum angle found.
*/
static WlzGMVertex *WlzMatchICPLoopTMaxCurv2D(WlzGMLoopT *gLT,
					WlzDVertex2 *gNr, int minEdg,
                                        double *dstAngle)
{
  int		idN,
		len = 0,
		okFlg = 1;
  double	cosAng,
  		mCosAng;
  WlzGMEdgeT	*cET,
  		*fET,
		*lET,
  		*tET0,
		*mET,
		*nET,
		*pET;
  WlzGMVertex	*mV = NULL;
  WlzDVertex2	n0,
  		n1;


  /* Find the first and last edge topology elements between which to search
   * for the vertex at the centre of a maximum curvature region. */
  if(gLT != gLT->next)
  {
    /* This loopT is not an open loop. */
    fET = gLT->edgeT;
    /* Count the edge topology elements. */
    tET0 = fET;
    while(((nET = tET0->next) != tET0->opp) && (nET != fET))
    {
      ++len;
      tET0 = nET;
    }
    if(len <= (2 * minEdg))
    {
      okFlg = 0;
    }
  }
  else
  {
    /* This loopT is an open loop, find it's end point. */
    tET0 = gLT->edgeT;
    while(((pET = tET0->prev) != tET0->opp) && (pET != gLT->edgeT))
    {
      tET0 = pET;
    }
    /* Go forwards minEdg edge topology elements to find the first
     * edge topology element. */
    idN = minEdg;
    while(((nET = tET0->next) != tET0->opp) && (idN-- > 0))
    {
      tET0 = nET;
    }
    if(idN == 0)
    {
      okFlg = 0;
    }
    else
    {
      fET = tET0;
      /* Find the other end point. */
      while(((nET = tET0->next) != tET0->opp) && (nET != pET))
      {
	tET0 = nET;
      }
      if(nET == pET)
      {
	okFlg = 0;
      }
    }
  }
  if(okFlg)
  {
    /* Go back minEdg edge topology elements to find the last edge topology
     * element. */
    idN = minEdg;
    while(((pET = tET0->prev) != fET) && (idN-- > 0))
    {
      tET0 = pET;
    }
    if(idN < -1)
    {
      okFlg = 0;
    }
  }
  if(okFlg)
  {
    lET = tET0;
    len -= minEdg;
  }
  /* Having found the first and last edge topology elements between which to
   * search, now do the search for the vertex at the point of maximal
   * curvature. Using the cosine of the angle between the normals which
   * the modulus of which will be minimal at the position of highest
   * curvature.
   */
  if(okFlg)
  {
    mCosAng = 2.0;
    cET = fET->next;
    while(cET != lET)
    {
      n0 = *(gNr + cET->prev->vertexT->diskT->vertex->idx);
      n1 = *(gNr + cET->next->vertexT->diskT->vertex->idx);
      cosAng = fabs(WLZ_VTX_2_DOT(n0, n1));
      if(cosAng < mCosAng)
      {
	mET = cET;
        mCosAng = cosAng;
      }
      cET = cET->next;
    }
    mV = mET->vertexT->diskT->vertex;
    if(dstAngle)
    {
      if(mCosAng >= 1.0)
      {
        *dstAngle = 0.0;
      }
      else
      {
        *dstAngle = acos(mCosAng);
      }
    }
  }
  return(mV);
}

/*!
* \return       Vertex at point of minimum curvature which avoids the
*		neighbourhoods of the loopT ends vertices (if they exist)
*		and the vertex at the position of highest curvature.
* \ingroup	WlzTransform
* \brief	Finding the vertex at the position of minimum curvature
*		is not a well posed problem. This function attempts to
*		find a suitable vertex by considering the loop to have
*		four regions, each with an equal number of vertices.
*		First  a search is made for the vertex at the position of
*		highest curvature, if this is found in region \f$i\f$
*		then the lowest curvature vertex in region
*		\f$(i + 2) mod 4\f$ is chosen as the vertex at the position
*		of lowest curvature in the loopT.
*		See WlzMatchICPLoopTMaxCurv2D() for a description of how
*		curvature is established.
* \param	gLT			Given loop topology element.
* \param	gNr			Given vertex normals.
* \param	minEdg			Minimum number of edges to
*					skip at the terminal ends of an
*					loop topology element that is
*					an open edge chain.
* \param	dstAngle		Destination pointer for minimum
*					angle between normals.
*/
static WlzGMVertex *WlzMatchICPLoopTMinCurv2D(WlzGMLoopT *gLT,
					WlzDVertex2 *gNr, int minEdg,
                                        double *dstAngle)
{
  int		idN,
		mOff,
		len,
		regLn,
		mReg,
		okFlg = 1;
  double	cosAng,
  		mCosAng;
  WlzGMEdgeT	*cET,
  		*fET,
		*lET,
  		*tET0,
		*mET,
		*nET,
		*pET;
  WlzGMVertex	*mV = NULL;
  WlzDVertex2	n0,
  		n1;


  /* Find the first and last edge topology elements between which to search
   * for the vertex at the centre of a minimum curvature region. */
  if(gLT != gLT->next)
  {
    /* This loopT is not an open loop. */
    len = 0;
    fET = gLT->edgeT;
    /* Count the edge topology elements. */
    tET0 = fET;
    while(((nET = tET0->next) != tET0->opp) && (nET != fET))
    {
      ++len;
      tET0 = nET;
    }
    if(len <= (2 * minEdg))
    {
      okFlg = 0;
    }
  }
  else
  {
    /* This loopT is an open loop, find it's end point. */
    tET0 = gLT->edgeT;
    while(((pET = tET0->prev) != tET0->opp) && (pET != gLT->edgeT))
    {
      tET0 = pET;
    }
    /* Go forwards minEdg edge topology elements to find the first
     * edge topology element. */
    idN = minEdg;
    while(((nET = tET0->next) != tET0->opp) && (idN-- > 0))
    {
      tET0 = nET;
    }
    if(idN < -1)
    {
      okFlg = 0;
    }
    else
    {
      len = 0;
      fET = tET0;
      /* Find the other end point, counting the edge topology elements while
       * traveling. */
      while(((nET = tET0->next) != tET0->opp) && (nET != pET))
      {
	++len;
	tET0 = nET;
      }
      if(nET == pET)
      {
	okFlg = 0;
      }
    }
  }
  if(okFlg)
  {
    /* Go back minEdg edge topology elements to find the last edge topology
     * element. */
    idN = minEdg;
    while(((pET = tET0->prev) != fET) && (idN-- > 0))
    {
      tET0 = pET;
    }
    if(idN < -1)
    {
      okFlg = 0;
    }
  }
  if(okFlg)
  {
    lET = tET0;
    len -= minEdg;
  }
  /* Having found the first and last edge topology elements between which to
   * search, now do the search for the vertex at the point of maximal
   * curvature. Using the cosine of the angle between the normals which
   * the modulus of which will be minimal at the position of highest
   * curvature.
   */
  if(okFlg)
  {
    mCosAng = 2.0;
    mOff = idN = 0;
    cET = fET->next;
    while(cET != lET)
    {
      n0 = *(gNr + cET->prev->vertexT->diskT->vertex->idx);
      n1 = *(gNr + cET->next->vertexT->diskT->vertex->idx);
      cosAng = fabs(WLZ_VTX_2_DOT(n0, n1));
      if(cosAng < mCosAng)
      {
	mET = cET;
	mOff = idN;
        mCosAng = cosAng;
      }
      cET = cET->next;
      ++idN;
    }
    /* Compute which of the four regions the maximum curvature vertex is
     * in. */
    regLn = len / 4;
    mReg = ((mOff / regLn) + 2) % 4;
    if(mReg > 0)
    {
      /* Search forward toward the end because this is not the first region. */
      mOff = mReg * regLn;
      /* Find start of region. */
      idN = 0;
      cET = fET->next;
      while(idN < mOff)
      {
	cET = cET->next;
	++idN;
      }
      /* Find minimum curvature in this region. */
      idN  = 0;
      mCosAng = 0.0;
      while(idN < regLn)
      {
	n0 = *(gNr + cET->prev->vertexT->diskT->vertex->idx);
	n1 = *(gNr + cET->next->vertexT->diskT->vertex->idx);
	cosAng = fabs(WLZ_VTX_2_DOT(n0, n1));
	if(cosAng > mCosAng)
	{
	  mET = cET;
	  mCosAng = cosAng;
	}
	cET = cET->next;
	++idN;
      }
    }
    else
    {
      /* Find end of first region. */
      idN = 0;
      cET = fET->next;
      while(idN < regLn)
      {
	cET = cET->next;
	++idN;
      }
      /* Search back for minimum curvature from the end of the first region. */
      mCosAng = 0.0;
      while(idN > 0)
      {
	n0 = *(gNr + cET->prev->vertexT->diskT->vertex->idx);
	n1 = *(gNr + cET->next->vertexT->diskT->vertex->idx);
	cosAng = fabs(WLZ_VTX_2_DOT(n0, n1));
	if(cosAng > mCosAng)
	{
	  mET = cET;
	  mCosAng = cosAng;
	}
	cET = cET->prev;
	--idN;
      }
    }
    mV = mET->vertexT->diskT->vertex;
    if(dstAngle)
    {
      if(mCosAng >= 1.0)
      {
        *dstAngle = 0.0;
      }
      else
      {
        *dstAngle = acos(mCosAng);
      }
    }
  }
  return(mV);
}

/*!
* \return				Error code.
* \ingroup	WlzTransform
* \brief	Removes all vertices in the workspace's vertex flags
*		buffer from the source model.
* \param	cbData			Callback data structure in which
*					the list of new shells is being built.
* \param	sGM			Source model.
* \param	iBuf			Buffer containing the indicies of
*					vertices to be removed from the
*					model.
* \param	nD			Number of vertices to be removed.
*/
static WlzErrorNum WlzMatchICPRemoveVertices(WlzMatchICPCbData *cbData,
				WlzGMModel *sGM, int *iBuf, int nD)
{
  int		idD;
  AlcVector   	*vec;
  WlzGMVertex	*sV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Remove all these high connectivity vertices from the model creating
   * new shells. Make sure that the transform of each of these new shells
   * is a copy of the transform of the shell with the high connectivity
   * vertex. */
  idD = 0;
  vec = sGM->res.vertex.vec;
  while((errNum == WLZ_ERR_NONE) && (idD < nD))
  {
    sV = (WlzGMVertex *)AlcVectorItemGet(vec, *(iBuf + idD));
    if(sV && (sV->idx >= 0))
    {
      errNum = WlzGMModelDeleteV(sGM, sV);
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = cbData->errNum;
      }
    }
    ++idD;
  }
  return(errNum);
}

/*!
* \return	Weight value in the range [0.0-1.0].
* \ingroup      WlzTransform
* \brief	Compute an additional weighting for the matched pair of
* 		2D verticies.
*		Each pair of matched verticies has an additional weighting
*		computed which relates to the sensitivity of the match to
*		perturbations.
*		Weight is :
*		  \li 0.0 	If the vertex and perturbed vertex are in
*				different shells.
*		  \li 1.0       If the vertex and perturbed vertex are in
*				the same shell but not the same loopTs.
*		  \li 2.0       If the vertex and perturbed vertex share the
*				same loopT, but there is more than one loopT
*				per shell.
*		  \li [2.0-3.0] If the vertex and perturbed vertex share the
*				same loopT which is the only loopT of the
*				parent shell. In this case the weight is
*				computed by comparing Euclidean and Model
*				distances between the vertices.
*		The weight is finally normalized to the range [0-1.0].
* \param	curTr			Current affine transform.
* \param	tree			Given kD-tree populated by the
*					target vertices such that the nodes of
*					the tree have the sameindicies as the
*					given target vertices.
* \param	tVx			The target vertices.
* \param	sVx			The source vertices.
* \param	tMVx			The matched target vertex.
* \param	sMVx			The matched source vertex.
* \param	tGM			Target geometric model.
* \param	sGM			Source geometric model.
* \param	maxDisp			Maximum displacement for source
*					verticies.
* \param	nScatter		Number of verticies scattered about
*					the target vertex to test it's
*					sensitivity.
*/
static double	WlzMatchICPWeightMatches2D(WlzAffineTransform *curTr,
				       	AlcKDTTree *tree,
				       	WlzDVertex2 *tVx, WlzDVertex2 *sVx,
					WlzDVertex2 tMVx, WlzDVertex2 sMVx,
				       	WlzGMModel *tGM, WlzGMModel *sGM,
				       	double maxDisp, int nScatter)
{
  /* TODO CHECK THIS FUNCTION! */
  int		idN;
  double	distE,
   		distM,
		wgt;
  WlzDVertex2	diff,
   		disp,
		sMTVx,
   		tMVx0,
 		sMTVx0,
   		sMVx0;
  WlzGMLoopT	*tMLT,
  		*tMLT0;
  WlzGMShell	*tMS,
  		*tMS0;
  WlzGMVertex	*tMV,
  		*tMV0;
  AlcKDTNode	*tNode;
  double	vxD[2];
  const double	delta = 0.000001;

  if(nScatter > 0)
  {
    wgt = 0.0;
    tMV = WlzGMModelMatchVertexG2D(tGM, tMVx);
    sMTVx = WlzAffineTransformVertexD2(curTr, sMVx, NULL);
    tMLT = tMV->diskT->vertexT->parent->parent;
    tMS = tMLT->parent;
    for(idN = 0; idN < nScatter; ++idN)
    {
      /* Compute a new source vertex with a random displacement
      * (distance < maxDist) from the source vertex. */
      disp.vtX = ((AlgRandUniform() * 2.0) - 1.0) * ALG_M_SQRT1_2 * maxDisp;
      disp.vtY = ((AlgRandUniform() * 2.0) - 1.0) * ALG_M_SQRT1_2 * maxDisp;
      sMVx0.vtX = sMVx.vtX + disp.vtX;
      sMVx0.vtY = sMVx.vtY + disp.vtY;
      /* Transfrom the source vertex using the current affine transform. */
      sMTVx0 = WlzAffineTransformVertexD2(curTr, sMVx0, NULL);
      /* Find nearest neighbour of the transformed source vertex in the
      * target model. */
      vxD[0] = sMTVx0.vtX;
      vxD[1] = sMTVx0.vtY;
      if((tNode = AlcKDTGetNN(tree, vxD, DBL_MAX, NULL, NULL)) != NULL)
      {
	tMVx0 = *(tVx + tNode->idx);
	tMV0 = WlzGMModelMatchVertexG2D(tGM, tMVx0);
	if(tMV0)
	{
	  tMLT0 = tMV->diskT->vertexT->parent->parent;
	  tMS0 = tMLT0->parent;
	  if(tMS == tMS0)
	  {
	    /* Same shell. */
	    wgt += 1.0;
	    if(tMLT == tMLT0)
	    {
	      /* Same loopT */
	      wgt += 1.0;
	      if(tMLT == tMLT->next)
	      {
		if(tMV == tMV0)
		{
		  /* Same vertex! */
		  wgt += 1.0;
		}
		else
		{
		  /* Compute the Euclidean distance between the transformed
		   * source vertex and the transformed displaced source
		   * vertex. */
		  WLZ_VTX_2_SUB(diff, sMTVx, sMTVx0);
		  distE = WLZ_VTX_2_LENGTH(diff);
		  if(distE > delta)
		  {
		    /* Compute the distance between the two verticies where the
		     * path between them is constrained to the edges of the
		     * common shell, this will be DBL_MAX if the two verticies
		     * are in different shells. */
		    distM = WlzGMVertexShellDist(tMV, tMV0, NULL);
		    /* Compare this distance with the displacement of the
		     * transformed source vertex and then compute the new
		     * weight. */
		    if(distM > 0.0)
		    {
		      wgt += (distE > distM)?  distM / distE: distE / distM;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    wgt /= (3.0 * nScatter);
  }
  return(wgt);
}

/*!
* \return	Weight value in the range [0.0-1.0].
* \ingroup      WlzTransform
* \brief	Compute an additional weighting for the matched pair of
* 		3D verticies.
*		Each pair of matched verticies has an additional weighting
*		computed which relates to the sensitivity of the match to
*		perturbations.
* \param	curTr			Current affine transform.
* \param	tree			Given kD-tree populated by the
*					target vertices such that the nodes of
*					the tree have the sameindicies as the
*					given target vertices.
* \param	tVx			The target vertices.
* \param	sVx			The source vertices.
* \param	tMVx			The matched target vertex.
* \param	sMVx			The matched source vertex.
* \param	tGM			Target geometric model.
* \param	sGM			Source geometric model.
* \param	maxDisp			Maximum displacement for source
*					verticies.
* \param	nScatter		Number of verticies scattered about
*					the target vertex to test it's
*					sensitivity.
*/
static double	WlzMatchICPWeightMatches3D(WlzAffineTransform *curTr,
				       	AlcKDTTree *tree,
				       	WlzDVertex3 *tVx, WlzDVertex3 *sVx,
					WlzDVertex3 tMVx, WlzDVertex3 sMVx,
				       	WlzGMModel *tGM, WlzGMModel *sGM,
				       	double maxDisp, int nScatter)
{
  /* TODO */
}

/*!
* \return				Comparison value.
* \ingroup	WlzTransform
* \brief	Compare match shell sizes for qsort. Sorted shells will
*		have the largest shells first and the smallest shells
*		last.
* \param	val0			Used to pass first match shell.
* \param	val1			Used to pass second match shell.
*/
static int	WlzMatchICPShellSzCmp(const void *val0, const void *val1)
{
  int 		cmp;

  cmp = ((WlzMatchICPShell *)val1)->size - ((WlzMatchICPShell *)val0)->size;
  return(cmp);
}

/*!
* \return	Comparison value.
* \ingroup	WlzTransform
* \brief	Compare double values for AlgHeapSortIdx. 
* \param	data			Double data for softing.
* \param	idx			Indicies to be sorted.
* \param	id0			Index of index of first datum.
* \param	id1			Index of index of second datum.
*/
static int	WlzMatchICPDblSortFnD(void *data, int *idx, int id0, int id1)
{
  int		cmp;
  double	diff;

  diff = *((double *)data + *(idx + id0)) - *((double *)data + *(idx + id1));
  if(diff > DBL_EPSILON)
  {
    cmp = 1;
  }
  else if(diff < -(DBL_EPSILON))
  {
    cmp = -1;
  }
  else
  {
    cmp = 0.0;
  }
  return(cmp);
}

/*!
* \return	Comparison value.
* \ingroup	WlzTransform
* \brief	Compare the target points of matched points, sorting
*		in increasing vtY then vtX.
* \param	p0		Pointer to first matched points.
* \param	p1		Pointer to second matched points.
*/
static int	WlzMatchICPTPPairSortFnD(void *p0, void *p1)
{
  int		cmp = 0;
  WlzMatchICPTPPair2D *tPP0,
  		*tPP1;
  
  tPP0 = (WlzMatchICPTPPair2D *)p0;
  tPP1 = (WlzMatchICPTPPair2D *)p1;
  if(tPP0->tVx.vtY > (tPP1->tVx.vtY + DBL_EPSILON))
  {
    cmp = 1;
  }
  else if(tPP0->tVx.vtY < (tPP1->tVx.vtY - DBL_EPSILON))
  {
    cmp = -1;
  }
  else
  {
    if(tPP0->tVx.vtX > (tPP1->tVx.vtX + DBL_EPSILON))
    {
      cmp = 1;
    }
    else if(tPP0->tVx.vtX < (tPP1->tVx.vtX - DBL_EPSILON))
    {
      cmp = -1;
    }
  }
  return(cmp);
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Callback function which is called by the Woolz geometric model
*		functions when a new element has been created or is about to
*		be free'd. This function builds a list of new shells and
*		keeps count of both the number of new and freed shells.
* \param	model			The geometric model.
* \param	elm			The element that's been created or
*					is about to be destroyed.
* \param	reason			The reason for the callback.
* \param	data			Used to pass the callback data
*					structure with the list of new shells.
*/
static void	WlzMatchICPShellCb(WlzGMModel *model, WlzGMElemP elm,
			           WlzGMCbReason reason, void *data)
{
  WlzMatchICPShellListElm *tLElm,
  			  *nLElm;
  WlzMatchICPCbData *nSCb;

  if(elm.core && (elm.core->type == WLZ_GMELM_SHELL))
  {
    nSCb = (WlzMatchICPCbData *)data;
    if(reason == WLZ_GMCB_NEW)
    {
      if((nLElm = WlzMatchICPShellListElmNew()) == NULL)
      {
        nSCb->errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	++(nSCb->nNS);
	nLElm->mShell.shell = elm.shell;
	WlzMatchICPShellListElmAppend(nSCb->list, nLElm);
      }
    }
    else
    {
      /* Check that the shell that's about to be deleted isn't in the list
       * that's being built. If it is then remove it before it's deleted. */
      tLElm = nSCb->list->head;
      while(tLElm)
      {
        nLElm = tLElm->next;
	if(elm.shell == tLElm->mShell.shell)
	{
	  /* Need to remove the list entry before the shell is deleted. */
          WlzMatchICPShellListElmUnlink(nSCb->list, tLElm);
	  WlzMatchICPShellListElmFree(tLElm);
	}
	tLElm = nLElm;
      }
      ++(nSCb->nFS);
    }
  }
}

/*!
* \return				New shell match list or NULL on error.
* \ingroup	WlzTransform
* \brief	Allocate a new match list.
* \param	void
*/
static WlzMatchICPShellList *WlzMatchICPShellListNew(void)
{
  WlzMatchICPShellList *nList;

  nList  = (WlzMatchICPShellList *)
           AlcCalloc(1, sizeof(WlzMatchICPShellList));
  return(nList);
}

/*!
* \return	New shell match list element or NULL on error.
* \ingroup	WlzTransform
* \brief	Allocate a new match list element.
* \param	void
*/
static WlzMatchICPShellListElm *WlzMatchICPShellListElmNew(void)
{
  WlzMatchICPShellListElm *nLElm;

  nLElm  = (WlzMatchICPShellListElm *)
           AlcCalloc(1, sizeof(WlzMatchICPShellListElm));
  return(nLElm);
}

/*!
* \return	void.
* \ingroup	WlzTransform
* \brief	Free's a unlinked match list element.
* \param	void
*/
static void	WlzMatchICPShellListElmFree(WlzMatchICPShellListElm *elm)
{
  AlcFree(elm);
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Free's a whole list of match list elements along with the
*		list itself. Transforms and shells held in the list are NOT
*		freed.
* \param	void
*/
static void	WlzMatchICPShellListFree(WlzMatchICPShellList *list)
{
  WlzMatchICPShellListElm *elm0,
  		*elm1;

  elm0 = list->head;
  while(elm0)
  {
    elm1 = elm0->next;
    AlcFree(elm0);
    elm0 = elm1;
  }
  AlcFree(list);
}


/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Inserts a new element onto the head of the match shell list.
* \param	list			Match shell list.
* \param	nElm			New match shell list element.
*/
static void	WlzMatchICPShellListElmInsert(WlzMatchICPShellList *list,
					WlzMatchICPShellListElm *nElm)
{
  if(list && nElm)
  {
    if(list->head)
    {
      list->head->prev = nElm;
      nElm->next = list->head;
      nElm->prev = NULL;
      list->head = nElm;
    }
    else
    {
      list->head = list->tail = nElm;
    }
  }
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Appends a new element onto the end of the match shell list.
* \param	list			Match shell list.
* \param	nElm			New match shell list element.
*/
static void	WlzMatchICPShellListElmAppend(WlzMatchICPShellList *list,
					WlzMatchICPShellListElm *nElm)
{
  if(list && nElm)
  {
    if(list->tail)
    {
      list->tail->next = nElm;
      nElm->prev = list->tail;
      nElm->next = NULL;
      list->tail = nElm;
    }
    else
    {
      list->head = list->tail = nElm;
    }
  }
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Unlinks the given element from the list of shells.
* \param	list			List of match shells.
* \param 	elm			Element of the list of match shells
*					to be unlinked.
*/
static void	WlzMatchICPShellListElmUnlink(WlzMatchICPShellList *list,
					  WlzMatchICPShellListElm *elm)
{
  if(list && elm)
  {
    if(elm->prev)
    {
      elm->prev->next = elm->next;
    }
    else
    {
      list->head = elm->next;
    }
    if(elm->next)
    {
      elm->next->prev = elm->prev;
    }
    else
    {
      list->tail = elm->prev;
    }
    elm->prev = elm->next = NULL;
  }
}
