#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        WlzGeoModel.c
* Date:         February 2000
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Basic operators for handling both manifold and
*		non-manifold geometric models (GM) within Woolz.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 16-01-01 bill Add WlzGMVertexNormal3D().
* 12-12-00 bill	Fixed linked list bug in WlzGMModelAddVertex().
* 27-11-00 bill Add shell geometry updates during 2D model creation,
*		write WlzGMShellMergeG(). Fix bug allowing insertion
*		of duplicate simplicies in WlzGMModelConstructSimplex2D().
*		Fix bug in WlzGMModelDeleteS().
* 21-11-00 bill	Fix shell leak in WlzGMModelJoinL2D().
* 16-11-00 bill Add WlzGMModelDeleteS().
* 12-10-00 bill Fix bug counting shells in WlzGMModelConstructSimplex3D().
* 10-10-00 bill Add WlzGMModelCopy().
* 25-08-00 bill Fix bugs in 3D model construction. Remove element's
*		data and flags.
* 10-08-00 bill Do some optimization so that when joining shells, the
*		shell who's geometry has the greatest bounding box
*		volume is retained.
************************************************************************/
#include <Wlz.h>
#include <limits.h>

static WlzErrorNum 	WlzGMShellMergeG(
			  WlzGMShell *shell0,
			  WlzGMShell *shell1);
static WlzErrorNum 	WlzGMShellTestOutVTK(
			  WlzGMShell *shell,
			  int *vIdxTb,
			  char *lFlg,
			  FILE *fP);
static WlzErrorNum 	WlzGMLoopTTestOutputVTK(
			  WlzGMLoopT *fLT,
			  int *vIdxTb,
			  FILE *fP);
static WlzErrorNum 	WlzGMShellTestOutPS(
			  WlzGMShell *shell,
			  char *lFlg,
			  char *eFlg,
			  FILE *fP,
			  WlzDVertex2 offset,
			  WlzDVertex2 scale);
static WlzErrorNum	WlzGMLoopTTestOutPS(
			  WlzGMLoopT *loopT,
			  char *eFlg,
			  FILE *fP,
			  WlzDVertex2 offset,
			  WlzDVertex2 scale);
static WlzErrorNum	WlzGMEdgeTTestOutPS(
			  WlzGMEdgeT *edgeT,
			  FILE *fP,
			  WlzDVertex2 offset,
			  WlzDVertex2 scale);
static void		WlzGMTestOutLinePS(
			  FILE *fP,
			  WlzDVertex2 pos0,
			  WlzDVertex2 pos1,
			  WlzDVertex2 offset,
			  WlzDVertex2 scale);
static unsigned int 	WlzGMHashPos2D(
			  WlzDVertex2 pos);
static unsigned int 	WlzGMHashPos3D(
			  WlzDVertex3 pos);
static void		WlzGMModelMatchEdgeTG2D(
			  WlzGMModel *model,
			  WlzGMEdgeT **matchET,
			  WlzDVertex2 *pos);
static WlzErrorNum      WlzGMModelConstructNewS3D(
                          WlzGMModel *model,
                          WlzDVertex3 *pos);
static WlzErrorNum      WlzGMModelExtend1V0E1S3D(
                          WlzGMModel *model,
                          WlzGMVertex *eV,
                          WlzDVertex3 pos0,
                          WlzDVertex3 pos1);
static WlzErrorNum      WlzGMModelExtend2V1E1S3D(
                          WlzGMModel *model,
                          WlzGMEdge *eE,
                          WlzDVertex3 pos2);
static WlzErrorNum      WlzGMModelExtend2V0E1S3D(
                          WlzGMModel *model,
                          WlzGMVertex *eV0,
                          WlzGMVertex *eV1,
                          WlzDVertex3 pos2);
static WlzErrorNum      WlzGMModelJoin2V0E0S3D(
                          WlzGMModel *model,
                          WlzGMVertex *eV0,
                          WlzGMVertex *eV1,
                          WlzDVertex3 pos2);
static WlzErrorNum      WlzGMModelJoin3V0E3S3D(
                          WlzGMModel *model,
                          WlzGMVertex **eV);
static WlzErrorNum      WlzGMModelJoin3V0E2S3D(
                          WlzGMModel *model,
                          WlzGMVertex **eV);
static WlzErrorNum      WlzGMModelExtend3V0E1S3D(
                          WlzGMModel *model,
                          WlzGMVertex **eV);
static WlzErrorNum      WlzGMModelExtend3V1E1S3D(
                          WlzGMModel *model,
                          WlzGMEdge *eE,
                          WlzGMVertex *eV);
static WlzErrorNum      WlzGMModelJoin3V1E2S3D(
                          WlzGMModel *model,
                          WlzGMEdge *eE,
                          WlzGMVertex *eV);
static WlzErrorNum      WlzGMModelExtend3V2E1S3D(
                          WlzGMModel *model,
                          WlzGMEdge *eE0,
                          WlzGMEdge *eE1);
static WlzErrorNum      WlzGMModelConstructNewS2D(
                          WlzGMModel *model,
                          WlzDVertex2 *pos);
static WlzErrorNum      WlzGMModelExtendL2D(
                          WlzGMModel *model,
                          WlzGMEdgeT *eET,
                          WlzDVertex2 nPos);
static WlzErrorNum      WlzGMModelConstructSplitL2D(
                          WlzGMModel *model,
                          WlzGMEdgeT *eET0,
                          WlzGMEdgeT *eET1);
static WlzErrorNum      WlzGMModelJoinL2D(
                          WlzGMModel *model,
                          WlzGMEdgeT *eET0,
                          WlzGMEdgeT *eET1);
static void             WlzGMLoopTSetT(
                          WlzGMLoopT *eLT);
static void             WlzGMModelAddVertex(
                          WlzGMModel *model,
                          WlzGMVertex *nV);
static void		WlzGMElmMarkFree(
			  int *idxP);

/* Creation  of geometric modeling elements */

/************************************************************************
* Function:	WlzGMModelNew
* Returns:	WlzGMModel *:		New empty model.
* Purpose:	Creates an empty non-manifold geometry model.
* Global refs:	-
* Parameters:	WlzGMModelType modType:	Type of model to create.
*		int blkSz:		Resource block size, used for
*					allocating storage for model
*					elements. A default size is
*					used if <= 0.
*		int vHTSz:		Vertex matching hash table size,
*					A default size is used if <= 0.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMModel	*WlzGMModelNew(WlzGMModelType modType,
			       int blkSz, int vHTSz, WlzErrorNum *dstErr)
{
  WlzGMModel	*model = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  unsigned int	vertexGSz,
		shellGSz;
  const unsigned int defBlkSz = 1024;
  const unsigned int defVHTSz = 1024 * 1024;

  if(blkSz <= 0)
  {
    blkSz = defBlkSz;
  }
  if(vHTSz <= 0)
  {
    vHTSz = defVHTSz;
  }
  switch(modType)
  {
    case WLZ_GMMOD_2I:
      vertexGSz = sizeof(WlzGMVertexG2I);
      shellGSz = sizeof(WlzGMShellG2I);
      break;
    case WLZ_GMMOD_2D:
      vertexGSz = sizeof(WlzGMVertexG2D);
      shellGSz = sizeof(WlzGMShellG2D);
      break;
    case WLZ_GMMOD_3I:
      vertexGSz = sizeof(WlzGMVertexG3I);
      shellGSz = sizeof(WlzGMShellG3I);
      break;
    case WLZ_GMMOD_3D:
      vertexGSz = sizeof(WlzGMVertexG3D);
      shellGSz = sizeof(WlzGMShellG3D);
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  /* Create the new model. All elements of the model including it's
   * resources are set to zero. */
  if((errNum == WLZ_ERR_NONE) &&
     ((model = AlcCalloc(1, sizeof(WlzGMModel))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  /* Create initial model resources */
  if(errNum == WLZ_ERR_NONE)
  {
    model->type = modType;
    if(((model->res.vertex.vec = AlcVectorNew(1, sizeof(WlzGMVertex),
    					      blkSz, NULL)) == NULL) ||
       ((model->res.vertexT.vec = AlcVectorNew(1, sizeof(WlzGMVertexT),
    					       blkSz, NULL)) == NULL) ||
       ((model->res.vertexG.vec = AlcVectorNew(1, vertexGSz,
    					       blkSz, NULL)) == NULL) ||
       ((model->res.diskT.vec = AlcVectorNew(1, sizeof(WlzGMDiskT),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.edge.vec = AlcVectorNew(1, sizeof(WlzGMEdge),
    					    blkSz, NULL)) == NULL) ||
       ((model->res.edgeT.vec = AlcVectorNew(1, sizeof(WlzGMEdgeT),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.loop.vec = AlcVectorNew(1, sizeof(WlzGMLoop),
    					    blkSz, NULL)) == NULL) ||
       ((model->res.loopT.vec = AlcVectorNew(1, sizeof(WlzGMLoopT),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.shell.vec = AlcVectorNew(1, sizeof(WlzGMShell),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.shellG.vec = AlcVectorNew(1, shellGSz,
    					      blkSz, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Create vertex hash table. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((model->vertexHT = (WlzGMVertex **)
    			  AlcCalloc(vHTSz, sizeof(WlzGMVertex *))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      model->vertexHTSz = vHTSz;
    }
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFree(model);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(model);
}

/************************************************************************
* Function:	WlzGMModelNewS
* Returns:	WlzGMShell *:		New empty shell.
* Purpose:	Creates an empty shell with a shell geometry element
*		for the model.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The geometric model.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMShell	*WlzGMModelNewS(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resS,
  		*resSG;
  WlzGMShell	*shell = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new shell and shell geometry element from the shell */
    resS = &(model->res.shell);
    resSG = &(model->res.shellG);
    if(((shell = (WlzGMShell *)
     	        (AlcVectorExtendAndGet(resS->vec, resS->numIdx))) == NULL) ||
       ((shell->geo.core = (WlzGMCore *)
    			   (AlcVectorExtendAndGet(resSG->vec,
					 	  resSG->numIdx))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resS->numElm);
    ++(resSG->numElm);
    shell->type = WLZ_GMELM_SHELL;
    shell->idx = (resS->numIdx)++;
    shell->geo.core->type = WlzGMModelGetSGeomType(model);
    shell->geo.core->idx = (resSG->numIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeS(model, shell);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(shell);
}

/************************************************************************
* Function:	WlzGMModelNewL
* Returns:	WlzGMLoop *:		New loop.
* Purpose:	Creates an empty loop with a loop geometry element.
*		The loop geometry element only has it's index set
*		to a meaningful value.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Parent model.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMLoop	*WlzGMModelNewL(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resL;
  WlzGMLoop	*loop = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new loop and loop geometry element from the model */
    resL = &(model->res.loop);
    if((loop = (WlzGMLoop *)
     	       (AlcVectorExtendAndGet(resL->vec, resL->numIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resL->numElm);
    loop->type = WLZ_GMELM_LOOP;
    loop->idx = (resL->numIdx)++;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(loop);
}

/************************************************************************
* Function:	WlzGMModelNewLT
* Returns:	WlzGMLoopT *:		New loop.
* Purpose:	Creates a loop topology element.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Parent model.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMLoopT	*WlzGMModelNewLT(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resLT;
  WlzGMLoopT	*loopT = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new loop topology element from the model */
    resLT = &(model->res.loopT);
    if((loopT = (WlzGMLoopT *)
     	        (AlcVectorExtendAndGet(resLT->vec, resLT->numIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resLT->numElm);
    loopT->type = WLZ_GMELM_LOOP_T;
    loopT->idx = (resLT->numIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeLT(model, loopT);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(loopT);
}

/************************************************************************
* Function:	WlzGMModelNewDT
* Returns:	WlzGMDiskT:		New disk topology element.
* Purpose:	Creates a new disk topology element.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzErrorNum *dstErr: 	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMDiskT     *WlzGMModelNewDT(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resDT;
  WlzGMDiskT	*diskT = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new disk topology element from the model */
    resDT = &(model->res.diskT);
    if((diskT = (WlzGMDiskT *)
     	        (AlcVectorExtendAndGet(resDT->vec, resDT->numIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resDT->numElm);
    diskT->type = WLZ_GMELM_DISK_T;
    diskT->idx = (resDT->numIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeDT(model, diskT);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(diskT);
}

/************************************************************************
* Function:	WlzGMModelNewE
* Returns:	WlzGMEdge:		New edge.
* Purpose:	Creates a new edge and an edge geometry element.
*		The edge geometry element only has it's index set
*		to a meaningful value.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzErrorNum *dstErr: 	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMEdge      *WlzGMModelNewE(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resE;
  WlzGMEdge	*edge = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new edge and edge geometry element from the model */
    resE = &(model->res.edge);
    if((edge = (WlzGMEdge *)
     	       (AlcVectorExtendAndGet(resE->vec, resE->numIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resE->numElm);
    edge->type = WLZ_GMELM_EDGE;
    edge->idx = (model->res.edge.numIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeE(model, edge);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(edge);
}

/************************************************************************
* Function:	WlzGMModelNewET
* Returns:	WlzGMEdgeT:		New edge topology element.
* Purpose:	Creates a new edge topology element.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzErrorNum *dstErr: 	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMEdgeT     *WlzGMModelNewET(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resET;
  WlzGMEdgeT	*edgeT = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new edge topology element from the model */
    resET = &(model->res.edgeT);
    if((edgeT = (WlzGMEdgeT *)
     	        (AlcVectorExtendAndGet(resET->vec, resET->numIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resET->numElm);
    edgeT->type = WLZ_GMELM_EDGE_T;
    edgeT->idx = (resET->numIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeET(model, edgeT);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(edgeT);
}

/************************************************************************
* Function:	WlzGMModelNewV
* Returns:	WlzGMVertex:		New vertex.
* Purpose:	Creates a new vertex and a vertex geometry element.
*		The vertex geometry element only has it's index set
*		to a meaningful value.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzErrorNum *dstErr: 	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMVertex      *WlzGMModelNewV(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resV,
  		*resVG;
  WlzGMVertex	*vertex = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new vertex and vertex geometry element from the model */
    resV = &(model->res.vertex);
    resVG = &(model->res.vertexG);
    if(((vertex = (WlzGMVertex *)
     	        (AlcVectorExtendAndGet(resV->vec, resV->numIdx))) == NULL) ||
       ((vertex->geo.core = (WlzGMCore *)
    			    (AlcVectorExtendAndGet(resVG->vec,
					 	   resVG->numIdx))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resV->numElm);
    ++(resVG->numElm);
    vertex->type = WLZ_GMELM_VERTEX;
    vertex->idx = (resV->numIdx)++;
    vertex->geo.core->type = WlzGMModelGetVGeomType(model);
    vertex->geo.core->idx = (resVG->numIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeV(model, vertex);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vertex);
}

/************************************************************************
* Function:	WlzGMModelNewVT
* Returns:	WlzGMVertexT:		New vertex topology element.
* Purpose:	Creates a new vertex topology element.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzErrorNum *dstErr: 	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMVertexT      *WlzGMModelNewVT(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resVT;
  WlzGMVertexT	*vertexT = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new vertex topology element from the model */
    resVT = &(model->res.vertexT);
    if((vertexT = (WlzGMVertexT *)
     	          (AlcVectorExtendAndGet(resVT->vec, resVT->numIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resVT->numElm);
    vertexT->type = WLZ_GMELM_VERTEX_T;
    vertexT->idx = (resVT->numIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeVT(model, vertexT);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vertexT);
}

/************************************************************************
* Function:     WlzGMModelCopy
* Returns:      WlzGMModel *:		New copy of the given model.
* Purpose:	Copies the given model. All unused elements are squeezed
*		out.
* Global refs:	-
* Parameters:	WlzGMModel *gM:	Given model.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMModel 	*WlzGMModelCopy(WlzGMModel *gM, WlzErrorNum *dstErr)
{
  int		gIdx,
  		nIdx,
  		nBkSz,
		nHTSz;
  WlzGMElemP	gElmP,
  		nElmP;
  WlzGMResIdxTb *resIdxTb = NULL;
  const int	minBkSz = 1024,
  		minHTSz = 1024;
  WlzGMModel 	*nM = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(gM == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Make a new model. */
    if((nBkSz = gM->res.vertex.numElm / 16) < minBkSz)
    {
      nBkSz = minBkSz;
    }
    if((nHTSz = gM->res.vertex.numElm / 8) < minHTSz)
    {
      nHTSz = minHTSz;
    }
    nM = WlzGMModelNew(gM->type, nBkSz, nHTSz, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute element index tables to be used in setting up the elements
     * of the new model while sqeezing out the unused elements. */
    resIdxTb = WlzGMModelResIdx(gM,
				WLZ_GMELMFLG_VERTEX |
				WLZ_GMELMFLG_VERTEX_G |
				WLZ_GMELMFLG_VERTEX_T |
				WLZ_GMELMFLG_DISK_T |
				WLZ_GMELMFLG_EDGE |
				WLZ_GMELMFLG_EDGE_T |
				WLZ_GMELMFLG_LOOP |
				WLZ_GMELMFLG_LOOP_T |
				WLZ_GMELMFLG_SHELL |
				WLZ_GMELMFLG_SHELL_G,
    				&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy verticies. */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.vertex.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.vertex.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->vertex.idxLut + gIdx);
        if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.vertex.vec,
							 nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = WLZ_GMELM_VERTEX;
	  nElmP.core->idx = nIdx;
	  nElmP.vertex->diskT = (WlzGMDiskT *)
	  			AlcVectorExtendAndGet(nM->res.diskT.vec,
						 *(resIdxTb->diskT.idxLut + 
						   gElmP.vertex->diskT->idx));
	  nElmP.vertex->geo.core = (WlzGMCore *)
	  			   AlcVectorExtendAndGet(nM->res.vertexG.vec,
						 *(resIdxTb->vertexG.idxLut + 
						 gElmP.vertex->geo.core->idx));
	  if(gElmP.vertex->next)
	  {
	    nElmP.vertex->next = (WlzGMVertex *)
				 AlcVectorExtendAndGet(nM->res.vertex.vec,
						  *(resIdxTb->vertex.idxLut +
						    gElmP.vertex->next->idx));
	  }
	  else
	  {
	    nElmP.vertex->next = NULL;
	  }
	  if((nElmP.vertex->diskT == NULL) ||
	     (nElmP.vertex->geo.core == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy vertex geometries */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.vertexG.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.vertexG.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->vertexG.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.vertexG.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  switch(nElmP.core->type)
	  {
	    case WLZ_GMELM_VERTEX_G2I:
	      nElmP.vertexG2I->vtx = gElmP.vertexG2I->vtx;
	      break;
	    case WLZ_GMELM_VERTEX_G2D:
	      nElmP.vertexG2D->vtx = gElmP.vertexG2D->vtx;
	      break;
	    case WLZ_GMELM_VERTEX_G3I:
	      nElmP.vertexG3I->vtx = gElmP.vertexG3I->vtx;
	      break;
	    case WLZ_GMELM_VERTEX_G3D:
	      nElmP.vertexG3D->vtx = gElmP.vertexG3D->vtx;
	      break;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy vertexT's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.vertexT.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.vertexT.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->vertexT.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.vertexT.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  nElmP.vertexT->next = (WlzGMVertexT *)
	  			AlcVectorExtendAndGet(nM->res.vertexT.vec,
						 *(resIdxTb->vertexT.idxLut +
						   gElmP.vertexT->next->idx));
	  nElmP.vertexT->prev = (WlzGMVertexT *)
	  			AlcVectorExtendAndGet(nM->res.vertexT.vec,
						 *(resIdxTb->vertexT.idxLut +
						   gElmP.vertexT->prev->idx));
	  nElmP.vertexT->diskT = (WlzGMDiskT *)
	  			 AlcVectorExtendAndGet(nM->res.diskT.vec,
						 *(resIdxTb->diskT.idxLut +
						   gElmP.vertexT->diskT->idx));
	  nElmP.vertexT->parent = (WlzGMEdgeT *)
	  			  AlcVectorExtendAndGet(nM->res.edgeT.vec,
						 *(resIdxTb->edgeT.idxLut +
						  gElmP.vertexT->parent->idx));
	  if((nElmP.vertexT->next == NULL) ||
	     (nElmP.vertexT->prev == NULL) ||
	     (nElmP.vertexT->diskT == NULL) ||
	     (nElmP.vertexT->parent == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy diskT's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.diskT.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.diskT.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->diskT.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.diskT.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  nElmP.diskT->next = (WlzGMDiskT *)
	  		      AlcVectorExtendAndGet(nM->res.diskT.vec,
					       *(resIdxTb->diskT.idxLut +
						 gElmP.diskT->next->idx));
	  nElmP.diskT->prev = (WlzGMDiskT *)
			      AlcVectorExtendAndGet(nM->res.diskT.vec,
					       *(resIdxTb->diskT.idxLut +
					         gElmP.diskT->prev->idx));
	  nElmP.diskT->vertex = (WlzGMVertex *)
	  			 AlcVectorExtendAndGet(nM->res.vertex.vec,
						 *(resIdxTb->vertex.idxLut +
						   gElmP.diskT->vertex->idx));
	  nElmP.diskT->vertexT = (WlzGMVertexT *)
	  			  AlcVectorExtendAndGet(nM->res.vertexT.vec,
						 *(resIdxTb->vertexT.idxLut +
						  gElmP.diskT->vertexT->idx));
	  if((nElmP.diskT->next == NULL) ||
	     (nElmP.diskT->prev == NULL) ||
	     (nElmP.diskT->vertex == NULL) ||
	     (nElmP.diskT->vertexT == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy edge's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.edge.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.edge.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->edge.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.edge.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  nElmP.edge->edgeT = (WlzGMEdgeT *)
	  		      AlcVectorExtendAndGet(nM->res.edgeT.vec,
					       *(resIdxTb->edgeT.idxLut +
						 gElmP.edge->edgeT->idx));
	  if(nElmP.edge->edgeT == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy edgeT's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.edgeT.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.edgeT.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->edgeT.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.edgeT.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  nElmP.edgeT->next = (WlzGMEdgeT *)
	  		       AlcVectorExtendAndGet(nM->res.edgeT.vec,
					       *(resIdxTb->edgeT.idxLut +
						 gElmP.edgeT->next->idx));
	  nElmP.edgeT->prev = (WlzGMEdgeT *)
	  		       AlcVectorExtendAndGet(nM->res.edgeT.vec,
					       *(resIdxTb->edgeT.idxLut +
						 gElmP.edgeT->prev->idx));
	  nElmP.edgeT->opp = (WlzGMEdgeT *)
	  		     AlcVectorExtendAndGet(nM->res.edgeT.vec,
					      *(resIdxTb->edgeT.idxLut +
					        gElmP.edgeT->opp->idx));
	  nElmP.edgeT->rad = (WlzGMEdgeT *)
	  		       AlcVectorExtendAndGet(nM->res.edgeT.vec,
					       *(resIdxTb->edgeT.idxLut +
						 gElmP.edgeT->rad->idx));
	  nElmP.edgeT->edge = (WlzGMEdge *)
	  		      AlcVectorExtendAndGet(nM->res.edge.vec,
					       *(resIdxTb->edge.idxLut +
						 gElmP.edgeT->edge->idx));
	  nElmP.edgeT->vertexT = (WlzGMVertexT *)
	  		         AlcVectorExtendAndGet(nM->res.vertexT.vec,
					          *(resIdxTb->vertexT.idxLut +
						   gElmP.edgeT->vertexT->idx));
	  nElmP.edgeT->parent = (WlzGMLoopT *)
	  		        AlcVectorExtendAndGet(nM->res.loopT.vec,
					       *(resIdxTb->loopT.idxLut +
						 gElmP.edgeT->parent->idx));

	  if((nElmP.edgeT->next == NULL) || 
	     (nElmP.edgeT->prev == NULL) || 
	     (nElmP.edgeT->opp == NULL) || 
	     (nElmP.edgeT->rad == NULL) || 
	     (nElmP.edgeT->edge == NULL) || 
	     (nElmP.edgeT->vertexT == NULL) || 
	     (nElmP.edgeT->parent == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy loop's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.loop.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.loop.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->loop.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.loop.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  nElmP.loop->loopT = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec,
					       *(resIdxTb->loopT.idxLut +
						 gElmP.loop->loopT->idx));
	  if(nElmP.loop->loopT == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy loopT's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.loopT.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.loopT.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->loopT.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.loopT.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  nElmP.loopT->next = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec,
					       *(resIdxTb->loopT.idxLut +
						 gElmP.loopT->next->idx));
	  nElmP.loopT->prev = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec,
					       *(resIdxTb->loopT.idxLut +
						 gElmP.loopT->prev->idx));
	  nElmP.loopT->opp = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec,
					       *(resIdxTb->loopT.idxLut +
						 gElmP.loopT->opp->idx));
	  nElmP.loopT->loop = (WlzGMLoop *)
	  		      AlcVectorExtendAndGet(nM->res.loop.vec,
					       *(resIdxTb->loop.idxLut +
						 gElmP.loopT->loop->idx));
	  nElmP.loopT->edgeT = (WlzGMEdgeT *)
	  		      AlcVectorExtendAndGet(nM->res.edgeT.vec,
					       *(resIdxTb->edgeT.idxLut +
						 gElmP.loopT->edgeT->idx));
	  nElmP.loopT->parent = (WlzGMShell *)
				AlcVectorExtendAndGet(nM->res.shell.vec,
						 *(resIdxTb->shell.idxLut +
						   gElmP.loopT->parent->idx));
	  if((nElmP.loopT->next == NULL) ||
	     (nElmP.loopT->prev == NULL) ||
	     (nElmP.loopT->opp == NULL) ||
	     (nElmP.loopT->loop == NULL) ||
	     (nElmP.loopT->edgeT == NULL) ||
	     (nElmP.loopT->parent == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy shell's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.shell.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.shell.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->shell.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.shell.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  nElmP.shell->next = (WlzGMShell *)
	  		      AlcVectorExtendAndGet(nM->res.shell.vec,
					       *(resIdxTb->shell.idxLut +
						 gElmP.shell->next->idx));
	  nElmP.shell->prev = (WlzGMShell *)
	  		      AlcVectorExtendAndGet(nM->res.shell.vec,
					       *(resIdxTb->shell.idxLut +
						 gElmP.shell->prev->idx));
	  nElmP.shell->geo.core = (WlzGMCore *)
	  		      AlcVectorExtendAndGet(nM->res.shellG.vec,
					       *(resIdxTb->shellG.idxLut +
						 gElmP.shell->geo.core->idx));
	  nElmP.shell->child = (WlzGMLoopT *)
	  		       AlcVectorExtendAndGet(nM->res.loopT.vec,
					       *(resIdxTb->loopT.idxLut +
						 gElmP.shell->child->idx));
	  nElmP.shell->parent = nM;
	  if((nElmP.shell->next == NULL) ||
	     (nElmP.shell->prev == NULL) ||
	     (nElmP.shell->geo.core == NULL) ||
	     (nElmP.shell->child == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy shellG's */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.shellG.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.shellG.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->shellG.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.shellG.vec,
						       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  switch(nElmP.core->type)
	  {
	    case WLZ_GMELM_SHELL_G2I:
	      nElmP.shellG2I->bBox = gElmP.shellG2I->bBox;
	      break;
	    case WLZ_GMELM_SHELL_G2D:
	      nElmP.shellG2D->bBox = gElmP.shellG2D->bBox;
	      break;
	    case WLZ_GMELM_SHELL_G3I:
	      nElmP.shellG3I->bBox = gElmP.shellG3I->bBox;
	      break;
	    case WLZ_GMELM_SHELL_G3D:
	      nElmP.shellG3D->bBox = gElmP.shellG3D->bBox;
	      break;
	  }
        }
      }
      ++gIdx;
    }
  }
  if(resIdxTb)
  {
    WlzGMModelResIdxFree(resIdxTb);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nM->child = (WlzGMShell *)AlcVectorItemGet(nM->res.shell.vec,
					       *(resIdxTb->shell.idxLut +
						 gM->child->idx));
    nM->res.vertex.numElm = nM->res.vertex.numIdx = gM->res.vertex.numElm;
    nM->res.vertexG.numElm = nM->res.vertexG.numIdx = gM->res.vertexG.numElm;
    nM->res.vertexT.numElm = nM->res.vertexT.numIdx = gM->res.vertexT.numElm;
    nM->res.diskT.numElm = nM->res.diskT.numIdx = gM->res.diskT.numElm;
    nM->res.edge.numElm = nM->res.edge.numIdx = gM->res.edge.numElm;
    nM->res.edgeT.numElm = nM->res.edgeT.numIdx = gM->res.edgeT.numElm;
    nM->res.loop.numElm = nM->res.loop.numIdx = gM->res.loop.numElm;
    nM->res.loopT.numElm = nM->res.loopT.numIdx = gM->res.loopT.numElm;
    nM->res.shell.numElm = nM->res.shell.numIdx = gM->res.shell.numElm;
    nM->res.shellG.numElm = nM->res.shellG.numIdx = gM->res.shellG.numElm;
    errNum = WlzGMModelRehashVHT(nM, 0);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFree(nM);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nM);
}

/* Freeing  of geometric modeling elements */

/************************************************************************
* Function:	WlzGMModelFree
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Free's a geometric model and it's elements.
* Global refs:	-
* Parameters:	WlzGMModel *shell:	The shell to free.
************************************************************************/
WlzErrorNum	WlzGMModelFree(WlzGMModel *model)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )AlcVectorFree(model->res.vertex.vec);
    (void )AlcVectorFree(model->res.vertexT.vec);
    (void )AlcVectorFree(model->res.vertexG.vec);
    (void )AlcVectorFree(model->res.diskT.vec);
    (void )AlcVectorFree(model->res.edge.vec);
    (void )AlcVectorFree(model->res.edgeT.vec);
    (void )AlcVectorFree(model->res.loop.vec);
    (void )AlcVectorFree(model->res.loopT.vec);
    (void )AlcVectorFree(model->res.shell.vec);
    (void )AlcVectorFree(model->res.shellG.vec);
    if(model->vertexHT)
    {
      AlcFree(model->vertexHT);
    }
    AlcFree((void *)model);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMElmMarkFree
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks an element as free by setting it's index to
*		a -ve value, the actual choice of value aids debugging.
* Global refs:	-
* Parameters:	int *idxP:		Pointer to elements index.
************************************************************************/
static void	WlzGMElmMarkFree(int *idxP)
{
  if(*idxP >= 0)
  {
    *idxP = -(*idxP + 1);
  }
}

/************************************************************************
* Function:	WlzGMModelFreeS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a shell and it's geometry as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMShell *shell:	Shell to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeS(WlzGMModel *model, WlzGMShell *shell)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (shell == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the loop so just mark it and it's geometry element
     * as invalid by making the indicies < 0. */
    WlzGMElmMarkFree(&(shell->idx));
    if(shell->geo.core != NULL)
    {
      WlzGMElmMarkFree(&(shell->geo.core->idx));
      --(model->res.shellG.numElm);
    }
    --(model->res.shell.numElm);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeL
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a loop as invalid and suitable for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMLoop *loop:	Loop to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeL(WlzGMModel *model, WlzGMLoop *loop)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (loop == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the loop so just mark it as invalid by making
     * the indicies < 0. */
    WlzGMElmMarkFree(&(loop->idx));
    --(model->res.loop.numElm);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeLT
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a loop topology as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMLoopT *loopT:	LoopT to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeLT(WlzGMModel *model, WlzGMLoopT *loopT)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (loopT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the loopT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(loopT->idx));
    --(model->res.loopT.numElm);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeDT
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks an disk topology as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMDiskT *diskT:	DiskT to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeDT(WlzGMModel *model, WlzGMDiskT *diskT)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (diskT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the diskT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(diskT->idx));
    --(model->res.diskT.numElm);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeE
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks an edge as invalid and suitable for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMEdge *edge:	Edge to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeE(WlzGMModel *model, WlzGMEdge *edge)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (edge == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the edge so just mark it as invalid by making
     * the indicies < 0. */
    WlzGMElmMarkFree(&(edge->idx));
    --(model->res.edge.numElm);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeET
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks an edge topology as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMEdgeT *edgeT:	EdgeT to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeET(WlzGMModel *model, WlzGMEdgeT *edgeT)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (edgeT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the edgeT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(edgeT->idx));
    --(model->res.edgeT.numElm);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeV
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a vertex and it's geometry as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMVertex *vertex:	Vertex to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeV(WlzGMModel *model, WlzGMVertex *vertex)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (vertex == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the vertex so just mark it and it's geometry element
     * as invalid by making the indicies < 0. */
    WlzGMElmMarkFree(&(vertex->idx));
    if(vertex->geo.core != NULL)
    {
      WlzGMElmMarkFree(&(vertex->geo.core->idx));
      --(model->res.vertexG.numElm);
    }
    --(model->res.vertex.numElm);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeVT
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a vertex topology as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMVertexT *vertexT:	VertexT to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeVT(WlzGMModel *model, WlzGMVertexT *vertexT)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (vertexT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the vertexT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(vertexT->idx));
    --(model->res.vertexT.numElm);
  }
  return(errNum);
}

/* Deletion of geometric modeling elements along with children */

/************************************************************************
* Function:	WlzGMModelDeleteS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Deletes a shell by unlinking it and then freeing it and
*		all it's children.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMShell *dS:		Shell to delete.
************************************************************************/
WlzErrorNum	WlzGMModelDeleteS(WlzGMModel *model, WlzGMShell *dS)
{
  WlzGMVertexT	*tVT;
  WlzGMEdgeT	*fET,
  		*nET,
		*tET;
  WlzGMLoopT	*fLT,
  		*nLT,
  		*tLT;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (dS == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Unlink the shell. */
    WlzGMShellUnlink(dS);
    /* For each loopT. */
    fLT = tLT = dS->child;
    do
    {
      /* For each edgeT. */
      fET = tET = tLT->edgeT;
      do
      {
	tVT = tET->vertexT;
	if(tVT->diskT->vertex->idx >= 0)
	{
	  (void )WlzGMModelFreeV(model, tVT->diskT->vertex);
	}
	if(tVT->diskT->idx >= 0)
	{
	  (void )WlzGMModelFreeDT(model, tVT->diskT);
	}
	if(tVT->idx >= 0)
	{
	  (void )WlzGMModelFreeVT(model, tVT);
	}
	if(tET->rad->idx >= 0)
	{
	  (void )WlzGMModelFreeET(model, tET->rad);
	}
	if(tET->edge->idx >= 0)
	{
	  (void )WlzGMModelFreeE(model, tET->edge);
        }
	nET = tET->next;
	if(tET->idx >= 0)
	{
	  (void )WlzGMModelFreeET(model, tET);
	}
	tET = nET;
      } while(tET != fET);
      /* Free the loopT and loop. */
      nLT = tLT->next;
      if(tLT->loop->idx >= 0)
      {
        (void )WlzGMModelFreeL(model, tLT->loop);
      }
      (void )WlzGMModelFreeLT(model, tLT);
      tLT = nLT;
    } while(tLT != fLT);
    (void )WlzGMModelFreeS(model, dS);
  }
  return(errNum);
}

/* Model access and testing. */

/************************************************************************
* Function:	WlzGMModelTypeValid
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Checks the model type is valid.
* Global refs:	-
* Parameters:	WlzGMModelType type:	Model type to check.
************************************************************************/
WlzErrorNum 	WlzGMModelTypeValid(WlzGMModelType type)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(type)
  {
    case WLZ_GMMOD_2I: /* FALLTHROUGH */
    case WLZ_GMMOD_2D: /* FALLTHROUGH */
    case WLZ_GMMOD_3I: /* FALLTHROUGH */
    case WLZ_GMMOD_3D:
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelTestOutVTK
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 3D model to the specified file using the
*		VTK polydata format.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Given model.
*		FILE *fP:		File ptr, may be NULL if
*					no output required.
************************************************************************/
WlzErrorNum   	WlzGMModelTestOutVTK(WlzGMModel *model, FILE *fP)
{
  int		idI,
  		idV,
		cnt,
		vCntMax;
  char		*lFlg = NULL;
  int		*vIdxTb = NULL;
  WlzDVertex3	vtx;
  AlcVector	*vec;
  WlzGMVertex	*vP;
  WlzGMShell	*fS,
  		*tS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(model->type)
    {
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate shell, loop and edge flags. */
    if((lFlg = (char *)AlcCalloc(model->res.loop.numIdx,
    				 sizeof(char))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate a vertex index table. */
    vCntMax = model->res.vertex.numIdx;
    if((vIdxTb = (int *)AlcCalloc(vCntMax, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the file header. */
    (void )fprintf(fP,
    		   "# vtk DataFile Version 1.0\n"
    		   "WlzGeoModel test output\n"
		   "ASCII\n"
		   "DATASET POLYDATA\n"
		   "POINTS %d float\n",
		   model->res.vertex.numElm);
    /* Output the verticies while building the index table. */
    idI = idV = 0;
    vec = model->res.vertex.vec;
    for(idI = 0; idI < model->res.vertex.numIdx; ++idI)
    {
      vP = (WlzGMVertex *)AlcVectorItemGet(vec, idI);
      if(vP->idx >= 0)
      {
	*(vIdxTb + vP->idx) = idV++;
	(void )WlzGMVertexGetG3D(vP, &vtx);
	(void )fprintf(fP, "%g %g %g\n", vtx.vtX, vtx.vtY, vtx.vtZ);
      }
    }
    if((tS = fS = model->child) != NULL)
    {
      (void )fprintf(fP,
      		     "POLYGONS %d %d\n",
		     model->res.loop.numElm, 4 * model->res.loop.numElm);
      do
      {
	if(tS->idx >= 0)
	{
	  errNum = WlzGMShellTestOutVTK(tS, vIdxTb, lFlg, fP);
	}
	tS = tS->next;
      } while((errNum == WLZ_ERR_NONE) && (tS->idx != fS->idx));
    }
  }
  if(lFlg)
  {
    AlcFree(lFlg);
  }
  if(vIdxTb)
  {
    AlcFree(vIdxTb);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellTestOutVTK
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 3D shell to the specified file using the
*		VTK polydata format.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell.
*		int *vIdxTb:		Vertex index table.
*		char *lFlg:		Loop flags.
*		FILE *fP:		File ptr, may be NULL if
*					no output required.
************************************************************************/
static WlzErrorNum WlzGMShellTestOutVTK(WlzGMShell *shell, int *vIdxTb,
					char *lFlg, FILE *fP)
{
  WlzGMLoopT	*tLT;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if(shell == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(shell->type != WLZ_GMELM_SHELL)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((tLT = shell->child)  != NULL)
  {
    do
    {
      if((tLT->idx >= 0) && (*(lFlg + tLT->loop->idx) == 0))
      {
        errNum = WlzGMLoopTTestOutputVTK(tLT, vIdxTb, fP);
	*(lFlg + tLT->loop->idx) = 1;
      }
      tLT = tLT->next;
    } while((errNum == WLZ_ERR_NONE) && (tLT->idx != shell->child->idx));
  }
  return(errNum);
}


/************************************************************************
* Function:	WlzGMLoopTTestOutputVTK
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 3D loopT to the specified file using the
*		VTK polydata format.
* Global refs:	-
* Parameters:	WlzGMLoopT *fLT:	Given loopT.
*		int *vIdxTb:		Vertex index table.
*		FILE *fP:		File ptr, may be NULL if
*					no output required.
************************************************************************/
static WlzErrorNum WlzGMLoopTTestOutputVTK(WlzGMLoopT *fLT, int *vIdxTb,
					   FILE *fP)
{
  WlzGMEdgeT	*tET;
  WlzGMVertex	*vBuf[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fLT == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(fLT->type != WLZ_GMELM_LOOP_T)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    /* Get the verticies around this loop and check that the loop is around
     * a triangular face. */
    if(fLT->edgeT == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else if((tET = fLT->edgeT)->type != WLZ_GMELM_EDGE_T)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      vBuf[0] = tET->vertexT->diskT->vertex;
      tET = tET->next;
      if(tET->idx == fLT->edgeT->idx)
      {
	errNum = WLZ_ERR_DOMAIN_DATA;
      }
      else
      {
	vBuf[1] = tET->vertexT->diskT->vertex;
	tET = tET->next;
	if(tET->idx == fLT->edgeT->idx)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else
	{
	  vBuf[2] = tET->vertexT->diskT->vertex;
	  tET = tET->next;
	  if(tET->idx != fLT->edgeT->idx)
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Output the face enclosed by this loop and mark it as having been
       * output. */
      (void )fprintf(fP, "3 %d %d %d\n",
		     *(vIdxTb + vBuf[0]->idx), *(vIdxTb + vBuf[1]->idx),
		     *(vIdxTb + vBuf[2]->idx));
    }
  }
  return(errNum);
}


/************************************************************************
* Function:	WlzGMModelTestOutPS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 2D model to the specified file as Postscript.
*		To make shell identification easier Postscript colors
*		can be cycled.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Given model.
*		FILE *fP:		File ptr, may be NULL if
*					no output required.
*		WlzDVertex2 offset:	Postscript offset.
*		WlzDVertex2 scale:	Postscript scale.
*		int nCol:		Number of colours to cycle
*					through.
************************************************************************/
WlzErrorNum   	WlzGMModelTestOutPS(WlzGMModel *model, FILE *fP,
				    WlzDVertex2 offset, WlzDVertex2 scale,
				    int nCol)
{
  int		colIdx  = 0;
  char		*lFlg = NULL,
		*eFlg = NULL;
  WlzGMShell	*shell0,
  		*shell1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	maxCol = 7;			  /* Colours to cycle though */
  const int	colR[] = {0,   255, 0,   255, 0,   200, 0},
		colG[] = {0,   0,   255, 200, 0,   0,   255},
		colB[] = {0,   0,   0,   0,   255, 255, 200};

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((model->type != WLZ_GMMOD_2I) && (model->type != WLZ_GMMOD_2D))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    /* Allocate shell, loop and edge flags. */
    if(((lFlg = (char *)AlcCalloc(model->res.loop.numIdx,
    				  sizeof(char))) == NULL) ||
       ((eFlg = (char *)AlcCalloc(model->res.edge.numIdx,
    				  sizeof(char))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((shell0 = model->child) != NULL)
    {
      if(nCol > maxCol)
      {
        nCol = maxCol;
      }
      shell1 = shell0;
      do
      {
	if(shell1->idx >= 0)
	{
	  /* Set colour. */
	  if(fP && (nCol > 1))
	  {
	    (void )fprintf(fP,
			   "%d %d %d setrgbcolor\n",
			   colR[colIdx], colG[colIdx], colB[colIdx]);
	    colIdx = ++colIdx % nCol;
	  }
	  errNum = WlzGMShellTestOutPS(shell1, lFlg, eFlg, fP, offset, scale);
	  shell1 = shell1->next;
	}
      }
      while((errNum == WLZ_ERR_NONE) && (shell1->idx != shell0->idx));
    }
  }
  if(lFlg)
  {
    AlcFree(lFlg);
  }
  if(eFlg)
  {
    AlcFree(eFlg);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellTestOutPS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 2D shell to the specified file as Postscript.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell.
*		char *lFlg:		Loop flags.
*		char *eFlg:		Edge flags.
*		FILE *fP:		Output file ptr, may be NULL if
*					no output required.
*		WlzDVertex2 offset:	Postscript offset.
*		WlzDVertex2 scale:	Postscript scale.
************************************************************************/
static WlzErrorNum WlzGMShellTestOutPS(WlzGMShell *shell,
				       char *lFlg, char *eFlg,
				       FILE *fP,
				       WlzDVertex2 offset, WlzDVertex2 scale)
{
  WlzGMLoopT	*loopT0,
  		*loopT1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(shell == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((loopT0 = shell->child) == NULL)
  {
    errNum = WLZ_ERR_NONE;	         /* TODO Is an empty shell an error? */
  }
  else
  {
    loopT1 = loopT0;
    do
    {
      if(*(lFlg + loopT1->loop->idx) == 0)
      {
	errNum = WlzGMLoopTTestOutPS(loopT1, eFlg, fP, offset, scale);
	*(lFlg + loopT1->loop->idx) = 1;
      }
      loopT1 = loopT1->next;
    }
    while((errNum == WLZ_ERR_NONE) && (loopT1->idx != loopT0->idx));
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMLoopTTestOutPS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 2D loop topology element to the specified
*		file as Postscript.
* Global refs:	-
* Parameters:	WlzGMLoopT *loopT:	Given loop topology element.
*		char *eFlg:		Edge flags.
*		FILE *fP:		Output file ptr, may be NULL if
*					no output required.
*		WlzDVertex2 offset:	Postscript offset.
*		WlzDVertex2 scale:	Postscript scale.
************************************************************************/
static WlzErrorNum WlzGMLoopTTestOutPS(WlzGMLoopT *loopT,
				       char *eFlg, FILE *fP,
			               WlzDVertex2 offset, WlzDVertex2 scale)
{
  WlzGMEdgeT	*edgeT0,
  		*edgeT1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((loopT == NULL) || ((edgeT0 = loopT->edgeT) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    edgeT1 = edgeT0;
    do
    {
      if((edgeT1->idx >= 0) && (*(eFlg + edgeT1->edge->idx) == 0))
      {
        errNum = WlzGMEdgeTTestOutPS(edgeT1, fP, offset, scale);
      }
      *(eFlg + edgeT1->edge->idx) = 1;
      edgeT1 = edgeT1->next;
    } 
    while((errNum == WLZ_ERR_NONE) && (edgeT1->idx != edgeT0->idx));
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMEdgeTTestOutPS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 2D edge topology element to the specified
*		file as Postscript.
* Global refs:	-
* Parameters:	WlzGMEdgeT *edgeT:	Given edge topology element.
*		FILE *fP:		Output file ptr, may be NULL if
*					no output required.
*		WlzDVertex2 offset:	Postscript offset.
*		WlzDVertex2 scale:	Postscript scale.
************************************************************************/
static WlzErrorNum WlzGMEdgeTTestOutPS(WlzGMEdgeT *edgeT, FILE *fP,
			               WlzDVertex2 offset, WlzDVertex2 scale)
{
  WlzDVertex2	pos0,
  		pos1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((edgeT == NULL) || (edgeT->opp == NULL) ||
     (edgeT->vertexT == NULL) || (edgeT->opp->vertexT == NULL) ||
     (edgeT->vertexT->diskT->vertex == NULL) ||
     (edgeT->opp->vertexT->diskT->vertex == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    if(((errNum = WlzGMVertexGetG2D(edgeT->vertexT->diskT->vertex,
    				    &pos0)) == WLZ_ERR_NONE) &&
       ((errNum = WlzGMVertexGetG2D(edgeT->next->vertexT->diskT->vertex,
    				    &pos1)) == WLZ_ERR_NONE) &&
       fP)
    {
      WlzGMTestOutLinePS(fP,  pos0, pos1, offset, scale);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMTestOutLinePS
* Returns:	void
* Purpose:	Outputs a 2D line segment, as defined by two endpoint
*		coordinates, to  the specifiedfile as Postscript.
* Global refs:	-
* Parameters:	FILE *fP:		Output file ptr.
************************************************************************/
static void	WlzGMTestOutLinePS(FILE *fP,
				   WlzDVertex2 pos0, WlzDVertex2 pos1,
				   WlzDVertex2 offset, WlzDVertex2 scale)
{
  pos0.vtX = offset.vtX + (scale.vtX * pos0.vtX);
  pos0.vtY = offset.vtY + (scale.vtY * pos0.vtY);
  pos1.vtX = offset.vtX + (scale.vtX * pos1.vtX);
  pos1.vtY = offset.vtY + (scale.vtY * pos1.vtY);
  (void )fprintf(fP,
  		 "%f %f moveto "
		 "%f %f lineto "
		 "stroke\n",
		 pos0.vtX, pos0.vtY,
		 pos1.vtX, pos1.vtY);
}

/* Geometry access, update and testing */

/************************************************************************
* Function:	WlzGMModelGetSGeomType
* Returns:	WlzGMElemType:		Shell's geometry type.
* Purpose:	Gets the shell's geometry type from the model's type.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model.
************************************************************************/
WlzGMElemType 	WlzGMModelGetSGeomType(WlzGMModel *model)
{
  WlzGMElemType	sGType;

  switch(model->type)
  {
    case WLZ_GMMOD_2I:
      sGType = WLZ_GMELM_SHELL_G2I;
      break;
    case WLZ_GMMOD_2D:
      sGType = WLZ_GMELM_SHELL_G2D;
      break;
    case WLZ_GMMOD_3I:
      sGType = WLZ_GMELM_SHELL_G3I;
      break;
    case WLZ_GMMOD_3D:
      sGType = WLZ_GMELM_SHELL_G3D;
      break;
  }
  return(sGType);
}

/************************************************************************
* Function:	WlzGMModelGetVGeomType
* Returns:	WlzGMElemType:		Shell's geometry type.
* Purpose:	Gets the verticies geometry type from the model's type.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model.
************************************************************************/
WlzGMElemType 	WlzGMModelGetVGeomType(WlzGMModel *model)
{
  WlzGMElemType	vGType;

  switch(model->type)
  {
    case WLZ_GMMOD_2I:
      vGType = WLZ_GMELM_VERTEX_G2I;
      break;
    case WLZ_GMMOD_2D:
      vGType = WLZ_GMELM_VERTEX_G2D;
      break;
    case WLZ_GMMOD_3I:
      vGType = WLZ_GMELM_VERTEX_G3I;
      break;
    case WLZ_GMMOD_3D:
      vGType = WLZ_GMELM_VERTEX_G3D;
      break;
  }
  return(vGType);
}

/************************************************************************
* Function:	WlzGMShellGetGBB3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Gets a shell's geometry using the given destination
*		pointer to a double precision bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDBox3 *bBox:		Given destination pointer to a
*					double precision bounding box.
************************************************************************/
WlzErrorNum	WlzGMShellGetGBB3D(WlzGMShell *shell, WlzDBox3 *bBox)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(bBox == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
        bBox->xMin = shell->geo.sg2I->bBox.xMin;
        bBox->yMin = shell->geo.sg2I->bBox.yMin;
	bBox->zMin = 0.0;
        bBox->xMax = shell->geo.sg2I->bBox.xMax;
        bBox->yMax = shell->geo.sg2I->bBox.yMax;
	bBox->zMax = 0.0;
	break;
      case WLZ_GMELM_SHELL_G2D:
        bBox->xMin = shell->geo.sg2D->bBox.xMin;
        bBox->yMin = shell->geo.sg2D->bBox.yMin;
	bBox->zMin = 0.0;
        bBox->xMax = shell->geo.sg2D->bBox.xMax;
        bBox->yMax = shell->geo.sg2D->bBox.yMax;
	bBox->zMax = 0.0;
	break;
      case WLZ_GMELM_SHELL_G3I:
        bBox->xMin = shell->geo.sg3I->bBox.xMin;
        bBox->yMin = shell->geo.sg3I->bBox.yMin;
	bBox->zMin = shell->geo.sg3I->bBox.zMin;
        bBox->xMax = shell->geo.sg3I->bBox.xMax;
        bBox->yMax = shell->geo.sg3I->bBox.yMax;
	bBox->zMax = shell->geo.sg3I->bBox.zMax;
	break;
      case WLZ_GMELM_SHELL_G3D:
        *bBox = shell->geo.sg3D->bBox;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellGetGBBV3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Gets the volume of the shell's geometry's bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDBox3 *bBox:		Given destination pointer for
*					the volume.
************************************************************************/
WlzErrorNum	WlzGMShellGetGBBV3D(WlzGMShell *shell, double *vol)
{
  WlzDBox3	bBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(vol == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((errNum = WlzGMShellGetGBB3D(shell, &bBox)) == WLZ_ERR_NONE)
  {
    *vol = (bBox.xMax - bBox.xMin) * (bBox.yMax - bBox.yMin) *
    	   (bBox.zMax - bBox.zMin);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellGetGBB2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Gets a shell's geometry using the given destination
*		pointer to a double precision bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDBox2 *bBox:		Given destination pointer to a
*					double precision bounding box.
************************************************************************/
WlzErrorNum	WlzGMShellGetGBB2D(WlzGMShell *shell, WlzDBox2 *bBox)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(bBox == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
        bBox->xMin = shell->geo.sg2I->bBox.xMin;
        bBox->yMin = shell->geo.sg2I->bBox.yMin;
        bBox->xMax = shell->geo.sg2I->bBox.xMax;
        bBox->yMax = shell->geo.sg2I->bBox.yMax;
	break;
      case WLZ_GMELM_SHELL_G2D:
        *bBox = shell->geo.sg2D->bBox;
	break;
      case WLZ_GMELM_SHELL_G3I:
        bBox->xMin = shell->geo.sg3I->bBox.xMin;
        bBox->yMin = shell->geo.sg3I->bBox.yMin;
        bBox->xMax = shell->geo.sg3I->bBox.xMax;
        bBox->yMax = shell->geo.sg3I->bBox.yMax;
	break;
      case WLZ_GMELM_SHELL_G3D:
	bBox->xMin = shell->geo.sg3D->bBox.xMin;
	bBox->yMin = shell->geo.sg3D->bBox.yMin;
	bBox->xMax = shell->geo.sg3D->bBox.xMax;
	bBox->yMax = shell->geo.sg3D->bBox.yMax;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellSetGBB3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a shell's geometry using the given double
*		precision bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDBox3 bBox:		Given double precision bounding
*					box.
************************************************************************/
WlzErrorNum	WlzGMShellSetGBB3D(WlzGMShell *shell, WlzDBox3 bBox)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
        shell->geo.sg2I->bBox.xMin = WLZ_NINT(bBox.xMin);
        shell->geo.sg2I->bBox.yMin = WLZ_NINT(bBox.yMin);
        shell->geo.sg2I->bBox.xMax = WLZ_NINT(bBox.xMax);
        shell->geo.sg2I->bBox.yMax = WLZ_NINT(bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G2D:
        shell->geo.sg2D->bBox.xMin = bBox.xMin;
	shell->geo.sg2D->bBox.yMin = bBox.yMin;
	shell->geo.sg2D->bBox.xMax = bBox.xMax;
	shell->geo.sg2D->bBox.yMax = bBox.yMax;
	break;
      case WLZ_GMELM_SHELL_G3I:
        shell->geo.sg3I->bBox.xMin = WLZ_NINT(bBox.xMin);
        shell->geo.sg3I->bBox.yMin = WLZ_NINT(bBox.yMin);
        shell->geo.sg3I->bBox.zMin = WLZ_NINT(bBox.zMin);
        shell->geo.sg3I->bBox.xMax = WLZ_NINT(bBox.xMax);
        shell->geo.sg3I->bBox.yMax = WLZ_NINT(bBox.yMax);
        shell->geo.sg3I->bBox.zMax = WLZ_NINT(bBox.zMax);
	break;
      case WLZ_GMELM_SHELL_G3D:
        shell->geo.sg3D->bBox = bBox;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellSetGBB2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a shell's geometry using the given double
*		precision bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDBox2 bBox:		Given double precision bounding
*					box.
************************************************************************/
WlzErrorNum	WlzGMShellSetGBB2D(WlzGMShell *shell, WlzDBox2 bBox)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
        shell->geo.sg2I->bBox.xMin = WLZ_NINT(bBox.xMin);
        shell->geo.sg2I->bBox.yMin = WLZ_NINT(bBox.yMin);
        shell->geo.sg2I->bBox.xMax = WLZ_NINT(bBox.xMax);
        shell->geo.sg2I->bBox.yMax = WLZ_NINT(bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G2D:
        shell->geo.sg2D->bBox = bBox;
	break;
      case WLZ_GMELM_SHELL_G3I:
        shell->geo.sg3I->bBox.xMin = WLZ_NINT(bBox.xMin);
        shell->geo.sg3I->bBox.yMin = WLZ_NINT(bBox.yMin);
        shell->geo.sg3I->bBox.zMin = 0;
        shell->geo.sg3I->bBox.xMax = WLZ_NINT(bBox.xMax);
        shell->geo.sg3I->bBox.yMax = WLZ_NINT(bBox.yMax);
        shell->geo.sg3I->bBox.zMax = 0;
	break;
      case WLZ_GMELM_SHELL_G3D:
        shell->geo.sg3D->bBox.xMin = bBox.xMin;
        shell->geo.sg3D->bBox.yMin = bBox.yMin;
        shell->geo.sg3D->bBox.zMin = 0.0;
        shell->geo.sg3D->bBox.xMax = bBox.xMax;
        shell->geo.sg3D->bBox.yMax = bBox.yMax;
        shell->geo.sg3D->bBox.zMax = 0.0;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellSetG3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a shell's geometry using the pair of double
*		precision points.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		int nPnt:		Number of points.
*		WlzDVertex3 *pos:	Array of points.
************************************************************************/
WlzErrorNum	WlzGMShellSetG3D(WlzGMShell *shell,
				 int nPnt, WlzDVertex3 *pos)
{
  WlzDBox3	bBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nPnt > 0) && pos)
  {
    bBox.xMin = bBox.xMax = pos->vtX;
    bBox.yMin = bBox.yMax = pos->vtY;
    bBox.zMin = bBox.zMax = pos->vtZ;
    while(--nPnt > 0)
    {
      ++pos;
      if(pos->vtX > bBox.xMax)
      {
	bBox.xMax = pos->vtX;
      }
      else if(pos->vtX < bBox.xMin)
      {
	bBox.xMin = pos->vtX;
      }
      if(pos->vtY > bBox.yMax)
      {
	bBox.yMax = pos->vtY;
      }
      else if(pos->vtY < bBox.xMin)
      {
	bBox.yMin = pos->vtY;
      }
      if(pos->vtZ > bBox.zMax)
      {
	bBox.zMax = pos->vtZ;
      }
      else if(pos->vtZ < bBox.zMin)
      {
	bBox.zMin = pos->vtZ;
      }
    }
    errNum = WlzGMShellSetGBB3D(shell, bBox);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellSetG2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a shell's geometry using the pair of double
*		precision points.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		int nPnt:		Number of points.
*		WlzDVertex2 *pos:	Array of points.
************************************************************************/
WlzErrorNum	WlzGMShellSetG2D(WlzGMShell *shell,
				 int nPnt, WlzDVertex2 *pos)
{
  WlzDBox2	bBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nPnt > 0) && pos)
  {
    bBox.xMin = bBox.xMax = pos->vtX;
    bBox.yMin = bBox.yMax = pos->vtY;
    while(--nPnt > 0)
    {
      ++pos;
      if(pos->vtX > bBox.xMax)
      {
	bBox.xMax = pos->vtX;
      }
      else if(pos->vtX < bBox.xMin)
      {
	bBox.xMin = pos->vtX;
      }
      if(pos->vtY > bBox.yMax)
      {
	bBox.yMax = pos->vtY;
      }
      else if(pos->vtY < bBox.xMin)
      {
	bBox.yMin = pos->vtY;
      }
    }
    errNum = WlzGMShellSetGBB2D(shell, bBox);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellUpdateG3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Updates a shell's geometry using the given double
*		precision position.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDVertex3 pos:	Given position.
************************************************************************/
WlzErrorNum	WlzGMShellUpdateG3D(WlzGMShell *shell, WlzDVertex3 pos)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	if((tI0 = pos.vtX) < shell->geo.sg2I->bBox.xMin)
	{
	  shell->geo.sg2I->bBox.xMin = tI0;
	}
	else if(tI0 > shell->geo.sg2I->bBox.xMax)
	{
	  shell->geo.sg2I->bBox.xMax = tI0;
	}
	if((tI0 = pos.vtY) < shell->geo.sg2I->bBox.yMin)
	{
	  shell->geo.sg2I->bBox.yMin = tI0;
	}
	else if(tI0 > shell->geo.sg2I->bBox.yMax)
	{
	  shell->geo.sg2I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G2D:
	if(pos.vtX < shell->geo.sg2D->bBox.xMin)
	{
	  shell->geo.sg2D->bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > shell->geo.sg2D->bBox.xMax)
	{
	  shell->geo.sg2D->bBox.xMax = pos.vtX;
	}
	if(pos.vtY < shell->geo.sg2D->bBox.yMin)
	{
	  shell->geo.sg2D->bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > shell->geo.sg2D->bBox.yMax)
	{
	  shell->geo.sg2D->bBox.yMax = pos.vtY;
	}
	break;
      case WLZ_GMELM_SHELL_G3I:
	if((tI0 = pos.vtX) < shell->geo.sg3I->bBox.xMin)
	{
	  shell->geo.sg3I->bBox.xMin = tI0;
	}
	else if(tI0 > shell->geo.sg3I->bBox.xMax)
	{
	  shell->geo.sg3I->bBox.xMax = tI0;
	}
	if((tI0 = pos.vtY) < shell->geo.sg3I->bBox.yMin)
	{
	  shell->geo.sg3I->bBox.yMin = tI0;
	}
	else if(tI0 > shell->geo.sg3I->bBox.yMax)
	{
	  shell->geo.sg3I->bBox.yMax = tI0;
	}
	if((tI0 = pos.vtZ) < shell->geo.sg3I->bBox.zMin)
	{
	  shell->geo.sg3I->bBox.zMin = tI0;
	}
	else if(tI0 > shell->geo.sg3I->bBox.zMax)
	{
	  shell->geo.sg3I->bBox.zMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G3D:
	if(pos.vtX < shell->geo.sg3D->bBox.xMin)
	{
	  shell->geo.sg3D->bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > shell->geo.sg3D->bBox.xMax)
	{
	  shell->geo.sg3D->bBox.xMax = pos.vtX;
	}
	if(pos.vtY < shell->geo.sg3D->bBox.yMin)
	{
	  shell->geo.sg3D->bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > shell->geo.sg3D->bBox.yMax)
	{
	  shell->geo.sg3D->bBox.yMax = pos.vtY;
	}
	if(pos.vtZ < shell->geo.sg3D->bBox.zMin)
	{
	  shell->geo.sg3D->bBox.zMin = pos.vtZ;
	}
	else if(pos.vtZ > shell->geo.sg3D->bBox.zMax)
	{
	  shell->geo.sg3D->bBox.zMax = pos.vtZ;
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
* Function:	WlzGMShellUpdateG2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Updates a shell's geometry using the given double
*		precision position.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDVertex2 pos:	Given position.
************************************************************************/
WlzErrorNum	WlzGMShellUpdateG2D(WlzGMShell *shell, WlzDVertex2 pos)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	if((tI0 = pos.vtX) < shell->geo.sg2I->bBox.xMin)
	{
	  shell->geo.sg2I->bBox.xMin = tI0;
	}
	else if(tI0 > shell->geo.sg2I->bBox.xMax)
	{
	  shell->geo.sg2I->bBox.xMax = tI0;
	}
	if((tI0 = pos.vtY) < shell->geo.sg2I->bBox.yMin)
	{
	  shell->geo.sg2I->bBox.yMin = tI0;
	}
	else if(tI0 > shell->geo.sg2I->bBox.yMax)
	{
	  shell->geo.sg2I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G2D:
	if(pos.vtX < shell->geo.sg2D->bBox.xMin)
	{
	  shell->geo.sg2D->bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > shell->geo.sg2D->bBox.xMax)
	{
	  shell->geo.sg2D->bBox.xMax = pos.vtX;
	}
	if(pos.vtY < shell->geo.sg2D->bBox.yMin)
	{
	  shell->geo.sg2D->bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > shell->geo.sg2D->bBox.yMax)
	{
	  shell->geo.sg2D->bBox.yMax = pos.vtY;
	}
	break;
      case WLZ_GMELM_SHELL_G3I:
	if((tI0 = pos.vtX) < shell->geo.sg3I->bBox.xMin)
	{
	  shell->geo.sg3I->bBox.xMin = tI0;
	}
	else if(tI0 > shell->geo.sg3I->bBox.xMax)
	{
	  shell->geo.sg3I->bBox.xMax = tI0;
	}
	if((tI0 = pos.vtY) < shell->geo.sg3I->bBox.yMin)
	{
	  shell->geo.sg3I->bBox.yMin = tI0;
	}
	else if(tI0 > shell->geo.sg3I->bBox.yMax)
	{
	  shell->geo.sg3I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G3D:
	if(pos.vtX < shell->geo.sg3D->bBox.xMin)
	{
	  shell->geo.sg3D->bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > shell->geo.sg3D->bBox.xMax)
	{
	  shell->geo.sg3D->bBox.xMax = pos.vtX;
	}
	if(pos.vtY < shell->geo.sg3D->bBox.yMin)
	{
	  shell->geo.sg3D->bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > shell->geo.sg3D->bBox.yMax)
	{
	  shell->geo.sg3D->bBox.yMax = pos.vtY;
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
* Function:	WlzGMShellUpdateGBB3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Updates a shell's geometry using the given double
*		precision position bounding box. The updated geometry
*		is the union of the existing geometry and the given
*		bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDBox3 bBox:		Given bounding box.
************************************************************************/
WlzErrorNum	WlzGMShellUpdateGBB3D(WlzGMShell *shell, WlzDBox3 bBox)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	if((tI0 = WLZ_NINT(bBox.xMin)) < shell->geo.sg2I->bBox.xMin)
	{
	  shell->geo.sg2I->bBox.xMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMin)) < shell->geo.sg2I->bBox.yMin)
	{
	  shell->geo.sg2I->bBox.yMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.xMax)) > shell->geo.sg2I->bBox.xMax)
	{
	  shell->geo.sg2I->bBox.xMax = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMax)) > shell->geo.sg2I->bBox.yMax)
	{
	  shell->geo.sg2I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G2D:
	if(bBox.xMin < shell->geo.sg2D->bBox.xMin)
	{
	  shell->geo.sg2D->bBox.xMin = bBox.xMin;
	}
	if(bBox.yMin < shell->geo.sg2D->bBox.yMin)
	{
	  shell->geo.sg2D->bBox.yMin = bBox.yMin;
	}
	if(bBox.xMax > shell->geo.sg2D->bBox.xMax)
	{
	  shell->geo.sg2D->bBox.xMax = bBox.xMax;
	}
	if(bBox.yMax > shell->geo.sg2D->bBox.yMax)
	{
	  shell->geo.sg2I->bBox.yMax = bBox.yMax;
	}
	break;
      case WLZ_GMELM_SHELL_G3I:
	if((tI0 = WLZ_NINT(bBox.xMin)) < shell->geo.sg3I->bBox.xMin)
	{
	  shell->geo.sg3I->bBox.xMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMin)) < shell->geo.sg3I->bBox.yMin)
	{
	  shell->geo.sg3I->bBox.yMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.zMin)) < shell->geo.sg3I->bBox.zMin)
	{
	  shell->geo.sg3I->bBox.zMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.xMax)) > shell->geo.sg3I->bBox.xMax)
	{
	  shell->geo.sg3I->bBox.xMax = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMax)) > shell->geo.sg3I->bBox.yMax)
	{
	  shell->geo.sg3I->bBox.yMax = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.zMax)) > shell->geo.sg3I->bBox.zMax)
	{
	  shell->geo.sg3I->bBox.zMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G3D:
	if(bBox.xMin < shell->geo.sg3D->bBox.xMin)
	{
	  shell->geo.sg3D->bBox.xMin = bBox.xMin;
	}
	if(bBox.yMin < shell->geo.sg3D->bBox.yMin)
	{
	  shell->geo.sg3D->bBox.yMin = bBox.yMin;
	}
	if(bBox.zMin < shell->geo.sg3D->bBox.zMin)
	{
	  shell->geo.sg3D->bBox.zMin = bBox.zMin;
	}
	if(bBox.xMax > shell->geo.sg3D->bBox.xMax)
	{
	  shell->geo.sg3D->bBox.xMax = bBox.xMax;
	}
	if(bBox.yMax > shell->geo.sg3D->bBox.yMax)
	{
	  shell->geo.sg3D->bBox.yMax = bBox.yMax;
	}
	if(bBox.zMax > shell->geo.sg3D->bBox.zMax)
	{
	  shell->geo.sg3D->bBox.zMax = bBox.zMax;
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
* Function:	WlzGMShellUpdateGBB2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Updates a shell's geometry using the given double
*		precision position bounding box. The updated geometry
*		is the union of the existing geometry and the given
*		bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDBox2 bBox:		Given bounding box.
************************************************************************/
WlzErrorNum	WlzGMShellUpdateGBB2D(WlzGMShell *shell, WlzDBox2 bBox)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	if((tI0 = WLZ_NINT(bBox.xMin)) < shell->geo.sg2I->bBox.xMin)
	{
	  shell->geo.sg2I->bBox.xMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMin)) < shell->geo.sg2I->bBox.yMin)
	{
	  shell->geo.sg2I->bBox.yMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.xMax)) > shell->geo.sg2I->bBox.xMax)
	{
	  shell->geo.sg2I->bBox.xMax = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMax)) > shell->geo.sg2I->bBox.yMax)
	{
	  shell->geo.sg2I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G2D:
	if(bBox.xMin < shell->geo.sg2D->bBox.xMin)
	{
	  shell->geo.sg2D->bBox.xMin = bBox.xMin;
	}
	if(bBox.yMin < shell->geo.sg2D->bBox.yMin)
	{
	  shell->geo.sg2D->bBox.yMin = bBox.yMin;
	}
	if(bBox.xMax > shell->geo.sg2D->bBox.xMax)
	{
	  shell->geo.sg2D->bBox.xMax = bBox.xMax;
	}
	if(bBox.yMax > shell->geo.sg2D->bBox.yMax)
	{
	  shell->geo.sg2I->bBox.yMax = bBox.yMax;
	}
	break;
      case WLZ_GMELM_SHELL_G3I:
	if((tI0 = WLZ_NINT(bBox.xMin)) < shell->geo.sg3I->bBox.xMin)
	{
	  shell->geo.sg3I->bBox.xMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMin)) < shell->geo.sg3I->bBox.yMin)
	{
	  shell->geo.sg3I->bBox.yMin = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.xMax)) > shell->geo.sg3I->bBox.xMax)
	{
	  shell->geo.sg3I->bBox.xMax = tI0;
	}
	if((tI0 = WLZ_NINT(bBox.yMax)) > shell->geo.sg3I->bBox.yMax)
	{
	  shell->geo.sg3I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_SHELL_G3D:
	if(bBox.xMin < shell->geo.sg3D->bBox.xMin)
	{
	  shell->geo.sg3D->bBox.xMin = bBox.xMin;
	}
	if(bBox.yMin < shell->geo.sg3D->bBox.yMin)
	{
	  shell->geo.sg3D->bBox.yMin = bBox.yMin;
	}
	if(bBox.xMax > shell->geo.sg3D->bBox.xMax)
	{
	  shell->geo.sg3D->bBox.xMax = bBox.xMax;
	}
	if(bBox.yMax > shell->geo.sg3D->bBox.yMax)
	{
	  shell->geo.sg3D->bBox.yMax = bBox.yMax;
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
* Function:	WlzGMShellComputeGBB
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Recomputes the shell's geometry by walking through
*		it's children.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
************************************************************************/
WlzErrorNum	WlzGMShellComputeGBB(WlzGMShell *shell)
{
  int		tI0,
  		fETIdx,
  		fLTIdx;
  double	tD0;
  WlzGMEdgeT	*tET;
  WlzGMLoopT	*tLT;
  WlzGMVertexGU	vGU;
  WlzGMShellGU	sGU;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || ((sGU = shell->geo).core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(shell->child)
  {
    switch(sGU.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	vGU = shell->child->edgeT->vertexT->diskT->vertex->geo;
	sGU.sg2I->bBox.xMin = sGU.sg2I->bBox.xMax = vGU.vg2I->vtx.vtX;
	sGU.sg2I->bBox.yMin = sGU.sg2I->bBox.yMax = vGU.vg2I->vtx.vtY;
	/* For each loopT */
	tLT = shell->child;
	fLTIdx = tLT->idx;
	do
	{
	  /* For each edgeT */
	  tET = tLT->edgeT;
	  fETIdx = tET->idx;
	  do
	  {
	    /* Update the new bounding box. */
	    vGU = tET->vertexT->diskT->vertex->geo;
	    if((tI0 = vGU.vg2I->vtx.vtX) < sGU.sg2I->bBox.xMin)
	    {
	      sGU.sg2I->bBox.xMin = tI0;
	    }
	    else if(tI0 > sGU.sg2I->bBox.xMax)
	    {
	      sGU.sg2I->bBox.xMax = tI0;
	    }
	    if((tI0 = vGU.vg2I->vtx.vtY) < sGU.sg2I->bBox.yMin)
	    {
	      sGU.sg2I->bBox.yMin = tI0;
	    }
	    else if(tI0 > sGU.sg2I->bBox.yMax)
	    {
	      sGU.sg2I->bBox.yMax = tI0;
	    }
	    /* Next edge topology element. */
	    tET = tET->next;
	  } while(tET->idx != fETIdx);
	  tLT = tLT->next;
	} while(tLT->idx != fLTIdx);
	break;
      case WLZ_GMELM_SHELL_G2D:
	vGU = shell->child->edgeT->vertexT->diskT->vertex->geo;
	sGU.sg2D->bBox.xMin = sGU.sg2D->bBox.xMax = vGU.vg2D->vtx.vtX;
	sGU.sg2D->bBox.yMin = sGU.sg2D->bBox.yMax = vGU.vg2D->vtx.vtY;
	/* For each loopT */
	tLT = shell->child;
	fLTIdx = tLT->idx;
	do
	{
	  /* For each edgeT */
	  tET = tLT->edgeT;
	  fETIdx = tET->idx;
	  do
	  {
	    /* Update the new bounding box. */
	    vGU = tET->vertexT->diskT->vertex->geo;
	    if((tD0 = vGU.vg2D->vtx.vtX) < sGU.sg2D->bBox.xMin)
	    {
	      sGU.sg2D->bBox.xMin = tD0;
	    }
	    else if(tD0 > sGU.sg2D->bBox.xMax)
	    {
	      sGU.sg2D->bBox.xMax = tD0;
	    }
	    if((tD0 = vGU.vg2D->vtx.vtY) < sGU.sg2D->bBox.yMin)
	    {
	      sGU.sg2D->bBox.yMin = tD0;
	    }
	    else if(tD0 > sGU.sg2D->bBox.yMax)
	    {
	      sGU.sg2D->bBox.yMax = tD0;
	    }
	    /* Next edge topology element. */
	    tET = tET->next;
	  } while(tET->idx != fETIdx);
	  tLT = tLT->next;
	} while(tLT->idx != fLTIdx);
	break;
      case WLZ_GMELM_SHELL_G3I:
	vGU = shell->child->edgeT->vertexT->diskT->vertex->geo;
	sGU.sg3I->bBox.xMin = sGU.sg3I->bBox.xMax = vGU.vg3I->vtx.vtX;
	sGU.sg3I->bBox.yMin = sGU.sg3I->bBox.yMax = vGU.vg3I->vtx.vtY;
	sGU.sg3I->bBox.zMin = sGU.sg3I->bBox.zMax = vGU.vg3I->vtx.vtZ;
	/* For each loopT */
	tLT = shell->child;
	fLTIdx = tLT->idx;
	do
	{
	  /* For each edgeT */
	  tET = tLT->edgeT;
	  fETIdx = tET->idx;
	  do
	  {
	    /* Update the new bounding box. */
	    vGU = tET->vertexT->diskT->vertex->geo;
	    if((tI0 = vGU.vg3I->vtx.vtX) < sGU.sg3I->bBox.xMin)
	    {
	      sGU.sg3I->bBox.xMin = tI0;
	    }
	    else if(tI0 > sGU.sg3I->bBox.xMax)
	    {
	      sGU.sg3I->bBox.xMax = tI0;
	    }
	    if((tI0 = vGU.vg3I->vtx.vtY) < sGU.sg3I->bBox.yMin)
	    {
	      sGU.sg3I->bBox.yMin = tI0;
	    }
	    else if(tI0 > sGU.sg3I->bBox.yMax)
	    {
	      sGU.sg3I->bBox.yMax = tI0;
	    }
	    if((tI0 = vGU.vg3I->vtx.vtZ) < sGU.sg3I->bBox.zMin)
	    {
	      sGU.sg3I->bBox.zMin = tI0;
	    }
	    else if(tI0 > sGU.sg3I->bBox.zMax)
	    {
	      sGU.sg3I->bBox.zMax = tI0;
	    }
	    /* Next edge topology element. */
	    tET = tET->next;
	  } while(tET->idx != fETIdx);
	  tLT = tLT->next;
	} while(tLT->idx != fLTIdx);
	break;
      case WLZ_GMELM_SHELL_G3D:
	vGU = shell->child->edgeT->vertexT->diskT->vertex->geo;
	sGU.sg3D->bBox.xMin = sGU.sg3D->bBox.xMax = vGU.vg3D->vtx.vtX;
	sGU.sg3D->bBox.yMin = sGU.sg3D->bBox.yMax = vGU.vg3D->vtx.vtY;
	sGU.sg3D->bBox.zMin = sGU.sg3D->bBox.zMax = vGU.vg3D->vtx.vtZ;
	/* For each loopT */
	tLT = shell->child;
	fLTIdx = tLT->idx;
	do
	{
	  /* For each edgeT */
	  tET = tLT->edgeT;
	  fETIdx = tET->idx;
	  do
	  {
	    /* Update the new bounding box. */
	    vGU = tET->vertexT->diskT->vertex->geo;
	    if((tD0 = vGU.vg3D->vtx.vtX) < sGU.sg3D->bBox.xMin)
	    {
	      sGU.sg3D->bBox.xMin = tD0;
	    }
	    else if(tD0 > sGU.sg3D->bBox.xMax)
	    {
	      sGU.sg3D->bBox.xMax = tD0;
	    }
	    if((tD0 = vGU.vg3D->vtx.vtY) < sGU.sg3D->bBox.yMin)
	    {
	      sGU.sg3D->bBox.yMin = tD0;
	    }
	    else if(tD0 > sGU.sg3D->bBox.yMax)
	    {
	      sGU.sg3D->bBox.yMax = tD0;
	    }
	    if((tD0 = vGU.vg3D->vtx.vtZ) < sGU.sg3D->bBox.zMin)
	    {
	      sGU.sg3D->bBox.zMin = tD0;
	    }
	    else if(tD0 > sGU.sg3D->bBox.zMax)
	    {
	      sGU.sg3D->bBox.zMax = tD0;
	    }
	    /* Next edge topology element. */
	    tET = tET->next;
	  } while(tET->idx != fETIdx);
	  tLT = tLT->next;
	} while(tLT->idx != fLTIdx);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellMergeG
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Merges the geometries of the two shells so that the
*		first shell's geometry is set to the union of the
*		two shell's bounding boxes.
* Global refs:	-
* Parameters:	WlzGMShell *shell0:	First shell with geometry to
*					be both used and set.
*		WlzGMShell *shell1:	Second shell with geometry to
*					be used.
************************************************************************/
static WlzErrorNum WlzGMShellMergeG(WlzGMShell *shell0, WlzGMShell *shell1)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell0 == NULL) || (shell0->geo.core == NULL) ||
     (shell1 == NULL) || (shell1->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(shell0->geo.core->type != shell1->geo.core->type)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  {
    switch(shell0->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	if(shell1->geo.sg2I->bBox.xMin < shell0->geo.sg2I->bBox.xMin)
	{
	  shell0->geo.sg2I->bBox.xMin = shell1->geo.sg2I->bBox.xMin;
	}
	if(shell1->geo.sg2I->bBox.yMin < shell0->geo.sg2I->bBox.yMin)
	{
	  shell0->geo.sg2I->bBox.yMin = shell1->geo.sg2I->bBox.yMin;
	}
	if(shell1->geo.sg2I->bBox.xMax > shell0->geo.sg2I->bBox.xMax)
	{
	  shell0->geo.sg2I->bBox.xMax = shell1->geo.sg2I->bBox.xMax;
	}
	if(shell1->geo.sg2I->bBox.yMax > shell0->geo.sg2I->bBox.yMax)
	{
	  shell0->geo.sg2I->bBox.yMax = shell1->geo.sg2I->bBox.yMax;
	}
	break;
      case WLZ_GMELM_SHELL_G2D:
	if(shell1->geo.sg2D->bBox.xMin < shell0->geo.sg2D->bBox.xMin)
	{
	  shell0->geo.sg2D->bBox.xMin = shell1->geo.sg2D->bBox.xMin;
	}
	if(shell1->geo.sg2D->bBox.yMin < shell0->geo.sg2D->bBox.yMin)
	{
	  shell0->geo.sg2D->bBox.yMin = shell1->geo.sg2D->bBox.yMin;
	}
	if(shell1->geo.sg2D->bBox.xMax > shell0->geo.sg2D->bBox.xMax)
	{
	  shell0->geo.sg2D->bBox.xMax = shell1->geo.sg2D->bBox.xMax;
	}
	if(shell1->geo.sg2D->bBox.yMax > shell0->geo.sg2D->bBox.yMax)
	{
	  shell0->geo.sg2D->bBox.yMax = shell1->geo.sg2D->bBox.yMax;
	}
	break;
      case WLZ_GMELM_SHELL_G3I:
	if(shell1->geo.sg3I->bBox.xMin < shell0->geo.sg3I->bBox.xMin)
	{
	  shell0->geo.sg3I->bBox.xMin = shell1->geo.sg3I->bBox.xMin;
	}
	if(shell1->geo.sg3I->bBox.yMin < shell0->geo.sg3I->bBox.yMin)
	{
	  shell0->geo.sg3I->bBox.yMin = shell1->geo.sg3I->bBox.yMin;
	}
	if(shell1->geo.sg3I->bBox.zMin < shell0->geo.sg3I->bBox.zMin)
	{
	  shell0->geo.sg3I->bBox.zMin = shell1->geo.sg3I->bBox.zMin;
	}
	if(shell1->geo.sg3I->bBox.xMax > shell0->geo.sg3I->bBox.xMax)
	{
	  shell0->geo.sg3I->bBox.xMax = shell1->geo.sg3I->bBox.xMax;
	}
	if(shell1->geo.sg3I->bBox.yMax > shell0->geo.sg3I->bBox.yMax)
	{
	  shell0->geo.sg3I->bBox.yMax = shell1->geo.sg3I->bBox.yMax;
	}
	if(shell1->geo.sg3I->bBox.zMax > shell0->geo.sg3I->bBox.zMax)
	{
	  shell0->geo.sg3I->bBox.zMax = shell1->geo.sg3I->bBox.zMax;
	}
	break;
      case WLZ_GMELM_SHELL_G3D:
	if(shell1->geo.sg3D->bBox.xMin < shell0->geo.sg3D->bBox.xMin)
	{
	  shell0->geo.sg3D->bBox.xMin = shell1->geo.sg3D->bBox.xMin;
	}
	if(shell1->geo.sg3D->bBox.yMin < shell0->geo.sg3D->bBox.yMin)
	{
	  shell0->geo.sg3D->bBox.yMin = shell1->geo.sg3D->bBox.yMin;
	}
	if(shell1->geo.sg3D->bBox.zMin < shell0->geo.sg3D->bBox.zMin)
	{
	  shell0->geo.sg3D->bBox.zMin = shell1->geo.sg3D->bBox.zMin;
	}
	if(shell1->geo.sg3D->bBox.xMax > shell0->geo.sg3D->bBox.xMax)
	{
	  shell0->geo.sg3D->bBox.xMax = shell1->geo.sg3D->bBox.xMax;
	}
	if(shell1->geo.sg3D->bBox.yMax > shell0->geo.sg3D->bBox.yMax)
	{
	  shell0->geo.sg3D->bBox.yMax = shell1->geo.sg3D->bBox.yMax;
	}
	if(shell1->geo.sg3D->bBox.zMax > shell0->geo.sg3D->bBox.zMax)
	{
	  shell0->geo.sg3D->bBox.zMax = shell1->geo.sg3D->bBox.zMax;
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
* Function:	WlzGMShellGInBB3D
* Returns:	int:			Non zero if point inside shell's
*					bounding box.
* Purpose:	Checks to see if the given double precision position
*		is within the shell's bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDVertex3 pos:	Given double precision position.
************************************************************************/
int		WlzGMShellGInBB3D(WlzGMShell *shell, WlzDVertex3 pos)
{
  int		inside = 0;

  if((shell != NULL) && (shell->geo.core != NULL))
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	inside = (pos.vtX >= shell->geo.sg2I->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg2I->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg2I->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg2I->bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G2D:
	inside = (pos.vtX >= shell->geo.sg2D->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg2D->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg2D->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg2D->bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G3I:
	inside = (pos.vtX >= shell->geo.sg3I->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg3I->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg3I->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg3I->bBox.yMax) &&
	         (pos.vtZ >= shell->geo.sg3I->bBox.zMin) &&
		 (pos.vtZ <= shell->geo.sg3I->bBox.zMax);
	break;
      case WLZ_GMELM_SHELL_G3D:
	inside = (pos.vtX >= shell->geo.sg3D->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg3D->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg3D->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg3D->bBox.yMax) &&
	         (pos.vtZ >= shell->geo.sg3D->bBox.zMin) &&
		 (pos.vtZ <= shell->geo.sg3D->bBox.zMax);
	break;
      default:
	break;
    }
  }
  return(inside);
}

/************************************************************************
* Function:	WlzGMShellGInBB2D
* Returns:	int:			Non zero if point inside shell's
*					bounding box.
* Purpose:	Checks to see if the given double precision position
*		is within the shell's bounding box.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDVertex2 pos:	Given double precision position.
************************************************************************/
int		WlzGMShellGInBB2D(WlzGMShell *shell, WlzDVertex2 pos)
{
  int		inside = 0;

  if((shell != NULL) && (shell->geo.core != NULL))
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G2I:
	inside = (pos.vtX >= shell->geo.sg2I->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg2I->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg2I->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg2I->bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G2D:
	inside = (pos.vtX >= shell->geo.sg2D->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg2D->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg2D->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg2D->bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G3I:
	inside = (pos.vtX >= shell->geo.sg3I->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg3I->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg3I->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg3I->bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G3D:
	inside = (pos.vtX >= shell->geo.sg3D->bBox.xMin) &&
		 (pos.vtX <= shell->geo.sg3D->bBox.xMax) &&
	         (pos.vtY >= shell->geo.sg3D->bBox.yMin) &&
		 (pos.vtY <= shell->geo.sg3D->bBox.yMax);
	break;
      default:
	break;
    }
  }
  return(inside);
}

/************************************************************************
* Function:	WlzGMVertexSetG3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a veticies geometry using the given double
*		precision position.
* Global refs:	-
* Parameters:	WlzGMVertex *vertex:	Given vertex with geometry to
*					be set.
*		WlzDVertex3 pos:	Given position.
************************************************************************/
WlzErrorNum	WlzGMVertexSetG3D(WlzGMVertex *vertex, WlzDVertex3 pos)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((vertex == NULL) || (vertex->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
        vertex->geo.vg2I->vtx.vtX = WLZ_NINT(pos.vtX);
        vertex->geo.vg2I->vtx.vtY = WLZ_NINT(pos.vtY);
	break;
      case WLZ_GMELM_VERTEX_G2D:
        vertex->geo.vg2D->vtx.vtX = pos.vtX;
	vertex->geo.vg2D->vtx.vtY = pos.vtY;
	break;
      case WLZ_GMELM_VERTEX_G3I:
        vertex->geo.vg3I->vtx.vtX = WLZ_NINT(pos.vtX);
	vertex->geo.vg3I->vtx.vtY = WLZ_NINT(pos.vtY);
	vertex->geo.vg3I->vtx.vtZ = WLZ_NINT(pos.vtZ);
	break;
      case WLZ_GMELM_VERTEX_G3D:
        vertex->geo.vg3D->vtx = pos;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMVertexSetG2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a veticies geometry using the given double
*		precision position.
* Global refs:	-
* Parameters:	WlzGMVertex *vertex:	Given vertex with geometry to
*					be set.
*		WlzDVertex2 pos:	Given position.
************************************************************************/
WlzErrorNum	WlzGMVertexSetG2D(WlzGMVertex *vertex, WlzDVertex2 pos)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((vertex == NULL) || (vertex->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
        vertex->geo.vg2I->vtx.vtX = WLZ_NINT(pos.vtX);
        vertex->geo.vg2I->vtx.vtY = WLZ_NINT(pos.vtY);
	break;
      case WLZ_GMELM_VERTEX_G2D:
        vertex->geo.vg2D->vtx = pos;
	break;
      case WLZ_GMELM_VERTEX_G3I:
        vertex->geo.vg3I->vtx.vtX = WLZ_NINT(pos.vtX);
	vertex->geo.vg3I->vtx.vtY = WLZ_NINT(pos.vtY);
	vertex->geo.vg3I->vtx.vtZ = 0;
	break;
      case WLZ_GMELM_VERTEX_G3D:
        vertex->geo.vg3D->vtx.vtX = pos.vtX;
	vertex->geo.vg3D->vtx.vtY = pos.vtY;
	vertex->geo.vg3D->vtx.vtZ = 0.0;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMVertexGetG3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Gets a veticies geometry into the given double
*		precision position destination pointer.
* Global refs:	-
* Parameters:	WlzGMVertex *vertex:	Given vertex with geometry to
*					be set.
*		WlzDVertex3 *dstPos:	Given position destination
*					pointer, may NOT be NULL.
************************************************************************/
WlzErrorNum	WlzGMVertexGetG3D(WlzGMVertex *vertex, WlzDVertex3 *dstPos)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((vertex == NULL) || (vertex->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
	dstPos->vtX = vertex->geo.vg2I->vtx.vtX;
	dstPos->vtY = vertex->geo.vg2I->vtx.vtY;
	dstPos->vtZ = 0;
	break;
      case WLZ_GMELM_VERTEX_G2D:
	dstPos->vtX = vertex->geo.vg2D->vtx.vtX;
	dstPos->vtY = vertex->geo.vg2D->vtx.vtY;
	dstPos->vtZ = 0.0;
	break;
      case WLZ_GMELM_VERTEX_G3I:
	dstPos->vtX = vertex->geo.vg3I->vtx.vtX;
	dstPos->vtY = vertex->geo.vg3I->vtx.vtY;
	dstPos->vtZ = vertex->geo.vg3I->vtx.vtZ;
	break;
      case WLZ_GMELM_VERTEX_G3D:
	*dstPos = vertex->geo.vg3D->vtx;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMVertexGetG2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Gets a veticies geometry into the given double
*		precision position destination pointer.
* Global refs:	-
* Parameters:	WlzGMVertex *vertex:	Given vertex with geometry to
*					be set.
*		WlzDVertex2 *dstPos:	Given position destination
*					pointer, may NOT be NULL.
************************************************************************/
WlzErrorNum	WlzGMVertexGetG2D(WlzGMVertex *vertex, WlzDVertex2 *dstPos)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((vertex == NULL) || (vertex->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
	dstPos->vtX = vertex->geo.vg2I->vtx.vtX;
	dstPos->vtY = vertex->geo.vg2I->vtx.vtY;
	break;
      case WLZ_GMELM_VERTEX_G2D:
	*dstPos = vertex->geo.vg2D->vtx;
	break;
      case WLZ_GMELM_VERTEX_G3I:
	dstPos->vtX = vertex->geo.vg3I->vtx.vtX;
	dstPos->vtY = vertex->geo.vg3I->vtx.vtY;
	break;
      case WLZ_GMELM_VERTEX_G3D:
	dstPos->vtX = vertex->geo.vg3D->vtx.vtX;
	dstPos->vtY = vertex->geo.vg3D->vtx.vtY;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMVertexCmp3D
* Returns:	WlzDVertex3:			{vertex x|y|z - given
						  x|y|z}.
* Purpose:	Compares column coordinates of the given vertex and
*		3D double precision position.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Shell with resources.
*		WlzGMVertex *vertex:	Vertex to free.
************************************************************************/
WlzDVertex3	WlzGMVertexCmp3D(WlzGMVertex *vertex, WlzDVertex3 pos)
{
  WlzDVertex3	cmp;


  cmp.vtZ = cmp.vtY = cmp.vtX = 0.0;
  if(vertex && vertex->geo.core)
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
        cmp.vtX = vertex->geo.vg2I->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg2I->vtx.vtY - pos.vtY;
        break;
      case WLZ_GMELM_VERTEX_G2D:
        cmp.vtX = vertex->geo.vg2D->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg2D->vtx.vtY - pos.vtY;
        break;
      case WLZ_GMELM_VERTEX_G3I:
        cmp.vtX = vertex->geo.vg3I->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg3I->vtx.vtY - pos.vtY;
        cmp.vtZ = vertex->geo.vg3I->vtx.vtZ - pos.vtZ;
        break;
      case WLZ_GMELM_VERTEX_G3D:
        cmp.vtX = vertex->geo.vg3D->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg3D->vtx.vtY - pos.vtY;
        cmp.vtZ = vertex->geo.vg3D->vtx.vtZ - pos.vtZ;
        break;
      default:
        break;
    }
  }
  return(cmp);
}

/************************************************************************
* Function:	WlzGMVertexCmp2D
* Returns:	WlzDVertex2:			{vertex x|y - given
						  x|y}.
* Purpose:	Compares column coordinates of the given vertex and
*		2D double precision position.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Shell with resources.
*		WlzGMVertex *vertex:	Vertex to free.
************************************************************************/
WlzDVertex2	WlzGMVertexCmp2D(WlzGMVertex *vertex, WlzDVertex2 pos)
{
  WlzDVertex2	cmp;


  cmp.vtY = cmp.vtX = 0.0;
  if(vertex && vertex->geo.core)
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
        cmp.vtX = vertex->geo.vg2I->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg2I->vtx.vtY - pos.vtY;
        break;
      case WLZ_GMELM_VERTEX_G2D:
        cmp.vtX = vertex->geo.vg2D->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg2D->vtx.vtY - pos.vtY;
        break;
      case WLZ_GMELM_VERTEX_G3I:
        cmp.vtX = vertex->geo.vg3I->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg3I->vtx.vtY - pos.vtY;
        break;
      case WLZ_GMELM_VERTEX_G3D:
        cmp.vtX = vertex->geo.vg3D->vtx.vtX - pos.vtX;
        cmp.vtY = vertex->geo.vg3D->vtx.vtY - pos.vtY;
        break;
      default:
        break;
    }
  }
  return(cmp);
}

/************************************************************************
* Function:	WlzGMVertexCmpSign3D
* Returns:	int:			Signed value for sorting.
						  x|y|z}.
* Purpose:	Compares column coordinates of the given vertex and
*		3D double precision position to find a signed value
*		for sorting.
* Global refs:	-
* Parameters:	WlzGMVertex *vertex:	Vertex to free.
*		WlzDVertex3 pos:	3D double precision position.
************************************************************************/
int		WlzGMVertexCmpSign3D(WlzGMVertex *vertex, WlzDVertex3 pos)
{
  int		cmp;
  WlzDVertex3	cmp3D;

  cmp3D = WlzGMVertexCmp3D(vertex, pos);
  if(cmp3D.vtZ < -WLZ_GM_TOLERANCE)
  {
    cmp = -1;
  }
  else if(cmp3D.vtZ > WLZ_GM_TOLERANCE)
  {
    cmp = 1;
  }
  else
  {
    if(cmp3D.vtY < -WLZ_GM_TOLERANCE)
    {
      cmp = -1;
    }
    else if(cmp3D.vtY > WLZ_GM_TOLERANCE)
    {
      cmp = 1;
    }
    else
    {
      if(cmp3D.vtX < -WLZ_GM_TOLERANCE)
      {
	cmp = -1;
      }
      else if(cmp3D.vtX > WLZ_GM_TOLERANCE)
      {
	cmp = 1;
      }
      else
      {
	cmp = 0;
      }
    }
  }
  return(cmp);
}

/************************************************************************
* Function:	WlzGMVertexCmpSign2D
* Returns:	int:			Signed value for sorting.
						  x|y|z}.
* Purpose:	Compares column coordinates of the given vertex and
*		2D double precision position to find a signed value
*		for sorting.
* Global refs:	-
* Parameters:	WlzGMVertex *vertex:	Vertex to free.
*		WlzDVertex2 pos:	2D double precision position.
************************************************************************/
int		WlzGMVertexCmpSign2D(WlzGMVertex *vertex, WlzDVertex2 pos)
{
  int		cmp;
  WlzDVertex2	cmp2D;

  cmp2D = WlzGMVertexCmp2D(vertex, pos);
  if(cmp2D.vtY < -WLZ_GM_TOLERANCE)
  {
    cmp = -1;
  }
  else if(cmp2D.vtY > WLZ_GM_TOLERANCE)
  {
    cmp = 1;
  }
  else
  {
    if(cmp2D.vtX < -WLZ_GM_TOLERANCE)
    {
      cmp = -1;
    }
    else if(cmp2D.vtX > WLZ_GM_TOLERANCE)
    {
      cmp = 1;
    }
    else
    {
      cmp = 0;
    }
  }
  return(cmp);
}

/************************************************************************
* Function:	WlzGMVertexDistSq2D
* Returns:	double:			Square of distance, -1.0 on
*					error.
* Purpose:	Calculates the square of the Euclidean distance
*		between the given vertex and the given 3D double
*		precision position.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Shell with resources.
*		WlzGMVertex *vertex:	Vertex to free.
************************************************************************/
double		WlzGMVertexDistSq3D(WlzGMVertex *vertex, WlzDVertex3 pos)
{
  double	tD0,
		tD1,
		tD2,
		dstSq = -1.0;

  if(vertex && vertex->geo.core)
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
        tD0 = pos.vtX - vertex->geo.vg2I->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg2I->vtx.vtY;
	dstSq = (tD0 * tD0) + (tD1 * tD1);
        break;
      case WLZ_GMELM_VERTEX_G2D:
        tD0 = pos.vtX - vertex->geo.vg2D->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg2D->vtx.vtY;
	dstSq = (tD0 * tD0) + (tD1 * tD1);
        break;
      case WLZ_GMELM_VERTEX_G3I:
        tD0 = pos.vtX - vertex->geo.vg3I->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg3I->vtx.vtY;
	tD2 = pos.vtZ - vertex->geo.vg3I->vtx.vtZ;
	dstSq = (tD0 * tD0) + (tD1 * tD1) + (tD2 * tD2);
        break;
      case WLZ_GMELM_VERTEX_G3D:
        tD0 = pos.vtX - vertex->geo.vg3D->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg3D->vtx.vtY;
	tD2 = pos.vtZ - vertex->geo.vg3D->vtx.vtZ;
	dstSq = (tD0 * tD0) + (tD1 * tD1) + (tD2 * tD2);
        break;
      default:
        break;
    }
  }
  return(dstSq);
}

/************************************************************************
* Function:	WlzGMVertexNormal3D
* Returns:	WlzDVertex3		Value of normal.
* Purpose:	Computes the value of the normal at the given vertex
*		which lies within the given model.
*		This function requires a buffer in which to store the
*		verticies found on the loops surrounding the given
*		vertex. For efficiency this can/should be reused
*		between calls of this function.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The given model.
*		WlzGMVertex *gV:	Given vertex in the model
*		int *sVBufSz:		Ptr to the number WlzGMVertex's
*					that can be held in *sVBuf.
*		WlzGMVertex ***sVBuf:	Ptr to an allocated buffer for
*					verticies, may NOT be NULL
*					although the buffer it points
*					to may be. The buffer should
*					be free'd using AlcFree when
*					it is no longer needed.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzDVertex3	WlzGMVertexNormal3D(WlzGMModel *model, WlzGMVertex *gV,
				    int *sVBufSz, WlzGMVertex ***sVBuf,
				    WlzErrorNum *dstErr)
{
  int		sIdx,
  		sCnt,
  		manifold = 1;
  double	tD0;
  WlzGMVertexT	*vT0,
  		*vT1;
  WlzGMEdgeT	*eT1;
  WlzDVertex3	nrm,
  		cNrm,
  		sNrm;
  WlzDVertex3	sVG[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nrm.vtX = 0.0;
  nrm.vtY = 0.0;
  nrm.vtZ = 0.0;
  if(gV->diskT->idx != gV->diskT->next->idx)
  {
    /* The surface in the imediate neighbourhood of the vertex is
     * not manifold.*/
    manifold = 0;
  }
  else
  {
    /* Find the ordered ring of verticies which surround the current
     * vertex, checking that each of the edges is part of a manifold
     * surface. */

    sCnt = 0;
    vT1 = vT0 = gV->diskT->vertexT;
    eT1 = vT1->parent;
    do
    {
      if((eT1->idx == eT1->rad->idx) ||
	  (eT1->idx != eT1->rad->rad->idx))
      {
	/* More than two loops are joined by the edge, so the surface
	 * in the imediate neighbourhood of the vertex is not manifold.
	 */
	manifold = 0;
      }
      else
      {
	if((*sVBufSz <= sCnt) || (*sVBuf == NULL))
	{
	  *sVBufSz = (*sVBufSz <= 0)? 64: *sVBufSz * 2;
	  if((*sVBuf = (WlzGMVertex **)
		      AlcRealloc(*sVBuf,
		      *sVBufSz * sizeof(WlzGMVertex *))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *(*sVBuf + sCnt++) = eT1->opp->vertexT->diskT->vertex;
	}
	/* Find the next edgeT directed from gV. */
	eT1 = eT1->prev->opp->rad;
	vT1 = eT1->vertexT;
      }
    }
    while((errNum == WLZ_ERR_NONE) && manifold &&
	(vT1->idx != vT0->idx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(manifold)
    {
      /* Compute the mean of the normals of the loops around the
       * current vertex. */
      sNrm.vtX = 0.0;
      sNrm.vtY = 0.0;
      sNrm.vtZ = 0.0;
      (void )WlzGMVertexGetG3D(gV, sVG + 0);
      (void )WlzGMVertexGetG3D(*(*sVBuf + 0), sVG + 2);
      for(sIdx = 0; sIdx < sCnt; ++sIdx)
      {
	*(sVG + 1) = *(sVG + 2);
	(void )WlzGMVertexGetG3D(*(*sVBuf + ((sIdx + 1) % sCnt)), sVG + 2);
	cNrm = WlzGeomTriangleNormal(*(sVG + 0), *(sVG + 1), *(sVG + 2));
	WLZ_VTX_3_ADD(sNrm, sNrm, cNrm);
      }
      tD0 = 1.0 / (double)sCnt;
      WLZ_VTX_3_SCALE(nrm, sNrm, tD0);
    }
    /* else normal undefined. */
  }
  if(dstErr)
  {
    *dstErr == errNum;
  }
  return(nrm);
}

/************************************************************************
* Function:	WlzGMVertexDistSq2D
* Returns:	double:			Square of distance, -1.0 on
*					error.
* Purpose:	Calculates the square of the Euclidean distance
*		between the given vertex and the given 2D double
*		precision position.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Shell with resources.
*		WlzGMVertex *vertex:	Vertex to free.
************************************************************************/
double		WlzGMVertexDistSq2D(WlzGMVertex *vertex, WlzDVertex2 pos)
{
  double	tD0,
		tD1,
		dstSq = -1.0;

  if(vertex && vertex->geo.core)
  {
    switch(vertex->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
        tD0 = pos.vtX - vertex->geo.vg2I->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg2I->vtx.vtY;
	dstSq = (tD0 * tD0) + (tD1 * tD1);
        break;
      case WLZ_GMELM_VERTEX_G2D:
        tD0 = pos.vtX - vertex->geo.vg2D->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg2D->vtx.vtY;
	dstSq = (tD0 * tD0) + (tD1 * tD1);
        break;
      case WLZ_GMELM_VERTEX_G3I:
        tD0 = pos.vtX - vertex->geo.vg3I->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg3I->vtx.vtY;
	dstSq = (tD0 * tD0) + (tD1 * tD1);
        break;
      case WLZ_GMELM_VERTEX_G3D:
        tD0 = pos.vtX - vertex->geo.vg3D->vtx.vtX;
	tD1 = pos.vtY - vertex->geo.vg3D->vtx.vtY;
	dstSq = (tD0 * tD0) + (tD1 * tD1);
        break;
      default:
        break;
    }
  }
  return(dstSq);
}

/************************************************************************
* Function:	WlzGMHashPos3D
* Returns:	unsigned int:		Hash value.
* Purpose:	Computes a hash value from a given 3D double precision
*		position.
* Global refs:	-
* Parameters:	WlzDVertex3 pos:	Given position.
************************************************************************/
static unsigned int WlzGMHashPos3D(WlzDVertex3 pos)
{
  unsigned int	hashVal;
  double	fF,
  		fI;
  const unsigned int pX = 399989, /* These are just different 6 digit primes */
  		pY = 599999,
		pZ = 999983;
  
  fF = modf(pos.vtX, &fI);
  fF = floor(fF / WLZ_GM_TOLERANCE) * WLZ_GM_TOLERANCE;
  fI *= pX;
  fF *= pY * pZ;
  hashVal = ((long long )fI + (long long )fF) & UINT_MAX;
  fF = modf(pos.vtY, &fI);
  fF = floor(fF / WLZ_GM_TOLERANCE) * WLZ_GM_TOLERANCE;
  fI *= pY;
  fF *= pZ * pX;
  hashVal ^= ((long long )fI + (long long )fF) & UINT_MAX;
  fF = modf(pos.vtZ, &fI);
  fF = floor(fF / WLZ_GM_TOLERANCE) * WLZ_GM_TOLERANCE;
  fI *= pZ;
  fF *= pX * pY;
  hashVal ^= ((long long )fI + (long long )fF) & UINT_MAX;
  return(hashVal);
}

/************************************************************************
* Function:	WlzGMHashPos2D
* Returns:	unsigned int:		Hash value.
* Purpose:	Computes a hash value from a given 2D double precision
*		position.
* Global refs:	-
* Parameters:	WlzDVertex2 pos:	Given position.
************************************************************************/
static unsigned int WlzGMHashPos2D(WlzDVertex2 pos)
{
  unsigned int	hashVal;
  double	fF,
  		fI;
  const unsigned int pX = 399989, /* These are just different 6 digit primes */
  		pY = 599999,
		pZ = 999983;
  
  fF = modf(pos.vtX, &fI);
  fF = floor(fF / WLZ_GM_TOLERANCE) * WLZ_GM_TOLERANCE;
  fI *= pX;
  fF *= pY * pZ;
  hashVal = ((long long )fI + (long long )fF) & UINT_MAX;
  fF = modf(pos.vtY, &fI);
  fF = floor(fF / WLZ_GM_TOLERANCE) * WLZ_GM_TOLERANCE;
  fI *= pY;
  fF *= pZ * pX;
  hashVal ^= ((long long )fI + (long long )fF) & UINT_MAX;
  return(hashVal);
}

/************************************************************************
* Function:	WlzGMModelMatchVertexG3D
* Returns:	WlzGMVertex *:		Matched vertex, or NULL if no
*					vertex with the given geometry.
* Purpose:	Attempts to find a vertex which matches the given
*		double precision 3D position.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzDVertex3 *pos:	Position to match.
************************************************************************/
WlzGMVertex	*WlzGMModelMatchVertexG3D(WlzGMModel *model, WlzDVertex3 gPos)
{
  int		cmp;
  unsigned int	hVal;
  WlzGMVertex	*tV,
  		*mV = NULL;

  hVal = WlzGMHashPos3D(gPos);
  if((tV = *(model->vertexHT + (hVal % model->vertexHTSz))) != NULL)
  {
    /* Have found the hash table head, a linked list of verticies with this
     * hash value. Now need to search for a matching vertex in the linked 
     * list. */
    do
    {
      if((cmp = WlzGMVertexCmpSign3D(tV, gPos)) == 0)
      {
        mV = tV;
      }
      else
      {
        tV = tV->next;
      }
    }
    while((cmp < 0) && tV);
  }
  return(mV);
}

/************************************************************************
* Function:	WlzGMModelMatchVertexG2D
* Returns:	WlzGMVertex *:		Matched vertex, or NULL if no
*					vertex with the given geometry.
* Purpose:	Attempts to find a vertex which matches the given
*		double precision 2D position.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzDVertex2 *pos:	Position to match.
************************************************************************/
WlzGMVertex	*WlzGMModelMatchVertexG2D(WlzGMModel *model, WlzDVertex2 gPos)
{
  int		cmp;
  unsigned int	hVal;
  WlzGMVertex	*tV,
  		*mV = NULL;

  hVal = WlzGMHashPos2D(gPos);
  if((tV = *(model->vertexHT + (hVal % model->vertexHTSz))) != NULL)
  {
    /* Have found the hash table head, a linked list of verticies with this
     * hash value. Now need to search for a matching vertex in the linked 
     * list. */
    do
    {
      if((cmp = WlzGMVertexCmpSign2D(tV, gPos)) == 0)
      {
        mV = tV;
      }
      else
      {
        tV = tV->next;
      }
    }
    while((cmp < 0) && tV);
  }
  return(mV);
}

/* Topology query */

/************************************************************************
* Function:	WlzGMModelFindNMEdges
* Returns:	WlzGMEdge *:		Ptr to an array of non-manifold
*					edge ptrs, NULL if none exist or
*					on error.
* Purpose:	Finds a loop topology element in common for the two
*		edge topology elements.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The given model.
*		int *dstNMCnt:		Destination pointer for number
*					of non-manifold edges found.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMEdge	**WlzGMModelFindNMEdges(WlzGMModel *model, int *dstNMCnt,
				        WlzErrorNum *dstErr)
{
  int		tEI,
  		nmEI;
  int		*nmEIP;
  WlzGMEdge	*tE;
  WlzGMEdge	**nmE = NULL;
  WlzGMEdgeT	*tET;
  AlcVector	*nmEIVec = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nmEIVec = AlcVectorNew(1, sizeof(int), 1024, NULL)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(model->type)
    {
      case WLZ_GMMOD_2I:
      case WLZ_GMMOD_2D:
        break;
      case WLZ_GMMOD_3I:
      case WLZ_GMMOD_3D:
	tEI = 0;
	nmEI = 0;
	while((errNum == WLZ_ERR_NONE) && (tEI < model->res.edge.numIdx))
	{
	  tE = (WlzGMEdge *)AlcVectorItemGet(model->res.edge.vec, tEI);
	  if(tE->idx >= 0)
	  {
	    tET = tE->edgeT;
	    if((tET->rad->idx == tET->idx) ||
	       (tET->rad->rad->idx != tET->idx))
	    {
	      if((nmEIP = (int *)
			  (AlcVectorExtendAndGet(nmEIVec, nmEI))) == NULL)
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		*nmEIP = tE->idx;
	      }
	      ++nmEI;
	    }
	  }
	  ++tEI;
	}
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (nmEI > 0))
  {
    if((nmE = (WlzGMEdge **)AlcMalloc(sizeof(WlzGMEdge *) * nmEI)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      for(tEI = 0; tEI < nmEI; ++tEI)
      {
        *(nmE + tEI) = (WlzGMEdge *)AlcVectorItemGet(model->res.edge.vec,
						     tEI);
      }
    }
  }
  if(dstNMCnt)
  {
    *dstNMCnt = (errNum == WLZ_ERR_NONE)? nmEI: 0;
  }
  if(nmEIVec)
  {
    (void )AlcVectorFree(nmEIVec);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nmE);
}

/************************************************************************
* Function:	WlzGMEdgeTCommonLoopT
* Returns:	WlzGMLoopT:		Common loop topology element,
*					NULL if it doesn't exist.
* Purpose:	Finds a loop topology element in common for the two
*		edge topology elements.
* Global refs:	-
* Parameters:	WlzGMEdgeT *eT0:	First edge topology element.
*		WlzGMEdgeT *eT1:	Second edge topology element.
************************************************************************/
WlzGMLoopT	*WlzGMEdgeTCommonLoopT(WlzGMEdgeT *eT0, WlzGMEdgeT *eT1)
{
  WlzGMLoopT	*loopT = NULL;

  if(eT0->parent->idx == eT1->parent->idx)
  {
    loopT = eT0->parent;
  }
  return(loopT);
}

/************************************************************************
* Function:	WlzGMVertexCommonEdge
* Returns:	WlzGMEdge:		Common edge, NULL if it doesn't
					exist.
* Purpose:	Finds the edge common to the two given verticies.
* Global refs:	-
* Parameters:	WlzGMVertex *eV0:	First vertex element.
*		WlzGMVertex *eV1:	Second vertex element.
************************************************************************/
WlzGMEdge	*WlzGMVertexCommonEdge(WlzGMVertex *eV0, WlzGMVertex *eV1)
{
  WlzGMVertexT	*vT0,
		*vT1,
  		*eVT0,
  		*eVT1;
  WlzGMEdge	*e0,
		*edge = NULL;

  vT0 = eVT0 = eV0->diskT->vertexT;
  vT1 = eVT1 = eV1->diskT->vertexT;
  do
  {
    e0 = vT0->parent->edge;
    do
    {
      if(e0->idx == vT1->parent->edge->idx)
      {
	edge = e0;
      }
    } while((edge == NULL) && ((vT1 = vT1->next)->idx != eVT1->idx));
  } while((edge == NULL) && ((vT0 = vT0->next)->idx != eVT0->idx));
  return(edge);
}

/************************************************************************
* Function:	WlzGMVertexCommonShell
* Returns:	WlzGMShell:		Common shell, NULL if it doesn't
					exist.
* Purpose:	Finds the shell common to the two given verticies.
* Global refs:	-
* Parameters:	WlzGMVertex *eV0:	First vertex element.
*		WlzGMVertex *eV1:	Second vertex element.
************************************************************************/
WlzGMShell	*WlzGMVertexCommonShell(WlzGMVertex *eV0, WlzGMVertex *eV1)
{
  WlzGMShell	*s0,
  		*s1,
		*shell = NULL;

  if(((s0 = WlzGMVertexGetShell(eV0)) != NULL) &&
     ((s1 = WlzGMVertexGetShell(eV1)) != NULL) &&
     (s0->idx == s1->idx))
  {
    shell = s0;
  }
  return(shell);
}

/************************************************************************
* Function:	WlzGMVertexGetShell
* Returns:	WlzGMShell:		The parent shell.
* Purpose:	Finds the parent shell of the given vertex.
* Global refs:	-
* Parameters:	WlzGMVertex *eV:	Given vertex element.
************************************************************************/
WlzGMShell	*WlzGMVertexGetShell(WlzGMVertex *eV)
{
  WlzGMShell	*pS = NULL;

  if(eV && (eV->type == WLZ_GMELM_VERTEX))
  {
    pS = eV->diskT->vertexT->parent->parent->parent;
  }
  return(pS);
}

/************************************************************************
* Function:	WlzGMEdgeGetShell
* Returns:	WlzGMShell:		The parent shell.
* Purpose:	Finds the parent shell of the given edge.
* Global refs:	-
* Parameters:	WlzGMEdge *eE:		Given edge element.
************************************************************************/
WlzGMShell	*WlzGMEdgeGetShell(WlzGMEdge *eE)
{
  WlzGMShell	*pS = NULL;

  if(eE && (eE->type == WLZ_GMELM_EDGE))
  {
    pS = eE->edgeT->parent->parent;
  }
  return(pS);
}

/************************************************************************
* Function:	WlzGMEdgeCommonVertex
* Returns:	WlzGMVertex:		The common vertex, NULL if it
*					doesn't exist.
* Purpose:	Finds the common vertex of the two given edges.
* Global refs:	-
* Parameters:	WlzGMEdge *eE0:		First edge element.
*		WlzGMEdge *eE1:		Second edge element.
************************************************************************/
WlzGMVertex	*WlzGMEdgeCommonVertex(WlzGMEdge *eE0, WlzGMEdge *eE1)
{
  WlzGMVertex	*tV,
  		*cV = NULL;

  if(((tV = eE0->edgeT->vertexT->diskT->vertex)->idx ==
      eE1->edgeT->vertexT->diskT->vertex->idx) ||
     (tV->idx == eE1->edgeT->opp->vertexT->diskT->vertex->idx))
  {
    cV = tV;
  }
  else if(((tV = eE0->edgeT->opp->vertexT->diskT->vertex)->idx ==
           eE1->edgeT->vertexT->diskT->vertex->idx) ||
          (tV->idx == eE1->edgeT->opp->vertexT->diskT->vertex->idx))
  {
    cV = tV;
  }
  return(cV);
}

/* Model list management */

/************************************************************************
* Function:	WlzGMVertexTAppend
* Returns:	void
* Purpose:	Append new vertex topology element onto a doubly
*		linked list of vertex topology element's, knowing that
*		neither is NULL.
*		Vertex topology elements are maintained in an unordered
*		doubly linked list.
* Global refs:	-
* Parameters:	WlzGMVertexT *eVT:	Existing vertexT in list of
*					vertexT's.
*		WlzGMVertexT *nVT:	New vertex topology element to
*					append to list after existing
*					element.
************************************************************************/
void	   	WlzGMVertexTAppend(WlzGMVertexT *eVT, WlzGMVertexT *nVT)
{
  nVT->next = eVT->next;
  nVT->prev = eVT;
  nVT->next->prev = nVT;
  nVT->prev->next = nVT;
}

/************************************************************************
* Function:	WlzGMDiskTAppend
* Returns:	void
* Purpose:	Append new edge topology onto a doubly linked list
*		of disk topology element's, knowing that neither is NULL.
*		Disk topology elements are maintained in an unordered
*		doubly linked list.
* Global refs:	-
* Parameters:	WlzGMDiskT *eDT:	Existing disk topology element in
*					list.
*		WlzGMDiskT *nDT:	New disk topology element to
*					append to list after existing
*					element.
************************************************************************/
void	   	WlzGMDiskTAppend(WlzGMDiskT *eDT, WlzGMDiskT *nDT)
{
  nDT->next = eDT->next;
  nDT->prev = eDT;
  nDT->next->prev = nDT;
  nDT->prev->next = nDT;
}

/************************************************************************
* Function:	WlzGMEdgeTAppend
* Returns:	void
* Purpose:	Append new edge topology onto a doubly linked list
*		of edge topology element's, knowing that neither is NULL.
*		Edge topology elements are maintained in an ordered
*		doubly linked list, CCW around the inside of loops and
*		CW around the outside of 2D loops.
* Global refs:	-
* Parameters:	WlzGMEdgeT *eET:	Existing edge topology element in
*					list.
*		WlzGMEdgeT *nET:	New edge topology element to
*					append to list after existing
*					element.
************************************************************************/
void	   	WlzGMEdgeTAppend(WlzGMEdgeT *eET, WlzGMEdgeT *nET)
{
  nET->next = eET->next;
  nET->prev = eET;
  nET->next->prev = nET;
  nET->prev->next = nET;
}

/************************************************************************
* Function:	WlzGMEdgeTInsert
* Returns:	void
* Purpose:	Insert new edge topology into a doubly linked list
*		of edge topology element's, knowing that neither is NULL.
*		Edge topology elements are maintained in an ordered
*		doubly linked list, CCW around the inside of loops and
*		CW around the outside of 2D loops.
* Global refs:	-
* Parameters:	WlzGMEdgeT *eET:	Existing edge topology element in
*					list.
*		WlzGMEdgeT *nET:	New edge topology element to
*					insert in list before existing
*					element.
************************************************************************/
void	   	WlzGMEdgeTInsert(WlzGMEdgeT *eET, WlzGMEdgeT *nET)
{
  nET->next = eET;
  nET->prev = eET->prev;
  nET->next->prev = nET;
  nET->prev->next = nET;
}

/************************************************************************
* Function:	WlzGMEdgeTInsertRadial
* Returns:	void
* Purpose:	Insert the given new edge topology element into a
*		radially sorted cyclic list of edge topology element's.
*		In 2D the radial edgeT is always the edgeT itself and
*		this function should not be called.
*		In 3D the radial edge is a similarly directed edge
*		topology element to the given edge topology element
*		who's off edge vertex is the next CCW off edge vertex
*		when viewed along the edge.
*                                                                 
*                                      O <-                       
*                                     /    \                      
*                                    /      \ Next radial edge    
*                                   /        |                    
*                                  /         |                    
*                                 X----------O                    
*                           Given edge into    Off edge vertex
*                           screen/page        of new simplex
*                                                                 
*		Given a directed edge which forms part of a loop in
*		a new triangle simplex. Define points around the
*		the simplex:
*                                                                 
*                         O p1                                     
*                       ^ |\                                      
*                       | | \                                     
*                       | |  \                                    
*                       | |   \                                   
*                       | |    \                                  
*                   nET | |     O p2                               
*                       | |    /                                  
*                       | |   /                                   
*                       | |  /                                    
*                       | | /                                     
*                       | |/                                      
*                         O p0                                    
*
* 		To find the next CCW simplex.
*		Find the normal vector:
*		  n0 = (p0 - p1) / |p0 - p1|
*		Find the normal vector perpendicular to the new
*		(triangular) simplex:
*		        n0 x (p2 - p0)
*		  n1 = ----------------
*		       |n0 x (p2 - p0)|
*		Find the normalized vector perpendicular to n0 and n1:
*		  n2 = n0 x n1
*		Because n1 and n2 are perpendicular axes in the plane
*		perpendicular to p1 - p0 we can compute the angle of
*		any vector projected onto the plane as:
*		             -1 / vi . n2 \
*		  theta = tan  | --------- |, 0 <= theta <= 2PI
*			        \ vi . n1 /
*		where
*		  vi = p2i - p0
*		and p2i are the p2 verticies on the candidate loop/face.
*		Use the projections to find the previous radial edge,
*		then insert the edge in the cyclic radial list.
* Global refs:	-
* Parameters:	WlzGMEdgeT *eET:	Existing edge topology element in
*					list. which MUST have all it's
*					fields set EXCEPT for the radial
*					edge link.
*		WlzGMEdgeT *nET:	New edge topology element to
*					insert in list before existing
*					element.
************************************************************************/
void	   	WlzGMEdgeTInsertRadial(WlzGMEdgeT *nET)
{
  double	tD0;
  WlzDVertex2	nGrd,
  		pGrd,
		tGrd;
  WlzDVertex3	p0,
  		p1,
		p2,
		n0,
		n1,
		n2,
		v2;
  WlzGMEdgeT	*pET,
		*fET,
		*tET;

  /* Get an edgeT of this edge with same direction as the given edgeT. */
  fET = nET->edge->edgeT;
  if(fET->vertexT->diskT->vertex->idx != nET->vertexT->diskT->vertex->idx)
  {
    /* fET is antiparrallel to the new edgeT so use it's opposite. */
    fET = fET->opp;
  }
  if(fET->rad->idx == fET->idx)
  {
    /* Edge not shared. */
    pET = fET;
  }
  else
  {
    /* Get the geometry of the new edgeT's simplex. */
    (void )WlzGMVertexGetG3D(nET->vertexT->diskT->vertex, &p0);
    (void )WlzGMVertexGetG3D(nET->next->vertexT->diskT->vertex, &p1);
    (void )WlzGMVertexGetG3D(nET->prev->vertexT->diskT->vertex, &p2);
    /* Compute the normal vectors: n0, n1 and n2. */
    WLZ_VTX_3_SUB(n0, p0, p1);
    tD0 = 1.0 / WLZ_VTX_3_LENGTH(n0);
    WLZ_VTX_3_SCALE(n0, n0, tD0);
    WLZ_VTX_3_SUB(v2, p2, p0);
    WLZ_VTX_3_CROSS(n1, n0, v2);
    tD0 = 1.0 / WLZ_VTX_3_LENGTH(n1);
    WLZ_VTX_3_SCALE(n1, n1, tD0);
    WLZ_VTX_3_CROSS(n2, n0, n1);
    /* Compute the angle from the origin to the point projected onto the
     * plane (defined by n1 and n2) by p2. */
    nGrd.vtX = WLZ_VTX_3_DOT(v2, n1);
    nGrd.vtY = WLZ_VTX_3_DOT(v2, n2);
    /* Search for the edgeT that should be the previous radial edge to the
     * new edgeT by searching for the edgeT which forms part of a simplex
     * with the smallest CW angle when viewed along the directed edgeT. */
    pET = tET = fET;
    (void )WlzGMVertexGetG3D(pET->prev->vertexT->diskT->vertex, &p2);
    WLZ_VTX_3_SUB(v2, p2, p0);
    pGrd.vtX = WLZ_VTX_3_DOT(v2, n1);
    pGrd.vtY = WLZ_VTX_3_DOT(v2, n2);
    while((tET = tET->rad)->idx != fET->idx)
    {
      (void )WlzGMVertexGetG3D(tET->prev->vertexT->diskT->vertex, &p2);
      WLZ_VTX_3_SUB(v2, p2, p0);
      tGrd.vtX = WLZ_VTX_3_DOT(v2, n1);
      tGrd.vtY = WLZ_VTX_3_DOT(v2, n2);
      /* Check to see if the test angle (gradient tGrd), is CCW wrt
       * the previous angle (gradient pGrd) and CW wrt the new angle
       * (gradient nGrd). If so the test angle and edgeT become the new
       * previous values. Do this by testing the orientation of the
       * triangle. */
      if(WlzGeomTriangleSnArea2(nGrd, pGrd, tGrd) > 0)
      {
	pET = tET;
        pGrd = tGrd;
      }
    }
  }
  /* Insert nET after pET and before pET->rad. */
  nET->rad = pET->rad;
  pET->rad = nET;
}

/************************************************************************
* Function:	WlzGMLoopTAppend
* Returns:	void
* Purpose:	Append new loop topology onto a doubly linked list
*		of loop topology element's, knowing that neither is NULL.
*		Loop topology elements are maintained in an unordered
*		doubly linked list.
* Global refs:	-
* Parameters:	WlzGMLoopT *eLT:	Existing loop topology element in
*					list.
*		WlzGMLoopT *nLT:	New loop topology element to
*					append to list after existing
*					element.
************************************************************************/
void	   	WlzGMLoopTAppend(WlzGMLoopT *eLT, WlzGMLoopT *nLT)
{
  nLT->next = eLT->next;
  nLT->prev = eLT;
  nLT->next->prev = nLT;
  nLT->prev->next = nLT;
}

/************************************************************************
* Function:	WlzGMLoopTUnlink
* Returns:	void
* Purpose:	Unlinks the given loop topology element from a doubly
*		linked list of loop topology element's, knowing that
*		it is not NULL. If this is the only loopT in parent
*		the parent's list of loopT's is set to NULL.
* Global refs:	-
* Parameters:	WlzGMLoopT *dLT:	Loop topology element to be
*					unlinked.
************************************************************************/
void		WlzGMLoopTUnlink(WlzGMLoopT *dLT)
{
  dLT->prev->next = dLT->next;
  dLT->next->prev = dLT->prev;
  if(dLT->idx == dLT->next->idx)
  {
    dLT->next = dLT->prev = NULL;
  }
  if(dLT->parent->child->idx == dLT->idx)
  {
    dLT->parent->child = dLT->next;
  }
}

/************************************************************************
* Function:	WlzGMShellAppend
* Returns:	void
* Purpose:	Append new shell onto a doubly linked list of
*		shells, knowing that neither is NULL.
*		A model's shells are maintained in an unordered
*		doubly linked list.
* Global refs:	-
* Parameters:	WlzGMShell *eS:		Existing shell in list of shells.
*		WlzGMShell *nS:		New shell to append to list
*					after existing shell.
************************************************************************/
void	   	WlzGMShellAppend(WlzGMShell *eS, WlzGMShell *nS)
{
  nS->next = eS->next;
  nS->prev = eS;
  nS->next->prev = nS;
  nS->prev->next = nS;
}

/************************************************************************
* Function:	WlzGMShellUnlink
* Returns:	void
* Purpose:	Unlinks the given shell from a doubly linked list of
*		shells kept by the model, knowing that the given shell
*		ptr is not NULL.
* Global refs:	-
* Parameters:	WlzGMShell *dS:		Shell to be unlinked.
************************************************************************/
void		WlzGMShellUnlink(WlzGMShell *dS)
{
  if(dS->idx == dS->parent->child->idx)
  {
    /* Change child held by parent model */
    if(dS->idx == dS->next->idx)
    {
      /* The model only has one shell and this is it. */
      dS->parent = NULL;
    }
    else
    {
      dS->parent->child = dS->next;
    }
  }
  if(dS->next != NULL)
  {
    /* Unlink the shell from the doubly linked list of shells. */
    dS->prev->next = dS->next;
    dS->next->prev = dS->prev;
  }
  dS->next = dS->prev = NULL;
}

/************************************************************************
* Function:	WlzGMShellJoinAndUnlink
* Returns:	void
* Purpose:	Joins the shell to be unlinked onto the shell that's to
*		be extended and then unlinks (deletes) it.
* Global refs:	-
* Parameters:	WlzGMShell *eS:		Shell to be extended.
*		WlzGMShell *dS:		Shell to be joined and then
*					unlinked.
************************************************************************/
void		WlzGMShellJoinAndUnlink(WlzGMShell *eShell, WlzGMShell *dShell)
{
  WlzGMLoopT	*eLT,
  		*dLT,
		*tLT;
  WlzDBox3	bBox;

  if(eShell && dShell &&
     ((eLT = eShell->child) != NULL) && ((dLT = dShell->child) != NULL))
  {
    /* Update extended shells geometry. */
    (void )WlzGMShellGetGBB3D(dShell, &bBox);
    (void )WlzGMShellUpdateGBB3D(eShell, bBox);
    /* Update extended shells topology. */
    tLT = dLT;
    do
    {
      tLT->parent = eShell;
      tLT = tLT->next;
    } while(tLT->idx != dLT->idx);
    tLT = eLT->prev;
    tLT->next = dLT;
    eLT->prev = dLT->prev;
    tLT->next->prev = tLT;
    eLT->prev->next = eLT;
    /* Unlink the now childless shell. */
    WlzGMShellUnlink(dShell);
  }
}

/************************************************************************
* Function:	WlzGMModelResIdx
* Returns:	WlzGMResIdxTb *:	New resource index table, NULL
*					on error.
* Purpose:	Makes an index look up table data structure for the
*		given model.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model with the vertex
*					resources.
*		unsigned int eMsk:	Element mask with bits set for
*					elemet resources to index.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMResIdxTb	*WlzGMModelResIdx(WlzGMModel *model, unsigned int eMsk,
				  WlzErrorNum *dstErr)
{
  int		iCnt,
		vCnt;
  unsigned int 	eMskTst;
  int		*iLut;
  AlcVector	*vec;
  WlzGMElemP	eP;
  WlzGMResIdxTb	*resIdxTb = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(model->type)
    {
      case WLZ_GMMOD_2I:
      case WLZ_GMMOD_2D:
      case WLZ_GMMOD_3I:
      case WLZ_GMMOD_3D:
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    eMskTst = eMsk;
    if((resIdxTb = (WlzGMResIdxTb *)
    		   AlcCalloc(1, sizeof(WlzGMResIdxTb))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Set initial element counts and allocate the index luts. */
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_VERTEX))
  {
    eMskTst &= ~(WLZ_GMELMFLG_VERTEX);
    vec = model->res.vertex.vec;
    resIdxTb->vertex.idxCnt = model->res.vertex.numIdx;
    if((resIdxTb->vertex.idxLut = (int *)
	    	AlcCalloc(resIdxTb->vertex.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_VERTEX_T))
  {
    eMskTst &= ~(WLZ_GMELMFLG_VERTEX_T);
    vec = model->res.vertexT.vec;
    resIdxTb->vertexT.idxCnt = model->res.vertexT.numIdx;
    if((resIdxTb->vertexT.idxLut = (int *)
    		AlcCalloc(resIdxTb->vertexT.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_VERTEX_G))
  {
    eMskTst &= ~(WLZ_GMELMFLG_VERTEX_G);
    vec = model->res.vertexG.vec;
    resIdxTb->vertexG.idxCnt = model->res.vertexG.numIdx;
    if((resIdxTb->vertexG.idxLut = (int *)
    		AlcCalloc(resIdxTb->vertexG.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_DISK_T))
  {
    eMskTst &= ~(WLZ_GMELMFLG_DISK_T);
    vec = model->res.diskT.vec;
    resIdxTb->diskT.idxCnt = model->res.diskT.numIdx;
    if((resIdxTb->diskT.idxLut = (int *)
	    	AlcCalloc(resIdxTb->diskT.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_EDGE))
  {
    eMskTst &= ~(WLZ_GMELMFLG_EDGE);
    vec = model->res.edge.vec;
    resIdxTb->edge.idxCnt = model->res.edge.numIdx;
    if((resIdxTb->edge.idxLut = (int *)
	    	AlcCalloc(resIdxTb->edge.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_EDGE_T))
  {
    eMskTst &= ~(WLZ_GMELMFLG_EDGE_T);
    vec = model->res.edgeT.vec;
    resIdxTb->edgeT.idxCnt = model->res.edgeT.numIdx;
    if((resIdxTb->edgeT.idxLut = (int *)
	    	AlcCalloc(resIdxTb->edgeT.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_LOOP))
  {
    eMskTst &= ~(WLZ_GMELMFLG_LOOP);
    vec = model->res.loop.vec;
    resIdxTb->loop.idxCnt = model->res.loop.numIdx;
    if((resIdxTb->loop.idxLut = (int *)
	    	AlcCalloc(resIdxTb->loop.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_LOOP_T))
  {
    eMskTst &= ~(WLZ_GMELMFLG_LOOP_T);
    vec = model->res.loopT.vec;
    resIdxTb->loopT.idxCnt = model->res.loopT.numIdx;
    if((resIdxTb->loopT.idxLut = (int *)
	    	AlcCalloc(resIdxTb->loopT.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_SHELL))
  {
    eMskTst &= ~(WLZ_GMELMFLG_SHELL);
    vec = model->res.shell.vec;
    resIdxTb->shell.idxCnt = model->res.shell.numIdx;
    if((resIdxTb->shell.idxLut = (int *)
	    	AlcCalloc(resIdxTb->shell.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_SHELL_G))
  {
    eMskTst &= ~(WLZ_GMELMFLG_SHELL_G);
    vec = model->res.shellG.vec;
    resIdxTb->shellG.idxCnt = model->res.shellG.numIdx;
    if((resIdxTb->shellG.idxLut = (int *)
	    	AlcCalloc(resIdxTb->shellG.idxCnt, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (eMskTst != 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(eMsk & WLZ_GMELMFLG_VERTEX)
    {
      /* Compute vertex index lut. */
      vCnt = iCnt = 0;
      vec = model->res.vertex.vec;
      iLut = resIdxTb->vertex.idxLut;
      while(vCnt < resIdxTb->vertex.idxCnt)
      {
	eP.vertex = (WlzGMVertex *)AlcVectorItemGet(vec, vCnt++);
	if(eP.vertex->idx >= 0)
	{
	  *(iLut + eP.vertex->idx) = iCnt++;
	}
      }
      resIdxTb->vertex.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_VERTEX_T)
    {
      /* Compute vertexT index lut. */
      vCnt = iCnt = 0;
      vec = model->res.vertexT.vec;
      iLut = resIdxTb->vertexT.idxLut;
      while(vCnt < resIdxTb->vertexT.idxCnt)
      {
	eP.vertexT = (WlzGMVertexT *)AlcVectorItemGet(vec, vCnt++);
	if(eP.vertexT->idx >= 0)
	{
	  *(iLut + eP.vertexT->idx) = iCnt++;
	}
      }
      resIdxTb->vertexT.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_VERTEX_G)
    {
      /* Compute vertexG index lut. */
      vCnt = iCnt = 0;
      vec = model->res.vertexG.vec;
      iLut = resIdxTb->vertexG.idxLut;
      switch(model->type)
      {
	case WLZ_GMMOD_2I:
	  while(vCnt < resIdxTb->vertexG.idxCnt)
	  {
	    eP.vertexG2I = (WlzGMVertexG2I *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.vertexG2I->idx >= 0)
	    {
	      *(iLut + eP.vertexG2I->idx) = iCnt++;
	    }
	  }
	  break;
	case WLZ_GMMOD_2D:
	  while(vCnt < resIdxTb->vertexG.idxCnt)
	  {
	    eP.vertexG2D = (WlzGMVertexG2D *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.vertexG2D->idx >= 0)
	    {
	      *(iLut + eP.vertexG2D->idx) = iCnt++;
	    }
	  }
	  break;
	case WLZ_GMMOD_3I:
	  while(vCnt < resIdxTb->vertexG.idxCnt)
	  {
	    eP.vertexG3I = (WlzGMVertexG3I *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.vertexG3I->idx >= 0)
	    {
	      *(iLut + eP.vertexG3I->idx) = iCnt++;
	    }
	  }
	  break;
	case WLZ_GMMOD_3D:
	  while(vCnt < resIdxTb->vertexG.idxCnt)
	  {
	    eP.vertexG3D = (WlzGMVertexG3D *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.vertexG3D->idx >= 0)
	    {
	      *(iLut + eP.vertexG3D->idx) = iCnt++;
	    }
	  }
	  break;
      }
      resIdxTb->vertexG.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_DISK_T)
    {
      /* Compute diskT index lut. */
      vCnt = iCnt = 0;
      vec = model->res.diskT.vec;
      iLut = resIdxTb->diskT.idxLut;
      while(vCnt < resIdxTb->diskT.idxCnt)
      {
	eP.diskT = (WlzGMDiskT *)AlcVectorItemGet(vec, vCnt++);
	if(eP.diskT->idx >= 0)
	{
	  *(iLut + eP.diskT->idx) = iCnt++;
	}
      }
      resIdxTb->diskT.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_EDGE)
    {
      /* Compute edge index lut. */
      vCnt = iCnt = 0;
      vec = model->res.edge.vec;
      iLut = resIdxTb->edge.idxLut;
      while(vCnt < resIdxTb->edge.idxCnt)
      {
	eP.edge = (WlzGMEdge *)AlcVectorItemGet(vec, vCnt++);
	if(eP.edge->idx >= 0)
	{
	  *(iLut + eP.edge->idx) = iCnt++;
	}
      }
      resIdxTb->edge.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_EDGE_T)
    {
      /* Compute edgeT index lut. */
      iCnt = vCnt = 0;
      vec = model->res.edgeT.vec;
      iLut = resIdxTb->edgeT.idxLut;
      while(vCnt < resIdxTb->edgeT.idxCnt)
      {
	eP.edgeT = (WlzGMEdgeT *)AlcVectorItemGet(vec, vCnt++);
	if(eP.edgeT->idx >= 0)
	{
	  *(iLut + eP.edgeT->idx) = iCnt++;
	}
      }
      resIdxTb->edgeT.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_LOOP)
    {
      /* Compute loop index lut. */
      iCnt = vCnt = 0;
      vec = model->res.loop.vec;
      iLut = resIdxTb->loop.idxLut;
      while(vCnt < resIdxTb->loop.idxCnt)
      {
	eP.loop = (WlzGMLoop *)AlcVectorItemGet(vec, vCnt++);
	if(eP.loop->idx >= 0)
	{
	  *(iLut + eP.loop->idx) = iCnt++;
	}
      }
      resIdxTb->loop.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_LOOP_T)
    {
      /* Compute loopT index lut. */
      iCnt = vCnt = 0;
      vec = model->res.loopT.vec;
      iLut = resIdxTb->loopT.idxLut;
      while(vCnt < resIdxTb->loopT.idxCnt)
      {
	eP.loopT = (WlzGMLoopT *)AlcVectorItemGet(vec, vCnt++);
	if(eP.loopT->idx >= 0)
	{
	  *(iLut + eP.loopT->idx) = iCnt++;
	}
      }
      resIdxTb->loopT.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_SHELL)
    {
      /* Compute shell index lut. */
      iCnt = vCnt = 0;
      vec = model->res.shell.vec;
      iLut = resIdxTb->shell.idxLut;
      while(vCnt < resIdxTb->shell.idxCnt)
      {
	eP.shell = (WlzGMShell *)AlcVectorItemGet(vec, vCnt++);
	if(eP.shell->idx >= 0)
	{
	  *(iLut + eP.shell->idx) = iCnt++;
	}
      }
      resIdxTb->shell.idxCnt = iCnt;
    }
    if(eMsk & WLZ_GMELMFLG_SHELL_G)
    {
      /* Compute shellG index lut. */
      iCnt = vCnt = 0;
      vec = model->res.shellG.vec;
      iLut = resIdxTb->shellG.idxLut;
      switch(model->type)
      {
	case WLZ_GMMOD_2I:
	  while(vCnt < resIdxTb->shellG.idxCnt)
	  {
	    eP.shellG2I = (WlzGMShellG2I *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.shellG2I->idx >= 0)
	    {
	      *(iLut + eP.shellG2I->idx) = iCnt++;
	    }
	  }
	  break;
	case WLZ_GMMOD_2D:
	  while(vCnt < resIdxTb->shellG.idxCnt)
	  {
	    eP.shellG2D = (WlzGMShellG2D *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.shellG2D->idx >= 0)
	    {
	      *(iLut + eP.shellG2D->idx) = iCnt++;
	    }
	  }
	  break;
	case WLZ_GMMOD_3I:
	  while(vCnt < resIdxTb->shellG.idxCnt)
	  {
	    eP.shellG3I = (WlzGMShellG3I *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.shellG3I->idx >= 0)
	    {
	      *(iLut + eP.shellG3I->idx) = iCnt++;
	    }
	  }
	  break;
	case WLZ_GMMOD_3D:
	  while(vCnt < resIdxTb->shellG.idxCnt)
	  {
	    eP.shellG3D = (WlzGMShellG3D *)AlcVectorItemGet(vec, vCnt++);
	    if(eP.shellG3D->idx >= 0)
	    {
	      *(iLut + eP.shellG3D->idx) = iCnt++;
	    }
	  }
	  break;
      }
      resIdxTb->shellG.idxCnt = iCnt;
    }
  }
  if((errNum != WLZ_ERR_NONE) && resIdxTb)
  {
    WlzGMModelResIdxFree(resIdxTb);
    resIdxTb = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(resIdxTb);
}

/************************************************************************
* Function:	WlzGMModelResIdxFree
* Returns:	void
* Purpose:	Frees a GM index look up table data structure.
* Global refs:	-
* Parameters:	WlzGMResIdxTb *resIdxTb: Given index lut data structure.
************************************************************************/
void		WlzGMModelResIdxFree(WlzGMResIdxTb *resIdxTb)
{
  if(resIdxTb)
  {
    if(resIdxTb->vertex.idxLut)
    {
      AlcFree(resIdxTb->vertex.idxLut);
    }
    if(resIdxTb->vertexT.idxLut)
    {
      AlcFree(resIdxTb->vertexT.idxLut);
    }
    if(resIdxTb->vertexG.idxLut)
    {
      AlcFree(resIdxTb->vertexG.idxLut);
    }
    if(resIdxTb->diskT.idxLut)
    {
      AlcFree(resIdxTb->diskT.idxLut);
    }
    if(resIdxTb->edge.idxLut)
    {
      AlcFree(resIdxTb->edge.idxLut);
    }
    if(resIdxTb->edgeT.idxLut)
    {
      AlcFree(resIdxTb->edgeT.idxLut);
    }
    if(resIdxTb->loop.idxLut)
    {
      AlcFree(resIdxTb->loop.idxLut);
    }
    if(resIdxTb->shell.idxLut)
    {
      AlcFree(resIdxTb->shell.idxLut);
    }
    if(resIdxTb->shellG.idxLut)
    {
      AlcFree(resIdxTb->shellG.idxLut);
    }
  }
}

/************************************************************************
* Function:	WlzGMModelRehashVHT
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Rehash the vertex matching hash table.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model.
*		int vHTSz:		New vertex matching hash table
*					size, no change if <= 0.
************************************************************************/
WlzErrorNum 	WlzGMModelRehashVHT(WlzGMModel *model, int vHTSz)
{
  int		idV,
  		vCnt;
  WlzGMVertex	*vertex;
  WlzGMVertex	**newVHT,
  		**oldVHT;
  AlcVector	*vec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    if(vHTSz > 0)
    {
      if((newVHT = AlcCalloc(vHTSz, sizeof(WlzGMVertex *))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	model->vertexHTSz = vHTSz;
	oldVHT = model->vertexHT;
	model->vertexHT = newVHT;
	AlcFree(oldVHT);
      }
    }
    else
    {
      (void )memset(model->vertexHT, 0,
      		    model->vertexHTSz * sizeof(WlzGMVertex *));
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idV = 0;
    vec = model->res.vertex.vec;
    vCnt = model->res.vertex.numIdx;
    while(idV < vCnt)
    {
      vertex = (WlzGMVertex *)AlcVectorItemGet(vec, idV);
      if(vertex->idx >= 0)
      {
	vertex->next = NULL;
	WlzGMModelAddVertex(model, vertex);
      }
      ++idV;
    }
  }
  return(errNum);
}

/* Model construction */

/************************************************************************
* Function:	WlzGMModelConstructNewS3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a new shell, face, loop and 3 edges (with
*		topological elements) and adds them to the model. Given
*		3 3D double precision points which are known not to
*		be in the model.
*		None of the verticies already exists within the shell.
*		Need to create a new shell (nShell) with: 1 loop (nL),
*		2 loop topology element (nLT0, nLT1), 3 edges (*nE),
*		6 edge topology elements (*nET0, *nET1), 3 disk
*		topology elements *(nDT), 3 verticies (*nV) with
*		geometry (*pos) and 6 vertex topology elements
*		(*nVT0, *nVT1).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   nV0 @- - - - - - - - - nE0 - - - - - - - -@ nV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                    nE2 ----\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE1
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            nV2
*
*                   DT0 = {nVT00, nVT10},
*                   DT1 = {nVT01, nVT11},
*                   DT2 = {nVT02, nVT12},
*                   LT0 = {nET00, nET01, nET02},
*                   LT1 = {nET10, nET11, nET12}
*
*                                                                      
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzDVertex3 *pos:	Positions of the 3 points.
************************************************************************/
static WlzErrorNum WlzGMModelConstructNewS3D(WlzGMModel *model,
					     WlzDVertex3 *pos)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell,
  		*nShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertex	*nV[3];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nShell = WlzGMModelNewS(model, &errNum)) != NULL) &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[2] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[2] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nV[0] = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nV[1] = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nV[2] = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* New verticies. */
    for(idx = 0; idx < 3; ++idx)
    {
      (void )WlzGMVertexSetG3D(nV[idx], *(pos + idx));
      nV[idx]->diskT = nDT[idx];
      WlzGMModelAddVertex(model, nV[idx]);
    }
    /* New vertex topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->next = nVT0[idx]->prev = nVT1[idx];
      nVT1[idx]->next = nVT1[idx]->prev = nVT0[idx];
      nVT0[idx]->diskT = nVT1[idx]->diskT = nDT[idx];
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nDT[idx]->next = nDT[idx]->prev = nDT[idx];
      nDT[idx]->vertex = nV[idx];
      nDT[idx]->vertexT = nVT0[idx];
    }
    /* New edges. */
    for(idx = 0; idx < 3; ++idx)
    {
      nE[idx]->edgeT = nET0[idx];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->edge = nE[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->edge = nE[pIdx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    nLT[0]->next = nLT[0]->prev = nLT[1];
    nLT[1]->next = nLT[1]->prev = nLT[0];
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = nShell;
    nLT[1]->parent = nShell;
    /* New shell in model */
    (void )WlzGMShellSetG3D(nShell, 3, pos);
    nShell->child = nLT[0];
    nShell->parent = model;
    nShell->next = nShell->prev = nShell;
    /* Is this the first shell of the model? */
    if((eShell = model->child) == NULL)
    {
      /* First shell in model */
      model->child = nShell;
    }
    else
    {
      /* Append new shell onto end of model's list of shells. */
      WlzGMShellAppend(eShell, nShell);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelExtend1V0E1S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Extends an existing shell within the model by adding a:
*		loop, 3 edges and 2 verticies (with topological elements).
*		The 2 given 3D double precision points which are known
*		not to be in the model.
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology elements
*		(nDT[0,1,2]), 2 verticies (nV[0,1]) with geometry
*		(pos0, pos1) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   nV0 @- - - - - - - - - nE0 - - - - - - - -@ nV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                    nE2 ----\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE1
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            eV
*
*                   DT0 = {nVT00, nVT10},
*                   DT1 = {nVT01, nVT11},
*                   DT2 = {nVT02, nVT12},
*                   LT0 = {nET00, nET01, nET02},
*                   LT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMVertex *eV:	Shared vertex.
*		WlzDVertex3 pos0:	Position of the first point.
*		WlzDVertex3 pos1:	Position of the second point.
************************************************************************/
static WlzErrorNum WlzGMModelExtend1V0E1S3D(WlzGMModel *model, WlzGMVertex *eV,
					    WlzDVertex3 pos0, WlzDVertex3 pos1)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertex	*nV[2];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[2] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[2] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nV[0] = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nV[1] = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    pos[0] = pos0;
    pos[1] = pos1;
    (void )WlzGMVertexGetG3D(eV, pos + 2);
    /* Get shell that's to be extended. */
    eShell = WlzGMVertexGetShell(eV);
    /* New verticies. */
    for(idx = 0; idx < 2; ++idx)
    {
      (void )WlzGMVertexSetG3D(nV[idx], pos[idx]);
      nV[idx]->diskT = nDT[idx];
      WlzGMModelAddVertex(model, nV[idx]);
    }
    /* New vertex topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->next = nVT0[idx]->prev = nVT1[idx];
      nVT1[idx]->next = nVT1[idx]->prev = nVT0[idx];
      nVT0[idx]->diskT = nVT1[idx]->diskT = nDT[idx];
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology elements */
    for(idx = 0; idx < 2; ++idx)
    {
      nDT[idx]->next = nDT[idx]->prev = nDT[idx];
      nDT[idx]->vertex = nV[idx];
      nDT[idx]->vertexT = nVT0[idx];
    }
    WlzGMDiskTAppend(eV->diskT, nDT[2]);
    nDT[2]->vertex = eV;
    nDT[2]->vertexT = nVT0[2];
    /* New edges. */
    for(idx = 0; idx < 3; ++idx)
    {
      nE[idx]->edgeT = nET0[idx];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->edge = nE[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->edge = nE[pIdx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update shell */
    (void )WlzGMShellUpdateG3D(eShell, pos0);
    (void )WlzGMShellUpdateG3D(eShell, pos1);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelExtend2V1E1S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Extends an existing shell within the model by adding a:
*		loop, 2 edges and 1 vertex (with topological elements).
*		The given 3D double precision point is known not to be
*		in the model.
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 2 edges (nE[0,1]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 1 disk topology element
*		(nDT), 1 vertex (nV) with geometry (pos0) and
*		6 vertex topology elements (nVT[0,1,2] and nVT1[0,1,2]).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - - eE- - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE1 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE0
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            nV
*
*                   eDT0 += {nVT00, nVT10}
*                   eDT1 += {nVT01, nVT11}
*                   nDT = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMEdge *eE:		Shared edge.
*		WlzDVertex3 nPos:	Position of new vertex.
************************************************************************/
static WlzErrorNum WlzGMModelExtend2V1E1S3D(WlzGMModel *model, WlzGMEdge *eE,
					    WlzDVertex3 nPos)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[2];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT;
  WlzGMVertex	*nV;
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nV = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shell that's to be extended. */
    eShell = WlzGMEdgeGetShell(eE);
    /* Get geometry of all 3 verticies. */
    (void )WlzGMVertexGetG3D(eE->edgeT->vertexT->diskT->vertex, pos + 0);
    (void )WlzGMVertexGetG3D(eE->edgeT->opp->vertexT->diskT->vertex, pos + 1);
    pos[2] = nPos;
    /* New vertex. */
    (void )WlzGMVertexSetG3D(nV, nPos);
    nV->diskT = nDT;
    WlzGMModelAddVertex(model, nV);
    /* New vertex topology elements. */
    WlzGMVertexTAppend(eE->edgeT->vertexT, nVT0[0]);
    WlzGMVertexTAppend(eE->edgeT->vertexT, nVT1[0]);
    nVT0[0]->diskT = nVT1[0]->diskT = eE->edgeT->vertexT->diskT;
    WlzGMVertexTAppend(eE->edgeT->opp->vertexT, nVT0[1]);
    WlzGMVertexTAppend(eE->edgeT->opp->vertexT, nVT1[1]);
    nVT0[1]->diskT = nVT1[1]->diskT = eE->edgeT->opp->vertexT->diskT;
    nVT0[2]->next = nVT0[2]->prev = nVT1[2];
    nVT1[2]->next = nVT1[2]->prev = nVT0[2];
    nVT0[2]->diskT = nVT1[2]->diskT = nDT;
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology element. */
    nDT->next = nDT->prev = nDT;
    nDT->vertex = nV;
    nDT->vertexT = nVT0[2];
    /* New edges. */
    for(idx = 0; idx < 2; ++idx)
    {
      nE[idx]->edgeT = nET0[idx + 1];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    nET1[1]->edge = nET0[0]->edge = eE;
    nET1[2]->edge = nET0[1]->edge = nE[0];
    nET1[0]->edge = nET0[2]->edge = nE[1];
    /* Need to set radial edges for nET0[0] and nET1[1]. */
    WlzGMEdgeTInsertRadial(nET0[0]);
    WlzGMEdgeTInsertRadial(nET1[1]);
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update shell */
    (void )WlzGMShellUpdateG3D(eShell, pos[0]);
    (void )WlzGMShellUpdateG3D(eShell, pos[1]);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelExtend2V0E1S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Extends an existing shell within the model by adding a:
*		loop, 3 edges and 1 vertex (with topological elements).
*		Both the given 3D double precision points are known not
*		to be in the model.
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]), 1 vertex (nV) with geometry (pos0) and
*		6 vertex topology elements (nVT[0,1,2] and nVT1[0,1,2]).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - nE0 - - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE2 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE1
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            nV
*
*                   nDT0 = {nVT00, nVT10}
*                   nDT1 = {nVT01, nVT11}
*		    nDT2 = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMVertex *eV0:	First shared vertex.
*		WlzGMVertex *eV1:	Second shared vertex.
*		WlzDVertex3 nPos:	Position of new vertex.
************************************************************************/
static WlzErrorNum WlzGMModelExtend2V0E1S3D(WlzGMModel *model,
					    WlzGMVertex *eV0, WlzGMVertex *eV1,
					    WlzDVertex3 nPos)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertex	*nV;
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[2] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[2] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nV = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shell that's to be extended. */
    eShell = WlzGMVertexGetShell(eV0);
    /* Get geometry of the 2 existing verticies and the new vertex. */
    (void )WlzGMVertexGetG3D(eV0, pos + 0);
    (void )WlzGMVertexGetG3D(eV1, pos + 1);
    pos[2] = nPos;
    /* New vertex. */
    (void )WlzGMVertexSetG3D(nV, nPos);
    nV->diskT = nDT[2];
    WlzGMModelAddVertex(model, nV);
    /* New vertex topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->next = nVT0[idx]->prev = nVT1[idx];
      nVT1[idx]->next = nVT1[idx]->prev = nVT0[idx];
      nVT0[idx]->diskT = nVT1[idx]->diskT = nDT[idx];
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology elements. */
    WlzGMDiskTAppend(eV0->diskT, nDT[0]);
    nDT[0]->vertex = eV0;
    nDT[0]->vertexT = nVT0[0];
    WlzGMDiskTAppend(eV1->diskT, nDT[1]);
    nDT[1]->vertex = eV1;
    nDT[1]->vertexT = nVT0[1];
    nDT[2]->next = nDT[2]->prev = nDT[2];
    nDT[2]->vertex = nV;
    nDT[2]->vertexT = nVT0[2];
    /* New edges. */
    for(idx = 0; idx < 3; ++idx)
    {
      nE[idx]->edgeT = nET0[idx];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->edge = nE[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->edge = nE[pIdx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update shell */
    (void )WlzGMShellUpdateG3D(eShell, nPos);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelJoin2V0E0S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Joins two existing shells within the model by adding:
*		1 loop, 3 edges and 1 vertex (with topological elements).
*		Both the given 3D double precision points are known not
*		to be in the model.
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]), 1 vertex (nV) with geometry (pos0) and
*		6 vertex topology elements (nVT[0,1,2] and nVT1[0,1,2]).
*		Need to delete 1 shell (dShell).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - nE0 - - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE2 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE1
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            nV
*
*                   nDT0 = {nVT00, nVT10}
*                   nDT1 = {nVT01, nVT11}
*		    nDT2 = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMVertex *eV0:	First shared vertex.
*		WlzGMVertex *eV1:	Second shared vertex.
*		WlzDVertex3 nPos:	Position of new vertex.
************************************************************************/
static WlzErrorNum WlzGMModelJoin2V0E0S3D(WlzGMModel *model,
					  WlzGMVertex *eV0, WlzGMVertex *eV1,
					  WlzDVertex3 nPos)
{
  int		idx,
  		nIdx,
		pIdx;
  double	eVol,
  		dVol;
  WlzGMShell	*eShell,
  		*dShell,
		*tShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertex	*nV,
  		*tV;
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[2] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[2] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nV = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shell that's to be extended and the shell thats to be joined
     * with it (and then deleated). */
    eShell = WlzGMVertexGetShell(eV0);
    dShell = WlzGMVertexGetShell(eV1);
    /* Optimization: Get the shell geometries and compare their volume.
     * Then make sure that the shell with the smaller volume (and hopefully
     * the smaller number of loops) is the one that's to be deleted. */
    (void )WlzGMShellGetGBBV3D(eShell, &eVol);
    (void )WlzGMShellGetGBBV3D(dShell, &dVol);
    if(eVol < dVol)
    {
      tShell = eShell; eShell = dShell; dShell = tShell;
      tV = eV0; eV0 = eV1; eV1 = tV;
    }
    /* Get geometry of the 2 existing verticies and the new vertex. */
    (void )WlzGMVertexGetG3D(eV0, pos + 0);
    (void )WlzGMVertexGetG3D(eV1, pos + 1);
    pos[2] = nPos;
    /* New vertex. */
    (void )WlzGMVertexSetG3D(nV, nPos);
    nV->diskT = nDT[2];
    WlzGMModelAddVertex(model, nV);
    /* New vertex topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->next = nVT0[idx]->prev = nVT1[idx];
      nVT1[idx]->next = nVT1[idx]->prev = nVT0[idx];
      nVT0[idx]->diskT = nVT1[idx]->diskT = nDT[idx];
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology elements. */
    WlzGMDiskTAppend(eV0->diskT, nDT[0]);
    nDT[0]->vertex = eV0;
    nDT[0]->vertexT = nVT0[0];
    WlzGMDiskTAppend(eV1->diskT, nDT[1]);
    nDT[1]->vertex = eV1;
    nDT[1]->vertexT = nVT0[1];
    nDT[2]->next = nDT[2]->prev = nDT[2];
    nDT[2]->vertex = nV;
    nDT[2]->vertexT = nVT0[2];
    /* New edges. */
    for(idx = 0; idx < 3; ++idx)
    {
      nE[idx]->edgeT = nET0[idx];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->edge = nE[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->edge = nE[pIdx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update the extended shell and delete the joined shell. */
    (void )WlzGMShellUpdateG3D(eShell, nPos);
    if(model->child->idx == dShell->idx)
    {
      model->child = eShell;
    }
    WlzGMShellJoinAndUnlink(eShell, dShell);
    WlzGMModelFreeS(model, dShell);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelJoin3V0E3S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Joins three existing shells within the model by adding:
*		1 loop and 3 edges (with topological elements).
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*		Need to delete 2 shells (dShell[0,1]).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - nE2 - - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE2 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE1
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            eV2
*
*                   nDT0 = {nVT00, nVT10}
*                   nDT1 = {nVT01, nVT11}
*		    nDT2 = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMVertex *gEV:	The three given shared verticies.
************************************************************************/
static WlzErrorNum WlzGMModelJoin3V0E3S3D(WlzGMModel *model, WlzGMVertex **gEV)
{
  int		idx,
  		nIdx,
		pIdx;
  double	eVol;
  double	dVol[2];
  WlzGMShell	*eShell,
  		*tShell;
  WlzGMShell	*dShell[2];
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertex	*tV;
  WlzGMVertex	*eV[3];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  eV[0] = gEV[0]; eV[1] = gEV[1]; eV[2] = gEV[2];
  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[2] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[2] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shell that's to be extended and the shells that are to be joined
     * with it (and then deleated). */
    eShell = WlzGMVertexGetShell(eV[0]);
    dShell[0] = WlzGMVertexGetShell(eV[1]);
    dShell[1] = WlzGMVertexGetShell(eV[2]);
    /* Optimization: Get the shell geometries and compare their volume.
     * Then make sure that the shell with the greatest volume (and hopefully
     * the greatest number of loops) is the one that's to be kept. */
    (void )WlzGMShellGetGBBV3D(eShell, &eVol);
    (void )WlzGMShellGetGBBV3D(dShell[0], &dVol[0]);
    (void )WlzGMShellGetGBBV3D(dShell[1], &dVol[1]);
    if((eVol < dVol[0]) || (eVol < dVol[1]))
    {
      if(dVol[0] > dVol[1])
      {
        tShell = eShell; eShell = dShell[0]; dShell[0] = tShell;
	tV = eV[0]; eV[0] = eV[1]; eV[1] = tV;
      }
      else
      {
        tShell = eShell; eShell = dShell[1]; dShell[1] = tShell;
	tV = eV[0]; eV[0] = eV[2]; eV[2] = tV;
      }
    }
    /* Get geometry of the 3 existing verticies. */
    for(idx = 0; idx < 3; ++idx)
    {
      (void )WlzGMVertexGetG3D(*(eV + idx), pos + idx);
    }
    /* New vertex topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->next = nVT0[idx]->prev = nVT1[idx];
      nVT1[idx]->next = nVT1[idx]->prev = nVT0[idx];
      nVT0[idx]->diskT = nVT1[idx]->diskT = nDT[idx];
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      WlzGMDiskTAppend(eV[idx]->diskT, nDT[idx]);
      nDT[idx]->vertex = eV[idx];
      nDT[idx]->vertexT = nVT0[idx];
    }
    /* New edges. */
    for(idx = 0; idx < 3; ++idx)
    {
      nE[idx]->edgeT = nET0[idx];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->edge = nE[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->edge = nE[pIdx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update the extended shell and delete the joined shells. */
    for(idx = 0; idx < 2; ++idx)
    {
      if(model->child->idx == dShell[idx]->idx)
      {
	model->child = eShell;
      }
      WlzGMShellJoinAndUnlink(eShell, dShell[idx]);
      WlzGMModelFreeS(model, dShell[idx]);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelJoin3V0E2S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Joins two existing shells within the model by adding:
*		1 loop and 3 edges (with topological elements).
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*		Need to delete 1 shell (dShell).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - nE0 - - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE2 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE1
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            eV2
*
*                   nDT0 = {nVT00, nVT10}
*                   nDT1 = {nVT01, nVT11}
*		    nDT2 = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMVertex *eV:	The three given shared verticies.
************************************************************************/
static WlzErrorNum WlzGMModelJoin3V0E2S3D(WlzGMModel *model, WlzGMVertex **eV)
{
  int		idx,
  		nIdx,
		pIdx;
  double	sVol[3];
  WlzGMShell	*dShell,
  		*eShell;
  WlzGMShell	*tShell[3];
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertex   *ev[3];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[2] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[2] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shells and find which is to be extended and which is to be joined
     * with it (and then deleated). Know that 2 of the shells are the same.
     * Compute volumes of the shells bounding boxes, compare the volumes,
     * delete the shell with the smallest volume and keep the shell with
     * the largest volume. This is an optimazation, hopefuly the shell with
     * greatest volume will have the greatest number of loops. */
    for(idx = 0; idx < 3; ++idx)
    {
      tShell[idx] = WlzGMVertexGetShell(*(eV + idx));
      (void )WlzGMShellGetGBBV3D(tShell[idx],sVol + idx);
    }
    if(sVol[0] > sVol[1])
    {
      if(sVol[0] > sVol[2])
      {
        idx = 0;
      }
      else
      {
        idx = 2;
      }
    }
    else
    {
      if(sVol[1] > sVol[2])
      {
        idx = 1;
      }
      else
      {
        idx = 2;
      }
    }
    eShell = tShell[idx];
    dShell = (tShell[idx]->idx == tShell[(idx + 1) % 3]->idx)?
             tShell[(idx + 2) % 3]: tShell[(idx + 1) % 3];
    if(tShell[0]->idx == tShell[1]->idx)
    {
      eShell = tShell[0];
      dShell = tShell[2];
    }
    else if(tShell[1]->idx == tShell[2]->idx)
    {
      eShell = tShell[1];
      dShell = tShell[0];
    }
    else /* tShell[2]->idx == tShell[0]->idx) */
    {
      eShell = tShell[2];
      dShell = tShell[1];
    }
    /* New vertex topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->next = nVT0[idx]->prev = nVT1[idx];
      nVT1[idx]->next = nVT1[idx]->prev = nVT0[idx];
      nVT0[idx]->diskT = nVT1[idx]->diskT = nDT[idx];
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      WlzGMDiskTAppend(eV[idx]->diskT, nDT[idx]);
      nDT[idx]->vertex = eV[idx];
      nDT[idx]->vertexT = nVT0[idx];
    }
    /* New edges. */
    for(idx = 0; idx < 3; ++idx)
    {
      nE[idx]->edgeT = nET0[idx];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->edge = nE[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->edge = nE[pIdx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update the extended shell and delete the joined shell. */
    if(model->child->idx == dShell->idx)
    {
      model->child = eShell;
    }
    WlzGMShellJoinAndUnlink(eShell, dShell);
    WlzGMModelFreeS(model, dShell);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelExtend3V0E1S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Extends an existing shell within the model by adding:
*		1 loop and 3 edges (with topological elements).
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - nE0 - - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE2 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE1
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            eV2
*
*                   nDT0 = {nVT00, nVT10}
*                   nDT1 = {nVT01, nVT11}
*		    nDT2 = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMVertex *eV:	Thre three shared vertex.
************************************************************************/
static WlzErrorNum WlzGMModelExtend3V0E1S3D(WlzGMModel *model,
					    WlzGMVertex **eV)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[2] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[2] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shell. */
    eShell = WlzGMVertexGetShell(*eV);
    /* Get geometry of the 3 existing verticies. */
    for(idx = 0; idx < 3; ++idx)
    {
      (void )WlzGMVertexGetG3D(*(eV + idx), pos + idx);
    }
    /* New vertex topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->next = nVT0[idx]->prev = nVT1[idx];
      nVT1[idx]->next = nVT1[idx]->prev = nVT0[idx];
      nVT0[idx]->diskT = nVT1[idx]->diskT = nDT[idx];
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      WlzGMDiskTAppend(eV[idx]->diskT, nDT[idx]);
      nDT[idx]->vertex = eV[idx];
      nDT[idx]->vertexT = nVT0[idx];
    }
    /* New edges. */
    for(idx = 0; idx < 3; ++idx)
    {
      nE[idx]->edgeT = nET0[idx];
    }
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->edge = nE[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->edge = nE[pIdx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* No change to shell. */
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelExtend3V1E1S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Extends an existing shell within the model by adding:
*		1 loop and 2 edges (with topological elements).
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 2 edges (nE[0,1]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - eE  - - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE1 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE0
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            eV2 = sV
*
*		    nDT = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMEdge *eE:		The shared edge.
*		WlzGMVertex *sV:	The shared vertex (not on the
*					shared edge).
************************************************************************/
static WlzErrorNum WlzGMModelExtend3V1E1S3D(WlzGMModel *model,
					    WlzGMEdge *eE, WlzGMVertex *sV)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[2];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT;
  WlzGMVertex	*eV[3];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shell that's to be extended. */
    eShell = WlzGMEdgeGetShell(eE);
    /* Get verticies. */
    eV[0] = eE->edgeT->vertexT->diskT->vertex;
    eV[1] = eE->edgeT->opp->vertexT->diskT->vertex;
    eV[2] = sV;
    /* Get geometry of all 3 verticies. */
    for(idx = 0; idx < 3; ++idx)
    {
      (void )WlzGMVertexGetG3D(*(eV + idx), pos + idx);
    }
    /* New vertex topology elements. */
    WlzGMVertexTAppend(eE->edgeT->vertexT, nVT0[0]);
    WlzGMVertexTAppend(eE->edgeT->vertexT, nVT1[0]);
    nVT0[0]->diskT = nVT1[0]->diskT = eE->edgeT->vertexT->diskT;
    WlzGMVertexTAppend(eE->edgeT->opp->vertexT, nVT0[1]);
    WlzGMVertexTAppend(eE->edgeT->opp->vertexT, nVT1[1]);
    nVT0[1]->diskT = nVT1[1]->diskT = eE->edgeT->opp->vertexT->diskT;
    nVT0[2]->next = nVT0[2]->prev = nVT1[2];
    nVT1[2]->next = nVT1[2]->prev = nVT0[2];
    nVT0[2]->diskT = nVT1[2]->diskT = nDT;
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology element. */
    nDT->next = nDT->prev = nDT;
    nDT->vertex = eV[2];
    nDT->vertexT = nVT0[2];
    /* New edges. */
    nE[0]->edgeT = nET0[1];
    nE[1]->edgeT = nET0[2];
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    nET1[1]->edge = nET0[0]->edge = eE;
    nET1[2]->edge = nET0[1]->edge = nE[0];
    nET1[0]->edge = nET0[2]->edge = nE[1];
    /* Need to set radial edges for nET0[0] and nET1[1]. */
    WlzGMEdgeTInsertRadial(nET0[0]);
    WlzGMEdgeTInsertRadial(nET1[1]);
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update shell */
    (void )WlzGMShellUpdateG3D(eShell, pos[2]);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelJoin3V1E2S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Joins two existing shells within the model by adding:
*		1 loop and 2 edges (with topological elements).
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 2 edges (nE[0,1]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*		Need to delete 1 shell (dShell).
*
*                           nVT00        nET00
*                         O ------------------------------->
*                   eV0 @- - - - - - - - eE  - - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  nE1 ------\ \                       / / /
*                           \   \ nET10               /   /
*                            \ \ \                   / /---- nE0
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @
*                                            eV2 = sV
*
*                   nDT = {nVT02, nVT12}
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMEdge *eE:		The shared edge.
*		WlzGMVertex *sV:	The shared vertex (not on the
*					shared edge).
************************************************************************/
static WlzErrorNum WlzGMModelJoin3V1E2S3D(WlzGMModel *model,
				          WlzGMEdge *eE, WlzGMVertex *sV)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*dShell,
  		*eShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[2];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT;
  WlzGMVertex	*eV[3];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE[0] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nE[1] = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shells that are to be extended and deleted. */
    eShell = WlzGMEdgeGetShell(eE);
    dShell = WlzGMVertexGetShell(sV);
    /* Get verticies. */
    eV[0] = eE->edgeT->vertexT->diskT->vertex;
    eV[1] = eE->edgeT->opp->vertexT->diskT->vertex;
    eV[2] = sV;
    /* Get geometry of all 3 verticies. */
    for(idx = 0; idx < 3; ++idx)
    {
      (void )WlzGMVertexGetG3D(*(eV + idx), pos + idx);
    }
    /* New vertex topology elements. */
    WlzGMVertexTAppend(eE->edgeT->vertexT, nVT0[0]);
    WlzGMVertexTAppend(eE->edgeT->vertexT, nVT1[0]);
    nVT0[0]->diskT = nVT1[0]->diskT = eE->edgeT->vertexT->diskT;
    WlzGMVertexTAppend(eE->edgeT->opp->vertexT, nVT0[1]);
    WlzGMVertexTAppend(eE->edgeT->opp->vertexT, nVT1[1]);
    nVT0[1]->diskT = nVT1[1]->diskT = eE->edgeT->opp->vertexT->diskT;
    nVT0[2]->next = nVT0[2]->prev = nVT1[2];
    nVT1[2]->next = nVT1[2]->prev = nVT0[2];
    nVT0[2]->diskT = nVT1[2]->diskT = nDT;
    for(idx = 0; idx < 3; ++idx)
    {
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New disk topology element. */
    nDT->next = nDT->prev = nDT;
    nDT->vertex = eV[2];
    nDT->vertexT = nVT0[2];
    /* New edges. */
    nE[0]->edgeT = nET0[1];
    nE[1]->edgeT = nET0[2];
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    nET1[1]->edge = nET0[0]->edge = eE;
    nET1[2]->edge = nET0[1]->edge = nE[0];
    nET1[0]->edge = nET0[2]->edge = nE[1];
    /* Need to set radial edges for nET0[0], nET1[1]. */
    WlzGMEdgeTInsertRadial(nET0[0]);
    WlzGMEdgeTInsertRadial(nET1[1]);
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update extended shell and delete the joined shell. */
    (void )WlzGMShellUpdateG3D(eShell, pos[2]);
    if(model->child->idx == dShell->idx)
    {
      model->child = eShell;
    }
    WlzGMShellJoinAndUnlink(eShell, dShell);
    WlzGMModelFreeS(model, dShell);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelExtend3V2E1S3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Extends a shell within the model by adding:
*		1 loop and 1 edge (with topological elements).
*		Need to create: 1 loop (nL), 2 loop topology elements
*		(nLT[0,1]), 1 edge (nE), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]) and 6 vertex topology elements
*		(nVT[0,1,2] and nVT1[0,1,2]).
*
*                           nVT02        nET02
*                         O ------------------------------->
*                   eV2 @- - - - - - - - nE  - - - - - - - - -@ eV0
*                        \   <--------------------------- O  /
*                       ^   O nVT12      nET10     nVT10  ^   O nVT00
*                        \ \ \                           / / /
*                         \   \                         /   /
*                  eE1 ------\ \                       / / /
*                           \   \ nET12               /   /
*                            \ \ \                   / /---- eE0
*                             \   \                 /   /
*                        nET01 \ \ \         nET11 / / / nET00
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT11   /   /
*                                  \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT12 O   v O   v
*                                        \   /
*                                          @
*                                            eV1 
*
*                   nLT0 = {nET00, nET01, nET02}
*                   nLT1 = {nET10, nET11, nET12}
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMEdge *eE:		The shared edges.
************************************************************************/
static WlzErrorNum WlzGMModelExtend3V2E1S3D(WlzGMModel *model,
					    WlzGMEdge *eE0, WlzGMEdge *eE1)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*eDT[2];
  WlzGMVertex	*eV[3];
  WlzGMVertexT	*eVT[3],
  		*nVT0[3],
  		*nVT1[3];
  WlzDVertex3	pos[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET0[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1[2] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nVT0[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT0[2] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[1] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1[2] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Get shell that's to be extended. */
    eShell = WlzGMEdgeGetShell(eE0);
    /* Get verticies. */
    eV[1] = WlzGMEdgeCommonVertex(eE0, eE1);
    eV[0] = (eVT[0] = eE0->edgeT->vertexT)->diskT->vertex;
    if(eV[0]->idx == eV[1]->idx)
    {
      eV[0] = (eVT[0] = eE0->edgeT->opp->vertexT)->diskT->vertex;
      eVT[1] = eE0->edgeT->vertexT;
    }
    else
    {
      eVT[1] = eE0->edgeT->opp->vertexT;
    }
    eV[2] = (eVT[2] = eE1->edgeT->vertexT)->diskT->vertex;
    if(eV[2]->idx == eV[1]->idx)
    {
      eV[2] = (eVT[2] = eE1->edgeT->opp->vertexT)->diskT->vertex;
    }
    /* Get geometry of all 3 verticies. */
    for(idx = 0; idx < 3; ++idx)
    {
      (void )WlzGMVertexGetG3D(*(eV + idx), pos + idx);
    }
    /* New vertex topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      WlzGMVertexTAppend(eVT[idx], nVT0[idx]);
      WlzGMVertexTAppend(eVT[idx], nVT1[idx]);
      nVT0[idx]->diskT = nVT1[idx]->diskT = eVT[idx]->diskT;
      nVT0[idx]->parent = nET0[idx];
      nVT1[idx]->parent = nET1[idx];
    }
    /* New edge. */
    nE->edgeT = nET0[2];
    /* New edge topology elements */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nET0[idx]->next = nET0[nIdx];
      nET0[idx]->prev = nET0[pIdx];
      nET0[idx]->opp = nET1[nIdx];
      nET0[idx]->rad = nET0[idx];
      nET0[idx]->vertexT = nVT0[idx];
      nET0[idx]->parent = nLT[0];
      nET1[idx]->next = nET1[pIdx];   /* Previous because reverse direction. */
      nET1[idx]->prev = nET1[nIdx];       /* Next because reverse direction. */
      nET1[idx]->opp = nET0[pIdx];
      nET1[idx]->rad = nET1[idx];
      nET1[idx]->vertexT = nVT1[idx];
      nET1[idx]->parent = nLT[1];
    }
    nET1[1]->edge = nET0[0]->edge = eE0;
    nET1[2]->edge = nET0[1]->edge = eE1;
    nET1[0]->edge = nET0[2]->edge = nE;
    /* Need to set radial edges. */
    WlzGMEdgeTInsertRadial(nET0[0]);
    WlzGMEdgeTInsertRadial(nET1[1]);
    WlzGMEdgeTInsertRadial(nET0[1]);
    WlzGMEdgeTInsertRadial(nET1[2]);
    /* New loop */
    nL->loopT = nLT[0];
    /* New loop topology elements. */
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->loop = nL;
    nLT[1]->loop = nL;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* No change to shell */
    /* No new disk topology elements, but may need to delete one. */
    eDT[0] = nET0[1]->rad->vertexT->diskT;
    eDT[1] = nET1[1]->rad->vertexT->diskT;
    if(eDT[0]->idx != eDT[1]->idx)
    {
      eVT[0] = eVT[1] = eDT[1]->vertexT;
      do
      {
        eVT[2] = eVT[1]->next;
	WlzGMVertexTAppend(eDT[0]->vertexT, eVT[1]);
	eVT[1]->diskT = eDT[0];
	eVT[1] = eVT[2];
      } while(eVT[1]->idx != eVT[0]->idx);
      if(eV[1]->diskT->idx == eDT[1]->idx)
      {
        eV[1]->diskT = eDT[0];
      }
      WlzGMModelFreeDT(model, eDT[1]);
    }
  }
  return(errNum);
}


/************************************************************************
* Function:	WlzGMModelConstructNewS2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a new shell, loop and edge (with topological
*		elements) and adds them to the model. Given a pair
*		of 2D double precision edge end points which are known
*		not to be in the model.
*		Neither of the verticies already exists within the shell.
*		Need to create a new shell with: 1 loop (nL), 1 loop
*		topology element (nLT), 1 edge (nE) 2 edge
*		topology elements (nET[0,1]), 2 disk topology elements,
*		(nDT[0,1]), 2 verticies (nV[0,1] with geometry
*		pos[0,1]) and 2 vertex topology elements (nVT[0,1]).
*
*		     / nV0
*		    /
*		   |  / nVT0    / nET0          / nE
*		   | /         /               /
*		   | O -----------------------/---->
*		    @- - - - - - - - - - - - - - - -@
*		     <---------------------------- O|
*		                      /           / |
*		                nET1 /      nVT1 /  |
*		                                    /
*		                               nV1 /
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzDVertex2 *pos:	Ptr to 1st then 2nd point.
************************************************************************/
static WlzErrorNum WlzGMModelConstructNewS2D(WlzGMModel *model,
					     WlzDVertex2 *pos)
{
  int		idx;
  WlzGMShell	*eShell,
  		*nShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET[2];
  WlzGMDiskT	*nDT[2];
  WlzGMVertex	*nV[2];
  WlzGMVertexT	*nVT[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nShell = WlzGMModelNewS(model, &errNum)) != NULL) &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET[0] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET[1] = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT[0] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nDT[1] = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nV[0] = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nV[1] = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT[0] = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT[1] = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    for(idx = 0; idx < 2; ++idx)
    {
      /* New vertcies */
      nV[idx]->diskT = nDT[idx];
      (void )WlzGMVertexSetG2D(nV[idx], *(pos + idx));
      WlzGMModelAddVertex(model, nV[idx]);
      /* New vertex topology elements */
      nVT[idx]->next = nVT[idx]->prev = nVT[idx];
      nVT[idx]->diskT = nDT[idx];
      nVT[idx]->parent = nET[idx];
      /* New disk topology elements */
      nDT[idx]->next = nDT[idx]->prev = nDT[idx];
      nDT[idx]->vertex = nV[idx];
      nDT[idx]->vertexT = nVT[idx];
    }
    /* New Edges and their topology elements */
    nE->edgeT = nET[0];
    nET[0]->edge = nET[1]->edge = nE;
    nET[0]->rad = nET[0]->next = nET[0]->prev = nET[0]->opp = nET[1];
    nET[0]->vertexT = nVT[0];
    nET[0]->parent = nLT;
    nET[1]->rad = nET[1]->next = nET[1]->prev = nET[1]->opp = nET[0];
    nET[1]->vertexT = nVT[1];
    nET[1]->parent = nLT;
    /* New loop and loop topology elements */
    nL->loopT = nLT;
    nLT->next = nLT->prev = nLT->opp = nLT;
    nLT->loop = nL;
    nLT->edgeT = nET[0];
    nLT->parent = nShell;
    /* New shell in model */
    (void )WlzGMShellSetG2D(nShell, 2, pos);
    nShell->child = nLT;
    nShell->parent = model;
    nShell->next = nShell->prev = nShell;
    /* Is this the first shell of the model? */
    if((eShell = model->child) == NULL)
    {
      /* First shell in model */
      model->child = nShell;
    }
    else
    {
      /* Append new shell onto end of model's list of shells. */
      WlzGMShellAppend(eShell, nShell);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelExtendL2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a new edge (with topological elements) and
*		adds them to the model. Given an existing edge topology
*		element at one end of the new edge and a 2D double
*		precision end points at the other end of the new edge,
*		where the new end point is known not to be in the model.
*		Only one vertex already exists within the model so no new
*		shell or loop is required, instead an existing loop is
*		extended using: 1 edge (nE), 2 edge topology elements
*		(nET0, nET1), 1 disk topology element (nDT), 1 vertex
*		(nV0 with geometry (nPos) and 2 vertex topology elements
*		(nVT0, nVT1).
*
*		    / nV0                    nE
*		   /                             /
*		   |   /nVT0    / nET0          /         /eET0
*		   |  /        /               /         /
*		   |  O ----------------------/--->   O----
*		    @- - - - - - - - - - - - - - - -@- - - -
*		     <--------------------------- O |<------
*		                      /          /  |    \
*		                nET1 /      nVT1/   |     \eET1
*		                                   /
*		                               nV1/
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMEdgeT *eET0:	Next edge topology element for the
*					new edge topology element directed
*					away from the new end point.
*		WlzDVertex2 nPos:	The new end point.
************************************************************************/
static WlzErrorNum WlzGMModelExtendL2D(WlzGMModel *model, WlzGMEdgeT *eET0,
				       WlzDVertex2 nPos)
{
  WlzGMVertex	*nV;
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMDiskT	*nDT;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*eET1,
		*nET0,
  		*nET1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model && eET0 &&
     ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nDT = WlzGMModelNewDT(model, &errNum)) != NULL) &&
     ((nV = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1 = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    eET1 = eET0->prev;
    /* New vertex */
    nV->diskT = nDT;
    (void )WlzGMVertexSetG2D(nV, nPos);
    WlzGMModelAddVertex(model, nV);
    /* Set vertex topology elements */
    nVT0->prev = nVT0->next = nVT0;
    nVT0->diskT = nDT;
    nVT0->parent = nET0;
    WlzGMVertexTAppend(eET0->vertexT, nVT1);
    nVT1->diskT = eET0->vertexT->diskT;
    nVT1->parent = nET1;
    /* New disk topology element */
    nDT->next = nDT->prev = nDT;
    nDT->vertex = nV;
    nDT->vertexT = nVT0;
    /* Edges and their topology elements */
    nE->edgeT = nET0;
    nET0->prev = nET1;
    nET0->next = eET0;
    eET0->prev = nET0;
    nET1->prev = eET1;
    eET1->next = nET1;
    nET1->next = nET0;
    nET0->opp = nET1;
    nET0->rad = nET0;
    nET0->edge = nE;
    nET0->vertexT = nVT0;
    nET0->parent = eET0->parent;
    nET1->opp = nET0;
    nET1->rad = nET1;
    nET1->edge = nE;
    nET1->vertexT = nVT1;
    nET1->parent = eET1->parent;
    /* Update the shell's geometry. */
    (void )WlzGMShellUpdateG2D(nET0->parent->parent, nPos);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelConstructSplitL2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	First of all checks that the edge segment doesn't
*		already exist! If it doesn't then a new edge is
*		constructed which splits and existing loop
*		within the model. Given the loop topology element
*		to be split and the existing edge topology elements at
*		new edges end points.
*		Split a loop within a shell by adding a new edge between the
*		two matched verticies.
*		Create: 1 edge (nE), 2 edge topology elements (nET0,
*		nET1), 2 vertex topology elements (nVT0, nVT1), a
*		loop (nL) and a loop topology element (nLT).
*
*		      O ------------->     O --------------->
*		    @- - - - - - - - - -@- - - - - - - - - - -@
*		    |  <------------ O  |  <--------------- O | O
*		  ^   O        /      ^   O                 ^
*		  | |    eET0 /       | |  \nVT1            | | |
*		  |   |               |   |                 |   |
*		  | | |        nET0 \ | | |/nE              | | |
*		  |   |              \|   /                 |   |
*		  | | |    eL         | |/| /nET1   nL      | | |
*		  |   |               |   |/                |   |
*		  | | |         nVT0 \  | |     / eET1        | |
*		  |   V               O   V    /            O   V
*		    |  O ------------>  |  O -------------->  |
*		    @- - - - - - - - - -@- - - - - - - - - - -@
*		      <-------------- O    <--------------- O
*
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMEdgeT *eET0:	Next edge topology element for the
*					first end point.
*		WlzGMEdgeT *eET1:	Next edge topology element for the
*					other end point.
************************************************************************/
static WlzErrorNum WlzGMModelConstructSplitL2D(WlzGMModel *model,
					       WlzGMEdgeT *eET0,
					       WlzGMEdgeT *eET1)
{
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0,
  		*nET1;
  WlzGMLoop	*nL;
  WlzGMLoopT	*eLT,
  		*nLT;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check for duplicate edge. */
  if(eET0->prev->idx != eET1->prev->opp->idx)
  {
    if(model && eET0 &&
       ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
       ((nET0 = WlzGMModelNewET(model, &errNum)) != NULL) &&
       ((nET1 = WlzGMModelNewET(model, &errNum)) != NULL) &&
       ((nVT0 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
       ((nVT1 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
       ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
       ((nLT = WlzGMModelNewLT(model, &errNum)) != NULL))
    {
      /* Vertex topology elements */
      WlzGMVertexTAppend(eET1->vertexT, nVT0);
      nVT0->diskT = eET1->vertexT->diskT;
      nVT0->parent = nET0;
      WlzGMVertexTAppend(eET0->vertexT, nVT1);
      nVT1->diskT = eET0->vertexT->diskT;
      nVT1->parent = nET1;
      /* Edges and their topology elements */
      nE->edgeT = nET0;
      nET0->next = eET0;
      nET0->prev = eET1->prev;
      nET0->opp = nET1;
      nET0->rad = nET0;
      nET0->edge = nE;
      nET0->vertexT = nVT0;
      nET0->parent = eET0->parent;
      nET1->next = eET1;
      nET1->prev = eET0->prev;
      nET1->opp = nET0;
      nET1->rad = nET1;
      nET1->edge = nE;
      nET1->vertexT = nVT1;
      nET1->parent = eET1->parent;
      nET0->next->prev = nET0;
      nET0->prev->next = nET0;
      nET1->next->prev = nET1;
      nET1->prev->next = nET1;
      /* Loop and loop topology element */
      eLT = eET0->parent;
      eLT->edgeT = nET0;
      nL->loopT = nLT;
      WlzGMLoopTAppend(eET0->parent, nLT);	       /* Add loopT to shell */
      nLT->opp = nLT;
      nLT->loop = nL;
      nLT->parent = eLT->parent;
      /* Set/update the loops by walking around the edge topology elements
       * setting their parents. */
      eLT->edgeT = nET0;
      nLT->edgeT = nET1;
      WlzGMLoopTSetT(eLT);
      WlzGMLoopTSetT(nLT);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelJoinL2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a new edge which joins two (different)
*		existing loops within the model.
*		Join two loops within two possibly different shells by
*		adding a new edge between the two matched verticies.
*		Create: 1 edge (nE), 2 edge topology elements (nET0,
*		nET1) and 2 vertex topology elements (nVT0, nVT1).
*		Delete: 1 loop (dL), 1 loop topology element (dLT)
*		and 1 shell (dShell). If edges were (they're not in 2D)
*		allowed to cross then it would be possible to have a
*		pair of verticies within the same shell but NOT within
*		a common loop. Check for this disallowed case and return
*		an error if it's found.
*
*	                                 eV1\
*	                                     \
*	                               nET0\  \
*	                       nE\          \  \          /eET1
*	                          \   /nVT0  \  \        /
*	             O -------->   \ O --------->\   O -------->
*	          @- - - - - - - -@- - - - - - - -@- - - - - - - -@
*	             <-------- O  \  <--------- O    <-------- O
*	                  /        \       \    \
*	             eET0/          \eV0    \    \nVT1
*	                                     \
*	                                      \nET1
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzGMEdgeT *eET0:	Next edge topology element for the
*					first end point.
*		WlzGMEdgeT *eET1:	Next edge topology element for the
*					other end point.
************************************************************************/
static WlzErrorNum WlzGMModelJoinL2D(WlzGMModel *model,
				     WlzGMEdgeT *eET0, WlzGMEdgeT *eET1)
{
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0,
  		*nET1;
  WlzGMLoop	*dL;
  WlzGMLoopT	*eLT,
  		*dLT;
  WlzGMShell	*dS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model && eET0 &&
     ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nVT0 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1 = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* Vertex topology elements */
    WlzGMVertexTAppend(eET0->vertexT, nVT0);
    nVT0->diskT = eET0->vertexT->diskT; /* eV0 */
    nVT0->parent = nET0;
    WlzGMVertexTAppend(eET1->vertexT, nVT1);
    nVT1->diskT = eET1->vertexT->diskT; /* eV1 */
    nVT1->parent = nET1;
    /* Edges and their topology elements */
    nE->edgeT = nET0;
    nET0->next = eET1;
    nET0->prev = eET0->prev;
    nET0->opp = nET1;
    nET0->rad = nET0;
    nET0->edge = nE;
    nET0->vertexT = nVT0;
    nET0->parent = eET0->parent;
    nET1->next = eET0;
    nET1->prev = eET1->prev;
    nET1->opp = nET0;
    nET1->rad = nET1;
    nET1->edge = nE;
    nET1->vertexT = nVT1;
    nET1->parent = eET0->parent;
    nET0->next->prev = nET0;
    nET0->prev->next = nET0;
    nET1->next->prev = nET1;
    nET1->prev->next = nET1;
    /* Loop and loop topology elements */
    eLT = eET0->parent;
    dLT = eET1->parent;
    dL = dLT->loop;
    dS = dLT->parent;
    /* Look to see if loops share the same shell. */
    if((dS = dLT->parent)->idx == eLT->parent->idx)
    {
      dS = NULL;
    }
    else
    {
      /* Merge the two shell's geometries. */
      (void )WlzGMShellMergeG(eLT->parent, dS);
    }
    /* Set/update the retained loop by walking around the edge topology
     * elements setting their parents. */
    WlzGMLoopTSetT(eLT);
    /* Unlink and free the deleted loop and loop topology element. */
    WlzGMLoopTUnlink(dLT);
    (void )WlzGMModelFreeLT(model, dLT);
    (void )WlzGMModelFreeL(model, dL);
    /* If deleted loop had a different parent shell then unlink it and
     * free the deleted loops parent shell if it's childless. */
    if(dS)
    {
      WlzGMShellUnlink(dS);
      (void )WlzGMModelFreeS(model, dS);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelConstructSimplex2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a 2D simplex (edge) defined by two double
*		precision end points. Either of the two points may
*		already exist within the model.
*		See WlzGMShellMatchVtxG2D() for the meaning of the
*		backwards, forwards and distance search parameters.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzDVertex2 *pos:	Pointer to first then second
*					positions.
************************************************************************/
WlzErrorNum	WlzGMModelConstructSimplex2D(WlzGMModel *model,
					     WlzDVertex2 *pos)
{
  unsigned int	matchCode;
  WlzGMEdgeT	*matchEdgeT[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

#ifdef WLZ_GEOMODEL_DEBUG
  (void )fprintf(stderr, "%g %g %g %g\n",
		 (pos + 0)->vtX, (pos + 0)->vtY,
		 (pos + 1)->vtX, (pos + 1)->vtY);
  (void )fflush(stderr);
#endif /* WLZ_GEOMODEL_DEBUG */
  WlzGMModelMatchEdgeTG2D(model, matchEdgeT, pos);
  matchCode = ((matchEdgeT[1] != NULL) << 1) | (matchEdgeT[0] != NULL);
  switch(matchCode)
  {
    case 0:
      /* Edge makes a new shell. */
      errNum = WlzGMModelConstructNewS2D(model, pos);
      break;
    case 1:
      /* Edge extends a loop and shell. */
      errNum = WlzGMModelExtendL2D(model, matchEdgeT[0], *(pos + 1));
      break;
    case 2:
      /* Edge extends a loop and shell. */
      errNum = WlzGMModelExtendL2D(model, matchEdgeT[1], *(pos + 0));
      break;
    case 3:
      /* Both verticies already exist within the model.
         Don't insert this simplex if it is already in the model. */
      if(matchEdgeT[0]->edge->idx != matchEdgeT[1]->edge->idx)
      {
	/* Both verticies already exist within the model so adding another
	 * edge between the two verticies will change the number of loops
	 * and shells. */
	if(WlzGMEdgeTCommonLoopT(matchEdgeT[0], matchEdgeT[1]) != NULL)
	{
	  /* Both of the verticies lie on a common loop then that loop needs
	   * to be SPLIT by a new edge. */
	  errNum = WlzGMModelConstructSplitL2D(model,
						matchEdgeT[0], matchEdgeT[1]);
	}
	else
	{
	  /* The two verticies do NOT share a common loop so loops will be
	   * JOINed by a new edge, destroying a loop.  */
	  errNum = WlzGMModelJoinL2D(model, matchEdgeT[0], matchEdgeT[1]);
	  
	}
      }
      break;
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelConstructSimplex3D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a 3D simplex (triangle) defined by three
*		double precision verticies, any of which may already
*		exist within the model.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model to add the segment to.
*		WlzDVertex3 *pos:	Pointer to triangle verticies,
************************************************************************/
WlzErrorNum	WlzGMModelConstructSimplex3D(WlzGMModel *model,
				       	     WlzDVertex3 *pos)
{
  int		idx0,
  		idx1,
		idx2,
		edgeCnt,
		shellCnt;
  WlzGMEdge	*cE[3];
  unsigned int	matchCode;
  WlzGMVertex	*matchV[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int 	bitCntTb[8] =	     /* Number of verticies matched for a given
  				      * match code index */
  {
    0, /* 0 */
    1, /* 1 */
    1, /* 2 */
    2, /* 3 */
    1, /* 4 */
    2, /* 5 */
    2, /* 6 */
    3  /* 7 */
  },
  		firstCtgBitTb[8] = /* First of contigous bits set, bits wrap
				    * around */
  {
    0, /* 0 */
    0, /* 1 */
    1, /* 2 */
    0, /* 3 */
    2, /* 4 */
    2, /* 5 */
    1, /* 6 */
    0  /* 7 */
  };
  const int 	matchVtxIdx[8][3] =  /* Indicies of the given simplex verticies
  				      * that have been matched for a given
				      * match code index */
  {
    {0, 0, 0}, /* 0 */
    {0, 1, 2}, /* 1 */
    {1, 2, 0}, /* 2 */
    {0, 1, 2}, /* 3 */
    {2, 0, 1}, /* 4 */
    {2, 0, 1}, /* 5 */
    {1, 2, 0}, /* 6 */
    {0, 1, 2}, /* 7 */
  };

  for(idx0 = 0; idx0 < 3; ++idx0)
  {
    matchV[idx0] = WlzGMModelMatchVertexG3D(model, *(pos + idx0));
  }
  matchCode = ((matchV[2] != NULL) << 2) |
  	      ((matchV[1] != NULL) << 1) | (matchV[0] != NULL);
  switch(bitCntTb[matchCode])
  {
    case 0:
      /* No verticies matched, construct a new shell. */
      errNum = WlzGMModelConstructNewS3D(model, pos);
      break;
    case 1:
      /* Single vertex matched, extend existing shell. */
      idx0 = matchVtxIdx[matchCode][0];
      idx1 = matchVtxIdx[matchCode][1];
      idx2 = matchVtxIdx[matchCode][2];
      errNum = WlzGMModelExtend1V0E1S3D(model, matchV[idx0],
				        *(pos + idx1), *(pos + idx2));
      break;
    case 2:
      /* Two verticies matched, check for existing edges and shells. */
      idx0 = matchVtxIdx[matchCode][0];
      idx1 = matchVtxIdx[matchCode][1];
      idx2 = matchVtxIdx[matchCode][2];
      if(WlzGMVertexCommonShell(matchV[idx0], matchV[idx1]) != NULL)
      {
        if((cE[0] = WlzGMVertexCommonEdge(matchV[idx0],
				          matchV[idx1])) != NULL)
	{
	  /* Matched 2 verticies connected by an edge. */
	  errNum = WlzGMModelExtend2V1E1S3D(model, cE[0], *(pos + idx2));
	}
	else
	{
	  /* Matched 2 verticies in same shell but not connected by an edge. */
	  errNum = WlzGMModelExtend2V0E1S3D(model, matchV[idx0], matchV[idx1],
					    *(pos + idx2));
	}
      }
      else
      {
	/* Matched 2 verticies in different shells. */
	errNum = WlzGMModelJoin2V0E0S3D(model, matchV[idx0], matchV[idx1],
					*(pos + idx2));
      }
      break;
    case 3:
      /* Three verticies matched, count existing edges and shells which use
       * these verticies. */
      edgeCnt = ((cE[0] = WlzGMVertexCommonEdge(matchV[0],
      						matchV[1])) != NULL) +
		((cE[1] = WlzGMVertexCommonEdge(matchV[1],
						matchV[2])) != NULL) +
		((cE[2] = WlzGMVertexCommonEdge(matchV[2],
						matchV[0])) != NULL);
      idx0 = (WlzGMVertexGetShell(matchV[0]))->idx;
      idx1 = (WlzGMVertexGetShell(matchV[1]))->idx;
      idx2 = (WlzGMVertexGetShell(matchV[2]))->idx;
      if((idx0 != idx1) && (idx1 != idx2) && (idx2 != idx0))
      {
        shellCnt = 3;
      }
      else if((idx0 == idx1) && (idx1 == idx2) && (idx2 == idx0))
      {
        shellCnt = 1;
      }
      else
      {
        shellCnt = 2;
      }
      switch(edgeCnt)
      {
        case 0:
	  switch(shellCnt)
	  {
	    case 1:
	      /* All 3 matched verticies share the same shell but there
	       * are no common edges between the 3 verticies. */
	      errNum = WlzGMModelExtend3V0E1S3D(model, matchV);
	      break;
	    case 2: /* 3V0E2S */
	      /* Two of the 3 matched verticies share the same shell but
	       * there is no common edge between any of the matched
	       * verticies. */
	      errNum = WlzGMModelJoin3V0E2S3D(model, matchV);
	      break;
	    case 3: /* 3V0E3S */
	      /* All 3 verticies are in different shells and there are
	       * no common edges between any of the matched verticies. */
	      errNum = WlzGMModelJoin3V0E3S3D(model, matchV);
	      break;
	  }
	  break;
	case 1:
	  idx2 = (cE[0] != NULL) | ((cE[1] != NULL) << 1) |
	  	 ((cE[2] != NULL) << 2);
	  idx0 = firstCtgBitTb[idx2];
	  idx1 = (idx0 + 2) % 3;
	  if(shellCnt == 1)
	  {
	    /* Two of the 3 matched verticies are connected by a single edge,
	     * all of the 3 are in a single shell. */
	     errNum = WlzGMModelExtend3V1E1S3D(model, cE[idx0], matchV[idx1]);
	  }
	  else
	  {
	    /* Two of the 3 matched verticies are connected by a single edge,
	     * the other vertex is in a different shell. */
	     errNum = WlzGMModelJoin3V1E2S3D(model, cE[idx0], matchV[idx1]);
	  }
	  break;
	case 2:
	  /* All 3 verticies share the same shell, and two pairs of verticies
	   * are connected by edges. */
	  idx2 = (cE[0] != NULL) | ((cE[1] != NULL) << 1) |
	  	 ((cE[2] != NULL) << 2);
	  idx0 = firstCtgBitTb[idx2];
	  idx1 = (idx0 + 1) % 3;
	  errNum = WlzGMModelExtend3V2E1S3D(model, cE[idx0], cE[idx1]);
	  break;
	case 3:
          /* All three verticies and all three edges are within the model,
	   * but is the loop formed by these verticies within it?
	   * TODO Check and insert loop if required. */
	  break;
      }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMLoopTSetT
* Returns:	void
* Purpose:	Walks around a loop topology element's child edge
*		topology elements setting their parent.
* Global refs:	-
* Parameters:	WlzGMLoopT *gLT:	The loop topology element.
************************************************************************/
static void	WlzGMLoopTSetT(WlzGMLoopT *gLT)
{
  WlzGMEdgeT	*fET,
  		*tET;
  
  if(gLT && ((tET = fET = gLT->edgeT) != NULL))
  {
    do
    {
      tET->parent = gLT;
    } while((tET = tET->next)->idx != fET->idx);
  }
}

/************************************************************************
* Function:	WlzGMModelAddVertex
* Returns:	void
* Purpose:	Adds a new vertex into the models vertex hash table.
*		The verticies geometry must have been set.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model.
*		WlzGMVertex *nV:	New vertex to insert into the
*					models hash table.
************************************************************************/
static void	WlzGMModelAddVertex(WlzGMModel *model, WlzGMVertex *nV)
{
  WlzDVertex3	nPos;
  unsigned int 	hVal;
  WlzGMVertex  	**hdV;
  WlzGMVertex	*tV,
		*pV;

  (void )WlzGMVertexGetG3D(nV, &nPos);
  hVal = WlzGMHashPos3D(nPos);
  hdV = model->vertexHT + (hVal % model->vertexHTSz);
  if(*hdV == NULL)
  {
    /* No verticies at this position in the hash table yet. */
    *hdV = nV;
    nV->next = NULL;
  }
  else
  {
    /* Need to search for position in the linked list. */
    tV = *hdV;
    pV = NULL;
    while(tV && (WlzGMVertexCmpSign3D(tV, nPos) < 0))
    {
      pV = tV;
      tV = tV->next;
    }
    if(pV == NULL)
    {
      nV->next = *hdV;
      *hdV = nV;
    }
    else if(tV == NULL)
    {
      nV->next = NULL;
      pV->next = nV;
    }
    else
    {
      nV->next = pV->next;
      pV->next = nV;
    }
  }
}

/************************************************************************
* Function:	WlzGMModelMatchEdgeTG2D
* Returns:	void
* Purpose:	Attempts to find verticies which match the two given
*		double precision positions.
*		The linking/matching workspace makes the matching much
*		more efficient, but it's only valid if the verticies
*		are maintained with indicies which increase with 
*		increasing row,column.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMEdgeT **matchET:	Array for return of matched
*					edge topology element pointers.
*		WlzDVertex2 *pos:	Pointer to first then second
*					positions.
************************************************************************/
static void	WlzGMModelMatchEdgeTG2D(WlzGMModel *model,
					WlzGMEdgeT **matchET,
				        WlzDVertex2 *pos)
{
  int		idD;
  double	tD0,
  		tD1,
		trX,
		trY,
		trCos,
		trSin,
		trScale;
  WlzDVertex2	tPosB,
  		tPosN,
		mPos,
  		oPos,
		tPos;
  WlzGMEdgeT	*eTB,
  		*eTL,
		*eTN,
		*eTS;
  WlzGMVertex	*vertex;
  WlzGMVertexT	*vertexT;
  WlzGMVertex	*matchV[2];
  const double	tol = 1.0e-06;

  *(matchET + 0) = *(matchET + 1) = NULL;
  /* First find any verticies that have geometry element equal one of
   * the given positions. */
  matchV[0] = WlzGMModelMatchVertexG2D(model, *(pos + 0));
  matchV[1] = WlzGMModelMatchVertexG2D(model, *(pos + 1));
  /* Now find the edge topology element(s) that are directed away from these
   * verticies. */
  for(idD = 0; idD < 2; ++idD)			 /* For each matched vertex. */
  {
    /* Check to see if this vertex has more than one vertex topology element.
     * If it does then more than one edge shares this vertex and the search
     * for the next edge topology elements becomes more complicated. */
    if((vertex = *(matchV + idD)) != NULL)
    {
      vertexT = vertex->diskT->vertexT;
      eTS = vertexT->parent;
      if(vertexT->idx == vertexT->next->idx)
      {
	/* The simple case just one edge uses this vertex. */
	*(matchET + idD) = eTS;
      }
      else
      {
	/* More than one edge uses this vertex. Oh dear! Now need to search
	 * for next edge, where the next edge is directed towards *(pos + idD)
	 * and away from *(pos + idD? 0: 1).
	 */
	mPos = *(pos + idD);
	oPos = *(pos + !idD);
	/* Compute the rigid body affine transform which takes the
	 * line segment from the matched vertex at mPos to the new vertex at
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
	  /* Search for next edge by walking around the matched vertex CCW,
	   * using the vertex topology elements, comparing the relative
	   * polar angles of the line segments between the verticies.
	   * The next edge topology element is the one with the smallest
	   * CCW polar angle relative to the new, as yet unmade, edge. */
	  eTB = eTS;	  /* Initialize the search, best (eTB) = start (eTS) */
	  (void )WlzGMVertexGetG2D(eTB->vertexT->diskT->vertex, &tPos);
	  tPosB.vtX = trScale *
		      ((tPos.vtX * trCos) - (tPos.vtY * trSin) + trX);
	  tPosB.vtY = trScale *
		      ((tPos.vtX * trSin) - (tPos.vtY * trCos) + trY);
	  /* Find next edge topology element (eTN) directed from the vertex.
	   * Because the edge topology elements are unordered need to
	   * search them all. */
	  eTN = eTB->opp->next;
	  /* For each edge topology element directed away from the vertex */
	  while(eTS->idx != eTN->idx)
	  {
	    (void )WlzGMVertexGetG2D(eTN->vertexT->diskT->vertex, &tPos);
	    tPosN.vtX = trScale *
			((tPos.vtX * trCos) - (tPos.vtY * trSin) + trX);
	    tPosN.vtY = trScale *
			((tPos.vtX * trSin) - (tPos.vtY * trCos) + trY);
	    /* Compare the best and next transformed vertex positions */
	    eTL = eTN;
	    if(WlzGeomCmpAngle(tPosN, tPosB) > 0)
	    {
	      /* Update the best (minimum angle) so far. */
	      tPosB = tPosN;
	      eTB = eTN;
	    }
	    eTN = eTL->opp->next;
	  }
	  *(matchET + idD) = eTB;
	}
      }
    }
  }
}
