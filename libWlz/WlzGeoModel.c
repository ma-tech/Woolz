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
* Purpose:      Basic operators for handling non-manifold geometric
*		models within Woolz.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>



static void	WlzGMModelMatchVtxG2D(
			  WlzGMModel *model,
			  WlzGMEdgeT **matchET,
			  WlzDVertex2 searchDist,
			  int *searchIdx,
			  WlzDVertex2 *pos);
static void	WlzGMModelMatchVtxGIdx2D(
			  WlzGMModel *model,
			  WlzGMVertex **matchV,
			  WlzDVertex2 searchDist,
			  int *searchIdx,
			  WlzDVertex2 *pos);
static void	WlzGMModelMatchVtxGLst2D(
			  WlzGMModel *model,
			  WlzGMVertex **matchV,
			  WlzDVertex2 *pos);
static void	WlzGMTestOutLinePS(
			  FILE *fP,
			  WlzDVertex2 pos0,
			  WlzDVertex2 pos1,
			  WlzDVertex2 offset,
			  WlzDVertex2 scale);

/* Creation  of geometric modeling elements */

/************************************************************************
* Function:	WlzGMModelNew
* Returns:	WlzGMModel *:		New empty model.
* Purpose:	Creates an empty non-manifold geometry model.
* Global refs:	-
* Parameters:	WlzGMModelType modType:	Type of model to create.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMModel	*WlzGMModelNew(WlzGMModelType modType, WlzErrorNum *dstErr)
{
  WlzGMModel	*model = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  unsigned int	vertexGSz,
		loopGSz,
		shellGSz;
  const unsigned int blkSz = 1024;

  switch(modType)
  {
    case WLZ_GMMOD_2I:
      vertexGSz = sizeof(WlzGMVertexG2I);
      loopGSz = sizeof(WlzGMLoopG2I);
      shellGSz = sizeof(WlzGMShellG2I);
      break;
    case WLZ_GMMOD_2D:
      vertexGSz = sizeof(WlzGMVertexG2D);
      loopGSz = sizeof(WlzGMLoopG2D);
      shellGSz = sizeof(WlzGMShellG2D);
      break;
    case WLZ_GMMOD_3I:
      vertexGSz = sizeof(WlzGMVertexG3I);
      loopGSz = sizeof(WlzGMLoopG3I);
      shellGSz = sizeof(WlzGMShellG3I);
      break;
    case WLZ_GMMOD_3D:
      vertexGSz = sizeof(WlzGMVertexG3D);
      loopGSz = sizeof(WlzGMLoopG3D);
      shellGSz = sizeof(WlzGMShellG3D);
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  /* Create the new model */
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
       ((model->res.edge.vec = AlcVectorNew(1, sizeof(WlzGMEdge),
    					    blkSz, NULL)) == NULL) ||
       ((model->res.edgeT.vec = AlcVectorNew(1, sizeof(WlzGMEdgeT),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.edgeG.vec = AlcVectorNew(1, sizeof(WlzGMEdgeG),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.face.vec = AlcVectorNew(1, sizeof(WlzGMFace),
    					    blkSz, NULL)) == NULL) ||
       ((model->res.faceT.vec = AlcVectorNew(1, sizeof(WlzGMFaceT),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.faceG.vec = AlcVectorNew(1, sizeof(WlzGMFaceG),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.loop.vec = AlcVectorNew(1, sizeof(WlzGMLoop),
    					    blkSz, NULL)) == NULL) ||
       ((model->res.loopT.vec = AlcVectorNew(1, sizeof(WlzGMLoopT),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.loopG.vec = AlcVectorNew(1, loopGSz,
    					     blkSz, NULL)) == NULL) ||
       ((model->res.shell.vec = AlcVectorNew(1, sizeof(WlzGMShell),
    					     blkSz, NULL)) == NULL) ||
       ((model->res.shellG.vec = AlcVectorNew(1, shellGSz,
    					      blkSz, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      model->res.vertex.nxtIdx = 1;
      model->res.vertexT.nxtIdx = 1;
      model->res.vertexG.nxtIdx = 1;
      model->res.edge.nxtIdx = 1;
      model->res.edgeT.nxtIdx = 1;
      model->res.edgeG.nxtIdx = 1;
      model->res.face.nxtIdx = 1;
      model->res.faceT.nxtIdx = 1;
      model->res.faceG.nxtIdx = 1;
      model->res.loop.nxtIdx = 1;
      model->res.loopT.nxtIdx = 1;
      model->res.loopG.nxtIdx = 1;
      model->res.shell.nxtIdx = 1;
      model->res.shellG.nxtIdx = 1;
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
     	        (AlcVectorExtendAndGet(resS->vec, resS->nxtIdx))) == NULL) ||
       ((shell->geo.core = (WlzGMCore *)
    			   (AlcVectorExtendAndGet(resSG->vec,
					 	  resSG->nxtIdx))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    shell->type = WLZ_GMELM_SHELL;
    shell->idx = (resS->nxtIdx)++;
    shell->geo.core->type = WlzGMModelGetSGeomType(model);
    shell->geo.core->idx = (resSG->nxtIdx)++;
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
* Function:	WlzGMModelNewF
* Returns:	WlzGMFace *:		New face.
* Purpose:	Creates an empty face with a face geometry element.
*		The face geometry element only has it's index set
*		to a meaningful value.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Parent model.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzGMFace	*WlzGMModelNewF(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource *resF,
  		*resFG;
  WlzGMFace	*face = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(model->type != WLZ_GMELM_SHELL)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    /* Get a new face and face geometry element from the model */
    resF = &(model->res.face);
    resFG = &(model->res.faceG);
    if(((face = (WlzGMFace *)
     	        (AlcVectorExtendAndGet(resF->vec, resF->nxtIdx))) == NULL) ||
       ((face->geo.core = (WlzGMCore *)
    			  (AlcVectorExtendAndGet(resFG->vec,
					         resFG->nxtIdx))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    face->type = WLZ_GMELM_FACE;
    face->idx = (resF->nxtIdx)++;
    face->geo.core->type = WLZ_GMELM_FACE_G;
    face->geo.core->idx = (resFG->nxtIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeF(model, face);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(face);
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
  WlzGMResource	*resL,
  		*resLG;
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
    resLG = &(model->res.loopG);
    if(((loop = (WlzGMLoop *)
     	        (AlcVectorExtendAndGet(resL->vec, resL->nxtIdx))) == NULL) ||
       ((loop->geo.core = (WlzGMCore *)
    			   (AlcVectorExtendAndGet(resLG->vec,
					   	  resLG->nxtIdx))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    loop->type = WLZ_GMELM_LOOP;
    loop->idx = (resL->nxtIdx)++;
    loop->geo.core->type = WlzGMModelGetLGeomType(model);
    loop->geo.core->idx = (resLG->nxtIdx)++;
  }
  /* Clear up on error */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFreeL(model, loop);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(loop);
}

/************************************************************************
* Function:	WlzGMModelNewLT
* Returns:	WlzGMLoop *:		New loop.
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
    /* Get a new loop and loop geometry element from the model */
    resLT = &(model->res.loopT);
    if((loopT = (WlzGMLoopT *)
     	        (AlcVectorExtendAndGet(resLT->vec, resLT->nxtIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    loopT->type = WLZ_GMELM_LOOP_T;
    loopT->idx = (resLT->nxtIdx)++;
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
  WlzGMResource	*resE,
  		*resEG;
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
    resEG = &(model->res.edgeG);
    if(((edge = (WlzGMEdge *)
     	        (AlcVectorExtendAndGet(resE->vec, resE->nxtIdx))) == NULL) ||
       ((edge->geo.core = (WlzGMCore *)
			  (AlcVectorExtendAndGet(resEG->vec,
					   	 resEG->nxtIdx))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    edge->type = WLZ_GMELM_EDGE;
    edge->idx = (model->res.edge.nxtIdx)++;
    edge->geo.core->type = WLZ_GMELM_EDGE_G;
    edge->geo.core->idx = (resEG->nxtIdx)++;
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
     	        (AlcVectorExtendAndGet(resET->vec, resET->nxtIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    edgeT->type = WLZ_GMELM_EDGE_T;
    edgeT->idx = (resET->nxtIdx)++;
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
     	        (AlcVectorExtendAndGet(resV->vec, resV->nxtIdx))) == NULL) ||
       ((vertex->geo.core = (WlzGMCore *)
    			    (AlcVectorExtendAndGet(resVG->vec,
					 	   resVG->nxtIdx))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vertex->type = WLZ_GMELM_VERTEX;
    vertex->idx = (resV->nxtIdx)++;
    vertex->geo.core->type = WlzGMModelGetVGeomType(model);
    vertex->geo.core->idx = (resVG->nxtIdx)++;
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
     	          (AlcVectorExtendAndGet(resVT->vec, resVT->nxtIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vertexT->type = WLZ_GMELM_VERTEX_T;
    vertexT->idx = (resVT->nxtIdx)++;
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
    (void )AlcVectorFree(model->res.edge.vec);
    (void )AlcVectorFree(model->res.edgeT.vec);
    (void )AlcVectorFree(model->res.edgeG.vec);
    (void )AlcVectorFree(model->res.face.vec);
    (void )AlcVectorFree(model->res.faceT.vec);
    (void )AlcVectorFree(model->res.faceG.vec);
    (void )AlcVectorFree(model->res.loop.vec);
    (void )AlcVectorFree(model->res.loopT.vec);
    (void )AlcVectorFree(model->res.loopG.vec);
    (void )AlcVectorFree(model->res.shell.vec);
    (void )AlcVectorFree(model->res.shellG.vec);
    AlcFree((void *)model);
  }
  return(errNum);
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
     * as invalid by making the indicies = 0. */
    shell->idx = 0;
    if(shell->geo.core != NULL)
    {
      shell->geo.core->idx = 0;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeL
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a loop and it's geometry as invalid and suitable
*		for reclaiming.
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
    /* Can't really free the loop so just mark it and it's geometry element
     * as invalid by making the indicies = 0. */
    loop->idx = 0;
    if(loop->geo.core != NULL)
    {
      loop->geo.core->idx = 0;
    }
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
     * the index = 0. */
    loopT->idx = 0;
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeF
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a face and it's geometry as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMFace *face:	Face to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeF(WlzGMModel *model, WlzGMFace *face)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (face == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the face so just mark it and it's geometry element
     * as invalid by making the indicies = 0. */
    face->idx = 0;
    if(face->geo.core != NULL)
    {
      face->geo.core->idx = 0;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeFT
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks a face topology as invalid and suitable
*		for reclaiming.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMFaceT *faceT:	FaceT to free.
************************************************************************/
WlzErrorNum	WlzGMModelFreeFT(WlzGMModel *model, WlzGMFaceT *faceT)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (faceT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Can't really free the faceT so just mark it as invalid by making
     * the index = 0. */
    faceT->idx = 0;
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelFreeE
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Marks an edge and it's geometry as invalid and suitable
*		for reclaiming.
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
    /* Can't really free the edge so just mark it and it's geometry element
     * as invalid by making the indicies = 0. */
    edge->idx = 0;
    if(edge->geo.core != NULL)
    {
      edge->geo.core->idx = 0;
    }
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
     * the index = 0. */
    edgeT->idx = 0;
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
     * as invalid by making the indicies = 0. */
    vertex->idx = 0;
    if(vertex->geo.core != NULL)
    {
      vertex->geo.core->idx = 0;
    }
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
     * the index = 0. */
    vertexT->idx = 0;
  }
  return(errNum);
}

/* Model access and testing */

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
* Function:	WlzGMModelTestOutPS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 2D model to the specified file as Postscript.
*		To make shell identification easier Postscript colors
*		can be cycled.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell.
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
  WlzGMShell	*shell0,
  		*shell1;
  int		colIdx  = 0;
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
    if((shell0 = model->child) != NULL)
    {
      if(nCol > maxCol)
      {
        nCol = maxCol;
      }
      shell1 = shell0;
      do
      {
        /* Set colour. */
	if(fP && (nCol > 1))
	{
	  (void )fprintf(fP,
	  	         "%d %d %d setrgbcolor\n",
			 colR[colIdx], colG[colIdx], colB[colIdx]);
	  colIdx = ++colIdx % nCol;
	}
	errNum = WlzGMShellTestOutPS(shell1, fP, offset, scale);
	if((shell1 = shell1->next) == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
      }
      while((errNum == WLZ_ERR_NONE) && (shell1->idx != shell0->idx));
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMShellTestOutPS
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Outputs a 2D shell to the specified file as Postscript.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell.
*		FILE *fP:		Output file ptr, may be NULL if
*					no output required.
*		WlzDVertex2 offset:	Postscript offset.
*		WlzDVertex2 scale:	Postscript scale.
************************************************************************/
WlzErrorNum   	WlzGMShellTestOutPS(WlzGMShell *shell, FILE *fP,
				WlzDVertex2 offset, WlzDVertex2 scale)
{
  WlzGMLoopT	*loopT0,
  		*loopT1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(shell == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((loopT0 = shell->child.loopT) == NULL)
  {
    errNum = WLZ_ERR_NONE;	         /* TODO Is an empty shell an error? */
  }
  else
  {
    loopT1 = loopT0;
    do
    {
      if((errNum = WlzGMLoopTTestOutPS(loopT1, fP,
      				       offset, scale)) == WLZ_ERR_NONE)
      {
        loopT1 = loopT1->next;
      }
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
*		FILE *fP:		Output file ptr, may be NULL if
*					no output required.
*		WlzDVertex2 offset:	Postscript offset.
*		WlzDVertex2 scale:	Postscript scale.
************************************************************************/
WlzErrorNum   	WlzGMLoopTTestOutPS(WlzGMLoopT *loopT, FILE *fP,
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
      if((errNum = WlzGMEdgeTTestOutPS(edgeT1, fP,
      				       offset, scale)) == WLZ_ERR_NONE)
      {
        edgeT1 = edgeT1->next;
      }
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
WlzErrorNum   	WlzGMEdgeTTestOutPS(WlzGMEdgeT *edgeT, FILE *fP,
			            WlzDVertex2 offset, WlzDVertex2 scale)
{
  WlzDVertex2	pos0,
  		pos1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((edgeT == NULL) || (edgeT->opp == NULL) ||
     (edgeT->vertexT == NULL) || (edgeT->opp->vertexT == NULL) ||
     (edgeT->vertexT->vertex == NULL) || (edgeT->opp->vertexT->vertex == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    if(((errNum = WlzGMVertexGetG2D(edgeT->vertexT->vertex,
    				    &pos0)) == WLZ_ERR_NONE) &&
       ((errNum = WlzGMVertexGetG2D(edgeT->next->vertexT->vertex,
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
* Function:	WlzGMModelGetLGeomType
* Returns:	WlzGMElemType:		Shell's geometry type.
* Purpose:	Gets the loop's geometry type from the model's type.
* Global refs:	-
* Parameters:	WlzGMModel *model:	The model.
************************************************************************/
WlzGMElemType 	WlzGMModelGetLGeomType(WlzGMModel *model)
{
  WlzGMElemType	lGType;

  switch(model->type)
  {
    case WLZ_GMMOD_2I:
      lGType = WLZ_GMELM_LOOP_G2I;
      break;
    case WLZ_GMMOD_2D:
      lGType = WLZ_GMELM_LOOP_G2D;
      break;
    case WLZ_GMMOD_3I:
      lGType = WLZ_GMELM_LOOP_G3I;
      break;
    case WLZ_GMMOD_3D:
      lGType = WLZ_GMELM_LOOP_G3D;
      break;
  }
  return(lGType);
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
* Function:	WlzGMShellSetG2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a shell's geometry using the pair of double
*		precision points.
* Global refs:	-
* Parameters:	WlzGMShell *shell:	Given shell with geometry to
*					be set.
*		WlzDVertex2 pos0:	First point position.
*		WlzDVertex2 pos1:	Second point position.
************************************************************************/
WlzErrorNum	WlzGMShellSetG2D(WlzGMShell *shell,
				  WlzDVertex2 pos0, WlzDVertex2 pos1)
{
  WlzDBox2	bBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(pos0.vtX > pos1.vtX)
  {
    bBox.xMin = pos1.vtX;
    bBox.xMax = pos0.vtX;
  }
  else
  {
    bBox.xMin = pos0.vtX;
    bBox.xMax = pos1.vtX;
  }
  if(pos0.vtY > pos1.vtY)
  {
    bBox.yMin = pos1.vtY;
    bBox.yMax = pos0.vtY;
  }
  else
  {
    bBox.yMin = pos0.vtY;
    bBox.yMax = pos1.vtY;
  }
  (void )WlzGMShellSetGBB2D(shell, bBox);
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
	dstPos->vtX = vertex->geo.vg2D->vtx.vtX;
	dstPos->vtY = vertex->geo.vg2D->vtx.vtY;
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
* Function:	WlzGMLoopSetGBB2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a loop's geometry using the given double
*		precision bounding box.
* Global refs:	-
* Parameters:	WlzGMLoop *loop:	Given loop with geometry to
*					be set.
*		WlzDBox2 bBox:		Given double precision bounding
*					box.
************************************************************************/
WlzErrorNum	WlzGMLoopSetGBB2D(WlzGMLoop *loop, WlzDBox2 bBox)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((loop == NULL) || (loop->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(loop->geo.core->type)
    {
      case WLZ_GMELM_LOOP_G2I:
        loop->geo.lg2I->bBox.xMin = WLZ_NINT(bBox.xMin);
        loop->geo.lg2I->bBox.yMin = WLZ_NINT(bBox.yMin);
        loop->geo.lg2I->bBox.xMax = WLZ_NINT(bBox.xMax);
        loop->geo.lg2I->bBox.yMax = WLZ_NINT(bBox.yMax);
	break;
      case WLZ_GMELM_LOOP_G2D:
        loop->geo.lg2D->bBox = bBox;
	break;
      case WLZ_GMELM_LOOP_G3I:
        loop->geo.lg3I->bBox.xMin = WLZ_NINT(bBox.xMin);
        loop->geo.lg3I->bBox.yMin = WLZ_NINT(bBox.yMin);
        loop->geo.lg3I->bBox.zMin = 0;
        loop->geo.lg3I->bBox.xMax = WLZ_NINT(bBox.xMax);
        loop->geo.lg3I->bBox.yMax = WLZ_NINT(bBox.yMax);
        loop->geo.lg3I->bBox.zMax = 0;
	break;
      case WLZ_GMELM_LOOP_G3D:
        loop->geo.lg3D->bBox.xMin = bBox.xMin;
        loop->geo.lg3D->bBox.yMin = bBox.yMin;
        loop->geo.lg3D->bBox.zMin = 0.0;
        loop->geo.lg3D->bBox.xMax = bBox.xMax;
        loop->geo.lg3D->bBox.yMax = bBox.yMax;
        loop->geo.lg3D->bBox.zMax = 0.0;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMLoopSetG2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Sets a loop's geometry using the pair of double
*		precision points.
* Global refs:	-
* Parameters:	WlzGMLoop *loop:	Given loop with geometry to
*					be set.
*		WlzDVertex2 pos0:	First point position.
*		WlzDVertex2 pos1:	Second point position.
************************************************************************/
WlzErrorNum	WlzGMLoopSetG2D(WlzGMLoop *loop,
				  WlzDVertex2 pos0, WlzDVertex2 pos1)
{
  WlzDBox2	bBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(pos0.vtX > pos1.vtX)
  {
    bBox.xMin = pos1.vtX;
    bBox.xMax = pos0.vtX;
  }
  else
  {
    bBox.xMin = pos0.vtX;
    bBox.xMax = pos1.vtX;
  }
  if(pos0.vtY > pos1.vtY)
  {
    bBox.yMin = pos1.vtY;
    bBox.yMax = pos0.vtY;
  }
  else
  {
    bBox.yMin = pos0.vtY;
    bBox.yMax = pos1.vtY;
  }
  (void )WlzGMLoopSetGBB2D(loop, bBox);
  return(errNum);
}

/************************************************************************
* Function:	WlzGMLoopUpdateG2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Updates a loop's geometry using the given double
*		precision position.
* Global refs:	-
* Parameters:	WlzGMLoop *loop:	Given loop with geometry to
*					be set.
*		WlzDVertex2 pos:	Given position.
************************************************************************/
WlzErrorNum	WlzGMLoopUpdateG2D(WlzGMLoop *loop, WlzDVertex2 pos)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((loop == NULL) || (loop->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(loop->geo.core->type)
    {
      case WLZ_GMELM_LOOP_G2I:
	if((tI0 = pos.vtX) < loop->geo.lg2I->bBox.xMin)
	{
	  loop->geo.lg2I->bBox.xMin = tI0;
	}
	else if(tI0 > loop->geo.lg2I->bBox.xMax)
	{
	  loop->geo.lg2I->bBox.xMax = tI0;
	}
	if((tI0 = pos.vtY) < loop->geo.lg2I->bBox.yMin)
	{
	  loop->geo.lg2I->bBox.yMin = tI0;
	}
	else if(tI0 > loop->geo.lg2I->bBox.yMax)
	{
	  loop->geo.lg2I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_LOOP_G2D:
	if(pos.vtX < loop->geo.lg2D->bBox.xMin)
	{
	  loop->geo.lg2D->bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > loop->geo.lg2D->bBox.xMax)
	{
	  loop->geo.lg2D->bBox.xMax = pos.vtX;
	}
	if(pos.vtY < loop->geo.lg2D->bBox.yMin)
	{
	  loop->geo.lg2D->bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > loop->geo.lg2D->bBox.yMax)
	{
	  loop->geo.lg2D->bBox.yMax = pos.vtY;
	}
	break;
      case WLZ_GMELM_LOOP_G3I:
	if((tI0 = pos.vtX) < loop->geo.lg3I->bBox.xMin)
	{
	  loop->geo.lg3I->bBox.xMin = tI0;
	}
	else if(tI0 > loop->geo.lg3I->bBox.xMax)
	{
	  loop->geo.lg3I->bBox.xMax = tI0;
	}
	if((tI0 = pos.vtY) < loop->geo.lg3I->bBox.yMin)
	{
	  loop->geo.lg3I->bBox.yMin = tI0;
	}
	else if(tI0 > loop->geo.lg3I->bBox.yMax)
	{
	  loop->geo.lg3I->bBox.yMax = tI0;
	}
	break;
      case WLZ_GMELM_LOOP_G3D:
	if(pos.vtX < loop->geo.lg3D->bBox.xMin)
	{
	  loop->geo.lg3D->bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > loop->geo.lg3D->bBox.xMax)
	{
	  loop->geo.lg3D->bBox.xMax = pos.vtX;
	}
	if(pos.vtY < loop->geo.lg3D->bBox.yMin)
	{
	  loop->geo.lg3D->bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > loop->geo.lg3D->bBox.yMax)
	{
	  loop->geo.lg3D->bBox.yMax = pos.vtY;
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
* Function:	WlzGMLoopGInBB2D
* Returns:	int:			Non zero if point inside loop's
*					bounding box.
* Purpose:	Checks to see if the given double precision position
*		is within the loop's bounding box.
* Global refs:	-
* Parameters:	WlzGMLoop *loop:	Given loop with geometry to
*					be set.
*		WlzDVertex2 pos:	Given double precision position.
************************************************************************/
int		WlzGMLoopGInBB2D(WlzGMLoop *loop, WlzDVertex2 pos)
{
  int		inside = 0;

  if((loop != NULL) && (loop->geo.core != NULL))
  {
    switch(loop->geo.core->type)
    {
      case WLZ_GMELM_LOOP_G2I:
	inside = (pos.vtX >= loop->geo.lg2I->bBox.xMin) &&
		 (pos.vtX <= loop->geo.lg2I->bBox.xMax) &&
	         (pos.vtY >= loop->geo.lg2I->bBox.yMin) &&
		 (pos.vtY <= loop->geo.lg2I->bBox.yMax);
	break;
      case WLZ_GMELM_LOOP_G2D:
	inside = (pos.vtX >= loop->geo.lg2D->bBox.xMin) &&
		 (pos.vtX <= loop->geo.lg2D->bBox.xMax) &&
	         (pos.vtY >= loop->geo.lg2D->bBox.yMin) &&
		 (pos.vtY <= loop->geo.lg2D->bBox.yMax);
	break;
      case WLZ_GMELM_LOOP_G3I:
	inside = (pos.vtX >= loop->geo.lg3I->bBox.xMin) &&
		 (pos.vtX <= loop->geo.lg3I->bBox.xMax) &&
	         (pos.vtY >= loop->geo.lg3I->bBox.yMin) &&
		 (pos.vtY <= loop->geo.lg3I->bBox.yMax);
	break;
      case WLZ_GMELM_LOOP_G3D:
	inside = (pos.vtX >= loop->geo.lg3D->bBox.xMin) &&
		 (pos.vtX <= loop->geo.lg3D->bBox.xMax) &&
	         (pos.vtY >= loop->geo.lg3D->bBox.yMin) &&
		 (pos.vtY <= loop->geo.lg3D->bBox.yMax);
	break;
      default:
	break;
    }
  }
  return(inside);
}


/* Topology query */

/************************************************************************
* Function:	WlzGMEdgeTCommonParent
* Returns:	WlzGMLoopT:		Common loop topology element,
					NULL if it doesn't exist.
* Purpose:	Finds a loop topology element in common for the two
*		edge topology elements.
*		edge topology elements are maintained in an ordered
*		doubly linked list.
* Global refs:	-
* Parameters:	WlzGMEdgeT *eT0:	First edge topology element.
*		WlzGMEdgeT *eT1:	Second edge topology element.
************************************************************************/
WlzGMLoopT	*WlzGMEdgeTCommonParent(WlzGMEdgeT *eT0, WlzGMEdgeT *eT1)
{
  WlzGMLoopT	*loopT = NULL;

  if((eT0->parent.core->type == WLZ_GMELM_LOOP_T) &&
     (eT0->parent.core->type == eT1->parent.core->type) &&
     (eT0->parent.loopT->idx == eT1->parent.loopT->idx))
  {
    loopT = eT0->parent.loopT;
  }
  return(loopT);
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
* Function:	WlzGMEdgeTAppend
* Returns:	void
* Purpose:	Append new edge topology onto a doubly linked list
*		of edge topology element's, knowing that neither is NULL.
*		Edge topology elements are maintained in an ordered
*		doubly linked list, CCW around the inside of loops and
*		CW around the outside of loops.
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
*		CW around the outside of loops.
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
  switch(dLT->parent.core->type)
  {
    case WLZ_GMELM_FACE_T:
      if(dLT->parent.faceT->loopT->idx == dLT->idx)
      {
        dLT->parent.faceT->loopT = dLT->next;
      }
      break;
    case WLZ_GMELM_SHELL:
      if(dLT->parent.shell->child.loopT->idx == dLT->idx)
      {
        dLT->parent.shell->child.loopT  = dLT->next;
      }
      break;
    default:
      /* Should never get to here! */
      break;
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
  /* Unlink the shell from the doubly linked list of shells. */
  dS->prev->next = dS->next;
  dS->next->prev = dS->prev;
  if(dS->idx == dS->next->idx)
  {
    dS->next = dS->prev = NULL;
  }
  if(dS->idx == dS->parent->child->idx)
  {
    /* Change child held by parent model */
    dS->parent->child = dS->next;
  }
}

/* Model construction */

/************************************************************************
* Function:	WlzGMModelConstructSLE2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a new shell, loop and edge (with topological
*		elements) and adds them to the model. Given a pair
*		of 2D double precision edge end points which are known
*		not to be in the model.
*		Neither of the verticies already exists within the shell.
*		Need to create a new shell with: 1 loop (nL), 1 loop
*		topology element (nLT), 1 edge (nE) 2 edg
*		topology elements (nET0, nET1), 2 verticies (nV0
*		with geometry *(pos + 0), nV1 with geometry
*		*(pos + 1)), 2 vertex topology elements (nVT0, nVT1).
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
*		WlzDVertex2 pos0:	First end point.
*		WlzDVertex2 pos1:	The other end point.
************************************************************************/
WlzErrorNum	WlzGMModelConstructSLE2D(WlzGMModel *model,
				   	 WlzDVertex2 nPos0, WlzDVertex2 nPos1)
{
  WlzGMShell	*eShell,
  		*nShell;
  WlzGMLoop	*nL;
  WlzGMLoopT	*nLT;
  WlzGMVertex	*nV0,
  		*nV1;
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0,
  		*nET1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nShell = WlzGMModelNewS(model, &errNum)) != NULL) &&
     ((nL = WlzGMModelNewL(model, &errNum)) != NULL) &&
     ((nLT = WlzGMModelNewLT(model, &errNum)) != NULL) &&
     ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nV0 = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nV1 = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1 = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    /* New vertcies and their topology elements */
    nV0->top = nVT0;
    (void )WlzGMVertexSetG2D(nV0, nPos0);
    nV1->top = nVT1;
    (void )WlzGMVertexSetG2D(nV1, nPos1);
    nVT0->vertex = nV0;
    nVT0->next = nVT0->prev = nVT0;
    nVT0->parent.edgeT = nET0;
    nVT1->vertex = nV1;
    nVT1->next = nVT1->prev = nVT1;
    nVT1->parent.edgeT = nET1;
    /* New Edges and their topology elements */
    nE->top = nET0;
    nET0->edge = nET1->edge = nE;
    nET0->rad = nET0->next = nET0->prev = nET0->opp = nET1;
    nET0->vertexT = nVT0;
    nET0->parent.loopT = nLT;
    nET1->rad = nET1->next = nET1->prev = nET1->opp = nET0;
    nET1->vertexT = nVT1;
    nET1->parent.loopT = nLT;
    /* New loop and loop topology elements */
    (void )WlzGMLoopSetG2D(nL, nPos0, nPos1);
    nL->top = nLT;
    nLT->next = nLT->prev = nLT->opp = nLT;
    nLT->loop = nL;
    nLT->edgeT = nET0;
    nLT->parent.shell = nShell;
    /* New shell in model */
    (void )WlzGMShellSetG2D(nShell, nPos0, nPos1);
    nShell->child.loopT = nLT;
    nShell->parent = model;
    /* Is this the first shell of the model? */
    if((eShell = model->child) == NULL)
    {
      /* First shell in model */
      model->child = nShell;
      nShell->next = nShell->prev = nShell;
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
* Function:	WlzGMModelConstructE2D
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Constructs a new edge (with topological elements) and
*		adds them to the model. Given an existing edge topology
*		element at one end of the new edge and a 2D double
*		precision end points at the other end of the new edge,
*		where the new end point is known not to be in the model.
*		Only one vertex already exists within the model so no new
*		shell or loop is required, instead an existing loop is
*		extended using: 1 edge (nE), 2 edge topology elements
*		(nET0, nET1), 1 vertex (nV0 with geometry (nPos)
*		and 2 vertex topology elements (nVT0, nVT1).
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
WlzErrorNum	WlzGMModelConstructE2D(WlzGMModel *model,
				       WlzGMEdgeT *eET0,
				       WlzDVertex2 nPos)
{
  WlzGMVertex	*nV;
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*eET1,
		*nET0,
  		*nET1;
  WlzGMLoopT	*eLT;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model && eET0 &&
     ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
     ((nET0 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nET1 = WlzGMModelNewET(model, &errNum)) != NULL) &&
     ((nV = WlzGMModelNewV(model, &errNum)) != NULL) &&
     ((nVT0 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
     ((nVT1 = WlzGMModelNewVT(model, &errNum)) != NULL))
  {
    eET1 = eET0->prev;
    /* Vertcies and their topology elements */
    nV->top = nVT0;
    (void )WlzGMVertexSetG2D(nV, nPos);
    nVT0->prev = nVT0->next = nVT0;
    nVT0->vertex = nV;
    nVT0->parent.edgeT = nET0;
    WlzGMVertexTAppend(eET0->vertexT, nVT1);
    nVT1->vertex = eET0->vertexT->vertex;
    nVT1->parent.edgeT = nET1;
    /* Edges and their topology elements */
    nE->top = nET0;
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
    /* Just need to update loop and shell geometry */
    eLT = eET0->parent.loopT;
    (void )WlzGMLoopUpdateG2D(eLT->loop, nPos);
    (void )WlzGMShellUpdateG2D(eLT->parent.shell, nPos);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelConstructESplitL2D
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
WlzErrorNum	WlzGMModelConstructESplitL2D(WlzGMModel *model,
					     WlzGMEdgeT *eET0, WlzGMEdgeT *eET1)
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
      nVT0->vertex = eET1->vertexT->vertex;
      nVT0->parent.edgeT = nET0;
      WlzGMVertexTAppend(eET0->vertexT, nVT1);
      nVT1->vertex = eET0->vertexT->vertex;
      nVT1->parent.edgeT = nET1;
      /* Edges and their topology elements */
      nE->top = nET0;
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
      eLT = eET0->parent.loopT;
      eLT->edgeT = nET0;
      nL->top = nLT;
      WlzGMLoopTAppend(eET0->parent.loopT, nLT);       /* Add loopT to shell */
      nLT->opp = nLT;
      nLT->loop = nL;
      nLT->parent = eLT->parent;
      /* Set/update the loops by walking around the edge topology elements
       * their parents and updating the geometry. */
      eLT->edgeT = nET0;
      nLT->edgeT = nET1;
      WlzGMLoopSetGT(eLT);
      WlzGMLoopSetGT(nLT);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMModelConstructEJoinL2D
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
WlzErrorNum	WlzGMModelConstructEJoinL2D(WlzGMModel *model,
				      WlzGMEdgeT *eET0, WlzGMEdgeT *eET1)
{
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0,
  		*nET1;
  WlzGMLoop	*eL,
  		*dL;
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
    nVT0->vertex = eET0->vertexT->vertex; /* eV0 */
    nVT0->parent.edgeT = nET0;
    WlzGMVertexTAppend(eET1->vertexT, nVT1);
    nVT1->vertex = eET1->vertexT->vertex; /* eV1 */
    nVT1->parent.edgeT = nET1;
    /* Edges and their topology elements */
    nE->top = nET0;
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
    eLT = eET0->parent.loopT;
    dLT = eET1->parent.loopT;
    eL = eLT->loop;
    dL = dLT->loop;
    dS = dLT->parent.shell;
    /* Look to see if loops share the same shell. */
    if((dS = dLT->parent.shell)->idx == eLT->parent.shell->idx)
    {
      dS = NULL;
    }
    /* Set/update the retained loop by walking around the edge topology
     * elements setting their parents and updating the geometry. */
    WlzGMLoopSetGT(eLT);
    /* Unlink and free the deleted loop and loop topology element. */
    WlzGMLoopTUnlink(dLT);
    (void )WlzGMModelFreeLT(model, dLT);
    (void )WlzGMModelFreeL(model, dL);
    /* If deleted loop had a different parent shell then unlink it and
     * free the deleted loops parent shell if it's childless. */
    if(dS)
    {
      WlzGMShellUnlink(dS);
      if(dS->child.core == NULL)
      {
        (void )WlzGMModelFreeS(model, dS);
      }
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
*		WlzDVertex2 searchDist:	Maximum search distance ignored,
*					if search indicies pointer is NULL.
*		int *searchIdx:		Last line vertex search index,
*					ignored if NULL.
*		WlzDVertex2 *pos:	Pointer to first then second
*					positions.
************************************************************************/
WlzErrorNum	WlzGMModelConstructSimplex2D(WlzGMModel *model,
					     WlzDVertex2 searchDist,
				       	     int *searchIdx, WlzDVertex2 *pos)
{
  unsigned int	matchCode;
  WlzGMLoopT	*eLT;
  WlzGMEdgeT	*matchEdgeT[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WlzGMModelMatchVtxG2D(model, matchEdgeT, searchDist, searchIdx, pos);
  matchCode = ((matchEdgeT[1] != NULL) << 1) | (matchEdgeT[0] != NULL);
  switch(matchCode)
  {
    case 0:
      /* Edge makes a new shell. */
      errNum = WlzGMModelConstructSLE2D(model, *(pos + 0), *(pos + 1));
      break;
    case 1:
      /* Edge extends a loop and shell. */
      errNum = WlzGMModelConstructE2D(model, matchEdgeT[0], *(pos + 1));
      break;
    case 2:
      /* Edge extends a loop and shell. */
      errNum = WlzGMModelConstructE2D(model, matchEdgeT[1], *(pos + 0));
      break;
    case 3:
      /* Both verticies already exist within the model so adding another
       * edge between the two verticies will change the number of loops
       * and shells. */
      if((eLT = WlzGMEdgeTCommonParent(matchEdgeT[0], matchEdgeT[1])) != NULL)
      {
        /* Both of the verticies lie on a common loop then that loop needs
	 * to be SPLIT by a new edge. */
	errNum = WlzGMModelConstructESplitL2D(model,
					      matchEdgeT[0], matchEdgeT[1]);
      }
      else
      {
        /* The two verticies do NOT share a common loop so loops will be
         * JOINed by a new edge, destroying a loop.  */
	errNum = WlzGMModelConstructEJoinL2D(model,
				      	     matchEdgeT[0], matchEdgeT[1]);
        
      }
      break;
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzGMLoopSetGT
* Returns:	void
* Purpose:	Walks around a loop topology element's child edge
*		topology elements setting their parent and the loop's
*		geometry.
*		It is not legal for a model to have inconsistant
*		geometries and this function relies on the loop 
*		geometry being of the same type as the vertex geometry.
* Global refs:	-
* Parameters:	WlzGMLoopT *gLT:	The loop topology element.
************************************************************************/
void		WlzGMLoopSetGT(WlzGMLoopT *gLT)
{
  int		tI0;
  double	tD0;
  WlzGMVertex	*tV;
  WlzGMEdgeT	*fET,
  		*tET;
  WlzGMLoopGU	gLG;
  
  if(gLT && ((fET = gLT->edgeT) != NULL))
  {
    tET = fET;
    gLG = gLT->loop->geo;
    tET->parent.loopT = gLT;
    tV = tET->vertexT->vertex;
    switch(tV->geo.core->type)
    {
      case WLZ_GMELM_VERTEX_G2I:
	gLG.lg2I->bBox.xMin = gLG.lg2I->bBox.xMax =  tV->geo.vg2I->vtx.vtX;
	gLG.lg2I->bBox.yMin = gLG.lg2I->bBox.yMax =  tV->geo.vg2I->vtx.vtY;
	while((tET = tET->next)->idx != fET->idx)
	{
	  tET->parent.loopT = gLT;
	  tV = tET->vertexT->vertex;
	  if((tI0 = tV->geo.vg2I->vtx.vtX) > gLG.lg2I->bBox.xMax)
	  {
	    gLG.lg2I->bBox.xMax = tI0;
	  }
	  else if(tI0 < gLG.lg2I->bBox.xMin)
	  {
	    gLG.lg2I->bBox.xMin = tI0;
	  }
	  if((tI0 = tV->geo.vg2I->vtx.vtY) > gLG.lg2I->bBox.yMax)
	  {
	    gLG.lg2I->bBox.yMax = tI0;
	  }
	  else if(tI0 < gLG.lg2I->bBox.yMin)
	  {
	    gLG.lg2I->bBox.yMin = tI0;
	  }
	}
	break;
      case WLZ_GMELM_VERTEX_G2D:
	gLG.lg2D->bBox.xMin = gLG.lg2D->bBox.xMax =  tV->geo.vg2D->vtx.vtX;
	gLG.lg2D->bBox.yMin = gLG.lg2D->bBox.yMax =  tV->geo.vg2D->vtx.vtY;
	while((tET = tET->next)->idx != fET->idx)
	{
	  tET->parent.loopT = gLT;
	  tV = tET->vertexT->vertex;
	  if((tD0 = tV->geo.vg2D->vtx.vtX) > gLG.lg2D->bBox.xMax)
	  {
	    gLG.lg2D->bBox.xMax = tD0;
	  }
	  else if(tD0 < gLG.lg2D->bBox.xMin)
	  {
	    gLG.lg2D->bBox.xMin = tD0;
	  }
	  if((tD0 = tV->geo.vg2D->vtx.vtY) > gLG.lg2D->bBox.yMax)
	  {
	    gLG.lg2D->bBox.yMax = tD0;
	  }
	  else if(tD0 < gLG.lg2D->bBox.yMin)
	  {
	    gLG.lg2D->bBox.yMin = tD0;
	  }
	}
	break;
      case WLZ_GMELM_VERTEX_G3I:
	gLG.lg3I->bBox.xMin = gLG.lg3I->bBox.xMax =  tV->geo.vg3I->vtx.vtX;
	gLG.lg3I->bBox.yMin = gLG.lg3I->bBox.yMax =  tV->geo.vg3I->vtx.vtY;
	gLG.lg3I->bBox.zMin = gLG.lg3I->bBox.zMax =  tV->geo.vg3I->vtx.vtZ;
	while((tET = tET->next)->idx != fET->idx)
	{
	  tET->parent.loopT = gLT;
	  tV = tET->vertexT->vertex;
	  if((tI0 = tV->geo.vg3I->vtx.vtX) > gLG.lg3I->bBox.xMax)
	  {
	    gLG.lg3I->bBox.xMax = tI0;
	  }
	  else if(tI0 < gLG.lg3I->bBox.xMin)
	  {
	    gLG.lg3I->bBox.xMin = tI0;
	  }
	  if((tI0 = tV->geo.vg3I->vtx.vtY) > gLG.lg3I->bBox.yMax)
	  {
	    gLG.lg3I->bBox.yMax = tI0;
	  }
	  else if(tI0 < gLG.lg3I->bBox.yMin)
	  {
	    gLG.lg3I->bBox.yMin = tI0;
	  }
	  if((tI0 = tV->geo.vg3I->vtx.vtZ) > gLG.lg3I->bBox.zMax)
	  {
	    gLG.lg3I->bBox.zMax = tI0;
	  }
	  else if(tI0 < gLG.lg3I->bBox.zMin)
	  {
	    gLG.lg3I->bBox.zMin = tI0;
	  }
	}
	break;
      case WLZ_GMELM_VERTEX_G3D:
	gLG.lg3D->bBox.xMin = gLG.lg3D->bBox.xMax =  tV->geo.vg3D->vtx.vtX;
	gLG.lg3D->bBox.yMin = gLG.lg3D->bBox.yMax =  tV->geo.vg3D->vtx.vtY;
	gLG.lg3D->bBox.zMin = gLG.lg3D->bBox.zMax =  tV->geo.vg3D->vtx.vtZ;
	while((tET = tET->next)->idx != fET->idx)
	{
	  tET->parent.loopT = gLT;
	  tV = tET->vertexT->vertex;
	  if((tD0 = tV->geo.vg3D->vtx.vtX) > gLG.lg3D->bBox.xMax)
	  {
	    gLG.lg3D->bBox.xMax = tD0;
	  }
	  else if(tD0 < gLG.lg3D->bBox.xMin)
	  {
	    gLG.lg3D->bBox.xMin = tD0;
	  }
	  if((tD0 = tV->geo.vg3D->vtx.vtY) > gLG.lg3D->bBox.yMax)
	  {
	    gLG.lg3D->bBox.yMax = tD0;
	  }
	  else if(tD0 < gLG.lg3D->bBox.yMin)
	  {
	    gLG.lg3D->bBox.yMin = tD0;
	  }
	  if((tD0 = tV->geo.vg3D->vtx.vtZ) > gLG.lg3D->bBox.zMax)
	  {
	    gLG.lg3D->bBox.zMax = tD0;
	  }
	  else if(tD0 < gLG.lg3D->bBox.zMin)
	  {
	    gLG.lg3D->bBox.zMin = tD0;
	  }
	}
	break;
    }
  }
}

/************************************************************************
* Function:	WlzGMModelMatchVtxG2D
* Returns:	void
* Purpose:	Attempts to find verticies which match the two given
*		double precision positions.
*		The forwards and backwards  vertex search indicies and
*		the double match distance are only valid if the
*		verticies are maintained with indicies which increase
*		with increasing row,column. Using these values makes
*		the search for matching verticies MUCH more efficient.
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMedgeT **matchET:	Array for return of matched
*					edge topology element pointers.
*		WlzDVertex2 searchDist:	Maximum search distance ignored,
*					if search indicies pointer is NULL.
*		int *searchIdx:		Last line vertex search index,
*					ignored if NULL.
*		WlzDVertex2 *pos:	Pointer to first then second
*					positions.
************************************************************************/
static void	WlzGMModelMatchVtxG2D(WlzGMModel *model,
				       WlzGMEdgeT **matchET,
				       WlzDVertex2 searchDist,
				       int *searchIdx,
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
  if(searchIdx)
  {
    WlzGMModelMatchVtxGIdx2D(model, matchV, searchDist, searchIdx, pos);
  }
  else
  {
    WlzGMModelMatchVtxGLst2D(model, matchV, pos);
  }
  /* Now find the edge topology element(s) that are directed away from these
   * verticies. */
  for(idD = 0; idD < 2; ++idD)			 /* For each matched vertex. */
  {
    /* Check to see if this vertex has more than one vertex topology element.
     * If it does then more than one edge shares this vertex and the search
     * for the next edge topology elements becomes more complicated. */
    if((vertex = *(matchV + idD)) != NULL)
    {
      vertexT = vertex->top;
      eTS = vertexT->parent.edgeT;
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
	  (void )WlzGMVertexGetG2D(eTB->vertexT->vertex, &tPos);
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
	    (void )WlzGMVertexGetG2D(eTN->vertexT->vertex, &tPos);
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

/************************************************************************
* Function:	WlzGMModelMatchVtxGIdx2D
* Returns:	void
* Purpose:	Attempts to find verticies which match the two given
*		double precision positions using the given search
*		indicies. This algorithm will only work if the verticies
*		are allocated (and thus have increasing indicies) in
*		WLZ_RASTERDIR_ILIC order.
*		Verticies are matched by searching form the second of
*		the two search indicies (*(searchIdx + 1)) backwards
*		until either both verticies have been matched of
*		one of the verticies is out of range (searchDistSq).
*		The search then looks from the last vertex matched on
*		the previous vertex line forwards
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMVertex *matchV:	Array for return of matched
*					vertex pointers.
*		WlzDVertex2 searchDist:	Maximum search distance ignored,
*					if search indicies pointer is NULL.
*		int *searchIdx:		Last line vertex search index,
*					ignored if NULL.
*		WlzDVertex2 *pos:	Pointer to first then second
*					positions.
************************************************************************/
static void	WlzGMModelMatchVtxGIdx2D(WlzGMModel *model,
					 WlzGMVertex **matchV,
					 WlzDVertex2 searchDist,
					 int *searchIdx,
					 WlzDVertex2 *pos)
{
  int		idD,
  		idS,
		nxtIdx;
  int		doSearch[2];
  WlzDVertex2	cmp;
  AlcVector	*vecV;
  WlzGMVertex	*tV0;

  *matchV = *(matchV + 1) = NULL;
  if(model->res.vertex.nxtIdx > 1)       /* Don't look if this is the first! */
  {
    doSearch[0] = 1;
    doSearch[1] = 1;
    vecV = model->res.vertex.vec;
    /* Search from the last vertex on the current line back for the given
     * search distance, but stop searching if both verticies are matched. */
    idS = nxtIdx = model->res.vertex.nxtIdx - 1;
    do
    {
      idD = 0;
      tV0 = (WlzGMVertex  *)AlcVectorItemGet(vecV, idS);
      if(tV0->idx > 0)
      {
	while(idD < 2)
	{
	  if(doSearch[idD])
	  {
	    cmp = WlzGMVertexCmp2D(tV0, *(pos + idD));
	    doSearch[idD] = (cmp.vtY > -searchDist.vtY);
	    if(doSearch[idD])
	    {
	      if(((cmp.vtX * cmp.vtX) +
	          (cmp.vtY * cmp.vtY)) <= WLZ_GM_TOLERANCE)
	      {
		doSearch[idD] = 0;
		*(matchV + idD) = tV0;
	      }
	    }
	  }
	  ++idD;
	}
      }
      --idS;
    } while((idS > 0) && (doSearch[0] || doSearch[1]));
    /* Search forwards from the last vertex matched on the last line of
     * verticies for the given search distance for matching, but stop
     * searching if both verticies are matched. */
    doSearch[0] = *(matchV + 0) == NULL;
    doSearch[1] = *(matchV + 1) == NULL;
    if((doSearch[0] || doSearch[1]) && (*searchIdx > 1))
    {
      idS = *searchIdx;
      do
      {
	idD = 0;
	tV0 = (WlzGMVertex  *)AlcVectorItemGet(vecV, idS);
	if(tV0->idx > 0)
	{
	  while(idD < 2)
	  { 
	    if(doSearch[idD])
	    {
	      cmp = WlzGMVertexCmp2D(tV0, *(pos + idD));
	      doSearch[idD] = (cmp.vtX < searchDist.vtX) &&
			      (cmp.vtY < searchDist.vtY);
	      if(doSearch[idD])
	      {
		if(((cmp.vtX * cmp.vtX) +
		    (cmp.vtY * cmp.vtY)) <= WLZ_GM_TOLERANCE)
		{
		  doSearch[idD] = 0;
		  *(matchV + idD) = tV0;
		  *searchIdx = idS;      /* Update last matched on last line */
		}
	      }
	    }
	    ++idD;
	  }
	}
	++idS;
      } while((idS < nxtIdx) && (doSearch[0] || doSearch[1]));
    }
  }
}

/************************************************************************
* Function:	WlzGMModelMatchVtxGLst2D
* Returns:	void
* Purpose:	Attempts to find verticies which match the two given
*		double precision positions using a hierarchical
*		exhaustive search.
*		This search doesn't depend on the verticies being in
*		any particular order and is slower than the index
*		search method.
*		The basic algorithm is:
*		  for each shell of the model
*		    if either vertex is inside shell's bounding box
*		      for each loop of the shell
*			if either vertex is inside loop's bounding box
*			  check each edgeT's vertexT's vertex geometry
*			end if
*		      end for
*		    end if
*                 end for
* Global refs:	-
* Parameters:	WlzGMModel *model:	Model with resources.
*		WlzGMVertex **matchV:	Array for return of matched
*					vertex pointers.
*		WlzDVertex2 *pos:	Pointer to first then second
*					point positions.
************************************************************************/
static void	WlzGMModelMatchVtxGLst2D(WlzGMModel *model,
					  WlzGMVertex **matchV,
					  WlzDVertex2 *pos)
{
  double	tD0;
  unsigned int	searchMsk,
  		loopMsk,
  		shellMsk;
  WlzGMVertex	*tstV;
  WlzGMEdgeT	*tstET,
  		*firstET;
  WlzGMLoop	*tstL;
  WlzGMLoopT	*tstLT,
  		*firstLT;
  WlzGMShell	*tstS,
  		*firstS;

  *matchV = *(matchV + 1) = NULL;
  if(model->res.vertex.nxtIdx > 1)       /* Don't look if this is the first! */
  {
    searchMsk = 0x03;
    *(matchV + 0) = *(matchV + 1) = NULL;
    if(model && ((firstS = model->child) != NULL))
    {
      /* For each shell of model */
      tstS = firstS;
      do
      {
	shellMsk = (((WlzGMShellGInBB2D(tstS, *(pos + 1)) != 0) << 1) ||
		     (WlzGMShellGInBB2D(tstS, *(pos + 1)) != 0)) & searchMsk;
	if((shellMsk != 0) && ((firstLT = tstS->child.loopT) != NULL))
	{
	  /* For each loop of the shell */
	  tstLT = firstLT;
	  do
	  {
	    tstL = tstLT->loop;
	    loopMsk = (((WlzGMLoopGInBB2D(tstL, *(pos + 1)) != 0) << 1) ||
			(WlzGMLoopGInBB2D(tstL, *(pos + 1)) != 0)) & searchMsk;
	    if((loopMsk != 0) && ((firstET = tstLT->edgeT) != NULL))
	    {
	      /* For each edge of the loop */
	      tstET = firstET;
	      do
	      {
		/* Test vertex */
		tstV = tstET->vertexT->vertex;
		if((shellMsk & 1) != 0)
		{
		  tD0 = WlzGMVertexDistSq2D(tstV, *(pos + 0));
		  if((tD0 <= WLZ_GM_TOLERANCE) && (tD0 >= -(WLZ_GM_TOLERANCE)))
		  {
		    *(matchV + 0) = tstV;
		    loopMsk &= ~1;
		    shellMsk &= ~1;
		    searchMsk &= ~1;
		  }
		}
		if((shellMsk & 2) != 0)
		{
		  tD0 = WlzGMVertexDistSq2D(tstV, *(pos + 1));
		  if((tD0 <= WLZ_GM_TOLERANCE) && (tD0 >= -(WLZ_GM_TOLERANCE)))
		  {
		    *(matchV + 1) = tstV;
		    loopMsk &= ~2;
		    shellMsk &= ~2;
		    searchMsk &= ~2;
		  }
		}
		tstET = tstET->next;
	      }
	      while(loopMsk && (tstET->idx != firstET->idx));
	    }
	    tstLT = tstLT->next;
	  } while(shellMsk && (tstLT->idx != firstLT->idx));
	}
	tstS = tstS->next;
      } while(searchMsk && (tstS->idx != firstS->idx));
    }
  }
}

