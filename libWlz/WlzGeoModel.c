#pragma ident "MRC HGU $Id$"
/*!
* \file		WlzGeoModel.c
* \author	Bill Hill
* \date		February 2000
* \version	$Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Basic operators for manipulating Woolz geometric models.
*		These can be either planar graphs or 3D surfaces, with
*		the surfaces being either manifold or non-manifold.
* \ingroup      WlzGeoModel
* \todo
*		\li Functions for finding the minimum distance between two
*		verticies are only partialy implemented.
* 		\li The element deletion functions are only partly written
*		and only partly work for 2D models.
*		\li Radial edges are not checked for in WlzGMModelDeleteS().
*		\li Unlinking edge topology elements using WlzGMEdgeTUnlink()
*		may not be valid if edge topology elements aren't unlinked
*		in pairs for terminal edge elements.
*		\li 3D simplex construction is not complete. If all three
*		verticies and all three edges are within the model, need to
*		check if the three verticies are in a common loop and then if
*		ther're not add a new loop, 2 loopT's, 6 edgeT's and 6
*		vertexT's.
*		\li WlzGMLoopTFindShell() code for loopTs connected by a
*		vertex not implemented yet.
* \bug          None known.
*/
#include <Wlz.h>
#include <float.h>
#include <limits.h>
#include <string.h>

static WlzGMShell 	*WlzGMLoopTFindShell(
			  WlzGMLoopT *gLT);
static void	 	WlzGMShellMergeG(
			  WlzGMShell *shell0,
			  WlzGMShell *shell1);
static void	 	WlzGMShellMergeT(
			  WlzGMShell *s0,
			  WlzGMShell *s1);
static void		WlzGMModelMatchEdgeTG2D(
			  WlzGMModel *model,
			  WlzGMEdgeT **matchET,
			  WlzDVertex2 *pos);
static unsigned int 	WlzGMHashPos2D(
			  WlzDVertex2 pos);
static unsigned int 	WlzGMHashPos3D(
			  WlzDVertex3 pos);
static double		WlzGMVertexShellDist2(
			  WlzGMVertex *v0,
			  WlzGMVertex *v1,
			  WlzErrorNum *dstErr);
static double		WlzGMVertexShellDist2CommonLT(
			  WlzGMVertexT *gVT0,
			  WlzGMVertexT *gVt1);
static WlzErrorNum 	WlzGMModelDeleteE3D(
			  WlzGMModel *model,
			  WlzGMEdge *dE);
static WlzErrorNum 	WlzGMModelDeleteE2D(
			  WlzGMModel *model,
			  WlzGMEdge *dE);
static WlzErrorNum 	WlzGMModelDeleteE2D0V1L(
			  WlzGMModel *model,
			  WlzGMEdge *dE);
static WlzErrorNum 	WlzGMModelDeleteE2D1V1L(
			  WlzGMModel *model,
			  WlzGMEdge *dE);
static WlzErrorNum 	WlzGMModelDeleteE2D2V1L(
			  WlzGMModel *model,
			  WlzGMEdge *dE);
static WlzErrorNum 	WlzGMModelDeleteE2D2V2L(
			  WlzGMModel *model,
			  WlzGMEdge *dE);
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
static WlzErrorNum 	WlzGMModelExtend3V3E1S3D(
			  WlzGMModel *model,
			  WlzGMEdge **cE);
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
static void		WlzGMLoopTSetAdjT(
			  WlzGMLoopT *gLT,
			  WlzGMShell *gS);
static void		WlzGMShellSetT(
			  WlzGMShell *gS);
static void             WlzGMLoopTSetT(
                          WlzGMLoopT *eLT);
static void             WlzGMModelAddVertex(
                          WlzGMModel *model,
                          WlzGMVertex *nV);
static void		WlzGMModelRemVertex(
			  WlzGMModel *model,
			  WlzGMVertex *dV);
static void		WlzGMElmMarkFree(
			  int *idxP);

/* Resource callback function list manipulation. */

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Add a resource allocation callback to the models resources.
* \param	model			The model.
* \param	fn			The callback function.
* \param	data			The callback data to be passed on
*					to the callback function.
*/
WlzErrorNum	WlzGMModelAddResCb(WlzGMModel *model, WlzGMCbFn fn,
				   void *data)
{
  WlzGMCbEntry	*newCbE = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model && fn)
  {
    if((newCbE = (WlzGMCbEntry *)AlcMalloc(sizeof(WlzGMCbEntry))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      newCbE->fn = fn;
      newCbE->data = data;
      newCbE->next = model->res.callbacks;
      model->res.callbacks = newCbE;
    }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Removes a resource allocation callback from the models
*		resources. Both the callback function and callback data
*		must have the same values as when the callback was added
*		since they are used to find the corresponding callback
*		entry.
* \param	model			The model.
* \param	fn			The callback function.
* \param	data			The callback data.
*/
void		WlzGMModelRemResCb(WlzGMModel *model, WlzGMCbFn fn,
				   void *data)
{
  WlzGMCbEntry	*cBE0,
  		*cBE1;

  if(model && fn && ((cBE1 = model->res.callbacks) != NULL))
  {
    cBE0 = NULL;
    while((fn != cBE1->fn) && (data != cBE1->data) && cBE1->next)
    {
      cBE0 = cBE1;
      cBE1 = cBE1->next;
    }
    if((fn == cBE1->fn) && (data == cBE1->data))
    {
      if(cBE0)
      {
        cBE0->next = cBE1->next;
      }
      else
      {
        model->res.callbacks = NULL;
      }
      AlcFree(cBE1);
    }
  }
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Calls all resource allocation callbacks passing on the
*		given model, element and reason along with the data from
*		the callback entry.
* \param	model			The model.
* \param	elm			The callback function.
* \param	reason			The callback data.
*/
static void	WlzGMModelCallResCb(WlzGMModel *model, WlzGMElemP elm,
				    WlzGMCbReason reason)
{
  WlzGMCbEntry	*cBE;

  if(model && elm.core && ((cBE = model->res.callbacks) != NULL))
  {
    do
    {
      (*(cBE->fn))(model, elm, reason, cBE->data);
    } while((cBE = cBE->next) != NULL);
  }
}

/* Creation  of geometric modeling elements */

/*!
* \return				New empty model.
* \ingroup      WlzGeoModel
* \brief	Creates an empty non-manifold geometry model.
* \param	modType			Type of model to create.
* \param	blkSz			Resource block size, used for
*					allocating storage for model
*					elements. A default size is
*					used if <= 0.
* \param	vHTSz			Vertex matching hash table size,
*                                       A default size is used if <= 0.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMModel	*WlzGMModelNew(WlzGMModelType modType,
			       int blkSz, int vHTSz, WlzErrorNum *dstErr)
{
  int		dim = 0;
  WlzGMModel	*model = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  unsigned int	vertexGSz,
		shellGSz;
  const unsigned int defBlkSz = 1024;
  const unsigned int defVHTSz2D = 16 * 1024,
  		     defVHTSz3D = 1024 * 1024;

  if(blkSz <= 0)
  {
    blkSz = defBlkSz;
  }
  switch(modType)
  {
    case WLZ_GMMOD_2I:
      dim = 2;
      vertexGSz = sizeof(WlzGMVertexG2I);
      shellGSz = sizeof(WlzGMShellG2I);
      break;
    case WLZ_GMMOD_2D:
      dim = 2;
      vertexGSz = sizeof(WlzGMVertexG2D);
      shellGSz = sizeof(WlzGMShellG2D);
      break;
    case WLZ_GMMOD_3I:
      dim = 3;
      vertexGSz = sizeof(WlzGMVertexG3I);
      shellGSz = sizeof(WlzGMShellG3I);
      break;
    case WLZ_GMMOD_3D:
      dim = 3;
      vertexGSz = sizeof(WlzGMVertexG3D);
      shellGSz = sizeof(WlzGMShellG3D);
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if((errNum == WLZ_ERR_NONE) && (vHTSz <= 0))
  {
    vHTSz = (dim == 3)? defVHTSz3D: defVHTSz2D;
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
       ((model->res.face.vec = AlcVectorNew(1, sizeof(WlzGMFace),
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

/*!
* \return				New empty shell.
* \ingroup      WlzGeoModel
* \brief	Creates an empty shell with a shell geometry element
*               for the model.
* \param	model			The geometric model.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMShell	*WlzGMModelNewS(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resS,
  		*resSG;
  WlzGMShell	*shell = NULL;
  WlzGMElemP	elm;
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
    elm.shell = shell;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
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

/*!
* \return				New face.
* \ingroup      WlzGeoModel
* \brief	Creates an empty face.
* \param	model			Parent model.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMFace	*WlzGMModelNewF(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resF;
  WlzGMFace	*face = NULL;
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Get a new face element from the model */
    resF = &(model->res.face);
    if((face = (WlzGMFace *)
     	       (AlcVectorExtendAndGet(resF->vec, resF->numIdx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ++(resF->numElm);
    face->type = WLZ_GMELM_FACE;
    face->idx = (resF->numIdx)++;
    elm.face = face;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(face);
}

/*!
* \return				New loop topology element.
* \ingroup      WlzGeoModel
* \brief	Creates a loop topology element.
* \param	model			Parent model.
* \param	dstErr			Destination error pointer, may
*					be null.
*/
WlzGMLoopT	*WlzGMModelNewLT(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resLT;
  WlzGMLoopT	*loopT = NULL;
  WlzGMElemP	elm;
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
    elm.loopT = loopT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(loopT);
}

/*!
* \return				New disk topology element.
* \ingroup      WlzGeoModel
* \brief	Creates a new disk topology element.
* \param	model			Model with resources.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMDiskT     *WlzGMModelNewDT(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resDT;
  WlzGMDiskT	*diskT = NULL;
  WlzGMElemP	elm;
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
    elm.diskT = diskT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
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

/*!
* \return				New edge.
* \ingroup      WlzGeoModel
* \brief	Creates a new edge.
*               The edge geometry element only has it's index set
*               to a meaningful value.
* \param	model			Model with resources.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMEdge      *WlzGMModelNewE(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resE;
  WlzGMEdge	*edge = NULL;
  WlzGMElemP	elm;
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
    elm.edge = edge;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(edge);
}

/*!
* \return				New edge topology element.
* \ingroup      WlzGeoModel
* \brief	Creates a new edge topology element.
* \param	model			Model with resources.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMEdgeT     *WlzGMModelNewET(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resET;
  WlzGMEdgeT	*edgeT = NULL;
  WlzGMElemP	elm;
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
    elm.edgeT = edgeT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(edgeT);
}

/*!
* \return				New vertex.
* \ingroup      WlzGeoModel
* \brief	Creates a new vertex and a vertex geometry element.
*               The vertex geometry element only has it's index set
*               to a meaningful value.
* \param	model			Model with resources.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMVertex      *WlzGMModelNewV(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resV,
  		*resVG;
  WlzGMVertex	*vertex = NULL;
  WlzGMElemP	elm;
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
    elm.vertex = vertex;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
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

/*!
* \return				New vertex topology element.
* \ingroup      WlzGeoModel
* \brief	Creates a new vertex topology element.
* \param	model			Model with resources.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMVertexT      *WlzGMModelNewVT(WlzGMModel *model, WlzErrorNum *dstErr)
{
  WlzGMResource	*resVT;
  WlzGMVertexT	*vertexT = NULL;
  WlzGMElemP	elm;
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
    elm.vertexT = vertexT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_NEW);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vertexT);
}

/*!
* \return				New copy of the given model.
* \ingroup      WlzGeoModel
* \brief	Copies the given model.
*		Because all unused elements are squeezed out, the indicies of
*		the elements may not be equal in the two models.
* \param	gM			Given model.
* \param	dstErr			Destination error pointer, may
*					be null.
*/
WlzGMModel 	*WlzGMModelCopy(WlzGMModel *gM, WlzErrorNum *dstErr)
{
  int		dim,
  		gIdx,
  		nIdx,
		tIdx,
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
    dim = WlzGMModelGetDimension(gM, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
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
				WLZ_GMELMFLG_FACE |
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
        if((nElmP.core = (WlzGMCore *)
			 AlcVectorExtendAndGet(nM->res.vertex.vec,
					       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = WLZ_GMELM_VERTEX;
	  nElmP.core->idx = nIdx;
	  tIdx = *(resIdxTb->diskT.idxLut + gElmP.vertex->diskT->idx);
	  nElmP.vertex->diskT = (WlzGMDiskT *)
	  			AlcVectorExtendAndGet(nM->res.diskT.vec, tIdx);
	  tIdx = *(resIdxTb->vertexG.idxLut + gElmP.vertex->geo.core->idx);
	  nElmP.vertex->geo.core = (WlzGMCore *)
	  			   AlcVectorExtendAndGet(nM->res.vertexG.vec,
						         tIdx);
	  if((gElmP.vertex->next) && (gElmP.vertex->next->idx >= 0))
	  {
	    tIdx = *(resIdxTb->vertex.idxLut + gElmP.vertex->next->idx);
	    nElmP.vertex->next = (WlzGMVertex *)
				 AlcVectorExtendAndGet(nM->res.vertex.vec,
						       tIdx);
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
	if((nElmP.core = (WlzGMCore *)
			 AlcVectorExtendAndGet(nM->res.vertexG.vec,
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
	if((nElmP.core = (WlzGMCore *)
			 AlcVectorExtendAndGet(nM->res.vertexT.vec,
			 		       nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  tIdx = *(resIdxTb->vertexT.idxLut + gElmP.vertexT->next->idx);
	  nElmP.vertexT->next = (WlzGMVertexT *)
	  			AlcVectorExtendAndGet(nM->res.vertexT.vec,
						      tIdx);
          tIdx = *(resIdxTb->vertexT.idxLut + gElmP.vertexT->prev->idx);
	  nElmP.vertexT->prev = (WlzGMVertexT *)
	  			AlcVectorExtendAndGet(nM->res.vertexT.vec,
						      tIdx);
	  tIdx = *(resIdxTb->diskT.idxLut + gElmP.vertexT->diskT->idx);
	  nElmP.vertexT->diskT = (WlzGMDiskT *)
	  			 AlcVectorExtendAndGet(nM->res.diskT.vec,
						       tIdx);
	  tIdx = *(resIdxTb->edgeT.idxLut + gElmP.vertexT->parent->idx);
	  nElmP.vertexT->parent = (WlzGMEdgeT *)
	  			  AlcVectorExtendAndGet(nM->res.edgeT.vec,
						 	tIdx);
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
	  tIdx = *(resIdxTb->diskT.idxLut + gElmP.diskT->next->idx);
	  nElmP.diskT->next = (WlzGMDiskT *)
	  		      AlcVectorExtendAndGet(nM->res.diskT.vec, tIdx);
	  tIdx = *(resIdxTb->diskT.idxLut + gElmP.diskT->prev->idx);
	  nElmP.diskT->prev = (WlzGMDiskT *)
			      AlcVectorExtendAndGet(nM->res.diskT.vec, tIdx);
	  tIdx = *(resIdxTb->vertex.idxLut + gElmP.diskT->vertex->idx);
	  nElmP.diskT->vertex = (WlzGMVertex *)
	  			 AlcVectorExtendAndGet(nM->res.vertex.vec,
				 		       tIdx);
	  tIdx = *(resIdxTb->vertexT.idxLut + gElmP.diskT->vertexT->idx);
	  nElmP.diskT->vertexT = (WlzGMVertexT *)
	  			  AlcVectorExtendAndGet(nM->res.vertexT.vec,
						        tIdx);
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
	  tIdx = *(resIdxTb->edgeT.idxLut + gElmP.edge->edgeT->idx);
	  nElmP.edge->edgeT = (WlzGMEdgeT *)
	  		      AlcVectorExtendAndGet(nM->res.edgeT.vec, tIdx);
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
	  tIdx = *(resIdxTb->edgeT.idxLut + gElmP.edgeT->next->idx);
	  nElmP.edgeT->next = (WlzGMEdgeT *)
	  		       AlcVectorExtendAndGet(nM->res.edgeT.vec, tIdx);
	  tIdx  = *(resIdxTb->edgeT.idxLut + gElmP.edgeT->prev->idx);
	  nElmP.edgeT->prev = (WlzGMEdgeT *)
	  		       AlcVectorExtendAndGet(nM->res.edgeT.vec, tIdx);
	  tIdx = *(resIdxTb->edgeT.idxLut + gElmP.edgeT->opp->idx);
	  nElmP.edgeT->opp = (WlzGMEdgeT *)
	  		     AlcVectorExtendAndGet(nM->res.edgeT.vec, tIdx);
	  tIdx = *(resIdxTb->edgeT.idxLut + gElmP.edgeT->rad->idx);
	  nElmP.edgeT->rad = (WlzGMEdgeT *)
	  		       AlcVectorExtendAndGet(nM->res.edgeT.vec, tIdx);
	  tIdx = *(resIdxTb->edge.idxLut + gElmP.edgeT->edge->idx);
	  nElmP.edgeT->edge = (WlzGMEdge *)
	  		      AlcVectorExtendAndGet(nM->res.edge.vec, tIdx);
	  tIdx = *(resIdxTb->vertexT.idxLut + gElmP.edgeT->vertexT->idx);
	  nElmP.edgeT->vertexT = (WlzGMVertexT *)
	  		         AlcVectorExtendAndGet(nM->res.vertexT.vec,
						       tIdx);
	  tIdx = *(resIdxTb->loopT.idxLut + gElmP.edgeT->parent->idx);
	  nElmP.edgeT->parent = (WlzGMLoopT *)
	  		        AlcVectorExtendAndGet(nM->res.loopT.vec, tIdx);

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
  if((errNum == WLZ_ERR_NONE) && (dim == 3))
  {
    /* Copy faces */
    gIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (gIdx < gM->res.face.numIdx))
    {
      gElmP.core = (WlzGMCore *)AlcVectorItemGet(gM->res.face.vec, gIdx);
      if((gElmP.core != NULL) && (gElmP.core->idx >= 0))
      {
	nIdx = *(resIdxTb->face.idxLut + gIdx);
	if((nElmP.core = (WlzGMCore *)AlcVectorExtendAndGet(nM->res.face.vec,
						            nIdx)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nElmP.core->type = gElmP.core->type;
	  nElmP.core->idx = nIdx;
	  tIdx = *(resIdxTb->loopT.idxLut + gElmP.face->loopT->idx);
	  nElmP.face->loopT = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec, tIdx);
	  if(nElmP.face->loopT == NULL)
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
	  tIdx = *(resIdxTb->loopT.idxLut + gElmP.loopT->next->idx);
	  nElmP.loopT->next = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec, tIdx);
	  tIdx = *(resIdxTb->loopT.idxLut + gElmP.loopT->prev->idx);
	  nElmP.loopT->prev = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec, tIdx);
	  tIdx = *(resIdxTb->loopT.idxLut + gElmP.loopT->opp->idx);
	  nElmP.loopT->opp = (WlzGMLoopT *)
	  		      AlcVectorExtendAndGet(nM->res.loopT.vec, tIdx);
	  if(dim == 2)
	  {
	    nElmP.loopT->face = NULL;
	  }
	  else
	  {
	    tIdx = *(resIdxTb->face.idxLut + gElmP.loopT->face->idx);
	    nElmP.loopT->face = (WlzGMFace *)
				AlcVectorExtendAndGet(nM->res.face.vec, tIdx);
	  }
	  tIdx = *(resIdxTb->edgeT.idxLut + gElmP.loopT->edgeT->idx);
	  nElmP.loopT->edgeT = (WlzGMEdgeT *)
	  		      AlcVectorExtendAndGet(nM->res.edgeT.vec, tIdx);
	  tIdx = *(resIdxTb->shell.idxLut + gElmP.loopT->parent->idx);
	  nElmP.loopT->parent = (WlzGMShell *)
				AlcVectorExtendAndGet(nM->res.shell.vec, tIdx);
	  if((nElmP.loopT->next == NULL) ||
	     (nElmP.loopT->prev == NULL) ||
	     (nElmP.loopT->opp == NULL) ||
	     ((dim != 2) && (nElmP.loopT->face == NULL)) ||
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
	  tIdx = *(resIdxTb->shell.idxLut + gElmP.shell->next->idx);
	  nElmP.shell->next = (WlzGMShell *)
	  		      AlcVectorExtendAndGet(nM->res.shell.vec, tIdx);
          tIdx = *(resIdxTb->shell.idxLut + gElmP.shell->prev->idx);
	  nElmP.shell->prev = (WlzGMShell *)
	  		      AlcVectorExtendAndGet(nM->res.shell.vec, tIdx);
          tIdx = *(resIdxTb->shellG.idxLut + gElmP.shell->geo.core->idx);
	  nElmP.shell->geo.core = (WlzGMCore *)
	  		      AlcVectorExtendAndGet(nM->res.shellG.vec, tIdx);
	  tIdx = *(resIdxTb->loopT.idxLut + gElmP.shell->child->idx);
	  nElmP.shell->child = (WlzGMLoopT *)
	  		       AlcVectorExtendAndGet(nM->res.loopT.vec, tIdx);
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
    nM->res.face.numElm = nM->res.face.numIdx = gM->res.face.numElm;
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Free's a geometric model and it's elements.
* \param	model			The shell to free.
*/
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
    (void )AlcVectorFree(model->res.face.vec);
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks an element as free by setting it's index to
*               a -ve value, the actual choice of value aids debugging.
* \param	idxP			Pointer to elements index.
*/
static void	WlzGMElmMarkFree(int *idxP)
{
  if(*idxP >= 0)
  {
    *idxP = -(*idxP + 1);
  }
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks a shell and it's geometry as invalid and suitable
*               for reclaiming.
* \param	model			Model with resources.
* \param	shell			Shell to free.
*/
WlzErrorNum	WlzGMModelFreeS(WlzGMModel *model, WlzGMShell *shell)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (shell == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (shell->idx >= 0))
  {
    elm.shell = shell;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
    /* Can't really free the shell so just mark it and it's geometry element
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks a face as invalid and suitable for reclaiming.
* \param	model			Model with resources.
* \param	face			Face to free.
*/
WlzErrorNum	WlzGMModelFreeF(WlzGMModel *model, WlzGMFace *face)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (face == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (face->idx >= 0))
  {
    elm.face = face;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
    /* Can't really free the face so just mark it as invalid by making
     * the indicies < 0. */
    WlzGMElmMarkFree(&(face->idx));
    --(model->res.face.numElm);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks a loop topology as invalid and suitable
*               for reclaiming.
* \param	model			Model with resources.
* \param	loopT			LoopT to free.
*/
WlzErrorNum	WlzGMModelFreeLT(WlzGMModel *model, WlzGMLoopT *loopT)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (loopT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (loopT->idx >= 0))
  {
    elm.loopT = loopT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
    /* Can't really free the loopT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(loopT->idx));
    --(model->res.loopT.numElm);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks an disk topology as invalid and suitable
*               for reclaiming.
* \param	model			Model with resources.
* \param	diskT			DiskT to free.
*/
WlzErrorNum	WlzGMModelFreeDT(WlzGMModel *model, WlzGMDiskT *diskT)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (diskT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (diskT->idx >= 0))
  {
    elm.diskT = diskT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
    /* Can't really free the diskT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(diskT->idx));
    --(model->res.diskT.numElm);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks an edge as invalid and suitable for reclaiming.
* \param	model			Model with resources.
* \param	edge			Edge to free.
*/
WlzErrorNum	WlzGMModelFreeE(WlzGMModel *model, WlzGMEdge *edge)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (edge == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (edge->idx >= 0))
  {
    elm.edge = edge;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
    /* Can't really free the edge so just mark it as invalid by making
     * the indicies < 0. */
    WlzGMElmMarkFree(&(edge->idx));
    --(model->res.edge.numElm);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks an edge topology as invalid and suitable for
*		reclaiming.
* \param	model			Model with resources.
* \param	edgeT			EdgeT to free.
*/
WlzErrorNum	WlzGMModelFreeET(WlzGMModel *model, WlzGMEdgeT *edgeT)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (edgeT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (edgeT->idx >= 0))
  {
    elm.edgeT = edgeT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
    /* Can't really free the edgeT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(edgeT->idx));
    --(model->res.edgeT.numElm);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks a vertex and it's geometry as invalid and suitable
*               for reclaiming. Watch out that vertices remain in the hash
*               table after being free'd.
* \param	model			Model with resources.
* \param	vertex			Vertex to free.
*/
WlzErrorNum	WlzGMModelFreeV(WlzGMModel *model, WlzGMVertex *vertex)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (vertex == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (vertex->idx >= 0))
  {
    elm.vertex = vertex;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Marks a vertex topology as invalid and suitable for
*		reclaiming.
* \param	model			Model with resources.
* \param	vertexT			VertexT to free.
*/
WlzErrorNum	WlzGMModelFreeVT(WlzGMModel *model, WlzGMVertexT *vertexT)
{
  WlzGMElemP	elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model == NULL) || (vertexT == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMModelTypeValid(model->type);
  }
  if((errNum == WLZ_ERR_NONE) && (vertexT->idx >= 0))
  {
    elm.vertexT = vertexT;
    WlzGMModelCallResCb(model, elm, WLZ_GMCB_FREE);
    /* Can't really free the vertexT so just mark it as invalid by making
     * the index < 0. */
    WlzGMElmMarkFree(&(vertexT->idx));
    --(model->res.vertexT.numElm);
  }
  return(errNum);
}

/* Deletion of geometric modeling elements along with children */

/*!
* \return	Woolz error code.
* \ingroup      WlzGeoModel
* \brief        Deletes a shell by unlinking it and then freeing it and
*               all it's children.
* \param        model                   Model with resources.
* \param        dS                      Shell to delete.
*/
WlzErrorNum	WlzGMModelDeleteS(WlzGMModel *model, WlzGMShell *dS)
{
  WlzGMVertexT	*cVT;
  WlzGMEdgeT	*cET,
  		*fET;
  WlzGMLoopT	*cLT,
		*fLT;

  if(model && dS && (dS->idx >= 0))
  {
    WlzGMShellUnlink(dS);
    if((cLT = fLT = dS->child) != NULL)
    {
      do 				     /* For each loopT of the shell. */
      {
	cLT = cLT->next;
	if(cLT->idx >= 0)
	{
	  cET = fET = cLT->edgeT;
	  do 				     /* For each edgeT of the loopT. */
	  {
	    /* Don't need to follow the opposite or radial edgeT's because
	     * each of these will be part of some loopT within this shell.
	     * All the element freeing functions check that the element's
	     * is valid, so no need to check for NULL element pointers or
	     * if already freed. */
	    cET = cET->next;
	    cVT = cET->vertexT;
	    (void )WlzGMModelFreeV(model,  cVT->diskT->vertex);
	    (void )WlzGMModelFreeDT(model,  cVT->diskT);
	    (void )WlzGMModelFreeVT(model, cVT);
	    (void )WlzGMModelFreeE(model, cET->edge);
	    (void )WlzGMModelFreeET(model, cET);
	  } while(cET != fET);
	  (void )WlzGMModelFreeF(model, cLT->face);
	  (void )WlzGMModelFreeLT(model, cLT);
	}
      } while(cLT != fLT);
    }
    WlzGMModelFreeS(model, dS);
  }
  return(WLZ_ERR_NONE);
}

/*!
* \return	Woolz error code.
* \brief	Deletes a face along with all the elements which depend
*		on it. This may change the geometry and/or the topology
*		of the parent shell.
* \todo		This function is yet to be written.
* \param	model			Model with resources.
* \param	dF			Face to delete.
*/
WlzErrorNum	WlzGMModelDeleteF(WlzGMModel *model, WlzGMFace *dF)
{
  /* TODO */
  return(WLZ_ERR_UNIMPLEMENTED);
}

/*!
* \return	Woolz error code.
* \brief	Deletes a edge along with all the elements which depend
*		on it. This may change the geometry and/or the topology
*		of the parent shell. On return all shell geometries are
*		valid.
* \todo		Deletion of edges has only been implemented for 2D models.
* \param	model			Model with resources.
* \param	dE			Edge to delete.
*/
WlzErrorNum	WlzGMModelDeleteE(WlzGMModel *model, WlzGMEdge *dE)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model && dE && (dE->idx >= 0))
  {
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	errNum = WlzGMModelDeleteE2D(model, dE);
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	errNum = WlzGMModelDeleteE3D(model, dE);
	break;
      default:
	errNum = WLZ_ERR_GMELM_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Deletes a edge along with all the elements which depend
*		on it. This may change the geometry and/or the topology
*		of the parent shell. Special case for 2D models.
*		The topology of the edge to be deleted is used to classify
*		the problem into one of four cases. This classification
*		is based on the number of the edge's verticies which connect
*		it to the rest of the shell and on the number of loop
*		topology elements that run along the edge. On return all
*		shell geometries are valid.
* \param	model			Model with resources.
* \param	dE			Edge to delete.
*/
static WlzErrorNum WlzGMModelDeleteE2D(WlzGMModel *model, WlzGMEdge *dE)
{
  int		nV,  /* Number of verticies on the edge connecting it to the
  	              * rest of the shell. */
  		nLT, /* Non-zero if the edge has two distinct edgeTs running
		      * along it. */
                con;
  WlzGMEdgeT	*tET0;
  WlzGMVertexT	*tVT0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tET0 = dE->edgeT;
  tVT0 = tET0->vertexT;
  nV = (tVT0 != tVT0->next);
  tVT0 = tET0->opp->vertexT;
  nV += (tVT0 != tVT0->next);
  nLT = tET0->parent != tET0->opp->parent;
  con = (nLT * 4) + nV;
  switch(con)
  {
    case 0:	/* nV = 0, nLT = 0: Edge is isolated in it's own shell. */
      errNum = WlzGMModelDeleteE2D0V1L(model, dE);
      break;
    case 1:	/* nV = 1, nLT = 0: Edge is a hair extending from it's shell,
    		 * attached to the shell at a single vertex. */
      errNum = WlzGMModelDeleteE2D1V1L(model, dE);
      break;
    case 2:	/* nV = 2, nLT = 0: Edge is a bridge which connects two
    		 * otherwise isolated shells. */
      errNum = WlzGMModelDeleteE2D2V1L(model, dE);
      break;
    case 6:	/* nV = 2, nLT = 1: Edge joins two loopTs within the shell. */
      errNum = WlzGMModelDeleteE2D2V2L(model, dE);
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Deletes a edge along with all the elements which depend
*		on it. This may change the geometry and/or the topology
*		of the parent shell. Special case for 3D models.
* \todo		This function is yet to be written.
* \param	model			Model with resources.
* \param	dE			Edge to delete.
*/
static WlzErrorNum WlzGMModelDeleteE3D(WlzGMModel *model, WlzGMEdge *dE)
{
  return(WLZ_ERR_UNIMPLEMENTED);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Deletes a vertex along with all the elements which depend
*		on it. All elements which depend on the vertex are unlinked
*		and then freed. If the vertex's parents were to depend
*		only on the vertex then they would be free'd too as
*		the Woolz geometric models can not hold an isolated
*		vertex.
*		The basic algorithm used is: Build a collection of edges
*		which use the vertex and then delete all edges in the
*		collection. If there are no edges (in the collection) then
*		just unlink and free the vertex.
*		The geometries of the existing and new shells will be
*		valid on return.
* \todo		This only works for 2D models, it needs to be extend to
*		3D models.
* \param	model			Model with resources.
* \param	dV			The vertex to delete.
*/
WlzErrorNum	WlzGMModelDeleteV(WlzGMModel *model, WlzGMVertex *dV)
{
  int		eCnt,
  		eMax;
  WlzGMVertexT	*tVT0,
  		*tVT1;
  WlzGMDiskT	*tDT0,
  		*tDT1;
  WlzGMEdge	**edgeCol = NULL;
  WlzGMShell	*pS;
  const int	eStp = 16;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model && dV && (dV->idx >= 0))
  {
    /* Allocate initial storage for the edge collection. With luck eStp
     * will have been chosen to be some small value that's > than the the
     * number of edges incident with most verticies and the edge collection
     * won't need reallocating. */
    eCnt = eMax = 0;
    /* Build a collection of all the edges that will be destroyed by deleting
     * the vertex. */
    tDT1 = tDT0 = dV->diskT;
    do
    {
      tDT1 = tDT1->next;
      tVT1 = tVT0 = tDT1->vertexT;
      do
      {
	tVT1 = tVT1->next;
	if(eCnt >= eMax)
	{
	  eMax = (eMax == 0)? eStp: eMax * 2;
	  if((edgeCol = (WlzGMEdge **)AlcRealloc(edgeCol,
					sizeof(WlzGMEdge *) * eMax)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *(edgeCol + eCnt++) = tVT1->parent->edge;
	}
      } while((errNum == WLZ_ERR_NONE) && (tVT1 != tVT0));
    } while((errNum == WLZ_ERR_NONE) && (tDT1 != tDT0));
    if(errNum == WLZ_ERR_NONE)
    {
      if(eCnt == 0)
      {
	/* Unlink the lone vertex (removing it from the hash table)
	 * and then free it. Currently this isn't needed because lone
	 * verticies aren't allowed in the models. */
      }
      else
      {
	/* Delete all the edges in the collection. */
	while(eCnt > 0)
	{
	  errNum = WlzGMModelDeleteE(model, *(edgeCol + --eCnt));
	}
      }
    }
    if(edgeCol)
    {
      AlcFree(edgeCol);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Deletes an isolated edge. Since the edge is in it's own
*		shell this is simply a matter of deleting the shell.
* \param	model			Model with resources.
* \param	dE			Edge to delete.
*/
static WlzErrorNum WlzGMModelDeleteE2D0V1L(WlzGMModel *model, WlzGMEdge *dE)
{
  WlzErrorNum	errNum;

  errNum = WlzGMModelDeleteS(model, dE->edgeT->parent->parent);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Deletes an edge connected to the rest of it's shell through
*		a single vertex of the edge.
*		On the shell has had it's geometry set.
* \verbatim
*
*                                              tVT0   tET0        
*                 -->   o------------------>   o------------------>   
*                                                      dE              
*                 ----@----------------------@----------------------@  
*                                                                      
*                 --o   <------------------o   <------------------o   
*                                                     tET1
*
* \endverbatim
* \param	model			Model with resources.
* \param	dE			Edge to delete.
*/
static WlzErrorNum WlzGMModelDeleteE2D1V1L(WlzGMModel *model, WlzGMEdge *dE)
{
  WlzGMLoopT	*cLT;
  WlzGMEdgeT	*tET0,
  		*tET1;
  WlzGMVertexT	*tVT0;
  WlzDVertex2	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tET0 = dE->edgeT;
  if(tET0->next != tET0->opp)
  {
    tET0 = tET0->opp;
  }
  tET1 = tET0->opp;
  tVT0 = tET0->vertexT;
  cLT = tET0->parent;
  if((cLT->edgeT == tET0) || (cLT->edgeT = tET1))
  {
    cLT->edgeT = tET0->prev;
  }
  tET0->prev->next = tET1->next;
  tET1->next->prev = tET0->prev;
  errNum = WlzGMVertexGetG2D(tVT0->diskT->vertex, &pos);
  if(errNum == WLZ_ERR_NONE)
  {
    WlzGMVertexTUnlink(tVT0);
    WlzGMModelFreeV(model, tET1->vertexT->diskT->vertex);
    WlzGMModelFreeDT(model, tET1->vertexT->diskT);
    WlzGMModelFreeVT(model, tET1->vertexT);
    WlzGMModelFreeE(model, dE);
    WlzGMModelFreeET(model, tET0);
    WlzGMModelFreeET(model, tET1);
    WlzGMModelFreeVT(model, tVT0);
    errNum = WlzGMShellDndateG2D(cLT->parent, pos);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Deletes an edge which connects two otherwise seperate shells
*		and in doing so creates a new shell.
*		On return both the current and new shell have had their
*		geometries set.
*		This is the most expensive of the edge deletions because
*		it changes both the topology and geometry of the model.
* \verbatim
*
*                 pET0          tVT0    tET0        
*                 ---------->   o------------------>   o------------
*                                        dE                          
*                 ------------@----------------------@--------------
*                                                                      
*                 ----------o   <------------------o   <------------
*                                       tET1       tVT1     pET1
*
* \endverbatim
* \param	model			Model with resources.
* \param	dE			Edge to delete.
*/
static WlzErrorNum WlzGMModelDeleteE2D2V1L(WlzGMModel *model, WlzGMEdge *dE)
{
  WlzGMEdgeT	*tET0,
  		*tET1,
		*pET0,
		*pET1,
		*nET0,
		*nET1;
  WlzGMVertexT	*tVT0,
  		*tVT1;
  WlzGMLoopT	*cLT,
  		*nLT;
  WlzGMShell	*cS,
  		*nS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tET0 = dE->edgeT;
  tET1 = tET0->opp;
  pET0 = tET0->prev;
  pET1 = tET1->prev;
  nET0 = tET0->next;
  nET1 = tET1->next;
  tVT0 = tET0->vertexT;
  tVT1 = tET1->vertexT;
  nET0->prev = pET1;
  pET1->next = nET0;
  nET1->prev = pET0;
  pET0->next = nET1;
  WlzGMVertexTUnlink(tVT0);
  WlzGMVertexTUnlink(tVT1);
  cLT = pET0->parent;
  nLT = WlzGMModelNewLT(model, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    cS = cLT->parent;
    nS = WlzGMModelNewS(model, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cLT->edgeT = pET0;
    nLT->edgeT = pET1;
    nLT->parent = nS;
    nLT->next = nLT->prev = nLT->opp = nLT;
    nS->child = nLT;
    nS->parent = model;
    nS->next = nS->prev = nS;
    WlzGMShellAppend(cS, nS);
    WlzGMLoopTSetT(nLT); /* TODO Check this fn. */
    WlzGMModelFreeE(model, dE);
    WlzGMModelFreeET(model, tET0);
    WlzGMModelFreeET(model, tET1);
    WlzGMModelFreeVT(model, tVT0);
    WlzGMModelFreeVT(model, tVT1);
    errNum = WlzGMShellComputeGBB(cS);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMShellComputeGBB(nS);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Deletes an edge which seperates two loop topology elements.
*		Because the edge is internal to the shell then the geometry
*		of the shell is not changed.
*
*                        o--------------->   o--------------->       
*                 -----@-------------------@-------------------@----
*                      | <---------------o | <---------------o |    
*                      |                 ^ | o tVT1 pET1       |    
*                      |                 | | |                 |    
*                      |                 | | |                 |    
*                      |            tET0 | | | tET1            |    
*                      |       tLT0      | | |      tLT1       |    
*                      |                 | |--- dE             |    
*                      |                 | | |                 |    
*                      |       pET0 tVT0 o | V                 |    
*                      | o---------------> | o---------------> |    
*                 -----@-------------------@-------------------@----
*                        <---------------o   <---------------o       
*                                                                      
*
* \param	model			Model with resources.
* \param	dE			Edge to delete.
*/
static WlzErrorNum WlzGMModelDeleteE2D2V2L(WlzGMModel *model, WlzGMEdge *dE)
{
  WlzGMEdgeT	*tET0,
  		*tET1,
		*pET0,
		*pET1,
		*nET0,
		*nET1;
  WlzGMVertexT	*tVT0,
  		*tVT1;
  WlzGMLoopT	*tLT0,
  		*tLT1;

  tET0 = dE->edgeT;
  tET1 = tET0->opp;
  pET0 = tET0->prev;
  pET1 = tET1->prev;
  nET0 = tET0->next;
  nET1 = tET1->next;
  tVT0 = tET0->vertexT;
  tVT1 = tET1->vertexT;
  nET0->prev = pET1;
  pET1->next = nET0;
  nET1->prev = pET0;
  pET0->next = nET1;
  WlzGMVertexTUnlink(tVT0);
  WlzGMVertexTUnlink(tVT1);
  tLT0 = pET0->parent;
  tLT1 = pET1->parent;
  if(tLT0->edgeT == tET0)
  {
    tLT0->edgeT = pET0;
  }
  WlzGMLoopTUnlink(tLT1);
  WlzGMLoopTSetT(tLT0);
  WlzGMModelFreeE(model, dE);
  WlzGMModelFreeET(model, tET0);
  WlzGMModelFreeET(model, tET1);
  WlzGMModelFreeVT(model, tVT0);
  WlzGMModelFreeVT(model, tVT1);
  WlzGMModelFreeLT(model, tLT1);
  return(WLZ_ERR_NONE);
}

/* Model access and testing. */

/*!
* \return				The model's dimension, either 2 or 3.
* \ingroup      WlzGeoModel
* \brief	Gets the models dimension from the model's type.
* \param	model			The model.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
int	 	WlzGMModelGetDimension(WlzGMModel *model, WlzErrorNum *dstErr)
{
  int		dim = 0;

  switch(model->type)
  {
    case WLZ_GMMOD_2I:
    case WLZ_GMMOD_2D:
      dim = 2;
      break;
    case WLZ_GMMOD_3I:
    case WLZ_GMMOD_3D:
      dim = 3;
      break;
  }
  if(dstErr)
  {
    *dstErr = (dim == 0)? WLZ_ERR_DOMAIN_TYPE: WLZ_ERR_NONE;
  }
  return(dim);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Checks the model type is valid.
* \param	type			Model type to check.
*/
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

/* Geometry access, update and testing */

/*!
* \return				Shell's geometry type.
* \ingroup      WlzGeoModel
* \brief	Gets the shell's geometry type from the model's type.
* \param	model			The model.
*/
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

/*!
* \return				Shell's geometry type.
* \ingroup      WlzGeoModel
* \brief	Gets the verticies geometry type from the model's type.
* \param	model			The model.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Gets a shell's geometry using the given destination
*               pointer to a double precision bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	bBox			Given destination pointer to a
*                                       double precision bounding box.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Gets the volume of the shell's geometry's bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	vol			Given destination pointer for
*                                       the volume.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Gets a shell's geometry using the given destination
*               pointer to a double precision bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	bBox			Given destination pointer to a
*                                       double precision bounding box.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Sets a shell's geometry using the given double
*               precision bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	bBox			Given double precision bounding
*                                       box.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Sets a shell's geometry using the given double
*               precision bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	bBox			Given double precision bounding
*                                       box.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Sets a shell's geometry using the pair of double
*               precision points.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	nPnt			Number of points.
* \param	pos			Array of points.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Sets a shell's geometry using the pair of double
*               precision points.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	nPnt			Number of points.
* \param	pos			Array of points.
*/
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

/*!
* \return	Woolz eror code.
* \ingroup      WlzGeoModel
* \brief	Sets the geometry of each of the shells in the given model.
* \param	model			Given geometric model.
*/
WlzErrorNum	WlzGMModelSetSG(WlzGMModel *model)
{
  int		idx,
  		vCnt;
  WlzBoxP	sGP;
  WlzVertexP	vGP;
  WlzGMShell	*fS,
  		*nS;
  WlzGMVertex	*nV;
  AlcVector	*vVec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* First pass: For each shell of the model just set the geometry to
     * enclose a single vertex of the shell. */
    nS = fS = model->child;
    switch(model->type)
    {
      case WLZ_GMMOD_2I:
        do
	{
          sGP.i2 = &(nS->geo.sg2I->bBox);
          nV = nS->child->edgeT->vertexT->diskT->vertex;
	  vGP.i2 = &(nV->geo.vg2I->vtx);
	  sGP.i2->xMin = sGP.i2->xMax = vGP.i2->vtX;
	  sGP.i2->yMin = sGP.i2->yMax = vGP.i2->vtY;
	  nS = nS->next;
	} while(nS != fS);
	break;
      case WLZ_GMMOD_2D:
        do
	{
          sGP.d2 = &(nS->geo.sg2D->bBox);
          nV = nS->child->edgeT->vertexT->diskT->vertex;
	  vGP.d2 = &(nV->geo.vg2D->vtx);
	  sGP.d2->xMin = sGP.d2->xMax = vGP.d2->vtX;
	  sGP.d2->yMin = sGP.d2->yMax = vGP.d2->vtY;
	  nS = nS->next;
	} while(nS != fS);
	break;
      case WLZ_GMMOD_3I:
        do
	{
          sGP.i3 = &(nS->geo.sg3I->bBox);
          nV = nS->child->edgeT->vertexT->diskT->vertex;
	  vGP.i3 = &(nV->geo.vg3I->vtx);
	  sGP.i3->xMin = sGP.i3->xMax = vGP.i3->vtX;
	  sGP.i3->yMin = sGP.i3->yMax = vGP.i3->vtY;
	  sGP.i3->zMin = sGP.i3->zMax = vGP.i3->vtZ;
	  nS = nS->next;
	} while(nS != fS);
	break;
      case WLZ_GMMOD_3D:
        do
	{
          sGP.d3 = &(nS->geo.sg3D->bBox);
          nV = nS->child->edgeT->vertexT->diskT->vertex;
	  vGP.d3 = &(nV->geo.vg3D->vtx);
	  sGP.d3->xMin = sGP.d3->xMax = vGP.d3->vtX;
	  sGP.d3->yMin = sGP.d3->yMax = vGP.d3->vtY;
	  sGP.d3->zMin = sGP.d3->zMax = vGP.d3->vtZ;
	  nS = nS->next;
	} while(nS != fS);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  /* Second pass: Now for each vertex of the model get the parent shell and
   * update it's geometry. */
  if(errNum == WLZ_ERR_NONE)
  {
    vCnt = model->res.vertex.numElm;
    vVec = model->res.vertex.vec;
    switch(model->type)
    {
      case WLZ_GMMOD_2I:
	for(idx = 0; idx < vCnt; ++idx)
	{
	  nV = (WlzGMVertex *)AlcVectorItemGet(vVec, idx);
	  nS = nV->diskT->vertexT->parent->parent->parent;
	  sGP.i2 = &(nS->geo.sg2I->bBox);
	  vGP.i2 = &(nV->geo.vg2I->vtx);
	  if(vGP.i2->vtX < sGP.i2->xMin)
	  {
	    sGP.i2->xMin = vGP.i2->vtX;
	  }
	  else if(vGP.i2->vtX > sGP.i2->xMax)
	  {
	    sGP.i2->xMax = vGP.i2->vtX;
	  }
	  if(vGP.i2->vtY < sGP.i2->yMin)
	  {
	    sGP.i2->yMin = vGP.i2->vtY;
	  }
	  else if(vGP.i2->vtY > sGP.i2->yMax)
	  {
	    sGP.i2->yMax = vGP.i2->vtY;
	  }
	}
	break;
      case WLZ_GMMOD_2D:
	for(idx = 0; idx < vCnt; ++idx)
	{
	  nV = (WlzGMVertex *)AlcVectorItemGet(vVec, idx);
	  nS = nV->diskT->vertexT->parent->parent->parent;
	  sGP.d2 = &(nS->geo.sg2D->bBox);
	  vGP.d2 = &(nV->geo.vg2D->vtx);
	  if(vGP.d2->vtX < sGP.d2->xMin)
	  {
	    sGP.d2->xMin = vGP.d2->vtX;
	  }
	  else if(vGP.d2->vtX > sGP.d2->xMax)
	  {
	    sGP.d2->xMax = vGP.d2->vtX;
	  }
	  if(vGP.d2->vtY < sGP.d2->yMin)
	  {
	    sGP.d2->yMin = vGP.d2->vtY;
	  }
	  else if(vGP.d2->vtY > sGP.d2->yMax)
	  {
	    sGP.d2->yMax = vGP.d2->vtY;
	  }
        }
	break;
      case WLZ_GMMOD_3I:
	for(idx = 0; idx < vCnt; ++idx)
	{
	  nV = (WlzGMVertex *)AlcVectorItemGet(vVec, idx);
	  nS = nV->diskT->vertexT->parent->parent->parent;
	  sGP.i3 = &(nS->geo.sg3I->bBox);
	  vGP.i3 = &(nV->geo.vg3I->vtx);
	  if(vGP.i3->vtX < sGP.i3->xMin)
	  {
	    sGP.i3->xMin = vGP.i3->vtX;
	  }
	  else if(vGP.i3->vtX > sGP.i3->xMax)
	  {
	    sGP.i3->xMax = vGP.i3->vtX;
	  }
	  if(vGP.i3->vtY < sGP.i3->yMin)
	  {
	    sGP.i3->yMin = vGP.i3->vtY;
	  }
	  else if(vGP.i3->vtY > sGP.i3->yMax)
	  {
	    sGP.i3->yMax = vGP.i3->vtY;
	  }
	  if(vGP.i3->vtZ < sGP.i3->zMin)
	  {
	    sGP.i3->zMin = vGP.i3->vtZ;
	  }
	  else if(vGP.i3->vtZ > sGP.i3->zMax)
	  {
	    sGP.i3->zMax = vGP.i3->vtZ;
	  }
	}
	break;
      case WLZ_GMMOD_3D:
	for(idx = 0; idx < vCnt; ++idx)
	{
	  nV = (WlzGMVertex *)AlcVectorItemGet(vVec, idx);
	  nS = nV->diskT->vertexT->parent->parent->parent;
	  sGP.d3 = &(nS->geo.sg3D->bBox);
	  vGP.d3 = &(nV->geo.vg3D->vtx);
	  if(vGP.d3->vtX < sGP.d3->xMin)
	  {
	    sGP.d3->xMin = vGP.d3->vtX;
	  }
	  else if(vGP.d3->vtX > sGP.d3->xMax)
	  {
	    sGP.d3->xMax = vGP.d3->vtX;
	  }
	  if(vGP.d3->vtY < sGP.d3->yMin)
	  {
	    sGP.d3->yMin = vGP.d3->vtY;
	  }
	  else if(vGP.d3->vtY > sGP.d3->yMax)
	  {
	    sGP.d3->yMax = vGP.d3->vtY;
	  }
	  if(vGP.d3->vtZ < sGP.d3->zMin)
	  {
	    sGP.d3->zMin = vGP.d3->vtZ;
	  }
	  else if(vGP.d3->vtZ > sGP.d3->zMax)
	  {
	    sGP.d3->zMax = vGP.d3->vtZ;
	  }
	}
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Updates a shell's geometry using the given double
*               precision position.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	pos			Given position.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Updates a shell's geometry using the new given double
*               precision position.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	pos			Given new position.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Updates a shell's geometry using the given double
*               precision position bounding box. The updated geometry
*               is the union of the existing geometry and the given
*               bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	bBox			Given bounding box.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Updates a shell's geometry using the given double
*               precision position bounding box. The updated geometry
*               is the union of the existing geometry and the given
*               bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	bBox			Given bounding box.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Updates the shell's geometry for the removed vertex
*		with the given position.
* \param	shell			Given shell with geometry.
* \param	pos			The removed vertex position.
*/
WlzErrorNum	WlzGMShellDndateG2D(WlzGMShell *shell, WlzDVertex2 pos)
{
  double	tD0;
  int		recompFlg[4];
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
	recompFlg[0] = (pos.vtX == shell->geo.sg2I->bBox.xMin) ||
		       (pos.vtX == shell->geo.sg2I->bBox.xMax) ||
		       (pos.vtY == shell->geo.sg2I->bBox.yMin) ||
		       (pos.vtY == shell->geo.sg2I->bBox.yMax);
	break;
      case WLZ_GMELM_SHELL_G2D:
	tD0 = pos.vtX - shell->geo.sg2D->bBox.xMin;
	recompFlg[0] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtX - shell->geo.sg2D->bBox.xMax;
	recompFlg[1] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtY - shell->geo.sg2D->bBox.yMin;
	recompFlg[2] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtY - shell->geo.sg2D->bBox.yMax;
	recompFlg[3] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	recompFlg[0] = recompFlg[0] || recompFlg[1] ||
		       recompFlg[2] || recompFlg[3];
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && recompFlg[0])
  {
    errNum = WlzGMShellComputeGBB(shell);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Updates the shell's geometry for the removed vertex
*		with the given position.
* \param	shell			Given shell with geometry.
* \param	pos			The removed vertex position.
*/
WlzErrorNum	WlzGMShellDndateG3D(WlzGMShell *shell, WlzDVertex3 pos)
{
  double	tD0;
  int		recompFlg[6];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((shell == NULL) || (shell->geo.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(shell->geo.core->type)
    {
      case WLZ_GMELM_SHELL_G3I:
	recompFlg[0] = (pos.vtX == shell->geo.sg3I->bBox.xMin) ||
		       (pos.vtX == shell->geo.sg3I->bBox.xMax) ||
		       (pos.vtY == shell->geo.sg3I->bBox.yMin) ||
		       (pos.vtY == shell->geo.sg3I->bBox.yMax) ||
		       (pos.vtZ == shell->geo.sg3I->bBox.zMin) ||
		       (pos.vtZ == shell->geo.sg3I->bBox.zMax);
	break;
      case WLZ_GMELM_SHELL_G3D:
	tD0 = pos.vtX - shell->geo.sg3D->bBox.xMin;
	recompFlg[0] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtX - shell->geo.sg3D->bBox.xMax;
	recompFlg[1] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtY - shell->geo.sg3D->bBox.yMin;
	recompFlg[2] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtY - shell->geo.sg3D->bBox.yMax;
	recompFlg[3] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtZ - shell->geo.sg3D->bBox.zMin;
	recompFlg[4] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	tD0 = pos.vtZ - shell->geo.sg3D->bBox.zMax;
	recompFlg[5] = (tD0 * tD0) < WLZ_GM_TOLERANCE;
	recompFlg[0] = recompFlg[0] || recompFlg[1] ||
		       recompFlg[2] || recompFlg[3] ||
		       recompFlg[4] || recompFlg[5] ;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && recompFlg[0])
  {
    errNum = WlzGMShellComputeGBB(shell);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Recomputes the shell's geometry by walking through
*               it's children.
* \param	shell			Given shell with geometry to
*                                       be set.
*/
WlzErrorNum	WlzGMShellComputeGBB(WlzGMShell *shell)
{
  int		tI0;
  double	tD0;
  WlzGMEdgeT	*fET,
  		*tET;
  WlzGMLoopT	*fLT,
  		*tLT;
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
	fLT = tLT = shell->child;
	do
	{
	  /* For each edgeT */
	  fET = tET = tLT->edgeT;
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
	  } while(tET != fET);
	  tLT = tLT->next;
	} while(tLT != fLT);
	break;
      case WLZ_GMELM_SHELL_G2D:
	vGU = shell->child->edgeT->vertexT->diskT->vertex->geo;
	sGU.sg2D->bBox.xMin = sGU.sg2D->bBox.xMax = vGU.vg2D->vtx.vtX;
	sGU.sg2D->bBox.yMin = sGU.sg2D->bBox.yMax = vGU.vg2D->vtx.vtY;
	/* For each loopT */
	fLT = tLT = shell->child;
	do
	{
	  /* For each edgeT */
	  fET = tET = tLT->edgeT;
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
	  } while(tET != fET);
	  tLT = tLT->next;
	} while(tLT != fLT);
	break;
      case WLZ_GMELM_SHELL_G3I:
	vGU = shell->child->edgeT->vertexT->diskT->vertex->geo;
	sGU.sg3I->bBox.xMin = sGU.sg3I->bBox.xMax = vGU.vg3I->vtx.vtX;
	sGU.sg3I->bBox.yMin = sGU.sg3I->bBox.yMax = vGU.vg3I->vtx.vtY;
	sGU.sg3I->bBox.zMin = sGU.sg3I->bBox.zMax = vGU.vg3I->vtx.vtZ;
	/* For each loopT */
	fLT = tLT = shell->child;
	do
	{
	  /* For each edgeT */
	  fET = tET = tLT->edgeT;
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
	  } while(tET != fET);
	  tLT = tLT->next;
	} while(tLT != fLT);
	break;
      case WLZ_GMELM_SHELL_G3D:
	vGU = shell->child->edgeT->vertexT->diskT->vertex->geo;
	sGU.sg3D->bBox.xMin = sGU.sg3D->bBox.xMax = vGU.vg3D->vtx.vtX;
	sGU.sg3D->bBox.yMin = sGU.sg3D->bBox.yMax = vGU.vg3D->vtx.vtY;
	sGU.sg3D->bBox.zMin = sGU.sg3D->bBox.zMax = vGU.vg3D->vtx.vtZ;
	/* For each loopT */
	fLT = tLT = shell->child;
	do
	{
	  /* For each edgeT */
	  fET = tET = tLT->edgeT;
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
	  } while(tET != fET);
	  tLT = tLT->next;
	} while(tLT != fLT);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Merges the geometries of the two shells so that the
*               first shell's geometry is set to the union of the
*               two shell's bounding boxes.
*		Both shells are assumed valid.
* \param	s0			First shell with geometry to
*                                       be both used and set.
* \param	s1			Second shell with geometry to
*                                       be used.
*/
static void	WlzGMShellMergeG(WlzGMShell *s0, WlzGMShell *s1)
{
  switch(s0->geo.core->type)
  {
    case WLZ_GMELM_SHELL_G2I:
      if(s1->geo.sg2I->bBox.xMin < s0->geo.sg2I->bBox.xMin)
      {
	s0->geo.sg2I->bBox.xMin = s1->geo.sg2I->bBox.xMin;
      }
      if(s1->geo.sg2I->bBox.yMin < s0->geo.sg2I->bBox.yMin)
      {
	s0->geo.sg2I->bBox.yMin = s1->geo.sg2I->bBox.yMin;
      }
      if(s1->geo.sg2I->bBox.xMax > s0->geo.sg2I->bBox.xMax)
      {
	s0->geo.sg2I->bBox.xMax = s1->geo.sg2I->bBox.xMax;
      }
      if(s1->geo.sg2I->bBox.yMax > s0->geo.sg2I->bBox.yMax)
      {
	s0->geo.sg2I->bBox.yMax = s1->geo.sg2I->bBox.yMax;
      }
      break;
    case WLZ_GMELM_SHELL_G2D:
      if(s1->geo.sg2D->bBox.xMin < s0->geo.sg2D->bBox.xMin)
      {
	s0->geo.sg2D->bBox.xMin = s1->geo.sg2D->bBox.xMin;
      }
      if(s1->geo.sg2D->bBox.yMin < s0->geo.sg2D->bBox.yMin)
      {
	s0->geo.sg2D->bBox.yMin = s1->geo.sg2D->bBox.yMin;
      }
      if(s1->geo.sg2D->bBox.xMax > s0->geo.sg2D->bBox.xMax)
      {
	s0->geo.sg2D->bBox.xMax = s1->geo.sg2D->bBox.xMax;
      }
      if(s1->geo.sg2D->bBox.yMax > s0->geo.sg2D->bBox.yMax)
      {
	s0->geo.sg2D->bBox.yMax = s1->geo.sg2D->bBox.yMax;
      }
      break;
    case WLZ_GMELM_SHELL_G3I:
      if(s1->geo.sg3I->bBox.xMin < s0->geo.sg3I->bBox.xMin)
      {
	s0->geo.sg3I->bBox.xMin = s1->geo.sg3I->bBox.xMin;
      }
      if(s1->geo.sg3I->bBox.yMin < s0->geo.sg3I->bBox.yMin)
      {
	s0->geo.sg3I->bBox.yMin = s1->geo.sg3I->bBox.yMin;
      }
      if(s1->geo.sg3I->bBox.zMin < s0->geo.sg3I->bBox.zMin)
      {
	s0->geo.sg3I->bBox.zMin = s1->geo.sg3I->bBox.zMin;
      }
      if(s1->geo.sg3I->bBox.xMax > s0->geo.sg3I->bBox.xMax)
      {
	s0->geo.sg3I->bBox.xMax = s1->geo.sg3I->bBox.xMax;
      }
      if(s1->geo.sg3I->bBox.yMax > s0->geo.sg3I->bBox.yMax)
      {
	s0->geo.sg3I->bBox.yMax = s1->geo.sg3I->bBox.yMax;
      }
      if(s1->geo.sg3I->bBox.zMax > s0->geo.sg3I->bBox.zMax)
      {
	s0->geo.sg3I->bBox.zMax = s1->geo.sg3I->bBox.zMax;
      }
      break;
    case WLZ_GMELM_SHELL_G3D:
      if(s1->geo.sg3D->bBox.xMin < s0->geo.sg3D->bBox.xMin)
      {
	s0->geo.sg3D->bBox.xMin = s1->geo.sg3D->bBox.xMin;
      }
      if(s1->geo.sg3D->bBox.yMin < s0->geo.sg3D->bBox.yMin)
      {
	s0->geo.sg3D->bBox.yMin = s1->geo.sg3D->bBox.yMin;
      }
      if(s1->geo.sg3D->bBox.zMin < s0->geo.sg3D->bBox.zMin)
      {
	s0->geo.sg3D->bBox.zMin = s1->geo.sg3D->bBox.zMin;
      }
      if(s1->geo.sg3D->bBox.xMax > s0->geo.sg3D->bBox.xMax)
      {
	s0->geo.sg3D->bBox.xMax = s1->geo.sg3D->bBox.xMax;
      }
      if(s1->geo.sg3D->bBox.yMax > s0->geo.sg3D->bBox.yMax)
      {
	s0->geo.sg3D->bBox.yMax = s1->geo.sg3D->bBox.yMax;
      }
      if(s1->geo.sg3D->bBox.zMax > s0->geo.sg3D->bBox.zMax)
      {
	s0->geo.sg3D->bBox.zMax = s1->geo.sg3D->bBox.zMax;
      }
      break;
  }
}

/*!
* \return				Non zero if point inside shell's
*                                       bounding box.
* \ingroup      WlzGeoModel
* \brief	Checks to see if the given double precision position
*               is within the shell's bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	pos			Given double precision position.
*/
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

/*!
* \return				Non zero if point inside shell's
*                                       bounding box.
* \ingroup      WlzGeoModel
* \brief	Checks to see if the given double precision position
*               is within the shell's bounding box.
* \param	shell			Given shell with geometry to
*                                       be set.
* \param	pos			Given double precision position.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Sets a veticies geometry using the given double
*               precision position.
* \param	vertex			Given vertex with geometry to
*                                       be set.
* \param	pos			Given position.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Sets a veticies geometry using the given double
*               precision position.
* \param	vertex			Given vertex with geometry to
*                                       be set.
* \param	pos			Given position.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Gets a veticies geometry into the given double
*               precision position destination pointer.
* \param	vertex			Given vertex with geometry to
*                                       be set.
* \param	dstPos			Given position destination
*                                       pointer, may NOT be NULL.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Gets a veticies geometry into the given double
*               precision position destination pointer.
* \param	vertex			Given vertex with geometry to
*                                       be set.
* \param	dstPos			Given position destination
*                                       pointer, may NOT be NULL.
*/
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

/*!
* \return				Position of vertex \f$-\f$ given
*					position.
* \ingroup      WlzGeoModel
* \brief	Compares the position of the given vertex with the
*		given 3D double precision position.
* \param	vertex			Given vertex.
* \param	pos			Given position.
*/
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

/*!
* \return                               Position of vertex \f$-\f$ given
*                                       position.
* \ingroup      WlzGeoModel
* \brief        Compares the position of the given vertex with the
*               given 2D double precision position.
* \param        vertex                  Given vertex.
* \param        pos                     Given position.
*/
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

/*!
* \return				The sign of the vertex position:
*					-1, 0 or +1.
* \ingroup      WlzGeoModel
* \brief	Compares the coordinates of the given vertex and
*		3D double precision position to find a signed value
*               for sorting.
* \param	vertex			Given vertex.
* \param	pos			Given position.
*/
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

/*!
* \return				The sign of the vertex position
*					- the given position: -1, 0 or +1.
* \ingroup      WlzGeoModel
* \brief	Compares the coordinates of the given vertex and
*		2D double precision position to find a signed value
*               for sorting.
* \param	vertex			Given vertex.
* \param	pos			Given position.
*/
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

/*!
* \return				Square of distance, -1.0 on
*                                       error.
* \ingroup      WlzGeoModel
* \brief	Calculates the square of the Euclidean distance
*               between the given vertex and the given 3D double
*               precision position.
* \param	vertex			Given vertex.
* \param	pos			Given position.
*/
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

/*!
* \return				Value of normal.
* \ingroup      WlzGeoModel
* \brief	Computes the value of the normal at the given vertex
*               which lies within the given model.
*               This function requires a buffer in which to store the
*               vertices found on the loops surrounding the given
*               vertex. For efficiency this can/should be reused
*               between calls of this function.
* \param	model			The given model.
* \param	gV			Given vertex in the model.
* \param	sVBufSz			Ptr to the number WlzGMVertex's
*                                       that can be held in *sVBuf.
* \param	sVBuf			Ptr to an allocated buffer for
*                                       vertices, may NOT be NULL
*                                       although the buffer it points
*                                       to may be. The buffer should
*                                       be free'd using AlcFree when
*                                       it is no longer needed.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
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
  if(gV->diskT != gV->diskT->next)
  {
    /* The surface in the imediate neighbourhood of the vertex is
     * not manifold.*/
    manifold = 0;
  }
  else
  {
    /* Find the ordered ring of vertices which surround the current
     * vertex, checking that each of the edges is part of a manifold
     * surface. */

    sCnt = 0;
    vT1 = vT0 = gV->diskT->vertexT;
    eT1 = vT1->parent;
    do
    {
      if((eT1 == eT1->rad) || (eT1 != eT1->rad->rad))
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
    while((errNum == WLZ_ERR_NONE) && manifold && (vT1 != vT0));
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
    *dstErr = errNum;
  }
  return(nrm);
}

/*!
* \return				Square of distance, -1.0 on
*                                       error.
* \ingroup      WlzGeoModel
* \brief	Calculates the square of the Euclidean distance
*               between the given vertex and the given 2D double
*               precision position.
* \param	vertex			Given vertex.
* \param	pos			Given position.
*/
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

/*!
* \return				Hash value.
* \ingroup      WlzGeoModel
* \brief	Computes a hash value from a given 3D double precision
*               position.
* \param	pos			Given position.
*/
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

/*!
* \return				Hash value.
* \ingroup      WlzGeoModel
* \brief	Computes a hash value from a given 2D double precision
*               position.
* \param	pos			Given position.
*/
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

/*!
* \return				Matched vertex, or NULL if no
*                                       vertex with the given geometry.
* \ingroup      WlzGeoModel
* \brief	Attempts to find a vertex which matches the given
*               double precision 3D position.
* \param	model			Model with resources.
* \param	gPos			Position to match.
*/
WlzGMVertex	*WlzGMModelMatchVertexG3D(WlzGMModel *model, WlzDVertex3 gPos)
{
  int		cmp;
  unsigned int	hVal;
  WlzGMVertex	*tV,
  		*mV = NULL;

  hVal = WlzGMHashPos3D(gPos);
  if((tV = *(model->vertexHT + (hVal % model->vertexHTSz))) != NULL)
  {
    /* Have found the hash table head, a linked list of vertices with this
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

/*!
* \return				Matched vertex, or NULL if no
*                                       vertex with the given geometry.
* \ingroup      WlzGeoModel
* \brief	Attempts to find a vertex which matches the given
*               double precision 2D position.
* \param	model			Model with resources.
* \param	gPos			Position to match.
*/
WlzGMVertex	*WlzGMModelMatchVertexG2D(WlzGMModel *model, WlzDVertex2 gPos)
{
  int		cmp;
  unsigned int	hVal;
  WlzGMVertex	*tV,
  		*mV = NULL;

  hVal = WlzGMHashPos2D(gPos);
  if((tV = *(model->vertexHT + (hVal % model->vertexHTSz))) != NULL)
  {
    /* Have found the hash table head, a linked list of vertices with this
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

/*!
* \return	Minimum constrained distance between the vertices or a
*		negative value if the two vertices are in different shells.
* \ingroup      WlzGeoModel
* \brief	Computes the minimum distance between the two vertices
*		where the path between them is constrained to the common
*		shell, this will be negative if the two vertices are
*		in different shells.
* \param	v0			First vertex.
* \param	v1			Second vertex.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
double		WlzGMVertexShellDist(WlzGMVertex *v0, WlzGMVertex *v1,
				     WlzErrorNum *dstErr)
{
  int 		dim;
  double	dist = -1.0;
  WlzGMShell	*s0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((v0 == NULL) || (v1 == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((s0 = WlzGMVertexCommonShell(v0, v1)) != NULL)
  {
    dim = WlzGMModelGetDimension(s0->parent, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(dim)
      {
        case 2:
	  dist = WlzGMVertexShellDist2(v0, v1, &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_UNIMPLEMENTED;
	  break;
      }
    }
  }
  if(dstErr)
  {
    errNum = *dstErr;
  }
  return(dist);
}

/*!
* \return	Minimum constrained distance between the vertices.
* \ingroup      WlzGeoModel
* \brief	Computes the minimum distance between the two vertices
*		where the path between them is constrained to the edges
*		of the common shell and the shell is the child of a 2D
*		model.
* \param	v0			First vertex.
* \param	v1			Second vertex.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static double	WlzGMVertexShellDist2(WlzGMVertex *v0, WlzGMVertex *v1,
				      WlzErrorNum *dstErr)
{
  double	dist = DBL_MAX;
  WlzGMVertexT	*vT0,
  		*vT1,
		*fVT0,
  		*fVT1;
  WlzGMLoopT	*lT0,
  		*lT1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check for a common loop topology element. */
  vT0 = fVT0 = v0->diskT->vertexT;
  vT1 = fVT1 = v1->diskT->vertexT;
  do
  {
    lT0 = vT0->parent->parent;
    do
    {
      lT1 = vT1->parent->parent;
      if(lT0 != lT1)
      {
        vT1 = vT1->next;
      }
    } while((vT1 != fVT1) && (lT0 != lT1));
    if(lT0 != lT1)
    {
      vT0 = vT0->next;
    }
  } while((vT0 != fVT0) && (lT0 != lT1));
  if(lT0 == lT1)
  {
    dist = WlzGMVertexShellDist2CommonLT(vT0, vT1);
  }
  else
  {
    /* Vertices are in different loops.
     * This code hasn't been implemented yet!
     */
    errNum = WLZ_ERR_UNIMPLEMENTED;
  }
  if(dstErr)
  {
    errNum = *dstErr;
  }
  return(dist);
}

/*!
* \return	Minimum constrained distance between the vertices of
*		the vertex topology elements.
* \ingroup      WlzGeoModel
* \brief	Computes the minimum distance between the vertices of
*		the two vertex topology elements which share a common
*		loop topology element a 2D model.
* \param	gVT0			First vertex topology element.
* \param	gVT1			Second vertex topology element.
*/
static double	WlzGMVertexShellDist2CommonLT(WlzGMVertexT *gVT0,
					      WlzGMVertexT *gVT1)
{
  double	dist = 0.0,
  		distN = 0.0,
  		distP = 0.0;
  WlzDVertex2	posN[2],
  		posP[2],
		diff;
  WlzGMEdgeT	*eTN,
  		*eTP;
  WlzGMVertexT	*vTN,
  		*vTP;

  if(gVT0->diskT->vertex != gVT1->diskT->vertex)
  {
    /* Move around the loopT both ways, using the next and previous vertexT's,
     * computing distance between the two verticies */
    eTN = eTP = gVT0->parent;
    (void )WlzGMVertexGetG2D(eTN->vertexT->diskT->vertex, posN + 1);
    posP[1] = posN[1];
    do
    {
      posN[0] = posN[1];
      posP[0] = posP[1];
      eTN = eTN->next;
      eTP = eTP->prev;
      vTN = eTN->vertexT;
      vTP = eTP->vertexT;
      (void )WlzGMVertexGetG2D(vTN->diskT->vertex, posN + 1);
      (void )WlzGMVertexGetG2D(vTP->diskT->vertex, posP + 1);
      WLZ_VTX_2_SUB(diff, posN[0], posN[1]);
      distN += WLZ_VTX_2_LENGTH(diff);
      WLZ_VTX_2_SUB(diff, posP[0], posP[1]);
      distP += WLZ_VTX_2_LENGTH(diff);
    } while((vTN != gVT1) && (vTP != gVT1));
    dist = (vTN == gVT1)? distN: distP;
  }
  return(dist);
}

/* Topology update, access and checking. */
/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Verifies the given model. If the invalid element destination
*		pointer is given then it MAY be set to point to an invalid
*		element of the model.
*		This function was written as an aid to debugging.
* \param	model			Given model.
* \param	dstElmP			Destination pointer for an
*					invalid element pointermay be NULL.
*/
WlzErrorNum	WlzGMVerifyModel(WlzGMModel *model, WlzGMElemP *dstElmP)
{
  int		maxS;
  WlzGMShell	*cS,
  		*fS;
  WlzGMElemP	elmP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmP.core = NULL;
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
    if((fS = model->child) != NULL)
    {
      if(fS->type != WLZ_GMELM_SHELL)
      {
	elmP.shell = fS;
        errNum = WLZ_ERR_GMELM_TYPE;
      }
      else
      {
	cS = fS;
	maxS = model->res.shell.numElm;
	do
	{
	  if(cS == NULL)
	  {
	    errNum = WLZ_ERR_GMELM_NULL;
	  }
	  else
	  {
	    errNum = WlzGMVerifyShell(cS, &elmP);
	    cS = cS->next;
	  }
	} while((errNum == WLZ_ERR_NONE) && (cS != fS) && (maxS-- > 0));
	if(errNum == WLZ_ERR_NONE)
	{
	  if(cS != fS)
	  {
	    elmP.shell = cS;
	    errNum = WLZ_ERR_GMELM_DATA;
	  }
	}
      }
    }
  }
  if((errNum != WLZ_ERR_NONE) && dstElmP)
  {
    dstElmP->core = elmP.core;
  }
  return(errNum);
}

/*!
* \return
* \ingroup      WlzGeoModel
* \brief	Verifies the given shell. If the invalid element destination
*               pointer is given then it MAY be set to point to an invalid
*               element of the model.
*               This function was written as an aid to debugging.
* \param	shell			Given shell.
* \param	dstElmP			Destination pointer for an
*					invalid element pointer, may be NULL.
*/
WlzErrorNum	WlzGMVerifyShell(WlzGMShell *shell, WlzGMElemP *dstElmP)
{
  int		maxLT;
  WlzGMLoopT	*cLT,
  		*fLT;
  WlzGMElemP	elmP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmP.core = NULL;
  if(shell == NULL)
  {
    errNum = WLZ_ERR_GMELM_NULL;
  }
  else if(shell->type != WLZ_GMELM_SHELL)
  {
    elmP.shell = shell;
    errNum = WLZ_ERR_GMELM_TYPE;
  }
  else if(shell->parent == NULL)
  {
    elmP.shell = shell;
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(shell->idx < 0)
  {
    elmP.shell = shell;
    errNum = WLZ_ERR_GMELM_DATA;
  }
  else if ((shell->next == NULL) || (shell->prev == NULL))
  {
    elmP.shell = shell;
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((shell->next->prev != shell) || (shell->prev->next != shell))
  {
    errNum = WLZ_ERR_GMELM_NULL;
  }
  else if((fLT = shell->child) == NULL)
  {
    errNum = WLZ_ERR_GMELM_NULL;
  }
  else
  {
    cLT = fLT;
    maxLT = shell->parent->res.loopT.numElm;
    do
    {
      if(cLT == NULL)
      {
        errNum = WLZ_ERR_GMELM_NULL;
      }
      else
      {
        errNum = WlzGMVerifyLoopT(cLT, &elmP);
	cLT = cLT->next;
      }
    } while((errNum == WLZ_ERR_NONE) && (cLT != fLT) && (maxLT-- > 0));
    if(errNum == WLZ_ERR_NONE)
    {
      if(cLT != fLT)
      {
        elmP.loopT = cLT;
	errNum = WLZ_ERR_GMELM_DATA;
      }
    }
  }
  if((errNum != WLZ_ERR_NONE) && dstElmP)
  {
    dstElmP->core = elmP.core;
  }
  return(errNum);
}

/*!
* \return
* \ingroup      WlzGeoModel
* \brief	Verifies the given loop topology element. If the invalid
*		element destination pointer is given then it MAY be set
*		to point to an invalid element of the model.
*               This function was written as an aid to debugging.
* \param	loopT			Given loop topology element.
* \param	dstElmP			Destination pointer for an
*					invalid element pointer, may be NULL.
*/
WlzErrorNum	WlzGMVerifyLoopT(WlzGMLoopT *loopT, WlzGMElemP *dstElmP)
{
  int		maxET;
  WlzGMDiskT	*dT;
  WlzGMEdgeT	*cET,
  		*fET;
  WlzGMElemP	elmP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmP.core = NULL;
  if((loopT == NULL) || (loopT->parent == NULL))
  {
    errNum = WLZ_ERR_GMELM_NULL;
  }
  else if(loopT->parent->parent == NULL)
  {
    elmP.shell = loopT->parent;
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(loopT->type != WLZ_GMELM_LOOP_T)
  {
    elmP.loopT = loopT;
    errNum = WLZ_ERR_GMELM_TYPE;
  }
  else if(loopT->parent->type != WLZ_GMELM_SHELL)
  {
    elmP.shell = loopT->parent;
    errNum = WLZ_ERR_GMELM_TYPE;
  }
  else if(loopT->idx < 0)
  {
    elmP.loopT = loopT;
    errNum = WLZ_ERR_GMELM_DATA;
  }
  else if((loopT->next == NULL) || (loopT->prev == NULL))
  {
    elmP.loopT = loopT;
    errNum = WLZ_ERR_GMELM_NULL;
  }
  else if((loopT->next->prev != loopT) || (loopT->prev->next != loopT))
  {
    elmP.loopT = loopT;
    errNum = WLZ_ERR_GMELM_DATA;
  }
  else if((fET = loopT->edgeT) == NULL)
  {
    errNum = WLZ_ERR_GMELM_NULL;
  }
  else
  {
    cET = fET;
    maxET = loopT->parent->parent->res.edgeT.numElm;
    do
    {
      if((cET == NULL) ||
	 (cET->next == NULL) || (cET->prev == NULL) ||
	 (cET->opp == NULL) || (cET->rad == NULL) ||
	 (cET->edge == NULL) || (cET->vertexT == NULL) ||
	 (cET->vertexT->diskT == NULL))
      {
	elmP.edgeT = cET;
	errNum = WLZ_ERR_GMELM_NULL;
      }
      else if((cET->idx < 0) ||
              (cET->next->idx < 0) || (cET->prev->idx < 0) || 
              (cET->opp->idx < 0) || (cET->rad->idx < 0) ||
              (cET->edge->idx < 0) || (cET->vertexT->idx < 0) ||
	      (cET->vertexT->diskT->idx < 0))
      {
        elmP.edgeT = cET;
	errNum = WLZ_ERR_GMELM_DATA;
      }
      else if((cET->type != WLZ_GMELM_EDGE_T) ||
              (cET->next->type != WLZ_GMELM_EDGE_T) ||
	      (cET->prev->type != WLZ_GMELM_EDGE_T) || 
              (cET->opp->type != WLZ_GMELM_EDGE_T) ||
	      (cET->rad->type != WLZ_GMELM_EDGE_T) ||
              (cET->edge->type != WLZ_GMELM_EDGE) ||
	      (cET->vertexT->type != WLZ_GMELM_VERTEX_T) ||
	      (cET->vertexT->diskT->type != WLZ_GMELM_DISK_T))
      {
        elmP.edgeT = cET;
	errNum = WLZ_ERR_GMELM_TYPE;
      }
      else
      {
	dT = cET->vertexT->diskT;
	if((dT->next == NULL) || (dT->prev == NULL) ||
	   (dT->vertex == NULL) || (dT->vertexT == NULL))
	{
	  elmP.diskT = dT;
	  errNum = WLZ_ERR_GMELM_NULL;
	}
	else if((dT->next->idx < 0) || (dT->prev->idx < 0) ||
	        (dT->vertex->idx < 0) || (dT->vertexT->idx < 0))
	{
	  elmP.diskT = dT;
	  errNum = WLZ_ERR_GMELM_DATA;
	}
	else if((dT->next->type != WLZ_GMELM_DISK_T) ||
		(dT->prev->type != WLZ_GMELM_DISK_T) ||
	        (dT->vertex->type != WLZ_GMELM_VERTEX) ||
	        (dT->vertexT->type != WLZ_GMELM_VERTEX_T))
	{
	  elmP.diskT = dT;
	  errNum = WLZ_ERR_GMELM_TYPE;
	}
	else if((dT->vertexT->diskT != dT) ||
	        (dT->next->prev != dT) || (dT->prev->next != dT))
	{
	  elmP.diskT = dT;
	  errNum = WLZ_ERR_GMELM_DATA;
	}
	cET = cET->next;
      }
    } while((errNum == WLZ_ERR_NONE) && (cET != fET) && (maxET-- > 0));
    if(errNum == WLZ_ERR_NONE)
    {
      if(cET != fET)
      {
        elmP.edgeT = cET;
	errNum = WLZ_ERR_GMELM_DATA;
      }
    }
  }
  if((errNum != WLZ_ERR_NONE) && dstElmP)
  {
    dstElmP->core = elmP.core;
  }
  return(errNum);
}

/*!
* \return				Ptr to an array of non-manifold
*                                       edge ptrs, NULL if none exist or
*                                       on error.
* \ingroup      WlzGeoModel
* \brief	Finds a loop topology element in common for the two
*               edge topology elements.
* \param	model			The given model.
* \param	dstNMCnt		Destination pointer for number
*                                       of non-manifold edges found.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
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
	    if((tET->rad == tET) ||
	       (tET->rad->rad != tET))
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

/*!
* \return				Common loop topology element,
*                                       NULL if it doesn't exist.
* \ingroup      WlzGeoModel
* \brief	Finds a loop topology element in common for the two
*               edge topology elements.
* \param	eT0			First edge topology element.
* \param	eT1			Second edge topology element.
*/
WlzGMLoopT	*WlzGMEdgeTCommonLoopT(WlzGMEdgeT *eT0, WlzGMEdgeT *eT1)
{
  WlzGMLoopT	*loopT = NULL;

  if(eT0->parent == eT1->parent)
  {
    loopT = eT0->parent;
  }
  return(loopT);
}

/*!
* \return				Common edge, NULL if it doesn't
                                        exist.
* \ingroup      WlzGeoModel
* \brief	Finds the edge common to the two given verticies.
* \param	eV0			First vertex element.
* \param	eV1			Second vertex element.
*/
WlzGMEdge	*WlzGMVertexCommonEdge(WlzGMVertex *eV0, WlzGMVertex *eV1)
{
  WlzGMEdge	*e0,
  		*e1,
		*eC = NULL;
  WlzGMVertexT	*vT0,
  		*vT1;
  WlzGMDiskT	*dT0,
  		*dT1;

  dT1 = eV1->diskT;
  do
  {
    vT1 = dT1->vertexT;
    do
    {
      e1 = vT1->parent->edge;
      dT0 = eV0->diskT;
      do
      {
	vT0 = dT0->vertexT;
	do
	{
	  e0 = vT0->parent->edge;
	  if(e0 == e1)
	  {
	      eC = e0;
	      goto MATCH_FOUND;
	  }
	  vT0 = vT0->next;
	} while(vT0 != dT0->vertexT);
	dT0 = dT0->next;
      } while(dT0 != eV0->diskT);
      vT1 = vT1->next;
    } while(vT1 != dT1->vertexT);
    dT1 = dT1->next;
  } while(dT1 != eV1->diskT);
  MATCH_FOUND:
  return(eC);
}

/*!
* \return				Common shell, NULL if it doesn't
                                        exist.
* \ingroup      WlzGeoModel
* \brief	Finds the shell common to the two given verticies.
* \param	eV0			First vertex element.
* \param	eV1			Second vertex element.
*/
WlzGMShell	*WlzGMVertexCommonShell(WlzGMVertex *eV0, WlzGMVertex *eV1)
{
  WlzGMShell	*s0,
  		*s1,
		*shell = NULL;

  if(((s0 = WlzGMVertexGetShell(eV0)) != NULL) &&
     ((s1 = WlzGMVertexGetShell(eV1)) != NULL) &&
     (s0 == s1))
  {
    shell = s0;
  }
  return(shell);
}

/*!
* \return				The parent shell.
* \ingroup      WlzGeoModel
* \brief	Finds the parent shell of the given vertex.
* \param	eV			Given vertex element.
*/
WlzGMShell	*WlzGMVertexGetShell(WlzGMVertex *eV)
{
  WlzGMShell	*pS = NULL;

  if(eV && (eV->type == WLZ_GMELM_VERTEX))
  {
    pS = eV->diskT->vertexT->parent->parent->parent;
  }
  return(pS);
}

/*!
* \return				The parent shell.
* \ingroup      WlzGeoModel
* \brief	Finds the parent shell of the given edge.
* \param	eE			Given edge element.
*/
WlzGMShell	*WlzGMEdgeGetShell(WlzGMEdge *eE)
{
  WlzGMShell	*pS = NULL;

  if(eE && (eE->type == WLZ_GMELM_EDGE))
  {
    pS = eE->edgeT->parent->parent;
  }
  return(pS);
}

/*!
* \return	The common vertex, NULL if it doesn't exist.
* \ingroup      WlzGeoModel
* \brief	Finds the common vertex of the two given edges.
* \param	eE0			First edge element.
* \param	eE1			Second edge element.
*/
WlzGMVertex	*WlzGMEdgeCommonVertex(WlzGMEdge *eE0, WlzGMEdge *eE1)
{
  WlzGMVertex	*tV,
  		*cV = NULL;

  if(((tV = eE0->edgeT->vertexT->diskT->vertex) ==
      eE1->edgeT->vertexT->diskT->vertex) ||
     (tV == eE1->edgeT->opp->vertexT->diskT->vertex))
  {
    cV = tV;
  }
  else if(((tV = eE0->edgeT->opp->vertexT->diskT->vertex) ==
           eE1->edgeT->vertexT->diskT->vertex) ||
          (tV == eE1->edgeT->opp->vertexT->diskT->vertex))
  {
    cV = tV;
  }
  return(cV);
}

/*!
* \return	The common vertex, NULL if it doesn't exist.
* \ingroup      WlzGeoModel
* \brief	Finds the common vertex of the two given edges and sets
*		the disk topology element destination pointers if there
*		is a common vertex.
*		NULL*
* \param	eE0			First edge element.
* \param	eE1			Second edge element.
* \param	dstDT0			First diskT destination pointer,
					must not be NULL.
* \param	dstDT1			Second diskT destination pointer,
					must not be NULL.
*/
WlzGMVertex	*WlzGMEdgeCommonVertexGetDiskTs(WlzGMEdge *eE0, WlzGMEdge *eE1,
				      WlzGMDiskT **dstDT0, WlzGMDiskT **dstDT1)
{
  WlzGMDiskT	*tDT00,
  		*tDT01,
		*tDT10,
		*tDT11;
  WlzGMVertex	*cV = NULL;

  tDT00 = eE0->edgeT->vertexT->diskT;
  tDT10 = eE1->edgeT->vertexT->diskT;
  if(tDT00->vertex == tDT10->vertex)
  {
    cV = tDT00->vertex;
    *dstDT0 = tDT00;
    *dstDT1 = tDT10;
  }
  else
  {
    tDT11 = eE1->edgeT->opp->vertexT->diskT;
    if(tDT00->vertex == tDT11->vertex)
    {
      cV = tDT00->vertex;
      *dstDT0 = tDT00;
      *dstDT1 = tDT11;
    }
    else
    {
      tDT01 = eE0->edgeT->opp->vertexT->diskT;
      if(tDT01->vertex == tDT10->vertex)
      {
	cV = tDT01->vertex;
	*dstDT0 = tDT01;
	*dstDT1 = tDT10;
      }
      else if(tDT01->vertex == tDT11->vertex)
      {
	cV = tDT01->vertex;
	*dstDT0 = tDT01;
	*dstDT1 = tDT11;
      }
    }
  }
  return(cV);
}

/*!
* \return				The common disk topology element,
*					NULL if it doesn't exist.
* \ingroup      WlzGeoModel
* \brief	Finds the common disk topology element of the two given
*		edges.
* \param	eE0			First edge element.
* \param	eE1			Second edge element.
*/
WlzGMDiskT	*WlzGMEdgeCommonDiskT(WlzGMEdge *eE0, WlzGMEdge *eE1)
{
  WlzGMDiskT	*tDT,
  		*cDT = NULL;

  if(((tDT = eE0->edgeT->vertexT->diskT) ==
      eE1->edgeT->vertexT->diskT) ||
     (tDT == eE1->edgeT->opp->vertexT->diskT))
  {
    cDT = tDT;
  }
  else if(((tDT = eE0->edgeT->opp->vertexT->diskT) ==
           eE1->edgeT->vertexT->diskT) ||
          (tDT == eE1->edgeT->opp->vertexT->diskT))
  {
    cDT = tDT;
  }
  return(cDT);
}

/*!
* \return	The common face or NULL if it doesn't exist.
* \brief	Finds the common face which is shared by the two given
*		edges.
* \param	eE0			First edge element.
* \param	eE1			Second edge element.
*/
WlzGMFace	*WlzGMEdgeCommonFace(WlzGMEdge *eE0, WlzGMEdge *eE1)
{
  WlzGMEdgeT	*cET0,
  		*cET1;
  WlzGMFace	*cF = NULL;

  cET0 = eE0->edgeT;
  do
  {
    cET1 = eE1->edgeT;
    do
    {
      if(cET0->parent->face == cET1->parent->face)
      {
        cF = cET0->parent->face;
	goto MATCH_FOUND;
      }
      cET1 = cET1->rad;
    } while(cET1 != eE1->edgeT);
    cET0 = cET0->rad;
  } while(cET0 != eE0->edgeT);
  MATCH_FOUND:
  return(cF);
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Merges the topologies of the two shells so that the
*               first shell 'absorbs' the second shell.
*		Both shells are assumed valid.
* \param	s0			First shell which will absorb the
*					loop topology elements of the second
*					shell.
* \param	s1			Second shell which will give it's
*					loop topology elements to the first.
*/
static void	WlzGMShellMergeT(WlzGMShell *s0, WlzGMShell *s1)
{
  WlzGMLoopT	*lt0,
  		*lt1;

  lt0 = s0->child;
  while((lt1 = s1->child) != NULL)
  {
    WlzGMLoopTUnlink(lt1);
    WlzGMLoopTAppend(lt0, lt1);
    lt1->parent = s0;
  }
}

/* Model list management */

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Append new vertex topology element onto a doubly
*               linked list of vertex topology element's, knowing that
*               neither is NULL.
*               Vertex topology elements are maintained in an unordered
*               doubly linked list.
* \param	eVT			Existing vertexT in list of
*                                       vertexT's.
* \param	nVT			New vertex topology element to
*                                       append to list after existing
*                                       element.
*/
void	   	WlzGMVertexTAppend(WlzGMVertexT *eVT, WlzGMVertexT *nVT)
{
  nVT->next = eVT->next;
  nVT->prev = eVT;
  nVT->next->prev = nVT;
  nVT->prev->next = nVT;
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Append new edge topology onto a doubly linked list
*               of disk topology element's, knowing that neither is NULL.
*               Disk topology elements are maintained in an unordered
*               doubly linked list.
* \param	eDT			Existing disk topology element in
*                                       list.
* \param	nDT			New disk topology element to
*                                       append to list after existing
*                                       element.
*/
void	   	WlzGMDiskTAppend(WlzGMDiskT *eDT, WlzGMDiskT *nDT)
{
  nDT->next = eDT->next;
  nDT->prev = eDT;
  nDT->next->prev = nDT;
  nDT->prev->next = nDT;
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Append new edge topology onto a doubly linked list
*               of edge topology element's, knowing that neither is NULL.
*               Edge topology elements are maintained in an ordered
*               doubly linked list, CCW around the inside of loops and
*               CW around the outside of 2D loops.
* \param	eET			Existing edge topology element in
*                                       list.
* \param	nET			New edge topology element to
*                                       append to list after existing
*                                       element.
*/
void	   	WlzGMEdgeTAppend(WlzGMEdgeT *eET, WlzGMEdgeT *nET)
{
  nET->next = eET->next;
  nET->prev = eET;
  nET->next->prev = nET;
  nET->prev->next = nET;
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Insert new edge topology into a doubly linked list
*               of edge topology element's, knowing that neither is NULL.
*               Edge topology elements are maintained in an ordered
*               doubly linked list, CCW around the inside of loops and
*               CW around the outside of 2D loops.
* \param	eET			Existing edge topology element in
*                                       list.
* \param	nET			New edge topology element to
*                                       insert in list before existing
*                                       element.
*/
void	   	WlzGMEdgeTInsert(WlzGMEdgeT *eET, WlzGMEdgeT *nET)
{
  nET->next = eET;
  nET->prev = eET->prev;
  nET->next->prev = nET;
  nET->prev->next = nET;
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Inserts the given new edge topology element into a
*		radially sorted cyclic list of edge topology element's.
*
*       	Inserts the given new edge topology element into a
*		radially sorted cyclic list of edge topology element's.
*		In 2D the radial edgeT is always the edgeT itself and
*		this function should not be called.
*		In 3D the radial edge is a similarly directed edge
*		topology element to the given edge topology element
*		who's off edge vertex is the next CCW off edge vertex
*		when viewed along the edge.
* \verbatim
*                                      O <-                       
*                                     /    \                      
*                                    /      \ Next radial edge    
*                                   /        |                    
*                                  /         |                    
*                                 X----------O                    
*                           Given edge into    Off edge vertex
*                           screen/page        of new simplex
* \endverbatim
*                                                                 
*		Given a directed edge which forms part of a loop
*		topology element in a new triangle simplex. Define
*		points around the the simplex:
* \verbatim
*                                O p1                                     
*                              ^ |\                                      
*                              | | \                                     
*                              | |  \                                    
*                              | |   \                                   
*                              | |    \                                  
*                          nET | |     O p2                               
*                              | |    /                                  
*                              | |   /                                   
*                              | |  /                                    
*                              | | /                                     
*                              | |/                                      
*                                O p0                                    
* \endverbatim
* 		To find the next CCW simplex.
*		Find the normal vector:
*		  \f[ n_0 = \frac{(p_0 - p_1)}{|p_0 - p_1|} \f]
*		Find the normal vector perpendicular to the new
*		(triangular) simplex:
*		  \f[ n_1 = \frac{n_0(p_2 - p_0)}{|n_0(p_2 - p_0)}\f]
*		Find the normalized vector perpendicular to \f$n_0\f$
*		and \f$n_1\f$:
*		  \f[n_2 = n_0 n_1 \f]
*		Because \f$n_1\f$ and \f$n_2\f$ are perpendicular axes
*		in the plane perpendicular to \f$p_1 - p_0\f$ we can
*		compute the angle of any vector projected onto the
*		plane as:
* \f[ \theta = \tan^{-1}{\frac{v_i n_2}{v_i n_1}}, 0 \leq \theta \leq 2\pi \f]
*		where
*		  \f[v_i = p_{2i} - p_0\f]
*		and \f$p_{2i}\f$ are the \f$p_2\f$ verticies on the
*		candidate face.
*		Use the projections to find the previous radial edge,
*		then insert the edge in the cyclic radial list.
* \param	nET 			New edge topology element to
*					insert in list. The new edge 
*					topology element MUST have all
*					it's fields set EXCEPT for the
*					radial edge link.
*/
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
  if(fET->vertexT->diskT->vertex != nET->vertexT->diskT->vertex)
  {
    /* fET is antiparrallel to the new edgeT so use it's opposite. */
    fET = fET->opp;
  }
  if(fET->rad == fET)
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
    while((tET = tET->rad) != fET)
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

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Append new loop topology onto a doubly linked list
*               of loop topology element's, knowing that neither is NULL.
*               Loop topology elements are maintained in an unordered
*               doubly linked list.
* \param	eLT			Existing loop topology element in
*                                       list.
* \param	nLT			New loop topology element to
*                                       append to list after existing
*                                       element.
*/
void	   	WlzGMLoopTAppend(WlzGMLoopT *eLT, WlzGMLoopT *nLT)
{
  nLT->next = eLT->next;
  nLT->prev = eLT;
  nLT->next->prev = nLT;
  nLT->prev->next = nLT;
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Unlinks the given loop topology element from a doubly
*               linked list of loop topology element's, knowing that
*               it is not NULL. If this is the only loop topology
*		element in parent the parent's list of loop topology
*		elements is set to NULL. Other elements which point
*		to this loop topology element are not modified.
* \param	dLT			Loop topology element to be
*                                       unlinked.
*/
void		WlzGMLoopTUnlink(WlzGMLoopT *dLT)
{
  WlzGMLoopT	*pLT,
  		*nLT;

  pLT = dLT->prev;
  nLT = dLT->next;
  pLT->next = dLT->next;
  nLT->prev = dLT->prev;
  if(dLT->parent->child == dLT)
  {
    dLT->parent->child = (dLT == nLT)? NULL: nLT;
  }
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Unlinks the given edge topology element from it's loop
*		topology element and it's opposite edge topology element.
*		The opposite edge topology element MUST be unlinked before
*		the model is valid again.
* \param	dET			Edge topology element to be
*                                       unlinked.
*/
void		WlzGMEdgeTUnlink(WlzGMEdgeT *dET)
{
    WlzGMEdgeT	*cET,
    		*pET;

  if(dET != dET->opp)
  {
    dET->prev->next = dET->opp->next;  
    dET->next->prev = dET->opp->prev;
    dET->opp->opp = dET->opp;
  }
  else
  {
    /* No longer have an opposite edgeT as it's already been unlinked. */
    /* TODO: Fix this code, it doesn't always work! */
    pET = dET->prev;
    while((cET = pET->opp->prev) != dET->prev)
    {
      pET = cET;
    }
    dET->prev->next = pET->opp;
    pET = dET->next;
    while((cET = pET->opp->next) != dET->next)
    {
      pET = cET;
    }
    dET->next->prev = pET->opp;
  }
  if(dET == dET->next)
  {
    dET->next = NULL;
  }
  if(dET == dET->parent->edgeT)
  {
    dET->parent->edgeT = (dET->next->idx >= 0)? dET->next: dET->prev;
  }
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Unlinks the given vertex topology element from a doubly
*               linked list of vertex topology element's, knowing that
*               it is not NULL. Since the parent of a vertex topology
*		element is always an edge topology element and an
*		edge topology element only ever uses a single vertex
*		topology element, set the parent's vertex topology
*		element pointer to NULL.
* \param	dVT			Vertex topology element to be
*                                       unlinked.
*/
void		WlzGMVertexTUnlink(WlzGMVertexT *dVT)
{
  WlzGMVertexT	*pVT,
  		*nVT;

  pVT = dVT->prev;
  nVT = dVT->next;
  pVT->next = dVT->next;  
  nVT->prev = dVT->prev;
  dVT->parent->vertexT = NULL;
  if(dVT->diskT->vertexT == dVT)
  {
    dVT->diskT->vertexT = (dVT == nVT)? NULL: nVT;
  }
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Unlinks the given disk topology element from a doubly
*               linked list of disk topology element's, knowing that
*               it is not NULL. If this is the only disk topology
*		element of the vertex the vertex's list of disk topology
*		elements is set to NULL. Other elements which point
*		to this disk topology element are not modified.
* \param	dDT			Disk topology element to be
*                                       unlinked.
*/
void		WlzGMDiskTUnlink(WlzGMDiskT *dDT)
{
  WlzGMDiskT	*pDT,
  		*nDT;

  pDT = dDT->prev;
  nDT = dDT->next;
  pDT->next = dDT->next;
  nDT->prev = dDT->prev;
  if(dDT->vertex->diskT == dDT)
  {
    dDT->vertex->diskT = (dDT == nDT)? NULL: nDT;
  }
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Joins the given pair of disk topology elements
*		on return the first of the diskTs has all the
*		vertexTs of both and the second of the diskTs
*		should be free'd without being unlinked.
* \param	gDT0			First disk topology element.
* \param	gDT1			Second disk topology element
*					which should be free'd on return.
*/
void		WlzGMDiskTJoin(WlzGMDiskT *gDT0, WlzGMDiskT *gDT1)
{
  WlzGMVertexT	*nVT0,
  		*pVT1;

  /* Set the diskT pointers of the second diskT's vertexT's. */
  nVT0 = gDT1->vertexT;
  do
  {
    nVT0->diskT = gDT0;
    nVT0 = nVT0->next;
  }
  while(nVT0 != gDT1->vertexT);
  /* Unlink the second diskT. */
  WlzGMDiskTUnlink(gDT1);
  /* Joint the lists of vertexTs. */
  nVT0 = gDT0->vertexT->next;
  pVT1 = gDT1->vertexT->prev;
  gDT0->vertexT->next = gDT1->vertexT;
  gDT1->vertexT->prev = gDT0->vertexT;
  nVT0->prev = pVT1;
  pVT1->next = nVT0;
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Append new shell onto a doubly linked list of
*               shells, knowing that neither is NULL.
*               A model's shells are maintained in an unordered
*               doubly linked list.
* \param	eS			Existing shell in list of shells.
* \param	nS			New shell to append to list
*                                       after existing shell.
*/
void	   	WlzGMShellAppend(WlzGMShell *eS, WlzGMShell *nS)
{
  nS->next = eS->next;
  nS->prev = eS;
  nS->next->prev = nS;
  nS->prev->next = nS;
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Unlinks the given shell from a doubly linked list of
*               shells kept by the model, knowing that the given shell
*               ptr is not NULL.
* \param	dS			Shell to be unlinked.
*/
void		WlzGMShellUnlink(WlzGMShell *dS)
{
  if(dS == dS->parent->child)
  {
    /* Change child held by parent model */
    if(dS == dS->next)
    {
      /* The model only has one shell and this is it. */
      dS->parent->child = NULL;
    }
    else
    {
      dS->parent->child = dS->next;
    }
  }
  /* Unlink the shell from the doubly linked list of shells. */
  dS->prev->next = dS->next;
  dS->next->prev = dS->prev;
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Joins the shell to be unlinked onto the shell that's to
*               be extended and then unlinks (but does not free) it.
*		Both shells must be valid.
* \param	eShell			Shell to be extended.
* \param	dShell			Shell to be joined and then
*                                       unlinked.
*/
void		WlzGMShellJoinAndUnlink(WlzGMShell *eShell, WlzGMShell *dShell)
{
  /* Update extended shell's geometry. */
  WlzGMShellMergeG(eShell, dShell);
  /* Update extended shell's topology. */
  WlzGMShellMergeT(eShell, dShell);
  /* Unlink the now childless shell. */
  WlzGMShellUnlink(dShell);
}

/*!
* \return				New resource index table, NULL
*                                       on error.
* \ingroup      WlzGeoModel
* \brief	Makes an index look up table data structure for the
*               given model.
* \param	model			The model with the vertex
*                                       resources.
* \param	eMsk			Element mask with bits set for
*                                       elemet resources to index.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzGMResIdxTb	*WlzGMModelResIdx(WlzGMModel *model, unsigned int eMsk,
				  WlzErrorNum *dstErr)
{
  int		dim,
  		iCnt,
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
    dim = WlzGMModelGetDimension(model, &errNum);
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
  if((errNum == WLZ_ERR_NONE) && (eMsk & WLZ_GMELMFLG_FACE))
  {
    eMskTst &= ~(WLZ_GMELMFLG_FACE);
    if(dim == 3)
    {
      vec = model->res.face.vec;
      resIdxTb->face.idxCnt = model->res.face.numIdx;
      if((resIdxTb->face.idxLut = (int *)
	    AlcCalloc(resIdxTb->face.idxCnt, sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
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
    if((eMsk & WLZ_GMELMFLG_FACE) && (dim == 3))
    {
      /* Compute face index lut. */
      iCnt = vCnt = 0;
      vec = model->res.face.vec;
      iLut = resIdxTb->face.idxLut;
      while(vCnt < resIdxTb->face.idxCnt)
      {
	eP.face = (WlzGMFace *)AlcVectorItemGet(vec, vCnt++);
	if(eP.face->idx >= 0)
	{
	  *(iLut + eP.face->idx) = iCnt++;
	}
      }
      resIdxTb->face.idxCnt = iCnt;
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

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Frees a GM index look up table data structure.
* \param	resIdxTb		Given index lut data structure.
*/
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
    if(resIdxTb->face.idxLut)
    {
      AlcFree(resIdxTb->face.idxLut);
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Rehash the vertex matching hash table.
* \param	model			The model.
* \param	vHTSz			New vertex matching hash table
*                                       size, no change if <= 0.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Constructs a new shell, face and 3 edges (with
*		topological elements) and adds them to the model. Given
*		3 3D double precision points which are known not to
*		be in the model.
* 
*   	        Constructs a new shell, face and 3 edges (with
*		topological elements) and adds them to the model. Given
*		3 3D double precision points which are known not to
*		be in the model.
*		None of the verticies already exists within the shell.
*		Need to create a new shell (nShell) with: 1 face (nF),
*		2 loop topology element (nLT0, nLT1), 3 edges (*nE),
*		6 edge topology elements (*nET0, *nET1), 3 disk
*		topology elements *(nDT), 3 verticies (*nV) with
*		geometry (*pos) and 6 vertex topology elements
*		(*nVT0, *nVT1).
* \verbatim
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
* \endverbatim
*                                                                      
* \param	model			The model to add the segment to.
* \param	pos			Positions of the 3 points.
*/
static WlzErrorNum WlzGMModelConstructNewS3D(WlzGMModel *model,
					     WlzDVertex3 *pos)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell,
  		*nShell;
  WlzGMFace	*nF;
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
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    nLT[0]->next = nLT[0]->prev = nLT[1];
    nLT[1]->next = nLT[1]->prev = nLT[0];
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Extends an existing shell within the model by adding a:
*		face, 3 edges and 2 verticies (with topological elements).
*
*       	Extends an existing shell within the model by adding a:
*		loop, 3 edges and 2 verticies (with topological elements).
*		The 2 given 3D double precision points which are known
*		not to be in the model.
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology elements
*		(nDT[0,1,2]), 2 verticies (nV[0,1]) with geometry
*		(pos0, pos1) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eV			Shared vertex.
* \param	pos0			Position of the first point.
* \param	pos1			Position of the second point.
*/
static WlzErrorNum WlzGMModelExtend1V0E1S3D(WlzGMModel *model, WlzGMVertex *eV,
					    WlzDVertex3 pos0, WlzDVertex3 pos1)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMFace	*nF;
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
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Extends an existing shell within the model by adding a:
*		face, 2 edges and 1 vertex (with topological elements).
*
*       	Extends an existing shell within the model by adding a:
*		face, 2 edges and 1 vertex (with topological elements).
*		The given 3D double precision point is known not to be
*		in the model.
*		Need to create: 1 face (nF), 2 loop topology elements
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eE			Shared edge.
* \param	nPos			Position of new vertex.
*/
static WlzErrorNum WlzGMModelExtend2V1E1S3D(WlzGMModel *model, WlzGMEdge *eE,
					    WlzDVertex3 nPos)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMFace	*nF;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[2];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT;
  WlzGMVertex	*nV;
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update shell */
    (void )WlzGMShellUpdateG3D(eShell, nPos);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Extends an existing shell within the model by adding a:
*		face, 3 edges and 1 vertex (with topological elements).
* 
*        	Extends an existing shell within the model by adding a:
*		face, 3 edges and 1 vertex (with topological elements).
*		Both the given 3D double precision points are known not
*		to be in the model.
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]), 1 vertex (nV) with geometry (pos0) and
*		6 vertex topology elements (nVT[0,1,2] and nVT1[0,1,2]).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eV0			First shared vertex.
* \param	eV1			Second shared vertex.
* \param	nPos			Position of new vertex.
*/
static WlzErrorNum WlzGMModelExtend2V0E1S3D(WlzGMModel *model,
					    WlzGMVertex *eV0, WlzGMVertex *eV1,
					    WlzDVertex3 nPos)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMFace	*nF;
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
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update shell */
    (void )WlzGMShellUpdateG3D(eShell, nPos);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Joins two existing shells within the model by adding:
*		1 face, 3 edges and 1 vertex (with topological elements).
*
*		Joins two existing shells within the model by adding:
*		1 face, 3 edges and 1 vertex (with topological elements).
*		Both the given 3D double precision points are known not
*		to be in the model.
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]), 1 vertex (nV) with geometry (pos0) and
*		6 vertex topology elements (nVT[0,1,2] and nVT1[0,1,2]).
*		Need to delete 1 shell (dShell).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eV0			First shared vertex.
* \param	eV1			Second shared vertex.
* \param	nPos			Position of new vertex.
*/
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
  WlzGMFace	*nF;
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
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update the extended shell and delete the joined shell. */
    (void )WlzGMShellUpdateG3D(eShell, nPos);
    if(model->child == dShell)
    {
      model->child = eShell;
    }
    WlzGMShellJoinAndUnlink(eShell, dShell);
    WlzGMModelFreeS(model, dShell);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Joins three existing shells within the model by adding:
*		1 face and 3 edges (with topological elements).
*
*		Joins three existing shells within the model by adding:
*		1 face and 3 edges (with topological elements).
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*		Need to delete 2 shells (dShell[0,1]).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	gEV			The three given shared verticies.
*/
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
  WlzGMFace	*nF;
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
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update the extended shell and delete the joined shells. */
    for(idx = 0; idx < 2; ++idx)
    {
      if(model->child == dShell[idx])
      {
	model->child = eShell;
      }
      WlzGMShellJoinAndUnlink(eShell, dShell[idx]);
      WlzGMModelFreeS(model, dShell[idx]);
    }
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Joins two existing shells within the model by adding:
*		1 face and 3 edges (with topological elements).
*
*		Joins two existing shells within the model by adding:
*		1 face and 3 edges (with topological elements).
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*		Need to delete 1 shell (dShell).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eV			The three given shared verticies.
*/
static WlzErrorNum WlzGMModelJoin3V0E2S3D(WlzGMModel *model, WlzGMVertex **eV)
{
  int		idx,
  		nIdx,
		pIdx;
  double	sVol[3];
  WlzGMShell	*dShell,
  		*eShell;
  WlzGMShell	*tShell[3];
  WlzGMFace	*nF;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[3];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT[3];
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    dShell = (tShell[idx] == tShell[(idx + 1) % 3])?
             tShell[(idx + 2) % 3]: tShell[(idx + 1) % 3];
    if(tShell[0] == tShell[1])
    {
      eShell = tShell[0];
      dShell = tShell[2];
    }
    else if(tShell[1] == tShell[2])
    {
      eShell = tShell[1];
      dShell = tShell[0];
    }
    else /* tShell[2] == tShell[0] */
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update the extended shell and delete the joined shell. */
    if(model->child == dShell)
    {
      model->child = eShell;
    }
    WlzGMShellJoinAndUnlink(eShell, dShell);
    WlzGMModelFreeS(model, dShell);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Extends an existing shell within the model by adding:
*		1 face and 3 edges (with topological elements).
*
*		Extends an existing shell within the model by adding:
*		1 face and 3 edges (with topological elements).
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 3 edges (nE[0,1,2]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eV			The three shared vertex.
*/
static WlzErrorNum WlzGMModelExtend3V0E1S3D(WlzGMModel *model,
					    WlzGMVertex **eV)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMFace	*nF;
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
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* No change to shell. */
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Extends an existing shell within the model by adding:
*		1 face and 2 edges (with topological elements).
*
*		Extends an existing shell within the model by adding:
*		1 face and 2 edges (with topological elements).
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 2 edges (nE[0,1]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eE			The shared edge.
* \param	sV			The shared vertex (not on the
*                                       shared edge).
*/
static WlzErrorNum WlzGMModelExtend3V1E1S3D(WlzGMModel *model,
					    WlzGMEdge *eE, WlzGMVertex *sV)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMFace	*nF;
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
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    WlzGMDiskTAppend(eV[2]->diskT, nDT);
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update shell */
    (void )WlzGMShellUpdateG3D(eShell, pos[2]);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Joins two existing shells within the model by adding:
*		1 face and 2 edges (with topological elements).
*
*		Joins two existing shells within the model by adding:
*		1 face and 2 edges (with topological elements).
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 2 edges (nE[0,1]), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]), 3 disk topology element
*		(nDT[0,1,2]) and 6 vertex topology elements (nVT[0,1,2]
*		and nVT1[0,1,2]).
*		Need to delete 1 shell (dShell).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eE			The shared edge.
* \param	sV			The shared vertex (not on the
*                                       shared edge).
*/
static WlzErrorNum WlzGMModelJoin3V1E2S3D(WlzGMModel *model,
				          WlzGMEdge *eE, WlzGMVertex *sV)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*dShell,
  		*eShell;
  WlzGMFace	*nF;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE[2];
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*nDT;
  WlzGMVertex	*eV[3];
  WlzGMVertexT	*tVT;
  WlzGMVertexT	*nVT0[3],
  		*nVT1[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* New vertex topology elements. */
    tVT = eE->edgeT->vertexT;
    WlzGMVertexTAppend(tVT, nVT0[0]);
    WlzGMVertexTAppend(tVT, nVT1[0]);
    nVT0[0]->diskT = nVT1[0]->diskT = tVT->diskT;
    tVT = eE->edgeT->opp->vertexT;
    WlzGMVertexTAppend(tVT, nVT0[1]);
    WlzGMVertexTAppend(tVT, nVT1[1]);
    nVT0[1]->diskT = nVT1[1]->diskT = tVT->diskT;
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
    WlzGMDiskTAppend(eV[2]->diskT, nDT);
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
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nF;
    nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = eShell;
    nLT[1]->parent = eShell;
    /* Update extended shell and delete the joined shell. */
    if(model->child == dShell)
    {
      model->child = eShell;
    }
    WlzGMShellJoinAndUnlink(eShell, dShell);
    WlzGMModelFreeS(model, dShell);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Extends a shell within the model by adding:
*		1 face and 1 edge (with topological elements).
* 
*		Extends a shell within the model by adding:
*		1 face and 1 edge (with topological elements).
*		Need to create: 1 face (nF), 2 loop topology elements
*		(nLT[0,1]), 1 edge (nE), 6 edge topology elements
*		(nET0[0,1,2], nET1[0,1,2]) and 6 vertex topology elements
*		(nVT[0,1,2] and nVT1[0,1,2]).
* \endverbatim
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
* \endverbatim
* \param	model			The model.
* \param	eE0			First shared edge.
* \param	eE1			Second shared edge.
*/
static WlzErrorNum WlzGMModelExtend3V2E1S3D(WlzGMModel *model,
					    WlzGMEdge *eE0, WlzGMEdge *eE1)
{
  int		idx,
  		nIdx,
		pIdx;
  WlzGMShell	*eShell;
  WlzGMFace	*nF;
  WlzGMLoopT	*nLT[2];
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0[3],
  		*nET1[3];
  WlzGMDiskT	*eDT[3];
  WlzGMVertex	*eV[3];
  WlzGMVertexT	*eVT[3],
  		*nVT0[3],
  		*nVT1[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
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
    /* Get verticies and disk topology elements. */
    eV[1] = WlzGMEdgeCommonVertexGetDiskTs(eE0, eE1, eDT + 0, eDT + 1);
    if(eDT[0] != eDT[1])
    {
      /* Join the disk topology elements. */
      WlzGMDiskTJoin(eDT[1], eDT[0]);
      WlzGMModelFreeDT(model, eDT[0]);
    }
    eDT[0] = eE0->edgeT->vertexT->diskT;
    eV[0] = eDT[0]->vertex;
    if(eV[0] == eV[1])
    {
      eDT[0] = eE0->edgeT->opp->vertexT->diskT;
      eV[0] = eDT[0]->vertex;
    }
    eDT[2] = eE1->edgeT->vertexT->diskT;
    eV[2] = eDT[2]->vertex;
    if(eV[2] == eV[1])
    {
      eDT[2] = eE1->edgeT->opp->vertexT->diskT;
      eV[2] = eDT[2]->vertex;
    }
    /* Set up new edge and vertex topology elements. */
    for(idx = 0; idx < 3; ++idx)
    {
      nIdx = (idx + 1) % 3;
      pIdx = (idx + 3 - 1) % 3;
      nVT0[idx]->parent = nET0[idx];
      nVT0[idx]->diskT = eDT[idx];
      WlzGMVertexTAppend(eDT[idx]->vertexT, nVT0[idx]);
      nVT1[idx]->parent = nET1[idx];
      nVT1[idx]->diskT = eDT[idx];
      WlzGMVertexTAppend(eDT[idx]->vertexT, nVT1[idx]);
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
    /* New and existing edges. */
    nE->edgeT = nET0[2];
    nET1[1]->edge = nET0[0]->edge = eE0;
    nET1[2]->edge = nET0[1]->edge = eE1;
    nET1[0]->edge = nET0[2]->edge = nE;
    /* Set radial edges. */
    WlzGMEdgeTInsertRadial(nET0[0]);
    WlzGMEdgeTInsertRadial(nET1[1]);
    WlzGMEdgeTInsertRadial(nET0[1]);
    WlzGMEdgeTInsertRadial(nET1[2]);
    /* New face and loop topology elements. */
    nF->loopT = nLT[0];
    nLT[0]->opp = nLT[1];
    nLT[1]->opp = nLT[0];
    nLT[0]->face = nLT[1]->face = nF;
    nLT[0]->edgeT = nET0[0];
    nLT[1]->edgeT = nET1[0];
    nLT[0]->parent = nLT[1]->parent = eShell;
    WlzGMLoopTAppend(eShell->child, nLT[0]);
    WlzGMLoopTAppend(eShell->child, nLT[1]);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzGeoModel
* \brief	If it doesn't already exist in the model, constructs
*		a new 3D simplex splices it into the given model.
*
*		Creates a new simplex has: A face (nF), 2 loop topology
*		elements (nLT0, nLT1), 6 edge topology elements (*nET0,
*		*nET1) and 6 vertex topology elements (*nVT0, *nVT1).
*		No disk topology elements, vertices or edges are
*		created by this function, nor is a shell.
* 
* \verbatim
*                           nVT00        nET00
*                  eDT0   O ------------------------------->        eDT1
*                   eV0 @- - - - - - - - - eE0 - - - - - - - -@ eV1
*                        \   <--------------------------- O  /
*                       ^   O nVT10      nET11     nVT11  ^   O nVT01
*                        \ \ \                           / / /
*                         \   \                         /   /
*                    e2 -----\ \                       / /----- e1
*                           \   \ nET10               /   /
*                            \ \ \        lt1        / / /      
*                             \   \                 /   /
*                        nET02 \ \ \         nET12 / / / nET01
*                               \   \             /   /
*                                \ \ \           / / /
*                                 \   \ nVT12   /   /
*                      lt0         \ \ \    |  / / /
*                                   \   \   | /   /
*                                    \ \ \  |/ / /
*                               nVT02 O   v O   v
*                                        \   /
*                                          @  eV2
*                                            eDT2
* \endverbatim
*                                             
*                   LT0 = {nET00, nET01, nET02},
*                   LT1 = {nET10, nET11, nET12}
*
* \param	model			The geometric model.
* \param	eE			Array of three existing edges to be
*					shared by the new facet.
*/
static WlzErrorNum WlzGMModelExtend3V3E1S3D(WlzGMModel *model, WlzGMEdge **eE)
{
  int		idx,
		pIdx,
		nIdx;
  WlzGMVertexT	*nVT0[3],
		*nVT1[3];
  WlzGMVertex	*eV;
  WlzGMDiskT	*eDT[2];
  WlzGMEdgeT	*nET0[3],
		*nET1[3];
  WlzGMLoopT	*nLT[2];
  WlzGMFace	*nF = NULL;
  WlzGMShell	*eShell;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check simplex is not already model. */
  if(WlzGMEdgeCommonFace(eE[0], eE[1]) == NULL)
  {
    /* Make a new simplex. */
    if(((nF = WlzGMModelNewF(model, &errNum)) != NULL) &&
       ((nLT[0] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
       ((nLT[1] = WlzGMModelNewLT(model, &errNum)) != NULL) &&
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
      eShell = WlzGMEdgeGetShell(eE[0]);
      /* For each vertex. */
      for(idx = 0; idx < 3; ++idx)
      {
	nIdx = (idx + 1) % 3;
	pIdx = (idx + 3 - 1) % 3;
	eV = WlzGMEdgeCommonVertexGetDiskTs(eE[idx], eE[pIdx],
						 eDT + 0, eDT + 1);
	if(eDT[0] != eDT[1])
	{
	  /* Join the disk topology elements. */
	  WlzGMDiskTJoin(eDT[0], eDT[1]);
	  WlzGMModelFreeDT(model, eDT[1]);
	}
	/* New vertex topology elements */
	nVT0[idx]->parent = nET0[idx];
	nVT0[idx]->diskT = eDT[0];
        WlzGMVertexTAppend(eDT[0]->vertexT, nVT0[idx]);
	nVT1[idx]->parent = nET1[idx];
	nVT1[idx]->diskT = eDT[0];
        WlzGMVertexTAppend(eDT[0]->vertexT, nVT1[idx]);
	/* New edge topology elements */
	nET0[idx]->next = nET0[nIdx];
	nET0[idx]->prev = nET0[pIdx];
	nET0[idx]->opp = nET1[nIdx];
	nET0[idx]->edge = eE[idx];
	nET0[idx]->vertexT = nVT0[idx];
	nET0[idx]->parent = nLT[0];
	nET1[idx]->next = nET1[pIdx]; /* Previous because reverse direction. */
	nET1[idx]->prev = nET1[nIdx];     /* Next because reverse direction. */
	nET1[idx]->opp = nET0[pIdx];
	nET1[idx]->edge = eE[pIdx];
	nET1[idx]->vertexT = nVT1[idx];
	nET1[idx]->parent = nLT[1];
      }
      for(idx = 0; idx < 3; ++idx)
      {
        WlzGMEdgeTInsertRadial(nET0[idx]);
        WlzGMEdgeTInsertRadial(nET1[idx]);
      }
      /* New loop topology elements. */
      nLT[0]->opp = nLT[1];
      nLT[0]->face = nF;
      nLT[0]->edgeT = nET0[0];
      nLT[0]->parent = eShell;
      WlzGMLoopTAppend(eShell->child, nLT[0]);
      nLT[1]->opp = nLT[0];
      nLT[1]->face = nF;
      nLT[1]->edgeT = nET1[0];
      nLT[1]->parent = eShell;
      WlzGMLoopTAppend(eShell->child, nLT[1]);
      /* New face. */
      nF->loopT = nLT[0];
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzGeoModel

/*!
* \return				Woolz error code.
* \brief	Constructs a new shell, edge (with topological
*		elements) and adds them to the model.
*
*		Constructs a new shell, edge (with topological
*		elements) and adds them to the model.
*		Given a pair of 2D double precision edge end points which
*		are known not to be in the model.
*		Neither of the verticies already exists within the shell.
*		Need to create a new shell with: 1 loop
*		topology element (nLT), 1 edge (nE) 2 edge
*		topology elements (nET[0,1]), 2 disk topology elements,
*		(nDT[0,1]), 2 verticies (nV[0,1] with geometry
*		pos[0,1]) and 2 vertex topology elements (nVT[0,1]).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	pos			Ptr to 1st then 2nd point.
*/
static WlzErrorNum WlzGMModelConstructNewS2D(WlzGMModel *model,
					     WlzDVertex2 *pos)
{
  int		idx;
  WlzGMShell	*eShell,
  		*nShell;
  WlzGMLoopT	*nLT;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET[2];
  WlzGMDiskT	*nDT[2];
  WlzGMVertex	*nV[2];
  WlzGMVertexT	*nVT[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model &&
     ((nShell = WlzGMModelNewS(model, &errNum)) != NULL) &&
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
    nET[0]->rad = nET[0];
    nET[0]->next = nET[0]->prev = nET[0]->opp = nET[1];
    nET[0]->vertexT = nVT[0];
    nET[0]->parent = nLT;
    nET[1]->rad = nET[1];
    nET[1]->next = nET[1]->prev = nET[1]->opp = nET[0];
    nET[1]->vertexT = nVT[1];
    nET[1]->parent = nLT;
    /* New loop topology elements */
    nLT->next = nLT->prev = nLT->opp = nLT;
    nLT->face = NULL;
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Constructs a new edge (with topological elements) and
*		adds them to the model.
*
*		Constructs a new edge (with topological elements) and
*		adds them to the model.
*		Given an existing edge topology element at one end
*		of the new edge and a 2D double precision end points
*		at the other end of the new edge, where the new end
*		point is known not to be in the model.
*		Only one vertex already exists within the model so no new
*		shell is required, instead an existing shell is
*		extended using: 1 edge (nE), 2 edge topology elements
*		(nET0, nET1), 1 disk topology element (nDT), 1 vertex
*		(nV0 with geometry (nPos) and 2 vertex topology elements
*		(nVT0, nVT1).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eET0			Next edge topology element for the
*                                       new edge topology element directed
*                                       away from the new end point.
* \param	nPos			The new end point.
*/
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Splits a loop by adding a new edge.
*
*		Splits a loop by adding a new edge.
*		First of all checks that the edge segment doesn't
*		already exist! If it doesn't then a new edge is
*		constructed which splits and existing loop
*		within the model. Given the loop topology element
*		to be split and the existing edge topology elements at
*		new edges end points.
*		Split a loop within a shell by adding a new edge between the
*		two matched verticies.
*		Create: 1 edge (nE), 2 edge topology elements (nET0,
*		nET1), 2 vertex topology elements (nVT0, nVT1) and
*		a loop topology element (nLT).
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eET0			Next edge topology element for the
*                                       first end point.
* \param	eET1			Next edge topology element for the
*                                       other end point.
*/
static WlzErrorNum WlzGMModelConstructSplitL2D(WlzGMModel *model,
					       WlzGMEdgeT *eET0,
					       WlzGMEdgeT *eET1)
{
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0,
  		*nET1;
  WlzGMLoopT	*eLT,
  		*nLT;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check for duplicate edge. */
  if(eET0->prev != eET1->prev->opp)
  {
    if(model && eET0 &&
       ((nE = WlzGMModelNewE(model, &errNum)) != NULL) &&
       ((nET0 = WlzGMModelNewET(model, &errNum)) != NULL) &&
       ((nET1 = WlzGMModelNewET(model, &errNum)) != NULL) &&
       ((nVT0 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
       ((nVT1 = WlzGMModelNewVT(model, &errNum)) != NULL) &&
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
      /* Loop topology element */
      eLT = eET0->parent;
      eLT->edgeT = nET0;
      WlzGMLoopTAppend(eET0->parent, nLT);	       /* Add loopT to shell */
      nLT->opp = nLT;
      nLT->face = NULL;
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Constructs a new edge which joins two (different)
*		existing loops within the model.
*
*		brief	Constructs a new edge which joins two (different)
*		existing loops within the model.
*		Join two loops within two possibly different shells by
*		adding a new edge between the two matched verticies.
*		Create: 1 edge (nE), 2 edge topology elements (nET0,
*		nET1) and 2 vertex topology elements (nVT0, nVT1).
*		Delete: 1 loop topology element (dLT)
*		and 1 shell (dShell). If edges were (they're not)
*		allowed to cross then it would be possible to have a
*		pair of verticies within the same shell but NOT within
*		a common loop. Check for this disallowed case and return
*		an error if it's found.
* \verbatim
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
* \endverbatim
* \param	model			The model to add the segment to.
* \param	eET0			Next edge topology element for the
*                                       first end point.
* \param	eET1			Next edge topology element for the
*                                       other end point.
*/
static WlzErrorNum WlzGMModelJoinL2D(WlzGMModel *model,
				     WlzGMEdgeT *eET0, WlzGMEdgeT *eET1)
{
  WlzGMVertexT	*nVT0,
  		*nVT1;
  WlzGMEdge	*nE;
  WlzGMEdgeT	*nET0,
  		*nET1;
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
    /* Loop topology elements */
    eLT = eET0->parent;
    dLT = eET1->parent;
    dS = dLT->parent;
    /* Look to see if loops share the same shell. */
    if((dS = dLT->parent) == eLT->parent)
    {
      dS = NULL;
    }
    else
    {
      /* Merge the two shell's geometries and topologies. */
      WlzGMShellMergeG(eLT->parent, dS);
      WlzGMShellMergeT(eLT->parent, dS);
    }
    /* Set/update the retained loop topology element by walking around the
     * edge topology  elements setting their parents. */
    WlzGMLoopTSetT(eLT);
    /* Unlink and free the deleted loop topology element. */
    WlzGMLoopTUnlink(dLT);
    (void )WlzGMModelFreeLT(model, dLT);
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

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Constructs a 2D simplex (edge) defined by two double
*               precision end points. Either of the two points may
*               already exist within the model.
*               See WlzGMShellMatchVtxG2D() for the meaning of the
*               backwards, forwards and distance search parameters.
* \param	model			The model to add the segment to.
* \param	pos			Pointer to first then second
*                                       positions.
*/
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
      if((matchEdgeT[0]->edge != matchEdgeT[1]->edge) &&
         (matchEdgeT[0]->prev->edge != matchEdgeT[1]->prev->edge))
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
	   * joined by a new edge, destroying a loop.  */
	  errNum = WlzGMModelJoinL2D(model, matchEdgeT[0], matchEdgeT[1]);
	  
	}
      }
      break;
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Constructs a 3D simplex (triangle) defined by three
*               double precision verticies, any of which may already
*               exist within the model.
* \param	model			The model to add the segment to.
* \param	pos			Pointer to triangle verticies.
*/
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
	   * need to check if the three verticies are in a common loop and then
	   * if ther're not add a new loop, 2 loopT's, 6 edgeT's and 6
	   * vertexT's. */
	  errNum = WlzGMModelExtend3V3E1S3D(model, cE);
	  break;
      }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Walks around a loop topology element's child edge
*               topology elements setting their parent and making
*		sure that the loop topology element's adjacent to
*		it have the correct parent too.
* \param	gLT			The loop topology element.
*/
static void	WlzGMLoopTSetT(WlzGMLoopT *gLT)
{
  WlzGMShell	*gS;
  WlzGMEdgeT	*fET,
  		*tET;
  
  if(gLT && ((tET = fET = gLT->edgeT) != NULL))
  {
    /* Set each of the edgeT's to have the given loopT as it's parent. */
    gS = gLT->parent;
    do
    {
      tET->parent = gLT;
    } while((tET = tET->next) != fET);
    WlzGMLoopTSetAdjT(gLT, gLT->parent);
  }
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Moves all loop topology elements connected to the given
*		loop topology element from their current shell to the
*		iparent shell of the given loop topology element.
* \param	gLT			The loop topology element.
* \param 	gS			The shell to which all adjacent
*					loop topology elements are to be
*					moved.
*/
static void	WlzGMShellSetT(WlzGMShell *gS)
{
  /* TODO Check this fn. */
  if(gS && (gS->child != NULL))
  {
    WlzGMLoopTSetAdjT(gS->child, gS);
  }
}

/*!
* \return	<void>
* \ingroup      WlzGeoModel
* \brief	Recursive function which moves loop topology elements
*		from their current shell to the given shell if they are
*		adjacent to the given loop topology element.
* \param	gLT			The loop topology element.
* \param 	gS			The shell to which all adjacent
*					loop topology elements are to be
*					moved.
*/
static void	WlzGMLoopTSetAdjT(WlzGMLoopT *gLT, WlzGMShell *gS)
{
  WlzGMLoopT	*cLT;
  WlzGMEdgeT	*cET,
  		*fET;

  /* TODO Check this fn. */
  if(gLT->parent != gS)
  {
    WlzGMLoopTUnlink(gLT);
    WlzGMLoopTAppend(gS->child, gLT);
    gLT->parent = gS;
  }
  cET = fET = gLT->edgeT;
  do
  {
    cLT = cET->parent;
    if(cLT->parent != gS)
    {
      WlzGMLoopTSetAdjT(cLT, gS);
    }
    cET = cET->next;
  } while(cET != fET);
}

/*!
* \return				The parent shell or NULL if no
*					parent shell was found.
* \ingroup      WlzGeoModel
* \brief	Walks around a loop topology element looking for another
*		loop topology element connected to the given one. If
*		such a loop topology element exists this function returns
*		it's shell.
* \param	gLT			The given loopT.
*/
static WlzGMShell *WlzGMLoopTFindShell(WlzGMLoopT *gLT)
{
  WlzGMEdgeT	*fET,
  		*tET;
  WlzGMShell	*fS = NULL;

  tET = fET = gLT->edgeT;
  do
  {
    tET = tET->next;
    if(tET->opp->parent != gLT)
    {
      /* LoopTs connected by an edge. */
      fS = tET->opp->parent->parent;
    }
    else if(tET->vertexT->diskT != tET->vertexT->diskT->next)
    {
      /* loopTs connected by a vertex. */
      /* TODO not implemented yet, needed for 3D. */
    }
  } while((fS == NULL) && (tET != fET));
  return(fS);
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Removes a vertex from the models vertex hash table.
*               The vertex's geometry must have been set.
* \param	model			The model.
* \param	dV			The vertex to remove from the
*                                       models hash table.
*/
static void	WlzGMModelRemVertex(WlzGMModel *model, WlzGMVertex *dV)
{
  WlzDVertex3	nPos;
  unsigned int 	hVal;
  WlzGMVertex  	**hdV;
  WlzGMVertex	*tV,
		*pV;

  (void )WlzGMVertexGetG3D(dV, &nPos);
  hVal = WlzGMHashPos3D(nPos);
  hdV = model->vertexHT + (hVal % model->vertexHTSz);
  if(*hdV)
  {
    /* Need to search for position in the linked list. */
    tV = *hdV;
    pV = NULL;
    while(tV && (WlzGMVertexCmpSign3D(tV, nPos) < 0))
    {
      pV = tV;
      tV = tV->next;
    }
    /* tV is now dV. */
    if(pV == NULL)
    {
      *hdV = dV->next;
    }
    else
    {
      pV->next = dV->next;
    }
    dV->next = NULL;
  }
}

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Adds a new vertex into the models vertex hash table.
*               The vertex's geometry must have been set.
* \param	model			The model.
* \param	nV			New vertex to insert into the
*                                       models hash table.
*/
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
    /* No vertex at this position in the hash table yet. */
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

/*!
* \return				<void>
* \ingroup      WlzGeoModel
* \brief	Attempts to find verticies which match the two given
*               double precision positions.
*               The linking/matching workspace makes the matching much
*               quite efficient, but it relies on the verticies being
*               maintained in a hash table which has indicies which
*		increase with increasing row and then increasing column.
* \param	model			Model with resources.
* \param	matchET			Array for return of matched
*                                       edge topology element pointers.
* \param	pos			Pointer to first then second
*                                       positions.
*/
static void	WlzGMModelMatchEdgeTG2D(WlzGMModel *model,
					WlzGMEdgeT **matchET,
				        WlzDVertex2 *pos)
{
  int		idD;
  double	tD0,
		len,
		trCos,
		trSin;
  WlzDVertex2	tr,
		rPos,
  		tPosB,
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
      if(vertexT == vertexT->next)
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
	 *   {p_x}' = p_x - m_x
	 *   {p_y}' = p_y - m_y
	 *   {p_x}'' = {p_x}'\cos(\theta) - {p_y}'\sin(\theta)
	 *   {p_y}'' = {p_x}'\sin(\theta) + {p_y}'\cos(\theta)
	 *   \sin(\theta) = \frac{{o_x}'}{|o'|}
	 *   \cos(\theta) = \frac{{o_y}'}{|o'|}
	 */
	tr.vtX = mPos.vtX;
	tr.vtY = mPos.vtY;
	rPos.vtX = oPos.vtX - tr.vtX;
	rPos.vtY = oPos.vtY - tr.vtY;
	tD0 = (rPos.vtX * rPos.vtX) + (rPos.vtY * rPos.vtY);
	if(tD0 > tol)
	{
	  /* Compute the transform rotation elements. */
	  len = sqrt(tD0);
	  trCos = rPos.vtX / len;
	  trSin = rPos.vtY / len;
	  /* Search for next edge by walking around the matched vertex CCW,
	   * using the vertex topology elements, comparing the relative
	   * angles of the line segments between the verticies.
	   * The next edge topology element is the one with the smallest
	   * CCW polar angle relative to the new, as yet unmade, edge. */
	  eTB = eTS;	  /* Initialize the search, best (eTB) = start (eTS) */
	  (void )WlzGMVertexGetG2D(eTB->opp->vertexT->diskT->vertex, &tPos);
	  rPos.vtX = tPos.vtX - tr.vtX;
	  rPos.vtY = tPos.vtY - tr.vtY;
	  tPosB.vtX =  (rPos.vtX * trCos) + (rPos.vtY * trSin);
	  tPosB.vtY = -(rPos.vtX * trSin) + (rPos.vtY * trCos);
	  /* Find next edge topology element (eTN) directed from the vertex.
	   * Because the edge topology elements are unordered need to
	   * search them all. */
	  eTN = eTB->opp->next;
	  /* For each edge topology element directed away from the vertex */
	  while(eTS != eTN)
	  {
	    (void )WlzGMVertexGetG2D(eTN->opp->vertexT->diskT->vertex, &tPos);
	    rPos.vtX = tPos.vtX - tr.vtX;
	    rPos.vtY = tPos.vtY - tr.vtY;
	    tPosN.vtX =  (rPos.vtX * trCos) + (rPos.vtY * trSin);
	    tPosN.vtY = -(rPos.vtX * trSin) + (rPos.vtY * trCos);
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

/* Model Features */

/*!
* \return				Number of simplicies in given shell.
* \ingroup      WlzGeoModel
* \brief	Counts the number of simplicies in the given shell.
*               For 2D models the simplicies are edges and for 3D models
*               they are loops.
* \param	shell			The given shell.
*/
int		WlzGMShellSimplexCnt(WlzGMShell *shell)
{
  int		dim,
  		cnt = 0;
  WlzGMEdgeT	*cET,
  		*fET;
  WlzGMLoopT	*cLT,
  		*fLT;

  if((shell != NULL) && (shell->idx >= 0) && (shell->parent != NULL) &&
     ((dim = WlzGMModelGetDimension(shell->parent, NULL)) != 0) &&
     ((cLT = fLT = shell->child) != NULL))
  {
    if(dim == 2)
    {
      /* For each loop topology element. */
      do
      {
	/* For each edge topology element. */
	cET = fET = cLT->edgeT;
	do
	{
	  cnt += cET->edge->edgeT == cET;
	  cET = cET->next;
	} while(cET != fET);
	cLT = cLT->next;
      } while(cLT != fLT);
    }
    else /* dim == 3 */
    {
      /* For each loop topology element. */
      do
      {
	cnt += cLT->face->loopT == cLT;
	cLT = cLT->next;
      } while(cLT != fLT);
    }
  }
  return(cnt);
}
