#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzGeoModelFilters.c
* \author       Bill Hill
* \date         November 2000
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Filters for Woolz geometric models (GM's).
* \ingroup      WlzGeoModel
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>
#include <limits.h>
#include <float.h>

static void			WlzGMFilterGeomLPL2D(
				  WlzGMModel *model,
				  WlzDVertex2 *vGIn,
				  WlzDVertex2 *vGOut,
				  double lambda,
				  int nonMan);
static void			WlzGMFilterGeomLPL3D(
				  WlzGMModel *model,
				  WlzDVertex3 *vGIn,
				  WlzDVertex3 *vGOut,
				  double lambda,
				  int nonMan);
static WlzDVertex2 		WlzGMFilterGeomLPL2Delta(
				  WlzGMModel *model,
				  WlzGMVertex *gV,
				  WlzDVertex2 *vGIn,
				  int nonMan);
static WlzDVertex3 		WlzGMFilterGeomLPL3Delta(
				  WlzGMModel *model,
				  WlzGMVertex *gV,
				  WlzDVertex3 *vGIn,
				  int nonMan);
static WlzErrorNum 		WlzGMFilterGeomSetVertices(
				  WlzGMModel *model,
				  WlzVertexP vtxBuf);

/*!
* \return	Woolz error code.
* \brief	Flips the orientation of edges in 2D models and
*		faces in 3D models. This is done in place.
* \param	model			Given model.
*/
WlzErrorNum	WlzGMFilterFlipOrient(WlzGMModel *model)
{
  int		eIdx,
  		eCnt;
  AlcVector	*eVec;
  WlzGMEdge	*edge;
  WlzGMFace	*face;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(model->type)   
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N:  
	eIdx = 0;
	eCnt = model->res.edge.numIdx;
	eVec = model->res.edge.vec;
	while(eIdx < eCnt)
	{
	  edge = (WlzGMEdge *)AlcVectorItemGet(eVec, eIdx);
	  if(edge && (edge->idx >= 0))
	  {
	    edge->edgeT = edge->edgeT->opp;
	  }
	  ++eIdx;
	}
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:  
	eIdx = 0;
	eCnt = model->res.face.numIdx;
	eVec = model->res.face.vec;
	while(eIdx < eCnt)
	{
	  face = (WlzGMFace *)AlcVectorItemGet(eVec, eIdx);
	  if(face && (face->idx >= 0))
	  {
	    face->loopT = face->loopT->opp;
	  }
	  ++eIdx;
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
* \return	Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Removes small shells from the given geometric model.
* \param	model			Given model.
* \param	minSpx			Minimum number of simplicies
*                                       (edges in a 2D model or faces
*                                       in a 3D model) in a shell
*                                       for it to stay in the model.
*/
WlzErrorNum	WlzGMFilterRmSmShells(WlzGMModel *model, int minSpx)
{
  WlzGMShell	*cS,
  		*fS,
		*nS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((fS = model->child) != NULL)
  {
    nS = fS->next;
    do
    {
      cS = nS;
      nS = nS->next;
      if(WlzGMShellSimplexCnt(cS) < minSpx)
      {
         WlzGMModelDeleteS(model, cS);
      }
    } while(cS != fS);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Performs a low pass filtering operation on the geometry of 
*		the model, given the band pass and band stop filter
*		parameters. These parameters specify the transfer function
*		of the low pass filter.
*                \li \f$k_{PB}\f$         The band pass frequency parameter 
*                                       which is the interface between the
*					pass and transition bands of the
*					filter.
*                \li \f$k_{SB}\f$     	The stop band frequency parameter
*					which is the interface between the
*					transition and stop bands of the
*					filter.
*                \li \f$\delta_{PB}\f$ 	The band pass delta parameter
*					which is the maximum deviation
*					from a flat response in the band
*					pass region.
*                \li \f$\delta_{SB}\f$ 	The stop band delta parameter
*					which is the maximum deviation
*					from a flat response in the stop
*					band region.
* \note		See WlzGMFilterGeomLPParam() and WlzGMFilterGeomLPLM().
* \param	model			Given model.
* \param	kPB			The band pass frequency parameter.
* \param	kSB			The band stop frequency parameter.
* \param	dPB			The pass band maximum deviation.
* \param	dSB			The stop band maximum deviation.
* \param	maxItr			Minimum number of itterations.
* \param	nonMan			If non-zero allows non manifold
* 					vertices to be filtered.
*/
WlzErrorNum	WlzGMFilterGeomLP(WlzGMModel *model, double kPB, double kSB,
				  double dPB, double dSB, int maxItr,
				  int nonMan)
{
  int		nItr;
  double	lambda,
  		mu;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    errNum = WlzGMFilterGeomLPParam(&lambda, &mu, &nItr, kPB, kSB, dPB, dSB);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nItr > maxItr)
    {
      nItr = maxItr;
    }
    errNum = WlzGMFilterGeomLPLM(model, lambda, mu, nItr, nonMan);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Computes the values of \f$\lambda\f$ and \f$\mu\f$ from the
*		given low pass filter parameters \f$k_{PB}\f$,  \f$k_{SB}\f$,
*		\f$\delta_{PB}\f$ and \f$\delta_{SB}\f$.
*		Valid filter parameters require:
*		\f$0 < k_{PB} < 1\f$, \f$k_{PB} < k_{SB} < 2\f$,
*		\f$0 < \delta_{PB} < 1\f$ and \f$0 < \delta_{SB} < 1\f$.
*		If the parameters are outside of these ranges an error is
*		returned, but these constrants are not sufficient for a good
*		filter. Design tips are don't make the transition band too
*		narrow (keep \f$k_{PB}\f$ away from \f$k_{SB}\f$) and
*		don't make \f$\delta_{PB}\f$ or \f$\delta_{SB}\f$ too small
*		(try values of 0.1 for \f$\delta_{PB}\f$ and
*		\f$\delta_{SB}\f$). Good values for \f$k_{PB}\f$ are
*		probably in the range \f$[0.01-0.1]\f$.
*		The values of \f$\lambda\f$ and \f$\mu\f$ in the transfer
*		function of the filter
*		\f[
		    f(k) = ((1 - {\lambda}k)(1 - {\mu}k))^N
		\f]
*		are constrained
*		by:
*		\f[
		    0 < N
		\f]
*		\f[
		    0 < \lambda < -\mu < 1
		\f]
*		\f[
		    \lambda < \frac{1}{k_{SB}}
		\f]
*		\f[
		    \frac{1}{\lambda} + \frac{1}{\mu} = k_{PB}
		\f]
*		\f[
		    {\left(\frac{{(\lambda - \mu)}^2}
		    		{-4\lambda\mu}\right)}^N < 1 + \delta_{PB}
		\f]
*		\f[
                    {\left(\frac{1 - \lambda k_{SB}}
		                {1 - \mu k_{SB}}\right)}^N < \delta_{SB}
		\f]
*		Following Taubin's fairing paper an additional constraint
*		is imposed \f$f(1) = -f(2)\f$, which allows a simple
*		algebraic solution for \f$\lambda\f$, \f$\mu\f$. The
*		value of \f$N\f$ is then computed by taking the maximum
*		of \f$N_{PB}\f$ and \f$N_{SB}\f$.
*		So the algorithm used to compute the values of \f$\lambda\f$,
*		\f$\mu\f$ and \f$N\f$ is:
*		\f[
		    \lambda = \frac{1}{3 k_{PB} - 5}
		              (k_{PB} - \sqrt{k_{PB}^2 - 6 k_{PB} + 10})
		\f]
*		\f[
		    \mu = \frac{1}{3 k_{PB} - 5}
		          (k_{PB} + \sqrt{k_{PB}^2 - 6 k_{PB} + 10})
		\f]
*		\f[
		    N_{PB} = \frac{\ln{(1 + \delta_{PB})}}
				  {\ln{({(\lambda - \mu)}^2 (-4 \lambda \mu))}}
		\f]
*		\f[
		    N_{SB} = \frac{\ln{(\delta_{SB})}}
				  {\ln{((1 - \lambda k_{SB})(1 - \mu k_{SB}))}}
		\f]
*		\f[
		    N = \max{(N_{PB}, N_{SB})}
		\f]
* \note		See WlzGMFilterGeomLP() and WlzGMFilterGeomLPLM().
* \param	dstLambda		Destination pointer for the positive
*					scale factor, must NOT be NULL.
* \param	dstMu			Destination pointer for the negative
*					scale factor, must NOT be NULL.
* \param	dstNItr			Destination pointer for the number
*					of itterations, may be NULL.
* \param	kPB			The band pass frequency parameter.
* \param	kSB			The band stop frequency parameter.
* \param	dPB			The pass band maximum deviation.
* \param	dSB			The stop band maximum deviation.
*/
WlzErrorNum	WlzGMFilterGeomLPParam(double *dstLambda, double *dstMu,
				       int *dstNItr,
				       double kPB, double kSB,
				       double dPB, double dSB)
{
  int		nItr;
  double	tD0,
  		tD1,
		lambda,
  		mu,
		nPB,
		nSB;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((kPB < DBL_EPSILON) || (kPB > 1.0) || (kPB > kSB) ||
     (kSB > (2.0 - DBL_EPSILON)))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    tD0 = sqrt(kPB * (kPB - 6.0) + 10.0);
    tD1 = 1.0 / ((3 * kPB) - 5.0);
    lambda = tD1 * (kPB - tD0);
    mu = tD1 * (kPB + tD0);
    tD0 = lambda - mu;
    nPB = log(1.0 + dPB) / log(-4.0 * tD0 * tD0 * lambda * mu);
    nSB = log(dSB) / log((1.0 - (lambda * kSB)) * (1.0 - (mu * kSB)));
    tD0 = (nPB > nSB)? nPB: nSB;
    nItr = (int )ceil(tD0);
    *dstLambda = lambda;
    *dstMu = mu;
    if(dstNItr)
    {
      *dstNItr = nItr;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeoModel
* \brief	Given values of \f$\lambda\f$, \f$\mu\f$ and the number of
*		itterations \f$N\f$, applies a low pass filter to the geometry
*		of the given model.
*
*		This filter is based on two papers by Taubin:
*		\li Gabriel Taubin. Curve and Surface Smoothing without
*		    Shrinkage. International Conference on Computer
*		    Vision (ICCV'95) Conference Procedings, 827-857, 1995.
*		\li Gabriel Taubin. A Signal Processing Approach to Fair
*		    Surface Design. Computer Graphics, 29, 351-358, 1995.
*
*		The filter repeatedly applies a pair of filters to the
*		geometry of the model vertices, of the form:
*		\f[
		{v_i}' = v_i + \lambda \Delta v_i
		\f]
*		and
*		\f[
		{v_i}'' = {v_i}' + \mu \Delta {v_i}'
		\f]
*		Where \f$v_i\f$ is the geometry of the i'th vertex and
*		\f$Delta v_i\f$ is given by
*		\f[
		\Delta v_i = \sum_{j \in i^*}{w_{ij}(v_j - v_i)}
		\f]
*		this is just the weighted sum of weighted differences between
*		the vertex and it's first order neighbourhood, with \f$i^*\f$
*		being the number of neighbours. The weights \f$w_{ij}\f$ are
*		computed using
*		\f[
		w_{ij} = \frac{1}{i^*}
		\f]
*		If the topology of a vertex is manifold, then Taubin's filter
*		is applied as set out in the two papers. But if the
*		topology is non-manifold different stratagies are applied
*		to prevent the shrinkage (in non-manifold regions) which
*		Taubin's filter would cause.
*		\li For 2D models, non-manifold vertices, with either a
*		    single incident edge or more than two incident edges
*		    have their geometry preserved.
*		\li For 3D models, vertices on non-manifold edges only
*		    those verticies along the edge are considered to be
*		    neibours and vertices which are themselves non-manifold
*		    have their geometry preserved.
*		The implementation of the filter makes use of a pair of
*		buffers which contain the vertex positions during the
*		itterations. Only when the filter has completed the
*		itterations are the vertex hash table and shell geometries
*		are recomputed.
* \note		See WlzGMFilterGeomLP() and WlzGMFilterGeomLPParam().
* \param	model			The geometric model.
* \param	lambda			Positive filter parameter.
* \param	mu			Negative filter parameter.
* \param	nItr			Number of itterations.
* \param	nonMan			If non-zero allows non manifold
					vertices to be filtered.
*/
WlzErrorNum 	WlzGMFilterGeomLPLM(WlzGMModel *model, double lambda, double mu,
				    int nItr, int nonMan)
{
  int 		itr,
  		dim,
		nVtx;
  WlzVertexP	vtxBuf[2];
  WlzVertexType	vtxType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vtxBuf[0].v = vtxBuf[1].v = NULL;
  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((lambda < DBL_EPSILON) || (lambda > (1.0 - DBL_EPSILON)) ||
          (-mu < DBL_EPSILON) || (-mu > (1.0 - DBL_EPSILON)) ||
          (nItr < 1))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    dim = WlzGMModelGetDimension(model, &errNum);
  }
  /* Allocate a pair of buffers for the vertex positions, with the first
   * buffer filled and the second just allocated. */
  if(errNum == WLZ_ERR_NONE)
  {
    vtxBuf[0] = WlzDVerticesFromGM(model, &nVtx, &vtxType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(vtxType)
    {
      case WLZ_VERTEX_D2:
        if((vtxBuf[1].d2 = (WlzDVertex2 *)
		           AlcMalloc(sizeof(WlzDVertex2) * nVtx)) == NULL)
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_VERTEX_D3:
        if((vtxBuf[1].d3 = (WlzDVertex3 *)
		           AlcMalloc(sizeof(WlzDVertex3) * nVtx)) == NULL)
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_DATA;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Perform the itterations. */
    if(dim == 2)
    {
      for(itr = 0; itr < nItr; ++itr)
      {
        WlzGMFilterGeomLPL2D(model, vtxBuf[0].d2, vtxBuf[1].d2,
			     lambda, nonMan);
        WlzGMFilterGeomLPL2D(model, vtxBuf[1].d2, vtxBuf[0].d2,
			     mu, nonMan);
      }
    }
    else /* dim == 3 */
    {
      for(itr = 0; itr < nItr; ++itr)
      {
        WlzGMFilterGeomLPL3D(model, vtxBuf[0].d3, vtxBuf[1].d3,
			     lambda, nonMan);
        WlzGMFilterGeomLPL3D(model, vtxBuf[1].d3, vtxBuf[0].d3,
			     mu, nonMan);
      }
    }
    errNum = WlzGMFilterGeomSetVertices(model, vtxBuf[0]);
  }
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzGeoModel
* \brief	Filters the geometry of the verticies in a 2D model using
*		the given input and output buffers for the vertex geometries.
* \note		See WlzGMFilterGeomLPLM().
* \param	model			The given model.
* \param	vGIn			Input vertex geometries.
* \param	vGOut			Output vertex geometries.
* \param	lambda			The filter parameter.
* \param	nonMan			If non-zero allows non manifold
					vertices to be filtered.
*/
static void	WlzGMFilterGeomLPL2D(WlzGMModel *model,
				     WlzDVertex2 *vGIn, WlzDVertex2 *vGOut,
				     double lambda, int nonMan)
{
  int		idx,
  		cnt;
  WlzDVertex2	tV0,
  		tV1;
  WlzGMVertex	*cV;
  AlcVector	*vec;

  cnt = model->res.vertex.numIdx;
  vec = model->res.vertex.vec;
  for(idx = 0; idx < cnt; ++idx)
  {
    cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
    if(cV->idx >= 0)
    {
      tV0 = *(vGIn + idx);
      tV1 = WlzGMFilterGeomLPL2Delta(model, cV, vGIn, nonMan);
      tV0.vtX += lambda * tV1.vtX;
      tV0.vtY += lambda * tV1.vtY;
      *(vGOut + idx) = tV0;
    }
  }
}

/*!
* \return	void
* \ingroup	WlzGeoModel
* \brief	Filters the geometry of the verticies in a 3D model using
*		the given input and output buffers for the vertex geometries.
* \note		See WlzGMFilterGeomLPLM().
* \param	model			The given model.
* \param	vGIn			Input vertex geometries.
* \param	vGOut			Output vertex geometries.
* \param	lambda			The filter parameter.
* \param	nonMan			If non-zero allows non manifold
					vertices to be filtered.
*/
static void	WlzGMFilterGeomLPL3D(WlzGMModel *model,
				     WlzDVertex3 *vGIn, WlzDVertex3 *vGOut,
				     double lambda, int nonMan)
{
  int		idx,
  		cnt;
  WlzDVertex3	tV0,
  		tV1;
  WlzGMVertex	*cV;
  AlcVector	*vec;

  cnt = model->res.vertex.numIdx;
  vec = model->res.vertex.vec;
  for(idx = 0; idx < cnt; ++idx)
  {
    cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
    if(cV->idx >= 0)
    {
      tV0 = *(vGIn + idx);
      tV1 = WlzGMFilterGeomLPL3Delta(model, cV, vGIn, nonMan);
      tV0.vtX += lambda * tV1.vtX;
      tV0.vtY += lambda * tV1.vtY;
      tV0.vtZ += lambda * tV1.vtZ;
      *(vGOut + idx) = tV0;
    }
  }
}

/*!
* \return	Vertex displacement.
* \ingroup	WlzGeoModel
* \brief	Computes the displacement of the given vertex, for use
*		by WlzGMFilterGeomLPL2D().
*		With the given vertex \f$v_0\f$ which has \f$i^*\f$
*		neighbours \f$\{v_j\} \: j\in[1 \ldots i^*]\f$ the
*		displacement is \f$ \sum_{j \in i^*}{(v_j - v_i)}\f$
*		if \f$i^* = 2\f$ and \f$0\f$ otherwise.
* \param	model			The given model.
* \param	cV			Given vertex.
* \param	vGIn			Buffer of vertex geometries.
* \param	nonMan			If non-zero allows non manifold
					vertices to be filtered.
*/
static WlzDVertex2 WlzGMFilterGeomLPL2Delta(WlzGMModel *model, WlzGMVertex *gV,
					    WlzDVertex2 *vGIn, int nonMan)
{
  WlzGMEdgeT	*fET,
  		*nET0,
		*nET1;
  WlzDVertex2	*tVGP;
  WlzDVertex2	sVG;

  /* Find the neighbouring verticies and sum their positions. */
  fET = gV->diskT->vertexT->parent->next;
  nET0 = fET->next->opp;
  nET1 = nET0->next->opp;
  if((nET1 == fET) || nonMan)
  {
    sVG = *(vGIn + nET0->vertexT->diskT->vertex->idx);
    tVGP = vGIn + nET1->vertexT->diskT->vertex->idx;
    sVG.vtX = 0.5 * (sVG.vtX + tVGP->vtX);
    sVG.vtY = 0.5 * (sVG.vtY + tVGP->vtY);
    tVGP = vGIn + gV->idx;
    sVG.vtX -= tVGP->vtX;
    sVG.vtY -= tVGP->vtY;
  }
  else
  {
    sVG.vtX = sVG.vtY = 0.0;
  }
  return(sVG);
}

/*!
* \return	Vertex displacement.
* \ingroup	WlzGeoModel
* \brief	Computes the displacement of the given vertex, for use
*		by WlzGMFilterGeomLPL3D().
*		If the model is manifold in the region of the vertex 
*		then the displacement of \f$v_0\f$ with \f$i^*\f$
*		neighbours \f$\{v_j\} \: j\in[1 \ldots i^*]\f$ is just
*		\f$ \sum_{j \in i^*}{(v_j - v_i)}\f$. If the model is
*		not manifold then the displacement is \f$0\f$.
* \param	model			The given model.
* \param	cV			Given vertex.
* \param	vGIn			Buffer of vertex geometries.
* \param	nonMan			If non-zero allows non manifold
					vertices to be filtered.
*/
static WlzDVertex3 WlzGMFilterGeomLPL3Delta(WlzGMModel *model, WlzGMVertex *gV,
					    WlzDVertex3 *vGIn, int nonMan)
{
  int		nN = 0,
  		manifold = 1;
  double	tD0;
  WlzGMVertex	*nV;
  WlzGMEdgeT	*fET,
  		*nET;
  WlzDVertex3	*tVGP;
  WlzDVertex3	nVG,
  		sVG;

  sVG.vtX = sVG.vtY = sVG.vtZ = 0.0;
  manifold = gV->diskT == gV->diskT->next;
  fET = gV->diskT->vertexT->parent->opp;
  if(manifold || nonMan)
  {
    nET = fET;
    do
    {
      ++nN;
      if(!nonMan)
      {
        manifold = (nET != nET->rad) && (nET == nET->rad->rad);
      }
      nV = nET->vertexT->diskT->vertex;
      nVG = *(vGIn + nV->idx);
      WLZ_VTX_3_ADD(sVG, sVG, nVG);
      nET = nET->rad->opp->prev;
    } while((nET != fET) && manifold);
    if(manifold || nonMan)
    {
      tD0 = 1.0 / nN;
      tVGP = vGIn + gV->idx;
      sVG.vtX = (tD0 * sVG.vtX) - tVGP->vtX;
      sVG.vtY = (tD0 * sVG.vtY) - tVGP->vtY;
      sVG.vtZ = (tD0 * sVG.vtZ) - tVGP->vtZ;
    }
    else
    {
      sVG.vtX = sVG.vtY = sVG.vtZ = 0.0;
    }
  }
  return(sVG);
}

/*!
* \return
* \brief	Sets the geometry of the vertices in the geometric model,
*		rebuilds the model's vertex hash table and sets the
*		geometry of model's shells.hells.
*		The given vertices must be double precision (either 2D or
*		3D) and must correspond to the indices of the vertices in
*		the model.
* \param	model			The geometric model.
* \param	vtxBuf			The buffer with vertex geometry values
*					that are to be set in the model.
*/
static WlzErrorNum WlzGMFilterGeomSetVertices(WlzGMModel *model,
					      WlzVertexP vtxBuf)
{
  int		idx,
  		cnt;
  AlcVector	*vec;
  WlzGMVertex	*cV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  
  cnt = model->res.vertex.numIdx;
  vec = model->res.vertex.vec;
  switch(model->type)
  {
    case WLZ_GMMOD_2I:
      for(idx = 0; idx < cnt; ++idx)
      {
        cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
	cV->geo.vg2I->vtx.vtX = WLZ_NINT((vtxBuf.i2 + idx)->vtX);
	cV->geo.vg2I->vtx.vtY = WLZ_NINT((vtxBuf.i2 + idx)->vtY);
      }
      break;
    case WLZ_GMMOD_2D:
      for(idx = 0; idx < cnt; ++idx)
      {
        cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
	cV->geo.vg2D->vtx = *(vtxBuf.d2 + idx);
      }
      break;
    case WLZ_GMMOD_2N:
      for(idx = 0; idx < cnt; ++idx)
      {
        cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
	cV->geo.vg2N->vtx = *(vtxBuf.d2 + idx);
      }
      break;
    case WLZ_GMMOD_3I:
      for(idx = 0; idx < cnt; ++idx)
      {
        cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
	cV->geo.vg3I->vtx.vtX = WLZ_NINT((vtxBuf.i3 + idx)->vtX);
	cV->geo.vg3I->vtx.vtY = WLZ_NINT((vtxBuf.i3 + idx)->vtY);
	cV->geo.vg3I->vtx.vtZ = WLZ_NINT((vtxBuf.i3 + idx)->vtZ);
      }
      break;
    case WLZ_GMMOD_3D:
      for(idx = 0; idx < cnt; ++idx)
      {
        cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
	cV->geo.vg3D->vtx = *(vtxBuf.d3 + idx);
      }
    case WLZ_GMMOD_3N:
      for(idx = 0; idx < cnt; ++idx)
      {
        cV = (WlzGMVertex *)AlcVectorItemGet(vec, idx);
	cV->geo.vg3N->vtx = *(vtxBuf.d3 + idx);
      }
      break;
    default:
      errNum = WLZ_ERR_GMELM_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMModelRehashVHT(model, 0);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMModelSetSG(model);
  }
  return(errNum);
}
