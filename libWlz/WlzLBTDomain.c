#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzLBTDomain.c
* \author       Bill Hill
* \date         December 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for creating and manipulating linear binary
*		tree domains.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

static int			WlzPartialItv2DCmpFn(
				  const void *vc,
				  const void *v0,
				  const void *v1);
static int			WlzLBTDomain2DNodeCmpFn(
				  const void *lPtr,
				  const void *ptr0,
				  const void *ptr1);
static int			WlzLBTNodeIdxFromKeys2D(
				  WlzLBTDomain2D *lDom,
				  int idN,
				  unsigned *keys,
				  int *dstFound);
static int			WlzLBTBndEdgNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static int			WlzLBTMinLogSzEdgeNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static int			WlzLBTBndEdgNbrDirIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir);
static int			WlzLBTMinLogSzEdgeDirNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir);
static int			WlzLBTMaxLogSzEdgeDirNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir,
				  int *dstSz);
static int			WlzLBTIdxCmpFn(
				  void *datum0,
				  void *datum1);
static int			WlzLBTQueueUnlink(
				  AlcCPQQueue *pQ,
				  AlcHashTable *hT);
static unsigned 		WlzLBTIdxHashFn(
				  void *datum);
static void			WlzLBTCondenseNodes2D(
				  WlzLBTDomain2D *lDom);
static void			WlzLBTSetNodeIndexObj2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static WlzErrorNum 		WlzLBTQueueInsert(
				  AlcCPQQueue *pQ,
				  AlcHashTable *hT,
				  int sz,
				  int idx);

/*!
* \return	New 2D linear binary tree domain.
* \ingroup	WlzAllocation
* \brief	Creates a linear binary tree domain without creating any
*		nodes, leaving the nodes pointer NULL. Only the type,
*		bounding box and depth are set. The depth \f$d\f$ is set
*		such that
*		\f[
		  2^d >= \max(k_l - k_1, l_l - l_1) + 1
		\f]
*		where \f$k_1\f$, \f$k_l\f$, \f$l_1\f$ and \f$l_l\f$
*		are the first column, last column, first line and
*		last line respectively.
* \param	type			Type of domain, which must be
*					WLZ_LBTDOMAIN_2D.
* \param	l1			First line.
* \param	ll			last line.
* \param	k1			First column.
* \param	kl			Last column.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLBTDomain2D	*WlzMakeLBTDomain2D(WlzObjectType type,
				    int l1,
				    int ll,
				    int k1,
				    int kl,
				    WlzErrorNum *dstErr)
{
  unsigned int	sz;
  WlzLBTDomain2D *dom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((type != WLZ_LBTDOMAIN_2D) || (ll < l1) || (kl < k1))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if((dom = (WlzLBTDomain2D *)AlcCalloc(sizeof(WlzLBTDomain2D), 1)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom->type = type;
    dom->line1 = l1;
    dom->lastln = ll;
    dom->kol1 = k1;
    dom->lastkl = kl;
    /* Compute the depth. */
    sz = WLZ_MAX(kl - k1, ll - l1) + 1;
    dom->depth = AlgBitNextPowerOfTwo(NULL, sz);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dom);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Frees the given 2D linear binary tree domain.
* \param	lDom			Given LBT domain.
*/
WlzErrorNum	WlzFreeLBTDomain2D(WlzLBTDomain2D *lDom)
{
  WlzDomain	dom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  dom.l2 = lDom;
  errNum = WlzFreeDomain(dom);
  return(errNum);
}

/*!
* \return	New 2D linear binary tree domain.
* \ingroup	WlzDomainOps
* \brief	Creates a new 2D linear binary tree domain from the given
*		domain.
* \param	dom			Given domain, which must be 2D.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLBTDomain2D	*WlzLBTDomain2DFromDomain(WlzDomain dom, WlzErrorNum *dstErr)
{
  WlzLBTDomain2D *lDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(dom.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(dom.core->type)
    {
      case WLZ_INTERVALDOMAIN_INTVL:
	lDom = WlzLBTDomain2DFromIDomain(dom.i, &errNum);
        break;
      case WLZ_INTERVALDOMAIN_RECT:
	lDom = WlzLBTDomain2DFromIDomain(dom.i, &errNum);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lDom);
}

/*!
* \return	New 2D interval domain.
* \ingroup	WlzDomainOps
* \brief	Creates a new interval domain from the given 2D linear
* 		binary tree domain.
* \param	lDom			Given 2D linear binary tree domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIntervalDomain *WlzLBTDomainToIDomain(WlzLBTDomain2D *lDom,
					 WlzErrorNum *dstErr)
{
  int		idI,
  		idN,
		nItv,
		nSz,
  		nPItv;
  WlzIVertex2	nPos;
  WlzLBTNode2D	*nod;
  WlzPartialItv2D *pItv0,
  		*pItvTb = NULL;
  WlzIntervalDomain *iDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  if(lDom == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(lDom->type != WLZ_LBTDOMAIN_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    /* Count number of node partial intervals. */
    nPItv = 0;
    nod = lDom->nodes;
    for(idN = 0; idN < lDom->nNodes; ++idN)
    {
      nPItv += WlzLBTNodeSz2D(nod++);
    }
    /* Allocate node partial intervals. */
    if((pItvTb = (WlzPartialItv2D *)AlcCalloc(sizeof(WlzPartialItv2D),
					      nPItv)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Set node partial intervals. */
    idI = 0;
    pItv0 = pItvTb;
    nod = lDom->nodes;
    for(idN = 0; idN < lDom->nNodes; ++idN)
    {
      nSz = WlzLBTNodeSz2D(nod); 
      WlzLBTKeyToPos2I(nod->keys, &nPos);
      for(idI = 0; idI < nSz; ++idI)
      {
        pItv0->ileft = nPos.vtX;
	pItv0->iright = nPos.vtX + nSz - 1;
	pItv0->ln = nPos.vtY + idI;
	++pItv0;
      }
      ++nod;
    }
    iDom = WlzIDomainFromPItv2D(lDom->line1, lDom->lastln,
    				lDom->kol1, lDom->lastkl,
				nPItv, pItvTb,
				&errNum);
  }
  /* Free storage. */
  AlcFree(pItvTb);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(iDom);
}

/*!
* \return 	Woolz interval domain or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Allocates and computes an interval domain from the given
*		table of partial intervals. The given partial intervals
*               must fit within the given domain bounding box.
* \param	line1			First line.
* \param	lastln			Last line.
* \param	kol1			First column.
* \param	lastkl			last column.
* \param	nPItv			Number of partial intervals.
* \param	pItv			Array of partial intervals.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzIntervalDomain *WlzIDomainFromPItv2D(int line1, int lastln,
					int kol1, int lastkl,
					int nPItv, WlzPartialItv2D *pItv,
					WlzErrorNum *dstErr)
{
  int		idI,
  		nItv;
  WlzPartialItv2D *pItv0,
  		*pItv1;
  WlzIntervalLine *itvLn;
  WlzInterval	*itv0,
		*itv1,
  		*itvTb = NULL;
  WlzIntervalDomain *iDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((line1 > lastln) || (kol1 > lastkl) || (nPItv < 1) || (pItv == NULL))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(nPItv > 1)
  {
    /* Sort the partial intervals. */
    AlgQSort(pItv, nPItv, sizeof(WlzPartialItv2D), NULL, WlzPartialItv2DCmpFn);
    /* Condense the partial intervals into intervals. */
    idI = 0;
    nItv = 0;
    pItv0 = pItv;
    pItv1 = pItv + 1;
    while(idI < nPItv)
    {
      if((pItv0->ln == pItv1->ln) && ((pItv0->iright + 1) == pItv1->ileft))
      {
        pItv0->iright = pItv1->iright;
      }
      else
      {
	*++pItv0 = *pItv1;
	++nItv;
      }
      ++idI;
      ++pItv1;
    }
    /* Allocate an interval domain, interval lines and intervals. */
    if(((iDom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				      line1, lastln, kol1, lastkl,
				      &errNum)) == NULL) ||
       ((itvTb = (WlzInterval *)AlcCalloc(sizeof(WlzInterval),
       					  nItv)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    iDom->freeptr = AlcFreeStackPush(iDom->freeptr, (void *)itvTb, NULL);
    /* Fill in the interval domain. */
    idI = 0;
    pItv1 = pItv0 = pItv;
    itv1 = itv0 = itvTb;
    while(idI < nItv)
    {
      itvLn = iDom->intvlines + pItv0->ln - iDom->line1;
      while((pItv1->ln == pItv0->ln) && (idI < nItv))
      {
	itv1->ileft = pItv1->ileft;
	itv1->iright = pItv1->iright;
	++(itvLn->nintvs);
	++itv1;
        ++pItv1;
	++idI;
      }
      itvLn->intvs = itv0;
      pItv0 = pItv1;
      itv0 = itv1;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeIntervalDomain(iDom);
    iDom = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(iDom);
}

/*!
* \return	New 2D linear binary tree domain.
* \ingroup	WlzDomainOps
* \brief	Creates a new 2D linear binary tree domain from the given
*		interval domain.
* \todo		There are two possible algorithms here:
*		Create all possible nodes for in Morton order condensing
*		nodes on the fly or the simpler and slower algorithm
*		create all possible nodes in scan order and then
*		sort to get the Morton order. Currently the latter
*		algorithm is used but this should be replaced with the former
*		when there's working code.
*		time than 
*		
* \param	iDom			Given domain, which must be an
*					interval domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLBTDomain2D	*WlzLBTDomain2DFromIDomain(WlzIntervalDomain *iDom,
					WlzErrorNum *dstErr)
{
  int		idI,
  		idN;
  WlzObject	*obj = NULL;
  WlzLBTNode2D	*nod;
  WlzLBTDomain2D *lDom = NULL;
  WlzDomain	dom;
  WlzValues	nullVal;
  WlzIVertex2	pos;
  WlzIntervalWSpace iWsp;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.i = iDom;
  nullVal.core = NULL;
  lDom = WlzMakeLBTDomain2D(WLZ_LBTDOMAIN_2D,
  			    iDom->line1, iDom->lastln,
			    iDom->kol1, iDom->lastkl, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom, nullVal, NULL, NULL, &errNum);
  }
  /* Compute area for the initial (and maximum) number of nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitRasterScan(obj, &iWsp, WLZ_RASTERDIR_ILIC);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    while((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE )
    {
      lDom->maxNodes += iWsp.rgtpos - iWsp.lftpos + 1;
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  /* Allocate nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((lDom->nodes = (WlzLBTNode2D *)AlcCalloc(sizeof(WlzLBTNode2D),
    						lDom->maxNodes)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      lDom->nNodes = lDom->maxNodes;
      lDom->freeptr = AlcFreeStackPush(lDom->freeptr, lDom->nodes, &alcErr);
      if(alcErr != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  /* Set keys of all initial nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitRasterScan(obj, &iWsp, WLZ_RASTERDIR_ILIC);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nod = lDom->nodes;
    while((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE )
    {
      pos.vtY = iWsp.linpos;
      pos.vtX = iWsp.lftpos;
      for(idI = 0; idI < iWsp.colrmn; ++idI)
      {
        WlzLBTPosToKey2D(pos, nod->keys);
	++nod;
	++(pos.vtX);
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  (void )WlzFreeObj(obj);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Condense the nodes to get the linear binary tree. */
    WlzLBTCondenseNodes2D(lDom);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lDom);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Balances the given LBT domain so that the neighbouring
*		nodes of each node are either of the same size or differ
*		in size by a ratio of 2:1. The function also enforces
*		maximum node size for all nodes and boundary nodes.
*		The neighbour finding algorithm used is quick and
*		simple but it requires an object in which the values
*		are set to the corresponding LBT domain indices.
*		For efficiency an associated interval domain may be given,
*		if the associated domain pointer is NULL then a domain
*		will be computed.
* \param	lDom			Given LBT domain.
* \param	iObj			Index object for finding neighbours
*					of nodes.
* \param	maxSz			Maximum node size.
*/
WlzErrorNum	WlzLBTBalanceDomain2D(WlzLBTDomain2D *lDom,
				      WlzObject *iObj,
				      int maxSz,
				      int maxBndSz)
{
  int		idN,
  		idM,
		idP,
		flg,
		sz0,
		sz1,
		sz2,
		sInc;
  int		idNN[4];
  WlzLBTNode2D	*nod[4];
  AlcCPQQueue	*pQ = NULL;
  AlcHashTable	*hT = NULL;
  WlzGreyValueWSpace *iGVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzDirection dirTab[4] = {
				   WLZ_DIRECTION_IC,
				   WLZ_DIRECTION_IL,
				   WLZ_DIRECTION_DC,
				   WLZ_DIRECTION_DL
  				 };

  maxSz = (maxSz < 1)? 1: AlgBitNextPowerOfTwo(NULL, maxSz) + 1;
  maxBndSz = (maxBndSz < 1)? 1: AlgBitNextPowerOfTwo(NULL, maxBndSz) + 1;
  iGVWSp = WlzGreyValueMakeWSp(iObj, &errNum);
  /* Create a priority queue and a hash table to hold the nodes
   * to be split. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((pQ = AlcCPQQueueNew(NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sz0 = (lDom->nNodes / 8);
    if((sz0 = (lDom->nNodes / 8)) < 1024)
    {
      sz0 = 1024;
    }
    if((hT = AlcHashTableNew(sz0, WlzLBTIdxCmpFn,
    			     WlzLBTIdxHashFn, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Put all nodes which need to be split into a priority queue using
   * the size as the priority. Nodes need to be split if they have a
   * size (cell side length) > twice that of the smallest neighbour,
   * a size greater than the maximum size or a size greater than the
   * maximum size for a boundary node. */
  idN = 0;
  while((errNum == WLZ_ERR_NONE) && (idN < lDom->nNodes))
  {
    flg = 0;
    nod[0] = lDom->nodes + idN;
    sz0 = WlzLBTNodeLogSz2D(nod[0]);
    if(sz0 > maxSz)
    {
      flg = 1;
    }
    else if(sz0 > 1)
    {
      sz1 = WlzLBTMinLogSzEdgeNbrIdx2D(lDom, iGVWSp, idN);
      if((sz1 >= 0) && (sz0 > sz1 + 1))
      {
        flg = 1;
      }
    }
    if(0 == flg)
    {
      flg = (sz0 > maxBndSz) &&
	    (WlzLBTBndEdgNbrIdx2D(lDom, iGVWSp, idN) != 0);
    }
    if(flg)
    {
      errNum = WlzLBTQueueInsert(pQ, hT, sz0, idN);
    }
    ++idN;
  }
  /* Pull the nodes from the priority queue and split them, returning
   * child nodes which have a size > the maximum node size or > twice
   * that of the smallest neighbour to the priority queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    while((idN = WlzLBTQueueUnlink(pQ, hT)) >= 0)
    {
      /* Split the node. */
      idNN[0] = idN;
      idNN[1] = lDom->nNodes;
      idNN[2] = lDom->nNodes + 1;
      idNN[3] = lDom->nNodes + 2;
      lDom->nNodes += 3;
      if(lDom->nNodes > lDom->maxNodes)
      {
        /* Reallocate the nodes. */
	sz0 = (lDom->maxNodes > 0)?
	      ((int )floor(lDom->nNodes / lDom->maxNodes) + 1) *
	      lDom->maxNodes: 1024;
	if((lDom->nodes = (WlzLBTNode2D *)
			  AlcRealloc(lDom->nodes,
			  	     sz0 * sizeof(WlzLBTNode2D))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  lDom->maxNodes = sz0;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	nod[0] = lDom->nodes + idNN[0];
	nod[1] = lDom->nodes + idNN[1];
	nod[2] = lDom->nodes + idNN[2];
	nod[3] = lDom->nodes + idNN[3];
	sInc = WlzLBTNodeSz2D(nod[0]) / 2;
	nod[1]->keys[0] = nod[0]->keys[0] + sInc;
	nod[1]->keys[1] = nod[0]->keys[1];
	nod[2]->keys[0] = nod[0]->keys[0];
	nod[2]->keys[1] = nod[0]->keys[1] + sInc;
	nod[3]->keys[0] = nod[0]->keys[0] + sInc;
	nod[3]->keys[1] = nod[0]->keys[1] + sInc;
	nod[0]->keys[2] >>= 1;
	nod[3]->keys[2] = nod[2]->keys[2] = nod[1]->keys[2] = nod[0]->keys[2];
	/* Set indices in the index object. */
	for(idN = 0; idN < 4; ++idN)
	{
	  WlzLBTSetNodeIndexObj2D(lDom, iGVWSp, idNN[idN]);
	}
	/* Test each child node to see if either the child node or it's
	 * neighbours need to be inserted into the queue. */
	idN = 0;
	sz0 = WlzLBTNodeLogSz2D(nod[0]);
	while((errNum == WLZ_ERR_NONE) && (idN < 4))
	{
	  /* Check nodes and reinsert any that need to be split back into
	   * the queue. */
	  flg = 0;
	  if(sz0 > maxSz)
	  {
	    flg = 1;
	  }
	  else if(sz0 > 1)
	  {
	    sz1 = WlzLBTMinLogSzEdgeNbrIdx2D(lDom, iGVWSp, idNN[idN]);
	    if((sz1 >= 0) && (sz0 > sz1 + 1))
	    {
	      flg = 1;
	    }
	  }
	  if(0 == flg)
	  {
	    flg = (sz0 > maxBndSz) &&
		  (WlzLBTBndEdgNbrIdx2D(lDom, iGVWSp, idNN[idN]) != 0);
	  }
	  if(flg)
	  {
	    errNum = WlzLBTQueueInsert(pQ, hT, sz0, idNN[idN]);
	  }
	  /* Check for any neighbours with node size > twice that of
	   * the node. */
	  idM = 0;
	  while((errNum == WLZ_ERR_NONE) && (idM < 4))
	  {
	    idP = WlzLBTMaxLogSzEdgeDirNbrIdx2D(lDom, iGVWSp, idNN[idN],
						dirTab[idM], &sz2);
	    if(sz2 > sz0 + 1)
	    {
	      errNum = WlzLBTQueueInsert(pQ, hT, sz2, idP);
	    }
	    ++idM;
	  }
	  ++idN;
	}
      }
    }
  }
  /* Free temporary storage. */
  WlzGreyValueFreeWSp(iGVWSp);
  (void )AlcCPQQueueFree(pQ);
  (void )AlcHashTableFree(hT);
  /* Sort the nodes so that all the split nodes are in the appropriate
   * place. */
  if(errNum == WLZ_ERR_NONE)
  {
    AlgQSort(lDom->nodes, lDom->nNodes, sizeof(WlzLBTNode2D), lDom,
             WlzLBTDomain2DNodeCmpFn);
  }
  return(errNum);
}

/*!
* \return	New object with node indices or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Creates a new 2D domain object with integer values 
*		which are the indices of the nodes containing the
*		pixels of the domain.
* \param	lDom			Given LBT domain.
* \param	iDom			Corresponding interval domain,
					may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzLBTMakeNodeIndexObj2D(WlzLBTDomain2D *lDom,
				       WlzIntervalDomain *iDom,
				       WlzErrorNum *dstErr)
{
  int 		idN;
  WlzDomain	tDom;
  WlzValues	tVal;
  WlzPixelV	iBkg;
  WlzObjectType iValTblType;
  WlzGreyValueWSpace *iGVWSp = NULL;
  WlzObject	*iObj = NULL;
  WlzIBox2	nBox;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = -1;
  tDom.core = NULL;
  tVal.core = NULL;
  if(iDom)
  {
    switch(iDom->type)
    {
      case WLZ_INTERVALDOMAIN_INTVL:
      case WLZ_INTERVALDOMAIN_RECT:
	tDom.i = iDom;
        break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  else
  {
    tDom.i = WlzLBTDomainToIDomain(lDom, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    iObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, tDom, tVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    iValTblType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
    iBkg.type = WLZ_GREY_INT;
    iBkg.v.ubv = 0;
    tVal.v = WlzNewValueTb(iObj, iValTblType, iBkg, &errNum);
    iObj->values = WlzAssignValues(tVal, NULL);
    iGVWSp = WlzGreyValueMakeWSp(iObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < lDom->nNodes; ++idN)
    {
      WlzLBTSetNodeIndexObj2D(lDom, iGVWSp, idN);
    }
    (void )WlzSetBackground(iObj, bgdV);
  }
  WlzGreyValueFreeWSp(iGVWSp);
  if(errNum != WLZ_ERR_NONE)
  {
    if(iObj)
    {
      (void )WlzFreeObj(iObj);
      iObj = NULL;
    }
    else
    {
      if((iDom == NULL) && (tDom.core != NULL))
      {
	(void )WlzFreeIntervalDomain(tDom.i);
      }
      if(tVal.core)
      {
        (void )WlzFreeValueTb(tVal.v);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(iObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Sets all index values in an existing LBT node
*		index object.
* \param	lDom
* \param	iObj
*/
WlzErrorNum	WlzLBTIndexObjSetAllNodes2D(WlzLBTDomain2D *lDom,
				        WlzObject *iObj)
{
  int		idN;
  WlzPixelV	bgdV;
  WlzGreyValueWSpace *iGVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.v.inv = -1;
  iGVWSp = WlzGreyValueMakeWSp(iObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < lDom->nNodes; ++idN)
    {
      WlzLBTSetNodeIndexObj2D(lDom, iGVWSp, idN);
    }
    (void )WlzSetBackground(iObj, bgdV);
  }
  WlzGreyValueFreeWSp(iGVWSp);
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Classifies the given LBT node by it's connectivity and
*		returns it's class and the counter-clockwise rotation
*		of the basic class pattern in multiples of 90 degrees.
* \param	lDom			Linear binary tree domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of the LBT node.
* \param	dstCls			Destination pointer for the class.
* \param	dstRot			Destination pointer for the rotation.
*/
void		WlzLBTClassifyNode2D(WlzLBTDomain2D *lDom,
				     WlzGreyValueWSpace *iGVWSp,
				     int idN, WlzLBTNodeClass2D *dstCls,
				     int *dstRot)
{
  int		idM,
		nNbr,
		nSz;
  unsigned	msk;
  WlzIBox2	nBB;
  const WlzDirection dirTab[4] = {
				  WLZ_DIRECTION_IC,
				  WLZ_DIRECTION_IL,
				  WLZ_DIRECTION_DC,
				  WLZ_DIRECTION_DL
			 	};
  const int	 rotTab[16] = {
				0, 		/*  0  o-o
						       | |
						       o-o */
				3,		/*  1  o-o
						       | o
						       o-o */
				0,		/*  2  ooo
						       | |
						       o-o */
				0,		/*  3  ooo
						       | o
						       o-o */
				1,		/*  4  o-o
						       o |
						       o-o */
				3,		/*  5  o-o
						       o o
						       o-o */
				1,		/*  6  ooo
						       o |
						       o-o */
				1,		/*  7  ooo
						       o o
						       o-o */
				2,		/*  8  o-o
						       | |
						       ooo */
				3,		/*  9  o-o
						       | o
						       ooo */
				0,		/* 10  ooo
						       | |
						       ooo */
				0,		/* 11  ooo
						       | o
						       ooo */
				2,		/* 12  o-o
						       o |
						       ooo */
				3,		/* 13  o-o
						       o o
						       ooo */
				2,		/* 14  ooo
						       o |
						       ooo */
				0};		/* 15  ooo
						       o o
						       ooo */
  const WlzLBTNodeClass2D clsTab[16] = {
				0, 		/*  0  o-o
						       | |
						       o-o */
				1,		/*  1  o-o
						       | o
						       o-o */
				1,		/*  2  ooo
						       | |
						       o-o */
				2,		/*  3  ooo
						       | o
						       o-o */
				1,		/*  4  o-o
						       o |
						       o-o */
				3,		/*  5  o-o
						       o o
						       o-o */
				2,		/*  6  ooo
						       o |
						       o-o */
				4,		/*  7  ooo
						       o o
						       o-o */
				1,		/*  8  o-o
						       | |
						       ooo */
				2,		/*  9  o-o
						       | o
						       ooo */
				3,		/* 10  ooo
						       | |
						       ooo */
				4,		/* 11  ooo
						       | o
						       ooo */
				2,		/* 12  o-o
						       o |
						       ooo */
				4,		/* 13  o-o
						       o o
						       ooo */
				4,		/* 14  ooo
						       o |
						       ooo */
				5};		/* 15  ooo
						       o o
						       ooo */


  nSz = WlzLBTNodeSz2D(lDom->nodes + idN);
  if(nSz == 1)
  {
    *dstRot = 0;
    *dstCls = 0;
  }
  else
  {
    WlzLBTKeyToBox2I((lDom->nodes + idN)->keys, &nBB);
    /* WLZ_DIRECTION_IC */
    WlzGreyValueGet(iGVWSp, 0, nBB.yMin, nBB.xMax + 1);
    idM = iGVWSp->gVal[0].inv;
    if(idM < 0)
    {
      WlzGreyValueGet(iGVWSp, 0, nBB.yMax, nBB.xMax + 1);
      idM = iGVWSp->gVal[0].inv;
    }
    nNbr = (idM >= 0) && (WlzLBTNodeSz2D(lDom->nodes + idM) < nSz);
    msk = nNbr;
    /* WLZ_DIRECTION_IL */
    WlzGreyValueGet(iGVWSp, 0, nBB.yMax + 1, nBB.xMin);
    idM = iGVWSp->gVal[0].inv;
    if(idM < 0)
    {
      WlzGreyValueGet(iGVWSp, 0, nBB.yMax + 1, nBB.xMax);
      idM = iGVWSp->gVal[0].inv;
    }
    nNbr = (idM >= 0) && (WlzLBTNodeSz2D(lDom->nodes + idM) < nSz);
    msk |= nNbr << 1;
    /* WLZ_DIRECTION_DC */
    WlzGreyValueGet(iGVWSp, 0, nBB.yMin, nBB.xMin - 1);
    idM = iGVWSp->gVal[0].inv;
    if(idM < 0)
    {
      WlzGreyValueGet(iGVWSp, 0, nBB.yMax, nBB.xMin - 1);
      idM = iGVWSp->gVal[0].inv;
    }
    nNbr = (idM >= 0) && (WlzLBTNodeSz2D(lDom->nodes + idM) < nSz);
    msk |= nNbr << 2;
    /* WLZ_DIRECTION_DL */
    WlzGreyValueGet(iGVWSp, 0, nBB.yMin - 1, nBB.xMin);
    idM = iGVWSp->gVal[0].inv;
    if(idM < 0)
    {
      WlzGreyValueGet(iGVWSp, 0, nBB.yMin - 1, nBB.xMax);
      idM = iGVWSp->gVal[0].inv;
    }
    nNbr = (idM >= 0) && (WlzLBTNodeSz2D(lDom->nodes + idM) < nSz);
    msk |= nNbr << 3;
    /* Look up the class and rotation. */
    *dstCls = clsTab[msk];
    *dstRot = rotTab[msk];
  }
}

/*!
* \return	Returns signed comparison for AlgQSort().
* \ingroup	WlzDomainOps
* \brief	Compares to partial intervals.
* \param	vc			Unused client data.
* \param	v0			First partial interval.
* \param	v1			Second partial interval.
*/
static int	WlzPartialItv2DCmpFn(const void *vc,
				     const void *v0, const void *v1)
{
  int		cmp;
  WlzPartialItv2D *pItv0,
  		*pItv1;

  pItv0 = (WlzPartialItv2D *)v0;
  pItv1 = (WlzPartialItv2D *)v1;
  if((cmp = pItv0->ln - pItv1->ln) == 0)
  {
    cmp = pItv0->ileft - pItv1->ileft;
  }
  return(cmp);
}

/*!
* \return	result of comparison.
* \ingroup	WlzDomainOps
* \brief	Compares to int's for the hash table.
* \param	datum0			Pointer to first int.
* \param	datum1			Pointer to second int.
*/
static int	WlzLBTIdxCmpFn(void *datum0, void *datum1)
{
  int		cmp;

  cmp = (int )datum0 - (int )datum1;
  return(cmp);
}

/*!
* \return	Hash value.
* \ingroup	WlzDomainOps
* \brief	Simple hash function for hashing node index.
* \param	datum
*/
static unsigned WlzLBTIdxHashFn(void *datum)
{
  unsigned	hVal;
  const unsigned p0 = 399989, /* These are three different 6 digit primes */
  		 p1 = 599999,
  		 p2 = 999983;
  hVal = (((unsigned )datum + p0) * p1) % p2;
  return(hVal);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Inserts a node's index into the pripority queue using the
*		size as the priority. The hash table is used to ensure that
*		nodes are only in the priority queue once.
* \param	pQ			Priority queue.
* \param	hT			Hash table.
* \param	sz			Log node size.
* \param	idx			Node index.
*/
static WlzErrorNum WlzLBTQueueInsert(AlcCPQQueue *pQ, AlcHashTable *hT,
				     int sz, int idx)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(AlcHashItemGet(hT, (void *)idx, NULL) == NULL)
  {
    if((AlcCPQEntryInsert(pQ, (float )sz, (void *)idx) != ALC_ER_NONE) ||
       (AlcHashTableEntryInsert(hT, (void *)idx, (void *)idx,
                                NULL) != ALC_ER_NONE))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  return(errNum);
}

/*!
* \return	Index of top priority node or -1 if queue empty.
* \ingroup	WlzDomainOps
* \brief	Unlinks the top priority queue entry and removes it from
*		the hash table.
* \param	pQ			Priority queue.
* \param	hT			Hash table.
*/
static int	WlzLBTQueueUnlink(AlcCPQQueue *pQ, AlcHashTable *hT)
{
  int		idx = -1;
  AlcHashItem	*hItem;
  AlcCPQItem	*pItem;

  if((pItem = AlcCPQItemUnlink(pQ)) != NULL)
  {
    idx = (int )(pItem->entry);
    AlcCPQItemFree(pQ, pItem);
    if((hItem = AlcHashItemGet(hT, (void *)(idx), NULL)) != NULL)
    {
      (void )AlcHashItemUnlink(hT, hItem, 1);
    }
  }
  return(idx);
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets values within a nodes bounding box in an index
*		object. All parameters are assumed valid.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
*/
static void	WlzLBTSetNodeIndexObj2D(WlzLBTDomain2D *lDom,
				      WlzGreyValueWSpace *iGVWSp, int idN)
{
  int		pX,
  		pY;
  WlzIBox2	nBox;
  WlzLBTNode2D	*nod;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox2I(nod->keys, &nBox);
  for(pY = nBox.yMin; pY <= nBox.yMax; ++pY)
  {
    for(pX = nBox.xMin; pX <= nBox.xMax; ++pX)
    {
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      *(iGVWSp->gPtr[0].inp) = idN;
    }
  }
}

/*!
* \return	Non-zero if the node is on the boundary of the domain.
* \ingroup	WlzDomainOps
* \brief	Looks for a non-existant neighbour to the given node.
* \param        lDom                    Given LBT domain.
* \param        iGVWSp                  Grey workspace for index object.
* \param        idN                     Index of node in the LBT domain.
*/
static int	WlzLBTBndEdgNbrIdx2D(WlzLBTDomain2D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN)
{
  int		isBnd;

  isBnd = WlzLBTBndEdgNbrDirIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_IC);
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_IL);
  }
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_DC);
  }
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_DL);
  }
  return(isBnd);
}

/*!
* \return	Log of the size of the smallest neighbouring node or
		-ve if there is no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds the size of the smallest neighbouring node.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
*/
static int	WlzLBTMinLogSzEdgeNbrIdx2D(WlzLBTDomain2D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN)
{
  int		sz,
  		minSz;

  minSz = WlzLBTMinLogSzEdgeDirNbrIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_IC);
  sz = WlzLBTMinLogSzEdgeDirNbrIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_IL);
  if((minSz < 0) || ((sz >= 0) && (sz < minSz)))
  {
    minSz = sz;
  }
  sz = WlzLBTMinLogSzEdgeDirNbrIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_DC);
  if((minSz < 0) || ((sz >= 0) && (sz < minSz)))
  {
    minSz = sz;
  }
  sz = WlzLBTMinLogSzEdgeDirNbrIdx2D(lDom, iGVWSp, idN, WLZ_DIRECTION_DL);
  if((minSz < 0) || ((sz >= 0) && (sz < minSz)))
  {
    minSz = sz;
  }
  return(minSz);
}

/*!
* \return	Non-zero if the node's neighbour in the given diection is
*		outside the domain.
* \ingroup	WlzDomainOps
* \brief	Looks for a non-existant neighbour to the given node
*		in the given direction.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
* \param	dir			Given direction.
*/
static int	WlzLBTBndEdgNbrDirIdx2D(WlzLBTDomain2D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN, WlzDirection dir)
{
  int		pX,
  		pY,
		isBnd = 0;
  WlzLBTNode2D	*nod;
  WlzIBox2	nBB;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox2I(nod->keys, &nBB);
  switch(dir)
  {
    case WLZ_DIRECTION_IC: /* FALLTHROUGH */
    case WLZ_DIRECTION_DC:
      pX = (dir == WLZ_DIRECTION_IC)? nBB.xMax + 1: nBB.xMin - 1;
      pY = nBB.yMin;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      isBnd = iGVWSp->gVal[0].inv < 0;
      while(!isBnd && (++pY <= nBB.yMax))
      {
        WlzGreyValueGet(iGVWSp, 0, pY, pX);
	isBnd = iGVWSp->gVal[0].inv < 0;
      }
      break;
    case WLZ_DIRECTION_IL: /* FALLTHROUGH */
    case WLZ_DIRECTION_DL:
      pX = nBB.xMin;
      pY = (dir == WLZ_DIRECTION_IL)? nBB.yMax + 1: nBB.yMin - 1;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      isBnd = iGVWSp->gVal[0].inv < 0;
      while(!isBnd && (++pX <= nBB.xMax))
      {
        WlzGreyValueGet(iGVWSp, 0, pY, pX);
	isBnd = iGVWSp->gVal[0].inv < 0;
      }
      break;
  }
  return(isBnd);
}

/*!
* \return	Log of the size of the smallest neighbouring node or
		-ve if there is no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds the size of the smallest neighbouring node in the
*		given direction.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
* \param	dir			Given direction.
*/
static int	WlzLBTMinLogSzEdgeDirNbrIdx2D(WlzLBTDomain2D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN, WlzDirection dir)
{
  int		id0,
  		id1,
		pX,
  		pY,
		sz,
		minSz = -1;
  WlzLBTNode2D	*nod;
  WlzIBox2	nBB;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox2I(nod->keys, &nBB);
  switch(dir)
  {
    case WLZ_DIRECTION_IC: /* FALLTHROUGH */
    case WLZ_DIRECTION_DC:
      pX = (dir == WLZ_DIRECTION_IC)? nBB.xMax + 1: nBB.xMin - 1;
      pY = nBB.yMin;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      id1 = iGVWSp->gVal[0].inv;
      minSz = sz = (id1 >= 0)? WlzLBTNodeLogSz2D(lDom->nodes + id1): -1;
      while(++pY <= nBB.yMax)
      {
	id0 = id1;
        WlzGreyValueGet(iGVWSp, 0, pY, pX);
	id1 = iGVWSp->gVal[0].inv;
	if((id1 != id0) && (id1 >= 0))
	{
	  sz = WlzLBTNodeLogSz2D(lDom->nodes + id1);
	  if((sz >= 0) && ((minSz < 0) || (sz < minSz)))
	  {
	    minSz = sz;
	  }
	}
      }
      break;
    case WLZ_DIRECTION_IL: /* FALLTHROUGH */
    case WLZ_DIRECTION_DL:
      pX = nBB.xMin;
      pY = (dir == WLZ_DIRECTION_IL)? nBB.yMax + 1: nBB.yMin - 1;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      id1 = iGVWSp->gVal[0].inv;
      minSz = sz = (id1 >= 0)? WlzLBTNodeLogSz2D(lDom->nodes + id1): -1;
      while(++pX <= nBB.xMax)
      {
	id0 = id1;
        WlzGreyValueGet(iGVWSp, 0, pY, pX);
	id1 = iGVWSp->gVal[0].inv;
	if((id1 != id0) && (id1 >= 0))
	{
	  sz = WlzLBTNodeLogSz2D(lDom->nodes + id1);
	  if((sz >= 0) && ((minSz < 0) || (sz < minSz)))
	  {
	    minSz = sz;
	  }
	}
      }
      break;
  }
  return(minSz);
}

/*!
* \return	Number of neighbours in given direction.
* \ingroup	WlzDomainOps
* \brief	Counts the number of neighbours of the given node in
*		the given direction.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
* \param	dir			Given direction.
*/
int		WlzLBTCountNodNbrDir2D(WlzLBTDomain2D *lDom,
				       WlzGreyValueWSpace *iGVWSp,
				       int idN, WlzDirection dir)
{
  int		id0,
  		id1,
		pX,
  		pY,
		nCnt = 0;
  WlzLBTNode2D	*nod;
  WlzIBox2	nBB;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox2I(nod->keys, &nBB);
  switch(dir)
  {
    case WLZ_DIRECTION_IC: /* FALLTHROUGH */
    case WLZ_DIRECTION_DC:
      pX = (dir == WLZ_DIRECTION_IC)? nBB.xMax + 1: nBB.xMin - 1;
      pY = nBB.yMin;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      id0 = iGVWSp->gVal[0].inv;
      nCnt = id0 >= 0;
      while(++pY <= nBB.yMax)
      {
	WlzGreyValueGet(iGVWSp, 0, pY, pX);
        id1 = iGVWSp->gVal[0].inv;
	nCnt += id0 != id1;
	id0 = id1;
      }
      break;
    case WLZ_DIRECTION_IL: /* FALLTHROUGH */
    case WLZ_DIRECTION_DL:
      pX = nBB.xMin;
      pY = (dir == WLZ_DIRECTION_IL)? nBB.yMax + 1: nBB.yMin - 1;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      id0 = iGVWSp->gVal[0].inv;
      nCnt = id0 >= 0;
      while(++pX <= nBB.xMax)
      {
	WlzGreyValueGet(iGVWSp, 0, pY, pX);
        id1 = iGVWSp->gVal[0].inv;
	nCnt += id0 != id1;
	id0 = id1;
      }
      break;
  }
  return(nCnt);
}

/*!
* \return	Index of maximum sized neighbour or -ve if no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds a maximum sized neighbour in the given direction
*		then returns it's index and size.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
* \param	dir			Given direction.
* \param	dstSz			Destination pointer for log size
*					of neighbour.
*/
static int	WlzLBTMaxLogSzEdgeDirNbrIdx2D(WlzLBTDomain2D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN, WlzDirection dir,
					int *dstSz)
{
  int		id0,
  		id1,
		pX,
  		pY,
		sz,
		idM = -1,
  		szM = -1;
  WlzLBTNode2D	*nod;
  WlzIBox2	nBB;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox2I(nod->keys, &nBB);
  switch(dir)
  {
    case WLZ_DIRECTION_IC: /* FALLTHROUGH */
    case WLZ_DIRECTION_DC:
      pX = (dir == WLZ_DIRECTION_IC)? nBB.xMax + 1: nBB.xMin - 1;
      pY = nBB.yMin;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      idM = id1 = iGVWSp->gVal[0].inv;
      szM = sz = (id1 >= 0)? WlzLBTNodeLogSz2D(lDom->nodes + id1): -1;
      while(++pY <= nBB.yMax)
      {
	id0 = id1;
        WlzGreyValueGet(iGVWSp, 0, pY, pX);
	id1 = iGVWSp->gVal[0].inv;
	if((id1 != id0) && (id1 >= 0))
	{
	  sz = WlzLBTNodeLogSz2D(lDom->nodes + id1);
	  if((sz >= 0) && ((szM < 0) || (sz > szM)))
	  {
	    szM = sz;
	    idM = id1;
	  }
	}
      }
      break;
    case WLZ_DIRECTION_IL: /* FALLTHROUGH */
    case WLZ_DIRECTION_DL:
      pX = nBB.xMin;
      pY = (dir == WLZ_DIRECTION_IL)? nBB.yMax + 1: nBB.yMin - 1;
      WlzGreyValueGet(iGVWSp, 0, pY, pX);
      idM = id1 = iGVWSp->gVal[0].inv;
      szM = sz = (id1 >= 0)? WlzLBTNodeLogSz2D(lDom->nodes + id1): -1;
      while(++pX <= nBB.xMax)
      {
	id0 = id1;
        WlzGreyValueGet(iGVWSp, 0, pY, pX);
	id1 = iGVWSp->gVal[0].inv;
	if((id1 != id0) && (id1 >= 0))
	{
	  sz = WlzLBTNodeLogSz2D(lDom->nodes + id1);
	  if((sz >= 0) && ((szM < 0) || (sz > szM)))
	  {
	    szM = sz;
	    idM = id1;
	  }
	}
      }
      break;
  }
  *dstSz = szM;
  return(idM);
}

/*!
* \return	Node size: 0 on error else 1, 2, 4, ....
* \ingroup	WlzDomainOps
* \brief	Computes the size of a 2D linear binary tree node.
* \param	nod			Given node.
*/
int		WlzLBTNodeSz2D(WlzLBTNode2D *nod)
{
  int		key,
  		sz = 0;

  if(nod)
  {
    key = nod->keys[2];
    sz = 1;
    while(key)
    {
      key >>= 1;
      sz <<= 1;
    }
  }
  return(sz);
}

/*!
* \return	Node size: -1 on error else 0, 1, 2, ....
* \ingroup	WlzDomainOps
* \brief	Computes the log (base 2) of the size of a 2D linear binary
*		tree node.
* \param	nod			Given node.
*/
int		WlzLBTNodeLogSz2D(WlzLBTNode2D *nod)
{
  int		key,
  		sz = -1;

  if(nod)
  {
    key = nod->keys[2];
    sz = 0;
    while(key)
    {
      key >>= 1;
      ++sz;
    }
  }
  return(sz);
}

/*!
* \return
* \ingroup
* \brief	Finds the index of the node which matches the given keys,
*		starting the search from the node with index \f$n\f$.
*		The index of the \f$n\f$ node is used for efficiency and
*		should be as close as possible to the node that will be
*		matched, but can be any integer value.
*		In many cases there will not be a node in the domain which
*		matches the given keys and in these cases the matched node
*		should share as many higher digits as possible with the
*		given key. For example if the given keys are {4, 0, 0}
*		(digits 0100) and the minimum and maximum keys found are
*		{3, 3, 0} (digits 0033) and {5, 2, 0} (digits 0121)
*		respectively the matched node would have keys {5, 2, 0}.
* \param	lDom			Given 2D linear binary tree domain.
* \param	idN			Index of the \f$n\f$'th node.
* \param	gKeys			Direction \f$d\f$.
* \param	dstFound		Destination pointer for flag which is
*					set non zero if an exact match is
*					found. May be NULL.
*/
int		WlzLBTNodeIdxFromKeys2D(WlzLBTDomain2D *lDom, int idN,
				      unsigned *gKeys, int *dstFound)
{
  int		idMin,
  		idMax,
		cmp,
		found = 0;
  unsigned	k0,
  		k1,
		dG,
		dMax,
		dMin,
		dMsk;
  unsigned	*cKeys,
  		*minKeys,
		*maxKeys;

  idMin = 0;
  idMax = lDom->nNodes;
  cKeys = (lDom->nodes + idN)->keys;
  while((found == 0) && ((idMax - idMin) > 1))
  {
    cmp = WlzLBTDomain2DNodeCmpFn(lDom, gKeys, cKeys);
    if(cmp < 0)
    {
      idMax = idN;
      idN = (idN + idMin) / 2;
      found = idMin == idMax;
    }
    else if (cmp > 0)
    {
      idMin = idN;
      idN = (idN + idMax) / 2;
      found = idMin == idMax;
    }
    else
    {
      found = 1;
    }
  }
  if(!found)
  {
    /* No node has been found with an exact match of the given keys so
     * need to decide which of the nodes lies in the same cell as the
     * node withe the given keys. */
    dMsk = 1 << lDom->depth;
    minKeys = (lDom->nodes + idMin)->keys;
    maxKeys = (lDom->nodes + idMax)->keys;
    do
    {

      k0 = (gKeys[0] & dMsk) != 0; 
      k1 = (gKeys[1] & dMsk) != 0; 
      dG = (k1 << 1) | k0;
      k0 = (minKeys[0] & dMsk) != 0; 
      k1 = (minKeys[1] & dMsk) != 0; 
      dMin = (k1 << 1) | k0;
      k0 = (maxKeys[0] & dMsk) != 0; 
      k1 = (maxKeys[1] & dMsk) != 0; 
      dMax = (k1 << 1) | k0;
      idN = (dG == dMin)? idMin: idMax;
      dMsk >>= 1;
    } while(dMsk && (dMin == dMax));
  }
  if(dstFound)
  {
    *dstFound = found;
  }
  return(idN);
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets the value of the LBT key for the given position,
*		where the position is relative to the first line and
*		column of the domain.
*		The key is encoded by interleaving the bits of the
*		column coordinate, line coordinate and term value
*		in that order.
*		The keys are ordered so that bit 0 is least significant.
* \param	pos			Position.
* \param	keys			Keys to be set in line, column and
*					term order.
*/
void		WlzLBTPosToKey2D(WlzIVertex2 pos, unsigned *keys)
{
  keys[0] = pos.vtX;
  keys[1] = pos.vtY;
  keys[2] = 0;
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Gets an array of WLZ_LBTDOMAIN_MAXDIGITS location digits
*		from the given keys.
* \param	keys			Keys to be read.
* \param	digits			Array of WLZ_LBTDOMAIN_MAXDIGITS
*					digits.
*/
void		WlzLBTGetKeyDigits2D(unsigned *keys, UBYTE *digits)
{
  int		idD;
  unsigned	k0,
  		k1,
		k2,
		dMsk;

  idD = 0;
  dMsk = 1;
  while(idD < WLZ_LBTDOMAIN_MAXDIGITS)
  {
    k0 = (keys[0] & dMsk) != 0; 
    k1 = (keys[1] & dMsk) != 0; 
    k2 = (keys[2] & dMsk) != 0; 
    digits[idD] = k0 + (k1 << 1) + (k2 << 2);
    ++idD;
    dMsk <<= 1;
  }
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets position (relative to the first line and column of the
* 		domain) which corresponds to the given LBT key.
* \param	key			Keys in line, column and term order.
* \param	pos			Position to be set.
*/
void		WlzLBTKeyToPos2I(unsigned *key, WlzIVertex2 *pos)
{
  pos->vtX = key[0];
  pos->vtY = key[1];
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets bounding box which corresponds to the given LBT key.
* \param	key			Keys in line, column and term order.
* \param	box			Position to be set.
*/
void		WlzLBTKeyToBox2I(unsigned *key, WlzIBox2 *box)
{
  unsigned int	k,
  		sz = 1;

  box->xMin = key[0];
  box->yMin = key[1];
  if(key[2] != 0)
  {
    k = key[2];
    while(k)
    {
      sz <<= 1;
      k >>= 1;
    }
  }
  box->xMax = box->xMin + sz - 1;
  box->yMax = box->yMin + sz - 1;
}

/*!
* \return	Sort value - less than, equal to, or greater than zero.
* \ingroup	WlzDomainOps
* \brief	Sorts the keys into accending order by node value.
* \param	ptrC			Used to pass the domain pointer.
* \param	ptr0			Used to pass first node.
* \param	ptr1			Used to pass second node.
*/
static int	WlzLBTDomain2DNodeCmpFn(const void *ptrC,
				        const void *ptr0, const void *ptr1)
{
  int           cmp;
  unsigned      k0,
                k1,
                dMsk,
                tMsk;
  WlzLBTDomain2D *lDom;
  WlzLBTNode2D  *nod0,
                *nod1;

  cmp = 0;
  nod0 = (WlzLBTNode2D *)ptr0;
  nod1 = (WlzLBTNode2D *)ptr1;
  lDom = (WlzLBTDomain2D *)ptrC;
  dMsk = 1 << (lDom->depth - 1);
  tMsk = nod0->keys[2] | nod1->keys[2];
  while((cmp == 0) && (dMsk != 0) && ((tMsk & dMsk) == 0))
  {
    k0 = nod0->keys[1] & dMsk;
    k1 = nod1->keys[1] & dMsk;
    if((cmp = k0 - k1) == 0)
    {
      k0 = nod0->keys[0] & dMsk;
      k1 = nod1->keys[0] & dMsk;
      cmp = k0 - k1;
    }
    dMsk >>= 1;
  }
  return(cmp);
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Condenses the sorted 2D LBT domain's nodes.
* 		Runs through the nodes looking for keys in which the key
*		values differ only in less significant digits and the
*		less significant digits form a complete set, eg:
*		321 323 330 331 332 333. The these keys are modified
*		using the term code, eg: 33X.
*		When all such sequences have been marked identical keys
*		are removed leaving the linear quadtree nodes.
   * where those which follow 
* \param	lDom
*/
static void	 WlzLBTCondenseNodes2D(WlzLBTDomain2D *lDom)
{
  int		idC,
  		idD,
		idP,
		idM,
		idN,
		dCnt,
  		nCnt,
		allPrvD;
  unsigned int	dMsk,
  		pDig,
  		cDig;
  WlzLBTNode2D	*cNod,
  		*pNod,
		*tNod;
  int		keyCnt[WLZ_LBTDOMAIN_MAXDIGITS];

  if(lDom->nNodes > 1)
  {
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, lDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
    /* Sort the nodes by key values. */
    AlgQSort(lDom->nodes, lDom->nNodes, sizeof(WlzLBTNode2D), lDom,
	     WlzLBTDomain2DNodeCmpFn);
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Sorted nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, lDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
    /* First pass: Mark key digits using the term digit when all previous
     * digits increment. */
    dMsk = 1;
    cNod = lDom->nodes;
    allPrvD = 1;
    for(idD = 0; idD < lDom->depth; ++idD)
    {
      allPrvD = allPrvD &&
                ((cNod->keys[0] & dMsk) == 0) &&
                ((cNod->keys[1] & dMsk) == 0);
      keyCnt[idD] = allPrvD;
      dMsk <<= 1;
    }
    for(idN = 1; idN < lDom->nNodes; ++idN)
    {
      idD = 0;
      dMsk = 1;
      dCnt = 4;
      allPrvD = 1;
      pNod = cNod++;
      while(allPrvD && (idD < lDom->depth))
      {
	pDig = ((pNod->keys[0] & dMsk) != 0) |
	       (((pNod->keys[1] & dMsk) != 0) << 1);
	cDig = ((cNod->keys[0] & dMsk) != 0) |
	       (((cNod->keys[1] & dMsk) != 0) << 1);
	if(idD == 0)
	{
	  if(cDig == 0)
	  {
	    keyCnt[idD] = 1;
	  }
	  else if(keyCnt[idD] && (cDig == pDig + 1))
	  {
	    ++keyCnt[idD];
	  }
	  else
	  {
	    keyCnt[idD] = 0;
	    allPrvD = 0;
	  }
	}
	else
	{
	  if((cDig == 0) && (keyCnt[idD - 1] == 1))
	  {
	    keyCnt[idD] = 1;
	  }
	  else if(keyCnt[idD] && ((cDig == pDig) || (cDig == pDig + 1)))
	  {
	    ++keyCnt[idD];
	  }
	  else
	  {
	    keyCnt[idD] = 0;
	    allPrvD = 0;
	  }
	}
	if(allPrvD && (keyCnt[idD] == dCnt))
	{
	  keyCnt[idD] = 0;
	  /* Mark key digits as term. */
	  tNod = cNod;
	  for(idM = 0; idM < dCnt; ++idM)
	  {
	    tNod->keys[0] &= ~dMsk;
	    tNod->keys[1] &= ~dMsk;
	    tNod->keys[2] |= dMsk;
	    --tNod;
	  }
	}
	++idD;
	dMsk <<= 1;
	dCnt *= 4;
      }
      while(idD < lDom->depth)
      {
        keyCnt[idD] = 0;
	++idD;
      }
    }
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Marked nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, lDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
    /* Second pass: Remove nodes with duplicate keys. */
    nCnt = 1;
    idP = 0;
    idC = 1;
    while(idC < lDom->nNodes)
    {
      pNod = lDom->nodes + idP;
      cNod = lDom->nodes + idC;
      if((cNod->keys[0] != pNod->keys[0]) ||
	 (cNod->keys[1] != pNod->keys[1]))
      {
	if(++idP != idC)
	{
	  pNod = lDom->nodes + idP;
	  pNod->keys[0] = cNod->keys[0];
	  pNod->keys[1] = cNod->keys[1];
	  pNod->keys[2] = cNod->keys[2];
	}
	++nCnt;
      }
      ++idC;
    }
    lDom->nNodes = nCnt;
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Condensed nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, lDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Outputs the nodes of the LBT ask text for testing.
* \param	fP		Output file.
* \param	lDom		The LBT domain to output.
*/
WlzErrorNum	WlzLBTTestOutputNodesTxt(FILE *fP, WlzLBTDomain2D *lDom)
{
  int		idD,
  		idN;
  WlzIBox2	nBB;
  WlzLBTNode2D	*nod;
  UBYTE		digits[WLZ_LBTDOMAIN_MAXDIGITS];
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(lDom == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(lDom->type != WLZ_LBTDOMAIN_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    (void )fprintf(fP, "type = %d\n", lDom->type);
    (void )fprintf(fP, "linkcount = %d\n", lDom->linkcount);
    (void )fprintf(fP, "freeptr = 0x%lx\n", (unsigned long )(lDom->freeptr));
    (void )fprintf(fP, "line1 %d\n", lDom->line1);
    (void )fprintf(fP, "lastln %d\n", lDom->lastln);
    (void )fprintf(fP, "kol1 %d\n", lDom->kol1);
    (void )fprintf(fP, "lastkl %d\n", lDom->lastkl);
    (void )fprintf(fP, "depth %d\n", lDom->depth);
    (void )fprintf(fP, "nNodes %d\n", lDom->nNodes);
    (void )fprintf(fP, "maxNodes %d\n", lDom->maxNodes);
    (void )fprintf(fP, "nodes:\n");
    nod = lDom->nodes;
    idN = 0;
    while(idN < lDom->nNodes)
    {
      (void )fprintf(fP, "%8d  ", idN);
      WlzLBTGetKeyDigits2D(nod->keys, digits);
      for(idD = 29; idD >= 0; --idD)
      {
        (void )fprintf(fP, "%d", digits[idD]);
      }
      WlzLBTKeyToBox2I(nod->keys, &nBB);
      (void )fprintf(fP, "  %d,%d,%d,%d\n",
                     nBB.xMin, nBB.yMin, nBB.xMax, nBB.yMax);
      ++idN;
      ++nod;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Outputs the nodes of the LBT as VTK polydata for testing.
* \param	fP		Output file.
* \param	lDom		The LBT domain to output.
*/
WlzErrorNum	WlzLBTTestOutputNodesVtk(FILE *fP, WlzLBTDomain2D *lDom)
{
  int		idN,
  		nOff;
  WlzLBTNode2D	*nod;
  WlzIBox2	*nodBB,
  		*nodeBB = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(lDom == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(lDom->type != WLZ_LBTDOMAIN_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(lDom->nNodes > 0)
  {
    if((nodeBB = (WlzIBox2 *)
                 AlcMalloc(sizeof(WlzIBox2) * lDom->nNodes)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      nodBB = nodeBB;
      nod = lDom->nodes;
      for(idN = 0; idN < lDom->nNodes; ++idN)
      {
        WlzLBTKeyToBox2I(nod->keys, nodBB);
	nodBB->xMax += 1;
	nodBB->yMax += 1;
        ++nod;
	++nodBB;
      }
      (void )fprintf(fP, "# vtk DataFile Version 1.0\n"
      		         "WlzLBTDomain2D test output\n"
			 "ASCII\n"
			 "DATASET POLYDATA\n"
			 "POINTS %d float\n",
			 lDom->nNodes * 4);
      nodBB = nodeBB;
      for(idN = 0; idN < lDom->nNodes; ++idN)
      {
        (void )fprintf(fP, "%d %d 0\n"
                           "%d %d 0\n"
                           "%d %d 0\n"
                           "%d %d 0\n",
			   nodBB->xMin, nodBB->yMin,
			   nodBB->xMax, nodBB->yMin,
			   nodBB->xMax, nodBB->yMax,
			   nodBB->xMin, nodBB->yMax);
        ++nodBB;
      }
      (void )fprintf(fP, "LINES %d %d\n",
                         lDom->nNodes * 4, lDom->nNodes * 12);
      for(idN = 0; idN < lDom->nNodes; ++idN)
      {
	nOff = idN * 4;
        (void )fprintf(fP, "2 %d %d\n"
	                   "2 %d %d\n"
			   "2 %d %d\n"
			   "2 %d %d\n",
			   nOff + 0, nOff + 1,
			   nOff + 1, nOff + 2,
			   nOff + 2, nOff + 3,
			   nOff + 3, nOff + 0);
      }
      AlcFree(nodeBB);
    }
  }
  return(errNum);
}

#ifdef WLZ_LBTDOMAIN_TEST_1

typedef enum _WlzLBTTestOutType
{
  WLZ_LBTDOMAIN_TEST_OUT_IDOM,
  WLZ_LBTDOMAIN_TEST_OUT_TXT,
  WLZ_LBTDOMAIN_TEST_OUT_VTK
} WlzLBTTestOutType;

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		balance = 0,
  		option,
  		ok = 1,
		maxNodSz = INT_MAX,
		maxBndNodSz = INT_MAX;
		usage = 0,
  		txtOut = 0,
  		vtkOut = 1;
  WlzLBTTestOutType outType = WLZ_LBTDOMAIN_TEST_OUT_VTK;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  FILE		*fP = NULL;
  WlzObject	*idObj = NULL,
  		*inObj = NULL,
  		*outObj = NULL;
  WlzDomain	outDom;
  WlzValues	nullVal;
  WlzLBTDomain2D *lDom = NULL;
  char		*inFileStr,
  		*outFileStr;
  static char	optList[] = "bitvho:",
  		outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  nullVal.core = NULL;
  outFileStr = outFileStrDef;
  inFileStr = inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        balance = 1;
	break;
      case 'i':
        outType = WLZ_LBTDOMAIN_TEST_OUT_IDOM;
	break;
      case 't':
        outType = WLZ_LBTDOMAIN_TEST_OUT_TXT;
	vtkOut = 0;
        break;
      case 'v':
	outType = WLZ_LBTDOMAIN_TEST_OUT_VTK;
        txtOut = 0;
        break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if((inFileStr == NULL) || (*inFileStr == '\0') ||
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
    }
    else
    {
      inFileStr = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inFileStr == NULL) || (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")? fopen(inFileStr, "r"):
                                       stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    if(inObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(inObj->type != WLZ_2D_DOMAINOBJ)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if(inObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: invalid object read from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
  }
  if(ok)
  {
    lDom = WlzLBTDomain2DFromDomain(inObj->domain, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to compute LBT domain (%d).\n",
		     *argv, errMsg);
    }
  }
  if(ok && balance)
  {
   idObj = WlzLBTMakeNodeIndexObj2D(lDom, inObj->domain.i, &errNum);
   if(errNum == WLZ_ERR_NONE)
   {
      errNum = WlzLBTBalanceDomain2D(lDom, idObj, maxNodSz, maxBndNodSz);
   }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
  }
  if(ok)
  {
    switch(outType)
    {
      case WLZ_LBTDOMAIN_TEST_OUT_IDOM:
        outDom.i = WlzLBTDomainToIDomain(lDom, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, outDom, nullVal,
	                       NULL, NULL, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzWriteObj(fP, outObj);
	}
	if(outObj)
	{
	  (void )WlzFreeObj(outObj);
	}
	else if(outDom.i)
	{
	  (void )WlzFreeIntervalDomain(outDom.i);
	}
	break;
      case WLZ_LBTDOMAIN_TEST_OUT_TXT:
        errNum = WlzLBTTestOutputNodesTxt(fP, lDom);
	break;
      case WLZ_LBTDOMAIN_TEST_OUT_VTK:
        errNum = WlzLBTTestOutputNodesVtk(fP, lDom);
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(idObj);
  (void )WlzFreeLBTDomain2D(lDom);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<output object>] [-h] [-b] [-t|v] [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output object file name.\n"
    "  -b  Produce a balanced tree.\n"
    "  -i  Woolz interval domain output.\n"
    "  -t  Text output.\n"
    "  -v  VTK output.\n");
  }
  return(!ok);
}
#endif /* WLZ_LBTDOMAIN_TEST_1 */
