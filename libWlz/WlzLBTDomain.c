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
* \ingroup
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

/* #define WLZ_LBTDOMAIN_DEBUG HACK */
#define WLZ_LBTDOMAIN_TEST_1 /* HACK */


enum _TODO /* HACK Enums for existing WlzType.h structures. */
{
  WLZ_LBTDOMAIN_2D
} TODO;

/*!
* \def		WLZ_LBTDOMAIN_MAXDIGITS
* \ingroup	WlzType
* \brief	The maximum number of linear binary tree key digits,
*		which must be less than the number of bits in an int.
*/
#define WLZ_LBTDOMAIN_MAXDIGITS (30)


/*!
* \struct	_WlzLBTNode2D
* \ingroup	WlzType
* \brief	A 2D linear binary tree node for spatial domain representation.
*		Typedef: ::WlzLBTNode.
*/
typedef struct _WlzLBTNode
{
  unsigned		keys[3];	/*!< A single location key which
  					     uses bit interleaving with
					     the bits interleaved between
					     the elements for column, line
					     and term in that order.
					     Each of the array elements must
					     have at least
					     WLZ_LBTDOMAIN_MAXDIGITS bits. */
} WlzLBTNode2D;

/*!
* \struct	_WlzLBTDomain2D
* \ingroup	WlzType
* \brief	A linear binary tree spatial domain representation.
*		Typedef: ::WlzLBTDomain.
*/
typedef struct _WlzLBTDomain2D
{
  WlzObjectType     	type;		/*!< From WlzCoreDomain. */
  int               	linkcount;	/*!< From WlzCoreDomain. */
  void              	*freeptr;	/*!< From WlzCoreDomain. */
  int               	line1;      	/*!< First line coordinate. */
  int               	lastln;		/*!< Last line coordinate. */
  int               	kol1;          	/*!< First column line
  					     coordinate. */
  int               	lastkl;        	/*!< Last column  line
  					     coordinate. */
  int		    	depth;		/*!< LBT depth. */
  int			nNodes;		/*!< Number of nodes in the
  					     tree. */
  int			maxNodes;	/*!< Number of nodes
  					     allocated. */
  WlzLBTNode2D 	    	*nodes;      	/*!< Array of nodes sorted by their
  					     location key. */
} WlzLBTDomain2D;

/* HACK For WlzProto.h. */
WlzLBTDomain2D			*WlzMakeLBTDomain2D(
				  WlzObjectType type,
				  int l1,
				  int ll,
				  int k1,
				  int kl,
				  WlzErrorNum *dstErr);
WlzLBTDomain2D			*WlzLBTDomain2DFromDomain(
				  WlzDomain dom,
				  WlzErrorNum *dstErr);
WlzLBTDomain2D			*WlzLBTDomain2DFromIntvlDomain(
				  WlzIntervalDomain *iDom,
				  WlzErrorNum *dstErr);
WlzErrorNum			WlzFreeLBTDomain2D(
				  WlzLBTDomain2D *lDom);
WlzErrorNum			WlzLBTTestOutputNodesTxt(
				  FILE *fP,
				  WlzLBTDomain2D *lDom);
void				WlzLBTPosToKey2D(
				  WlzIVertex2 pos,
				  unsigned *keys);
void				WlzLBTGetKeyDigits2D(
				  unsigned *keys,
				  UBYTE *digits);
void				WlzLBTKeyToPos2I(
				  unsigned *key,
				  WlzIVertex2 *pos);
void				WlzLBTKeyToBox2I(
				  int depth,
				  unsigned *key,
				  WlzIBox2 *box);

static int			WlzLBTDomain2DNodeSortFn(
				  const void *ptr0,
				  const void *ptr1);
static void			WlzLBTCondenseNodes2D(
				  WlzLBTDomain2D *lDom);


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
  
  /* TODO use dom.l = lDom instead of core! */
  dom.core = (WlzCoreDomain *)lDom;
  errNum = WlzFreeDomain(dom);
  return(errNum);
}

/*!
* \return	New 2D linear binary tree domain.
* \ingroup	DomainOps
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
	lDom = WlzLBTDomain2DFromIntvlDomain(dom.i, &errNum);
        break;
      case WLZ_INTERVALDOMAIN_RECT:
	/* TODO there must be a fast way of generating the LBT from a
	 * rectangular interval domain. */
	lDom = WlzLBTDomain2DFromIntvlDomain(dom.i, &errNum);
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
* \return	New 2D linear binary tree domain.
* \ingroup	DomainOps
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
WlzLBTDomain2D	*WlzLBTDomain2DFromIntvlDomain(WlzIntervalDomain *iDom,
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
* \return	void
* \ingroup	DomainOps
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
* \ingroup	DomainOps
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
* \ingroup	DomainOps
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
* \ingroup	DomainOps
* \brief	Sets bounding box which corresponds to the given LBT key.
* \param	depth			Depth of the LBT domain.
* \param	key			Keys in line, column and term order.
* \param	box			Position to be set.
*/
void		WlzLBTKeyToBox2I(int depth, unsigned *key, WlzIBox2 *box)
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
  box->xMax = box->xMin + sz;
  box->yMax = box->yMin + sz;
}

/*!
* \return	Sort value - less than, equal to, or greater than zero.
* \ingroup	DomainOps
* \brief	Sorts the keys into accending order by node value.
* \param	ptr0
* \param	ptr1
*/
static int	WlzLBTDomain2DNodeSortFn(const void *ptr0, const void *ptr1)
{
  int		d0,
  		d1,
		cmp,
		idD;
  unsigned 	k0,
  		k1,
		k2,
		dMsk;
  WlzLBTNode2D	*nod0,
  		*nod1;

  cmp = 0;
  nod0 = (WlzLBTNode2D *)ptr0;
  nod1 = (WlzLBTNode2D *)ptr1;
  idD = WLZ_LBTDOMAIN_MAXDIGITS - 1;
  dMsk = 1 << idD;
  while((cmp == 0) && (idD >= 0))
  {
    k0 = (nod0->keys[0] & dMsk) != 0;
    k1 = (nod0->keys[1] & dMsk) != 0;
    k2 = (nod0->keys[2] & dMsk) != 0;
    d0 = (k2 << 2) + (k1 << 1) + k0;
    k0 = (nod1->keys[0] & dMsk) != 0;
    k1 = (nod1->keys[1] & dMsk) != 0;
    k2 = (nod1->keys[2] & dMsk) != 0;
    d1 = (k2 << 2) + (k1 << 1) + k0;
    cmp = d0 - d1;
    --idD;
    dMsk >>= 1;
  }
  return(cmp);
}

/*!
* \return	void
* \ingroup	DomainOps
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
    qsort(lDom->nodes, lDom->nNodes, sizeof(WlzLBTNode2D),
	  WlzLBTDomain2DNodeSortFn);
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
* \ingroup	DomainOps
* \brief	Outputs the nodes of the LBT ask text for testing.
* \param	fP		Output file.
* \param	lDom		The LBT domain to output.
*/
WlzErrorNum	WlzLBTTestOutputNodesTxt(FILE *fP, WlzLBTDomain2D *lDom)
{
  int		idD,
  		idN;
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
      (void )fprintf(fP, "\n");
      ++idN;
      ++nod;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	DomainOps
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
        WlzLBTKeyToBox2I(lDom->depth, nod->keys, nodBB);
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
extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		option,
  		ok = 1,
		usage = 0,
  		txtOut = 0,
  		vtkOut = 1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  WlzLBTDomain2D *lDom = NULL;
  char		*inFileStr,
  		*outFileStr;
  static char	optList[] = "tvho:",
  		outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inFileStr = inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 't':
        txtOut = 1;
	vtkOut = 0;
        break;
      case 'v':
	vtkOut = 1;
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
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL) ||
       (txtOut &&
        ((errNum = WlzLBTTestOutputNodesTxt(fP, lDom)) != WLZ_ERR_NONE)) ||
       (vtkOut &&
        ((errNum = WlzLBTTestOutputNodesVtk(fP, lDom)) != WLZ_ERR_NONE)))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Write to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
  }
  if(lDom)
  {
    (void )WlzFreeLBTDomain2D(lDom);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<output object>] [-h] [-t|v] [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output object file name.\n"
    "  -t  Text output.\n"
    "  -v  VTK output.\n");
  }
  return(!ok);
}
#endif /* WLZ_LBTDOMAIN_TEST_1 */
