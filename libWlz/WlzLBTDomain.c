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
* \enum		_WlzDirection
* \ingroup	WlzType
* \brief	Basic directions.
*		TODO Add to WlzType.h
*/
typedef enum _WlzDirection
{
  WLZ_DIRECTION_IC,		/*!< Increasing columns (++x, east). */
  WLZ_DIRECTION_IL,		/*!< Increasing lines (++y, south). */
  WLZ_DIRECTION_IP,		/*!< Increasing planes (++z, up). */
  WLZ_DIRECTION_DC,		/*!< Decreasing columns (--x, west). */
  WLZ_DIRECTION_DL,		/*!< Decreasing lines (--y, north). */
  WLZ_DIRECTION_DP		/*!< Decreasing planes (--z, down). */
} WlzDirection;

/*!
* \def		WLZ_LBTDOMAIN_MAXDIGITS
* \ingroup	WlzType
* \brief	The maximum number of linear binary tree key digits,
*		which must be less than the number of bits in an int.
*/
#define WLZ_LBTDOMAIN_MAXDIGITS (30)
/*!
* \def		WLZ_LBTDOMAIN_SPLIT
* \ingroup	WlzType
* \brief	Used to mark a linear binary tree split node.
*		A split node has:
*		  \li key[0] = index of the first node in the split.
*		  \li key[1] = index of the last node in the split.
*		  \li key[n] = WLZ_LBTDOMAIN_SPLIT (n = dimension).
*/
#define WLZ_LBTDOMAIN_SPLIT	(1<<((WLZ_LBTDOMAIN_MAXDIGITS)+1))


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
    
/*!
* \struct	_WlzLBTIdxElm
* \ingroup	DomainOps
* \brief	Element for use in doubley linked list of indices.
*/
typedef struct	_WlzLBTIdxElm
{
  int		idx;
  struct _WlzLBTIdxElm *next;
  struct _WlzLBTIdxElm *prev;
} WlzLBTIdxElm;

/*!
* \struct	_WlzLBTBalWsp2D
* \ingroup      DomainOps
* \brief 	Workspace used in balancing 2D linear binary tree domains.
*/
typedef struct _WlzLBTBalWsp2D
{
  WlzLBTIdxElm	*table[WLZ_LBTDOMAIN_MAXDIGITS]; /*!< Table of nodes (using
  					     the node indices) which are to
					     be split. */
  WlzLBTIdxElm	*pool;			/*!< Pool of index elements available
  					     for use. */
  WlzBlockStack	*res;			/*!< Block stack used to allocate the
  					     index elements.
} WlzLBTBalWsp2D;

/*!
* \struct	_WlzPartialItv2D
* \ingroup	DomainOps
* \brief	Data structure that can be used to hold partial intervals.
* 		These can then be sorted and condensed to find the intervals
*		for an interval domain.
*/

typedef struct _WlzPartialItv2D
{
  int			ileft;
  int			iright;
  int			ln;
} WlzPartialItv2D;

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
WlzLBTDomain2D			*WlzLBTDomain2DFromIDomain(
				  WlzIntervalDomain *iDom,
				  WlzErrorNum *dstErr);
WlzIntervalDomain 		*WlzLBTDomainToIDomain(
				  WlzLBTDomain2D *lDom,
				  WlzErrorNum *dstErr);
WlzIntervalDomain 		*WlzIDomainFromPItv2D(
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  int nPItv,
				  WlzPartialItv2D *pItv,
				  WlzErrorNum *dstErr);
WlzErrorNum			WlzLBTBalanceDomain2D(
				  WlzLBTDomain2D *sDom,
				  WlzIntervalDomain *iDom);
WlzErrorNum			WlzFreeLBTDomain2D(
				  WlzLBTDomain2D *lDom);
WlzErrorNum			WlzLBTTestOutputNodesTxt(
				  FILE *fP,
				  WlzLBTDomain2D *lDom);
int				WlzLBTNodeSz2D(
				  WlzLBTNode2D *nod);
int				WlzLBTNodeEdgeNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  int idN,
				  WlzDirection dir);
void				WlzLBTNodeEdgeNbr2D(
				  WlzLBTDomain2D *lDom,
				  int idN,
				  WlzDirection dir,
				  unsigned *keys);
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

static int			WlzPartialItv2DCmpFn(
				  const void *v0,
				  const void *v1);
static int			WlzLBTDomain2DNodeCmpFn(
				  const void *ptr0,
				  const void *ptr1);
static int			WlzLBTNodeIdxFromKeys2D(
				  WlzLBTDomain2D *lDom,
				  int idN,
				  unsigned *keys,
				  int *dstFound);
static void			WlzLBTCondenseNodes2D(
				  WlzLBTDomain2D *lDom);
static void			WlzLBTNodeEdgeNbr2DE(
				  WlzLBTDomain2D *lDom,
				  int idI,
				  WlzLBTNode2D *sNod,
				  unsigned *keys);
static void			WlzLBTNodeEdgeNbr2DS(
				  WlzLBTDomain2D *lDom,
				  int idI,
				  WlzLBTNode2D *sNod,
				  unsigned *keys);
static void			WlzLBTNodeEdgeNbr2DW(
				  WlzLBTDomain2D *lDom,
				  int idI,
				  WlzLBTNode2D *sNod,
				  unsigned *keys);
static void			WlzLBTNodeEdgeNbr2DN(
				  WlzLBTDomain2D *lDom,
				  int idI,
				  WlzLBTNode2D *sNod,
				  unsigned *keys);


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
	lDom = WlzLBTDomain2DFromIDomain(dom.i, &errNum);
        break;
      case WLZ_INTERVALDOMAIN_RECT:
	/* TODO there must be a fast way of generating the LBT from a
	 * rectangular interval domain. */
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
* \ingroup	DomainOps
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
* \return	Returns signed comparison for qsort().
* \ingroup	DomainOps
* \brief	Compares to partial intervals.
* \param	v0			First partial interval.
* \param	v1			Second partial interval.
*/
static int	WlzPartialItv2DCmpFn(const void *v0, const void *v1)
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
* \return 	Woolz interval domain or NULL on error.
* \ingroup	DomainOps
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
    qsort(pItv, nPItv, sizeof(WlzPartialItv2D), WlzPartialItv2DCmpFn);
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
* \ingroup	DomainOps
* \brief	Balances the given LBT domain so that the neighbouring
*		nodes of each node are either of the same size or differ
*		in size by a ratio of 2:1.
*		The neighbour finding algorithm used is quick and
*		simple but it requires an object in which the values
*		are set to the corresponding LBT domain indices.
*		For efficiency an associated interval domain may be given,
*		if the associated domain pointer is NULL then a domain
*		will be computed.
*		however it will still be possible to free it using
*		WlzFreeLBTDomain2D().
* \param	lDom			Given LBT domain.
* \param	aDom			Given assosciate interval domain
*					for neighbour finding.
*/
WlzErrorNum	WlzLBTBalanceDomain2D(WlzLBTDomain2D *lDom,
				      WlzIntervalDomain *iDom)
{
  int		idD,
  		idN,
  		idM,
		idS,
		sz0,
		sz1;
  int		idNN[4];
  WlzLBTNode2D	*nod[4];
  WlzLBTNode2D	*nod[4];
  WlzLBTBalWsp2D balWsp;
  WlzDomain	tDom;
  WlzObject	*iObj = NULL;
  WlzValues	tVal;
  WlzPixelV	iBkg;
  WlzObjectType iValTblType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzDirection dirLUT[4] =
                  {
		    WLZ_DIRECTION_IC,
		    WLZ_DIRECTION_IL,
		    WLZ_DIRECTION_DC,
		    WLZ_DIRECTION_DL
		  };

  /* Compute an index object which has values set to the LBT domain
   * node indices. This is used for neighbour finding. */
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
  }
  /* Put indices of all nodes which have a size (cell side length) >
   * twice that of the smallest neighbour into the split table.
   * The split table is indexed by log (base 2) of the node size with
   * a sorted (by node index) doublely linked list of node indices at each
   * split table entry. */
  (void )memset(balWsp, 0, sizeof(WlzLBTBalWsp2D));
  idN = 0;
  while((errNum == WLZ_ERR_NONE) && (idN < lDom->nNodes))
  {
    nod[0] = lDom->nodes + idN;
    sz0 = WlzLBTNodeSzLog2D(nod[0]);
    for(idD = 0; idD < 4; ++idD)
    {
      idM = WlzLBTNodeEdgeNbrIdx2D(lDom, idN, iObj, dirLUT[idD],
      				   WLZ_LBTDOMAIN_NBR_SMALLEST);
      if(idM > 0)
      {
        nod[1] = lDom->nodes + idM;
        sz1 = WlzLBTNodeSzLog2D(nod[1]);
	if(sz0 > sz1 + 1)
	{
	  /* Add child node to table of nodes to be split. */
	  errNum = WlzLBTSplitTblAddNod2D(&balWsp, sz0, lDom->nodes, idN);
	  break;
	}
      }
    }
    ++idN;
  }
  /* Split all nodes in the split table to get a balanced tree. */
  idS = WLZ_LBTDOMAIN_MAXDIGITS - 1;
  while((errNum == WLZ_ERR_NONE) && (idS > 1))
  {
    while(balanced.table[idS])
    {
      /* Get node to split. */
      idN = WlzLBTSplitTblGetNod2D(&balWsp, idS);
      /* Split the node. */
      idNN[0] = idN;
      idNN[1] = lDom->nNodes;
      idNN[2] = lDom->nNodes + 1;
      idNN[3] = lDom->nNodes + 2;
      lDom->nodes += 3;
      if(lDom->nodes > lDom->maxNodes)
      {
        /* Reallocte the nodes. */
	errNum = WLZ_ERR_UNIMPLEMENTED; /* TODO */
      }
      if(errNum == WLZ_ERR_NONE)
      {
	nod[0] = lDom->nodes + idNN[0];
	nod[1] = lDom->nodes + idNN[1];
	nod[2] = lDom->nodes + idNN[2];
	nod[3] = lDom->nodes + idNN[3];
	nod[1]->keys[0] = nod[0].keys[0] + sInc;
	nod[1]->keys[1] = nod[0].keys[1];
	nod[2]->keys[0] = nod[0].keys[0];
	nod[2]->keys[1] = nod[0].keys[1] + sInc;
	nod[3]->keys[0] = nod[0].keys[0] + sInc;
	nod[3]->keys[1] = nod[0].keys[1] + sInc;
	nod[3]->keys[2] = nod[2]->keys[2] = nod[1]->keys[2] =
			  nod[0]->keys[2] >>= 1;
      }
      /* Test each of the child node to see if it needs spliting and if
       * so add it to the split table in the balance workspace. */
      idN = 0;
      while((errNum == WLZ_ERR_NONE) && (idN < 4))
      {
	nod[0] = lDom->nodes + idNN[idN];
	sz0 = WlzLBTNodeSzLog2D(nod[0]);
	for(idD = 0; idD < 4; ++idD)
	{
	  idM = WlzLBTNodeEdgeNbrIdx2D(lDom, idNN[idN], iObj, dirLUT[idD],
	  			       WLZ_LBTDOMAIN_NBR_SMALLEST);
	  if(idM > 0)
	  {
            nod[1] = lDom->nodes + idM;
	    sz1 = WlzLBTNodeSzLog2D(nod[1]);
	    if(sz0 > sz1 + 1)
	    {
	      /* Add child node to table of nodes still to be split. */
	      errNum = WlzLBTSplitTblAddNod2D(&balWsp, sz0, lDom->nodes, idN);
	      break;
	    }
	  }
	}
        ++idN;
      }
    }
    ++idS;
  }
  (void )WlzFreeObj(iObj);
  /* Sort the nodes so that all the split nodes are in the appropriate
   * place. */
  if(errNum == WLZ_ERR_NONE)
  {
    qsort(lDom->nodes, lDom->nNodes, sizeof(WlzLBTNode2D),
          WlzLBTDomain2DNodeCmpFn);
  }
  return(errNum);
}

/*!
* \return	Node size: 0 on error else 1, 2, 4, ....
* \ingroup	DomainOps
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
* \ingroup	DomainOps
* \brief	Computes the log (base 2) of the size of a 2D linear binary
*		tree node.
* \param	nod			Given node.
*/
int		WlzLBTNodeSzLog2D(WlzLBTNode2D *nod)
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
* \return	Index of neighbor or < 0 if the neighbor is not in the
*		domain.
* \ingroup	DomainOps
* \brief	Finds the index of the neighbor of the \f$n\f$'th node in
*		direction \f$d\f$. For efficiency the validity of the given
*		parameters are not checked.
*		Finding the neighbor is a three stage process: First the
*		keys of the a neighbour in the appropriate direction are
*		computed then  the node which shares the most decomposition
*		cells with the computed keys is found. The node then found
*		may not be adjacent to the given node: Either because there
*		is no adjacent node (boundary) or because it is not the
*		correct node in the adjacent cell. In the third stage
*		these final checks are made.
* \param	lDom			Given 2D linear binary tree domain.
* \param	idN			Index of the \f$n\f$'th node.
* \param	dir			Direction \f$d\f$.
*/
int		WlzLBTNodeEdgeNbrIdx2D(WlzLBTDomain2D *lDom,
				       int idN, WlzDirection dir)
{
  int		idM = -1;
  unsigned 	keys[3];

  WlzLBTNodeEdgeNbr2D(lDom, idN, dir, keys);
  idM = WlzLBTNodeIdxFromKeys2D(lDom, idN, keys, NULL);
  /* HACK TODO idM might not be the appropriate node need to use direction. */
  return(idM);
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
  idMax = lDom->lastNode;
  cKeys = (lDom->nodes + idN)->keys;
  while((found == 0) && ((idMax - idMin) > 1))
  {
    cmp = WlzLBTDomain2DNodeCmpFn(gKeys, cKeys);
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
    if(!found)
    {
      cKeys = (lDom->nodes + idN)->keys;
      if((*(cKeys + 2) & WLZ_LBTDOMAIN_SPLIT) != 0)
      {
	/* Make sure that the target node is in the split before going in. */
        cmp = WlzLBTDomain2DNodeCmpFn(cKeys, gKeys);
	if(cmp > 0)
	{
	  idN = idN + 1;
	}
	else if(cmp < 0)
	{
	  idN = idN - 1;
	}
	else
	{
	  /* Update the values of idMin and idMax for the split. */
	  idMin = cKeys[0];
	  idMax = cKeys[1];
	  idN = (idMin + idMax) / 2;
	}
      }
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
* \ingroup	DomainOps
* \brief	Computes the keys of a neighbor of the \f$n\f$'th node in
*		direction \f$d\f$. For efficiency the validity of the given
*		parameters are not checked.
*		The keys computed are for a neighbor of the same size as
*		the given \f$n\f$'th node.
* \param	lDom			Given 2D linear binary tree domain.
* \param	idN			Index of the \f$n\f$'th node.
* \param	dir			Direction \f$d\f$.
* \param	keys			Keys to be set in line, column and
*					term order.
*/
void		WlzLBTNodeEdgeNbr2D(WlzLBTDomain2D *lDom,
				    int idN, WlzDirection dir,
				    unsigned *keys)
{
  int		idI = 0;
  unsigned	mskI = 1;
  WlzLBTNode2D	*sNod;

  /* Find first non term digit. */
  sNod = lDom->nodes + idN;
  while((idI < lDom->depth) && ((sNod->keys[2] & mskI) != 0))
  {
    mskI <<= 1;
    ++idI;
  }
  /* Find neighbor using first non term digit. */
  keys[2] = 0;
  switch(dir)
  {
    case WLZ_DIRECTION_IC:
      WlzLBTNodeEdgeNbr2DE(lDom, idI, sNod, keys);
      break;
    case WLZ_DIRECTION_IL:
      WlzLBTNodeEdgeNbr2DS(lDom, idI, sNod, keys);
      break;
    case WLZ_DIRECTION_DC:
      WlzLBTNodeEdgeNbr2DW(lDom, idI, sNod, keys);
      break;
    case WLZ_DIRECTION_DL:
      WlzLBTNodeEdgeNbr2DN(lDom, idI, sNod, keys);
      break;
  }
}

/*!
* \return	void
* \ingroup	DomainOps
* \brief	Computes the keys of a neighbor of the given node in
		an east direction. See WlzLBTNodeEdgeNbr2DS().
* \param	lDom			Given 2D linear binary tree.
* \param	idI			Index of first non term bit in the
*					location keys.
* \param	sNod			Given node.
* \param	keys			Keys to be set in line, column and
*					term order.
*/
static void	WlzLBTNodeEdgeNbr2DE(WlzLBTDomain2D *lDom, int idI,
				     WlzLBTNode2D *sNod, unsigned *keys)
{
  int		idJ;
  unsigned 	mskI,
  		mskJ,
		mskD;

  keys[1] = sNod->keys[1];
  /* Is source node at west of block? */
  mskI = 1 << idI;
  if((sNod->keys[0] & mskI) == 0)
  {
    /* Source node at west of parent block so neighbor is at the east of the
     * same block. */
    keys[0] = sNod->keys[0] | mskI;
  }
  else
  {
    /* Source node is to the east of parent block. */
    keys[0] = sNod->keys[0] ^ mskI;
    while(idI < lDom->depth)
    {
      mskI = 1 << idI;
      mskJ = 1 << (idI - 1);
      if((sNod->keys[0] & mskJ) == 0)
      {
	mskD = 1 << lDom->depth;
	while(mskI < mskD)
	{
	  keys[0] = (keys[0] & ~mskI) | ((sNod->keys[0] & mskI) ^ mskI);
          mskI <<= 1;
	}
	break;
      }
      else
      {
	keys[0] |= (sNod->keys[0] & mskI) ^ mskI;
	++idI;
      }
    }
  }
}

/*!
* \return	void
* \ingroup	DomainOps
* \brief	Computes the keys of a neighbor of the given node in
		an south direction.
*		The algorithm is based on the adjacency algorithm of
*		Gargantini, I. An Effective Way to Represent Quadtrees.
*		Communications of the ACM 25, 12 (Dec 1982), 905-910.
* \param	lDom			Given 2D linear binary tree.
* \param	idI			Index of first non term bit in the
*					location keys.
* \param	sNod			Given node.
* \param	keys			Keys to be set in line, column and
*					term order.
*/
static void	WlzLBTNodeEdgeNbr2DS(WlzLBTDomain2D *lDom, int idI,
				     WlzLBTNode2D *sNod, unsigned *keys)
{
  int		idJ;
  unsigned 	mskI,
  		mskJ,
		mskD;

  keys[0] = sNod->keys[0];
  /* Is source node at north of block? */
  mskI = 1 << idI;
  if((sNod->keys[1] & mskI) == 0)
  {
    /* Source node at north of parent block so neighbor is at the south of the
     * same block. */
    keys[1] = sNod->keys[1] | mskI;
  }
  else
  {
    /* Source node is to the south of parent block. */
    keys[1] = sNod->keys[1] ^ mskI;
    while(idI < lDom->depth)
    {
      mskI = 1 << idI;
      mskJ = 1 << (idI - 1);
      if((sNod->keys[1] & mskJ) == 0)
      {
	/* It's all the same from here up. */
	mskD = 1 << lDom->depth;
	while(mskI < mskD)
	{
	  keys[1] = (keys[1] & ~mskI) | ((sNod->keys[1] & mskI) ^ mskI);
          mskI <<= 1;
	}
	break;
      }
      else
      {
	/* North if south, south if north. */
	keys[1] |= (sNod->keys[1] & mskI) ^ mskI;
	++idI;
      }
    }
  }
}

/*!
* \return	void
* \ingroup	DomainOps
* \brief	Computes the keys of a neighbor of the given node in
		an west direction. See WlzLBTNodeEdgeNbr2DS().
* \param	lDom			Given 2D linear binary tree.
* \param	idI			Index of first non term bit in the
*					location keys.
* \param	sNod			Given node.
* \param	keys			Keys to be set in line, column and
*					term order.
*/
static void	WlzLBTNodeEdgeNbr2DW(WlzLBTDomain2D *lDom, int idI,
				     WlzLBTNode2D *sNod, unsigned *keys)
{
  int		idJ;
  unsigned 	mskI,
  		mskJ,
		mskD;

  keys[1] = sNod->keys[1];
  /* Is source node at east of block? */
  mskI = 1 << idI;
  if((sNod->keys[0] & mskI) != 0)
  {
    /* Source node at east of parent block so neighbor is at the west of the
     * same block. */
    keys[1] = sNod->keys[1] & ~mskI;
  }
  else
  {
    /* Source node is to the west of parent block. */
    keys[0] = sNod->keys[0] ^ mskI;

    while(idI < lDom->depth)
    {
      mskI = 1 << idI;
      mskJ = 1 << (idI - 1);
      if((sNod->keys[0] & mskJ) != 0)
      {
	mskD = 1 << lDom->depth;
	while(mskI < mskD)
	{
	  keys[0] = (keys[0] & ~mskI) | ((sNod->keys[0] & mskI) ^ mskI);
          mskI <<= 1;
	}
	break;
      }
      else
      {
	keys[0] |= (sNod->keys[0] & mskI) ^ mskI;
	++idI;
      }
    }
  }
}

/*!
* \return	void
* \ingroup	DomainOps
* \brief	Computes the keys of a neighbor of the given node in
		an north direction. See WlzLBTNodeEdgeNbr2DS().
* \param	lDom			Given 2D linear binary tree.
* \param	idI			Index of first non term bit in the
*					location keys.
* \param	sNod			Given node.
* \param	keys			Keys to be set in line, column and
*					term order.
*/
static void	WlzLBTNodeEdgeNbr2DN(WlzLBTDomain2D *lDom, int idI,
				     WlzLBTNode2D *sNod, unsigned *keys)
{
  int		idJ;
  unsigned 	mskI,
  		mskJ,
		mskD;

  keys[0] = sNod->keys[0];
  /* Is source node at south of block? */
  mskI = 1 << idI;
  if((sNod->keys[1] & mskI) != 0)
  {
    /* Source node at south of parent block so neighbor is at the north of the
     * same block. */
    keys[0] = sNod->keys[0] & ~mskI;
  }
  else
  {
    /* Source node is to the south of parent block. */
    keys[1] = sNod->keys[1] ^ mskI;
    while(idI < lDom->depth)
    {
      mskI = 1 << idI;
      mskJ = 1 << (idI - 1);
      if((sNod->keys[1] & mskJ) != 0)
      {
	mskD = 1 << lDom->depth;
	while(mskI < mskD)
	{
	  keys[1] = (keys[1] & ~mskI) | ((sNod->keys[1] & mskI) ^ mskI);
          mskI <<= 1;
	}
	break;
      }
      else
      {
	keys[1] |= (sNod->keys[1] & mskI) ^ mskI;
	++idI;
      }
    }
  }
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
static int	WlzLBTDomain2DNodeCmpFn(const void *ptr0, const void *ptr1)
{
  int		d0,
  		d1,
		cmp,
		idD;
  unsigned 	k0,
  		k1,
		dMsk,
		tMsk;
  WlzLBTNode2D	*nod0,
  		*nod1;

  cmp = 0;
  nod0 = (WlzLBTNode2D *)ptr0;
  nod1 = (WlzLBTNode2D *)ptr1;
  idD = WLZ_LBTDOMAIN_MAXDIGITS - 1; /* TODO Inefficient, use lDom->depth */
  dMsk = 1 << idD;
  tMsk = nod0->keys[2] | nod1->keys[2];
  while((cmp == 0) && (idD >= 0) && ((tMsk & dMsk) == 0))
  {
    k0 = (nod0->keys[0] & dMsk) != 0;
    k1 = (nod0->keys[1] & dMsk) != 0;
    d0 = (k1 << 1) + k0;
    k0 = (nod1->keys[0] & dMsk) != 0;
    k1 = (nod1->keys[1] & dMsk) != 0;
    d1 = (k1 << 1) + k0;
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
		usage = 0,
  		txtOut = 0,
  		vtkOut = 1;
  WlzLBTTestOutType outType = WLZ_LBTDOMAIN_TEST_OUT_VTK;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
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
    errNum = WlzLBTBalanceDomain2D(lDom, inObj->domain.i);
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
  if(lDom)
  {
    (void )WlzFreeLBTDomain2D(lDom);
  }
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
