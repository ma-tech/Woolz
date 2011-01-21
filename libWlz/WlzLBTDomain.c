#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzLBTDomain_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzLBTDomain.c
* \author       Bill Hill
* \date         December 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	Functions for creating and manipulating linear binary
* 		tree domains.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>
#include <limits.h>
#include <float.h>

static int			WlzPartialItv2DCmpFn(
				  const void *vc,
				  const void *v0,
				  const void *v1);
static int			WlzPartialItv3DCmpFn(
				  const void *vc,
				  const void *v0,
				  const void *v1);
static int			WlzLBTDomain3DNodeCmpFn(
				  const void *lPtr,
				  const void *ptr0,
				  const void *ptr1);
static int			WlzLBTDomain2DNodeCmpFn(
				  const void *lPtr,
				  const void *ptr0,
				  const void *ptr1);
#ifdef WLZ_UNUSED_FUNCTIONS
static int			WlzLBTNodeIdxFromKeys2D(
				  WlzLBTDomain2D *lDom,
				  int idN,
				  unsigned *keys,
				  int *dstFound);
#endif /* WLZ_UNUSED_FUNCTIONS */
static int			WlzLBTBndEdgNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static int			WlzLBTBndEdgNbrIdx3D(
				  WlzLBTDomain3D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static int			WlzLBTMinLogSzEdgeNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static int			WlzLBTMinLogSzEdgeNbrIdx3D(
				  WlzLBTDomain3D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static int			WlzLBTBndEdgNbrDirIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir);
static int			WlzLBTBndEdgNbrDirIdx3D(
				  WlzLBTDomain3D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir);
static int			WlzLBTMinLogSzEdgeDirNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir);
static int			WlzLBTMinLogSzEdgeDirNbrIdx3D(
				  WlzLBTDomain3D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir);
static int			WlzLBTMaxLogSzEdgeDirNbrIdx2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir,
				  int *dstSz);
static int			WlzLBTMaxLogSzEdgeDirNbrIdx3D(
				  WlzLBTDomain3D *lDom,
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
#ifdef WLZ_UNUSED_FUNCTIONS
static int			WlzLBNodeAtBoundary3D(
				  WlzPlaneDomain *pDom,
				  WlzLBTNode3D *nod);
static int			WlzLBNodeAtBoundary2D(
				  WlzIntervalDomain *iDom,
				  WlzLBTNode2D *nod);
#endif /* WLZ_UNUSED_FUNCTIONS */
static int 			WlzLBTClassifyRotFromEdgeMask2D(
				  unsigned msk);
static unsigned 		WlzLBTIdxHashFn(
				  void *datum);
static void			WlzLBTCondenseNodes3D(
				  WlzLBTDomain3D *lDom,
				  int bndCndFlg);
static void			WlzLBTCondenseNodes2D(
				  WlzLBTDomain2D *lDom,
				  int bndCndFlg);
static void			WlzLBTSetNodeIndexObj2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static void			WlzLBTSetNodeIndexObj3D(
				  WlzLBTDomain3D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN);
static WlzErrorNum 		WlzLBTQueueInsert(
				  AlcCPQQueue *pQ,
				  AlcHashTable *hT,
				  int sz,
				  int idx);
static WlzLBTNodeClass2D 	WlzLBTClassifyNodeFromEdgeMask2D(
				  unsigned msk);

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
* \return	New 3D linear binary tree domain.
* \ingroup	WlzAllocation
* \brief	Creates a linear binary tree domain without creating any
*		nodes, leaving the nodes pointer NULL. Only the type,
*		bounding box and depth are set. The depth \f$d\f$ is set
*		such that
*		\f[
		  2^d >= \max(p_l - p_1, k_l - k_1, l_l - l_1) + 1
		\f]
*		where \f$k_1\f$, \f$k_l\f$, \f$l_1\f$ and \f$l_l\f$
*		are the first column, last column, first line and
*		last line respectively.
* \param	type			Type of domain, which must be
*					WLZ_LBTDOMAIN_3D.
* \param	p1			First plane.
* \param	pl			Last plane.
* \param	l1			First line.
* \param	ll			Last line.
* \param	k1			First column.
* \param	kl			Last column.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLBTDomain3D	*WlzMakeLBTDomain3D(WlzObjectType type,
				    int p1,
				    int pl,
				    int l1,
				    int ll,
				    int k1,
				    int kl,
				    WlzErrorNum *dstErr)
{
  unsigned int	sz;
  WlzLBTDomain3D *dom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((type != WLZ_LBTDOMAIN_3D) || (pl < p1) || (ll < l1) || (kl < k1))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if((dom = (WlzLBTDomain3D *)AlcCalloc(sizeof(WlzLBTDomain3D), 1)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom->type = type;
    dom->plane1 = p1;
    dom->lastpl = pl;
    dom->line1 = l1;
    dom->lastln = ll;
    dom->kol1 = k1;
    dom->lastkl = kl;
    /* Compute the depth. */
    sz = WLZ_MAX(kl - k1, ll - l1);
    sz = WLZ_MAX(sz, pl - p1) + 1;
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
* \brief	Frees the given 3D linear binary tree domain.
* \param	lDom			Given LBT domain.
*/
WlzErrorNum	WlzFreeLBTDomain3D(WlzLBTDomain3D *lDom)
{
  WlzDomain	dom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  dom.l3 = lDom;
  errNum = WlzFreeDomain(dom);
  return(errNum);
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
* \return	New 2D or 3D linear binary tree domain.
* \ingroup	WlzDomainOps
* \brief	Creates a new 2D or 3D linear binary tree domain from
		the given object's domain. The given object must be
		either a 2D or 3D domain object.
*		domain.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDomain	WlzLBTDomainFromObj(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzDomain	dom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	dom.l2 = WlzLBTDomain2DFromDomain(obj->domain, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
	dom.l3 = WlzLBTDomain3DFromDomain(obj->domain, &errNum);
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
  return(dom);
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
* \return	New 2D linear binary tree domain.
* \ingroup	WlzDomainOps
* \brief	Creates a new 2D linear binary tree domain from the given
*		domain.
* \param	dom			Given domain, which must be 2D.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLBTDomain3D	*WlzLBTDomain3DFromDomain(WlzDomain dom, WlzErrorNum *dstErr)
{
  WlzLBTDomain3D *lDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(dom.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(dom.core->type)
    {
      case WLZ_PLANEDOMAIN_DOMAIN:
	lDom = WlzLBTDomain3DFromPDomain(dom.p, &errNum);
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
* \return	New 3D plane domain.
* \ingroup	WlzDomainOps
* \brief	Creates a new plane domain from the given 3D linear
* 		binary tree domain.
* \param	lDom			Given 3D linear binary tree domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPlaneDomain *WlzLBTDomainToPDomain(WlzLBTDomain3D *lDom,
					 WlzErrorNum *dstErr)
{
  int		idI,
		idJ,
  		idN,
		nSz,
  		nPItv;
  WlzIVertex3	nPos;
  WlzLBTNode3D	*nod;
  WlzPartialItv3D *pItv0,
  		*pItvTb = NULL;
  WlzPlaneDomain *pDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  if(lDom == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(lDom->type != WLZ_LBTDOMAIN_3D)
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
      nSz = WlzLBTNodeSz3D(nod++);
      nPItv += nSz * nSz;
    }
    /* Allocate node partial intervals. */
    if((pItvTb = (WlzPartialItv3D *)AlcCalloc(sizeof(WlzPartialItv3D),
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
      nSz = WlzLBTNodeSz3D(nod); 
      WlzLBTKeyToPos3I(nod->keys, &nPos);
      for(idJ = 0; idJ < nSz; ++idJ)
      {
	for(idI = 0; idI < nSz; ++idI)
	{
	  pItv0->ileft = nPos.vtX;
	  pItv0->iright = nPos.vtX + nSz - 1;
	  pItv0->ln = nPos.vtY + idI;
	  pItv0->pl = nPos.vtZ + idJ;
	  ++pItv0;
	}
      }
      ++nod;
    }
    pDom = WlzPDomainFromPItv3D(lDom->plane1, lDom->lastpl,
   			        lDom->line1, lDom->lastln,
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
  return(pDom);
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
* \return 	Woolz interval domain or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Allocates and computes a plane domain from the given
*		table of partial intervals. The given partial intervals
*               must fit within the given domain bounding box.
* \param	plane1			First plane.
* \param	lastpl			Last plane.
* \param	line1			First line.
* \param	lastln			Last line.
* \param	kol1			First column.
* \param	lastkl			last column.
* \param	nPItv			Number of partial intervals.
* \param	pItv			Array of partial intervals.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPlaneDomain	*WlzPDomainFromPItv3D(int plane1, int lastpl,
				      int line1, int lastln,
				      int kol1, int lastkl,
				      int nPItv, WlzPartialItv3D *pItv,
				      WlzErrorNum *dstErr)
{
  int		idI,
		idP,
  		nItv,
		nItv2;
  WlzIBox2	bBox2;
  WlzDomain	tDom;
  WlzPartialItv3D *pItv0,
  		*pItv1,
		*pItv2;
  WlzIntervalLine *itvLn;
  WlzInterval	*itv0,
		*itv1,
  		*itvTb = NULL;
  WlzIntervalDomain *iDom;
  WlzPlaneDomain *pDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((plane1 > lastpl) || (line1 > lastln) || (kol1 > lastkl) ||
     (nPItv < 1) || (pItv == NULL))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(nPItv > 1)
  {
    /* Sort the partial intervals. */
    AlgQSort(pItv, nPItv, sizeof(WlzPartialItv3D), NULL, WlzPartialItv3DCmpFn);
    /* Condense the partial intervals into intervals. */
    idI = 0;
    nItv = 0;
    pItv0 = pItv;
    pItv1 = pItv + 1;
    while(idI < nPItv)
    {
      if((pItv0->pl == pItv1->pl) && (pItv0->ln == pItv1->ln) &&
         ((pItv0->iright + 1) == pItv1->ileft))
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
  }
  /* Create the plane domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    pDom = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
                              plane1, lastpl, line1, lastln, kol1, lastkl,
			      &errNum);
  }
  /* For each plane of the LBT domain create and set the interval domain. */
  idP = plane1;
  pItv0 = pItv;
  while((errNum == WLZ_ERR_NONE) && (idP <= lastpl))
  {
    pItv1 = pItv0;
    /* Count intervals in this plane and find their 2D bounding box. */
    nItv2 = 0;
    bBox2.xMin = pItv1->ileft;
    bBox2.xMax = pItv1->iright;
    bBox2.yMax = bBox2.yMin = pItv1->ln;
    while(pItv1->pl == idP)
    {
      if(bBox2.xMin > pItv1->ileft)
      {
        bBox2.xMin = pItv1->ileft;
      }
      if(bBox2.xMax < pItv1->iright)
      {
        bBox2.xMax = pItv1->iright;
      }
      if(bBox2.yMin > pItv1->ln)
      {
        bBox2.yMin = pItv1->ln;
      }
      else if(bBox2.yMax < pItv1->ln)
      {
        bBox2.yMax = pItv1->ln;
      }
      ++nItv2;
      ++pItv1;
    }
    iDom = NULL;
    if(nItv2 > 0)
    {
      /* Allocate an interval domain, interval lines and intervals
       * then set the intervals. */
      if(((iDom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					bBox2.yMin, bBox2.yMax,
					bBox2.xMin, bBox2.yMax,
					&errNum)) == NULL) ||
	 ((itvTb = (WlzInterval *)AlcCalloc(sizeof(WlzInterval),
					    nItv2)) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	iDom->freeptr = AlcFreeStackPush(iDom->freeptr, (void *)itvTb, NULL);
	/* Fill in the interval domain. */
	idI = 0;
	pItv2 = pItv1 = pItv0;
	itv1 = itv0 = itvTb;
	while(idI < nItv2)
	{
	  itvLn = iDom->intvlines + pItv1->ln - iDom->line1;
	  while((pItv2->ln == pItv1->ln) && (idI < nItv2))
	  {
	    itv1->ileft = pItv2->ileft;
	    itv1->iright = pItv2->iright;
	    ++(itvLn->nintvs);
	    ++itv1;
	    ++pItv2;
	    ++idI;
	  }
	  itvLn->intvs = itv0;
	  pItv1 = pItv2;
	  itv0 = itv1;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      tDom.i = iDom;
      *(pDom->domains + idP - plane1) = WlzAssignDomain(tDom, NULL);
    }
    pItv0 += nItv2;
    ++idP;
  }

  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreePlaneDomain(pDom);
    pDom = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pDom);
}

/*!
* \return	New 3D linear binary tree domain.
* \ingroup	WlzDomainOps
* \brief	Creates a new 3D linear binary tree domain from the given
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
* \param	pDom			Given plane domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzLBTDomain3D	*WlzLBTDomain3DFromPDomain(WlzPlaneDomain *pDom,
					WlzErrorNum *dstErr)
{
  int		idI,
		idP,
		pCnt;
  WlzObject	*obj = NULL,
  		*eObj = NULL;
  WlzLBTNode3D	*nod;
  WlzLBTDomain3D *lDom = NULL;
  WlzDomain	dom;
  WlzDomain	*dom2;
  WlzValues	nullVal;
  WlzIVertex3	pos;
  WlzIntervalWSpace iWsp;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.p = pDom;
  nullVal.core = NULL;
  lDom = WlzMakeLBTDomain3D(WLZ_LBTDOMAIN_3D,
			    pDom->plane1, pDom->lastpl,
			    pDom->line1, pDom->lastln,
			    pDom->kol1, pDom->lastkl, &errNum);
  /* Compute volume for the initial (and maximum) number of nodes. */
  idP = 0;
  pCnt = pDom->lastpl - pDom->plane1 + 1;
  while((errNum == WLZ_ERR_NONE) && (idP < pCnt))
  {
    if(((dom2 = pDom->domains + idP) != NULL) &&
       (dom2->core != NULL) &&
       (dom2->core->type != WLZ_EMPTY_DOMAIN))
    {
      obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dom2, nullVal, NULL, NULL,
			&errNum);
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
      (void )WlzFreeObj(obj);
    }
    ++idP;
  }
  /* Allocate nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((lDom->nodes = (WlzLBTNode3D *)AlcCalloc(sizeof(WlzLBTNode3D),
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
  /* Create an eroded domain to allow easy testing for boundary voxels. */
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, nullVal, NULL, NULL,
		      &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      eObj = WlzErosion(obj, WLZ_26_CONNECTED, &errNum);
    }
    (void )WlzFreeObj(obj);
  }
  /* Set keys of all initial nodes. */
  idP = 0;
  nod = lDom->nodes;
  while((errNum == WLZ_ERR_NONE) && (idP < pCnt))
  {
    if(((dom2 = pDom->domains + idP) != NULL) &&
       (dom2->core->type != WLZ_EMPTY_DOMAIN))
    {
      obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dom2, nullVal, NULL, NULL,
                          &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	pos.vtZ = pDom->plane1 + idP;
	errNum = WlzInitRasterScan(obj, &iWsp, WLZ_RASTERDIR_ILIC);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	while((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE )
	{
	  pos.vtY = iWsp.linpos;
	  pos.vtX = iWsp.lftpos;
	  for(idI = 0; idI < iWsp.colrmn; ++idI)
	  {
	    if(WlzInsideDomain(eObj, pos.vtZ, pos.vtY, pos.vtX, NULL) == 0)
	    {
	      nod->flags = WLZ_LBT_NODE_FLAG_BOUNDARY;
	    }
	    WlzLBTPosToKey3D(pos, nod->keys);
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
    }
    ++idP;
  }
  (void )WlzFreeObj(eObj);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Condense the nodes to get the linear binary tree. */
    WlzLBTCondenseNodes3D(lDom, 0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lDom);
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
  int		idI;
  WlzObject	*obj = NULL,
  		*eObj = NULL;
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
  /* Create an eroded domain to allow easy testing for boundary voxels. */
  if(errNum == WLZ_ERR_NONE)
  {
    eObj = WlzErosion(obj, WLZ_8_CONNECTED, &errNum);
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
	if(WlzInsideDomain(eObj, 0.0, pos.vtY, pos.vtX, NULL) == 0)
	{
	  nod->flags = WLZ_LBT_NODE_FLAG_BOUNDARY;
	}
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
  (void )WlzFreeObj(eObj);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Condense the nodes to get the linear binary tree. */
    WlzLBTCondenseNodes2D(lDom, 0);
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
* \param        maxBndSz		Maximum boundary node size
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
* \param        maxBndSz		Maximum boundary node size
*/
WlzErrorNum	WlzLBTBalanceDomain3D(WlzLBTDomain3D *lDom,
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
  int		idNN[8];
  WlzLBTNode3D	*nod[8];
  AlcCPQQueue	*pQ = NULL;
  AlcHashTable	*hT = NULL;
  WlzGreyValueWSpace *iGVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzDirection dirTab[6] = {
				   WLZ_DIRECTION_IC,
				   WLZ_DIRECTION_IL,
				   WLZ_DIRECTION_IP,
				   WLZ_DIRECTION_DC,
				   WLZ_DIRECTION_DL,
				   WLZ_DIRECTION_DP
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
    sz0 = WlzLBTNodeLogSz3D(nod[0]);
    if(sz0 > 0)
    {
      if(sz0 > maxSz)
      {
	flg = 1;
      }
      else if(sz0 > 1)
      {
	sz1 = WlzLBTMinLogSzEdgeNbrIdx3D(lDom, iGVWSp, idN);
	if((sz1 >= 0) && (sz0 > sz1 + 1))
	{
	  flg = 1;
	}
      }
      if(0 == flg)
      {
	flg = (sz0 > maxBndSz) &&
	      (WlzLBTBndEdgNbrIdx3D(lDom, iGVWSp, idN) != 0);
      }
      if(flg)
      {
	errNum = WlzLBTQueueInsert(pQ, hT, sz0, idN);
      }
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
      for(idM = 0; idM < 7; ++idM)
      {
        idNN[idM + 1] = lDom->nNodes + idM;
      }
      lDom->nNodes += 7;
      if(lDom->nNodes > lDom->maxNodes)
      {
        /* Reallocate the nodes. */
	sz0 = (lDom->maxNodes > 0)?
	      ((int )floor(lDom->nNodes / lDom->maxNodes) + 1) *
	      lDom->maxNodes: 1024;
	if((lDom->nodes = (WlzLBTNode3D *)
			  AlcRealloc(lDom->nodes,
			  	     sz0 * sizeof(WlzLBTNode3D))) == NULL)
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
        for(idM = 0; idM < 8; ++idM)
	{
	  nod[idM] = lDom->nodes + idNN[idM];
	}
	sInc = WlzLBTNodeSz3D(nod[0]) / 2;
	nod[1]->keys[0] = nod[0]->keys[0] + sInc;
	nod[1]->keys[1] = nod[0]->keys[1];
	nod[1]->keys[2] = nod[0]->keys[2];
	nod[2]->keys[0] = nod[0]->keys[0];
	nod[2]->keys[1] = nod[0]->keys[1] + sInc;
	nod[2]->keys[2] = nod[0]->keys[2];
	nod[3]->keys[0] = nod[0]->keys[0] + sInc;
	nod[3]->keys[1] = nod[0]->keys[1] + sInc;
	nod[3]->keys[2] = nod[0]->keys[2];
	nod[4]->keys[0] = nod[0]->keys[0];
	nod[4]->keys[1] = nod[0]->keys[1];
	nod[4]->keys[2] = nod[0]->keys[2] + sInc;
	nod[5]->keys[0] = nod[0]->keys[0] + sInc;
	nod[5]->keys[1] = nod[0]->keys[1];
	nod[5]->keys[2] = nod[0]->keys[2] + sInc;
	nod[6]->keys[0] = nod[0]->keys[0];
	nod[6]->keys[1] = nod[0]->keys[1] + sInc;
	nod[6]->keys[2] = nod[0]->keys[2] + sInc;
	nod[7]->keys[0] = nod[0]->keys[0] + sInc;
	nod[7]->keys[1] = nod[0]->keys[1] + sInc;
	nod[7]->keys[2] = nod[0]->keys[2] + sInc;
	nod[0]->keys[3] >>= 1;
	nod[7]->keys[3] = nod[6]->keys[3] =
	nod[5]->keys[3] = nod[4]->keys[3] =
	nod[3]->keys[3] = nod[2]->keys[3] =
	nod[1]->keys[3] = nod[0]->keys[3];
	/* Set indices in the index object. */
	for(idN = 0; idN < 8; ++idN)
	{
	  WlzLBTSetNodeIndexObj3D(lDom, iGVWSp, idNN[idN]);
	}
	/* Test each child node to see if either the child node or it's
	 * neighbours need to be inserted into the queue. */
	idN = 0;
	sz0 = WlzLBTNodeLogSz3D(nod[0]);
	while((errNum == WLZ_ERR_NONE) && (idN < 8))
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
	    sz1 = WlzLBTMinLogSzEdgeNbrIdx3D(lDom, iGVWSp, idNN[idN]);
	    if((sz1 >= 0) && (sz0 > sz1 + 1))
	    {
	      flg = 1;
	    }
	  }
	  if(0 == flg)
	  {
	    flg = (sz0 > maxBndSz) &&
		  (WlzLBTBndEdgNbrIdx3D(lDom, iGVWSp, idNN[idN]) != 0);
	  }
	  if(flg)
	  {
	    errNum = WlzLBTQueueInsert(pQ, hT, sz0, idNN[idN]);
	  }
	  /* Check for any neighbours with node size > twice that of
	   * the node. */
	  idM = 0;
	  while((errNum == WLZ_ERR_NONE) && (idM < 8))
	  {
	    idP = WlzLBTMaxLogSzEdgeDirNbrIdx3D(lDom, iGVWSp, idNN[idN],
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
    AlgQSort(lDom->nodes, lDom->nNodes, sizeof(WlzLBTNode3D), lDom,
             WlzLBTDomain3DNodeCmpFn);
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
* \return	New object with node indices or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Creates a new 3D domain object with integer values 
*		which are the indices of the nodes containing the
*		pixels of the domain.
* \param	lDom			Given LBT domain.
* \param	pDom			Corresponding interval domain,
					may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzLBTMakeNodeIndexObj3D(WlzLBTDomain3D *lDom,
				       WlzPlaneDomain *pDom,
				       WlzErrorNum *dstErr)
{
  int 		idN;
  WlzDomain	tDom;
  WlzValues	tVal;
  WlzPixelV	iBkg;
  WlzObjectType iValTblType;
  WlzGreyValueWSpace *iGVWSp = NULL;
  WlzObject	*iObj = NULL;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = -1;
  tDom.core = NULL;
  tVal.core = NULL;
  if(pDom)
  {
    tDom.p = pDom;
  }
  else
  {
    tDom.p = WlzLBTDomainToPDomain(lDom, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    iObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, tDom, tVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    iValTblType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
    iBkg.type = WLZ_GREY_INT;
    iBkg.v.ubv = 0;
    tVal.vox = WlzNewValuesVox(iObj, iValTblType, iBkg, &errNum);
    iObj->values = WlzAssignValues(tVal, NULL);
    iGVWSp = WlzGreyValueMakeWSp(iObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < lDom->nNodes; ++idN)
    {
      WlzLBTSetNodeIndexObj3D(lDom, iGVWSp, idN);
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
      if((pDom == NULL) && (tDom.core != NULL))
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
* \brief	Sets all index values in an existing 2D LBT node
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

  bgdV.type = WLZ_GREY_INT;
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
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Sets all index values in an existing 3D LBT node
*		index object.
* \param	lDom
* \param	iObj
*/
WlzErrorNum	WlzLBTIndexObjSetAllNodes3D(WlzLBTDomain3D *lDom,
				        WlzObject *iObj)
{
  int		idN;
  WlzPixelV	bgdV;
  WlzGreyValueWSpace *iGVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = -1;
  iGVWSp = WlzGreyValueMakeWSp(iObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < lDom->nNodes; ++idN)
    {
      WlzLBTSetNodeIndexObj3D(lDom, iGVWSp, idN);
    }
    (void )WlzSetBackground(iObj, bgdV);
  }
  WlzGreyValueFreeWSp(iGVWSp);
  return(errNum);
}

/*!
* \return	Counter-clockwise rotation in 90 degree steps,
*		with range [0-3].
* \ingroup	WlzDomainOps
* \brief	Looks up the rotation of a 2D node in steps of 90
*		degrees from the given edge mask.
* \param	msk			Given edge mask.
*/
static int 	WlzLBTClassifyRotFromEdgeMask2D(unsigned msk)
{
  int		rot;
  const int	rotTab[16] = {
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
  rot = rotTab[msk & 0xf];
  return(rot);
}

/*!
* \return	Node classification.
* \ingroup	WlzDomainOps
* \brief	Looks up the 2D node (or face of 3D node) class from the
*		given edge mask.
* \param	msk			Given edge mask.
*/
static WlzLBTNodeClass2D WlzLBTClassifyNodeFromEdgeMask2D(unsigned msk)
{
  WlzLBTNodeClass2D cls;
  const WlzLBTNodeClass2D clsTab[16] = {
				WLZ_LBT_NODE_CLASS_2D_0,	/*  0  o-o
						                       | |
						                       o-o */
				WLZ_LBT_NODE_CLASS_2D_1,	/*  1  o-o
						       		       | o
						       		       o-o */
				WLZ_LBT_NODE_CLASS_2D_1,	/*  2  ooo
						       		       | |
						       		       o-o */
				WLZ_LBT_NODE_CLASS_2D_2,	/*  3  ooo
						       		       | o
						       		       o-o */
				WLZ_LBT_NODE_CLASS_2D_1,	/*  4  o-o
						       		       o |
						       		       o-o */
				WLZ_LBT_NODE_CLASS_2D_3,	/*  5  o-o
						       		       o o
						       		       o-o */
				WLZ_LBT_NODE_CLASS_2D_2,	/*  6  ooo
						       		       o |
						       		       o-o */
				WLZ_LBT_NODE_CLASS_2D_4,	/*  7  ooo
						       		       o o
						       		       o-o */
				WLZ_LBT_NODE_CLASS_2D_1,	/*  8  o-o
						       		       | |
						       		       ooo */
				WLZ_LBT_NODE_CLASS_2D_2,	/*  9  o-o
						       		       | o
						       		       ooo */
				WLZ_LBT_NODE_CLASS_2D_3,	/* 10  ooo
						       		       | |
						       		       ooo */
				WLZ_LBT_NODE_CLASS_2D_4,	/* 11  ooo
						       		       | o
						       		       ooo */
				WLZ_LBT_NODE_CLASS_2D_2,	/* 12  o-o
						       		       o |
						       		       ooo */
				WLZ_LBT_NODE_CLASS_2D_4,	/* 13  o-o
						       		       o o
						       		       ooo */
				WLZ_LBT_NODE_CLASS_2D_4,	/* 14  ooo
						       		       o |
						       		       ooo */
				WLZ_LBT_NODE_CLASS_2D_5};	/* 15  ooo
						       		       o o
						       		       ooo */
  cls = clsTab[msk & 0xf];
  return(cls);
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

  nSz = WlzLBTNodeSz2D(lDom->nodes + idN);
  if(nSz == 1)
  {
    *dstRot = 0;
    *dstCls = WLZ_LBT_NODE_CLASS_2D_0;
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
    *dstCls = WlzLBTClassifyNodeFromEdgeMask2D(msk);
    *dstRot = WlzLBTClassifyRotFromEdgeMask2D(msk);
  }
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Classifies the given LBT node's face by it's connectivity
*		and returns it's class and the counter-clockwise rotation
*		of the basic class pattern, when viewed from outside of the
*		node's cube looking at the face, in multiples of 90 degrees.
*		This function uses the face index and orientation system
*		defined in WlzCMeshFromBalLBTDom3D().
* \param	lDom			Linear binary tree domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of the LBT node.
* \param	idF			Index of the face.
* \param	vtx			Coordinates of the 6 vertices
*					at the origin of faces of the
					LBT node's bounding box.
* \param	dstCls			Destination pointer for the class.
* \param	dstRot			Destination pointer for the rotation.
*/
void		WlzLBTClassifyNodeFace3D(WlzLBTDomain3D *lDom,
				     WlzGreyValueWSpace *iGVWSp,
				     int idN, int idF,
				     WlzDVertex3 *vtx,
				     WlzLBTNodeClass2D *dstCls,
				     int *dstRot)
{
  /* TODO Check this function! */
  int		idE,
		idO,
		idNN,
		nSz,
		nNSz,
		sSz;
  unsigned	msk;
  WlzDVertex3	off,
  		pos;
  const int 	offTbl0[6][3] = /* Offsets from 1st vertex of cube faces
  				 * to neighbours. These are small
				 * offsets to ensure the correct
				 * neighbouring face is located. */
     {
       { 0, -1,  0},
       { 0,  0, -1},
       { 1,  0,  0},
       { 0,  0,  1},
       {-1,  0,  0},
       { 0,  1,  0}
     };
  const	int	offTbl1[6][2][3] = /* Offsets from 1st vertex of cube faces
				    * to neighbours. These are multiples
				    * of the current cube size, with
				    * all as defined in the function
				    * WlzCMeshFromBalLBTDom3D() for o1 and
				    * o2, o3 = -o1 and o4 = -o2. */
     {
       /*     o1            o2   */
       {{ 1,  0,  0}, { 0,  0,  1}},
       {{-1,  0,  0}, { 0,  1,  0}},
       {{ 0,  0, -1}, { 0,  1,  0}},
       {{ 1,  0,  0}, { 0,  1,  0}},
       {{ 0,  0,  1}, { 0,  1,  0}},
       {{ 0,  0,  1}, { 1,  0,  0}}
     };

  msk = 0;
  nSz = WlzLBTNodeSz3D(lDom->nodes + idN);
  if(nSz == 1)
  {
    *dstRot = 0;
    *dstCls = WLZ_LBT_NODE_CLASS_2D_0;
  }
  else
  {
    /* First check the size of the neighbour which shares the current face.
     * If the neighbour is smaller than the current node or the node is on
     * the outer boundary, then set the flags for all it's edges, ie set
     * the mask value to 0xf. */
    off.vtX = offTbl0[idF][0];
    off.vtY = offTbl0[idF][1];
    off.vtZ = offTbl0[idF][2];
    WLZ_VTX_3_ADD(pos, vtx[idF], off);
    WlzGreyValueGet(iGVWSp, pos.vtZ, pos.vtY, pos.vtX);
    if(iGVWSp->bkdFlag)
    {
      msk = 0xf;
    }
    else
    {
      idNN = *(iGVWSp->gPtr[0].inp);
      if(idNN >= 0)
      {
        nNSz = WlzLBTNodeSz3D(lDom->nodes + idNN);
	if(nNSz < nSz)
	{
	  msk = 0xf;
	}
      }
    }
    /* If mask not set because of smaller node abutting face then check
     * the neighbours which share the edges of the given face. */
    if(msk == 0)
    {
      for(idE = 0; idE < 4; ++idE)
      {
	idO = (idE % 2) + 1;
	sSz = (idE > 1)? -nSz: nSz;
        off.vtX = offTbl1[idF][idO][0] * sSz;
	off.vtY = offTbl1[idF][idO][1] * sSz;
	off.vtZ = offTbl1[idF][idO][2] * sSz;
	WLZ_VTX_3_ADD(pos, vtx[idF], off);
	WlzGreyValueGet(iGVWSp, pos.vtZ, pos.vtY, pos.vtX);
	if(iGVWSp->bkdFlag)
	{
	  msk |= 1 << idE;
	}
	else
	{
	  idNN = *(iGVWSp->gPtr[0].inp);
	  if(idNN >= 0)
	  {
	    nNSz = WlzLBTNodeSz3D(lDom->nodes + idNN);
	    if(nNSz < nSz)
	    {
	      msk |= 1 << idE;
	    }
	  }
	}
      }
    }
    /* Look up the class and rotation. */
    *dstCls = WlzLBTClassifyNodeFromEdgeMask2D(msk);
    *dstRot = WlzLBTClassifyRotFromEdgeMask2D(msk);
  }
}

/*!
* \return	Returns signed comparison for AlgQSort().
* \ingroup	WlzDomainOps
* \brief	Compares two 2D partial intervals.
* \param	dummy			Unused client data.
* \param	v0			First partial interval.
* \param	v1			Second partial interval.
*/
static int	WlzPartialItv2DCmpFn(const void *dummy,
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
* \return	Returns signed comparison for AlgQSort().
* \ingroup	WlzDomainOps
* \brief	Compares two 3D partial intervals.
* \param	dummy			Unused client data.
* \param	v0			First partial interval.
* \param	v1			Second partial interval.
*/
static int	WlzPartialItv3DCmpFn(const void *dummy,
				     const void *v0, const void *v1)
{
  int		cmp;
  WlzPartialItv3D *pItv0,
  		*pItv1;

  pItv0 = (WlzPartialItv3D *)v0;
  pItv1 = (WlzPartialItv3D *)v1;
  if((cmp = pItv0->pl - pItv1->pl) == 0)
  {
    if((cmp = pItv0->ln - pItv1->ln) == 0)
    {
      cmp = pItv0->ileft - pItv1->ileft;
    }
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
  int		i0,
  		i1,
		cmp;

  i0 = (char *)datum0 - (char *)0;
  i1 = (char *)datum1 - (char *)0;
  cmp = i0 - i1;
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
  unsigned	dVal,
  		hVal;
  const unsigned p0 = 399989, /* These are three different 6 digit primes */
  		 p1 = 599999,
  		 p2 = 999983;

  dVal = (char *)datum - (char *)0;
  hVal = ((dVal + p0) * p1) % p2;
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
  char		*v;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

#ifdef WLZ_LBTDOMAIN_DEBUG
  (void )fprintf(stderr, "WlzLBTQueueInsert idx = %d sz = %d\n", idx, sz);
#endif /* WLZ_LBTDOMAIN_DEBUG */ 
  v = (char *)0 + idx;
  if(AlcHashItemGet(hT, v, NULL) == NULL)
  {
    if((AlcCPQEntryInsert(pQ, (float )sz, v) != ALC_ER_NONE) ||
       (AlcHashTableEntryInsert(hT, v, v, NULL) != ALC_ER_NONE))
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
    char	*v;

    v = pItem->entry;
    idx = v - (char *)0;
    AlcCPQItemFree(pQ, pItem);
    if((hItem = AlcHashItemGet(hT, v, NULL)) != NULL)
    {
      (void )AlcHashItemUnlink(hT, hItem, 1);
#ifdef WLZ_LBTDOMAIN_DEBUG
      (void )fprintf(stderr, "WlzLBTQueueUnlink idx = %d\n", idx);
#endif /* WLZ_LBTDOMAIN_DEBUG */ 
    }
  }
  return(idx);
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets values within a node's bounding box in an index
*		object for a 2D LBT. All parameters are assumed valid.
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
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets values within a node's bounding box in an index
*		object for a 3D LBT. All parameters are assumed valid.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
*/
static void	WlzLBTSetNodeIndexObj3D(WlzLBTDomain3D *lDom,
				      WlzGreyValueWSpace *iGVWSp, int idN)
{
  int		pX,
  		pY,
		pZ;
  WlzIBox3	nBox;
  WlzLBTNode3D	*nod;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox3I(nod->keys, &nBox);
  for(pZ = nBox.zMin; pZ <= nBox.zMax; ++pZ)
  {
    for(pY = nBox.yMin; pY <= nBox.yMax; ++pY)
    {
      for(pX = nBox.xMin; pX <= nBox.xMax; ++pX)
      {
	WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	*(iGVWSp->gPtr[0].inp) = idN;
      }
    }
  }
}

/*!
* \return	Non-zero if the node is on the boundary of the domain.
* \ingroup	WlzDomainOps
* \brief	Looks for a non-existant neighbour to the given node in the
*		2D LBT.
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
* \return	Non-zero if the node is on the boundary of the domain.
* \ingroup	WlzDomainOps
* \brief	Looks for a non-existant neighbour to the given node in the
*		3D LBT.
* \param        lDom                    Given LBT domain.
* \param        iGVWSp                  Grey workspace for index object.
* \param        idN                     Index of node in the LBT domain.
*/
static int	WlzLBTBndEdgNbrIdx3D(WlzLBTDomain3D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN)
{
  int		isBnd;

  isBnd = WlzLBTBndEdgNbrDirIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_IC);
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_IL);
  }
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_IP);
  }
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_DC);
  }
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_DL);
  }
  if(isBnd == 0)
  {
    isBnd = WlzLBTBndEdgNbrDirIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_DP);
  }
  return(isBnd);
}

/*!
* \return	Log of the size of the smallest neighbouring node or
		-ve if there is no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds the size of the smallest neighbouring node in a 2D LBT.
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
* \return	Log of the size of the smallest neighbouring node or
		-ve if there is no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds the size of the smallest neighbouring node in a 3D LBT.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
*/
static int	WlzLBTMinLogSzEdgeNbrIdx3D(WlzLBTDomain3D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN)
{
  int		sz,
  		minSz;

  minSz = WlzLBTMinLogSzEdgeDirNbrIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_IC);
  sz = WlzLBTMinLogSzEdgeDirNbrIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_IL);
  if((minSz < 0) || ((sz >= 0) && (sz < minSz)))
  {
    minSz = sz;
  }
  sz = WlzLBTMinLogSzEdgeDirNbrIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_IP);
  if((minSz < 0) || ((sz >= 0) && (sz < minSz)))
  {
    minSz = sz;
  }
  sz = WlzLBTMinLogSzEdgeDirNbrIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_DC);
  if((minSz < 0) || ((sz >= 0) && (sz < minSz)))
  {
    minSz = sz;
  }
  sz = WlzLBTMinLogSzEdgeDirNbrIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_DL);
  if((minSz < 0) || ((sz >= 0) && (sz < minSz)))
  {
    minSz = sz;
  }
  sz = WlzLBTMinLogSzEdgeDirNbrIdx3D(lDom, iGVWSp, idN, WLZ_DIRECTION_DP);
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
* \brief	Looks for a non-existant neighbour to the given 2D LBT node
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
    default:
      break;
  }
  return(isBnd);
}

/*!
* \return	Non-zero if the node's neighbour in the given diection is
*		outside the domain.
* \ingroup	WlzDomainOps
* \brief	Looks for a non-existant neighbour to the given 3D LBT node
*		in the given direction.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
* \param	dir			Given direction.
*/
static int	WlzLBTBndEdgNbrDirIdx3D(WlzLBTDomain3D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN, WlzDirection dir)
{
  int		pX,
  		pY,
		pZ,
		isBnd = 0;
  WlzLBTNode3D	*nod;
  WlzIBox3	nBB;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox3I(nod->keys, &nBB);
  switch(dir)
  {
    case WLZ_DIRECTION_IC: /* FALLTHROUGH */
    case WLZ_DIRECTION_DC:
      pX = (dir == WLZ_DIRECTION_IC)? nBB.xMax + 1: nBB.xMin - 1;
      pY = nBB.yMin;
      pZ = nBB.zMin;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      isBnd = iGVWSp->gVal[0].inv < 0;
      while(!isBnd && (++pZ <= nBB.zMax))
      {
	while(!isBnd && (++pY <= nBB.yMax))
	{
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  isBnd = iGVWSp->gVal[0].inv < 0;
	}
      }
      break;
    case WLZ_DIRECTION_IL: /* FALLTHROUGH */
    case WLZ_DIRECTION_DL:
      pX = nBB.xMin;
      pY = (dir == WLZ_DIRECTION_IL)? nBB.yMax + 1: nBB.yMin - 1;
      pZ = nBB.zMin;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      isBnd = iGVWSp->gVal[0].inv < 0;
      while(!isBnd && (++pZ <= nBB.zMax))
      {
	while(!isBnd && (++pX <= nBB.xMax))
	{
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  isBnd = iGVWSp->gVal[0].inv < 0;
	}
      }
      break;
    case WLZ_DIRECTION_IP: /* FALLTHROUGH */
    case WLZ_DIRECTION_DP:
      pX = nBB.xMin;
      pY = nBB.yMin;
      pZ = (dir == WLZ_DIRECTION_IP)? nBB.zMax + 1: nBB.zMin - 1;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      isBnd = iGVWSp->gVal[0].inv < 0;
      while(!isBnd && (++pY <= nBB.yMax))
      {
	while(!isBnd && (++pX <= nBB.xMax))
	{
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  isBnd = iGVWSp->gVal[0].inv < 0;
	}
      }
      break;
  }
  return(isBnd);
}


/*!
* \return	Log of the size of the smallest neighbouring node or
*		-ve if there is no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds the size of the smallest neighbouring node in the
*		given direction in a 2D LBT.
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
    default:
      break;
  }
  return(minSz);
}

/*!
* \return	Log of the size of the smallest neighbouring node or
*		-ve if there is no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds the size of the smallest neighbouring node in the
*		given direction in a 3D LBT.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
* \param	dir			Given direction.
*/
static int	WlzLBTMinLogSzEdgeDirNbrIdx3D(WlzLBTDomain3D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN, WlzDirection dir)
{
  int		id0,
  		id1,
		pX,
  		pY,
		pZ,
		sz,
		minSz = -1;
  WlzLBTNode3D	*nod;
  WlzIBox3	nBB;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox3I(nod->keys, &nBB);
  switch(dir)
  {
    case WLZ_DIRECTION_IC: /* FALLTHROUGH */
    case WLZ_DIRECTION_DC:
      pX = (dir == WLZ_DIRECTION_IC)? nBB.xMax + 1: nBB.xMin - 1;
      pY = nBB.yMin;
      pZ = nBB.zMin;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      id1 = iGVWSp->gVal[0].inv;
      minSz = sz = (id1 >= 0)? WlzLBTNodeLogSz3D(lDom->nodes + id1): -1;
      while(++pZ <= nBB.zMax)
      {
	while(++pY <= nBB.yMax)
	{
	  id0 = id1;
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  id1 = iGVWSp->gVal[0].inv;
	  if((id1 != id0) && (id1 >= 0))
	  {
	    sz = WlzLBTNodeLogSz3D(lDom->nodes + id1);
	    if((sz >= 0) && ((minSz < 0) || (sz < minSz)))
	    {
	      minSz = sz;
	    }
	  }
	}
      }
      break;
    case WLZ_DIRECTION_IL: /* FALLTHROUGH */
    case WLZ_DIRECTION_DL:
      pX = nBB.xMin;
      pY = (dir == WLZ_DIRECTION_IL)? nBB.yMax + 1: nBB.yMin - 1;
      pZ = nBB.zMin;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      id1 = iGVWSp->gVal[0].inv;
      minSz = sz = (id1 >= 0)? WlzLBTNodeLogSz3D(lDom->nodes + id1): -1;
      while(++pZ <= nBB.zMax)
      {
	while(++pX <= nBB.xMax)
	{
	  id0 = id1;
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  id1 = iGVWSp->gVal[0].inv;
	  if((id1 != id0) && (id1 >= 0))
	  {
	    sz = WlzLBTNodeLogSz3D(lDom->nodes + id1);
	    if((sz >= 0) && ((minSz < 0) || (sz < minSz)))
	    {
	      minSz = sz;
	    }
	  }
	}
      }
      break;
    case WLZ_DIRECTION_IP: /* FALLTHROUGH */
    case WLZ_DIRECTION_DP:
      pX = nBB.xMin;
      pY = nBB.yMin;
      pZ = (dir == WLZ_DIRECTION_IP)? nBB.zMax + 1: nBB.zMin - 1;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      id1 = iGVWSp->gVal[0].inv;
      minSz = sz = (id1 >= 0)? WlzLBTNodeLogSz3D(lDom->nodes + id1): -1;
      while(++pY <= nBB.yMax)
      {
	while(++pX <= nBB.xMax)
	{
	  id0 = id1;
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  id1 = iGVWSp->gVal[0].inv;
	  if((id1 != id0) && (id1 >= 0))
	  {
	    sz = WlzLBTNodeLogSz3D(lDom->nodes + id1);
	    if((sz >= 0) && ((minSz < 0) || (sz < minSz)))
	    {
	      minSz = sz;
	    }
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
    default:
      break;
  }
  return(nCnt);
}

/*!
* \return	Index of maximum sized neighbour or -ve if no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds a maximum sized neighbour in the given direction
*		in the 2D LBT and then returns it's index and size.
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
    default:
      break;
  }
  *dstSz = szM;
  return(idM);
}

/*!
* \return	Index of maximum sized neighbour or -ve if no neighbour.
* \ingroup	WlzDomainOps
* \brief	Finds a maximum sized neighbour in the given direction
*		in the 3D LBT and then returns it's index and size.
* \param	lDom			Given LBT domain.
* \param	iGVWSp			Grey workspace for index object.
* \param	idN			Index of node in the LBT domain.
* \param	dir			Given direction.
* \param	dstSz			Destination pointer for log size
*					of neighbour.
*/
static int	WlzLBTMaxLogSzEdgeDirNbrIdx3D(WlzLBTDomain3D *lDom,
					WlzGreyValueWSpace *iGVWSp,
					int idN, WlzDirection dir,
					int *dstSz)
{
  int		id0,
  		id1,
		pX,
  		pY,
		pZ,
		sz,
		idM = -1,
  		szM = -1;
  WlzLBTNode3D	*nod;
  WlzIBox3	nBB;

  nod = lDom->nodes + idN;
  WlzLBTKeyToBox3I(nod->keys, &nBB);
  switch(dir)
  {
    case WLZ_DIRECTION_IC: /* FALLTHROUGH */
    case WLZ_DIRECTION_DC:
      pX = (dir == WLZ_DIRECTION_IC)? nBB.xMax + 1: nBB.xMin - 1;
      pY = nBB.yMin;
      pZ = nBB.zMin;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      idM = id1 = iGVWSp->gVal[0].inv;
      szM = sz = (id1 >= 0)? WlzLBTNodeLogSz3D(lDom->nodes + id1): -1;
      while(++pZ <= nBB.zMax)
      {
	while(++pY <= nBB.yMax)
	{
	  id0 = id1;
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  id1 = iGVWSp->gVal[0].inv;
	  if((id1 != id0) && (id1 >= 0))
	  {
	    sz = WlzLBTNodeLogSz3D(lDom->nodes + id1);
	    if((sz >= 0) && ((szM < 0) || (sz > szM)))
	    {
	      szM = sz;
	      idM = id1;
	    }
	  }
	}
      }
      break;
    case WLZ_DIRECTION_IL: /* FALLTHROUGH */
    case WLZ_DIRECTION_DL:
      pX = nBB.xMin;
      pY = (dir == WLZ_DIRECTION_IL)? nBB.yMax + 1: nBB.yMin - 1;
      pZ = nBB.zMin;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      idM = id1 = iGVWSp->gVal[0].inv;
      szM = sz = (id1 >= 0)? WlzLBTNodeLogSz3D(lDom->nodes + id1): -1;
      while(++pZ <= nBB.zMax)
      {
	while(++pX <= nBB.xMax)
	{
	  id0 = id1;
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  id1 = iGVWSp->gVal[0].inv;
	  if((id1 != id0) && (id1 >= 0))
	  {
	    sz = WlzLBTNodeLogSz3D(lDom->nodes + id1);
	    if((sz >= 0) && ((szM < 0) || (sz > szM)))
	    {
	      szM = sz;
	      idM = id1;
	    }
	  }
	}
      }
      break;
    case WLZ_DIRECTION_IP: /* FALLTHROUGH */
    case WLZ_DIRECTION_DP:
      pX = nBB.xMin;
      pY = nBB.yMin;
      pZ = (dir == WLZ_DIRECTION_IP)? nBB.zMax + 1: nBB.zMin - 1;
      WlzGreyValueGet(iGVWSp, pZ, pY, pX);
      idM = id1 = iGVWSp->gVal[0].inv;
      szM = sz = (id1 >= 0)? WlzLBTNodeLogSz3D(lDom->nodes + id1): -1;
      while(++pY <= nBB.yMax)
      {
	while(++pX <= nBB.xMax)
	{
	  id0 = id1;
	  WlzGreyValueGet(iGVWSp, pZ, pY, pX);
	  id1 = iGVWSp->gVal[0].inv;
	  if((id1 != id0) && (id1 >= 0))
	  {
	    sz = WlzLBTNodeLogSz3D(lDom->nodes + id1);
	    if((sz >= 0) && ((szM < 0) || (sz > szM)))
	    {
	      szM = sz;
	      idM = id1;
	    }
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
* \return	Node size: 0 on error else 1, 2, 4, ....
* \ingroup	WlzDomainOps
* \brief	Computes the size of a 3D linear binary tree node.
* \param	nod			Given node.
*/
int		WlzLBTNodeSz3D(WlzLBTNode3D *nod)
{
  int		key,
  		sz = 0;

  if(nod)
  {
    key = nod->keys[3];
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
* \brief	Computes the log (base 2) of the size of a 3D linear binary
*		tree node.
* \param	nod			Given node.
*/
int		WlzLBTNodeLogSz3D(WlzLBTNode3D *nod)
{
  int		key,
  		sz = -1;

  if(nod)
  {
    key = nod->keys[3];
    sz = 0;
    while(key)
    {
      key >>= 1;
      ++sz;
    }
  }
  return(sz);
}

#ifdef WLZ_UNUSED_FUNCTIONS
/*!
* \return
* \ingroup      WlzDomainOps
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
static int	WlzLBTNodeIdxFromKeys2D(WlzLBTDomain2D *lDom, int idN,
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
#endif /* WLZ_UNUSED_FUNCTIONS */

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets the value of the LBT key for the given position,
*		where the position is relative to the first line and
*		column of the domain.
*		The key is encoded by interleaving the bits of the
*		column coordinate, line coordinate, plane coordinate
*		and term value in that order.
*		The keys are ordered so that bit 0 is least significant.
* \param	pos			Position.
* \param	keys			Keys to be set in line, column and
*					term order.
*/
void		WlzLBTPosToKey3D(WlzIVertex3 pos, unsigned *keys)
{
  keys[0] = pos.vtX;
  keys[1] = pos.vtY;
  keys[2] = pos.vtZ;
  keys[3] = 0;
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
*		from the given 3D keys.
* \param	keys			Keys to be read.
* \param	digits			Array of WLZ_LBTDOMAIN_MAXDIGITS
*					digits.
*/
void		WlzLBTGetKeyDigits3D(unsigned *keys, WlzUByte *digits)
{
  int		idD;
  unsigned	k0,
  		k1,
		k2,
		k3,
		dMsk;

  idD = 0;
  dMsk = 1;
  while(idD < WLZ_LBTDOMAIN_MAXDIGITS)
  {
    k0 = (keys[0] & dMsk) != 0; 
    k1 = (keys[1] & dMsk) != 0; 
    k2 = (keys[2] & dMsk) != 0; 
    k3 = (keys[3] & dMsk) != 0; 
    digits[idD] = (WlzUByte )(k0 + (k1 << 1) + (k2 << 2) + (k3 << 3));
    ++idD;
    dMsk <<= 1;
  }
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Gets an array of WLZ_LBTDOMAIN_MAXDIGITS location digits
*		from the given 2D keys.
* \param	keys			Keys to be read.
* \param	digits			Array of WLZ_LBTDOMAIN_MAXDIGITS
*					digits.
*/
void		WlzLBTGetKeyDigits2D(unsigned *keys, WlzUByte *digits)
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
    digits[idD] = (WlzUByte )(k0 + (k1 << 1) + (k2 << 2));
    ++idD;
    dMsk <<= 1;
  }
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets position (relative to the first line and column of the
* 		domain) which corresponds to the given 3D LBT key.
* \param	key			Keys in line, column and term order.
* \param	pos			Position to be set.
*/
void		WlzLBTKeyToPos3I(unsigned *key, WlzIVertex3 *pos)
{
  pos->vtX = key[0];
  pos->vtY = key[1];
  pos->vtZ = key[2];
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets position (relative to the first line and column of the
* 		domain) which corresponds to the given 3D LBT key.
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
* \brief	Sets bounding box which corresponds to the given 3D LBT key.
* \param	key			Keys in line, column and term order.
* \param	box			Position to be set.
*/
void		WlzLBTKeyToBox3I(unsigned *key, WlzIBox3 *box)
{
  unsigned int	k,
  		sz = 1;

  box->xMin = key[0];
  box->yMin = key[1];
  box->zMin = key[2];
  if(key[3] != 0)
  {
    k = key[3];
    while(k)
    {
      sz <<= 1;
      k >>= 1;
    }
    box->xMax = box->xMin + sz - 1;
    box->yMax = box->yMin + sz - 1;
    box->zMax = box->zMin + sz - 1;
  }
  else
  {
    box->xMax = box->xMin;
    box->yMax = box->yMin;
    box->zMax = box->zMin;
  }
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Sets bounding box which corresponds to the given 2D LBT key.
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
    box->xMax = box->xMin + sz - 1;
    box->yMax = box->yMin + sz - 1;
  }
  else
  {
    box->xMax = box->xMin;
    box->yMax = box->yMin;
  }
}

/*!
* \return	Sort value - less than, equal to, or greater than zero.
* \ingroup	WlzDomainOps
* \brief	Sorts the 3D keys into accending order by node value.
* \param	ptrC			Used to pass the domain pointer.
* \param	ptr0			Used to pass first node.
* \param	ptr1			Used to pass second node.
*/
static int	WlzLBTDomain3DNodeCmpFn(const void *ptrC,
				        const void *ptr0, const void *ptr1)
{
  int           cmp;
  unsigned      k0,
                k1,
                dMsk,
                tMsk;
  WlzLBTDomain3D *lDom;
  WlzLBTNode3D  *nod0,
                *nod1;

  cmp = 0;
  nod0 = (WlzLBTNode3D *)ptr0;
  nod1 = (WlzLBTNode3D *)ptr1;
  lDom = (WlzLBTDomain3D *)ptrC;
  dMsk = 1 << (lDom->depth - 1);
  tMsk = nod0->keys[3] | nod1->keys[3];
  while((cmp == 0) && (dMsk != 0) && ((tMsk & dMsk) == 0))
  {
    k0 = nod0->keys[2] & dMsk;
    k1 = nod1->keys[2] & dMsk;
    if((cmp = k0 - k1) == 0)
    {
      k0 = nod0->keys[1] & dMsk;
      k1 = nod1->keys[1] & dMsk;
      if((cmp = k0 - k1) == 0)
      {
	k0 = nod0->keys[0] & dMsk;
	k1 = nod1->keys[0] & dMsk;
	cmp = k0 - k1;
      }
    }
    dMsk >>= 1;
  }
  return(cmp);
}

/*!
* \return	Sort value - less than, equal to, or greater than zero.
* \ingroup	WlzDomainOps
* \brief	Sorts the 2D keys into accending order by node value.
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
* \brief	Condenses the sorted 3D LBT domain's nodes.
* 		Runs through the nodes looking for keys in which the key
*		values differ only in less significant digits and the
*		less significant digits form a complete set, eg:
*		330 331 332 333 ... 337. The these keys are modified
*		using the term code, eg: 33X.
*		When all such sequences have been marked identical keys
*		are removed leaving the linear quadtree nodes.
* \param	lDom			The LBT domain.
* \param	bndCndFlg		If non-zero then condensation of
*					boundary nodes is allowed.
*/
static void	 WlzLBTCondenseNodes3D(WlzLBTDomain3D *lDom,
				       int bndCndFlg)
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
#ifdef WLZ_LBTDOMAIN_DEBUG
  WlzDomain	tDom;
#endif
  WlzLBTNode3D	*cNod,
  		*pNod,
		*tNod;
  int		keyCnt[WLZ_LBTDOMAIN_MAXDIGITS];

  if(lDom->nNodes > 1)
  {
#ifdef WLZ_LBTDOMAIN_DEBUG
    tDom.l3 = lDom;
    (void )fprintf(stderr, "Nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
    /* Sort the nodes by key values. */
    AlgQSort(lDom->nodes, lDom->nNodes, sizeof(WlzLBTNode3D), lDom,
	     WlzLBTDomain3DNodeCmpFn);
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Sorted nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
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
                ((cNod->keys[1] & dMsk) == 0) &&
                ((cNod->keys[2] & dMsk) == 0);
      keyCnt[idD] = allPrvD;
      dMsk <<= 1;
    }
    for(idN = 1; idN < lDom->nNodes; ++idN)
    {
      dMsk = 1;
      dCnt = 8;
      allPrvD = 1;
      pNod = cNod++;
      /* Avoid condensing nodes which are on the boundary of the domain. */
      if((bndCndFlg == 0) && ((pNod->flags & WLZ_LBT_NODE_FLAG_BOUNDARY) != 0))
      {
	for(idD = 0; idD < lDom->depth; ++idD)
	{
	  allPrvD = allPrvD &&
		    ((cNod->keys[0] & dMsk) == 0) &&
		    ((cNod->keys[1] & dMsk) == 0) &&
		    ((cNod->keys[2] & dMsk) == 0);
	  keyCnt[idD] = allPrvD;
	  dMsk <<= 1;
	}
      }
      else
      {
	idD = 0;
	while(allPrvD && (idD < lDom->depth))
	{
	  /* Compute previous and current key digit for the current depth.
	   * These digits are in the range [0-3]. */
	  pDig =  ((pNod->keys[0] & dMsk) != 0) |
		 (((pNod->keys[1] & dMsk) != 0) << 1) |
		 (((pNod->keys[2] & dMsk) != 0) << 2);
	  cDig =  ((cNod->keys[0] & dMsk) != 0) |
		 (((cNod->keys[1] & dMsk) != 0) << 1) |
		 (((cNod->keys[2] & dMsk) != 0) << 2);
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
	      tNod->keys[2] &= ~dMsk;
	      tNod->keys[3] |= dMsk;
	      --tNod;
	    }
	  }
	  ++idD;
	  dMsk <<= 1;
	  dCnt *= 8;
	}
	while(idD < lDom->depth)
	{
	  keyCnt[idD++] = 0;
	}
      }
    }
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Marked nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
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
	 (cNod->keys[1] != pNod->keys[1]) ||
	 (cNod->keys[2] != pNod->keys[2]))
      {
	if(++idP != idC)
	{
	  pNod = lDom->nodes + idP;
	  pNod->keys[0] = cNod->keys[0];
	  pNod->keys[1] = cNod->keys[1];
	  pNod->keys[2] = cNod->keys[2];
	  pNod->keys[3] = cNod->keys[3];
	}
	++nCnt;
      }
      ++idC;
    }
    lDom->nNodes = nCnt;
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Condensed nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
  }
}

/*!
* \return	void
* \ingroup	WlzDomainOps
* \brief	Condenses the sorted 2D LBT domain's nodes.
* 		Runs through the nodes looking for keys in which the key
*		values differ only in less significant digits and the
*		less significant digits form a complete set, eg:
*		330 331 332 333. The these keys are modified using the
*		term code, eg: 33X.
*		When all such sequences have been marked identical keys
*		are removed leaving the linear quadtree nodes.
* \param	lDom			The LBT domain.
* \param	bndCndFlg		If non-zero then condensation of
*					boundary nodes is allowed.
*/
static void	 WlzLBTCondenseNodes2D(WlzLBTDomain2D *lDom,
				       int bndCndFlg)
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
#ifdef WLZ_LBTDOMAIN_DEBUG
  WlzDomain	tDom;
#endif
  WlzLBTNode2D	*cNod,
  		*pNod,
		*tNod;
  int		keyCnt[WLZ_LBTDOMAIN_MAXDIGITS];

  if(lDom->nNodes > 1)
  {
#ifdef WLZ_LBTDOMAIN_DEBUG
    tDom.l2 = lDom;
    (void )fprintf(stderr, "Nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
    /* Sort the nodes by key values. */
    AlgQSort(lDom->nodes, lDom->nNodes, sizeof(WlzLBTNode2D), lDom,
	     WlzLBTDomain2DNodeCmpFn);
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Sorted nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
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
      dMsk = 1;
      dCnt = 4;
      allPrvD = 1;
      pNod = cNod++;
      /* Avoid condensing nodes which are on the boundary of the domain. */
      if((bndCndFlg == 0) && ((pNod->flags & WLZ_LBT_NODE_FLAG_BOUNDARY) != 0))
      {
	for(idD = 0; idD < lDom->depth; ++idD)
	{
	  allPrvD = allPrvD &&
		    ((cNod->keys[0] & dMsk) == 0) &&
		    ((cNod->keys[1] & dMsk) == 0);
	  keyCnt[idD] = allPrvD;
	  dMsk <<= 1;
	}
      }
      else
      {
	idD = 0;
	while(allPrvD && (idD < lDom->depth))
	{
	  /* Compute previous and current key digit for the current depth.
	   * These digits are in the range [0-3]. */
	  pDig =  ((pNod->keys[0] & dMsk) != 0) |
		 (((pNod->keys[1] & dMsk) != 0) << 1);
	  cDig =  ((cNod->keys[0] & dMsk) != 0) |
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
	  keyCnt[idD++] = 0;
	}
      }
    }
#ifdef WLZ_LBTDOMAIN_DEBUG
    (void )fprintf(stderr, "Marked nodes:\n");
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
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
    (void )WlzLBTTestOutputNodesTxt(stderr, tDom);
#endif /* WLZ_LBTDOMAIN_DEBUG */
  }
}

#ifdef WLZ_UNUSED_FUNCTIONS
/*!
* \return	Non zero value if the node is at the domain's boundary.
* \ingroup	WlzDomainOps
* \brief	Checks whether the 3D LBT node is at the boundary of
*		the interval domain.
* \param	iDom
* \param	nod
*/
static int	WlzLBNodeAtBoundary3D(WlzPlaneDomain *pDom,
				      WlzLBTNode3D *nod)
{
  int		bnd = 0;
  WlzIVertex3	inc,
  		pos0,
  		pos1;

  WlzLBTKeyToPos3I(nod->keys, &pos0);
  for(inc.vtZ = -1; inc.vtZ <= 1; ++inc.vtZ)
  {
    pos1.vtZ = pos0.vtZ + inc.vtZ;
    for(inc.vtY = -1; inc.vtY <= 1; ++inc.vtY)
    {
      pos1.vtY = pos0.vtY + inc.vtY;
      for(inc.vtX = -1; inc.vtX <= 1; ++inc.vtX)
      {
	pos1.vtX = pos0.vtX + inc.vtX;
	if(WlzInsideDomain3D(pDom, pos1.vtZ, pos1.vtY, pos1.vtX, NULL) == 0)
	{
	  bnd = 1;
	  goto done;
	}
      }
    }
  }
done:
  return(bnd);
}

/*!
* \return	Non zero value if the node is at the domain's boundary.
* \ingroup	WlzDomainOps
* \brief	Checks whether the 2D LBT node is at the boundary of
*		the interval domain.
* \param	iDom
* \param	nod
*/
static int	WlzLBNodeAtBoundary2D(WlzIntervalDomain *iDom,
				      WlzLBTNode2D *nod)
{
  int		bnd = 0;
  WlzIVertex2	inc,
  		pos0,
  		pos1;

  WlzLBTKeyToPos2I(nod->keys, &pos0);
  for(inc.vtY = -1; inc.vtY <= 1; ++inc.vtY)
  {
    pos1.vtY = pos0.vtY + inc.vtY;
    for(inc.vtX = -1; inc.vtX <= 1; ++inc.vtX)
    {
      pos1.vtX = pos0.vtX + inc.vtX;
      if(WlzInsideDomain2D(iDom, pos1.vtY, pos1.vtX, NULL) == 0)
      {
        bnd = 1;
	goto done;
      }
    }
  }
done:
  return(bnd);
}
#endif /* WLZ_UNUSED_FUNCTIONS */

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Outputs the nodes of the LBT ask text for testing.
* \param	fP		Output file.
* \param	dom		The LBT domain to output.
*/
WlzErrorNum	WlzLBTTestOutputNodesTxt(FILE *fP, WlzDomain dom)
{
  int		idD,
  		idN;
  WlzBox	nBB;
  WlzLBTNode2D	*nod2;
  WlzLBTNode3D	*nod3;
  WlzUByte	digits[WLZ_LBTDOMAIN_MAXDIGITS];
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(dom.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(dom.core->type)
    {
      case WLZ_LBTDOMAIN_2D:
	(void )fprintf(fP, "type = %d\n", (int )(dom.l2->type));
	(void )fprintf(fP, "linkcount = %d\n", dom.l2->linkcount);
	(void )fprintf(fP, "freeptr = %p\n", dom.l2->freeptr);
	(void )fprintf(fP, "line1 %d\n", dom.l2->line1);
	(void )fprintf(fP, "lastln %d\n", dom.l2->lastln);
	(void )fprintf(fP, "kol1 %d\n", dom.l2->kol1);
	(void )fprintf(fP, "lastkl %d\n", dom.l2->lastkl);
	(void )fprintf(fP, "depth %d\n", dom.l2->depth);
	(void )fprintf(fP, "nNodes %d\n", dom.l2->nNodes);
	(void )fprintf(fP, "maxNodes %d\n", dom.l2->maxNodes);
	(void )fprintf(fP, "nodes:\n");
	nod2 = dom.l2->nodes;
	idN = 0;
	while(idN < dom.l2->nNodes)
	{
	  (void )fprintf(fP, "%8d  ", idN);
	  WlzLBTGetKeyDigits2D(nod2->keys, digits);
	  for(idD = 29; idD >= 0; --idD)
	  {
	    (void )fprintf(fP, "%d", digits[idD]);
	  }
	  WlzLBTKeyToBox2I(nod2->keys, &(nBB.i2));
	  (void )fprintf(fP, "  %d,%d,%d,%d",
			 nBB.i2.xMin, nBB.i2.yMin, nBB.i2.xMax, nBB.i2.yMax);
	  (void )fprintf(fP, "  %d\n",
			 (nBB.i2.xMax - nBB.i2.xMin + 1) *
			 (nBB.i2.yMax - nBB.i2.yMin + 1));
	  ++idN;
	  ++nod2;
	}
        break;
      case WLZ_LBTDOMAIN_3D:
	(void )fprintf(fP, "type = %d\n", (int )(dom.l3->type));
	(void )fprintf(fP, "linkcount = %d\n", dom.l3->linkcount);
	(void )fprintf(fP, "freeptr = %pn", dom.l3->freeptr);
	(void )fprintf(fP, "plane1 %d\n", dom.l3->plane1);
	(void )fprintf(fP, "lastpl %d\n", dom.l3->lastpl);
	(void )fprintf(fP, "line1 %d\n", dom.l3->line1);
	(void )fprintf(fP, "lastln %d\n", dom.l3->lastln);
	(void )fprintf(fP, "kol1 %d\n", dom.l3->kol1);
	(void )fprintf(fP, "lastkl %d\n", dom.l3->lastkl);
	(void )fprintf(fP, "depth %d\n", dom.l3->depth);
	(void )fprintf(fP, "nNodes %d\n", dom.l3->nNodes);
	(void )fprintf(fP, "maxNodes %d\n", dom.l3->maxNodes);
	(void )fprintf(fP, "nodes:\n");
	nod3 = dom.l3->nodes;
	idN = 0;
	while(idN < dom.l3->nNodes)
	{
	  (void )fprintf(fP, "%8d  ", idN);
	  WlzLBTGetKeyDigits3D(nod3->keys, digits);
	  for(idD = 29; idD >= 0; --idD)
	  {
	    (void )fprintf(fP, "%d", digits[idD]);
	  }
	  WlzLBTKeyToBox3I(nod3->keys, &(nBB.i3));
	  (void )fprintf(fP, "  %d,%d,%d,%d,%d,%d",
			 nBB.i3.xMin, nBB.i3.yMin, nBB.i3.zMin,
			 nBB.i3.xMax, nBB.i3.yMax, nBB.i3.zMax);
	  (void )fprintf(fP, "  %d\n",
			 (nBB.i3.xMax - nBB.i3.xMin + 1) *
			 (nBB.i3.yMax - nBB.i3.yMin + 1) *
			 (nBB.i3.zMax - nBB.i3.zMin + 1));
	  ++idN;
	  ++nod3;
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
* \ingroup	WlzDomainOps
* \brief	Outputs the nodes of the LBT as VTK polydata for testing.
* \param	fP		Output file.
* \param	dom		The LBT domain to output.
*/
WlzErrorNum	WlzLBTTestOutputNodesVtk(FILE *fP, WlzDomain dom)
{
  int		idN,
  		nOff;
  WlzLBTNode2D	*nod2;
  WlzLBTNode3D	*nod3;
  WlzBox	nodBB;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(dom.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(dom.core->type)
    {
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      case WLZ_LBTDOMAIN_2D:
	if(dom.l2->nNodes > 0)
	{
	  (void )fprintf(fP, "# vtk DataFile Version 1.0\n"
			     "WlzLBTDomain 2D test output\n"
			     "ASCII\n"
			     "DATASET POLYDATA\n"
			     "POINTS %d float\n",
			     dom.l2->nNodes * 4);
	  nod2 = dom.l2->nodes;
	  for(idN = 0; idN < dom.l2->nNodes; ++idN)
	  {
	    WlzLBTKeyToBox2I(nod2->keys, &(nodBB.i2));
	    nodBB.i2.xMax += 1;
	    nodBB.i2.yMax += 1;
	    (void )fprintf(fP, "%d %d 0\n"
			       "%d %d 0\n"
			       "%d %d 0\n"
			       "%d %d 0\n",
			       nodBB.i2.xMin, nodBB.i2.yMin,
			       nodBB.i2.xMax, nodBB.i2.yMin,
			       nodBB.i2.xMax, nodBB.i2.yMax,
			       nodBB.i2.xMin, nodBB.i2.yMax);
	    ++nod2;
	  }
	  (void )fprintf(fP, "LINES %d %d\n",
			     dom.l2->nNodes * 4, dom.l2->nNodes * 12);
	  for(idN = 0; idN < dom.l2->nNodes; ++idN)
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
	}
	break;
      case WLZ_LBTDOMAIN_3D:
	if(dom.l3->nNodes > 0)
	{
	  (void )fprintf(fP, "# vtk DataFile Version 1.0\n"
			 "WlzLBTDomain 3D test output\n"
			 "ASCII\n"
			 "DATASET POLYDATA\n"
			 "POINTS %d float\n",
			 dom.l3->nNodes * 8);
	  nod3 = dom.l3->nodes;
	  for(idN = 0; idN < dom.l3->nNodes; ++idN)
	  {
	    WlzLBTKeyToBox3I(nod3->keys, &(nodBB.i3));
	    nodBB.i3.xMax += 1;
	    nodBB.i3.yMax += 1;
	    nodBB.i3.zMax += 1;
	    (void )fprintf(fP, "%d %d %d\n"
			       "%d %d %d\n"
			       "%d %d %d\n"
			       "%d %d %d\n"
			       "%d %d %d\n"
			       "%d %d %d\n"
			       "%d %d %d\n"
			       "%d %d %d\n",
		     nodBB.i3.xMin, nodBB.i3.yMin, nodBB.i3.zMin,
		     nodBB.i3.xMax, nodBB.i3.yMin, nodBB.i3.zMin,
		     nodBB.i3.xMax, nodBB.i3.yMax, nodBB.i3.zMin,
		     nodBB.i3.xMin, nodBB.i3.yMax, nodBB.i3.zMin,
		     nodBB.i3.xMin, nodBB.i3.yMin, nodBB.i3.zMax,
		     nodBB.i3.xMax, nodBB.i3.yMin, nodBB.i3.zMax,
		     nodBB.i3.xMax, nodBB.i3.yMax, nodBB.i3.zMax,
		     nodBB.i3.xMin, nodBB.i3.yMax, nodBB.i3.zMax);
	    ++nod3;
	  }
	  (void )fprintf(fP, "POLYGONS %d %d\n",
			 dom.l3->nNodes * 6, dom.l3->nNodes * 30);
	  for(idN = 0; idN < dom.l3->nNodes; ++idN)
	  {
	    nOff = idN * 8;
	    (void )fprintf(fP, "4 %d %d %d %d\n"
			       "4 %d %d %d %d\n"
			       "4 %d %d %d %d\n"
			       "4 %d %d %d %d\n"
			       "4 %d %d %d %d\n"
			       "4 %d %d %d %d\n",
			       nOff + 0, nOff + 1, nOff + 2, nOff + 3,
			       nOff + 1, nOff + 0, nOff + 4, nOff + 5,
			       nOff + 2, nOff + 1, nOff + 5, nOff + 6,
			       nOff + 3, nOff + 2, nOff + 6, nOff + 7,
			       nOff + 0, nOff + 3, nOff + 7, nOff + 4,
			       nOff + 7, nOff + 6, nOff + 5, nOff + 4);
	  }
	}
        break;
    }
  }
  return(errNum);
}

#ifdef WLZ_LBTDOMAIN_TEST_1

typedef enum _WlzLBTTestOutType
{
  WLZ_LBTDOMAIN_TEST_OUT_IDOM,
  WLZ_LBTDOMAIN_TEST_OUT_IDX,
  WLZ_LBTDOMAIN_TEST_OUT_TXT,
  WLZ_LBTDOMAIN_TEST_OUT_VTK
} WlzLBTTestOutType;

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		dim = 0,
  		balance = 0,
  		option,
  		ok = 1,
		maxNodSz = INT_MAX,
		maxBndNodSz = INT_MAX,
		usage = 0;
  WlzLBTTestOutType outType = WLZ_LBTDOMAIN_TEST_OUT_VTK;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  FILE		*fP = NULL;
  WlzObject	*idxObj = NULL,
  		*inObj = NULL,
  		*outObj = NULL;
  WlzDomain	lDom,
  		outDom;
  WlzValues	nullVal;
  char		*inFileStr,
  		*outFileStr;
  static char	optList[] = "bditvho:",
  		outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  lDom.core = NULL;
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
      case 'd':
        outType = WLZ_LBTDOMAIN_TEST_OUT_IDOM;
	break;
      case 'i':
        outType = WLZ_LBTDOMAIN_TEST_OUT_IDX;
	break;
      case 't':
        outType = WLZ_LBTDOMAIN_TEST_OUT_TXT;
        break;
      case 'v':
	outType = WLZ_LBTDOMAIN_TEST_OUT_VTK;
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
    else if(inObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      switch(inObj->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  dim = 2;
	  break;
        case WLZ_3D_DOMAINOBJ:
	  dim =3;
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
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
    lDom = WlzLBTDomainFromObj(inObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to compute LBT domain (%d).\n",
		     *argv, errMsg);
    }
  }
  if(ok && (balance || (outType == WLZ_LBTDOMAIN_TEST_OUT_IDX)))
  {
    switch(lDom.core->type)
    {
      case WLZ_LBTDOMAIN_2D:
        idxObj = WlzLBTMakeNodeIndexObj2D(lDom.l2, inObj->domain.i, &errNum);
        if(errNum == WLZ_ERR_NONE)
        {
	   errNum = WlzLBTBalanceDomain2D(lDom.l2, idxObj, maxNodSz,
	                                  maxBndNodSz);
        }
	break;
      case WLZ_LBTDOMAIN_3D:
        idxObj = WlzLBTMakeNodeIndexObj3D(lDom.l3, inObj->domain.p, &errNum);
        if(errNum == WLZ_ERR_NONE)
        {
	   errNum = WlzLBTBalanceDomain3D(lDom.l3, idxObj, maxNodSz,
	                                  maxBndNodSz);
        }
	break;
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
      case WLZ_LBTDOMAIN_TEST_OUT_IDX:
        errNum = WlzWriteObj(fP, idxObj);
	break;
      case WLZ_LBTDOMAIN_TEST_OUT_IDOM:
        switch(lDom.core->type)
	{
          case WLZ_LBTDOMAIN_2D:
            outDom.i = WlzLBTDomainToIDomain(lDom.l2, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      outObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, outDom, nullVal,
				   NULL, NULL, &errNum);
	    }
	    break;
          case WLZ_LBTDOMAIN_3D:
            outDom.p = WlzLBTDomainToPDomain(lDom.l3, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      outObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, outDom, nullVal,
				   NULL, NULL, &errNum);
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
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
  (void )WlzFreeObj(idxObj);
  (void )WlzFreeLBTDomain2D(lDom.l2);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<output object>] [-h] [-b] [-i|t|v] [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output file name.\n"
    "  -b  Produce a balanced tree.\n"
    "  -d  Woolz interval or plane domain output.\n"
    "  -i  Index object output.\n"
    "  -t  Text output.\n"
    "  -v  VTK output.\n");
  }
  return(!ok);
}
#endif /* WLZ_LBTDOMAIN_TEST_1 */
