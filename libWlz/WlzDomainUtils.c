#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzDomainUtils.c
* \author       Bill Hill, Richard Baldock
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Utility functions for Woolz domains.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

/*!
* \return	The type of the domain, WLZ_NULL if the domain is NULL.
* \ingroup	WlzDomainOps
* \brief	Gets the type of the given domain.
* \param	dom			Given domain.
*/
WlzObjectType	WlzDomainType(WlzDomain dom)
{
  WlzObjectType	type = WLZ_NULL;

  if(dom.core)
  {
    type = dom.core->type;
  }
  return(type);
}

/*!
* \return	<void>
* \ingroup	WlzDomainOps
* \brief	Sets bits in the given byte packed bit line which are within
* 		the given interval.
* \param	bitLnP			The bit line pointer.
* \param	iLft			Left coordinate of interval.
* \param	iRgt			Right coordinate of interval.
* \param	size			Number of bits in the bit line.
*/
void		WlzBitLnSetItv(UBYTE *bitLnP, int iLft, int iRgt, int size)
{
  int		bitIdx,
  		bitCnt,
		bytIdx,
		bytCnt,
		mskIdx,
		mskCnt;
  UBYTE		*bytP;
  const UBYTE	bitMsk[9] = 
  {
    0x00, 
    0x01,
    0x03,
    0x07,
    0x0f,
    0x1f,
    0x3f,
    0x7f,
    0xff
  };

  if((iLft < size) && (iRgt >= 0))
  {
    if(iLft < 0)
    {
      iLft = 0;
    }
    if(iRgt >= size)
    {
      iRgt = size - 1;
    }
    /* Set bits in this interval. */
    bitIdx = iLft;
    bitCnt = iRgt - iLft + 1;
    bytIdx = bitIdx >> 3;
    mskIdx = bitIdx & 7;
    mskCnt = bitCnt;
    if(mskCnt > 8)
    {
      mskCnt = 8;
    }
    bytP = bitLnP + bytIdx;
    /* Set bits in first byte. */
    *bytP++ |= 0xff & (bitMsk[mskCnt] << mskIdx);
    bitCnt -= 8 - mskIdx;
    if(bitCnt > 0)
    {
      /* Set bits in all mid interval bytes. */
      bytCnt = bitCnt >> 3;
      while(bytCnt-- > 0)
      {
	*bytP++ = 0xff;
      }
      /* Set bits in trailing byte. */
      bitCnt &= 7;
      if(bitCnt)
      {
	*bytP |= bitMsk[bitCnt];
      }
    }
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Adds an interval line to a interval domain given
*               an allocated interval domain a byte packed bitmask
*               for the interval line and a pool of available
*               intervals.
* \param	iDom			Given interval domain.
* \param	bitLn			Byte packed bitmap for line.
* \param	line			The line coordinate.
* \param	width			Width of the line, ie number of valid
* 					bits in bitmask.
* \param	iPool			Interval pool.
*/
WlzErrorNum 	WlzDynItvLnFromBitLn(WlzIntervalDomain *iDom,
				     UBYTE *bitLn, int line, int width,
				     WlzDynItvPool *iPool)
{
  int		iLft = 0,
  		iLen = 0,
		inDomFlg,
		bitIdx,
		bitMsk,
		bitIdxWas;
  UBYTE		*bytP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iDom == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(iDom->type != WLZ_INTERVALDOMAIN_INTVL)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((iPool == NULL) || (bitLn == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    bitMsk = 1;
    bitIdx = 0;
    inDomFlg = 0;
    bytP = bitLn;
    while((errNum == WLZ_ERR_NONE) && (bitIdx < width))
    {
      if(*bytP == 0)
      {
	if(inDomFlg)
	{
	  errNum = WlzDynItvAdd(iDom, iPool, line, iLft, iLen);
	  inDomFlg = 0;
	  iLen = 0;
	}
	do
	{
	  ++bytP;
	  bitIdx += 8;
	}
	while((*bytP == 0) && (bitIdx < width));
	if(bitIdx > width)
	{
	  bitIdx = width;
	}
      }
      else if(*bytP == 0xff)
      {
	bitIdxWas = bitIdx;
	if(inDomFlg == 0)
	{
	  inDomFlg = 1;
	  iLft = bitIdx;
	}
        do
	{
	  ++bytP;
	  bitIdx += 8;
	}
	while((*bytP == 0xff) && (bitIdx < width));
	if(bitIdx > width)
	{
	  bitIdx = width;
	}
	iLen += bitIdx - bitIdxWas;
      }
      else
      {
	if(inDomFlg)
	{
	  if((inDomFlg = *bytP & bitMsk) != 0)
	  {
	    ++iLen;
	  }
	  else
	  {
	    errNum = WlzDynItvAdd(iDom, iPool, line, iLft, iLen);
	    iLen = 0;
	  }
	}
	else
	{
	  if((inDomFlg = *bytP & bitMsk) != 0)
	  {
	    iLft = bitIdx;
	    iLen = 1;
	  }
	}
	if((bitMsk = bitMsk << 1) > (1<<7))
	{
	  bitMsk = 1;
	  ++bytP;
	}
	++bitIdx;
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && inDomFlg)
  {
    errNum = WlzDynItvAdd(iDom, iPool, line, iLft, iLen);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Adds an interval to a interval domain given an allocated
* 		interval domain, a pool of available intervals, the line, the
* 		intervals left most column and the intertvals width.
* \param	iDom			Given interval domain.
* \param	iPool			Interval pool.
* \param	line			The line.
* \param	iLft			Left most column of interval.
* \param	iLen			Width of the interval.
*/
WlzErrorNum 	WlzDynItvAdd(WlzIntervalDomain *iDom, WlzDynItvPool *iPool,
			     int line, int iLft, int iLen)
{
  WlzInterval	*itv;
  WlzIntervalLine *itvLn;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;
  AlcErrno	alcErr = ALC_ER_NONE;
#ifdef WLZ_DYNITV_TUNE_MALLOC
  static int 	tuneMalloc = 0;
#endif /* WLZ_DYNITV_TUNE_MALLOC */

  if(iDom == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(iDom->type != WLZ_INTERVALDOMAIN_INTVL)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(line < iDom->line1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(iPool == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    iDom->lastln = line;
    itvLn = iDom->intvlines + line - iDom->line1;
    if((iPool->itvBlock == NULL) ||
       (iPool->itvsInBlock - iPool->offset - 1) < 1)
    {
#ifdef WLZ_DYNITV_TUNE_MALLOC
      (void )fprintf(stderr, "WlzDynItvAdd(tuneMalloc): %d\n", ++tuneMalloc);
#endif /* WLZ_DYNITV_TUNE_MALLOC */
      iPool->offset = 0;
      if((iPool->itvBlock = AlcMalloc(iPool->itvsInBlock *
				     sizeof(WlzInterval))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	iDom->freeptr = AlcFreeStackPush(iDom->freeptr,
					 (void *)(iPool->itvBlock),
					 &alcErr);
        if(alcErr != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(itvLn->nintvs > 0)
	{
	  (void )memcpy((void *)(iPool->itvBlock), (void *)(itvLn->intvs),
			itvLn->nintvs * sizeof(WlzInterval));
	}
	itvLn->intvs = iPool->itvBlock;
	iPool->offset = itvLn->nintvs;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(itvLn->intvs == NULL)
    {
      itvLn->intvs = (WlzInterval *)(iPool->itvBlock) + iPool->offset;
    }
    if(itvLn->intvs)
    {
      itv = itvLn->intvs + itvLn->nintvs;
      itv->ileft = iLft;
      itv->iright = iLft + iLen - 1;
      ++(itvLn->nintvs);
      ++(iPool->offset);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Standardises an interval domain by ensuring that it
*		has a minaimal bounding box, striping away empty lines
*		as required. The domain is modified "in place".
* \param	idom			The given interval domain.
*/
WlzErrorNum WlzStandardIntervalDomain(WlzIntervalDomain *idom)
{
  WlzIntervalLine	*itvl;
  WlzInterval		*intv;
  int			i,j,k,r;

  /* check the domain */
  if( idom == NULL ){
    return WLZ_ERR_DOMAIN_NULL;
  }

  /* check domain type */
  if( idom->type != WLZ_INTERVALDOMAIN_INTVL ){
    return WLZ_ERR_NONE;
  }

  /* check interval line is non NULL */
  if( idom->intvlines == NULL ){
    return WLZ_ERR_INTERVALLINE_NULL;
  }

  /* find first non-empty line */
  while (idom->line1 < idom->lastln && idom->intvlines->nintvs == 0) {
    idom->line1++;
    idom->intvlines++;
  }

  /* check interval line has no intervals */
  if( idom->intvlines->nintvs == 0 ){
    return WLZ_ERR_EOO;
  }

  /* find last non-empy line */
  while (idom->intvlines[idom->lastln-idom->line1].nintvs == 0) {
    if (idom->line1 == idom->lastln) {
      idom->kol1 = idom->lastkl = 0;
      return( WLZ_ERR_NONE );
    }
    idom->lastln--;
  }

  /* find left-most end point of non-empty lines */
  itvl = idom->intvlines;
  k = itvl->intvs->ileft;	/* first line IS now non-empty */
  for (i=idom->line1; i<=idom->lastln; i++) {
    intv = itvl->intvs;
    for (j=0; j<itvl->nintvs; j++) {
      if (intv->ileft < k)
	k = intv->ileft;
      intv++;
    }
    itvl++;
  }

  /* update interval domain and interval end points, find right-most */
  idom->kol1 += k;
  r = 0;
  itvl = idom->intvlines;
  for (i=idom->line1; i<=idom->lastln; i++) {
    intv = itvl->intvs;
    for (j=0; j<itvl->nintvs; j++) {
      intv->ileft -= k;
      intv->iright -= k;
      if (intv->iright > r)
	r = intv->iright;
      intv++;
    }
    itvl++;
  }
  idom->lastkl = idom->kol1 + r;

  return( WLZ_ERR_NONE );
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Standardizes a plane domain and corresponding voxel-table
* 		(voxel-tables must have exactly matching valuetables) by
* 		stripping leading and trailing NULL domains and standardising
* 		each domain in turn. The bounding box is reset to be minimal.
*		Both the domain and values may be modified by this function.
*		Any plane which has an invalid interval domain is discarded.
* \param	pdom			Given plane domain, must NOT be NULL..
* \param	voxtb			Corresponding voxel table, maybe NULL.
*/
WlzErrorNum WlzStandardPlaneDomain(WlzPlaneDomain 	*pdom,
				   WlzVoxelValues	*voxtb)
{
  /* local variables */
  WlzObject 	tempobj;
  WlzDomain 	*domains,
  		*tDom;
  WlzValues	*tVal;
  int	 	line1, lastln, kol1, lastkl;
  int 		p, nplanes, np, firstplane, lastplane, bndFnd;

  /* check domain argument */
  if( pdom == NULL ){
    return( WLZ_ERR_DOMAIN_NULL );
  }
  
  /* set local variables */
  nplanes = pdom->lastpl - pdom->plane1 + 1;
  domains = pdom->domains;

  /* check line and column limits */
  bndFnd = 0;
  for(p = 0; p < nplanes; ++p, ++domains)
  {
    if((*domains).core)
    {
      if(WlzStandardIntervalDomain((*domains).i) == WLZ_ERR_NONE)
      {
	if(bndFnd == 0)
	{
	  np = p;
	  bndFnd = 1;
	  line1 = (*domains).i->line1;
	  lastln = (*domains).i->lastln;
	  kol1 = (*domains).i->kol1;
	  lastkl = (*domains).i->lastkl;
	}
	else
	{
	  if((*domains).i->line1 < line1)
	  {
	    line1 = (*domains).i->line1;
	  }
	  if((*domains).i->lastln > lastln)
	  {
	    lastln = (*domains).i->lastln;
	  }
	  if((*domains).i->kol1 < kol1)
	  {
	    kol1 = (*domains).i->kol1;
	  }
	  if((*domains).i->lastkl > lastkl)
	  {
	    lastkl = (*domains).i->lastkl;
	  }
	}
      }
      else
      {
        /* There's something wrong with this plane! Throw it away. */
	WlzFreeIntervalDomain(domains->i);
	domains->core = NULL;
	if(voxtb)
	{
	  WlzFreeValueTb((voxtb->values + p)->v);
	  (voxtb->values + p)->core = NULL;
	}
      }
    }
  }

  pdom->line1 = line1;
  pdom->lastln = lastln;
  pdom->kol1 = kol1;
  pdom->lastkl = lastkl;
  domains -= nplanes;

  /* strip leading and trailing NULL or empty domains */
  tempobj.type = WLZ_2D_DOMAINOBJ;
  tempobj.domain.i = NULL;
  tempobj.values.v = NULL;
  tempobj.plist = NULL;
  tempobj.assoc = NULL;

  /* find the last non-NULL or non-empty plane, np is the number
     of such planes */
  lastplane = 0;
  for(p=0, np=0; p < nplanes; p++){
    tempobj.domain.i = pdom->domains[p].i;
    if( pdom->domains[p].i != NULL && !WlzIsEmpty(&tempobj, NULL) ){
      lastplane = p;
      np++;
    }
  }
  for( p = nplanes - 1, firstplane = p; p >=0; p--){
    tempobj.domain.i = pdom->domains[p].i;
    if( pdom->domains[p].i != NULL && !WlzIsEmpty(&tempobj, NULL)  )
      firstplane = p;
  }
  if( np == 0 ){
    firstplane = 0;
    lastplane = 0;
  }
  if( firstplane > lastplane ){
    return( WLZ_ERR_PLANEDOMAIN_DATA );
  }

  /* free interval and valuedomains that become dereferenced */
  for(p=0; p < firstplane; p++){
    if( pdom->domains[p].i )
      WlzFreeIntervalDomain( pdom->domains[p].i );
    if( voxtb && voxtb->values[p].v )
      WlzFreeValueTb( voxtb->values[p].v );
  }
  for(p=lastplane+1; p < nplanes; p++){
    if( pdom->domains[p].i )
      WlzFreeIntervalDomain( pdom->domains[p].i );
    if( voxtb && voxtb->values[p].v )
      WlzFreeValueTb( voxtb->values[p].v );
  }

  /* shift the remaining planes */
  for(p=firstplane, np=0; p <= lastplane; p++, np++){
    pdom->domains[np] = pdom->domains[p];
    if( voxtb )
      voxtb->values[np] = voxtb->values[p];
  }
  pdom->plane1 += firstplane;
  pdom->lastpl  = pdom->plane1 + lastplane - firstplane;
  if( voxtb )
    {
      voxtb->plane1 = pdom->plane1;
      voxtb->lastpl = pdom->lastpl;
    }

  /* return */
  return( WLZ_ERR_NONE );
}
