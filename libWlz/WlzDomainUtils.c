#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzDomainUtils.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Utility functions for Woolz domains.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzStandardIntervalDomain				*
*   Date       : Mon Oct 21 18:38:20 1996				*
*************************************************************************
*   Synopsis   :standardise an interval domain - minimal bounding box,	*
*		strip leading and trailing empty lines. The domain	*
*		is modified "in place".					*
*   Returns    :WlzErrorNum: error return values: WLZ_ERR_NONE,		*
*		WLZ_ERR_DOMAIN_NULL						*
*   Parameters :WlzIntervalDomain *idom: the domain to be standardised	*
*   Global refs:None							*
************************************************************************/

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

  /* find first non-empty line */
  while (idom->line1 < idom->lastln && idom->intvlines->nintvs == 0) {
    idom->line1++;
    idom->intvlines++;
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

/************************************************************************
*   Function   : WlzStandardPlaneDomain					*
*   Date       : Mon Oct 21 17:40:00 1996				*
*************************************************************************
*   Synopsis   :Standardize a plane domain and corresponding voxel-table*
*		(voxel-tables must have exactly matching valuetables)	*
*		by stripping leading and trailing NULL domains and	*
*		standardising each domain in turn. The bounding box is	*
*		reset to be minimal.					*
*   Returns    :WlzErrorNum: error return values: WLZ_ERR_DOMAIN_NULL,		*
*		WLZ_ERR_NONE, WLZ_ERR_PLANEDOMAIN_DATA. 				*
*   Parameters :WlzPlaneDomain 	*pdom: the planedomain must be non-NULL	*
*		WlzVoxelValues	*voxtb: corresponding voxel table - may	*
*		br NULL.						*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzStandardPlaneDomain(WlzPlaneDomain 	*pdom,
				   WlzVoxelValues	*voxtb)
{
  /* local variables */
  WlzObject 	tempobj;
  WlzDomain 	*domains;
  int	 	line1, lastln, kol1, lastkl;
  int 		p, nplanes, np, firstplane, lastplane;

  /* check domain argument */
  if( pdom == NULL ){
    return( WLZ_ERR_DOMAIN_NULL );
  }
  
  /* set local variables */
  nplanes = pdom->lastpl - pdom->plane1 + 1;
  domains = pdom->domains;

  /* check for at least one intervaldomain
     the loop exits with the plane index of the first non-NULL domain.
     np is set to the value of the first non-NULL index.
     If this value is >= nplanes then all domains are NULL. */
  for (p=0; p < nplanes; p++, domains++){
    if( (*domains).i != NULL ){
      WlzStandardIntervalDomain((*domains).i);
      line1 = (*domains).i->line1;
      lastln = (*domains).i->lastln;
      kol1 = (*domains).i->kol1;
      lastkl = (*domains).i->lastkl;
      break;
    }
  }
  if( (np = p) >= nplanes ){
    return( WLZ_ERR_NONE );
  }

  /* check line and column limits */
  for(p=np; p < nplanes; p++, domains++){
    if ( (*domains).i != NULL ){
      WlzStandardIntervalDomain((*domains).i);
      if( (*domains).i->line1 < line1 )  { line1  = (*domains).i->line1;}
      if( (*domains).i->lastln > lastln ){ lastln = (*domains).i->lastln;}
      if( (*domains).i->kol1 < kol1 )    { kol1   = (*domains).i->kol1;}
      if( (*domains).i->lastkl > lastkl ){ lastkl = (*domains).i->lastkl;}
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
    if( pdom->domains[p].i != NULL && WlzArea(&tempobj, NULL) > 0 ){
      lastplane = p;
      np++;
    }
  }
  for( p = nplanes - 1, firstplane = p; p >=0; p--){
    tempobj.domain.i = pdom->domains[p].i;
    if( pdom->domains[p].i != NULL && WlzArea(&tempobj, NULL) > 0  )
      firstplane = p;
  }
  if( np == NULL ){
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
