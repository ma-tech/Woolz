#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzLabel.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Tue Aug 19 17:47:51 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzBinaryOps
* \brief        Segment a given woolz object into disconnected regions.
*               
* \todo         Make object array grow as required when requested, e.g.
 for maxNumObjs=0.
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
* 03-03-2K bill	Replace WlzPushFreePtr(), WlzPopFreePtr() and 
*		WlzFreeFreePtr() with AlcFreeStackPush(),
*		AlcFreeStackPop() and AlcFreeStackFree().
*/

#include <stdlib.h>
#include <Wlz.h>

/*
 * The allocation list contains the interval end-points of the intervals of
 * the current and previous lines.  Its length therefore needs to be twice the
 * maximum number of intervals in a line of the input object (see mintcount()).
 */
typedef struct {
  struct _WlzLLink	*a_link;
  WlzInterval		a_int;
} WlzLAllocBuf;

/*
 * A chain-link is used to hold either data (line number, interval
 * end points) or a pointer to another chain.  Do not confuse this
 * latter use with the pointer which points to the next link in
 * the same chain.
 */
typedef struct _WlzLLink {
  struct _WlzLLink *l_link;
  union {
    unsigned long line;
    struct _WlzLLink *u_link;
    WlzInterval intv;
  } l_u;
} WlzLLink;

static WlzLLink	*newchain1	(WlzLLink *chainbase,
				 unsigned long entry);
static WlzLLink *newchain2	(WlzLLink *chainbase,
				 int entry1,
				 int entry2);
static WlzLLink *addlink1	(WlzLLink *chainbase,
				 WlzLLink *chain,
				 unsigned long entry);
static WlzLLink *addlink2	(WlzLLink *chainbase,
				 WlzLLink *chain,
				 int entry1,
				 int entry2);
static WlzLLink *chainalloc	(int flag, int n);
static void 	chainFree	(void);
static WlzLLink *getlink	(WlzLLink *chainbase);
static int 	buckle 		(WlzLLink *baselink,
				 WlzLLink *chainbase,
				 int length);
static int	freechain	(WlzLLink *l);
static int 	join		(WlzLLink *list1,
				 WlzLLink *list2);
static int 	merge		(WlzLLink **ak1,
				 WlzLLink **ak2,
				 WlzLLink *freechain);
static int 	monotonejoin	(WlzLLink **amlist1,
				 WlzLLink *mlist2);
static int 	mintcount	(WlzIntervalDomain *idom,
				 int *maxinline);

extern WlzErrorNum	WlzLabel3d	(WlzObject	*obj,
					 int		*numobj,
					 WlzObject	**objlist,
					 int		nobj,
					 int		ignlns,
					 WlzConnectType connect);

/* function:     WlzLabel    */
/*! 
* \ingroup      WlzBinaryOps
* \brief        Segment a domain into connected parts. Connectivity
 is defined by the connect parameter and can be 4- or 8-connected for
 2D objects and 6-, 18- or 26-connected for 3D objects. Note this
 version requires that there is sufficient space in the objects array
 defined by maxNumObjs and this is not extended. This should be changed
 in future so that the array is extended as required.
*
* \return       Error number.
* \param    obj	input object to be segmented
* \param    mm	number of objects return
* \param    dstArrayObjs	object array return, allocated in the
 procedure.
* \param    maxNumObjs	maximum number of object to return (determines
 the size of the array)
* \param    ignlns	ignore objects with num lines <= ignlns
* \param    connect	connectivity to determine connected regions
* \par      Source:
*                WlzLabel.c
*/
WlzErrorNum WlzLabel(
  WlzObject	*obj,
  int		*mm,
  WlzObject	***dstArrayObjs,
  int		maxNumObjs,
  int		ignlns,
  WlzConnectType connect)
{ 
  WlzIntervalDomain 	*jdp;
  WlzRagRValues 	*jvp;
  WlzIntervalWSpace 	iwsp;
  WlzInterval 		*itvl;
  WlzInterval 		*jtvl;
  WlzLAllocBuf 		*crntal, *altemp;
  WlzLAllocBuf		*al, *crntst, *crlast;
  WlzLAllocBuf		*precal, *precst, *prlast;
  int			nints, mkl;
  int			maxinline,nob,line,ended,lend,chainlistsize;
  WlzLLink		*freechain, *alprec, *alloc, *link1, *link2;
  int			lftcrn, lftprc, rtcrn, rtprec;
  int			jrtcrn, mxkl, jl, jr;
  int			oll, ofl, jjj;
  WlzDomain		domain;
  WlzValues		values;
  WlzConnectType	jdqt;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzObject		**objlist;

  /* see HISTORY for comments from the FORTRAN version */

  /* now we allocate space for the objects */
  if( (objlist = (WlzObject **)AlcMalloc(sizeof(WlzObject *) * maxNumObjs))
     == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else {
    *dstArrayObjs = objlist;
  }
  

  /* check object note *mm is always set to zero on error
     return because the "Too many objects" error can return
     nobj valid objects therefore if *mm != 0 there are valid objects
     in objlist which must be freed */
  if( obj == NULL ){
    *mm = 0;
    return WLZ_ERR_OBJECT_NULL;
  }

  /* check types */
  switch( obj->type ){

  case WLZ_2D_DOMAINOBJ:
    if( obj->domain.core == NULL ){
      *mm = 0;
      return WLZ_ERR_DOMAIN_NULL;
    }

    switch( obj->domain.core->type ){

    case WLZ_INTERVALDOMAIN_INTVL:
      break;

    case WLZ_INTERVALDOMAIN_RECT:
      if( (obj->domain.i->lastln - obj->domain.i->line1) < ignlns ){
	*mm = 0;
	return( WLZ_ERR_NONE );
      }

      if( maxNumObjs < 1 ){
	*mm = 0;
	return( WLZ_ERR_INT_DATA );
      }

      objlist[0] = WlzAssignObject(
	WlzMakeMain(obj->type, obj->domain,
		    obj->values, NULL, NULL, &errNum), NULL);
      *mm = 1;
      return WLZ_ERR_NONE;
      
    default:
      *mm = 0;
      return WLZ_ERR_DOMAIN_TYPE;

    }
    break;

  case WLZ_3D_DOMAINOBJ:
    if( obj->domain.core == NULL ){
      *mm = 0;
      return WLZ_ERR_DOMAIN_NULL;
    }

    return WlzLabel3d(obj, mm, objlist, maxNumObjs, ignlns, connect);

  case WLZ_EMPTY_OBJ:
    *mm = 0;
    return WLZ_ERR_NONE;

  case WLZ_TRANS_OBJ:
    /* should be able to do this quite easily */
  default:
    *mm = 0;
    return WLZ_ERR_OBJECT_TYPE;

  }

  /* check connectivity parameter */
  switch( connect ){

  default:
    return WLZ_ERR_INT_DATA;

  case WLZ_8_CONNECTED:
    jdqt = 0;
    break;

  case WLZ_4_CONNECTED:
    jdqt = 1;
    break;

  }
    
  /*
   * Allocate and initialise the working spaces.
   *
   * The chain store size is not easily predictable.
   * Here we use ((twice the number of lines) plus (one third the
   * number of intervals)), very much smaller than the
   * previous allocation, but apparently adequate to
   * segment both metaphase spreads and single connected objects.
   * In any case, this version of label can if necessary allocate
   * more space when the original chain store fills up.
   */
  chainlistsize = mintcount(obj->domain.i, &maxinline) / 3;
  chainlistsize += (1 + obj->domain.i->lastln - obj->domain.i->line1)*2;
  if( (freechain = chainalloc(0,chainlistsize)) == NULL ){
    *mm = 0;
    return WLZ_ERR_MEM_ALLOC;
  }
  buckle(freechain, freechain, chainlistsize);
  /*
   * The allocation list size is determined as twice
   * (1 + maximum number of intervals in a line of input object).
   */
  if( (al = (WlzLAllocBuf *)
       AlcCalloc(2 * (maxinline+1), sizeof(WlzLAllocBuf))) == NULL ){
    chainFree();
    *mm = 0;
    return WLZ_ERR_MEM_ALLOC;
  }
  jvp = obj->values.v;

  /* initialise interval scanning */
  WlzInitRasterScan (obj, &iwsp, WLZ_RASTERDIR_ILIC);
  line = iwsp.linpos;
  ended = 0 ;
  lend = 0 ;
  WlzNextInterval(&iwsp);
  crntst = al + 1 + maxinline;
  nob = 0 ;
  prlast = al ;
  precst = al + 1 ;

  /* commence new line */
  while (!ended) {
    crntal = crntst-1 ;

    line++;
    if (lend) { 
      ended = 1 ;
      lend = 0 ;
    }
    else {
      /* read the intervals into the allocation buffer */
      while (iwsp.linpos <= line) {
	crntal++;
	crntal->a_int.ileft = iwsp.lftpos;
	crntal->a_int.iright = iwsp.rgtpos;
	if (WlzNextInterval(&iwsp) != 0) {
	  lend = 1 ;
	  break ;
	}
      }
    }
    crlast = crntal ;
    crntal = crntst ;
    precal = precst ;
    alloc = NULL ;

    /* test whether last interval in current line dealt with */
    while (crntal <= crlast) {
      lftcrn = crntal->a_int.ileft ;
      jrtcrn = crntal->a_int.iright ;
      rtcrn = jrtcrn+2 ;
      alprec = precal->a_link;
      lftprc = precal->a_int.ileft ;
      rtprec = precal->a_int.iright +2 ;

      /* test whether last interval in preceeding line dealt with */
      if (precal > prlast || rtcrn <= lftprc+(int )jdqt) { 
	/* is interval in current line already allocated */
	if (!alloc) { 
	  /* start a new object and allocate this interval to it */
	  /* interval list */
	  link1 = newchain2(freechain,lftcrn,jrtcrn);
	  /* line-of-intervals list */
	  link2 = newchain1(freechain, (unsigned long) link1);
	  /* "object" -- first line, last line, line-of-intervals pointer */
	  crntal->a_link  = newchain1(freechain,line) ;
	  crntal->a_link = addlink1 (freechain, crntal->a_link, (unsigned long)line);
	  crntal->a_link = addlink1 (freechain, crntal->a_link, (unsigned long) link2);
	}
	/* move on to next interval in current line */
	crntal++ ;
	alloc = NULL;
      }
      else {
	if (rtprec > lftcrn+ (int )jdqt) { 
	  /* case of overlapping intervals: */
	  /* is intvl in current line already allocated ? */
	  if (!alloc) { 
	    /* allocate this interval and add to object list */
	    alloc = alprec ;
	    link1 = alloc->l_link ;
	    crntal->a_link = alloc ;
	    /* test whether this line has already been started for this object */
	    if (link1->l_u.line != line) {
	      /* update last line */
	      link1->l_u.line = line;
	      /* add a link to the line list for this line */
	      link1 = newchain2(freechain,lftcrn,jrtcrn);
	      alloc->l_u.u_link = addlink1 (freechain, alloc->l_u.u_link, (unsigned long) link1);
	    }
	    else {
	      /* add interval to interval list for last line */
	      link1 = alloc->l_u.u_link;
	      link1->l_u.u_link = addlink2 (freechain, link1->l_u.u_link, lftcrn, jrtcrn);
	    }
	  }
	  else {
	    /* merge lists and reallocate intervals */
	    /* test whether both already allocated to same object */
	    if (alprec != alloc) { 
	      merge (&alprec, &alloc, freechain);
	      /* reallocate intervals in preceding line */
	      for (altemp = precst; altemp <= prlast; altemp++)	{
		if (altemp->a_link == alloc) { 
		  altemp->a_link = alprec ;	
		}	
	      }
	      /* reallocate intervals in current line */
	      for (altemp = crntst; altemp <= crntal; altemp++)	{
		if (altemp->a_link == alloc) { 
		  altemp->a_link = alprec ;	
		}	
	      }
	      alloc = alprec ;
	    }
	  }
	  if (rtcrn < rtprec) { 
	    /* move to next interval in this line */
	    crntal++;
	    alloc = NULL ;
	    continue;	/* the outer while loop */
	  }
	}
	/* move on to next interval in preceding line */
	precal++ ;
      }
    }

    /* all intervals in current line dealt with: */
    /* find and construct any finished objects */
    if (precst <= prlast) { 
      for (precal = precst; precal <= prlast; precal++) {
	alprec = precal->a_link ;
	/* has the object to which this interval was allocated been dealt with */
	if (alprec) { 
	  /* remove any later intervals allocated to the same object */
	  for (altemp = precal; altemp <= prlast; altemp++) {
	    if (altemp->a_link == alprec) { 
	      altemp->a_link  = NULL ;	
	    }	
	  }
	  /* test if this object has intervals in the current line */
	  /* and if so skip to end of outer loop */
	  if (crntst <= crlast) { 
	    for (altemp = crntst; altemp <= crlast; altemp++)	{
	      if (altemp->a_link == alprec)
		goto loopend;
	    }
	  }
	  /* construct object - first find line and column bounds */
	  link1 = alprec->l_link;
	  oll = (int) link1->l_u.line;
	  link1 = link1->l_link;
	  ofl = (int) link1->l_u.line;
	  link1 = alprec->l_u.u_link;
	  mkl = 0; /* Just to keep lint happy. */
	  mxkl = 0; /* Just to keep lint happy. */
	  nints = 0;
	  for (jjj=ofl; jjj<=oll; jjj++) {
	    link1 = link1->l_link;
	    link2 = link1->l_u.u_link;
	    do {
	      link2 = link2->l_link;
	      jl = link2->l_u.intv.ileft;
	      jr = link2->l_u.intv.iright;
	      if (nints == 0 || jl < mkl)
		mkl = jl ;	
	      if (nints == 0 || jr > mxkl)
		mxkl = jr ;	
	      nints++;
	      /* test for end of line */
	    } 
	    while (link2 != link1->l_u.u_link);
	  }
	  /* test whether object large enough, if not ignore it */
	  /* test for height or width less than threshold */
	  if (oll-ofl < ignlns || mxkl-mkl < ignlns) {
	    for (jjj = ofl; jjj <= oll; jjj++)	{
	      link1 = link1->l_link;
	      link2 = link1->l_u.u_link;
	      /* recover chain space */
	      join (freechain, link2) ;
	    }
	  }
	  else {

	    /* test for object array overflow */
	    if (nob >= maxNumObjs) {
	      AlcFree((void *) al);
	      chainFree();
	      *mm = nob;
	      return(WLZ_ERR_INT_DATA);
	    }
	    link1 = alprec->l_u.u_link;
	    /* set up domain and object, and update counts
	       need to test successful space allocation here */
	    jdp = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					ofl,oll,mkl,mxkl, &errNum);
	    domain.i = jdp;
	    values.v = jvp;
	    *objlist = WlzAssignObject(
	      WlzMakeMain(WLZ_2D_DOMAINOBJ,
			  domain,values,NULL,obj,&errNum), NULL);
	    objlist++;
	    nob++;
	    /* get the size correct here !!! */
	    itvl = (WlzInterval *) AlcMalloc (nints * sizeof(WlzInterval));
	    jdp->freeptr = AlcFreeStackPush(jdp->freeptr, (void *)itvl, NULL);
	    /* write intervals and interval pointers lists */
	    for (jjj=ofl; jjj<=oll; jjj++) {
	      jtvl = itvl;
	      nints = 0;
	      link1 = link1->l_link;
	      link2 = link1->l_u.u_link;
	      do {
		/*
		 * the following increment since the interval
		 * list is circular, monotone increasing,
		 * but link2 originally points at the
		 * rightmost interval (see comment at top).
		 */
		link2 = link2->l_link;
		itvl->ileft = link2->l_u.intv.ileft - mkl;
		itvl->iright = link2->l_u.intv.iright - mkl;
		itvl++;
		nints++;
		/* test for end of line */
	      } 
	      while (link2 != link1->l_u.u_link);
	      WlzMakeInterval(jjj,jdp,nints,jtvl);
	      join (freechain, link2) ;
	    }
	  }
	  /* return line list etc to free-list store */
	  join (freechain, alprec->l_u.u_link);
	  join (freechain, alprec) ;
	}
      loopend:			 
	;
      }
    }

    /* update pointers - swap the two halves of the allocation list */
    /* before getting the intervals in the next line */
    altemp = crntst ;
    crntst = precst ;
    precst = altemp ;
    prlast = crlast ;
  }

  *mm = nob ;
  AlcFree((void *) al);
  chainFree();

  return WLZ_ERR_NONE;
} 


/*
 * make a new chain: extract a link from the free-list pointed to by
 * "chainbase", link it to itself, add appropriate entries.
 * Newchain1 will handle single scalar and pointer entries so long
 * as the C compiler has pointers and integers the same length (this
 * is a language requirement).
 */
static WlzLLink *newchain1(WlzLLink *chainbase,
			      unsigned long entry)
{
  WlzLLink *newlink;

  newlink = getlink(chainbase);
  newlink->l_link = newlink;
  newlink->l_u.line = entry;
  return(newlink);
}

/* Newchain2 is for intervals */
static WlzLLink *newchain2(WlzLLink *chainbase,
			      int entry1,
			      int entry2)
{
  WlzLLink *newlink;
  newlink = getlink(chainbase);
  newlink->l_link = newlink;
  newlink->l_u.intv.ileft = entry1;
  newlink->l_u.intv.iright = entry2;
  return(newlink);
}

/*
 * Add to a chain: extract a link from the free-list pointed to by
 * "chainbase", link it to chain, add appropriate entries.
 * Addchain1 will handle single scalar and pointer entries so long
 * as the C compiler has pointers and integers the same length (this
 * is a language requirement).
 */
static WlzLLink *addlink1(WlzLLink *chainbase,
			     WlzLLink *chain,
			     unsigned long entry)
{
  WlzLLink *newlink;
  newlink = getlink(chainbase);
  newlink->l_link = chain->l_link;
  chain->l_link = newlink;
  newlink->l_u.line = entry;
  return(newlink);
}

/* Addchain2 handles intervals */
static WlzLLink *addlink2(WlzLLink *chainbase,
			     WlzLLink *chain,
			     int entry1,
			     int entry2)
{
  WlzLLink *newlink;
  newlink = getlink(chainbase);
  newlink->l_link = chain->l_link;
  chain->l_link = newlink;
  newlink->l_u.intv.ileft = entry1;
  newlink->l_u.intv.iright = entry2;
  return(newlink);
}

/*
 * allocate chain list and store information in static structure
 *
 * The size of the chain link is not easily predictable in advance.
 * We therefore permit extension of the chain list by Malloc() calls,
 * and require a structure to store a sensible chain list size
 * increment.  Allocated chunks of chain list are themselves chained
 * together (using the first chain location)
 * in order to permit space-freeing; the original base of the chain
 * and most recent chunk base are also stored in this structure.
 */
struct allocthings {
  int allocsize;
  WlzLLink *orig_base;
  WlzLLink *chunk_base;
};
static struct allocthings allocthings;

static WlzLLink *chainalloc(int flag, int n)
{
  WlzLLink *chain;
  /*
   * memory allocation size determined by first allocation
   */
  if (flag)
    n = allocthings.allocsize;
  else
    allocthings.allocsize = n;
  /*
   * allocate memory, chain to existing allocated memory
   */
  if( (chain = (WlzLLink *) AlcMalloc((n+1) * sizeof(WlzLLink))) == NULL ){
    return( NULL );
  }
  chain->l_link = NULL;
  if (flag) {
    allocthings.chunk_base->l_link = chain;
    allocthings.chunk_base = chain;
  }
  else {
    allocthings.orig_base = allocthings.chunk_base = chain;
  }
  /*
   * retain first link for chaining allocated memory chunks
   */
  return(chain+1);
}

static void chainFree(void)
{
  (void )freechain(allocthings.orig_base);
}

static int freechain(WlzLLink *l)
{
  if (l != NULL) {
    freechain(l->l_link);
    AlcFree((char *) l);
  }

  return( 0 );
}

/* get an unused link from the free chain store */
static WlzLLink *getlink(WlzLLink *chainbase)
{
  WlzLLink *newlink;
  newlink = chainbase->l_link;
  if (newlink ==chainbase) {
    newlink = chainalloc(1,0);
    buckle(chainbase,newlink,allocthings.allocsize);
    newlink = chainbase->l_link;
  }
  chainbase->l_link = newlink->l_link;
  return(newlink);
}

/* link up the chain store */
static int buckle (WlzLLink *baselink,
		   WlzLLink *chainbase,
		   int length)
{
  WlzLLink *lin, *linwas;
  chainbase->l_link = baselink;
  baselink->l_link = chainbase;
  lin = chainbase;
  for (; length > 1; length--) {
    linwas = lin++;
    lin->l_link = linwas->l_link;
    linwas->l_link = lin;
  }

  return( 0 );
}


/* merge chain list2 into chain (e.g. free-list) list1 */
static int join(WlzLLink *list1,
		WlzLLink *list2)
{ 
  WlzLLink *last1;
  if ( list1 != NULL && list2 != NULL) { 
    last1 = list1->l_link;
    list1->l_link = list2->l_link ;
    list2->l_link = last1 ;
  }

  return( 0 );
} 


/* merge two objects whose reference lists begin at k1,k2 into */
/* one object beginning at k1 */
static int merge(WlzLLink **ak1,
		 WlzLLink **ak2,
		 WlzLLink *freechain)
{ 
  WlzLLink *intvln1, *intvln2, *k1;
  WlzLLink *nk1, *nnk1, *iplp1;
  WlzLLink *k2, *nk2, *nnk2, *iplp2;
  int more, i, ll1, ll2;
  int fl1, fl2;
  do {
    k1 = *ak1 ;
    nk1 = k1->l_link;
    nnk1 = nk1->l_link;
    /* first line of k1 */
    fl1 = (int) nnk1->l_u.line;

    k2 = *ak2 ;
    nk2 = k2->l_link;
    nnk2 = nk2->l_link;
    /* first line of k2 */
    fl2 = (int) nnk2->l_u.line;

    /*
     * if k1's first line is greater than k2's first line then
     * the objects are interchanged so that k2 absorbs k1 -
     * this is why k1 and k2 are passed by address.
     */
    more = fl2-fl1 ;
    if (more<0) {
      *ak1 = k2 ;
      *ak2 = k1 ;
    }
    else {

      /* last line of k1 */
      ll1 = (int) nk1->l_u.line;
      /* interval pointer list pointer for first object */
      iplp1 = k1->l_u.u_link;
      /* last line of k2 */
      ll2 = (int) nk2->l_u.line;
      /* interval pointer list pointer for second object */
      iplp2 = k2->l_u.u_link;

      /*
       * add extra location to interval pointers
       * list of shorter object - of course, "last lines"
       * of objects being merged will differ by exactly one.
       */
      if (ll1 > ll2) { 
	iplp2 = addlink1 (freechain, iplp2, 0);
	k2->l_u.u_link = iplp2 ;
	nk2->l_u.line = ll1 ;
      } 
      else if (ll1 < ll2) {
	iplp1 = addlink1 (freechain, iplp1, 0);
	k1->l_u.u_link = iplp1 ;
	nk1->l_u.line = ll2 ;
      }
      intvln1 = iplp1->l_link;
      intvln2 = iplp2->l_link;

      for (i=0; i<more; i++)
	intvln1 = intvln1->l_link;

      /* chain lists for corresponding lines are 'monotone-joined'
	 until list of pointers to interval lists for object with
	 fewer lines is exhausted */
      do {
	monotonejoin (&(intvln1->l_u.u_link), intvln2->l_u.u_link);
	intvln1 = intvln1->l_link;
	intvln2 = intvln2->l_link;
      } 
      while (intvln1 != iplp1->l_link);
      /* return space to free list */
      join (freechain, k2) ;
      join (freechain, iplp2) ;
      return (0) ;
    }
  } 
  while (more < 0);

  return( 0 );
} 


/*
 * merge two monotone chainlists (of intervals) into one monotone chain.
 * The address of the combined list is returned via amlist1.
 * Magnitude of link determined by left end point of link interval
 */
static int monotonejoin(WlzLLink **amlist1,
			WlzLLink *mlist2)
{ 
  int mlist;
  WlzLLink *mlist1;
  WlzLLink *l1, *nl2, *l;
  mlist1 = *amlist1;
  /* if first chain is empty the join consists of the second only */
  if ( mlist1 == NULL) {
    *amlist1 = mlist2;
    return( 0 );
  }
  /* if second chain is empty the join consists of the first only */
  if ( mlist2 == NULL)
    return( 0 );
  /* if mlist1 consists of single link swap it to be second list */
  if ( mlist1 == mlist1->l_link) {
    l = mlist1 ;
    *amlist1 = mlist1 = mlist2 ;
    mlist2 = l ;
  }
  mlist = 1 ;
  /* 
   * Merge the lists monotonically.
   * When all links of mlist2 have been added to mlist1,
   * exit with amlist1 as address of combined list.
   * Note that the following code was modified to correct a bug 07-07-92,
   * in which the bizarre order of the circular lists was not fully
   * taken into account (see comment at top).  The solution adopted
   * here is robust rather than maximally efficient, in that each
   * element of list2 is compared with list1 "from scratch".
   */
  while (mlist) {
    l1 = mlist1;
    nl2 = mlist2->l_link;
    if (nl2 == nl2->l_link)
      mlist = 0 ;	
    /*
     * If 2nd list 2nd element (expected to be leftmost in 2nd list)
     * (and not distinguished from 1st element if length of list is one)
     * is less than 1st list 1st element (expected to be rightmost)
     * then search for the correct slot; otherwise the correct slot
     * is immediately after 1st list 1st element, and must be followed
     * by revision of 1st list pointer "mlist1" (see below).
     */
    if (nl2->l_u.intv.ileft < mlist1->l_u.intv.ileft) {
      for (;;) {
	l = l1->l_link;
	if (nl2->l_u.intv.ileft < l->l_u.intv.ileft)
	  break;
	l1 = l;
      }
    }
    /*
     * detach 2nd element from list2 and link into list1
     */
    mlist2->l_link = nl2->l_link ;
    nl2->l_link = l1->l_link ;
    l1->l_link = nl2 ;
    /*
     * Revise the list1 pointer if the added interval is rightmost
     */
    if (nl2->l_u.intv.ileft >= mlist1->l_u.intv.ileft)
      *amlist1 = mlist1 = nl2 ;
  }

  return( 0 );
}

/*
 * revised version of library intcount().
 * (i) only takes case of type 1 interval domain,
 * (ii) returns maximum-intervals-in-a-line
 */
static int mintcount(WlzIntervalDomain *idom,
		     int *maxinline)
{
  int	ll,intcount,max;
  WlzIntervalLine *intl;
  intl = idom->intvlines;
  max = intcount = 0;
  ll = idom->lastln - idom->line1;
  for ( ; ll >= 0; ll--) {
    if (intl->nintvs > max)
      max = intl->nintvs;
    intcount += (intl++)->nintvs;
  }
  *maxinline = max;
  return(intcount);
}
