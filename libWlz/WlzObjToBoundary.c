#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzObjToBoundary_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzObjToBoundary.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Computes a boundary list from a Woolz object.
* \ingroup	WlzBoundary
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/*
 * local structure definitions
 */
typedef struct {
  int line;
  int left;
  int right;
  WlzBoundList *lb;		/* boundary component assigned */
  WlzBoundList *rb;
} WlzBoundInterval;

typedef struct {
  WlzBoundInterval *bi;	/* interval point lies in */
  int interval_props;	/* at interval end ? and which ? */
  WlzBoundList *bl;	/* current boundary component */
  WlzBoundList *bwas;	/* previous boundary component of this point */
  int l;
  int k;
} WlzBoundPoint;
/*
 * constant definitions
 */
#define LEFT_END	1
#define RIGHT_END	2
#define UP		1
#define DOWN		2 

static WlzBoundInterval	*makeboundints	(WlzObject *obj,
					 int icount);
static int 		nextunprocessed	(WlzBoundInterval *bi,
					 WlzBoundInterval *bitop,
					 WlzBoundPoint *bp);
static int 		bnd_link	(WlzBoundList *bl,
					 WlzBoundList *bp);
static int 		newbpoint	(WlzBoundPoint *cp,
					 WlzBoundInterval *bi,
					 WlzIVertex2 **vtx,
					 int *nv,
					 int dir);
static int 		extranewbpoint	(WlzBoundPoint *cp,
					 WlzBoundInterval *bi,
					 WlzIVertex2 **vtx,
					 int *nv,
					 int dir);
static void 		vertcompress	(WlzIVertex2 **vtx,
					 int *nv);
static int 		wraparound	(WlzPolygonDomain *pdom,
					 int wrap);

/* function:     WlzObjToBoundary    */
/*! 
* \ingroup      WlzBoundary
* \brief        Compute the recursive boundary structure of input
 object "obj",with "wrap" points repeated in each boundary polygon
 to ensurewrap-around (usual value for drawing purposes : 1).
 All boundary polygon arcs lie in one of the 8 directions vertical,
 horizontal, diagonal.
*
* \return       Boundary object, NULL on error.
* \param    obj	Input domain object.
* \param    wrap	Number ov overlapping vertices within each polygon.
* \param    dstErr	Error return.
* \par      Source:
*                WlzObjToBoundary.c
*/
WlzObject *WlzObjToBoundary(
  WlzObject	*obj,
  int		wrap,
  WlzErrorNum	*dstErr)
{
  int 			prev_direc, line1, lastln;
  int 			icount, nv;
  WlzBoundInterval 	*bi, *bitop;
  WlzBoundInterval 	*ci, *ci_ov, *ci_next;
  static WlzBoundList 	uhole = {WLZ_BOUNDLIST_HOLE,0,NULL,
				 NULL,NULL,NULL,1,NULL};
  WlzBoundList 		*bnd;
  WlzBoundPoint 	cur_pt;
  WlzObject 		*bobj=NULL;
  WlzIVertex2 		*vtx, *pvtx;
  WlzDomain		domain, *rtnDomains, *domains;
  WlzValues		bvalues,
  			values;
  int			p;
  WlzObjectType		objType;
  WlzObject		*obj1, *obj2;
  char			*procStr="WlzObjToBoundary";
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /*
   * Initial checks.
   */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* switch on object type */
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      break;

    case WLZ_TRANS_OBJ:
      if( (bobj = WlzObjToBoundary(obj->values.obj, wrap, &errNum))
	 == NULL ){
	break;
      }
      bvalues.obj = bobj;
      return WlzMakeMain(obj->type, obj->domain, bvalues,
			 NULL, obj, dstErr);

    case WLZ_2D_POLYGON:
      if( (bnd = WlzMakeBoundList(WLZ_BOUNDLIST_PIECE, 0,
				  obj->domain.poly, &errNum)) == NULL){
	break;
      }
      domain.b = bnd;
      values.core = NULL;
      return WlzMakeMain(WLZ_BOUNDLIST, domain, values,
			 NULL, NULL, dstErr);

    case WLZ_BOUNDLIST:
    case WLZ_EMPTY_OBJ:
      return WlzMakeMain(obj->type, obj->domain, obj->values,
			 NULL, NULL, dstErr);

    case WLZ_3D_DOMAINOBJ:
      /* check the domain */
      if( obj->domain.p ){
	switch( obj->domain.p->type ){
	case WLZ_PLANEDOMAIN_DOMAIN:
	case WLZ_PLANEDOMAIN_POLYGON:
	  if((domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_BOUNDLIST,
					    obj->domain.p->plane1,
					    obj->domain.p->lastpl,
					    obj->domain.p->line1,
					    obj->domain.p->lastln,
					    obj->domain.p->kol1,
					    obj->domain.p->lastkl,
					    &errNum)) != NULL){
	    domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
	    domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
	    domain.p->voxel_size[2] = obj->domain.p->voxel_size[2];
	    values.core = NULL;
	    domains = obj->domain.p->domains;
	    rtnDomains = domain.p->domains;
	    objType = (obj->domain.p->type == WLZ_PLANEDOMAIN_DOMAIN) ?
	      WLZ_2D_DOMAINOBJ : WLZ_2D_POLYGON;
	    for(p=obj->domain.p->plane1; p <= obj->domain.p->lastpl;
		p++, domains++, rtnDomains++){
	      if( (*domains).core ){
		if( (obj1 = WlzMakeMain(objType, *domains, values,
					NULL, NULL, &errNum)) &&
		   (obj2 = WlzObjToBoundary(obj1, wrap, &errNum)) ){
		  *rtnDomains = WlzAssignDomain(obj2->domain, NULL);
		  WlzFreeObj(obj2);
		  WlzFreeObj(obj1);
		}
		else {
		  if( obj1 ){
		    WlzFreeObj(obj1);
		  }
		  WlzFreePlaneDomain(domain.p);
		  domain.p = NULL;
		  break;
		}
	      }
	    }
	    if( domain.p ){
	      return WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
				 NULL, NULL, dstErr);
	    }
	  }
	  break;
					       
	case WLZ_PLANEDOMAIN_BOUNDLIST:
	case WLZ_EMPTY_DOMAIN:
	  return WlzMakeMain(obj->type, obj->domain, obj->values,
			     NULL, NULL, dstErr);

	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* now check the domain pointer */
  if( (errNum == WLZ_ERR_NONE) && (obj->domain.core == NULL) ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }

  /*
   * copy object domain into more appropriate structure,
   * fix so that outside pieces are inserted in universal hole.
   * This should work fine for rectangular domains.
   */
  if( errNum == WLZ_ERR_NONE ){
    icount = WlzIntervalCount(obj->domain.i, &errNum);
    if( (bi = makeboundints(obj, icount)) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      bi->rb = &uhole;
      WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	      ("%s - uhole address planted %p\n",
	       procStr, bi->rb));
      bitop = bi + icount + 1;
    }
  }

  /*
   * allocate a vertex table.  Since each interval endpoint contributes
   * up to two points to a boundary
   * we may need up to 4*intcount vertices.
   */
  if((errNum == WLZ_ERR_NONE) &&
     ((pvtx = (WlzIVertex2 *)
       AlcMalloc(sizeof(WlzIVertex2)*icount*4)) == NULL) ){
    AlcFree((void *) bi);
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  /*
   * object line bounds
   */
  if( errNum == WLZ_ERR_NONE ){
    line1 = obj->domain.i->line1;
    lastln = obj->domain.i->lastln;
  }

  /*
   * Check for unprocessed interval end-point -
   * implies that nested boundary structure is incomplete
   */
  if( errNum == WLZ_ERR_NONE ){
    uhole.down = NULL;
    while (nextunprocessed(bi,bitop,&cur_pt) != 0) {
      /*
       * plant first boundary point in vertex list
       */
      vtx=pvtx;
      vtx->vtX = cur_pt.k;
      vtx->vtY = cur_pt.l;
      vtx++;
      nv=1;
      /*
       * Find if piece or hole;
       * Plant pointer back to containing element of hierarchy.
       * If piece, previous search direction is set to UP.
       * If hole, previous search direction is set to DOWN.
       * In either case, this will almost immediately be reversed
       * as the boundary is constructed.
       */
      if (cur_pt.interval_props == LEFT_END) {
	prev_direc = UP;
      }
      else {
	prev_direc = DOWN;
      }
      /*
       * search round the boundary for next endpoint
       * until we get back to the beginning.
       * While searching UP, we prefer to acquire interval left
       * end points, while searching DOWN, the preference is to
       * interval right end points.
       *
       * In the diagrams that follow, the current point is "O",
       * and the next interval end-point in the boundary is "X".
       * Additional points are "*" and "+".  The difference is
       * to do with whether the point is "essential" or is added
       * in order to restrict the boundary to the 8 directions.
       * If either type of additional point
       * is at an interval end-point IT MUST NOT BE RECORDED
       * AS SUCH in the bound_interval structure.
       * Other domain points "-".
       */
      while (cur_pt.bwas == NULL) {
	WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
		("%s - current point: (%3d,%3d), direction: %s,"
		 "int-prop: %s\n",
		 procStr, cur_pt.l, cur_pt.k,
		 (prev_direc == DOWN)?"DOWN":"UP",
		 (cur_pt.interval_props == LEFT_END)?"LEFT_END":"RIGHT_END"));
	ci = cur_pt.bi;
	if (prev_direc == DOWN) {
	  ci_ov = ci + 1;

	  /* this test was within an ifdef DEBUG should be removed or put
	     in to test for an error condition - test for now */
	  if (cur_pt.interval_props != RIGHT_END) {
	    AlcFree((void *) bi);
	    AlcFree((void *) pvtx);
	    return NULL;
	  }

	  /*
	   * No overlapping interval below
	   *
	   *     X----O
	   *
	   * reverse direction
	   */
	  /* last line ? */
	  if (cur_pt.l == lastln) {
	    newbpoint(&cur_pt,ci,&vtx,&nv,LEFT_END);
	    prev_direc = UP;
	    continue;	/* the outer while loop */
	  }
	  /* find next line */
	  while (ci_ov->line <= ci->line)
	    ci_ov++;
	  /* find first overlapping interval, if any */
	  while (ci_ov->line == ci->line+1 && ci_ov->right < ci->left-1)
	    ci_ov++;
	  if (ci_ov->line > ci->line+1 || ci_ov->left > ci->right+1) {
	    /* nothing overlaps */
	    newbpoint(&cur_pt,ci,&vtx,&nv,LEFT_END);
	    prev_direc = UP;
	    continue;
	  }
	  /*
	   * Adjacent or overlapping interval below,
	   * and another interval to right in this
	   * line which is also adjacent or overlaps.
	   *
	   *    ---O    X---    ---O    X---
	   *    ----*--+----    ----*--+
	   *
	   *    ---O    X---    ---O    X---
	   *        *--+----        *--+
	   *
	   * reverse direction
	   */
	  /* find rightmost overlapping interval */
	  while ((ci_ov+1)->line == ci->line+1 && (ci_ov+1)->left <= ci->right+1)
	    ci_ov++;
	  if (ci_ov->right > ci->right) {
	    /*
	     * suitable overlapping interval in next line.
	     * check if it overlaps the next
	     * interval in this line.
	     */
	    ci_next = ci + 1;
	    if (ci_next->line == ci->line && ci_next->left <= ci_ov->right+1) {
	      extranewbpoint(&cur_pt,ci_next,&vtx,&nv,LEFT_END);
	      prev_direc = UP;
	      continue;
	    }
	  }
	  /*
	   * Overlapping or adjacent interval below,
	   * Any further interval on this line not
	   * overlapping or adjacent.
	   *
	   *    ---O             --------O         +---O
	   *    ----+---X        -----X       ----X
	   *
	   *    ----O            ----O
	   *    ----X                 +---X
	   */
	
	  /* this was within an ifdef DEBUG - now an error condition test */
	  if (ci_ov->right < ci->left-1) {
	    WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
		    ("%s - Overlap logic chopped "
		     "(ci %d, %d-%d)(ov %d, %d-%d)\n",
		     ci->line, ci->left, ci->right,
		     ci_ov->line, ci_ov->left, ci_ov->right));
	    AlcFree((void *) bi);
	    AlcFree((void *) pvtx);
	    return NULL;
	  }
	  newbpoint(&cur_pt,ci_ov,&vtx,&nv,RIGHT_END);
	  continue;
	}
	else {
	  ci_ov = ci - 1;

	  /* this was within an ifdef DEBUG - now an error condition test */
	  if (cur_pt.interval_props != LEFT_END) {
	    AlcFree((void *) bi);
	    AlcFree((void *) pvtx);
	    return NULL;
	  }

	  /*
	   * No overlapping interval above
	   *
	   *     O----X
	   *
	   * reverse direction
	   */
	  /* first line ? */
	  if (cur_pt.l == line1) {
	    newbpoint(&cur_pt,ci,&vtx,&nv,RIGHT_END);
	    prev_direc = DOWN;
	    continue;	/* the outer while loop */
	  }
	  /* find next line */
	  while (ci_ov->line >= ci->line)
	    ci_ov--;
	  /* find first overlapping interval, if any */
	  while (ci_ov->line == ci->line-1 && ci_ov->left > ci->right+1)
	    ci_ov--;
	  if (ci_ov->line < ci->line-1 || ci_ov->right < ci->left-1) {
	    /* nothing overlaps */
	    newbpoint(&cur_pt,ci,&vtx,&nv,RIGHT_END);
	    prev_direc = DOWN;
	    continue;
	  }
	  /*
	   * Adjacent or overlapping interval above,
	   * and another interval to left in this
	   * line which is also adjacent or overlaps.
	   *
	   *    ----+--*----    ----+--*
	   *    ---X    O---    ---X    O---
	   *
	   *        +--*----        +--*
	   *    ---X    O---    ---X    O---
	   *
	   * reverse direction
	   */
	  /* find leftmost overlapping interval */
	  while ((ci_ov-1)->line == ci->line-1 && (ci_ov-1)->right >= ci->left-1)
	    ci_ov--;
	  if (ci_ov->left < ci->left) {
	    /*
	     * suitable overlapping interval in next line.
	     * check if it overlaps the next
	     * interval in this line.
	     */
	    ci_next = ci - 1;
	    if (ci_next->line == ci->line && ci_next->right >= ci_ov->left-1) {
	      extranewbpoint(&cur_pt,ci_next,&vtx,&nv,RIGHT_END);
	      prev_direc = DOWN;
	      continue;
	    }
	  }
	  /*
	   * Overlapping or adjacent interval above
	   * Any previous interval on this line not
	   * overlapping or adjacent.
	   *
	   *    X--+----             X----           X----
	   *        O---         O--------      O---+
	   *
	   *    X----            X---+
	   *    O----                 O---
	   */

	  /* this was within an ifdef DEBUG - now an error condition test */
	  if (ci_ov->left > ci->right+1) {
	    WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
		    ("%s - Overlap logic chopped "
		     "(ci %d, %d-%d)(ov %d, %d-%d)\n",
		     ci->line, ci->left, ci->right,
		     ci_ov->line, ci_ov->left, ci_ov->right));
	    AlcFree((void *) bi);
	    AlcFree((void *) pvtx);
	    return NULL;
	  }
	  newbpoint(&cur_pt,ci_ov,&vtx,&nv,LEFT_END);
	  continue;
	}
      }

      /* this was within an ifdef DEBUG - now an error condition test */
      if (cur_pt.bwas != cur_pt.bl) {
	AlcFree((void *) bi);
	AlcFree((void *) pvtx);
	errNum = WLZ_ERR_UNSPECIFIED;
	break;
      }

      /*
       * Attach the boundlist polygon for this component
       */
      cur_pt.bl->poly =
	WlzAssignPolygonDomain(
	  WlzMakePolygonDomain(WLZ_POLYGON_INT,nv,pvtx,nv+wrap,1, &errNum),
	  NULL);

      /*
       * fix the wrap
       */
      if (wrap > 1){
	wraparound(cur_pt.bl->poly, wrap);
      }

    }
  }

  /*
   * Make the object.  Remember that the first element of the
   * bound_interval array is the universal hole pointer, and
   * the top-level boundlist is at the second element.
   * Also remember that makemain() increments the "idom" linkcount.
   */
  if( errNum == WLZ_ERR_NONE ){
    (bi+1)->lb->linkcount=0;
    domain.b = (bi+1)->lb;
    values.core = NULL;
    bobj = WlzMakeMain(WLZ_BOUNDLIST, domain, values, NULL, NULL, &errNum);

    AlcFree((void *) pvtx);
    AlcFree((void *) bi);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return bobj;
}

/*
 * Make the bound_interval table.
 * Make an extra location at the start for the universal-hole pointer
 * which is also useful for storing line1-1.  Make an extra location
 * at the end for lastln+1.
 */
static WlzBoundInterval *makeboundints(WlzObject *obj,
				       int icount)
{
  WlzBoundInterval *bi, *bib;
  WlzIntervalWSpace iwsp;

  if( (bib = (WlzBoundInterval *) AlcCalloc(sizeof(WlzBoundInterval),
					    icount+2)) == NULL ){
    return NULL;
  }
  bi = bib;
  bi->line = obj->domain.i->line1-1;

  WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_ILIC);
  while( WlzNextInterval(&iwsp) == WLZ_ERR_NONE ){
    bi++;
    bi->left = iwsp.lftpos;
    bi->right = iwsp.rgtpos;
    bi->line = iwsp.linpos;
  }
  (++bi)->line = obj->domain.i->lastln+1;
  return(bib);
}


/*
 * Find the next unprocessed interval end-point, construct
 * a new boundary segment, and link appropriately to previous
 * segments.   The first interval is a special case since it
 * is ALWAYS part of the boundary of the first piece component.
 * To permit a uniform call and linking structure, a zero'th
 * interval has previously been set up to point to a "universal hole"
 * component, whose "down" pointer is set to the first boundary
 * piece by the first call to this function.
 */
static int nextunprocessed(WlzBoundInterval *bi,
			   WlzBoundInterval *bitop,
			   WlzBoundPoint *bp)
{
  WlzBoundList *bl;

  bp->bwas = NULL;
  /*
   * look for an unprocessed interval, find enclosing piece or
   * hole, check its down pointer.  If NULL, fill it in to point
   * to new segment.  If non-NULL, follow it down and then
   * sideways through the parallel segments.
   */
  while (++bi != bitop) {
    if (bi->lb == NULL) {
      /*
       * left end.  Previous interval right end
       * tags enclosing hole (possibly uhole).
       */
      bl = (WlzBoundList *)AlcCalloc(sizeof(WlzBoundList), 1);
      bl->type = WLZ_BOUNDLIST_PIECE;
      bl->linkcount = 1;
      bnd_link(bl,(bi-1)->rb);
      bi->lb = bl;
      bp->bi = bi;
      bp->interval_props = LEFT_END;
      bp->bl = bl;
      bp->k = bi->left;
      bp->l = bi->line;
      return 1;
    }
    if (bi->rb == NULL) {
      /*
       * right end.  Corresponding left end tags
       * enclosing piece.
       */
      bl = (WlzBoundList *) AlcCalloc(sizeof(WlzBoundList),1);
      bl->type = WLZ_BOUNDLIST_HOLE;
      bl->linkcount = 1;
      bnd_link(bl,bi->lb);
      bi->rb = bl;
      bp->bi = bi;
      bp->interval_props = RIGHT_END;
      bp->bl = bl;
      bp->k = bi->right;
      bp->l = bi->line;
      return 1;
    }
  }
  /*
   * No further unprocessed intervals
   */
  return 0;
}

static int bnd_link(WlzBoundList *bl,
		    WlzBoundList *bp)
{
  WlzBoundList *nbr;

  WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	  ("Bnd_link entry bp %p bp->type %s bl->type %s\n",
	   bp,
	   (bp->type == WLZ_BOUNDLIST_HOLE)?"HOLE":"PIECE",
	   (bl->type == WLZ_BOUNDLIST_HOLE)?"HOLE":"PIECE"));

  if( bp->type == bl->type ){
    nbr = bp;
    while (nbr->next != NULL){
      nbr = nbr->next;
    }
    nbr->next = bl;
  }
  else if( bp->down == NULL ){
    bp->down = bl;
  }
  else {
    WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	    ("Bnd_link down not NULL (%p)\n", bp->down));

    nbr = bp->down;
    while (nbr->next != NULL){
      nbr = nbr->next;
    }
    nbr->next = bl;
  }
  WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	  ("Bnd_link exit\n"));

  return 0;
}


static int newbpoint(WlzBoundPoint *cp,
		     WlzBoundInterval *bi,
		     WlzIVertex2 **vtx,
		     int *nv,
		     int dir)
{
  WlzIVertex2 *wtx = (*vtx)-1;
  int xdiff, ydiff;

  switch( dir ){
  case LEFT_END:
    cp->bwas = bi->lb;
    bi->lb = cp->bl;
    cp->k = bi->left;
    break;
  case RIGHT_END:
    cp->bwas = bi->rb;
    bi->rb = cp->bl;
    cp->k = bi->right;
    break;
  }
  cp->bi = bi;
  cp->l = bi->line;

  /*
   * Check whether an additional point must be inserted:
   * ydiff non-zero and xdiff > 1.
   */
  ydiff = cp->l - wtx->vtY;
  if( ydiff > 0 ){
    xdiff = cp->k - wtx->vtX;
    if (xdiff > 1) {
      (*vtx)->vtX = wtx->vtX + 1;
      (*vtx)->vtY = cp->l;
      (*nv)++;
      vertcompress(vtx, nv);
      (*vtx)++;
    }
    else if( xdiff < -1 ){
      (*vtx)->vtX = cp->k + 1;
      (*vtx)->vtY = wtx->vtY;
      /* don't need vertcompress() call here */
      (*vtx)++;
      (*nv)++;
    }
  }
  else if( ydiff < 0 ){
    xdiff = cp->k - wtx->vtX;
    if( xdiff > 1 ){
      (*vtx)->vtX = cp->k - 1;
      (*vtx)->vtY = wtx->vtY;
      /* don't need vertcompress() call here */
      (*vtx)++;
      (*nv)++;
    }
    else if( xdiff < -1 ){
      (*vtx)->vtX = wtx->vtX - 1;
      (*vtx)->vtY = cp->l;
      (*nv)++;
      vertcompress(vtx,nv);
      (*vtx)++;
    }
  }
  cp->interval_props = dir;
  (*vtx)->vtX = cp->k;
  (*vtx)->vtY = cp->l;
  (*nv)++;
  vertcompress(vtx,nv);
  (*vtx)++;
  return( 0 );
}

static int extranewbpoint(WlzBoundPoint *cp,
			  WlzBoundInterval *bi,
			  WlzIVertex2 **vtx,
			  int *nv,
			  int dir)
{
  switch( dir ){
  case LEFT_END:
    (*vtx)->vtX = cp->k + 1;
    (*vtx)->vtY = cp->l + 1;
    break;
  case RIGHT_END:
    (*vtx)->vtX = cp->k - 1;
    (*vtx)->vtY = cp->l - 1;
    break;
  }

  WLZ_DBG((WLZ_DBG_ALLOC|WLZ_DBG_LVL_1),
	  ("WlzObjToBoundary - extranewbpoint: (%d,%d), "
	   "direction %s)(extra)\n",
	   (*vtx)->vtY, (*vtx)->vtX,
	   (dir == UP)?"UP":"DOWN"));

  (*nv)++;
  vertcompress(vtx,nv);
  (*vtx)++;
  newbpoint(cp,bi,vtx,nv,dir);
  return 0;
}


/*
 * Compress colinear runs of vertices.  Since we only permit horizontal,
 * vertical, and diagonal directions, a check for non-reversal plus both
 * segments zero or both non-zero in each of X and Y is sufficient.
 * ENSURE where this function is called that vtx is the most recent vertex
 * and nv is the vertex count INCLUDING the most recent.
 */
static void vertcompress(WlzIVertex2 **vtx,
			 int *nv)
{
  WlzIVertex2	*v1, *v2, *v3;
  int		x1, x2, y1, y2;

  v1 = *vtx;
  v2 = v1-1;
  /*
   * two points may be identical if interval where up-down
   * direction reverses is only one point long
   */
  x1 = v1->vtX - v2->vtX;
  y1 = v1->vtY - v2->vtY;
  if (x1 != 0 || y1 != 0) {
    if (*nv < 3)
      return;
    v3 = v1-2;
    x2 = v2->vtX - v3->vtX;
    /* check if X reversal */
    if ((x1 <= 0 && x2 > 0) || (x1 >= 0 && x2 < 0) || (x1 != 0 && x2 == 0))
      return;
    y2 = v2->vtY - v3->vtY;
    /* check if Y reversal */
    if ((y1 <= 0 && y2 > 0) || (y1 >= 0 && y2 < 0) || (y1 != 0 && y2 == 0))
      return;
  }
  /* colinear or x1=y1=0  */
  v2->vtX = v1->vtX;
  v2->vtY = v1->vtY;
  (*vtx)--;
  (*nv)--;

  return;
}


static int wraparound(WlzPolygonDomain *pdom,
		      int wrap)
{
  WlzIVertex2 *vtx, *btx;

  btx = pdom->vtx;
  /*
   * There is a default wrap of 1, included in pdom->nvertices
   */
  wrap--;
  vtx = btx + pdom->nvertices;
  btx++;
  pdom->nvertices += wrap;
  for (; wrap>0; wrap--,btx++,vtx++) {
    vtx->vtX = btx->vtX;
    vtx->vtY = btx->vtY;
  }

  return 0;
}



