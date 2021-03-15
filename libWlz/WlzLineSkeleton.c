#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLineSkeleton_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzLineSkeleton.c
* \author       Bill Hill
* \date         July 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \brief	Functions to compute the line skeleton of an object's domain.
* \ingroup	WlzFeatures
* \todo 	The node priority queue is implemented as a simple doubly
* 		linked list that is sorted on insert. There are some hacks
* 		to guess nearby nodes and speed it up a little, but the
* 		queue is really inefficient and should be raplaced by
* 		something like AlcHeap().
*/

#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
* \struct 	_WlzLSDomEnt
* \ingroup	WlzFeatures
* \brief 	A queue entry holding a domain object (which still has a
* 		line skeleton segment to be determined) along with source
* 		and destination line skeleton segment end points. For use
* 		with AlcHeap().
* 		Typedef: ::WlzLSDomEnt
*/
typedef struct _WlzLSDom
{

  double        priority;	/*!< Heap priority, use domain size. */
  WlzObject	*obj;		/*!< Object holding domain. */
  WlzIVertex3	src;		/*!< Start of path in domain. */
  WlzIVertex3   dst;		/*!< End of path in domain. */
}  WlzLSDomEnt;

/*!
* \enum 	_WlzLSNodState
* \ingroup	WlzFeatures
* \brief	State of a node when computing a line skeleton segment.
* 		Typedef: ::WlzLSNodState
*/
typedef enum _WlzLSNodState
{
  WLZLS_NOD_NEW   	= 0,	/*!< New unused node. */
  WLZLS_NOD_OPEN	= 1,	/*!< Open node in frontier. */
  WLZLS_NOD_CLOSED	= 2,	/*!< Closed node which is no longer active. */
} WlzLSNodState;

/*!
* \struct	_WlzLSNod
* \ingroup	WlzFeatures
* \brief	A node entry in a node queue data structure used when
* 		computing a line skeleton segment.
* 		Typedef: ::WlzLSNod
*/
typedef struct _WlzLSNod
{
  struct _WlzLSNod *prevQ;	/*!< Previous node in the priority queue. */
  struct _WlzLSNod *nextQ;	/*!< Next node in the priority queue. */
  struct _WlzLSNod *nextT;	/*!< Next node in the hash table. */
  struct _WlzLSNod *from;	/*!< Node from which the frontier propagated to
  				     this node. */
  WlzIVertex3	pos;		/*!< Position of the node. */
  WlzLSNodState	state;		/*!< Node state. */
  double	cost;		/*!< Cost from source to this node. */
} WlzLSNod;

/*!
* \struct	_WlzLSNodQueue
* \ingroup	WlzFeatures
* \brief	A queue of line segment nodes together with a table
* 		for accessing them by coordinate value.
* 		The slots of the hash table are allocated at creation of
* 		the queue and the queue/table entries are allocated in
* 		blocks.
*/
typedef struct _WlzLSNodQueue
{
  int		nSlt;		/*!< Number of hash table slots. Must not
  				     be changed once allocated. */
  struct _WlzLSNod *top;	/*!< Top node in the queue with circular doubly
  				     linked nodes. Lowest cost nodes are at the
				     top of the queue. */
  struct _WlzLSNod **tbl;	/*!< Node hash table with singly linked nodes
  				     in slots and slots indexed by position
				     based hash key. */
  struct _WlzLSNod *freeQ;	/*!< Singly linked list of available allocated
  				     nodes. */
  size_t	blkSz;		/*!< Number of nodes to allocate in a block. */
  struct _AlcBlockStack *freeStack; /*!< Free stack used to efficiently
  				         allocate and free blocks of nodes. */
} WlzLSNodQueue;


static void			WlzLSNodQInsertNod(
				  WlzLSNodQueue *q,
				  WlzLSNod *m,
				  WlzLSNod *n);
static void			WlzLSNodQFree(
				  WlzLSNodQueue *q);
static int			WlzLSPosCmp(
				  WlzIVertex3 p0,
				  WlzIVertex3 p1);
static unsigned int		WlzLSNodQKey(
				  WlzLSNodQueue *q,
				  WlzIVertex3 p);
static WlzErrorNum		WlzLSSplitDom(
				  AlcHeap *objQ,
				  WlzObject *obj,
				  WlzObject *skel,
				  int bDstMax);
static WlzErrorNum 		WlzLSBuildPath(
				  WlzObjectType sObjType,
				  WlzGreyValueWSpace *sWSp,
				  WlzPoints *pts,
				  WlzObject *obj,
				  WlzIVertex3 ps,
				  WlzIVertex3 pd,
				  WlzGreyValueWSpace *bDstWSp,
				  WlzGreyValueWSpace *dDstWSp);
static WlzLSNod			*WlzLSNodQNewNod(
				  WlzLSNodQueue *q);
static WlzLSNod 		*WlzLSNodQPop(
				  WlzLSNodQueue *q);
static WlzLSNod			*WlzLSNodQMatch(
				  WlzLSNodQueue *q,
				  WlzIVertex3 p,
				  int createFlg);
static WlzLSNodQueue		*WlzLSNodQueueNew(
				  size_t nSlt,
				  size_t blkSz);

/*!
* \return	New Woolz object, the domain of which is the line skeleton
* 		or NULL on error.
* \ingroup	WlzFeatures
* \brief	Computes and returns the line skeleton of the given object's
* 		domain.
* \param	gObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzLineSkeleton(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr)
{
  int		bDstMax;
  WlzIVertex3	p0,
  		p1;
  WlzObject 	*bDst = NULL,
		*pDst = NULL,
		*p0Dmn = NULL;
  WlzGreyValueWSpace *bDstWSp = NULL,
		     *pDstWSp = NULL,
		     *sWSp = NULL;
  WlzLSDomEnt 	*qEnt;
  AlcHeap	*objQ = NULL;
  WlzObject	*skel = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Compute object in which values are distance from the object boundary. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*bDmn;

    bDmn = WlzAssignObject(WlzBoundaryDomain(gObj, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      bDst = WlzAssignObject(
             WlzDistanceTransform(gObj, bDmn, WLZ_OCTAGONAL_DISTANCE, 0.0, 0.0,
                                  &errNum), NULL);
    }
    (void )WlzFreeObj(bDmn);
    if(errNum == WLZ_ERR_NONE)
    {
      bDstWSp = WlzGreyValueMakeWSp(bDst, &errNum);
    }
  }
  /* Find a maximum valued location within the boundary distance object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzPixelV	maxDst;

    p0 = WlzGreyExtremumPos(bDst, 1, &maxDst, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzValueConvertPixel(&maxDst, maxDst, WLZ_GREY_INT);
      bDstMax = maxDst.v.inv;
      p0Dmn = WlzMakeSinglePixelObject(gObj->type, p0.vtX, p0.vtY, p0.vtZ,
      			&errNum);
    }
  }
  /* Compute object in which values are distance from the maximum valued
   * location within the boundary distance object. This is the initial seed
   * point. */
  if(errNum == WLZ_ERR_NONE)
  {
    pDst = WlzAssignObject(
           WlzDistanceTransform(gObj, p0Dmn, WLZ_OCTAGONAL_DISTANCE, 0.0, 0.0,
                                &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      pDstWSp = WlzGreyValueMakeWSp(pDst, &errNum);
    }
  }
  (void )WlzFreeObj(p0Dmn);
  /* Find point in the domain which is most distant from the initial seed point
   * this is the initial destination. */
  if(errNum == WLZ_ERR_NONE)
  {
    p1 = WlzGreyExtremumPos(pDst, 1, NULL, &errNum);
  }
  /* Create a heap into which to put domains for which the line skeleton is
   * still to be determined and initialise the queue with the given domain . */
  if(errNum == WLZ_ERR_NONE)
  {
    if((objQ = AlcHeapNew(sizeof(WlzLSDomEnt), 0, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzLSDomEnt ent;
      WlzValues	nullVal = {0};

      ent.src = p0;
      ent.dst = p1;
      ent.obj = WlzAssignObject(
      		WlzMakeMain(gObj->type, gObj->domain, nullVal,
		            NULL, NULL, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
        ent.priority = WlzVolume(ent.obj, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(AlcHeapInsertEnt(objQ, &ent) != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }
  /* Create a new object for the skeleton. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObjectType gTT;
    WlzPixelV	  zV;

    zV.type = WLZ_GREY_INT;
    zV.v.inv = 0;
    gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
    skel = WlzNewObjectValues(gObj, gTT, zV, 1, zV, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
     sWSp = WlzGreyValueMakeWSp(skel, &errNum);
    }
  }
  /* Construct a central line skeleton path within the domain from the seed
   * point. */
  while((errNum == WLZ_ERR_NONE) && ((qEnt = AlcHeapTop(objQ)) != NULL))
  {
    /* Compute new minimum cost path from destination to source. */
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzLSBuildPath(gObj->type, sWSp, NULL, qEnt->obj, qEnt->dst,
          qEnt->src, bDstWSp, pDstWSp);
    }
    /* Remove object covered by new path, split current object and add
     * parts to the queue. */
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzLSSplitDom(objQ, qEnt->obj, skel, bDstMax);
    }
  }
  AlcHeapFree(objQ);
  WlzGreyValueFreeWSp(bDstWSp);
  WlzGreyValueFreeWSp(pDstWSp);
  WlzGreyValueFreeWSp(sWSp);
  (void )WlzFreeObj(bDst);
  (void )WlzFreeObj(pDst);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(skel);
    skel = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(skel);
}

/*!
* \return	New Woolz object, the domain of which is the line skeleton
* 		segment or NULL on error.
* \ingroup	WlzFeatures
* \brief	Computes and returns the line skeleton segment between to
* 		given points in the given object's domain. This function
* 		is unlikely to be symetric so swapping the two points will
* 		probably give a different path.
* \param	gObj			Given object.
* \param	rObjType		Return object type, must be
* 					WLZ_2D_DOMAINOBJ, WLZ_3D_DOMAINOBJ
* 					or WLZ_POINTS.
* \param	p0			First point in object domain.
* \param	p1			Second point in object domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzLineSkeletonSegment(
				  WlzObject *gObj,
				  WlzObjectType rObjType,
				  WlzIVertex3 p0,
				  WlzIVertex3 p1,
				  WlzErrorNum *dstErr)
{
  WlzObject 	*bDst = NULL,
		*pDst = NULL;
  WlzPoints	*pts = NULL;
  WlzGreyValueWSpace *bDstWSp = NULL,
		     *pDstWSp = NULL,
		     *sWSp = NULL;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(rObjType)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	if(rObjType != gObj->type)
	{
	  errNum = WLZ_ERR_OBJECT_TYPE;
	}
        break;
      case WLZ_POINTS:
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(!WlzInsideDomain(gObj, p0.vtZ, p0.vtY, p0.vtX, NULL) ||
       !WlzInsideDomain(gObj, p1.vtZ, p1.vtY, p1.vtX, NULL))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Compute object in which values are distance from the object boundary. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*bDmn;

    bDmn = WlzAssignObject(WlzBoundaryDomain(gObj, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      bDst = WlzAssignObject(
             WlzDistanceTransform(gObj, bDmn, WLZ_OCTAGONAL_DISTANCE, 0.0, 0.0,
                                  &errNum), NULL);
    }
    (void )WlzFreeObj(bDmn);
    if(errNum == WLZ_ERR_NONE)
    {
      bDstWSp = WlzGreyValueMakeWSp(bDst, &errNum);
    }
  }
  /* Compute object in which values are distances from the first point. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*p0Dmn = NULL;

    p0Dmn = WlzMakeSinglePixelObject(gObj->type, p0.vtX, p0.vtY, p0.vtZ,
      			             &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      pDst = WlzAssignObject(
	     WlzDistanceTransform(gObj, p0Dmn, WLZ_OCTAGONAL_DISTANCE, 0.0, 0.0,
                                &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      pDstWSp = WlzGreyValueMakeWSp(pDst, &errNum);
    }
    (void )WlzFreeObj(p0Dmn);
  }
  /* Create a new object for the skeleton line segment. */
  if(errNum == WLZ_ERR_NONE)
  {
    switch(rObjType)
    {
      case WLZ_2D_DOMAINOBJ:  /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	{
	  WlzObjectType gTT;
	  WlzPixelV	  zV;

	  zV.type = WLZ_GREY_INT;
	  zV.v.inv = 0;
	  gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
	  rObj = WlzNewObjectValues(gObj, gTT, zV, 1, zV, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    sWSp = WlzGreyValueMakeWSp(rObj, &errNum);
	  }
	}
	break;
      case WLZ_POINTS:
        {
	  int		maxPts;
	  WlzObjectType pdt;
	  WlzVertexP	vxp = {0};

	  pdt = (gObj->type == WLZ_2D_DOMAINOBJ)? WLZ_POINTS_2I: WLZ_POINTS_3I;
	  /* Choose a maximum number of points equal to twice the max distance
	   * from p0 to p1. */
          WlzGreyValueGet(pDstWSp, p1.vtZ, p1.vtY, p1.vtX);
	  maxPts = 2 * pDstWSp->gVal[0].inv;
          pts = WlzMakePoints(pdt, 0, vxp, maxPts, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzDomain	pDom;
	    WlzValues	nullVal = {0};

	    pDom.pts = pts;
	    rObj = WlzMakeMain(WLZ_POINTS, pDom, nullVal, NULL, NULL, &errNum);
	    if(errNum != WLZ_ERR_NONE)
            {
	      (void )WlzFreeDomain(pDom);
	      pDom.core = NULL;
	    }
	  }
	}
	break;
      default:
        break;
    }
  }
  /* Construct a central line skeleton segment within the domain between
   * the two points. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzLSBuildPath(rObjType, sWSp, pts, gObj, p1, p0,
        bDstWSp, pDstWSp);
  }
  WlzGreyValueFreeWSp(bDstWSp);
  WlzGreyValueFreeWSp(pDstWSp);
  WlzGreyValueFreeWSp(sWSp);
  (void )WlzFreeObj(bDst);
  (void )WlzFreeObj(pDst);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Compute the line skeleton segment between the source and
* 		destination points within the given domain.
* \param	sObjType	Skeleton object type: Must be WLZ_2D_DOMAINOBJ,
* 				WLZ_3D_DOMAINOBJ or WLZ_POINTS.
* \param	sWSp		Skeleton workspace to which line segment
* 				should be added, may be NULL if sObjType
* 				is neither WLZ_2D_DOMAINOBJ or
* 				WLZ_3D_DOMAINOBJ.
* \param	pts		Points domain of the appropriate type and
* 				with sufficient points allocated, may be NULL if
* 				sObjType is not WLZ_POINTS.
* \param	obj		Given object with domain.
* \param	ps		Source point.
* \param	pd		Destination point.
* \param	bDstWSp		Workspace for distance from domain boundary.
* \param	dDstWSp		Workspace for distance from destination point.
*/
static WlzErrorNum 		WlzLSBuildPath(
				  WlzObjectType sObjType,
				  WlzGreyValueWSpace *sWSp,
				  WlzPoints *pts,
				  WlzObject *obj,
				  WlzIVertex3 ps,
				  WlzIVertex3 pd,
				  WlzGreyValueWSpace *bDstWSp,
				  WlzGreyValueWSpace *dDstWSp)
{
  const int	*nNbr;
#ifdef WLZ_LS_DEBUG
  int		debugCnt = 0;
#endif
  WlzLong	qsz = 0;
  WlzLSNod	*cNod,
  		*complete = NULL;
  WlzLSNodQueue	*nodQ = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	nNbr2[] = {9,17};
  const int	nNbr3[] = {0,26};
  const int	nbrTbl[26][3] =
  {
    {-1, -1, -1}, { 0, -1, -1}, { 1, -1, -1}, 
    { 1,  0, -1}, { 1,  1, -1}, { 0,  1, -1}, 
    {-1,  1, -1}, {-1,  0, -1}, { 0,  0, -1}, 
    {-1, -1,  0}, { 0, -1,  0}, { 1, -1,  0}, 
    { 1,  0,  0},               { 1,  1,  0}, 
    { 0,  1,  0}, {-1,  1,  0}, {-1,  0,  0}, 
    { 0,  0,  1}, {-1,  0,  1}, {-1,  1,  1}, 
    { 0,  1,  1}, { 1,  1,  1}, { 1,  0,  1}, 
    { 1, -1,  1}, { 0, -1,  1}, {-1, -1,  1}
  };

  nNbr = (obj->type == WLZ_2D_DOMAINOBJ)? nNbr2: nNbr3;
  /* Initialise node queue and get a first node. */
  qsz = WlzVolume(obj, &errNum) / 4;
  if(errNum == WLZ_ERR_NONE)
  {
    if(((nodQ = WlzLSNodQueueNew(qsz, qsz)) == NULL) ||
       ((cNod = WlzLSNodQNewNod(nodQ)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Set first node to start position and add to the open node queue. */
    cNod->pos = ps;
    cNod->from = NULL;
    cNod->cost = 0.0;
    WlzLSNodQInsertNod(nodQ, NULL, cNod);
  }
  while((errNum == WLZ_ERR_NONE) && !complete &&
        ((cNod = WlzLSNodQPop(nodQ)) != NULL))
  {
#ifdef WLZ_LS_DEBUG
    ++debugCnt;
#endif
    /* Check if we've arrived at the destination with the current node. */
    if((cNod->pos.vtX == pd.vtX) &&
       (cNod->pos.vtY == pd.vtY) &&
       (cNod->pos.vtZ == pd.vtZ))
    {
      complete = cNod;
    }
    else
    {
      int	idN;
      WlzLSNod	*lastNod = NULL;

      /* Examine each neighbour of the current node. */
      for(idN = nNbr[0]; idN < nNbr[1]; ++idN)
      {
	const int   *nP;
        WlzIVertex3 nPos;

	nP = nbrTbl[idN];
	nPos.vtX = cNod->pos.vtX + nP[0];
	nPos.vtY = cNod->pos.vtY + nP[1];
	nPos.vtZ = cNod->pos.vtZ + nP[2];
	if(WlzInsideDomain(obj, nPos.vtZ, nPos.vtY, nPos.vtX, NULL))
	{
	  WlzLSNod *nNod;

	  if((nNod = WlzLSNodQMatch(nodQ, nPos, 1)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	  if(nNod->state != WLZLS_NOD_CLOSED)
	  {
	    double	b,
	    		c;

	    WlzGreyValueGet(bDstWSp, nPos.vtZ, nPos.vtY, nPos.vtX);
	    WlzGreyValueGet(dDstWSp, nPos.vtZ, nPos.vtY, nPos.vtX);
	    b = bDstWSp->gVal[0].inv;
	    c = cNod->cost + (dDstWSp->gVal[0].inv) / (0.1 + (b * b));
	    if((nNod->state == WLZLS_NOD_NEW) ||
	       ((nNod->state == WLZLS_NOD_OPEN) && (c < nNod->cost)))
	    {
	      nNod->cost = c;
	      nNod->from = cNod;
	      WlzLSNodQInsertNod(nodQ, lastNod, nNod);
	    }
	    if(nNod->state == WLZLS_NOD_OPEN)
	    {
	      lastNod = nNod;
	    }
	  }
	}
      }
    }
    /* Mark this node closed. */
    cNod->state = WLZLS_NOD_CLOSED;
  }
  if(complete)
  {
    WlzLSNod 	*nod;

    nod = complete;
    switch(sObjType)
    {
      case WLZ_POINTS:
	{
	  int	idx = 0;
	  if(obj->type == WLZ_2D_DOMAINOBJ) /* 2D */
	  {
	    WlzIVertex2 *vxp;

	    vxp = pts->points.i2;
	    while(nod && (idx < pts->maxPoints))
	    {
	      vxp[idx].vtX = nod->pos.vtX;
	      vxp[idx].vtY = nod->pos.vtY;
	      nod = nod->from;
	      ++idx;
	    }
	  }
	  else /* 3D */
	  {
	    WlzIVertex3 *vxp;

	    vxp = pts->points.i3;
	    while(nod && (idx < pts->maxPoints))
	    {
	      vxp[idx] = nod->pos;
	      nod = nod->from;
	      ++idx;
	    }
	  }
	  pts->nPoints = idx;
	}
        break;
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	while(nod)
	{
	  WlzIVertex3 pos;

	  pos = nod->pos;
	  WlzGreyValueGet(sWSp, pos.vtZ, pos.vtY, pos.vtX);
	  WlzGreyValueGet(bDstWSp, pos.vtZ, pos.vtY, pos.vtX);
	  *(sWSp->gPtr[0].inp) = bDstWSp->gVal[0].inv;
	  nod = nod->from;
	}
        break;
      default:
        break;
    }
    /* Add line skeleton segment to the skeleton object. */

#ifdef WLZ_LS_DEBUG
    (void )fprintf(stderr, "WlzLSBuildPath() %d\n", debugCnt);
#endif
  }
  WlzLSNodQFree(nodQ);
  return(errNum);
}


/*!
* \return	New line skeleton node queue.
* \ingroup	WlzFeatures
* \brief	Create a new line skeleton node queue.
* \param        nSlt		Number of hash table slots, default if <= 0.
* \param	blkSz		Number of nodes to allocate in each block,
* 				default if <= 0.
*/
static WlzLSNodQueue		*WlzLSNodQueueNew(
				  size_t nSlt,
				  size_t blkSz)
{
  WlzLSNodQueue	*q;

  if((q = (WlzLSNodQueue *)AlcCalloc(1, sizeof(WlzLSNodQueue))) != NULL)
  {
    q->nSlt = (nSlt > 0)? nSlt: 1024;
    q->blkSz = (blkSz > 0)? blkSz: 1024;
    if((q->tbl = (WlzLSNod **)AlcCalloc(q->nSlt, sizeof(WlzLSNod *))) == NULL)
    {
      AlcFree(q);
      q = NULL;
    }
  }
  return(q);
}

/*!
* \ingroup	WlzFeatures
* \brief	Free the node queue and all it's entry nodes.
* \param	q		Given node queue.
*/
static void			WlzLSNodQFree(
				  WlzLSNodQueue *q)
{
  AlcFree(q->tbl);
  (void )AlcBlockStackFree(q->freeStack);
  AlcFree(q);
}

/*!
* \return	New node of NULL on error.
* \ingroup	WlzFeatures
* \brief	Get a new node.
* \param	q		Node queue.
*/
static WlzLSNod			*WlzLSNodQNewNod(
				  WlzLSNodQueue *q)
{
  WlzLSNod	*n = NULL;

  /* Make sure there are nodes avaliable on the free node queue,
   * if not then allocate some more. */
  if(q->freeQ == NULL)
  {
    /* Allocate a block of new nodes. */
    if((q->freeStack = AlcBlockStackNew((size_t )q->blkSz,
       				sizeof(WlzLSNod), q->freeStack,
				NULL)) != NULL)
    {
      int	idN;
      WlzLSNod 	*n0;

      /* Put the block of new nodes onto the free queue. */
      n0 = (WlzLSNod *)(q->freeStack->elements);
      q->freeQ = n0;
      for(idN = 1; idN < q->freeStack->maxElm; ++idN)
      {
	WlzLSNod	*n1;

	n1 = n0 + 1;
	n0->nextQ = n1;
	n0 = n1;
      }
      n0->nextQ = NULL;
    }
  }
  /* Get a node from the free queue. */
  if(q->freeQ)
  {
    n = q->freeQ;
    q->freeQ = n->nextQ;
  }
  if(n)
  {
    (void )memset(n, 0, sizeof(WlzLSNod));
  }
  return(n);
}

/*!
* \return	Node at the top of the queue or NULL if the queue is empty.
* \ingroup	WlzFeatures
* \brief	Pops the node with the lowest cost from the top of the
* 		queue.
* \param	q			Node queue.
*/
static WlzLSNod 		*WlzLSNodQPop(
				  WlzLSNodQueue *q)
{
  WlzLSNod	*n;

  n = q->top;
  if(n)
  {
    if(n == n->nextQ)
    {
      q->top = NULL;
    }
    else
    {
      q->top = n->nextQ;
      n->prevQ->nextQ = q->top;
      n->nextQ->prevQ = q->top;
      q->top->prevQ = n->prevQ;
      q->top->nextQ = n->nextQ->nextQ;
    }
  }
  return(n);
}

/*!
* \return	Matched (possibly new) node of NULL on error.
* \ingroup	WlzFeatures
* \brief	Finds the node with the matching position in the queue.
* 		If the create flag is set and a matching node is not found
* 		then one is created, this should then be inserted into the
* 		queue.
* \param	q			Node queue.
* \param	p			Position of node.
* \param	createFlg		Create a new node (if match fails)
* 					when non-zero.
*/
static WlzLSNod			*WlzLSNodQMatch(
				  WlzLSNodQueue *q,
				  WlzIVertex3 p,
				  int createFlg)
{
  int		c;
  unsigned int  k;
  WlzLSNod	*n;

  k = WlzLSNodQKey(q, p);
  n = q->tbl[k];
  while(n && ((c = WlzLSPosCmp(n->pos, p)) < 0))
  {
    n = n->nextT;
  }
  if(!n || c)
  {
    if((n = WlzLSNodQNewNod(q)) != NULL)
    {
      n->pos = p;
    }
  }
  return(n);
}

/*!
* \ingroup	WlzFeatures
* \brief	Either inserts or removes and then re-inserts the given node
* 		into the node queue. The given node will be placed into the
* 		queue's hash table and priority queue.
* \param	q		Node queue.
* \param	m		A nearby node to speed up search, may be NULL.
* \param	n		Node to insert/reinsert into the queue.
*/
static void			WlzLSNodQInsertNod(
				  WlzLSNodQueue *q,
				  WlzLSNod *m,
				  WlzLSNod *n)
{
  unsigned int	k;
  WlzLSNod 	**hd;
#ifdef WLZ_LS_DEBUG
  int 		debugQDepth = 0,
  		debugHDepth = 0;
#endif

  /* Look for the node in the queue's hash table and insert it if it isn't
   * already in it. */
  n->nextT = NULL;
  k = WlzLSNodQKey(q, n->pos);
  hd = q->tbl + k;
  if(*hd == NULL)
  {
    *hd = n;
  }
  else
  {
    int		c = -1;

    c = WlzLSPosCmp((*hd)->pos, n->pos);
    if(c > 0)
    {
      n->nextT = *hd;
      *hd = n;
    }
    else if(c < 0)
    {
      WlzLSNod 	*n0,
		*n1;

      n0 = *hd;
      n1 = n0->nextT;
      if(!n1)
      {
        n0->nextT = n;
      }
      else
      {
	while(((c = WlzLSPosCmp(n1->pos, n->pos)) < 0) && (n1->nextT != NULL))
	{
	  n0 = n1;
	  n1 = n1->nextT;
#ifdef WLZ_LS_DEBUG
          ++debugHDepth;
#endif
	}
	if(c)
	{
	  n0->nextT = n;
	  n->nextT = n1;
	}
      }
    }
  }
  /* If the node is already in the queue, remove it from the queue. */
  if(n->state == WLZLS_NOD_OPEN)
  {
    WlzLSNod  *n0,
	      *n1;

    n0 = n->prevQ;
    n1 = n->nextQ;
    if(n0 == n)
    {
      q->top = NULL;
    }
    else
    {
      n0->nextQ = n1;
      n1->prevQ = n0;
      if(q->top == n)
      {
	q->top = n1;
      }
      if(m == NULL)
      {
	/* If not given a nearby node, guess that priority won't change much. */
        m = n0;
      }
    }
  }
  /* Insert the node into the queue. */
  if(q->top == NULL)
  {
    q->top = n;
    n->prevQ = n->nextQ = n;
  }
  else
  {
    double	cost;
    WlzLSNod	*tn;
    
    tn = q->top;
    cost = n->cost;
    if(!(cost > tn->cost))
    {
      n->nextQ = tn;
      n->prevQ = tn->prevQ;
      q->top = n;
    }
    else if(!(cost < tn->prevQ->cost))
    {
      n->nextQ = tn;
      n->prevQ = tn->prevQ;
    }
    else
    {
      int	c;
      WlzLSNod	*lcn,
		*hcn;

      if((m == NULL) || (m->state != WLZLS_NOD_OPEN))
      {
	double	d0,
		d1;

	/* If not given a guess at a nearby node then use either the lowest cost
	 * (top) node or the highest cost node (top->prevQ) depending on which
	 * has the closest cost value. */
	d0 = cost - tn->cost;
	d1 = tn->prevQ->cost - cost;
	m = (d0 < d1)? tn: tn->prevQ;
      }
      lcn = hcn = m;
      if(cost < m->cost)
      {
	if(m == tn)
	{
	  hcn = hcn->nextQ;
	}
	else
	{
	  lcn = lcn->prevQ;
	}
	while((c = (hcn->cost > cost)) && (lcn != tn))
	{
	  hcn = lcn;
	  lcn = lcn->prevQ;
#ifdef WLZ_LS_DEBUG
	  ++debugQDepth;
#endif
	}
	if(c)
	{
	  n->prevQ = tn;
	}
	else
	{
	  n->prevQ = hcn;
	}
	n->nextQ = n->prevQ->nextQ;
      }
      else /* !(cost < m->cost) */
      {
	if(m == tn->prevQ)
	{
	  lcn = lcn->nextQ;
	}
	else
	{
	  hcn = hcn->nextQ;
	}
        while((c = (lcn->cost < cost)) && (hcn != tn->prevQ))
	{
	  lcn = hcn;
	  hcn = hcn->nextQ;
#ifdef WLZ_LS_DEBUG
	  ++debugQDepth;
#endif
	}
	if(c)
	{
	  n->nextQ = tn->prevQ;
	}
	else
	{
	  n->nextQ = lcn;
	}
	n->prevQ = n->nextQ->prevQ;
      }
    }
    n->nextQ->prevQ = n;
    n->prevQ->nextQ = n;
  }
  n->state = WLZLS_NOD_OPEN;
#ifdef WLZ_LS_DEBUG
  (void )fprintf(stderr, "WlzLSNodQInsertNod() %d %d\n", debugHDepth, debugQDepth);
#endif
}

/*!
* \return	Key for node queue hash table.
* \ingroup	WlzFeatures
* \brief	Computes a slot key for the node queue hash table using
* 		the nodes position.
* \param	q		Node queue.
* \param	p		Position of node in (or to be in) the queue.
*/
static unsigned int		WlzLSNodQKey(
				  WlzLSNodQueue *q,
				  WlzIVertex3 p)
{
  unsigned int	k;

  k = p.vtX + 0x9e3779b9;
  k ^= p.vtY + 0x9e3779b9 + (k << 6) + (k >> 2);
  k ^= p.vtZ + 0x9e3779b9 + (k << 6) + (k >> 2);
  k ^= 0x9e3779b9 + (k << 6) + (k >> 2);
  k = k % q->nSlt;
  return(k);
}

/*!
* \return	Sort order integer.
* \ingroup	WlzFeatures
* \brief	Computes a sort order integer for the given pair of vertices
* 		and can (sort of) be thought of as p0 - p1.
* \param	p0		Position of 1st node.
* \param	p1		Position of 2nd node.
*/
static int			WlzLSPosCmp(
				  WlzIVertex3 p0,
				  WlzIVertex3 p1)
{
  int		c = 0;

  if(!(c = p0.vtZ - p1.vtZ))
  {
    if(!(c = p0.vtY - p1.vtY))
    {
      c = p0.vtX - p1.vtX;
    }
  }
  return(c);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Removes the domain covered by the new path from the current
* 		domain, label it and where the resulting domains are
* 		significant add them to the domain queue.
* \todo		This function remains unimplemented.
* \param	objQ		Domain queue.
* \param	obj		Current object domain.
* \param	path		Given path.
* \param	bDstMax		Maximum distance of any point in the object
* 				domain from the boundary.
*/
static WlzErrorNum		WlzLSSplitDom(
				  AlcHeap *objQ,
				  WlzObject *obj,
				  WlzObject *path,
				  int bDstMax)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WLZ_ERR_UNIMPLEMENTED;
  return(errNum);
}

