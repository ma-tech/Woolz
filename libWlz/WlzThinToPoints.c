#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzThinToPoints_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzThinToPoints.c
* \author       Bill Hill
* \date         March 2016
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
* \brief	Functions for extracting point locations by thinnning.
* \ingroup	WlzFeatures
*/

#include <Wlz.h>

typedef struct _WlzThinToPointsHeapEntry
{
  double	priority;                    
  int		thr;
  WlzLong	size;
  WlzObject	*obj;
} WlzThinToPointsHeapEntry;

static WlzErrorNum		WlzThinToPointsPopHeap(
				  AlcHeap *heap,
				  WlzObject *gObj,
				  int gThin,
				  int gStart,
				  int gInc,
				  WlzThresholdType gHiLo);
static WlzErrorNum		WlzThinToPointsSplitHeap(
				  AlcHeap *heap,
				  int dim,
				  int gThin,
				  int gInc,
				  WlzThresholdType gHiLo);
static WlzPoints		*WlzThinToPointsFromHeap(
				  AlcHeap *heap,
				  int dim,
				  WlzErrorNum *dstErr);

/*!
* \return	New Woolz points object.
* \ingroup	WlzFeatures
* \brief	Extracts a points object in which the points are the
*  		centres of thinned regions of the input domain object.
*  		If the object has grey values and gThin is set then
*  		the object is erroded by successive thresholding while
*  		collecting all disconnected regions or if the object
*  		has no values the object is simply labeled to give a
*  		collection of disconnected regions. The disconnected
*  		regions are then eroded (morphologicaly) again
*  		collecting disconnected regions.
* \param	gObj			Given domain object, with or
* 					without values.
* \param	gThin			Use object grey values.
* \param	gStart			Initial threshold value for grey
* 					thinning.
* \param	gInc			Increment for grey thinning.
* \param	dstErr			Destination error pointer, may be NUll.
*/
WlzObject			*WlzThinToPoints(
				  WlzObject *gObj,
				  int gThin,
				  int gStart,
				  int gInc,
				  WlzErrorNum *dstErr)
{
  int		dim = 0;
  WlzObject	*pObj = NULL;
  AlcHeap	*heap = NULL;
  WlzThresholdType hiLo = WLZ_THRESH_HIGH;
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
      case WLZ_2D_DOMAINOBJ:
        dim = 2;
	break;
      case WLZ_3D_DOMAINOBJ:
	dim = 3;
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(gObj->values.core == NULL)
    {
      gThin = 0;
    }
    if(gThin)
    {
      if(gInc == 0)
      {
	errNum = WLZ_ERR_PARAM_DATA;
      }
      else if(gInc < 0)
      {
        hiLo = WLZ_THRESH_LOW;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((heap = AlcHeapNew(sizeof(WlzThinToPointsHeapEntry),
                           1024, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzThinToPointsPopHeap(heap, gObj, gThin, gStart, gInc, hiLo);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzThinToPointsSplitHeap(heap, dim, gThin, gInc, hiLo);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain	dom;

    dom.pts = WlzThinToPointsFromHeap(heap, dim, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValues	nullVal;

      nullVal.core = NULL;
      pObj = WlzMakeMain(WLZ_POINTS, dom, nullVal, NULL, NULL, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	pObj = NULL;
        (void )WlzFreeDomain(dom);
      }
    }
  }
  AlcHeapFree(heap);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pObj);
}

/*!
* \return	Wolz error code.
* \ingroup	WlzFeatures
* \brief	Given an empty heap, populate it with the labeled objects
* 		resulting from the given object. If the given object has
* 		grey values then gThin will be set and the objct will be
* 		thresholded at the initial grey threshold.
* \param	heap			Given heap to be populated.
* \param	gObj			Given object.
* \param	gThin			Thin using grey values if non-zero.
* \param	gStart			Initial threshold level if gThin set.
* \param	gInc			Threshold increment if gThin set.
* \param	gHiLo			Threshold high or low (only used if
* 					gThin set).
*/
static WlzErrorNum		WlzThinToPointsPopHeap(
				  AlcHeap *heap,
				  WlzObject *gObj,
				  int gThin,
				  int gStart,
				  int gInc,
				  WlzThresholdType gHiLo)
{
  int		nLObj = 0;
  WlzObject	*tObj = NULL;
  WlzObject	**lObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gThin)
  {
    WlzPixelV	thrV;

    thrV.type = WLZ_GREY_INT;
    thrV.v.inv = gStart;
    tObj = WlzAssignObject(WlzThreshold(gObj, thrV, gHiLo, &errNum), NULL);
  }
  else
  {
    tObj = WlzAssignObject(gObj, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzLong	maxObj;
    WlzConnectType con;

    maxObj = WlzVolume(tObj, &errNum);
    if(errNum == WLZ_ERR_NONE) 
    {
      if(maxObj == 0)
      {
        errNum = WLZ_ERR_DOMAIN_DATA;
      }
      else
      {
	if(tObj->type == WLZ_2D_DOMAINOBJ)
	{
	  con = WLZ_8_CONNECTED;
	}
	else
	{
	  con = WLZ_26_CONNECTED;
	}
	maxObj = 8 + (maxObj / 2);
	errNum = WlzLabel(tObj, &nLObj, &lObj, maxObj, 0, con);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    for(i = 0; i < nLObj; ++i)
    {
      WlzLong	sz;
      WlzThinToPointsHeapEntry ent;

      sz = WlzVolume(lObj[i], &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	if(sz > 0)
	{
	  ent.obj = lObj[i];
	  ent.thr = gStart;
	  ent.size = sz;
	  ent.priority = sz;
	  if(AlcHeapInsertEnt(heap, &ent))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	else
	{
	  (void )WlzFreeObj(lObj[i]);
	}
      }
    }
    AlcFree(lObj);
  }
  (void )WlzFreeObj(tObj);
  return(errNum);
}

/*!
* \return	Wolz error code.
* \ingroup	WlzFeatures
* \brief	Given a heap with an initial population of objects, continue
* 		to split the heap objects using their grey values (if gThin
* 		set) and/or erosion. When objects can be thinned no more
* 		before they disappear their heap priority is reduced below
* 		1.0.
* \param	heap			Given heap to be populated.
* \param	dim			Either 2 or 3 representing the object
* 					dimension.
* \param	gThin			Use grey values.
* \param	gInc			Threshold increment.
* \param	gHiLo			Threshold high or low.
*/
static WlzErrorNum		WlzThinToPointsSplitHeap(
				  AlcHeap *heap,
				  int dim,
				  int gThin,
				  int gInc,
				  WlzThresholdType gHiLo)
{
  int		more = 1;
  WlzConnectType con = WLZ_26_CONNECTED;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(dim == 2)
  {
    con = WLZ_8_CONNECTED;
  }
  do
  {
    WlzThinToPointsHeapEntry *top;

    top = AlcHeapTop(heap);
    if(top->priority < 1.0)
    {
      more = 0;
    }
    else
    {
      int	sz = 0;
      WlzObject	*nObj = NULL;
      WlzThinToPointsHeapEntry ent;

      ent = *top;
      AlcHeapEntFree(heap);
      if(ent.obj->values.core)
      {
	/* Object has grey values, threshold away. */
	WlzPixelV	thrV;

	ent.thr += gInc;
	thrV.v.inv = ent.thr;
	thrV.type = WLZ_GREY_INT;
	nObj = WlzAssignObject(
	       WlzThreshold(ent.obj, thrV, gHiLo, &errNum), NULL);
        
      }
      else
      {
	/* Object has no grey values, erode away. */

        nObj = WlzAssignObject(WlzErosion(ent.obj, con, &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	sz = WlzVolume(nObj, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(sz < 1)
	{
	  /* Object has been eroded/thresholded away. */
	  if(ent.obj->values.core)
	  {
	    /* Remove object values. */
	    WlzObject	*dObj;
	    WlzValues	nullVal;

	    nullVal.core = NULL;
	    dObj = WlzAssignObject(
	           WlzMakeMain(ent.obj->type, ent.obj->domain, nullVal,
		               NULL, NULL, &errNum), NULL);
	    (void )WlzFreeObj(ent.obj);
	    ent.obj = dObj;
	  }
	  else
	  {
	    /* If object had no values then mark it done by setting priority
	     * < 1.0. */
	    ent.priority = 0.5;
	  }
	  /* Put original object (with no values) back onto the heap. */
	  if(errNum == WLZ_ERR_NONE)
	  {
	    (void )AlcHeapInsertEnt(heap, &ent);
	  }
	}
	else if(sz == ent.size)
	{
	  /* Object domain unchanged, so just put it back. */
	  if(AlcHeapInsertEnt(heap, &ent))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	else
	{
	  /* Object remains after erosion/thresholding, but is reduced. */
	  int		maxObj,
			nLObj = 0;
	  WlzObject 	**lObj = NULL;

	  maxObj = 8 + (sz / 2);
	  errNum = WlzLabel(nObj, &nLObj, &lObj, maxObj, 0, con);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i;

	    WlzFreeObj(ent.obj);
	    for(i = 0; (errNum == WLZ_ERR_NONE) && (i < nLObj); ++i)
	    {
	      WlzLong	sz;

	      sz = WlzVolume(lObj[i], &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(sz > 0)
		{
		  ent.obj = lObj[i];
		  ent.size = sz;
		  ent.priority = sz;
		  if(AlcHeapInsertEnt(heap, &ent))
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		else
		{
		  (void )WlzFreeObj(lObj[i]);
		}
	      }
	    }
            AlcFree(lObj);
	  }
	}
      }
      (void )WlzFreeObj(nObj);
    }
  } while((errNum == WLZ_ERR_NONE) && more);
  return(errNum);
}

/*!
* \return	New points domain or NULL on error.
* \ingroup	WlzFeatures
* \brief	Creates a points domain with each point at the centroid
* 		of the objects on the given heap.
* \param	heap			Given heap.
* \param	dim			Dimension, either 2 or 3.
* \param	dstErr			Destination error point, may be NULL.
*/
static WlzPoints 		*WlzThinToPointsFromHeap(
				  AlcHeap *heap,
				  int dim,
				  WlzErrorNum *dstErr)
{
  int		nPts;
  WlzVertexP	nullVP;
  WlzPoints	*pts = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nullVP.v = NULL;
  nPts = heap->nEnt;
  pts = WlzMakePoints((dim == 2)? WLZ_POINTS_2I: WLZ_POINTS_3I,
		      0, nullVP, nPts, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    for(i = 0; i < nPts; ++i)
    {
      WlzThinToPointsHeapEntry *top;

      top = AlcHeapTop(heap);
      if((top == NULL) || (top->obj == NULL))
      {
        break;
      }
      if(dim == 2)
      {
	WlzIVertex2 *p;
	WlzDVertex2 c;

	p = &(pts->points.i2[i]);
	(void )WlzStandardIntervalDomain(top->obj->domain.i);
	c = WlzCentreOfMass2D(top->obj, 1, NULL, &errNum);
	WLZ_VTX_2_NINT(*p, c);
      }
      else
      {
	WlzIVertex3 *p;
	WlzDVertex3 c;

	p = &(pts->points.i3[i]);
	(void )WlzStandardPlaneDomain(top->obj->domain.p,
	                              top->obj->values.vox);
	c = WlzCentreOfMass3D(top->obj, 1, NULL, &errNum);
	WLZ_VTX_3_NINT(*p, c);
      }
      (void )WlzFreeObj(top->obj);
      AlcHeapEntFree(heap);
    }
    pts->nPoints = i;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pts);
}
