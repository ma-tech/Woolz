#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDomainFill3D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzDomainFill3D.c
* \author       Bill Hill
* \date         September 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief        Functions to fill holes in 3D voxel domain objects.
* \ingroup      WlzDomainOps
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Wlz.h>

static void			WlzDomFill3DDoBound2D(
				  WlzBoundList **bnd,
				  WlzDomain dom);

/*!
* \return	New Woolz object without holes or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Fills the holes in the given object's domain (which are by
* 		definition not connected to the outside). When the given
* 		object's domain has more than one component part, the
* 		object should first be labeled, this function should then be
* 		called for each of the labeled parts and then the union of
* 		the filled domains should be formed.
* \param	srcObj			Given 3D domain object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 			*WlzDomainFill3D(
  				  WlzObject *srcObj,
    				  WlzErrorNum *dstErr)
{
  int		nPln = 0;
  WlzObject	*bndObj = NULL,
  		*filObj = NULL,
  		*gvnObj = NULL,
		*sedObj = NULL,
		*shlObj = NULL;
  WlzPixelV	zeroV;
  WlzValues 	nullVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nullVal.core = NULL;
  zeroV.type = WLZ_GREY_UBYTE;
  zeroV.v.ubv = 0;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    gvnObj = WlzMakeMain(srcObj->type, srcObj->domain, nullVal,
    		         NULL, NULL, &errNum);
  }
  /* Erode the given object and then subtract this from the given object
   * to create a then shell 1 voxel thick just inside the given objects's
   * domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*difObj = NULL,
		*dilObj = NULL,
    		*erdObj = NULL;

    erdObj = WlzAssignObject(
    	   WlzErosion(gvnObj, WLZ_26_CONNECTED, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      difObj = WlzAssignObject(
      	       WlzDiffDomain(gvnObj, erdObj, &errNum), NULL);
    }
    /*
    if(errNum == WLZ_ERR_NONE)
    {
      dilObj = WlzAssignObject(
      	       WlzDilation(difObj, WLZ_26_CONNECTED, &errNum), NULL);
    }
    */
    if(errNum == WLZ_ERR_NONE)
    {
      WlzIBox3	clipBox;

      /* Clip the dilated shell domain to make sure it stays within the
       * bounding box of the given object then all planes will align. */
      clipBox.xMin = gvnObj->domain.p->kol1;
      clipBox.yMin = gvnObj->domain.p->line1;
      clipBox.zMin = gvnObj->domain.p->plane1;
      clipBox.xMax = gvnObj->domain.p->lastkl;
      clipBox.yMax = gvnObj->domain.p->lastln;
      clipBox.zMax = gvnObj->domain.p->lastpl;
      shlObj = WlzAssignObject(
	       /*
	       WlzClipObjToBox3D(dilObj, clipBox, &errNum), NULL);
	       */
	       WlzClipObjToBox3D(difObj, clipBox, &errNum), NULL);
    }
    (void )WlzFreeObj(difObj);
    /* (void )WlzFreeObj(dilObj); */
    (void )WlzFreeObj(erdObj);
  }
  /* Make sure that the bounding box of the thin shell domain fits it and
   * that it's first and last planes have interrvals. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(shlObj->domain.p, NULL);
  }
  /* Create a value table for the shell object with values set to zero. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	val;
    WlzObjectType tType;

    tType = WlzGreyTableType(WLZ_GREY_TAB_INTL, WLZ_GREY_UBYTE, NULL);
    val.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
				  shlObj->domain.p->plane1,
				  shlObj->domain.p->lastpl,
				  zeroV, NULL, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      int	p;

      nPln = shlObj->domain.p->lastpl - shlObj->domain.p->plane1 + 1;
      shlObj->values = WlzAssignValues(val, NULL);
#ifdef _OPENMP
#pragma omp parallel for shared(shlObj)
#endif
      for(p = 0; p < nPln; ++p)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  WlzDomain dom2;
	  WlzErrorNum errNum2;

	  dom2 = shlObj->domain.p->domains[p];
	  if(dom2.core)
	  {
	    WlzValues val2;
	    WlzObject *shlObj2;
	    shlObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2, nullVal, NULL, NULL,
				&errNum2);
	    if(errNum2 == WLZ_ERR_NONE)
	    {
	      val2.i = WlzMakeIntervalValues(tType, shlObj2, zeroV, &errNum2);
	      /* WlzMakeIntervalValues() sets all values to zero. */
	    }
	    if(errNum2 == WLZ_ERR_NONE)                          
	    {
	      shlObj->values.vox->values[p] = WlzAssignValues(val2, NULL);
	    }
	    (void )WlzFreeObj(shlObj2);
	    if(errNum2 == WLZ_ERR_NONE)
	    {
#ifdef _OPENMP
#pragma omp critical
	      {
#endif
		if((errNum == WLZ_ERR_NONE) && (errNum2 != WLZ_ERR_NONE))
		{
		  errNum = errNum2;
		}
#ifdef _OPENMP
	      }
#endif
	    }
	  }
	}
      }
    }
  }
  /* Compute the (plane-wise) boundary list for the given object. */
  if(errNum == WLZ_ERR_NONE)
  {
    bndObj = WlzObjToBoundary(gvnObj, 0, &errNum);
  }
  /* Sweep down through the boundary object setting the values of
   * those voxels in the shell object to a non zero value when they
   * correspond to top level boundaries. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		p;

#ifdef _OPENMP
#pragma omp parallel for shared(bndObj,shlObj)
#endif
    for(p = 0; p < nPln; ++p)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	WlzDomain bDom2;

	bDom2 = bndObj->domain.p->domains[p];
	if(bDom2.core)
	{
	  WlzDomain iDom2;
	  WlzValues iVal2;
	  WlzObject *iObj2 = NULL;
	  WlzGreyValueWSpace *gVWSp = NULL;
	  WlzErrorNum errNum2 = WLZ_ERR_NONE;

	  iDom2 = shlObj->domain.p->domains[p];
	  iVal2 = shlObj->values.vox->values[p];
	  iObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, iDom2, iVal2, NULL, NULL,
			      &errNum2);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    gVWSp = WlzGreyValueMakeWSp(iObj2, &errNum2);
	  }
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    WlzBoundList *bnd,
			 *bnd2;

	    bnd2 = bDom2.b;
	    for(bnd = bnd2; bnd != NULL; bnd = bnd->next)
	    {
	      if(bnd->poly != NULL)
	      {
		WlzPolygonDomain *ply;

		ply = bnd->poly;
		if(ply)
		{
		  int	i;
		  WlzIVertex2 *vtx;

		  vtx = ply->vtx;
		  for(i = 0; i < ply->nvertices; ++i)
		  {
		    WlzGreyValueGet(gVWSp, 0, vtx[i].vtY, vtx[i].vtX);
		    *(gVWSp->gPtr[0].ubp) = 255;
		  }
		}
	      }
	    }
	  }
	  else
	  {
#ifdef _OPENMP
#pragma omp critical
	    {
#endif
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = errNum2;
	      }
#ifdef _OPENMP
	    }
#endif
	  }
	  (void )WlzFreeObj(iObj2);
	  WlzGreyValueFreeWSp(gVWSp);
	}
      }
    }
  }
  /* Threshold the shell object, throwing away all but where the voxels
   * are set to create a seed domain. Then remove the value table from
   * the shell object and free it as it's no longer needed. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*tObj = NULL;
    WlzPixelV	tV;

    tV.type = WLZ_GREY_UBYTE;
    tV.v.ubv = 1;
    tObj = WlzAssignObject(
    	   WlzThreshold(shlObj, tV, WLZ_THRESH_HIGH, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      sedObj = WlzAssignObject(
      	       WlzMakeMain(tObj->type, tObj->domain, nullVal,
			   NULL, NULL, &errNum), NULL);
    }
    (void )WlzFreeObj(tObj); tObj = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      tObj = WlzAssignObject(
             WlzMakeMain(shlObj->type, shlObj->domain, nullVal,
			 NULL, NULL, &errNum), NULL);
    }
    (void )WlzFreeObj(shlObj); shlObj = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      shlObj = tObj;
      tObj = NULL;
    }
    (void )WlzFreeObj(tObj);
#ifdef WLZ_DOMOMAINFILL3D_DEBUG
    {
      FILE	*fP;
      fP = fopen("debug-shlObj-00.wlz", "w");
      (void )WlzWriteObj(fP, shlObj);
      (void )fclose(fP);
    }
#endif
  }
  /* Label the shell domain using 26-connectivity in 3D and then
   * keep only those component objects which intersect the seed domain.
   * Then free the shell and seed domains replacing the shell domain
   * with the union of the intersecting labeled component objects.
   * Finaly free the intersecting component objects, keeping only the
   * new shell domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i,
		j,
    		nCSObj = 0;
    WlzObject	**csObj = NULL;
    const int 	maxCSObj = 1000;

    errNum = WlzLabel(shlObj, &nCSObj, &csObj, maxCSObj, 1,
		      WLZ_26_CONNECTED);
    if(errNum == WLZ_ERR_NONE)
    {
      for(i = 0; i < nCSObj; ++i)
      {
        if(!WlzHasIntersection(csObj[i], sedObj, &errNum))
	{
	  (void )WlzFreeObj(csObj[i]);
	  csObj[i] = NULL;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Squeeze out any NULL objects reseting their number.*/
      for(i = 0, j = 0; i < nCSObj; ++i)
      {
        if(csObj[i])
	{
	  csObj[j++] = csObj[i];
	}
      }
      nCSObj = j;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzObject	*iObj = NULL,
      		*uObj = NULL;

      uObj = WlzUnionN(nCSObj, csObj, 0, &errNum);
      iObj = WlzAssignObject(
               WlzIntersect2(uObj, shlObj, &errNum),  NULL);
      (void )WlzFreeObj(uObj);
      (void )WlzFreeObj(shlObj);
      shlObj = WlzAssignObject(iObj, NULL);
#ifdef WLZ_DOMOMAINFILL3D_DEBUG
      {
	FILE	*fP;
	fP = fopen("debug-shlObj-01.wlz", "w");
	(void )WlzWriteObj(fP, shlObj);
	(void )fclose(fP);
      }
#endif
    }
    if(csObj)
    {
      for(i = 0; i < nCSObj; ++i)
      {
        (void )WlzFreeObj(csObj[i]);
      }
      (void )AlcFree(csObj);
    }
  }
  /* Sweep down through the boundary lists again removing boundaries that
   * do not intersect the new shell domain and creating a new filled
   * object from them. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		p;
    WlzDomain	filDom;

    filDom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN, 
		shlObj->domain.p->plane1, shlObj->domain.p->lastpl,
                shlObj->domain.p->line1, shlObj->domain.p->lastln,
                shlObj->domain.p->kol1, shlObj->domain.p->lastkl,
		&errNum);
#ifdef _OPENMP
#pragma omp parallel for shared(bndObj,shlObj)
#endif
    for(p = 0; p < nPln; ++p)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	WlzDomain bDom2;
	
	bDom2 = bndObj->domain.p->domains[p];
	if(bDom2.core)
	{
	  WlzDomain 	sDom2;
	  WlzObject	*fObj2 = NULL;
	  WlzErrorNum	errNum2 = WLZ_ERR_NONE;

	  sDom2 = shlObj->domain.p->domains[p];
	  WlzDomFill3DDoBound2D(&(bDom2.b), sDom2);
	  fObj2 = WlzBoundToObj(bDom2.b, WLZ_SIMPLE_FILL, &errNum2);
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    if(fObj2)
	    {
	      filDom.p->domains[p] = WlzAssignDomain(fObj2->domain, NULL);
	    }
	  }
	  else
	  {
#ifdef _OPENMP
#pragma omp critical
	    {
#endif
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = errNum2;
	      }
#ifdef _OPENMP
	    }
#endif
	  }
	  (void )WlzFreeObj(fObj2);
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzStandardPlaneDomain(filDom.p, NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      filObj = WlzMakeMain(srcObj->type, filDom, nullVal,
    		         NULL, NULL, &errNum);
    }
  }
  (void )WlzFreeObj(bndObj);
  (void )WlzFreeObj(gvnObj);
  (void )WlzFreeObj(shlObj);
  if((errNum != WLZ_ERR_NONE) && (filObj != NULL))
  {
    (void )WlzFreeObj(filObj);
    filObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(filObj);
}

/*!
* \ingroup	WlzDomainOps
* \brief	Remove boundaries that do not intersect the given pixel
* 		domain by just skiping them in the boundary linked lists.
* \param	bnd		Pointer to the given boundary list pointer.
* \param	dom		Given pixel domain.
*/
static void			WlzDomFill3DDoBound2D(
				  WlzBoundList **bnd,
				  WlzDomain dom)
{
  if(bnd && *bnd)
  {
    do
    {
      int	isn = 0;
      WlzPolygonDomain *ply;
      WlzBoundList **dwnBnd,
		   **nxtBnd;

      dwnBnd = &((*bnd)->down);
      nxtBnd = &((*bnd)->next);
      if((ply = (*bnd)->poly) != NULL)
      {
	int	i;
	WlzIVertex2 *vtx;

	vtx = ply->vtx;
	for(i = 0; i < ply->nvertices; ++i)
	{
	  isn = WlzInsideDomain2D(dom.i, vtx[i].vtY,  vtx[i].vtX, NULL);
	  if(isn)
	  {
	    break;
	  }
	}
      }
      if(isn)
      {
	WlzDomFill3DDoBound2D(dwnBnd, dom);
        bnd = nxtBnd;
      }
      else
      {
        *bnd = *nxtBnd;
      }
    }
    while(*bnd);
  }
}
