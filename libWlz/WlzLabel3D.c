#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLabel3D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzLabel3D.c
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
* \brief	3D labeling (segmention).
* \ingroup	WlzBinaryOps
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Wlz.h>

/*!
* \return	A compund array object containing the labeled object
* 		components of the given object.
* \ingroup	WlzBinaryOps
* \brief	Labels (segments) a 3D domain object into connected component
* 		objects using their connectivity. 
* \param	gObj		Given object to be labeled.
* \param	maxObj		Maximum number of objects to be found in any
* 				plane.
* \param	ignLn		Ignore objects within a plane which have
* 				\f$\geq\f$ the given number of lines.
* \param	con		The connectivity to use in 3D.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject			*WlzLabel3D(
				  WlzObject *gObj,
				  int maxObj,
				  int ignLn,
				  WlzConnectType con,
				  WlzErrorNum *dstErr)
{
  int		nPln,		/* Number of planes in the input object. */
  		nObjs = 0,	/* Number of 3D objects found by labeling. */
		totNFrg = 0,	/* Total number of 2D fragment objects found
				 * by labeling the planes. */
		maxFrgPln = 0;  /* Maximum number of 2D fragment objects in
				 * any plane. */
  size_t	objBufSz = 0;
  WlzObject	*lObj = NULL;	/* The return compound object. */
  WlzConnectType con2 = WLZ_0_CONNECTED,
  		 con3 = WLZ_0_CONNECTED;
  AlcUFTree	*uft = NULL;	/* Union-find tree for grouping the framents
  				 * into 3D objects. */
  WlzObject	**frgObjBuf = NULL; /* Buffer big enough for all the fragment
  				     * objects of any single plane for each
				     * thread. */
  WlzObject	***frgTbl = NULL; /* Objects in each plane, indexed by plane
  				   * then fragment within plane. */
  int		*nFrgTbl = NULL, /* Table of the number of fragments per plane,
  				  * indexed by plane. */
		*cNFrgTbl = NULL, /* Cumulative version of nFrgTbl, simplifies
			           * OpenMP code. */
  		*frgToObjTb = NULL; /* Look up table from fragment index to
				       3D labeled object index. */
  WlzPixelV	bgdV;
  WlzCompoundArray *objs = NULL;
  WlzValues	nulVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nulVal.core = NULL;
  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((gObj->values.core != NULL) &&
         (gObj->values.core->type != WLZ_VOXELVALUETABLE_GREY))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    switch(con)
    {
      case WLZ_4_CONNECTED:  /* FALLTHROUGH */
      case WLZ_6_CONNECTED:
        con2 = WLZ_4_CONNECTED;
        con3 = WLZ_0_CONNECTED;
	break;
      case WLZ_8_CONNECTED:  /* FALLTHROUGH */
      case WLZ_18_CONNECTED: /* FALLTHROUGH */
      case WLZ_26_CONNECTED:
        con2 = WLZ_8_CONNECTED;
        con3 = WLZ_8_CONNECTED;
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  /* Create a fragment objects table. */
  if(errNum == WLZ_ERR_NONE)
  {
    nPln = gObj->domain.p->lastpl - gObj->domain.p->plane1 + 1;
    if(((frgTbl = (WlzObject ***)
                  AlcCalloc(nPln, sizeof(WlzObject **))) == NULL) ||
       ((nFrgTbl = (int *)AlcCalloc(nPln, sizeof(int))) == NULL) ||
       ((cNFrgTbl = (int *)AlcCalloc(nPln, sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Label each plane populating the fragment objects table. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		p;
    WlzDomain	*doms;

    doms = gObj->domain.p->domains;
#ifdef _OPENMP
#pragma omp parallel for shared(nFrgTbl,frgTbl,doms,nulVal,maxObj,ignLn,con2)
#endif
    for(p = 0; p < nPln; ++p)
    {
      if(errNum == WLZ_ERR_NONE)
      {
        WlzDomain	dom2;
	WlzErrorNum	errNum2 = WLZ_ERR_NONE;

	dom2 = doms[p];
	if(dom2.core != NULL)
	{
	  WlzObject	*obj2;

	  obj2 = WlzAssignObject(
	         WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2, nulVal,
		 	     NULL, NULL, &errNum2), NULL);
	  if(obj2)
	  {
	    int		nFrg = 0;
	    WlzObject	**frgObj = NULL;

	    errNum2 = WlzLabel(obj2, &nFrg, &frgObj, maxObj, ignLn, con2);
	    WlzFreeObj(obj2);
	    if(errNum2 == WLZ_ERR_NONE)
	    {
	      if(frgObj != NULL)
	      {
	        if(nFrg > 0)
		{
		  nFrgTbl[p] = nFrg;
		  frgTbl[p] = frgObj;
	        }
		else
		{
		  AlcFree(frgObj);
		}
	      }
	    }
	  }
	}
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
  if(errNum == WLZ_ERR_NONE)
  {
    int		p;

    for(p = 0; p < nPln; ++p)
    {
      int	nFrg;

      nFrg = nFrgTbl[p];
      cNFrgTbl[p] = totNFrg;
      if(nFrg > maxFrgPln)
      {
	maxFrgPln = nFrg;
      }
      totNFrg += nFrg;
    }
  }
  if(errNum == WLZ_ERR_NONE) 
  {
    if(totNFrg < 1)
    {
      /* No fragments but given a domain object. */
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if(totNFrg == 1)
    {
      /* Just have a single fragment so no need for complex connectivity
       * code, just make a compound object with the single fragment. */
      int	p;
      WlzObject	*obj2 = NULL;
      WlzDomain dom;

      for(p = 0; p < nPln; ++p)
      {
	if(nFrgTbl[p])
	{
	  obj2 = *(frgTbl[p]);
	  break;
	}
      }
      dom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				 gObj->domain.p->plane1 + p,
				 gObj->domain.p->lastpl + p,
				 obj2->domain.i->line1,
				 obj2->domain.i->lastln,
				 obj2->domain.i->kol1,
				 obj2->domain.i->lastkl,
				 &errNum);
      objs = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, 1, NULL,
				  WLZ_3D_DOMAINOBJ, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        objs->n = 1;
	objs->o[0] = WlzAssignObject(
		     WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, nulVal, NULL, NULL,
				 &errNum), NULL);
      }
    }
    else
    {
      /* Have many fragments and so need to establish their connectivity. */
      /* Create a find-union tree for all the fragment objects. */
      if(errNum == WLZ_ERR_NONE)
      {
	if((uft = AlcUFTreeNew(totNFrg, totNFrg)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      /* Walk down through the planes checking for conectivity between
       * the fragments of the current plane and the previous plane updating
       * the union-find tree to build the connectivity. */
      if(errNum == WLZ_ERR_NONE)
      {
	int		p;		/* Current plane */

#ifdef _OPENMP
#pragma omp parallel for shared(nFrgTbl,cNFrgTbl,con3,uft)
#endif
	for(p = 1; p < nPln; ++p)
	{

	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	q,		/* Previous plane. */
		    nFrgP,	/* Number of fgragment objects on plane p */
		    nFrgQ; 	/* Number of fgragment objects on plane q */
	    WlzErrorNum errNum2 = WLZ_ERR_NONE;

	    q = p - 1;
	    nFrgP = nFrgTbl[p];
	    nFrgQ = nFrgTbl[q];
	    if((nFrgP > 0) && (nFrgQ > 0))
	    {
	      int	  fP;	/* Fragment within the current plane p */
	      WlzObject **frgP,	/* Array of fragment objects on plane p */
			**frgQ;	/* Array of fragment objects on plane q */

	      frgP = frgTbl[p];
	      frgQ = frgTbl[q];
	      for(fP = 0; (errNum2 == WLZ_ERR_NONE) && (fP < nFrgP); ++fP)
	      {
		int	fQ;
		WlzObject *fObjP;

		if(con3 == WLZ_0_CONNECTED)
		{
		  fObjP = WlzAssignObject(frgP[fP], NULL);
		}
		else
		{
		  fObjP = WlzDilation(frgP[fP], con3, &errNum2);
		}
		if((fObjP != NULL) &&
		   (fObjP->type == WLZ_2D_DOMAINOBJ) &&
		   (fObjP->domain.core != NULL))
		{
		  for(fQ = 0; (errNum2 == WLZ_ERR_NONE) && (fQ < nFrgQ); ++fQ)
		  {
		    WlzObject *fObjQ;

		     fObjQ = frgQ[fQ];
		    /* Avoid calling WlzHasIntersection() if possible, saves
		     * a significant amount of CPU time. */
		    if((fObjQ != NULL) &&
		       (fObjQ->type == WLZ_2D_DOMAINOBJ) &&
		       (fObjQ->domain.core != NULL) &&
		       (fObjP->domain.i->kol1 < fObjQ->domain.i->lastkl) &&
		       (fObjQ->domain.i->kol1 < fObjP->domain.i->lastkl) &&
		       (fObjP->domain.i->line1 < fObjQ->domain.i->lastln) &&
		       (fObjQ->domain.i->line1 < fObjP->domain.i->lastln) &&
		       WlzHasIntersection(fObjQ, fObjP, &errNum2))
		    {
#ifdef _OPENMP
#pragma omp critical
		      {
#endif
			AlcUFTreeUnion(uft, cNFrgTbl[q] + fQ, cNFrgTbl[p] + fP);
#ifdef _OPENMP
		      }
#endif
		    }
		  }
		}
		(void )WlzFreeObj(fObjP);
	      }
	    }
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
      /* The union-find tree now holds the connectivity of all the fragments.
       * Create a simple table which maps the fragment index to the 3D labeled
       * object index. */
      if(errNum == WLZ_ERR_NONE)
      {
	nObjs = uft->nCmp;
	if((frgToObjTb = (int *)AlcMalloc(totNFrg * sizeof(int))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  int	f,
		    o = 0;

	  /* First create mapping from the union-find tree components to 3D
	   * labeled object index. */
	  for(f = 0; f < totNFrg; ++f)
	  {
	    if(uft->pr[f] == f)
	    {
	      frgToObjTb[f] = o++;
	    }
	  }
	  /* Now assign the 3D labeled object index to each fragment. */
	  for(f = 0; f < totNFrg; ++f)
	  {
	    frgToObjTb[f] = frgToObjTb[AlcUFTreeFind(uft, f)];
	  }
	}
      }
      /* Finished with the union-find tree so free it. */
      AlcUFTreeFree(uft);
      /* Create a buffer large enough to hold all fragments on any plane (for
       * each thread if using OpenMP) and a compound array object for the
       * labeled 3D objects initialise them with the bounding box of the given
       * object. The plane domains will be fixed later. */
      if(errNum == WLZ_ERR_NONE)
      {
	objBufSz = maxFrgPln * sizeof(WlzObject *);
	if((frgObjBuf = (WlzObject **)AlcMalloc(objBufSz)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	objs = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, nObjs, NULL,
				    WLZ_3D_DOMAINOBJ, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	int		i;

	for(i = 0; (errNum == WLZ_ERR_NONE) && (i < nObjs); ++i)
	{
	  WlzDomain dom;

	  dom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				     gObj->domain.p->plane1,
				     gObj->domain.p->lastpl,
				     gObj->domain.p->line1,
				     gObj->domain.p->lastln,
				     gObj->domain.p->kol1,
				     gObj->domain.p->lastkl,
				     &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    objs->o[i] = WlzAssignObject(
			 WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, nulVal, NULL, NULL,
				     &errNum), NULL);
	  }
	}
      }
      /* Use the fragment to object index table to allocate the fragments
       * within each plane to the appropriate object.*/
      if(errNum == WLZ_ERR_NONE)
      {
	int		p;		/* Current plane */

	for(p = 0; (errNum == WLZ_ERR_NONE) && (p < nPln); ++p)
	{
	  int	fP,		/* Fragment within current plane p */
		  nFrgP,
		  frgIdxP; 	/* Index of the first fragment on plane p */
	  WlzObject **frgP;

	  frgP = frgTbl[p];
	  nFrgP = nFrgTbl[p];
	  frgIdxP = cNFrgTbl[p];
	  for(fP = 0; (errNum == WLZ_ERR_NONE) && (fP < nFrgP); ++fP)
	  {
	    if(frgP[fP])
	    {
	      /* Put all fragments on this plane which are part of the same
	       * labeled 3D object into a buffer. Form their union and the
	       * add the domain to the 3D object. */
	      int	  fP1,
		    objIdx,
		    frgBufIdx;
	      WlzObject *obj2 = NULL;

	      frgBufIdx = 0;
	      objIdx = frgToObjTb[frgIdxP + fP];
	      for(fP1 = fP; fP1 < nFrgP; ++fP1)
	      {
		if(objIdx == frgToObjTb[frgIdxP + fP1])
		{
		  frgObjBuf[frgBufIdx++] = frgP[fP1];
		  frgP[fP1] = NULL;
		}
	      }
	      if(frgBufIdx > 0)
	      {
		obj2 = WlzAssignObject(
		       WlzUnionN(frgBufIdx, frgObjBuf, 0, &errNum), NULL);
	      }
	      if(obj2)
	      {
		objs->o[objIdx]->domain.p->domains[p] =
		    WlzAssignDomain(obj2->domain, NULL);
		(void )WlzFreeObj(obj2);
	      }
	      for(fP1 = 0; fP1 < frgBufIdx; ++fP1)
	      {
		(void )WlzFreeObj(frgObjBuf[fP1]);
	      }
	    }
	  }
	}
      }
    }
  }
  /* For each of the labeled objects, standardise the plane domain, set the
   * voxel size, values and possibly background for the objects. */
  if((errNum == WLZ_ERR_NONE) && (totNFrg > 0) && (gObj->values.core != NULL))
  {
    bgdV = WlzGetBackground(gObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      int		i;
      WlzPlaneDomain *gPDom;

      gPDom = gObj->domain.p;
#ifdef _OPENMP
#pragma omp parallel for shared(objs)
#endif
      for(i = 0; i < nObjs; ++i)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  WlzPlaneDomain *nPDom;
	  WlzErrorNum    errNum2 = WLZ_ERR_NONE;
	  WlzValues	nVal;

	  nVal.core = NULL;
	  nPDom = objs->o[i]->domain.p;
	  nPDom->voxel_size[0] = gPDom->voxel_size[0];
	  nPDom->voxel_size[1] = gPDom->voxel_size[1];
	  nPDom->voxel_size[2] = gPDom->voxel_size[2];
	  if((errNum2 == WLZ_ERR_NONE) && (gObj->values.core != NULL))
	  {
	    int		p;
	    WlzValues	gVal;

	    gVal = gObj->values;
	    nVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					   nPDom->plane1, nPDom->lastpl,
					   bgdV, NULL, &errNum2);
	    if(errNum2 == WLZ_ERR_NONE)
	    {
	      objs->o[i]->values = WlzAssignValues(nVal, NULL);
	      for(p = nPDom->plane1; p <= nPDom->lastpl; ++p)
	      {
		int		gP,
				  nP;

		nP = p - nPDom->plane1;
		gP = p - gPDom->plane1;
		nVal.vox->values[nP] = WlzAssignValues(
				       gVal.vox->values[gP], NULL);
	      }
	    }
	  }
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    errNum2 = WlzStandardPlaneDomain(nPDom, nVal.vox);
	  }
	  if(errNum2 != WLZ_ERR_NONE)
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
  /* Free the buffers, array of fragments per plane and the fragment
   * table. */
  AlcFree(frgObjBuf);
  AlcFree(frgToObjTb);
  if(frgTbl)
  {
    int		p;

    for(p = 0; p < nPln; ++p)
    {
      AlcFree(frgTbl[p]);
    }
    AlcFree(frgTbl);
  }
  AlcFree(nFrgTbl);
  AlcFree(cNFrgTbl);
  if(errNum == WLZ_ERR_NONE)
  {
    lObj = (WlzObject *)objs;
  }
  else if(objs)
  {
    WlzFreeObj((WlzObject *)objs);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(lObj);
}

