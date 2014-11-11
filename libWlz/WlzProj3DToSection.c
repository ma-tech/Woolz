#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzProj3DToSection_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzProj3DToSection.c
* \author       Bill Hill
* \date         June 2011
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
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
* \brief 	Functions for projecting 3D domains onto sections.
* \ingroup	WlzTransform
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	New domain object.
* \ingroup	WlzTransform
* \brief	Projects the given 3D doimain onto a section using
* 		a given projection, set of mask domains and set of
* 		within section transforms and a section view.
* \param	gvnObj			Given 3D domain to be projected
* 					onto a section the section
* \param	nMask			Number of mask domain objects.
* \param	maskObj			Mask domain objects and mask objects
* 					that are NULL are assumed to represent
* 					the universal domain.
* \param	prjView			Projection view.
* \param	nPlnTr			Number of in section transforms,
* 					must be equal to the number of
* 					mask domains.
* \param	plnTrObj		Within-plane transforms any
* 				        within-plane transforms that are NULL
* 				        are assumed to represent an identity
* 				        transform.
* \param	secView			The section view transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzProj3DToSection(WlzObject *gvnObj,
				int nMask, WlzObject **maskObj,
				WlzThreeDViewStruct *prjView,
				int nPlnTr, WlzObject **plnTrObj,
				WlzThreeDViewStruct *secView,
				WlzErrorNum *dstErr)
{
  WlzObject	*plnObj = NULL,
  		*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check parameters. */
  if((gvnObj == NULL) || (prjView == NULL) || (secView == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((nMask < 1) || (nMask != nPlnTr) ||
          (maskObj == NULL || (plnTrObj == NULL)))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    int		idx;

    for(idx = 0; idx < nMask; ++idx)
    {
      if(maskObj[idx] != NULL)
      {
        if(maskObj[idx]->type != WLZ_3D_DOMAINOBJ)
	{
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
	}
	else if(maskObj[idx]->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	  break;
	}
      }
      if(plnTrObj[idx] != NULL)
      {
        if(plnTrObj[idx]->type != WLZ_MESH_TRANS)
	{
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
	}
	else if (plnTrObj[idx]->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx;

    for(idx = 0; idx < nMask; ++idx)
    {
      WlzObject *isnObj = NULL,
      		*prjObj = NULL,
		*prtObj = NULL;

      isnObj = WlzAssignObject(
               (maskObj[idx] == NULL)?
	       gvnObj:
               WlzIntersect2(gvnObj, maskObj[idx], &errNum), NULL);
      if((errNum == WLZ_ERR_NONE) && (WlzIsEmpty(isnObj, NULL) == 0))
      {
        WlzThreeDViewStruct *prjView1 = NULL;

        prjView1 = WlzMake3DViewStructCopy(prjView, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
          errNum = WlzInit3DViewStruct(prjView1, isnObj);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  prjObj = WlzProjectObjToPlane(isnObj, prjView1,
	  				WLZ_PROJECT_INT_MODE_NONE, 0, NULL,
					0.0, &errNum);
	}
        (void )WlzFree3DViewStruct(prjView1);
	if(errNum == WLZ_ERR_NONE)
	{
	  if(plnTrObj[idx] == NULL)
	  {
	    prtObj = prjObj;
	    prjObj = NULL;
	  }
	  else
	  {
	    prtObj = WlzMeshTransformObj(prjObj, plnTrObj[idx]->domain.mt,
					 WLZ_INTERPOLATION_NEAREST, &errNum);
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(plnObj == NULL)
	  {
	    plnObj = prtObj;
	    prtObj = NULL;
	  }
	  else
	  {
	    WlzObject *tmpObj;

	    tmpObj = WlzUnion2(plnObj, prtObj, &errNum);
	    (void )WlzFreeObj(plnObj);
	    plnObj = tmpObj;
	  }
	}
      }
      (void )WlzFreeObj(isnObj);
      (void )WlzFreeObj(prjObj);
      (void )WlzFreeObj(prtObj);
      if(errNum != WLZ_ERR_NONE)
      {
	(void )WlzFreeObj(plnObj);
	plnObj = NULL;
        break;
      }
    }
    if((plnObj != NULL) && (WlzIsEmpty(plnObj, NULL) != 0))
    {
      (void )WlzFreeObj(plnObj);
      plnObj = NULL;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (plnObj != NULL))
  {
    WlzThreeDViewStruct *secView1 = NULL;

    (void )WlzAssignObject(plnObj, NULL);
    secView1 = WlzMake3DViewStructCopy(secView, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzInit3DViewStruct(secView1, plnObj);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rtnObj = Wlz3DViewTransformObj(plnObj, secView1, &errNum);
    }
    (void )WlzFree3DViewStruct(secView1);
  }
  (void )WlzFreeObj(plnObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*
cat cs14-billeye-warped.wlz embryo_2_3D_body.wlz | \
WlzIntersect > j1.wlz
Wlz3DGetProjection -b embryo_2_WM_left_proj.bib j1.wlz > j2.wlz
WlzMeshTransformObj -m embryo_2_WM_left_body_tr.wlz j2.wlz | \
Wlz3DViewTransformObj -b embryo_2_WM_left.bib > cs14-billeye-warped-wm.wlz
*/
