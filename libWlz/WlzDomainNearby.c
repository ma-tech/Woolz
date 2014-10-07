#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDomainNearby_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzDomainNearby.c
* \author       Bill Hill
* \date         October 2014
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
* \brief	Functions for computing the portion of a domain which
* 		is nearby (all pixels/voxels less than a given distance
* 		from) given locations.
* \ingroup	WlzMorphologyOps
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	A new domain in which all pixels/voxels are nearby the given
* 		location(s). The object returned will either be of the same
* 		type as the reference object or an empty object if (all)
* 		the given location(s) are outside of the reference object's
* 		domain. On error a NULL pointer will be returned.
* \ingroup	WlzMorphologyOps
* \brief	
* \param	refObj		Reference object which must be of type
* 				WLZ_EMPTY_OBJ, WLZ_2D_DOMAINOBJ or
* 				WLZ_3D_DOMAINOBJ.
* \param	nPos		Number of locations.
* \param	pos		Locations which are either of type
* 				WlzDVertex2 if the given object is of type
* 				WLZ_2D_DOMAINOBJ or WlzDVertex3 if the
* 				given object is of type WLZ_3D_DOMAINOBJ.
* \param	dFn		The connectivity to use when establishing
* 				distance. This must be one of:
*				<ul>
*				<li>WLZ_OCTAGONAL_DISTANCE
*				<li>WLZ_4_DISTANCE
*				<li>WLZ_6_DISTANCE
*				<li>WLZ_8_DISTANCE
*				<li>WLZ_18_DISTANCE
*				<li>WLZ_26_DISTANCE
*				</ul>
*				Distances must be appropriate for the given
*				object type although all are acceptable
*				for an empty object.
* \param	dMax		Maximum distance for a pixel/voxel to be
* 				considered nearby.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject			*WlzDomainNearby(
				  WlzObject *refObj,
				  int nPos, WlzVertexP pos,
				  WlzDistanceType dFn,
				  double dMax,
				  WlzErrorNum *dstErr)
{
  int		dim = 0;
  WlzObject	*rtnObj = NULL;
  WlzObject	**parObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(refObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(refObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(refObj->type)
    {
      case WLZ_EMPTY_OBJ:
        dim = 0;
	break;
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
    switch(dFn)
    {
      case WLZ_OCTAGONAL_DISTANCE:
        break;
      case WLZ_4_DISTANCE:  /* FALLTHROUGH */
      case WLZ_8_DISTANCE:
	if(dim == 3)
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	}
        break;
      case WLZ_6_DISTANCE:  /* FALLTHROUGH */
      case WLZ_18_DISTANCE: /* FALLTHROUGH */
      case WLZ_26_DISTANCE:
	if(dim == 2)
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	}
        break;
      default:
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((parObj = (WlzObject **)
    		 AlcCalloc(nPos, sizeof(WlzObject *))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i,
    		in = 0;

    for(i = 0; (errNum == WLZ_ERR_NONE) && (i < nPos); ++i)
    {
      WlzIVertex3 p;

      switch(dim)
      {
        case 2:
	  p.vtZ = 0;
	  WLZ_VTX_2_NINT(p, pos.d2[i]);
	  in = WlzInsideDomain2D(refObj->domain.i,
	                         p.vtY, p.vtX, &errNum);
	  break;
	case 3:
	  WLZ_VTX_3_NINT(p, pos.d3[i]);
	  in = WlzInsideDomain3D(refObj->domain.p,
	                         p.vtZ, p.vtY, p.vtX, &errNum);
	  break;
	default:
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(!in)
	{
	  parObj[i] = WlzMakeEmpty(&errNum);
	}
	else
	{
	  WlzObject *dObj = NULL,
	  	    *pObj = NULL,
		    *tObj = NULL;

	  pObj = WlzAssignObject(
	         WlzMakeSinglePixelObject(refObj->type, p.vtX, p.vtY, p.vtZ,
		                          &errNum), NULL);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dObj = WlzAssignObject(
	           WlzDistanceTransform(refObj, pObj, dFn, 0.0, dMax,
		                        &errNum), NULL);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzPixelV thrV;

	    thrV.type = WLZ_GREY_INT;
	    thrV.v.inv = 1;
	    tObj = WlzAssignObject(
	           WlzThreshold(dObj, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    parObj[i] = WlzAssignObject(
	    	        WlzUnion2(tObj, pObj, &errNum), NULL);
	  }
	  (void )WlzFreeObj(dObj);
	  (void )WlzFreeObj(pObj);
	  (void )WlzFreeObj(tObj);
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rtnObj = WlzUnionN(nPos, parObj, 0, &errNum);
  }
  if(parObj)
  {
    int		i;

    for(i = 0; i < nPos; ++i)
    {
      (void )WlzFreeObj(parObj[i]);
    }
    AlcFree(parObj);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

