#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzOffsetDist_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzOffsetDist.c
* \author       Bill Hill
* \date         October 2015
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2015],
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
* \brief	Functions to compute offset distance objects.
* \ingroup	WlzFeatures
*/
#include <Wlz.h>

/*!
* \return	New object with minimum distances.
* \ingroup	WlzFeatures
* \brief	Computes an object with domain and values that are the
*               set of minimum distances between the two given objects.
*		An equidistant boundary is computed between the domains
*		of the two given objects, within maxDist of each object's
*		domain and within the convex hull of the union of the
*		two given object's domains. This function will probably
*		only be useful where one of the objects tends to track
*		the other.
* \param	obj0			First given object.
* \param	obj1			Second given object.
* \param	maxDist			Maximum distance for offset. This
* 					is used to compute a distance object,
* 					large distances will significantly
* 					increase the processing time.
* 					If zero the maximum distance will
* 					be determined by the convex hull.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzOffsetDist(
				  WlzObject *obj0,
				  WlzObject *obj1,
				  int maxDist,
				  WlzErrorNum *dstErr)
{
  WlzObject	*o[2];
  WlzObject	*eObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((o[0] = obj0) == NULL) ||
     ((o[1] = obj1) == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((o[0]->type != o[1]->type) ||
          ((o[0]->type != WLZ_2D_DOMAINOBJ) &&
           (o[0]->type != WLZ_3D_DOMAINOBJ)))
  {
      errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((o[0]->domain.core == NULL) ||
	  (o[1]->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  /* Compute distance transforms of the two given objects out to a given
   * maximum distance and then using these distances the equi-distant
   * domain between these two objects. The values of the eqi-distant object
   * are those of the distance between the objects.*/
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    WlzObject	*sObj = NULL,
    		*cObj = NULL;
    WlzObject	*dObj[2] = {NULL};

    /* Create structuring element with which to dilate the given object
     * domains(by maxDist). */
    if(maxDist > 0)
    {
      sObj = WlzAssignObject(
	     WlzMakeSphereObject(o[0]->type, maxDist, 0, 0, 0, &errNum), NULL);
    }
    /* Create domain or convex hull of the union of the two given object
     * domains. */
    if(errNum == WLZ_ERR_NONE)
    {
      WlzObject	*uObj = NULL,
      		*xObj = NULL;

      uObj = WlzAssignObject(
             WlzUnionN(2, o, 0, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
        xObj = WlzAssignObject(
	       WlzObjToConvexHull(uObj, &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        cObj = WlzAssignObject(
	       WlzConvexHullToObj(xObj, o[0]->type, &errNum), NULL);
      }
      (void )WlzFreeObj(xObj);
      (void )WlzFreeObj(uObj);
    }
    /* Dilate the two given objects and find the intersection of the
     * dilated domains with each other and the convex hull computed
     * above. Within his domain compute the distances. */
    if(errNum == WLZ_ERR_NONE)
    {
      for(i = 0; i < 2; ++i)
      {
	WlzObject *rObj = NULL;

	if(maxDist > 0)
	{
	  WlzObject *tObj;

	  tObj = WlzAssignObject(
		 WlzStructDilation(o[i], sObj, &errNum), NULL);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    rObj = WlzAssignObject(
		   WlzIntersect2(tObj, cObj, &errNum), NULL);
	  }
	  (void )WlzFreeObj(tObj);
	}
	else
	{
	  rObj = WlzAssignObject(cObj, NULL);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  dObj[i] = WlzAssignObject(
		    WlzDistanceTransform(rObj, o[!i],
					 WLZ_OCTAGONAL_DISTANCE,
					 0.0, maxDist, &errNum), NULL);
	}
        (void )WlzFreeObj(rObj);
        if(errNum == WLZ_ERR_NONE)
	{
	  WlzPixelV bgdV;

	  bgdV.type = WLZ_GREY_INT;
	  bgdV.v.inv = maxDist;
	  errNum = WlzSetBackground(dObj[i], bgdV);
	}
        if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
    /* Find the domain which is equi-distant from the two given objects,
     * within the xDist range and within the convex hull of the union of
     * the two given object's domains. */
    (void )WlzFreeObj(sObj); sObj = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      int	empty = 0;
      WlzLong	vol = 0;
      WlzObject *qObj = NULL,
      		*tObj = NULL;

      qObj = WlzAssignObject(
             WlzImageArithmetic(dObj[0], dObj[1], WLZ_BO_EQ, 0, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	WlzPixelV thrV;

	thrV.type = WLZ_GREY_INT;
	thrV.v.inv = 1;
        tObj = WlzAssignObject(
	       WlzThreshold(qObj, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
      }
      /* Check that the eqi-distant domain is of a reasonable size ie has
       * a area or volume greater than half the maximum distance. */
      if(errNum == WLZ_ERR_NONE)
      {
        vol = WlzVolume(tObj, &errNum);
	if((maxDist / 2) >= vol)
	{
	  empty = 1;
	}
      }
      if((errNum == WLZ_ERR_NONE) && !empty)
      {
	WlzObject *mObj;
	WlzPixelV tmpV;

	tmpV.type = WLZ_GREY_INT;
	tmpV.v.inv = 0;
        mObj = WlzAssignObject(
	       WlzGreyTemplate(dObj[0], tObj, tmpV, &errNum), NULL);
        if(errNum == WLZ_ERR_NONE)
	{
	  tmpV.v.inv = 1;
	  eObj = WlzThreshold(mObj, tmpV, WLZ_THRESH_HIGH, &errNum);
	}
	(void )WlzFreeObj(mObj);
      }
      (void )WlzFreeObj(tObj);
      (void )WlzFreeObj(qObj);
      if((errNum == WLZ_ERR_NONE) && !empty)
      {
	WlzLong vol;

	vol = WlzVolume(eObj, &errNum);
	if((maxDist / 2) >= vol)
	{
	  empty = 1;
	}
      }
      if((errNum == WLZ_ERR_NONE) && empty)
      {
        (void )WlzFreeObj(eObj);
	eObj = WlzMakeEmpty(&errNum);
      }
    }
    (void )WlzFreeObj(cObj);
    (void )WlzFreeObj(sObj);
    (void )WlzFreeObj(dObj[0]);
    (void )WlzFreeObj(dObj[1]);
  }
  return(eObj);
}
