#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzInteriority_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzInteriority.c
* \author       Bill Hill
* \date         January 2015
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
* \brief	Functions to compute an interiority scores.
* \ingroup	WlzFeatures
*/
#include <Wlz.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double			WlzInteriorityPrv(
  				  WlzObject *dis,
				  WlzObject *tstObj,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzInteriorityCompDisObj(
				  WlzObject *refObj,
				  WlzErrorNum *dstErr);

/*!
* \return	Array of interiority scores with a score for each of the
* 		test objects or NULL on error.
* \ingroup	WlzFeatures
* \brief	Computes an interiority score for each test object with
* 		respect to the reference object. See WlzInteriority().
		This function computes the distance transform for the
		reference domain once and then uses this for each of the
		test domains.
		If not NULL the returned array of scores should be freed
		using AlcFree().
* \param	refObj			Given reference object.
* \param 	nTstObj			Number of test objects.
* \param	tstObjs			Array of test objects.
* \param	dstErr			Destination error pointer, may be NULL.
*/
double				*WlzInteriorityN(
				  WlzObject *refObj,
				  int nTstObj,
				  WlzObject **tstObjs,
				  WlzErrorNum *dstErr)
{
  WlzObject	*dis = NULL;
  double	*scores = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check given objects. */
  if(refObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((refObj->type != WLZ_2D_DOMAINOBJ) &&
          (refObj->type != WLZ_3D_DOMAINOBJ) &&
	  (refObj->type != WLZ_EMPTY_OBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((refObj->type != WLZ_EMPTY_OBJ) &&
          (refObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(nTstObj < 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(tstObjs == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    int 	i;

    for(i = 0; (errNum == WLZ_ERR_NONE) && (i < nTstObj); ++i)
    {
      if(tstObjs[i] == NULL)
      {
        errNum = WLZ_ERR_OBJECT_NULL;
      }
      if(tstObjs[i]->type != WLZ_EMPTY_OBJ)
      {
	if(tstObjs[i]->type != refObj->type)
	{
	  errNum = WLZ_ERR_OBJECT_TYPE;
	}
	else if(tstObjs[i]->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
      }
    }
  }
  /* Allocate the scores array. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((scores = (double *)AlcCalloc(nTstObj, sizeof(double))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Compute the boundary distance object. */
  if((errNum == WLZ_ERR_NONE) && (refObj->type != WLZ_EMPTY_OBJ))
  {
    dis = WlzAssignObject(
          WlzInteriorityCompDisObj(refObj, &errNum), NULL);
  }
  /* For each test object compute it's interiority score with respect
   * to the distance object. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i = 0; i < nTstObj; ++i)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	WlzErrorNum errNum2 = WLZ_ERR_NONE;

	scores[i] = WlzInteriorityPrv(dis, tstObjs[i], &errNum2);
	if(errNum2 != WLZ_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp critical
	  {
#endif
            errNum = errNum2;
#ifdef _OPENMP
	  }
#endif
	}
      }
    }
  }

  (void )WlzFreeObj(dis);
  if(errNum != WLZ_ERR_NONE)
  {
    AlcFree(scores);
    scores = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(scores);
}

/*!
* \return	Interiority score for the test object.
* \ingroup	WlzFeatures
* \brief	Computes an interiority score for the test object with
* 		respect to the reference object. The interiority score
* 		is the mean distance of those elements of the intersection
* 		of the reference and test object from the boundary of the
* 		reference object, ie:
*		\f[
     		s = t \cap r
		b = t \neg t^-
		i = \frac{1}{|s|} \sum_{s}{D(b, s)}
		\f]
		where \f$t\f$ is a test object, \f$t^-\f$ is the eroded
		test object, \f$|x|\f$ is the cardinality of \f$x\f$ and
		\f$D\f$ is the distance operator.
* \param	refObj			Given reference object.
* \param	tstObj			Array of test objects.
* \param	dstErr			Destination error pointer, may be NULL.
*/
double				WlzInteriority(
				  WlzObject *refObj,
				  WlzObject *tstObj,
				  WlzErrorNum *dstErr)
{
  int 		empty = 0;
  WlzObject	*dis = NULL;
  double	score = 0.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check given objects. */
  if((refObj == NULL) || (tstObj== NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((refObj->type == WLZ_EMPTY_OBJ) || (tstObj->type == WLZ_EMPTY_OBJ))
  {
    empty = 1;
  }
  else if((refObj->type != WLZ_2D_DOMAINOBJ) &&
          (refObj->type != WLZ_3D_DOMAINOBJ) &&
	  (refObj->type != WLZ_EMPTY_OBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((refObj->domain.core == NULL) || (tstObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(!empty)
  {
    /* Compute the boundary distance object. */
    dis = WlzAssignObject(
             WlzInteriorityCompDisObj(refObj, &errNum), NULL);
  }
  /* Compute the interiority score with respect to the distance object. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(!empty)
    {
      score = WlzInteriorityPrv(dis, tstObj, &errNum);
    }
  }
  (void )WlzFreeObj(dis);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(score);
}

/*!
* \return	Interiority score for the test object.
* \ingroup	WlzFeatures
* \brief	Computes an interiority score for the test object with
* 		respect to the given distance object. See WlzInteriority().
* \param	disObj			Given distance object.
* \param	tstObjs			Array of test objects.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static double			WlzInteriorityPrv(
  				  WlzObject *disObj,
				  WlzObject *tstObj,
				  WlzErrorNum *dstErr)
{
  double	score = 0.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tstObj->type != WLZ_EMPTY_OBJ)
  {
    WlzObject	*isn = NULL;

    isn = WlzAssignObject(
    	  WlzIntersect2(disObj, tstObj, &errNum), NULL);
    if((errNum == WLZ_ERR_NONE) &&
       (isn != NULL) && (WlzIsEmpty(isn, NULL) == 0))
    {
      WlzPixelV	bgd;
      WlzObjectType gtt;
      WlzObject *ist = NULL,
      		*isv = NULL;

      bgd.v.inv = 0;
      bgd.type = WLZ_GREY_INT;
      gtt = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
      isv = WlzNewObjectValues(isn, gtt, bgd, 0, bgd, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        ist = WlzAssignObject(
	      WlzGreyTransfer(isv, disObj, &errNum), NULL);
      }
      (void )WlzFreeObj(isv);
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzGreyStats(ist, NULL, NULL, NULL, NULL, NULL, &score, NULL,
			    &errNum);
      }
      (void )WlzFreeObj(ist);
      
    }
    (void )WlzFreeObj(isn);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(score);
}

/*!
* \return	Distance object or null on error.
* \ingroup	WlzFeatures
* \brief	Computes a distance object where all distances are from
* 		the boundary of the given reference object.
* \param	refObj			Given reference object, known to
* 					be valid.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzInteriorityCompDisObj(
				  WlzObject *refObj,
				  WlzErrorNum *dstErr)
{
  WlzObject	*dis = NULL,
  		*bndObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bndObj = WlzAssignObject(WlzBoundaryDomain(refObj, &errNum), NULL);
  if(errNum == WLZ_ERR_NONE)
  {
    dis = WlzDistanceTransform(refObj, bndObj, WLZ_OCTAGONAL_DISTANCE,
			       0.0, 0.0, &errNum);
  }
  (void )WlzFreeObj(bndObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dis);
}
