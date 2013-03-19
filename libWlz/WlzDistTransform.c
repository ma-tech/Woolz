#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDistTransform_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzDistTransform.c
* \author       Konstantinos Liakos, Bill Hill
* \date         October 2004
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
* \brief	Distance transform functions which calculate the distance
* 		of every pixel/voxel in a foreground object from a reference
* 		object.
* \ingroup	WlzMorphologyOps
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static WlzObject 		*WlzDistSample(
				  WlzObject *obj,
				  int dim,
				  double scale,
    			          WlzErrorNum *dstErr);

/*!
* \return	Distance object which shares the given foreground object's
*		domain and has integer distance values, null on error.
* \ingroup	WlzMorphologyOps
* \brief	Computes the distance of every pixel/voxel in the foreground
* 		object from the reference object.
*
*		A distance transform maps all position within a  forground
*		domain to their distances from a reference domain.
*		The distance transforms implemented within this function
*		use efficient morphological primitives.
*		
*		Given two domains,
*		\f$\Omega_r\f$ the reference domain and \f$\Omega_f\f$
*		the domain specifying the region of interest,
*		a domain with a thin shell \f$\Omega_i\f$
*		is iteratively expanded from it's initial domain
*		corresponding to the reference domain \f$\Omega_r\f$.
*		At each iteration
*		\f$\Omega_i\f$ is dilated and clipped
*		by it's intersection with \f$\Omega_f\f$ until \f$\Omega_i\f$
*		becomes the null domain \f$\emptyset\f$.
*		At each iteration the current distance is recorded in a value
*		table which
*		covers the domain \f$\Omega_f\f$.
*
*		An octagonal distance scheme may be used in which
*		the distance metric is alternated between 4 and 8
*		connected for 2D and 6 and 26 connectivities in 3D.
*		See: G. Borgefors. "Distance Transformations in Arbitrary
*		Dimensions" CVGIP 27:321-345, 1984.
*
* 		An approximate Euclidean distance transform may be computed
* 		by: Scaling the given foreground and reference objects using
* 		the given approximation scale parameter, dilating the
* 		reference domain using a sphere with a radius having the same
* 		value as the scale parameter and then finaly sampling the
* 		scaled distances.
* \param	forObj			Foreground object.
* \param	refObj			Reference object.
* \param	dFn			Distance function which must be
*					appropriate to the dimension of
*					the foreground and reference objects.
* \param	dParam			Parameter required for distance
* 					function. Currently only
* 					WLZ_APX_EUCLIDEAN_DISTANCE requires a
* 					parameter. In this case the parameter
* 					is the approximation scale.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzDistanceTransform(WlzObject *forObj, WlzObject *refObj,
				   WlzDistanceType dFn, double dParam,
				   WlzErrorNum *dstErr)
{
  int 		idP,
		lastP,
		dim = 0,
  		notDone = 1;
  double	scale;
  WlzObject	*tmpObj,
		*sObj = NULL,
		*sForObj = NULL,
		*sRefObj = NULL,
		*dilObj = NULL,
  		*dstObj = NULL,
		*difObj = NULL,
  		*curItrObj = NULL;
  WlzObject 	*bothObj[2];
  WlzDomain	*difDoms;
  WlzPixelV 	dstV,
  		bgdV;
  WlzValues 	*difVals;
  WlzAffineTransform *tr = NULL;
  WlzConnectType con = WLZ_0_CONNECTED;
  WlzObjectType dstGType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzValues 	difVal,
  		dstVal,
		nullVal;
  /* By defining WLZ_DIST_TRANSFORM_ENV these normalization parameters
   * are read from the environment. This is useful for optimization. */
#ifndef WLZ_DIST_TRANSFORM_ENV
  const
#endif /* ! WLZ_DIST_TRANSFORM_ENV */
  /* These normalizarion factors have been choosen to minimize the sum of
   * squares of the deviation of the distance values from Euclidean values
   * over a radius 100 circle or sphere, where the distances are computed
   * from the circumference of the sphere towards it's centre. The values
   * were established by experiment. */
  double	nrmDist4 =  0.97,
		nrmDist6 =  0.91,
		nrmDist8 =  1.36,
		nrmDist18 = 1.34,
		nrmDist26 = 1.60;

#ifdef WLZ_DIST_TRANSFORM_ENV
  double	val;
  char		*envStr;

  if(((envStr = getenv("WLZ_DIST_TRANSFORM_NRMDIST4")) != NULL) &&
     (sscanf(envStr, "%lg", &val) == 1))
  {
    nrmDist4 = val;
  }
  if(((envStr = getenv("WLZ_DIST_TRANSFORM_NRMDIST6")) != NULL) &&
     (sscanf(envStr, "%lg", &val) == 1))
  {
    nrmDist6 = val;
  }
  if(((envStr = getenv("WLZ_DIST_TRANSFORM_NRMDIST8")) != NULL) &&
     (sscanf(envStr, "%lg", &val) == 1))
  {
    nrmDist8 = val;
  }
  if(((envStr = getenv("WLZ_DIST_TRANSFORM_NRMDIST18")) != NULL) &&
     (sscanf(envStr, "%lg", &val) == 1))
  {
    nrmDist18 = val;
  }
  if(((envStr = getenv("WLZ_DIST_TRANSFORM_NRMDIST26")) != NULL) &&
     (sscanf(envStr, "%lg", &val) == 1))
  {
    nrmDist26 = val;
  }
#endif /* WLZ_DIST_TRANSFORM_ENV */
  scale = dParam;
  nullVal.core = NULL;
  /* Check parameters. */
  if((forObj == NULL) || (refObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(((forObj->type != WLZ_2D_DOMAINOBJ) &&
           (forObj->type != WLZ_3D_DOMAINOBJ)) ||
          ((refObj->type != WLZ_POINTS) &&
	   (refObj->type != forObj->type)))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((forObj->domain.core == NULL) || (refObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bgdV.type = WLZ_GREY_INT;
    bgdV.v.inv = 0;
    dstV.type = WLZ_GREY_DOUBLE;
    dstV.v.dbv = 0.0;
    switch(forObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	switch(dFn)
	{
	  case WLZ_4_DISTANCE: /* FALLTHROUGH */
	  case WLZ_8_DISTANCE: /* FALLTHROUGH */
	  case WLZ_OCTAGONAL_DISTANCE: /* FALLTHROUGH */
	  case WLZ_APX_EUCLIDEAN_DISTANCE:
	    dim = 2;
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	switch(dFn)
	{
	  case WLZ_6_DISTANCE:  /* FALLTHROUGH */
	  case WLZ_18_DISTANCE: /* FALLTHROUGH */
	  case WLZ_26_DISTANCE: /* FALLTHROUGH */
	  case WLZ_OCTAGONAL_DISTANCE: /* FALLTHROUGH */
	  case WLZ_APX_EUCLIDEAN_DISTANCE:
	    dim = 3;
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
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
      case WLZ_4_DISTANCE:
        con = WLZ_4_CONNECTED;
	break;
      case WLZ_6_DISTANCE:
        con = WLZ_6_CONNECTED;
	break;
      case WLZ_8_DISTANCE:
        con = WLZ_8_CONNECTED;
	break;
      case WLZ_18_DISTANCE:
        con = WLZ_18_CONNECTED;
	break;
      case WLZ_26_DISTANCE:
        con = WLZ_26_CONNECTED;
	break;
      case WLZ_OCTAGONAL_DISTANCE:
        con = (dim == 2)? WLZ_8_CONNECTED: WLZ_26_CONNECTED;
	break;
      case WLZ_APX_EUCLIDEAN_DISTANCE:
        con = (dim == 2)? WLZ_8_CONNECTED: WLZ_26_CONNECTED;
	if(scale < 1.0)
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	}
	break;
      case WLZ_EUCLIDEAN_DISTANCE:
	errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  /* Create scaled domains and a sphere domain for structual erosion if the
   * distance function is approximate Euclidean. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(dFn == WLZ_APX_EUCLIDEAN_DISTANCE)
    {
      tr = (dim == 2)?
	   WlzAffineTransformFromScale(WLZ_TRANSFORM_2D_AFFINE,
	                               scale, scale, 0.0, &errNum):
	   WlzAffineTransformFromScale(WLZ_TRANSFORM_3D_AFFINE,
	                               scale, scale, scale, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	tmpObj = WlzMakeMain(forObj->type, forObj->domain, nullVal,
			     NULL, NULL, &errNum);
	if(tmpObj)
	{
	  sForObj = WlzAssignObject(
	            WlzAffineTransformObj(tmpObj, tr,
		                          WLZ_INTERPOLATION_NEAREST,
		    			  &errNum), NULL);
	  (void )WlzFreeObj(tmpObj);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(refObj->type == WLZ_POINTS)
	{
	  sRefObj = WlzPointsToDomObj(refObj->domain.pts, scale, &errNum);
	}
	else /* type == WLZ_2D_DOMAINOBJ || type == WLZ_3D_DOMAINOBJ */
	{
	  tmpObj = WlzMakeMain(refObj->type, refObj->domain, nullVal,
			       NULL, NULL, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    sRefObj = WlzAssignObject(
		      WlzAffineTransformObj(tmpObj, tr,
					    WLZ_INTERPOLATION_NEAREST,
					    &errNum), NULL);
	  }
	}
	(void )WlzFreeObj(tmpObj);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	sObj = WlzAssignObject(
	       WlzMakeSphereObject(forObj->type, scale,
	                           0.0, 0.0, 0.0, &errNum), NULL);
      }
      (void )WlzFreeAffineTransform(tr);
    }
    else
    {
      sForObj = WlzAssignObject(
	        WlzMakeMain(forObj->type, forObj->domain, nullVal,
	  		    NULL, NULL, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	if(refObj->type == WLZ_POINTS)
	{
	  sRefObj = WlzPointsToDomObj(refObj->domain.pts, 1.0, &errNum);
	}
	else
	{
	  sRefObj = WlzAssignObject(
		    WlzMakeMain(refObj->type, refObj->domain, nullVal,
				NULL, NULL, &errNum), NULL);
	}
      }
    }
  }
  /* Create new values for the computed distances. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstGType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
    if(dim == 2)
    {
      dstVal.v = WlzNewValueTb(sForObj, dstGType, bgdV, &errNum);
    }
    else
    {
      dstVal.vox = WlzNewValuesVox(sForObj, dstGType, bgdV, &errNum);
    }
  }
  /* Create a distance object using the foreground object's domain and
   * the new distance values. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzMakeMain(sForObj->type, sForObj->domain, dstVal,
			 NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bothObj[0] = sForObj;
    errNum = WlzGreySetValue(dstObj, dstV);
  }
  /* Dilate the reference object while setting the distances in each
   * dilated shell. */
  while((errNum == WLZ_ERR_NONE) && notDone)
  {
    if(dFn == WLZ_APX_EUCLIDEAN_DISTANCE)
    {
      dstV.v.dbv += 1.0;
    }
    else
    {
      switch(con)
      {
	case WLZ_4_CONNECTED:
	  dstV.v.dbv += nrmDist4;
	  break;
	case WLZ_6_CONNECTED:
	  dstV.v.dbv += nrmDist6;
	  break;
	case WLZ_8_CONNECTED:
	  dstV.v.dbv += nrmDist8;
	  break;
	case WLZ_18_CONNECTED:
	  dstV.v.dbv += nrmDist18;
	  break;
	case WLZ_26_CONNECTED:
	  dstV.v.dbv += nrmDist26;
	  break;
        default:
	  errNum = WLZ_ERR_CONNECTIVITY_TYPE;
	  break;
      }
    }
    if(dFn == WLZ_APX_EUCLIDEAN_DISTANCE)
    {
      dilObj = WlzStructDilation(sRefObj, sObj, &errNum);
    }
    else
    {
      dilObj = WlzDilation(sRefObj, con, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(sForObj->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  curItrObj = WlzAssignObject(
	              WlzIntersect2(dilObj, sForObj, &errNum), NULL);
	  break;
        case WLZ_3D_DOMAINOBJ:
	  bothObj[1] = dilObj;
	  curItrObj = WlzAssignObject(
	              WlzIntersectN(2, bothObj, 1, &errNum), NULL);
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    (void)WlzFreeObj(dilObj);
    /* Create difference object for the expanding shell. */
    if(errNum == WLZ_ERR_NONE)
    {
      difObj = WlzDiffDomain(curItrObj, sRefObj, &errNum);
    }
    if((difObj == NULL) || WlzIsEmpty(difObj, &errNum))
    {
      notDone = 0;
    }
    else
    {
      /* Assign the distance object's values to the difference object
       * and set all it's values to the current distance. */
      if(errNum == WLZ_ERR_NONE)
      {
	switch(sForObj->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	    difObj->values = WlzAssignValues(dstObj->values, NULL);
	    errNum = WlzGreySetValue(difObj, dstV);
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    /* 3D is more complex than 2D: Need to create a temporary
	     * voxel valuetable and assign the individual 2D values. */
	    difVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					     difObj->domain.p->plane1,
					     difObj->domain.p->lastpl,
					     bgdV, NULL, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      difObj->values = WlzAssignValues(difVal, NULL);
	      difDoms = difObj->domain.p->domains;
	      difVals = difObj->values.vox->values;
	      idP = difObj->domain.p->plane1;
	      lastP = difObj->domain.p->lastpl;
	      while(idP <= lastP)
	      {
		if((*difDoms).core)
		{
		  dstVal = dstObj->values.vox->values[idP - 
						      dstObj->domain.p->plane1];
		  *difVals = WlzAssignValues(dstVal, NULL);
		}
		++idP;
		++difDoms;
		++difVals;
	      }
	      if(difObj->domain.p->lastpl > difObj->domain.p->plane1)
	      {
		errNum = WlzGreySetValue(difObj, dstV);
	      }
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
      (void )WlzFreeObj(sRefObj);
      sRefObj = WlzAssignObject(curItrObj, NULL);
      (void )WlzFreeObj(curItrObj);
    }
    (void )WlzFreeObj(difObj); difObj = NULL;
    if(dFn == WLZ_OCTAGONAL_DISTANCE)
    {
      /* Alternate connectivities for octagonal distance. */
      if(dim == 2)
      {
	con = (con == WLZ_4_CONNECTED)? WLZ_8_CONNECTED: WLZ_4_CONNECTED;
      }
      else /* dim == 3 */
      {
	con = (con == WLZ_6_CONNECTED)? WLZ_26_CONNECTED: WLZ_6_CONNECTED;
      }
    }
  }
  (void )WlzFreeObj(sObj);
  (void )WlzFreeObj(sForObj);
  (void )WlzFreeObj(sRefObj);
  (void )WlzFreeObj(curItrObj);
  if((errNum == WLZ_ERR_NONE) && (dFn == WLZ_APX_EUCLIDEAN_DISTANCE))
  {
    tmpObj = WlzDistSample(dstObj, dim, scale, &errNum);
    (void )WlzFreeObj(dstObj);
    dstObj = tmpObj;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(dstObj); dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	Sampled object.
* \ingroup	WlzMorphologyOps
* \brief	Samples the scaled distance object, interpolating the distance
* 		values.
* \param	obj			Given object to be sampled.
* \param	dim			Dimension either 2D or 3D.
* \param	scale			Scale factor.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzDistSample(WlzObject *obj, int dim, double scale,
    			        WlzErrorNum *dstErr)
{
  WlzObject	*sObj = NULL;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(scale < 1.0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    scale = 1.0 / scale;
    tr = (dim == 2)?
	 WlzAffineTransformFromScale(WLZ_TRANSFORM_2D_AFFINE,
				     scale, scale, 0.0, &errNum):
	 WlzAffineTransformFromScale(WLZ_TRANSFORM_3D_AFFINE,
				     scale, scale, scale, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sObj = WlzAffineTransformObj(obj, tr, WLZ_INTERPOLATION_LINEAR,
				 &errNum);
  }
  (void )WlzFreeAffineTransform(tr);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(sObj);
}
