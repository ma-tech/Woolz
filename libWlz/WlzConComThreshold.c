#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConComThreshold_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzConComThreshold.c
* \author       Bill Hill
* \date         March 2014
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
* \brief	Functions to perform connected component thresholding.
* \ingroup	WlzThreshold
*/
#include <stdlib.h>
#include <Wlz.h>

static int 			WlzConComThreshComp(
				  WlzObject *pObj,
				  WlzGreyValueWSpace *gVWSp,
				  WlzIVertex3 pos,
				  WlzThresholdType rHiLo,
	      			  float *dstVal,
				  WlzErrorNum *dstErr);
static int 			WlzConComThreshComp2D(
				  WlzObject *pObj,
				  WlzGreyValueWSpace *gVWSp,
				  WlzIVertex3 off,
	      			  float *dstMin,
				  float *dstMax,
				  float *dstSum,
				  WlzErrorNum *dstErr);
static int 			WlzConComThreshComp3D(
				  WlzObject *pObj,
				  WlzGreyValueWSpace *gVWSp,
				  WlzIVertex3 off,
	      			  float *dstMin,
				  float *dstMax,
				  float *dstSum,
				  WlzErrorNum *dstErr);
/*!
* \return	New (thresholded) Woolz object or NULL on error.
* \ingroup	WlzThreshold
* \brief	Performs connected component thresholding of the given object
* 		using the given 2D seed points.
* 		See also WlzConComThreshold3D(), WlzConComThreshold() and
* 		WlzThreshold().
* \param	gObj			Given object to threshold, this must
* 					be either an empty object or a 2D
* 					domain object with values.
* \param	nSeeds			Number of seed positions supplied.
* \param	seeds			The given seed positions.
* \param	tHiLo			Threshold mode parameter.
* \param	xtr			Percent extra to add to threshold
* 					value.
* \param	rad			Radius of region around seed in which
* 					to establish threshold value.
* \param	rHiLo			How to choose threshold from the
* 					region around each seed.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 			*WlzConComThreshold2D(
				  WlzObject *gObj,
		                  int nSeeds,
				  WlzIVertex2 *seeds,
				  WlzThresholdType tHiLo,
				  int xtr,
				  double rad,
				  WlzThresholdType rHiLo,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzVertexP	seedsP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  seedsP.i2 = seeds;
  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_EMPTY_OBJ: /* FALLTHROUGH */
      case WLZ_2D_DOMAINOBJ:
        rObj = WlzConComThreshold(gObj, nSeeds, WLZ_VERTEX_I2, seedsP,
	                          tHiLo, xtr, rad, rHiLo, &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New (thresholded) Woolz object or NULL on error.
* \ingroup	WlzThreshold
* \brief	Performs connected component thresholding of the given object
* 		using the given 3D seed points.
* 		See also WlzConComThreshold2D(), WlzConComThreshold() and
* 		WlzThreshold().
* \param	gObj			Given object to threshold, this must
* 					be an empty object or either a 2 or
* 					3D domain object with values.
* \param	nSeeds			Number of seed positions supplied.
* \param	seeds			The given seed positions.
* \param	tHiLo			Threshold mode parameter.
* \param	xtr			Percent extra to add to threshold
* 					value.
* \param	rad			Radius of region around seed in which
* 					to establish threshold value.
* \param	rHiLo			How to choose threshold from the
* 					region around each seed.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 			*WlzConComThreshold3D(
				  WlzObject *gObj,
		                  int nSeeds,
				  WlzIVertex3 *seeds,
				  WlzThresholdType tHiLo,
				  int xtr,
				  double rad,
				  WlzThresholdType rHiLo,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzVertexP	seedsP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  seedsP.i3 = seeds;
  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_EMPTY_OBJ: /* FALLTHROUGH */
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
        rObj = WlzConComThreshold(gObj, nSeeds, WLZ_VERTEX_I3, seedsP,
				  tHiLo, xtr, rad, rHiLo, &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New (thresholded) Woolz object or NULL on error.
* \ingroup	WlzThreshold
* \brief	Performs connected component thresholding of the given object
* 		using the given coordinates of seed points within the given
* 		object.
* 		See also WlzConComThreshold2D(), WlzConComThreshold3D() and
* 		WlzThreshold().
* \param	gObj			Given object to threshold, this may
* 					be an empty object or either a 2 or
* 					3D domain object with values.
* \param	nSeeds			Number of seed positions supplied.
* \param	seedType		
* \param	seeds			The given seed positions. If the given
* 					object is a 3D domain object the the
* 					seeds must be of type WlzIVertex3
* 					otherwise they may be WlzIVertex2 or
* 					WlzIVertex3.
* \param	tHiLo			Threshold mode parameter with
* 					  WLZ_THRESH_HIGH >= threshold value,
* 					  WLZ_THRESH_EQUAL == threshold value
* 					  and
* 					  WLZ_THRESH_LOW  <= threshold value.
* 					  Watch out this is slightly different
* 					  from WlzThreshold().
* \param	xtr			Percent extra to add to threshold value
* 					where
* 					\[
					v_{thresh} = v_{seed}
					             \frac{100 \pm xtr}{100}
 					\]
					with the sign being:
					\f$-\f$ for WLZ_THRESH_HIGH and
					\f$+\f$ for WLZ_THRESH_LOW. If the
					value of tHiLo is WLZ_THRESH_EQUAL
					then this parameter is ignored.
* \param	rad			Radius of region around seed in which
* 					to establish threshold value, 0.0
* 					implies a single voxel.
* \param	rHiLo			How to choose threshold from the
* 					region around each seed:
* 					WLZ_THRESH_HIGH - highest,
* 					WLZ_THRESH_EQUAL - mean and
* 					WLZ_THRESH_LOW - lowest.
* 					Unused if region radius is 0.0.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 			*WlzConComThreshold(
				  WlzObject *gObj,
		                  int nSeeds,
				  WlzVertexType seedType,
				  WlzVertexP seeds,
				  WlzThresholdType tHiLo,
				  int xtr,
				  double rad,
				  WlzThresholdType rHiLo,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL,			           /* Return object. */
  		*pObj = NULL;	   /* Region object, needs shifting to seed. */
  WlzConnectType con = WLZ_0_CONNECTED;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minLn = 1,		/* Minimum number of lines in a x-y
  					   plane. */
  		maxConCom = 10000;	/* Maximum number of components that
  					   can be labeled. */
  const double	minRad = 0.5;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(nSeeds <= 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_EMPTY_OBJ: /* FALLTHROUGH */
      case WLZ_2D_DOMAINOBJ:
	con = WLZ_8_CONNECTED;
        if((seedType != WLZ_VERTEX_I2) && (seedType != WLZ_VERTEX_I3))
	{
	  errNum = WLZ_ERR_PARAM_TYPE;
	}
	else if((rad > minRad) && (gObj->type == WLZ_2D_DOMAINOBJ))
	{
	  pObj = WlzAssignObject(
	         WlzMakeSphereObject(WLZ_2D_DOMAINOBJ,
		                     rad, 0, 0, 0, &errNum), NULL);
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	con = WLZ_26_CONNECTED;
        if(seedType != WLZ_VERTEX_I3)
	{
	  errNum = WLZ_ERR_PARAM_TYPE;
	}
	else if(rad > minRad)
	{
	  pObj = WlzAssignObject(
	         WlzMakeSphereObject(WLZ_3D_DOMAINOBJ,
		                     rad, 0, 0, 0, &errNum), NULL);
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeEmpty(&errNum);
    if(gObj->type != WLZ_EMPTY_OBJ)
    {
      if(gObj->domain.core == NULL)
      {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if(gObj->values.core == NULL)
      {
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else
      {
	WlzGreyValueWSpace *gVWSp = NULL;

	gVWSp = WlzGreyValueMakeWSp(gObj, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  int	s;

	  for(s = 0; s < nSeeds; ++s)
	  {
	    int		ist = 0;   /* Seed or seed region intersects object. */
	    WlzPixelV   tVal;	           /* Threshold value for this seed. */
	    WlzIVertex3 pos;			    /* Position of this seed */
	    WlzObject 	*tObj = NULL; 		       /* Thresholded object */

	    if(seedType == WLZ_VERTEX_I2)
	    {
	      pos.vtX = seeds.i2[s].vtX;
	      pos.vtY = seeds.i2[s].vtY;
	      pos.vtZ = 0;
	    }
	    else /* seedType == WLZ_VERTEX_I3 */
	    {
	      pos = seeds.i3[s];
	    }
	    /* Compute threshold value for this seed. */
	    if(pObj)
	    {
	      tVal.type = WLZ_GREY_FLOAT;
	      ist = WlzConComThreshComp(pObj, gVWSp, pos, rHiLo,
	      			        &(tVal.v.flv), &errNum);
	    }
	    else
	    {
	      WlzGreyValueGet(gVWSp, pos.vtZ, pos.vtY, pos.vtX);
	      if(gVWSp->bkdFlag == 0)
	      {
		ist = 1;
		tVal.type = gVWSp->gType;
		tVal.v = gVWSp->gVal[0];
		(void )WlzValueConvertPixel(&(tVal), tVal, WLZ_GREY_FLOAT);
	      }
	    }
	    if((errNum == WLZ_ERR_NONE) && ist)
	    {
	      switch(tHiLo)
	      {
	        case WLZ_THRESH_HIGH:
		  tVal.v.flv = tVal.v.flv * (100.0f - xtr) / 100.0f;
		  break;
		case WLZ_THRESH_LOW:
		  /* Make sure the seed is in it if lower, hence add 1. */
		  tVal.v.flv = (tVal.v.flv + 1.0f) * (100.0f + xtr) / 100.0f;
		  break;
		case WLZ_THRESH_EQUAL:
		  break;
	      }
#define WLZ_CONCOMTHRESHOLD_DEBUG
#ifdef WLZ_CONCOMTHRESHOLD_DEBUG
	      (void )fprintf(stderr,
	                     "WLZ_CONCOMTHRESHOLD_DEBUG %d %d,%d,%d %c %f\n",
			     s, pos.vtX, pos.vtY, pos.vtZ,
			     (tHiLo == WLZ_THRESH_HIGH)?'H':
			     (tHiLo == WLZ_THRESH_LOW)?'L':'E',
			     tVal.v.flv);
#endif /* WLZ_CONCOMTHRESHOLD_DEBUG */
	      /* Threshold the object. */
	      tObj = WlzAssignObject(
		     WlzThreshold(gObj, tVal, tHiLo, &errNum), NULL);
	    }
	    if((errNum == WLZ_ERR_NONE) && tObj && tObj &&
	       (tObj->type != WLZ_EMPTY_OBJ))

	    {
	      int	nC = 0;
	      WlzObject **cAry = NULL;

	      /* Label thresholded object. */
	      errNum = WlzLabel(tObj, &nC, &cAry, maxConCom, minLn, con);
	      if((errNum == WLZ_ERR_NONE) && (nC > 0))
	      {
		int	c;

		/* For each component of the thresholded object, check if
		 * this seed is inside it. Can only be inside one. */
		for(c = 0; c < nC; ++c)
		{
		  int	in;

		  if(pObj)
		  {
		    WlzObject *qObj;

		    qObj = WlzAssignObject(
		    	   WlzShiftObject(pObj,
				pos.vtX, pos.vtY, pos.vtZ,
				&errNum), NULL);
                    if(errNum == WLZ_ERR_NONE)
		    {
		      WlzObject *iObj;

		      iObj = WlzAssignObject(
		             WlzIntersect2(cAry[c], qObj, &errNum), NULL);
		      in = WlzIsEmpty(iObj, NULL) == 0;
		      (void )WlzFreeObj(iObj);
		    }
		    (void )WlzFreeObj(qObj);
		  }
		  else
		  {
		    in = WlzInsideDomain(cAry[c], pos.vtZ, pos.vtY, pos.vtX,
					 &errNum);
		  }
		  if((errNum == WLZ_ERR_NONE) & in)
		  {
		    WlzObject *uObj;

		    uObj = WlzAssignObject(
			   WlzUnion2(rObj, cAry[c], &errNum), NULL);
		    (void )WlzFreeObj(rObj);
		    rObj = uObj;
		    break;               /* Break out of the component loop as
		    			    the seed can only be in one
					    component. */
		  }
		  if(errNum != WLZ_ERR_NONE)
		  {
		    break;
		  }
		}
		for(c = 0; c < nC; ++c)
		{
		  (void )WlzFreeObj(cAry[c]);
		}
		(void )AlcFree(cAry);
	      }
	    }
	    (void )WlzFreeObj(tObj);
	    if(errNum != WLZ_ERR_NONE)
	    {
	      break;
	    }
	  }
	}
	WlzGreyValueFreeWSp(gVWSp);
      }
    }
    if(rObj)
    {
      if(errNum == WLZ_ERR_NONE)
      {
        WlzObject *tObj;

	tObj = WlzMakeMain(rObj->type, rObj->domain, rObj->values,
	                   rObj->plist, rObj->assoc, &errNum);
        (void )WlzFreeObj(rObj);
	rObj = tObj;
      }
    }
  }
  (void )WlzFreeObj(pObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Volume (or area) of intersection.
* \ingroup	WlzThreshold
* \brief	Given a sampling region object centred at the origin and
* 		a valid grey value workspace for an object. This function
* 		finds the min, max or mean value of the intersection if one
* 		exists.
* \param	pObj		Sampling region object.
* \param	gVWSp		Grey value workspace.
* \param	pos		Position of the sample.
* \param	rHiLo		Used to indicate max, min and mean.
* \param	dstVal		Destination pointer for the value.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static int 			WlzConComThreshComp(
				  WlzObject *pObj,
				  WlzGreyValueWSpace *gVWSp,
				  WlzIVertex3 pos,
				  WlzThresholdType rHiLo,
	      			  float *dstVal,
				  WlzErrorNum *dstErr)
{
  int		vol = 0;
  float		min = 0.0,
  		max = 0.0,
		sum = 0.0,
		val = 0.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(pObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(pObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gVWSp == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(pObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
      {
	vol = WlzConComThreshComp2D(pObj, gVWSp, pos,
				    &min, &max, &sum, &errNum);
	break;
      }
      case WLZ_3D_DOMAINOBJ:
	vol = WlzConComThreshComp3D(pObj, gVWSp, pos,
				    &min, &max, &sum, &errNum);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (vol > 0))
  {
    switch(rHiLo)
    {
      case WLZ_THRESH_HIGH:
	val = max;
        break;
      case WLZ_THRESH_EQUAL:
	val = sum / vol;
        break;
      case WLZ_THRESH_LOW:
	val = min;
        break;
    }
  }
  *dstVal = val;
  *dstErr = errNum;
  return(vol);
}

/*!
* \return	Volume (or area) of intersection.
* \ingroup	WlzThreshold
* \brief	Given a 2D sampling region object centred at the origin and
* 		a valid grey value workspace for an object. This function
* 		finds the min, max and sum value of the intersection if one
* 		exists.
* 		This is fairly inefficient but the overall cost should
* 		be insignificant compared to that of the threshold and
* 		label functions.
* \param	pObj		Sampling region object.
* \param	gVWSp		Grey value workspace.
* \param	off		Offset of the sample region.
* \param	dstMin		Destination pointer for the minimum value.
* \param	dstMax		Destination pointer for the maximum value.
* \param	dstSum		Destination pointer for the sum value.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static int 			WlzConComThreshComp2D(
				  WlzObject *pObj,
				  WlzGreyValueWSpace *gVWSp,
				  WlzIVertex3 off,
	      			  float *dstMin,
				  float *dstMax,
				  float *dstSum,
				  WlzErrorNum *dstErr)
{
  int		vol = 0;
  float		min = 0,
  		max = 0,
		sum = 0;
  WlzIVertex3 	pos;
  WlzIntervalWSpace iWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  pos.vtZ = off.vtZ;
  errNum = WlzInitRasterScan(pObj, &iWSp, WLZ_RASTERDIR_ILIC);
  while((errNum == WLZ_ERR_NONE) && 
        (errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE)
  {
    int		i;
    
    pos.vtY = iWSp.linpos + off.vtY;
    for(i = iWSp.lftpos; i <= iWSp.rgtpos; ++i)
    {
      pos.vtX = i + off.vtX;
      WlzGreyValueGet(gVWSp, pos.vtZ, pos.vtY, pos.vtX);
      if(gVWSp->bkdFlag == 0)
      {
	float	   fv;
        WlzPixelV  val;

        val.type = gVWSp->gType;
	val.v = gVWSp->gVal[0];
	(void )WlzValueConvertPixel(&val, val, WLZ_GREY_FLOAT);
	fv = val.v.flv;
	if(vol++ == 0)
	{
	  min = max = sum = fv;
	}
	else
	{
	  sum += fv;
	  if(min > fv)
	  {
	    min = fv;
	  }
	  else if(max < fv)
	  {
	    max = fv;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  *dstMin = min;
  *dstMax = max;
  *dstSum = sum;
  *dstErr = errNum;
  return(vol);
}

/*!
* \return	Volume (or area) of intersection.
* \ingroup	WlzThreshold
* \brief	Given a 3D sampling region object centred at the origin and
* 		a valid grey value workspace for an object. This function
* 		finds the min, max and sum value of the intersection if one
* 		exists.
* \param	pObj		Sampling region object.
* \param	gVWSp		Grey value workspace.
* \param	off		Offset of the sample region.
* \param	dstMin		Destination pointer for the minimum value.
* \param	dstMax		Destination pointer for the maximum value.
* \param	dstSum		Destination pointer for the sum value.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static int 			WlzConComThreshComp3D(
				  WlzObject *pObj,
				  WlzGreyValueWSpace *gVWSp,
				  WlzIVertex3 off,
	      			  float *dstMin,
				  float *dstMax,
				  float *dstSum,
				  WlzErrorNum *dstErr)
{
  int		p,
  		pln1,
  		lpln,
  		vol = 0;
  float		min = 0,
  		max = 0,
		sum = 0;
  WlzIVertex3 	pos;
  WlzDomain	*doms;
  WlzValues	nulVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  pos.vtX = off.vtX;
  pos.vtY = off.vtY;
  nulVal.core = NULL;
  doms = pObj->domain.p->domains;
  pln1 = pObj->domain.p->plane1;
  lpln = pObj->domain.p->lastpl;
  for(p = pln1; p < lpln; ++p)
  {
    int		pIdx;
    WlzObject	*obj2 = NULL;

    pIdx = pln1 - pObj->domain.p->plane1;
    if(doms[pIdx].core)
    {
      pos.vtZ = off.vtZ + p;
      obj2 = WlzAssignObject(
             WlzMakeMain(WLZ_2D_DOMAINOBJ, doms[pIdx], nulVal,
                         NULL, NULL, &errNum), NULL);
      if(obj2)
      {
	float	pMin,
		pMax,
		pSum,
		pVol;

	pVol = WlzConComThreshComp2D(obj2, gVWSp, pos, &pMin, &pMax, &pSum,
				     &errNum);
	if(pVol > 0)
	{
	  if(vol == 0)
	  {
	    min = pMin;
	    max = pMax;
	  }
	  else
	  {
	    if(pMin < min)
	    {
	      min = pMin;
	    }
	    if(pMax > max)
	    {
	      max = pMax;
	    }
	  }
	  vol += pVol;
	  sum += pSum;
	}
	(void )WlzFreeObj(obj2);
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  *dstMin = min;
  *dstMax = max;
  *dstSum = sum;
  *dstErr = errNum;
  return(vol);
}
