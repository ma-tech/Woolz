#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDistAllNearest_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzDistAllNearest.c
* \author       Bill Hill
* \date         April 2021
* \version      $Id$
* \par
* Copyright (C), [2021],
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
* \brief	Functions to compute nearest points in a reference domain at
* 		all points in a foreground domain.
* \ingroup	WlzMorphologyOps
*/

#include <stddef.h>
#include <limits.h>                             
#include <float.h>                             
#include <Wlz.h>                             

static WlzConnectType		WlzDANAlternateCon(
				  WlzConnectType con);
static WlzValues		WlzDANNewValues(
				  WlzObject *obj,
				  WlzObjectType gTT,
				  WlzPixelV bgdV,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzDANSetDstValues(
				  WlzGreyValueWSpace **gVWSp,
				  WlzObject *uknObj,
				  WlzObject *kwnObj,
				  WlzConnectType con);
static WlzErrorNum		WlzDANSetValues(
				  WlzObject *vObj,
				  WlzObject *mObj,
				  WlzPixelV mV);
static WlzErrorNum		WlzDANInitRefValues(
				  WlzGreyValueWSpace **gVWSp,
				  WlzObject *refObj);
static void 			WlzDANFreeValues(
				  int dim,
				  WlzValues val);
/*!
* \return	A Woolz compound array object or NULL on error.
* \ingroup	WlzMorphologyOps
* \brief	Computes the nearest points in a reference domain at all points
* 		in a foreground domain and (optionally) the distance of the
* 		reference domain points from their nearest points within the
* 		reference domain.
* 		The returned compound array object has integer reference point
* 		coordinates with column (x), line (y) and for 3D objects plane
* 		(z) having indices 0, 1 and 2 respectively.
* 		Distances are correct for the given distance function and are
* 		not normalised to remove over estimates with respect to
* 		constrained Euclidean distance.
* 		The intersection of the foreground and reference object
* 		is be used to ensure that the reference object is within
* 		the foreground object.
* \param        gForObj         Given foreground object.
* \param        gRefObj         Given reference object.
* \param	dFn		Distance function which must be appropriate to
* 				the dimension of the foreground and reference
* 				objects.
* \param	dstDstObj	Destination distance object pointer for the
* 				(float) distance values, may be NULL.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzCompoundArray		*WlzDistAllNearest(
				  WlzObject *gForObj,
				  WlzObject *gRefObj,
				  WlzDistanceType dFn,
				  WlzObject **dstDstObj,
				  WlzErrorNum *dstErr)
{
  int		dim = 0;
  WlzObject	*forObj = NULL,
  		*refObj = NULL,
		*dstObj = NULL;
  WlzConnectType con = WLZ_0_CONNECTED;
  WlzCompoundArray *nrpObj = NULL;
  WlzGreyValueWSpace *nrpGVWSp[4] = {0}; 
  WlzValues	nulVal = {0};
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  /* Check parameters. */
  if((gForObj == NULL) || (gRefObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gRefObj->type != gForObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else
  {
    switch(gForObj->type)
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
    if((gForObj->domain.core == NULL) || (gRefObj->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
  }
  /* Establish connectivity. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(dim == 2)
    {
      switch(dFn)
      {
	case WLZ_4_DISTANCE:
	  con = WLZ_4_CONNECTED;
	  break;
	case WLZ_8_DISTANCE:
	  con = WLZ_8_CONNECTED;
	  break;
	case WLZ_OCTAGONAL_DISTANCE:
	  con = WLZ_8_CONNECTED;
	  break; 
	default: 
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
      }
    }
    else
    {
      switch(dFn)
      {
	case WLZ_6_DISTANCE:
	  con = WLZ_6_CONNECTED;
	  break;
	case WLZ_18_DISTANCE:
	  con = WLZ_18_CONNECTED;
	  break;
	case WLZ_26_DISTANCE:
	  con = WLZ_26_CONNECTED;
	  break;
	case WLZ_OCTAGONAL_DISTANCE:
	  con = WLZ_26_CONNECTED;
	  break;
	default:
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
      }
    }
  }
  /* Create domain objects without values from the given foreground and
   * intersection of the foreground and reference objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    forObj = WlzMakeMain(gForObj->type, gForObj->domain, nulVal,
		         NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    refObj = WlzIntersect2(forObj, gRefObj, &errNum);
  }
  /* Create compound array for the integer coordinate values to be computed
   * and a float valued distance object all covering the foreground domain.
   * Within these objects set values inside the reference object to zero.
   * Initialise grey value workspaces for these objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzPixelV	bgdV;
    WlzObjectType gTT;

    bgdV.type = WLZ_GREY_INT;
    bgdV.v.inv = 0;
    gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT, NULL);
    nrpObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, dim, NULL,
                                  forObj->type, &errNum);
    for(int i = 0; (errNum == WLZ_ERR_NONE) && (i < dim); ++i)
    {
      WlzValues	val;

      val = WlzDANNewValues(forObj, gTT, bgdV, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	nrpObj->o[i] = WlzAssignObject(
                       WlzMakeMain(forObj->type, forObj->domain, val,
		                   NULL, NULL, &errNum), NULL);
        if(errNum != WLZ_ERR_NONE)
	{
	  WlzDANFreeValues(dim, val);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzDANSetValues(nrpObj->o[i], refObj, bgdV);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        nrpGVWSp[i] = WlzGreyValueMakeWSp(nrpObj->o[i], &errNum);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValues val = {0};

      bgdV.type = WLZ_GREY_FLOAT;
      bgdV.v.flv = 0.0;
      gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_FLOAT, NULL);
      val = WlzDANNewValues(forObj, gTT, bgdV, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	dstObj = WlzMakeMain(forObj->type, forObj->domain, val,
	                     NULL, NULL, &errNum);
        if(errNum != WLZ_ERR_NONE)
	{
	  WlzDANFreeValues(dim, val);
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzPixelV iV;

      /* Initialise all distances outside of the reference object to be
       * the maximum float and those within it to be zero. */
      iV.type = WLZ_GREY_FLOAT;
      iV.v.flv = FLT_MAX;
      errNum = WlzGreySetValue(dstObj, iV);
      if(errNum == WLZ_ERR_NONE)
      {
	iV.v.flv = 0.0f;
        errNum = WlzDANSetValues(dstObj, refObj, iV);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nrpGVWSp[3] = WlzGreyValueMakeWSp(dstObj, &errNum);
    }
  }
  /* Set cooordinates in (known) reference object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject       *tObj;

    tObj = WlzIntersect2(refObj, forObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzDANInitRefValues(nrpGVWSp, tObj);
    }
    (void )WlzFreeObj(tObj);
  }
  /* Propagate distances and nearest reference location out in shells.
   * Each unknown shell is given by:
   *
   *   U_i = K^{+}_i - K_i \cap R
   *   K_{i + 1} = K_i \cup U_{i}
   *
   * where U, K and R are the unknown, known and reference domains.
   * Operators ^{+}, -, \cap and \cup are dilation, difference in domain,
   * intersection and union. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*kwnObj = NULL,   /* Known distances and points. */
    		*uknObj = NULL;   /* Shell with distances and points to be
		                     computed. */

    kwnObj = WlzMakeMain(refObj->type, refObj->domain, nulVal,
	                 NULL, NULL, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      do
      {
        WlzObject	*tObj[2] = {0};

	if(dFn == WLZ_OCTAGONAL_DISTANCE)
	{
	  con = WlzDANAlternateCon(con);
	}
	WlzFreeObj(uknObj); uknObj = NULL;
	tObj[0] = WlzDilation(kwnObj, con, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  tObj[1] = WlzDiffDomain(tObj[0], kwnObj, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  uknObj = WlzIntersect2(tObj[1], forObj, &errNum);
	}
	(void )WlzFreeObj(tObj[0]);
	(void )WlzFreeObj(tObj[1]);
	if(errNum == WLZ_ERR_NONE)
	{
	  /* Check for all reference object covered, if so set uknObj to NULL
	   * to end the propagation. */
	  if(uknObj && WlzIsEmpty(uknObj, NULL))
	  {
	    (void )WlzFreeObj(uknObj);
	    uknObj = NULL;
	  }
	}
	if(uknObj)
	{
	  if(errNum == WLZ_ERR_NONE)
	  {
	    /* Set distance and coordinate object's values within the shell. */
	    errNum = WlzDANSetDstValues(nrpGVWSp, uknObj, kwnObj, con);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzObject *tObj;

	    tObj = WlzUnion2(kwnObj, uknObj, &errNum);
	    (void )WlzFreeObj(kwnObj);
	    kwnObj = tObj;
	  }
	}
      } while((errNum == WLZ_ERR_NONE) && (uknObj != NULL));
    }
    (void )WlzFreeObj(kwnObj);
    (void )WlzFreeObj(uknObj);
  }
  (void )WlzFreeObj(forObj);
  (void )WlzFreeObj(refObj);
  for(int i = 0; i < 4; ++i)
  {
    WlzGreyValueFreeWSp(nrpGVWSp[i]);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstDstObj)
    {
      *dstDstObj = dstObj;
      dstObj = NULL;
    }
  }
  else
  {
    (void )WlzFreeObj((WlzObject *)nrpObj);
    nrpObj = NULL;
  }
  (void )WlzFreeObj(dstObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nrpObj);
}

/*!
* \return	Wool error code.
* \ingroup	WlzMorphologyOps
* \brief	Sets initial coordinate values in the known object.
* \param	gVWSp			Array of grey value workspaces for
* 					coordinate values.
* \param	refObj			Object with reference domain.
*/
static WlzErrorNum		WlzDANInitRefValues(
				  WlzGreyValueWSpace **gVWSp,
				  WlzObject *refObj)
{
  int		dim;
  size_t	off[3];
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  off[0] = offsetof(WlzIVertex3, vtX); 
  off[1] = offsetof(WlzIVertex3, vtY); 
  off[2] = offsetof(WlzIVertex3, vtZ); 
  dim = (refObj->type == WLZ_2D_DOMAINOBJ)? 2: 3;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < dim; ++i)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      WlzErrorNum	errNum2 = WLZ_ERR_NONE;
      WlzIterateWSpace *iWSp;

      iWSp = WlzIterateInit(refObj, WLZ_RASTERDIR_IPILIC, 0, &errNum2);
      while(errNum2 == WLZ_ERR_NONE)
      {
	if((errNum2 = WlzIterate(iWSp)) == WLZ_ERR_NONE)
	{
	  char 	*base;

	  base = (char *)&(iWSp->pos);
	  WlzGreyValueGet(gVWSp[i],
	                  iWSp->pos.vtZ, iWSp->pos.vtY, iWSp->pos.vtX);
	  *(gVWSp[i]->gPtr[0].inp) = *(int *)(base + off[i]);
	}
      }
      if(errNum2 == WLZ_ERR_EOO)
      {
	errNum2 = WLZ_ERR_NONE;
      }
      WlzIterateWSpFree(iWSp);
#ifdef _OPENMP
#pragma omp critical (WlzDANInitRefValues)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = errNum2;
	}
      }
#else
      errNum = errNum2;
#endif
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMorphologyOps
* \brief	Given a valid grey value workspaces for integer cooordinates
* 		and float distances. For use by WlzDistAllNearest().
* \param	gVWSp			Array of grey value workspaces with 
* 					components 0 = x coordinate, 1 = y
* 					coordinate, (for 3D) 2 = z coordinate
* 					and 3 = distance.
* \param	uknObj			Current shell object with unknown values
* 					of coordinate or distance. Distance
* 					values are intialised to FLT_MAX.
* \param	kwnObj			Known object with all known coordinates
* 					and distances.
* \param	con			Connectivity to use.
*/
static WlzErrorNum		WlzDANSetDstValues(
				  WlzGreyValueWSpace **gVWSp,
				  WlzObject *uknObj,
				  WlzObject *kwnObj,
				  WlzConnectType con)
{
  int		dim,
  		nNbr;
  const double  *nbrDst;
  const WlzIVertex3 *nbrOff;
  WlzIterateWSpace *iWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double  nbrDst2[8] =    {1.0, 1.0, 1.0, 1.0,
				 ALG_M_SQRT2, ALG_M_SQRT2,
				 ALG_M_SQRT2, ALG_M_SQRT2};
  const double  nbrDst3[26] =     {1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				   ALG_M_SQRT2, ALG_M_SQRT2,
				   ALG_M_SQRT2, ALG_M_SQRT2,
				   ALG_M_SQRT2, ALG_M_SQRT2,
				   ALG_M_SQRT2, ALG_M_SQRT2,
				   ALG_M_SQRT2, ALG_M_SQRT2,
				   ALG_M_SQRT2, ALG_M_SQRT2,
				   ALG_M_SQRT3, ALG_M_SQRT3,
				   ALG_M_SQRT3, ALG_M_SQRT3,
				   ALG_M_SQRT3, ALG_M_SQRT3,
				   ALG_M_SQRT3, ALG_M_SQRT3};
  const WlzIVertex3 nbrOff2[8] = {{-1,  0,  0},  /* 0 */
				   { 1,  0,  0},  /* 1 */
				   { 0, -1,  0},  /* 2 */
				   { 0,  1,  0},  /* 3 */
				   {-1, -1,  0},  /* 4 */
				   { 1, -1,  0},  /* 5 */
				   {-1,  1,  0},  /* 6 */
				   { 1,  1,  0}}; /* 7 */
  const WlzIVertex3 nbrOff3[26] = {{-1,  0,  0},  /*  0 */
				   { 1,  0,  0},  /*  1 */
				   { 0, -1,  0},  /*  2 */
				   { 0,  1,  0},  /*  3 */
				   { 0,  0, -1},  /*  4 */
				   { 0,  0,  1},  /*  5 */
				   {-1,  0, -1},  /*  6 */
				   { 1,  0, -1},  /*  7 */
				   { 0, -1, -1},  /*  8 */
				   { 0,  1, -1},  /*  9 */
				   {-1,  0,  1},  /* 10 */
				   { 1,  0,  1},  /* 11 */
				   { 0, -1,  1},  /* 12 */
				   { 0,  1,  1},  /* 13 */
				   {-1, -1,  0},  /* 14 */
				   { 1, -1,  0},  /* 15 */
				   {-1,  1,  0},  /* 16 */
				   { 1,  1,  0},  /* 17 */
				   {-1, -1, -1},  /* 18 */
				   { 1, -1, -1},  /* 19 */
				   {-1,  1, -1},  /* 20 */
				   { 1,  1, -1},  /* 21 */
				   {-1, -1,  1},  /* 22 */
				   { 1, -1,  1},  /* 23 */
				   {-1,  1,  1},  /* 24 */
				   { 1,  1,  1}}; /* 25 */

  switch(con)
  {
    case WLZ_4_CONNECTED:
      nNbr = 4;
      nbrDst = nbrDst2;
      nbrOff = nbrOff2;
      break;
    case WLZ_8_CONNECTED:
      nNbr = 8;
      nbrDst = nbrDst2;
      nbrOff = nbrOff2;
      break;
    case WLZ_6_CONNECTED:
      nNbr = 6;
      nbrDst = nbrDst3;
      nbrOff = nbrOff3;
      break;
    case WLZ_18_CONNECTED:
      nNbr = 18;
      nbrDst = nbrDst3;
      nbrOff = nbrOff3;
      break;
    case WLZ_26_CONNECTED:
      nNbr = 26;
      nbrDst = nbrDst3;
      nbrOff = nbrOff3;
      break;
    default:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
  }
  dim = (uknObj->type == WLZ_2D_DOMAINOBJ)? 2: 3;
  iWSp = WlzIterateInit(uknObj, WLZ_RASTERDIR_IPILIC, 0, &errNum);
  while(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzIterate(iWSp);
    if(errNum == WLZ_ERR_NONE)
    {
      int	nIn = 0,
      		nMin = -1;
      double	dMin,
      		dSum = 0.0;
      float	*dUknP;
      WlzIVertex3 pos;

      WlzGreyValueGet(gVWSp[3],
                      iWSp->pos.vtZ, iWSp->pos.vtY, iWSp->pos.vtX);
      dUknP = gVWSp[3]->gPtr[0].flp;
      dMin = *dUknP;
      for(int i = 0; i < nNbr; ++i)
      {

	WLZ_VTX_3_ADD(pos, iWSp->pos, nbrOff[i]);
	if(WlzInsideDomain(kwnObj, pos.vtZ, pos.vtY, pos.vtX, NULL))
	{
	  double    d;

	  WlzGreyValueGet(gVWSp[3], pos.vtZ, pos.vtY, pos.vtX);
	  d = gVWSp[3]->gVal[0].flv + nbrDst[i];
	  ++nIn;
	  dSum += d;
	  if(d < dMin)
	  {
	    dMin = d;
	    nMin = i;
	  }
	}
      }
      if(nMin >= 0)
      {
        *dUknP = dSum / nIn;
	WLZ_VTX_3_ADD(pos, iWSp->pos, nbrOff[nMin]);
	for(int i = 0; i < dim; ++i)
	{
	  int	v;

	  WlzGreyValueGet(gVWSp[i], pos.vtZ, pos.vtY, pos.vtX);
	  v = gVWSp[i]->gVal[0].inv;
	  WlzGreyValueGet(gVWSp[i],
			  iWSp->pos.vtZ, iWSp->pos.vtY, iWSp->pos.vtX);
	  *(gVWSp[i]->gPtr[0].inp) = v;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  WlzIterateWSpFree(iWSp);
  return(errNum);
}

/*!
* \return	Next connectivity.
* \ingroup	WlzMorphologyOps
* \brief	Alternates connectivity for when using octagonal distances.
* \param	con			Current connectivity.
*/
static WlzConnectType		WlzDANAlternateCon(
				  WlzConnectType con)
{
  switch(con)
  {
    case WLZ_4_CONNECTED:
      con = WLZ_8_CONNECTED;
      break;
    case WLZ_6_CONNECTED:
      con = WLZ_26_CONNECTED;
      break;
    case WLZ_8_CONNECTED:
      con = WLZ_4_CONNECTED;
      break;
    case WLZ_26_CONNECTED:
      con = WLZ_6_CONNECTED;
      break;
    default:
      break;
  }
  return(con);
}

/*!
* \return	Woolz values union.
* \ingroup	WlzMorphologyOps
* \brief	Creates a new value table, simply a wrapper for WlzNewValueTb()
* 		and WlzNewValuesVox().
* \param	obj			Given object with valid type and domain.
* \param	gTT			Required grey table type.
* \param	bgdV			Background value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzValues		WlzDANNewValues(
				  WlzObject *obj,
				  WlzObjectType gTT,
				  WlzPixelV bgdV,
				  WlzErrorNum *dstErr)
{
  WlzValues	val = {0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        val.v = WlzNewValueTb(obj, gTT, bgdV, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
        val.vox = WlzNewValuesVox(obj, gTT, bgdV, &errNum);
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
  return(val);
}

/*!
* \ingroup	WlzMorphologyOps
* \brief	Frees the given value table, simply a wrapper for
* 		WlzFreeValueTb() and WlzFreeVoxelValueTb().
* \param	dim			Dimension, 1 for pixel value table and
* 					2 for voxel value table.
* \param	val			Values union.
*/
static void 			WlzDANFreeValues(
				  int dim,
				  WlzValues val)
{
  if(dim == 2)
  {
    (void )WlzFreeValueTb(val.v);
  }
  else
  {
    (void )WlzFreeVoxelValueTb(val.vox);
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMorphologyOps
* \brief	Sets values in one object to a constant value using the domain
* 		of a second object.
* \param	vObj			Object with values to be set in place.
* \param	mObj			Mask object with domain.
* \param	mV			Mask value to be set.
*/
static WlzErrorNum		WlzDANSetValues(
				  WlzObject *vObj,
				  WlzObject *mObj,
				  WlzPixelV mV)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(vObj->type != mObj->type)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((vObj->domain.core == NULL) || (mObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(vObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(mObj->type)    
    {
      case WLZ_2D_DOMAINOBJ:
        {
	  WlzObject   *tObj;

	  tObj = WlzMakeMain(vObj->type, mObj->domain, vObj->values,
	                     NULL, NULL, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzGreySetValue(tObj, mV);
	  }
	  (void )WlzFreeObj(tObj);
	}
	break;
      case WLZ_3D_DOMAINOBJ:
        {
	  WlzValues   tVal;
	  WlzObject   *tObj = NULL;

	  tVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
	                                 mObj->domain.p->plane1,
					 mObj->domain.p->lastpl,
					 mV, NULL, &errNum);
	  if(errNum == WLZ_ERR_NONE)
          {
	    tObj = WlzMakeMain(mObj->type, mObj->domain, tVal,
	    		       NULL, NULL, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
          {
	    int		mp1,
	    		mlp,
			vp1;

            mp1 = mObj->domain.p->plane1;
	    mlp = mObj->domain.p->lastpl;
	    vp1 = vObj->domain.p->plane1;
	    for(int p = mp1; p <= mlp; ++p)
	    {
	      WlzDomain	*pD;
	      WlzValues *pV;

	      pD = mObj->domain.p->domains + p - mp1;
	      if((*pD).core)
	      {
		pV = vObj->values.vox->values + p - vp1;
	      }
	      *(tObj->values.vox->values + p - mp1) =
	          WlzAssignValues(*pV, NULL);
	    }
	    errNum = WlzGreySetValue(tObj, mV);
	  }
	  if(tObj)
	  {
	    (void )WlzFreeObj(tObj);
	  }
	  else if(tVal.core)
	  {
	    (void )WlzFreeVoxelValueTb(tVal.vox);
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

