#pragma ident "MRC HGU $Id$"
/*!
* \file        	WlzRegICP.c
* \author       Bill Hill
* \date         November 2000
* \version      $Id$
* \note
*               Copyright:
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par  Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions for the registration of two objects using an
*               implementation of the iterative closest point
*               algorithm.
*               See
*                 * Besel P. J. and McKay N. D. A method for
*                   registration of 3-D shapes. PAMI, 14(2): 239-256,
*                   1992.
*                 * Zhang Z. Iterative point matching for
*                   registration of free-form curves and surfaces.
*                   International Journal of Computer Vision,
*                   13(2):119-152, 1994.
* \ingroup      WlzTransform
* \todo         -
* \bug
*/

#include <float.h>
#include <Wlz.h>

/*!
* \struct	_WlzRegICPWSp
* \ingroup      WlzTransform
* \brief	A workspace data structure for use in the ICP functions.
*/
struct _WlzRegICPWSp
{
  /* Itteration control. */
  int		itr;		/*!< Current iteration */
  double	sumDistCur;	/*!< Current sum of distances between NN */
  double	sumDistLst;	/*!< Last sum of distances between NN */
  /* Nearest neighbour search. */
  double	maxDist;	/*!< Maximum distance to consider for a NN */
  AlcKDTTree	*tTree;		/*!< kD-tree */
  int		*sNN;		/*!< Indicies of NN to source verticies */
  double	*dist;		/*!< NN distances */
  /* Vericies and normals. */
  WlzVertexType vType;		/*!< Type of verticies WLZ_VERTEX_D2 or
  				     WLZ_VERTEX_D3 */
  int		nT;		/*!< Number of target verticies/normals */
  int		nS;		/*!< Number of source verticies/normals */
  int		nMatch;		/*!< Minimum of nT and nS */
  WlzVertexP	gTVx;		/*!< Given target verticies */
  WlzVertexP	gSVx;		/*!< Given source verticies */
  WlzVertexP	gTNr;		/*!< Given target normals */
  WlzVertexP	gSNr;		/*!< Given source normals */
  WlzVertexP    tSVx;		/*!< Transformed source verticies. */
  WlzVertexP    tSNr;		/*!< Transformed source normals. */
  WlzVertexP    nNTVx;		/*!< NN ordered target verticies. */
  /* Match weighting */
  double	kVx;		/*!< Coefficient for vertex weighting */
  double	kNr;		/*!< Coefficient for normal weighting */
  double	*wgtVx;		/*!< Weights for vertex matches */
  double	*wgtNr;		/*!< Weights for normal matches */
  /* Affine transform. */
  WlzAffineTransform *tr;	/*!< Current affine transform */
};
typedef struct _WlzRegICPWSp WlzRegICPWSp;

static void     		WlzRegICPTrans(
				  WlzRegICPWSp *wSp);
static void			WlzRegICPFindNN(
				  WlzRegICPWSp *wSp);
static void			WlzRegICPWeight(
				  WlzRegICPWSp *wSp);
static int			WlzRegICPItr(
				  WlzRegICPWSp *wSp,
			          WlzTransformType trType,
				  int maxItr,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzRegICPCompTransform(
				  WlzRegICPWSp *wSp,
				  WlzTransformType trType,
				  int *dstConv);
static WlzErrorNum 		WlzRegICPCheckVerticies(
				  WlzVertexP *vData,
				  int *vCnt,
				  WlzVertexType *vType);
static WlzErrorNum 		WlzRegICPBuildTree(
				  WlzRegICPWSp *wSp);

/*!
* \ingroup	WlzTransform
* \return				Affine transform which brings
*					the two objects into register.
* \brief	Registers the two given objects using the iterative
* 		closest point algorithm. Maximal gradient contours
*		and gradients are extracted from the two given domain
*		objects with values. An affine transform is computed,
*		which when applied to the source object takes it
*		into register with the target object.
* \param	tObj			The target object.
* \param	sObj			The source object to be
*					registered with target object.
* \param	initTr			 Initial affine transform
*					to be applied to the source
*					object prior to using the ICP
*					algorithm. May be NULL.
* \param	trType			Required transform type.
* \param	ctrLo			Low threshold contour gradient
*					value.
* \param	ctrHi			High threshold contour gradient
*					value.
* \param	ctrWth			Contour filter width.
* \param	dstConv			Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
* \param	dstItr			Destination ptr for the number
*					of iterations, may be NULL.
* \param	maxItr			Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzAffineTransform *WlzRegICPObjsGrd(WlzObject *tObj, WlzObject *sObj,
				     WlzAffineTransform *initTr,
				     WlzTransformType trType,
				     double ctrLo, double ctrHi, double ctrWth,
				     int *dstConv, int *dstItr, int maxItr,
				     WlzErrorNum *dstErr)
{
  int		idN,
  		idM,
		conv,
		itr;
  double	tD0;
  int		vCnt[2];
  WlzVertexType	vType;
  WlzDomain	rDom;
  WlzValues	rVal;
  WlzObject	*cObj;
  WlzObject	*gObj[2];
  WlzVertexP	nData[2],
  		vData[2];
  WlzDVertex2	tDV2;
  WlzDVertex3	tDV3;
  WlzGMVertex	*mVtx;
  WlzGMModel	*model;
  WlzAffineTransform *regTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
#ifdef WLZ_REGICP_DEBUG
  FILE		*fP = NULL;
  char		dbgStr[256];
#endif /* WLZ_REGICP_DEBUG */

  gObj[0] = tObj;
  gObj[1] = sObj;
  rDom.core = NULL;
  rVal.core = NULL;
  if((tObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(tObj->type != sObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  /* Create maximal gradient contours and extract the maximal gradient
   * verticies and their gradient normals. */
  /* For the target object. */
  for(idN = 0; idN < 2; ++idN)
  {
    cObj = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      rDom.ctr = WlzContourObjGrd(gObj[idN], NULL, ctrLo, ctrHi, ctrWth,
				  &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      cObj = WlzMakeMain(WLZ_CONTOUR, rDom, rVal, NULL, NULL, &errNum);
    }
#ifdef WLZ_REGICP_DEBUG
    if(errNum == WLZ_ERR_NONE)
    {
      (void )sprintf(dbgStr, "dbg-cObj%d.wlz", idN);
      if((fP = fopen(dbgStr, "w")) != NULL)
      {
	(void )WlzWriteObj(fP, cObj);
        (void )fclose(fP);
	fP = NULL;
      }
    }
#endif /* WLZ_REGICP_DEBUG */
    /* Get an array of position verticies and normals from the model,
     * don't use the gradient normals because they are eratic. */
    if(errNum == WLZ_ERR_NONE)
    {
      vData[idN] = WlzVerticiesFromObj(cObj, nData + idN, vCnt + idN, &vType,
				       &errNum);
    }
    (void )WlzFreeObj(cObj);
#ifdef WLZ_REGICP_DEBUG
    if(errNum == WLZ_ERR_NONE)
    {
      (void )sprintf(dbgStr, "dbg-vn%d.num", idN);
      if((fP = fopen(dbgStr, "w")) != NULL)
      {
	switch(vType)
	{
	  case WLZ_VERTEX_D2:
	    for(idM = 0; idM < vCnt[idN]; ++idM)
	    {
	      (void )fprintf(fP, "%d %g %g %g %g\n",
	               idM,
		       (vData[idN].d2 + idM)->vtX, (vData[idN].d2 + idM)->vtY,
		       (nData[idN].d2 + idM)->vtX, (nData[idN].d2 + idM)->vtY);
	    }
	    break;
	  case WLZ_VERTEX_D3:
	    for(idM = 0; idM < vCnt[idN]; ++idM)
	    {
	      (void )fprintf(fP, "%d %g %g %g %g  %g %g\n",
	               idM,
		       (vData[idN].d3 + idM)->vtX, (vData[idN].d3 + idM)->vtY,
		       (vData[idN].d3 + idM)->vtZ,
		       (nData[idN].d3 + idM)->vtX, (nData[idN].d3 + idM)->vtY,
		       (nData[idN].d3 + idM)->vtZ);
	    }
	    break;
	}
        (void )fclose(fP);
	fP = NULL;
      }
    }
#endif /* WLZ_REGICP_DEBUG */
  }
  /* Register the position verticies and normals. */
  if(errNum == WLZ_ERR_NONE)
  {
    regTr = WlzRegICPVerticies(vData[0], nData[0], vCnt[0],
			       vData[1], nData[1], vCnt[1],
			       vType, initTr,
			       trType, &conv, &itr, maxItr, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstConv)
    {
      *dstConv = conv;
    }
    if(dstItr)
    {
      *dstItr = itr;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/*!
* \ingroup	WlzTransform
* \return				Affine transform which brings
*					the two objects into register.
* \brief	Registers the two given objects using the iterative
* 		closest point algorithm. An affine transform is
*		computed, which when applied to the source object
*		takes it into register with the target object.
* \param	tObj			The target object.
* \param	sObj			The source object to be
*					registered with target object.
* \param	initTr			Initial affine transform
*					to be applied to the source
*					object prior to using the ICP
*					algorithm. May be NULL.
* \param	trType			Required transform type.
* \param	dstConv			Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
* \param	dstItr			Destination ptr for the number
*					of iterations, may be NULL.
* \param	maxItr			Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzAffineTransform *WlzRegICPObjs(WlzObject *tObj, WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv, int *dstItr, int maxItr,
				  WlzErrorNum *dstErr)
{
  int		idx,
		conv,
		itr;
  int		vCnt[2];
  WlzVertexType	vType[2];
  WlzObject	*objs[2];
  WlzVertexP	nData[2],
  		vData[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzAffineTransform *regTr = NULL;

  nData[0].v = nData[1].v = NULL;
  vData[0].v = vData[1].v = NULL;
  objs[0] = tObj; objs[1] = sObj;
  if((tObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      vData[idx] = WlzVerticiesFromObj(objs[idx], nData + idx, vCnt + idx,
      				       vType + idx, &errNum);
      ++idx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check type of verticies and convert to double if required. */
    errNum = WlzRegICPCheckVerticies(vData, vCnt, vType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
     regTr = WlzRegICPVerticies(vData[0], nData[0], vCnt[0],
     				vData[1], nData[1], vCnt[1],
     				vType[0], initTr,
				trType, &conv, &itr, maxItr, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstConv)
    {
      *dstConv = conv;
    }
    if(dstItr)
    {
      *dstItr = itr;
    }
  }
  for(idx = 0; idx < 2; ++idx)
  {
    if(vData[idx].v)
    {
      AlcFree(vData[idx].v);
    }
    if(nData[idx].v)
    {
      AlcFree(nData[idx].v);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/*!
* \ingroup	WlzTransform
* \return				Woolz error code.
* \brief	Checks the two sets of verticies and promotes the
*		vertex type to double (2 or 3D) if not already double
*		verticies.
* \param	vData			The sets of verticies.
* \param	vCnt			The numbers of verticies.
* \param	vType			The types of verticies.
*/
static WlzErrorNum WlzRegICPCheckVerticies(WlzVertexP *vData, int *vCnt,
					   WlzVertexType *vType)
{
  int		idx,
  		copiedFlg;
  int		dim[2],
  		type[2];
  WlzVertexP	tVP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find vertex set dimensions. */
  idx = 0;
  do
  {
    switch(vType[idx])
    {
      case WLZ_VERTEX_I2: /* FALLTROUGH */
      case WLZ_VERTEX_F2: /* FALLTROUGH */
      case WLZ_VERTEX_D2:
	dim[idx] = 2;
	break;
      case WLZ_VERTEX_I3: /* FALLTROUGH */
      case WLZ_VERTEX_F3: /* FALLTROUGH */
      case WLZ_VERTEX_D3:
	dim[idx] = 3;
	break;
      default:
	errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  } while((++idx < 2) && (errNum == WLZ_ERR_NONE));
  if((errNum == WLZ_ERR_NONE) && (dim[0] != dim[1]))
  {
    /* This is unlikly to work, should it be allowed?  */
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find vertex set primitive types. */
    idx = 0;
    do
    {
      switch(vType[idx])
      {
	case WLZ_VERTEX_I2: /* FALLTROUGH */
	case WLZ_VERTEX_I3:
	  type[idx] = 0;
	case WLZ_VERTEX_F2: /* FALLTROUGH */
	case WLZ_VERTEX_F3:
	  type[idx] = 1;
	case WLZ_VERTEX_D2:
	case WLZ_VERTEX_D3:
	  type[idx] = 2;
	  break;
      }
    } while(++idx < 2);
    idx = 0;
    do
    {
      if((type[0] == 0) || (type[0] == 1) ||
         (type[1] == 0) || (type[1] == 1))
      {
        if(dim[idx] == 2)
	{
	  if((tVP.d2 = (WlzDVertex2 *)
	                AlcMalloc(vCnt[idx] * sizeof(WlzDVertex2))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    copiedFlg = 0;
	    switch(vType[idx])
	    {
	      case WLZ_VERTEX_I2:
		WlzValueCopyIVertexToDVertex(tVP.d2, vData[idx].i2,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	      case WLZ_VERTEX_F2:
		WlzValueCopyFVertexToDVertex(tVP.d2, vData[idx].f2,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	    }
	    if(copiedFlg)
	    {
	      vType[idx] = WLZ_VERTEX_D2;
	      AlcFree(vData[idx].v);
	      vData[idx].d2 = tVP.d2;
	    }
	  }
	}
	else /* dim[idx] == 3 */
	{
	  if((tVP.d3 = (WlzDVertex3 *)
	                AlcMalloc(vCnt[idx] * sizeof(WlzDVertex3))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    copiedFlg = 0;
	    switch(vType[idx])
	    {
	      case WLZ_VERTEX_I3:
		WlzValueCopyIVertexToDVertex3(tVP.d3, vData[idx].i3,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	      case WLZ_VERTEX_F3:
		WlzValueCopyFVertexToDVertex3(tVP.d3, vData[idx].f3,
					      vCnt[idx]);
		copiedFlg = 1;
		break;
	    }
	    if(copiedFlg)
	    {
	      vType[idx] = WLZ_VERTEX_D3;
	      AlcFree(vData[idx].v);
	      vData[idx].d3 = tVP.d3;
	    }
	  }
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (++idx < 2));
  }
  return(errNum);
}

/*!
* \ingroup	WlzTransform
* \return				Affine transform which brings
*					the two sets of verticies into
*					register.
* \brief	Registers the two given sets of verticies using the
*		iterative closest point algorithm. An affine transform
*		is computed, which when applied to the source verticies
*		takes it into register with the target verticies.
*		The verticies and their normals are known to be either
*		WlzDVertex2 or WlzDVertex3.
* \param	tVx			Target verticies.
* \param	tNr			Target normals, may be NULL.
* \param	tCnt			Number of target verticies.
* \param	sVx			Source verticies.
* \param	sNr			Source normals, may be NULL.
* \param	sCnt			Number of source verticies.
* \param	vType			Type of the verticies.
* \param	initTr			Initial affine transform
*					to be applied to the source
*					object prior to using the ICP
*					algorithm. May be NULL.
* \param	trType			Required transform type.
* \param	dstConv			Destination ptr for the
*					convergence flag (non zero
*					on convergence), may be NULL.
* \param	dstItr			Destination ptr for the number
*					of iterations, may be NULL.
* \param	maxItr			Maximum number of iterations,
*					if <= 0 then infinite iterations
*					are allowed.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzAffineTransform	*WlzRegICPVerticies(WlzVertexP tVx, WlzVertexP tNr,
					    int tCnt,
					    WlzVertexP sVx, WlzVertexP sNr,
					    int sCnt,
					    WlzVertexType vType,
					    WlzAffineTransform *initTr,
					    WlzTransformType trType,
					    int *dstConv, int *dstItr,
					    int maxItr,
					    WlzErrorNum *dstErr)
{

  int		conv,
		maxCnt;
  WlzAffineTransform *regTr = NULL;
  WlzRegICPWSp	wSp;
  AlcKDTTree	*tTree = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  /* Setup workspace. */
  wSp.itr = 0;
  wSp.sumDistCur = DBL_MAX;
  wSp.sumDistLst = DBL_MAX;
  wSp.maxDist = DBL_MAX;
  wSp.tTree = NULL;
  wSp.sNN = NULL;
  wSp.dist = NULL;
  wSp.vType = vType;
  wSp.nT = tCnt;
  wSp.nS = sCnt;
  wSp.nMatch = WLZ_MIN(tCnt, sCnt);
  wSp.gTVx = tVx;
  wSp.gSVx = sVx;
  wSp.gTNr = tNr;
  wSp.gSNr = sNr;
  wSp.tSVx.v = NULL;
  wSp.tSNr.v = NULL;
  wSp.nNTVx.v = NULL;
  wSp.kVx = 0.5;
  wSp.kNr = 0.5;
  wSp.wgtVx = NULL;
  wSp.wgtNr = NULL;
  wSp.tr = NULL;
  maxCnt = WLZ_MAX(tCnt, sCnt);
  if(((wSp.sNN = (int *)AlcMalloc(sizeof(int) * maxCnt)) == NULL) ||
     ((wSp.dist = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL) ||
     ((wSp.wgtVx = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL) ||
     ((wSp.wgtNr = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(initTr)
    {
      wSp.tr = WlzAffineTransformCopy(initTr, &errNum);
    }
    else
    {
      if(vType == WLZ_VERTEX_D2)
      {
	wSp.tr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum);
      }
      else /* vType == WLZ_VERTEX_D3 */
      {
        wSp.tr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, &errNum);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(vType == WLZ_VERTEX_D2)
    {
      if(((wSp.tSVx.d2 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
      						maxCnt)) == NULL) ||
         ((wSp.nNTVx.d2 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
      						maxCnt)) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      if((errNum == WLZ_ERR_NONE) && (tNr.v != NULL))
      {
	if((wSp.tSNr.d2 = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
						maxCnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
    else /* vType == WLZ_VERTEX_D3 */
    {
      if(((wSp.tSVx.d3 = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) *
      						maxCnt)) == NULL) ||
         ((wSp.nNTVx.d3 = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) *
      						maxCnt)) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      if((errNum == WLZ_ERR_NONE) && (tNr.v != NULL))
      {
	if((wSp.tSNr.d3 = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) *
						maxCnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Build k-D tree for the target verticies. */
    errNum = WlzRegICPBuildTree(&wSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(trType)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
      case WLZ_TRANSFORM_3D_AFFINE:
        conv = WlzRegICPItr(&wSp, (trType == WLZ_TRANSFORM_2D_AFFINE)?
			          WLZ_TRANSFORM_2D_REG: WLZ_TRANSFORM_3D_REG,
			    maxItr, &errNum);
	if(conv && (errNum == WLZ_ERR_NONE))
	{
	  wSp.sumDistCur = DBL_MAX;
	  wSp.sumDistLst = DBL_MAX;
          conv = WlzRegICPItr(&wSp, trType, maxItr, &errNum);
	}
        break;
      case WLZ_TRANSFORM_2D_REG:
      case WLZ_TRANSFORM_3D_REG:
        conv = WlzRegICPItr(&wSp, trType, maxItr, &errNum);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstConv)
    {
      *dstConv = conv;
    }
    if(dstItr)
    {
      *dstItr = wSp.itr;
    }
  }
  AlcFree(wSp.sNN);
  AlcFree(wSp.dist);
  AlcFree(wSp.wgtVx);
  AlcFree(wSp.wgtNr);
  AlcFree(wSp.tSVx.v);
  AlcFree(wSp.tSNr.v);
  AlcFree(wSp.nNTVx.v);
  if(wSp.tr)
  {
    if((errNum == WLZ_ERR_NONE) && conv)
    {
      regTr = wSp.tr;
    }
    else
    {
      (void )WlzFreeAffineTransform(wSp.tr);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/*!
* \ingroup	WlzTransform
* \return				Woolz error code
* \brief	Allocates and populates a k-D tree from the given verticies.
* 		The verticies are either WlzDVertex2 orWlzDVertex3
* \param	wSp			ICP registration workspace.
*/
static WlzErrorNum WlzRegICPBuildTree(WlzRegICPWSp *wSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Use wSp->sNN as a temporay buffer of shuffling indcies. */
  wSp->tTree = WlzVerticiesBuildTree(wSp->vType, wSp->nT, wSp->gTVx,
  				     wSp->sNN, &errNum);
  return(errNum);
}

/*!
* \ingroup	WlzTransform
* \return				Nonzero if the iteration has
* 					converged.
* \brief	The iterative loop of the ICP, which iterates to find the
* 		registration transform.
* \param	wSp			ICP registration workspace.
* \param	trType			Type of affine transform.
* \param	maxItr			Maximum number of iterations.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static int	WlzRegICPItr(WlzRegICPWSp *wSp,
			     WlzTransformType trType, int maxItr,
			     WlzErrorNum *dstErr)
{
  int		conv = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  do
  {
    wSp->sumDistLst = wSp->sumDistCur;
    /* Apply current transform to the verticies and normals. */
    WlzRegICPTrans(wSp);
    /* Find closest points. */
    WlzRegICPFindNN(wSp);
    /* Compute weightings of the matches. */
    WlzRegICPWeight(wSp);
    /* Compute registration transform and check for convergence. */
    errNum = WlzRegICPCompTransform(wSp, trType, &conv);
  } while((errNum == WLZ_ERR_NONE) && (wSp->itr++ < maxItr) && (conv == 0));
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(conv);
}

/*!
* \ingroup	WlzTransform
* \return       <void>
* \brief        Transforms the source verticies and normals using the
*               current affine transform.
* \param        wSp                     ICP registration workspace.
*/
static void     WlzRegICPTrans(WlzRegICPWSp *wSp)
{
  int		idx;

  if(wSp->vType == WLZ_VERTEX_D2)
  {
    for(idx = 0; idx < wSp->nS; ++idx)
    {
      *(wSp->tSVx.d2 + idx) = WlzAffineTransformVertexD2(wSp->tr,
      						*(wSp->gSVx.d2 + idx), NULL);
    }
    if(wSp->gSNr.v)
    {
      for(idx = 0; idx < wSp->nS; ++idx)
      {
        *(wSp->tSNr.d2 + idx) = WlzAffineTransformNormalD2(wSp->tr,
						*(wSp->gSNr.d2 + idx), NULL);
      }
    }
  }
  else /* wSp->vType == WLZ_VERTEX_D3 */
  {
    for(idx = 0; idx < wSp->nS; ++idx)
    {
      *(wSp->tSVx.d3 + idx) = WlzAffineTransformVertexD3(wSp->tr,
      						*(wSp->gSVx.d3 + idx), NULL);
    }
    if(wSp->gSNr.v)
    {
      for(idx = 0; idx < wSp->nS; ++idx)
      {
        *(wSp->tSNr.d3 + idx) = WlzAffineTransformNormalD3(wSp->tr,
						*(wSp->gSNr.d3 + idx), NULL);
      }
    }
  }
}

/*!
* \ingroup	WlzTransform
* \return	<void>
* \brief	Finds nearest neighbour matches in the target tree for
*		the source verticies, sets the nearest neighbour
*		indicies and permutes the NN ordered target verticies in
*		the workspace.
* \param	wSp			ICP registration workspace.
*/
static void	WlzRegICPFindNN(WlzRegICPWSp *wSp)
{
  int		idx;
  WlzDVertex2	tVD2;
  WlzDVertex3	tVD3;
  double	datD[3];
  AlcKDTNode	*node;

  wSp->sumDistCur = 0.0;
  for(idx = 0; idx < wSp->nMatch; ++idx)
  {
    if(wSp->vType == WLZ_VERTEX_D2)
    {
      tVD2 = *(wSp->tSVx.d2 + idx);
      datD[0] = tVD2.vtX;
      datD[1] = tVD2.vtY;
    }
    else /* wSp->vType == WLZ_VERTEX_D3 */
    {
      tVD3 = *(wSp->tSVx.d3 + idx);
      datD[0] = tVD3.vtX;
      datD[1] = tVD3.vtY;
      datD[2] = tVD3.vtZ;
    }
    node = AlcKDTGetNN(wSp->tTree, datD, wSp->maxDist, wSp->dist + idx, NULL);
    *(wSp->sNN + idx) = node->idx;
    if(wSp->vType == WLZ_VERTEX_D2)
    {
      *(wSp->nNTVx.d2 + idx) = *(wSp->gTVx.d2 + node->idx);
    }
    else /* wSp->vType == WLZ_VERTEX_D3 */
    {
      *(wSp->nNTVx.d3 + idx) = *(wSp->gTVx.d3 + node->idx);
    }
    wSp->sumDistCur += *(wSp->dist + idx);
  }
}

/*!
* \ingroup	WlzTransform
* \return	<void>
* \brief	Weights the matched verticies by combining weightings
*		for the vertex position and normal matches.
*		  
* \param	wSp			ICP registration workspace.
*/
static void	WlzRegICPWeight(WlzRegICPWSp *wSp)
/* #define WLZ_REGICP_SAMEWEIGHT */
#ifdef WLZ_REGICP_SAMEWEIGHT
{
  int		idx;

  for(idx = 0; idx < wSp->nMatch; ++idx)
  {
    *(wSp->wgtVx + idx) = wSp->kVx;
  }
  if(wSp->gSNr.v)
  {
    for(idx = 0; idx < wSp->nMatch; ++idx)
    {
      *(wSp->wgtNr +idx) = wSp->kNr;
    }
  }
}
#else /* WLZ_REGICP_SAMEWEIGHT */
{
  int		idx;
  double	tD0,
  		wVx,
		wNr,
		maxDist;
  WlzVertex	sV,
  		tV;

  /* Find the maximum distance. */
  maxDist = 1.0;
  for(idx = 0; idx < wSp->nMatch; ++idx)
  {
    if((tD0 = *(wSp->dist + idx)) > maxDist)
    {
      maxDist = tD0;
    }
  }
  /* Compute weights. */
  for(idx = 0; idx < wSp->nMatch; ++idx)
  {
    tD0 = (maxDist - *(wSp->dist + idx)) / maxDist;
    wVx = (1.0 - wSp->kVx) + (wSp->kVx * tD0);
    if(wSp->gSNr.v)
    {
      if(wSp->vType == WLZ_VERTEX_D2)
      {
        tV.d2 = *(wSp->gTNr.d2 + *(wSp->sNN + idx));
	sV.d2 = *(wSp->tSNr.d2 + idx);
	tD0 = WLZ_VTX_2_DOT(sV.d2, tV.d2);
	wNr = (1.0 - wSp->kNr) + (wSp->kNr * tD0 * tD0);
      }
      else /* wSp->vType == WLZ_VERTEX_D3 */
      {
        tV.d3 = *(wSp->gTNr.d3 + *(wSp->sNN + idx));
	sV.d3 = *(wSp->tSNr.d3 + idx);
	tD0 = WLZ_VTX_3_DOT(sV.d3, tV.d3);
	wNr = (1.0 - wSp->kNr) + (wSp->kNr * tD0 * tD0);
      }
    }
    *(wSp->wgtVx + idx) = wVx * wNr;
  }
}
#endif /* WLZ_REGICP_SAMEWEIGHT */

/*!
* \ingroup	WlzTransform
* \return		 		Woolz error code.
* \brief	Computes an affine transform from matched verticies
*		and the weights.
* \param	wSp			ICP registration workspace.
* \param	trType	 		Required transform type.
* \param	dstConv			Dst ptr for convergence flag,
*					non-zero if converged ie the
*					current affine transform is an
*					identity transform. May NOT be
*					NULL.
*/
static WlzErrorNum WlzRegICPCompTransform(WlzRegICPWSp *wSp,
				          WlzTransformType trType,
					  int *dstConv)
{
  int		idx,
		conv = 0;
  WlzAffineTransform *curTr = NULL,
  		 *newTr = NULL;
  WlzVertexP	nullP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	convThr = 0.001;
#ifdef WLZ_REGICP_DEBUG
  FILE		*fP = NULL;
#endif /* WLZ_REGICP_DEBUG */

  nullP.v = NULL;
#ifdef WLZ_REGICP_DEBUG
  if(wSp->vType == WLZ_VERTEX_D2)
  {
    if((fP = fopen("dbg-match.num", "w")) != NULL)
    {
      for(idx = 0; idx < wSp->nMatch; ++idx)
      {
        (void )fprintf(fP,
		"%g %d %g %g %g %g %4d %g %g %g %g\n",
		       *(wSp->wgtVx + idx),
		       *(wSp->sNN + idx),
		       (wSp->nNTVx.d2 + idx)->vtX,
		       (wSp->nNTVx.d2 + idx)->vtY,
		       (wSp->gTNr.d2 + *(wSp->sNN + idx))->vtX,
		       (wSp->gTNr.d2 + *(wSp->sNN + idx))->vtY,
		       idx,
		       (wSp->tSVx.d2 + idx)->vtX,
		       (wSp->tSVx.d2 + idx)->vtY,
		       (wSp->tSNr.d2 + idx)->vtX,
		       (wSp->tSNr.d2 + idx)->vtY);
      }
      (void )fclose(fP);
      fP = NULL;
    }
  }
  else /* wSp->vType == WLZ_VERTEX_D3 */
  {
    if((fP = fopen("dbg-match.num", "w")) != NULL)
    {
      for(idx = 0; idx < wSp->nMatch; ++idx)
      {
        (void )fprintf(fP,
		"%g %d %g %g %g %g %d %g %g %g %g\n",
		       *(wSp->wgtVx + idx),
		       *(wSp->sNN + idx),
		       (wSp->nNTVx.d3 + idx)->vtX,
		       (wSp->nNTVx.d3 + idx)->vtY,
		       (wSp->nNTVx.d3 + idx)->vtZ,
		       (wSp->gTNr.d3 + *(wSp->sNN + idx))->vtX,
		       (wSp->gTNr.d3 + *(wSp->sNN + idx))->vtY,
		       (wSp->gTNr.d3 + *(wSp->sNN + idx))->vtZ,
		       idx,
		       (wSp->tSVx.d3 + idx)->vtX,
		       (wSp->tSVx.d3 + idx)->vtY,
		       (wSp->tSVx.d3 + idx)->vtZ,
		       (wSp->tSNr.d3 + idx)->vtX,
		       (wSp->tSNr.d3 + idx)->vtY,
		       (wSp->tSNr.d3 + idx)->vtZ);
      }
      (void )fclose(fP);
      fP = NULL;
    }
  }
#endif /* WLZ_REGICP_DEBUG */
  /* Compute new affine trasform. */
  switch(trType)
  {
    case WLZ_TRANSFORM_2D_REG:
    case WLZ_TRANSFORM_3D_REG:
      newTr = WlzAffineTransformLSq2(wSp->vType, wSp->nMatch,
      				     wSp->wgtVx, wSp->tSVx, wSp->nNTVx,
				     0, NULL, nullP, nullP,
				     trType, &errNum);
      break;
    case WLZ_TRANSFORM_2D_AFFINE:
    case WLZ_TRANSFORM_3D_AFFINE:
      newTr = WlzAffineTransformLSqWgt(wSp->vType, wSp->nMatch, 
      				       wSp->wgtVx, wSp->tSVx, wSp->nNTVx,
				       trType, &errNum);
      break;
    default:
      errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef WLZ_REGICP_DEBUG
      (void )fprintf(stderr, "WlzRegICP newTr->mat = \n");
      (void )AlcDouble2WriteAsci(stderr, newTr->mat, 4, 4);
#endif /* WLZ_REGICP_DEBUG */
    conv = ((wSp->sumDistLst - wSp->sumDistCur) / wSp->sumDistLst) < convThr;
    if(conv == 0)
    {
      if(wSp->tr == NULL)
      {
	wSp->tr = newTr;
	newTr = NULL;
      }
      else
      {
	curTr = WlzAffineTransformProduct(wSp->tr, newTr, &errNum);
	WlzFreeAffineTransform(wSp->tr);
	wSp->tr = curTr;
	curTr = NULL;
      }
    }
#ifdef WLZ_REGICP_DEBUG
    (void )fprintf(stderr, "WlzRegICP conv = %d\n", conv);
    (void )fprintf(stderr, "WlzRegICP wSp->itr = %d\n", wSp->itr);
    (void )fprintf(stderr, "WlzRegICP wSp->tr->mat = \n");
    (void )AlcDouble2WriteAsci(stderr, newTr->mat, 4, 4);
    (void )fprintf(stderr, "\n");
#endif /* WLZ_REGICP_DEBUG */
  }
  *dstConv = conv;
  if(curTr)
  {
    WlzFreeAffineTransform(curTr);
  }
  if(newTr)
  {
    WlzFreeAffineTransform(newTr);
  }
  return(errNum);
}

/* #define WLZ_REGICP_TEST */
#ifdef WLZ_REGICP_TEST

/* Test main() for WlzRegICPObjs(). */

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
		grdFlg = 0,
  		option,
  		ok = 1,
		usage = 0;
  FILE		*fP = NULL;
  char		*outObjFileStr;
  char		*inObjFileStr[2];
  WlzValues	nullVal;
  WlzDomain	outDom;
  WlzObject	*outObj;
  WlzObject	*inObj[2];
  WlzTransformType trType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "aghro:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  inObj[0] = inObj[1] = NULL;
  nullVal.core = NULL;
  outDom.core = NULL;
  outObj = NULL;
  outObjFileStr = (char *)outFileStrDef;
  inObjFileStr[0] = inObjFileStr[1] = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        trType = WLZ_TRANSFORM_2D_AFFINE;
	break;
      case 'g':
        grdFlg = 1;
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'r':
        trType = WLZ_TRANSFORM_2D_REG;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 2) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inObjFileStr[0] = *(argv + optind);
        inObjFileStr[1] = *(argv + optind + 1);
      }
    }
  }
  if(ok)
  {
    idx = 0;
    while(ok && (idx < 2))
    {
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP, &errNum),
	  				 NULL)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read object from file %s\n",
		       *argv, inObjFileStr);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
      }
      ++idx;
    }
  }
  if(ok)
  {
    if(inObj[0]->type != inObj[1]->type)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
      ok = 0;
    }
    else
    {
      switch(inObj[0]->type)
      {
      	case WLZ_2D_DOMAINOBJ:
	  break;
      	case WLZ_3D_DOMAINOBJ:
	  switch(trType)
	  {
	    case WLZ_TRANSFORM_2D_REG:
	      trType = WLZ_TRANSFORM_3D_REG;
	      break;
	    case WLZ_TRANSFORM_2D_AFFINE:
	      trType = WLZ_TRANSFORM_3D_AFFINE;
	      break;
	  }
	  break;
	case WLZ_CONTOUR:
	  if((inObj[0]->domain.core == NULL) ||
	     (inObj[0]->domain.ctr->model == NULL))
	  {
	    errNum = WLZ_ERR_DOMAIN_NULL;
	  }
	  else
	  {
	    switch(inObj[0]->domain.ctr->model->type)
	    {
	      case WLZ_GMMOD_2I: /* FALLTHROUGH */
	      case WLZ_GMMOD_2D:
		break;
	      case WLZ_GMMOD_3I: /* FALLTHROUGH */
	      case WLZ_GMMOD_3D:
	        switch(trType)
		{
		  case WLZ_TRANSFORM_2D_REG:
		    trType = WLZ_TRANSFORM_3D_REG;
		    break;
		  case WLZ_TRANSFORM_2D_AFFINE:
		    trType = WLZ_TRANSFORM_3D_AFFINE;
		    break;
		}
		break;
	      default:
		errNum = WLZ_ERR_DOMAIN_TYPE;
		ok = 0;
		break;
	    }
	  }
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  ok = 0;
	  break;
      }
    }
  }
  if(ok)
  {
    if(grdFlg)
    {
      outDom.t = WlzRegICPObjsGrd(inObj[0], inObj[1], NULL,
				       trType, 50.0, 50.0, 1.6,
				       NULL, NULL, 200, &errNum);
    }
    else
    {
      outDom.t = WlzRegICPObjs(inObj[0], inObj[1], NULL,
			       trType, NULL, NULL, 100, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to register contours.\n",
      		     argv[0]);
    }
  }
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_AFFINE_TRANS, outDom, nullVal, NULL, NULL,
    			 &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to make affine transform object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output affine transform object "
		     "to file %s (%s).\n",
		     *argv, outObjFileStr, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  for(idx = 0; idx < 2; ++idx)
  {
    if(inObj[idx])
    {
      (void )WlzFreeObj(inObj[idx]);
    }
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(outDom.core)
  {
    (void )WlzFreeAffineTransform(outDom.t);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s",
      *argv,
      " [-o<output object>] [-a] [-h] [-g] [-r]\n"
      "          [<input object 0>] [<input object 1>]\n"
      "Options:\n"
      "  -a  Compute affine transform.\n"
      "  -g  Use maximal gradient surfaces.\n"
      "  -h  Prints this usage information.\n"
      "  -r  Compute registration transform, affine but no scale or shear.\n"
      "  -o  Output object file name.\n"
      "Computes a registration transform which registers the second of the\n"
      "given objects to the first using an ICP algorithm.\n"
      "list from the given input object.\n");
  }
  return(!ok);
}

#endif /* WLZ_REGICP_TEST */
