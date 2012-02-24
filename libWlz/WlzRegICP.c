#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRegICP_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzRegICP.c
* \author       Bill Hill
* \date         November 2000
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
* \brief	Functions for the registration of two objects using an
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
*/

#include <float.h>
#include <Wlz.h>

/*!
* \struct	_WlzRegICPWSp
* \ingroup      WlzTransform
* \brief	A workspace data structure for use in the ICP functions.
*/
typedef struct _WlzRegICPWSp
{
  /* Itteration control. */
  int		itr;		/*!< Current iteration */
  double	delta;		/*!< Convergence metric tollerance. */
  double	curMetric;	/*!< Current sum of distances between NN */
  double	prvMetric;	/*!< Last sum of distances between NN */
  /* Nearest neighbour search. */
  double	maxDist;	/*!< Maximum distance to consider for a NN */
  AlcKDTTree	*tTree;		/*!< kD-tree */
  int		*sNN;		/*!< Indicies of NN to source vertices */
  double	*dist;		/*!< NN distances */
  /* Vericies and normals. */
  WlzVertexType vType;		/*!< Type of vertices WLZ_VERTEX_D2 or
  				     WLZ_VERTEX_D3 */
  int		sgnNrm;		/*!< Non zero if normals have reliably signed
  				     components. */
  int		nT;		/*!< Number of target vertices/normals */
  int		nS;		/*!< Number of source vertices/normals */
  int		nMatch;		/*!< Minimum of nT and nS */
  WlzVertexP	gTVx;		/*!< Given target vertices */
  WlzVertexP	gSVx;		/*!< Given source vertices */
  WlzVertexP	gTNr;		/*!< Given target normals */
  WlzVertexP	gSNr;		/*!< Given source normals */
  WlzVertexP    tSVx;		/*!< Transformed source vertices. */
  WlzVertexP    tSNr;		/*!< Transformed source normals. */
  WlzVertexP    nNTVx;		/*!< NN ordered target vertices. */
  /* Match weighting */
  double	*wgt;		/*!< Weights for matches */
  /* Affine transform. */
  WlzAffineTransform *prvTr;	/*!< Previous affine transform */
  WlzAffineTransform *curTr;	/*!< Current affine transform */
}  WlzRegICPWSp;

static void     		WlzRegICPTrans(
				  WlzRegICPWSp *wSp);
static void			WlzRegICPFindNN(
				  WlzRegICPWSp *wSp);
static int			WlzRegICPItr(
				  WlzRegICPWSp *wSp,
			          WlzTransformType trType,
				  int maxItr,
				  double minDistWgt,
				  WlzErrorNum *dstErr);
static double			WlzRegICPWeight(
				  WlzRegICPWSp *wSp,
				  double minDistWgt);
static WlzErrorNum		WlzRegICPCompTransform(
				  WlzRegICPWSp *wSp,
				  WlzTransformType trType);
static WlzErrorNum 		WlzRegICPCheckVertices(
				  WlzVertexP *vData,
				  int *vCnt,
				  WlzVertexType *vType);
static WlzErrorNum 		WlzRegICPBuildTree(
				  WlzRegICPWSp *wSp);
static WlzAffineTransform 	*WlzRegICPTreeAndVerticesSimple(
				  AlcKDTTree *tree,
				  WlzTransformType trType,
				  WlzVertexType vType,
				  int sgnNrm,
				  int nT,
				  WlzVertexP tVx,
				  WlzVertexP tNr,
				  int nS,
				  int *sIdx,
				  WlzVertexP sVx,
				  WlzVertexP sNr,
				  WlzVertexP tVxBuf,
				  WlzVertexP sVxBuf,
				  double *wgtBuf,
				  int maxItr,
				  WlzAffineTransform *initTr,
				  double *prvMetric,
				  double *curMetric,
				  int *dstConv,
				  WlzRegICPUsrWgtFn usrWgtFn,
				  void *usrWgtData,
				  double delta,
				  double minDistWgt,
				  WlzErrorNum *dstErr);

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
* \param	delta			Tolerance for mean value of
*					registration metric.
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzAffineTransform *WlzRegICPObjsGrd(WlzObject *tObj, WlzObject *sObj,
				     WlzAffineTransform *initTr,
				     WlzTransformType trType,
				     double ctrLo, double ctrHi, double ctrWth,
				     int *dstConv, int *dstItr, int maxItr,
				     double delta, double minDistWgt,
				     WlzErrorNum *dstErr)
{
  int		idN,
		conv,
		itr;
  int		vCnt[2];
  WlzVertexType	vType;
  WlzDomain	rDom;
  WlzValues	rVal;
  WlzObject	*cObj;
  WlzObject	*gObj[2];
  WlzVertexP	nData[2],
  		vData[2];
  WlzAffineTransform *regTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
#ifdef WLZ_REGICP_DEBUG
  int		idM;
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
   * vertices and their gradient normals. */
  /* For the target object. */
  for(idN = 0; idN < 2; ++idN)
  {
    cObj = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      rDom.ctr = WlzContourObjGrd(gObj[idN], ctrLo, ctrHi, ctrWth,
				  1, &errNum);
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
    /* Get an array of position vertices and normals from the model,
     * don't use the gradient normals because they are eratic. */
    if(errNum == WLZ_ERR_NONE)
    {
      vData[idN] = WlzVerticesFromObj(cObj, nData + idN, vCnt + idN, &vType,
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
  /* Register the position vertices and normals. */
  if(errNum == WLZ_ERR_NONE)
  {
    regTr = WlzRegICPVertices(vData[0], nData[0], vCnt[0],
			       vData[1], nData[1], vCnt[1],
			       vType, 1, initTr,
			       trType, &conv, &itr, maxItr,
			       delta, minDistWgt, &errNum);
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
* \return				Affine transform which brings
*					the two objects into register.
* \ingroup	WlzTransform
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
* \param	delta			Tolerance for mean value of
*					registration metric.
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzAffineTransform *WlzRegICPObjs(WlzObject *tObj, WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv, int *dstItr, int maxItr,
				  double delta, double minDistWgt,
				  WlzErrorNum *dstErr)
{
  int		idx,
		conv,
		itr,
		sgnNrm;
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
    sgnNrm = 0;
    if((tObj->type == WLZ_CONTOUR) && (sObj->type == WLZ_CONTOUR) &&
       (sObj->domain.ctr->model != NULL) && (tObj->domain.ctr->model != NULL))
    {
      switch(tObj->domain.ctr->model->type)
      {
        case WLZ_GMMOD_2N:
	case WLZ_GMMOD_3N:
	  sgnNrm = 1;
	  break;
        default:
	  errNum = WLZ_ERR_PARAM_TYPE;
	  break;
      }
    }
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      vData[idx] = WlzVerticesFromObj(objs[idx], nData + idx, vCnt + idx,
      				       vType + idx, &errNum);
      ++idx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check type of vertices and convert to double if required. */
    errNum = WlzRegICPCheckVertices(vData, vCnt, vType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
     regTr = WlzRegICPVertices(vData[0], nData[0], vCnt[0],
     				vData[1], nData[1], vCnt[1],
     				vType[0], sgnNrm, initTr,
				trType, &conv, &itr, maxItr,
				delta, minDistWgt, &errNum);
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
* \return				3D Woolz domain object.
* \ingroup	WlzTransform
* \brief	Computes a 3D domain object in which the double
*		precission values are the weighted sums of distances
*		as used by WlzRegICPObjs().
*		Thehese values are computed for the given range
*		of translations and rotations.
* \param	tObj			2D target object.
* \param	sObj			2D source object to be
*					registered with target object.
* \param	initTr			Initial affine transform
*					to be applied to the source
*					object prior to using the ICP
*					algorithm. May be NULL.
* \param	xMin			Minimum column translation.
* \param	xMax			Maximum column translation.
* \param	xStep			Increment for column translation.
* \param	yMin			Minimum line translation.
* \param	yMax			Maximum line translation.
* \param	yStep			Increment for line translation.
* \param	rMin			Minimum rotation (radians).
* \param	rMax			Maximum rotation (radians).
* \param	rStep			Increment for rotation (radians).
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzObject	*WlzRegICPObjWSD2D(WlzObject *tObj, WlzObject *sObj,
				   WlzAffineTransform *initTr,
				   double xMin, double xMax, double xStep,
				   double yMin, double yMax, double yStep,
				   double rMin, double rMax, double rStep,
				   double minDistWgt,
				   WlzErrorNum *dstErr)
{
  int		idx,
		sgnNrm;
  int		vCnt[2];
  WlzVertexType	vType[2];
  WlzObject	*objs[2];
  WlzVertexP	nData[2],
  		vData[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*mObj = NULL;

  nData[0].v = nData[1].v = NULL;
  vData[0].v = vData[1].v = NULL;
  objs[0] = tObj; objs[1] = sObj;
  if((tObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    sgnNrm = 0;
    if((tObj->type == WLZ_CONTOUR) && (sObj->type == WLZ_CONTOUR) &&
       (sObj->domain.ctr->model != NULL) && (tObj->domain.ctr->model != NULL))
    {
      switch(tObj->domain.ctr->model->type)
      {
        case WLZ_GMMOD_2N:
	case WLZ_GMMOD_3N:
	  sgnNrm = 1;
	  break;
        default:
	  errNum = WLZ_ERR_PARAM_TYPE;
	  break;
      }
    }
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      vData[idx] = WlzVerticesFromObj(objs[idx], nData + idx, vCnt + idx,
      				       vType + idx, &errNum);
      ++idx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check type of vertices and convert to double if required. */
    errNum = WlzRegICPCheckVertices(vData, vCnt, vType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
     mObj = WlzRegICPVerticesWSD2D(vData[0], nData[0], vCnt[0],
				   vData[1], nData[1], vCnt[1],
				   vType[0], sgnNrm, initTr,
				   xMin, xMax, xStep,
				   yMin, yMax, yStep,
				   rMin, rMax, rStep,
				   minDistWgt,
				   &errNum);
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
  return(mObj);
}

/*!
* \return				Woolz error code.
* \ingroup	WlzTransform
* \brief	Checks the two sets of vertices and promotes the
*		vertex type to double (2 or 3D) if not already double
*		vertices.
* \param	vData			The sets of vertices.
* \param	vCnt			The numbers of vertices.
* \param	vType			The types of vertices.
*/
static WlzErrorNum WlzRegICPCheckVertices(WlzVertexP *vData, int *vCnt,
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
	      default:
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
	      default:
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
* \return				Affine transform which brings
*					the two sets of vertices into
*					register.
* \ingroup	WlzTransform
* \brief	Registers the two given sets of vertices using the
*		iterative closest point algorithm. An affine transform
*		is computed, which when applied to the source vertices
*		takes it into register with the target vertices.
*		The vertices and their normals are known to be either
*		WlzDVertex2 or WlzDVertex3.
* \param	tVx			Target vertices.
* \param	tNr			Target normals, may be NULL.
* \param	tCnt			Number of target vertices.
* \param	sVx			Source vertices.
* \param	sNr			Source normals, may be NULL.
* \param	sCnt			Number of source vertices.
* \param	vType			Type of the vertices.
* \param	sgnNrm			Non zero if the normals have reliably
*					signed components.
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
* \param	delta			Tolerance for mean value of
*					registration metric.
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzAffineTransform	*WlzRegICPVertices(WlzVertexP tVx, WlzVertexP tNr,
					    int tCnt,
					    WlzVertexP sVx, WlzVertexP sNr,
					    int sCnt,
					    WlzVertexType vType, int sgnNrm,
					    WlzAffineTransform *initTr,
					    WlzTransformType trType,
					    int *dstConv, int *dstItr,
					    int maxItr,
				     	    double delta,
					    double minDistWgt,
					    WlzErrorNum *dstErr)
{
  int		conv,
		maxCnt;
  WlzAffineTransform *regTr = NULL;
  WlzRegICPWSp	wSp;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  /* Setup workspace. */
  wSp.itr = 0;
  wSp.delta = delta;
  wSp.curMetric = DBL_MAX;
  wSp.prvMetric = DBL_MAX;
  wSp.maxDist = DBL_MAX;
  wSp.tTree = NULL;
  wSp.sNN = NULL;
  wSp.dist = NULL;
  wSp.vType = vType;
  wSp.sgnNrm = sgnNrm;
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
  wSp.wgt = NULL;
  wSp.prvTr = NULL;
  wSp.curTr = NULL;
  maxCnt = WLZ_MAX(tCnt, sCnt);
  if(((wSp.sNN = (int *)AlcMalloc(sizeof(int) * maxCnt)) == NULL) ||
     ((wSp.dist = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL) ||
     ((wSp.wgt = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(initTr)
    {
      wSp.curTr = WlzAffineTransformCopy(initTr, &errNum);
    }
    else
    {
      if(vType == WLZ_VERTEX_D2)
      {
	wSp.curTr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum);
      }
      else /* vType == WLZ_VERTEX_D3 */
      {
        wSp.curTr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, &errNum);
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
    /* Build k-D tree for the target vertices. */
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
				   maxItr, minDistWgt, &errNum);
	if(conv && (errNum == WLZ_ERR_NONE))
	{
          conv = WlzRegICPItr(&wSp, trType, maxItr, minDistWgt, &errNum);
	}
        break;
      case WLZ_TRANSFORM_2D_REG:
      case WLZ_TRANSFORM_3D_REG:
        conv = WlzRegICPItr(&wSp, trType, maxItr, minDistWgt, &errNum);
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
  AlcFree(wSp.wgt);
  AlcFree(wSp.tSVx.v);
  AlcFree(wSp.tSNr.v);
  AlcFree(wSp.nNTVx.v);
  (void )WlzFreeAffineTransform(wSp.prvTr);
  if(wSp.curTr)
  {
    if((errNum == WLZ_ERR_NONE) && conv)
    {
      regTr = wSp.curTr;
    }
    else
    {
      (void )WlzFreeAffineTransform(wSp.curTr);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(regTr);
}

/*!
* \return				3D Woolz domain object.
* \ingroup	WlzTransform
* \brief	Computes a 3D domain object in which the double
*		precission values are the weighted sums of distances
*		as used by WlzRegICPObjs().
*		Thehese values are computed for the given range
*		of translations and rotations.
* \param	tVx			Target vertices.
* \param	tNr			Target normals, may be NULL.
* \param	tCnt			Number of target vertices.
* \param	sVx			Source vertices.
* \param	sNr			Source normals, may be NULL.
* \param	sCnt			Number of source vertices.
* \param	vType			Type of the vertices.
* \param	sgnNrm			Non zero if the normals have reliably
*					signed components.
* \param	initTr			Initial affine transform
*					to be applied to the source
*					object prior to using the ICP
*					algorithm. May be NULL.
* \param	xMin			Minimum column translation.
* \param	xMax			Maximum column translation.
* \param	xStep			Increment for column translation.
* \param	yMin			Minimum line translation.
* \param	yMax			Maximum line translation.
* \param	yStep			Increment for line translation.
* \param	rMin			Minimum rotation (radians).
* \param	rMax			Maximum rotation (radians).
* \param	rStep			Increment for rotation (radians).
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzObject	*WlzRegICPVerticesWSD2D(WlzVertexP tVx, WlzVertexP tNr,
				int tCnt,
				WlzVertexP sVx, WlzVertexP sNr, int sCnt,
				WlzVertexType vType, int sgnNrm,
				WlzAffineTransform *initTr,
				double xMin, double xMax, double xStep,
				double yMin, double yMax, double yStep,
				double rMin, double rMax, double rStep,
				double minDistWgt,
				WlzErrorNum *dstErr)
{
  int		idX,
  		idY,
		idZ,
		maxCnt;
  WlzObject	*mObj = NULL;
  WlzRegICPWSp	wSp;
  WlzIVertex3	mOrg,
  		mSz;
  WlzAffineTransform *tTr = NULL;
  double	***mAry = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  /* Setup workspace. */
  wSp.itr = 0;
  wSp.curMetric = DBL_MAX;
  wSp.prvMetric = DBL_MAX;
  wSp.maxDist = DBL_MAX;
  wSp.tTree = NULL;
  wSp.sNN = NULL;
  wSp.dist = NULL;
  wSp.vType = vType;
  wSp.sgnNrm = sgnNrm;
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
  wSp.wgt = NULL;
  wSp.curTr = NULL;
  maxCnt = WLZ_MAX(tCnt, sCnt);
  if((fabs(xStep) < DBL_EPSILON) ||
     (fabs(yStep) < DBL_EPSILON) ||
     (fabs(rStep) < DBL_EPSILON))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((wSp.sNN = (int *)AlcMalloc(sizeof(int) * maxCnt)) == NULL) ||
       ((wSp.dist = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL) ||
       ((wSp.wgt = (double *)AlcMalloc(sizeof(double) * maxCnt)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
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
  if(errNum == WLZ_ERR_NONE)
  {
    mOrg.vtX = WLZ_NINT(xMin / xStep);
    mOrg.vtY = WLZ_NINT(yMin / yStep);
    mOrg.vtZ = WLZ_NINT(rMin / rStep);
    mSz.vtX = WLZ_NINT(xMax / xStep) - mOrg.vtX;
    mSz.vtY = WLZ_NINT(yMax / yStep) - mOrg.vtY;
    mSz.vtZ = WLZ_NINT(rMax / rStep) - mOrg.vtZ;
    if(AlcDouble3Calloc(&mAry, (size_t )mSz.vtZ,
    		        (size_t )mSz.vtY, (size_t )mSz.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = WlzFromArray3D((void ***)mAry,
                          mSz, mOrg, WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE,
    			  0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Build k-D tree for the target vertices. */
    errNum = WlzRegICPBuildTree(&wSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idZ = 0; (idZ < mSz.vtZ) && (errNum == WLZ_ERR_NONE); ++idZ)
    {
      for(idY = 0; (idY < mSz.vtY) && (errNum == WLZ_ERR_NONE); ++idY)
      {
	for(idX = 0; (idX < mSz.vtX) && (errNum == WLZ_ERR_NONE); ++idX)
	{
	  tTr = WlzAssignAffineTransform(
	        WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
	  				(mOrg.vtX + idX) * xStep,
	  				(mOrg.vtY + idY) * yStep,
					0.0,
					1.0,
	  				(mOrg.vtZ + idZ) * rStep,
					0.0, 0.0, 0.0, 0.0, 0,
					&errNum), NULL);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    (void )WlzFreeAffineTransform(wSp.curTr);
	    if(initTr)
	    {
	      wSp.curTr = WlzAssignAffineTransform(
	      	       WlzAffineTransformProduct(tTr, initTr, &errNum), NULL);
	    }
	    else
	    {
	      wSp.curTr = WlzAssignAffineTransform(tTr, NULL);
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzRegICPTrans(&wSp);
	    WlzRegICPFindNN(&wSp);
	    /* Compute weighted sum of distances for mObj. */
	    mAry[idZ][idY][idX] = WlzRegICPWeight(&wSp, minDistWgt);
	  }
	}
      }
    }
  }
  AlcFree(wSp.sNN);
  AlcFree(wSp.dist);
  AlcFree(wSp.wgt);
  AlcFree(wSp.tSVx.v);
  AlcFree(wSp.tSNr.v);
  AlcFree(wSp.nNTVx.v);
  (void )WlzFreeAffineTransform(wSp.curTr);
  if(errNum != WLZ_ERR_NONE)
  {
    if(mObj)
    {
      (void )WlzFreeObj(mObj);
    }
    else if(mAry)
    {
      (void )Alc3Free((void ***)mAry);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return				Woolz error code
* \ingroup	WlzTransform
* \brief	Allocates and populates a k-D tree from the given vertices.
* 		The vertices are either WlzDVertex2 orWlzDVertex3
* \param	wSp			ICP registration workspace.
*/
static WlzErrorNum WlzRegICPBuildTree(WlzRegICPWSp *wSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Use wSp->sNN as a temporay buffer of shuffling indcies. */
  wSp->tTree = WlzVerticesBuildTree(wSp->vType, wSp->nT, wSp->gTVx,
  				     wSp->sNN, &errNum);
  return(errNum);
}

/*!
* \return				Nonzero if the iteration has
* 					converged.
* \ingroup	WlzTransform
* \brief	The iterative loop of the ICP, which iterates to find the
* 		registration transform.
* \param	wSp			ICP registration workspace.
* \param	trType			Type of affine transform.
* \param	maxItr			Maximum number of iterations.
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static int	WlzRegICPItr(WlzRegICPWSp *wSp,
			     WlzTransformType trType, int maxItr,
			     double minDistWgt, WlzErrorNum *dstErr)
{
  int		conv = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  do
  {
#ifdef WLZ_REGICP_DEBUG
    (void )fprintf(stderr, "WlzRegICP wSp->itr = %d\n", wSp->itr);
#endif /* WLZ_REGICP_DEBUG */
    wSp->prvMetric = wSp->curMetric;
    /* Apply current transform to the vertices and normals. */
    WlzRegICPTrans(wSp);
    /* Find closest points. */
    WlzRegICPFindNN(wSp);
    /* Compute weightings of the matches and current metric. */
    wSp->curMetric = WlzRegICPWeight(wSp, minDistWgt);
    /* Check for convergence. */
    conv = (wSp->prvMetric - wSp->curMetric) < wSp->delta;
#ifdef WLZ_REGICP_DEBUG
    (void )fprintf(stderr, "WlzRegICP conv = %s\n", (conv)? "TRUE": "FALSE");
#endif /* WLZ_REGICP_DEBUG */
    if(conv)
    {
      (void )WlzFreeAffineTransform(wSp->curTr);
      wSp->curTr = wSp->prvTr;
      wSp->prvTr = NULL;
    }
    else
    {
      /* Compute registration transform. */
      errNum = WlzRegICPCompTransform(wSp, trType);
    }
  } while((errNum == WLZ_ERR_NONE) && (wSp->itr++ < maxItr) && (conv == 0));
  if(wSp->itr > maxItr)
  {
    errNum = WLZ_ERR_ALG_CONVERGENCE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(conv);
}

/*!
* \return       void
* \ingroup	WlzTransform
* \brief        Transforms the source vertices and normals using the
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
      *(wSp->tSVx.d2 + idx) = WlzAffineTransformVertexD2(wSp->curTr,
      						*(wSp->gSVx.d2 + idx), NULL);
    }
    if(wSp->gSNr.v)
    {
      for(idx = 0; idx < wSp->nS; ++idx)
      {
        *(wSp->tSNr.d2 + idx) = WlzAffineTransformNormalD2(wSp->curTr,
						*(wSp->gSNr.d2 + idx), NULL);
      }
    }
  }
  else /* wSp->vType == WLZ_VERTEX_D3 */
  {
    for(idx = 0; idx < wSp->nS; ++idx)
    {
      *(wSp->tSVx.d3 + idx) = WlzAffineTransformVertexD3(wSp->curTr,
      						*(wSp->gSVx.d3 + idx), NULL);
    }
    if(wSp->gSNr.v)
    {
      for(idx = 0; idx < wSp->nS; ++idx)
      {
        *(wSp->tSNr.d3 + idx) = WlzAffineTransformNormalD3(wSp->curTr,
						*(wSp->gSNr.d3 + idx), NULL);
      }
    }
  }
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Finds nearest neighbour matches in the target tree for
*		the source vertices, sets the nearest neighbour
*		indicies and permutes the NN ordered target vertices in
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
  }
}

/*!
* \return	Mean of the weighted distance measures.
* \ingroup	WlzTransform
* \brief	Weights the matched vertices by combining weightings
*		for the vertex position and normal matches.
*
*		the weight function callbacks are not used in this
*		function, instead the weights \f$w\f$ are simply computed
*		as the product of seperate vertex distance and normal
*		product weights: \f$w_v\f$ and \f$w_n\f$. With:
*		\f[
		w_v = \left\{
		      \begin{array}{ll}
		      W_{min} & \textrm{if maximum distance} \\
		      1.0 & \textrm{if minimum distance}
		      \end{array} \right.
		\f]
*		and
*		\f[
		w_n = \left\{
		      \begin{array}{ll}
		      \acute{w_n} & \textrm{for reliably signed normals} \\
		      (\acute{w_n})^2 & \textrm{for unreliably signed normals}
		      \end{array} \right.
		\f]
*		here the normal weight \f$\acute{w_n}\f$ is given by
*		the scalar product of the normals clamped to the
*		range 0.0 - 1.0, ie if \f$n \cdot n < 0\f$ then
*		\f$\acute{w_n} = 0\f$.
*
*		The minimum distance weighting \f$W_{min}\f$ can
*		be used to control the spatial localisation of the
*		registration, with \f$W_{min} = 0\f$ giving good
*		localisation and \f$W_{min} = 0.25\f$ giving a more
*		global registration.
* \param	wSp			ICP registration workspace.
* \param	minVxWgt		Minimum distance weighting
* 					\f$W_{min}\f$, range [0-1].
*/
static double	WlzRegICPWeight(WlzRegICPWSp *wSp, double minVxWgt)
{
  int		idx;
  double	tD0,
		w0,
		w1,
		w2,
  		wVx,
		wNr,
		minDist,
		maxDist,
		meanSumWgt = 0.0;
  WlzVertex	sV,
  		tV;

  /* Find the maximum and minimum distances. */
  minDist = maxDist = *(wSp->dist + 0);
  for(idx = 1; idx < wSp->nMatch; ++idx)
  {
    if((tD0 = *(wSp->dist + idx)) < minDist)
    {
      minDist = tD0;
    }
    else if(tD0 > maxDist)
    {
      maxDist = tD0;
    }
  }
  /* Compute weights. */
  w0 = maxDist - minDist;
  w1 = 1.0 - minVxWgt;
  w2 = (w0 > DBL_EPSILON)? w1 / w0: 1.0;
  for(idx = 0; idx < wSp->nMatch; ++idx)
  {
    /* Use linear weighting for distance such that:
     *   w = minVxWgt, d = maxDist
     * and
     *   w = 1.0, d = minDist
     */
    if(w0 > DBL_EPSILON)
    {
      wVx = 1.0 - w2 * (*(wSp->dist + idx) - minDist);
    }
    else
    {
      wVx = 1.0;
    }
    if(wSp->gSNr.v)
    {
      if(wSp->vType == WLZ_VERTEX_D2)
      {
        tV.d2 = *(wSp->gTNr.d2 + *(wSp->sNN + idx));
	sV.d2 = *(wSp->tSNr.d2 + idx);
	wNr = WLZ_VTX_2_DOT(sV.d2, tV.d2);
      }
      else /* wSp->vType == WLZ_VERTEX_D3 */
      {
        tV.d3 = *(wSp->gTNr.d3 + *(wSp->sNN + idx));
	sV.d3 = *(wSp->tSNr.d3 + idx);
	wNr = WLZ_VTX_3_DOT(sV.d3, tV.d3);
      }
      if(wSp->sgnNrm && (wNr < 0.0))
      {
        wNr = 0.0;
      }
    }
    if(wSp->sgnNrm)
    {
      tD0 = wVx * wNr;
    }
    else
    {
      tD0 = wVx * wNr * wNr;
    }
    meanSumWgt += tD0 * *(wSp->dist + idx);
    *(wSp->wgt + idx) = tD0;
  }
  meanSumWgt /= wSp->nMatch;
  return(meanSumWgt);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes an affine transform from matched vertices
*		and the weights.
* \param	wSp			ICP registration workspace.
* \param	trType	 		Required transform type.
*/
static WlzErrorNum WlzRegICPCompTransform(WlzRegICPWSp *wSp,
				          WlzTransformType trType)
{
  WlzAffineTransform *newTr = NULL;
  WlzVertexP	nullP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
#ifdef WLZ_REGICP_DEBUG
  int		idx;
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
		       *(wSp->wgt + idx),
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
		       *(wSp->wgt + idx),
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
  if(errNum == WLZ_ERR_NONE)
  {
    newTr = WlzAffineTransformLSq(wSp->vType, wSp->nMatch, wSp->nNTVx,
				  wSp->nMatch, wSp->tSVx,
				  wSp->nMatch, wSp->wgt,
				  trType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef WLZ_REGICP_DEBUG
    (void )fprintf(stderr, "WlzRegICP wSp->prvMetric = %g\n", wSp->prvMetric);
    (void )fprintf(stderr, "WlzRegICP wSp->curMetric = %g\n", wSp->curMetric);
    (void )fprintf(stderr, "WlzRegICP newTr->mat = \n");
    (void )AlcDouble2WriteAsci(stderr, newTr->mat, 4, 4);
#endif /* WLZ_REGICP_DEBUG */
    (void )WlzFreeAffineTransform(wSp->prvTr);
    wSp->prvTr = wSp->curTr;
    wSp->curTr = WlzAffineTransformProduct(newTr, wSp->prvTr, &errNum);
#ifdef WLZ_REGICP_DEBUG
    (void )fprintf(stderr, "WlzRegICP wSp->curTr->mat = \n");
    (void )AlcDouble2WriteAsci(stderr, wSp->curTr->mat, 4, 4);
    (void )fprintf(stderr, "\n");
#endif /* WLZ_REGICP_DEBUG */
  }
  (void )WlzFreeAffineTransform(newTr);
  return(errNum);
}

/*!
* \return	Affine transform found.
* \ingroup	WlzTransform
* \brief	Registers the given vertices using the already built
*		kD-tree and the given buffers.
*		This function will attempt to find a rigid body registration
*		before attempting a general affine registration.
* \param	tree			Given kD-tree populated by the
*					target vertices such that the
*					nodes of the tree have the same
*					indicies as the given target vertices
*					and normals.
* \param	trType			The required type of transform,
*					must be either WLZ_TRANSFORM_2D_REG,
*					or WLZ_TRANSFORM_2D_AFFINE.
* \param	vType			Vertex type.
* \param	sgnNrm			Non zero if sign of normal components
* 					is meaningful.
* \param	nT			Number of target vertices.
* \param	tVx			The target vertices.
* \param	tNr			The target normals.
* \param        nS			Number of source vertices.
* \param	sIdx			Indicies of the source
*					vertices/normals.
* \param        sVx 			The source vertices.
* \param	sNr			The source normals.
* \param	tVxBuf			A buffer with room for at least
*					nS vertices.
* \param	sVxBuf			A buffer with room for at least
*					nS vertices.
* \param	wgtBuf			A buffer with room for at least
*					nS doubles.
* \param	maxItr			Maximum number of iterations.
* \param	initTr			Initial affine transform.
* \param	dstConv			Destination pointer for a
*					convergence flag which is set to
*					a non zero value if the registration
*					converges.
* \param	usrWgtFn		User supplied weight function, may be
* 					NULL.
* \param	usrWgtData		User supplied weight data, may be NULL.
* \param	delta			Tolerance for mean value of
*					registration metric.
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzAffineTransform *WlzRegICPTreeAndVertices(AlcKDTTree *tree,
				WlzTransformType trType,
				WlzVertexType vType, int sgnNrm,
				int nT, WlzVertexP tVx, WlzVertexP tNr,
				int nS, int *sIdx,
				WlzVertexP sVx, WlzVertexP sNr,
				WlzVertexP tVxBuf, WlzVertexP sVxBuf,
				double *wgtBuf, int maxItr,
				WlzAffineTransform *initTr, int *dstConv,
				WlzRegICPUsrWgtFn usrWgtFn, void *usrWgtData,
				double delta, double minDistWgt,
				WlzErrorNum *dstErr)
{
  int		conv = 0;
  double	prvMetric = DBL_MAX,
  		curMetric = DBL_MAX;
  WlzTransformType trType0;
  WlzAffineTransform *newTr0 = NULL,
  		*newTr1 = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(trType)
  {
    case WLZ_TRANSFORM_2D_REG: /* FALLTHROUGH */
    case WLZ_TRANSFORM_2D_AFFINE:
      trType0 = WLZ_TRANSFORM_2D_REG;
      break;
    case WLZ_TRANSFORM_3D_REG: /* FALLTHROUGH */
    case WLZ_TRANSFORM_3D_AFFINE:
      trType0 = WLZ_TRANSFORM_3D_REG;
      break;
    default:
      errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newTr0 = WlzRegICPTreeAndVerticesSimple(tree, trType0, vType, sgnNrm,
    					 nT, tVx, tNr, nS, sIdx, sVx, sNr,
					 tVxBuf, sVxBuf, wgtBuf,
					 maxItr, initTr,
					 &prvMetric, &curMetric, &conv,
					 usrWgtFn, usrWgtData,
					 delta, minDistWgt,
					 &errNum);
#ifdef WLZ_REGICP_DEBUG
    (void )fprintf(stderr, "WlzRegICP newTr0->mat = \n");
    (void )AlcDouble2WriteAsci(stderr, newTr0->mat, 4, 4);
#endif /* WLZ_REGICP_DEBUG */
  }
  if((errNum == WLZ_ERR_NONE) && conv && (trType != trType0))
  {
    newTr1 = WlzRegICPTreeAndVerticesSimple(tree, trType, vType, sgnNrm,
					 nT, tVx, tNr, nS, sIdx, sVx, sNr,
					 tVxBuf, sVxBuf, wgtBuf,
					 maxItr, newTr0,
					 &prvMetric, &curMetric, &conv,
					 usrWgtFn, usrWgtData,
					 delta, minDistWgt,
					 &errNum);
    (void )WlzFreeAffineTransform(newTr0);
    newTr0 = newTr1;
#ifdef WLZ_REGICP_DEBUG
      (void )fprintf(stderr, "WlzRegICP newTr0->mat = \n");
      (void )AlcDouble2WriteAsci(stderr, newTr0->mat, 4, 4);
#endif /* WLZ_REGICP_DEBUG */
  }
  if(dstConv)
  {
    *dstConv = conv;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newTr0);
}

/*!
* \return				Affine transform found.
* \ingroup	WlzTransform
* \brief	Registers the given vertices using the already built
*		kD-tree and the given buffers. Unlike
*		WlzRegICPTreeAndVertices() this function does not
*		perform a registration transform before a general
*		affine transform.
* \param	tree			Given kD-tree populated by the
*					target vertices such that the
*					nodes of the tree have the same
*					indicies as the given target vertices
*					and normals.
* \param	trType			The required type of transform,
*					must be either WLZ_TRANSFORM_2D_REG,
*					or WLZ_TRANSFORM_2D_AFFINE.
* \param	vType			Type of vertex.
* \param	sgnNrm			Non zero if sign of normal components
* 					is meaningful.
* \param	nT			Number of target vertices.
* \param	tVx			The target vertices.
* \param	tNr			The target normals.
* \param        nS			Number of source vertices.
* \param	sIdx			Indicies of the source
*					vertices/normals.
* \param        sVx 			The source vertices.
* \param	sNr			The source normals.
* \param	tVxBuf			A buffer with room for at least
*					nS vertices. Used for target vertices.
* \param	sTVxBuf			A buffer with room for at least
*					nS vertices. Used for transformed
*					source vertices.
* \param	wgtBuf			A buffer with room for at least
*					nS doubles.
* \param	maxItr			Maximum number of iterations.
* \param	initTr			Initial affine transform.
* \param	gPrvMetric		Previous value of registration metric.
* \param	gCurMetric		Current value of registration metric.
* \param	dstConv			Destination pointer for a
*					convergence flag which is set to
*					a non zero value if the registration
*					converges.
* \param	usrWgtFn		User supplied weight function, may be
* 					NULL.
* \param	usrWgtData		User supplied weight data, may be NULL.
* \param	delta			Tolerance for mean value of
*					registration metric.
* \param	minDistWgt		Minimum distance weighting.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
static WlzAffineTransform *WlzRegICPTreeAndVerticesSimple(AlcKDTTree *tree,
				WlzTransformType trType,
				WlzVertexType vType, int sgnNrm,
				int nT, WlzVertexP tVx, WlzVertexP tNr,
				int nS, int *sIdx,
				WlzVertexP sVx, WlzVertexP sNr,
				WlzVertexP tVxBuf, WlzVertexP sTVxBuf,
				double *wgtBuf, int maxItr,
				WlzAffineTransform *initTr,
				double *gPrvMetric, double *gCurMetric,
				int *dstConv,
				WlzRegICPUsrWgtFn usrWgtFn, void *usrWgtData,
				double delta, double minDistWgt,
				WlzErrorNum *dstErr)
{
  int		idS,
  		idM,
		idV,
		itr = 0,
		conv = 0;
  AlcKDTNode	*tNode;
  WlzAffineTransform *invTr = NULL,
		*prvTr = NULL,
  		*curTr = NULL,
  		*newTr = NULL;
  WlzVertexP	nullP;
  WlzVertex	dV,
  		sV,
		sN,
		sTV,
		sTN,
  		tV,
		tN;
  double	wgt0,
		wgt1,
		wgt2,
		wMaxDist,
		wMinDist,
		prvMetric,
		curMetric,
		dist,
  		wNr,
		wVx;
  double	vxD[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
 
  nullP.v = NULL;
  curTr = (initTr == NULL)?
  	  WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum):
	  WlzAffineTransformCopy(initTr, &errNum);
  curMetric = *gCurMetric;
  if(errNum == WLZ_ERR_NONE)
  {
    invTr = WlzAffineTransformInverse(curTr, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    do
    {
      idM = 0;
      wMinDist = DBL_MAX;
      wMaxDist = 0.0;
      prvMetric = curMetric;
      curMetric = 0.0;
      /* Populate the buffers with source vertices, nearest neighbours
       * in the target tree and scalar product of vertex normals. While
       * doinf this also find maximum and minimum source - target vertex
       * distances */
      for(idS = 0; idS < nS; ++idS)
      {
	idV = *(sIdx + idS);
	if(vType == WLZ_VERTEX_D2)
	{
	  sV.d2 = *(sVx.d2 + idV);
	  sN.d2 = *(sNr.d2 + idV);
	  sTV.d2 = WlzAffineTransformVertexD2(curTr, sV.d2, NULL);
	  sTN.d2 = WlzAffineTransformNormalD2(curTr, sN.d2, NULL);
	  *(sTVxBuf.d2 + idM) = sTV.d2;
	  vxD[0] = sTV.d2.vtX;
	  vxD[1] = sTV.d2.vtY;
	}
	else /* vType == WLZ_VERTEX_D3 */
	{
	  sV.d3 = *(sVx.d3 + idV);
	  sN.d3 = *(sNr.d3 + idV);
	  sTV.d3 = WlzAffineTransformVertexD3(curTr, sV.d3, NULL);
	  sTN.d3 = WlzAffineTransformNormalD3(curTr, sN.d3, NULL);
	  *(sTVxBuf.d3 + idM) = sTV.d3;
	  vxD[0] = sTV.d3.vtX;
	  vxD[1] = sTV.d3.vtY;
	  vxD[2] = sTV.d3.vtZ;
	}
	if((tNode = AlcKDTGetNN(tree, vxD, DBL_MAX, &dist, NULL)) != NULL)
	{
	  if(wMinDist > dist)
	  {
	    wMinDist = dist;
	  }
	  else if(wMaxDist < dist)
	  {
	    wMaxDist = dist;
	  }
	  if(vType == WLZ_VERTEX_D2)
	  {
	    tV.d2 = *(tVx.d2 + tNode->idx);
	    tN.d2 = *(tNr.d2 + tNode->idx);
	    *(tVxBuf.d2 + idM) = tV.d2;
	    *(wgtBuf + idM) = WLZ_VTX_2_DOT(sTN.d2, tN.d2);
	  }
	  else /* vType == WLZ_VERTEX_D3 */
	  {
	    tV.d3 = *(tVx.d3 + tNode->idx);
	    tN.d3 = *(tNr.d3 + tNode->idx);
	    *(tVxBuf.d3 + idM) = tV.d3;
	    *(wgtBuf + idM) = WLZ_VTX_3_DOT(sTN.d3, tN.d3);
	  }
	  ++idM;
	}
      }
      if(idM == 0)
      {
        errNum = WLZ_ERR_ALG_CONVERGENCE;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* Compute weightings. */
	wgt0 = wMaxDist - wMinDist;
	wgt1 = 1.0 - minDistWgt;
	wgt2 = (wgt0 > DBL_EPSILON)? wgt1 / wgt0: 1.0;
	for(idS = 0; idS < idM; ++idS)
	{
	  if(vType == WLZ_VERTEX_D2)
	  {
	    sTV.d2 = *(sTVxBuf.d2 + idS);
	    tV.d2 = *(tVxBuf.d2 + idS);
	    WLZ_VTX_2_SUB(dV.d2, tV.d2, sTV.d2);
	    dist = WLZ_VTX_2_LENGTH(dV.d2);
	  }
	  else /* vType == WLZ_VERTEX_D3 */
	  {
	    sTV.d3 = *(sTVxBuf.d3 + idS);
	    tV.d3 = *(tVxBuf.d3 + idS);
	    WLZ_VTX_3_SUB(dV.d3, tV.d3, sTV.d3);
	    dist = WLZ_VTX_3_LENGTH(dV.d3);
	  }
	  if(wgt0 > DBL_EPSILON)
	  {
	    wVx = 1.0 - wgt2 * (dist - wMinDist);
	  }
	  else
	  {
	    wVx = 1.0;
	  }
	  wNr = *(wgtBuf + idS);
	  if(sgnNrm && (wNr < 0.0))
	  {
	    wNr = 0.0;
	  }
	  if(usrWgtFn)
	  {
	    if(vType == WLZ_VERTEX_D2)
	    {
	      sV.d2 = WlzAffineTransformVertexD2(invTr, sTV.d2, NULL);
	    }
	    else /* vType == WLZ_VERTEX_D3 */
	    {
	      sV.d3 = WlzAffineTransformVertexD3(invTr, sTV.d3, NULL);
	    }
	    *(wgtBuf + idS) = (*usrWgtFn)(vType, curTr, tree, tVx, sVx,
	    				  tV, sV, wVx, wNr, usrWgtData);
	  }
	  else
	  {
	    *(wgtBuf + idS) = wVx * wNr;
	  }
	  curMetric += dist * *(wgtBuf + idS);
	}
	curMetric /= idM;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        /* Test for convergence. */
	conv = (prvMetric - curMetric) < delta;
	if(conv)
	{
	  if(prvTr)
	  {
	    (void )WlzFreeAffineTransform(curTr);
	    curTr = prvTr;
	    prvTr = NULL;
	  }
	}
	else
	{
	  /* Compute the transform. */
	  newTr = WlzAffineTransformLSq(vType, idM, tVxBuf,
					idM, sTVxBuf,
					idM, wgtBuf,
					trType, &errNum);
	  if(errNum == WLZ_ERR_PARAM_DATA)
	  {
	    /* If the registration failed because of the given vertex
	     * position and normal data then set the iteration count
	     * higher than the limit. This will be set as a convergence
	     * failure later. */
	    itr = maxItr;
	  }
	  else if(errNum == WLZ_ERR_NONE)
	  {
	    if(curTr == NULL)
	    {
	      curTr = newTr;
	    }
	    else
	    {
	      (void )WlzFreeAffineTransform(prvTr);
	      prvTr = curTr;
	      curTr = WlzAffineTransformProduct(newTr, prvTr, &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      (void )WlzFreeAffineTransform(invTr);
	      invTr = WlzAffineTransformInverse(curTr, &errNum);
	    }
	  }
	  (void )WlzFreeAffineTransform(newTr);
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (itr++ < maxItr) && (conv == 0));
  }
  if(itr >= maxItr)
  {
    errNum = WLZ_ERR_ALG_CONVERGENCE;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeAffineTransform(curTr);
    curTr = NULL;
  }
  if(errNum == WLZ_ERR_ALG_SINGULAR)
  {
    errNum = WLZ_ERR_NONE;
    conv = 0;
  }
  (void )WlzFreeAffineTransform(invTr);
  (void )WlzFreeAffineTransform(prvTr);
  *gPrvMetric = prvMetric;
  *gCurMetric = curMetric;
  if(dstConv)
  {
    *dstConv = conv;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(curTr);
}
