#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDistMetric_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzDistMetric.c
* \author       Bill Hill
* \date         August 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Functions to compute the Hausdorff distance, mean nearest
* 		neighbour and the median nearest neighbour distances
* 		between two datasets.
* \ingroup	WlzFeatures
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>
#include <float.h>
#include <limits.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes any combination of the Hausdorff, mean
*		nearest neighbour, median nearest neighbour and
*		minimum nearest neighbour distances
*		between the vertices of the given geometric models.
*		See WlzDistMetricVertex2D() for details of the metrics.
* \param	model0			First geometric model.
* \param	model1			Second geometric model.
* \param	dstDistH		Destination pointer for the directed
*					Hausdorff distance, may be NULL.
* \param	dstDistM		Destination pointer for the directed
*					mean nearest neighbour distance, may
*					be NULL.
* \param	dstDistN		Destination pointer for the directed
*					median nearest neighbour distance, may
*					be NULL.
* \param	dstDistI		Destination pointer for the minimum
*					nearest neighbour distance, may
*					be NULL.
*/
WlzErrorNum 	WlzDistMetricGM(WlzGMModel *model0, WlzGMModel *model1,
			        double *dstDistH, double *dstDistM,
				double *dstDistN, double *dstDistI)
{
  int		nV[2];
  WlzVertexP	vP[2];
  WlzVertexType vType[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model0 == NULL) || (model1 == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(model0->type != model1->type)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    vP[0] = WlzVerticesFromGM(model0, NULL, NULL, nV + 0, vType + 0, &errNum); 
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vP[1] = WlzVerticesFromGM(model1, NULL, NULL, nV + 1, vType + 1, &errNum); 
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(vType[0] != vType[1])
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(vType[0])
    {
      case WLZ_VERTEX_D2:
	errNum = WlzDistMetricVertex2D(nV[0], vP[0].d2, nV[1], vP[1].d2,
				       dstDistH, dstDistM, dstDistN, dstDistI);
	break;
      case WLZ_VERTEX_D3:
	errNum = WlzDistMetricVertex3D(nV[0], vP[0].d3, nV[1], vP[1].d3,
				       dstDistH, dstDistM, dstDistN, dstDistI);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_DATA;
	break;
    }
  }
  AlcFree(vP[0].v);
  AlcFree(vP[1].v);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes any combination of the directed Hausdorff, mean
*		nearest neighbour, median nearest neighbour and minimum
*		nearest neighbour distances
*		between the vertices of the given geometric models.
*		See WlzDistMetricDirVertex2D() for details of the metrics.
* \param	model0			First geometric model.
* \param	model1			Second geometric model.
* \param	dstDistH		Destination pointer for the directed
*					Hausdorff distance, may be NULL.
* \param	dstDistM		Destination pointer for the directed
*					mean nearest neighbour distance, may
*					be NULL.
* \param	dstDistN		Destination pointer for the directed
*					median nearest neighbour distance, may
*					be NULL.
* \param	dstDistI		Destination pointer for the minimum
*					nearest neighbour distance, may
*					be NULL.
*/
WlzErrorNum 	WlzDistMetricDirGM(WlzGMModel *model0, WlzGMModel *model1,
			        double *dstDistH, double *dstDistM,
				double *dstDistN, double *dstDistI)
{
  int		nV[2];
  WlzVertexP	vP[2];
  WlzVertexType vType[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((model0 == NULL) || (model1 == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(model0->type != model1->type)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    vP[0] = WlzVerticesFromGM(model0, NULL, NULL, nV + 0, vType + 0, &errNum); 
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vP[1] = WlzVerticesFromGM(model1, NULL, NULL, nV + 1, vType + 1, &errNum); 
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(vType[0] != vType[1])
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(vType[0])
    {
      case WLZ_VERTEX_D2:
	errNum = WlzDistMetricDirVertex2D(nV[0], vP[0].d2, nV[1], vP[1].d2,
				          dstDistH, dstDistM, dstDistN,
					  dstDistI);
	break;
      case WLZ_VERTEX_D3:
	errNum = WlzDistMetricDirVertex3D(nV[0], vP[0].d3, nV[1], vP[1].d3,
				          dstDistH, dstDistM, dstDistN,
					  dstDistI);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_DATA;
	break;
    }
  }
  AlcFree(vP[0].v);
  AlcFree(vP[1].v);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes any combination of the Hausdorff, mean
*		nearest neighbour, median nearest neighbour and
*		minimum nearest neighbour distances
*		between the given sets of vertices.
*		Each of these distance measures is the maximum of the two
*		possible directed measures:
*		\f[
		D = \max{(d(A,B), d(B,A))}
		\f]
*		Where \f$D\f$ is a non-directed distance metric and \f$d\f$
*		is the associated directed distance metric. \f$A\f$ and
*		\f$B\f$ are the two datasets.
* \param	n0			Number of vertices in the first
*					array.
* \param	vx0			First array of vertices.
* \param	n1			Number of vertices in the second
*					array.
* \param	vx1			Second array of vertices.
* \param	dstDistH		Destination pointer for the directed
*					Hausdorff distance, may be NULL.
* \param	dstDistM		Destination pointer for the directed
*					mean nearest neighbour distance, may
*					be NULL.
* \param	dstDistN		Destination pointer for the directed
*					median nearest neighbour distance, may
*					be NULL.
* \param	dstDistI		Destination pointer for the minimum
*					nearest neighbour distance, may
*					be NULL.
*/
WlzErrorNum	WlzDistMetricVertex2D(int n0, WlzDVertex2 *vx0,
				      int n1, WlzDVertex2 *vx1,
				      double *dstDistH, double *dstDistM,
				      double *dstDistN, double *dstDistI)
{
  double 	distH[2],
  		distM[2],
  		distN[2],
		distI[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzDistMetricDirVertex2D(n0, vx0, n1, vx1,
  				    (dstDistH)? distH + 0: NULL,
				    (dstDistM)? distM + 0: NULL,
				    (dstDistN)? distN + 0: NULL,
				    (dstDistI)? distI + 0: NULL);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzDistMetricDirVertex2D(n1, vx1, n0, vx0,
				      (dstDistH)? distH + 1: NULL,
				      (dstDistM)? distM + 1: NULL,
				      (dstDistN)? distN + 1: NULL,
				      NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstDistH)
    {
      *dstDistH  = WLZ_MAX(distH[0], distH[1]);
    }
    if(dstDistM)
    {
      *dstDistM  = WLZ_MAX(distM[0], distM[1]);
    }
    if(dstDistN)
    {
      *dstDistN  = WLZ_MAX(distN[0], distN[1]);
    }
    if(dstDistI)
    {
      *dstDistI  = distI[0];
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes any combination of the Hausdorff, mean
*		nearest neighbour, median nearest neighbour and
*		minimum nearest neighbour distances
*		between the given sets of vertices.
*		See WlzDistMetricVertex2D() for an explaination of the
*		distance metrics.
* \param	n0			Number of vertices in the first
*					array.
* \param	vx0			First array of vertices.
* \param	n1			Number of vertices in the second
*					array.
* \param	vx1			Second array of vertices.
* \param	dstDistH		Destination pointer for the directed
*					Hausdorff distance, may be NULL.
* \param	dstDistM		Destination pointer for the directed
*					mean nearest neighbour distance, may
*					be NULL.
* \param	dstDistN		Destination pointer for the directed
*					median nearest neighbour distance, may
*					be NULL.
* \param	dstDistI		Destination pointer for the minimum
*					nearest neighbour distance, may
*					be NULL.
*/
WlzErrorNum	WlzDistMetricVertex3D(int n0, WlzDVertex3 *vx0,
				      int n1, WlzDVertex3 *vx1,
				      double *dstDistH, double *dstDistM,
				      double *dstDistN, double *dstDistI)
{
  double 	distH[2],
  		distM[2],
		distN[2],
		distI[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzDistMetricDirVertex3D(n0, vx0, n1, vx1,
  				    (dstDistH)? distH + 0: NULL,
				    (dstDistM)? distM + 0: NULL,
				    (dstDistN)? distN + 0: NULL,
				    (dstDistI)? distI + 0: NULL);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzDistMetricDirVertex3D(n1, vx1, n0, vx0,
				      (dstDistH)? distH + 1: NULL,
				      (dstDistM)? distM + 1: NULL,
				      (dstDistN)? distN + 1: NULL,
				      NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstDistH)
    {
      *dstDistH  = WLZ_MAX(distH[0], distH[1]);
    }
    if(dstDistM)
    {
      *dstDistM  = WLZ_MAX(distM[0], distM[1]);
    }
    if(dstDistN)
    {
      *dstDistN  = WLZ_MAX(distN[0], distN[1]);
    }
    if(dstDistI)
    {
      *dstDistI  = distI[0];
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes any combination of the directed Hausdorff, mean
*		nearest neighbour, median nearest neighbour
*		and minimum mean neighbour distances
*		between the given sets of vertices.
*		The directed Hausdorff distance metric is:
*		\f[
                \max_{a \in A}{\min_{b \in B}{\|a - b\|}}
		\f]
*		The directed mean nearest neighbour distance metric is:
*               \f[
                \frac{1}{n_a} \sum_{a \in A}{\min_{b \in B}{\|a - b\|}}
                \f]
*		Likewise the directed median nearest neighbour distance metric
*		is:
*		\f[
                \mathrm{median}_{a \in A}{\min_{b \in B}{\|a - b\|}}
		\f]
*		Where \f$A\f$ and \f$B\f$ are the two datasets.
* \param	n0			Number of vertices in the first
*					array.
* \param	vx0			First array of vertices.
* \param	n1			Number of vertices in the second
*					array.
* \param	vx1			Second array of vertices.
* \param	dstDistH		Destination pointer for the directed
*					Hausdorff distance, may be NULL.
* \param	dstDistM		Destination pointer for the directed
*					mean nearest neighbour distance, may
*					be NULL.
* \param	dstDistN		Destination pointer for the directed
*					median nearest neighbour distance, may
*					be NULL.
* \param	dstDistI		Destination pointer for the minimum
*					nearest neighbour distance, may
*					be NULL.
*/
WlzErrorNum 	WlzDistMetricDirVertex2D(int n0, WlzDVertex2 *vx0,
				         int n1, WlzDVertex2 *vx1,
				         double *dstDistH, double *dstDistM,
				         double *dstDistN, double *dstDistI)
{
  int		id0,
  		cCnt;
  double	cDist,
		mDist,
  		sDist,
		iDist;
  int		*iWSp = NULL;
  double	*nnDist = NULL;
  double	vxD2[2];
  WlzVertexP	tVP;
  AlcKDTNode    *tNode;
  AlcKDTTree    *tTree = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((vx0 == NULL) || (vx1 == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((n0 <= 0) || (n1 <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((iWSp = (int *)AlcMalloc(n1 * sizeof(int))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    if(dstDistN)
    {
      if((nnDist = (double *)AlcMalloc(n0 * sizeof(double))) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tVP.d2 = vx1;
    tTree = WlzVerticesBuildTree(WLZ_VERTEX_D2, n1, tVP, iWSp, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cCnt = 0;
    iDist = DBL_MAX;
    sDist = mDist = 0.0;
    for(id0 = 0; id0 < n0; ++id0)
    {
      vxD2[0] = (vx0 + id0)->vtX;
      vxD2[1] = (vx0 + id0)->vtY;
      tNode = AlcKDTGetNN(tTree, vxD2, DBL_MAX, &cDist, NULL);
      if(tNode)
      {
	sDist += cDist;
	if(cDist > mDist)
	{
	  mDist = cDist;
	}
	if(nnDist)
	{
	  *(nnDist + cCnt) = cDist;
	}
	if(cDist < iDist)
	{
	  iDist = cDist;
	}
        ++cCnt;
      }
    }
    if(cCnt <= 0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstDistH)
    {
      *dstDistH = mDist;
    }
    if(dstDistM)
    {
      *dstDistM = sDist / cCnt;
    }
    if(dstDistN)
    {
      id0 = cCnt / 2;
      AlgRankSelectD(nnDist, cCnt, id0);
      *dstDistN = *(nnDist + id0);
    }
    if(dstDistI)
    {
      *dstDistI = iDist;
    }
  }
  AlcFree(iWSp);
  AlcFree(nnDist);
  AlcKDTTreeFree(tTree);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Computes any combination of the directed Hausdorff, mean
*		nearest neighbour, median nearest neighbour and
*		minimum nearest neighbour distances
*		between the given sets of vertices.
*		See WlzDistMetricDirVertex2D() for an explaination of the
*		distance metrics.
* \param	n0			Number of vertices in the first
*					array.
* \param	vx0			First array of vertices.
* \param	n1			Number of vertices in the second
*					array.
* \param	vx1			Second array of vertices.
* \param	dstDistH		Destination pointer for the directed
*					Hausdorff distance, may be NULL.
* \param	dstDistM		Destination pointer for the directed
*					mean nearest neighbour distance, may
*					be NULL.
* \param	dstDistN		Destination pointer for the directed
*					median nearest neighbour distance, may
*					be NULL.
* \param	dstDistI		Destination pointer for the minimum
*					nearest neighbour distance, may
*					be NULL.
*/
WlzErrorNum	WlzDistMetricDirVertex3D(int n0, WlzDVertex3 *vx0,
				      int n1, WlzDVertex3 *vx1,
				      double *dstDistH, double *dstDistM,
				      double *dstDistN, double *dstDistI)
{
  int		id0,
  		cCnt;
  double	cDist,
		mDist,
  		sDist,
		iDist;
  int		*iWSp = NULL;
  double	*nnDist = NULL;
  double	vxD3[3];
  WlzVertexP	tVP;
  AlcKDTNode    *tNode;
  AlcKDTTree    *tTree = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((vx0 == NULL) || (vx1 == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((n0 <= 0) || (n1 <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((iWSp = (int *)AlcMalloc(n1 * sizeof(int))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    if(dstDistN)
    {
      if((nnDist = (double *)AlcMalloc(n0 * sizeof(double))) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tVP.d3 = vx1;
    tTree = WlzVerticesBuildTree(WLZ_VERTEX_D3, n1, tVP, iWSp, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cCnt = 0;
    iDist = DBL_MAX;
    sDist = mDist = 0.0;
    for(id0 = 0; id0 < n0; ++id0)
    {
      vxD3[0] = (vx0 + id0)->vtX;
      vxD3[1] = (vx0 + id0)->vtY;
      vxD3[2] = (vx0 + id0)->vtZ;
      tNode = AlcKDTGetNN(tTree, vxD3, DBL_MAX, &cDist, NULL);
      if(tNode)
      {
	sDist += cDist;
	if(cDist > mDist)
	{
	  mDist = cDist;
	}
	if(nnDist)
	{
	  *(nnDist + cCnt) = cDist;
	}
	if(cDist < iDist)
	{
	  iDist = cDist;
	}
        ++cCnt;
      }
    }
    if(cCnt <= 0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstDistH)
    {
      *dstDistH = mDist;
    }
    if(dstDistM)
    {
      *dstDistM = sDist / cCnt;
    }
    if(dstDistN)
    {
      id0 = cCnt / 2;
      AlgRankSelectD(nnDist, cCnt, id0);
      *dstDistN = *(nnDist + id0);
    }
    if(dstDistI)
    {
      *dstDistI = iDist;
    }
  }
  AlcFree(iWSp);
  AlcFree(nnDist);
  AlcKDTTreeFree(tTree);
  return(errNum);
}

#ifdef WLZ_DISTMETRIC_TEST_MAIN
int		main(int argc, char *argv[])
{
  double	distH,
  		distM,
		distN,
		distI;
  static WlzDVertex2 vx0[4],
  		vx1[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
                               
  vx0[0].vtX = 0; vx0[0].vtY = 2;
  vx0[1].vtX = 3; vx0[1].vtY = 2;
  vx0[2].vtX = 3; vx0[2].vtY = 5;
  vx0[3].vtX = 0; vx0[3].vtY = 5;
  vx1[0].vtX = 4; vx1[0].vtY = 1;
  vx1[1].vtX = 7; vx1[1].vtY = 1;
  vx1[2].vtX = 7; vx1[2].vtY = 3;
  vx1[3].vtX = 4; vx1[3].vtY = 3;
  errNum = WlzDistMetricVertex2D(4, vx0, 4, vx1,
  				 &distH, &distM, &distN, &distI);
  if(errNum == WLZ_ERR_NONE)
  {
    (void )printf("H %1.8f, M %1.8f, N %1.8f, I %1.8f\n",
    		  distH, distM, distN, distI);
    (void )printf("Results should be:\n"
                  "H 4.47213595, M 3.06138078, N 4.12310563, I 1.41421356\n");
  }
  return(errNum);
}
#endif /* WLZ_DISTMETRIC_TEST_MAIN */
