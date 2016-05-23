#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFitPlane_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzFitPlane.c
* \author       Bill Hill
* \date         May 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief 	Functions for computing best fit planes.
* \ingroup	WlzTransform
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Fits a plane to the given vertices using singular value
* 		decomposition. The best fit plane is returned (provided
* 		there is no error) as a unit normal \f$n\f$ and the
* 		centroid of the given vertices ( \f$r_0\f$ ) which in a
* 		point in the plane, with the plane equation
* 		\f[
		\vec{n}\cdot(\vec{r} - \vec{r_0}) = \vec{0}
 		\f]
* \param	vtxType			Given vertex type.
* \param	nVtx			Number of given vertices.
* \param	vtx			Given vertices.
* \param	dstPinP			Destination pointer for the return
* 					of the median point (\f$r_0\f$) in
* 					the plane.
* \param	dstNrm			Destination pointer for the unit
*   					normal to the plane (\f$n\f$).
*/
WlzErrorNum WlzFitPlaneSVD(WlzVertexType vtxType, int nVtx, WlzVertexP vtx,
			   WlzDVertex3 *dstPinP, WlzDVertex3 *dstNrm)
{
  AlgMatrix     aMat,
  		vMat;
  WlzDVertex3	centroid,
                normal;
  double	*sVec = NULL;
  AlgError	algErr = ALG_ERR_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  aMat.core = NULL;
  vMat.core = NULL;
  if(nVtx < 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(vtx.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((sVec = (double *)AlcCalloc(sizeof(double), 3)) == NULL) ||
       ((aMat.rect = AlgMatrixRectNew(nVtx, 3, NULL)) == NULL) ||
       ((vMat.rect = AlgMatrixRectNew(3, 3, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Compute a Nx3 matrix from the given vertices then subtract their
   * centroid. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    double	*row;
    
    WLZ_VTX_3_ZERO(centroid);
    switch(vtxType)
    {
      case WLZ_VERTEX_I3:
	for(i = 0; i < nVtx; ++i)
	{
	  WlzIVertex3 *v;

	  v = vtx.i3 + i;
          row = aMat.rect->array[i];
	  centroid.vtX += row[0] = v->vtX;
	  centroid.vtY += row[1] = v->vtY;
	  centroid.vtZ += row[2] = v->vtZ;
	}
	break;
      case WLZ_VERTEX_F3:
	for(i = 0; i < nVtx; ++i)
	{
	  WlzFVertex3 *v;

	  v = vtx.f3 + i;
          row = aMat.rect->array[i];
	  centroid.vtX += row[0] = v->vtX;
	  centroid.vtY += row[1] = v->vtY;
	  centroid.vtZ += row[2] = v->vtZ;
	}
	break;
      case WLZ_VERTEX_D3:
	for(i = 0; i < nVtx; ++i)
	{
	  WlzDVertex3 *v;

	  v = vtx.d3 + i;
          row = aMat.rect->array[i];
	  centroid.vtX += row[0] = v->vtX;
	  centroid.vtY += row[1] = v->vtY;
	  centroid.vtZ += row[2] = v->vtZ;
	}
	break;
      default:
	errNum = WLZ_ERR_PARAM_TYPE;
        break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WLZ_VTX_3_SCALE(centroid, centroid, (1.0 / nVtx));
      for(i = 0; i < nVtx; ++i)
      {
	row = aMat.rect->array[i];
	row[0] -= centroid.vtX;
	row[1] -= centroid.vtY;
	row[2] -= centroid.vtZ;
      }
    }
  }
  /* Compute SVD. */
  if(errNum == WLZ_ERR_NONE)
  {
    algErr = AlgMatrixSVDecomp(aMat, sVec, vMat);
    errNum = WlzErrorFromAlg(algErr);
  }
  /* Find the vector corresponding to the least singular value. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i,
    		m;
    double	len;
    double	**ary;
    const double eps = 1.0e-06;

    m  = 0;
    for(i = 1; i < 3; ++i)
    {
      if(sVec[i] < sVec[m])
      {
        m = i;
      }
    }
    ary = vMat.rect->array;
#ifdef WLZ_FITPLANESVD_DEBUG
    (void )fprintf(stderr,
                   "WlzFitPlaneSVD() singular values = {%lg, %lg, %lg}\n",
		   sVec[0], sVec[1], sVec[2]);
    (void )fprintf(stderr,
                   "WlzFitPlaneSVD() index of min singular value = %d\n",
		   m);
    (void )fprintf(stderr,
                   "WlzFitPlaneSVD() singular vectors = \n"
		   "{\n");
    for(i = 0; i < 3; ++i)
    {
      (void )fprintf(stderr,
                     " {%lg, %lg, %lg}\n",
		     ary[i][0], ary[i][1], ary[i][2]);
    }
    (void )fprintf(stderr,
		   "}\n");
#endif /* WLZ_FITPLANESVD_DEBUG */
    normal.vtX = ary[0][m];
    normal.vtY = ary[1][m];
    normal.vtZ = ary[2][m];
    len = WLZ_VTX_3_LENGTH(normal);
    if(len < eps)
    {
      errNum = WLZ_ERR_ALG_SINGULAR;
    }
    else
    {
      WLZ_VTX_3_SCALE(normal, normal, (1.0 / len));
    }
  }
  AlgMatrixRectFree(aMat.rect);
  AlgMatrixRectFree(vMat.rect);
  AlcFree(sVec);
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstNrm)
    {
      *dstNrm = normal;
    }
    if(dstPinP)
    {
      *dstPinP = centroid;
    }
  }
  return(errNum);
}
