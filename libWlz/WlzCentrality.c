#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCentrality_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCentrality.c
* \author       Bill Hill
* \date         April 2008
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
* \brief	Functions for computing  the centrality of a feature
* 		domain with respect to a boundary domain. See
* 		WlzCentrality().
* \ingroup	WlzFeatures
*/

#include <math.h>
#include <float.h>
#include <Wlz.h>


static double			WlzCentrality2D(
				  WlzObject *fObj,
				  WlzObject *bObj,
				  int nRay,
				  int binFlg,
				  double *dstMaxR,
				  WlzErrorNum *dstErr);
static void			WlzCentralityCompPolarTbl2D(
				  WlzDVertex2 cmV,
				  int pTblSz,
				  double *pTbl,
				  WlzObject *bndObj);
static void			WlzCentralityCompPolarTblPly2D(
				  WlzDVertex2 cmV,
				  int pTblSz,
				  double *pTbl,
				  WlzPolygonDomain *ply);
static void			WlzCentralityCompPolarTblBnd2D(
				  WlzDVertex2 cmV,
				  int pTblSz,
				  double *pTbl,
				  WlzBoundList  *bndLst);
static void			WlzCentralityUpdate2D(
				  double *fNum,
				  double *fDnm,
				  WlzDVertex2 cmV,
				  int pTblSz,
				  double *pTbl,
				  double mass,
				  WlzIVertex2 pos);

/*!
* \return	Centrality feature value in range [0.0-1.0].
* \ingroup	WlzFeatures
* \brief	Computes the centrality of a feature domain with respect
* 		to a boundary domain. The boundary domain must enclose
* 		the feature domain.
* 		
* 		Rays are projected from the centre of mass of the boundary
* 		domain to it's boundary. The maximum value of this distance
* 		for a given angle \f$i\f$ being \f$R_i\f$. At each point
* 		in the feature domain intersected by the ray the currrent
* 		value of the centrality is updated. The distance from the
* 		centre of mass to each feature point is defined to be
* 		\f$r_{i,j}\f$ and the mass at this point \f$m_{i,j}\f$.
* 		The centrality of the feature domain \f$\Omega_f\f$ with
* 		respect to the border domain \f$Omega_b\f$ is defined to be
* 		\f[
		c = \frac{\sum_{i,j}{m_{i,j}(R_i - r_{i,j})}}
		         {\sum_{i,j}{m_{i,j}R}}
                \f]
* \param	fObj			Feature domain object, \f$Omega_f\f$.
* \param	bObj			Boundary domain object, \f$Omega_b\f$.
* \param	nRay			Number of equally spaced rays
* 					projected from the centre of mass.
* \param	binFlg			Treat as binary object if non-zero,
* 					with all masses having value 1.0.
* \param	dstMaxR			Destination pointer for maximum
* 					boundary radius, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
double		WlzCentrality(WlzObject *fObj, WlzObject *bObj,
			      int nRay, int binFlg,
			      double *dstMaxR,
			      WlzErrorNum *dstErr)
{
  double	cent = 0.0,
  		maxR = 0.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fObj == NULL) || (bObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(fObj->type != bObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((fObj->domain.core == NULL) || (bObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((binFlg == 0) && (fObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(fObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	cent = WlzCentrality2D(fObj, bObj, nRay, binFlg, &maxR, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      default:
        errNum = WLZ_ERR_NONE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (dstMaxR != NULL))
  {
    *dstMaxR = maxR;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cent);
}

/*!
* \return	Centrality feature value in range [0.0-1.0].
* \ingroup	WlzFeatures
* \brief	Computes the centrality of a feature domain with respect
* 		to a boundary domain, where the domains are 2D. See
* 		WlzCentrality().
* \param	fObj			Feature domain object.
* \param	bObj			Boundary domain object.
* \param	nRay			Number of equally spaced rays projected
* 					from the centre of mass.
* \param	binFlg			Treat as binary object if non-zero.
* \param	dstMaxR			Destination pointer for maximum
* 					boundary radius, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static double	WlzCentrality2D(WlzObject *fObj, WlzObject *bObj,
				int nRay, int binFlg,
				double *dstMaxR, WlzErrorNum *dstErr)
{
  int		idx;
  double	cent = 0.0,
  		fNum = 0.0,
		fDnm = 0.0,
		maxR = 0.0;
  WlzIVertex2	pos;
  WlzDVertex2	cmV;
  double	*pTbl = NULL;
  WlzObject	*tmpObj,
		*bndObj = NULL;
  WlzGreyP	gPix;
  WlzGreyWSpace gWsp;
  WlzIntervalWSpace iWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Get centre of mass of the boundary object. */
  cmV = WlzCentreOfMass2D(bObj, 1, NULL, &errNum);
  /* Get boundary of the boundary domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    bndObj = WlzObjToBoundary(bObj, 0, &errNum);
  }
  /* Make sure the boundary is 8-connected. */
  if(errNum == WLZ_ERR_NONE)
  {
    tmpObj = WlzBoundTo8Bound(bndObj, &errNum);
    (void )WlzFreeObj(bndObj);
    bndObj = tmpObj;
  }
  /* Compute polar maximum boundary table. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((pTbl = (double *)AlcMalloc(nRay * sizeof(double))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCentralityCompPolarTbl2D(cmV, nRay, pTbl, bndObj);
    if(dstMaxR != NULL)
    {
      for(idx = 0; idx < nRay; ++idx)
      {
        if(pTbl[idx] > maxR)
	{
	  maxR = pTbl[idx];
	}
      }
    }
  }
  (void )WlzFreeObj(bndObj);
  /* Scan through feature domain adding to feature values. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(binFlg != 0)
    {
      errNum = WlzInitRasterScan(fObj, &iWsp, WLZ_RASTERDIR_ILIC);
      while((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE)
      {
	pos.vtY = iWsp.linpos;
	for(pos.vtX = iWsp.lftpos; pos.vtX <= iWsp.rgtpos; ++pos.vtX)
	{
	  WlzCentralityUpdate2D(&fNum, &fDnm, cmV, nRay, pTbl, 1.0, pos);
	}
      }
    }
    else
    {
      errNum = WlzInitGreyScan(fObj, &iWsp, &gWsp);
      while((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE)
      {
	gPix = gWsp.u_grintptr;
	pos.vtY = iWsp.linpos;
	for(pos.vtX = iWsp.lftpos; pos.vtX <= iWsp.rgtpos; ++pos.vtX)
	{
	  switch(gWsp.pixeltype)
	  {
	    case WLZ_GREY_INT:
	      WlzCentralityUpdate2D(&fNum, &fDnm, cmV, nRay, pTbl,
	                            *(gPix.inp)++, pos);
	      break;
	    case WLZ_GREY_SHORT:
	      WlzCentralityUpdate2D(&fNum, &fDnm, cmV, nRay, pTbl,
	                            *(gPix.shp)++, pos);
	      break;
	    case WLZ_GREY_UBYTE:
	      WlzCentralityUpdate2D(&fNum, &fDnm, cmV, nRay, pTbl,
	                            *(gPix.ubp)++, pos);
	      break;
	    case WLZ_GREY_FLOAT:
	      WlzCentralityUpdate2D(&fNum, &fDnm, cmV, nRay, pTbl,
	                            *(gPix.flp)++, pos);
	      break;
	    case WLZ_GREY_DOUBLE:
	      WlzCentralityUpdate2D(&fNum, &fDnm, cmV, nRay, pTbl,
	                            *(gPix.dbp)++, pos);
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cent = (fabs(fDnm) > DBL_EPSILON)? fNum / fDnm: DBL_MAX;
    if(dstMaxR != NULL)
    {
      *dstMaxR = maxR;
    }
  }
  AlcFree(pTbl);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cent);
}

/*!
* \ingroup	WlzFeatures
* \brief	Computes the polar table with quantized angles as the index
* 		and maximum radius from the boundary object's centre of mass
* 		to it's boundary as the value.
* \param	cmV			Position of the centre of mass.
* \param	pTblSz			Number of angle increments and
* 					table size.
* \param	pTbl			The polar table which is to have it's
* 					values computed.
* \param	bndObj			The boundary domain object.
*/
static void	WlzCentralityCompPolarTbl2D(WlzDVertex2 cmV, 
			        int pTblSz, double *pTbl, WlzObject *bndObj)
{
  int		idx;
  double	rad;

  for(idx = 0; idx < pTblSz; ++idx)
  {
    pTbl[idx] = -1.0;
  }
  WlzCentralityCompPolarTblBnd2D(cmV, pTblSz, pTbl, bndObj->domain.b);
  /* Fill in any gaps in the table. */
  idx = 0;
  rad = -1.0;
  while(idx < pTblSz)
  {
    if(pTbl[idx] > 0.0)
    {
      rad = pTbl[idx];
    }
    else
    {
      pTbl[idx] = rad;
    }
    ++idx;
  }
  idx = 0;
  while((pTbl[idx] < 0.0) && (idx < pTblSz))
  {
    pTbl[idx] = rad;
    ++idx;
  }
}

/*!
* \ingroup	WlzFeatures
* \brief	Computes values of the polar table for a single boundary
* 		polyline.
* \param	cmV			Position of the centre of mass.
* \param	pTblSz			Number of angle increments in the
* 					polar angle/distance table.
* \param	pTbl			Polar angle/distance table.
* \param	ply			Given polyline.
*/
static void	WlzCentralityCompPolarTblPly2D(WlzDVertex2 cmV, 
			        int pTblSz, double *pTbl,
				WlzPolygonDomain *ply)
{
  int		idA,
  		idV;
  double	ang,
  		rad,
		pTblSz1;
  WlzVertexP	vtxP;
  WlzDVertex2	tmpV;

  vtxP.i2 = ply->vtx;
  pTblSz1 = pTblSz - 0.1;
  switch(ply->type)
  {
    case WLZ_POLYGON_INT:
      for(idV = 0; idV < ply->nvertices; ++idV)
      {
	tmpV.vtX = vtxP.i2[idV].vtX;
	tmpV.vtY = vtxP.i2[idV].vtY;
	ang = WlzGeomPolar2D(cmV, tmpV, &rad);
	idA = (int )(floor((ang * pTblSz1) / (2.0 * ALG_M_PI)));
	if(rad > pTbl[idA])
	{
	  pTbl[idA] = rad;
	}
      }
      break;
    case WLZ_POLYGON_FLOAT:
      for(idV = 0; idV < ply->nvertices; ++idV)
      {
	tmpV.vtX = vtxP.f2[idV].vtX;
	tmpV.vtY = vtxP.f2[idV].vtY;
	ang = WlzGeomPolar2D(cmV, tmpV, &rad);
	idA = (int )(floor((ang * pTblSz1) / (2.0 * ALG_M_PI)));
	if(rad > pTbl[idA])
	{
	  pTbl[idA] = rad;
	}
      }
      break;
    case WLZ_POLYGON_DOUBLE:
      for(idV = 0; idV < ply->nvertices; ++idV)
      {
	ang = WlzGeomPolar2D(cmV, vtxP.d2[idV], &rad);
	idA = (int )(floor((ang * pTblSz1) / (2.0 * ALG_M_PI)));
	if(rad > pTbl[idA])
	{
	  pTbl[idA] = rad;
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup	WlzFeatures
* \brief	Computes values of the polar table for a boundary
* 		list. This is a recursive function.
* \param	cmV			Position of the centre of mass.
* \param	pTblSz			Number of angle increments in the
* 					polar angle/distance table.
* \param	pTbl			Polar angle/distance table.
* \param	bndLst			Given boundary list.
*/
static void	WlzCentralityCompPolarTblBnd2D(WlzDVertex2 cmV, 
			        int pTblSz, double *pTbl,
				WlzBoundList  *bndLst)
{
  if(bndLst)
  {
    WlzCentralityCompPolarTblPly2D(cmV, pTblSz, pTbl, bndLst->poly);
    WlzCentralityCompPolarTblBnd2D(cmV, pTblSz, pTbl, bndLst->next);
    WlzCentralityCompPolarTblBnd2D(cmV, pTblSz, pTbl, bndLst->down);
  }
}

/*!
* \ingroup	WlzFeatures
* \brief	Updates the numererator and denominator sums for the given
* 		position. See WlzCentrality().
* \param	fNum			Pointer to current numererator
* 					for update.
* \param	fDnm			Pointer to current denominator
* 					for update.
* \param	cmV			Position of the centre of mass.
* \param	pTblSz			Number of angle increments in the
* 					polar angle/distance table.
* \param	pTbl			Polar angle/distance table.
* \param	mass			Mass at given position.
* \param	pos			Given position.
*/
static void	WlzCentralityUpdate2D(double *fNum, double *fDnm,
				WlzDVertex2 cmV,
				int pTblSz, double *pTbl,
				double mass, WlzIVertex2 pos)
{
  int		idA;
  double	ang,
  		rad;
  WlzDVertex2	posD;

  posD.vtX = pos.vtX;
  posD.vtY = pos.vtY;
  ang = WlzGeomPolar2D(cmV, posD, &rad);
  idA = (int )(floor((ang * (pTblSz - 1.0)) / (2.0 * ALG_M_PI)));
  *fNum += mass * (pTbl[idA] - rad);
  *fDnm += mass * pTbl[idA];
}

