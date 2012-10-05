#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGeoModelCellGridWSp_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGeoModelCellGridWSp.c
* \author       Bill Hill
* \date         January 2010
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
* \brief	Functions to create and free geometric models cell grid
* 		work spaces.
* \ingroup	WlzGeoModel
*/

#include <Wlz.h>
#include <limits.h>
#include <float.h>

/*!
* \return	New cell grid for the model.
* \ingroup	WlzGeoModel
* \brief	Makes a new cell grid for the given elelemnt type in the
* 		given 3D model. See WlzGeoModelGridWSpSet3D().
* \param	model			Given model.
* \param	elemType		Geometric model element type for the
* 					grid of cells.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGMGridWSp3D 	*WlzGeoModelGridWSpNew3D(WlzGMModel *model,
                                         WlzGMElemType elemType,
                                         WlzErrorNum *dstErr)
{
  WlzGMGridWSp3D *grid = NULL;
  WlzGMResource	*res;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(model->type)
    {
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
        break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    res = WlzGMModelGetRes(model, elemType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(res->numElm < 1)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Allocate the cell grid data structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((grid = (WlzGMGridWSp3D *)
               AlcCalloc(1, sizeof(WlzGMGridWSp3D))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    grid->elemType = elemType;
  }
  /* Populate the cells of the cell grid data structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGeoModelGridWSpSet3D(grid, model, elemType);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGeoModelGridFree3D(grid);
    grid = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(grid);
}

/*!
* \return	New cell grid for the model.
* \ingroup	WlzGeoModel
* \brief	Populates a 3D geometric model grid workspace from the given
* 		3D model.
* 		Currently this function has only been implimented for when
* 		geometric models elements are faces all other types will
* 		result in and unimplimented error code.
* \param	grid			Given grid to populate using model.
* \param	model			Given model.
* \param	elemType		Geometric model element type for the
* 					grid of cells.
*/
WlzErrorNum 	WlzGeoModelGridWSpSet3D(WlzGMGridWSp3D *grid,
                                        WlzGMModel *model,
                                        WlzGMElemType elemType)
{
  WlzDBox3	mBox;
  WlzGMResource	*res = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(grid == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(elemType)
    {
      case WLZ_GMELM_FACE:
        break;
      default:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    res = WlzGMModelGetRes(model, elemType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(res->numElm < 1)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  else
  {
    switch(model->type)
    {
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
        break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mBox = WlzBoundingBoxGModel3D(model, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzIVertex3	newNCells;
    WlzDVertex3	mSz;

    mSz.vtX = mBox.xMax - mBox.xMin;
    mSz.vtY = mBox.yMax - mBox.yMin;
    mSz.vtZ = mBox.zMax - mBox.zMin;
    grid->org.vtX = mBox.xMin;
    grid->org.vtY = mBox.yMin;
    grid->org.vtZ = mBox.zMin;
    grid->cellSz = ALG_MAX3(mSz.vtX, mSz.vtY, mSz.vtZ);
    grid->cellSz = grid->cellSz / cbrt(fabs(res->numElm) + 1.0);
    newNCells.vtX = (int )ceil(mSz.vtX / grid->cellSz) + 1;
    newNCells.vtY = (int )ceil(mSz.vtY / grid->cellSz) + 1;
    newNCells.vtZ = (int )ceil(mSz.vtZ / grid->cellSz) + 1;
    if((grid->cells == NULL) ||
       (grid->nCells.vtX != newNCells.vtX) ||
       (grid->nCells.vtY != newNCells.vtY) ||
       (grid->nCells.vtZ != newNCells.vtZ))
    {
      grid->nCells.vtX = newNCells.vtX;
      grid->nCells.vtY = newNCells.vtY;
      grid->nCells.vtZ = newNCells.vtZ;
      (void )Alc3Free((void ***)(grid->cells));
      if(AlcPtr3Malloc((void *****)&(grid->cells), grid->nCells.vtZ,
                       grid->nCells.vtY, grid->nCells.vtX) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    grid->cellVecMax = 0;
    if(grid->cellVec == NULL)
    {
      size_t cellVecBlkSz;

      cellVecBlkSz = res->numElm / 16;
      if(cellVecBlkSz < 1024)
      {
	cellVecBlkSz = 1024;
      }
      if((grid->cellVec = AlcVectorNew(cellVecBlkSz, sizeof(WlzGMGridWSpCell3D),
				       cellVecBlkSz, NULL)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  /* Populate the cells of the cell grid data structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    int       idE,
    	      idX,
	      idY,
	      idZ;
    AlcVector *eVec;

    /* First clear the cell lists. */
    for(idZ = 0; idZ < grid->nCells.vtZ; ++idZ)
    {
      for(idY = 0; idY < grid->nCells.vtY; ++idY)
      {
	WlzGMGridWSpCell3D **cellP;

	cellP = *(*(grid->cells + idZ) + idY);
	for(idX = 0; idX < grid->nCells.vtX; ++idX)
	{
	  cellP[idX] = NULL;
	}
      }
    }
    eVec = res->vec;
    /* For each element of the model. */
    for(idE = 0; idE < res->numIdx; ++idE)
    {
      WlzGMElemP elemP;

      elemP.core = (WlzGMCore *)AlcVectorItemGet(eVec, idE);
      if(elemP.core->idx >= 0)
      {
	WlzDVertex3 fVtx[3];
	WlzIBox3 cBox;
	WlzDBox3 fBox;

	/* Get the element vertices and bounding box. */
	(void )WlzGMFaceGetG3D(elemP.face, fVtx + 0, fVtx + 1, fVtx + 2);
	fBox = WlzBoundingBoxVtx3D(3, fVtx, NULL);
	/* Compute range of grid cells which intersect the face's bounding
	 * box. */
	cBox = WlzGeoModelGridCellsInDBox(grid, fBox);
	/* For each cell that the bonding box intersects. */
	for(idZ = cBox.zMin; idZ <= cBox.zMax; ++idZ)
	{
	  WlzDVertex3 cellMin,
	  	      cellMax;

	  cellMin.vtZ = grid->org.vtZ + grid->cellSz * idZ;
	  cellMax.vtZ = grid->org.vtZ + grid->cellSz * (idZ + 1);
	  for(idY = cBox.yMin; idY <= cBox.yMax; ++idY)
	  {
	    cellMin.vtY = grid->org.vtY + grid->cellSz * idY;
	    cellMax.vtY = grid->org.vtY + grid->cellSz * (idY + 1);
	    for(idX = cBox.xMin; idX <= cBox.xMax; ++idX)
	    {
	      cellMin.vtX = grid->org.vtX + grid->cellSz * idX;
	      cellMax.vtX = grid->org.vtX + grid->cellSz * (idX + 1);
	      /* Check for an intersection between the face and the cell. */
	      if(WlzGeomTriangleAABBIntersect3D(fVtx[0], fVtx[1], fVtx[2],
	                                        cellMin, cellMax, 0) != 0)
	      {
		/* Add the face to the cells list. */
		WlzGMGridWSpCell3D *cell;
		WlzGMGridWSpCell3D **headCell;

	        cell = (WlzGMGridWSpCell3D *)
		       AlcVectorExtendAndGet(grid->cellVec, grid->cellVecMax);
                if(cell == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		  goto RETURN;
		}
		else
		{
		  ++(grid->cellVecMax);
		  cell->elem = elemP;
		  headCell = *(*(grid->cells + idZ) + idY) + idZ;
		  cell->next = *headCell;
		  *headCell = cell;
		}
	      }
	    }
	  }
	}
      }
    }
  }
RETURN:
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeoModel
* \brief	Frees the given geometric model cell grid.
* \param	grid			Given grid to free.
*/
WlzErrorNum	WlzGeoModelGridFree3D(WlzGMGridWSp3D *grid)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(grid != NULL)
  {
    (void )Alc3Free((void ***)(grid->cells));
    (void )AlcVectorFree(grid->cellVec);
  }
  else
  {
    (void )AlcFree(grid);
  }
  return(errNum);
}

/*!
* \return	Box of grid cells which enclose the given box.
* \ingroup	WlzGeoModel
* \brief	Computes the box of grid cells which encloses the given
* 		double box.
* \param	grid			Given grid of cells.
* \param	dBox			Given double box.
*/
WlzIBox3	WlzGeoModelGridCellsInDBox(WlzGMGridWSp3D *grid,
				           WlzDBox3 dBox)
{
  WlzIBox3	cBox;
  const double	tol = ALG_DBL_TOLLERANCE;

  if((cBox.xMin = (int )floor((dBox.xMin - grid->org.vtX) /
                              grid->cellSz - tol)) < 0)
  {
    cBox.xMin = 0;
  }
  else if(cBox.xMin >= grid->nCells.vtX)
  {
    cBox.xMin = grid->nCells.vtX - 1;
  }
  if((cBox.xMax = (int )ceil((dBox.xMax - grid->org.vtX) /
                             grid->cellSz + tol)) < 0)
  {
    cBox.xMax = 0;
  }
  else if(cBox.xMax >= grid->nCells.vtX)
  {
    cBox.xMax = grid->nCells.vtX - 1;
  }
  if((cBox.yMin = (int )floor((dBox.yMin - grid->org.vtY) /
                              grid->cellSz - tol)) < 0)
  {
    cBox.yMin = 0;
  }
  else if(cBox.yMin >= grid->nCells.vtY)
  {
    cBox.yMin = grid->nCells.vtY - 1;
  }
  if((cBox.yMax = (int )ceil((dBox.yMax - grid->org.vtY) /
                              grid->cellSz + tol)) < 0)
  {
    cBox.yMax = grid->nCells.vtY - 1;
  }
  else if(cBox.yMax >= grid->nCells.vtY)
  {
    cBox.yMax = grid->nCells.vtY - 1;
  }
  if((cBox.zMin = (int )floor((dBox.zMin - grid->org.vtZ) /
                              grid->cellSz - tol)) < 0)
  {
    cBox.zMin = 0;
  }
  else if(cBox.zMin >= grid->nCells.vtZ)
  {
    cBox.zMin = grid->nCells.vtZ - 1;
  }
  if((cBox.zMax = (int )ceil((dBox.zMax - grid->org.vtZ) /
                             grid->cellSz + tol)) < 0)
  {
    cBox.zMax = grid->nCells.vtZ - 1;
  }
  else if(cBox.zMax >= grid->nCells.vtZ)
  {
    cBox.zMax = grid->nCells.vtZ - 1;
  }
  return(cBox);
}
