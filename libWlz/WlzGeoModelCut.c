#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGeoModelCut_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzGeoModelCut.c
* \author       Bill Hill
* \date         January 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
* \brief	Functions to cut geometric models (GM's).
* \ingroup	WlzGeoModel
*/

#include <Wlz.h>
#include <limits.h>
#include <float.h>

static WlzErrorNum 		WlzGMModelCut2D(
				  WlzGMModel *cut,
				  WlzGMModel *given,
				  WlzGMModel *knife);
static WlzErrorNum 		WlzGMModelCutDom2D(
				  WlzGMModel *cut,
				  WlzGMModel *given,
				  WlzObject *knife);
static WlzErrorNum 		WlzGMModelCut3D(
				  WlzGMModel *cut,
				  WlzGMModel *given,
				  WlzGMModel *knife);
static WlzErrorNum 		WlzGMModelCutDom3D(
				  WlzGMModel *cut,
				  WlzGMModel *given,
				  WlzObject *knife);
/*!
* \return	A new geometric model.
* \ingroup	WlzGeoModel
* \brief	Creates a new geometric model from the given model such
* 		that no elements of the given model which intersect
* 		with elements in the knife model are included.
* 		While this code seems to work for simple models errors
* 		have been seen with complex models. The code for this
* 		function and those static functions that it calls should
* 		be considered "under construction".
* \param	given			Given model.
* \param	knife			Knife model.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGMModel	*WlzGMModelCut(WlzGMModel *given, WlzGMModel *knife,
			       WlzErrorNum *dstErr)
{
  int		dim = 2,
  		nBkSz,
  		nHTSz;
  WlzGMModel	*cut = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const int	minBkSz = 1024,
  		minHTSz = 1024;

  if((given == NULL) || (knife == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(given->type != knife->type)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    switch(given->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N:
        dim = 2;
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
	dim = 3;
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make a new model. */
    if((nBkSz = given->res.vertex.numElm / 16) < minBkSz)
    {
      nBkSz = minBkSz;
    }
    if((nHTSz = given->res.vertex.numElm) < minHTSz)
    {
      nHTSz = minHTSz;
    }
    cut = WlzGMModelNew(given->type, nBkSz, nHTSz, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dim == 2)
    {
      errNum = WlzGMModelCut2D(cut, given, knife);
    }
    else /* dim == 3 */
    {
      errNum = WlzGMModelCut3D(cut, given, knife);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzGMModelFree(cut);
      cut = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cut);
}

/*!
* \return	A new geometric model.
* \ingroup	WlzGeoModel
* \brief	Creates a new geometric model from the given model such
* 		that no elements of the given model which intersect
* 		with elements in the knife domain are included.
* \param	given			Given model.
* \param	knife			Knife domain object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGMModel	*WlzGMModelCutDom(WlzGMModel *given, WlzObject *knife,
			       WlzErrorNum *dstErr)
{
  int		dim = 2,
  		nBkSz,
  		nHTSz;
  WlzGMModel	*cut = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const int	minBkSz = 1024,
  		minHTSz = 1024;

  if(knife == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((knife->type != WLZ_2D_DOMAINOBJ) &&
          (knife->type != WLZ_3D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((given == NULL) || (knife->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(given->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N:
	if(knife->type != WLZ_2D_DOMAINOBJ)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
        dim = 2;
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
	if(knife->type != WLZ_3D_DOMAINOBJ)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	dim = 3;
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make a new model. */
    if((nBkSz = given->res.vertex.numElm / 16) < minBkSz)
    {
      nBkSz = minBkSz;
    }
    if((nHTSz = given->res.vertex.numElm) < minHTSz)
    {
      nHTSz = minHTSz;
    }
    cut = WlzGMModelNew(given->type, nBkSz, nHTSz, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dim == 2)
    {
      errNum = WlzGMModelCutDom2D(cut, given, knife);
    }
    else /* dim == 3 */
    {
      errNum = WlzGMModelCutDom3D(cut, given, knife);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzGMModelFree(cut);
      cut = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cut);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeoModel
* \brief	Copies elements from the given model to the cut model on
* 		the condition that the elements do not intersects any
* 		elements of the knide model.
* 		All models are assumed to be valid 2D models.
* \param	given			Given model.
* \param	knife			Knife model.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzErrorNum WlzGMModelCut2D(WlzGMModel *cut, WlzGMModel *given,
				   WlzGMModel *knife)
{
  WlzErrorNum   errNum = WLZ_ERR_UNIMPLEMENTED;

  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeoModel
* \brief	Copies elements from the given model to the cut model on
* 		the condition that the elements do not lie within the knife
* 		domain.
* \param	given			Given model.
* \param	knife			Knife object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzErrorNum WlzGMModelCutDom2D(WlzGMModel *cut, WlzGMModel *given,
				      WlzObject *knife)
{
  WlzErrorNum   errNum = WLZ_ERR_UNIMPLEMENTED;

  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeoModel
* \brief	Copies elements from the given model to the cut model on
* 		the condition that the elements do not intersects any
* 		elements of the knide model.
* 		All models are assumed to be valid 3D models.
* \param	given			Given model.
* \param	knife			Knife model.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzErrorNum WlzGMModelCut3D(WlzGMModel *cut, WlzGMModel *given,
				   WlzGMModel *knife)
{
  WlzGMGridWSp3D *kGrid = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((cut == NULL) || (knife == NULL) || (given == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    int cutDim,
    	givenDim,
    	knifeDim;

    cutDim = WlzGMModelGetDimension(cut, NULL);
    givenDim = WlzGMModelGetDimension(given, NULL);
    knifeDim = WlzGMModelGetDimension(knife, NULL);
    if((cutDim != 3) || (givenDim != 3) || (knifeDim != 3))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  /* Create a cell grid workspace for the knife model. */
  if(errNum == WLZ_ERR_NONE)
  {
    kGrid = WlzGeoModelGridWSpNew3D(knife, WLZ_GMELM_FACE, &errNum);
  }
  /* For each face of the given model, if it does not intersect a face of
   * the knife model then add it to the new cut model. */
  if(errNum == WLZ_ERR_NONE)
  {
    int idF,
        idX,
	idY,
	idZ,
    	gMaxF;
    AlcVector *gFVec;

    gMaxF = given->res.face.numIdx;
    gFVec = given->res.face.vec;
    for(idF = 0; idF < gMaxF; ++idF)
    {
      WlzGMFace *gFace;

      gFace = (WlzGMFace *)AlcVectorItemGet(gFVec, idF);
      if(gFace->idx >= 0)
      {
        int	add = 1;
	WlzDVertex3 gVtx[3];
	WlzIBox3 cBox;
	WlzDBox3 fBox;

	(void )WlzGMFaceGetG3D(gFace, gVtx + 0, gVtx + 1, gVtx + 2);
	fBox = WlzBoundingBoxVtx3D(3, gVtx, NULL);
	cBox = WlzGeoModelGridCellsInDBox(kGrid, fBox);
	/* For each cell that the bonding box intersects. */
        for(idZ = cBox.zMin; idZ <= cBox.zMax; ++idZ)
        {
          WlzDVertex3 kCellMin,
                      kCellMax;

          kCellMin.vtZ = kGrid->org.vtZ + kGrid->cellSz * idZ;
          kCellMax.vtZ = kGrid->org.vtZ + kGrid->cellSz * (idZ + 1);
          for(idY = cBox.yMin; idY <= cBox.yMax; ++idY)
          {
            kCellMin.vtY = kGrid->org.vtY + kGrid->cellSz * idY;
            kCellMax.vtY = kGrid->org.vtY + kGrid->cellSz * (idY + 1);
            for(idX = cBox.xMin; idX <= cBox.xMax; ++idX)
            {
	      WlzDVertex3 kVtx[3];
	      WlzGMGridWSpCell3D *kCell;

              kCellMin.vtX = kGrid->org.vtX + kGrid->cellSz * idX;
              kCellMax.vtX = kGrid->org.vtX + kGrid->cellSz * (idX + 1);
              /* Test if the face of the given model intersects the faces
	       * of the knife grid cells. */
	      kCell = *(*(*(kGrid->cells + idZ) + idY) + idX);
	      while(kCell != NULL)
	      {
		WlzGMFace *kFace;

		kFace = kCell->elem.face;
		(void )WlzGMFaceGetG3D(kFace, kVtx + 0, kVtx + 1, kVtx + 2);
		if(WlzGeomTriangleTriangleIntersect3D(
		                  gVtx[0], gVtx[1], gVtx[2],
		                  kVtx[0], kVtx[1], kVtx[2]) != 0)
		{
		  add = 0;
		  idX = cBox.xMax + 1;            /* Break from these loops. */
		  idY = cBox.yMax + 1;
		  idZ = cBox.zMax + 1;
		  break;
		}
		if(errNum != WLZ_ERR_NONE)
		{
		  goto RETURN;
		}
	        kCell = kCell->next;
	      }
            }
	  }
	}
	if(add)
	{
	  /* Given model face does not intersect the knife model face
	   * so just add it to the cut model. */
	  errNum = WlzGMModelConstructSimplex3D(cut, gVtx);
	}
      }
    }
  }
RETURN:
  if(kGrid)
  {
    (void )WlzGeoModelGridFree3D(kGrid);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeoModel
* \brief	Copies elements from the given model to the cut model on
* 		the condition that the elements do not lie within the
* 		knide domain.
* 		It is assumed that the model is a valid 3D model and that
* 		the knife is a valid 3D (plane) domain object.
* \param	given			Given model.
* \param	knife			Knife object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzErrorNum WlzGMModelCutDom3D(WlzGMModel *cut, WlzGMModel *given,
				      WlzObject *knife)
{
  WlzPlaneDomain *kDom = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(knife == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((cut == NULL) || (given == NULL) ||
          ((kDom = knife->domain.p) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    int idF,
    	gMaxF;
    AlcVector *gFVec;

    /* For each face of the given model, if it does not intersect a face of
     * the knife model then add it to the new cut model. */
    gMaxF = given->res.face.numIdx;
    gFVec = given->res.face.vec;
    for(idF = 0; idF < gMaxF; ++idF)
    {
      WlzGMFace *gFace;

      gFace = (WlzGMFace *)AlcVectorItemGet(gFVec, idF);
      if(gFace->idx >= 0)
      {
        int	add;
	WlzDVertex3 v[3];

	(void )WlzGMFaceGetG3D(gFace, v + 0, v + 1, v + 2);
	add = (WlzInsideDomain3D(kDom, v[0].vtZ, v[0].vtY, v[0].vtX,
	                         NULL) == 0) &&
	      (WlzInsideDomain3D(kDom, v[1].vtZ, v[1].vtY, v[1].vtX, 
				 NULL) == 0) &&
	      (WlzInsideDomain3D(kDom, v[2].vtZ, v[2].vtY, v[2].vtX, 
				 NULL) == 0);
	if(add)
	{
	  /* Given model face does not intersect the knife domain so add it
	   * to the cut model. */
	  errNum = WlzGMModelConstructSimplex3D(cut, v);
	}
      }
    }
  }
  return(errNum);
}
