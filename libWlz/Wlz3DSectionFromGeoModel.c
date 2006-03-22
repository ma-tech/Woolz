#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _Wlz3DSectionFromGeoModel_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/Wlz3DSectionFromGeoModel.c
* \author       Bill Hill
* \date         September 2001
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
* \brief	Functions to cut a 2D geometric model from a 3D
* 		geometric model.
* \ingroup	WlzSectionTransform
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

/*!
* \return				New 2D geometric model or NULL
*					on error.
* \ingroup      WlzSectionTransform
* \brief	Get a 2D geometric model which is the intersection a
*		plane, specified by the 3D view structure, and the given
*		3D geometric model.
* \param	gModel			Given 3D geometric model.
* \param	view			Given view structure which specifies
* 					the plane.
* \param	dstErr			Destination pointer for an error
*					code, may be NULL.
*/
WlzGMModel	*WlzGetSectionFromGMModel(WlzGMModel *gModel,
				WlzThreeDViewStruct *view, WlzErrorNum *dstErr)
{
  int		iFlg;
  WlzGMEdgeT	*cET;
  WlzGMLoopT	*cLT,
  		*fLT;
  WlzGMShell	*cS,
  		*fS;
  WlzGMModel	*nModel = NULL;
  WlzDBox3	sG;
  double	pln[4];
  WlzDVertex2	eG2D[2];
  WlzDVertex3	eG3D[2],
  		fG[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Create a new 2D model. */
  nModel = WlzGMModelNew(WLZ_GMMOD_2D, 0, 0, &errNum);
  /* For each shell in the model look to see if the shell's bounding box
   * intersects the plane of the view structure. */
  if(errNum == WLZ_ERR_NONE)
  {
    cS = fS = gModel->child;
    do
    {
      errNum = WlzGMShellGetGBB3D(cS, &sG);
      /* Check for intersection. */
      if(errNum == WLZ_ERR_NONE)
      {
        iFlg = Wlz3DViewIntersectAABB(view, sG);
      }
      if(iFlg && (errNum == WLZ_ERR_NONE))
      {
	Wlz3DViewGetPlaneEqn(view, pln + 0, pln + 1, pln + 2, pln + 3);
        /* For each loop topology element of the shell. Avoiding double
	 * counting by only considering the first loop topology element
	 * of each face. */
        cLT = fLT = cS->child;
	do
	{
	  if(cLT == cLT->face->loopT)
	  {
            /* Get the geometry of the faces verticies. */
	    cET = cLT->edgeT;
	    if(cET->prev != cET->next->next)
	    {
	      errNum = WLZ_ERR_DOMAIN_DATA;
	    }
	    else
	    {
	      (void )WlzGMVertexGetG3D(cET->vertexT->diskT->vertex,
	      			       fG + 0);
	      (void )WlzGMVertexGetG3D(cET->next->vertexT->diskT->vertex,
	      			       fG + 1);
	      (void )WlzGMVertexGetG3D(cET->prev->vertexT->diskT->vertex,
	      			       fG + 2);
              /* Compute the intersection of the face with the plane. */
	      iFlg = WlzGeomPlaneTriangleIntersect(pln[0], pln[1],
	      				pln[2], pln[3], fG[0], fG[1], fG[2],
	      				eG3D + 0, eG3D + 1);
              /* If there is an intersection then add the edge segment to
	       * the new model. */
	      if(iFlg == 2)
	      {
		eG2D[0].vtX = eG3D[0].vtX; eG2D[0].vtY = eG3D[0].vtY;
		eG2D[1].vtX = eG3D[1].vtX; eG2D[1].vtY = eG3D[1].vtY;
	        errNum = WlzGMModelConstructSimplex2D(nModel, eG2D);
	      }
	    }
	  }
	  cLT = cLT->next;
	} while((errNum == WLZ_ERR_NONE) && (cLT != fLT));
      }
      cS = cS->next;
    } while((errNum == WLZ_ERR_NONE) && (cS != fS));
  }
  if((errNum != WLZ_ERR_NONE) && nModel)
  {
    (void )WlzGMModelFree(nModel);
    nModel = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nModel);
}
