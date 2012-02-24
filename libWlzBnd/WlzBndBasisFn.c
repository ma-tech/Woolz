#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBndBasisFn_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzBnd/WlzBndBasisFn.c
* \author       Guangjie Feng
* \date         August 2003
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
* \brief   	Binding for basis functions.
* \ingroup	LibWlzBnd
*/
#include <WlzBnd.h>


 WlzBasisFn *WlzBndBasisFnMQ3DFromCPts(int arraySizeVec0,
     WlzDVertex3 *arrayVec0,
               int arraySizeVec1,
               WlzDVertex3 *arrayVec1,
               double delta,
               WlzErrorNum *dstErr)
   {
      return WlzBasisFnMQ3DFromCPts(arraySizeVec0,
                  arrayVec0,
                  arrayVec1,
                  delta,
		  NULL, NULL,
                  dstErr);

   }

   WlzDVertex3   WlzBndBasisFnTransformVertexD3(WlzBasisFnTransform *basisTr,
                  WlzDVertex3 srcVx,
                  WlzErrorNum *dstErr)
   {
   WlzDVertex3   dstVx;
   WlzErrorNum   errNum = WLZ_ERR_NONE;

   dstVx.vtX = 0;
   dstVx.vtY = 0;
   dstVx.vtZ = 0;

   if(basisTr == NULL)
   {
      errNum = WLZ_ERR_OBJECT_NULL;
   }
   else
   {
      switch(basisTr->basisFn->type)
      {
      case WLZ_FN_BASIS_3DGAUSS:
      /*dstVx = WlzBasisFnValueGauss2D(basisTr->basisFn, srcVx);
      dstVx.vtX += srcVx.vtX;
      dstVx.vtY += srcVx.vtY;*/
      break;
      case WLZ_FN_BASIS_3DPOLY:
      /*dstVx = WlzBasisFnValuePoly2D(basisTr->basisFn, srcVx);
      dstVx.vtX += srcVx.vtX;
      dstVx.vtY += srcVx.vtY;*/
      break;
      case WLZ_FN_BASIS_3DMQ:
      dstVx = WlzBasisFnValueMQ3D(basisTr->basisFn, srcVx);
      dstVx.vtX += srcVx.vtX;
      dstVx.vtY += srcVx.vtY;
      dstVx.vtZ += srcVx.vtZ;
      break;
      case WLZ_FN_BASIS_3DTPS:
      /*dstVx = WlzBasisFnValueTPS2D(basisTr->basisFn, srcVx);
      dstVx.vtX += srcVx.vtX;
      dstVx.vtY += srcVx.vtY;*/
      break;
      case WLZ_FN_BASIS_3DCONF_POLY:
      /*dstVx = WlzBasisFnValueConf2D(basisTr->basisFn, srcVx);
      dstVx.vtX += srcVx.vtX;
      dstVx.vtY += srcVx.vtY;*/
      break;
      default:
      errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
      }
   }
   if(dstErr)
   {
      *dstErr = errNum;
   }
   return(dstVx);
   }


   WlzBasisFnTransform *WlzBndBasisFnTrFromCPts3(WlzFnType type, int order,
                int arraySizeVec0, WlzDVertex3 *arrayVec0, int arraySizeVec1, WlzDVertex3 *arrayVec1, WlzErrorNum *dstErr)
    {
        return  WlzBasisFnTrFromCPts3D(type, order, arraySizeVec0, arrayVec0, arraySizeVec1, arrayVec1, NULL, dstErr);
    }

