#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBndBasisFn_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzBndBasisFn.c
* \author       Guangjie Feng
* \date         August 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief   Main (top-level) Woolz binding header file which includes
*      all other header files required by the Woolz binding
*      library.
* \todo         -
* \bug          None known.
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

