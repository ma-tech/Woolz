#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzBndFunction.c
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
        return  WlzBasisFnTrFromCPts3(type, order, arraySizeVec0, arrayVec0, arraySizeVec1, arrayVec1, dstErr);
    }

	WlzErrorNum WlzSetVoxelSize(WlzObject *obj, double x, double y, double z){
		WlzErrorNum errNum = WLZ_ERR_NONE;
		
		if (obj->type == WLZ_3D_DOMAINOBJ){
			obj->domain.p->voxel_size[0] = x;
			obj->domain.p->voxel_size[1] = y;
			obj->domain.p->voxel_size[2] = z;
		} else {
			errNum = WLZ_ERR_OBJECT_TYPE;
		}

		return errNum;
	}

	
	WlzObjectType WlzGetObjectType(WlzObject *obj){
		return obj->type;
	}

	WlzErrorNum	WlzExplode(int *dstExpObjCount,
			     WlzObject ***dstExpObjVecP,
				 WlzObject *srcObj){

		WlzErrorNum errNum = WLZ_ERR_NONE;
		WlzCompoundArray *cmpObj = NULL;

		if (srcObj->type == WLZ_COMPOUND_ARR_1
			|| srcObj->type == WLZ_COMPOUND_ARR_2){
			cmpObj = (WlzCompoundArray *) srcObj;
			*dstExpObjVecP = cmpObj->o;
			*dstExpObjCount = cmpObj->n;
		} else{
			errNum = WLZ_ERR_OBJECT_NULL;
		}

		return(errNum);
	}

    const char *WlzGetPropName(WlzObject *obj){
		WlzErrorNum errNum = WLZ_ERR_NONE;
		WlzProperty prop;
		const char *name = NULL;
		char *dst = NULL;

		if (obj->plist){
			prop = WlzGetProperty(obj->plist->list, WLZ_PROPERTY_NAME, &errNum);
			if (prop.core && (errNum == WLZ_ERR_NONE)){
				switch (prop.core->type){
				case WLZ_PROPERTY_NAME:
					dst = prop.name->name;
					break;
				case WLZ_PROPERTY_GREY:
					dst = prop.greyV->name;
					break;
				}
			}
		}

		name = AlcStrDup((const) dst);

		return name;
	}

	WlzObject *WlzGetContourObj(WlzObject *inObj){
		int	flip = 1,
		nrm = 0,
		nItr = 10,
		nonMan = 0,
		unitVoxelSz = 0,
		filterGeom = 1,
		setbackVz = 0;
		double lambda = 0,
  		mu = 0,
  		ctrVal = 100,
  		ctrWth = 1.0,
		filterPB = 0.1,
		filterSB = 1.1,
		xsize = 1,
		ysize = 1,
		zsize = 1;

		WlzDomain ctrDom;
		WlzValues	dumVal;
		WlzContourMethod ctrMtd = WLZ_CONTOUR_MTD_BND;
		const double	filterDPB = 0.25,
  		filterDSB = 0.10;
        const char	*errMsgStr;
		WlzObject *outObj = NULL;
		WlzErrorNum   errNum = WLZ_ERR_NONE;

		ctrDom.core = NULL;
		dumVal.core = NULL;

		if (inObj && (inObj->type == WLZ_3D_DOMAINOBJ) &&
			(inObj->domain.core) 
			&& (inObj->domain.core->type == WLZ_2D_DOMAINOBJ)){
				xsize = inObj->domain.p->voxel_size[0];
				ysize = inObj->domain.p->voxel_size[1];
				zsize = inObj->domain.p->voxel_size[2];
				
				inObj->domain.p->voxel_size[0] = 1.0;
				inObj->domain.p->voxel_size[1] = 1.0;
				inObj->domain.p->voxel_size[2] = 1.0;
			} else {
			errNum = WLZ_ERR_OBJECT_TYPE;
			}

			if (!errNum){
				ctrDom.ctr = WlzContourObj(inObj, ctrMtd,ctrVal, ctrWth, 
					nrm, &errNum);
				if (!errNum && filterGeom){
					errNum = WlzGMFilterGeomLPParam(&lambda, &mu, &nItr, 
						filterPB, filterSB, filterDPB, filterDSB);
					if(!errNum){
						errNum = WlzGMFilterGeomLPLM(ctrDom.ctr->model,
							lambda, mu, nItr, nonMan);
						if (!errNum){
							outObj = WlzMakeMain(WLZ_CONTOUR, ctrDom, 
								dumVal, NULL, NULL, &errNum);
							if (setbackVz){
								inObj->domain.p->voxel_size[0] = xsize;
								inObj->domain.p->voxel_size[1] = ysize;
								inObj->domain.p->voxel_size[2] = zsize;
							}
							if (!errNum && flip && ctrDom.core && ctrDom.ctr->model){
								errNum  = WlzGMFilterFlipOrient(ctrDom.ctr->model);
							}
						}
					}
				}
			}
		(void)WlzStringFromErrorNum(errNum, &errMsgStr);
		return outObj;
	}

	int WlzDestroyObj(WlzObject *obj){
		int success = 0;
		WlzErrorNum   errNum = WLZ_ERR_NONE;
		
		//obj->linkcount = 1;
		errNum = WlzFreeObj(obj) ;
		//obj = NULL;
		if (errNum == WLZ_ERR_NONE){
			success = 1;
		}

		return success;
	}

