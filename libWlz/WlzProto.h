#ifndef WLZ_PROTO_H
#define WLZ_PROTO_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzProto_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzProto.h
* \author       Bill Hill
* \date         March 1999
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
* \brief	Defines the Woolz function prototypes.
*		
*		To allow bindings to other languages (ie Java) to be
*               generated automaticaly from this file, the following
*               conventions MUST be followed:
*		<ul>
*		  <li>
*                 The preprocessor identifier #WLZ_EXT_BIND controls
*                 whether code is bound.
*		  </li>
*		  <li>
*                 All prototypes to be bound must provide both
*                 parameter types and parameter identifiers.
*		  </li>
*		  <li>
*                 If a function returns a WlzErrorNum type it is
*                 assumed to be the error code.
*		  </li>
*		  <li>
*                 If a parameter is used to return the error code then
*                 it must be declared as 'WlzErrorNum *dstErr' and it
*                 may be a good idea to make it the last parameter.
*		  </li>
*		  <li>
*                 If either a function return type or parameter is
*                 'char' and it is not a pointer then it is assummed
*                 to be a byte.
*		  </li>
*		  <li>
*                 All Woolz types and functions start with 'Wlz'.
*		  </li>
*		  <li>
*                 All enum's are typedef'd.
*		  </li>
*		  <li>
*                 All structs are typedef'd.
*		  </li>
*		  <li>
*                 If a parameter identifier starts with 'dst' and is
*                 a pointer then it is assumed to be a destination
*                 pointer.
*		  </li>
*		  <li>
*                 All function return types or parameters of type
*                 'char *' are assumed to be strings unless the
*                 parameter identifier starts with 'dst'.
*		  </li>
*		  <li>
*                 All arrays must either start with 'array' or
*                 'dstArray', with 'dstArray' being used to identify
*                 arrays which are allocated within the function being
*                 called. An array pointer must be preceded by it's
*                 size, the type should be int, WlzIVertex2, ... and
*                 the identifier should be the same as array's but with
*                 'size' prepended as is the example 'sizeArrayName'.
*                 Arrays can not be pointers to void and wrappers may
*                 be needed to avoid this.
*		  </li>
*		  <li>
*		  Because of a parsing bug WlzUByte should be expanded to
*		  unsigned char.
*		</ul>
* \ingroup	Wlz
*/

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#ifdef _WIN32
#define fabsf(x) fabs(x)
#define _setmode setmode
#endif

/************************************************************************
* Wlz2DContains.c							*
************************************************************************/
extern WlzObject 		*Wlz2DContains(
				  WlzObject *obj,
				  double x,
				  double y,
				  WlzErrorNum *dstErr);
/************************************************************************
* Wlz3DProjection.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject 		*WlzGetProjectionFromObject(
				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  Wlz3DProjectionIntFn intFunc,
				  void	*intFuncData,
				  WlzErrorNum *dstErr);
#endif

/************************************************************************
* Wlz3DSection.c							*
************************************************************************/
extern WlzObject 		*WlzGetSectionFromObject(
				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  WlzInterpolationType	interp,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzGetMaskedSectionFromObject(
				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  WlzInterpolationType	interp,
				  WlzErrorNum *dstErr);

/************************************************************************
* Wlz3DSubSection.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject 		*WlzGetSubSectionFromObject(
				  WlzObject 	*obj,
				  WlzObject	*subDomain,
				  WlzThreeDViewStruct *viewStr,
				  WlzInterpolationType	interp,
				  WlzObject	**maskRtn,
				  WlzErrorNum *dstErr);
#endif

/************************************************************************
* Wlz3DSectionSegmentObject.c						*
************************************************************************/
extern WlzErrorNum 		Wlz3DSectionSegmentObject(
				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  int *dstArraySizeObjs,
				  WlzObject ***dstArrayObjs);

/************************************************************************
* Wlz3DViewStructUtils.c						*
************************************************************************/
extern WlzThreeDViewStruct 	*WlzMake3DViewStruct(
				  WlzObjectType type,
				  WlzErrorNum *dstErr);
extern WlzThreeDViewStruct	*WlzMake3DViewStructCopy(
				  WlzThreeDViewStruct *given,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzFree3DViewStruct(
				  WlzThreeDViewStruct *viewStr);
extern WlzErrorNum 		WlzInit3DViewStruct(
				  WlzThreeDViewStruct *viewStr,
				  WlzObject *obj);
extern WlzErrorNum 		WlzInit3DViewStructAffineTransform(
                                  WlzThreeDViewStruct	*viewStr);
extern WlzErrorNum 		Wlz3DViewStructSetupTransformLuts(
                                  WlzThreeDViewStruct	*viewStr);
extern WlzErrorNum 		Wlz3DViewStructTransformBB(
                                  WlzObject	*obj,
                                  WlzThreeDViewStruct	*viewStr);
extern WlzErrorNum 		Wlz3DSectionTransformVtx(
				  WlzDVertex3 *vtx,
				  WlzThreeDViewStruct *viewStr);
extern WlzErrorNum 		Wlz3DSectionTransformVtxR(
		                  WlzThreeDViewStruct *viewStr,
				  WlzDVertex3 vtx,
				  WlzDVertex3 *dstVtx);
extern WlzErrorNum 		Wlz3DSectionTransformInvVtx(
				  WlzDVertex3 *vtx,
				  WlzThreeDViewStruct *viewStr);
extern WlzErrorNum 		Wlz3DSectionTransformInvVtxR(
                                  WlzThreeDViewStruct *viewStr,
				  WlzDVertex3 vtx,
				  WlzDVertex3 *dstVtx);
extern WlzErrorNum 		Wlz3DSectionIncrementDistance(
				  WlzThreeDViewStruct *viewStr,
				  double incr);
extern WlzDVertex2 		Wlz3DViewGetIntersectionPoint(
				  WlzThreeDViewStruct *vS1,
				  WlzThreeDViewStruct *vS2,
				  WlzErrorNum *dstErr);
extern double 			Wlz3DViewGetIntersectionAngle(
				  WlzThreeDViewStruct *vS1,
				  WlzThreeDViewStruct *vS2,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern int 			Wlz3DViewGetBoundingBoxIntersection(
				  WlzThreeDViewStruct *viewStr,
				  WlzDVertex3 *rtnVtxs,
				  WlzErrorNum *dstErr);
extern int 			Wlz3DViewGetGivenBBIntersection(
                                  WlzThreeDViewStruct *viewStr,
				  WlzDVertex3 bbMin,
				  WlzDVertex3 bbMax,
				  WlzDVertex3 *rtnVtxs,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern int			Wlz3DViewGetBoundingBoxIntersectionA(
				  WlzThreeDViewStruct *viewStr,
				  int *dstSizeArrayVtxs,
				  WlzDVertex3 **dstArrayVtxs,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		Wlz3DViewGetFixed(
				  WlzThreeDViewStruct *vs,
				  double *dstX,
				  double *dstY,
				  double *dstZ);
extern WlzErrorNum		Wlz3DViewSetFixed(
				  WlzThreeDViewStruct *vs,	
				  double x,
				  double y,
				  double z);
extern WlzErrorNum		Wlz3DViewGetTheta(
				  WlzThreeDViewStruct *vs,
				  double *dstVal);
extern WlzErrorNum		Wlz3DViewSetTheta(
				  WlzThreeDViewStruct *vs,
				  double val);
extern WlzErrorNum		Wlz3DViewGetPhi(
				  WlzThreeDViewStruct *vs,
				  double *dstVal);
extern WlzErrorNum		Wlz3DViewSetPhi(
				  WlzThreeDViewStruct *vs,
				  double val);
extern WlzErrorNum		Wlz3DViewGetZeta(
				  WlzThreeDViewStruct *vs,
				  double *dstVal);
extern WlzErrorNum		Wlz3DViewSetZeta(
				  WlzThreeDViewStruct *vs,
				  double val);
extern WlzErrorNum		Wlz3DViewGetDist(
				  WlzThreeDViewStruct *vs,
				  double *dstVal);
extern WlzErrorNum		Wlz3DViewSetDist(
				  WlzThreeDViewStruct *vs,
				  double val);
extern WlzErrorNum		Wlz3DViewGetScale(
				  WlzThreeDViewStruct *vs,
				  double *dstVal);
extern WlzErrorNum		Wlz3DViewSetScale(
				  WlzThreeDViewStruct *vs,
				  double val);
extern WlzErrorNum		Wlz3DViewGetViewMode(
				  WlzThreeDViewStruct *vs,
				  WlzThreeDViewMode *dstVal);
extern WlzErrorNum		Wlz3DViewSetViewMode(
				  WlzThreeDViewStruct *vs,
				  WlzThreeDViewMode val);
extern WlzErrorNum		Wlz3DViewGetUp(
				  WlzThreeDViewStruct *vs,
				  double *dstX,
				  double *dstY,
				  double *dstZ);
extern WlzErrorNum		Wlz3DViewSetUp(
				  WlzThreeDViewStruct *vs,
				  double x,
				  double y,
				  double z);
extern WlzErrorNum		Wlz3DViewGetFixed2(
				  WlzThreeDViewStruct *vs,
				  double *dstX,
				  double *dstY,
				  double *dstZ);
extern WlzErrorNum		Wlz3DViewSetFixed2(
				  WlzThreeDViewStruct *vs,
				  double x,
				  double y,
				  double z);
extern WlzErrorNum		Wlz3DViewGetFixedLineAngle(
				  WlzThreeDViewStruct *vs,
				  double *dstVal);
extern WlzErrorNum		Wlz3DViewSetFixedLineAngle(
				  WlzThreeDViewStruct *vs,
				  double val);
extern WlzErrorNum		Wlz3DViewGetMaxvals(
				  WlzThreeDViewStruct *vs,
				  double *dstX,
				  double *dstY,
				  double *dstZ);
extern WlzErrorNum		Wlz3DViewGetMinvals(
				  WlzThreeDViewStruct *vs,
				  double *dstX,
				  double *dstY,
				  double *dstZ);
extern void 			Wlz3DViewGetPlaneEqn(
				  WlzThreeDViewStruct *view,
				  double *dstA,
				  double *dstB,
				  double *dstC,
				  double *dstD);
extern int             		Wlz3DViewIntersectAABB(
				  WlzThreeDViewStruct *view,
                                  WlzDBox3 box);

/************************************************************************
* WlzAutoCor.c
************************************************************************/
extern WlzObject		*WlzAutoCor(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGetSectionFromGMModel.c
************************************************************************/
extern WlzGMModel		*WlzGetSectionFromGMModel(
				  WlzGMModel *gModel,
				  WlzThreeDViewStruct *view,
				  WlzErrorNum *dstErr);

/************************************************************************
* Wlz3DViewTransformObj.c						*
************************************************************************/
extern WlzObject 		*Wlz3DViewTransformObj(
				  WlzObject *srcObj,
				  WlzThreeDViewStruct *viewStr,
				  WlzErrorNum *dstErr);
extern WlzObject 		*Wlz3DViewTransformBitmap(
                                  int arraySizeBitData,
				  unsigned char *arrayBitData,
				  int width,
				  int height,
				  int x_offset,
				  int y_offset,
				  double x,
				  double y,
				  double z,
				  double theta,
				  double phi,
				  double distance,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzAffineTransform.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzAffineTransform	*WlzAffineTransformFromMatrix(
				  WlzTransformType type,
				  double **arrayMat,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzAffineTransform	*WlzAffineTransformFromTranslation(
				  WlzTransformType type,
				  double tX,
				  double tY,
				  double tZ,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformFromScale(
				  WlzTransformType type,
				  double sX,
				  double sY,
				  double sZ,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformFromRotation(
				  WlzTransformType type,
				  double rX,
				  double rY,
				  double rZ,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformFromPrimVal(
				  WlzTransformType type,
				  double trX,
				  double trY,
				  double trZ,
				  double trScale,
				  double trTheta,
				  double trPhi,
				  double trAlpha,
				  double trPsi,
				  double trXsi,
				  int trInvert,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzAffineTransform	*WlzAffineTransformFromPrim(
				  WlzAffineTransformPrim prim);
#endif /* WLZ_EXT_BIND */
extern WlzAffineTransform	*WlzAffineTransformFromSpin(
				  double spX,
				  double spY,
				  double spTheta,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformFromSpinSqueeze(
				  double spX,
				  double spY,
				  double spTheta,
				  double sqX,
				  double sqY,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformProduct(
				  WlzAffineTransform *trans0,
				  WlzAffineTransform *trans1,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformInverse(
				  WlzAffineTransform *trans,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformCopy(
				  WlzAffineTransform *trans,
				  WlzErrorNum *dstErr);
extern int			WlzAffineTransformDimension(
				  WlzAffineTransform *trans,
				  WlzErrorNum *dstErr);
extern int			WlzAffineTransformIsTranslate(
				  WlzAffineTransform *trans,
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern int			WlzAffineTransformIsIdentity(
				  WlzAffineTransform *trans,
				  WlzErrorNum *dstErr);
extern int			WlzAffineTransformIsIdentityTol(
				  WlzAffineTransform *trans,
				  double tolTn,
				  double tolTx,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzAffineTransformObj(
				  WlzObject *srcObj,
				  WlzAffineTransform *trans,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzAffineTransformObjCb(
    				  WlzObject *srcObj,
				  WlzAffineTransform *trans,
				  WlzInterpolationType interp,
				  void *cbData,
				  WlzAffineTransformCbFn cbFn,
				  WlzErrorNum *dstErr);
#endif
extern WlzDVertex2     		WlzAffineTransformVertexD2(
				  WlzAffineTransform *trans,
				  WlzDVertex2 srcVtx,
				  WlzErrorNum *dstErr);
extern WlzDVertex3     		WlzAffineTransformVertexD3(
				  WlzAffineTransform *trans,
				  WlzDVertex3 srcVtx,
				  WlzErrorNum *dstErr);
extern WlzFVertex2     		WlzAffineTransformVertexF2(
				  WlzAffineTransform *trans,
				  WlzFVertex2 srcVtx,
				  WlzErrorNum *dstErr);
extern WlzFVertex3     		WlzAffineTransformVertexF3(
				  WlzAffineTransform *trans,
				  WlzFVertex3 srcVtx,
				  WlzErrorNum *dstErr);
extern WlzIVertex2     		WlzAffineTransformVertexI2(
				  WlzAffineTransform *trans,
				  WlzIVertex2 srcVtx,
				  WlzErrorNum *dstErr);
extern WlzIVertex3     		WlzAffineTransformVertexI3(
				  WlzAffineTransform *trans,
				  WlzIVertex3 srcVtx,
				  WlzErrorNum *dstErr);
extern WlzDVertex2		WlzAffineTransformNormalD2(
				  WlzAffineTransform *trans,
				  WlzDVertex2 srcNrm,
				  WlzErrorNum *dstErr);
extern WlzDVertex3		WlzAffineTransformNormalD3(
				  WlzAffineTransform *trans,
				  WlzDVertex3 srcNrm,
				  WlzErrorNum *dstErr);
extern WlzIBox2        		WlzAffineTransformBBoxI2(
				  WlzAffineTransform *tr,
				  WlzIBox2 srcBox,
				  WlzErrorNum *dstErr);
extern WlzDBox2        		WlzAffineTransformBBoxD2(
				  WlzAffineTransform *tr,
				  WlzDBox2 srcBox,
				  WlzErrorNum *dstErr);
extern WlzIBox3        		WlzAffineTransformBBoxI3(
				  WlzAffineTransform *tr,
				  WlzIBox3 srcBox,
				  WlzErrorNum *dstErr);
extern WlzDBox3        		WlzAffineTransformBBoxD3(
				  WlzAffineTransform *tr,
				  WlzDBox3 srcBox,
				  WlzErrorNum *dstErr);
extern WlzContour      		*WlzAffineTransformContour(
				  WlzContour *srcCtr,
				  WlzAffineTransform *tr,
				  int newModFlg,
				  WlzErrorNum *dstErr);
extern WlzGMModel      		*WlzAffineTransformGMModel(
				  WlzGMModel *srcM,
				  WlzAffineTransform *tr,
				  int newModFlg,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzAffineTransformGMShell(
				  WlzGMShell *shell,
				  WlzAffineTransform *tr);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum	     	WlzAffineTransformPrimGet(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim *prim);
extern WlzErrorNum 		WlzAffineTransformMatrixSet(
				  WlzAffineTransform *trans,
				  double **arrayMat);
#endif /* WLZ_EXT_BIND */
extern WlzErrorNum		WlzAffineTransformTranslationSet(
				  WlzAffineTransform *tr,
				  double tX,
				  double tY,
				  double tZ);
extern WlzErrorNum		WlzAffineTransformScaleSet(
				  WlzAffineTransform *tr,
				  double sX,
				  double sY,
				  double sZ);
extern WlzErrorNum		WlzAffineTransformRotationSet(
				  WlzAffineTransform *tr,
				  double rX,
				  double rY,
				  double rZ);
extern WlzErrorNum 		WlzAffineTransformPrimValSet(
				WlzAffineTransform *trans,
				  double trX,
				  double trY,
				  double trZ,
				  double trScale,
				  double trTheta,
				  double trPhi,
				  double trAlpha,
				  double trPsi,
				  double trXsi,
				  int trInvert);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum 		WlzAffineTransformPrimSet(
				  WlzAffineTransform *tr,
				  WlzAffineTransformPrim prim);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzAffineTransformLSq.c						*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzAffineTransform 	*WlzAffineTransformLSq(
				  WlzVertexType vType,
				  int sizeArrayVT,
				  WlzVertexP arrayVT,
				  int sizeArrayVS,
				  WlzVertexP arrayVS,
				  int sizeArrayVW,
				  double *arrayVW,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzAffineTransform 	*WlzAffineTransformLSq2D(
				  int sizeArrayVT,
				  WlzDVertex2 *arrayVT,
				  int sizeArrayVS,
				  WlzDVertex2 *arrayVS,
				  int sizeArrayVW,
				  double *arrayVW,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSq3D(
				  int sizeArrayVT,
				  WlzDVertex3 *arrayVT,
				  int sizeArrayVS,
				  WlzDVertex3 *arrayVS,
				  int sizeArrayVW,
				  double *arrayVW,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzAffineTransform 	*WlzAffineTransformLSqTrans2D(
				  WlzDVertex2 *vT,
				  WlzDVertex2 *vS,
				  double *vW,
				  int nV,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform  	*WlzAffineTransformLSqTrans3D(
				  WlzDVertex3 *vT,
				  WlzDVertex3 *vS,
				  double *vW,
				  int nV, WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSqGen2D(
				  WlzDVertex2 *vT,
				  WlzDVertex2 *vS,
				  double *vW,
				  int nV,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSqGen3D(
				  WlzDVertex3 *vT,
				  WlzDVertex3 *vS,
				  double *vW,
				  int nV,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSqScale2D(
				  WlzDVertex2 *vT,
				  WlzDVertex2 *vS,
				  double *vW,
				  int nVtx,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSqReg2D(
				  WlzDVertex2 *vT,
				  WlzDVertex2 *vS,
				  double *vW,
				  int nVtx,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform  	*WlzAffineTransformLSqReg3D(
				  WlzDVertex3 *vT,
				  WlzDVertex3 *vS,
				  double *vW,
				  int nV,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSqRegWlz2D(
				  WlzDVertex2 *vT,
				  WlzDVertex2 *vS,
				  int nV,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSqDQ2D(
				  int nV,
				  double *vW,
				  WlzDVertex2 *vT,
				  WlzDVertex2 *vS,
				  int nN,
				  double *nW,
				  WlzDVertex2 *nT,
				  WlzDVertex2 *nS,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSqDQ3D(
				  int nV,
				  double *vW,
				  WlzDVertex3 *vT,
				  WlzDVertex3 *vS,
				  int nN,
				  double *nW,
				  WlzDVertex3 *nT,
				  WlzDVertex3 *nS,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzArea.c								*
************************************************************************/
extern int 			WlzArea(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);


/************************************************************************
* WlzArray.c								*
************************************************************************/
extern WlzErrorNum		WlzToBArray2D(
				  WlzIVertex2 *dstSizeArrayDat,
				  unsigned char ***dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex2 origin,
				  WlzIVertex2 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToIArray2D(
				  WlzIVertex2 *dstSizeArrayDat,
				  int ***dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex2 origin,
				  WlzIVertex2 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToSArray2D(
				  WlzIVertex2 *dstSizeArrayDat,
				  short ***dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex2 origin,
				  WlzIVertex2 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToUArray2D(
				  WlzIVertex2 *dstSizeArrayDat,
				  unsigned char ***dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex2 origin,
				  WlzIVertex2 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToFArray2D(
				  WlzIVertex2 *dstSizeArrayDat,
				  float ***dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex2 origin,
				  WlzIVertex2 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToDArray2D(
				  WlzIVertex2 *dstSizeArrayDat,
				  double ***dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex2 origin,
				  WlzIVertex2 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToBArray3D(
				  WlzIVertex3 *dstSizeArrayDat,
				  unsigned char ****dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex3 origin,
				  WlzIVertex3 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToIArray3D(
				  WlzIVertex3 *dstSizeArrayDat,
				  int ****dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex3 origin,
				  WlzIVertex3 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToSArray3D(
				  WlzIVertex3 *dstSizeArrayDat,
				  short ****dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex3 origin,
				  WlzIVertex3 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToUArray3D(
				  WlzIVertex3 *dstSizeArrayDat,
				  unsigned char ****dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex3 origin,
				  WlzIVertex3 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToFArray3D(
				  WlzIVertex3 *dstSizeArrayDat,
				  float ****dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex3 origin,
				  WlzIVertex3 size,
				  int noiseFlag);
extern WlzErrorNum		WlzToDArray3D(
				  WlzIVertex3 *dstSizeArrayDat,
				  double ****dstArrayDat,
				  WlzObject *srcObj,
				  WlzIVertex3 origin,
				  WlzIVertex3 size,
				  int noiseFlag);
extern WlzObject		*WlzFromBArray1D(
				  WlzIVertex2 arraySizeDat,
				  unsigned char *arrayDat,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromIArray2D(
				  WlzIVertex2 arraySizeDat,
				  int **arrayDat,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromSArray2D(
				  WlzIVertex2 arraySizeDat,
				  short **arrayDat,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromUArray2D(
				  WlzIVertex2 arraySizeDat,
				  unsigned char **arrayDat,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromBArray2D(
				  WlzIVertex2 arraySizeDat,
				  unsigned char **arrayDat,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromFArray2D(
				  WlzIVertex2 arraySizeDat,
				  float **arrayDat,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromDArray2D(
				  WlzIVertex2 arraySizeDat,
				  double **arrayDat,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromBArray3D(
				  WlzIVertex3 arraySizeDat,
				  WlzUByte ***arrayDat,
				  WlzIVertex3 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromIArray3D(
				  WlzIVertex3 arraySizeDat,
				  int ***arrayDat,
				  WlzIVertex3 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromSArray3D(
				  WlzIVertex3 arraySizeDat,
				  short ***arrayDat,
				  WlzIVertex3 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromUArray3D(
				  WlzIVertex3 arraySizeDat,
				  unsigned char ***arrayDat,
				  WlzIVertex3 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromFArray3D(
				  WlzIVertex3 arraySizeDat,
				  float ***arrayDat,
				  WlzIVertex3 arrayOrigin,
				  WlzErrorNum *dstErrNum);
extern WlzObject		*WlzFromDArray3D(
				  WlzIVertex3 arraySizeDat,
				  double ***arrayDat,
				  WlzIVertex3 arrayOrigin,
				  WlzErrorNum *dstErrNum);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzToArray1D(
				  WlzGreyP gP,
				  WlzGreyType gType,
				  WlzIBox3 gBufBox,
				  int gOffset,
				  WlzObject *obj);
extern WlzErrorNum 		WlzToArray2D(
				  void ***dstP,
				  WlzObject *srcObj,
				  WlzIVertex2 size,
				  WlzIVertex2 origin,
				  int noiseFlg,
				  WlzGreyType dstGreyType);
extern WlzErrorNum 		WlzToArray3D(
				  void ****dstP,
				  WlzObject *srcObj,
				  WlzIVertex3 size,
				  WlzIVertex3 origin,
				  int noiseFlg,
				  WlzGreyType dstGreyType);
extern WlzObject		*WlzFromArray1D(
				  WlzObjectType oType,
				  WlzIVertex3 sz,
				  WlzIVertex3 org,
				  WlzGreyType gType,
				  WlzGreyP gDat,
				  int noCopy,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzFromArray2D(
				  void **arrayP,
				   WlzIVertex2 arraySize,
				   WlzIVertex2 arrayOrigin,
				   WlzGreyType dstGreyType,
				   WlzGreyType srcGreyType,
				   double valOffset,
				   double valScale,
				   int clampFlg,
				   int noCopyFlg,
				   WlzErrorNum *dstErr);
extern WlzObject 		*WlzFromArray3D(
				  void ***arrayP,
				 WlzIVertex3 arraySize,
				 WlzIVertex3 arrayOrigin,
				 WlzGreyType dstGreyType,
				 WlzGreyType srcGreyType,
				 double valOffset,
				 double valScale,
				 int clampFlg,
				 int noCopyFlg,
				 WlzErrorNum *dstErr);
extern int 			WlzArrayStats3D(
				  void ***arrayP,
				  WlzIVertex3 arraySize,
				  WlzGreyType greyType,
				  double *dstMin,
				  double *dstMax,
				  double *dstSum,
				  double *dstSumSq,
				  double *dstMean,
				  double *dstStdDev);
extern int			WlzArrayStats2D(
				  void **arrayP,
				  WlzIVertex2 arraySize,
				  WlzGreyType greyType,
				  double *dstMin,
				  double *dstMax,
				  double *dstSum,
				  double *dstSumSq,
				  double *dstMean,
				  double *dstStdDev);
extern int			WlzArrayStats1D(
				  void *arrayP,
				  int arraySize,
				  WlzGreyType greyType,
				  double *dstMin,
				  double *dstMax,
				  double *dstSum,
				  double *dstSumSq,
				  double *dstMean,
				  double *dstStdDev);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzAssign.c								*
************************************************************************/
extern WlzObject		*WlzAssignObject(
				  WlzObject *object,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzDomain		WlzAssignDomain(
				  WlzDomain domain,
				  WlzErrorNum *dstErr);
extern WlzValues 		WlzAssignValues(
				  WlzValues values,
				  WlzErrorNum *dstErr);
extern WlzProperty		WlzAssignProperty(
				  WlzProperty property,
				  WlzErrorNum *dstErr);
extern WlzTransform		WlzAssignTransform(
				  WlzTransform t,
				  WlzErrorNum *dstErr);
extern WlzPropertyList          *WlzAssignPropertyList(
                                  WlzPropertyList *plist,
                                  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAssignAffineTransform(
				  WlzAffineTransform *,
				  WlzErrorNum *dstErr);
extern WlzThreeDViewStruct	*WlzAssign3DViewStruct(
                                  WlzThreeDViewStruct *viewStr,
                                  WlzErrorNum	*dstErr);
extern WlzBoundList 		*WlzAssignBoundList(
				  WlzBoundList *blist,
				  WlzErrorNum *dstErr);
extern WlzPolygonDomain 	*WlzAssignPolygonDomain(
				  WlzPolygonDomain *poly,
				  WlzErrorNum *dstErr);
extern WlzGMModel		*WlzAssignGMModel(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern int 			WlzUnlink(
				  int *linkcount,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzBackground.c							*
************************************************************************/
extern WlzPixelV 		WlzGetBackground(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzSetBackground(
				  WlzObject *obj,
				  WlzPixelV bckgrnd);

/************************************************************************
* WlzBasisFn.c								*
************************************************************************/
extern WlzErrorNum		WlzBasisFnFree(
				  WlzBasisFn *basisFn);
extern WlzDVertex2		WlzBasisFnValueGauss2D(
				  WlzBasisFn *basisFn,
				  WlzDVertex2 srcVx);
extern WlzDVertex2		WlzBasisFnValueIMQ2D(
				  WlzBasisFn *basisFn,
				  WlzDVertex2 srcVx);
extern WlzDVertex3		WlzBasisFnValueIMQ3D(
				  WlzBasisFn *basisFn,
				  WlzDVertex3 srcVx);
extern WlzDVertex2		WlzBasisFnValuePoly2D(
				  WlzBasisFn *basisFn,
				  WlzDVertex2 srcVx);
extern WlzDVertex2		WlzBasisFnValueMQ2D(
				  WlzBasisFn *basisFn,
				  WlzDVertex2 srcVx);
extern WlzDVertex3		WlzBasisFnValueMQ3D(
				  WlzBasisFn *basisFn,
				  WlzDVertex3 srcVx);
extern WlzDVertex2		WlzBasisFnValueTPS2D(
				  WlzBasisFn *basisFn,
				  WlzDVertex2 srcVx);
extern WlzDVertex2		WlzBasisFnValueConf2D(
				  WlzBasisFn *basisFn,
				  WlzDVertex2 srcVx);
WlzDVertex3     		WlzBasisFnValueMOS3D(
				  WlzBasisFn *basisFn,
				  WlzDVertex3 srcVx);
extern double          		WlzBasisFnValueScalarMOS3D(
				  WlzBasisFn *basisFn,
				  WlzDVertex3 srcVx);
extern double          		WlzBasisFnValueMOSPhi(
				  double r,
				  double delta,
				  double tau);
#ifndef WLZ_EXT_BIND
extern WlzBasisFn		*WlzBasisFnGauss2DFromCPts(
				  int nPts,
				  WlzDVertex2 *dPts,
				  WlzDVertex2 *sPts,
				  double delta,
				  WlzBasisFn *prvBasisFn,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr);
extern WlzBasisFn		*WlzBasisFnPoly2DFromCPts(
				  int nPts,
				  int order,
				  WlzDVertex2 *dPts,
				  WlzDVertex2 *sPts,
				  WlzErrorNum *dstErr);
extern WlzBasisFn 		*WlzBasisFnConf2DFromCPts(
				  int nPts,
				  int order,
				  WlzDVertex2 *dPts,
				  WlzDVertex2 *sPts,
				  WlzErrorNum *dstErr);
extern WlzBasisFn		*WlzBasisFnIMQ2DFromCPts(
				  int nPts,
				  WlzDVertex2 *dPts,
				  WlzDVertex2 *sPts,
				  double delta,
				  WlzBasisFn *prvBasisFn,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr);
extern WlzBasisFn		*WlzBasisFnIMQ3DFromCPts(
				  int nPts,
				  WlzDVertex3 *dPts,
				  WlzDVertex3 *sPts,
				  double delta,
				  WlzBasisFn *prvBasisFn,
				  WlzCMesh3D *mesh,
				  WlzErrorNum *dstErr);
extern WlzBasisFn		*WlzBasisFnMQ2DFromCPts(
				  int nPts,
				  WlzDVertex2 *dPts,
				  WlzDVertex2 *sPts,
				  double delta,
				  WlzBasisFn *prvBasisFn,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr);
extern WlzBasisFn		*WlzBasisFnMQ3DFromCPts(
				  int nPts,
				  WlzDVertex3 *dPts,
				  WlzDVertex3 *sPts,
				  double delta,
				  WlzBasisFn *prvBasisFn,
				  WlzCMesh3D *mesh,
				  WlzErrorNum *dstErr);
extern WlzBasisFn		*WlzBasisFnTPS2DFromCPts(
				  int nPts,
				  WlzDVertex2 *dPts,
				  WlzDVertex2 *sPts,
				  WlzBasisFn *prvBasisFn,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr);
extern WlzBasisFn 		*WlzBasisFnMOS3DFromCPts(
				  int nPts,
				  WlzDVertex2 *dPts,
				  WlzDVertex2 *sPts,
				  double *alpha,
				  double *param,
				  WlzErrorNum *dstErr);
extern WlzBasisFn 		*WlzBasisFnScalarMOS3DFromCPts(
				  int nPts,
				  WlzDVertex3 *cPts,
				  double *cVal,
				  double *alpha,
				  double *param,
				  WlzErrorNum *dstErr);
#endif

/************************************************************************
* WlzBasisFnTransform.c							*
************************************************************************/
extern WlzBasisFnTransform	*WlzMakeBasisFnTransform(
				  WlzErrorNum *dstErr);
extern WlzBasisFnTransform	*WlzBasisFnTrFromCPts2D(
				  WlzFnType type,
				  int order,
				  int sizeArrayDPts,
				  WlzDVertex2 *arrayDPts,
				  int sizeArraySPts,
				  WlzDVertex2 *arraySPts,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr);
extern WlzBasisFnTransform	*WlzBasisFnTrFromCPts2DParam(
				  WlzFnType type,
				  int order,
				  int sizeArrayDPts,
				  WlzDVertex2 *arrayDPts,
				  int sizeArraySPts,
				  WlzDVertex2 *arraySPts,
				  WlzCMesh2D *mesh,
				  int sizeArrayParam,
				  double *arrayParam,
				  WlzErrorNum *dstErr);
extern WlzBasisFnTransform	*WlzBasisFnTrFromCPts3DParam(
				  WlzFnType type,
				  int order,
				  int sizeArrayDPts,
				  WlzDVertex3 *arrayDPts,
				  int sizeArraySPts,
				  WlzDVertex3 *arraySPts,
				  WlzCMesh3D *mesh,
				  int sizeArrayParam,
				  double *arrayParam,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzBasisFnTPS2DChangeCPts(
				  WlzBasisFnTransform *basisTr,
				  int nDPts,
				  WlzDVertex2 *dPts,
				  int nSPts,
				  WlzDVertex2 *sPts,
				  WlzObject *cObj);
extern WlzErrorNum		WlzBasisFnTPS2DChangeCPtsParam(
				  WlzBasisFnTransform *basisTr,
				  int nDPts,
				  WlzDVertex2 *dPts,
				  int nSPts,
				  WlzDVertex2 *sPts,
				  WlzObject *cObj,
				  int nParam,
				  double *param);
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzBasisFnMakeCMeshTr(
				  WlzBasisFnTransform *basisTr,
				  WlzCMeshP mesh,
				  WlzErrorNum *dstErr);
#endif
extern WlzObject		*WlzBasisFnMakeCMeshTr2D(
				  WlzBasisFnTransform *basisTr,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzBasisFnMakeCMeshTr3D(
				  WlzBasisFnTransform *basisTr,
				  WlzCMesh3D *mesh,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzBasisFnInvertMakeCMeshTr(
				  WlzBasisFnTransform *basisTr,
				  WlzCMeshP mesh,
				  WlzErrorNum *dstErr);
#endif
extern WlzObject		*WlzBasisFnInvertMakeCMeshTr2D(
				  WlzBasisFnTransform *basisTr,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzBasisFnInvertMakeCMeshTr3D(
				  WlzBasisFnTransform *basisTr,
				  WlzCMesh3D *mesh,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzBasisFnSetMesh(
				  WlzMeshTransform *mesh,
				  WlzBasisFnTransform *basisTr);
extern WlzErrorNum		WlzBasisFnSetCMesh(
				  WlzObject *mObj,
				  WlzBasisFnTransform *basisTr);
extern WlzErrorNum		WlzBasisFnSetCMesh2D(
				  WlzObject *mObj,
				  WlzBasisFnTransform *basisTr);
extern WlzErrorNum		WlzBasisFnSetCMesh3D(
				  WlzObject *mObj,
				  WlzBasisFnTransform *basisTr);
extern WlzErrorNum 		WlzBasisFnFreeTransform(
				  WlzBasisFnTransform *basisTr);
extern WlzObject		*WlzBasisFnTransformObj(
				  WlzObject *srcObj,
				  WlzBasisFnTransform *basisTr,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
extern WlzPolygonDomain		*WlzBasisFnTransformPoly2(
				  WlzPolygonDomain *srcPoly,
				  WlzBasisFnTransform *basisTr,
				  int newPoly,
				  WlzErrorNum *dstErr);
extern WlzBoundList		*WlzBasisFnTransformBoundList(
				  WlzBoundList *srcBnd,
				  WlzBasisFnTransform *basisTr,
				  int newBnd,
				  WlzErrorNum *dstErr);
extern WlzContour 		*WlzBasisFnTransformContour(
				  WlzContour *srcCtr,
				  WlzBasisFnTransform *basisTr,
				  int newCtr,
				  WlzErrorNum *dstErr);
extern WlzGMModel		*WlzBasisFnTransformGMModel(WlzGMModel *srcM,
				  WlzBasisFnTransform *basisTr,
				  int newModel,
				  WlzErrorNum *dstErr);
extern WlzIVertex2		WlzBasisFnTransformVertexI(
				  WlzBasisFnTransform *basisTr,
				  WlzIVertex2 srcVxF,
				  WlzErrorNum *dstErr);
extern WlzFVertex2		WlzBasisFnTransformVertexF(
				  WlzBasisFnTransform *basisTr,
			          WlzFVertex2 srcVxF,
				  WlzErrorNum *dstErr);
extern WlzDVertex2		WlzBasisFnTransformVertexD(
				  WlzBasisFnTransform *basisTr,
			          WlzDVertex2 srcVxF,
				  WlzErrorNum *dstErr);
extern WlzDVertex2     		WlzBasisFnTransformNormalD(
				  WlzBasisFnTransform *basisTr,
				  WlzDVertex2 srcVx,
				  WlzDVertex2 srcNr,
				  WlzDVertex2 *dstVx,
				  WlzErrorNum *dstErr);
extern WlzBasisFnTransform 	*WlzBasisFnTrFromCPts3D(
				  WlzFnType type,
				  int order,
				  int nDPts,
				  WlzDVertex3 *dPts,
				  int nSPts,
				  WlzDVertex3 *sPts,
				  WlzCMesh3D *mesh,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzBoundaryUtils.c							*
************************************************************************/
extern WlzErrorNum 		WlzBoundaryToPolyObjArray(
				  WlzObject *bndObj,
				   int *dstNumObjs,
				   WlzObject ***dstObjArray);
extern WlzErrorNum		WlzBoundObjToPolyDomArray(
				  WlzObject *bndObj,
				  int *dstSizeArrayPoly,
				  WlzPolygonDomain ***dstArrayPoly);
extern int			WlzBoundPolyCount(
				  WlzBoundList *bnd,
				  WlzErrorNum *dstErr);
extern int			WlzBoundVtxCount(
				  WlzBoundList *bnd,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzBoundingBox.c							*
************************************************************************/
extern WlzIBox2 		WlzBoundingBox2I(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzDBox2 		WlzBoundingBox2D(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzIBox3 		WlzBoundingBox3I(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzDBox3 		WlzBoundingBox3D(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzDBox3			WlzBoundingBoxVtx3D(
				  int sizeArrayVtx,
				  WlzDVertex3 *arrayVtx,
				  WlzErrorNum *dstErr);
extern WlzIBox3			WlzBoundingBoxVtx3I(
				  int sizeArrayVtx,
				  WlzIVertex3 *arrayVtx,
				  WlzErrorNum *dstErr);
extern WlzIBox2        		WlzBoundingBox2DTo2I(
				  WlzDBox2 bBox2D);
extern WlzDBox2        		WlzBoundingBox2ITo2D(
				  WlzIBox2 bBox2I);
extern WlzIBox3        		WlzBoundingBox3DTo3I(
				  WlzDBox3 bBox3D);
extern WlzDBox3        		WlzBoundingBox3ITo3D(
				  WlzIBox3 bBox3I);
extern WlzDBox3			WlzBoundingBox3FTo3D(
				  WlzFBox3 bBox3F);
extern WlzFBox3			WlzBoundingBox3DTo3F(
				  WlzDBox3 bBox3D);
extern WlzIBox2			WlzBoundingBoxUnion2I(
				  WlzIBox2 box0,
				  WlzIBox2 box1);
extern WlzFBox2			WlzBoundingBoxUnion2F(
				  WlzFBox2 box0,
				  WlzFBox2 box1);
extern WlzDBox2			WlzBoundingBoxUnion2D(
				  WlzDBox2 box0,
				  WlzDBox2 box1);
extern WlzIBox3			WlzBoundingBoxUnion3I(
				  WlzIBox3 box0,
				  WlzIBox3 box1);
extern WlzFBox3			WlzBoundingBoxUnion3F(
				  WlzFBox3 box0,
				  WlzFBox3 box1);
extern WlzDBox3			WlzBoundingBoxUnion3D(
				  WlzDBox3 box0,
				  WlzDBox3 box1);
extern WlzDBox2			WlzBoundingBoxGModel2D(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzDBox3			WlzBoundingBoxGModel3D(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzBoundToObj.c							*
************************************************************************/
extern WlzObject		*WlzBoundToObj(
				  WlzBoundList *bound,
				  WlzPolyFillMode fillMode,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzBoundaryToObj(
				  WlzObject *boundary,
				  WlzPolyFillMode fillMode,
				  WlzErrorNum *dstErr);

/************************************************************************
 * WlzCannyDeriche.c							*
 ************************************************************************/
extern WlzObject		*WlzCannyDeriche(
				  WlzObject **dstGObj,
				  WlzObject *srcObj,
				  double alpha,
				  double mult,
				  WlzPixelV pMinGrdV,
				  WlzPixelV sMinGrdV,
				  WlzErrorNum *dstErr);

/************************************************************************
 * WlzCentrality.c							*
 ************************************************************************/
 extern double          	WlzCentrality(
 				  WlzObject *fObj,
				  WlzObject *bObj,
				  int nRay,
				  int binFlg,
				  double *dstMaxR,
				  WlzErrorNum *dstErr);

/************************************************************************
 * WlzCentreOfMass.c							*
 ************************************************************************/
extern WlzDVertex2 		WlzCentreOfMass2D(
				  WlzObject *srcObj,
				  int binObjFlg,
				  double *dstMass,
				  WlzErrorNum *dstErr);
extern WlzDVertex3 		WlzCentreOfMass3D(
				  WlzObject *srcObj,
				  int binObjFlg,
				  double *dstMass,
				  WlzErrorNum *dstErr);
extern WlzDVertex2		WlzCentreOfMassVtx2D(
				  int nVtx,
				  WlzDVertex2 *vtx);
extern WlzDVertex3		WlzCentreOfMassVtx3D(
				  int nVtx,
				  WlzDVertex3 *vtx);

/************************************************************************
* WlzClipObjToBox.c							*
************************************************************************/
extern WlzObject		*WlzClipObjToBox2D(
				  WlzObject *srcObj,
				  WlzIBox2 clipBox,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzClipObjToBox3D(
				  WlzObject *srcObj,
				  WlzIBox3 clipBox,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCMeshCurvature.c							*
************************************************************************/
extern WlzObject		*WlzCMeshCurvToImage(
				  WlzObject *inObj,
				  double scale,
				  int meanCrv,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
extern WlzObject       		*WlzCMeshComputeCurvatures(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzObject       		*WlzCMeshComputeNormalsElm(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzObject       		*WlzCMeshComputeNormalsNod(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshComputeCurvaturesFromNodNorm(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCMeshIntersect.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzCMeshIntersect(
				  WlzObject *obj0,
				  WlzObject *obj1,
				  int dsp0,
				  int **dstNodTab,
				  WlzErrorNum *dstErr);
extern WlzCMesh2D		*WlzCMeshIntersect2Mesh2D(
				  WlzMeshTransform *tr0,
				  WlzMeshTransform *tr1,
				  int dsp0,
				  int **dstNodTab,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzObject		*WlzCMeshIntersectDom2D5(
				  WlzObject *sObj,
				  WlzObject *cObj,
				  double delta,
				  double scale,
				  WlzErrorNum *dstErr);
extern WlzDVertex2		WlzCMeshClosePointDom2D5(
				  WlzObject *sObj,
				  WlzObject *cObj,
				  double scale,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCMeshScan.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzCMeshNod2D		*WlzCMeshNextNod2D(
				  WlzCMesh2D *mesh,
				  int *idx,
				  int all);
extern WlzCMeshNod3D		*WlzCMeshNextNod3D(
				  WlzCMesh3D *mesh,
				  int *idx,
				  int all);
extern WlzCMeshElm2D		*WlzCMeshNextElm2D(
				  WlzCMesh2D *mesh,
				  int *idx,
				  int all);
extern WlzCMeshElm3D		*WlzCMeshNextElm3D(
				  WlzCMesh3D *mesh,
				  int *idx,
				  int all);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzCMeshSurfMap.c                                                   *
************************************************************************/
extern WlzObject		*WlzCMeshCompSurfMap(
				  WlzObject *inObj,
				  int nDV,
				  WlzDVertex3 *dV,
				  int nSV,
				  WlzDVertex3 *sV,
				  WlzErrorNum *dstErr);
extern WlzObject       		*WlzCMeshCompSurfMapIdx(
				  WlzCMesh2D5 *mesh,
				  int nP,
				  WlzDVertex3 *dPV,
				  int *pIdx,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshToContour(
				  WlzObject *mObj,
				  double disp,
				  WlzErrorNum *dstErr);
extern WlzGMModel		*WlzCMeshToGMModel(
				  WlzObject *mObj,
				  double disp,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshFlatten2D5(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzCMeshCompSurfMapConformalIdx(
				  WlzCMesh2D5 *mesh,
				  int nPN,
				  WlzDVertex3 *dPV,
				  int *pIdx,
				  double minA,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzCMeshFMar.c.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum     		WlzCMeshFMarNodes2D(
				  WlzCMesh2D *mesh,
				  double *distances,
				  int sizeArraySeedPos,
				  WlzDVertex2 *arraySeedPos);
extern WlzErrorNum     		WlzCMeshFMarNodes3D(
				  WlzCMesh3D *mesh,
				  double *distances,
				  int sizeArraySeedPos,
				  WlzDVertex3 *arraySeedPos);
#endif /* WLZ_EXT_BIND */
extern WlzObject		*WlzCMeshDistance2D(
				  WlzObject *mObj,
				  WlzObjectType rObjType,
				  int sizeArraySeeds,
				  WlzDVertex2 *arraySeeds,
                                  WlzInterpolationType itp,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshDistance3D(
				  WlzObject *mObj,
				  WlzObjectType rObjType,
				  int sizeArraySeeds,
				  WlzDVertex3 *arraySeeds,
                                  WlzInterpolationType itp,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCMeshTransform.c							*
************************************************************************/
extern WlzObject		*WlzCMeshTransformInvert(
				  WlzObject *mObj,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshTransformFromObj(
				  WlzObject *srcObj,
				  WlzMeshGenMethod method,
				  double minDist,
				  double maxDist,
				  WlzObject **dstDilObj,
				  int delOut,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshTransformObj(
				  WlzObject *srcObj,
				  WlzObject *mObj,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
extern WlzCompoundArray		*WlzCMeshTransformManyObjAsIdx(
				  WlzCompoundArray *srcObj,
				  WlzObject *mObj,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzCMeshTransformVtxAry2I(
				  WlzObject *mObj,
				  int sizeArrayVtx,
				  WlzIVertex2 *arrayVtx);
extern WlzErrorNum		WlzCMeshTransformVtxAry3I(
				  WlzObject *mObj,
				  int sizeArrayVtx,
				  WlzIVertex3 *arrayVtx);
extern WlzErrorNum		WlzCMeshTransformVtxAry2F(
				  WlzObject *mObj,
				  int sizeArrayVtx,
				  WlzFVertex2 *arrayVtx);
extern WlzErrorNum		WlzCMeshTransformVtxAry3F(
				  WlzObject *mObj,
				  int sizeArrayVtx,
				  WlzFVertex3 *arrayVtx);
extern WlzErrorNum		WlzCMeshTransformVtxAry2D(
				  WlzObject *mObj,
				  int sizeArrayVtx,
				  WlzDVertex2 *arrayVtx);
extern WlzErrorNum		WlzCMeshTransformVtxAry2D5(
				  WlzObject *mObj,
				  int sizeArrayVtx,
				  WlzDVertex3 *arrayVtx);
extern WlzErrorNum		WlzCMeshTransformVtxAry3D(
				  WlzObject *mObj,
				  int sizeArrayVtx,
				  WlzDVertex3 *arrayVtx);
extern WlzErrorNum		WlzScaleCMeshValue(
                                  double scale,
                                  WlzObject *obj);
extern WlzObject*		WlzCopyScaleCMeshValue(
                                  double scale,
                                  WlzObject *obj,
                                  WlzErrorNum* errNumP);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzCMeshGetNodesAndEdges(
				  WlzObject *mObj,
				  int *dstSizeArrayNod,
				  WlzDVertex2 **dstArrayNod,
				  int *dstSizeArrayDsp,
				  WlzDVertex2 **dstArrayDsp,
				  int *dstSizeArrayEdg,
				  int **dstArrayEdg);
extern WlzObject		*WlzCMeshToDomObj(
				  WlzObject *mObj,
				  int trans,
				  double scale,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshToDomObjValues(
				  WlzObject *dObj,
				  WlzObject *mObj,
                                  WlzInterpolationType itp,
				  WlzErrorNum *dstErr);
extern WlzDBox2			WlzCMeshTransformGetBBox2D(
				  WlzObject *mObj,
				  int trans,
				  WlzErrorNum *dstErr);
extern WlzDBox3			WlzCMeshTransformGetBBox3D(
				  WlzObject *mObj,
				  int trans,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern	WlzErrorNum		WlzCMeshAffineProduct(
				  WlzObject *trC,
				  WlzAffineTransform *trA,
                                  int order);
extern WlzObject       		*WlzCMeshMeshMeshProduct(
				  WlzMeshTransform *tr0,
                                  WlzMeshTransform *tr1,
				  WlzErrorNum *dstErr);
extern WlzObject       		*WlzCMeshMeshProduct(
				  WlzObject *tr0,
                                  WlzMeshTransform *tr1,
				  int order,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshProduct(
				  WlzObject *tr0,
				  WlzObject *tr1,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCMeshUtils.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern double			WlzCMeshElmMaxSqEdgLen2D(
				  WlzCMeshElm2D *elm);
extern double			WlzCMeshElmMaxSqEdgLen2D5(
				  WlzCMeshElm2D5 *elm);
extern double			WlzCMeshElmMaxSqEdgLen3D(
				  WlzCMeshElm3D *elm);
extern void                     WlzCMeshUpdateMaxSqEdgLen2D(
                                  WlzCMesh2D *mesh);
extern void                     WlzCMeshUpdateMaxSqEdgLen2D5(
                                  WlzCMesh2D5 *mesh);
extern void                     WlzCMeshUpdateMaxSqEdgLen3D(
                                  WlzCMesh3D *mesh);
extern void			WlzCMeshUpdateBBox2D(
				  WlzCMesh2D *mesh);
extern void			WlzCMeshUpdateBBox2D5(
				  WlzCMesh2D5 *mesh);
extern void			WlzCMeshUpdateBBox3D(
				  WlzCMesh3D *mesh);
extern void			WlzCMeshSetNodFlags(
				  WlzCMeshP mesh,
				  unsigned int flags);
extern void			WlzCMeshSetNodFlags2D(
				  WlzCMesh2D *mesh,
				  unsigned int flags);
extern void			WlzCMeshSetNodFlags3D(
				  WlzCMesh3D *mesh,
				  unsigned int flags);
extern void			WlzCMeshClearNodFlags(
				  WlzCMeshP mesh,
				  unsigned int flags);
extern void			WlzCMeshClearNodFlags2D(
				  WlzCMesh2D *mesh,
				  unsigned int flags);
extern void			WlzCMeshClearNodFlags3D(
				  WlzCMesh3D *mesh,
				  unsigned int flags);
extern void			WlzCMeshClearElmFlags(
				  WlzCMeshP mesh,
				  unsigned int flags);
extern void			WlzCMeshClearElmFlags2D(
				  WlzCMesh2D *mesh,
				  unsigned int flags);
extern void			WlzCMeshClearElmFlags3D(
				  WlzCMesh3D *mesh,
				  unsigned int flags);
extern void			WlzCMeshElmGetNodes2D(
				  WlzCMeshElm2D *elm,
				  WlzCMeshNod2D **dstNod0,
				  WlzCMeshNod2D **dstNod1,
				  WlzCMeshNod2D **dstNod2);
extern void			WlzCMeshElmGetNodes2D5(
				  WlzCMeshElm2D5 *elm,
				  WlzCMeshNod2D5 **dstNod0,
				  WlzCMeshNod2D5 **dstNod1,
				  WlzCMeshNod2D5 **dstNod2);
extern void			WlzCMeshElmGetNodes3D(
				  WlzCMeshElm3D *elm,
				  WlzCMeshNod3D **dstNod0,
				  WlzCMeshNod3D **dstNod1,
				  WlzCMeshNod3D **dstNod2,
				  WlzCMeshNod3D **dstNod3);
extern void            		WlzCMeshSetVertices2D(
				  WlzCMesh2D *mesh,
				  WlzDVertex2 *vtxBuf,
                                  int update);
extern void            		WlzCMeshSetVertices3D(
				  WlzCMesh3D *mesh,
				  WlzDVertex3 *vtxBuf,
                                  int update);
extern void			WlzCMeshGetCellStats(
				  WlzCMeshP mesh,
				  int *dstMinNodPerCell,
				  int *dstMaxNodPerCell,
				  double *dstMeanNodPerCell,
				  int *dstMinElmPerCell,
				  int *dstMaxElmPerCell,
				  double *dstMeanElmPerCell);
extern void			WlzCMeshGetCellStats2D(
				  WlzCMesh2D *mesh,
				  int *dstMinNodPerCell,
				  int *dstMaxNodPerCell,
				  double *dstMeanNodPerCell,
				  int *dstMinElmPerCell,
				  int *dstMaxElmPerCell,
				  double *dstMeanElmPerCell);
extern void			WlzCMeshGetCellStats3D(
				  WlzCMesh3D *mesh,
				  int *dstMinNodPerCell,
				  int *dstMaxNodPerCell,
				  double *dstMeanNodPerCell,
				  int *dstMinElmPerCell,
				  int *dstMaxElmPerCell,
				  double *dstMeanElmPerCell);
extern int			WlzCMeshCountBoundNodes(
				  WlzCMeshP mesh);
extern int			WlzCMeshCountBoundNodes2D(
				  WlzCMesh2D *mesh);
extern int			WlzCMeshCountBoundNodes2D5(
				  WlzCMesh2D5 *mesh);
extern int			WlzCMeshCountBoundNodes3D(
				  WlzCMesh3D *mesh);
extern int			WlzCMeshSetBoundNodFlags(
				  WlzCMeshP mesh);
extern int			WlzCMeshSetBoundNodFlags2D(
				  WlzCMesh2D *mesh);
extern int			WlzCMeshSetBoundNodFlags2D5(
				  WlzCMesh2D5 *mesh);
extern int			WlzCMeshSetBoundNodFlags3D(
				  WlzCMesh3D *mesh);
extern int			WlzCMeshSetBoundElmFlags(
				  WlzCMeshP mesh);
extern int			WlzCMeshSetBoundElmFlags2D(
				  WlzCMesh2D *mesh);
extern int			WlzCMeshSetBoundElmFlags3D(
				  WlzCMesh3D *mesh);
extern int			WlzCMeshNodIsBoundary2D(
				  WlzCMeshNod2D *nod);
extern int			WlzCMeshNodIsBoundary2D5(
				  WlzCMeshNod2D5 *nod);
extern int			WlzCMeshNodIsBoundary3D(
				  WlzCMeshNod3D *nod);
extern int			WlzCMeshElmIsBoundary2D(
				  WlzCMeshElm2D *elm);
extern int			WlzCMeshElmIsBoundary3D(
				  WlzCMeshElm3D *elm);
extern int 			WlzCMeshSetNodIdxTbl2D(
				  WlzCMesh2D *mesh,
				  int *idxTB);
extern int			WlzCMeshSetElmIdxTbl2D(
			 	  WlzCMesh2D *mesh,
				  int *idxTb);
extern int			WlzCMeshSetElmIdxTbl2D5(
			 	  WlzCMesh2D5 *mesh,
				  int *idxTb);
extern int			WlzCMeshSetElmIdxTbl3D(
			 	  WlzCMesh3D *mesh,
				  int *idxTb);
extern int 			WlzCMeshSetNodIdxTbl2D5(
				  WlzCMesh2D5 *mesh,
				  int *idxTB);
extern int 			WlzCMeshSetNodIdxTbl3D(
				  WlzCMesh3D *mesh,
				  int *idxTB);
extern int			*WlzCMeshMakeNodIdxTbl2D(
				  WlzCMesh2D *mesh,
				  int *dstNNod);
extern int			*WlzCMeshMakeNodIdxTbl2D5(
				  WlzCMesh2D5 *mesh,
				  int *dstNNod);
extern int			*WlzCMeshMakeNodIdxTbl3D(
				  WlzCMesh3D *mesh,
				  int *dstNNod);
extern int			*WlzCMeshMakeElmIdxTbl2D(
				  WlzCMesh2D *mesh,
				  int *dstNElm);
extern int			*WlzCMeshMakeElmIdxTbl2D5(
				  WlzCMesh2D5 *mesh,
				  int *dstNElm);
extern int			*WlzCMeshMakeElmIdxTbl3D(
				  WlzCMesh3D *mesh,
				  int *dstNElm);
extern int			WlzCMeshElmRingNodIndices2D(
				  WlzCMeshElm2D *elm,
				  int *maxIdxBuf,
				  int **idxBuf,
				  WlzErrorNum *dstErr);
extern int			WlzCMeshNodRingElmIndices2D5(
				  WlzCMeshNod2D5 *nod, 
				  int *maxIdxBuf,
				  int **idxBuf, 
				  WlzErrorNum *dstErr);
extern int			WlzCMeshNodRingNodIndices2D(
				  WlzCMeshNod2D *nod, 
				  int *maxIdxBuf,
				  int **idxBuf, 
				  WlzErrorNum *dstErr);
extern int			WlzCMeshNodRingNodIndices2D5(
				  WlzCMeshNod2D5 *nod, 
				  int *maxIdxBuf,
				  int **idxBuf, 
				  WlzErrorNum *dstErr);
extern int			WlzCMeshNodRingNodIndices3D(
				  WlzCMeshNod3D *nod, 
				  int *maxIdxBuf,
				  int **idxBuf, 
				  WlzErrorNum *dstErr);
extern double                   WlzCMeshElmSnArea22D(
                                  WlzCMeshElm2D *elm);
extern double			WlzCMeshElmSqArea22D5(
				  WlzCMeshElm2D5 *elm);
extern double                   WlzCMeshElmSnVolume63D(
                                  WlzCMeshElm3D *elm);
extern WlzDBox2			WlzCMeshElmBBox2D(
				  WlzCMeshElm2D *elm);
extern WlzDBox3			WlzCMeshElmBBox2D5(
				  WlzCMeshElm2D5 *elm);
extern WlzDBox3			WlzCMeshElmBBox3D(
				  WlzCMeshElm3D *elm);
extern WlzErrorNum		WlzCMeshLaplacianSmooth(
				  WlzCMeshP mesh,
				  int itr,
				  double alpha,
				  int doBnd,
				  int update);
extern WlzErrorNum		WlzCMeshGetBoundNodes(
				  WlzCMeshP mesh,
				  int *dstSizeArrayNod,
				  int **dstArrayNod);
extern WlzErrorNum		WlzCMeshGetBoundNodes2D(
				  WlzCMesh2D *mesh,
				  int *dstSizeArrayNod,
				  int **dstArrayNod);
extern WlzErrorNum		WlzCMeshGetBoundNodes2D5(
				  WlzCMesh2D5 *mesh,
				  int *dstSizeArrayNod,
				  int **dstArrayNod);
extern WlzErrorNum		WlzCMeshGetBoundNodes3D(
				  WlzCMesh3D *mesh,
				  int *dstSizeArrayNod,
				  int **dstArrayNod);
extern WlzErrorNum		WlzCMeshLaplacianSmooth2D(
				  WlzCMesh2D *mesh,
				  int itr,
				  double alpha,
				  int doBnd,
				  int update);
extern WlzErrorNum		WlzCMeshLaplacianSmooth3D(
				  WlzCMesh3D *mesh,
				  int itr,
				  double alpha,
				  int doBnd,
				  int update);
extern WlzErrorNum		WlzCMeshLPFilter(
				  WlzCMeshP mesh,
				  double kPB,
				  double kSB,
				  double dPB,
				  double dSB,
				  int maxItr,
				  int doBnd,
				  int update);
extern WlzErrorNum		WlzCMeshLPFilterLM(
				  WlzCMeshP mesh,
				  double lambda,
				  double mu,
				  int nItr,
				  int doBnd,
				  int update);
extern WlzErrorNum		WlzCMeshSetVertices(
				  WlzCMeshP mesh,
				  WlzVertexP vtxBuf,
				  int update);
extern WlzErrorNum		WlzCMeshVerify(
				  WlzCMeshP mesh,
				  void **dstElm,
				  int allErr,
				  FILE *fP);
extern WlzErrorNum 		WlzCMeshVerify2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshElm2D **dstElm,
				 int allErr,
				 FILE *fP);
extern WlzErrorNum 		WlzCMeshVerify2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshElm2D5 **dstElm,
				 int allErr,
				 FILE *fP);
extern WlzErrorNum 		WlzCMeshVerify3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshElm3D **dstElm,
				  int allErr,
				  FILE *fP);
extern WlzErrorNum		WlzCMeshCmpElmFeat(
				  WlzCMeshP mesh,
				  int *dstNElm,
				  int **dstIdx,
				  double **dstVol,
				  double **dstMinLen,
				  double **dstMaxLen);
extern WlzErrorNum		WlzCMeshCmpElmFeat2D(
				  WlzCMesh2D *mesh,
				  int *dstNElm,
				  int **dstIdx,
				  double **dstVol,
				  double **dstMinLen,
				  double **dstMaxLen);
extern WlzErrorNum		WlzCMeshCmpElmFeat3D(
				  WlzCMesh3D *mesh,
				  int *dstNElm,
				  int **dstIdx,
				  double **dstVol,
				  double **dstMinLen,
				  double **dstMaxLen);
extern WlzErrorNum		WlzCMeshFixNegativeElms(
				  WlzCMeshP mesh);
extern WlzErrorNum		WlzCMeshFixNegativeElms2D(
				  WlzCMesh2D *mesh);
extern WlzErrorNum		WlzCMeshFixNegativeElms3D(
				  WlzCMesh3D *mesh);
extern WlzCMeshP		WlzCMeshCopy(
				  WlzCMeshP mesh,
				  int squuze,
				  size_t datSz,
				  AlcVector **newDat,
				  AlcVector *gvnDat,
				  WlzErrorNum *dstErr);
extern WlzCMesh2D		*WlzCMeshCopy2D(
				  WlzCMesh2D *mesh,
				  int squuze,
				  size_t datSz,
				  AlcVector **newDat,
				  AlcVector *gvnDat,
				  WlzErrorNum *dstErr);
extern WlzCMesh2D5		*WlzCMeshCopy2D5(
				  WlzCMesh2D5 *mesh,
				  int squuze,
				  size_t datSz,
				  AlcVector **newDat,
				  AlcVector *gvnDat,
				  WlzErrorNum *dstErr);
extern WlzCMesh3D		*WlzCMeshCopy3D(
				  WlzCMesh3D *mesh,
				  int squuze,
				  size_t datSz,
				  AlcVector **newDat,
				  AlcVector *gvnDat,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzCompDispIncGrey.c						        *
************************************************************************/
extern WlzCompoundArray		*WlzCompDispIncGrey(
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCompoundArrayToScalar.c						*
************************************************************************/
extern WlzObject		*WlzCompoundArrayToScalar(
				  WlzCompoundArray *cpd,
				  WlzBinaryOperatorType bop,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCompThresh.c							*
************************************************************************/
extern WlzErrorNum		WlzCompThreshold(
				  double *dstThrVal,
				  WlzObject *histObj,
				  WlzCompThreshType method,
				  double extraFrac);
extern WlzErrorNum		WlzCompThresholdVT(
				  WlzObject *hObj,
				  WlzCompThreshType method,
				  double param0,
				  double param1,
				  double extraFrac,
				  WlzPixelV *dstTV,
				  WlzThresholdType *dstTType);

/************************************************************************
* WlzConstruct3D.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzConstruct3DObjFromFile(
				  int sizeArrayFileStr,
				  char **arrayFileStr,
				  int plane1,
				  float xSz,
				  float ySz,
				  float zSz,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzObject		*WlzConstruct3DObjFromObj(
				  int sizeArrayObjs,
				  WlzObject **arrayObjs,
				  int plane1,
				  float xSz,
				  float ySz,
				  float zSz,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzContour.c
************************************************************************/
extern WlzContour		*WlzContourObj(
				  WlzObject *srcObj,
				  WlzContourMethod ctrMtd,
				  double ctrVal,
				  double ctrWth,
				  int nrmFlg,
				  WlzErrorNum *dstErr);
extern WlzContour		*WlzContourObjGrd(
				  WlzObject *srcObj,
				  double ctrLo,
				  double ctrHi,
				  double ctrWth,
				  int nrmFlg,
				  WlzErrorNum *dstErr);
extern WlzContour 		*WlzContourGrdObj2D(
				  WlzObject *srcObj,
				  WlzObject *gXObj,
				  WlzObject *gYObj,
				  double grdLo,
				  double grdHi,
				  double ftrPrm,
				  int nrmFlg,
				  WlzErrorNum *dstErr);
extern WlzContour		*WlzContourRBFBndObj3D(
				  WlzObject *srcObj,
				  int bErosion,
				  int bDilation,
				  int sDilation,
				  int sFac,
				  int oFac,
				  double sAlpha,
				  double oAlpha,
				  double delta,
				  double tau,
				  double samFac,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzContour		*WlzContourFromPoints(
				  WlzObject *dObj,
				  WlzVertexType vtxType,
				  int nSPts,
				  WlzVertexP sPts,
				  double sAlpha,
				  int nIPts,
				  WlzVertexP iPts,
				  double iDist,
				  double iAlpha,
				  int nOPts,
				  WlzVertexP oPts,
				  double oDist,
				  double oAlpha,
				  double delta,
				  double tau,
				  double samFac,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzConvertPix.c							*
************************************************************************/
extern WlzObject 		*WlzConvertPix(
				  WlzObject *obj,
				  WlzGreyType nPixType,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzConvertVtx(
                                  WlzObject	*obj,
				  WlzVertexType	newVtxType,
				  WlzErrorNum	*dstErr);
extern WlzBoundList 		*WlzConvertBoundType(
                                  WlzBoundList	*bound,
				  WlzObjectType	type,
				  WlzErrorNum	*dstErr);
extern WlzPolygonDomain 	*WlzConvertPolyType(
                                  WlzPolygonDomain	*pdom,
				  WlzObjectType		type,
				  WlzErrorNum		*dstErr);

/************************************************************************
* WlzConvexHull.c							*
************************************************************************/
extern WlzObject		*WlzObjToConvexHull(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzObjToConvexPolygon(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzConvexHullClarkson.c                                               *
************************************************************************/
#ifndef WLZ_EXT_BIND
extern int			WlzConvHullClarkson2D(
				  WlzDVertex2 *vtx,
				  int n,
				  int **dstIdx,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzConvolve.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzConvolveObj(
				  WlzObject *inObj,
				  WlzConvolution *conv,
				  int newObjFlg,
				  WlzErrorNum	*dstErr);
extern int			WlzConvolveSeqParFn(
				  WlzSeqParWSpace *spWSpace,
				  void *spData);
extern int			WlzConvolutionSum(
				  WlzConvolution *conv);
extern int			WlzConvolutionNormalise(
				  WlzConvolution *conv);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzCopy.c
************************************************************************/
extern WlzObject		*WlzCopyObject(
				  WlzObject *inObj,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzCopyObjectGreyValues(
				  WlzObject *dObj,
				  WlzObject *sObj);
#ifndef WLZ_EXT_BIND
extern WlzDomain		WlzCopyDomain(
				  WlzObjectType inObjType,
				  WlzDomain inDom,
				  WlzErrorNum *dstErr);
extern WlzValues		WlzCopyValues(
				  WlzObjectType inObjType,
				  WlzValues inVal,
				  WlzDomain inDom,
				  WlzErrorNum *dstErr);
extern WlzPropertyList		*WlzCopyPropertyList(
				  WlzPropertyList *gList,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzCCor.c
************************************************************************/
extern double			WlzCCorS2D(
				  WlzObject *obj0,
				  WlzObject *obj1,
				  int unionFlg,
				  int normFlg,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCutObjToBox.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject 		*WlzCutObjToValBox2D(
				  WlzObject *srcObj,
				  WlzIBox2 cutBox,
				  WlzGreyType rtnGreyType,
				  void *valP,
				  int bgdNoise,
				  double bgdMu,
				  double bgdSigma,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCutObjToValBox3D(
				  WlzObject *srcObj,
				  WlzIBox3 cutBox,
				  WlzGreyType rtnGreyType,
				  void *valP,
				  int bgdNoise,
				  double bgdMu,
				  double bgdSigma,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzObject		*WlzCutObjToBox2D(
				  WlzObject *srcObj,
				  WlzIBox2 cutBox,
				  WlzGreyType rtnGreyType,
				  int bgdNoise,
				  double bgdMu,
				  double bgdSigma,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCutObjToBox3D(
				  WlzObject *srcObj,
				  WlzIBox3 cutBox,
				  WlzGreyType rtnGreyType,
				  int bgdNoise,
				  double bgdMu,
				  double bgdSigma,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzDiffDomain.c							*
************************************************************************/
extern WlzObject		*WlzDiffDomain(
				  WlzObject *obj1,
				  WlzObject *obj2,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzDiffDomain3d.c							*
************************************************************************/
extern WlzObject 		*WlzDiffDomain3d(
				  WlzObject *obj1,
				  WlzObject *obj2,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzDilation.c								*
************************************************************************/
extern WlzObject		*WlzDilation(
				  WlzObject *obj,
				  WlzConnectType connectivity,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzDistMetric.c							*
************************************************************************/
extern WlzErrorNum     		WlzDistMetricGM(
				  WlzGMModel *model0, 
				  WlzGMModel *model1,
				  double *dstDistH,
				  double *dstDistM,
				  double *dstDistN,
				  double *dstDistI);
extern WlzErrorNum     		WlzDistMetricDirGM(
				  WlzGMModel *model0,
				  WlzGMModel *model1,
				  double *dstDistH,
				  double *dstDistM,
				  double *dstDistN,
				  double *dstDistI);
extern WlzErrorNum     		WlzDistMetricVertex2D(
				  int n0,
				  WlzDVertex2 *vx0,
				  int n1,
				  WlzDVertex2 *vx1,
				  double *dstDistH,
				  double *dstDistM,
				  double *dstDistN,
				  double *dstDistI);
extern WlzErrorNum     		WlzDistMetricVertex3D(
				  int n0,
				  WlzDVertex3 *vx0,
				  int n1,
				  WlzDVertex3 *vx1,
				  double *dstDistH,
				  double *dstDistM,
				  double *dstDistN,
				  double *dstDistI);
extern WlzErrorNum     		WlzDistMetricDirVertex2D(
				  int n0,
				  WlzDVertex2 *vx0,
				  int n1,
				  WlzDVertex2 *vx1,
				  double *dstDistH,
				  double *dstDistM,
				  double *dstDistN,
				  double *dstDistI);
extern WlzErrorNum     		WlzDistMetricDirVertex3D(
				  int n0,
				  WlzDVertex3 *vx0,
				  int n1,
				  WlzDVertex3 *vx1,
				  double *dstDistH,
				  double *dstDistM,
				  double *dstDistN,
				  double *dstDistI);


/************************************************************************
* WlzDistTransform.c                                                    *
************************************************************************/
extern WlzObject 		*WlzDistanceTransform(
				  WlzObject *forObj,
				  WlzObject *refObj,
				  WlzDistanceType dFn,
				  double dParam,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzDomainFill.c							*
************************************************************************/
extern WlzObject 		*WlzDomainFill(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzDomainUtils.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObjectType		WlzDomainType(
				  WlzDomain dom);
extern void			WlzBitLnSetItv(
				  WlzUByte *bitLn,
				  int iLft,
				  int iRgt,
				  int size);
extern WlzErrorNum		WlzDynItvAdd(
				  WlzIntervalDomain *iDom,
				  WlzDynItvPool *iPool,
				  int line,
				  int iLft,
				  int iLen);
extern WlzErrorNum		WlzDynItvLnFromBitLn(
				  WlzIntervalDomain *iDom,
				  WlzUByte *bitLn,
				  int line,
				  int width,
				  WlzDynItvPool *iPool);
extern WlzErrorNum		WlzStandardPlaneDomain(
				  WlzPlaneDomain *pdom,
			 	  WlzVoxelValues *voxtb);
extern WlzErrorNum		WlzStandardIntervalDomain(
				  WlzIntervalDomain *idom);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzErosion.c								*
************************************************************************/
extern WlzObject		*WlzErosion(
				  WlzObject *obj,
				  WlzConnectType connectivity,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzErosion4.c								*
************************************************************************/
extern WlzObject		*WlzErosion4(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzError.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum 		WlzErrorFromAlg(
				  AlgError algErr);
#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzExplode3D.c							*
************************************************************************/
extern WlzErrorNum 		WlzExplode3D(
				  int *dstArraySizeExpObj,
				  WlzObject ***dstArrayExpObj,
				  WlzObject *srcObj);

/************************************************************************
* WlzFacts.c								*
************************************************************************/
extern WlzErrorNum		WlzObjectFacts(
				  WlzObject *obj,
				  FILE *factsFile,
				  char **dstStr,
				  int verbose);

/************************************************************************
* WlzFillBlankPlanes.c							*
************************************************************************/
extern WlzErrorNum		WlzFillBlankPlanes(
				  WlzObject *obj,
				  int min_domain);

/************************************************************************
* WlzFilterNObjValues.c
************************************************************************/
extern WlzObject		*WlzFilterNObjValues(
				  WlzObject *rObj,
				  int sizeArrayObjs,
				  WlzObject **arrayObjs,
				  int fn,
				  double rank,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzFreeSpace.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern void			*WlzPushFreePtr(
				  void *stack,
				  void *data,
				  WlzErrorNum *dstErr);
extern void			*WlzPopFreePtr(
				  void *stack,
				  void **data,
				  WlzErrorNum *dstErr);
#endif /* !WLZ_EXT_BIND */
extern WlzErrorNum		WlzFreeObj(
				  WlzObject *obj);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzFreeDomain(
				  WlzDomain domain);
extern WlzErrorNum		WlzFreeIntervalDomain(
				  WlzIntervalDomain *idom);
extern WlzErrorNum		WlzFreePlaneDomain(
				  WlzPlaneDomain *pdom);
extern WlzErrorNum		WlzFreeValues(
				  WlzValues values);
extern WlzErrorNum		WlzFreeValueTb(
				  WlzRagRValues *vdom);
extern WlzErrorNum		WlzFreeVoxelValueTb(
				  WlzVoxelValues *voxtab);
extern WlzErrorNum		WlzFreePolyDmn(
				  WlzPolygonDomain *poly);
extern WlzErrorNum		WlzFreeBoundList(
				  WlzBoundList *blist);
extern WlzErrorNum		WlzFreeConvHull(
				  WlzConvHullValues *convh);
extern WlzErrorNum		WlzFreeHistogramDomain(
				  WlzHistogramDomain *hist);
/*extern WlzErrorNum		WlzFreeWarpTrans(
				  WlzWarpTrans *obj);*/
/*extern WlzErrorNum		WlzFreeFMatchObj(
				  WlzFMatchObj *obj);*/
extern WlzErrorNum		WlzFree3DWarpTrans(
				  Wlz3DWarpTrans *obj);
extern WlzErrorNum		WlzFreeContour(
				  WlzContour *ctr);
extern WlzErrorNum 		WlzFreeIndexedValues(
				  WlzIndexedValues *ixv);
#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzFourier.c
************************************************************************/
extern WlzObject       		*WlzFourierTransformObj(
				  WlzObject *iObj,
				  int fwd,
                                  WlzErrorNum *dstErr);

/************************************************************************
* WlzDrawDomain.c
************************************************************************/
extern WlzObject		*WlzDrawDomainObj(
				  WlzDVertex2 org,
				  WlzThreeDViewStruct *view,
				  int keep2D,
				  char *cmdStr,
				  int *dstErrIdx,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGauss.c								*
************************************************************************/
extern WlzObject		*WlzGauss2(
				  WlzObject *obj,
				  double wx,
				  double wy,
				  int x_deriv,
				  int y_deriv,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		Wlz1DConv(
				  WlzSepTransWSpace *stwspc,
				  void *params); 
#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzGaussNoise.c							*
************************************************************************/
extern WlzErrorNum		WlzGaussNoise(
				  WlzObject *rtnObj,
				  WlzPixelV tmplVal);

/************************************************************************
* WlzGeometryTrackUpAndDown_s.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzDVertex3  		*WlzGeometryTrackUpAndDown_s(
				  int numberOfPixelsZ,
				  int startfile,
				  int numOfTrackupDown,
				  double distantForInAndOutPGuid,
				  double distantForInAndOutP,
			          unsigned char **TwoDImageFilesNameList,
				  int numOf2DWlzFiles,
				  int downOrUp,
				  int sectionLength_N,
				  int subSubSectionLength_L,
				  int numberOfSampleP_k,
				  char *surfacePointFileName,
				  char *surfaceInPointFileName,
				  char *surfaceOutPointFileName,
				  int startShell,
				  int endShell,
				  int startSection,
				  int endSection,
				  double minDis,
				  WlzErrorNum *dstErr);

#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzGeoModel.c
************************************************************************/
#ifndef WLZ_EXT_BIND
/* Resource callback function list manipulation. */
extern WlzErrorNum		WlzGMModelAddResCb(
				  WlzGMModel *model,
				  WlzGMCbFn fn,
				  void *data);
extern void			WlzGMModelRemResCb(
				  WlzGMModel *model,
				  WlzGMCbFn fn,
				  void *data);
/* Creation  of geometric modeling elements */
extern WlzGMModel		*WlzGMModelNew(
				  WlzGMModelType modType,
				  int blkSz,
				  int vHTSz,
				  WlzErrorNum *dstErr);
extern WlzGMShell		*WlzGMModelNewS(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzGMFace		*WlzGMModelNewF(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzGMLoopT		*WlzGMModelNewLT(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzGMEdge     	        *WlzGMModelNewE(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzGMEdgeT	     	*WlzGMModelNewET(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzGMVertex     		*WlzGMModelNewV(
			 	 WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzGMVertexT	    	*WlzGMModelNewVT(
				  WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzGMModel		*WlzGMModelNewFromS(
				  WlzGMShell *gS,
				  WlzErrorNum *dstErr);
extern WlzGMModel		*WlzGMModelCopy(
				  WlzGMModel *gM,
				  WlzErrorNum *dstErr);
extern void			WlzGMModelAddVertexToHT(
				  WlzGMModel *model,
				  WlzGMVertex *nV);
/* Freeing  of geometric modeling elements */
extern WlzErrorNum		WlzGMModelFree(
				  WlzGMModel *model);
extern WlzErrorNum     	 	WlzGMModelFreeS(
				  WlzGMModel *model,
				  WlzGMShell *shell);
extern WlzErrorNum		WlzGMModelFreeF(
				  WlzGMModel *model,
				  WlzGMFace *loop);
extern WlzErrorNum     		 WlzGMModelFreeLT(
				  WlzGMModel *model,
			  	  WlzGMLoopT *loopT);
extern WlzErrorNum      	WlzGMModelFreeE(
			  	  WlzGMModel *model,
			  	  WlzGMEdge *edge);
extern WlzErrorNum      	WlzGMModelFreeET(
			  	  WlzGMModel *model,
			  	  WlzGMEdgeT *edgeT);
extern WlzErrorNum		WlzGMModelFreeDT(
			  	  WlzGMModel *model,
			  	  WlzGMDiskT *diskT);
extern WlzErrorNum      	WlzGMModelFreeV(
			  	  WlzGMModel *model,
			  	  WlzGMVertex *vertex);
extern WlzErrorNum      	WlzGMModelFreeVT(
			  	  WlzGMModel *model,
			  	  WlzGMVertexT *vertexT);
/* Deletion of geometric modeling elements along with children and any parents
 * that depend solely on the element being deleted. */
extern WlzErrorNum     		WlzGMModelDeleteV(
			  	  WlzGMModel *model,
			   	  WlzGMVertex *dV);
extern WlzErrorNum     		WlzGMModelDeleteE(
			  	  WlzGMModel *model,
			  	  WlzGMEdge *dE);
extern WlzErrorNum     		WlzGMModelDeleteF(
			  	  WlzGMModel *model,
			  	  WlzGMFace *dL);
extern WlzErrorNum		WlzGMModelDeleteS(
			  	  WlzGMModel *model,
			  	  WlzGMShell *shell);
/* Searching */
extern WlzGMVertex		*WlzGMModelMatchVertexG3D(
			  	  WlzGMModel *model,
			  	  WlzDVertex3 gPos);
extern WlzGMVertex		*WlzGMModelMatchVertexG2D(
			  	  WlzGMModel *model,
			  	  WlzDVertex2 gPos);
/* Geometry access and query functions */
extern WlzGMElemType 		WlzGMModelGetSGeomType(
			  	  WlzGMModel *model);
extern WlzGMElemType 		WlzGMModelGetVGeomType(
			  	  WlzGMModel *model);
extern WlzErrorNum		WlzGMShellSetGBB3D(
			  	  WlzGMShell *shell,
			  	  WlzDBox3 bBox);
extern WlzErrorNum		WlzGMShellGetGBB3D(
			  	  WlzGMShell *shell,
			  	  WlzDBox3 *bBox);
extern WlzErrorNum		WlzGMShellGetGBBV3D(
			  	  WlzGMShell *shell,
			  	  double *vol);
extern WlzErrorNum		WlzGMShellGetGBB3D(
                          	  WlzGMShell *shell,
			  	  WlzDBox3 *bBox);
extern WlzErrorNum		WlzGMShellGetGBB2D(
                          	  WlzGMShell *shell,
			  	  WlzDBox2 *bBox);
extern WlzErrorNum		WlzGMShellSetGBB2D(
			  	  WlzGMShell *shell,
			  	  WlzDBox2 bBox);
extern WlzErrorNum		WlzGMShellSetGBB2D(
			  	  WlzGMShell *shell,
			  	  WlzDBox2 bBox);
extern WlzErrorNum		WlzGMShellSetG3D(
			  	  WlzGMShell *shell,
			  	  int nPnt,
			  	  WlzDVertex3 *pos);
extern WlzErrorNum		WlzGMShellSetG2D(
			  	  WlzGMShell *shell,
			  	  int nPnt,
			  	  WlzDVertex2 *pos);
extern int			WlzGMShellGInBB3D(
			  	  WlzGMShell *shell,
			  	  WlzDVertex3 pos);
extern int			WlzGMShellGInBB2D(
			  	  WlzGMShell *shell,
			  	  WlzDVertex2 pos);
extern WlzErrorNum		WlzGMModelSetSG(
			  	  WlzGMModel *model);
extern WlzErrorNum		WlzGMShellUpdateG3D(
			  	  WlzGMShell *shell,
			  	  WlzDVertex3 pos);
extern WlzErrorNum		WlzGMShellUpdateG2D(
			  	  WlzGMShell *shell,
			  	  WlzDVertex2 pos);
extern WlzErrorNum		WlzGMShellUpdateGBB3D(
			  	  WlzGMShell *shell,
			  	  WlzDBox3 bBox);
extern WlzErrorNum		WlzGMShellUpdateGBB2D(
			  	  WlzGMShell *shell,
			  	  WlzDBox2 bBox);
extern WlzErrorNum		WlzGMVertexSetG3D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex3 pos);
extern WlzErrorNum		WlzGMVertexSetG3N(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex3 pos,
			  	  WlzDVertex3 nrm);
extern WlzErrorNum		WlzGMVertexSetG2D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex2 pos);
extern WlzErrorNum		WlzGMVertexSetG2N(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex2 pos,
			  	  WlzDVertex2 nrm);
extern WlzErrorNum		WlzGMShellDndateG2D(
			  	  WlzGMShell *shell,
			  	  WlzDVertex2 pos);
extern WlzErrorNum		WlzGMShellDndateG3D(
			  	  WlzGMShell *shell,
			  	  WlzDVertex3 pos);
extern WlzErrorNum		WlzGMShellComputeGBB(
			  	  WlzGMShell *shell);
extern WlzErrorNum		WlzGMVertexGetG3N(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex3 *dstPos,
			  	  WlzDVertex3 *dstNrm);
extern WlzErrorNum		WlzGMVertexGetG2N(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex2 *dstPos,
			  	  WlzDVertex2 *dstNrm);
extern WlzErrorNum		WlzGMVertexGetG3D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex3 *dstPos);
extern WlzErrorNum		WlzGMVertexGetG2D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex2 *dstPos);
extern WlzDVertex3		WlzGMVertexCmp3D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex3 pos);
extern WlzDVertex2		WlzGMVertexCmp2D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex2 pos);
extern int			WlzGMVertexCmpSign3D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex3 pos);
extern int			WlzGMVertexCmpSign2D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex2 pos);
extern double			WlzGMVertexDistSq3D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex3 pos);
extern double			WlzGMVertexDistSq2D(
			  	  WlzGMVertex *vertex,
			  	  WlzDVertex2 pos);
extern WlzErrorNum		WlzGMFaceGetG3D(
				  WlzGMFace *face,
				  WlzDVertex3 *dstPos0,
				  WlzDVertex3 *dstPos1,
				  WlzDVertex3 *dstPos2);
extern double			WlzGMVertexShellDist(
			  	  WlzGMVertex *v0,
			  	  WlzGMVertex *v1,
			  	  double maxDist,
			  	  WlzErrorNum *dstErr);
extern WlzDVertex3		WlzGMVertexNormal3D(
			  	  WlzGMModel *model,
			  	  WlzGMVertex *gV,
			  	  int *sVBufSz,
			  	  WlzGMVertex ***sVBuf,
			  	  WlzErrorNum *dstErr);
extern WlzGMVertex 		*WlzGMModelLoopTMaxMinCurv2D(
			  	  WlzGMLoopT *gLT,
			  	  int minLen,
			  	  int lnLen,
			  	  WlzBinaryOperatorType mOrM,
			  	  double *dstAlg);
/* Model access and testing */
extern WlzErrorNum 		WlzGMModelTypeValid(
				  WlzGMModelType type);
extern int			WlzGMModelGetDimension(
			  	  WlzGMModel *model,
			  	  WlzErrorNum *dstErr);
extern WlzGMResource		*WlzGMModelGetRes(
				  WlzGMModel *model,
				  WlzGMElemType elemType,
				  WlzErrorNum *dstErr);
/* Topology validity checks (useful for debugging) */
extern WlzErrorNum		WlzGMVerifyModel(
			  	  WlzGMModel *model,
			  	  WlzGMElemP *dstElmP);
extern WlzErrorNum		WlzGMVerifyShell(
			  	  WlzGMShell *shell,
			  	  WlzGMElemP *dstElmP);
extern WlzErrorNum		WlzGMVerifyLoopT(
			  	  WlzGMLoopT *loopT,
			  	  WlzGMElemP *dstElmP);
/* Topology query */
extern WlzGMEdge		**WlzGMModelFindNMEdges(
			  	  WlzGMModel *model,
			  	  int *dstNMCnt,
			  	  WlzErrorNum *dstErr);
extern WlzGMLoopT		*WlzGMEdgeTCommonLoopT(
			  	  WlzGMEdgeT *eT0,
			  	  WlzGMEdgeT *eT1);
extern WlzGMVertex		*WlzGMEdgeCommonVertex(
			  	  WlzGMEdge *eE0,
			  	  WlzGMEdge *eE1);
extern WlzGMVertex		*WlzGMEdgeCommonVertexGetDiskTs(
			  	  WlzGMEdge *eE0,
			  	  WlzGMEdge *eE1,
			  	  WlzGMDiskT **dstDT0,
			  	  WlzGMDiskT **dstDT1);
extern WlzGMDiskT		*WlzGMEdgeCommonDiskT(
			  	  WlzGMEdge *eE0,
			  	  WlzGMEdge *eE1);
extern WlzGMShell		*WlzGMEdgeGetShell(
			  	  WlzGMEdge *eE);
extern WlzGMFace		*WlzGMEdgeCommonFace(
			  	  WlzGMEdge *eE0,
			  	  WlzGMEdge *eE1);
extern WlzGMEdge		*WlzGMVertexCommonEdge(
			  	  WlzGMVertex *eV0,
			  	  WlzGMVertex *eV1);
extern WlzGMShell		*WlzGMVertexCommonShell(
			  	  WlzGMVertex *eV0,
			  	  WlzGMVertex *eV1);
extern WlzGMShell		*WlzGMVertexGetShell(
			  	  WlzGMVertex *eV);
/* Model list management */
extern void	   		WlzGMVertexTAppend(
			  	  WlzGMVertexT *eVT,
			  	  WlzGMVertexT *nVT);
extern void			WlzGMDiskTAppend(
			  	  WlzGMDiskT *eDT,
			  	  WlzGMDiskT *nDT);
extern void			WlzGMEdgeTAppend(
			  	  WlzGMEdgeT *eET,
			  	  WlzGMEdgeT *nET);
extern void			WlzGMEdgeTInsert(
			  	  WlzGMEdgeT *eET,
			  	  WlzGMEdgeT *nET);
extern void			WlzGMEdgeTInsertRadial(
			  	  WlzGMEdgeT *nET);
extern void	   		WlzGMLoopTAppend(
			  	  WlzGMLoopT *pS,
			  	  WlzGMLoopT *nLT);
extern void	   		WlzGMShellAppend(
			  	  WlzGMShell *pS,
			  	  WlzGMShell *newS);
extern void			WlzGMEdgeTUnlink(
			  	  WlzGMEdgeT *dLT);
extern void			WlzGMVertexTUnlink(
			  	  WlzGMVertexT *dLT);
extern void			WlzGMDiskTUnlink(
			  	  WlzGMDiskT *dLT);
extern void			WlzGMDiskTJoin(
			  	  WlzGMDiskT *gDT0,
			  	  WlzGMDiskT *gDT1);
extern void			WlzGMLoopTUnlink(
			  	  WlzGMLoopT *dLT);
extern void			WlzGMShellUnlink(
			  	  WlzGMShell *dS);
extern void			WlzGMShellJoinAndUnlink(
			  	  WlzGMShell *eShell,
			  	  WlzGMShell *dShell);
extern WlzGMResIdxTb		*WlzGMModelResIdx(
			  	  WlzGMModel *model,
			  	  unsigned int eMsk,
			  	  WlzErrorNum *dstErr);
extern void			WlzGMModelResIdxFree(
			  	  WlzGMResIdxTb *resIdxTb);
extern void			WlzGMModelRemVertex(
			  	  WlzGMModel *model,
			  	  WlzGMVertex *dV);
extern WlzErrorNum		WlzGMModelRehashVHT(
			  	  WlzGMModel *model,
			  	  int vHTSz);
/* Model construction */
extern WlzErrorNum		WlzGMModelConstructS(
			  	  WlzGMModel *cM,
			  	  WlzGMShell *gS);
extern WlzErrorNum		WlzGMModelConstructSimplex3D(
			  	  WlzGMModel *model,
			  	  WlzDVertex3 *pos);
extern WlzErrorNum		WlzGMModelConstructSimplex3N(
			  	  WlzGMModel *model,
			  	  WlzDVertex3 *pos,
			  	  WlzDVertex3 *nrm);
extern WlzErrorNum		WlzGMModelConstructSimplex2D(
			  	  WlzGMModel *model,
			  	  WlzDVertex2 *pos);
extern WlzErrorNum		WlzGMModelConstructSimplex2N(
			  	  WlzGMModel *model,
			  	  WlzDVertex2 *pos,
			  	  WlzDVertex2 *nrm);
/* Model Features */
extern int			WlzGMShellSimplexCnt(
			  	  WlzGMShell *gShell);
#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzGeoModelCellGridWSp.c
************************************************************************/
extern WlzGMGridWSp3D		*WlzGeoModelGridWSpNew3D(
				  WlzGMModel *model,
				  WlzGMElemType elemType,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzGeoModelGridWSpSet3D(
				  WlzGMGridWSp3D *grid,
				  WlzGMModel *model,
				  WlzGMElemType elemType);
extern WlzErrorNum		WlzGeoModelGridFree3D(
				  WlzGMGridWSp3D *grid);
extern WlzIBox3			WlzGeoModelGridCellsInDBox(
				  WlzGMGridWSp3D *grid,
				  WlzDBox3 dBox);

/************************************************************************
* WlzGeoModelCut.c
************************************************************************/
extern WlzGMModel		*WlzGMModelCut(WlzGMModel *given,
				  WlzGMModel *knife,
				  WlzErrorNum *dstErr);
extern WlzGMModel		*WlzGMModelCutDom(WlzGMModel *given,
				  WlzObject *knife,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGeoModelFilters.c
************************************************************************/
extern WlzErrorNum		WlzGMFilterFlipOrient(
				  WlzGMModel *model);
extern WlzErrorNum		WlzGMFilterRmSmShells(
				  WlzGMModel *model,
				  int maxElm);
extern WlzErrorNum		WlzGMFilterGeomLP(
				  WlzGMModel *model,
				  double kPB,
				  double kSB,
				  double dPB,
				  double dSB,
				  int maxItr,
				  int nonMan);
extern WlzErrorNum     		WlzGMFilterGeomLPParam(
				  double *dstLambda,
				  double *dstMu,
				  int *dstNItr,
				  double kPB,
				  double kSB,
				  double dPB,
				  double dSB);
extern WlzErrorNum     		WlzGMFilterGeomLPLM(
				  WlzGMModel *model,
				  double lambda,
				  double mu,
				  int nItr,
				  int nonMan);

/************************************************************************
* WlzGeoModelStats.c
************************************************************************/
extern int			WlzGMModelSpxStats(
				  WlzGMModel *model,
				  double *dstMin,
				  double *dstMax,
				  double *dstSum,
				  double *dstSumSq,
				  double *dstMean,
				  double *dstStdDev,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGeometry.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern int			WlzGeomTriangleCircumcentre(
				  WlzDVertex2 *dstCcVx,
				  WlzDVertex2 vx0,
				  WlzDVertex2 vx1,
				  WlzDVertex2 vx2);
extern int			WlzGeomRectFromWideLine(
				  WlzDVertex2 s,
				  WlzDVertex2 t,
				  double width,
				  WlzDVertex2 *v);
#endif /* !WLZ_EXT_BIND */
extern WlzDVertex2     		WlzGeomTriangleCen2D(
				  WlzDVertex2 v0,
				  WlzDVertex2 v1,
                                  WlzDVertex2 v2);
extern WlzDVertex3     		WlzGeomTriangleCen3D(
				  WlzDVertex3 v0,
				  WlzDVertex3 v1,
                                  WlzDVertex3 v2);
extern WlzDVertex3     		WlzGeomTetrahedronCen3D(
				  WlzDVertex3 v0,
				  WlzDVertex3 v1,
                                  WlzDVertex3 v2,
				  WlzDVertex3 v3);
extern int			WlzGeomVxInTriangle2D(
				  WlzDVertex2 p0,
				  WlzDVertex2 p1,
			          WlzDVertex2 p2,
				  WlzDVertex2 pP);
extern int			WlzGeomVxInTriangle3D(
 				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  WlzDVertex3 p2,
				  WlzDVertex3 pG,
				  double ppMax);
extern int			WlzGeomVxInTetrahedron(
				  WlzDVertex3 vx0,
				  WlzDVertex3 vx1,
			          WlzDVertex3 vx2,
			          WlzDVertex3 vx3,
				  WlzDVertex3 vxP);
extern int 			WlzGeomInTriangleCircumcircle(
				  WlzDVertex2 vx0,
				  WlzDVertex2 vx1,
			          WlzDVertex2 vx2,
				  WlzDVertex2 gVx);
extern double			WlzGeomTriangleSnArea2(
				  WlzDVertex2 vx0,
				  WlzDVertex2 vx1,
				  WlzDVertex2 vx2);
extern double			WlzGeomTetraSnVolume6(
				  WlzDVertex3 vx0,
				  WlzDVertex3 vx1,
				  WlzDVertex3 vx2,
				  WlzDVertex3 vx3);
extern double			WlzGeomTriangleArea2Sq3(
				  WlzDVertex3 vx0,
				  WlzDVertex3 vx1,
				  WlzDVertex3 vx2);
extern double			WlzGeomTetraInSphereDiam(
				  WlzDVertex3 vx0,
				  WlzDVertex3 vx1,
				  WlzDVertex3 vx2,
				  WlzDVertex3 vx3);
extern double			WlzGeomTetraInSphereRegDiam(
				  double side);
extern int			WlzGeomLineSegmentsIntersect(
				  WlzDVertex2 p0,
				  WlzDVertex2 p1,
				  WlzDVertex2 q0,
				  WlzDVertex2 q1,
				  WlzDVertex2 *dstN);
extern int			WlzGeomCmpAngle(
				  WlzDVertex2 p0,
				  WlzDVertex2 p1);
extern int             		WlzGeomVtxEqual2D(
				  WlzDVertex2 pos0,
				  WlzDVertex2 pos1,
                                  double tolSq);
extern int             		WlzGeomVtxEqual3D(
				  WlzDVertex3 pos0,
				  WlzDVertex3 pos1,
                                  double tol);
extern void			WlzGeomVtxSortRadial(
				  int nV,
				  WlzDVertex3 *vP,
				  int *idxBuf,
				  WlzDVertex2 *wV,
				  WlzDVertex3 rV);
extern WlzDVertex3		WlzGeomTriangleNormal(
				  WlzDVertex3 v0,
				  WlzDVertex3 v1,
				  WlzDVertex3 v2);
extern int			WlzGeomTriangleAABBIntersect2D(
				  WlzDVertex2 t0,
				  WlzDVertex2 t1,
				  WlzDVertex2 t2,
				  WlzDVertex2 b0,
				  WlzDVertex2 b1,
				  int bbInt);
extern int			WlzGeomTriangleAABBIntersect3D(
				  WlzDVertex3 t0,
				  WlzDVertex3 t1,
				  WlzDVertex3 t2,
				  WlzDVertex3 b0,
				  WlzDVertex3 b1,
				  int bbInt);
extern int			WlzGeomTetrahedronAABBIntersect3D(
				  WlzDVertex3 t0,
				  WlzDVertex3 t1,
				  WlzDVertex3 t2,
				  WlzDVertex3 t3,
				  WlzDVertex3 b0,
				  WlzDVertex3 b1,
				  int bbInt);
extern int			WlzGeomPlaneAABBIntersect(
				  double a,
				  double b,
				  double c,
				  double d,
				  WlzDBox3 box);
extern int			WlzGeomTriangleTriangleIntersect2D(
				  WlzDVertex2 s0,
				  WlzDVertex2 s1,
				  WlzDVertex2 s2,
				  WlzDVertex2 t0,
				  WlzDVertex2 t1,
				  WlzDVertex2 t2);
extern int			WlzGeomTriangleTriangleIntersect3D(
				  WlzDVertex3 s0,
				  WlzDVertex3 s1,
				  WlzDVertex3 s2,
				  WlzDVertex3 t0,
				  WlzDVertex3 t1,
				  WlzDVertex3 t2);
extern int			WlzGeomPlaneLineIntersect(
				  double a,
				  double b,
				  double c,
				  double d,
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  WlzDVertex3 *dstIsn);
extern double			WlzGeomEllipseVxDistSq(
				  WlzDVertex2 centre,
				  WlzDVertex2 sAx,
				  WlzDVertex2 gPnt);
extern int			WlzGeomPlaneTriangleIntersect(
                                  double a, double b,
				  double c, double d,
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  WlzDVertex3 p2,
				  WlzDVertex3 *dstIsn0,
				  WlzDVertex3 *dstIsn1);
extern unsigned int    		WlzGeomHashVtx2D(
				  WlzDVertex2 pos,
				  double tol);
extern unsigned int    		WlzGeomHashVtx3D(
				  WlzDVertex3 pos,
				  double tol);
extern int			WlzGeomCmpVtx2D(
				  WlzDVertex2 pos0,
				  WlzDVertex2 pos1,
				  double tol);
extern int			WlzGeomCmpVtx3D(
				  WlzDVertex3 pos0,
				  WlzDVertex3 pos1,
				  double tol);
extern WlzDVertex2		WlzGeomUnitVector2D(
				  WlzDVertex2 vec);
extern WlzDVertex3		WlzGeomUnitVector3D(
				  WlzDVertex3 vec);
extern WlzDVertex2		WlzGeomUnitVector2D5(
				  WlzDVertex2 pos1,
				  WlzDVertex2 pos0);
extern WlzDVertex3		WlzGeomUnitVector3D2(
				  WlzDVertex3 pos1,
				  WlzDVertex3 pos0);
extern int             		WlzGeomVertexInDiamCircle(
				  WlzDVertex2 lPos0,
				  WlzDVertex2 lPos1,
				  WlzDVertex2 pos);
extern int			WlzGeomItrSpiralRing(
				  int step);
extern int			WlzGeomItrSpiral2I(
				  int step,
				  int *dstX,
				  int *dstY);
extern int			WlzGeomItrSpiralShell(
				  int step);
extern int			WlzGeomItrSpiral3I(
				  int step,
				  int *dstX,
				  int *dstY,
				  int *dstZ);
extern double          		WlzGeomDist2D(
				  WlzDVertex2 v0,
				  WlzDVertex2 v1);
extern double          		WlzGeomDist3D(
				  WlzDVertex3 v0,
				  WlzDVertex3 v1);
extern double          		WlzGeomDistSq2D(
				  WlzDVertex2 v0,
				  WlzDVertex2 v1);
extern double          		WlzGeomDistSq3D(
				  WlzDVertex3 v0,
				  WlzDVertex3 v1);
extern int			WlzGeomTriangleAffineSolve(
				  double *xTr,
				  double *yTr,
				  double dd,
				  WlzDVertex2 *sVx,
				  WlzDVertex2 *dVx,
				  double thresh);
extern int			WlzGeomTetraAffineSolve(
				  double *tr,
				  WlzDVertex3 *sVx,
				  WlzDVertex3 *dVx,
				  double thresh);
extern WlzDVertex2		WlzGeomObjLineSegIntersect2D(
				  WlzObject *obj,
				  WlzDVertex2 p0,
				  WlzDVertex2 p1,
				  double tol,
				  int inside,
				  int method,
				  int *dstStat);
extern WlzDVertex3		WlzGeomObjLineSegIntersect3D(
				  WlzObject *obj,
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  double tol,
				  int inside,
				  int method,
				  int *dstStat);
extern double                   WlzGeomPolar2D(
                                  WlzDVertex2 org,
                                  WlzDVertex2 dest,
                                  double *dstRad);
extern double			WlzGeomCos3V(
				  WlzDVertex2 v0,
				  WlzDVertex2 v1,
				  WlzDVertex2 v2);
extern int             		WlzGeomVtxOnLineSegment2D(
				  WlzDVertex2 tst,
                                  WlzDVertex2 seg0,
				  WlzDVertex2 seg1,
                                  double tol);
extern int             		WlzGeomVtxOnLineSegment3D(
				  WlzDVertex3 tst,
                                  WlzDVertex3 seg0,
				  WlzDVertex3 seg1,
                                  WlzDVertex3 *dstN);
extern double			WlzGeomArcLength2D(
				  WlzDVertex2 a,
				  WlzDVertex2 b,
				  WlzDVertex2 c);
extern WlzDVertex3		WlzGeomLinePlaneIntersection(
				  WlzDVertex3 v,
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  WlzDVertex3 p2,
				  WlzDVertex3 p3,
				  int *dstPar);
extern int			WlzGeomLineTriangleIntersect3D(
				  WlzDVertex3 org,
				  WlzDVertex3 dir,
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  WlzDVertex3 p2,
				  int *dstPar,
				  double *dstT,
				  double *dstU,
				  double *dstV);
extern int			WlzGeomLineLineSegmentIntersect3D(
				  WlzDVertex3 r0,
				  WlzDVertex3 rD,
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  WlzDVertex3 *dstN);
extern int			WlzGeomVtxOnLine3D(
				  WlzDVertex3 p0,
				  WlzDVertex3 r0,
				  WlzDVertex3 rD);
extern double			WlzGeomInterpolateTri2D(
				  WlzDVertex2 p0,
				  WlzDVertex2 p1,
                                  WlzDVertex2 p2,
				  double v0,
				  double v1,
				  double v2,
				  WlzDVertex2 pX);
extern double			WlzGeomInterpolatePoly2D(
				  int n,
				  WlzDVertex2 *p,
				  double *v,
				  double *w,
				  WlzDVertex2 q,
				  WlzErrorNum *dstErr);
extern double			WlzGeomInterpolateConvexPoly2D(
				  int n,
				  WlzDVertex2 *p,
				  double *v,
				  double *w,
				  WlzDVertex2 q);
extern double			WlzGeomInterpolateTet3D(
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
                                  WlzDVertex3 p2,
                                  WlzDVertex3 p3,
				  double v0,
				  double v1,
				  double v2,
				  double v3,
				  WlzDVertex3 pX);
extern void			WlzGeomMap3DTriangleTo2D(
				  WlzDVertex3 p0,
				  WlzDVertex3 p1,
				  WlzDVertex3 p2,
				  WlzDVertex2 *dstQ1,
				  WlzDVertex2 *dstQ2);
extern WlzErrorNum     		WlzGeomCurvature(
				  int dstArraySizeCrv,
				  double *dstArrayCrv,
				  WlzDVertex3 nrm,
                                 int sizeArrayVtx,
				 WlzDVertex3 *arrayVtx);
extern WlzErrorNum		WlzGeometryLSqOPlane(
				  WlzDVertex3 *dstNrm,
				  WlzDVertex3 *dstCen,
				  int sizeArrayVtx,
				  WlzDVertex3 *arrayVtx);
extern double			WlzGeomTriangleVtxDistSq3D(
				  WlzDVertex3 *dstPT,
				  int *dstZT,
				  int *dstIT,
				  double *dstL0,
				  double *dstL1,
				  WlzDVertex3 vT,
				  WlzDVertex3 v0,
				  WlzDVertex3 v1,
				  WlzDVertex3 v2);
extern double          		WlzGeomTriangleVtxDistSq2D(
				  WlzDVertex2 *dstU,
				  int *dstEI,
				  WlzDVertex2 vT,
				  WlzDVertex2 v0,
				  WlzDVertex2 v1,
				  WlzDVertex2 v2);
extern double			WlzGeomTetrahedronVtxDistSq3D(
				  WlzDVertex3 *dstU,
				  int *dstFI,
				  WlzDVertex3 vT,
				  WlzDVertex3 v0,
				  WlzDVertex3 v1,
				  WlzDVertex3 v2,
				  WlzDVertex3 v3);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzGeomPolyTriangulate2D(
				  int sizeArrayPVtx,
				  WlzDVertex2 *arrayPVtx,
				  int *dstSizeArrayTri,
				  int **dstArrayTri);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzGreyCrossing.c							*
************************************************************************/
extern WlzObject		*WlzGreyCrossing(
				  WlzObject *inObj,
				  int newObjFlg,
				  int cVal,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyDitherObj.c							*
************************************************************************/
extern WlzObject 		*WlzGreyDitherObj(WlzObject *o,
				   unsigned int  destBits,
				   WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyGradient.c							*
************************************************************************/
extern WlzObject 		*WlzGreyGradient(
				  WlzObject **dstGrdZ,
				  WlzObject **dstGrdY,
				  WlzObject **dstGrdX,
				  WlzObject *srcObj,
				  WlzRsvFilter *flt,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyInvertMinMax.c							*
************************************************************************/
extern WlzErrorNum 		WlzGreyInvertMinMax(
				  WlzObject *obj,
				  WlzPixelV min,
				  WlzPixelV max);

/************************************************************************
* WlzGreyMask.c								*
************************************************************************/
extern WlzObject 		*WlzGreyMask(
				  WlzObject *obj,
				  WlzObject *mask,
				  WlzPixelV maskVal,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyModGradient.c							*
************************************************************************/
extern WlzObject 		*WlzGreyModGradient(WlzObject *obj,
				   double width,
				   WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyNormalise.c							*
************************************************************************/
extern WlzErrorNum 		WlzGreyNormalise(
				  WlzObject *obj,
				  int dither);

/************************************************************************
* WlzGreyRange.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum 		WlzGreyRange(
				  WlzObject *obj,
				  WlzPixelV *dstMin,
				  WlzPixelV *dstMax);
#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzGreyScan.c								*
************************************************************************/
extern WlzErrorNum 		WlzInitGreyScan(
				  WlzObject *obj,
				  WlzIntervalWSpace *iwsp,
				  WlzGreyWSpace *gwsp);
extern WlzErrorNum 		WlzInitGreyRasterScan(
				  WlzObject *obj,
				  WlzIntervalWSpace *iwsp,
				  WlzGreyWSpace *gwsp,
				  WlzRasterDir raster,
				  int tranpl);
extern WlzErrorNum 		WlzInitGreyWSpace(
				  WlzObject *obj,
				  WlzIntervalWSpace *iwsp,
				  WlzGreyWSpace *gwsp,
				  int tranpl);
extern WlzErrorNum 		WlzNextGreyInterval(
				  WlzIntervalWSpace *iwsp);
extern WlzErrorNum 		WlzGreyInterval(
				  WlzIntervalWSpace *iwsp);

/************************************************************************
* WlzGreySetIncValues.c  						*
************************************************************************/
extern WlzObject		*WlzGreyNewHilbertRankValues(
				  WlzObject *in,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzGreySetHilbertRankValues(
				  WlzObject *in,
				  unsigned int *val);

/************************************************************************
* WlzGreySetIncValues.c  						*
************************************************************************/
extern WlzObject		*WlzGreyNewIncValues(
				  WlzObject *in,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzGreySetIncValues(
				  WlzObject *obj,
				  int *val);

/************************************************************************
* WlzGreySetRange.c							*
************************************************************************/
extern WlzErrorNum		WlzGreySetRange(
				  WlzObject *obj,
				  WlzPixelV min,
				  WlzPixelV max,
				  WlzPixelV Min,
				  WlzPixelV Max,
				  int dither);

/************************************************************************
* WlzGreySetRangeLut.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzGreySetRangeLut(
				  WlzObject *obj,
				  WlzPixelV min,
				  WlzPixelV max,
				  WlzPixelP lut);
#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzGreySetValue.c							*
************************************************************************/
extern WlzErrorNum 		WlzGreySetValue(
				  WlzObject *obj,
				  WlzPixelV tmplVal);

/************************************************************************
* WlzGreyStats.c							*
************************************************************************/
#ifdef WLZ_EXT_BIND
extern int          		WlzGreySize(
				  WlzGreyType gType);
#else /* WLZ_EXT_BIND */
extern size_t          		WlzGreySize(
				  WlzGreyType gType);
#endif /* WLZ_EXT_BIND */
extern int 			WlzGreyStats(
				  WlzObject *obj,
				  WlzGreyType *dstGType,
				  double *dstMin,
				  double *dstMax,
				  double *dstSum,
				  double *dstSumSq,
				  double *dstMean,
				  double *dstStdDev,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyTemplate.c							*
************************************************************************/
extern WlzObject 		*WlzGreyTemplate(
				  WlzObject *obj,
				  WlzObject *tmpl,
				  WlzPixelV tmplVal,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyTransfer.c							*
************************************************************************/
extern WlzObject 		*WlzGreyTransfer(
				  WlzObject *obj,
				  WlzObject *tmpl,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzGreyValue.c							*
************************************************************************/
extern WlzGreyValueWSpace 	*WlzGreyValueMakeWSp(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern void 			WlzGreyValueFreeWSp(
				  WlzGreyValueWSpace *gVWSp);
extern void            		WlzGreyValueGet(
				  WlzGreyValueWSpace *gVWSp,
				  double plane,
				  double line,
				  double kol);
extern void            		WlzGreyValueGetCon(
				  WlzGreyValueWSpace *gVWSp,
				  double plane,
				  double line,
				  double kol);
extern WlzGreyType     		WlzGreyValueGetGreyType(
				  WlzGreyValueWSpace *gVWSp,
				  WlzErrorNum *dstErr);
extern int			WlzGreyValueGetI(
				  WlzGreyValueWSpace *gVWSp,
				  double plane,
				  double line,
				  double kol);
extern double			WlzGreyValueGetD(
				  WlzGreyValueWSpace *gVWSp,
				  double plane,
				  double line,
				  double kol);
extern void	                WlzGreyValueGetDir(
				  WlzGreyValueWSpace *gVWSp,
				  int plane,
				  int line,
				  int kol);


/************************************************************************
 * WlzGreyVariance.c
 ************************************************************************/
extern WlzObject       		*WlzGreyVariance(
				  WlzObject *inObj,
				  int newObjFlag,
                                  int kernelSz,
				  double scale,
				  WlzErrorNum *dstErr);


/************************************************************************
 * WlzHistogram.c
 ************************************************************************/
extern WlzObject 		*WlzHistogramObj(
				  WlzObject *srcObj,
				  int nBins,
				  double binOrigin,
				  double binSize,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzHistogramCopy(
				  WlzObject *srcHistObj,
				  WlzObjectType rtnType,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzHistogramRebin(
				  WlzObject *srcHistObj,
				  WlzObjectType rtnType,
				  int rtnMaxBins,
				  int rtnNBins,
				  double rtnOrigin,
				  double rtnBinSize,
				  WlzErrorNum *dstErr);
extern double			WlzHistogramDistance(
				  WlzObject *histObj0,
				  WlzObject *histObj1,
				  WlzErrorNum *dstErrNum);
extern double 			WlzHistogramBinSum(
				  WlzHistogramDomain *histDom);
extern int 			WlzHistogramBinMax(
				  WlzHistogramDomain *histDom);
extern WlzErrorNum 		WlzHistogramCummulative(
				  WlzObject *srcHist);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzHistogramConvolve(
				  WlzObject *histObj,
				  int krnSz,
				  double *krn);
#endif /* WLZ_EXT_BIND */
extern WlzErrorNum		WlzHistogramCnvGauss(
				  WlzObject *histObj,
				  double sigma,
				  int deriv);
extern WlzErrorNum		WlzHistogramRsvGauss(
				  WlzObject *histObj,
				  double sigma,
				  int deriv);
extern WlzErrorNum		WlzHistogramRsvFilter(
				  WlzObject *histObj,
				  WlzRsvFilter *flt);
extern WlzErrorNum 		WlzHistogramSmooth(
				  WlzObject *histObj,
				  int width);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzHistogramFindPeaks(
				  WlzObject *histObj,
				  double sigma,
				  double thresh,
				  int *dstSizeArrayPk,
				  int **dstArrayPk,
				  WlzHistFeature feat);
extern WlzErrorNum		WlzHistogramFitPeaks(
				  WlzObject *histObj,
				  int numDbn,
				  double smooth,
				  double thresh,
				  double tol,
				  int *dstSizeArrayMu,
				  double **dstArrayMu,
				  int *dstSizeArraySigma,
				  double **dstArraySigma,
				  int *dstSizeArrayAlpha,
				  double **dstArrayAlpha,
				  double *dstLL);
#endif /* WLZ_EXT_BIND */
extern WlzErrorNum 		WlzHistogramNorm(
				  WlzObject *histObj,
				  double maxVal);
extern WlzErrorNum 		WlzHistogramMapValues(
				  WlzObject *srcObj,
				  WlzObject *mapHistObj,
				  int dither);
extern WlzErrorNum 		WlzHistogramMatchObj(
				  WlzObject *srcObj,
				  WlzObject *targetHist,
				  int independentPlanes,
				  int smoothing,
				  double minDist,
				  double maxDist,
				  int dither);
extern WlzErrorNum 		WlzHistogramEqualiseObj(
				  WlzObject *srcObj,
				  int smoothing,
				  int dither);

/************************************************************************
* WlzHyThreshold.c							*
************************************************************************/
extern WlzObject		*WlzHyThreshold(
				  WlzObject *srcObj,
				  WlzPixelV pThrV,
				  WlzPixelV sThrV,
				  WlzThresholdType hilo,
				  WlzConnectType con,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzImageArithmetic.c							*
************************************************************************/
extern WlzObject 		*WlzImageArithmetic(
				  WlzObject *obj0,
				  WlzObject *obj1,
			          WlzBinaryOperatorType op,
			          int overwrite,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzIndexObj.c
************************************************************************/
extern WlzObject		*WlzIndexObjFromCompound(
				  WlzCompoundArray *cObj,
				  WlzErrorNum *dstErr);
extern WlzCompoundArray 	*WlzIndexObjToCompound(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzInsideDomain.c							*
************************************************************************/
extern int 			WlzInsideDomain(
				  WlzObject *obj,
				  double plane,
				  double line,
				  double kol,
				  WlzErrorNum *dstErr);
extern int			WlzInsideDomain2D(
    				  WlzIntervalDomain *iDom,
				  int line,
				  int kol,
				  WlzErrorNum *dstErr);
extern int			WlzInsideDomain3D(
    				  WlzPlaneDomain *pDom,
				  int plane,
				  int line,
				  int kol,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzIntersect2.c							*
************************************************************************/
extern WlzObject 		*WlzIntersect2(
				  WlzObject *obj1,
				  WlzObject *obj2,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzIntersectN.c							*
************************************************************************/
extern WlzObject 		*WlzIntersectN(
				  int sizeArrayObjs,
				  WlzObject **arrayObjs,
				  int uvt,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzHasIntersect.c							*
************************************************************************/
extern int 			WlzHasIntersection(
				  WlzObject *obj1,
				  WlzObject *obj2,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzIntervalCount.c							*
************************************************************************/
extern int			WlzIDomMaxItvLn(
				  WlzIntervalDomain *iDom);
extern int 			WlzIntervalCount(
				  WlzIntervalDomain *idom,
			    	  WlzErrorNum *dstErr);

/************************************************************************
* WlzIntervalDomScan.c							*
************************************************************************/
extern WlzErrorNum 		WlzInitRasterScan(
				  WlzObject *obj,
				  WlzIntervalWSpace *iwsp,
				  WlzRasterDir raster);
extern WlzErrorNum 		WlzNextInterval(
				  WlzIntervalWSpace *iwsp);
extern WlzErrorNum 		WlzNextLine(
				  WlzIntervalWSpace *iwsp,
			          int step);
extern WlzErrorNum 		WlzInitLineScan(
				  WlzObject *obj,
				  WlzIntervalWSpace *iwsp,
				  WlzRasterDir raster,
				  int scale,
				  int firstline);

/************************************************************************
* WlzIntRescaleObj.c							*
************************************************************************/
extern WlzObject 		*WlzIntRescaleObj(WlzObject *obj,
				  int scale,
				  int expand,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzIterate.c								*
************************************************************************/
extern WlzRasterDir		WlzRasterDir2D(
				  WlzRasterDir dir);
extern void			WlzIterateWSpFree(
				  WlzIterateWSpace *itWSp);
extern WlzIterateWSpace		*WlzIterateInit(
				  WlzObject *obj,
				  WlzRasterDir dir,
				  int grey,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzIterate(
				  WlzIterateWSpace *itWSp);

/************************************************************************
* WlzKrig.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern void			WlzKrigSetModelFn(
				  WlzKrigModelFn *fn,
				  WlzKrigModelFnType type,
				  double c0,
				  double c1,
				  double a);
extern double          		WlzKrigModelFnNugget(
				  WlzKrigModelFn *f,
				  double h);
extern double          		WlzKrigModelFnLinear(
				  WlzKrigModelFn *f,
				  double h);
extern double          		WlzKrigModelFnSpherical(
				  WlzKrigModelFn *f,
				  double h);
extern double          		WlzKrigModelFnExponential(
				  WlzKrigModelFn *f,
				  double h);
extern double          		WlzKrigModelFnGaussian(
				  WlzKrigModelFn *f,
				  double h);
extern double          		WlzKrigModelFnQuadratic(
				  WlzKrigModelFn *f,
				  double h);
extern WlzErrorNum     		WlzKrigOSetModelSV2D(
				  AlgMatrix modelSV,
				  WlzKrigModelFn *modelFn,
				  int n,
				  WlzDVertex2 *nbr,
				  int *wSp);
extern WlzErrorNum     		WlzKrigOSetModelSV3D(
				  AlgMatrix modelSV,
				  WlzKrigModelFn *modelFn,
				  int n,
				  WlzDVertex3 *nbr,
				  int *wSp);
extern WlzErrorNum     		WlzKrigOSetPosSV2D(
				  double *posSV,
				  WlzKrigModelFn *modelFn,
				  int n,
				  WlzDVertex2 *nbr,
				  WlzDVertex2 pos);
extern WlzErrorNum     		WlzKrigOSetPosSV3D(
				  double *posSV,
				  WlzKrigModelFn *modelFn,
				  int n,
				  WlzDVertex3 *nbr,
				  WlzDVertex3 pos);
extern WlzErrorNum     		WlzKrigOWeightsSolve(
				  AlgMatrix modelSV,
				  double *posSV,
				  int *wSp,
				  double eps);
extern WlzErrorNum		WlzKrigReallocBuffers2D(
				  WlzDVertex2 **dstNbrPosBuf,
				  double **dstPosSV,
				  int **dstWSp,
				  AlgMatrix *dstModelSV,
				  int *dstMaxNbrIdxBuf,
				  int nNbrC,
				  int nNbrL);
extern WlzErrorNum		WlzKrigReallocBuffers3D(
				  WlzDVertex3 **dstNbrPosBuf,
				  double **dstPosSV,
				  int **dstWSp,
				  AlgMatrix *dstModelSV,
				  int *dstMaxNbrIdxBuf,
				  int nNbrC,
				  int nNbrL);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzLBTDomain.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzLBTDomain3D           *WlzMakeLBTDomain3D(
                                  WlzObjectType type,
                                  int p1,
                                  int lp,
                                  int l1,
                                  int ll,
                                  int k1,
                                  int lk,
                                  WlzErrorNum *dstErr);
extern WlzLBTDomain2D           *WlzMakeLBTDomain2D(
                                  WlzObjectType type,
                                  int l1,
                                  int ll,
                                  int k1,
                                  int kl,
                                  WlzErrorNum *dstErr);
extern WlzDomain		WlzLBTDomainFromObj(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzLBTDomain2D           *WlzLBTDomain2DFromDomain(
                                  WlzDomain dom,
                                  WlzErrorNum *dstErr);
extern WlzLBTDomain3D		*WlzLBTDomain3DFromDomain(
				  WlzDomain dom,
				  WlzErrorNum *dstErr);
extern WlzLBTDomain2D           *WlzLBTDomain2DFromIDomain(
                                  WlzIntervalDomain *iDom,
                                  WlzErrorNum *dstErr);
extern WlzLBTDomain3D		*WlzLBTDomain3DFromPDomain(
				  WlzPlaneDomain *pDom,
				  WlzErrorNum *dstErr);
extern WlzIntervalDomain        *WlzLBTDomainToIDomain(
                                  WlzLBTDomain2D *lDom,
                                  WlzErrorNum *dstErr);
extern WlzPlaneDomain		*WlzLBTDomainToPDomain(
				  WlzLBTDomain3D *lDom,
				  WlzErrorNum *dstErr);
extern WlzIntervalDomain        *WlzIDomainFromPItv2D(
                                  int line1,
                                  int lastln,
                                  int kol1,
                                  int lastkl,
                                  int nPItv,
                                  WlzPartialItv2D *pItv,
                                  WlzErrorNum *dstErr);
extern WlzPlaneDomain           *WlzPDomainFromPItv3D(
                                  int plane1,
                                  int lastpl,
                                  int line1,
                                  int lastln,
                                  int kol1,
                                  int lastkl,
                                  int nPItv,
                                  WlzPartialItv3D *pItv,
                                  WlzErrorNum *dstErr);
extern WlzErrorNum              WlzLBTBalanceDomain2D(
                                  WlzLBTDomain2D *sDom,
                                  WlzObject *idDom,
				  int maxSz,
				  int maxBndNdSz);
extern WlzErrorNum              WlzLBTBalanceDomain3D(
				  WlzLBTDomain3D *lDom,
				  WlzObject *iObj,
				  int maxSz,
				  int maxBndSz);
extern WlzObject                *WlzLBTMakeNodeIndexObj2D(
                                  WlzLBTDomain2D *lDom,
                                  WlzIntervalDomain *iDom,
                                  WlzErrorNum *dstErr);
extern WlzObject		*WlzLBTMakeNodeIndexObj3D(
				  WlzLBTDomain3D *lDom,
				  WlzPlaneDomain *pDom,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzLBTIndexObjSetAllNodes2D(
				  WlzLBTDomain2D *lDom,
				  WlzObject *iObj);
extern WlzErrorNum		WlzLBTIndexObjSetAllNodes3D(
				  WlzLBTDomain3D *lDom,
				  WlzObject *iObj);
extern WlzErrorNum              WlzFreeLBTDomain3D(
                                  WlzLBTDomain3D *lDom);
extern WlzErrorNum              WlzFreeLBTDomain2D(
                                  WlzLBTDomain2D *lDom);
extern WlzErrorNum              WlzLBTTestOutputNodesTxt(
                                  FILE *fP,
                                  WlzDomain dom);
extern WlzErrorNum              WlzLBTTestOutputNodesVtk(
                                  FILE *fP,
                                  WlzDomain dom);
extern int                      WlzLBTNodeSz2D(
                                  WlzLBTNode2D *nod);
extern int			WlzLBTNodeSz3D(
				  WlzLBTNode3D *nod);
extern int                    	WlzLBTNodeLogSz3D(
                                  WlzLBTNode3D *nod);
extern int                    	WlzLBTNodeLogSz2D(
                                  WlzLBTNode2D *nod);
extern int			WlzLBTCountNodNbrDir2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzDirection dir);
extern void                     WlzLBTPosToKey3D(
                                  WlzIVertex3 pos,
                                  unsigned *keys);
extern void                     WlzLBTPosToKey2D(
                                  WlzIVertex2 pos,
                                  unsigned *keys);
extern void                     WlzLBTGetKeyDigits3D(
                                  unsigned *keys,
                                  WlzUByte *digits);
extern void                     WlzLBTGetKeyDigits2D(
                                  unsigned *keys,
                                  WlzUByte *digits);
extern void                     WlzLBTKeyToPos3I(
                                  unsigned *key,
                                  WlzIVertex3 *pos);
extern void                     WlzLBTKeyToPos2I(
                                  unsigned *key,
                                  WlzIVertex2 *pos);
extern void                     WlzLBTKeyToBox3I(
                                  unsigned *key,
                                  WlzIBox3 *box);
extern void                     WlzLBTKeyToBox2I(
                                  unsigned *key,
                                  WlzIBox2 *box);
extern void			WlzLBTClassifyNode2D(
				  WlzLBTDomain2D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  WlzLBTNodeClass2D *cls,
				  int *dstRot);
extern void			WlzLBTClassifyNodeFace3D(WlzLBTDomain3D *lDom,
				  WlzGreyValueWSpace *iGVWSp,
				  int idN,
				  int idF,
				  WlzDVertex3 *vtx,
				  WlzLBTNodeClass2D *dstCls,
				  int *dstRot);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* "THE BEAST" - WlzLabel.c						*
************************************************************************/
extern WlzErrorNum 		WlzLabel(
				  WlzObject *obj,
				  int *dstArraySizeObjs,
				  WlzObject ***dstArrayObjs,
				  int maxNumObjs,
				  int ignlns,
				  WlzConnectType connect);

/************************************************************************
* WlzLaplacian.c							*
************************************************************************/
extern WlzObject       		*WlzLaplacian(
				  WlzObject *srcObj,
				  int kSize,
				  int mkNewObjFlg,
				  int modFlg,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzLineArea.c								*
************************************************************************/
extern int			WlzLineArea(
				  WlzObject *obj,
		       		  WlzErrorNum *dstErr);

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzLUT.c								*
************************************************************************/
extern WlzObject       		*WlzMakeLUTObject(
				  WlzGreyType vType,
				  int bin1,
				  int lastbin,
				  WlzErrorNum *dstErr);
extern WlzLUTDomain  		*WlzMakeLUTDomain(
				  int bin1,
				  int lastbin,
				  WlzErrorNum *dstErr);
extern WlzLUTValues  		*WlzMakeLUTValues(
				  WlzGreyType vType,
				  int maxVal,
				  WlzErrorNum *dstErr);
extern WlzErrorNum     		WlzFreeLUTValues(
				  WlzValues val);
extern WlzObject       		*WlzLUTGreyTransformNew(
				  WlzGreyTransformType gTrType,
				  WlzGreyType gType,
				  int il,
				  int iu,
				  WlzGreyV gol,
				  WlzGreyV gou,
				  double p0,
				  double p1,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzLUTGreyTransformSet(
				  WlzObject *obj,
				  WlzGreyTransformType gTrType,
				  WlzGreyType gType,
				  int il,
				  int iu,
				  WlzGreyV gol,
				  WlzGreyV gou,
				  double p0,
				  double p1);
extern WlzObject       		*WlzLUTMergeToRGBA(
				  WlzObject *r,
				  WlzObject *g,
				  WlzObject *b,
				  WlzObject *a,
				  WlzErrorNum *dstErr);
extern WlzObject       		*WlzLUTTransformObj(
				  WlzObject *vObj,
				  WlzObject *tObj,
				  WlzGreyType rGType,
				  int inplace,
				  int dither,
				  WlzErrorNum *dstErr);
extern WlzPixelV		WlzLUTTransformPixelValue(
				  WlzObject *tObj,
				  WlzPixelV gPix,
				  WlzGreyType rType,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzLUTTransformGreyValues(
				  WlzObject *tObj,
				  WlzGreyP rP,
				  WlzGreyType rType,
				  WlzGreyP gP,
				  WlzGreyType gType,
				  int nVal,
				  int dither);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzMakeAffineTransform.c						*
************************************************************************/
extern WlzAffineTransform 	*WlzMakeAffineTransform(
				  WlzTransformType type,
				  WlzErrorNum *dstErr);
extern WlzErrorNum    		WlzFreeAffineTransform(
				  WlzAffineTransform *trans);

/************************************************************************
* WlzMakeCompound.c							*
************************************************************************/
extern WlzCompoundArray 	*WlzMakeCompoundArray(
				  WlzObjectType type,
				  int mode,
				  int sizeArrayObjs,
				  WlzObject **arrayObjs,
				  WlzObjectType otype,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzMakeIntervalValues.c						*
************************************************************************/
extern WlzIntervalValues 	*WlzMakeIntervalValues(
				  WlzObjectType type,
				  WlzObject *obj,
				  WlzPixelV bckgrnd,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzMakeProperties.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzPropertyList		*WlzMakePropertyList(
				  WlzErrorNum *dstErr);
extern WlzSimpleProperty	*WlzMakeSimpleProperty(
				  int size,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzFreeSimpleProperty(
				  WlzSimpleProperty *prop);
extern WlzEMAPProperty          *WlzMakeEMAPProperty(
                                  WlzEMAPPropertyType type,
				  char *modelUID,
				  char *anatomyUID,
				  char *targetUID,
				  char *targetVersion,
				  char *stage,
				  char *subStage,
				  char *modelName,
				  char *version,
				  char *fileName,
				  char *comment,
				  WlzErrorNum *dstErr);
extern WlzNameProperty 		*WlzMakeNameProperty(
				  char *name,
				  WlzErrorNum *dstErr);
extern WlzGreyProperty 		*WlzMakeGreyProperty(
				  char *name,
				  WlzPixelV val,
				  WlzErrorNum *dstErr);
extern WlzTextProperty 		*WlzMakeTextProperty(
				  char *name,
				  char *text,
				  WlzErrorNum *dstErr);
extern WlzErrorNum              WlzChangeEMAPProperty(
                                  WlzEMAPProperty *prop,
                                  WlzEMAPPropertyType type,
				  char *modelUID,
				  char *anatomyUID,
				  char *targetUID,
				  char *targetVersion,
				  char *stage,
				  char *subStage,
				  char *modelName,
				  char *version,
				  char *fileName,
				  char *comment);
extern WlzErrorNum		WlzFreeEMAPProperty(
				  WlzEMAPProperty *prop);
extern WlzErrorNum		WlzFreeProperty(
				  WlzProperty prop);
extern WlzErrorNum		WlzFreePropertyList(
				  WlzPropertyList *plist);
extern void			WlzFreePropertyListEntry(
				  void *prop);
extern WlzProperty 		WlzGetProperty(
                                  AlcDLPList *plist,
				  WlzObjectType type,
				  WlzErrorNum *dstErr);
extern WlzErrorNum              WlzRemoveProperty(
                                  AlcDLPList *plist,
				  WlzProperty prop);
#endif /* WLZ_EXT_BIND */

/************************************************************************
 * WlzMakeStructs.c
 ************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject       		*WlzMarkerLattice(
			 	  WlzObject *gObj,
				  WlzMarkerType mType,
                                  int  mSz,
				  int mSep,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzMakeMarkers(
				  WlzVertexType vType,
				  int nVtx,
				  WlzVertexP vtx,
				  WlzMarkerType mType,
				  int mSz,
				  WlzErrorNum *dstErr);
extern WlzPoints		*WlzMakePoints(
    				  WlzObjectType type,
				  int nVtx,
				  WlzVertexP vtxP,
				  int maxVtx,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzPolygonDomain		*WlzMakePolygonDomain(
				  WlzObjectType type,
				  int sizeArrayVertices,
				  WlzIVertex2 *arrayVertices,
				  int maxVertices,
				  int copy,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzObject 		*WlzMakeMain(
				  WlzObjectType type,
				  WlzDomain domain,
				  WlzValues values,
				  WlzPropertyList *prop,
				  WlzObject *assoc,
				  WlzErrorNum *dstErr);
extern WlzIntervalDomain	*WlzMakeIntervalDomain(
				  WlzObjectType type,
				  int l1,
				  int ll,
				  int k1,
				  int kl,
				  WlzErrorNum *dstErr);

extern WlzErrorNum		WlzMakeInterval(
				  int l,
				  WlzIntervalDomain *idom,
				  int nints,
				  WlzInterval *itvl); 
extern WlzErrorNum		WlzMakeValueLine(
				  WlzRagRValues *vtb,
				  int line,
				  int k1,
				  int kl,
				  int *greyptr); 
extern WlzPlaneDomain		*WlzMakePlaneDomain(
				  WlzObjectType type,
				  int p1,
				  int pl,
				  int l1,
				  int ll,
				  int k1,
				  int kl,
				  WlzErrorNum *dstErr);
extern WlzRagRValues 		*WlzMakeValueTb(
				  WlzObjectType type,
				  int l1,
				  int ll,
				  int k1,
				  WlzPixelV bckgrnd,
				  WlzObject *orig,
				  WlzErrorNum *dstErr);
extern WlzRectValues		*WlzMakeRectValueTb(
				  WlzObjectType type,
				  int l1,
				  int ll,
				  int k1,
				  int width,
				  WlzPixelV bckgrnd,
				  int *greyptr,
				  WlzErrorNum *dstErr);
extern WlzVoxelValues		*WlzMakeVoxelValueTb(
				  WlzObjectType type,
				  int p1,
				  int pl,
				  WlzPixelV bckgrnd,
				  WlzObject *original,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzNewObjectValues(
				  WlzObject *sObj,
				  WlzObjectType tType,
				  WlzPixelV bgdV,
				  int setFG,
				  WlzPixelV fgdV,
				  WlzErrorNum *dstErr);
extern WlzVoxelValues		*WlzNewValuesVox(
				  WlzObject *sObj,
				  WlzObjectType gTType,
				  WlzPixelV bgdV,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzMakeRect(
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  WlzGreyType pixeltype,
				  int *greyptr,
				  WlzPixelV bckgrnd,
				  WlzPropertyList *plist,
				  WlzObject *assoc_obj,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzMakeEmpty(
				  WlzErrorNum *dstErr);
extern WlzBoundList		*WlzMakeBoundList(
				  WlzObjectType type,
				  int wrap,
				  WlzPolygonDomain *poly,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzMakeCuboid(
				  int plane1,
				  int lastpl,
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  WlzGreyType pixType,
				  WlzPixelV bgdV,
				  WlzPropertyList *plist,
				  WlzObject *assocObj,
				  WlzErrorNum *dstErr);
extern WlzHistogramDomain	*WlzMakeHistogramDomain(
				  WlzObjectType type,
				  int maxBins,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzMakeHistogram(
				  WlzObjectType type,
				  int maxBins,
				  WlzErrorNum *dstErr);
extern WlzIVertex2		*WlzMakeIVertex(
				  int nverts,
				  WlzErrorNum *dstErr);
extern WlzObject              	*WlzNewGrey(
				  WlzObject *iobj,
				  WlzErrorNum *dstErr);
extern WlzRagRValues		*WlzNewValueTb(
				  WlzObject *obj,
				  WlzObjectType type,
				  WlzPixelV bckgrnd,
				  WlzErrorNum *dstErr);
extern WlzIntervalDomain	*WlzNewIDomain(
				  WlzObjectType type,
				  WlzIntervalDomain *idom,
				  WlzErrorNum *dstErr);
extern WlzContour		*WlzMakeContour(
				  WlzErrorNum *dstErr);
extern WlzIndexedValues		*WlzMakeIndexedValues(
				  WlzObject *dObj,
				  int rank,
				  int *dim,
				  WlzGreyType vType,
				  WlzValueAttach attach,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzExtendIndexedValues(
				  WlzIndexedValues *ixv,
				  size_t n);
extern WlzErrorNum		WlzIndexedValuesSet(
				  WlzObject *obj,
				  size_t cnt,
				  void *val);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzMatchICP.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzMatchICPObjs(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  int *dstNMatch,
				  WlzVertexP *dstTMatch,
				  WlzVertexP *dstSMatch,
				  int maxItr,
				  int minSpx,
				  int minSegSpx,
				  int brkFlg,
				  double maxDisp,
				  double maxAng,
				  double maxDeform,
				  int matchImpNN,
				  double matchImpThr,
				  double delta);
extern WlzErrorNum  		WlzMatchICPCtr(
				  WlzContour *tCtr,
				  WlzContour *sCtr,
                                  WlzAffineTransform *initTr,
                                  int maxItr,
				  int minSpx,
				  int minSegSpx,
                                  int *dstNMatch,
				  WlzVertexP *dstTMatch,
                                  WlzVertexP *dstSMatch,
				  int brkFlg,
                                  double  maxDisp,
				  double maxAng,
                                  double maxDeform,
                                  int matchImpNN,
				  double matchImpThr,
                                  WlzRegICPUsrWgtFn usrWgtFn,
                                  void *usrWgtData,
				  double delta);
extern double          		WlzMatchICPWeightMatches(
				  WlzVertexType vType,
				  WlzAffineTransform *curTr,
				  AlcKDTTree *tree,
				  WlzVertexP tVx,
				  WlzVertexP sVx,
				  WlzVertex tMVx,
				  WlzVertex sMVx,
				  double wVx,
				  double wNr,
				  void *data);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzMeshGen.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzCMeshElmFuse2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshElm2D *gElm);
extern WlzObject		*WlzCMeshComputeNormalsIxv2D5(
				  WlzObject *gObj,
				  int nrmFlg,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzCMeshExtract2D(
				  WlzObject *gObj,
				  int applyDsp,
				  WlzErrorNum *dstErr);
extern WlzCMesh2D               *WlzCMeshNew2D(
                                  WlzErrorNum *dstErr);
extern WlzCMesh2D5               *WlzCMeshNew2D5(
                                  WlzErrorNum *dstErr);
extern WlzCMesh3D               *WlzCMeshNew3D(
                                  WlzErrorNum *dstErr);
extern WlzCMesh2D               *WlzCMeshFromBalLBTDom2D(
                                  WlzLBTDomain2D *lDom,
                                  WlzObject *iObj,
                                  WlzErrorNum *dstErr);
extern WlzCMesh3D               *WlzCMeshFromBalLBTDom3D(
                                  WlzLBTDomain3D *lDom,
                                  WlzObject *iObj,
                                  WlzErrorNum *dstErr);
extern WlzCMeshP		WlzCMeshFromObj(
				  WlzObject *obj,
				  double minElmSz,
				  double maxElmSz,
				  WlzObject **dstDilObj,
				  int conform,
				  WlzErrorNum *dstErr);
extern WlzCMesh2D    		*WlzCMeshFromObj2D(
                                  WlzObject *obj,
                                  double minElmSz,
                                  double maxElmSz,
				  WlzObject **dstDilObj,
				  int conform,
                                  WlzErrorNum *dstErr);
extern WlzCMesh3D               *WlzCMeshFromObj3D(
                                  WlzObject *obj,
				  double minElmSz,
				  double maxElmSz,
				  WlzObject **dstDilObj,
				  int conform,
                                  WlzErrorNum *dstErr);
extern WlzCMesh2D5		*WlzCMeshFromGM(WlzGMModel *model,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzCMeshSetElm2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshElm2D *elm,
				  WlzCMeshNod2D *nod0,
				  WlzCMeshNod2D *nod1,
				  WlzCMeshNod2D *nod2,
				  int allowFlip);
extern WlzErrorNum		WlzCMeshSetElm2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshElm2D5 *elm,
				  WlzCMeshNod2D5 *nod0,
				  WlzCMeshNod2D5 *nod1,
				  WlzCMeshNod2D5 *nod2,
				  int allowFlip);
extern WlzErrorNum		WlzCMeshSetElm3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshElm3D *elm,
				  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  int allowFlip);
extern WlzErrorNum		WlzCMeshAddNewNodCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddNewNodCb2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddNewNodCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddNewElmCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddNewElmCb2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddNewElmCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelNodCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelNodCb2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelNodCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelEdgCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelEdgCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelElmCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelElmCb2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshAddDelElmCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemNewNodCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemNewNodCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemNewElmCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemNewElmCb2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemNewElmCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelNodCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelNodCb2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelNodCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelEdgCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelEdgCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelElmCb2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelElmCb2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzErrorNum		WlzCMeshRemDelElmCb3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshCbFn fn,
				  void *data);
extern WlzCMeshNod2D            *WlzCMeshNewNod2D(
                                  WlzCMesh2D *mesh,
                                  WlzDVertex2 pos,
                                  WlzErrorNum *dstErr);
extern WlzCMeshNod2D5            *WlzCMeshNewNod2D5(
                                  WlzCMesh2D5 *mesh,
                                  WlzDVertex3 pos,
                                  WlzErrorNum *dstErr);
extern WlzCMeshNod3D            *WlzCMeshNewNod3D(
                                  WlzCMesh3D *mesh,
                                  WlzDVertex3 pos,
                                  WlzErrorNum *dstErr);
extern WlzCMeshNod2D		*WlzCMeshAllocNod2D(
				  WlzCMesh2D *mesh);
extern WlzCMeshNod2D5		*WlzCMeshAllocNod2D5(
				  WlzCMesh2D5 *mesh);
extern WlzCMeshNod3D		*WlzCMeshAllocNod3D(
				  WlzCMesh3D *mesh);
extern WlzCMeshNod2D            *WlzCMeshMatchNod2D(
                                  WlzCMesh2D *mesh,
                                  WlzDVertex2 pos);
extern WlzCMeshNod2D5           *WlzCMeshMatchNod2D5(
                                  WlzCMesh2D5 *mesh,
                                  WlzDVertex3 pos);
extern WlzCMeshElm2D		*WlzCMeshAllocElm2D(
				  WlzCMesh2D *mesh);
extern WlzCMeshElm2D5		*WlzCMeshAllocElm2D5(
				  WlzCMesh2D5 *mesh);
extern WlzCMeshElm3D		*WlzCMeshAllocElm3D(
				  WlzCMesh3D *mesh);
extern WlzCMeshElm2D            *WlzCMeshNewElm2D(
                                  WlzCMesh2D *mesh,
                                  WlzCMeshNod2D *nod0,
                                  WlzCMeshNod2D *nod1,
                                  WlzCMeshNod2D *nod2,
				  int allowFlip,
                                  WlzErrorNum *dstErr);
extern WlzCMeshElm2D5 		*WlzCMeshNewElm2D5(
				  WlzCMesh2D5 *mesh,
				   WlzCMeshNod2D5 *nod0,
				   WlzCMeshNod2D5 *nod1,
				   WlzCMeshNod2D5 *nod2,
				  int allowFlip,
				   WlzErrorNum *dstErr);
extern WlzCMeshElm3D            *WlzCMeshNewElm3D(
                                  WlzCMesh3D *mesh,
                                  WlzCMeshNod3D *nod0,
                                  WlzCMeshNod3D *nod1,
                                  WlzCMeshNod3D *nod2,
                                  WlzCMeshNod3D *nod3,
				  int allowFlip,
                                  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzCMeshFree(
				  WlzCMeshP mesh);
extern WlzErrorNum              WlzCMeshFree2D(
                                  WlzCMesh2D *mesh);
extern WlzErrorNum              WlzCMeshFree2D5(
                                  WlzCMesh2D5 *mesh);
extern WlzErrorNum              WlzCMeshFree3D(
                                  WlzCMesh3D *mesh);
extern WlzErrorNum		WlzCMeshDelElm2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshElm2D *elm);
extern WlzErrorNum		WlzCMeshDelElm3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshElm3D *elm);
extern WlzErrorNum		WlzCMeshDelNod2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshNod2D *nod);
extern WlzErrorNum		WlzCMeshDelNod2D5(
				  WlzCMesh2D5 *mesh,
				  WlzCMeshNod2D5 *nod);
extern WlzErrorNum		WlzCMeshBoundConform2D(
				  WlzCMesh2D *mesh,
				  WlzObject *obj,
				  int method,
				  double tol);
extern WlzErrorNum		WlzCMeshBoundConform3D(
				  WlzCMesh3D *mesh,
				  WlzObject *obj,
				  int method,
				  double tol);
extern double			WlzCMeshElmMinEdgLnSq2D(
				  WlzCMeshElm2D *elm);
extern double			WlzCMeshElmMinEdgLnSq3D(
				  WlzCMeshElm3D *elm);
extern WlzErrorNum              WlzCMeshAffineTransformMesh2D(
                                  WlzCMesh2D *mesh,
                                  WlzAffineTransform *tr);
extern WlzErrorNum              WlzCMeshAffineTransformMesh2D5(
                                  WlzCMesh2D5 *mesh,
                                  WlzAffineTransform *tr);
extern WlzErrorNum              WlzCMeshAffineTransformMesh3D(
                                  WlzCMesh3D *mesh,
                                  WlzAffineTransform *tr);
extern WlzErrorNum		WlzCMeshReassignGridCells2D(
				  WlzCMesh2D *mesh,
				  int newNumNod);
extern WlzErrorNum		WlzCMeshReassignGridCells2D5(
				  WlzCMesh2D5 *mesh,
				  int newNumNod);
extern WlzErrorNum		WlzCMeshReassignGridCells3D(
				  WlzCMesh3D *mesh,
				  int newNumNod);
extern void                     WlzCMeshNodFree2D(
                                  WlzCMesh2D *mesh,
                                  WlzCMeshNod2D *nod);
extern void                     WlzCMeshNodFree2D5(
                                  WlzCMesh2D5 *mesh,
                                  WlzCMeshNod2D5 *nod);
extern void                     WlzCMeshNodFree3D(
                                  WlzCMesh3D *mesh,
                                  WlzCMeshNod3D *nod);
extern void                     WlzCMeshElmFree2D(
                                  WlzCMesh2D *mesh,
                                  WlzCMeshElm2D *elm);
extern void                     WlzCMeshElmFree2D5(
                                  WlzCMesh2D5 *mesh,
                                  WlzCMeshElm2D5 *elm);
extern void                     WlzCMeshElmFree3D(
                                  WlzCMesh3D *mesh,
                                  WlzCMeshElm3D *elm);
extern void			WlzCMeshDelUnusedNodes3D(
				  WlzCMesh3D *mesh);
extern void			WlzCMeshDelUnusedNodes2D5(
				  WlzCMesh2D5 *mesh);
extern void			WlzCMeshDelUnusedNodes2D(
				  WlzCMesh2D *mesh);
extern int                      WlzCMeshMatchNNod2D(
                                  WlzCMesh2D *mesh,
                                  int nNod,
                                  WlzDVertex2 *nPos,
                                  WlzCMeshNod2D **mNod);
extern int			WlzCMeshMatchNNod2D5(
				  WlzCMesh2D5 *mesh,
				  int nNod,
				  WlzDVertex3 *nPos,
				  WlzCMeshNod2D5 **mNod);
extern int                      WlzCMeshMatchNNod3D(
                                  WlzCMesh3D *mesh,
                                  int nNod,
                                  WlzDVertex3 *nPos,
                                  WlzCMeshNod3D **mNod);
extern int			WlzCMeshMatchNNodIdx2D(
				  WlzCMesh2D *mesh,
				  int nNod,
				  WlzDVertex2 *nPos,
				  int *mIdx);
extern int			WlzCMeshMatchNNodIdx2D5(
				  WlzCMesh2D5 *mesh,
				  int nNod,
				  WlzDVertex3 *nPos,
				  int *mIdx);
extern int			WlzCMeshMatchNNodIdx3D(
				  WlzCMesh3D *mesh,
				  int nNod,
				  WlzDVertex3 *nPos,
				  int *mIdx);
extern int                      WlzCMeshLocateNod2D(
                                  WlzCMesh2D *mesh,
                                  WlzDVertex2 nPos,
                                  WlzCMeshNod2D **dstNod);
extern int                      WlzCMeshLocateNod2D5(
                                  WlzCMesh2D5 *mesh,
                                  WlzDVertex3 nPos,
                                  WlzCMeshNod2D5 **dstNod);
extern int                      WlzCMeshLocateNod3D(
                                  WlzCMesh3D *mesh,
                                  WlzDVertex3 nPos,
                                  WlzIVertex3 *dstGPos,
                                  WlzCMeshNod3D **dstPrev,
                                  WlzCMeshNod3D **dstNod);
extern int             		WlzCMeshClosestNod2D(
				  WlzCMesh2D *mesh,
				  WlzDVertex2 pos);
extern int             		WlzCMeshClosestNod2D5(
				  WlzCMesh2D5 *mesh,
				  WlzDVertex3 pos);
extern int             		WlzCMeshClosestNod3D(
				  WlzCMesh3D *mesh,
				  WlzDVertex3 pos);
extern int			WlzCMeshElmEnclosingPos(
				  WlzCMeshP mesh,
				  int lastElmIdx,
				  double pX,
				  double pY,
				  double pZ,
				  int exhaustive,
				  int *dstCloseNod);
extern int			WlzCMeshElmEnclosingPos2D(
				  WlzCMesh2D *mesh,
				  int lastElmIdx,
				  double pX,
				  double pY,
				  int exhaustive,
				  int *dstCloseNod);
extern int			WlzCMeshElmEnclosingPos2D5(
				  WlzCMesh2D5 *mesh,
				  int lastElmIdx,
				  double pX,
				  double pY,
				  double pZ,
				  int exhaustive,
				  int *dstCloseNod);
extern int			WlzCMeshElmEnclosingPos3D(
				  WlzCMesh3D *mesh,
				  int lastElmIdx,
				  double pX,
				  double pY,
				  double pZ,
				  int exhaustive,
				  int *dstCloseNod);
extern int			WlzCMeshElmEnclosesPos2D(
				  WlzCMeshElm2D *elm,
				  WlzDVertex2 gPos);
extern int			WlzCMeshElmEnclosesPos2D5(
				  WlzCMeshElm2D5 *elm,
				  WlzDVertex3 gPos);
extern int			WlzCMeshElmEnclosesPos3D(
				  WlzCMeshElm3D *elm,
				  WlzDVertex3 gPos);
extern void			WlzCMeshDbgOutVTK(
				  FILE *fP,
				  WlzCMeshP mesh);
extern void			WlzCMeshDbgOutVTK2D(FILE *fP,
				  WlzCMesh2D *mesh);
extern void			WlzCMeshDbgOutVTK3D(FILE *fP,
				  WlzCMesh3D *mesh);
extern int			WlzCMeshElmClosestPosIn(
				  WlzCMeshP mesh,
				  WlzVertexP dstPos,
				  WlzVertex pos,
				  double dMax);
extern int			WlzCMeshElmClosestPosIn2D(
				  WlzCMesh2D *mesh,
                                  WlzDVertex2 *dstPos,
				  WlzDVertex2 pos,
				  double dMax);
extern int			WlzCMeshElmClosestPosIn2D5(
				  WlzCMesh2D5 *mesh,
                                  WlzDVertex3 *dstPos,
				  WlzDVertex3 pos,
				  double dMax);
extern int			WlzCMeshElmClosestPosIn3D(
				  WlzCMesh3D *mesh,
                                  WlzDVertex3 *dstPos,
				  WlzDVertex3 pos,
				  double dMax);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzMeshTransform.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzMeshFreeTransform(
				  WlzMeshTransform *mesh);
extern WlzMeshTransform 	*WlzMeshFromObj(
				  WlzObject *inObj,
				  WlzMeshGenMethod method,
				  double minDist,
				  double maxDist,
				  WlzErrorNum *dstErr);
extern WlzMeshTransform 	*WlzMeshTransformNew(
				  unsigned int maxElem,
				  unsigned int maxNode,
				  WlzErrorNum *dstErr);
extern WlzMeshTransform 	*WlzMeshTransformAdapt(
				  WlzMeshTransform *mesh,
				  double minArea,
				  WlzErrorNum *dstErr);
extern WlzMeshTransform 	*WlzMeshTransformCopy(
				  WlzMeshTransform *mesh,
			          WlzErrorNum *dstErr);
extern WlzErrorNum		WlzMeshTransformVerify(
				  WlzMeshTransform *mesh,
				  int dispFlg,
				  int *badElm,
				  WlzMeshError *dstErrMsk);
extern WlzObject		*WlzMeshTransformObj(
				  WlzObject *srcObj,
				  WlzMeshTransform *mesh,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
extern WlzMeshTransform 	*WlzMeshTransformFromCPts(
				  WlzObject *obj,
				  WlzFnType basisFnType,
				  int polyOrder,
				  int nDPts,
				  WlzDVertex2 *dPts,
				  int nSPts,
				  WlzDVertex2 *sPts,
				  WlzMeshGenMethod meshGenMtd,
				  double meshMinDist,
				  double meshMaxDist,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzMeshAffineProduct(
				  WlzMeshTransform *mTr,
				  WlzAffineTransform *aTr,
				  int order);
extern WlzDVertex2 		WlzMeshTransformVtx(
				  WlzDVertex2 vtx,
				  WlzMeshTransform *mesh,
				  WlzErrorNum *dstErr);
extern double 			WlzClassValCon4(
                                  double *gVals,
				  double xOffset,
				  double yOffset);
extern WlzMeshTransform  	*WlzMeshFromObjBox(
				  WlzObject *srcObj,
				  WlzIBox2 *dstBox,
				  int boxDilation,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzMeshUtils.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern int			WlzMeshElemNodeIdxFromVx(
				  WlzMeshTransform *mesh,
				  WlzMeshElem *elm,
				  WlzDVertex2 gVx);
extern int			WlzMeshElemNodeIdxFromNodeIdx(
				  int *nodes,
				  int mNodId);
extern int			WlzMeshElemNbrIdxFromNodes(
				  WlzMeshElem *elm,
				  int nodId0,
				  int nodId1);
extern int			WlzMeshElemFindVx(
				  WlzMeshTransform *mesh,
				  WlzDVertex2 gVx,
				  int startElm,
				  int *lastElm,
				  int *existsFlg,
				  WlzErrorNum *dstErr);
extern int             		WlzMeshClosestNod2D(
				  WlzMeshTransform *mesh,
				  WlzDVertex2 pos);
extern double          		WlzMeshMaxEdgeLenSq(
				  WlzMeshTransform *mesh,
                                  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzMeshDomainAdd(
				  WlzMeshTransform *mesh,
				  WlzObject *obj,
				  double minDist,
				  WlzDVertex2 scale);
extern WlzErrorNum		WlzMeshVxVecAdd(
				  WlzMeshTransform *mesh,
				  WlzDVertex2 *vxVec,
				  int nVx,
				  double minDistSq,
				  unsigned int nodeFlgs);
extern WlzErrorNum		WlzMeshExpand(
				  WlzMeshTransform *mesh,
				  int nElem,
				  int nNodes);
extern WlzErrorNum		WlzMeshElemVerify(
				  WlzMeshTransform *mesh,
				  int dispFlg,
				  WlzMeshElem *elm,
				  WlzMeshError *dstErrMsk);
extern WlzErrorNum		WlzMeshElemSplit(
				  WlzMeshTransform *mesh,
				  int sElmIdx);
extern WlzErrorNum		WlzMeshNodeAdd(
				  WlzMeshTransform *mesh,
				  int startElm,
				  WlzDVertex2 gVx,
				  unsigned int nodeFlgs);
extern WlzErrorNum		WlzMeshNodeDelIdx(
				  WlzMeshTransform *mesh,
				  int startElm,
				  int *nodIdx,
				  int nNod);
extern WlzErrorNum		WlzMeshSqueeze(
				  WlzMeshTransform *mesh);
extern WlzErrorNum		WlzMeshGetNodesAndEdges(
				  WlzMeshTransform *mesh,
				  int *dstSizeArrayNod,
				  WlzDVertex2 **dstArrayNod,
				  int *dstSizeArrayDsp,
				  WlzDVertex2 **dstArrayDsp,
				  int *dstSizeArrayEdg,
				  int **dstArrayEdg);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzNMSuppress.c							*
************************************************************************/
extern WlzObject		*WlzNMSuppress(
				  WlzObject *grdM,
				  WlzObject *grdZ,
				  WlzObject *grdY,
				  WlzObject *grdX,
			          WlzPixelV minThrV,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzObjToBoundary.c							*
************************************************************************/
extern WlzObject 		*WlzObjToBoundary(
				  WlzObject *obj,
				  int wrap,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzOccupancy.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern	WlzErrorNum		Wlz3DSectionOcc(
				  WlzObject *obj,
				  WlzThreeDViewStruct *vs,
				  double sep,
				  double *dstFirst,
				  double *dstLast,
				  int *dstArraySizeOcc,
				  int **dstArrayOcc);
#endif /* WLZ_EXT_BIND */
extern WlzObject		*WlzDomainOccupancy(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzPoints.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzPointsToDomObj(
    				  WlzPoints *pnt,
				  double scale,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzPolarSample.c							*
************************************************************************/
extern WlzObject 		*WlzPolarSample(
				  WlzObject *srcObj,
				  WlzIVertex2 org,
				  double angleInc,
				  double distInc,
				  int nLines,
				  int outFlg,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzPolyToObj.c							*
************************************************************************/
extern WlzObject 		*WlzPolyToObj(
				  WlzPolygonDomain *pgdm,
				  WlzPolyFillMode fillMode,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzPolygonToObj(
				  WlzObject *polygon,
				  WlzPolyFillMode fillMode,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzPolyTo8Polygon(
				  WlzPolygonDomain *pgdm,
				  int wrap,
				  WlzErrorNum *dstErr);
extern WlzObject                *WlzBoundTo8Bound(
                                  WlzObject *gvnObj,
                                  WlzErrorNum *dstErr);
extern int			WlzPolyCrossings(
				  WlzIVertex2 vtx,
				  WlzPolygonDomain *pgdm,
				  WlzErrorNum *dstErr);
extern int			WlzInsidePolyEO(
				  WlzIVertex2 vtx,
				  WlzPolygonDomain *pgdm,
				  WlzErrorNum *dstErr);
extern int			WlzPolyCrossingsD(
				  WlzDVertex2 vtx,
				  WlzPolygonDomain *pgdm,
				  WlzErrorNum *dstErr);
extern int			WlzInsidePolyEOD(
				  WlzDVertex2 vtx,
				  WlzPolygonDomain *pgdm,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzPolyEquispace.c							*
************************************************************************/
extern WlzPolygonDomain		*WlzPolyEquispace(
                                  WlzPolygonDomain *poly,
				  int wrap,
				  double spacing,
				  int keepOrigVtxs,
				  WlzErrorNum *dstErr);
extern double 			WlzPolyLength(
                                  WlzPolygonDomain *poly,
				  int wrap,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzPolyDecimate.c							*
************************************************************************/
extern WlzPolygonDomain  	*WlzPolyDecimate(
                                  WlzPolygonDomain *poly,
				  int wrap,
				  double maxDist,
				  WlzErrorNum *dstErr);
extern WlzBoundList 		*WlzBoundDecimate(
                                  WlzBoundList *bound,
				  double maxDist,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzPolyReverse.c							*
************************************************************************/
extern WlzPolygonDomain 	*WlzPolyReverse(
                                  WlzPolygonDomain *poly,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzPolySmooth.c							*
************************************************************************/
extern WlzPolygonDomain		*WlzPolySmooth(
                                  WlzPolygonDomain *poly,
				  int wrap,
				  int iterations,
				  WlzErrorNum *dstErr);
extern WlzBoundList 		*WlzBoundSmooth(
                                  WlzBoundList *bound,
				  int iterations,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzPrinicipalAngle.c
************************************************************************/
extern double 			WlzPrincipalAngle(
				  WlzObject *srcObj,
				  WlzDVertex2 cMass,
				  int binObjFlg,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzProj3DToSection.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject       		*WlzProj3DToSection(
				  WlzObject *gvnObj,
				  int nMask,
				  WlzObject **maskObj,
				  WlzThreeDViewStruct *prjView,
				  int nPlnTr,
				  WlzObject **plnTrObj,
				  WlzThreeDViewStruct *secView,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzRank.c
************************************************************************/
extern WlzErrorNum		WlzRankFilter(
				  WlzObject *gObj,
				  int fSz,
				  double rank);

/************************************************************************
* WlzRaster.c
************************************************************************/
extern WlzObject		*WlzRasterObj(
				  WlzObject *gObj,
			          WlzErrorNum *dstErr);
extern void            		WlzRasterLineSetItv2D(
				  WlzIntervalDomain *iDom,
                                  WlzIVertex2 v0,
				  WlzIVertex2 v1);

/************************************************************************
* WlzReadObj.c								*
************************************************************************/
extern WlzObjectType		WlzReadObjType(
				  FILE *fp,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzReadObj(
				  FILE *fP,
			          WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzMeshTransform3D 	*WlzReadMeshTransform3D(
				  FILE *fP,
				  WlzErrorNum *dstErr);
#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzRegCCor.c
************************************************************************/
extern WlzAffineTransform 	*WlzRegCCorObjs(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  WlzDVertex2 maxTran,
				  double maxRot,
				  int maxItr,
				  WlzWindowFnType winFn,
				  int noise,
				  int inv,
				  int *dstConv,
				  double *dstCCor,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRegICP.c
************************************************************************/
extern WlzObject       		*WlzRegICPObjWSD2D(WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  double xMin,
				  double xMax,
				  double xStep,
				  double yMin,
				  double yMax,
				  double yStep,
				  double rMin,
				  double rMax,
				  double rStep,
				  double minDistWgt,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzRegICPVerticesWSD2D(WlzVertexP tVx,
				  WlzVertexP tNr,
				  int tCnt,
				  WlzVertexP sVx,
				  WlzVertexP sNr,
				  int sCnt,
				  WlzVertexType vType,
				  int sgnNrm,
				  WlzAffineTransform *initTr,
				  double xMin,
				  double xMax,
				  double xStep,
				  double yMin,
				  double yMax,
				  double yStep,
				  double rMin,
				  double rMax,
				  double rStep,
				  double minDistWgt,
				  WlzErrorNum *dstErr);
#endif /* !WLZ_EXT_BIND */
extern WlzAffineTransform	*WlzRegICPObjs(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv,
				  int *dstItr,
				  int maxItr,
				  double delta,
				  double minDistWgt,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzRegICPObjsGrd(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  double ctrLo,
				  double ctrHi,
				  double ctrWth,
				  int *dstConv,
				  int *dstItr,
				  int maxItr,
				  double delta,
				  double minDistWgt,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzAffineTransform	*WlzRegICPVertices(
				  WlzVertexP tVx,
				  WlzVertexP tNr,
				  int tCnt,
				  WlzVertexP sVx,
				  WlzVertexP sNr,
				  int sCnt,
				  WlzVertexType vType,
				  int sgnNrm,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv,
				  int *dstItr,
				  int maxItr,
				  double delta,
				  double minDistWgt,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzRegICPTreeAndVertices(
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
				  int *dstConv,
				  WlzRegICPUsrWgtFn usrWgtFn,
				  void *usrWgtData,
				  double delta,
				  double minDistWgt,
				  WlzErrorNum *dstErr); 
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzRsvFilter.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum		WlzRsvFilterBuffer(
				  WlzRsvFilter *ftr,
				  double *datBuf,
				  double *tmpBuf0,
				  double *tmpBuf1,
				  int bufSz);
#endif /* WLZ_EXT_BIND */
extern void			WlzRsvFilterFreeFilter(
				  WlzRsvFilter *ftr);
extern WlzRsvFilter		*WlzRsvFilterMakeFilter(
				  WlzRsvFilterName name,
				  double prm,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzRsvFilterObj(
				  WlzObject *srcObj,
				  WlzRsvFilter *ftr,
				  int actionMsk,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzSampleObj.c							*
************************************************************************/
extern WlzObject 		*WlzSampleObj(
				  WlzObject *srcObj,
				  WlzIVertex3 samFac,
				  WlzSampleFn samFn,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzSampleObjPoint2D(
				  WlzObject *srcObj,
				  WlzIVertex2 samFac,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzSampleObjPoint3D(
				  WlzObject *srcObj,
				  WlzIVertex3 samFac,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzSampleValuesAndCoords.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum     		WlzSampleValuesAndCoords(
				  WlzObject *obj,
				  WlzGreyType *dstGType,
				  int *dstNVal,
				  WlzGreyP *dstValP,
				  WlzVertexP *dstCoords,
				  WlzSampleFn samFn,
				  int samFac);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzScalarArithmeticOp.c
************************************************************************/
extern WlzErrorNum     		WlzGreyIncValues2D(
				  WlzObject *obj);
extern WlzErrorNum		WlzGreyIncValuesInDomain(
				  WlzObject *gObj,
				  WlzObject *dObj);
extern WlzObject		*WlzScalarAdd(
				  WlzObject *o1,
				  WlzPixelV pval,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzScalarSubtract(
				  WlzObject *o1,
				  WlzPixelV pval,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzScalarMultiply(
				  WlzObject *o1,
				  WlzPixelV pval,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzScalarDivide(
				  WlzObject *o1,
				  WlzPixelV pval,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzScalarBinaryOp2(
				  WlzObject *o1,
				  WlzPixelV pval,
				  WlzBinaryOperatorType op,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzScalarMulAdd(
				  WlzObject *obj,
				  WlzPixelV m,
				  WlzPixelV a,
				  WlzGreyType gType,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzScalarBinaryOp.c
************************************************************************/
extern WlzErrorNum 		WlzScalarBinaryOp(
				  WlzObject *o1,
				  WlzPixelV pval,
				  WlzObject *o3,
				  WlzBinaryOperatorType op);

/************************************************************************
* WlzScalarFeatures2D.c
************************************************************************/
extern WlzErrorNum     		WlzScalarFeatures2D(
				  WlzObject *obj,
				  int *dstSizeArrayFeat,
				  WlzIVertex2 **dstArrayFeat,
				  WlzScalarFeatureType fType,
				  WlzThresholdType thrHL,
				  WlzPixelV thrV,
				  double filterV,
				  double minDist);

/************************************************************************
WlzScalarFn.c
************************************************************************/
extern WlzObject 		*WlzScalarFn(
				  WlzObject *inObj,
				  WlzFnType fn,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzSepTrans.c
************************************************************************/
#ifndef WLZ_EXT_BIND
typedef WlzErrorNum 		(*WlzIntervalConvFunc)(
				  WlzSepTransWSpace *stwspc,
				  void *params);
extern WlzObject 		*WlzSepTrans(
				  WlzObject *obj,
				  WlzIntervalConvFunc x_fun,
				  void *x_params,
				  WlzIntervalConvFunc y_fun,
				  void *y_params,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzSeqPar.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject 		*WlzSeqPar(
				  WlzObject *inObj,
				  int newObjFlg,
				  int sequentialFlg,
				  WlzRasterDir rasterDir,
				  int bdrSz,
				  int bkgVal,
				  void *transformData,
				  int (*transformFn)(WlzSeqParWSpace *, void *),
				  WlzErrorNum	*dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzShadeCorrect.c								*
************************************************************************/
extern WlzObject		*WlzShadeCorrect(
				  WlzObject *srcObj,
				  WlzObject *shdObj,
				  double nrmVal,
				  int inPlace,
				  WlzErrorNum *dstErr);

extern WlzObject		*WlzShadeCorrectBFDF(
				  WlzObject *srcObj,
				  WlzObject *shdBFObj,
				  WlzObject *shdDFObj,
				  double nrmVal,
				  int inPlace,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzShift.c								*
************************************************************************/
extern WlzObject		*WlzShiftObject(
				  WlzObject *inObj,
				  int xShift,
				  int yShift,
				  int zShift,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern WlzDomain	 	WlzShiftDomain(
				  WlzObjectType inObjType,
				  WlzDomain inDom,
				  int xShift,
				  int yShift,
				  int zShift,
				  WlzErrorNum *dstErr);
extern WlzValues	 	WlzShiftValues(
				  WlzObjectType inObjType,
				  WlzValues inVal,
				  WlzDomain inDom,
				  int xShift,
				  int yShift,
				  int zShift,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzSkeleton.c								*
************************************************************************/
extern WlzObject 		*WlzSkeleton(
				  WlzObject *srcObj,
				  int smoothpasses,
				  WlzConnectType minCon,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzSobel.c								*
************************************************************************/
extern WlzObject 		*WlzSobel(
				  WlzObject *srcObj,
				  int hFlg,
				  int vFlg,
				  WlzErrorNum	*dstErr);

/************************************************************************
* WlzSnapFit.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum     		WlzSnapFit(
				  WlzObject *tObj,
				   WlzObject *sObj,
				   WlzAffineTransform *tr,
				   WlzVertexType *vType,
				   int *dstNVtx,
				   WlzVertexP *dstTVtxP,
				   WlzVertexP *dstSVtxP,
				   double maxCDist,
				   double minTDist,
				   double minSDist);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzSplitObj.c        							*
************************************************************************/
extern WlzErrorNum		WlzSplitObj(
				  WlzObject *refObj,
				  WlzObject *ppObj,
				  int bWidth,
				  double bgdFrac,
				  double sigma,
				  WlzCompThreshType compThrMethod,
				  int nReqComp,
				  int *dstSizeArrayComp,
				  WlzObject ***dstArrayComp);
extern WlzErrorNum		WlzSplitMontageObj(
				  WlzObject *mObj,
				  WlzPixelV gapV,
				  double tol,
				  int bWidth,
				  int minArea,
				  int maxComp,
				  int *dstNComp,
				  WlzObject ***dstComp);

/************************************************************************
* WlzStdStructElements.c						*
************************************************************************/
extern WlzObject 		*WlzMakeSpecialStructElement(
				  WlzSpecialStructElmType eType,
				  int elmIndex,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzMakeSinglePixelObject(
				  WlzObjectType spOType,
				  int k,
				  int l,
				  int p,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzMakeCircleObject(
				  double radius,
				  double x,
				  double y,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzMakeSphereObject(
				  WlzObjectType sOType,
				  double radius,
				  double x,
				  double y,
				  double z,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzMakeStdStructElement(
				  WlzObjectType ssTType,
				  WlzDistanceType dType,
				  double radius,
				  WlzErrorNum *dstErr);
extern WlzObject                *WlzMakeRectangleObject(
				  double radiusX,
				  double radiusY,
				  double x,
				  double y,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzMakeQuadrilateral(
				  double x0,
				  double y0,
				  double x1,
				  double y1,
				  double x2,
				  double y2,
				  double x3,
				  double y3,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzMakeCuboidObject(
				  WlzObjectType sOType,
				  double radiusX,
				  double radiusY,
				  double radiusZ,
				  double x,
				  double y,
				  double z,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzStringTypes.c							*
************************************************************************/
extern const char 		*WlzStringFromObjType(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern const char		*WlzStringFromObjTypeValue(
				  WlzObjectType objType,
				  WlzErrorNum *dstErr);
extern WlzObjectType 		WlzStringToObjType(
				  const char *oTypeStr,
				  WlzErrorNum *dstErr);
extern const char 		*WlzStringFromObjDomainType(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzObjectType 		WlzStringToObjDomainType(
				  const char *oDomTypeStr,
				  WlzErrorNum *dstErr);
extern const char 		*WlzStringFromObjValuesType(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzObjectType 		WlzStringToObjValuesType(
				  const char *oValTypeStr,
				  WlzErrorNum *dstErr);
extern const char 		*WlzStringFromGreyType(
				  WlzGreyType gType,
				  WlzErrorNum *dstErr);
extern const char      		*WlzStringFromScalarFeatureType(
				  WlzScalarFeatureType fType,
				  WlzErrorNum *dstErr);
extern WlzScalarFeatureType 	WlzStringToScalarFeatureType(
				  const char *fTypeStr,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
extern const char		*WlzStringFromPropertyType(
				  WlzProperty prop,
				  WlzErrorNum *dstErr);
#endif /* !WLZ_EXT_BIND */
extern WlzObjectType		WlzStringToPropertyType(
				  const char *pStr,
				  WlzErrorNum *dstErr);
extern const char      		*WlzStringFromEMAPPropertyType(
				  WlzEMAPProperty *eProp,
				  WlzErrorNum *dstErr);
extern WlzEMAPPropertyType 	WlzStringToEMAPPropertyType(
				  const char *pStr,
				  WlzErrorNum *dstErr);
extern WlzTransformType 	WlzStringToTransformType(
				  const char *tStr,
				  WlzErrorNum *dstErr);
extern const char 		*WlzStringFromTransformType(
				  WlzTransformType tType,
				  WlzErrorNum *dstErr);
extern WlzMeshGenMethod		WlzStringToMeshGenMethod(
				  const char *tStr,
				  WlzErrorNum *dstErr);
extern const char		*WlzStringFromMeshGenMethod(
				  WlzMeshGenMethod mtd,
				  WlzErrorNum *dstErr);
extern WlzFnType		WlzStringToFnType(
				  const char *tStr,
				  WlzErrorNum *dstErr);
extern const char		*WlzStringFromFnType(
				  WlzFnType fn,
				  WlzErrorNum *dstErr);
extern WlzGMModelType		WlzStringToGMModelType(
				  const char *tStr,
				  WlzErrorNum *dstErr);
extern const char		*WlzStringFromGMModelType(
				  WlzGMModelType mType,
				  WlzErrorNum *dstErr);
extern WlzGreyType 		WlzStringToGreyType(
				  const char *gStr,
				 WlzErrorNum *dstErr);
extern const char      		*WlzStringFromInterpolationType(
				  WlzInterpolationType iType,
                                  WlzErrorNum *dstErr);
extern WlzInterpolationType 	WlzStringToInterpolationType(
				  const char *iStr,
                                  WlzErrorNum *dstErr);
extern const char		*WlzStringFromValueAttachType(
				  WlzValueAttach at,
				  WlzErrorNum *dstErr);
extern WlzValueAttach 		WlzStringToValuesAttachType(
				  const char *aStr,
                                  WlzErrorNum *dstErr);
extern const char		*WlzStringFromMarkerType(
				  WlzMarkerType mtype,
				  WlzErrorNum *dstErr);
extern WlzMarkerType		WlzStringToMarkerType(
				  const char *mStr,
				  WlzErrorNum *dstErr);
extern const char      		*WlzStringFromThreeDViewMode(
				  WlzThreeDViewMode vMode,
                                  WlzErrorNum *dstErr);
extern WlzThreeDViewMode 	WlzStringToThreeDViewMode(
				  const char *vStr,
                                  WlzErrorNum *dstErr);
extern const char 		*WlzStringFromErrorNum(
				  WlzErrorNum gvnErr,
			          const char **dstMsgStr);
extern WlzErrorNum  		WlzStringToErrorNum(
				  const char *gStr);

/************************************************************************
* WlzStringUtils.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern int 			WlzStringMatchValue(
				  int *dstValue,
				  const char *targetStr,
				  const char *testStr,
				  ...);
extern int 			WlzValueMatchString(
				  char **dstStr,
				  int	targetVal,
				  const char *testStr,
				  ...);
extern char 			*WlzStringWhiteSpSkip(
				  char *str);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzStructDilation.c							*
************************************************************************/
extern WlzObject 		*WlzStructDilation(
				  WlzObject *obj,
				  WlzObject *structElm,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzStructErosion.c							*
************************************************************************/
extern WlzObject 		*WlzStructErosion(WlzObject *obj,
				  WlzObject *structElm,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRGBAConvert.c							*
************************************************************************/
extern WlzCompoundArray		*WlzRGBAToCompound(
				  WlzObject *obj,
				  WlzRGBAColorSpace colSpc,
				  WlzErrorNum *dstErr);
  
  
extern WlzObject 		*WlzCompoundToRGBA(
				  WlzCompoundArray *cmpnd,
				  WlzRGBAColorSpace colSpc,
				  WlzErrorNum *dstErr);
                               
extern WlzObject		*WlzRGBAToModulus(
                                  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzIndexToRGBA(
				  WlzObject *obj,
				  unsigned char colormap[3][256],
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzRGBAToChannel(
				  WlzObject *obj,
				  WlzRGBAColorChannel chan,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRGBAImageArithmetic.c			       			*
************************************************************************/
extern WlzObject 		*WlzRGBAImageArithmetic(
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzBinaryOperatorType	op,
				  int overwrite,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRGBAScalarBinaryOp.c			       			*
************************************************************************/
extern WlzObject 		*WlzRGBAScalarBinaryOp(
				  WlzObject		*o1,
				  WlzPixelV		pval,
				  WlzBinaryOperatorType op,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRGBARange.c							*
************************************************************************/
extern WlzErrorNum 		WlzRGBAModulusRange(
				  WlzObject *obj,
				  double *min,
				  double *max);

/************************************************************************
* WlzRGBAPixelUtils.c							*
************************************************************************/
extern double			WlzRGBAPixelValue(
				  WlzPixelV pixVal,
				  WlzRGBAColorChannel chan,
				  WlzErrorNum *dstErr);
extern void 			WlzRGBAConvertRGBToHSV_UBYTENormalised(
				  int *col);

/************************************************************************
* WlzRGBAMultiThreshold.c						*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject 		*WlzRGBAMultiThreshold(
				  WlzObject *obj,
				  WlzPixelV lowVal,
				  WlzPixelV highVal,
				  WlzUInt combineMode,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzObject 		*WlzRGBASliceThreshold(
				  WlzObject *obj,
				  WlzPixelV lowVal,
				  WlzPixelV highVal,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzRGBABoxThreshold(
				  WlzObject *obj,
				  WlzPixelV lowVal,
				  WlzPixelV highVal,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzRGBAEllipsoidThreshold(
				  WlzObject *obj,
				  WlzPixelV lowVal,
				  WlzPixelV highVal,
				  double eccentricity,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRGBAGreyStats.c							*
************************************************************************/
extern int 			WlzRGBAGreyStats(
				  WlzObject *srcObj,
				  WlzRGBAColorSpace colSpc,
				  WlzGreyType *dstGType,
				  double *dstMin,
				  double *dstMax,
				  double *dstSum,
				  double *dstSumSq,
				  double *dstMean,
				  double *dstStdDev,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRGBAModGradient.c							*
************************************************************************/
extern WlzObject 		*WlzRGBAModGradient(
				  WlzObject *obj,
				  double width,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzRGBChanRatio.c							*
************************************************************************/
extern WlzObject		*WlzRGBChanRatio(WlzObject *rgbObj,
				  WlzRGBAColorChannel num,
				  WlzRGBAColorChannel den,
				  WlzRGBAColorChannel mul,
				  int useMul,
				  int norm,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzCbThreshold.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObject		*WlzCbThreshold(
				  WlzObject *obj,
				  WlzThreshCbFn threshCb,
				  void *clientData,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzThreshold.c							*
************************************************************************/
extern WlzObject 		*WlzThreshold(
				  WlzObject *obj,
				  WlzPixelV threshV,
				  WlzThresholdType highlow,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzTiledValues.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzTiledValues		*WlzMakeTiledValues(
				  int dim,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzFreeTiledValues(
				  WlzTiledValues *tVal);
extern WlzErrorNum		WlzMakeTiledValuesTiles(
				  WlzTiledValues *tVal);
extern WlzObject		*WlzMakeTiledValuesFromObj(
				  WlzObject *gObj,
				  size_t tileSz,
				  int copyValues,
				  WlzGreyType gType,
				  WlzPixelV bgdV,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzTransform.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzEmptyTransform	*WlzMakeEmptyTransform(
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzFreeTransform(
				  WlzTransform tr);
extern WlzErrorNum		WlzFreeEmptyTransform(
				  WlzEmptyTransform *tr);
extern WlzTransform		WlzTransformProduct(
				  WlzTransform tr0,
				  WlzTransform tr1,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzTransposeObj.c							*
************************************************************************/
extern WlzObject 		*WlzTransposeObj(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzUnion2.c								*
************************************************************************/
extern WlzObject		*WlzUnion2(
				  WlzObject *obj1,
				  WlzObject *obj2,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzUnionN.c
************************************************************************/
extern WlzObject 		*WlzUnionN(
				  int sizeArrayObjs,
				  WlzObject **arrayObjs,
				  int uvt,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzValueTableUtils.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzObjectType 		WlzGreyTableType(
				  WlzObjectType table_type,
				  WlzGreyType grey_type,
				  WlzErrorNum *dstErr);
extern WlzGreyType 		WlzGreyTableTypeToGreyType(
				  WlzObjectType objType,
				  WlzErrorNum *dstErr);
extern WlzObjectType 		WlzGreyTableTypeToTableType(
				  WlzObjectType objType,
				  WlzErrorNum *dstErr);
extern WlzObjectType		WlzGreyTableIsTiled(
				  WlzObjectType gTabType);
extern WlzGreyType		WlzGreyTypeFromObj(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzDVertex3		WlzVozelSz(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern void			*WlzIndexedValueGet(
				  WlzIndexedValues *ixv,
				  int idx);
extern void			*WlzIndexedValueExtGet(
				  WlzIndexedValues *ixv,
				  int idx);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzValueUtils.c							*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern void			WlzValueSetInt(
				  int *vec,
				  int value,
				  size_t count);
extern void			WlzValueSetShort(
				  short *vec,
				  short value,
				  size_t count);
extern void			WlzValueSetUByte(
				  WlzUByte *vec,
				  WlzUByte value,
				  size_t count);
extern void			WlzValueSetFloat(
				  float *vec,
				  float value,
				  size_t count);
extern void			WlzValueSetDouble(
				  double *vec,
				  double value,
				  size_t count);
extern void			WlzValueSetRGBA(
				  WlzUInt *vec,
				  WlzUInt value,
				  size_t count);
extern void			WlzValueSetDVertex(
				  WlzDVertex2 *vec,
				  WlzDVertex2 value,
				  size_t count);
extern void			WlzValueSetFVertex(
				  WlzFVertex2 *vec,
				  WlzFVertex2 value,
				  size_t count);
extern void			WlzValueSetIVertex(
				  WlzIVertex2 *vec,
				  WlzIVertex2 value,
				  size_t count);
extern void			WlzValueSetGrey(
				  WlzGreyP vec,
				  size_t vecOff,
				  WlzGreyV value,
				  WlzGreyType gType,
				  size_t count);
extern void			WlzValueClampGreyIntoGrey(
				  WlzGreyP dst,
				  size_t dstOff,
				  WlzGreyType dstType,
				  WlzGreyP src,
				  size_t srcOff,
				  WlzGreyType srcType,
				  size_t count);
extern void			WlzValueClampIntIntoShort(
				  short *dst,
				  int *src,
				  size_t count);
extern void			WlzValueClampIntIntoUByte(
				  WlzUByte *dst,
				  int *src,
				  size_t count);
extern void			WlzValueClampShortIntoUByte(
				  WlzUByte *dst,
				  short *src,
				  size_t count);
extern void			WlzValueClampFloatIntoInt(
				  int *dst,
				  float *src,
				  size_t count);
extern void			WlzValueClampFloatIntoShort(
				  short *dst,
				  float *src,
				  size_t count);
extern void			WlzValueClampFloatIntoUByte(
				  WlzUByte *dst,
				  float *src,
				  size_t count);
extern void			WlzValueClampDoubleIntoInt(
				  int *dst,
				  double *src,
				  size_t count);
extern void			WlzValueClampDoubleIntoShort(
 				  short *dst,
				  double *src,
				  size_t count);
extern void			WlzValueClampDoubleIntoUByte(
				  WlzUByte *dst,
				  double *src,
				  size_t count);
extern void			WlzValueClampDoubleIntoFloat(
				  float *dst,
				  double *src,
				  size_t count);
extern void			WlzValueClampDoubleIntoRGBA(
				  WlzUInt *dst,
				  double *src,
				  size_t count);
extern void			WlzValueClampIntIntoRGBA(
				  WlzUInt *dst,
				  int *src,
				  size_t count);
extern void			WlzValueClampShortIntoRGBA(
				  WlzUInt *dst,
				  short *src,
				  size_t count);
extern void			WlzValueClampFloatIntoRGBA(
				  WlzUInt *dst,
				  float *src,
				  size_t count);
extern void			WlzValueClampDoubleIntoRGBA(
				  WlzUInt *dst,
				  double *src,
				  size_t count);
extern void			WlzValueClampIntToShort(
				  int *vec,
				  size_t count);
extern void			WlzValueClampIntToUByte(
				  int *vec,
				  size_t count);
extern void			WlzValueClampShortToUByte(
				  short *vec,
				  size_t count);
extern void			WlzValueClampDoubleToInt(
				  double *vec,
				  size_t count);
extern void			WlzValueClampDoubleToShort(
				  double *vec,
				  size_t count);
extern void			WlzValueClampDoubleToUByte(
				  double *vec,
				  size_t count);
extern void			WlzValueClampDoubleToRGBA(
				  double *vec,
				  size_t count);
extern void			WlzValueClampDoubleToFloat(
				  double *vec,
				  size_t count);
extern void			WlzValueClampFloatToInt(
				  float *vec,
				  size_t count);
extern void			WlzValueClampFloatToShort(
				  float *vec,
				  size_t count);
extern void			WlzValueClampFloatToUByte(
				  float *vec,
				  size_t count);
extern void			WlzValueCopyIntToInt(
				  int *dst,
				  int *src,
				  size_t count);
extern void			WlzValueCopyIntToShort(
				  short *dst,
				  int *src,
				  size_t count);
extern void			WlzValueCopyIntToUByte(
				  WlzUByte *dst,
				  int *src,
				  size_t count);
extern void			WlzValueCopyIntToFloat(
				  float *dst,
				  int *src,
				  size_t count);
extern void			WlzValueCopyIntToDouble(
				  double *dst,
				  int *src,
				  size_t count);
extern void			WlzValueCopyIntToRGBA(
				  WlzUInt *dst,
				  int *src,
				  size_t count);
extern void			WlzValueCopyShortToInt(
				  int *dst,
				  short *src,
				  size_t count);
extern void			WlzValueCopyShortToShort(
				  short *dst,
				  short *src,
				  size_t count);
extern void			WlzValueCopyShortToUByte(
				  WlzUByte *dst,
				  short *src,
				  size_t count);
extern void			WlzValueCopyShortToFloat(
				  float *dst,
				  short *src,
				  size_t count);
extern void			WlzValueCopyShortToDouble(
				  double *dst,
				  short *src,
				  size_t count);
extern void			WlzValueCopyShortToRGBA(
				  WlzUInt *dst,
				  short *src,
				  size_t count);
extern void			WlzValueCopyUByteToInt(
				  int *dst,
				  WlzUByte *src,
				  size_t count);
extern void			WlzValueCopyUByteToShort(
				  short *dst,
				  WlzUByte *src,
				  size_t count);
extern void			WlzValueCopyUByteToUByte(
				  WlzUByte *dst,
				  WlzUByte *src,
				  size_t count);
extern void			WlzValueCopyUByteToFloat(
				  float *dst,
				  WlzUByte *src,
				  size_t count);
extern void			WlzValueCopyUByteToDouble(
				  double *dst,
				  WlzUByte *src,
				  size_t count);
extern void			WlzValueCopyUByteToRGBA(
				  WlzUInt *dst,
				  WlzUByte *src,
				  size_t count);
extern void			WlzValueCopyFloatToInt(
				  int *dst,
				  float *src,
				  size_t count);
extern void			WlzValueCopyFloatToShort(
				  short *dst,
				  float *src,
				  size_t count);
extern void			WlzValueCopyFloatToUByte(
				  WlzUByte *dst,
				  float *src,
				  size_t count);
extern void			WlzValueCopyFloatToFloat(
				  float *dst,
				  float *src,
				  size_t count);
extern void			WlzValueCopyFloatToDouble(
				  double *dst,
				  float *src,
				  size_t count);
extern void			WlzValueCopyFloatToRGBA(
				  WlzUInt *dst,
				  float *src,
				  size_t count);
extern void			WlzValueCopyDoubleToInt(
				  int *dst,
				  double *src,
				  size_t count);
extern void			WlzValueCopyDoubleToShort(
				  short *dst,
				  double *src,
				  size_t count);
extern void			WlzValueCopyDoubleToUByte(
				  WlzUByte *dst,
				  double *src,
				  size_t count);
extern void			WlzValueCopyDoubleToFloat(
				  float *dst,
				  double *src,
				  size_t count);
extern void			WlzValueCopyDoubleToDouble(
				  double *dst,
				  double *src,
				  size_t count);
extern void			WlzValueCopyDoubleToRGBA(
				  WlzUInt *dst,
				  double *src,
				  size_t count);
extern void			WlzValueCopyRGBAToInt(
				  int *dst,
				  WlzUInt *src,
				  size_t count);
extern void			WlzValueCopyRGBAToShort(
				  short *dst,
				  WlzUInt *src,
				  size_t count);
extern void			WlzValueCopyRGBAToUByte(
				  WlzUByte *dst,
				  WlzUInt *src,
				  size_t count);
extern void			WlzValueCopyRGBAToFloat(
				  float *dst,
				  WlzUInt *src,
				  size_t count);
extern void			WlzValueCopyRGBAToDouble(
				  double *dst,
				  WlzUInt *src,
				  size_t count);
extern void			WlzValueCopyRGBAToRGBA(
				  WlzUInt *dst,
				  WlzUInt *src,
				  size_t count);
extern void			WlzValueCopyGreyToGrey(
				  WlzGreyP dst,
				  size_t dstOff,
				  WlzGreyType dstType,
				  WlzGreyP src,
				  size_t srcOff,
				  WlzGreyType srcType,
				  size_t count);
extern void			WlzValueCopyDVertexToDVertex(
				  WlzDVertex2 *dst,
				  WlzDVertex2 *src,
				  size_t count);
extern void			WlzValueCopyDVertexToFVertex(
				  WlzFVertex2 *dst,
				  WlzDVertex2 *src,
				  size_t count);
extern void			WlzValueCopyDVertexToIVertex(
				  WlzIVertex2 *dst,
				  WlzDVertex2 *src,
				  size_t count);
extern void			WlzValueCopyFVertexToDVertex(
				  WlzDVertex2 *dst,
				  WlzFVertex2 *src,
				  size_t count);
extern void			WlzValueCopyFVertexToFVertex(
				  WlzFVertex2 *dst,
				  WlzFVertex2 *src,
				  size_t count);
extern void			WlzValueCopyFVertexToIVertex(
				  WlzIVertex2 *dst,
				  WlzFVertex2 *src,
				  size_t count);
extern void			WlzValueCopyIVertexToDVertex(
				  WlzDVertex2 *dst,
				  WlzIVertex2 *src,
				  size_t count);
extern void			WlzValueCopyIVertexToFVertex(
				  WlzFVertex2 *dst,
				  WlzIVertex2 *src,
				  size_t count);
extern void			WlzValueCopyIVertexToIVertex(
				  WlzIVertex2 *dst,
				  WlzIVertex2 *src,
				  size_t count);
extern void			WlzValueCopyDVertexToDVertex3(
				  WlzDVertex3 *dst,
				  WlzDVertex3 *src,
				  size_t count);
extern void			WlzValueCopyDVertexToFVertex3(
				  WlzFVertex3 *dst,
				  WlzDVertex3 *src,
				  size_t count);
extern void			WlzValueCopyDVertexToIVertex3(
				  WlzIVertex3 *dst,
				  WlzDVertex3 *src,
				  size_t count);
extern void			WlzValueCopyFVertexToDVertex3(
				  WlzDVertex3 *dst,
				  WlzFVertex3 *src,
				  size_t count);
extern void			WlzValueCopyFVertexToFVertex3(
				  WlzFVertex3 *dst,
				  WlzFVertex3 *src,
				  size_t count);
extern void			WlzValueCopyFVertexToIVertex3(
				  WlzIVertex3 *dst,
				  WlzFVertex3 *src,
				  size_t count);
extern void			WlzValueCopyIVertexToDVertex3(
				  WlzDVertex3 *dst,
				  WlzIVertex3 *src,
				  size_t count);
extern void			WlzValueCopyIVertexToFVertex3(
				  WlzFVertex3 *dst,
				  WlzIVertex3 *src,
				  size_t count);
extern void			WlzValueCopyIVertexToIVertex3(
				  WlzIVertex3 *dst,
				  WlzIVertex3 *src,
				  size_t count);
extern void			WlzIndexedValueBufWeight(
				  WlzGreyP gP,
				  int gO,
				  WlzIndexedValues *ixv,
				  double *wgt,
				  int nIdx,
				  int *idx);
extern WlzErrorNum 		WlzValueConvertPixel(
				  WlzPixelV *dstPix,
				  WlzPixelV srcPix,
				  WlzGreyType dstType);
extern int			WlzValueMedianInt(
				  int *vec,
				  size_t count);
extern double			WlzValueMedianDouble(
				  double *vec,
				  size_t count);
extern int			WlzValueDitherI(
				  int p0,
				  int p1,
				  int p2);
extern size_t			WlzIndexedValueSize(
				  WlzIndexedValues *ixv,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzVerifyObj.c							*
************************************************************************/
extern WlzErrorNum 		WlzVerifyObject(
				  WlzObject *obj,
				  int fix);
#ifndef WLZ_EXT_BIND
extern WlzErrorNum 		WlzVerifyIntervalDomain(
				  WlzDomain dom,
				  int fix);
extern WlzErrorNum 		WlzVerifyIntervalLine(
				  WlzIntervalLine *intvline,
				  int fix);
extern WlzErrorNum 		WlzVerifyInterval(
				  WlzInterval *intv,
				  int fix);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzVersion.c 	 							*
************************************************************************/
extern const char		*WlzVersion();

/************************************************************************
* WlzVerticies.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzVertexP		WlzVerticesFromObj(
				  WlzObject *obj,
				  WlzVertexP *dstNr,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
extern WlzVertexP		WlzVerticesFromObjBnd(
				  WlzObject *obj,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzErrorNum		WlzVerticesFromObjBnd2I(
				  WlzObject *obj,
				  int *dstSizeArrayVtx,
				  WlzIVertex2 **dstArrayVtx);
extern WlzErrorNum		WlzVerticesFromObjBnd3I(
				  WlzObject *obj,
				  int *dstSizeArrayVtx,
				  WlzIVertex3 **dstArrayVtx);
extern WlzErrorNum		WlzVerticesFromObj2I(
				  WlzObject *obj,
				  int *dstSizeArrayVtx,
				  WlzIVertex2 **dstArrayVtx);
extern WlzErrorNum		WlzVerticesFromObj3I(
				  WlzObject *obj,
				  int *dstSizeArrayVtx,
				  WlzIVertex3 **dstArrayVtx);
#ifndef WLZ_EXT_BIND
extern WlzVertexP      		WlzDVerticesFromCMesh(
				  WlzCMeshP mesh,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  int skip,
				  WlzErrorNum *dstErr);
extern WlzVertexP		WlzDVerticesFromGM(
				  WlzGMModel *model,
				  int *dstCnt,
				  WlzVertexType *dstType,
				  WlzErrorNum *dstErr);
extern WlzVertexP               WlzVerticesFromGM(
                                  WlzGMModel *model,
                                  WlzVertexP *dstNr,
                                  int **dstVId,
                                  int *dstCnt,
                                  WlzVertexType *dstType,
                                  WlzErrorNum *dstErr);
extern AlcKDTTree      		*WlzVerticesBuildTree(
				  WlzVertexType vType,
				  int nV,
				  WlzVertexP vtx,
				  int *shfBuf,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */


/************************************************************************
* WlzVolume.c								*
************************************************************************/
extern int 			WlzVolume(
				  WlzObject *obj,
		     		  WlzErrorNum *dstErr);

/************************************************************************
* WlzEmpty.c								*
************************************************************************/
extern int 			WlzIsEmpty(
				  WlzObject *obj,
		     		  WlzErrorNum *dstErr);

/************************************************************************
* WlzWindow.c								*
************************************************************************/
extern WlzObject 		*WlzWindow(
				  WlzObject *srcObj,
				  WlzWindowFnType winFn,
				  WlzIVertex2 origin,
				  WlzIVertex2 radius,
				  WlzErrorNum	*dstErr);
extern WlzWindowFnType 		WlzWindowFnValue(
				  const char *winFnName);
extern const char 		*WlzWindowFnName(
				  WlzWindowFnType winFnValue);

/************************************************************************
* WlzWriteObj.c								*
************************************************************************/
extern WlzErrorNum 		WlzWriteObj(
				  FILE *fp,
			          WlzObject *obj);

#ifndef WLZ_EXT_BIND
extern WlzErrorNum  		WlzWriteMeshTransform3D(
				  FILE *fp,
			          WlzMeshTransform3D *obj);
#endif /* !WLZ_EXT_BIND */
				  
/************************************************************************
* WlzMwrAngle.c								*
************************************************************************/
extern double 			WlzMwrAngle(
				  WlzObject *cvh,
			  	  WlzErrorNum *dstErr);

/************************************************************************
* Wlz3DWarpMQ_S.c					         	*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzErrorNum 		WlzTetrahedronProducerFromCube(
				  int neighbourLeft,
				  int neighbourRight,
				  int neighbourBack,
				  int neighbourFront,
				  int neighbourDown,
				  int neighbourUp,
				  int nxmax,
				  int nymax,
				  int nzmax,
				  int nx,
				  int ny,
				  int nz,
				  int indexOfNextFirst,
				  WlzMeshElem3D *elements);
extern WlzMeshTransform3D 	*WlzTetrahedronMeshFromObj(
				  WlzObject *wObjC,
				  const WlzDBox3 bBoxS,
				  const int numOfElemAlonX,
				  const int numOfElemAlonY,
				  const int numOfElemAlonZ,
				  WlzErrorNum *errNum);
extern WlzErrorNum  		read_WlzTiePoints(
				  FILE *fp,
				  int *nTiePP,
				  WlzDVertex3 **vxVec0,
				  WlzDVertex3 **vxVec1,
				  const int ReadDisplacement);
extern 	WlzErrorNum 		WlzTetrahedronProducer(
				  int nxmax,
				  int nymax,
				  int nzmax,
				  double xmin,
				  double xmax,
				  double ymin ,
				  double ymax,
				  double zmin,
				  double zmax,
				  char *strchar,
				  int nCutBelow,
				  double zConst,
				  int nCutUp,
				  int outputMesh3D,
				  int outputTranMesh3D,
				  int outputCutedPlaneVTK,
				  int outputCutedPlaneTetraVTK,
				  int outputOriginalPlaneVTK,
				  int ReadDisplacement,
				  WlzIBox3 bBox0);
extern WlzErrorNum 		WlzGetTransformedMesh(
				  WlzMeshTransform3D *wmt3D,
				  WlzBasisFnTransform* basisTr);
extern WlzObject       		*WlzMeshTransformObj_3D(
				  WlzObject *srcObj,
				  WlzMeshTransform3D *wmt3D,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
extern  void 			WlzEffWriteMeshTransform3DWithDisplacementVTK(
				  FILE *fp,
               			  WlzMeshTransform3D *wmt3D);
extern  void 		     WlzEffWriteMeshTransform3DWithoutDisplacementVTK(
				  FILE *fp,
               			  WlzMeshTransform3D *wmt3D);
extern WlzMeshTransform2D5  	*Wlz2D5TransformFromCut3Dmesh(
				  double zConst,
				  WlzMeshTransform3D *wmt3D,
				  WlzErrorNum *disErr);
extern void			WlzEffWriteOriginalPlaneVTKByDis(
				  FILE *fp,
				  WlzMeshTransform2D5 *wmt2D5);
extern void			WlzEffWriteOriginalPlaneVTKByPos(
				  FILE *fp,
				  WlzMeshTransform2D5 *wmt2D5);
extern int   			WlzIsoIntersectWithTetrahadronIndex(
				  double zConst,
				  const WlzMeshTransform3D *wmt3D,
				  int *intersectIndex,
				  WlzDVertex3 *planepoints,
				  int **linkList,
				  int *noRedundancyCutingNum,
				  int *numOfTrangularElem);
extern void 			WlzMakeAffine3D4pointsTrFn(
				  WlzDVertex3 sr1,
				  WlzDVertex3 sr2,
				  WlzDVertex3 sr3,
				  WlzDVertex3 sr4,
				  WlzDVertex3 targ1,
				  WlzDVertex3 targ2,
				  WlzDVertex3 targ3,
				  WlzDVertex3 targ4,
				  double **Affine3D4pointsTrFun);
extern void 			WlzMakeAffine3D4pointsTrFnByGauss(
				  WlzDVertex3 sr1,
				  WlzDVertex3 sr2,
				  WlzDVertex3 sr3,
				  WlzDVertex3 sr4,
				  WlzDVertex3 targ1,
				  WlzDVertex3 targ2,
				  WlzDVertex3 targ3,
				  WlzDVertex3 targ4,
				  double **Affine3D4pointsTrFun);
					      
/************************************************************************
* WlzDistanceMap_S.c					         	*
************************************************************************/
extern WlzObject		*WlzGreyValueMixing_s(
				  WlzObject *srcObj,
				  WlzObject *tarObj,
				  double xmiddle,
				  WlzErrorNum *dstErr);

#endif /* !WLZ_EXT_BIND */

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
}
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#endif	/* !WLZ_PROTO_H Don't put anything after this line */
