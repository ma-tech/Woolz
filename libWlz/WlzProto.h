#ifndef WLZ_PROTO_H
#define WLZ_PROTO_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzProto.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Defines the Woolz function prototypes.
*		To allow bindings to other languages (ie Java) to be
*		automaticly generated from this file, the following
*		conventions MUST be followed:
*		* The preprocessor identifier #WLZ_EXT_BIND controls
*		  whether code is bound.
*		* All prototypes to be bound must provide both
*		  parameter types and parameter identifiers.
*		* If a function returns a WlzErrorNum type it is
*		  assumed to be the error code.
*		* If a parameter is used to return the error code then
*		  it must be declared as 'WlzErrorNum *dstErr' and it
*		  may be a good idea to make it the last parameter.
*		* If either a function return type or parameter is
*		  'char' and it is not a pointer then it is assummed
*		  to be a byte.
*		* All Woolz types and functions start with 'Wlz'.
*		* All enum's are typedef'd.
*		* All structs are typedef'd.
*		* If a parameter identifier starts with 'dst' and is
*		  a pointer then it is assumed to be a destination
*		  pointer.
*		* All function return types or parameters of type
*		  'char *' are assumed to be strings unless the
*		  parameter identifier starts with 'dst'.
*		* All arrays must either start with 'array' or
*		  'dstArray', with 'dstArray' being used to identify
*		  arrays which are allocated within the	function being
*		  called. An array pointer must be followed by it's
*		  size, the type should be int, WlzIVertex2, ... and
*		  the identifier should be the same as array's but with
*		  'size' prepended as is the example 'sizeArrayName'.
*		  Arrays can not be pointers to void and wrappers may
*		  be needed to avoid this.
*		* UBYTE should be expanded to unsigned char - (why?) RAB
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 19-09-02 nick Add WlzJavaUtils.c and moved WlzObjGetType() to it.
* 12-09-02 nick Add WlzObjGetType().
* 21-08-02 nick Add Wlz2DContains().
* 26-04-01 bill Add WlzContourObjGrd().
* 23-04-01 bill Add WlzRegICPVertices().
* 23-03-01 bill Add WlzGreyTypeFromObj().
* 14-03-01 bill Add WlzRasterObj().
* 01-02-01 bill Add WlzRegCCorObjs().
* 26-01-01 bill Change parameters of WlzGeomVtxSortRadial().
* 23-01-01 bill Add WlzShadeCorrect().
* 16-01-01 bill	Add WlzGeomTriangleNormal() and WlzGMVertexNormal3D().
* 22-12-00 bill Add initial affine transform parameter to WlzRegICPObjs().
*		Add normals to WlzVerticesFromObj().
* 01-12-00 RAB	Add WlzHasIntersection().
* 01-12-00 RAB	Add get and set functions for 3D view structure params.
* 29-11-00 bill Rename WlzAffineTransformLSq to WlzAffineTransformLSq2D
*		add WlzAffineTransformLSq3D and a new WlzAffineTransformLSq.
* 28-11-00 bill	Add WlzGMShellComputeGBB() and WlzGMFilterRmSmShells().
* 16-11-00 bill Add WlzGMModelDeleteS().
* 06-11-00 bill Add WlzVerticesFromObj() and WlzRegICPObjs().
* 03-11-00 bill Add WlzValueCopy[IFD]VertexTo[IFD]Vertex3().
* 10-10-00 bill Add WlzAffineTransformContour(),
*               WlzAffineTransformGMModel(),
*               WlzAffineTransformBBoxI2(),
*               WlzAffineTransformBBoxD2(),
*               WlzAffineTransformBBoxI3() and
*               WlzAffineTransformBBoxD3().
* 28-09-00 bill	Add WlzAffineTransformDimension(),
* 		WlzAffineTransformTranslationSet(),
*		WlzAffineTransformFromTranslation(),
*		WlzAffineTransformRotationSet(),
*		WlzAffineTransformFromScale(),
*		WlzAffineTransformScaleSet() and
*		WlzAffineTransformFromRotation().
* 10-08-00 bill	Modify WlzMakeContour().
*		Add WlzStringFromGMModelType(), WlzStringToGMModelType()
*		and WlzAssignGMModel().
* 17-03-00 bill Strip out all the WlzEdge functions and add
* 		functions for Woolz geometric models.
* 04-02-00 bill Add WlzHistogramFitPeaks().
* 02-02-00 bill Added WlzRsvFilterBuffer(), WlzHistogramRsvGauss(),
*		WlzHistogramRsvFilter() and modified WlzHistogramSmooth().
*		Add WlzHistogramConvolve() and WlzHistogramCnvGauss().
*		Add WlzHistogramFindPeaks().
* 29-11-00 bill Add all the WlzContour and WlzEdge functions.
* 31-08-00 bill Modify array size parameters for JavaWoolz.
*		Allowed more functions to be bound using JavaWoolz.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************************************************************
* Wlz2DContains.c							*
************************************************************************/
extern WlzObject 		*Wlz2DContains(
				  WlzObject *obj,
				  double x,
				  double y,
				  WlzErrorNum *dstErr);
/************************************************************************
* Wlz3DSection.c							*
************************************************************************/
extern WlzObject 		*WlzGetSectionFromObject(
				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzGetMaskedSectionFromObject(
				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  WlzErrorNum *dstErr);

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
extern WlzThreeDViewStruct 	*WlzRead3DViewStruct(
				  FILE *fp,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzWrite3DViewStruct(
				  FILE *fp,
				  WlzThreeDViewStruct *viewstr);
extern WlzThreeDViewStruct 	*WlzMake3DViewStruct(
				  WlzObjectType type,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzFree3DViewStruct(
				  WlzThreeDViewStruct *viewStr);
extern WlzErrorNum 		WlzInit3DViewStruct(
				  WlzThreeDViewStruct *viewStr,
				  WlzObject *obj);
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
                                  WlzThreeDViewStruct	*viewStr,
				  WlzDVertex3		bbMin,
				  WlzDVertex3		bbMax,
				  WlzDVertex3		*rtnVtxs,
				  WlzErrorNum		*dstErr);
#endif /* WLZ_EXT_BIND */
extern int		Wlz3DViewGetBoundingBoxIntersectionA(
                          WlzThreeDViewStruct	*viewStr,
			  int *dstSizeArrayVtxs,
			  WlzDVertex3 **dstArrayVtxs,
			  WlzErrorNum *dstErr);
extern WlzErrorNum Wlz3DViewGetFixed(WlzThreeDViewStruct	*vs,
				     double			*dstX,
				     double			*dstY,
				     double			*dstZ);
extern WlzErrorNum Wlz3DViewSetFixed(WlzThreeDViewStruct	*vs,	
				     double			x,
				     double	       		y,
				     double			z);
extern WlzErrorNum Wlz3DViewGetTheta(WlzThreeDViewStruct	*vs,
				     double			*dstVal);
extern WlzErrorNum Wlz3DViewSetTheta(WlzThreeDViewStruct	*vs,
				     double			val);
extern WlzErrorNum Wlz3DViewGetPhi(WlzThreeDViewStruct	*vs,
				   double		*dstVal);
extern WlzErrorNum Wlz3DViewSetPhi(WlzThreeDViewStruct	*vs,
				   double		val);
extern WlzErrorNum Wlz3DViewGetZeta(WlzThreeDViewStruct	*vs,
				    double		*dstVal);
extern WlzErrorNum Wlz3DViewSetZeta(WlzThreeDViewStruct	*vs,
				    double		val);
extern WlzErrorNum Wlz3DViewGetDist(WlzThreeDViewStruct	*vs,
				    double		*dstVal);
extern WlzErrorNum Wlz3DViewSetDist(WlzThreeDViewStruct	*vs,
				    double		val);
extern WlzErrorNum Wlz3DViewGetScale(WlzThreeDViewStruct	*vs,
				     double			*dstVal);
extern WlzErrorNum Wlz3DViewSetScale(WlzThreeDViewStruct	*vs,
				     double			val);
extern WlzErrorNum Wlz3DViewGetViewMode(WlzThreeDViewStruct	*vs,
					WlzThreeDViewMode	*dstVal);
extern WlzErrorNum Wlz3DViewSetViewMode(WlzThreeDViewStruct	*vs,
					WlzThreeDViewMode	val);
extern WlzErrorNum Wlz3DViewGetUp(WlzThreeDViewStruct	*vs,
				  double		*dstX,
				  double		*dstY,
				  double		*dstZ);
extern WlzErrorNum Wlz3DViewSetUp(WlzThreeDViewStruct	*vs,
				  double		x,
				  double	        y,
				  double		z);
extern WlzErrorNum Wlz3DViewGetFixed2(WlzThreeDViewStruct	*vs,
				      double			*dstX,
				      double			*dstY,
				      double			*dstZ);
extern WlzErrorNum Wlz3DViewSetFixed2(WlzThreeDViewStruct	*vs,
				      double			x,
				      double	      		y,
				      double			z);
extern WlzErrorNum Wlz3DViewGetFixedLineAngle(WlzThreeDViewStruct	*vs,
					      double		*dstVal);
extern WlzErrorNum Wlz3DViewSetFixedLineAngle(WlzThreeDViewStruct	*vs,
					      double		val);
extern WlzErrorNum Wlz3DViewGetMaxvals(WlzThreeDViewStruct	*vs,
				       double			*dstX,
				       double			*dstY,
				       double			*dstZ);
extern WlzErrorNum Wlz3DViewGetMinvals(WlzThreeDViewStruct	*vs,
				       double			*dstX,
				       double			*dstY,
				       double			*dstZ);
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
				  unsigned char *bitData,
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
extern WlzObject 		*WlzAffineTransformObj(
				  WlzObject *srcObj,
				  WlzAffineTransform *trans,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
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
extern WlzAffineTransform	*WlzAffineTransformLSqWgt(
				  WlzVertexType vtxType,
				  int nVtx,
				  double *wgt,
				  WlzVertexP vtxVec0,
				  WlzVertexP vtxVec1,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSq2(
				  WlzVertexType vtxType,
				  int nV,
				  double *vW,
				  WlzVertexP v0,
				  WlzVertexP v1,
				  int nN,
				  double *nW,
				  WlzVertexP n0,
				  WlzVertexP n1,
				  WlzTransformType trType,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAffineTransformLSq(
				  WlzVertexType vtxType,
				  int arraySizeVec0,
				  WlzVertexP arrayVec0,
				  int arraySizeVec1,
				  WlzVertexP arrayVec1,
				  WlzTransformType tType,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */
extern WlzAffineTransform 	*WlzAffineTransformLSq2D(
				  int arraySizeVec0,
				  WlzDVertex2 *arrayVec0,
				  int arraySizeVec1,
				  WlzDVertex2 *arrayVec1,
				  WlzTransformType tType,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzAffineTransformLSq3D(
				  int arraySizeVec0,
				  WlzDVertex3 *arrayVec0,
				  int arraySizeVec1,
				  WlzDVertex3 *arrayVec1,
				  WlzTransformType tType,
				  WlzErrorNum *dstErr);

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
				  WlzErrorNum *dstErr );
extern AlcDLPList               *WlzAssignPropertyList(
                                  AlcDLPList	*plist,
                                  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzAssignAffineTransform(
				  WlzAffineTransform *,
				  WlzErrorNum *dstErr);
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
* WlzBasisFnTransform.c							*
************************************************************************/
extern WlzBasisFnTransform	*WlzBasisFnTrFromCPts(
				  WlzBasisFnType type,
				  int order,
				  int sizeArrayDPts,
				  WlzDVertex2 *arrayDPts,
				  int sizeArraySPts,
				  WlzDVertex2 *arraySPts,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzBasisFnSetMesh(
				  WlzMeshTransform *mesh,
				  WlzBasisFnTransform *basis);
extern WlzErrorNum 		WlzBasisFnFreeTransform(
				  WlzBasisFnTransform *basis);
extern WlzObject		*WlzBasisFnTransformObj(
				  WlzObject *srcObj,
				  WlzBasisFnTransform *basis,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
extern WlzIVertex2		WlzBasisFnTransformVertexI(
				  WlzBasisFnTransform *basis,
				  WlzIVertex2 srcVxF,
				  WlzErrorNum *dstErr);
extern WlzFVertex2		WlzBasisFnTransformVertexF(
				  WlzBasisFnTransform *basis,
			          WlzFVertex2 srcVxF,
				  WlzErrorNum *dstErr);
extern WlzDVertex2		WlzBasisFnTransformVertexD(
				  WlzBasisFnTransform *basis,
			          WlzDVertex2 srcVxF,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzBoundaryUtils.c							*
************************************************************************/
extern WlzErrorNum		WlzBoundObjToPolyDomArray(
				  WlzObject *bndObj,
				  int *dstSizeArrayPoly,
				  WlzPolygonDomain ***dstArrayPoly);
extern int			WlzBoundPolyCount(
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
extern WlzIBox2        		WlzBoundingBox2DTo2I(
				  WlzDBox2 bBox2D);
extern WlzDBox2        		WlzBoundingBox2ITo2D(
				  WlzIBox2 bBox2I);
extern WlzIBox3        		WlzBoundingBox3DTo3I(
				  WlzDBox3 bBox3D);
extern WlzDBox3        		WlzBoundingBox3ITo3D(
				  WlzIBox3 bBox3I);

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
* WlzCompThresh.c							*
************************************************************************/
extern WlzErrorNum		WlzCompThreshold(
				  double *dstThrVal,
				  WlzObject *histObj,
				  WlzCompThreshType method,
				  double extraFrac);

/************************************************************************
* WlzContour.c
************************************************************************/
extern WlzContour		*WlzContourObj(
				  WlzObject *srcObj,
				  WlzContourMethod ctrMtd,
				  double ctrVal,
				  double ctrWth,
				  WlzErrorNum *dstErr);
extern WlzContour		*WlzContourObjGrd(
				  WlzObject *srcObj,
				  double ctrLo,
				  double ctrHi,
				  double ctrWth,
				  WlzErrorNum *dstErr);

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

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzConvolve.c								*
************************************************************************/
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
extern WlzSimpleProperty	*WlzCopySimpleProperty(
				  WlzSimpleProperty *inPLst,
				  WlzErrorNum *dstErr);
#endif /* WLZ_EXT_BIND */

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
* WlzDilation.c								*
************************************************************************/
extern WlzObject		*WlzDilation(
				  WlzObject *obj,
				  WlzConnectType connectivity,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzDomainFill.c							*
************************************************************************/
extern WlzObject 		*WlzDomainFill(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzDomainUtils.c							*
************************************************************************/
extern WlzObjectType		WlzDomainType(
				  WlzDomain dom);
extern void			WlzBitLnSetItv(
				  UBYTE *bitLn,
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
				  UBYTE *bitLn,
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

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzError.c								*
************************************************************************/
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
#endif /* !WLZ_EXT_BIND */

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
* WlzGeoModel.c
************************************************************************/
#ifndef WLZ_EXT_BIND
/* Resource callback function list manipulation. */
extern WlzErrorNum	WlzGMModelAddResCb(
			  WlzGMModel *model,
			  WlzGMCbFn fn,
			  void *data);
extern void		WlzGMModelRemResCb(
			  WlzGMModel *model,
			  WlzGMCbFn fn,
			  void *data);
/* Creation  of geometric modeling elements */
extern WlzGMModel	*WlzGMModelNew(
			  WlzGMModelType modType,
			  int blkSz,
			  int vHTSz,
			  WlzErrorNum *dstErr);
extern WlzGMShell	*WlzGMModelNewS(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
extern WlzGMFace	*WlzGMModelNewF(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
extern WlzGMLoopT	*WlzGMModelNewLT(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
extern WlzGMEdge      *WlzGMModelNewE(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
extern WlzGMEdgeT     *WlzGMModelNewET(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
extern WlzGMVertex     *WlzGMModelNewV(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
extern WlzGMVertexT    *WlzGMModelNewVT(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
extern WlzGMModel	*WlzGMModelCopy(
			  WlzGMModel *gM,
			  WlzErrorNum *dstErr);
/* Freeing  of geometric modeling elements */
extern WlzErrorNum	WlzGMModelFree(
			  WlzGMModel *model);
extern WlzErrorNum      WlzGMModelFreeS(
			  WlzGMModel *model,
			  WlzGMShell *shell);
extern WlzErrorNum      WlzGMModelFreeF(
			  WlzGMModel *model,
			  WlzGMFace *loop);
extern WlzErrorNum      WlzGMModelFreeLT(
			  WlzGMModel *model,
			  WlzGMLoopT *loopT);
extern WlzErrorNum      WlzGMModelFreeE(
			  WlzGMModel *model,
			  WlzGMEdge *edge);
extern WlzErrorNum      WlzGMModelFreeET(
			  WlzGMModel *model,
			  WlzGMEdgeT *edgeT);
extern WlzErrorNum	WlzGMModelFreeDT(
			  WlzGMModel *model,
			  WlzGMDiskT *diskT);
extern WlzErrorNum      WlzGMModelFreeV(
			  WlzGMModel *model,
			  WlzGMVertex *vertex);
extern WlzErrorNum      WlzGMModelFreeVT(
			  WlzGMModel *model,
			  WlzGMVertexT *vertexT);
/* Deletion of geometric modeling elements along with children and any parents
 * that depend solely on the element being deleted. */
extern WlzErrorNum     	WlzGMModelDeleteV(
			  WlzGMModel *model,
			   WlzGMVertex *dV);
extern WlzErrorNum     	WlzGMModelDeleteE(
			  WlzGMModel *model,
			  WlzGMEdge *dE);
extern WlzErrorNum     	WlzGMModelDeleteF(
			  WlzGMModel *model,
			  WlzGMFace *dL);
extern WlzErrorNum	WlzGMModelDeleteS(
			  WlzGMModel *model,
			  WlzGMShell *shell);
/* Searching */
extern WlzGMVertex	*WlzGMModelMatchVertexG3D(
			  WlzGMModel *model,
			  WlzDVertex3 gPos);
extern WlzGMVertex	*WlzGMModelMatchVertexG2D(
			  WlzGMModel *model,
			  WlzDVertex2 gPos);
/* Geometry access and query functions */
extern WlzGMElemType 	WlzGMModelGetSGeomType(
			  WlzGMModel *model);
extern WlzGMElemType 	WlzGMModelGetVGeomType(
			  WlzGMModel *model);
extern WlzErrorNum	WlzGMShellSetGBB3D(
			  WlzGMShell *shell,
			  WlzDBox3 bBox);
extern WlzErrorNum	WlzGMShellGetGBB3D(
			  WlzGMShell *shell,
			  WlzDBox3 *bBox);
extern WlzErrorNum	WlzGMShellGetGBBV3D(
			  WlzGMShell *shell,
			  double *vol);
extern WlzErrorNum	WlzGMShellGetGBB3D(
                          WlzGMShell *shell,
			  WlzDBox3 *bBox);
extern WlzErrorNum	WlzGMShellGetGBB2D(
                          WlzGMShell *shell,
			  WlzDBox2 *bBox);
extern WlzErrorNum	WlzGMShellSetGBB2D(
			  WlzGMShell *shell,
			  WlzDBox2 bBox);
extern WlzErrorNum	WlzGMShellSetGBB2D(
			  WlzGMShell *shell,
			  WlzDBox2 bBox);
extern WlzErrorNum	WlzGMShellSetG3D(
			  WlzGMShell *shell,
			  int nPnt,
			  WlzDVertex3 *pos);
extern WlzErrorNum	WlzGMShellSetG2D(
			  WlzGMShell *shell,
			  int nPnt,
			  WlzDVertex2 *pos);
extern int		WlzGMShellGInBB3D(
			  WlzGMShell *shell,
			  WlzDVertex3 pos);
extern int		WlzGMShellGInBB2D(
			  WlzGMShell *shell,
			  WlzDVertex2 pos);
extern WlzErrorNum	WlzGMModelSetSG(
			  WlzGMModel *model);
extern WlzErrorNum	WlzGMShellUpdateG3D(
			  WlzGMShell *shell,
			  WlzDVertex3 pos);
extern WlzErrorNum	WlzGMShellUpdateG2D(
			  WlzGMShell *shell,
			  WlzDVertex2 pos);
extern WlzErrorNum	WlzGMShellUpdateGBB3D(
			  WlzGMShell *shell,
			  WlzDBox3 bBox);
extern WlzErrorNum	WlzGMShellUpdateGBB2D(
			  WlzGMShell *shell,
			  WlzDBox2 bBox);
extern WlzErrorNum	WlzGMVertexSetG3D(
			  WlzGMVertex *vertex,
			  WlzDVertex3 pos);
extern WlzErrorNum	WlzGMVertexSetG2D(
			  WlzGMVertex *vertex,
			  WlzDVertex2 pos);
extern WlzErrorNum	WlzGMShellDndateG2D(
			  WlzGMShell *shell,
			  WlzDVertex2 pos);
extern WlzErrorNum	WlzGMShellDndateG3D(
			  WlzGMShell *shell,
			  WlzDVertex3 pos);
extern WlzErrorNum	WlzGMShellComputeGBB(
			  WlzGMShell *shell);
extern WlzErrorNum	WlzGMVertexGetG3D(
			  WlzGMVertex *vertex,
			  WlzDVertex3 *dstPos);
extern WlzErrorNum	WlzGMVertexGetG2D(
			  WlzGMVertex *vertex,
			  WlzDVertex2 *dstPos);
extern WlzDVertex3	WlzGMVertexCmp3D(
			  WlzGMVertex *vertex,
			  WlzDVertex3 pos);
extern WlzDVertex2	WlzGMVertexCmp2D(
			  WlzGMVertex *vertex,
			  WlzDVertex2 pos);
extern int		WlzGMVertexCmpSign3D(
			  WlzGMVertex *vertex,
			  WlzDVertex3 pos);
extern int		WlzGMVertexCmpSign2D(
			  WlzGMVertex *vertex,
			  WlzDVertex2 pos);
extern double		WlzGMVertexDistSq3D(
			  WlzGMVertex *vertex,
			  WlzDVertex3 pos);
extern double		WlzGMVertexDistSq2D(
			  WlzGMVertex *vertex,
			  WlzDVertex2 pos);
extern double		WlzGMVertexShellDist(
			  WlzGMVertex *v0,
			  WlzGMVertex *v1,
			  WlzErrorNum *dstErr);
extern WlzDVertex3	WlzGMVertexNormal3D(
			  WlzGMModel *model,
			  WlzGMVertex *gV,
			  int *sVBufSz,
			  WlzGMVertex ***sVBuf,
			  WlzErrorNum *dstErr);
extern WlzGMVertex 	*WlzGMModelLoopTMaxMinCurv2D(
			  WlzGMLoopT *gLT,
			  int minLen,
			  int lnLen,
			  WlzBinaryOperatorType mOrM,
			  double *dstAlg);
/* Model access and testing */
extern WlzErrorNum 	WlzGMModelTypeValid(
			WlzGMModelType type);
extern int		WlzGMModelGetDimension(
			  WlzGMModel *model,
			  WlzErrorNum *dstErr);
/* Topology validity checks (useful for debugging) */
extern WlzErrorNum	WlzGMVerifyModel(
			  WlzGMModel *model,
			  WlzGMElemP *dstElmP);
extern WlzErrorNum	WlzGMVerifyShell(
			  WlzGMShell *shell,
			  WlzGMElemP *dstElmP);
extern WlzErrorNum	WlzGMVerifyLoopT(
			  WlzGMLoopT *loopT,
			  WlzGMElemP *dstElmP);
/* Topology query */
extern WlzGMEdge	**WlzGMModelFindNMEdges(
			  WlzGMModel *model,
			  int *dstNMCnt,
			  WlzErrorNum *dstErr);
extern WlzGMLoopT	*WlzGMEdgeTCommonLoopT(
			  WlzGMEdgeT *eT0,
			  WlzGMEdgeT *eT1);
extern WlzGMVertex	*WlzGMEdgeCommonVertex(
			  WlzGMEdge *eE0,
			  WlzGMEdge *eE1);
extern WlzGMVertex	*WlzGMEdgeCommonVertexGetDiskTs(
			  WlzGMEdge *eE0,
			  WlzGMEdge *eE1,
			  WlzGMDiskT **dstDT0,
			  WlzGMDiskT **dstDT1);
extern WlzGMDiskT	*WlzGMEdgeCommonDiskT(
			  WlzGMEdge *eE0,
			  WlzGMEdge *eE1);
extern WlzGMShell	*WlzGMEdgeGetShell(
			  WlzGMEdge *eE);
extern WlzGMFace	*WlzGMEdgeCommonFace(
			  WlzGMEdge *eE0,
			  WlzGMEdge *eE1);
extern WlzGMEdge	*WlzGMVertexCommonEdge(
			  WlzGMVertex *eV0,
			  WlzGMVertex *eV1);
extern WlzGMShell	*WlzGMVertexCommonShell(
			  WlzGMVertex *eV0,
			  WlzGMVertex *eV1);
extern WlzGMShell	*WlzGMVertexGetShell(
			  WlzGMVertex *eV);
/* Model list management */
extern void	   	WlzGMVertexTAppend(
			  WlzGMVertexT *eVT,
			  WlzGMVertexT *nVT);
extern void		WlzGMDiskTAppend(
			  WlzGMDiskT *eDT,
			  WlzGMDiskT *nDT);
extern void		WlzGMEdgeTAppend(
			  WlzGMEdgeT *eET,
			  WlzGMEdgeT *nET);
extern void		WlzGMEdgeTInsert(
			  WlzGMEdgeT *eET,
			  WlzGMEdgeT *nET);
extern void		WlzGMEdgeTInsertRadial(
			  WlzGMEdgeT *nET);
extern void	   	WlzGMLoopTAppend(
			  WlzGMLoopT *pS,
			  WlzGMLoopT *nLT);
extern void	   	WlzGMShellAppend(
			  WlzGMShell *pS,
			  WlzGMShell *newS);
extern void		WlzGMEdgeTUnlink(
			  WlzGMEdgeT *dLT);
extern void		WlzGMVertexTUnlink(
			  WlzGMVertexT *dLT);
extern void		WlzGMDiskTUnlink(
			  WlzGMDiskT *dLT);
extern void		WlzGMDiskTJoin(
			  WlzGMDiskT *gDT0,
			  WlzGMDiskT *gDT1);
extern void		WlzGMLoopTUnlink(
			  WlzGMLoopT *dLT);
extern void		WlzGMShellUnlink(
			  WlzGMShell *dS);
extern void		WlzGMShellJoinAndUnlink(
			  WlzGMShell *eShell,
			  WlzGMShell *dShell);
extern WlzGMResIdxTb	*WlzGMModelResIdx(
			  WlzGMModel *model,
			  unsigned int eMsk,
			  WlzErrorNum *dstErr);
extern void		WlzGMModelResIdxFree(
			  WlzGMResIdxTb *resIdxTb);
extern WlzErrorNum	WlzGMModelRehashVHT(
			  WlzGMModel *model,
			  int vHTSz);
/* Model construction */
extern WlzErrorNum	WlzGMModelConstructSimplex3D(
			  WlzGMModel *model,
			  WlzDVertex3 *pos);
extern WlzErrorNum	WlzGMModelConstructSimplex2D(
			  WlzGMModel *model,
			  WlzDVertex2 *pos);
/* Model Features */
extern int		WlzGMShellSimplexCnt(
			  WlzGMShell *gShell);

#endif /* !WLZ_EXT_BIND */

/************************************************************************
* WlzGeoModelFilters.c
************************************************************************/
extern WlzErrorNum		WlzGMFilterRmSmShells(
				  WlzGMModel *model,
				  int maxElm);
extern WlzErrorNum		WlzGMFilterGeomLP(
				  WlzGMModel *model,
				  double kPB,
				  double kSB,
				  double dPB,
				  double dSB,
				  int maxItr);
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
				  int nItr);

/************************************************************************
* WlzGeometry.c								*
************************************************************************/
#ifndef WLZ_EXT_BIND
extern int			WlzGeomTriangleCircumcentre(
				  WlzDVertex2 *dstCcVx,
				  WlzDVertex2 vx0,
				  WlzDVertex2 vx1,
				  WlzDVertex2 vx2);
#endif /* !WLZ_EXT_BIND */
extern int			WlzGeomVxInTriangle(
				  WlzDVertex2 vx0,
				  WlzDVertex2 vx1,
			          WlzDVertex2 vx2,
				  WlzDVertex2 vxP);
extern int 			WlzGeomInTriangleCircumcircle(
				  WlzDVertex2 vx0,
				  WlzDVertex2 vx1,
			          WlzDVertex2 vx2,
				  WlzDVertex2 gVx);
extern double			WlzGeomTriangleSnArea2(
				  WlzDVertex2 vx0,
				  WlzDVertex2 vx1,
				  WlzDVertex2 vx2);
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
extern int			WlzGeomPlaneAABBIntersect(
				  double a,
				  double b,
				  double c,
				  double d,
				  WlzDBox3 box);
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
				  WlzObject *obj);

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
* WlzGreySetRange.c							*
************************************************************************/
extern WlzErrorNum		WlzGreySetRange(
				  WlzObject *obj,
				  WlzPixelV min,
				  WlzPixelV max,
				  WlzPixelV Min,
				  WlzPixelV Max);

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
* WlzInsideDomain.c							*
************************************************************************/
extern int 			WlzInsideDomain(
				  WlzObject *obj,
				  double plane,
				  double line,
				  double kol,
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
extern int WlzHasIntersection(WlzObject	*obj1,
			      WlzObject	*obj2,
			      WlzErrorNum	*dstErr);
/************************************************************************
* WlzIntervalCount.c							*
************************************************************************/
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
* WlzJavaUtils.c							*
************************************************************************/
extern int 		        WlzObjGetType(WlzObject *obj);
extern int			WlzGetIntervalDomainLine1(WlzObject *obj);
extern int			WlzGetIntervalDomainLastLn(WlzObject *obj);
extern int			WlzGetIntervalDomainKol1(WlzObject *obj);
extern int			WlzGetIntervalDomainLastKl(WlzObject *obj);
extern int			WlzGetPlaneDomainPlane1(WlzObject *obj);
extern int			WlzGetPlaneDomainLastPl(WlzObject *obj);
extern int			WlzGetPlaneDomainLine1(WlzObject *obj);
extern int			WlzGetPlaneDomainLastLn(WlzObject *obj);
extern int			WlzGetPlaneDomainKol1(WlzObject *obj);
extern int			WlzGetPlaneDomainLastKl(WlzObject *obj);

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

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzMakeProperties.c							*
************************************************************************/
extern WlzSimpleProperty	*WlzMakeSimpleProperty(
				  int size,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzFreeSimpleProperty(
				  WlzSimpleProperty *prop);
extern WlzEMAPProperty          *WlzMakeEMAPProperty(
                                  WlzEMAPPropertyType	type,
				  int			theilerStage,
				  char			*modelName,
				  char			*version,
				  char			*fileName,
				  char			*comment,
				  WlzErrorNum *dstErr);
extern WlzErrorNum              WlzChangeEMAPProperty(
                                  WlzEMAPProperty	*prop,
                                  WlzEMAPPropertyType	type,
				  int			theilerStage,
				  char			*modelName,
				  char			*version,
				  char			*fileName,
				  char			*comment);
extern WlzErrorNum		WlzFreeEMAPProperty(
				  WlzEMAPProperty *prop);
extern WlzErrorNum		WlzFreeProperty(WlzProperty prop);
extern WlzErrorNum		WlzFreePropertyList(AlcDLPList	*plist);
extern void			WlzFreePropertyListEntry(void *prop);
extern WlzProperty 		WlzGetProperty(
                                  AlcDLPList		*plist,
				  WlzObjectType		type,
				  WlzErrorNum 		*dstErr);
extern WlzErrorNum              WlzRemoveProperty(
                                  AlcDLPList		*plist,
				  WlzProperty		prop);
#endif /* WLZ_EXT_BIND */

#ifndef WLZ_EXT_BIND
/************************************************************************
 * WlzMakeStructs.c
 ************************************************************************/
extern WlzObject 		*WlzMakeMain(
				  WlzObjectType type,
				  WlzDomain domain,
				  WlzValues values,
				  AlcDLPList *prop,
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
extern WlzPolygonDomain		*WlzMakePolyDmn(
				  WlzObjectType type,
				  WlzIVertex2 *vtxs,
				  int nvertices,
				  int max_nverts,
				  int copy,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzMakeRect(
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  WlzGreyType pixeltype,
				  int *greyptr,
				  WlzPixelV bckgrnd,
				  AlcDLPList *plist,
				  WlzObject *assoc_obj,
				  WlzErrorNum *dstErr);
extern WlzObject		*WlzMakeEmpty(
				  WlzErrorNum *dstErr);
extern WlzBoundList		*WlzMakeBoundList(
				  WlzObjectType type,
				  int wrap,
				  WlzPolygonDomain *poly,
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
				  extern WlzObject *WlzNewGrey(
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
#endif /* WLZ_EXT_BIND */

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzMatchICP.c
************************************************************************/
extern WlzErrorNum		WlzMatchICPObjs(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  int *dstNMatch,
				  WlzVertexP *dstTMatch,
				  WlzVertexP *dstSMatch,
				  int maxItr,
				  int minSpx,
				  int brkFlg,
				  double maxDisp,
				  double maxAng,
				  double maxDeform,
				  int matchImpNN,
				  double matchImpThr);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* WlzMeshTransform.c							*
************************************************************************/
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
extern WlzDVertex2 		WlzMeshTransformVtx(WlzDVertex2 vtx,
						    WlzMeshTransform *mesh,
						    WlzErrorNum *dstErr);
extern double 			WlzClassValCon4(
                                  double	*gVals,
				  double	xOffset,
				  double	yOffset);
#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzMeshUtils.c							*
************************************************************************/
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
				  int *existsFlg,
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
extern int			WlzPolyCrossings(WlzIVertex2	vtx,
						 WlzPolygonDomain *pgdm,
						 WlzErrorNum	  *dstErr);
extern int			WlzInsidePolyEO(WlzIVertex2	vtx,
						WlzPolygonDomain  *pgdm,
						WlzErrorNum	  *dstErr);

/************************************************************************
* WlzPolyEquispace.c							*
************************************************************************/
extern WlzPolygonDomain		*WlzPolyEquispace(
                                  WlzPolygonDomain	*poly,
				  int			wrap,
				  double		spacing,
				  int			keepOrigVtxs,
				  WlzErrorNum		*dstErr);
extern double 			WlzPolyLength(
                                  WlzPolygonDomain	*poly,
				  int			wrap,
				  WlzErrorNum		*dstErr);
/************************************************************************
* WlzPolyDecimate.c							*
************************************************************************/
extern WlzPolygonDomain  	*WlzPolyDecimate(
                                  WlzPolygonDomain	*poly,
				  int			wrap,
				  double		maxDist,
				  WlzErrorNum		*dstErr);
extern WlzBoundList 		*WlzBoundDecimate(
                                  WlzBoundList	*bound,
				  double 	maxDist,
				  WlzErrorNum	*dstErr);
/************************************************************************
* WlzPolyReverse.c							*
************************************************************************/
extern WlzPolygonDomain 	*WlzPolyReverse(
                                  WlzPolygonDomain 	*poly,
				  WlzErrorNum		*dstErr);
/************************************************************************
* WlzPolySmooth.c							*
************************************************************************/
extern WlzPolygonDomain		*WlzPolySmooth(
                                  WlzPolygonDomain	*poly,
				  int			wrap,
				  int 			iterations,
				  WlzErrorNum		*dstErr);
extern WlzBoundList 		*WlzBoundSmooth(
                                  WlzBoundList	*bound,
				  int 		iterations,
				  WlzErrorNum	*dstErr);
/************************************************************************
* WlzPrinicipalAngle.c
************************************************************************/
extern double 			WlzPrincipalAngle(
				  WlzObject *srcObj,
				  WlzDVertex2 cMass,
				  int binObjFlg,
				  WlzErrorNum *dstErr);

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

/************************************************************************
* WlzReadObj.c								*
************************************************************************/
extern WlzObject		*WlzReadObj(
				  FILE *fp,
			          WlzErrorNum *dstErr);

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
				  int *dstConv,
				  double *dstCCor,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzRegICP.c
************************************************************************/
extern WlzAffineTransform	*WlzRegICPObjs(
				  WlzObject *tObj,
				  WlzObject *sObj,
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv,
				  int *dstItr,
				  int maxItr,
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
				  WlzAffineTransform *initTr,
				  WlzTransformType trType,
				  int *dstConv,
				  int *dstItr,
				  int maxItr,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform	*WlzRegICPTreeAndVertices(
				  AlcKDTTree *tree,
				  WlzTransformType trType,
				  WlzVertexType vType,
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

/************************************************************************
* WlzScalarArithmeticOp.c						*
************************************************************************/
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

/************************************************************************
* WlzScalarBinaryOp.c							*
************************************************************************/
extern WlzErrorNum 		WlzScalarBinaryOp(
				  WlzObject *o1,
				  WlzPixelV pval,
				  WlzObject *o3,
				  WlzBinaryOperatorType op);

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzSepTrans.c								*
************************************************************************/
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

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzSeqPar.c								*
************************************************************************/
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
* WlzShift.c								*
************************************************************************/
extern WlzObject		*WlzShadeCorrect(
				  WlzObject *srcObj,
				  WlzObject *shdObj,
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

/************************************************************************
* WlzStringTypes.c							*
************************************************************************/
extern const char 		*WlzStringFromObjType(
				  WlzObject *obj,
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
extern WlzTransformType 	WlzStringToTransformType(
				  const char *tStr,
				  WlzErrorNum *dstErr);
extern const char 		*WlzStringFromTransformType(
				  WlzTransformType tType,
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
extern const char 		*WlzStringFromErrorNum(
				  WlzErrorNum gvnErr,
			          const char **dstMsgStr);
extern WlzErrorNum  		WlzStringToErrorNum(
				  const char *gStr);

#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzStringUtils.c							*
************************************************************************/
extern int 			WlzStringMatchValue(
				  int *dstValue,
				  const char *targetStr,
				  const char *testStr,
				  ...);
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
* WlzThreshold.c							*
************************************************************************/
extern WlzObject 		*WlzThreshold(
				  WlzObject *obj,
				  WlzPixelV threshV,
				  WlzThresholdType highlow,
				  WlzErrorNum *dstErr);

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
extern WlzGreyType		WlzGreyTypeFromObj(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
extern WlzDVertex3		WlzVozelSz(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
#ifndef WLZ_EXT_BIND
/************************************************************************
* WlzValueUtils.c							*
************************************************************************/
extern void			WlzValueSetInt(
				  int *vec,
				  int value,
				  int count);
extern void			WlzValueSetShort(
				  short *vec,
				  short value,
				  int count);
extern void			WlzValueSetUByte(
				  UBYTE *vec,
				  UBYTE value,
				  int count);
extern void			WlzValueSetFloat(
				  float *vec,
				  float value,
				  int count);
extern void			WlzValueSetDouble(
				  double *vec,
				  double value,
				  int count);
extern void			WlzValueSetDVertex(
				  WlzDVertex2 *vec,
				  WlzDVertex2 value,
				  int count);
extern void			WlzValueSetFVertex(
				  WlzFVertex2 *vec,
				  WlzFVertex2 value,
				  int count);
extern void			WlzValueSetIVertex(
				  WlzIVertex2 *vec,
				  WlzIVertex2 value,
				  int count);
extern void			WlzValueSetGrey(
				  WlzGreyP vec,
				  int vecOff,
				  WlzGreyV value,
				  WlzGreyType gType,
				  int count);
extern void			WlzValueClampGreyIntoGrey(
				  WlzGreyP dst,
				  int dstOff,
				  WlzGreyType dstType,
				  WlzGreyP src,
				  int srcOff,
				  WlzGreyType srcType,
				  int count);
extern void			WlzValueClampIntIntoShort(
				  short *dst,
				  int *src,
				  int count);
extern void			WlzValueClampIntIntoUByte(
				  UBYTE *dst,
				  int *src,
				  int count);
extern void			WlzValueClampShortIntoUByte(
				  UBYTE *dst,
				  short *src,
				  int count);
extern void			WlzValueClampFloatIntoInt(
				  int *dst,
				  float *src,
				  int count);
extern void			WlzValueClampFloatIntoShort(
				  short *dst,
				  float *src,
				  int count);
extern void			WlzValueClampFloatIntoUByte(
				  UBYTE *dst,
				  float *src,
				  int count);
extern void			WlzValueClampDoubleIntoInt(
				  int *dst,
				  double *src,
				  int count);
extern void			WlzValueClampDoubleIntoUByte(
				  UBYTE *dst,
				  double *src,
				  int count);
extern void			WlzValueClampDoubleIntoFloat(
				  float *dst,
				  double *src,
				  int count);
extern void			WlzValueClampIntToShort(
				  int *vec,
				  int count);
extern void			WlzValueClampIntToUByte(
				  int *vec,
				  int count);
extern void			WlzValueClampShortToUByte(
				  short *vec,
				  int count);
extern void			WlzValueClampDoubleToInt(
				  double *vec,
				  int count);
extern void			WlzValueClampDoubleToShort(
				  double *vec,
				  int count);
extern void			WlzValueClampDoubleToUByte(
				  double *vec,
				  int count);
extern void			WlzValueClampDoubleToFloat(
				  double *vec,
				  int count);
extern void			WlzValueClampFloatToInt(
				  float *vec,
				  int count);
extern void			WlzValueClampFloatToShort(
				  float *vec,
				  int count);
extern void			WlzValueClampFloatToUByte(
				  float *vec,
				  int count);
extern void			WlzValueMaskIntToShort(
				  int *vec,
				  int count);
extern void			WlzValueMaskIntToUByte(
				  int *vec,
				  int count);
extern void			WlzValueMaskShortToUByte(
				  short *vec,
				  int count);
extern void			WlzValueCopyIntToInt(
				  int *dst,
				  int *src,
				  int count);
extern void			WlzValueCopyIntToShort(
				  short *dst,
				  int *src,
				  int count);
extern void			WlzValueCopyIntToUByte(
				  UBYTE *dst,
				  int *src,
				  int count);
extern void			WlzValueCopyIntToFloat(
				  float *dst,
				  int *src,
				  int count);
extern void			WlzValueCopyIntToDouble(
				  double *dst,
				  int *src,
				  int count);
extern void			WlzValueCopyShortToInt(
				  int *dst,
				  short *src,
				  int count);
extern void			WlzValueCopyShortToShort(
				  short *dst,
				  short *src,
				  int count);
extern void			WlzValueCopyShortToUByte(
				  UBYTE *dst,
				  short *src,
				  int count);
extern void			WlzValueCopyShortToFloat(
				  float *dst,
				  short *src,
				  int count);
extern void			WlzValueCopyShortToDouble(
				  double *dst,
				  short *src,
				  int count);
extern void			WlzValueCopyUByteToInt(
				  int *dst,
				  UBYTE *src,
				  int count);
extern void			WlzValueCopyUByteToShort(
				  short *dst,
				  UBYTE *src,
				  int count);
extern void			WlzValueCopyUByteToUByte(
				  UBYTE *dst,
				  UBYTE *src,
				  int count);
extern void			WlzValueCopyUByteToFloat(
				  float *dst,
				  UBYTE *src,
				  int count);
extern void			WlzValueCopyUByteToDouble(
				  double *dst,
				  UBYTE *src,
				  int count);
extern void			WlzValueCopyFloatToInt(
				  int *dst,
				  float *src,
				  int count);
extern void			WlzValueCopyFloatToShort(
				  short *dst,
				  float *src,
				  int count);
extern void			WlzValueCopyFloatToUByte(
				  UBYTE *dst,
				  float *src,
				  int count);
extern void			WlzValueCopyFloatToFloat(
				  float *dst,
				  float *src,
				  int count);
extern void			WlzValueCopyFloatToDouble(
				  double *dst,
				  float *src,
				  int count);
extern void			WlzValueCopyDoubleToInt(
				  int *dst,
				  double *src,
				  int count);
extern void			WlzValueCopyDoubleToShort(
				  short *dst,
				  double *src,
				  int count);
extern void			WlzValueCopyDoubleToUByte(
				  UBYTE *dst,
				  double *src,
				  int count);
extern void			WlzValueCopyDoubleToFloat(
				  float *dst,
				  double *src,
				  int count);
extern void			WlzValueCopyDoubleToDouble(
				  double *dst,
				  double *src,
				  int count);
extern void			WlzValueCopyGreyToGrey(
				  WlzGreyP dst,
				  int dstOff,
				  WlzGreyType dstType,
				  WlzGreyP src,
				  int srcOff,
				  WlzGreyType srcType,
				  int count);
extern void			WlzValueCopyDVertexToDVertex(
				  WlzDVertex2 *dst,
				  WlzDVertex2 *src,
				  int count);
extern void			WlzValueCopyDVertexToFVertex(
				  WlzFVertex2 *dst,
				  WlzDVertex2 *src,
				  int count);
extern void			WlzValueCopyDVertexToIVertex(
				  WlzIVertex2 *dst,
				  WlzDVertex2 *src,
				  int count);
extern void			WlzValueCopyFVertexToDVertex(
				  WlzDVertex2 *dst,
				  WlzFVertex2 *src,
				  int count);
extern void			WlzValueCopyFVertexToFVertex(
				  WlzFVertex2 *dst,
				  WlzFVertex2 *src,
				  int count);
extern void			WlzValueCopyFVertexToIVertex(
				  WlzIVertex2 *dst,
				  WlzFVertex2 *src,
				  int count);
extern void			WlzValueCopyIVertexToDVertex(
				  WlzDVertex2 *dst,
				  WlzIVertex2 *src,
				  int count);
extern void			WlzValueCopyIVertexToFVertex(
				  WlzFVertex2 *dst,
				  WlzIVertex2 *src,
				  int count);
extern void			WlzValueCopyIVertexToIVertex(
				  WlzIVertex2 *dst,
				  WlzIVertex2 *src,
				  int count);
extern void			WlzValueCopyDVertexToDVertex3(
				  WlzDVertex3 *dst,
				  WlzDVertex3 *src,
				  int count);
extern void			WlzValueCopyDVertexToFVertex3(
				  WlzFVertex3 *dst,
				  WlzDVertex3 *src,
				  int count);
extern void			WlzValueCopyDVertexToIVertex3(
				  WlzIVertex3 *dst,
				  WlzDVertex3 *src,
				  int count);
extern void			WlzValueCopyFVertexToDVertex3(
				  WlzDVertex3 *dst,
				  WlzFVertex3 *src,
				  int count);
extern void			WlzValueCopyFVertexToFVertex3(
				  WlzFVertex3 *dst,
				  WlzFVertex3 *src,
				  int count);
extern void			WlzValueCopyFVertexToIVertex3(
				  WlzIVertex3 *dst,
				  WlzFVertex3 *src,
				  int count);
extern void			WlzValueCopyIVertexToDVertex3(
				  WlzDVertex3 *dst,
				  WlzIVertex3 *src,
				  int count);
extern void			WlzValueCopyIVertexToFVertex3(
				  WlzFVertex3 *dst,
				  WlzIVertex3 *src,
				  int count);
extern void			WlzValueCopyIVertexToIVertex3(
				  WlzIVertex3 *dst,
				  WlzIVertex3 *src,
				  int count);
extern WlzErrorNum 		WlzValueConvertPixel(
				  WlzPixelV *dstPix,
				  WlzPixelV srcPix,
				  WlzGreyType dstType);
extern int			WlzValueMedianInt(
				  int *vec,
				  int count);
extern double			WlzValueMedianDouble(
				  double *vec,
				  int count);
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
* WlzVerticies.c
************************************************************************/
#ifndef WLZ_EXT_BIND
extern WlzVertexP		WlzVerticesFromObj(
				  WlzObject *obj,
				  WlzVertexP *dstNr,
				  int *dstCnt,
				  WlzVertexType *dstType,
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
/************************************************************************
* WlzMwrAngle.c								*
************************************************************************/
extern double WlzMwrAngle(WlzObject *cvh,
			  WlzErrorNum *dstErr);


#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif	/* !WLZ_PROTO_H Don't put anything after this line */
