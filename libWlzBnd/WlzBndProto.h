#ifndef WLZBND_PROTO_H
#define WLZBND_PROTO_H
#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBndProto_h[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzBndProto.h
* \author       Bill Hill
* \date         August 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2006 Medical research Council, UK.
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
* \brief	Defines the Woolz binding library function proto
*		types.
* \ingroup	WlzBnd
* \todo         -
* \bug          None known.
*/


#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************************************************************
* Wlz2DContains.c
************************************************************************/
extern WlzErrorNum		WlzBndGreyInvert(
				  WlzObject *obj);

/************************************************************************
* WlzJavaUtils.c
************************************************************************/
extern int 			WlzBndObjGetType(
				  WlzObject *obj);
extern int 			WlzBndGetIntervalDomainLine1(
				  WlzObject *obj);
extern int 			WlzBndGetIntervalDomainLastLn(
				  WlzObject *obj);
extern int 			WlzBndGetIntervalDomainKol1(
				  WlzObject *obj);
extern int 			WlzBndGetIntervalDomainLastKl(
				  WlzObject *obj);
extern int 			WlzBndGetPlaneDomainPlane1(
				  WlzObject *obj);
extern int 			WlzBndGetPlaneDomainLastPl(
				  WlzObject *obj);
extern int 			WlzBndGetPlaneDomainLine1(
				  WlzObject *obj);
extern int 			WlzBndGetPlaneDomainLastLn(
				  WlzObject *obj);
extern int 			WlzBndGetPlaneDomainKol1(
				  WlzObject *obj);
extern int 			WlzBndGetPlaneDomainLastKl(
				  WlzObject *obj);
extern WlzObject 		*WlzBndSetGreyValues(
				  WlzObject *domainObj,
				  WlzObject *valuesObj,
				  WlzErrorNum *dstErr);
extern WlzAffineTransform 	*WlzBndGetAffineTransform(
				  WlzThreeDViewStruct *viewStr,
				  WlzErrorNum *dstErr);
extern void 			WlzBndGetTransformMatrix(
				  WlzIVertex2 *dstSizeArrayDat,
				  double  ***dstArrayDat,
				  WlzAffineTransform *trans,
				  WlzErrorNum *dstErr);

/************************************************************************
* WlzBnd3DWarp.c
************************************************************************/
extern WlzObject 		*WlzBnd3DWarpObj(
				  WlzObject *wObjS,
				  int arraySizeVec0,
				  WlzDVertex3 *arrayVec0,
				  int arraySizeVec1,
				  WlzDVertex3 *arrayVec1,
				  WlzErrorNum *dstErr);
extern WlzObject 		*WlzBnd3DWarpFile(
				  char *inFileStr,
				  char *TiePointsFileStr,
				  char *outFileStr,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzBndBasisFn.c
************************************************************************/
extern WlzBasisFn 		*WlzBndBasisFnMQ3DFromCPts(
				  int arraySizeVec0,
				  WlzDVertex3 *arrayVec0,
				  int arraySizeVec1,
				  WlzDVertex3 *arrayVec1,
				  double delta,
				  WlzErrorNum *dstErr);
extern WlzDVertex3		WlzBndBasisFnTransformVertexD3(
				  WlzBasisFnTransform *basisTr,
				  WlzDVertex3 srcVx,
				  WlzErrorNum *dstErr);
extern WlzBasisFnTransform 	*WlzBndBasisFnTrFromCPts3(
				  WlzFnType type,
				  int order,
				  int arraySizeVec0,
				  WlzDVertex3 *arrayVec0,
				  int arraySizeVec1,
				  WlzDVertex3 *arrayVec1,
				  WlzErrorNum *dstErr);
/************************************************************************
* WlzBndFunction.c
************************************************************************/
extern WlzErrorNum 		WlzSetVoxelSize(
				  WlzObject *obj,
				  double x,
				  double y,
				  double z);
extern WlzObjectType 		WlzGetObjectType(
				  WlzObject *obj);
extern WlzErrorNum 		WlzExplode(
				  int *dstExpObjCount,
				  WlzObject ***dstExpObjVecP,
				  WlzObject *srcObj);
extern WlzObject 		*WlzGetContourObj(
				  WlzObject *inObj);
extern const char 		*WlzGetPropName(
				  WlzObject *obj);
extern int 			WlzDestroyObj(
				  WlzObject *obj);
#endif	/* !WLZBND_PROTO_H Don't put anything after this line */
