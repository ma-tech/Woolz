		#ifndef WLZBND_PROTO_H
		#define WLZBND_PROTO_H
		#pragma ident "MRC HGU $Id$"
		/*!
		* \file         WlzBndProto.h
		* \author       Bill Hill
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
		* \brief	Defines the Woolz binding library function prototypes.
		* \note		See the rules in WlzType.h for prototypes.
		* \ingroup	WlzBnd
		* \todo         -
		* \bug          None known.
		*/
		#ifdef  __cplusplus
		extern "C" {
		#endif /* __cplusplus */

		/************************************************************************
		* Wlz2DContains.c							*
		************************************************************************/
		extern WlzErrorNum		WlzBndGreyInvert(
						WlzObject *obj);

		/************************************************************************
		* WlzJavaUtils.c							*
		************************************************************************/
		extern int 		        WlzBndObjGetType(WlzObject *obj);
		extern int			WlzBndGetIntervalDomainLine1(WlzObject *obj);
		extern int			WlzBndGetIntervalDomainLastLn(WlzObject *obj);
		extern int			WlzBndGetIntervalDomainKol1(WlzObject *obj);
		extern int			WlzBndGetIntervalDomainLastKl(WlzObject *obj);
		extern int			WlzBndGetPlaneDomainPlane1(WlzObject *obj);
		extern int			WlzBndGetPlaneDomainLastPl(WlzObject *obj);
		extern int			WlzBndGetPlaneDomainLine1(WlzObject *obj);
		extern int			WlzBndGetPlaneDomainLastLn(WlzObject *obj);
		extern int			WlzBndGetPlaneDomainKol1(WlzObject *obj);
		extern int			WlzBndGetPlaneDomainLastKl(WlzObject *obj);
		extern WlzObject *WlzBndSetGreyValues(WlzObject *domainObj, WlzObject *valuesObj,
								WlzErrorNum *dstErr);
		extern WlzAffineTransform *WlzBndGetAffineTransform(
								WlzThreeDViewStruct *viewStr,
								WlzErrorNum *dstErr);
		extern void WlzBndGetTransformMatrix(WlzIVertex2 *dstSizeArrayDat,
						     double  ***dstArrayDat,
						     WlzAffineTransform *trans,
						     WlzErrorNum *dstErr);

		extern WlzObject *WlzBnd3DWarpObj(WlzObject *wObjS, int arraySizeVec0, WlzDVertex3 *arrayVec0, int arraySizeVec1, WlzDVertex3 *arrayVec1, WlzErrorNum *dstErr);
		extern WlzObject *WlzBnd3DWarpFile(char *inFileStr, char *TiePointsFileStr, char *outFileStr, WlzErrorNum *dstErr);
		
		extern WlzBasisFn *WlzBndBasisFnMQ3DFromCPts(int arraySizeVec0,	WlzDVertex3 *arrayVec0,
						int arraySizeVec1, WlzDVertex3 *arrayVec1,
						double delta, WlzErrorNum *dstErr);
		extern WlzDVertex3	WlzBndBasisFnTransformVertexD3(WlzBasisFnTransform *basisTr,
							WlzDVertex3 srcVx, WlzErrorNum *dstErr);
                extern WlzBasisFnTransform *WlzBndBasisFnTrFromCPts3(WlzFnType type,
					  int order, int arraySizeVec0, WlzDVertex3 *arrayVec0,
					  int arraySizeVec1, WlzDVertex3 *arrayVec1,
					  WlzErrorNum *dstErr);
							
		#endif	/* !WLZBND_PROTO_H Don't put anything after this line */