#ifndef WLZ_MACRO_H
#define WLZ_MACRO_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzMacro.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      C pre-processor directives, eg macros.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 13-03-01 bill Added WLZ_SIGN(), WLZ_VTX_2_SIGN(), WLZ_VTX_3_SIGN(),
		WLZ_VTX_2_ABS() and WLZ_VTX_3_ABS().
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************************************************************
* Simple macros not always available elsewhere.				*
************************************************************************/
#define	WLZ_MAX(X,Y)	(((X)>=(Y))?(X):(Y))
#define	WLZ_MIN(X,Y)	(((X)<=(Y))?(X):(Y))
#define	WLZ_ABS(X)	(((X)>0)?(X):(-(X)))
#define	WLZ_NINT(X)	((int)(((X)<0)?((X)-(0.5)):((X)+(0.5))))
#define	WLZ_SIGN(X)	(((X)<0)?-1:((X)>0)?1:0)

/************************************************************************
* Pixel value clamping. Used in WlzGreyScan() and WlzConvertPix().	*
************************************************************************/
#define WLZ_CLAMP(v, min, max) (v<min ? min : v>max ? max : v)

/************************************************************************
* Math constants.							*
************************************************************************/
#define	WLZ_M_E		(2.7182818284590452354)
#define	WLZ_M_LOG2E	(1.4426950408889634074)
#define	WLZ_M_LOG10E	(0.43429448190325182765)
#define	WLZ_M_LN2	(0.69314718055994530942)
#define	WLZ_M_LN10	(2.30258509299404568402)
#define WLZ_M_PI	(3.14159265358979323846)
#define WLZ_M_PI_2	(1.57079632679489661923)
#define WLZ_M_PI_4	(0.78539816339744830961)
#define	WLZ_M_1_PI	(0.31830988618379067154)
#define	WLZ_M_2_PI	(0.63661977236758134308)
#define	WLZ_M_2_SQRTPI	(1.12837916709551257390)
#define WLZ_M_SQRT2	(1.41421356237309504880)
#define	WLZ_M_SQRT1_2	(0.70710678118654752440)
#define WLZ_M_SQRT3     (1.73205080756887729352)
#define WLZ_M_SQRT3_2   (0.86602540378443864676)

/************************************************************************
* Mesh tolerance control.						*
************************************************************************/
#define WLZ_MESH_TOLERANCE      (1.0E-04)
#define WLZ_MESH_TOLERANCE_SQ   (WLZ_MESH_TOLERANCE * WLZ_MESH_TOLERANCE)
#define WLZ_MESH_ELEM_AREA_TOLERANCE   (WLZ_M_SQRT3_2 * WLZ_MESH_TOLERANCE_SQ)

/************************************************************************
* Byte packed bitmap macros
************************************************************************/
#define WLZ_BIT_SET(A,B)	*((A)+((B)>>3))|=(1<<((B)&7))
#define WLZ_BIT_GET(A,B)	(*((A)+((B)>>3))&(1<<((B)&7)))

/************************************************************************
* Colour values macros
************************************************************************/
#define WLZ_RGBA_RED_GET(V) 	((((UINT) V) & 0xff)>>0)
#define WLZ_RGBA_GREEN_GET(V) 	((((UINT) V) & 0xff00)>>8)
#define WLZ_RGBA_BLUE_GET(V) 	((((UINT) V) & 0xff0000)>>16)
#define WLZ_RGBA_ALPHA_GET(V) 	((((UINT) V) & 0xff000000)>>24)

#define WLZ_RGBA_RED_SET(V,C)	(V = (((UINT) V)&0xffffff00) | \
                                               (((UINT) C)&0xff))
#define WLZ_RGBA_GREEN_SET(V,C)	(V = (((UINT) V)&0xffff00ff) | \
                                               ((((UINT) C)&0xff)<<8))
#define WLZ_RGBA_BLUE_SET(V,C)	(V = (((UINT) V)&0xff00ffff) | \
                                               ((((UINT) C)&0xff)<<16))
#define WLZ_RGBA_ALPHA_SET(V,C)	(V = (((UINT) V)&0x00ffffff) | \
                                               ((((UINT) C)&0xff)<<24))
#define WLZ_RGBA_RGBA_SET(V,R,G,B,A) (V = ((((UINT) R)&0xff) + \
                                           ((((UINT) G)&0xff)<<8) + \
                                           ((((UINT) B)&0xff)<<16) + \
                                           ((((UINT) A)&0xff)<<24))) 

#define WLZ_RGBA_MEAN(V) 	((WLZ_RGBA_RED_GET(V) + \
                                  WLZ_RGBA_GREEN_GET(V) + \
                                  WLZ_RGBA_BLUE_GET(V))/3.0)

#define WLZ_RGBA_MODULUS_2(V)(WLZ_RGBA_RED_GET(V)*WLZ_RGBA_RED_GET(V) + \
                              WLZ_RGBA_GREEN_GET(V)*WLZ_RGBA_GREEN_GET(V) + \
                              WLZ_RGBA_BLUE_GET(V)*WLZ_RGBA_BLUE_GET(V))

#define WLZ_RGBA_MODULUS(V) (sqrt((double) WLZ_RGBA_MODULUS_2(V)))

/************************************************************************
* Vertex macros.							*
************************************************************************/
/* WLZ_VTX_2_ABS: Absolute value of a Wlz[DFI]Vertex2 */
#define WLZ_VTX_2_ABS(A,X) \
		(A).vtX = WLZ_ABS((X).vtX), \
		(A).vtY = WLZ_ABS((X).vtY)
/* WLZ_VTX_3_ABS: Absolute value of a Wlz[DFI]Vertex3 */
#define WLZ_VTX_3_ABS(A,X) \
		(A).vtX = WLZ_ABS((X).vtX), \
		(A).vtY = WLZ_ABS((X).vtY), \
		(A).vtZ = WLZ_ABS((X).vtZ)
/* WLZ_VTX_2_SIGN: Sign of a Wlz[DFI]Vertex2 */
#define WLZ_VTX_2_SIGN(S,X) \
		(S).vtX = WLZ_SIGN((X).vtX), \
		(S).vtY = WLZ_SIGN((X).vtY)
/* WLZ_VTX_3_SIGN: Sign of a Wlz[DFI]Vertex3 */
#define WLZ_VTX_3_SIGN(S,X) \
		(S).vtX = WLZ_SIGN((X).vtX), \
		(S).vtY = WLZ_SIGN((X).vtY), \
		(S).vtZ = WLZ_SIGN((X).vtZ)
/* WLZ_VTX_2_SET: Set Wlz[DFI]Vertex2 */
#define WLZ_VTX_2_SET(U,X,Y) \
		(U).vtX = (X), \
		(U).vtY = (Y)

/* WLZ_VTX_3_SET: Set Wlz[DFI]Vertex3 */
#define WLZ_VTX_3_SET(U,X,Y,Z) \
		(U).vtX = (X), \
		(U).vtY = (Y), \
		(U).vtZ = (Z)

/* WLZ_VTX_2_ADD: Add two Wlz[DFI]Vertex2's */
#define WLZ_VTX_2_ADD(U,V,W) \
		(U).vtX = (V).vtX + (W).vtX, \
		(U).vtY = (V).vtY + (W).vtY

/* WLZ_VTX_3_ADD: Add two Wlz[DFI]Vertex3's */
#define WLZ_VTX_3_ADD(U,V,W) \
		(U).vtX = (V).vtX + (W).vtX, \
		(U).vtY = (V).vtY + (W).vtY, \
		(U).vtZ = (V).vtZ + (W).vtZ

/* WLZ_VTX_2_ADD: Add three Wlz[DFI]Vertex2's */
#define WLZ_VTX_2_ADD3(U,V,W,X) \
		(U).vtX = (V).vtX + (W).vtX + (X).vtX, \
		(U).vtY = (V).vtY + (W).vtY + (X).vtY

/* WLZ_VTX_3_ADD: Add three Wlz[DFI]Vertex3's */
#define WLZ_VTX_3_ADD3(U,V,W,X) \
		(U).vtX = (V).vtX + (W).vtX + (X).vtX, \
		(U).vtY = (V).vtY + (W).vtY + (X).vtY, \
		(U).vtZ = (V).vtZ + (W).vtZ + (X).vtZ

/* WLZ_VTX_2_SUB: Subtract two Wlz[DFI]Vertex2's */
#define WLZ_VTX_2_SUB(U,V,W) \
		(U).vtX = (V).vtX - (W).vtX, \
		(U).vtY = (V).vtY - (W).vtY

/* WLZ_VTX_3_SUB: Subtract two Wlz[DFI]Vertex3's */
#define WLZ_VTX_3_SUB(U,V,W) \
		(U).vtX = (V).vtX - (W).vtX, \
		(U).vtY = (V).vtY - (W).vtY, \
		(U).vtZ = (V).vtZ - (W).vtZ

/* WLZ_VTX_2_SCALE: Scale a Wlz[DFI]Vertex2 */
#define WLZ_VTX_2_SCALE(U,V,C) \
		(U).vtX = (V).vtX * (C), \
		(U).vtY = (V).vtY * (C)

/* WLZ_VTX_3_SCALE: Scale a Wlz[DFI]Vertex3 */
#define WLZ_VTX_3_SCALE(U,V,C) \
		(U).vtX = (V).vtX * (C), \
		(U).vtY = (V).vtY * (C), \
		(U).vtZ = (V).vtZ * (C)

/* WLZ_VTX_2_DOT: Dot (scalar) product of two Wlz[DFI]Vertex2's */
#define WLZ_VTX_2_DOT(V,W) \
		((V).vtX * ((W).vtX) + \
		((V).vtY *  (W).vtY))

/* WLZ_VTX_3_DOT: Dot (scalar) product of two Wlz[DFI]Vertex3s */
#define WLZ_VTX_3_DOT(V,W) \
		((V).vtX * ((W).vtX) + \
		((V).vtY *  (W).vtY) + \
		((V).vtZ *  (W).vtZ))

/* WLZ_VTX_3_CROSS: Cross (vector) product of two Wlz[DFI]Vertex3s */
#define WLZ_VTX_3_CROSS(U,V,W) \
		(U).vtX = ((V).vtY * (W).vtZ) - ((W).vtY * (V).vtZ), \
		(U).vtY = ((V).vtZ * (W).vtX) - ((W).vtZ * (V).vtX), \
		(U).vtZ = ((V).vtX * (W).vtY) - ((W).vtX * (V).vtY)

/* WLZ_VTX_2_LENGTH: Square of length of a Wlz[DFI]Vertex2 */
#define WLZ_VTX_2_SQRLEN(U) \
		(((U).vtX * (U).vtX) + \
		 ((U).vtY * (U).vtY))

/* WLZ_VTX_3_LENGTH: Square of length of a Wlz[DFI]Vertex3 */
#define WLZ_VTX_3_SQRLEN(U) \
		(((U).vtX * (U).vtX) + \
		 ((U).vtY * (U).vtY) + \
		 ((U).vtZ * (U).vtZ))

/* WLZ_VTX_2_LENGTH: Length of a Wlz[DFI]Vertex2 */
#define WLZ_VTX_2_LENGTH(U) \
		(sqrt(WLZ_VTX_2_SQRLEN(U)))

/* WLZ_VTX_3_LENGTH: Length of a Wlz[DFI]Vertex3 */
#define WLZ_VTX_3_LENGTH(U) \
		(sqrt(WLZ_VTX_3_SQRLEN(U)))

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif
