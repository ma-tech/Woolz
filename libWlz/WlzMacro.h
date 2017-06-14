#ifndef WLZMACRO_H
#define WLZMACRO_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMacro_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzMacro.h
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
* \brief	Woolz C pre-processor directives, eg macros.
* \ingroup	Wlz
*/

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#define WLZ_VERSION	PACKAGE_VERSION

/************************************************************************
* Simple macros not always available elsewhere.				*
************************************************************************/
#define	WLZ_MAX(X,Y)	(((X)>=(Y))?(X):(Y))
#define	WLZ_MIN(X,Y)	(((X)<=(Y))?(X):(Y))
#define	WLZ_ABS(X)	(((X)>0)?(X):(-(X)))
#define	WLZ_NINT(X)	((int)(((X)<0)?((X)-(0.5)):((X)+(0.5))))
#define	WLZ_SIGN(X)	(((X)<0)?-1:((X)>0)?1:0)

/*!
* \def		WLZ_SWAP(T,A,B)
* \brief	Swaps the values of the variables A and B using the
* 		temporary variable T.
*/
#define WLZ_SWAP(T,A,B)	{(T)=(A);(A)=(B);(B)=(T);}

/************************************************************************
* Pixel value clamping. Used in WlzGreyScan() and WlzConvertPix().	*
************************************************************************/
/*!
* \def		WLZ_CLAMP(V,N,X)
* \brief	Clamps the value of V such that V >= N and V <= X.
*/
#define WLZ_CLAMP(V,N,X) (((V)<(N))?(N):((V)>(X))?(X):(V))

/*!
* \def		WLZ_CLAMP_DOUBLE_TO_GREYP(U,V,G)
* \brief	Clamps a single double value from V setting U offset by
* 		F such that U >= MIN and U <= MAX where MIN and MAX are
* 		the minimum and maximum values for the grey type G.
* 		The only valid types for G are: WLZ_GREY_UBYTE,
* 		WLZ_GREY_SHORT, WLZ_GREY_INT, WLZ_GREY_FLOAT, and
* 		WLZ_GREY_DOUBLE; there is no check that G is none of
* 		these types.
*/
#define WLZ_CLAMP_DOUBLE_TO_GREYP(U,F,V,G) \
		{\
		  switch((G))\
		  {\
		    case WLZ_GREY_UBYTE:\
		      (U).inp[(F)]=WLZ_CLAMP((V),0,255);\
		      break;\
		    case WLZ_GREY_SHORT:\
		      (U).shp[(F)]=WLZ_CLAMP((V),SHRT_MIN,SHRT_MAX);\
		      break;\
		    case WLZ_GREY_INT:\
		      (U).inp[(F)]=WLZ_CLAMP((V),-FLT_MAX,FLT_MAX);\
		      break;\
		    case WLZ_GREY_FLOAT:\
		      (U).flp[(F)]=WLZ_CLAMP((V),-FLT_MAX,FLT_MAX);\
		      break;\
		    case WLZ_GREY_DOUBLE:\
		      (U).dbp[(F)]=(V);\
		      break;\
		    default:\
		      break;\
		  }\
		}\

/************************************************************************
* Math constants.							*
************************************************************************/
#define	WLZ_M_E		(2.7182818284590452354)0
#define	WLZ_M_LOG2E	(1.4426950408889634074)
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
#define WLZ_RGBA_RED_GET(V) 	((((WlzUInt )V) & 0xff)>>0)
#define WLZ_RGBA_GREEN_GET(V) 	((((WlzUInt )V) & 0xff00)>>8)
#define WLZ_RGBA_BLUE_GET(V) 	((((WlzUInt )V) & 0xff0000)>>16)
#define WLZ_RGBA_ALPHA_GET(V) 	((((WlzUInt )V) & 0xff000000)>>24)

#define WLZ_RGBA_RED_SET(V,C)	(V = (((WlzUInt )V)&0xffffff00) | \
                                               (((WlzUInt) C)&0xff))
#define WLZ_RGBA_GREEN_SET(V,C)	(V = (((WlzUInt )V)&0xffff00ff) | \
                                               ((((WlzUInt) C)&0xff)<<8))
#define WLZ_RGBA_BLUE_SET(V,C)	(V = (((WlzUInt )V)&0xff00ffff) | \
                                               ((((WlzUInt) C)&0xff)<<16))
#define WLZ_RGBA_ALPHA_SET(V,C)	(V = (((WlzUInt )V)&0x00ffffff) | \
                                               ((((WlzUInt) C)&0xff)<<24))
#define WLZ_RGBA_RGBA_SET(V,R,G,B,A) (V = ((((WlzUInt )R)&0xff) + \
                                           ((((WlzUInt )G)&0xff)<<8) + \
                                           ((((WlzUInt )B)&0xff)<<16) + \
                                           ((((WlzUInt )A)&0xff)<<24))) 

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
/* WLZ_VTX_2_COPY: Copy values of Wlz[DFI]Vertex2 */
#define WLZ_VTX_2_COPY(D,S) \
		(D).vtX = (S).vtX, \
		(D).vtY = (S).vtY
/* WLZ_VTX_3_COPY: Copy values of Wlz[DFI]Vertex3 */
#define WLZ_VTX_3_COPY(D,S) \
		(D).vtX = (S).vtX, \
		(D).vtY = (S).vtY, \
		(D).vtZ = (S).vtZ
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

/* WLZ_VTX_4_ADD: Add four Wlz[DFI]Vertex3's */
#define WLZ_VTX_3_ADD4(U,V,W,X,Y) \
		(U).vtX = (V).vtX + (W).vtX + (X).vtX + (Y).vtX, \
		(U).vtY = (V).vtY + (W).vtY + (X).vtY + (Y).vtY, \
		(U).vtZ = (V).vtZ + (W).vtZ + (X).vtZ + (Y).vtZ

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
/* WLZ_VTX_2_SCALE_ADD: Scale a Wlz[DFI]Vertex2 then add a Wlz[DFI]Vertex2. */
#define WLZ_VTX_2_SCALE_ADD(U,V,C,W) \
		(U).vtX = (V).vtX * (C) + (W.vtX), \
		(U).vtY = (V).vtY * (C) + (W.vtY)

/* WLZ_VTX_3_SCALE_ADD: Scale a Wlz[DFI]Vertex3 then add a Wlz[DFI]Vertex3. */
#define WLZ_VTX_3_SCALE_ADD(U,V,C,W) \
		(U).vtX = (V).vtX * (C) + (W.vtX), \
		(U).vtY = (V).vtY * (C) + (W.vtY), \
		(U).vtZ = (V).vtZ * (C) + (W.vtZ)

/* WLZ_VTX_2_DOT: Dot (scalar) product of two Wlz[DFI]Vertex2's */
#define WLZ_VTX_2_DOT(V,W) \
		(((V).vtX * (W).vtX) + \
		 ((V).vtY * (W).vtY))

/* WLZ_VTX_3_DOT: Dot (scalar) product of two Wlz[DFI]Vertex3s */
#define WLZ_VTX_3_DOT(V,W) \
		(((V).vtX * (W).vtX) + \
		 ((V).vtY * (W).vtY) + \
		 ((V).vtZ * (W).vtZ))

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

/* WLZ_VTX_2_ZERO: Set vector to zero. */
#define WLZ_VTX_2_ZERO(U) \
		(U).vtX = 0, \
		(U).vtY = 0

/* WLZ_VTX_3_ZERO: Set vector to zero. */
#define WLZ_VTX_3_ZERO(U) \
		(U).vtX = 0, \
		(U).vtY = 0, \
		(U).vtZ = 0

/* WLZ_VTX_2_NEGATE: Set vector to zero. */
#define WLZ_VTX_2_NEGATE(U,V) \
		(U).vtX = -((V).vtX), \
		(U).vtY = -((V).vtY)

/* WLZ_VTX_3_NEGATE: Set vector to zero. */
#define WLZ_VTX_3_NEGATE(U,V) \
		(U).vtX = -((V).vtX), \
		(U).vtY = -((V).vtY), \
		(U).vtZ = -((V).vtZ)

/* WLZ_VTX_2_EQUAL: Tests for equality of vertices. */
#define WLZ_VTX_2_EQUAL(U,V,T) \
		((fabs((U).vtX - (V).vtX) < (T)) && \
		 (fabs((U).vtY - (V).vtY) < (T)))

/* WLZ_VTX_3_EQUAL: Tests for equality of vertices. */
#define WLZ_VTX_3_EQUAL(U,V,T) \
		((fabs((U).vtX - (V).vtX) < (T)) && \
		 (fabs((U).vtY - (V).vtY) < (T)) && \
		 (fabs((U).vtZ - (V).vtZ) < (T)))

/* WLZ_VTX_2_FABS: Floating point absolute value. */
#define WLZ_VTX_2_FABS(U,V) \
              	((U).vtX = fabs((V).vtX), \
              	 (U).vtY = fabs((V).vtY))

/* WLZ_VTX_3_FABS: Floating point absolute value. */
#define WLZ_VTX_3_FABS(U,V) \
              	((U).vtX = fabs((V).vtX), \
              	 (U).vtY = fabs((V).vtY), \
              	 (U).vtZ = fabs((V).vtZ))

/* WLZ_VTX_2_NINT: Nearest integer position. */
#define WLZ_VTX_2_NINT(U,P) \
        	((U).vtX = WLZ_NINT((P).vtX), \
		 (U).vtY = WLZ_NINT((P).vtY))

/* WLZ_VTX_3_NINT: Nearest integer position. */
#define WLZ_VTX_3_NINT(U,P) \
        	((U).vtX = WLZ_NINT((P).vtX), \
		 (U).vtY = WLZ_NINT((P).vtY), \
		 (U).vtZ = WLZ_NINT((P).vtZ))

/************************************************************************
* CMesh node access
************************************************************************/
#define WLZ_CMESH_ELM2D_GET_NODE_0(e)	((e)->edu[0].nod)
#define WLZ_CMESH_ELM2D_GET_NODE_1(e)	((e)->edu[1].nod)
#define WLZ_CMESH_ELM2D_GET_NODE_2(e)	((e)->edu[2].nod)
#define WLZ_CMESH_ELM2D5_GET_NODE_0(e)	((e)->edu[0].nod)
#define WLZ_CMESH_ELM2D5_GET_NODE_1(e)	((e)->edu[1].nod)
#define WLZ_CMESH_ELM2D5_GET_NODE_2(e)	((e)->edu[2].nod)
#define WLZ_CMESH_ELM3D_GET_NODE_0(e)	((e)->face[0].edu[0].nod)
#define WLZ_CMESH_ELM3D_GET_NODE_1(e)	((e)->face[0].edu[1].nod)
#define WLZ_CMESH_ELM3D_GET_NODE_2(e)	((e)->face[0].edu[2].nod)
#define WLZ_CMESH_ELM3D_GET_NODE_3(e)	((e)->face[1].edu[1].nod)

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
}
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#endif /* WLZMACRO_H */
