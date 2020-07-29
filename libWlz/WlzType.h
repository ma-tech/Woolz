#ifndef WLZ_TYPE_H
#define WLZ_TYPE_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzType_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzType.h
* \author       Bill Hill
* \date         April 2001
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
* \brief	Defines the Woolz types. These are enumerations and
* 		structures which have been typedef'd.
* \ingroup	Wlz
*/

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

/*!
* \typedef	WlzUByte
* \ingroup	WlzType
* \brief	An eight bit unsigned integer.
*/
typedef unsigned char WlzUByte;

/*!
* \typedef	WlzUInt
* \ingroup	WlzType
* \brief        A 32 bit unsigned integer.
*/
typedef unsigned int  WlzUInt;

/*!
* \typedef	WlzLong
* \ingroup	WlzType
* \brief        A 64 bit integer.
*/
#ifdef WLZ_EXT_BIND
typedef long  WlzLong;
#else /* WLZ_EXT_BIND */
typedef long long  WlzLong;
#endif /* WLZ_EXT_BIND */

/*!
* \typedef	WlzULong
* \ingroup	WlzType
* \brief        A 64 bit unsigned integer.
*/
#ifdef WLZ_EXT_BIND
typedef unsigned long  WlzULong;
#else /* WLZ_EXT_BIND */
typedef unsigned long long  WlzULong;
#endif /* WLZ_EXT_BIND */

/*!
* \enum		_WlzGreyType
* \ingroup	WlzType
* \brief	The valid grey value types.
*		Typedef: ::WlzGreyType.
*/
typedef enum _WlzGreyType
{
  WLZ_GREY_LONG			= 0,	/*!< Signed WlzLong integer. */
  WLZ_GREY_INT			= 1,	/*!< Signed integer. */
  WLZ_GREY_SHORT		= 2,	/*!< Signed short. */
  WLZ_GREY_UBYTE		= 3,	/*!< Unsigned WlzUByte integer. */
  WLZ_GREY_FLOAT		= 4,	/*!< Single precision floating
  					     point. */
  WLZ_GREY_DOUBLE		= 5,	/*!< Double precision floating
  					     point. */
  WLZ_GREY_BIT			= 6,	/*!< Single bit. */
  WLZ_GREY_RGBA			= 7,	/*!< Eight bit red, green, blue and
  					     alpha components packed into an
					     unsigned 32 bit integer, with:
					     v = r|(g<<8)|(b<<16)|(a<<24),
					     where r, g, b and a are the
					     unsigned byte red, green, blue
					     and alpha components of the
					     value. */
  WLZ_GREY_ERROR			/*!< An invalid grey type used to
  					     return error conditions.
					     Always the last enumerator! */
} WlzGreyType;

/*!
* \def		WLZ_GREY_TABLE_TYPE
* \ingroup	WlzType
* \brief	For historical reasons a pixel/voxel value table
*               encodes both the grey type and the table type in
*               a single type. This macro achieves this with the
*               addition of array type.
* \par
* AR:		Array type, 0 for scalar values, 1 for any array
* 		type.
* TT:		Table type.
* GT:		Grey type.
*/
#define WLZ_GREY_TABLE_TYPE(AR,TT,GT)    ((100*(!!(AR)))+(10*(TT))+(GT))

/*!
* \def		WLZ_GREY_TABLE_TO_GREY_TYPE
* \ingroup	WlzType
* \brief	Grey type encoded in grey table type. See
* 		WLZ_GREY_TABLE_TYPE.
* \par
* GTT:		Grey table type.
*/
#define WLZ_GREY_TABLE_TO_GREY_TYPE(GTT)  ((GTT)%10)

/*!
* \def		WLZ_GREY_TABLE_TO_TABLE_TYPE
* \ingroup	WlzType
* \brief	Table type encoded in grey table type. See
* 		WLZ_GREY_TABLE_TYPE.
* \par
* GTT:		Grey table type.
*/
#define WLZ_GREY_TABLE_TO_TABLE_TYPE(GTT) (((GTT)%100)/10)

/*!
* \def		WLZ_GREY_TABLE_TO_RANK
* \ingroup	WlzType
* \brief	Rank encoded in grey table type. See
* 		WLZ_GREY_TABLE_TYPE.
* \par
* GTT:		Grey table type.
*/
#define WLZ_GREY_TABLE_TO_RANK(GTT)       ((GTT)/100)

/*!
* \enum		_WlzGreyTableType
* \ingroup	WlzType
* \brief	Woolz pixel/voxel value table types.
*/
typedef enum	_WlzGreyTableType
{
  WLZ_GREY_TAB_RAGR		= 0,	/*!< Base ragged rectangle grey value
  					     table. */
  WLZ_GREY_TAB_RECT		= 1,	/*!< Base rectangular grey value
  					     table. */
  WLZ_GREY_TAB_INTL		= 2,	/*!< Base interval grey value table. */
  WLZ_POINT_VALUES              = 3,    /*!< Values defined at points. */
  WLZ_INDEXED_VALUES		= 4,	/*!< Indexed value table. */
  WLZ_FEAT_TAB_RAGR		= 5,	/*!< Base ragged rectangle feature
  					     table. */
  WLZ_FEAT_TAB_RECT 		= 6,	/*!< Base rectangular feature table. */
  WLZ_GREY_TAB_TILED		= 7     /*!< Tiled grey value table. */
} WlzGreyTableType;

/*!
* \enum		_WlzObjectType
* \ingroup	WlzType
* \brief	The Woolz object types.
*		Many of the integer enumeration values are required for
*		historical reasons but, with the execption of WLZ_NULL,
*		the numerical value should never be used explicitly.
*
*		For historical reasons the base grey table types are used
*		to form explicit grey table types which include the grey type.
*		The functions or macros for extracting and synthesising these
*		types should be used. These are:
*               WLZ_GREY_TABLE_TYPE,
*               WLZ_GREY_TABLE_TO_GREY_TYPE,
*               WLZ_GREY_TABLE_TO_TABLE_TYPE,
*               WLZ_GREY_TABLE_TO_RANK,
*		WlzGreyTableType(),
*		WlzGreyTableTypeToGreyType(),
*		WlzGreyTableTypeToRank(),
*		WlzGreyTableTypeToTableType().
*		Typedef: ::WlzObjectType.
*/
typedef enum _WlzObjectType
{
  /**********************************************************************
  * Top level Woolz object types.
  **********************************************************************/
  WLZ_NULL			= 0,	/*!< The NULL object is an invalid
  					     object used to return an error
					     condition. */
  WLZ_2D_DOMAINOBJ		= 1,	/*!< 2D spatial domain object. */
  WLZ_3D_DOMAINOBJ		= 2,	/*!< 3D spatial domain object. */
  WLZ_TRANS_OBJ			= 3,	/*!< Object in which the domain is
  					     a transformation and the values
					     are an object. */
  WLZ_3D_WARP_TRANS		= 4,	/*!< 3D finite element warp
  					     transformation. */
  WLZ_2D_POLYGON		= 10,	/*!< 2D polygon. */
  WLZ_BOUNDLIST			= 11,	/*!< Boundary tree. */
  WLZ_CONV_HULL			= 12,	/*!< Convex hull. */
  WLZ_HISTOGRAM			= 13,	/*!< Histogram. */
  WLZ_3D_POLYGON		= 14,	/*!< 3D polygon. */
  WLZ_CONTOUR			= 15,	/*!< Contour in either 2D or 3D. */
  WLZ_CMESH_2D			= 16,   /*!< Constrained mesh with triangular
                                             elements in 2D. */
  WLZ_CMESH_3D			= 17,   /*!< Constrained mesh with tetrahedral
  					     elements in 3D. */
  WLZ_CMESH_2D5			= 18,   /*!< Constrained mesh with triangular
  					     elements in 3D. */
  WLZ_RECTANGLE			= 20,	/*!< Rectangle. */
  WLZ_POINTS			= 21,   /*!< Points. */
  WLZ_SPLINE			= 30,	/*!< Splines. */
  WLZ_CONVOLVE_INT		= 50,	/*!< Integer convolution. */
  WLZ_CONVOLVE_FLOAT		= 51,	/*!< Floating point convolution. */
  WLZ_AFFINE_TRANS		= 63,	/*!< Affine transform, either 2D or
  					    3D. */
  WLZ_WARP_TRANS		= 64,	/*!< A finite element warp
  					     transformation. */
  WLZ_FMATCHOBJ			= 65,	/*!< Matched features for finite
  					     element warps. */
  WLZ_TEXT			= 70,	/*!< Simple ascii text. */
  WLZ_COMPOUND_ARR_1		= 80,	/*!< Compound array of objects
  					     with same type. */
  WLZ_COMPOUND_ARR_2		= 81,	/*!< Compound array of objects
  					     with differnt type. */
  WLZ_COMPOUND_LIST_1		= 82,	/*!< Linked list of objects
  					     with same type. */
  WLZ_COMPOUND_LIST_2		= 83,	/*!< Linked list of objects
  					     with different type. */
  WLZ_LUT			= 90,	/*!< Simple look up table. */
  WLZ_PROPERTY_OBJ		= 110,	/*!< An object which only has a
  					     property list. */
  WLZ_EMPTY_OBJ			= 127,	/*!< Empty object: An object which
  					     occupies no space and has no
					     values. */
  WLZ_MESH_TRANS                = 128,  /*!< Mesh transform. */
  WLZ_CMESH_TRANS               = 129,  /*!< Constrained mesh transform,
  					     either 2D or 3D. */

  WLZ_EMPTY_DOMAIN,			/*!< Empty domain: A domain which
  					     occupies no space. */
  WLZ_EMPTY_VALUES,			/*!< Empty values: A values which
  					     has no values! */
  /**********************************************************************
  * Plane domain types.
  **********************************************************************/
  WLZ_INTERVALDOMAIN_INTVL	= 1,	/*!< Spatial domain defined by
  					     a collection of intervals. */
  WLZ_INTERVALDOMAIN_RECT	= 2,	/*!< Spatial domain defined by an
  					     axis aligned rectangle. */
  WLZ_LBTDOMAIN_2D		= 3,	/*!< 2D linear binary tree domain. */
  WLZ_LBTDOMAIN_3D		= 4,	/*!< 3D linear binary tree domain. */
  WLZ_PLANEDOMAIN_DOMAIN	= WLZ_2D_DOMAINOBJ, /*!< 3D spatial domain
  					     composed of 2D spatial domains. */
  WLZ_PLANEDOMAIN_POLYGON	= WLZ_2D_POLYGON, /*!< 3D polygon domain
  					     composed of 2D polygon domains. */
  WLZ_PLANEDOMAIN_BOUNDLIST	= WLZ_BOUNDLIST, /*!< 3D boundary domain
  					     composed of 2D boundary
					     domains. */
  WLZ_PLANEDOMAIN_CONV_HULL	= WLZ_CONV_HULL, /*!< 3D convex hull
  					     composed of 2D convex hulls. */
  WLZ_PLANEDOMAIN_HISTOGRAM	= WLZ_HISTOGRAM, /*!< 3D histogram domain
  					     composed of 2D histogram
					     domains. */
  WLZ_PLANEDOMAIN_AFFINE	= WLZ_AFFINE_TRANS, /*!< 3D affine domain
  					     composed of 2D affine transform
					     domains. */
  WLZ_PLANEDOMAIN_WARP		= WLZ_WARP_TRANS, /*!< 3D warp domain
  					     composed of 2D warp domains. */
  /**********************************************************************
  * Value table types.
  **********************************************************************/
  WLZ_VALUETABLE_RAGR_INT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RAGR, WLZ_GREY_INT),
  					/*!< Ragged rectangle int value
					     table. */
  WLZ_VALUETABLE_RAGR_SHORT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RAGR, WLZ_GREY_SHORT),
  					/*!< Ragged rectangle short value
					     table. */
  WLZ_VALUETABLE_RAGR_UBYTE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RAGR, WLZ_GREY_UBYTE),
  					/*!< Ragged rectangle unsigned byte
					     value table. */
  WLZ_VALUETABLE_RAGR_FLOAT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RAGR, WLZ_GREY_FLOAT),
  					/*!< Ragged rectangle single precision
					     floating point value table. */
  WLZ_VALUETABLE_RAGR_DOUBLE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE),
  					/*!< Ragged rectangle double precision
					     floating point value table. */
  WLZ_VALUETABLE_RAGR_BIT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RAGR, WLZ_GREY_BIT),
  					/*!< Ragged rectangle single bit (packed
					     in unsigned bytes) value table. */
  WLZ_VALUETABLE_RAGR_RGBA	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RAGR, WLZ_GREY_RGBA),
  					/*!< Ragged rectangle red, green, blue,
					     alpha value table. */
  WLZ_VALUETABLE_RECT_INT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RECT, WLZ_GREY_INT),
  					/*!< Rectangular int value table. */
  WLZ_VALUETABLE_RECT_SHORT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RECT, WLZ_GREY_SHORT),
  					/*!< Rectangular short value table. */
  WLZ_VALUETABLE_RECT_UBYTE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RECT, WLZ_GREY_UBYTE),
  					/*!< Rectangular unsigned byte value
					     table. */
  WLZ_VALUETABLE_RECT_FLOAT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RECT, WLZ_GREY_FLOAT),
  					/*!< Rectangular single precision
					     floating point value table. */
  WLZ_VALUETABLE_RECT_DOUBLE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RECT, WLZ_GREY_DOUBLE),
  					/*!< Rectangular double precision
					     floating point value table. */
  WLZ_VALUETABLE_RECT_BIT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RECT, WLZ_GREY_BIT),
  					/*!< Rectangular single bit (packed in
					     unsigned bytes) value table. */
  WLZ_VALUETABLE_RECT_RGBA	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_RECT, WLZ_GREY_RGBA),
  					/*!< Rectangular red, green, blue,
					     alpha value table.  */ 
  WLZ_VALUETABLE_INTL_INT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_INTL, WLZ_GREY_INT),
  					/*!< Interval coded int value table. */
  WLZ_VALUETABLE_INTL_SHORT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_INTL, WLZ_GREY_SHORT),
  					/*!< Interval coded short value
					     table. */
  WLZ_VALUETABLE_INTL_UBYTE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_INTL, WLZ_GREY_UBYTE),
  					/*!< Interval coded unsigned byte value
					     table. */
  WLZ_VALUETABLE_INTL_FLOAT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_INTL, WLZ_GREY_FLOAT),
  					/*!< Interval coded single precision
					     floating point value table. */
  WLZ_VALUETABLE_INTL_DOUBLE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_INTL, WLZ_GREY_DOUBLE),
  					/*!< Interval coded double precision
					     floating point value table. */
  WLZ_VALUETABLE_INTL_BIT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_INTL, WLZ_GREY_BIT),
  					/*!< Interval coded single bit (packed
					     unsigned bytes) value table. */
  WLZ_VALUETABLE_INTL_RGBA	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_INTL, WLZ_GREY_RGBA),
  					/*!< Interval coded red, green, blue,
					     alpha value table.  */
  WLZ_VALUETABLE_TILED_INT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_TILED, WLZ_GREY_INT),
  					/*!< Tiled int value table. */
  WLZ_VALUETABLE_TILED_SHORT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_TILED, WLZ_GREY_SHORT),
  					/*!< Tiled short value table. */
  WLZ_VALUETABLE_TILED_UBYTE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_TILED, WLZ_GREY_UBYTE),
  					/*!< Tiled unsigned byte value
					     table. */
  WLZ_VALUETABLE_TILED_FLOAT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_TILED, WLZ_GREY_FLOAT),
  					/*!< Tiled single precision
					     floating point value table. */
  WLZ_VALUETABLE_TILED_DOUBLE	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_TILED, WLZ_GREY_DOUBLE),
  					/*!< Tiled double precision
					     floating point value table. */
  WLZ_VALUETABLE_TILED_BIT	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_TILED, WLZ_GREY_BIT),
  					/*!< Tiled single bit (packed
					     unsigned bytes) value table. */
  WLZ_VALUETABLE_TILED_RGBA	= WLZ_GREY_TABLE_TYPE(
                                    0, WLZ_GREY_TAB_TILED, WLZ_GREY_RGBA),
  					/*!< Tiled red, green, blue, table. */
  WLZ_VALUETABLE_TILED_ARY_INT	= WLZ_GREY_TABLE_TYPE(
                                    1, WLZ_GREY_TAB_TILED, WLZ_GREY_INT),
  					/*!< Tiled int value table. */
  WLZ_VALUETABLE_TILED_ARY_SHORT = WLZ_GREY_TABLE_TYPE(
                                    1, WLZ_GREY_TAB_TILED, WLZ_GREY_SHORT),
  					/*!< Tiled short value table. */
  WLZ_VALUETABLE_TILED_ARY_UBYTE = WLZ_GREY_TABLE_TYPE(
                                    1, WLZ_GREY_TAB_TILED, WLZ_GREY_UBYTE),
  					/*!< Tiled unsigned byte value
					     table. */
  WLZ_VALUETABLE_TILED_ARY_FLOAT = WLZ_GREY_TABLE_TYPE(
                                    1, WLZ_GREY_TAB_TILED, WLZ_GREY_FLOAT),
  					/*!< Tiled single precision
					     floating point value table. */
  WLZ_VALUETABLE_TILED_ARY_DOUBLE = WLZ_GREY_TABLE_TYPE(
                                    1, WLZ_GREY_TAB_TILED, WLZ_GREY_DOUBLE),
  					/*!< Tiled double precision
					     floating point value table. */
  WLZ_VALUETABLE_TILED_ARY_BIT	= WLZ_GREY_TABLE_TYPE(
                                    1, WLZ_GREY_TAB_TILED, WLZ_GREY_BIT),
  					/*!< Tiled single bit (packed
					     unsigned bytes) value table. */
  WLZ_VALUETABLE_TILED_ARY_RGBA	= WLZ_GREY_TABLE_TYPE(
                                    1, WLZ_GREY_TAB_TILED, WLZ_GREY_RGBA),
  					/*!< Tiled red, green, blue, table. */
  WLZ_FEATVALUETABLE_RAGR	= 50,	/*!< Ragged rectangle features
  					     value table. */
  WLZ_FEATVALUETABLE_RECT	= 60,	/*!< Rectangular features
  					     value table. */
  WLZ_VOXELVALUETABLE_GREY	= 1,	/*!< Grey value voxel value table. */
  WLZ_VOXELVALUETABLE_CONV_HULL,	/*!< Convex hull voxel value table. */
  /**********************************************************************
  * Polygon domain types.					
  **********************************************************************/
  WLZ_POLYGON_INT		= 1,	/*!< Integer polygon domain. */
  WLZ_POLYGON_FLOAT		= 2,	/*!< Single precision floating point
  					     polygon domain. */
  WLZ_POLYGON_DOUBLE		= 3,	/*!< Double precision floating point
  					     polygon domain. */
  /**********************************************************************
  * Boundary list types.					
  **********************************************************************/
  WLZ_BOUNDLIST_PIECE		= 0,	/*!< Piece contains foreground. */
  WLZ_BOUNDLIST_HOLE		= 1,	/*!< Piece contains background. */
  /**********************************************************************
  * Convex hull types.						
  **********************************************************************/
  WLZ_CONVHULL_VALUES		= 1,	/*!< Convex hull values. */
  WLZ_CONVHULL_DOMAIN_2D	= 2,	/*!< 2D convex hull domain. */
  WLZ_CONVHULL_DOMAIN_3D	= 3,	/*!< 3D convex hull domain. */
  /**********************************************************************
  * Histogram domain types. WLZ_HISTOGRAMDOMAIN_OLD_INT and	
  * WLZ_HISTOGRAMDOMAIN_OLD_FLOAT exist only to allow old files to be
  * read, they should not be used anywhere else.		
  **********************************************************************/
  WLZ_HISTOGRAMDOMAIN_OLD_INT	= 1,	/*!< Historical compatability. */
  WLZ_HISTOGRAMDOMAIN_OLD_FLOAT	= 2,	/*!< Historical compatability. */
  WLZ_HISTOGRAMDOMAIN_INT	= 3,	/*!< Integer histogram domain. */
  WLZ_HISTOGRAMDOMAIN_FLOAT	= 4,	/*!< Floating point histogram
  					     domain. */
  /**********************************************************************
  * Rectangle object domain types.				
  **********************************************************************/
  WLZ_RECTANGLE_DOMAIN_INT	= 1,	/*!< Integer rectangle domain. */
  WLZ_RECTANGLE_DOMAIN_FLOAT	= 2,	/*!< Floating point rectangle
  					     domain. */
  /**********************************************************************
  * 3D view structure types (also used for object and domain type).
  **********************************************************************/
  WLZ_3D_VIEW_STRUCT		= 160,	/*!< 3D view structure. */
  /**********************************************************************
  * Property list types.					
  **********************************************************************/
  WLZ_PROPERTYLIST		= 170,	/*! A property list which contains
  					    a linked list of properties. */
  /**********************************************************************
  * Property types.					
  **********************************************************************/
  WLZ_PROPERTY_SIMPLE		= 180,	/*!< Simple property. */
  WLZ_PROPERTY_EMAP		= 181,	/*!< EMAP property. */
  WLZ_PROPERTY_NAME		= 182,	/*!< Ascii name property. */
  WLZ_PROPERTY_GREY		= 183,	/*!< Grey value property. */
  WLZ_PROPERTY_TEXT		= 184,	/*!< Text property. */
  /**********************************************************************
  * Points domain types.
  **********************************************************************/
  WLZ_POINTS_2I			= 21,	/*!< Integer 2D points domain. */
  WLZ_POINTS_2D			= 22,	/*!< Integer 2D points domain. */
  WLZ_POINTS_3I	                = 23,   /*!< Double precision floating point
					      3D points domain. */
  WLZ_POINTS_3D	                = 24,   /*!< Double precision floating point
					      3D points domain. */
  /**********************************************************************
  * Spline domain types.
  **********************************************************************/
  WLZ_BSPLINE_C2D              = 30,  /*!< 2D B-spline line curve domains. */
  WLZ_BSPLINE_C3D              = 31,  /*!< 3D B-spline line curve domains. */
  /**********************************************************************
  * WLZ_DUMMY_ENTRY is not an object type.			
  * Keep it the last enumerator!				
  **********************************************************************/
  WLZ_DUMMY_ENTRY			/*!< Not a Woolz object type.
  					     Keep it the last enumerator! */
} WlzObjectType;


/*! 
* \enum		_WlzEMAPPropertyType
* \ingroup      WlzProperty
* \brief        Sub types of EMAP properties
*		Typedef: ::WlzEMAPPropertyType.
*/
typedef enum _WlzEMAPPropertyType
{
  WLZ_EMAP_PROPERTY_GREY_MODEL = 1,	/*!< Reference Model */
  WLZ_EMAP_PROPERTY_GREY_OTHER,		/*!< Grey object e.g. derived from model */
  WLZ_EMAP_PROPERTY_DOMAIN_ANATOMY,	/*!< Anatomy Domain */
  WLZ_EMAP_PROPERTY_DOMAIN_OTHER,	/*!< Other Domain */
  WLZ_EMAP_PROPERTY_TRANSFORM,		/*!< woolz transform for EMAP models */
  WLZ_EMAP_PROPERTY_DUMMY		/*!< Dummy property */
} WlzEMAPPropertyType;

/*!
* \enum		_WlzRasterDir
* \ingroup	WlzAccess
* \brief	Raster scan directions as used by WlzIntervalWSpace
*               and WlzIterateWSpace. These are built using bit masks
*               with bits set for decreasing in each of the directions.
*		Typedef: ::WlzRasterDir.
*/
typedef enum _WlzRasterDir
{
  WLZ_RASTERDIR_DC    = (1),		/*!< Used to build directions . */
  WLZ_RASTERDIR_DL    = (1 << 1),	/*!< Used to build directions. */
  WLZ_RASTERDIR_DP    = (1 << 2),	/*!< Used to build directions. */
  WLZ_RASTERDIR_ILIC  = (0),		/*!< Increasing lines, increasing
  					     columns. */
  WLZ_RASTERDIR_ILDC  = WLZ_RASTERDIR_DC,
  					/*!< Increasing lines, decreasing
  					     columns. */
  WLZ_RASTERDIR_DLIC  = WLZ_RASTERDIR_DL,
  					/*!< Decreasing lines, increasing
  					     columns. */
  WLZ_RASTERDIR_DLDC  = (WLZ_RASTERDIR_DL | WLZ_RASTERDIR_DC),
                                      	/*!< Decreasing lines, decreasing
  					     columns. */
  WLZ_RASTERDIR_IPILIC  = WLZ_RASTERDIR_ILIC,
  					/*!< Increasing planes,
                                             increasing lines, increasing
  					     columns. */
  WLZ_RASTERDIR_IPILDC  = WLZ_RASTERDIR_ILDC,
  					/*!< Increasing planes,
                                             increasing lines, decreasing
  					     columns. */
  WLZ_RASTERDIR_IPDLIC  = WLZ_RASTERDIR_DLIC,
  					/*!< Increasing planes,
                                             decreasing lines, increasing
  					     columns. */
  WLZ_RASTERDIR_IPDLDC  = WLZ_RASTERDIR_DLDC,
  					/*!< Increasing planes,
                                             decreasing lines, decreasing
  					     columns. */
  WLZ_RASTERDIR_DPILIC  = (WLZ_RASTERDIR_DP | WLZ_RASTERDIR_ILIC),
  					/*!< Decreasing planes,
                                             increasing lines, increasing
  					     columns. */
  WLZ_RASTERDIR_DPILDC  = (WLZ_RASTERDIR_DP | WLZ_RASTERDIR_ILDC),
  					/*!< Decreasing planes,
                                             increasing lines, decreasing
  					     columns. */
  WLZ_RASTERDIR_DPDLIC  = (WLZ_RASTERDIR_DP | WLZ_RASTERDIR_DLIC),
  					/*!< Decreasing planes,
                                             decreasing lines, increasing
  					     columns. */
  WLZ_RASTERDIR_DPDLDC  = (WLZ_RASTERDIR_DP | WLZ_RASTERDIR_DLDC)
  					/*!< Decreasing planes,
                                             decreasing lines, decreasing
  					     columns. */
} WlzRasterDir;

/*!
* \enum         _WlzDirection
* \ingroup      WlzType
* \brief        Basic directions.
*/
typedef enum _WlzDirection
{
  WLZ_DIRECTION_IC,             /*!< Increasing columns, ++x, east. */
  WLZ_DIRECTION_IL,             /*!< Increasing lines, ++y, south. */
  WLZ_DIRECTION_IP,             /*!< Increasing planes, ++z, up. */
  WLZ_DIRECTION_DC,             /*!< Decreasing columns, --x, west. */
  WLZ_DIRECTION_DL,             /*!< Decreasing lines, --y, north. */
  WLZ_DIRECTION_DP              /*!< Decreasing planes, --z, down. */
} WlzDirection;

/*!
* \enum		_WlzTransformType
* \ingroup	WlzTransform
* \brief	Types of spatial transformation.
*		Typedef: ::WlzTransformType.
*/
typedef enum _WlzTransformType
{
  WLZ_TRANSFORM_EMPTY = 0,		/*!< Undefined transform. */
  WLZ_TRANSFORM_2D_AFFINE = 1,		/*!< General 2D affine transform. */
  WLZ_TRANSFORM_2D_REG,	      		/*!< 2D affine but only rotation
  					     and translation. */
  WLZ_TRANSFORM_2D_TRANS,       	/*!< 2D affine but only translation. */
  WLZ_TRANSFORM_2D_NOSHEAR,             /*!< 2D affine but no shear. */
  WLZ_TRANSFORM_3D_AFFINE,		/*!< General 3D affine transform. */
  WLZ_TRANSFORM_3D_REG,	      		/*!< 3D affine but only rotation
  					     and translation. */
  WLZ_TRANSFORM_3D_TRANS,       	/*!< 3D affine but only translation. */
  WLZ_TRANSFORM_3D_NOSHEAR,             /*!< 3D affine but no shear. */
  WLZ_TRANSFORM_2D_BASISFN,		/*!< 2D basis function transform. */
  WLZ_TRANSFORM_2D5_BASISFN,   		/*!< 2.5D (plane wise) basis function
  					     transform. */
  WLZ_TRANSFORM_3D_BASISFN,             /*!< 3D basis function transform. */
  WLZ_TRANSFORM_2D_MESH,                /*!< 2D triangular mesh transform. */
  WLZ_TRANSFORM_2D5_MESH,     		/*!< 2.5D (plane wise) triangular
  					     mesh transform. */
  WLZ_TRANSFORM_3D_MESH,		/*!< 3D tetrahedral mesh transform. */
  WLZ_TRANSFORM_2D_CMESH = WLZ_CMESH_2D, /*!< 2D conforming triangular mesh
  				              transform. */
  WLZ_TRANSFORM_2D5_CMESH = WLZ_CMESH_2D5, /*!< 3D conforming triangular mesh
  				                transform. */
  WLZ_TRANSFORM_3D_CMESH = WLZ_CMESH_3D	/*!< 3D conforming tetrahedral mesh
  					     transform. */
} WlzTransformType;

/*!
* \enum		_WlzMeshElemType
* \ingroup	WlzTransform
* \brief	Mesh transform element types.
*		Typedef: ::WlzMeshElemType.
*/
typedef enum _WlzMeshElemType
{
  WLZ_MESH_ELEM_TRILINEAR,
  WLZ_MESH_ELEM_TRIINCOMPRESSIBLE,
  WLZ_MESH_ELEM_TRICOMPRESSIBLE
} WlzMeshElemType;

/*!
* \enum		_WlzMeshElemFlags
* \ingroup	WlzTransform
* \brief	Mesh transform element flag bit masks.
*		Typedef: ::WlzMeshElemFlags.
*/
typedef enum _WlzMeshElemFlags
{
  WLZ_MESH_ELEM_FLAGS_NONE      = (0),
  WLZ_MESH_ELEM_FLAGS_NBR_0     = (1),	  /*!< Neighbour on side 0 exists */
  WLZ_MESH_ELEM_FLAGS_NBR_1     = (1<<1), /*!< Neighbour on side 1 exists */
  WLZ_MESH_ELEM_FLAGS_NBR_2     = (1<<2), /*!< Neighbour on side 2 exists */
  WLZ_MESH_ELEM_FLAGS_ZOMBIE	= (1<<3), /*!< Dead, awaiting replacement */
  WLZ_MESH_ELEM_FLAGS_REFINE	= (1<<4), /*!< Available for refinement */
  WLZ_MESH_ELEM_FLAGS_OUTSIDE	= (1<<5)  /*!< Outside object's domain */
} WlzMeshElemFlags;

/*!
* \enum		_WlzMeshNodeFlags
* \ingroup	WlzTransform
* \brief	Mesh transform node flag masks.
*		Typedef: ::WlzMeshNodeFlags.
*/
typedef enum _WlzMeshNodeFlags
{
  WLZ_MESH_NODE_FLAGS_NONE      = (0),
  WLZ_MESH_NODE_FLAGS_BBOX     	= (1),	  /*!< Created from bounding box. */
  WLZ_MESH_NODE_FLAGS_BLOCK	= (1<<1), /*!< Created by block fill. */
  WLZ_MESH_NODE_FLAGS_IDOM	= (1<<2), /*!< Created to fill interval
  					       domain */
  WLZ_MESH_NODE_FLAGS_POLY	= (1<<3), /*!< Created along polygon domain */
  WLZ_MESH_NODE_FLAGS_ZOMBIE	= (1<<4)  /*!< Dead, awaiting replacement */
} WlzMeshNodeFlags;

/*!
* \enum		_WlzMeshGenMethod
* \ingroup	WlzTransform
* \brief	Mesh generation methods.
*		Typedef: ::WlzMeshGenMethod.
*/
typedef enum _WlzMeshGenMethod
{
  WLZ_MESH_GENMETHOD_BLOCK,	       	/*!< Uniform (triangulated) block
  					     grid. */
  WLZ_MESH_GENMETHOD_GRADIENT, 		/*!< Triangulated grid based on image
  					     gradient. */
  WLZ_MESH_GENMETHOD_CONFORM		/*!< Mesh conforming to domain. */
} WlzMeshGenMethod;

/*!
* \enum		_WlzMeshError
* \ingroup	WlzTransform
* \brief	Mesh error bit masks.
*		Typedef: ::WlzMeshError.
*/
typedef enum _WlzMeshError
{
  WLZ_MESH_ERR_NONE		= (0),	  /*!< No error, mesh valid */
  WLZ_MESH_ERR_ELEM_CW		= (1),	  /*!< Element not CCW */
  WLZ_MESH_ERR_ELEM_INDEX	= (1<<1), /*!< Element index invalid */
  WLZ_MESH_ERR_ELEM_NODE	= (1<<2), /*!< Element node invalid */
  WLZ_MESH_ERR_ELEM_ZOMBIE	= (1<<3), /*!< Element is a zombie */
  WLZ_MESH_ERR_DELEM_CW		= (1<<4), /*!< Displaced element not CCW */
  WLZ_MESH_ERR_NELEM_INDEX	= (1<<5), /*!< Neighbour index invalid */
  WLZ_MESH_ERR_NELEM_NODE	= (1<<6), /*!< Neighbour node invalid */
  WLZ_MESH_ERR_NELEM_NOTNBR	= (1<<7), /*!< Neighbour not a neighbour */
  WLZ_MESH_ERR_NELEM_ZOMBIE	= (1<<8)  /*!< Neighbour is a zombie */
} WlzMeshError;

/*!
* \enum		_WlzConnectType
* \ingroup	WlzType
* \brief	Connectivity in a 2D or 3D digital space.
*		Typedef: ::WlzConnectType.
*/
typedef enum _WlzConnectType
{
  WLZ_0_CONNECTED		= 0,
  WLZ_8_CONNECTED,
  WLZ_4_CONNECTED,
  WLZ_6_CONNECTED,
  WLZ_18_CONNECTED,
  WLZ_26_CONNECTED
} WlzConnectType;

/*!
* \enum		_WlzDistanceType
* \ingroup	WlzType
* \brief	Distance metrics in a 2D or 3D digital space.
*		Typedef: ::WlzDistanceType.
*/
typedef enum _WlzDistanceType
{
  WLZ_8_DISTANCE		= WLZ_8_CONNECTED,
  WLZ_4_DISTANCE		= WLZ_4_CONNECTED,
  WLZ_6_DISTANCE		= WLZ_6_CONNECTED,
  WLZ_18_DISTANCE		= WLZ_18_CONNECTED,
  WLZ_26_DISTANCE		= WLZ_26_CONNECTED,
  WLZ_OCTAGONAL_DISTANCE,
  WLZ_EUCLIDEAN_DISTANCE,
  WLZ_APX_EUCLIDEAN_DISTANCE		/*! Approximate Euclidean. */
} WlzDistanceType;

/*!
* \enum 	_WlzRCCClassIdx
* \ingroup	WlzType
* \brief	Discrete Region Connected Calculus (RCC) clasification
* 		indices. The classifications indices are for bits set
* 		in a classification bitmask.
* 		See also WlzRCCClass() and WlzRegConCalcRCC().
*		Typedef: ::WlzRCCClassIdx.
*/
typedef enum _WlzRCCClassIdx
{
  WLZ_RCCIDX_DC		= (0),		/*!< Edge Connected. */
  WLZ_RCCIDX_EC		= (1),		/*!< Edge Connected. */
  WLZ_RCCIDX_EQ		= (2),		/*!< Equal. */
  WLZ_RCCIDX_PO		= (3),		/*!< Partial Overlap. */
  WLZ_RCCIDX_TPP	= (4),		/*!< Tangential Proper Part. */
  WLZ_RCCIDX_NTPP	= (5),		/*!< Non-Tangential Proper Part. */
  WLZ_RCCIDX_TPPI	= (6),		/*!< Tangential Proper Part Inverse. */
  WLZ_RCCIDX_NTPPI	= (7),		/*!< Non-Tangential Proper Part
  				     	     Inverse. */
  WLZ_RCCIDX_TSUR	= (8),		/*!< Tangential SURrounds. */
  WLZ_RCCIDX_TSURI	= (9),		/*!< Tangential SURrounds Inverse. */
  WLZ_RCCIDX_NTSUR	= (10),		/*!< Non-Tangential SURrounds. */
  WLZ_RCCIDX_NTSURI	= (11),		/*!< Non-Tangential SURrounds
                                             Inverse. */
  WLZ_RCCIDX_ENC	= (12),		/*!< ENCloses. */
  WLZ_RCCIDX_ENCI	= (13),		/*!< ENCloses Inverse. */
  WLZ_RCCIDX_OST	= (14),		/*!< OffSeT. */
  WLZ_RCCIDX_CNT	= (15)		/*!< Not a classification index, but
                                             the number of classification
					     indices, keep last. */
} WlzRCCClassIdx;

/*!
* \enum 	_WlzRCCClass
* \ingroup	WlzType
* \brief	A Discrete Region Connected Calculus (RCC) clasification
* 		of an ordered pair of spatial regions. The classifications
* 		are bit masks to allow for multiple classifications.
* 		See also WlzRegConCalcRCC().
*		Typedef: ::WlzRCCClass.
*/
typedef enum _WlzRCCClass
{
  WLZ_RCC_EMPTY	= 0, 			/*!< Empty:
  					     One or both of the domains is
					     empty. */
  WLZ_RCC_DC	= (1<<WLZ_RCCIDX_DC),	/*!< Disconnected. */
  WLZ_RCC_EC	= (1<<WLZ_RCCIDX_EC),	/*!< Edge Connected. */
  WLZ_RCC_EQ	= (1<<WLZ_RCCIDX_EQ),	/*!< Equal. */
  WLZ_RCC_PO	= (1<<WLZ_RCCIDX_PO),	/*!< Partial Overlap. */
  WLZ_RCC_TPP	= (1<<WLZ_RCCIDX_TPP),	/*!< Tangential Proper Part. */
  WLZ_RCC_NTPP	= (1<<WLZ_RCCIDX_NTPP),	/*!< Non-Tangential Proper Part. */
  WLZ_RCC_TPPI	= (1<<WLZ_RCCIDX_TPPI),	/*!< Tangential Proper Part inverse. */
  WLZ_RCC_NTPPI	= (1<<WLZ_RCCIDX_NTPPI), /*!< Non-Tangential Proper Part
  				     	     inverse. */
  WLZ_RCC_TSUR	= (1<<WLZ_RCCIDX_TSUR),	/*!< Tangential Surrounds:
  					     The first domain is surrounded by
					     and edge connected to the
					     second domain. */
  WLZ_RCC_TSURI	= (1<<WLZ_RCCIDX_TSURI), /*!< Tangential Surrounds inverse:
  					     The second domain is surrounded by
					     and edge connected to the first
					     domain. */
  WLZ_RCC_NTSUR	= (1<<WLZ_RCCIDX_NTSUR), /*!< Non-Tangential Surrounds:
  					     The first domain is surrounded by
					     but not edge connected to the
					     second domain. */
  WLZ_RCC_NTSURI = (1<<WLZ_RCCIDX_NTSURI), /*!< Non-Tangential Surrounds
        				     inverse:
  					     The second domain is surrounded by
					     but not edge connected to the
					     first domain. */
  WLZ_RCC_ENC	= (1<<WLZ_RCCIDX_ENC),	/*!< Encloses:
  					     The majority of the first domain
					     is within the convex hull of the
					     second domain. */
  WLZ_RCC_ENCI	= (1<<WLZ_RCCIDX_ENCI),	/*!< Encloses inverse:
  					     The majority of the second domain
					     is within the convex hull of the
					     first domain. */
  WLZ_RCC_OST	= (1<<WLZ_RCCIDX_OST),	/*!< Offset:
  					     There is a narrow distribution of
					     equi-distant distances between
					     the first and second domain. */
  WLZ_RCC_MSK    = ((1<<WLZ_RCCIDX_CNT)-1) /*!< Not a clasification but a bit
  					     mask for all the possible
					     classifications. */
} WlzRCCClass;

/*!
* \enum		_WlzSpecialStructElmType
* \ingroup	WlzMorphologyOps
* \brief	Special structuring elements for morphological operations.
*		Typedef: ::WlzSpecialStructElmType.
*/
typedef enum _WlzSpecialStructElmType
{
  WLZ_SPEC_STRUCT_ELM_H4,
  WLZ_SPEC_STRUCT_ELM_EX4,
  WLZ_SPEC_STRUCT_ELM_A8,
  WLZ_SPEC_STRUCT_ELM_H6,
  WLZ_SPEC_STRUCT_ELM_H5,
  WLZ_SPEC_STRUCT_ELM_H7,
  WLZ_SPEC_STRUCT_ELM_A3,
  WLZ_SPEC_STRUCT_ELM_E1,
  WLZ_SPEC_STRUCT_ELM_E2,
  WLZ_SPEC_STRUCT_ELM_V2
} WlzSpecialStructElmType;

/*!
* \enum		_WlzBinaryOperatorType
* \ingroup	WlzArithmetic
* \brief	 Binary operators.						
*		Typedef: ::WlzBinaryOperatorType.
*/
typedef enum _WlzBinaryOperatorType
{
  WLZ_BO_ADD		= 0,
  WLZ_BO_SUBTRACT,
  WLZ_BO_MULTIPLY,
  WLZ_BO_DIVIDE,
  WLZ_BO_MODULUS,
  WLZ_BO_EQ,
  WLZ_BO_NE,
  WLZ_BO_GT,
  WLZ_BO_GE,
  WLZ_BO_LT,
  WLZ_BO_LE,
  WLZ_BO_AND,
  WLZ_BO_OR,
  WLZ_BO_XOR,
  WLZ_BO_MAX,
  WLZ_BO_MIN,
  WLZ_BO_MAGNITUDE
} WlzBinaryOperatorType;

/*!
* \enum 	_WlzCompThreshType
* \ingroup	WlzThreshold
* \brief	Automatic threshold computation methods.			
*  		The histogram may need to be smoothed for these algorithms
*		to work.
		Typedef: ::WlzCompThreshType.
*/
typedef enum _WlzCompThreshType
{
  WLZ_COMPTHRESH_FOOT,	/*!< Can not determine threshold type.
  			   The threshold value is intercept of a line fitted
			   to the upper slope of the histogram main peak with
			   the abscissa.  */
  WLZ_COMPTHRESH_DEPTH, /*!< Can not determine threshold type.
  			   The threshold value is that point to the right of
  			   the histogram peak that is maximally distant from
			   the chord joining the peak and the histogram right
			   hand end point.
			   The histogram may need to be smoothed for this
			   algorithm to work.  */
  WLZ_COMPTHRESH_GRADIENT, /*!< Can not determine threshold type.
  			   The threshold value is the first point to the
			   right of the histogram main peak at which the
			   gradient falls to zero (cf finding a minimum).
			   To find the slope of the histogram at some
			   point a straight line is fitted through the
			   point \f$\pm\f$ a fixed number of points to
			   either side. */
  WLZ_COMPTHRESH_FRACMIN, /*!< The threshold type is determined from the
  			   position of the histogram's maximum and the
			   threshold value is the minimum following
                           the specified fraction of the values.  */
  WLZ_COMPTHRESH_SMOOTHSPLIT, /*!< The threshold value is found by
  			   heavily smoothing the histogram and looking
			   for the minimum. Successively lesser smoothing
			   values are then applied and at each iteration
			   the minimum closest to the previous is found. */
  WLZ_COMPTHRESH_OTSU	   /*!< The threshold value is found by using
  			   Otsu's method. This is a clustering-based algorithm
			   which computes an optimum threshold value that
			   separates the two classes of an (assumed)
			   bi-modal histogram. */
} WlzCompThreshType;

/*!
* \enum		_WlzInterpolationType
* \ingroup	WlzType
* \brief	Interpolation methods.
*		Typedef: ::WlzInterpolationType.
*/
typedef enum _WlzInterpolationType
{
  WLZ_INTERPOLATION_NEAREST     = 0,	/*!< Nearest neighbour. */
  WLZ_INTERPOLATION_LINEAR,		/*!< Linear or tri-linear. */
  WLZ_INTERPOLATION_CLASSIFY_1,		/*!< Classification by probability. */
  WLZ_INTERPOLATION_CALLBACK,		/*!< Callback function computes
					     each interpolated value. */
  WLZ_INTERPOLATION_ORDER_2,		/*!< Second order interpolation. */
  WLZ_INTERPOLATION_BARYCENTRIC,	/*!< Barycentric mesh interpolation. */
  WLZ_INTERPOLATION_KRIG 	        /*!< Kriging mesh interpolation. */
} WlzInterpolationType;

/*!
* \enum		_WlzThresholdType
* \ingroup	WlzThreshold
* \brief	Threshold value selection.
* 		Typedef: ::WlzThresholdType.
*/
typedef enum _WlzThresholdType
{
  WLZ_THRESH_LOW		= 0, /*!< Threshold <  given value */
  WLZ_THRESH_HIGH,		     /*!< Threshold >= given value */
  WLZ_THRESH_EQUAL		     /*!< Threshold == given value */
} WlzThresholdType;


/*!
* \enum		_WlzRGBAThresholdType
* \ingroup	WlzThreshold
* \brief	Colour threshold type selection.
* 		Typedef: ::WlzRGBAThresholdType.
*/
typedef enum _WlzRGBAThresholdType
{
  WLZ_RGBA_THRESH_NONE,
  WLZ_RGBA_THRESH_SINGLE,
  WLZ_RGBA_THRESH_MULTI,
  WLZ_RGBA_THRESH_PLANE,
  WLZ_RGBA_THRESH_SLICE,
  WLZ_RGBA_THRESH_BOX,
  WLZ_RGBA_THRESH_SPHERE
} WlzRGBAThresholdType;

/*!
* \enum		_WlzRGBAColorSpace
* \ingroup	WlzValueUtils
* \brief	Colour space (i.e. rgb, hsb, grey etc.) selection.
* 		Typedef: ::WlzRGBAColorSpace.
*/
typedef enum _WlzRGBAColorSpace
{
  WLZ_RGBA_SPACE_GREY,
  WLZ_RGBA_SPACE_RGB,
  WLZ_RGBA_SPACE_HSB,
  WLZ_RGBA_SPACE_CMY
} WlzRGBAColorSpace;

/*!
* \enum		_WlzRGBAColorChannel
* \ingroup	WlzValueUtils
* \brief	Colour channel selection.
* 		Typedef: ::WlzRGBAColorChannel.
*/
typedef enum _WlzRGBAColorChannel
{
  WLZ_RGBA_CHANNEL_GREY,
  WLZ_RGBA_CHANNEL_RED,
  WLZ_RGBA_CHANNEL_GREEN,
  WLZ_RGBA_CHANNEL_BLUE,
  WLZ_RGBA_CHANNEL_HUE,
  WLZ_RGBA_CHANNEL_SATURATION,
  WLZ_RGBA_CHANNEL_BRIGHTNESS,
  WLZ_RGBA_CHANNEL_CYAN,
  WLZ_RGBA_CHANNEL_MAGENTA,
  WLZ_RGBA_CHANNEL_YELLOW,
  WLZ_RGBA_CHANNEL_DUMMY
} WlzRGBAColorChannel;

/*!
* \enum		_WlzPolyFillMode	
* \ingroup	WlzPolyline
* \brief	Polygon fill modes.
*		Typedef: ::WlzPolyFillMode.
*/
typedef enum _WlzPolyFillMode
{
  WLZ_SIMPLE_FILL,	/*!< Fill all pixels with winding number > 0 */
  WLZ_EVEN_ODD_FILL,	/*!< Fill all pixels with odd winding number */
  WLZ_VERTEX_FILL	/*!< Fill all pixels lying under the polyline */
} WlzPolyFillMode;

/*!
* \enum		_WlzGreyTransformType	
* \ingroup	WlzTransform
* \brief	Grey-level transform types.
*		Typedef: ::WlzGreyTransformType.
*/
typedef enum _WlzGreyTransformType {
  WLZ_GREYTRANSFORMTYPE_IDENTITY,	/*!< No value change. */
  WLZ_GREYTRANSFORMTYPE_LINEAR,		/*!< Linear interpolation. */
  WLZ_GREYTRANSFORMTYPE_GAMMA,		/*!< Gamma function. */
  WLZ_GREYTRANSFORMTYPE_SIGMOID		/*!< Sigmoid function. */
} WlzGreyTransformType;

/*!
* \enum		_WlzThreeDStdViews
* \ingroup	WlzSectionTransform
* \brief	Standard 3D views.
		Typedef: ::WlzThreeDStdViews.
*/
typedef enum _WlzThreeDStdViews
{
  WLZ_X_Y_VIEW,
  WLZ_Y_Z_VIEW,
  WLZ_Z_X_VIEW,
  WLZ_ARBITRARY_VIEW
} WlzThreeDStdViews;

/*!
* \enum		_WlzThreeDViewMode
* \ingroup	WlzSectionTransform
* \brief	3D viewing modes which determine the angle at which the plane
*		is cut through a 3D volume.
*		Typedef: ::WlzThreeDViewMode.
*/
typedef enum _WlzThreeDViewMode
{
  WLZ_STATUE_MODE,			/*!< Corresponds to "walking around a
  					     statue" so that if the view is say
					     from the left-hand side then the
					     section will have the top of the
					     statue to the left. */
  WLZ_UP_IS_UP_MODE,			/*!< The projection of the vector up
  					     onto the section image will be
					     "up". This is ill-defined if the
					     viewing direction is very close to
					     "up". */
  WLZ_FIXED_LINE_MODE,
  WLZ_ZERO_ZETA_MODE,
  WLZ_ZETA_MODE
} WlzThreeDViewMode;

/*!
* \enum		_WlzWindowFnType
* \ingroup	WlzValuesFilters
* \brief	Types of window functions.
*		Window functions are used to weight the grey values of Woolz
*		domain objects according to the spatial distribution.
*		Typedef: WlzWindowFnType.
*
*               The window functions are defined as:
*                 \f[
                    blackman(r,R) = 0.42 +
                                    0.50 \cos{(\frac{2.00 \pi r}{R - 1.00})} +
                                    0.08 \cos{(\frac{4.00 \pi r}{R - 1.00})}
                  \f]
                  \f[
                    hamming(r,R) = 0.54 +
                                   0.46 \cos{(\frac{2.00 \pi r}{R - 1.00})}
                  \f]
                  \f[
                    hanning(r,R) = 0.50 +
                                   0.50 \cos{(\frac{2.00 \pi r}{R - 1.00})}
                  \f]
                  \f[
                    parzen(r,R) = 1.00 - 
                                  |frac{2.00 r - (R - 1.00)}{R + 1.00}|
                  \f]
                  \f[
                    welch(r,R) = 1.00 -
                                 {(\frac{2.00 r - (R - 1.00)}{R + 1.00})}^2
                  \f]
                  \f[
                    rectangle(r,R) = \left\{
                                     \begin{array}{ll}
                                     1.00 & |r_x| \leq R, |r_y| \leq R \\
                                     0.00 & |r_x| > R, |r_y| > R
                                     \end{array}
                                     \right.
                  \f]
*               Where \f$R\f$ is the window radius and \f$r\f$ the radial
*               distance from the centre of the window.
*/
typedef enum _WlzWindowFnType
{
  WLZ_WINDOWFN_NONE,			/*!< No window function. */
  WLZ_WINDOWFN_BLACKMAN,		/*!< Blackman window. */
  WLZ_WINDOWFN_HAMMING,			/*!< Hanning window. */
  WLZ_WINDOWFN_HANNING,
  WLZ_WINDOWFN_PARZEN,			/*!< Parzen window. */
  WLZ_WINDOWFN_RECTANGLE,		/*!< Rectangular window. */
  WLZ_WINDOWFN_WELCH,			/*!< Welch window. */
  WLZ_WINDOWFN_UNSPECIFIED		/*!< Not a window function also has
  					     a value equal to the number of
					     window functions defined because
					     it is the last in the enum. */
} WlzWindowFnType;

/*!
* \enum		_WlzSampleFn
* \ingroup	WlzTransform
* \brief	Sampling functions.						
* 		Typedef: ::WlzSampleFn.
*/
typedef enum _WlzSampleFn
{
  WLZ_SAMPLEFN_NONE, 			/*!< No sampling function */
  WLZ_SAMPLEFN_POINT,			/*!< Point sampling */
  WLZ_SAMPLEFN_MEAN,			/*!< Mean value sample of data */
  WLZ_SAMPLEFN_GAUSS,	                /*!< Gaussian weighted sample of data */
  WLZ_SAMPLEFN_MIN,			/*!< Minimum value sampling */
  WLZ_SAMPLEFN_MAX,			/*!< Maximum value sampling */
  WLZ_SAMPLEFN_MEDIAN			/*!< Median value sampling */
} WlzSampleFn;

/*!
* \enum		_WlzScalarFeatureType
* \ingroup	WlzFeatures
* \brief	Scalar features of objects.
*/
typedef enum _WlzScalarFeatureType
{
  WLZ_SCALARFEATURE_VALUE,		/*!< Grey value. */
  WLZ_SCALARFEATURE_GRADIENT		/*!< Gradient of grey values. */
} WlzScalarFeatureType;

/*!
* \enum		_WlzDGTensorFeatureType
* \ingroup	WlzFeatures
* \brief	Features of Jacobian deformation gradient tensors.
*/
typedef enum _WlzDGTensorFeatureType
{
  WLZ_DGTENSOR_FEATURE_NONE	= 0,	/*!< No feature, null case. */
  WLZ_DGTENSOR_FEATURE_DETJAC	= 1,	/*!< Determinant of Jacobian tensor. */
  WLZ_DGTENSOR_FEATURE_EIGENVEC	= 2, 	/*!< Eigen vectors. */
  WLZ_DGTENSOR_FEATURE_EIGENVAL	= 3, 	/*!< Eigen values. */
  WLZ_DGTENSOR_FEATURE_LIMIT    = 4  	/*!< Not a feature but used to
                                             itterate through the features,
					     keep it the largest value. */
} WlzDGTensorFeatureType;

/*!
* \def		WLZ_DGTENSOR_FEATURE_MASK
* \ingroup	WlzFeatures
* \brief	Bit mask generator for features of Jacobian deformation
* 		gradient tensors.
*/
#define WLZ_DGTENSOR_FEATURE_MASK(F)	(((F)>0)?(1<<((F)-1)):(0))

/*!
* \enum		_WlzVertexType
* \ingroup	WlzType
* \brief	2D and 3D vertex types.
*		Typedef: ::WlzVertexType.
*/
typedef enum _WlzVertexType
{
  WLZ_VERTEX_ERROR	= 0,		/*!< Used to indicate not a vertex
  					     type. */
  WLZ_VERTEX_I2		= 1,		/*!< 2D integer vertex. */
  WLZ_VERTEX_F2,			/*!< 2D single precision floating
  					     point vertex. */
  WLZ_VERTEX_D2,			/*!< 2D double precision floating
  					     point vertex. */
  WLZ_VERTEX_I3,			/*!< 3D integer vertex. */
  WLZ_VERTEX_F3,			/*!< 3D single precision floating
  					     point vertex. */
  WLZ_VERTEX_D3,			/*!< 3D double precision floating
  					     point vertex. */
  WLZ_VERTEX_L2,			/*!< 2D long integer vertex. */
  WLZ_VERTEX_L3				/*!< 3D long integer vertex. */
} WlzVertexType;

/*!
* \struct	_WlzLVertex2
* \ingroup	WlzType
* \brief	2D long integer vertex.
* 		Typedef: ::WlzIVertex2.
*/
typedef struct _WlzLVertex2
{
  WlzLong   	vtY;
  WlzLong   	vtX;
} WlzLVertex2;

/*!
* \struct	_WlzIVertex2
* \ingroup	WlzType
* \brief	2D integer vertex.
* 		Typedef: ::WlzIVertex2.
*/
typedef struct _WlzIVertex2
{
  int   	vtY;
  int   	vtX;
} WlzIVertex2;

/*!
* \struct	_WlzFVertex2
* \ingroup	WlzType
* \brief	2D single precision float point vertex.
* 		Typedef: ::WlzFVertex2.
*/
typedef struct _WlzFVertex2
{
  float 	vtY;
  float 	vtX;
} WlzFVertex2;

/*!
* \struct	_WlzDVertex2
* \ingroup	WlzType
* \brief	2D double precision float point vertex.
* 		Typedef: ::WlzDVertex2.
*/
typedef struct _WlzDVertex2
{
  double 	vtY;
  double 	vtX;
} WlzDVertex2;


/*!
* \struct	_WlzLVertex3
* \ingroup	WlzType
* \brief	3D long integer vertex.
* 		Typedef: ::WlzLVertex3.
*/
typedef struct _WlzLVertex3
{
  WlzLong	vtX;
  WlzLong	vtY;
  WlzLong	vtZ;
} WlzLVertex3;

/*!
* \struct	_WlzIVertex3
* \ingroup	WlzType
* \brief	3D integer vertex.
* 		Typedef: ::WlzIVertex3.
*/
typedef struct _WlzIVertex3
{
  int		vtX;
  int		vtY;
  int		vtZ;
} WlzIVertex3;

/*!
* \struct       _WlzFVertex3
* \ingroup	WlzType
* \brief	3D single precision float point vertex.
*               Typedef: ::WlzFVertex3.
*/
typedef struct _WlzFVertex3
{
  float		vtX;
  float		vtY;
  float		vtZ;
} WlzFVertex3;

/*!
* \struct       _WlzDVertex3
* \ingroup	WlzType
* \brief	3D double precision float point vertex.
*               Typedef: ::WlzDVertex3.
*/
typedef struct _WlzDVertex3
{
  double	vtX;
  double	vtY;
  double	vtZ;
} WlzDVertex3;

/*!
* \union	_WlzVertexP
* \ingroup	WlzType
* \brief	Union of vertex pointers.
*		Typedef: ::WlzVertexP.
*/
typedef union _WlzVertexP
{
  void		*v;
  WlzIVertex2	*i2;
  WlzLVertex2	*l2;
  WlzFVertex2	*f2;
  WlzDVertex2	*d2;
  WlzIVertex3	*i3;
  WlzLVertex3	*l3;
  WlzFVertex3	*f3;
  WlzDVertex3	*d3;
} WlzVertexP;

/*!
* \union	_WlzVertex
* \ingroup	WlzType
* \brief	Union of vertex values.
*		Typedef: ::WlzVertex.
*/
typedef union _WlzVertex
{
  WlzIVertex2	i2;
  WlzLVertex2	l2;
  WlzFVertex2	f2;
  WlzDVertex2	d2;
  WlzIVertex3	i3;
  WlzLVertex3	l3;
  WlzFVertex3	f3;
  WlzDVertex3	d3;
} WlzVertex;

/*!
* \struct	_WlzIBox2
* \ingroup	WlzType
* \brief	2D integer axis aligned rectangle (box).
* 		Typedef: ::WlzIBox2.
*/
typedef struct _WlzIBox2
{
  int		xMin;
  int		yMin;
  int		xMax;
  int		yMax;
} WlzIBox2;

/*!
* \struct       _WlzDBox2
* \ingroup	WlzType
* \brief	2D double precision floating point axis aligned rectangle
* 		(box).
*		Typedef: ::WlzDBox2.
*/
typedef struct _WlzDBox2
{
  double	xMin;
  double	yMin;
  double	xMax;
  double	yMax;
} WlzDBox2;

/*!
* \struct       _WlzFBox2
* \ingroup	WlzType
* \brief	2D single precision floating point axis aligned rectangle
* 		(box).
*		Typedef: ::WlzFBox2.
*/
typedef struct _WlzFBox2
{
  float		xMin;
  float		yMin;
  float		xMax;
  float		yMax;
} WlzFBox2;

/*!
* \struct	_WlzIBox3
* \ingroup	WlzType
* \brief	3D integer axis aligned cubiod (box).
*		Typedef: ::WlzIBox3.
*/
typedef struct _WlzIBox3
{
  int		xMin;
  int		yMin;
  int		zMin;
  int		xMax;
  int		yMax;
  int		zMax;
} WlzIBox3;

/*!
* \struct	_WlzDBox3
* \ingroup	WlzType
* \brief	3D double precision floating point axis aligned cubiod (box).
*		Typedef: ::WlzDBox3.
*/
typedef struct _WlzDBox3
{
  double	xMin;
  double	yMin;
  double	zMin;
  double	xMax;
  double	yMax;
  double	zMax;
} WlzDBox3;

/*!
* \struct	_WlzFBox3
* \ingroup	WlzType
* \brief	3D single precision floating point axis aligned cubiod (box).
*		Typedef: ::WlzFBox3.
*/
typedef struct _WlzFBox3
{
  float		xMin;
  float		yMin;
  float		zMin;
  float		xMax;
  float		yMax;
  float		zMax;
} WlzFBox3;

/*!
* \union	_WlzBoxP
* \ingroup	WlzType
* \brief	Union of axis aligned box pointers.
* 		Typedef: ::WlzBoxP.
*/
typedef union _WlzBoxP
{
  void		*v;
  WlzIBox2	*i2;
  WlzFBox2	*f2;
  WlzDBox2	*d2;
  WlzIBox3	*i3;
  WlzFBox3	*f3;
  WlzDBox3	*d3;
} WlzBoxP;

/*!
* \union	_WlzBox
* \ingroup	WlzType
* \brief	Union of axis aligned boxes.
*		Typedef: ::WlzBox.
*/
typedef union _WlzBox
{
  WlzIBox2	i2;
  WlzFBox2	f2;
  WlzDBox2	d2;
  WlzIBox3	i3;
  WlzFBox3	f3;
  WlzDBox3	d3;
} WlzBox;

/************************************************************************
* Grey values.
************************************************************************/

/*!
* \union	_WlzGreyP
* \ingroup	WlzType
* \brief	A union of pointers to grey values.
*		Typedef: ::WlzGreyP.
*/
typedef union _WlzGreyP
{
  void    	*v; 			/*!< Can save a cast when assigning. */
  WlzLong 	*lnp;
  int     	*inp;
  short   	*shp;
  WlzUByte 	*ubp;
  float   	*flp;
  double  	*dbp;
  WlzUInt 	*rgbp;
  char		**bytes;
  unsigned char **ubytes;
} WlzGreyP;

/*!
* \union	_WlzGreyV
* \ingroup	WlzType
* \brief	A union of grey values.
*		Typedef: ::WlzGreyV.
*/
typedef union _WlzGreyV
{
  void		*v;		        /*!< Can save a cast when assigning. */
  WlzLong 	lnv;
  int 		inv;
  short 	shv;
  WlzUByte 	ubv;
  float 	flv;
  double 	dbv;
  WlzUInt 	rgbv;
  char		bytes[8];
  unsigned char	ubytes[8];
} WlzGreyV;

/*!
* \struct	_WlzPixelV
* \ingroup	WlzType
* \brief	A typed grey value.
*		Typedef: ::WlzPixelV.
*/
typedef struct _WlzPixelV
{
  WlzGreyType	type;			/*!< Type of grey value. */
  WlzGreyV	v;			/*!< The grey value. */
} WlzPixelV;

/*!
* \struct	_WlzPixelP
* \ingroup	WlzType
* \brief	A typed grey value pointer.
*		Typedef: ::WlzPixelP.
*/
typedef struct _WlzPixelP
{
  WlzGreyType	type;			/*!< Type of grey value. */
  WlzGreyP	p;			/*!< Pointer to the grey value(s). */
} WlzPixelP;

/*!
* \struct	_WlzGreyTransformParam
* \ingroup	WlzTransform
* \brief	Grey-level transform parameters.
*		Typedef: ::WlzGreyTransformParam.
*/
typedef struct _WlzGreyTransformParam
{
  WlzGreyTransformType type; 		/*!< Grey transform type. */
  WlzPixelV	il;       		/*!< Input minimum grey value. */
  WlzPixelV	iu;       		/*!< Input maximum grey value. */
  WlzPixelV	ol;       		/*!< Output minimum grey value. */
  WlzPixelV	ou;       		/*!< Output maximum grey value. */
  double	p0;			/*!< First parameter, used for
                                             gamma (\f$\gamma\f$) or
                                             sigmoid (\f$\mu\f$). */
  double	p1;			/*!< Second parameter, used for
                                             sigmoid (\f$\sigma\f$). */
} WlzGreyTransformParam;

/************************************************************************
* Markers.
************************************************************************/
/*!
* \enum		_WlzMarkerType
* \ingroup	WlzType
* \brief	Basic markers.
* 		Typedef: ::WlzMarkerType.
*/
typedef enum _WlzMarkerType
{
  WLZ_MARKER_NONE       = 0,		/*!< No marker. */
  WLZ_MARKER_POINT,			/*!< Single point, pixel or voxel. */
  WLZ_MARKER_SPHERE			/*!< Circle or sphere. */
} WlzMarkerType;


/************************************************************************
* Data structures for geometric models.
************************************************************************/

/*!
* \def		WLZ_GM_TOLERANCE
* \ingroup	WlzGeoModel
* \brief	Tolerance for geometric queries in geometric models.
*/
#define	WLZ_GM_TOLERANCE	(1.0e-06)

/*!
* \def		WLZ_GM_TOLERANCE_SQ
* \ingroup	WlzGeoModel
* \brief	Square of tolerance for geometric queries in geometric models.
*/
#define	WLZ_GM_TOLERANCE_SQ	(WLZ_GM_TOLERANCE * WLZ_GM_TOLERANCE)

/*!
* \enum		_WlzGMModelType
* \ingroup	WlzGeoModel
* \brief	Types of geometric models.
*		Typedef: WlzGMModelType.
*/
typedef enum _WlzGMModelType
{
  WLZ_GMMOD_2I = 1,
  WLZ_GMMOD_2D,
  WLZ_GMMOD_3I,
  WLZ_GMMOD_3D,
  WLZ_GMMOD_2N,
  WLZ_GMMOD_3N
} WlzGMModelType;

/*!
* \enum		_WlzGMElemType
* \ingroup	WlzGeoModel
* \brief        Types of geometric model elements.
*               Typedef: ::WlzGMElemType.
*/
typedef enum _WlzGMElemType
{
  WLZ_GMELM_NONE = 0,
  WLZ_GMELM_VERTEX,
  WLZ_GMELM_VERTEX_G2I,
  WLZ_GMELM_VERTEX_G2D,
  WLZ_GMELM_VERTEX_G2N,
  WLZ_GMELM_VERTEX_G3I,
  WLZ_GMELM_VERTEX_G3D,
  WLZ_GMELM_VERTEX_G3N,
  WLZ_GMELM_VERTEX_T,
  WLZ_GMELM_DISK_T,
  WLZ_GMELM_EDGE,
  WLZ_GMELM_EDGE_T,
  WLZ_GMELM_FACE,
  WLZ_GMELM_LOOP_T,
  WLZ_GMELM_SHELL,
  WLZ_GMELM_SHELL_G2I,
  WLZ_GMELM_SHELL_G2D,
  WLZ_GMELM_SHELL_G3I,
  WLZ_GMELM_SHELL_G3D
} WlzGMElemType;

/*!
* \enum 	_WlzGMElemTypeFlags
* \ingroup	WlzGeoModel
* \brief	Bit masks for the types of geometric model elements.
*               Typedef: ::WlzGMElemTypeFlags.
*/
typedef enum _WlzGMElemTypeFlags
{
  WLZ_GMELMFLG_VERTEX =	    (1 << 0),
  WLZ_GMELMFLG_VERTEX_G =   (1 << 1),
  WLZ_GMELMFLG_VERTEX_T =   (1 << 2),
  WLZ_GMELMFLG_DISK_T =     (1 << 3),
  WLZ_GMELMFLG_EDGE =       (1 << 4),
  WLZ_GMELMFLG_EDGE_T =     (1 << 5),
  WLZ_GMELMFLG_FACE =       (1 << 6),
  WLZ_GMELMFLG_LOOP_T =     (1 << 7),
  WLZ_GMELMFLG_SHELL =      (1 << 8),
  WLZ_GMELMFLG_SHELL_G =    (1 << 9)
} WlzGMElemTypeFlags;

/*!
* \union	_WlzGMElemP
* \ingroup	WlzGeoModel
* \brief	A union of pointers to all the valid geometric model
*		elements.
*		Typedef: ::WlzGMElemP.
*/
typedef union _WlzGMElemP
{
  struct _WlzGMCore *core;
  struct _WlzGMVertex *vertex;
  struct _WlzGMVertexG2I *vertexG2I;
  struct _WlzGMVertexG2D *vertexG2D;
  struct _WlzGMVertexG2N *vertexG2N;
  struct _WlzGMVertexG3I *vertexG3I;
  struct _WlzGMVertexG3D *vertexG3D;
  struct _WlzGMVertexG3N *vertexG3N;
  struct _WlzGMVertexT *vertexT;
  struct _WlzGMDiskT	*diskT;
  struct _WlzGMEdge 	*edge;
  struct _WlzGMEdgeT 	*edgeT;
  struct _WlzGMFace 	*face;
  struct _WlzGMLoopT 	*loopT;
  struct _WlzGMShell 	*shell;
  struct _WlzGMShellG2I *shellG2I;
  struct _WlzGMShellG2D *shellG2D;
  struct _WlzGMShellG3I *shellG3I;
  struct _WlzGMShellG3D *shellG3D;
} WlzGMElemP;

/*!
* \struct	_WlzGMCore
* \ingroup	WlzGeoModel
* \brief	The core geometric model element from which all
* 		geometric modeling elements inherit the type and
*		index fields.
*		Typedef: ::WlzGMCore.
*/
typedef struct _WlzGMCore
{
  WlzGMElemType type;			/*!< Any WlzGMElemType */
  int		idx;			/*!< Unique identifier for the
  					     instance of this type in it's
					     model. */
} WlzGMCore;

/*!
* \struct	_WlzGMVertexG2I
* \ingroup	WlzGeoModel
* \brief	The position of a point in 2D integer space.
*		Typedef: ::WlzGMVertexG2I.
*/
typedef struct _WlzGMVertexG2I
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_G2I */
  int		idx;		    	/*!< Unique identifier for the
  					     vertex geometry. */
  WlzIVertex2	vtx;			/*!< Where the point lies in space. */
} WlzGMVertexG2I;

/*!
* \struct	_WlzGMVertexG2D
* \ingroup	WlzGeoModel
* \brief	The position of a point in 2D double precision space.
*		Typedef: ::WlzGMVertexG2D.
*/
typedef struct _WlzGMVertexG2D
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_G2D */
  int		idx;		    	/*!< Unique identifier for vertex
  					     geometry. */
  WlzDVertex2	vtx;			/*!< Where the point lies in space. */
} WlzGMVertexG2D;

/*!
* \struct	_WlzGMVertexG2N
* \ingroup	WlzGeoModel
* \brief	The position of a point in 2D double precision space
*		and the normal vector at that point.
*		Note that the data structure is the same as ::WlzGMVertexG2D
*		until the normal, this is important as it allows the type,
*		index and position to be established without knowing whether
*		the vertex geometry is ::WlzGMVertexG2D or ::WlzGMVertexG2N.
*		Typedef: ::WlzGMVertexG2N.
*/
typedef struct _WlzGMVertexG2N
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_G2N */
  int		idx;		        /*!< Unique identifier for vertex
  				             geometry */
  WlzDVertex2	vtx;			/*!< Where the point lies in space */
  WlzDVertex2	nrm;			/*!< Normal at the point. */
} WlzGMVertexG2N;

/*!
* \struct	_WlzGMVertexG3I
* \ingroup	WlzGeoModel
* \brief	The position of a point in 3D integer space.
*		Typedef: ::WlzGMVertexG3I.
*/
typedef struct _WlzGMVertexG3I
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_G3I */
  int		idx;		    	/*!< Unique identifier for vertex
  					     geometry. */
  WlzIVertex3	vtx;			/*!< Where the point lies in space. */
} WlzGMVertexG3I;

/*!
* \struct	_WlzGMVertexG3D
* \ingroup	WlzGeoModel
* \brief	The position of a point in 3D double precision space.
*		Typedef: ::WlzGMVertexG3D.
*/
typedef struct _WlzGMVertexG3D
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_G3D */
  int		idx;		        /*!< Unique identifier for vertex
  				             geometry */
  WlzDVertex3	vtx;			/*!< Where the point lies in space */
} WlzGMVertexG3D;

/*!
* \struct	_WlzGMVertexG3N
* \ingroup	WlzGeoModel
* \brief	The position of a point in 3D double precision space
*		and the normal vector at that point.
*		Note that the data structure is the same as ::WlzGMVertexG3D
*		until the normal, this is important as it allows the type,
*		index and position to be established without knowing whether
*		the vertex geometry is ::WlzGMVertexG3D or ::WlzGMVertexG3N.
*		Typedef: ::WlzGMVertexG3N.
*/
typedef struct _WlzGMVertexG3N
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_G3N */
  int		idx;		        /*!< Unique identifier for vertex
  				             geometry */
  WlzDVertex3	vtx;			/*!< Where the point lies in space */
  WlzDVertex3	nrm;			/*!< Normal at the point. */
} WlzGMVertexG3N;

/*!
* \union	_WlzGMVertexGU
* \ingroup	WlzGeoModel
* \brief	A union of pointers to the geometric properties of a point.
*		Typedef: ::WlzGMVertexGU.
*/
typedef union _WlzGMVertexGU
{
  WlzGMCore 	*core;
  WlzGMVertexG2I *vg2I;
  WlzGMVertexG2D *vg2D;
  WlzGMVertexG2N *vg2N;
  WlzGMVertexG3I *vg3I;
  WlzGMVertexG3D *vg3D;
  WlzGMVertexG3N *vg3N;
} WlzGMVertexGU;

/*!
* \struct	_WlzGMVertexT
* \ingroup	WlzGeoModel
* \brief	The topological properties of a point in space.
*		The ordering of the linked list of vertex topology elements
*		formed by the 'next' and 'prev' pointers is not significant.
*		Typedef: ::WlzGMVertexT.
*/
typedef struct _WlzGMVertexT
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_T */
  int	 	idx; 	    	        /*!< Unique identifier for vertex
  					     topology element. */
  struct _WlzGMVertexT *next;		/*!< Next vertexT in disk. */
  struct _WlzGMVertexT *prev;		/*!< Previous vertexT in disk. */
  struct _WlzGMDiskT *diskT;		/*!< The disk topology element that
  					     this vertex topology element
					     is in. */
  struct _WlzGMEdgeT *parent; 		/*!< Parent of this vertex topology
  					     element. */
} WlzGMVertexT;

/*!
* \struct	_WlzGMVertex
* \ingroup	WlzGeoModel
* \brief	A single point in space defined in terms of both it's
*		geometry and it's topology.
*		Typedef: ::WlzGMVertex.
*/
typedef struct _WlzGMVertex
{
  WlzGMElemType type; 			/*!< WLZ_GMELM_VERTEX */
  int		idx;			/*!< Unique identifier for vertex. */
  struct _WlzGMDiskT *diskT;		/*!< A disk topology element of this
  					     vertex, others can be found by
					     following the diskT's next/prev
					     fields. */
  WlzGMVertexGU geo;	 		/*!< Geometry of this vertex. */
  struct _WlzGMVertex *next;		/*!< Next in sorted list. */
} WlzGMVertex;

/*!
* \struct	_WlzGMDiskT
* \ingroup	WlzGeoModel
* \brief	A topological disk around a vertex. In 2D or 3D manifold
*		there is one disk per vertex. But in a 3D non-manifold shell
*		many sheets (manifold surfaces components) may be connected
*		at a single vertex, in which case there is one disk per sheet.
*		The disk encodes the radial order of the vertex topology
*		elements around the vertex.
*		Typedef: ::WlzGMDiskT.
*/
typedef struct _WlzGMDiskT
{
  WlzGMElemType type; 			/*!< WLZ_GMELM_DISK_T */
  int		idx;			/*!< Unique identifier for vertex. */
  struct _WlzGMDiskT *next;		/*!< Next diskT of vertex. */
  struct _WlzGMDiskT *prev;		/*!< Previous diskT of vertex. */
  WlzGMVertex	*vertex;		/*!< The vertex that this disk cycles
  					     around. */
  WlzGMVertexT	*vertexT;		/*!< A vertex topology element in this
  					     disk topology element. */
} WlzGMDiskT;

/*!
* \struct	_WlzGMEdgeT
* \ingroup	WlzGeoModel
* \brief	The topological properties of a directed edge.
*		Typedef: ::WlzGMEdgeT.
*/
typedef struct _WlzGMEdgeT
{
  WlzGMElemType type;			/*!< WLZ_GMELM_EDGE_T */
  int		idx;	      		/*!< Unique identifier for the edge
  					 *   topology element. */
  struct _WlzGMEdgeT *next;		/*!< Next edgeT in the parent. */
  struct _WlzGMEdgeT *prev;		/*!< Previous edgeT in the parent. */
  struct _WlzGMEdgeT *opp;		/*!< Opposite edge topology element. */
  struct _WlzGMEdgeT *rad;		/*!< The radial edge topology
  					     element. */
  struct _WlzGMEdge *edge;		/*!< The edge. */
  struct _WlzGMVertexT *vertexT;  	/*!< Vertex FROM which this edge
   					     topology element is directed. */
  struct _WlzGMLoopT *parent;		/*!< Parent of this edge topology
  					     element. */
} WlzGMEdgeT;

/*!
* \struct	_WlzGMEdge
* \ingroup	WlzGeoModel
* \brief	A line or curve between a pair of vertices.
*		Although this only has a topological component a geometric
*		component would allow curves to be represented.
*		Typedef: ::WlzGMEdge.
*/
typedef struct _WlzGMEdge
{
  WlzGMElemType type;	        	/*!< WLZ_GMELM_EDGE */
  int		idx;		        /*!< Unique identifier for edge. */
  WlzGMEdgeT	*edgeT;	       		/*!< One of the many edge topology
  					     elements from which the others
					     can be found by following it's
					     opp/rad fields. */
} WlzGMEdge;

/*!
* \struct	_WlzGMLoopT
* \ingroup	WlzGeoModel
* \brief	A directed loop or the topological properties of a loop.
*		Typedef: ::WlzGMLoopT.
*/
typedef struct _WlzGMLoopT
{
  WlzGMElemType type;			/*!< WLZ_GMELM_LOOP_T */
  int		idx;	      		/*!< Unique identifier for loop
  					     topology element. */
  struct _WlzGMLoopT *next;		/*!< The next loopT in the parent. */
  struct _WlzGMLoopT *prev;		/*!< The previous loopT in the
  					     parent. */
  struct _WlzGMLoopT *opp;	        /*!< The opposite loop topology
  					     element. */
  struct _WlzGMFace *face;		/*!< The face not used in 2D models. */
  WlzGMEdgeT *edgeT;			/*!< An edge topology element in
  					     this loop topology element, the
					     others can be found by walking
					     the edgeT's next/prev fields. */
  struct _WlzGMShell *parent;		/*!< Parent of this loopT. */
} WlzGMLoopT;

/*!
* \struct	_WlzGMFace
* \ingroup	WlzGeoModel
* \brief	A circuit of edges.
*		Typedef: ::WlzGMFace.
*/
typedef struct _WlzGMFace
{
  WlzGMElemType type;			/*!< WLZ_GMELM_FACE */
  int		idx;			/*!< Unique identifier for face. */
  WlzGMLoopT 	*loopT;			/*!< A directed loop topology element
  					     of the face. */
} WlzGMFace;

/*!
* \struct	_WlzGMShellG2I
* \ingroup	WlzGeoModel
* \brief	The geometric properties of a shell in 2D integer space.
*		Typedef: ::WlzGMShellG2I.
*/
typedef struct _WlzGMShellG2I
{
  WlzGMElemType type;			/*!< WLZ_GMELM_SHELL_G2I */
  int		idx;	     		/*!< Unique identifier for shell
  					     geometry element. */
  WlzIBox2	bBox;			/*!< The bounding box of the shell. */
} WlzGMShellG2I;

/*!
* \struct	_WlzGMShellG2D
* \ingroup	WlzGeoModel
* \brief        The geometric properties of a shell in 2D double precision
*		space.
*		Typedef: ::WlzGMShellG2D.
*/
typedef struct _WlzGMShellG2D
{
  WlzGMElemType type;			/*!< WLZ_GMELM_SHELL_G2D */
  int		idx;	     		/*!< Unique identifier for shell
  					     geometry element. */
  WlzDBox2	bBox;			/*!< Bounding box of the shell. */
} WlzGMShellG2D;

/*!
* \struct	_WlzGMShellG3I
* \ingroup	WlzGeoModel
* \brief	The geometric properties of a shell in 3D integer space.
*		Typedef: ::WlzGMShellG3I.
*/
typedef struct _WlzGMShellG3I
{
  WlzGMElemType type;			/*!< WLZ_GMELM_SHELL_G3I */
  int		idx;	     		/*!< Unique identifier for shell
  					     geometry element. */
  WlzIBox3	bBox;			/*!< Bounding box of the shell. */
} WlzGMShellG3I;

/*!
* \struct	_WlzGMShellG3D
* \ingroup	WlzGeoModel
* \brief	The geometric properties of a shell in 3D double precision.
*		space.
*		Typedef: ::WlzGMShellG3D.
*/
typedef struct _WlzGMShellG3D
{
  WlzGMElemType type;			/*!< WLZ_GMELM_SHELL_G3D */
  int		idx;	     		/*!< Unique identifier for shell
  					     geometry element. */
  WlzDBox3	bBox;			/*!< Bounding box of the shell. */
} WlzGMShellG3D;

/*!
* \union	_WlzGMShellGU
* \ingroup	WlzGeoModel
* \brief	A union of pointers to the geometric properties of a shell.
* 		Typedef: ::WlzGMShellGU.
*/
typedef union _WlzGMShellGU
{
  WlzGMCore 	*core;
  WlzGMShellG2I *sg2I;
  WlzGMShellG2D *sg2D;
  WlzGMShellG3I *sg3I;
  WlzGMShellG3D *sg3D;
} WlzGMShellGU;

/*!
* \struct	_WlzGMShell
* \ingroup	WlzGeoModel
* \brief	A shell which is a collection of connected geometric
*		modeling elements.
*		Typedef: ::WlzGMShell.
*/
typedef struct _WlzGMShell
{
  WlzGMElemType type;	      		/*!< WLZ_GMELM_SHELL */
  int		idx;			/*!< The shell's index. */
  struct _WlzGMShell *next;		/*!< the next shell in the model. */
  struct _WlzGMShell *prev;		/*!< The previous shell in the
  					     model. */
  WlzGMShellGU geo;			/*!< The shell's geometry. */
  WlzGMLoopT	*child;			/*!< A child loop topology element
  					     of the shell from which all
					     the others can be reached
					     by walking the next/prev
					     fields. */
  struct _WlzGMModel *parent;		/*!< The parent model of the shell. */
} WlzGMShell;

/*!
* \enum		_WlzGMCbReason
* \ingroup	WlzGeoModel
* \brief	The reason a callback function is called.
*		Typedef: ::WlzGMCbReason.
*/
typedef enum _WlzGMCbReason
{
  WLZ_GMCB_NEW,				/*!< New element has been created. */
  WLZ_GMCB_FREE				/*!< Existing element is about to
  					     be free'd. */
} WlzGMCbReason;

/*!
* \struct	_WlzGMResource
* \ingroup	WlzGeoModel
* \brief	A resource vector (extensible array) used for allocating
*		geometric modeling elements.
*		Typedef: ::WlzGMResource.
*/
typedef struct _WlzGMResource
{
  unsigned int  numElm;			/*!< Number of element type in model. */
  unsigned int	numIdx;			/*!< Number of elements/indicies which
   					     have been pulled from the vector,
  					     with elm->idx < numIdx for all
					     elements. */
  AlcVector	*vec;			/*!< Vector (extensible array) of
  					     elements. */
} WlzGMResource;

/*!
* \struct	_WlzGMModelR
* \ingroup	WlzGeoModel
* \brief	The resources used by a model.
*		Typedef: ::WlzGMModelR.
*/
typedef struct _WlzGMModelR
{
  struct _WlzGMCbEntry	*callbacks;		/*!< Linked list of functions which
  					    are called when new elements are
					    created or existing elements are
					    destroyed. */
  WlzGMResource vertex;			/*!< Vertex elements. */
  WlzGMResource vertexT;	        /*!< Vertex topology elements. */
  WlzGMResource vertexG;	     	/*!< Vertex geometry elements. */
  WlzGMResource diskT;			/*!< Disk geometry elements. */
  WlzGMResource edge;			/*!< Edge elements. */
  WlzGMResource edgeT;		        /*!< Edge topology element elements. */
  WlzGMResource face;			/*!< Face elements. */
  WlzGMResource loopT;		        /*!< Loop topology element elements. */
  WlzGMResource	shell;			/*!< Shell elements. */
  WlzGMResource	shellG;			/*!< Shell geometry elements. */
} WlzGMModelR;

/*!
* \struct	_WlzGMModel
* \ingroup	WlzGeoModel
* \brief	A geometric model which can represent both 2D graphs
*		and 3D surfaces, with the surfaces being either
*		manifold or non-manifold.
*		The geometric model inherits it's core fields from
*		the Woolz core domain.
*		Typedef: ::WlzGMModel.
*/
typedef struct _WlzGMModel
{
  WlzGMModelType type;			/*!< Type of model integer or double
  				         * precision, 2D or 3D */
  int		linkcount;		/*!< Core. */
  void		*freeptr;		/*!< Core. */
  WlzGMShell	*child;			/*!< A child shell of the model, others
  					     can be reached by walking the
					     next/prev fields. */
  int		vertexHTSz;		/*!< Vertex hash table size. */
  WlzGMVertex	**vertexHT;		/*!< Vertex hash table. */
  WlzGMModelR	res;			/*!< Model resources. */
} WlzGMModel;

#ifndef WLZ_EXT_BIND
/*!
* \typedef	WlzGMCbFn
* \ingroup	WlzGeoModel
* \brief	A pointer function to a function called when elements of a
*		Woolz geometric model are either created or deleted.
*/
typedef void	(*WlzGMCbFn)(WlzGMModel *, WlzGMElemP ,
			     WlzGMCbReason, void *);

/*!
* \struct	_WlzGMCbEntry
* \ingroup	WlzGeoModel
* \brief	
*/
typedef struct	_WlzGMCbEntry
{
  WlzGMCbFn	fn;
  void		*data;
  struct _WlzGMCbEntry *next;
} WlzGMCbEntry;

#endif /* WLZ_EXT_BIND */

/*!
* \struct	_WlzGMResIdx
* \ingroup	WlzGeoModel
* \brief	A resource index look up table (::WlzGMResIdxTb).
*		The array of indicies is a look up table from the indicies of a
*		GM to contigous indicies suitable for copying or outputing a
*		resource vector without holes.
*		Typedef: ::WlzGMResIdx.
*/
typedef struct _WlzGMResIdx
{
  int           idxCnt;         	/*!< Number of indicies in lut. */
  int           *idxLut;        	/*!< Index look up table. */
} WlzGMResIdx;

/*!
* \struct	_WlzGMResIdxTb
* \ingroup	WlzGeoModel
* \brief	Resource look up tables for all geometric elements in
*		a model.
*		Typedef: ::WlzGMResIdxTb.
*/
typedef struct _WlzGMResIdxTb
{
  WlzGMResIdx   vertex;
  WlzGMResIdx   vertexT;
  WlzGMResIdx   vertexG;
  WlzGMResIdx   diskT;
  WlzGMResIdx   edge;
  WlzGMResIdx   edgeT;
  WlzGMResIdx   face;
  WlzGMResIdx   loopT;
  WlzGMResIdx   shell;
  WlzGMResIdx   shellG;
} WlzGMResIdxTb;

/*!
* \struct      _WlzGMGridWSpCell3D
* \ingroup     WlzGeoModel
* \brief       A single cell entry in an axis aligned grid for a 3D model.
*              Typedef: ::WlzGMGridWSpCell3D
*/
typedef struct _WlzGMGridWSpCell3D
{
  WlzGMElemP    elem;                   /*! Geometric model element which
                                            intersects the cuboid of this
                                            cell. */
  struct _WlzGMGridWSpCell3D *next;     /*! The next cell in the linked list of
                                            cells. */
} WlzGMGridWSpCell3D;

/*!
* \struct      _WlzGMGridWSp3D
* \ingroup     WlzGeoModel
* \brief       An axis aligned grid of cuboid cells. This has an array (the
*              grid) of linked lists of cells, with the entries in each
*              list holding the faces of the 3D model which intersect the
*              cuboid of the cell.
*              Typedef: ::WlzGMGridWSp3D
*/
typedef struct _WlzGMGridWSp3D
{
  WlzGMModelType elemType;              /*! Element type that the grid is
                                            defined for. */
  WlzIVertex3   nCells;                 /*! Dimensions of the cell grid
                                            array in terms of the number of
                                            cells. */
  double        cellSz;                 /*! Each cell is an axis aligned
                                            square with this side length. */
  WlzDVertex3   org;                    /*! Model coordinate value at the cell
                                            origin. */
  struct _WlzGMGridWSpCell3D ****cells; /*! Array of linked list of cells. */
  int           cellVecMax;             /*! Number of cells allocated and the
                                            index into cell vector of next
                                            cell available. */
  AlcVector     *cellVec;               /*! Extensible vector of cells for
                                            adding to the array of linked lists
                                            of cells. */
} WlzGMGridWSp3D;

/************************************************************************
* Data structures for linear binary tree domains.
************************************************************************/

/*!
* \enum		_WlzLBTNodeClass2D
* \ingroup	WlzType
* \brief	Classification of a 2D LBT node a face of a 3D LBT node bycl
*		its connectivity with it's neighbouring nodes (2D) or faces
*		(3D)
*/
typedef enum _WlzLBTNodeClass2D
{
  WLZ_LBT_NODE_CLASS_2D_0, 	/*!< Node has no more than one neighbour on
                                     any edge.
\verbatim
                                  0
2D  +---+       3D               /
    |   |           z      5----/4
    |   |           ^  y   |   / |
    |   |           | %    |  1  |
    +---+           |/     |     |
                    +--->x 2-----3
\endverbatim */
  WLZ_LBT_NODE_CLASS_2D_1, 	/*!< Node has no more than one neighbour on all
  			  	     edges except one.
\verbatim
                                  0
2D  +-+-+       3D               /
    |   |           z      5--6-/4
    |   |           ^  y   |   / |
    |   |           | %    |  1  |
    +---+           |/     |     |
                    +--->x 2-----3
\endverbatim */
  WLZ_LBT_NODE_CLASS_2D_2,	/*!< Node only has more than one neighbour on
  				     adjacent edges.
\verbatim
                                  0
2D  +-+-+       3D               /
    |   |           z      5--6-/4
    |   +           ^  y   |   / |
    |   |           | %    |  1  7
    +---+           |/     |     |
                    +--->x 2-----3
\endverbatim */
  WLZ_LBT_NODE_CLASS_2D_3,	/*!< Node only has more than one neighbour on
  				     opposite edges.
\verbatim
                                  0
2D  +-+-+       3D               /
    |   |           z      5--6-/4
    |   |           ^  y   |   / |
    |   |           | %    |  1  |
    +-+-+           |/     |     |
                    +--->x 2--7--3
\endverbatim */

  WLZ_LBT_NODE_CLASS_2D_4,	/*!< Node has more than one neighbour on all
  				     edges except one.
\verbatim
                                  0
2D  +-+-+       3D               /
    |   |           z      5--6-/4
    |   +           ^  y   |   / |
    |   |           | %    |  1  8
    +-+-+           |/     |     |
                    +--->x 2--7--3
\endverbatim */
  WLZ_LBT_NODE_CLASS_2D_5	/*!< Node has more than one neighbour on all
  				     edges.
\verbatim
                                  0
2D  +-+-+       3D               /
    |   |           z      5--6-/4
    +   +           ^  y   |   / |
    |   |           | %    9  1  8
    +-+-+           |/     |     |
                    +--->x 2--7--3
\endverbatim */
} WlzLBTNodeClass2D;

/*!
* \def          WLZ_LBTDOMAIN_MAXDIGITS
* \ingroup      WlzType
* \brief        The maximum number of linear binary tree key digits,
*               which must be less than the number of bits in an int.
*/
#define WLZ_LBTDOMAIN_MAXDIGITS (30)

/*!
* \enum		_WlzLBTNodeFlags
* \ingroup	WlzGeoModel
* \brief	The reason a callback function is called.
*		Typedef: ::WlzGMCbReason.
*/
typedef enum _WlzLBTNodeFlags
{
  WLZ_LBT_NODE_FLAG_NONE	= 0,
  WLZ_LBT_NODE_FLAG_BOUNDARY	= 1     /*!< Node is adjacent to the objects
  					     boundary. */
} WlzLBTNodeFlags;

/*!
* \struct       _WlzLBTNode2D
* \ingroup      WlzType
* \brief        A 2D linear binary tree node for spatial domain representation.
*               Typedef: ::WlzLBTNode2D.
*/
typedef struct _WlzLBTNode2D
{
  unsigned		flags;		/*!< Bit flags fo the node. */
  unsigned              keys[3];        /*!< A single location key which
                                             uses bit interleaving with
                                             the bits interleaved between
                                             the elements for column, line
                                             and term in that order.
                                             Each of the array elements must
                                             have at least
                                             WLZ_LBTDOMAIN_MAXDIGITS bits. */
} WlzLBTNode2D;

/*!
* \struct       _WlzLBTNode3D
* \ingroup      WlzType
* \brief        A 3D linear binary tree node for spatial domain representation.
*               Typedef: ::WlzLBTNode3D.
*/
typedef struct _WlzLBTNode3D
{
  unsigned		flags;		/*!< Bit flags fo the node. */
  unsigned              keys[4];        /*!< A single location key which
                                             uses bit interleaving with
                                             the bits interleaved between
                                             the elements for column, line
                                             plane and term in that order.
                                             Each of the array elements must
                                             have at least
                                             WLZ_LBTDOMAIN_MAXDIGITS bits. */
} WlzLBTNode3D;

/*!
* \struct       _WlzLBTDomain2D
* \ingroup      WlzType
* \brief        A 2D linear binary tree spatial domain representation.
*               Typedef: ::WlzLBTDomain2D.
*/
typedef struct _WlzLBTDomain2D
{
  WlzObjectType         type;           /*!< From WlzCoreDomain. */
  int                   linkcount;      /*!< From WlzCoreDomain. */
  void                  *freeptr;       /*!< From WlzCoreDomain. */
  int                   line1;          /*!< First line coordinate. */
  int                   lastln;         /*!< Last line coordinate. */
  int                   kol1;           /*!< First column line
                                             coordinate. */
  int                   lastkl;         /*!< Last column  line
                                             coordinate. */
  int                   depth;          /*!< LBT depth. */
  int                   nNodes;         /*!< Number of nodes in the
                                             tree. */
  int                   maxNodes;       /*!< Number of nodes
                                             allocated. */
  WlzLBTNode2D          *nodes;         /*!< Array of nodes sorted by their
                                             location key. */
} WlzLBTDomain2D;

/*!
* \struct       _WlzLBTDomain3D
* \ingroup      WlzType
* \brief        A 3D linear binary tree spatial domain representation.
*               Typedef: ::WlzLBTDomain3D.
*/
typedef struct _WlzLBTDomain3D
{
  WlzObjectType         type;           /*!< From WlzCoreDomain. */
  int                   linkcount;      /*!< From WlzCoreDomain. */
  void                  *freeptr;       /*!< From WlzCoreDomain. */
  int               	plane1;     	/*!< First plane coordinate. */
  int               	lastpl; 	/*!< Last plane coordinate. */
  int                   line1;          /*!< First line coordinate. */
  int                   lastln;         /*!< Last line coordinate. */
  int                   kol1;           /*!< First column line
                                             coordinate. */
  int                   lastkl;         /*!< Last column  line
                                             coordinate. */
  int                   depth;          /*!< LBT depth. */
  int                   nNodes;         /*!< Number of nodes in the
                                             tree. */
  int                   maxNodes;       /*!< Number of nodes
                                             allocated. */
  WlzLBTNode3D          *nodes;         /*!< Array of nodes sorted by their
                                             location key. */
} WlzLBTDomain3D;

/************************************************************************
* Data structures for contours (both 2D and 3D).
************************************************************************/

/*!
* \enum		_WlzContourMethod
* \ingroup	WlzContour
* \brief	Contour generation methods.
		Typedef: ::WlzContourMethod.
*/
typedef enum _WlzContourMethod
{
  WLZ_CONTOUR_MTD_ISO,                  /*!< Iso-value. */
  WLZ_CONTOUR_MTD_GRD,            	/*!< Maximum gradient value. */
  WLZ_CONTOUR_MTD_BND,			/*!< Object boundary. */
  WLZ_CONTOUR_MTD_RBFBND		/*!< Object boundary established using
  					     a radial basis function. */
} WlzContourMethod;

/*!
* \struct	_WlzContour
* \ingroup	WlzContour
* \brief	A collection of 2D polylines or 3D surface elements
*		represented by a Woolz geometric model.
*		Typedef: ::WlzContour.
*/
typedef struct _WlzContour
{
  WlzObjectType type;                   /*!< WLZ_CONTOUR. */
  int		linkcount;		/*!< Core. */
  void		*freeptr;		/*!< Core. */
  WlzGMModel	*model;			/*!< The Woolz geometric model
  					     defining the contour. */
} WlzContour;


#ifndef WLZ_EXT_BIND
/*!
* \typedef      Wlz3DProjectionIntFn
* \ingroup      WlzFunction
* \brief        Callback function for the WlzGetProjectionFromObject().
*/
typedef WlzPixelV (*Wlz3DProjectionIntFn)(WlzPixelP, int, int, void *,
                                          WlzErrorNum *);
#endif

/*!
* \struct	_WlzLUTDomain
* \ingroup	WlzType
* \brief	A look up table domain.
* 		Typedef: ::WlzLUTDomain.
*/
typedef struct _WlzLUTDomain
{
  WlzObjectType type;		        /*!< WLZ_LUT. */
  int		linkcount;		/*!< Core. */
  void		*freeptr;		/*!< Core. */
  int		bin1;			/*!< Index of the first bin in the
                                             LUT. */
  int		lastbin;		/*!< Index of the last bin in the
                                             LUT. */
} WlzLUTDomain;

/*!
* \union	_WlzValues
* \ingroup	WlzType
* \brief	The union of Woolz values.
*		Typedef: ::WlzValues.
*/
typedef union _WlzValues
{
  struct _WlzCoreValues     *core;
  struct _WlzRagRValues     *v;
  struct _WlzRectValues     *r;
  struct _WlzIntervalValues *i;
  struct _WlzConvHullValues *c;
  struct _WlzVoxelValues    *vox;
  struct _WlzObject	    *obj;
  struct _WlzFeatValues     *fv;
  struct _WlzRectFeatValues *rfv;
  struct _WlzIndexedValues  *x;
  struct _WlzTiledValues    *t;
  struct _WlzLUTValues      *lut;
  struct _WlzPointValues    *pts;
} WlzValues;

/*!
* \union	_WlzDomain
* \ingroup	WlzType
* \brief	The union of Woolz domains.
*		Typedef: ::WlzDomain.
*/
typedef union _WlzDomain
{
  struct _WlzCoreDomain      *core;
  struct _WlzIntervalDomain  *i;
  struct _WlzPlaneDomain     *p;
  struct _WlzPolygonDomain   *poly;
  struct _WlzBoundList       *b;
  struct _WlzHistogramDomain *hist;
  struct _WlzRect	     *r;
  struct _WlzFRect           *fr;
  struct _WlzAffineTransform *t;
  struct _WlzWarpTrans       *wt;
  struct _WlzContour 	     *ctr;
  struct _WlzMeshTransform   *mt;
  struct _WlzLBTDomain2D     *l2;
  struct _WlzLBTDomain3D     *l3;
  struct _WlzCMesh2D	     *cm2;
  struct _WlzCMesh2D5	     *cm2d5;
  struct _WlzCMesh3D	     *cm3;
  struct _WlzPoints	     *pts;
  struct _WlzLUTDomain       *lut;
  struct _WlzConvHullDomain2 *cvh2;
  struct _WlzConvHullDomain3 *cvh3;
  struct _WlzThreeDViewStruct *vs3d;
  struct _WlzBSpline          *bs;
} WlzDomain;

/*!
* \struct	_WlzPropertyList
* \ingroup	WlzProperty
* \brief	A property list which has a type, link count and a
*		linked list of properties.
*/
typedef struct	_WlzPropertyList
{
  WlzObjectType	type;
  int		linkcount;
  AlcDLPList	*list;
} WlzPropertyList;

/*!
* \struct	_WlzCoreProperty
* \ingroup	WlzProperty
* \brief	Core property with sufficient to data to provide the type
*		and enough to allow the property to be freed.
*		Typedef: ::WlzCoreProperty.
*/
typedef struct _WlzCoreProperty
{
  WlzObjectType		type;			/*!< Type */
  int 			linkcount;		/*!< linkcount */
  void 			*freeptr;		/*!< free pointer */
} WlzCoreProperty;

/*!
* \struct	_WlzSimpleProperty
* \ingroup	WlzProperty
* \brief	A simple property to hold arbitrary length string data.
*		Read and writing then coercing to a structure with
*		numerical values will not be portable.
*		Typedef: ::WlzSimpleProperty.
*/
typedef struct _WlzSimpleProperty
{
  WlzObjectType		type;			/*!< Type */
  int 			linkcount;		/*!< linkcount */
  void 			*freeptr;		/*!< free pointer */
  unsigned long		size;			/*!< Data size of the property */
  void 			*prop;			/*!< Pointer to the property */
} WlzSimpleProperty;

/*!
* \def		EMAP_PROPERTY_MODELNAME_LENGTH
* \ingroup	WlzProperty
* \brief	Maximum length of the model name in an EMAP property.
*/
#define EMAP_PROPERTY_MODELNAME_LENGTH	32

/*!
* \def		EMAP_PROPERTY_UID_LENGTH
* \ingroup	WlzProperty
* \brief	Maximum length of the model or anatomy UID in an EMAP property.
*/
#define EMAP_PROPERTY_UID_LENGTH	16

/*!
* \def		EMAP_PROPERTY_VERSION_LENGTH
* \ingroup	WlzProperty
* \brief	Maximum length of the version string in an EMAP property.
*/
#define EMAP_PROPERTY_VERSION_LENGTH	16

/*!
* \def		EMAP_PROPERTY_AUTHORNAME_LENGTH
* \ingroup	WlzProperty
* \brief	Maximum length of the author strings in an EMAP property.
*/
#define	EMAP_PROPERTY_AUTHORNAME_LENGTH 64

/*!
* \def		EMAP_PROPERTY_MACHINENAME_LENGTH
* \ingroup	WlzProperty
* \brief	Maximum length of the machine name strings in an EMAP property.
*/
#define	EMAP_PROPERTY_MACHINENAME_LENGTH 64

/*!
* \def		EMAP_PROPERTY_STAGE_LENGTH
* \ingroup	WlzProperty
* \brief	Maximum length of the stage strings in an EMAP property.
*/
#define	EMAP_PROPERTY_STAGE_LENGTH 32

/*!
* \struct	_WlzEMAPProperty
* \ingroup	WlzProperty
* \brief	A property to hold EMAP information to attach to
*		the reference models, anatomy and GE domains. MAPaint
*		and atlas tools will propogate the information as 
*		required.
*		Typedef: ::WlzEMAPProperty.
*/
typedef struct _WlzEMAPProperty
{
  WlzObjectType		type;				/*!< Type */
  int 			linkcount;			/*!< linkcount */
  void 			*freeptr;			/*!< free pointer */
  WlzEMAPPropertyType	emapType;			/*!< EMAP property type */
  char			modelUID[EMAP_PROPERTY_UID_LENGTH]; /*!< model UID */
  char			anatomyUID[EMAP_PROPERTY_UID_LENGTH];/*!< anatomy UID */
  char			targetUID[EMAP_PROPERTY_UID_LENGTH];/*!< target model UID */
  char			targetVersion[EMAP_PROPERTY_VERSION_LENGTH]; /*!< target model version */
  char			stage[EMAP_PROPERTY_STAGE_LENGTH];/*!< Embryo  stage */
  char			subStage[EMAP_PROPERTY_STAGE_LENGTH];/*!< Embryo  sub-stage */
  char			modelName[EMAP_PROPERTY_MODELNAME_LENGTH];/*!< Volume model name */
  char			version[EMAP_PROPERTY_VERSION_LENGTH];/*!< Model version */
  char			*fileName;			/*!< Original filename (not very useful) */
  long			creationTime;			/*!< Original creation time */
  char			creationAuthor[EMAP_PROPERTY_AUTHORNAME_LENGTH];/*!< Creation author */
  char			creationMachineName[EMAP_PROPERTY_MACHINENAME_LENGTH];/*!< Original creation machine name */
  long			modificationTime;		/*!< Modification time */
  char			modificationAuthor[EMAP_PROPERTY_AUTHORNAME_LENGTH];/*!< Modification author */
  char			*comment;			/*!< Text comment
  							     string. */
} WlzEMAPProperty;

/*!
* \struct	_WlzNameProperty
* \ingroup	WlzProperty
* \brief	A simple null terminated ASCII string for the object's
* 		name.
*		Typedef: ::WlzNameProperty.
*/
typedef struct _WlzNameProperty
{
  WlzObjectType		type;			/*!< Type. */
  int			linkcount;		/*!< linkcount. */
  void			*freeptr;		/*!< Free pointer. */
  char			*name;			/*!< A simple ASCII name
						     string. */
} WlzNameProperty;

/*!
* \struct	_WlzGreyProperty
* \ingroup	WlzProperty
* \brief	A single grey value, which for example might represent
*		the preferred display colour of a binary domain.
*		Typedef: ::WlzGreyProperty.
*/
typedef struct _WlzGreyProperty
{
  WlzObjectType		type;			/*!< Type. */
  int			linkcount;		/*!< linkcount. */
  void			*freeptr;		/*!< Free pointer. */
  char			*name;			/*!< An associated name
  						     string which conveys
						     the meaning of the 
						     value. May be NULL */
  struct _WlzPixelV	value;			/*!< The pixel value which
  						     both encodes the grey type
						     and it's value. */
} WlzGreyProperty;

/*!
* \struct	_WlzTextProperty
* \ingroup	WlzProperty
* \brief	A pair of simple null terminated ASCII strings one for the
*               property name and one for it's value.
*		Typedef: ::WlzTextProperty.
*/
typedef struct _WlzTextProperty
{
  WlzObjectType		type;			/*!< Type. */
  int			linkcount;		/*!< linkcount. */
  void			*freeptr;		/*!< Free pointer. */
  char			*name;			/*!< Name string. */
  char			*text;			/*!< Text string. */
} WlzTextProperty;

/*!
* \union	_WlzProperty
* \ingroup	WlzProperty
* \brief	A union of pointers for properties.
*		Typedef: WlzProperty.
*/
typedef union _WlzProperty
{
  struct _WlzCoreProperty	*core;
  struct _WlzSimpleProperty	*simple;
  struct _WlzEMAPProperty	*emap;
  struct _WlzNameProperty	*name;
  struct _WlzGreyProperty	*greyV;
  struct _WlzTextProperty	*text;
} WlzProperty;

/************************************************************************
* The Woolz objects.
************************************************************************/
/*!
* \struct	_WlzCoreObject
* \ingroup	WlzType
* \brief	The core Woolz object type which can be used to determine
*		the type of a Woolz object. 
*		Typedef: ::WlzCoreObject.
*/
typedef struct _WlzCoreObject
{
  WlzObjectType	type;			/*!< The Woolz object type. */
  int		linkcount;		/*!< The link count: A counter
  					     for the number of references to
					     the object, which should only be
					     accessed through WlzUnlink(),
					     WlzAssignObject() and
					     WlzFreeObj(). */
} WlzCoreObject;

/*!
* \struct	_WlzObject
* \ingroup	WlzType
* \brief	The Woolz object.
*		Typedef: ::WlzObject.
*/
typedef struct _WlzObject
{
  WlzObjectType      type;		/*!< From WlzCoreObject. */
  int                linkcount;		/*!< From WlzCoreObject. */
  WlzDomain          domain;		/*!< The objects domain: It's
					     spatial extent or
  					     geometric properties. */
  WlzValues          values;		/*! The values defined within
  					    the object's domain. */
  WlzPropertyList    *plist;		/*! A list of the object's
  					    properties. */
  struct _WlzObject  *assoc;		/*! An object which is
  					    assosciated with this object. */
} WlzObject;

/*!
* \struct	_WlzCompoundArray
* \ingroup	WlzType
* \brief	A compound object implemented as either an array or
*		a linked list of other objects. There is a distinction between
*		an compound of the same type (e.g. resulting from a labelling)
*		and a compound of different types (e.g. resulting from a range
*		of image processes from a single original object).
*		Typedef: ::WlzCompoundArray.
*/
typedef struct _WlzCompoundArray
{
  WlzObjectType type;			/*!< From WlzCoreObject. */
  int           linkcount;		/*!< From WlzCoreObject. */
  WlzObjectType otype;		       	/*!< The permitted type if
  					    constrained. */
  int           n;		 	/*!< The  number of objects */
  WlzObject     **o;			/*!< The list of Woolz object
  					     pointers. */
  WlzPropertyList *plist;		/*! A list of the object's
  					    properties.  */
  WlzObject     *assoc;
} WlzCompoundArray;

/************************************************************************
* Domains.
************************************************************************/
/*!
* \struct	_WlzCoreDomain
* \ingroup	WlzType
* \brief	The core domain: All Woolz domains have all the fields
*		of the core domain in the same order and before any
*		others, so allowing a domain to be assigned, freed
*		and have it's type established.
		Typedef: ::WlzCoreDomain.
*/
typedef struct _WlzCoreDomain
{
  WlzObjectType   type;				/*!< The type of domain. */
  int             linkcount;			/*!< The link count: A counter
  						     for the number of
						     references to the domain,
						     which should only be
						     accessed through
						     WlzUnlink(),
						     WlzAssignDomain() and
						     WlzFreeDomain(). */
  void            *freeptr;			/*!< A stack with pointers
  						     that can be freed by
						     AlcFreeStackFree(). */
} WlzCoreDomain;

/*!
* \struct	_WlzIntervalDomain
* \ingroup	WlzType
* \brief	A 2D domain defining an arbitrary  region of space in 2D.
*		The domain may be of type WLZ_INTERVALDOMAIN_INTVL or
*		WLZ_INTERVALDOMAIN_RECT. If the domain is of type
*		WLZ_INTERVALDOMAIN_RECT then the intvlines field is
*		not used. For WLZ_INTERVALDOMAIN_INTVL domains the
*		intervals in a line must be contiguous.
*		Typedef: ::WlzIntervalDomain.
*/
typedef struct _WlzIntervalDomain
{
  WlzObjectType     type;			/*!< From WlzCoreDomain. */
  int               linkcount;			/*!< From WlzCoreDomain. */
  void              *freeptr;			/*!< From WlzCoreDomain. */
  int               line1;			/*!< First line coordinate. */
  int               lastln;			/*!< Last line coordinate. */
  int               kol1;			/*!< First column
  						     coordinate. */
  int               lastkl;			/*!< Last column coordinate. */
  struct _WlzIntervalLine *intvlines;   	/*!< Array of interval line
  						     structures. */
} WlzIntervalDomain;

/*!
* \struct	_WlzPlaneDomain
* \ingroup	WlzType
* \brief	A 3D domain defining an arbitrary  region of space in 3D.
*		The 3D plane domain composed of plane-wise array of 2D domains.
*		Typedef: ::WlzPlaneDomain.
*/
typedef struct _WlzPlaneDomain
{
  WlzObjectType     type;			/*!< From WlzCoreDomain. */
  int               linkcount;			/*!< From WlzCoreDomain. */
  void              *freeptr;			/*!< From WlzCoreDomain. */
  int               plane1;     		/*!< First plane coordinate. */
  int               lastpl;     		/*!< Last plane coordinate. */
  int               line1;      		/*!< First line coordinate. */
  int               lastln;			/*!< Last line coordinate. */
  int               kol1;       	    	/*!< First column line
  						     coordinate. */
  int               lastkl;     	     	/*!< Last column  line
  						     coordinate. */
  WlzDomain 	    *domains;       		/*!< Array of pointers to
  						     2D domains. */
  float 	    voxel_size[3];      	/*!< Array of nominal voxel
  						     dimensions. */
} WlzPlaneDomain;

/************************************************************************
* Intervals.						
************************************************************************/
/*!
* \struct	_WlzIntervalLine
* \ingroup	WlzType
* \brief	A line of intervals.
*		Typedef: ::WlzIntervalLine.
*/
typedef struct _WlzIntervalLine
{
  int                  nintvs;		/*!< Number of intervals on line. */
  struct _WlzInterval *intvs;		/*!< Array of intervals. */
} WlzIntervalLine;

/*!
* \struct	_WlzInterval
* \ingroup	WlzType
* \brief	A single interval.
*		Typedef: ::WlzInterval.
*/
typedef struct _WlzInterval
{
  int ileft;				/*!< Left most pixel of interval. */
  int iright;				/*!< Right most pixel of interval. */
} WlzInterval;

/*!
* \struct	_WlzDynItvPool
* \ingroup	WlzType
* \brief	Dynamic interval pool, for building interval domains.
* 		Typedef: ::WlzDynItvPool.
*/
typedef struct _WlzDynItvPool
{
  WlzInterval	*itvBlock;		/*!< Array of intervals. */
  int		itvsInBlock;	        /*!< Number of intervals in array. */
  int		offset;			/*!< Offset into array for next
  					     available interval. */
} WlzDynItvPool;

/*!
* \struct       _WlzPartialItv2D
* \ingroup      DomainOps
* \brief        Data structure that can be used to hold partial intervals.
*               These can then be sorted and condensed to find the intervals
*               for an interval domain.
*/
typedef struct _WlzPartialItv2D
{
  int                   ileft;
  int                   iright;
  int                   ln;
} WlzPartialItv2D;

/*!
* \struct       _WlzPartialItv3D
* \ingroup      DomainOps
* \brief        Data structure that can be used to hold partial intervals.
*               These can then be sorted and condensed to find the intervals
*               for a plane domain.
*/
typedef struct _WlzPartialItv3D
{
  int                   ileft;
  int                   iright;
  int                   ln;
  int			pl;
} WlzPartialItv3D;

/************************************************************************
* Value tables.
************************************************************************/

/*!
* \struct	_WlzCoreValues
* \ingroup	WlzType
* \brief	All Woolz value tables must have all the fields of the
*		core values, in the same order and before any others,
*		so allowing a values to be assigned, freed and have it's
*		type established.
*		Typedef: ::WlzCoreValues.
*/
typedef struct _WlzCoreValues
{
  WlzObjectType type;	
  int       linkcount;	
} WlzCoreValues;

/*!
* \struct	_WlzValueLine
* \ingroup	WlzType
* \brief	Grey values along a line.
*		Typedef: ::WlzValueLine.
*/
typedef struct _WlzValueLine
{
  int      vkol1;			/*!< Relative left end. */
  int      vlastkl;			/*!< Relative right end. */
  WlzGreyP values;		  	/*!< Array of values. */
} WlzValueLine;

/*!
* \struct	_WlzTiledValueBuffer
* \ingroup	WlzType
* \brief	Position of and data for locating and buffering any interval
* 		of values in either 2 or 3D tiled value table.
*/
typedef struct _WlzTiledValueBuffer
{
  int		pl;			/*!< Plane of interval (relative to 
  					     tiled value table. */
  int		ln;			/*!< Line of interval (relative to 
  					     tiled value table. */
  int		kl[2];			/*!< Left most then right most column
  					     of interval (relative to tiled
					     value table. */
  size_t	lo;			/*!< Partial line offset within a
  					     tile. */
  size_t	li;			/*!< Partial line index to a tile. */
  int		valid;			/*!< Non-zero if the line buffer has
  					     valid values. */
  int		mode;			/*!< Valid access modes for the
  					     tiled values. */
  WlzGreyType	gtype;			/*!< Grey type buffer allocaed for. */
  WlzGreyP	lnbuf;			/*!< Buffer large enough to hold any
  					     single line of values padded to
					     an integral number of tile widths.
					     The values are all relative to the
					     value table origin. */
} WlzTiledValueBuffer;

/*!
* \struct	_WlzRagRValues
* \ingroup	WlzType
* \brief	The ragged rectangle values table.
*		The type encodes both the type of value table and the type of
*		grey value.
*		Typedef: ::WlzRagRValues.
*/
typedef struct _WlzRagRValues
{
  WlzObjectType type;			/*!< From WlzCoreValues. */
  int       linkcount;			/*!< From WlzCoreValues. */
  void      *freeptr;			/*!< From WlzCoreValues. */
  WlzValues original_table; 		/*!< If non-NULL, the values table
  					     which owns the raw values we
					     are using. */
  int       line1;			/*!< First line. */
  int       lastln;			/*!< Last line. */
  int       kol1;			/*!< First column. */
  int       width;			/*!< Width. */
  WlzPixelV bckgrnd;			/*!< Background value for pixels not
  					     in object. */
  WlzValueLine *vtblines;	      	/*!< Array of value table line
  					     structures. */
} WlzRagRValues;

/*!
* \struct	_WlzRectValues
* \ingroup	WlzType
* \brief	The rectangle values table.
*		The type encodes both the type of value table and the type of
*		grey value.
*		Typedef: ::WlzRectValues.
*/
typedef struct _WlzRectValues
{
  WlzObjectType type;			/*!< From WlzCoreValues. */
  int 		linkcount;		/*!< From WlzCoreValues. */
  void 		*freeptr;		/*!< From WlzCoreValues. */
  WlzValues 	original_table; 	/*!< If non-NULL, the values table
  					     which owns the raw values we
					     are using. */
  int 		line1;			/*!< First line. */
  int 		lastln;			/*!< Last line. */
  int 		kol1;			/*!< First column. */
  int 		width;			/*!< Width. */
  WlzPixelV 	bckgrnd; 		/*!< Background value for points
   					     not in object. */
  WlzGreyP 	values;			/*!< Contiguous array of values. */
} WlzRectValues;

/*!
* \struct	_WlzValueIntervalLine
* \ingroup	WlzType
* \brief	One line's worth of grey value intervals.
*		Typedef: ::WlzValueIntervalLine.
*/
typedef struct _WlzValueIntervalLine
{
  int           nintvs;			/*!< Number of grey value intervals. */
  WlzValueLine *vtbint;			/*!< Pointer to grey value
  					     intervals. */
} WlzValueIntervalLine;

/*!
* \struct	_WlzIntervalValues
* \ingroup	WlzType
* \brief	An interval structured value table.
*		The type encodes both the type of value table and the type of
*		grey value.
*               Typedef: ::WlzIntervalValues.
*/
typedef struct _WlzIntervalValues
{
  WlzObjectType type;			/*!< From WlzCoreValues. */
  int      	linkcount;      	/*!< From WlzCoreValues. */
  void    	*freeptr;		/*!< From WlzCoreValues. */
  WlzValues 	original_table;		/*!< If non-NULL, the values table
   					     which owns the raw values we
  					     are using. */
  int   	line1;	      		/*!< First line. */
  int    	lastln;	      		/*!< Last line. */
  int     	kol1;	      		/*!< First column. */
  int      	width;	      		/*!< Width. */
  WlzPixelV 	bckgrnd;        	/*!< Background value for points
  					     not in object. */
  WlzValueIntervalLine *vil;   		/*!< Pointers to structures of grey
  					     table lines. */
} WlzIntervalValues;

/*!
* \struct	_WlzVoxelValues
* \ingroup	WlzType
* \brief	Voxel value table.
*		Typedef: ::WlzVoxelValues.
*/
typedef struct _WlzVoxelValues
{
  WlzObjectType	type;			/*!< From WlzCoreValues. */
  int           linkcount;		/*!< From WlzCoreValues. */
  void          *freeptr;		/*!< From WlzCoreValues. */
  WlzValues 	original_table;		/*!< If non-NULL, the values table
  					     which owns the raw values we
					     are using. */
  int           plane1;			/*!< First plane. */
  int           lastpl;			/*!< Last plane. */
  WlzPixelV     bckgrnd;		/*!< Background value for points
  					     not in object. */
  WlzValues     *values;   		/*!< Array of pointers to value
  					     tables. */
} WlzVoxelValues;

/*!
* \struct       _WlzValueAttach
* \ingroup      WlzType
* \brief        Specifies what values (for example thoose in an indexed
                value table) are attached to.
*               Enum: ::WlzValueAttach
*/
typedef enum _WlzValueAttach
{
  WLZ_VALUE_ATTACH_NONE,                /*!< No attachment. */
  WLZ_VALUE_ATTACH_NOD,                 /*!< Attached to mesh nodes. */
  WLZ_VALUE_ATTACH_ELM                 /*!< Attached to mesh elements. */
} WlzValueAttach;

/*!
* \struct       _WlzIndexedValues
* \ingroup      WlzType
* \brief        In indexed value table.
*               Typedef: ::WlzIndexedValues.
*/
typedef struct _WlzIndexedValues
{
  WlzObjectType type;                   /*!< From WlzCoreValues:
                                             WLZ_INDEXED_VALUES. */
  int           linkcount;              /*!< From WlzCoreValues. */
  void          *freeptr;		/*!< From WlzCoreValues, although
  					     this won't free the values
					     themselves. */
  int           rank;                   /*!< The rank of the individual values.
                                             Here the rank for a scalar is 0,
                                             for a 1D array it is 1 and for
                                             and for individual values that
                                             are nD arrays the rank is n. */
  int           *dim;                   /*!< The dimensions of individual
                                             indexed values. The dimensions
                                             are a 1D array with the number
                                             of entries equal to the rank.
					     A dimension array is only
					     allocated if the rank > 0,
					     for rank == 0 dim is NULL. */
  WlzGreyType   vType;                  /*!< The type of the data in the
                                             individual values. */
  WlzValueAttach attach;                /*!< Specifies what the values are
                                             attached to. */
  AlcVector     *values;                /*!< The indexed values. */
} WlzIndexedValues;

/*!
* \struct       _WlzTiledValues
* \ingroup      WlzType
* \brief        A  tiled value table for both 2 an 3D domain objects.
*               Typedef: ::WlzTiledValues.
*
* 		Individual pixel/voxel values may contain a single
* 		(scalar) value of a multidimensional array of values.
* 		When an array of values is used the values are held
* 		contiguously.
* 		
* 		The grey values are stored in square or cubic tiles
* 		with tileSz being the number of values in each tile
* 		irrespective of the grey value type. An index to the
* 		tiles is stored as a simple one dimensional array.
*
* 		To access a grey value at some position \f$(x,y,z)\f$
* 		which is known to be within the tiled value table:
* 		First the position relative to the first column (\f$x_0\f$),
* 		line (\f$y_0\f$) and plane (\f$z_0\f$)of the value table
* 		is computed:
*		\f$x_r = x - x_0, y_r = y - y_0, z_r = z - z_0\f$.
*		The tile index (\f$i\f$) and within tile offset (\f$o\f$)
*		are then computed:
*		\f{eqnarray*}
		i_x = x_r / w \\
		i_y = y_r / w \\
		i_z = z_r / w \\
		o_x = x_r - (w i_x) \\
		o_y = y_r - (w i_y) \\
		o_z = z_r - (w i_z) \\
		i = (((i_z ni_y) + i_y) ni_x) + i_x \\
		o = (((o_z nt_y) + o_y) nt_x) + o_x \\
		g = T \left[ ( I \left[ i \right] t_s) + o \right]
		\f}
*		Where \f$I\f$ is the index table and \f$T\f$ the start
*		of the tile data.
*
* 		The tile data may be memory mapped instead of read into
* 		memory in which case the file descriptor will have a
* 		non-negative value. This can be used to close the file.
*
* 		A memory mapped tiled values object can only have it's
* 		grey values changed if the file was opened for writing
* 		attempting to change the grey values of an object only
* 		opened for read will lead to a memory fault. The read/write
* 		status is respected by the interval scanning access function
* 		WlzNextGreyInterval(), but code which uses random access
* 		via WlzGreyValueGet() or sequential access via WlzIterate()
* 		should check this. Attempting to write to the memory mapped
* 		values of an object only opened for reading will give a
* 		segmentation fault. A tiled values object can be written
* 		to only if the file descriptor is invalid (< 0) or if
* 		the file was opened in write or append mode. The function
* 		WlzTiledValuesMode() may also be used to determine the
* 		appropriate access mode(s) for the values table.
*/
typedef struct _WlzTiledValues
{
  WlzObjectType type;                   /*!< From WlzCoreValues:
                                             built from WLZ_GREY_TAB_TILED
					     and the grey value type. */
  int           linkcount;              /*!< From WlzCoreValues. */
  void          *freeptr;		/*!< From WlzCoreValues, although
  					     this won't free the values
					     themselves. */
  WlzValues 	original_table;		/*!< If non-NULL, the values table
  					     which owns the raw values we
					     are using. */
  int		dim;			/*!< The dimension of the value table,
  					     ie 2 or 3D. */
  int		kol1;			/*!< First column. */
  int		lastkl;			/*!< Last column. */
  int		line1;			/*!< First line. */
  int		lastln;			/*!< Last line. */
  int		plane1;			/*!< First plane. */
  int		lastpl;			/*!< Last plane. */
  WlzPixelV     bckgrnd;		/*!< Background value. */
  unsigned int  vRank;                  /*!< The rank of the individual values.
                                             Here the rank for a scalar is 0,
                                             for a 1D array it is 1 and for
                                             and for individual values that
                                             are nD arrays the rank is n.
					     The rank will only ever be > 0
					     when the grey table type is
					     composed using rank > 0. */
  unsigned int  *vDim;                  /*!< The dimensions of individual
                                             values. The dimensions are a
					     1D array with the number of
					     entries equal to the rank.
					     A dimension array is only
					     allocated if the rank > 0,
					     for rank == 0 dim is NULL. */
  unsigned int	vpe;			/*!< Values per element which has the
  					     value
			       \f$ \prod_i^\textrm{vRank}{\textrm{vDim}[i]} \f$
					     as compulted by
					     WlzTiledValuesValPerElm()
					     (value 1 when vRank = 0) 
					     is included since it is frequently
					     used. */
  size_t	tileSz;		        /*!< The number of elements in each
  					     tile which may be less than
					     the number of bytes. */
  size_t        tileWidth;              /*!< Width of the tiles (\f$w\f$). */
  size_t	numTiles;		/*!< The total number of tiles. */
  int		*nIdx;			/*!< Number of index columns,
  					     lines, .... */
  unsigned int	*indices;		/*!< Table of tile indices. */
  int		fd;			/*!< File descriptor if tiles are
  					     memory mapped else -1. */
  long		tileOffset;             /*!< Offset from the start of the
  					     file to the tiles. This may be
					     set even if not memory mapped. */
  WlzGreyP 	tiles;			/*!< The tiles. */
} WlzTiledValues;

/*!
* \def		WLZ_TILEDVALUES_TILE_SIZE
* \ingroup	WlzType
* \brief	The default number of pixels/voxels in a tiled value
* 		table tile. Chosen so that tiles will occupy an integer
* 		number of disk block,
*/
#define WLZ_TILEDVALUES_TILE_SIZE	(4096)

/*!
* \struct	_WlzLUTValues
* \ingroup	WlzType
* \brief	Look up table values.
* 		Typedef: ::WlzLUTValues
*/
typedef struct _WlzLUTValues
{
  WlzObjectType	type;			/*!< WLZ_LUT. */
  int		linkcount;		/*!< From WlzCoreValues. */
  void		*freeptr;		/*!< From WlzCoreValues. */
  WlzGreyType	vType;			/*!< Type for the LUT values. */
  int		maxVal;			/*!< Number of values allocated. */
  WlzGreyP	val;			/*!< LUT values. */
} WlzLUTValues;

/*!
* \struct	_WlzPointValues
* \ingroup	WlzType
* \brief	Point values - values with arbitrary rank and dimension
* 		defined at points.
* 		Typedef: ::WlzPointValues
*/
typedef struct _WlzPointValues
{
  WlzObjectType	type;			/*!< WLZ_POINT_VALUES. */
  int		linkcount;		/*!< From WlzCoreValues. */
  void		*freeptr;		/*!< From WlzCoreValues. */
  int           rank;                   /*!< The rank of the individual values.
                                             Here the rank for a scalar is 0,
                                             for a 1D array it is 1 and for
                                             and for individual values that
                                             are nD arrays the rank is n. */
  int           *dim;                   /*!< The dimensions of individual
                                             indexed values. The dimensions
                                             are a 1D array with the number
                                             of entries equal to the rank.
					     A dimension array is only
					     allocated if the rank > 0,
					     for rank == 0 dim is NULL. */
  int		pSz;			/*!< Size of each point this is
                                             the product of the dimensions
					     times the value type size and
					     is frequently used for accessing
					     the values. */
  WlzGreyType   vType;                  /*!< The type of the data in the
                                             individual values. */
  size_t	maxPoints;		/*!< Number of points allocated. */
  WlzGreyP      values;                 /*!< The indexed values. */
} WlzPointValues;

/************************************************************************
* Point domains.						
************************************************************************/
/*!
* \struct	_WlzPoints
* \ingroup	WlzFeatures
* \brief	An array of either 2D or 3D points which may have
*       	either integral of floating point values. Possible
*       	types are: WLZ_POINTS_2I, WLZ_POINTS_2D, WLZ_POINTS_3I
*       	and WLZ_POINTS_3D.
*       	Typedef: ::WlzPoints
*/
typedef struct _WlzPoints
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int		linkcount;		/*!< From WlzCoreDomain. */
  void		*freeptr;		/*!< From WlzCoreDomain. */
  int		nPoints;		/*!< Number of points. */
  int		maxPoints;		/*!< The maximum number of points
					     for which space has been
					     allocated. */
  WlzVertexP	points;			/*!< Array of point vertices. */
} WlzPoints;


/************************************************************************
* Spline domains.
************************************************************************/
/*!
* \def		WLZ_BSPLINE_ORDER_MIN
* \ingroup	WlzType
* \brief	The minimum order of the B-spline in a WlzBSpline domain.
*/
#define WLZ_BSPLINE_ORDER_MIN	(1)

/*!
* \def		WLZ_BSPLINE_ORDER_MAX
* \ingroup	WlzType
* \brief	The maximum order of the B-spline in a WlzBSpline domain.
*/
#define WLZ_BSPLINE_ORDER_MAX	(5)

/*!
* \struct	_WlzBSpline
* \ingroup	WlzFeatures
* \brief	Spline based line curves in either 2 or 3D.
* 		Possible types are: WLZ_BSPLINE_C2D or WLZ_BSPLINE_C3D.
* 		Typedef: ::WlzBSpline
*/
typedef struct _WlzBSpline
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int		linkcount;		/*!< From WlzCoreDomain. */
  void		*freeptr;		/*!< From WlzCoreDomain. */
  int		order;			/*!< Order of the B-spline. */
  int		nKnots;			/*!< Number of knots (also number
  					     of coefficients per dimension). */
  int		maxKnots;		/*!< Knots space allocated. */
  double	*knots;			/*!< Array of knots (parametric). */
  double	*coefficients;		/*!< Array of B-spline coefficients. */
} WlzBSpline;

/************************************************************************
* Polygon domains.						
************************************************************************/
/*!
* \struct	_WlzPolygonDomain
* \ingroup	WlzPolyline
* \brief	A 2D polyline domain with possible types: WLZ_POLYGON_INT, 
*		WLZ_POLYGON_FLOAT  or WLZ_POLYGON_DOUBLE. 
*		Typedef: ::WlzPolygonDomain.
*/
typedef struct _WlzPolygonDomain
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int linkcount;			/*!< From WlzCoreDomain. */
  void *freeptr;			/*!< From WlzCoreDomain. */
  int nvertices;			/*!< Number of vertices. */
  int maxvertices;			/*!< The maximum number of vertices
  					     for which space has been
					     allocated. */
  WlzIVertex2 *vtx; 			/*!< Array of vertices.
  					     This may need casting according
					     to the type field. */
} WlzPolygonDomain;

/*!
* \struct	_WlzPolygonDomain3
* \ingroup	WlzPolyline
* \brief	A 2D polyline domain with possible types:WLZ_POLYGON_INT, 
*		WLZ_POLYGON_FLOAT  or WLZ_POLYGON_DOUBLE. 
*		Typedef: ::WlzPolygonDomain.
*/
typedef struct _WlzPolygonDomain3
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int 	linkcount;			/*!< From WlzCoreDomain. */
  void *freeptr;			/*!< From WlzCoreDomain. */
  int nvertices;			/*!< Number of vertices. */
  int maxvertices;	   		/*!< The maximum number of vertices
  					     for which space has been
					     allocated. */
  WlzIVertex2 *vtx; 			/*!< Array of vertices.
  					     This may need casting according
					     to the type field. */
} WlzPolygonDomain3;

/************************************************************************
* Boundary list.						
************************************************************************/
/*!
* \struct	_WlzBoundList
* \ingroup	WlzBoundary
* \brief	A complete list of a set of boundaries which is encoded
*		in tree form.
*/
typedef struct _WlzBoundList
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int linkcount;			/*!< From WlzCoreDomain. */
  void *freeptr;			/*!< From WlzCoreDomain. */
  struct _WlzBoundList *up;		/*!< The containing hole or piece,
				   	     NULL if the universal hole
					     (very top). */
  struct _WlzBoundList *next;		/*!< Next hole or piece at same level
  					     and lying within same piece or
					     hole, NULL if no more at this
					     level. */
  struct _WlzBoundList *down;		/*!< First enclosed structure, NULL if
  				             none. */
  int wrap;				/*!< Wrap number: The number of points
  					     of boundary included both at
					     start and end of polygon
					     representation. */
  WlzPolygonDomain *poly;	       /*!< The polygon representation of
  					    this boundary. */
} WlzBoundList;

/*!
* \struct       _WlzConvHullDomain2
* \ingroup	WlzConvexHull
* \brief        A 2D convex hull with counter clockwise ordered vertices
*               and segments implicitly defined by the polygon of the
*               ordered vertices.
*               Typedef: ::WlzConvHullDomain3
*/
typedef struct _WlzConvHullDomain2
{
  WlzObjectType         type;           /*!< From WlzCoreDomain
					     (WLZ_CONVHULL_DOMAIN_2D). */
  int                   linkcount;      /*!< From WlzCoreDomain. */
  void                  *freeptr;       /*!< From WlzCoreDomain. */
  int                   nVertices;      /*!< Number of vertices. */
  int                   maxVertices;    /*!< The maximum number of vertices
                                             for which space has been
                                             allocated. */
  WlzVertexType         vtxType;        /*!< Vertex type, either WLZ_VERTEX_I2
                                             or WLZ_VERTEX_D2. */
  WlzVertex             centroid;       /*!< Centroid of the convex hull. */
  WlzVertexP            vertices;       /*!< Array of vertices. */
} WlzConvHullDomain2;

/*!
* \struct       _WlzConvHullDomain3
* \ingroup	WlzConvexHull
* \brief        A 3D convex hull with coordinate vertices and faces defined
*               by vertex index triples.
*               Typedef: ::WlzConvHullDomain3
*/
typedef struct _WlzConvHullDomain3
{
  WlzObjectType         type;           /*!< From WlzCoreDomain
					     (WLZ_CONVHULL_DOMAIN_3D). */
  int                   linkcount;      /*!< From WlzCoreDomain. */
  void                  *freeptr;       /*!< From WlzCoreDomain. */
  int                   nVertices;      /*!< Number of vertices. */
  int                   maxVertices;    /*!< The maximum number of vertices
                                             for which space has been
                                             allocated. */
  int                   nFaces;         /*!< Number of faces. */
  int                   maxFaces;       /*!< The maximum number of faces
                                             for which space has been
                                             allocated. */
  WlzVertexType         vtxType;        /*!< Vertex type, either WLZ_VERTEX_I3
                                             or WLZ_VERTEX_D3. */
  WlzVertex             centroid;       /*!< Centroid of the convex hull. */
  WlzVertexP            vertices;       /*!< Array of vertices. */
  int                   *faces;         /*!< Array of face vertex indices,
                                             3 vertex indices per face. */
} WlzConvHullDomain3;

/************************************************************************
* Histograms.						
***********************************************************************/
/*!
* \struct	_WlzHistogramDomain
* \ingroup	WlzHistogram
* \brief	Histograms are Woolz domains and not values as might be
*		expected.
*		Typedef: ::WlzHistogramDomain.
*/
typedef struct _WlzHistogramDomain
{
  WlzObjectType	type; 			/*!< From WlzCoreDomain. */
  int		linkcount; 	      	/*!< From WlzCoreDomain. */
  void		*freeptr; 		/*!< From WlzCoreDomain. */
  int		maxBins;	       	/*!< Number of histogram bins
  					     allocated. */
  int		nBins;			/*!< Number of histogram bins used. */
  double	origin; 	 	/*!< Lowest grey value of first
  					     histogram bin. */
  double	binSize; 	        /*!< Grey value range for a histogram
  					    bin. */
  WlzGreyP	binValues; 		/*!< Histogram values:
  					     Int for WLZ_HISTOGRAMDOMAIN_INT or
  			      		     Double for
					     WLZ_HISTOGRAMDOMAIN_FLOAT. */
} WlzHistogramDomain;

/*!
* \enum		_WlzHistFeature
* \ingroup      WlzHistogram
* \brief	Features of histograms.
*		Typedef: ::WlzHistFeature.
*/
typedef enum _WlzHistFeature
{
  WLZ_HIST_FEATURE_NONE     = (0),	/*!< No feature. */
  WLZ_HIST_FEATURE_PEAK	    = (1<<0),   /*!< Histogram peak. */
  WLZ_HIST_FEATURE_TROUGH   = (1<<1)	/*!< Histogram trough. */
} WlzHistFeature;

/************************************************************************
* Rectangle domains.						
************************************************************************/
/*!
* \struct	_WlzRect
* \ingroup	WlzFeatures
* \brief	An integer rectangle domain.
*		Side from (l[0],k[0]) to (l[1],k[1]) is a long side.
*		The vertices are cyclic.
*		Typedef: ::WlzIRect.
*/
typedef struct  _WlzRect
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int linkcount;			/*!< From WlzCoreDomain. */
  void *freeptr;			/*!< From WlzCoreDomain. */
  int irk[4];				/*!< Column vertex coordinates. */
  int irl[4];				/*!< Line vertex coordinates. */
  float rangle;			 	/*!< Angle of long side to
  					     vertical (radians). */
} WlzIRect;

/*!
* \struct	_WlzFRect
* \ingroup	WlzFeatures
* \brief	A single precision floating point rectangle domain.
*		Side from (l[0],k[0]) to (l[1],k[1]) is a long side.
*		The vertices are cyclic.
*		Typedef: ::WlzFRect.
*/
typedef struct _WlzFRect
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int linkcount;			/*!< From WlzCoreDomain. */
  void *freeptr;			/*!< From WlzCoreDomain. */
  float frk[4];				/*!< Column vertex coordinates. */
  float frl[4];				/*!< Line vertex coordinates. */
  float rangle;			 	/*!< Angle of long side to vertical
  					     (radians). */
} WlzFRect;

/************************************************************************
* Convolution and other value filters.
************************************************************************/
/*!
* \struct	_WlzConvolution
* \ingroup	WlzValuesFilters
* \brief	A 2D space domain convolution mask.
*		To reduce computational cost at the expense of data storage
*		the complete convolution is used even if highly symmetrical.
*		Typedef: ::WlzConvolution.
*/
typedef struct	_WlzConvolution
{
  WlzObjectType type;			/*!< Identifies a convolution mask. */
  int linkcount;			/*!< Reference count. */
  int xsize, ysize;			/*!< Size of mask which must be odd. */
  int *cv;			      	/*!< The convolution mask with
  					     size\f$\times\f$size elements. */
  int divscale;   	     		/*!< Scale factor by which the
  					     convolution  is divided. */
  int offset;				/*!< Offset which is added to the
  					     scaled convolution. */
  int modflag;		    		/*!< If non-zero the absolute value
  					     of the scaled, offset convolution
					     is used. */
} WlzConvolution;

/*!
* \enum		_WlzRsvFilterActionMask
* \ingroup	WlzValuesFilters
* \brief	The action to be performed by a recursive filter.
*		These values are bit masks which may be combined.
*		Typedef: ::WlzRsvFilterActionMask.
*/
typedef enum _WlzRsvFilterActionMask
{
  WLZ_RSVFILTER_ACTION_NONE     = (0),	/*!< No filtering. */
  WLZ_RSVFILTER_ACTION_X        = (1<<0), /*!< Filter along lines. */
  WLZ_RSVFILTER_ACTION_Y        = (1<<1), /*!< Filter through columns. */
  WLZ_RSVFILTER_ACTION_Z        = (1<<2)  /*!< Filter through planes. */
} WlzRsvFilterActionMask;

/*!
* \enum		_WlzRsvFilterName
* \ingroup	WlzValuesFilters
* \brief	Recursive filter types that can be used to define a recursive
* 		filter along with a filter parameter (eg sigma for Gaussian).
*		Typedef: ::WlzRsvFilterName.
*/
typedef enum _WlzRsvFilterName
{
  WLZ_RSVFILTER_NAME_NONE,      	/*!< No name, application defined
  					     filter. */
  WLZ_RSVFILTER_NAME_DERICHE_0,         /*!< Deriche's smoothing operator. */
  WLZ_RSVFILTER_NAME_DERICHE_1,         /*!< Deriche's edge operator. */
  WLZ_RSVFILTER_NAME_DERICHE_2,         /*!< Deriche's 2nd edge operator. */
  WLZ_RSVFILTER_NAME_GAUSS_0,           /*!< Gaussian. */
  WLZ_RSVFILTER_NAME_GAUSS_1,           /*!< First derivative of Gaussian. */
  WLZ_RSVFILTER_NAME_GAUSS_2            /*!< Second derivative of Gaussian. */
} WlzRsvFilterName;

/*!
* \struct	_WlzRsvFilter
* \ingroup	WlzValuesFilters
* \brief	The parameters
*		\f$a_j\f$, \f$b_j\f$ and \f$c\f$ with \f$j\in[0\cdots2]\f$
*		which define a recursive filter:
*		\f[
		  y^{+}[i] = a_0 x[i + 0] + a_1 x[i - 1] -
		             b_0 y^{+}[i - 1] - b_1 y^{+}[i - 2]
		\f]
*		\f[
		  y^{-}[i] = a_2 x[i + 1] + a_3 x[i + 2] -
		             b_0 y^{-}[i + 1] - b_1 y^{-}[i + 2]
		\f]
*		\f[
		  y[i] = c (y^{+}[i] + y^{-}[i])
		\f]
*		Typedef: ::WlzRsvFilter.
*/ 
typedef struct _WlzRsvFilter
{
  WlzRsvFilterName	name;		/*!< The filter name. */
  double		a[4];		/*!< Feed forward coefficients. */
  double		b[2];		/*!< Feed back coefficients. */
  double		c;		/*!< Normalization parameter. */
} WlzRsvFilter;
 
/************************************************************************
* Conforming mesh data structures.
************************************************************************/

/*!
* \enum		_WlzCMeshElmFlags
* \ingroup	WlzMesh
* \brief	Conforming mesh element flags. These are bit masks which are
*		used in a conforming mesh's elements flags.
*		Typedef: ::WlzCMeshElmFlags.
*/
typedef enum _WlzCMeshElmFlags
{
  WLZ_CMESH_ELM_FLAG_NONE	= (0),
  WLZ_CMESH_ELM_FLAG_BOUNDARY	= (1),	/*!< Element intersects the boundary
  					     of the domain to which it
					     should conform. */
  WLZ_CMESH_ELM_FLAG_OUTSIDE 	= (1<<1), /*!< Element is outside the domain to
  					     which the mesh should
					     conform. */
  WLZ_CMESH_ELM_FLAG_KNOWN	= (1<<2), /*!< A property of the element is
  					     known. */
  WLZ_CMESH_ELM_FLAG_ALL	= (65535) /*!< All possible flags, 0xffff
                                               decimal representation required
					       for JavaWoolz. */
} WlzCMeshElmFlags;

/*!
* \enum		_WlzCMeshNodFlags
* \ingroup	WlzMesh
* \brief	Conforming mesh node flags. These are bit masks which are
*		used in a conforming mesh's node flags.
*		Typedef: ::WlzCMeshNodFlags.
*/
typedef enum _WlzCMeshNodFlags
{
  WLZ_CMESH_NOD_FLAG_NONE	= (0),
  WLZ_CMESH_NOD_FLAG_ADJUSTED	= (1),	   /*!< Node position adjusted. */
  WLZ_CMESH_NOD_FLAG_BOUNDARY   = (1<<1),  /*!< Node is on a boundary of the
                                                mesh. */
  WLZ_CMESH_NOD_FLAG_UPWIND     = (1<<2),   /*!< Node is upwind of the active
                                                region having already been
						processed. */
  WLZ_CMESH_NOD_FLAG_ACTIVE     = (1<<3),  /*!< Node is in the active region
  						and is being processed. */
  WLZ_CMESH_NOD_FLAG_KNOWN	= (1<<4),  /*!< Property associated with
                                                node is known. */
  WLZ_CMESH_NOD_FLAG_OUTSIDE 	= (1<<5),  /*!< Node is outside the domain to
  					        which the mesh should
						conform. */
  WLZ_CMESH_NOD_FLAG_ALL	= (65535) /*!< All possible flags, 0xffff
  					       decimal representation required
					       for JavaWoolz. */
} WlzCMeshNodFlags;

/*!
* \struct	_WlzCMeshEntCore
* \ingroup	WlzMesh
* \brief	A core node/element structure containing the initial fields
* 		common to all node and element structures.
* 		Typedef: ::WlzCMeshEntCore
*/
typedef struct 	_WlzCMeshEntCore
{
  int		idx;                    /*!< The node/element index. */
  unsigned int	flags;			/*!< Bitwise description of the
  					     node/element. */
} WlzCMeshEntCore;

/*!
* \struct       _WlzCMeshNod2D
* \ingroup      WlzMesh
* \brief        A node of a 2D mesh.
*               Typedef: ::WlzCMeshNod2D.
*/
typedef struct _WlzCMeshNod2D
{
  int           idx;                    /*!< The node index from 
  					     WlzCMeshEntCore. */
  unsigned int  flags;                  /*!< Bitwise description of node
  					     from WlzCMeshEntCore. */
  WlzDVertex2   pos;                    /*!< Node position. */
  struct _WlzCMeshEdgU2D *edu;          /*!< One of many edge uses which is
                                             directed from the node. A
                                             node is shared by many parents. */
  struct _WlzCMeshNod2D *next;          /*!< Next node in bucket. */
  void		*prop;      		/*!< Node properties. */
} WlzCMeshNod2D;

/*!
* \struct       _WlzCMeshNod2D5
* \ingroup      WlzMesh
* \brief        A node of a 2D5 mesh with a 3D position but 2D connectivity.
*               Typedef: ::WlzCMeshNod2D5.
*/
typedef struct _WlzCMeshNod2D5
{
  int           idx;                    /*!< The node index from
  					     WlzCMeshEntCore. */
  unsigned int  flags;                  /*!< Bitwise description of node
  					     from WlzCMeshEntCore. */
  WlzDVertex3   pos;                    /*!< Node position. */
  struct _WlzCMeshEdgU2D5 *edu;         /*!< One of many edge uses which is
                                             directed from the node. A
                                             node is shared by many parents. */
  struct _WlzCMeshNod2D5 *next;         /*!< Next node in bucket. */
  void		*prop;      		/*!< Node properties. */
} WlzCMeshNod2D5;

/*!
* \struct       _WlzCMeshNod3D
* \ingroup      WlzMesh
* \brief        A node of a 3D mesh.
*               Typedef: ::WlzCMeshNod3D.
*/
typedef struct _WlzCMeshNod3D
{
  int           idx;                    /*!< The node index from
  					     WlzCMeshEntCore. */
  unsigned int  flags;                  /*!< Bitwise description of node
  					     from WlzCMeshEntCore. */
  WlzDVertex3   pos;                    /*!< Node position. */
  struct _WlzCMeshEdgU3D *edu;          /*!< One of many edge uses which is
                                             directed from the node. A
                                             node is shared by many parents. */
  struct _WlzCMeshNod3D *next;          /*!< Next node in bucket. */
  void		*prop;      		/*!< Node properties. */
} WlzCMeshNod3D;

/*!
* \union       _WlzCMeshNodP
* \ingroup      WlzMesh
* \brief        A node pointer for a 2 or 3D mesh.
*               Typedef: ::WlzCMeshNodP.
*/
typedef union _WlzCMeshNodP
{
  void			  *v;		/*!< Generic pointer. */
  struct _WlzCMeshEntCore *core;        /*!< Core pointer. */
  struct _WlzCMeshNod2D   *n2;		/*!< 2D node pointer. */
  struct _WlzCMeshNod2D5  *n2d5;	/*!< 2D5 node pointer. */
  struct _WlzCMeshNod3D   *n3;		/*!< 3D node pointer. */
} WlzCMeshNodP;

/*!
* \struct       _WlzCMeshEdgU2D
* \ingroup      WlzMesh
* \brief        A 2D CCW directed (half) edge within the parent simplex.
*               Typedef: ::WlzCMeshEdgU2D.
*/
typedef struct _WlzCMeshEdgU2D
{
  struct _WlzCMeshNod2D *nod;           /*!< Node from which this edge is
                                             directed. */
  struct _WlzCMeshEdgU2D *next;         /*!< Next directed edge, previous
                                             can be found using next->next. */
  struct _WlzCMeshEdgU2D *opp;          /*!< Opposite directed edge. */
  struct _WlzCMeshEdgU2D *nnxt;         /*!< Next edge directed from the
                                             same node (un-ordered). */
  struct _WlzCMeshElm2D *elm;           /*!< Parent element. */
} WlzCMeshEdgU2D;

/*!
* \struct       _WlzCMeshEdgU2D5
* \ingroup      WlzMesh
* \brief        A 2D CCW directed (half) edge within the parent simplex.
*               Typedef: ::WlzCMeshEdgU2D5.
*/
typedef struct _WlzCMeshEdgU2D5
{
  struct _WlzCMeshNod2D5 *nod;           /*!< Node from which this edge is
                                              directed. */
  struct _WlzCMeshEdgU2D5 *next;         /*!< Next directed edge, previous
                                             can be found using next->next. */
  struct _WlzCMeshEdgU2D5 *opp;          /*!< Opposite directed edge. */
  struct _WlzCMeshEdgU2D5 *nnxt;         /*!< Next edge directed from the
                                             same node (un-ordered). */
  struct _WlzCMeshElm2D5 *elm;           /*!< Parent element. */
} WlzCMeshEdgU2D5;

/*!
* \struct       _WlzCMeshEdgU3D
* \ingroup      WlzMesh
* \brief        A 3D directed (half) edge within the parent face.
*               Typedef: ::WlzCMeshEdgU3D.
*/
typedef struct _WlzCMeshEdgU3D
{
  struct _WlzCMeshNod3D *nod;           /*!< Node from which this edge is
                                             directed. */
  struct _WlzCMeshEdgU3D *next;         /*!< Next directed edge, previous
                                             can be found using next->next. */
  struct _WlzCMeshEdgU3D *nnxt;         /*!< Next edge directed from the
                                             same node (un-ordered). */
  struct _WlzCMeshFace *face;           /*!< Parent face. */
} WlzCMeshEdgU3D;

/*!
* \union       _WlzCMeshEdgUP
* \ingroup      WlzMesh
* \brief        An edge use pointer for a 2 or 3D mesh.
*               Typedef: ::WlzCMeshEdgUP.
*/
typedef union _WlzCMeshEdgUP
{
  void			  *v;		/*!< Generic pointer. */
  struct _WlzCMeshEdgU2D  *e2;		/*!< 2D node pointer. */
  struct _WlzCMeshEdgU2D5 *e2d5;	/*!< 2D5 node pointer. */
  struct _WlzCMeshEdgU3D  *e3;		/*!< 3D node pointer. */
} WlzCMeshEdgUP;

/*!
* \struct       _WlzCMeshFace
* \ingroup      WlzMesh
* \brief        A directed face within the parent simplex.
*               Typedef: ::WlzCMeshFace.
*/
typedef struct _WlzCMeshFace
{
  struct _WlzCMeshEdgU3D edu[3];        /*!< Directed edges of the face. */
  struct _WlzCMeshFace *opp;            /*!< Opposite face on neighboring
                                             mesh element. */
  struct _WlzCMeshElm3D *elm;           /*!< Parent mesh element. */
} WlzCMeshFace;

/*!
* \struct       _WlzCMeshElm2D
* \ingroup      WlzMesh
* \brief        A single 2D triangular mesh element.
*               Typedef: ::WlzCMeshElm2D.
*/
typedef struct _WlzCMeshElm2D
{
  int           idx;                    /*!< The element index from
  					     WlzCMeshEntCore. */
  unsigned int  flags;                  /*!< Element flags from
  					     WlzCMeshEntCore. */
  struct _WlzCMeshEdgU2D edu[3];        /*!< Edges of the mesh element. */
  struct _WlzCMeshCellElm2D *cElm;	/*!< First cell element from which
                                             all other cell elements can be
					     reached using the next pointer. */
  void		*prop;      		/*!< Element properties. */
} WlzCMeshElm2D;

/*!
* \struct       _WlzCMeshElm2D5
* \ingroup      WlzMesh
* \brief        A single 3D triangular mesh element.
*               Typedef: ::WlzCMeshElm2D5.
*/
typedef struct _WlzCMeshElm2D5
{
  int           idx;                     /*!< The element index from
  					      WlzCMeshEntCore. */
  unsigned int  flags;                   /*!< Element flags from
  					      WlzCMeshEntCore. */
  struct _WlzCMeshEdgU2D5 edu[3];        /*!< Edges of the mesh element. */
  struct _WlzCMeshCellElm2D5 *cElm;	 /*!< First cell element from which
                                              all other cell elements can be
					      reached using the next
					      pointer. */
  void		*prop;      		 /*!< Element properties. */
} WlzCMeshElm2D5;

/*!
* \struct       _WlzCMeshElm3D
* \ingroup      WlzMesh
* \brief        A single 3D tetrahedral mesh element.
*               Typedef: ::WlzCMeshElm3D.
*		
*		The following diagram depicts a mesh element, viewed from
*		above with the apex (node n0) pointing towards the viewer.
*		From this view, face f3 is at the rear of the tetrahedron.
* \verbatim
                                        O n1                         
                                       /|\                   
                                      / | \                  
                                     /  |  \                 
                                    /   |   \                
                                   /    |    \    f3          
                                  /     |     \  /           
                                 /      |      \/            
                                /       |      .\            
                               /        |     .  \            
                              /         |    .    \                  
                             /          |   .      \                 
                            /           O           \                
                           /      f1   / \    f0     \               
                          /          /  n0 \          \              
                         /         /         \         \             
                        /        /             \        \            
                       /       /                 \       \            
                      /      /                     \      \            
                     /     /                         \     \            
                    /    /              f2             \    \            
                   /   /                                 \   \            
                  /  /                                     \  \           
                 / /                                         \ \           
             n3 O-----------------------------------------------O n2        
\endverbatim
*		The relationship between nodes and faces depicted in
*		this diagrap is used in constructing new mesh elements.
*
*		The tetrahedron can be opened a net and viewed looking
*		at the directed edges uses associated with the faces from
*		*outside* the element surface as shown below:
* \verbatim
                                    n1                                
   n0 O-----------------------------O-----------------------------0 n0  
       \    ---------------------> / \    -------------------->  /    
        \ %         e0            / % \ %         e1            /     
         \ \                   / / /   \ \                  /  /      
          \ \                 / / /   \ \ \                /  /       
           \ \               / / /     \ \ \              /  /        
            \ \     f0      / / /       \ \ \     f1     /  /        
             \ \ e2     e1 / / /         \ \ \ e1    e0 /  /         
              \ \         / / /           \ \ \        /  /           
               \ \       / / /     f3      \ \ \      /  /             
                \ \     / / /  e0       e1  \ \ \    /  /               
                 \ \   / / /                 \ \    /  /                 
                  \   / / /                   \ \  %  /                   
                   \ % /           e2          % \   /                    
                    \ /  <--------------------    \ /                       
                     O-----------------------------O                         
                   n2 \ %  -------------------->  / n3                
                       \ \         e1            /                    
                        \ \                   / /                     
                         \ \                 / /                      
                          \ \               / /                       
                           \ \     f2      / /                        
                            \ \           / /                         
                             \ \ e0   e2 / /                          
                              \ \       / /                           
                               \ \     / /                            
                                \ \   / /                            
                                 \   / /                              
                                  \ % /                               
                                   \ /                                
                                    O                                 
                                   n0                                 
\endverbatim
*		Each mesh element has four faces and each face has three
*		directed edges, with the edges directed counter clockwise
*		(when viewed from outside the mesh element) away from their
*		nodes.
*
*		<table>
*		  <caption>
*		  Indexing of faces, edge uses and nodes in 3D elements.
*		  </caption>
*                 <tr><td>Face</td> <td>Edge Use</td> <td>Node</td></tr>
*                 <tr><td>0</td>    <td>0</td>        <td>0</td></tr>
*                 <tr><td>0</td>    <td>1</td>        <td>1</td></tr>
*                 <tr><td>0</td>    <td>2</td>        <td>2</td></tr>
*                 <tr><td>1</td>    <td>0</td>        <td>0</td></tr>
*                 <tr><td>1</td>    <td>1</td>        <td>3</td></tr>
*                 <tr><td>1</td>    <td>2</td>        <td>1</td></tr>
*                 <tr><td>2</td>    <td>0</td>        <td>0</td></tr>
*                 <tr><td>2</td>    <td>1</td>        <td>2</td></tr>
*                 <tr><td>2</td>    <td>2</td>        <td>3</td></tr>
*                 <tr><td>3</td>    <td>0</td>        <td>1</td></tr>
*                 <tr><td>3</td>    <td>1</td>        <td>1</td></tr>
*                 <tr><td>3</td>    <td>2</td>        <td>3</td></tr>
*		</table>

*		The face opposite a node, or node opposite a face
*		can easily be found using: \f$f + n = 3\f$, where
*		\f$f\f$ is the face and \f$n\f$ is the node.
*
*/
typedef struct _WlzCMeshElm3D
{
  int           idx;                    /*!< The element index from
  					     WlzCMeshEntCore. */
  unsigned int  flags;                  /*!< Element flags from
  					     WlzCMeshEntCore. */
  struct _WlzCMeshFace face[4];         /*!< Faces of the mesh element. */
  struct _WlzCMeshCellElm3D *cElm;	/*!< First cell element from which
                                             all other cell elements can be
					     reached using the next pointer. */
  void		*prop;      		/*!< Element properties. */
} WlzCMeshElm3D;

/*!
* \union       _WlzCMeshElmP
* \ingroup      WlzMesh
* \brief        A element pointer for a 2 or 3D mesh.
*               Typedef: ::WlzCMeshElmP.
*/
typedef union _WlzCMeshElmP
{
  void			  *v;		/*!< Generic pointer. */
  struct _WlzCMeshEntCore *core;        /*!< Core pointer. */
  struct _WlzCMeshElm2D   *e2;		/*!< 2D element pointer. */
  struct _WlzCMeshElm2D5  *e2d5;	/*!< 2D5 element pointer. */
  struct _WlzCMeshElm3D   *e3;		/*!< 3D element pointer. */
} WlzCMeshElmP;

/*!
* \struct	_WlzCMeshCellElm2D
* \ingroup      WlzMesh
* \brief	Data structure which is used to link lists of 2D elements
* 		with the grid cells that they intersect.
*		Typedef: ::WlzCMeshCell2D.
*/
typedef struct _WlzCMeshCellElm2D
{
  struct _WlzCMeshElm2D *elm;		/*! The element. */
  struct _WlzCMeshCell2D *cell;		/*! The cell. */
  struct _WlzCMeshCellElm2D *next;	/*! Next element intersecting cell or
  					    next cell element in the free
					    list. */
  struct _WlzCMeshCellElm2D *nextCell;	/*! Next cell which this element
  					    intersects. */
} WlzCMeshCellElm2D;

/*!
* \struct	_WlzCMeshCellElm2D5
* \ingroup      WlzMesh
* \brief	Data structure which is used to link lists of 2D5 elements
* 		with the grid cells that they intersect.
*		Typedef: ::WlzCMeshCell2D5.
*/
typedef struct _WlzCMeshCellElm2D5
{
  struct _WlzCMeshElm2D5 *elm;		/*! The element. */
  struct _WlzCMeshCell2D5 *cell;	/*! The cell. */
  struct _WlzCMeshCellElm2D5 *next;	/*! Next element intersecting cell or
  					    next cell element in the free
					    list. */
  struct _WlzCMeshCellElm2D5 *nextCell;	/*! Next cell which this element
  					    intersects. */
} WlzCMeshCellElm2D5;

/*!
* \struct	_WlzCMeshCellElm3D
* \ingroup      WlzMesh
* \brief	Data structure which is used to link lists of 3D elements
* 		with the grid cells that they intersect.
*		Typedef: ::WlzCMeshCell3D.
*/
typedef struct _WlzCMeshCellElm3D
{
  struct _WlzCMeshElm3D *elm;		/*! The element. */
  struct _WlzCMeshCell3D *cell;		/*! The cell. */
  struct _WlzCMeshCellElm3D *next;	/*! Next element intersecting cell or
  				            next cell element in the free
					    list. */
  struct _WlzCMeshCellElm3D *nextCell;	/*! Next cell which this element
  					    intersects. */
} WlzCMeshCellElm3D;

/*!
* \struct 	_WlzCMeshCell2D
* \ingroup	WlzMesh
* \brief	A single cell of a spatial grid or array of 2D cells.
*		Typedef: ::WlzCMeshCell2D.
*/
typedef struct _WlzCMeshCell2D
{
  struct _WlzCMeshNod2D *nod;		/*! Head of a linked list of nodes
  					    which are located within the
					    cell. */
  struct _WlzCMeshCellElm2D *cElm;	/*! Cell element data structure for
  					    an element which intersects this
					    cell. */
} WlzCMeshCell2D;

/*!
* \struct 	_WlzCMeshCell2D5
* \ingroup	WlzMesh
* \brief	A single cell of a spatial grid or array of 2D5 cells.
*		Typedef: ::WlzCMeshCell2D5.
*/
typedef struct _WlzCMeshCell2D5
{
  struct _WlzCMeshNod2D5 *nod;		/*! Head of a linked list of nodes
  					    which are located within the
					    cell. */
  struct _WlzCMeshCellElm2D5 *cElm;	/*! Cell element data structure for
  					    an element which intersects this
					    cell. */
} WlzCMeshCell2D5;

/*!
* \struct 	_WlzCMeshCell3D
* \ingroup	WlzMesh
* \brief	A single cell of a spatial grid or array of 3D cells.
*		Typedef: ::WlzCMeshCell3D.
*/
typedef struct _WlzCMeshCell3D
{
  struct _WlzCMeshNod3D *nod;		/*! Head of a linked list of nodes
  					    which are located within the
					    cell. */
  struct _WlzCMeshCellElm3D *cElm;	/*! Cell element data structure for
  					    an element which intersects this
					    cell. */
} WlzCMeshCell3D;

/*!
* \struct       _WlzCMeshCellGrid2D
* \ingroup      WlzMesh
* \brief        A spatial grid or array of square 2D cells that are used for
* 		fast node and element location queries.
*               Typedef: ::WlzCMeshCellGrid2D.
*/
typedef struct _WlzCMeshCellGrid2D
{
  WlzIVertex2	nCells;			/*! Dimensions of the cell grid
  					    array in terms of the number of
					    cells. */
  double	cellSz;			/*! Each cell is an axis aligned
  					    square with this side length. */
  struct _WlzCMeshCell2D **cells;	/*! Array of cells. */
  WlzCMeshCellElm2D *freeCE;	        /*! List of free cell elements for
                                            re/use. */
  AlcBlockStack	*allCE;                 /*! Allocated cell elements. */
} WlzCMeshCellGrid2D;

/*!
* \struct       _WlzCMeshCellGrid2D5
* \ingroup      WlzMesh
* \brief        A spatial grid or array of cubiod 3D cells that are used for
* 		fast 2D5 node and element location queries.
*               Typedef: ::WlzCMeshCellGrid2D5.
*/
typedef struct _WlzCMeshCellGrid2D5
{
  WlzIVertex3	nCells;			/*! Dimensions of the cell grid
  					    array in terms of the number of
					    cells. */
  double	cellSz;			/*! Each cell is an axis aligned
  					    cube with this side length. */
  struct _WlzCMeshCell2D5 ***cells;	/*! Array of cells. */
  WlzCMeshCellElm2D5 *freeCE;	        /*! List of free cell elements for
                                            re/use. */
  AlcBlockStack	*allCE;                 /*! Allocated cell elements. */
} WlzCMeshCellGrid2D5;

/*!
* \struct       _WlzCMeshCellGrid3D
* \ingroup      WlzMesh
* \brief        A spatial grid or array of square 3D cells that are used for
* 		fast node and element location queries.
*               Typedef: ::WlzCMeshCellGrid3D.
*/
typedef struct _WlzCMeshCellGrid3D
{
  WlzIVertex3	nCells;			/*! Dimensions of the cell grid
  					    array in terms of the number of
					    cells. */
  double	cellSz;			/*! Each cell is an axis aligned
  					    cube with this side length. */
  struct _WlzCMeshCell3D ***cells;	/*! Array of cells. */
  WlzCMeshCellElm3D *freeCE;	        /*! List of free cell elements for
                                            re/use. */
  AlcBlockStack	*allCE;                 /*! Allocated cell elements. */
} WlzCMeshCellGrid3D;

#ifndef WLZ_EXT_BIND
/*!
* \typedef	WlzCMeshCbFn
* \ingroup	WlzMesh
* \brief	A pointer to a function called to make mesh entity
*               properties.
*		Parameters passed are: mesh, entity, data.
*/
typedef WlzErrorNum (*WlzCMeshCbFn)(void *, void *, void *);
#endif

/*!
* \struct	_WlzCMeshCbEntry
* \ingroup	WlzMesh
* \brief	Callback entry for list of callbacks.
*		Typedef: ::WlzCMeshCbEntry.
*/
typedef struct _WlzCMeshCbEntry
{
#ifndef WLZ_EXT_BIND
  WlzCMeshCbFn	fn;
#else
  void		*fn;
#endif
  void		*data;
  struct _WlzCMeshCbEntry *next;
} WlzCMeshCbEntry;

/*!
* \struct       _WlzCMeshEntRes
* \ingroup      WlzMesh
* \brief        Resources used for efficient allocation and recycling of
*               mesh entities.
*               Typedef: ::WlzCMeshEntRes.
*/
typedef struct _WlzCMeshEntRes
{
  unsigned int  numEnt;                 /*!< Number of valid entities in
                                             vector. */
  unsigned int  maxEnt;                 /*!< Space allocated in vector. */
  unsigned int  nextIdx;                /*!< Index of next free mesh entity
                                             in vector. */
  AlcVector     *vec;                   /*!< Vector (extensible array) of
                                             mesh entities. */
  WlzCMeshCbEntry *newEntCb;		/*!< Callbacks for new entities. */
  WlzCMeshCbEntry *delEntCb;		/*!< Callbacks for deleted entities. */
} WlzCMeshEntRes;

/*!
* \struct       _WlzCMeshRes
* \ingroup      WlzMesh
* \brief        Resources used for efficient allocation, recycling and
*               location of mesh elements and nodes.
*               Typedef: ::WlzCMeshRes.
*/
typedef struct _WlzCMeshRes
{
  struct _WlzCMeshEntRes nod;            /*!< Node resources. */
  struct _WlzCMeshEntRes elm;            /*!< Element resources. */
} WlzCMeshRes;

/*!
* \union	_WlzCMeshEntP
* \ingroup	WlzMesh
* \brief	Union of pointers to top level mesh entities.
* 		Typedef: ::WlzCMeshEntP
*/
typedef union _WlzCMeshEntP
{
  void			  *v;		/*!< Generic pointer. */
  struct _WlzCMeshEntCore *core;        /*!< Core pointer. */
  struct _WlzCMeshNod2D   *n2;		/*!< 2D node pointer. */
  struct _WlzCMeshNod2D   *n2d5;	/*!< 2D5 node pointer. */
  struct _WlzCMeshNod3D   *n3;		/*!< 3D node pointer. */
  struct _WlzCMeshElm2D   *e2;		/*!< 2D element pointer. */
  struct _WlzCMeshElm2D   *e2d5;	/*!< 2D5 element pointer. */
  struct _WlzCMeshElm3D   *e3;		/*!< 3D element pointer. */
} WlzCMeshEntP;

/*!
* \union	_WlzCMeshEntPP
* \ingroup	WlzMesh
* \brief	Union of second level pointers to top level mesh entities.
* 		Typedef: ::WlzCMeshEntP
*/
typedef union _WlzCMeshEntPP
{
  void			  **v;		/*!< Generic pointer. */
  struct _WlzCMeshEntCore **core;       /*!< Core pointer. */
  struct _WlzCMeshNod2D   **n2;		/*!< 2D node pointer. */
  struct _WlzCMeshNod2D   **n2d5;	/*!< 2D5 node pointer. */
  struct _WlzCMeshNod3D   **n3;		/*!< 3D node pointer. */
  struct _WlzCMeshElm2D	  **e2;		/*!< 2D element pointer. */
  struct _WlzCMeshElm2D	  **e2d5;	/*!< 2D5 element pointer. */
  struct _WlzCMeshElm3D	  **e3;		/*!< 3D element pointer. */
} WlzCMeshEntPP;

/*!
* \struct       _WlzCMesh2D
* \ingroup      WlzMesh
* \brief        A graph based mesh model for 2D boundary conforming
*               simplical meshes.
*               The mesh inherits it's core fields from the Woolz core
*               domain.
*               Typedef: ::WlzCMesh2D.
*/
typedef struct _WlzCMesh2D
{
  int           type;                   /*!< Type of mesh. */
  int           linkcount;              /*!< Core. */
  void          *freeptr;               /*!< Core. */
  double	maxSqEdgLen;		/*!< Maximum of squared edge lengths
  					     which can be used to restrict
					     geometric searches. This may not
					     be correct if nodes have been
					     deleted or modified so it should
					     not be relied upon for any more
					     than an upper limit. */
  WlzDBox2      bBox;                   /*!< Axis aligned bounding box of
                                             the mesh. */
  WlzCMeshCellGrid2D cGrid;		/*!< Cell grid for fast node and
  					     element location queries. */
  struct _WlzCMeshRes res;              /*!< Mesh resources. */

} WlzCMesh2D;

/*!
* \struct       _WlzCMesh2D5
* \ingroup      WlzMesh
* \brief        A graph based mesh model for 2D5 boundary conforming
*               simplical meshes.
*               The mesh inherits it's core fields from the Woolz core
*               domain.
*               Typedef: ::WlzCMesh2D5.
*/
typedef struct _WlzCMesh2D5
{
  int           type;                   /*!< Type of mesh. */
  int           linkcount;              /*!< Core. */
  void          *freeptr;               /*!< Core. */
  double	maxSqEdgLen;		/*!< Maximum of squared edge lengths
  					     which can be used to restrict
					     geometric searches. This may not
					     be correct if nodes have been
					     deleted or modified so it should
					     not be relied upon for any more
					     than an upper limit. */
  WlzDBox3      bBox;                   /*!< Axis aligned bounding box of
                                             the mesh. */
  WlzCMeshCellGrid2D5 cGrid;		/*!< Cell grid for fast node and
  					     element location queries. */
  struct _WlzCMeshRes res;              /*!< Mesh resources. */

} WlzCMesh2D5;

/*!
* \struct       _WlzCMesh3D
* \ingroup      WlzMesh
* \brief        A graph based mesh model for 3D boundary conforming
*               simplical meshes.
*               The mesh inherits it's core fields from the Woolz core
*               domain.
*               Typedef: ::WlzCMesh3D.
*/
typedef struct _WlzCMesh3D
{
  int           type;                   /*!< Type of mesh. */
  int           linkcount;              /*!< Core. */
  void          *freeptr;               /*!< Core. */
  double	maxSqEdgLen;		/*!< Maximum of squared edge lengths
  					     which can be used to restrict
					     geometric searches. This may not
					     be correct if nodes have been
					     deleted or modified so it should
					     not be relied upon for any more
					     than an upper limit. */
  WlzDBox3      bBox;                   /*!< Axis aligned bounding box of
                                             the mesh. */
  WlzCMeshCellGrid3D cGrid;		/*!< Cell grid for fast node and
  					     element location queries. */
  struct _WlzCMeshRes res;              /*!< Mesh resources. */

} WlzCMesh3D;

/*!
* \union	_WlzCMeshP
* \ingroup   	WlzMesh
* \brief	Union of 2D and 3D conforming simplical mesh pointers.
*/
typedef union _WlzCMeshP
{
  void		*v;
  WlzCoreDomain *core;
  WlzCMesh2D	*m2;
  WlzCMesh2D5	*m2d5;
  WlzCMesh3D	*m3;
} WlzCMeshP;

/************************************************************************
* Functions
************************************************************************/

/*!
* \enum         _WlzFnType
* \ingroup	WlzFunction
* \brief	The types of function.
* 		Typedef: ::WlzFnType.
*/
typedef enum _WlzFnType
{
  WLZ_FN_BASIS_2DGAUSS,        		/*!< Gaussian basis function. */
  WLZ_FN_BASIS_3DGAUSS,
  WLZ_FN_BASIS_2DIMQ,			/*!< Inverse multiquadric basis
  				             function. */
  WLZ_FN_BASIS_3DIMQ,
  WLZ_FN_BASIS_2DPOLY,                   /*!< Polynomial basis function. */
  WLZ_FN_BASIS_3DPOLY,
  WLZ_FN_BASIS_2DMQ,			/*!< Multiquadric basis function. */
  WLZ_FN_BASIS_3DMQ,
  WLZ_FN_BASIS_2DTPS,              	/*!< Thin plate spline basis fn. */
  WLZ_FN_BASIS_3DTPS,
  WLZ_FN_BASIS_2DCONF_POLY,		/*!< 2D Conformal polynomial basis
  					     function. */
  WLZ_FN_BASIS_3DCONF_POLY,
  WLZ_FN_BASIS_3DMOS,			/*!< 3D Multi-order spline. */
  WLZ_FN_BASIS_SCALAR_3DMOS,		/*!< 3D Multi-order spline with scalar
                                             values. */
  WLZ_FN_SCALAR_MOD,			/*!< Modulus (abs() or fabs()). */
  WLZ_FN_SCALAR_EXP,                    /*!< Exponential (exp()). */
  WLZ_FN_SCALAR_LOG,                    /*!< Logarithm (log()). */
  WLZ_FN_SCALAR_SQRT,		        /*!< Square root (x^-1/2). */
  WLZ_FN_SCALAR_INVSQRT,		/*!< Inverse square root (x^-1/2). */
  WLZ_FN_SCALAR_SQR,			/*!< Square (x * x). */
  WLZ_FN_COUNT				/*!< Not a function but the number
  					     of functions. Keep this the
					     last of the enums! */
} WlzFnType;

/*!
* \typedef	WlzBasisEvalFn
* \ingroup	WlzFunction
* \brief	An alternative basis function evaluation function that may
*		may be called.
*/
#ifdef WLZ_EXT_BIND
typedef void *WlzBasisEvalFn;
#else /* WLZ_EXT_BIND */
typedef double (*WlzBasisEvalFn)(void *, double);
#endif /* WLZ_EXT_BIND */

/*!
* \typedef	WlzBasisDistFn
* \ingroup	WlzFunction
* \brief	An alternative basis function distance function that may
*		may be called.
*/
#ifdef WLZ_EXT_BIND
typedef void *WlzBasisDistFn;
#else /* WLZ_EXT_BIND */
typedef double (*WlzBasisDistFn)(void *, int, WlzVertex, void *);
#endif /* WLZ_EXT_BIND */

/*!
* \struct	_WlzBasisFn
* \ingroup	WlzFunction
* \brief	A basis function.
*		Typedef: ::WlzBasisFn.
*/
typedef struct _WlzBasisFn
{
  WlzFnType     type;       		/*!< The transform basis function. */
  int           nPoly;          	/*!< Polynomial order + 1. */
  int           nBasis;             	/*!< Number of basis function
  					     coefficients. */
  int           nVtx;                   /*!< Number of control point
  					     vertices. */
  int		maxVx;			/*!< Maximum number of vertices space
  					     has been allocated for. */
  WlzVertexP    poly;          		/*!< Polynomial coefficients. */
  WlzVertexP    basis;         		/*!< Basis function coefficients. */
  WlzVertexP    vertices;     		/*!< Control point vertices, these are
  					     the destination vertices of the
					     control points used to define the
					     basis function transform. */
  WlzVertexP	sVertices;		/*!< Source vertices used to compute
  					     the transform. Not always used
					     and may be NULL. */
  void		*param;			/*!< Other parameters used by the
  					     basis function, e.g. delta in
					     the MQ and Gauss basis
					     functions. Must be allocated
					     in a single block to allow
					     the parameters to be freed
					     by AlcFree().
					     For the MQ this is actually
					     \f$(\delta r)^2\f$, where \f$r\f$
					     is the range of the landmarks. */
#ifdef WLZ_EXT_BIND
  void		*evalFn;
#else
  WlzBasisEvalFn evalFn;		/*!< An alternative basis function
  					     evaluation function that may
					     be called if non NULL. */
#endif
  WlzHistogramDomain *evalData;		/*!< Data passed to the alternative
  					     basis function evaluation
					     function if the function pointer
					     is non NULL. AlcFree() will be
					     called to free the data when the
					     basis function is free'd. */
#ifdef WLZ_EXT_BIND
  void		*distFn;
  void		*mesh;
#else
  WlzBasisDistFn distFn;		/*!< An alternative basis function
  					     distance function that may
					     be called if non NULL. */
  WlzCMeshP	mesh;			/*!< Associated mesh. */
#endif
  double 	**distMap;		/*!< Array of one dimensional arrays.
					     There is a one dimensional array
					     of mesh node distances for each
					     landmark pair. Where the mesh is
					     the mesh in the associated
					     transform. In each one dimensional                                              array the node distances are
					     indexed bu the node indices.
					     There is a one dimensional array
					     for each of the control points.
					     Essentialy these arrays cache
					     distances to avoid recomputation.
					     If a distance array is no longer
					     valid it should be freed and set
					     to NULL, it will then be
					     recomputed as needed.
					     Athough the number of control
					     points may vary the number of
					     mesh nodes must remain constant. */
} WlzBasisFn;

/*!
* \struct	_WlzThreshCbStr
* \ingroup	WlzType
* \brief	Callback structure from WlzCbThreshold()
*		Typedef: ::WlzThreshCbStr.
*/
typedef struct _WlzThreshCbStr
{
  WlzPixelP	pix;
  WlzIVertex3	pos;
} WlzThreshCbStr;

/*!
* \typedef	WlzThreshCbFn
* \ingroup	WlzFunction
* \brief	Callback function for the WlzCbThreshold()
*/
#ifdef WLZ_EXT_BIND
typedef void *WlzThreshCbFn;
#else /* WLZ_EXT_BIND */
typedef int (*WlzThreshCbFn)(WlzObject *, void *, WlzThreshCbStr *);
#endif /* WLZ_EXT_BIND */

/************************************************************************
* Transforms
************************************************************************/
/*!
* \union	_WlzTransform
* \ingroup	WlzTransform
* \brief	A union of all valid transforms.
*		Typedef: ::WlzTransform.
*/
typedef union _WlzTransform
{
  struct _WlzCoreTransform *core;	/*!< Core transform. */
  struct _WlzEmptyTransform *empty;	/*!< Empty (zero) transform. */
  struct _WlzAffineTransform *affine;	/*!< Affine transforms, 2D or 3D. */
  struct _WlzBasisFnTransform *basis;	/*!< Any basis function transform. */
  struct _WlzMeshTransform *mesh;	/*!< Any convex mesh transform. */
  struct _WlzObject *obj;               /*!< Some transforms are objects
  					     with a domain and values (eg
					     conforming mesh transforms). */
} WlzTransform;

/*!
* \struct	_WlzCoreTransform
* \ingroup	WlzTransform
* \brief	The core transform, with members common to all transforms.
*		Typedef: ::WlzCoreTransform.
*/
typedef struct _WlzCoreTransform
{
  WlzTransformType type;       		/*!< From WlzCoreDomain. */
  int           linkcount;      	/*!< From WlzCoreDomain. */
  void 		*freeptr;		/*!< From WlzCoreDomain. */
} WlzCoreTransform;

/*!
* \struct	_WlzEmptyTransform
* \ingroup	WlzTransform
* \brief	An empty transform, with members common to all transforms.
* 		An empty transform is a compact represetation of a zero
* 		transform avoiding the use of NULL which implies an error.
*		Typedef: ::WlzCoreTransform.
*/
typedef struct _WlzEmptyTransform
{
  WlzTransformType type;       		/*!< From WlzCoreDomain. */
  int           linkcount;      	/*!< From WlzCoreDomain. */
  void 		*freeptr;		/*!< From WlzCoreDomain. */
} WlzEmptyTransform;

/*!
* \struct	_WlzAffineTransform
* \ingroup	WlzTransform
* \brief	Either a 2D or 3D affine transform.
*		The homogeneous matrix (mat) is always allocated as a 4x4
*		AlcDouble2Alloc style array. It is used as a 3x3
*		matrix for 2D and as a 4x4 matrix for 3D affine transforms. 
*		Typedef: ::WlzAffineTransform.
*/
typedef struct _WlzAffineTransform
{
  WlzTransformType type;       		/*!< From WlzCoreDomain. */
  int           linkcount;      	/*!< From WlzCoreDomain. */
  void 		*freeptr;		/*!< From WlzCoreDomain. */
  double        **mat;			/*!< A 4x4 homogeneous matrix which is
  					     used as a 3x3 matrix for 2D
					     transforms and as a 4x4 matrix for
					     3D affine transforms. */
} WlzAffineTransform;

/*!
* \struct	_WlzAffineTransformPrim
* \ingroup	WlzTransform
* \brief	Affine tranform primitives.
*		Typedef: ::WlzAffineTransformPrim.
*/
typedef struct _WlzAffineTransformPrim
{
  double        tx,             	/*!< X translation. */
		ty,             	/*!< Y translation. */
		tz,             	/*!< Z translation. */
		scale,          	/*!< Scale transformation. */
		theta,          	/*!< Rotation about z-axis. */
		phi,            	/*!< Rotation about y-axis. */
		alpha,          	/*!< Shear strength. */
		psi,            	/*!< Shear angle in x-y plane. */
		xsi;			/*!< 3D shear angle. */
  int           invert;                 /*!< Non-zero if reflection about
  					     the y-axis. */
} WlzAffineTransformPrim;

/*!
* \struct	_WlzBasisFnTransform
* \ingroup	WlzTransform
* \brief	A basis function transform.
* 		The delta is used by the MQ and Gauss basis functions:
*		For the MQ basis function delta = \f$R^2\f$,
*		and for the Gaussian basis function delta = \f$1/s^2\f$.
*		Typedef: ::WlzBasisFnTransform.
*/
typedef struct _WlzBasisFnTransform
{
  WlzTransformType type;       		/*!< From WlzCoreDomain. */
  int           linkcount;      	/*!< From WlzCoreDomain. */
  void 		*freeptr;		/*!< From WlzCoreDomain. */
  WlzBasisFn	*basisFn;		/*!< The basis function for the
  					     transform. */
} WlzBasisFnTransform;

/*!
* \struct	_WlzMeshNode
* \ingroup	WlzTransform
* \brief	Defines a node within a mesh transform.
*		Typedef: ::WlzMeshNode.
*/
typedef struct  _WlzMeshNode
{
  unsigned int	flags;			/*!< Mesh node flags. */
  WlzDVertex2	position;		/*!< Node position. */
  WlzDVertex2	displacement;		/*!< Node displacement. */
} WlzMeshNode;

/*!
* \struct	_WlzMeshNode3D
* \ingroup	WlzTransform
* \brief	Defines a 3D node within a mesh transform.
*  added by J. Rao 10/09/2001
*/
struct  _WlzMeshNode3D
{
  unsigned int	flags;			/*!< Mesh node flags */
  WlzDVertex3	position;		/*!< Node position */
  WlzDVertex3	displacement;		/*!< Node displacement */
};
typedef struct _WlzMeshNode3D WlzMeshNode3D;



/*!
* \struct	_WlzMeshNode2D5
* \ingroup	WlzTransform
* \brief	Defines a 2D5 node within a mesh transform.
*  added by J. Rao 23/10/2001
*/
typedef struct  _WlzMeshNode2D5
{
  unsigned int	flags;			/*!< Mesh node flags */
  WlzDVertex2	position;		/*!< Node position */
  WlzDVertex3	displacement;		/*!< Node displacement */
} WlzMeshNode2D5;

/*!
* \struct	_WlzMeshElem
* \ingroup	WlzTransform
* \brief	Defines an triangular mesh element within a mesh transform.
* 		The nodes and neighbours are indexed such that:		
* 		Neighbour 0 shares nodes 1 and 2, neighbour 1 shares nodes 2
*		and 0 and neighbour 2 shares nodes 0 and 1. All the nodes
*		are stored in counter clockwise (CCW) order.
*		Typedef: ::WlzMeshElem.
*/
typedef struct _WlzMeshElem
{
  WlzMeshElemType type;         	/*!< Type of mesh element. */
  int           idx;            	/*!< Index of this element. */
  unsigned int  flags;          	/*!< Mesh element flags. */
  int           nodes[3];       	/*!< Node indicies (CCW order). */
  int           neighbours[3];          /*!< Indicies of neighbouring
  					     elements. */
  double        strainU[3];             /*!< Constants of strain energy
  					     function. */
  double        strainA[3];		/*!< Constants of strain energy
  					     function. */
} WlzMeshElem;

/*!
* \struct	_WlzMeshElem3D 
* \ingroup	WlzTransform
* \brief	Defines an tetrahedral mesh element within a mesh transform.
* 		The nodes and neighbours are indexed such that:		
* 		Neighbour 0 shares surface ( nodes 0, 1 and 2), neighbour 1 
*               shares surface (nodes 1, 3 and 2), neighbour 2 shares surface
*               (nodes 2, 3 and 0 ) and  neighbour 3 shares surface (nodes 0, 3 and 1 )
*		All the nodes stored in the following sequence:
*               0-1-2 formed a counter clockwise (CCW) order using the dircection of a 
*               surface outwards from the tetrahedron. 0-2-3 are also stored in  
*		counter clockwise (CCW) order.				
*               added by J. Rao     10/09/2001
*
*/
typedef struct _WlzMeshElem3D
{
  WlzMeshElemType type;         	/*!< Type of mesh element */
  int           idx;            	/*!< Index of this element */
  unsigned int  flags;          	/*!< Mesh element flags */
  int           nodes[4];       	/*!< Node indicies (CCW order) */
  int           neighbours[4];          /*!< Indicies of neighbouring
  					     elements */
} WlzMeshElem3D;

/*!
* \struct	_WlzMeshTransform
* \ingroup	WlzTransform
* \brief	A mesh convex transform.
*		Typedef: ::WlzMeshElem.
*/
typedef struct _WlzMeshTransform
{
  WlzTransformType type;       		/*!< From WlzCoreDomain. */
  int           linkcount;      	/*!< From WlzCoreDomain. */
  void 		*freeptr;		/*!< From WlzCoreDomain. */
  int           nElem;          	/*!< Number of elements. */
  int           nNodes;         	/*!< Number of vertex nodes. */
  int           maxElem;        	/*!< Space allocated for elements. */
  int           maxNodes;       	/*!< Space allocated for vertex 
  					     nodes. */
  WlzMeshElem   *elements;      	/*!< Mesh elements. */
  WlzMeshNode	*nodes;			/*!< Mesh nodes. */
} WlzMeshTransform;


#ifndef WLZ_EXT_BIND
/*!
* \struct	_WlzMeshTransform3D
* \ingroup	WlzTransform
* \brief	Defines a mesh transform.
*               added by J. Rao 10/09/2001        
*/
typedef struct _WlzMeshTransform3D
{
  WlzTransformType type;       		/*!< From the core domain. */
  int           linkcount;      	/*!< From the core domain. */
  void 		*freeptr;		/*!< From the core domain. */
  int           nElem;          	/*!< Number of elements */
  int           nNodes;         	/*!< Number of vertex nodes */
  int           maxElem;        	/*!< Space allocated for elements */
  int           maxNodes;       	/*!< Space allocated for vertex 
  					     nodes */
  WlzMeshElem3D         *elements;     	/*!< Mesh elements */
  WlzMeshNode3D 	*nodes;		/*!< Mesh nodes */
} WlzMeshTransform3D;

#endif /* WLZ_EXT_BIND */

/*!
* \struct	_WlzMeshTransform2D5
* \ingroup	WlzTransform
* \brief	Defines a mesh transform.
*               added by J. Rao 23/10/2001        
*/
typedef struct _WlzMeshTransform2D5
{
  WlzTransformType type;       		/*!< From the core domain. */
  int           linkcount;      	/*!< From the core domain. */
  void 		*freeptr;		/*!< From the core domain. */
  int           nElem;          	/*!< Number of elements */
  int           nNodes;         	/*!< Number of vertex nodes */
  int           maxElem;        	/*!< Space allocated for elements */
  int           maxNodes;       	/*!< Space allocated for vertex 
  					     nodes */
  double        zConst;                 /*!< z plane const */
  WlzMeshElem           *elements;     	/*!< Mesh elements */
  WlzMeshNode2D5 	*nodes;		/*!< Mesh nodes */
} WlzMeshTransform2D5;

/************************************************************************
* User weighting functions and callback data structures for ICP based
* registration and matching.
************************************************************************/
#ifndef WLZ_EXT_BIND
/*!
* \typedef	WlzRegICPUsrWgtFn
* \ingroup	WlzTransform
* \brief	A pointer to a function called for user code weighting
*		of the matched vertices.
*/
typedef double	(*WlzRegICPUsrWgtFn)(WlzVertexType,
			     WlzAffineTransform *,
			     AlcKDTTree *,
			     WlzVertexP, WlzVertexP, WlzVertex, WlzVertex,
			     double, double, void *);

/*!
* \typedef	WlzMatchICPWeightCbData
* \ingroup	WlzTransform
* \brief	A data structure for wieghting vertex matches within
* 		WlzMatchICPWeightMatches().
*		Typedef: ::WlzMatchICPWeightCbData.
*/
typedef struct _WlzMatchICPWeightCbData
{
  WlzGMModel	*tGM;
  WlzGMModel	*sGM;
  int		nScatter;
  double 	maxDisp;
} WlzMatchICPWeightCbData;

#endif /* WLZ_EXT_BIND */

/************************************************************************
* Sequential/local transformation workspace structure.		
************************************************************************/
/*!
* \struct	_WlzSeqParWSpace
* \ingroup	WlzValuesFilters
* \brief	
*		Typedef: ::WlzSeqParWSpace.
*/
typedef struct _WlzSeqParWSpace
{
  int **adrptr;
  int kdelta;
  int ldelta;
  int brdrsz;
} WlzSeqParWSpace;

/*!
* \struct	_Wlz1DConvMask
* \ingroup	WlzValuesFilters
* \brief
*		Typedef: ::Wlz1DConvMask.
*/
typedef struct _Wlz1DConvMask
{
  int mask_size;
  int *mask_values;
  int norm_factor;
} Wlz1DConvMask;

/*!
* \struct	_WlzSepTransWSpace
* \ingroup	WlzValuesFilters
* \brief
*		Typedef: ::WlzSepTransWSpace.
*/
typedef struct _WlzSepTransWSpace
{
  WlzPixelP	inbuf;
  WlzPixelP	outbuf;
  int		len;
  WlzPixelV	bckgrnd;
} WlzSepTransWSpace;

/************************************************************************
* Standard workspace structure for interval objects.		
************************************************************************/
/*!
* \struct	_WlzIntervalWSpace
* \ingroup	WlzAccess
* \brief	The standard workspace structure for interval objects.
*		Typedef: ::WlzIntervalWSpace.
*/
typedef struct _WlzIntervalWSpace
{
  WlzObject *objaddr;			/*!< The current object. */
  int dmntype;				/*!< Domain type. */
  int lineraster;			/*!< Line scan direction as follows:
  					     <ul>
					       <li>
						1 increasing rows.
					       </li>
					       <li>
						-1 decreasing rows.
					       </li>
					     </ul> */
  int colraster;			/*!< Column scan direction as follows:
  					     <ul>
					       <li>
						1 increasing columns.
					       </li>
					       <li>
						-1 decreasing columns.
					       </li>
					     </ul> */
  WlzIntervalDomain *intdmn;		/*!< Pointer to interval structure. */
  WlzIntervalLine *intvln;	        /*!< Pointer to current line of
  					     intervals. */
  WlzInterval	*intpos;	        /*!< Pointer to current interval - 
  					     in the case of
					     WLZ_INTERVALDOMAIN_RECT this
					     is set up to point to the column
					     bounds in the interval domain
					     structure. */
  int colpos;				/*!< Column position. */
  int colrmn;				/*!< Columns remaining. */
  int linbot;				/*!< First line. */
  int linpos;				/*!< Line position. */
  int linrmn;				/*!< Lines remaining. */
  int intrmn;				/*!< Intervals remaining in line. */
  int lftpos;				/*!< Left end of interval. */
  int rgtpos;				/*!< Right end of interval. */
  int nwlpos;	  			/*!< Non-zero if new line, counts
  					     line increment since the last
  		     			     interval. */
  int plnpos;                           /*!< Plane position, for 3D domains
  					     and value tables. */
  struct _WlzGreyWSpace *gryptr;	/*!< Pointer to grey value table
  					     workspace. */
} WlzIntervalWSpace;

/************************************************************************
* Standard workspace for grey value table manipulations 		
************************************************************************/
/*!
* \struct	_WlzGreyWSpace
* \ingroup      WlzAccess
* \brief	The standard workspace for grey value table manipulations.
*		Typedef: ::WlzGreyWSpace.
*/
typedef struct _WlzGreyWSpace
{
  int gvio;				/*!< Grey value I/O switch:
  					     <ul>
					       <li>
					       0 = input to object only
					       </li>
					       <li>
					       1 = output from object only
					       </li>
					     </ul>
					    Only relevant if tranpl set, as all
					    grey-tables are unpacked. */
  int tranpl;				/*!< If non-zero, transplant values
  					     to a buffer whose address is
					     u_grintptr. Direction of
					     transplant in gvio. */
  WlzGreyType pixeltype;		/*!< Grey type. */
  WlzObjectType gdomaintype;		/*!< Value table type. */
  WlzValues gtable;			/*!< Grey value table. */
  WlzValueLine *gline;	       		/*!< Pointer to current grey table
  					     line pointer. */
  WlzTiledValueBuffer *tvb;		/*!< Tiled values buffer. */
  WlzIntervalWSpace *intptr;	      	/*!< Pointer to interval table
  					     workspace. */
  WlzGreyP u_grintptr;	    		/*!< Pointer to interval grey table.
  					     Always points to lowest order
					     column, whatever the value of
					     raster. */
} WlzGreyWSpace;


/*!
* \struct	_WlzIterateWSpace
* \ingroup	WlzAccess
* \brief	A workspace structure for interval objects which allows
* 		iteration through an object's pixels/voxels.
*		Typedef: ::WlzIterateWSpace.
*/
typedef struct _WlzIterateWSpace
{
  WlzObject	*obj;			/*!< The object being iterated
  					     through. */
  WlzObject	*obj2D;                 /*!< Object for the current plane. */
  WlzIntervalWSpace *iWSp;		/*!< Interval workspace for the current
  					     2D object. */
  WlzGreyWSpace *gWSp;			/*!< Grey workspace for the current
  					     2D object. */
  WlzRasterDir	dir;			/*!< Scanning direction. */
  int 		grey;                   /*!< Non-zero if initialised for
  					     grey values. */
  int		itvPos;			/*!< Offset into the current
  					     interval. */
  int		plnIdx;			/*!< Offset into planes of a 3D object
  					     for the current plane. */
  int		plnRmn;			/*!< Number of planes remaining for
  					     the current 2D object when
					     iterating through a 3D object. */
  WlzIVertex3	pos;			/*!< Current position. */
  WlzGreyType	gType;			/*!< The current grey type. If
  					     WLZ_GREY_ERROR then the work space
					     has not been initialised for grey
					     data. */
  WlzGreyP	gP;			/*!< Pointer to current grey value
                                             This will be NULL if grey values
					     are not appopriate, ie gType
					     is WLZ_GREY_ERROR. */
} WlzIterateWSpace;

/*!
* \struct	_WlzGreyValueWSpace
* \ingroup	WlzAccess
* \brief	Workspace for random access grey value manipulations.
*		Typedef: ::WlzGreyValueWSpace.
*/
typedef struct _WlzGreyValueWSpace
{
  WlzObjectType objType;          	/*!< Type of object, either
  					     WLZ_2D_DOMAINOBJ or
                                             WLZ_3D_DOMAINOBJ. */
  WlzDomain     domain;             	/*!< The object's domain. */
  WlzValues     values;         	/*!< The object's values. */
  WlzObjectType gTabType;		/*!< The grey table type. */
  WlzObjectType	*gTabTypes3D;		/*!< Cache for efficiency with 2D
  					     value types. Not used for
					     tiled values. */
  WlzAffineTransform *invTrans;        	/*!< If the object is a WLZ_TRANS_OBJ
					     then used to cache the inverse
					     transform. */
  WlzIntervalDomain *iDom2D;       	/*!< Current/last plane or 2D object
  					     domain. */
  WlzValues     values2D;          	/*!< Current/last plane or 2D object
  					     values. Not used for tiled
					     values. */
  int           plane;                  /*!< Current/last plane position. */
  WlzGreyType   gType;          	/*!< Grey type. */
  WlzObjectType gTabType2D;           	/*!< Current/last plane or 2D grey
  					     table. Not used for tiled
					     values. */
  WlzGreyV      gBkd;              	/*!< Background grey value, always
  					     of gType. */
  WlzGreyP      gPtr[8];        	/*!< One, four or eight grey
  					     pointers. */
  WlzGreyV      gVal[8];        	/*!< One, four or eight grey
  					     values. */
  unsigned	bkdFlag;	  	/*!< Flag set if background used
  					     with a bitmask to indicate
					     which values are background.
					     Value is 0 if there are no
					     background values. */
} WlzGreyValueWSpace;

/************************************************************************
* File I/O flags
************************************************************************/
/*!
* \enum		_WlzIOFlags
* \ingroup	WlzIO
* \brief	Flags for Woolz file I/O.
*/
typedef enum _WlzIOFlags
{
  WLZ_IOFLAGS_NONE	= (0),		/*!< No flags set. */
  WLZ_IOFLAGS_READ	= (1),		/*!< Read flag bit. */
  WLZ_IOFLAGS_WRITE	= (1<<1)	/*!< Write flag bit. */
} WlzIOFlags;

/************************************************************************
* Transform callback functions
************************************************************************/

/*!
* \typedef	WlzAffineTransformCbFn
* \ingroup	WlzTransform
* \brief	Callback function for the WlzAffineTransformCb()
* 		which may be used for value interpolation of an
* 		interval of destination values.
*
* 		The parameters are:
*		\arg cbData Callback data passed to WlzAffineTransformCb().
*		\arg gWSp   Target grey workspace.
*		\arg gVWSp  Source grey value workspace.
*		\arg invTr  Affine transform from target to source.
*		\arg pln    Plane of interval in destination.
*		\arg ln     Line of interval in destination.
*/
#ifndef WLZ_EXT_BIND
typedef WlzErrorNum (*WlzAffineTransformCbFn)(void *cbData,
    				      struct _WlzGreyWSpace *gWSp,
                                      struct _WlzGreyValueWSpace *gVWSp,
				      struct _WlzAffineTransform *invTr,
				      int pln, int ln);
#endif

/************************************************************************
* Finite Element and Warping structures.			
************************************************************************/

/*!
* \enum		_WlzElementType
* \ingroup	WlzTransform
* \brief	The types of elements in a finite element warp mesh.
*		Typedef: ::WlzElementType.
*/
typedef enum _WlzElementType
{
    WLZ_LINEAR_RECT 		= 11,	/*!< Linear and rectangular. */
    WLZ_INCOMPRESSIBLE_RECT,		/*!< Incompessible and rectangular. */
    WLZ_COMPRESSIBLE_RECT,		/*!< Compressible and rectangular. */
    WLZ_LINEAR_TRI 		= 21,	/*!< Linear and triangular. */
    WLZ_INCOMPRESSIBLE_TRI,		/*!< Incompessible and triangular. */
    WLZ_COMPRESSIBLE_TRI		/*!< Compressible and triangular. */
} WlzElementType;


#define WLZ_LINEAR 1
#define WLZ_INCOMPRESSIBLE 2
#define WLZ_COMPRESSIBLE 3

#define WLZ_RECTANGULAR 1
#define WLZ_TRIANGULAR 2

/*!
* \struct	_WlzTElement
* \ingroup	WlzTransform
* \brief	Triangular finite element warping mesh element.
*		Typedef: ::WlzTElement.
*/
typedef struct _WlzTElement
{
    WlzElementType type;	/*!< Type of element (linear, compressible,
				     incompressible). */
    int n;			/*!< Global element number. */
    int nodes[3];		/*!< Global node numbers - in anti-clockwise
    				     order! */
    float u[3];			/*!< E = u[0] if type linear. */
    float a[3];			/*!< \f$\mu\f$ = a[0] if type linear. */
} WlzTElement;

/*!
* \struct       _WlzRElement
* \ingroup	WlzTransform
* \brief        Rectangular finite element warping mesh element.
*		Typedef: ::WlzRElement.
*/
typedef struct _WlzRElement
{
    WlzElementType type;	/*!< Type of element (linear, compressible,
				     incompressible). */
    int n;			/*!< Global element number. */
    int nodes[4];		/*!< Global node numbers - in anti-clockwise
    				     order! */
    float u[3];			/*!< E = u[0] if type linear. */
    float a[3];			/*!< \f$\mu\f$ = a[0] if type linear. */
} WlzRElement;

/*!
* \struct	_WlzWarpTrans
* \ingroup	WlzTransform
* \brief	Finite element warp transformation.
* 		Typedef: ::WlzWarpTrans.
*/
typedef struct _WlzWarpTrans
{
  int type;			/*!< From WlzCoreDomain. */
  int linkcount;		/*!< From WlzCoreDomain. */
  int nelts;			/*!< Number of elements */
  int nodes;			/*!< Number of nodes. */
  WlzDVertex2 *ncoords; 	/*!< Array of nodal coordinates. */
  WlzTElement *eltlist;		/*!< List of elements. */
  WlzDVertex2 *displacements;	/*!< Array of nodal displacements. */
  float imdisp;			/*!< Max displacement in warped image. */
  float iterdisp;		/*!< Max displacement during last iteration. */
} WlzWarpTrans;

/************************************************************************
* Feature matching structures for finite element warping.
************************************************************************/

/*!
* \enum		_WlzMatchType
* \ingroup	WlzRegistration
* \brief	Finite element warping match types.
* 		Typedef: ::WlzMatchType.
*/
typedef enum _WlzMatchType
{
    WLZ_DISCARD_POINT 	= -1,
    WLZ_NODE_ATTACH   	= 0,
    WLZ_ELEMENT_ATTACH 	= 1
} WlzMatchType;

#define WLZ_MAX_NODAL_DEGREE 20

/*!
* \struct	_WlzFeatureVector
* \ingroup	WlzRegistration
* \brief	Finite element warping feature vector.
*		Typedef: ::WlzFeatureVector.
*/
typedef struct _WlzFeatureVector
{
  int direction;
  float magnitude;
  float mean1;
  float mean2;
  float std1;
  float std2;
} WlzFeatureVector;


/*!
* \struct	_WlzFeatValueLine
* \ingroup	WlzRegistration
* \brief	A line of finite element warping feature vectors.
*		Typedef: ::WlzFeatValueLine.
*/
typedef struct _WlzFeatValueLine
{
  int	 vkol1;				/*!< Left end. */
  int	 vlastkl;			/*!< Right end. */
  WlzFeatureVector *values;		/*!< Array of feature vector values. */
} WlzFeatValueLine;


/*!
* \struct	_WlzFeatValues
* \ingroup      WlzRegistration
* \brief	A ragged rectangular feature value table.
*		Typedef: ::WlzFeatValues.
*/
typedef struct _WlzFeatValue
{
  WlzObjectType type;			/*!< From WlzCoreValues. */
  int		linkcount;		/*!< From WlzCoreValues. */
  void 		*freeptr;		/*!< From WlzCoreValues. */
  WlzValues 	original_table; 	/*!< If non-NULL, the values table
  					     which owns the raw values we
					     are using. */
  int       	line1;			/*!< First line. */
  int       	lastln;			/*!< Last line. */
  int       	kol1;			/*!< First column. */
  int       	width;			/*!< Width. */
  WlzFeatureVector backgrnd;		/*!< Background value for feature
  					     vectors not in object. */
  WlzFeatValueLine *vtblines;	      	/*!< Array of feature value table line
  					     structures. */
} WlzFeatValues;
    

/*!
* \struct	_WlzRectFeatValues
* \ingroup	WlzRegistration
* \brief	A rectangular feature value table.
*		Typedef: ::WlzRectFeatValues.
*/
typedef struct	_WlzRectFeatValues
{
  WlzObjectType	type;			/*!< From WlzCoreValues. */
  int		linkcount;		/*!< From WlzCoreValues. */
  void		*freeptr;		/*!< From WlzCoreValues. */
  WlzValues 	original_table;		/*!< If non-NULL, the values table
  					     which owns the raw values we
					     are using. */
  int		line1;			/*!< First line. */
  int		lastln;			/*!< Last line. */
  int		kol1;			/*!< First column. */
  int		width;			/*!< Width. */
  WlzFeatureVector backgrnd;		/*!< Background value for feature
  				             vectors not in object. */
  WlzFeatureVector *values;		/*!< Contiguous array of feature
  					     vector values. */
} WlzRectFeatValues;
    

/*!
* \struct	_WlzFMatchPoint
* \ingroup	WlzRegistration
* \brief	Finite element warping feature match point.
* 		Typedef: ::WlzFMatchPoint.
*/
typedef struct _WlzFMatchPoint
{
  WlzMatchType	type;			/*!< Match type. */
  int 		node;			/*!< Node or element to which point
  					     attached. */
  WlzFVertex2 	ptcoords;		/*!< Coordinate of interesting
  					     point. */
  int 		elements[WLZ_MAX_NODAL_DEGREE]; /*!< Array of elements in
  					    which to search for point. */
  WlzFeatureVector *features;		/*!< Pointer to features of point. */
} WlzFMatchPoint;

/*!
* \struct	_WlzFMatchObj
* \ingroup	WlzRegistration
* \brief	A finite element warping feature match, interesting points
*		object.
* 		Typedef: ::WlzFMatchObj.
*/
typedef struct _WlzFMatchObj
{
  WlzObjectType	type;			/*!< From WlzCoreObject. */
  int		linkcount;		/*!< From WlzCoreObject. */
  int 		nopts;			/*!< Number of interesting points. */
  WlzFMatchPoint *matchpts;		/*!< Array of interesting points. */
} WlzFMatchObj;

/*!
* \struct	_Wlz3DWarpTrans
* \ingroup	WlzRegistration
* \brief	A plane-wise 3D finite element warping transform.
*		Typedef: ::Wlz3DWarpTrans.
*/
typedef struct _Wlz3DWarpTrans
{
  WlzObjectType	type;			/*!< From WlzCoreObject. */
  int		linkcount;		/*!< From WlzCoreObject. */
  WlzPlaneDomain *pdom;			/*!< Plane domain. */
  WlzFMatchObj 	**intptdoms; 		/*!< Array of pointers to interesting
  					     point objects. */
  int 		iteration;		/*!< Current iteration. */
  int 		currentplane;		/*!< Current plane. */
  float 	maxdisp;		/*!< Maximum displacement. */
  WlzPropertyList *plist;		/*!< A list of the object's
  					     properties. */
  WlzObject 	*assoc;			/*!< Associated object. */
} Wlz3DWarpTrans;

/*!
* \enum		_Wlz3DViewStructInitMask
* \ingroup	WlzTransform
* \brief	Mesh transform element flag bit masks.
*		Typedef: ::WlzMeshElemFlags.
*/
typedef enum _Wlz3DViewStructInitMask
{
  WLZ_3DVIEWSTRUCT_INIT_NONE	= (0),
  WLZ_3DVIEWSTRUCT_INIT_TRANS	= (1),	  /*!< Initialised transform */
  WLZ_3DVIEWSTRUCT_INIT_BB	= (1<<1), /*!< Initialised bounding box */
  WLZ_3DVIEWSTRUCT_INIT_LUT     = (1<<2), /*!< Initialised look-up tables */
#ifndef WLZ_EXT_BIND
  WLZ_3DVIEWSTRUCT_INIT_ALL	= 0x7	  /*!< Initialised all convenience mask */
#else
  WLZ_3DVIEWSTRUCT_INIT_ALL	= (7)     /*!< Initialised all convenience mask */
#endif
} Wlz3DViewStructInitMask;

/************************************************************************
* 3D section structure.						
************************************************************************/
/*!
* \struct	_WlzThreeDViewStruct
* \ingroup	WlzSectionTransform
* \brief	Defines a planar section through a 3D volume.
*		Typedef: ::WlzThreeDViewStruct.
*/
typedef struct _WlzThreeDViewStruct
{
  WlzObjectType	type;			/*!< Identifies the 3D view data
  					     structure: WLZ_3D_VIEW_STRUCT. */
  int		linkcount;		/*!< Core. */
  void		*freeptr;		/*!< Core. */
  int		initialised;		/*!< Non zero if the 3D view structure
  					     has been initialized. */
  WlzDVertex3	fixed;			/*!< Fixed point. */
  double	theta;			/*!< Angle of rotation about the z-axis
  					     (radians). */
  double	phi;			/*!< Angle between the viewing
  					     direction and the original
					     z-axis (radians). */
  double	zeta;			
  double	dist;			/*!< Perpendicular distance from the
  					     fixed point to the view plane. */
  double	scale;			/*!< Overall scale parameter */
  double	voxelSize[3];		/*!< Voxel rescaling if required */
  int		voxelRescaleFlg;	/*!< Voxel rescaling mode */
  WlzInterpolationType interp;		/*!< use pixel interpolation */
  WlzThreeDViewMode view_mode;		/*!< Determines the angle at which the
  					     section cut. */
  WlzDVertex3	up;			/*!< Up vector. */
  WlzDVertex3	fixed_2;		/*!< Second fixed point. */
  double	fixed_line_angle;	/*!< Angle of fixed line. */
  WlzObject	*ref_obj;
  WlzDVertex3	minvals;
  WlzDVertex3	maxvals;
  double	*xp_to_x, *xp_to_y, *xp_to_z;
  double	*yp_to_x, *yp_to_y, *yp_to_z;
/*  double	**rotation;*/
  WlzAffineTransform	*trans;		/*!< Affine transform for given
					  parameters. Could include the
					  voxel size rescaling */
} WlzThreeDViewStruct;

/*!
* \typedef	WlzProjectIntMode
* \ingroup	WlzTransform
* \brief	3D to 2D projection integration modes.
*/
typedef enum	_WlzProjectIntMode
{
  WLZ_PROJECT_INT_MODE_NONE,		/*!< No integration. */
  WLZ_PROJECT_INT_MODE_DOMAIN,		/*!< Integration using constant
  					     density within a domain. */
  WLZ_PROJECT_INT_MODE_VALUES		/*!< Integration using voxel value
  					     dependant density. */
} WlzProjectIntMode;

/*!
* \enum         _WlzKrigModelFnType
* \ingroup      WlzType
* \brief        Enumerated values for kriging variogram model functions. See
* 		the functions for details.
*/
typedef enum _WlzKrigModelFnType
{
  WLZ_KRIG_MODELFN_INVALID      = 0,    /*!< Invalid function, may be used
                                             to indicate an error. */
  WLZ_KRIG_MODELFN_NUGGET,		/*!< Nugget model function, see
  					     WlzKrigModelFnNugget(). */
  WLZ_KRIG_MODELFN_LINEAR,		/*!< Linear model function, see
  					     WlzKrigModelFnNugget(). */
  WLZ_KRIG_MODELFN_SPHERICAL,           /*!< Spherical model function, see
                                             WlzKrigModelFnSpherical(). */
  WLZ_KRIG_MODELFN_EXPONENTIAL,         /*!< Exponential model function, see
                                             WlzKrigModelFnExponential(). */
  WLZ_KRIG_MODELFN_GAUSSIAN,            /*!< Gaussian model function, see
                                             WlzKrigModelFnGaussian(). */
  WLZ_KRIG_MODELFN_QUADRATIC            /*!< Quadratic model function, see
                                             WlzKrigModelFnQuadratic(). */
} WlzKrigModelFnType;

/*!
* \struct       _WlzKrigModelFn
* \ingroup      WlzType
* \brief        Parameters and function pointer for a kriging model function.
*/
typedef struct _WlzKrigModelFn
{
  WlzKrigModelFnType    type;           /*!< Kriging model function type. */
  double                c0;             /*!< Nugget (offset) parameter. */
  double                c1;             /*!< Sill (slope) parameter. */
  double                a;              /*!< Range parameter. */
#ifdef WLZ_EXT_BIND
  void                  *fn;
#else
  double                (*fn)(struct _WlzKrigModelFn *, double h);
#endif
                                        /*!< Function pointer. */
} WlzKrigModelFn;


#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
}
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#endif	/* !WLZ_TYPE_H Don't put anything after this line */
