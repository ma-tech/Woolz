#ifndef WLZ_TYPE_H
#define WLZ_TYPE_H
#pragma ident "MRC HGU $Id$"
/*!**********************************************************************
* \file         WlzType.h
* \author       Bill Hill
* \date         April 2001
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzType
* \brief        Defines the Woolz types. These are enumerations and
*		structures which have been typedef'd.
* \bug          None known.
* \note
* Maintenance log with most recent changes at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
* \typedef	UBYTE
* \ingroup	WlzType
* \brief	An eight bit unsigned integer.
*/
typedef unsigned char UBYTE;
/*!
* \typedef	UINT
* \ingroup	WlzType
* \brief        A 32 bit unsigned integer.
*/
typedef unsigned int  UINT;

/*!
* \enum		_WlzGreyType
* \ingroup	WlzType
* \brief	The valid grey value types.
*		Typedef: ::WlzGreyType.
*/
typedef enum _WlzGreyType
{
  WLZ_GREY_LONG			= 0,	/*!< Signed long integer. */
  WLZ_GREY_INT			= 1,	/*!< Signed integer. */
  WLZ_GREY_SHORT		= 2,	/*!< Signed short. */
  WLZ_GREY_UBYTE		= 3,	/*!< Unsigned byte. */
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
* \enum		_WlzObjectType
* \ingroup	WlzType
* \brief	The Woolz object types.
*		Many of the integer enumeration values are required for
*		historical reasons but, with the execption of WLZ_NULL,
*		the numerical value should never be used explicitly.
*
*		The base grey table types are used to form explicit
*		grey table types which include the grey type. The functions
*		for extracting and synthesising these types should be
*		used: These are WlzGreyTableType(),
*		WlzGreyTableTypeToGreyType() and
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
  WLZ_RECTANGLE			= 20,	/*!< Rectangle. */
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
  WLZ_PROPERTY_OBJ		= 110,	/*!< An object which only has a
  					     property list. */
  WLZ_EMPTY_OBJ			= 127,	/*!< Empty object: An object which
  					     occupies no space and has no
					     values. */
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
  WLZ_GREY_TAB_RAGR		= 0,	/*!< Base ragged rectangle grey value
  					     table. */
  WLZ_GREY_TAB_RECT		= 1,	/*!< Base rectangular grey value
  					     table. */
  WLZ_GREY_TAB_INTL		= 2,	/*!< Base interval grey value table. */
  WLZ_FEAT_TAB_RAGR		= 5,	/*!< Base ragged rectangle feature
  					     table. */
  WLZ_FEAT_TAB_RECT 		= 6,	/*!< Base rectangular feature table. */
  WLZ_VALUETABLE_RAGR_INT	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_INT),
  					/*!< Ragged rectangle int value
					     table. */
  WLZ_VALUETABLE_RAGR_SHORT	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_SHORT),
  					/*!< Ragged rectangle short value
					     table. */
  WLZ_VALUETABLE_RAGR_UBYTE	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_UBYTE),
  					/*!< Ragged rectangle unsigned byte
					     value table. */
  WLZ_VALUETABLE_RAGR_FLOAT	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_FLOAT),
  					/*!< Ragged rectangle single precision
					     floating point value table. */
  WLZ_VALUETABLE_RAGR_DOUBLE	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_DOUBLE),
  					/*!< Ragged rectangle double precision
					     floating point value table. */
  WLZ_VALUETABLE_RAGR_BIT	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_BIT),
  					/*!< Ragged rectangle single bit (packed
					     in unsigned bytes) value table. */
  WLZ_VALUETABLE_RAGR_RGBA	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_RGBA),
  					/*!< Ragged rectangle red, green, blue,
					     alpha value table. */
  WLZ_VALUETABLE_RECT_INT	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_INT),
  					/*!< Rectangular int value
					     table. */
  WLZ_VALUETABLE_RECT_SHORT	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_SHORT),
  					/*!< Rectangular short value
					     table. */
  WLZ_VALUETABLE_RECT_UBYTE	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_UBYTE),
  					/*!< Rectangular unsigned byte value
					     table. */
  WLZ_VALUETABLE_RECT_FLOAT	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_FLOAT),
  					/*!< Rectangular single precision
					     floating point value table. */
  WLZ_VALUETABLE_RECT_DOUBLE	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_DOUBLE),
  					/*!< Rectangular double precision
					     floating point value table. */
  WLZ_VALUETABLE_RECT_BIT	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_BIT),
  					/*!< Rectangular single bit (packed in
					     unsigned bytes) value table. */
  WLZ_VALUETABLE_RECT_RGBA	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_RGBA),
  					/*!< Rectangular red, green, blue,
					     alpha value table.  */ 
  WLZ_VALUETABLE_INTL_INT	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_INT),
  					/*!< Interval coded int value
					     table. */
  WLZ_VALUETABLE_INTL_SHORT	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_SHORT),
  					/*!< Interval coded short value
					     table. */
  WLZ_VALUETABLE_INTL_UBYTE	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_UBYTE),
  					/*!< Interval coded unsigned byte value
					     table. */
  WLZ_VALUETABLE_INTL_FLOAT	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_FLOAT),
  					/*!< Interval coded single precision
					     floating point value table. */
  WLZ_VALUETABLE_INTL_DOUBLE	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_DOUBLE),
  					/*!< Interval coded double precision
					     floating point value table. */
  WLZ_VALUETABLE_INTL_BIT	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_BIT),
  					/*!< Interval coded single bit (packed
					     unsigned bytes) value table. */
  WLZ_VALUETABLE_INTL_RGBA	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_RGBA),
  					/*!< Interval coded red, green, blue,
					     alpha value table.  */
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
  * 3D view structure types.					
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
  WLZ_EMAP_PROPERTY_ANATOMY_DOMAIN,	/*!< Anatomy Domain */
  WLZ_EMAP_PROPERTY_OTHER_DOMAIN,	/*!< Other Domain */
  WLZ_EMAP_PROPERTY_DUMMY		/*!< Dummy property */
} WlzEMAPPropertyType;

/*!
* \enum		_WlzRasterDir
* \ingroup	WlzAccess
* \brief	Raster scan directions as used by WlzGreyScan().
*		Typedef: ::WlzRasterDir.
*/
typedef enum _WlzRasterDir
{
  WLZ_RASTERDIR_ILIC  = 0,		/*!< Increasing lines, increasing
  					     columns */
  WLZ_RASTERDIR_ILDC  = 1,		/*!< Increasing lines, decreasing
  					     columns */
  WLZ_RASTERDIR_DLIC  = 2,		/*!< Decreasing lines, increasing
  					     columns */
  WLZ_RASTERDIR_DLDC  = 3		/*!< Decreasing lines, decreasing
  					     columns */
} WlzRasterDir;

/*!
* \enum		_WlzTransformType
* \ingroup	WlzTransform
* \brief	Types of spatial transformation.
*		Typedef: ::WlzTransformType.
*/
typedef enum _WlzTransformType
{
  WLZ_TRANSFORM_EMPTY = 0,		/*!< Undefined transform. */
  WLZ_TRANSFORM_2D_AFFINE = 1,		/*!< General 2D affine transform */
  WLZ_TRANSFORM_2D_REG,	      		/*!< 2D affine but only rotation
  					     and translation */
  WLZ_TRANSFORM_2D_TRANS,       	/*!< 2D affine but only translation */
  WLZ_TRANSFORM_2D_NOSHEAR,             /*!< 2D affine but no shear */
  WLZ_TRANSFORM_3D_AFFINE,		/*!< General 3D affine transform */
  WLZ_TRANSFORM_3D_REG,	      		/*!< 3D affine but only rotation
  					     and translation */
  WLZ_TRANSFORM_3D_TRANS,       	/*!< 3D affine but only translation */
  WLZ_TRANSFORM_3D_NOSHEAR,             /*!< 3D affine but no shear */
  WLZ_TRANSFORM_2D_BASISFN,		/*!< 2D basis function transform */
  WLZ_TRANSFORM_2D5_BASISFN,   		/*!< 2.5D (plane wise) basis function
  					     transform */
  WLZ_TRANSFORM_3D_BASISFN,             /*!< 3D basis function transform */
  WLZ_TRANSFORM_2D_MESH,                /*!< 2D triangular mesh transform */
  WLZ_TRANSFORM_2D5_MESH,     		/*!< 2.5D (plane wise) triangular
  					     mesh transform */
  WLZ_TRANSFORM_3D_MESH			/*!< 3D tetrahedral mesh transform */
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
  WLZ_MESH_ELEM_FLAGS_REFINE	= (1<<4)  /*!< Available for refinement */
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
  WLZ_MESH_GENMETHOD_GRADIENT  		/*!< Triangulated grid based on image
  					     gradient. */
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
  WLZ_EUCLIDEAN_DISTANCE
} WlzDistanceType;

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
  WLZ_COMPTHRESH_FOOT,	/*!< The threshold value is intercept of a line fitted
			 *  to the upper slope of the histogram main peak with
			 *  the abscissa.
			 */
  WLZ_COMPTHRESH_DEPTH, /*!< The threshold value is that point to the right of
  			 *  the histogram peak that is maximally distant from
			 *  the chord joining the peak and the histogram right
			 *  hand end point.
			 *  The histogram may need to be smoothed for this
			 *  algorithm to work.
			 */
  WLZ_COMPTHRESH_GRADIENT /*!< The threshold value is the first point to the
			 *  right of the histogram main peak at which the
			 *  gradient falls to zero (cf finding a minimum).
			 *  To find the slope of the histogram at some
			 *  point a straight line is fitted through the
			 *  point \f$\pm\f$ a fixed number of points to
			 * either side. */
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
  WLZ_INTERPOLATION_CLASSIFY_1		/*!< Classification by probability. */
} WlzInterpolationType;

/*!
* \enum		_WlzThresholdType
* \ingroup	WlzThreshold
* \brief	Threshold value selection.
* 		Typedef: ::WlzThresholdType.
*/
typedef enum _WlzThresholdType
{
  WLZ_THRESH_LOW		= 0, /*!< Threshold < thresh_value */
  WLZ_THRESH_HIGH		     /*!< Threshold >= thresh_value */
} WlzThresholdType;

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
  WLZ_GREYTRANSFORMTYPE_LINEAR,		/*!< linear interpolation */
  WLZ_GREYTRANSFORMTYPE_GAMMA,		/*!< gamma function */
  WLZ_GREYTRANSFORMTYPE_EXPONENTIAL,	/*!< exponential function */
  WLZ_GREYTRANSFORMTYPE_SIGMOID		/*!< sigmoid function */
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
* \enum		_WlzVertexType
* \ingroup	WlzType
* \brief	2D and 3D vertex types.
*		Typedef: ::WlzVertexType.
*/
typedef enum _WlzVertexType
{
  WLZ_VERTEX_I2		= 1,		/*!< 2D integer vertex. */
  WLZ_VERTEX_F2,			/*!< 2D single precision floating
  					     point vertex. */
  WLZ_VERTEX_D2,			/*!< 2D double precision floating
  					     point vertex. */
  WLZ_VERTEX_I3,			/*!< 3D integer vertex. */
  WLZ_VERTEX_F3,			/*!< 3D single precision floating
  					     point vertex. */
  WLZ_VERTEX_D3				/*!< 3D double precision floating
  					     point vertex. */
} WlzVertexType;

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
  WlzFVertex2	*f2;
  WlzDVertex2	*d2;
  WlzIVertex3	*i3;
  WlzFVertex3	*f3;
  WlzDVertex3	*d3;
} WlzVertexP;

/*!
* \union	_WlzVertex
* \ingroup	WlzType
* \brief	Union of vertex values.
*		Typedef: ::WlzVertexV.
*/
typedef union _WlzVertex
{
  WlzIVertex2	i2;
  WlzFVertex2	f2;
  WlzDVertex2	d2;
  WlzIVertex3	i3;
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
* \union	_WlzBoxP
* \ingroup	WlzType
* \brief	Union of axis aligned box pointers.
* 		Typedef: ::WlzBoxP.
*/
typedef union _WlzBoxP
{
  void		*v;
  WlzIBox2	*i2;
  WlzDBox2	*d2;
  WlzIBox3	*i3;
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
  WlzDBox2	d2;
  WlzIBox3	i3;
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
  long *lnp;
  int  *inp;
  short *shp;
  UBYTE *ubp;
  float *flp;
  double *dbp;
  UINT *rgbp;
} WlzGreyP;

/*!
* \union	_WlzGreyV
* \ingroup	WlzType
* \brief	A union of grey values.
*		Typedef: ::WlzGreyV.
*/
typedef union _WlzGreyV
{
  long lnv;
  int inv;
  short shv;
  UBYTE ubv;
  float flv;
  double dbv;
  UINT rgbv;
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
  WLZ_GMMOD_3D
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
  WLZ_GMELM_VERTEX_G3I,
  WLZ_GMELM_VERTEX_G3D,
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
  struct _WlzGMVertexG3I *vertexG3I;
  struct _WlzGMVertexG3D *vertexG3D;
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
* \brief	The geometric properties of point in 2D integer space.
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
* \brief	The geometric properties of point in 2D double precision space.
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
* \struct	_WlzGMVertexG3I
* \ingroup	WlzGeoModel
* \brief	The geometric properties of point in 3D integer space.
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
* \brief	The geometricproperties of point in 3D double precision space.
*		Typedef: ::WlzGMVertexG3I.
*/
typedef struct _WlzGMVertexG3D
{
  WlzGMElemType type;			/*!< WLZ_GMELM_VERTEX_G3D */
  int		idx;		        /*!< Unique identifier for vertex
  				             geometry */
  WlzDVertex3	vtx;			/*!< Where the point lies in space */
} WlzGMVertexG3D;

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
  WlzGMVertexG3I *vg3I;
  WlzGMVertexG3D *vg3D;
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
* \brief	A line or curve between a pair of verticies.
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
* \struct	_WlzGMShellR
* \ingroup	WlzGeoModel
* \brief	The resources used by a model.
*		Typedef: ::WlzGMModelR.
*/
typedef struct _WlzGMShellR
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
  WLZ_CONTOUR_MTD_BND			/*!< Object boundary. */
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
} WlzValues;

/*!
* \union	_WlzDomain
* \ingroup	WlzType
* \brief	The union of Woolz domains.
*		Typedef: ::_WlzDomain.
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
  int			theilerStage;			/*!< Theiler stage */
  char			modelName[EMAP_PROPERTY_MODELNAME_LENGTH];/*!< Volume model name */
  char			version[EMAP_PROPERTY_VERSION_LENGTH];/*!< Model version */
  char			*fileName;			/*!< Original filename (not very useful) */
  long			creationTime;			/*!< Original creation time */
  char			creationAuthor[EMAP_PROPERTY_AUTHORNAME_LENGTH];/*!< Creation author */
  char			creationMachineName[EMAP_PROPERTY_AUTHORNAME_LENGTH];/*!< Original creation machine name - useful for location */
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
  int           n;		 	/*!< The maximum number of objects
  					     (array length). */
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

/************************************************************************
* Grey value tables.
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
  void*         freeptr;		/*!< From WlzCoreValues. */
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

/************************************************************************
* Polygon domains.						
************************************************************************/
/*!
* \struct	_WlzPolygonDomain
* \ingroup	WlzPolyline
* \brief	A 2D polyline domain with possible types:WLZ_POLYGON_INT, 
*		WLZ_POLYGON_FLOAT  or WLZ_POLYGON_DOUBLE. 
*		Typedef: ::WlzPolygonDomain.
*/
typedef struct _WlzPolygonDomain
{
  WlzObjectType type;			/*!< From WlzCoreDomain. */
  int linkcount;			/*!< From WlzCoreDomain. */
  void *freeptr;			/*!< From WlzCoreDomain. */
  int nvertices;			/*!< Number of verticies. */
  int maxvertices;			/*!< The maximum number of verticies
  					     for which space has been
					     allocated. */
  WlzIVertex2 *vtx; 			/*!< Array of verticies.
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
  int nvertices;			/*!< Number of verticies. */
  int maxvertices;	   		/*!< The maximum number of verticies
  					     for which space has been
					     allocated. */
  WlzIVertex2 *vtx; 			/*!< Array of verticies.
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


/************************************************************************
* Chord and convex hull parameters.				
************************************************************************/
/*!
* \struct	_WlzChord
* \ingroup	WlzConvexHull
* \brief	A cord described by the equation:
*		\f[
		  c = a x - b y
		\f]
*		Typedef: ::WlzChord.
*/
typedef struct _WlzChord 
{
  int sig;			     	/*!< Non-zero if judged to be
  					     significant. */
  int acon;			 	/*!< Chord equation paramter,
  					     \f$a\f$. */
  int bcon;			    	/*!< Chord equation paramter,
  					     \f$b\f$. */
  int ccon;				/*!< Chord equation paramter,
  					     \f$c\f$. */
  double cl;				/*!< Eight \f$\times\f$ chord
  					     length. */
  int bl;			   	/*!< Line number of bay bottom or
  					     bulge top. */
  int bk;	                 	/*!, Column number of bay bottom or
  					     bulge top. */
  int barea;				/*!< Eight \f$\times\f$ bay or
  					     bulge area. */
  int bd;				/*!< Eight \f$\times\f$ bay maximum
  					     depth or bulge max height. */
} WlzChord;

/*!
* \struct	_WlzConvHullValues
* \ingroup	WlzConvexHull
* \brief	A 2D convex hull.
* 		Typedef: ::WlzConvHullValues.
*/
typedef struct _WlzConvHullValues
{
  WlzObjectType type;			/*!< From WlzCoreValues. */
  int           linkcount;		/*!< From WlzCoreValues. */
  void		*freeptr;		/*!< From WlzCoreValues. */
  WlzValues original_table; 	        /*!< If non-NULL, valuetable which
  					     owns the raw values we are
					     using. */
  int           nchords;		/*!< Number of chords. */
  int           nsigchords;		/*!< Number of significant chords. */
  int           mdlin;		  	/*!< Mid-line of enclosed
  					     originating object. */
  int           mdkol;			/*!< Mid-column of enclosed
  					     originating object. */
  WlzChord      *ch;
} WlzConvHullValues;

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
  WLZ_FN_BASIS_2DGAUSS,        		/*!< Gaussian basis function */
  WLZ_FN_BASIS_2DPOLY,                   /*!< Polynomial basis function */
  WLZ_FN_BASIS_2DMQ,			/*!< Multiquadric basis function */
  WLZ_FN_BASIS_2DTPS,              	/*!< Thin plate spline basis
  					     function */
  WLZ_FN_BASIS_2DCONF_POLY,		/*!< Conformal polynomial basis
  					     function */
  WLZ_FN_COUNT				/*!< Not a function but the number
  					     of functions. Keep this the
					     last in the enums! */
} WlzFnType;

/*!
* \struct	_WlzBasisFn
* \ingroup	WlzFunction
* \brief	A basis function.
*		Typedef: ::WlzAffineTransformPrim.
*/
typedef struct _WlzBasisFn
{
  WlzFnType     type;       		/*!< The transform basis function. */
  int           nPoly;          	/*!< Polynomial order + 1. */
  int           nBasis;             	/*!< Number of basis function
  					     coefficients. */
  int           nVtx;                   /*!< Number of control point
  					     verticies. */
  WlzVertexP    poly;          		/*!< Polynomial coefficients. */
  WlzVertexP    basis;         		/*!< Basis function coefficients. */
  WlzVertexP    vertices;     		/*!< Control point vertices. */
  void		*param;			/*!< Other parameters used by the
  					     basis function, e.g. delta in
					     the MQ and Gauss basis
					     functions. Must be allocated
					     in a single block to allow
					     the parameters to be freed
					     by AlcFree(). */
} WlzBasisFn;

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
  struct _WlzAffineTransform *affine;	/*!< Affine transforms, 2D or 3D. */
  struct _WlzBasisFnTransform *basis;	/*!< Any basis function transform. */
  struct _WlzMeshTransform *mesh;	/*!< Any mesh transform. */
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
  WlzBasisFn	*basisFn;			/*!< The basis function for the
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
* \struct	_WlzMeshTransform
* \ingroup	WlzTransform
* \brief	A mesh transform.
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

/************************************************************************
* User weighting function for ICP based registration.
************************************************************************/
#ifndef WLZ_EXT_BIND
/*!
* \typedef	WlzRegICPUsrWgtFn
* \ingroup	WlzTransform
* \brief	A pointer to a function called for user code weighting
*		of the matched verticies.
*/
typedef double	(*WlzRegICPUsrWgtFn)(WlzVertexType,
			     WlzAffineTransform *,
			     AlcKDTTree *,
			     WlzVertexP, WlzVertexP, WlzVertex, WlzVertex,
			     void *);
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
  WlzIntervalWSpace *intptr;	      	/*!< Pointer to interval table
  					     workspace. */
  WlzGreyP u_grintptr;	    		/*!< Pointer to interval grey table.
  					     Always points to lowest order
					     column, whatever the value of
					     raster. */
} WlzGreyWSpace;

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
  WlzDomain     domain;         	/*!< The object's domain. */
  WlzValues     values;         	/*!< The object's values. */
  WlzObjectType	*gTabTypes3D;		/*!< Cache for efficiency with 2D
  					     value types. */
  WlzAffineTransform *invTrans;        	/*!< If the object is a WLZ_TRANS_OBJ
					     then used to cache the inverse
					     transform. */
  WlzIntervalDomain *iDom2D;       	/*!< Current/last plane or 2D object
  					     domain. */
  WlzValues     values2D;          	/*!< Current/last plane or 2D object
  					     values. */
  int           plane;                  /*!< Current/last plane position. */
  WlzGreyType   gType;          	/*!< Grey type. */
  WlzObjectType gTabType2D;           	/*!< Current/last plane or 2D grey
  					     table. */
  WlzGreyV      gBkd;              	/*!< Background grey value, always
  					     of gType. */
  WlzGreyP      gPtr[8];        	/*!< One, four or eight grey
  					     pointers. */
  WlzGreyV      gVal[8];        	/*!< One, four or eight grey
  					     values. */
  int		bkdFlag;	  	/*!< Flag set to 1 if background used
  					     else 0. */
} WlzGreyValueWSpace;

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
  double	scale;
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
  double	**rotation;
} WlzThreeDViewStruct;

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif	/* !WLZ_TYPE_H Don't put anything after this line */
