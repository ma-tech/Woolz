#ifndef WLZ_TYPE_H
#define WLZ_TYPE_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzType.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Defines the Woolz types. These are typedef'd
*		primatives, enumerations and structures.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************************************************************
* UBYTE: An 8 bit unsigned integer.					* 
************************************************************************/
typedef unsigned char UBYTE;

/************************************************************************
* WlzGreyType: The valid grey value types.				*
************************************************************************/
typedef enum
{
  WLZ_GREY_LONG			= 0,
  WLZ_GREY_INT			= 1,
  WLZ_GREY_SHORT		= 2,
  WLZ_GREY_UBYTE		= 3,
  WLZ_GREY_FLOAT		= 4,
  WLZ_GREY_DOUBLE		= 5,
  /**********************************************************************
  * WLZ_GREY_ERROR is not a grey type. It is an invalid grey type!	*
  * Keep it the last enumerator!					*
  **********************************************************************/
  WLZ_GREY_ERROR
} WlzGreyType;

/************************************************************************
* WlzObjectType: The (top-level) Woolz object types.			*
* Many of the integer enumeration values are required for historical	*
* reasons but should never be used explicitly.				*
************************************************************************/
typedef enum
{
  WLZ_NULL			= 0,
  WLZ_2D_DOMAINOBJ		= 1,
  WLZ_3D_DOMAINOBJ		= 2,
  WLZ_TRANS_OBJ			= 3,
  WLZ_3D_WARP_TRANS		= 4,
  WLZ_2D_POLYGON		= 10,
  WLZ_BOUNDLIST			= 11,
  WLZ_CONV_HULL			= 12,
  WLZ_HISTOGRAM			= 13,
  WLZ_3D_POLYGON		= 14,
  WLZ_RECTANGLE			= 20,
  WLZ_VECTOR_INT		= 30,
  WLZ_VECTOR_FLOAT		= 31,
  WLZ_POINT_INT			= 40,
  WLZ_POINT_FLOAT		= 41,
  WLZ_CONVOLVE_INT		= 50,
  WLZ_CONVOLVE_FLOAT		= 51,
  WLZ_DISP_FRAME		= 60,
  WLZ_DISP_GRID			= 61,
  WLZ_DISP_FRAMEX		= 62,
  WLZ_AFFINE_TRANS		= 63,
  WLZ_WARP_TRANS		= 64,
  WLZ_FMATCHOBJ			= 65,
  WLZ_TEXT			= 70,
  WLZ_COMPOUND_ARR_1		= 80,
  WLZ_COMPOUND_ARR_2		= 81,
  WLZ_COMPOUND_LIST_1		= 82,
  WLZ_COMPOUND_LIST_2		= 83,
  /**********************************************************************
  * WLZ_PROPERTY_OBJ: An object with only a property list.		*
  **********************************************************************/
  WLZ_PROPERTY_OBJ		= 110,
  /**********************************************************************
  * WLZ_EMPTY_OBJ, WLZ_EMPTY_DOMAIN, WLZ_EMPTY_VALUES: Empty entities.	*
  **********************************************************************/
  WLZ_EMPTY_OBJ			= 127,
  WLZ_EMPTY_DOMAIN,
  WLZ_EMPTY_VALUES,
  /**********************************************************************
  * WLZ_INTERVALDOMAIN_INTVL, WLZ_INTERVALDOMAIN_RECT: Interval and	*
  * rectangular interval domain types.					*
  **********************************************************************/
  WLZ_INTERVALDOMAIN_INTVL	= 1,
  WLZ_INTERVALDOMAIN_RECT	= 2,
  /**********************************************************************
  * Planedomain types.							*
  **********************************************************************/
  WLZ_PLANEDOMAIN_DOMAIN	= WLZ_2D_DOMAINOBJ,
  WLZ_PLANEDOMAIN_POLYGON	= WLZ_2D_POLYGON,
  WLZ_PLANEDOMAIN_BOUNDLIST	= WLZ_BOUNDLIST,
  WLZ_PLANEDOMAIN_CONV_HULL	= WLZ_CONV_HULL,
  WLZ_PLANEDOMAIN_HISTOGRAM	= WLZ_HISTOGRAM,
  WLZ_PLANEDOMAIN_AFFINE	= WLZ_AFFINE_TRANS,
  WLZ_PLANEDOMAIN_WARP		= WLZ_WARP_TRANS,
  /**********************************************************************
  * Value table types. A value table type encodes both the type of	*
  * table and the type of grey value it contains.			*
  * There are functions provided for extracting and synthesising these	*
  * types which should be used outside of this header file.		*
  * WLZ_GREY_TAB_RAGR: Ragged rectangle grey value table types.		*
  * WLZ_GREY_TAB_RECT: Ragged rectangular grey value table types.	*
  * WLZ_GREY_TAB_INTL: Interval grey value table types.			*
  * WLZ_FEAT_TAB_RAGR: Ragged-rectangle feature valuetable types.	*
  * WLZ_FEATVALUETABLE_RECT: Rectangular feature value table types.	*
  **********************************************************************/
  WLZ_GREY_TAB_RAGR		= 0,
  WLZ_GREY_TAB_RECT		= 1,
  WLZ_GREY_TAB_INTL		= 2,
  WLZ_FEAT_TAB_RAGR		= 5,
  WLZ_FEAT_TAB_RECT 		= 6,
  WLZ_VALUETABLE_RAGR_INT	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_INT),
  WLZ_VALUETABLE_RAGR_SHORT	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_SHORT),
  WLZ_VALUETABLE_RAGR_UBYTE	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_UBYTE),
  WLZ_VALUETABLE_RAGR_FLOAT	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_FLOAT),
  WLZ_VALUETABLE_RAGR_DOUBLE	= ((10 * WLZ_GREY_TAB_RAGR) + WLZ_GREY_DOUBLE),
  WLZ_VALUETABLE_RECT_INT	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_INT),
  WLZ_VALUETABLE_RECT_SHORT	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_SHORT),
  WLZ_VALUETABLE_RECT_UBYTE	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_UBYTE),
  WLZ_VALUETABLE_RECT_FLOAT	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_FLOAT),
  WLZ_VALUETABLE_RECT_DOUBLE	= ((10 * WLZ_GREY_TAB_RECT) + WLZ_GREY_DOUBLE),
  WLZ_VALUETABLE_INTL_INT	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_INT),
  WLZ_VALUETABLE_INTL_SHORT	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_SHORT),
  WLZ_VALUETABLE_INTL_UBYTE	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_UBYTE),
  WLZ_VALUETABLE_INTL_FLOAT	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_FLOAT),
  WLZ_VALUETABLE_INTL_DOUBLE	= ((10 * WLZ_GREY_TAB_INTL) + WLZ_GREY_DOUBLE),
  WLZ_FEATVALUETABLE_RAGR	= 50,
  WLZ_FEATVALUETABLE_RECT	= 60,
  /**********************************************************************
  * Voxel value table types.						*
  **********************************************************************/
  WLZ_VOXELVALUETABLE_GREY	= 1,
  WLZ_VOXELVALUETABLE_CONV_HULL,
  /**********************************************************************
  * Polygon domain types.						*
  **********************************************************************/
  WLZ_POLYGON_INT		= 1,
  WLZ_POLYGON_FLOAT		= 2,
  WLZ_POLYGON_DOUBLE		= 3,
  /**********************************************************************
  * Boundary list types.						*
  **********************************************************************/
  WLZ_BOUNDLIST_PIECE		= 0,
  WLZ_BOUNDLIST_HOLE		= 1,
  /**********************************************************************
  * Convex hull types.							*
  **********************************************************************/
  WLZ_CONVHULL_VALUES		= 1,
  /**********************************************************************
  * Histogram domain types. WLZ_HISTOGRAMDOMAIN_OLD_INT and		*
  * WLZ_HISTOGRAMDOMAIN_OLD_FLOAT exist only to allow old files to be	*
  * read, they should not be used anywhere else.			*
  **********************************************************************/
  WLZ_HISTOGRAMDOMAIN_OLD_INT	= 1,
  WLZ_HISTOGRAMDOMAIN_OLD_FLOAT	= 2,
  WLZ_HISTOGRAMDOMAIN_INT	= 3,
  WLZ_HISTOGRAMDOMAIN_FLOAT	= 4,
  /**********************************************************************
  * Rectangle object domain types.					*
  **********************************************************************/
  WLZ_RECTANGLE_DOMAIN_INT	= 1,
  WLZ_RECTANGLE_DOMAIN_FLOAT	= 2,
  /**********************************************************************
  * Property list types.						*
  **********************************************************************/
  WLZ_PROPERTY_SIMPLE		= 1,
  /**********************************************************************
  * 3D view structure types.						*
  **********************************************************************/
  WLZ_3D_VIEW_STRUCT		= 160,
  /* leave this last in the list */
  /**********************************************************************
  * WLZ_DUMMY_ENTRY is not an object type.				*
  * Keep it the last enumerator!					*
  **********************************************************************/
  WLZ_DUMMY_ENTRY
} WlzObjectType;

/************************************************************************
*  Raster scan directions as used by WlzGreyScan().			*
************************************************************************/
typedef enum
{
  WLZ_RASTERDIR_ILIC  = 0,	/* Increasing lines, increasing columns */
  WLZ_RASTERDIR_ILDC  = 1,	/* Increasing lines, decreasing columns */
  WLZ_RASTERDIR_DLIC  = 2,	/* Decreasing lines, increasing columns */
  WLZ_RASTERDIR_DLDC  = 3	/* Decreasing lines, decreasing columns */
} WlzRasterDir;

/************************************************************************
* Transform types.							*
************************************************************************/
typedef enum
{
  WLZ_TRANSFORM_2D_AFFINE = 1,		      /* General 2D affine transform */
  WLZ_TRANSFORM_2D_REG,	      /* 2D affine but only rotation and translation */
  WLZ_TRANSFORM_2D_TRANS,       	   /* 2D affine but only translation */
  WLZ_TRANSFORM_2D_NOSHEAR,                        /* 2D affine but no shear */
  WLZ_TRANSFORM_3D_AFFINE,		      /* General 3D affine transform */
  WLZ_TRANSFORM_3D_REG,	      /* 3D affine but only rotation and translation */
  WLZ_TRANSFORM_3D_TRANS,       	   /* 3D affine but only translation */
  WLZ_TRANSFORM_3D_NOSHEAR,                        /* 3D affine but no shear */
  WLZ_TRANSFORM_2D_BASISFN,		      /* 2D basis function transform */
  WLZ_TRANSFORM_2D5_BASISFN,   /* 2.5D (plane wise) basis function transform */
  WLZ_TRANSFORM_3D_BASISFN,                   /* 3D basis function transform */
  WLZ_TRANSFORM_2D_MESH,                     /* 2D triangular mesh transform */
  WLZ_TRANSFORM_2D5_MESH,     /* 2.5D (plane wise) triangular mesh transform */
  WLZ_TRANSFORM_3D_MESH			    /* 3D tetrahedral mesh transform */
} WlzTransformType;

/************************************************************************
* Types of basis function for basis function transforms.		*
************************************************************************/
typedef enum
{
  WLZ_BASISFN_GAUSS,                              /* Gaussian basis function */
  WLZ_BASISFN_POLY,                             /* Polynomial basis function */
  WLZ_BASISFN_MQ,			      /* Multiquadric basis function */
  WLZ_BASISFN_TPS,                       /* Thin plate spline basis function */
  WLZ_BASISFN_CONF_POLY		      /* conformal polynomial basis function */
} WlzBasisFnType;

/************************************************************************
* Mesh transform element types.						*
************************************************************************/
typedef enum
{
  WLZ_MESH_ELEM_TRILINEAR,
  WLZ_MESH_ELEM_TRIINCOMPRESSIBLE,
  WLZ_MESH_ELEM_TRICOMPRESSIBLE
} WlzMeshElemType;

/************************************************************************
* Mesh transform element flag masks.					*
************************************************************************/
typedef enum
{
  WLZ_MESH_ELEM_FLAGS_NONE      = (0),
  WLZ_MESH_ELEM_FLAGS_NBR_0     = (1),         /* Neighbour on side 0 exists */
  WLZ_MESH_ELEM_FLAGS_NBR_1     = (1<<1),      /* Neighbour on side 1 exists */
  WLZ_MESH_ELEM_FLAGS_NBR_2     = (1<<2),      /* Neighbour on side 2 exists */
  WLZ_MESH_ELEM_FLAGS_ZOMBIE	= (1<<3),      /* Dead, awaiting replacement */
  WLZ_MESH_ELEM_FLAGS_REFINE	= (1<<4)         /* Available for refinement */
} WlzMeshElemFlags;

/************************************************************************
* Mesh transform node flag masks.					*
************************************************************************/
typedef enum
{
  WLZ_MESH_NODE_FLAGS_NONE      = (0),
  WLZ_MESH_NODE_FLAGS_BBOX     	= (1),	       /* Created from bounding box. */
  WLZ_MESH_NODE_FLAGS_BLOCK	= (1<<1),          /* Created by block fill. */
  WLZ_MESH_NODE_FLAGS_IDOM	= (1<<2), /* Created to fill interval domain */
  WLZ_MESH_NODE_FLAGS_POLY	= (1<<3),    /* Created along polygon domain */
  WLZ_MESH_NODE_FLAGS_ZOMBIE	= (1<<4)       /* Dead, awaiting replacement */
} WlzMeshNodeFlags;

/************************************************************************
* Mesh generation methods.						*
************************************************************************/
typedef enum
{
  WLZ_MESH_GENMETHOD_BLOCK,	       /* Uniform (triangulated) block grid. */
  WLZ_MESH_GENMETHOD_GRADIENT  /* Triangulated grid based on image gradient. */
} WlzMeshGenMethod;

/************************************************************************
* Mesh error masks.							*
************************************************************************/
typedef enum
{
  WLZ_MESH_ERR_NONE		= (0),		/* No error, mesh valid */
  WLZ_MESH_ERR_ELEM_CW		= (1),		/* Element not CCW */
  WLZ_MESH_ERR_ELEM_INDEX	= (1<<1),	/* Element index invalid */
  WLZ_MESH_ERR_ELEM_NODE	= (1<<2),	/* Element node invalid */
  WLZ_MESH_ERR_ELEM_ZOMBIE	= (1<<3),	/* Element is a zombie */
  WLZ_MESH_ERR_DELEM_CW		= (1<<4),	/* Displaced element not CCW */
  WLZ_MESH_ERR_NELEM_INDEX	= (1<<5),	/* Neighbour index invalid */
  WLZ_MESH_ERR_NELEM_NODE	= (1<<6),	/* Neighbour node invalid */
  WLZ_MESH_ERR_NELEM_NOTNBR	= (1<<7),	/* Neighbour not a neighbour */
  WLZ_MESH_ERR_NELEM_ZOMBIE	= (1<<8)	/* Neighbour is a zombie */
} WlzMeshError;

/************************************************************************
* Vertex data types.							* 
************************************************************************/
typedef enum
{
  WLZ_VERTEX_INT		= 1,
  WLZ_VERTEX_FLOAT,
  WLZ_VERTEX_DOUBLE
} WlzVertexType;

/************************************************************************
* Connectivity.								* 
************************************************************************/
typedef enum
{
  WLZ_8_CONNECTED		= 0,
  WLZ_4_CONNECTED,
  WLZ_6_CONNECTED,
  WLZ_18_CONNECTED,
  WLZ_26_CONNECTED
} WlzConnectType;

/************************************************************************
* Distance metrics.							* 
************************************************************************/
typedef enum
{
  WLZ_8_DISTANCE		= WLZ_8_CONNECTED,
  WLZ_4_DISTANCE		= WLZ_4_CONNECTED,
  WLZ_6_DISTANCE		= WLZ_6_CONNECTED,
  WLZ_18_DISTANCE		= WLZ_18_CONNECTED,
  WLZ_26_DISTANCE		= WLZ_26_CONNECTED,
  WLZ_OCTAGONAL_DISTANCE,
  WLZ_EUCLIDEAN_DISTANCE
} WlzDistanceType;

/************************************************************************
* Special structure elements.						*
************************************************************************/
typedef enum
{
  WLZ_SPEC_STRUCT_ELM_H4,
  WLZ_SPEC_STRUCT_ELM_EX4,
  WLZ_SPEC_STRUCT_ELM_A8,
  WLZ_SPEC_STRUCT_ELM_H6,
  WLZ_SPEC_STRUCT_ELM_H5,
  WLZ_SPEC_STRUCT_ELM_H7,
  WLZ_SPEC_STRUCT_ELM_A3,
  WLZ_SPEC_STRUCT_ELM_E1,
  WLZ_SPEC_STRUCT_ELM_E2
} WlzSpecialStructElmType;

/************************************************************************
*  Binary operators.							*
************************************************************************/
typedef enum
{
  WLZ_ADD		= 0,
  WLZ_SUBTRACT,
  WLZ_MULTIPLY,
  WLZ_DIVIDE,
  WLZ_MODULUS,
  WLZ_EQ,
  WLZ_NE,
  WLZ_GT,
  WLZ_GE,
  WLZ_LT,
  WLZ_LE,
  WLZ_AND,
  WLZ_OR,
  WLZ_XOR,
  WLZ_MAX,
  WLZ_MIN,
  WLZ_MAGNITUDE
} WlzBinaryOperatorType;

/************************************************************************
*  Automatic threshold computation methods.				*
************************************************************************/
typedef enum
{
  WLZ_COMPTHRESH_FOOT,
  WLZ_COMPTHRESH_DEPTH,
  WLZ_COMPTHRESH_GRADIENT,
  WLZ_COMPTHRESH_MINIMUM
} WlzCompThreshType;

/************************************************************************
* Interpolation methods.						*
************************************************************************/
typedef enum
{
  WLZ_INTERPOLATION_NEAREST     = 0,
  WLZ_INTERPOLATION_LINEAR
} WlzInterpolationType;

/************************************************************************
* Threshold value selection.						*
************************************************************************/
typedef enum
{
  WLZ_THRESH_LOW		= 0,  		/* Threshold < thresh_value */
  WLZ_THRESH_HIGH		     	       /* Threshold >= thresh_value */
} WlzThresholdType;

/************************************************************************
* Polygon fill modes.							* 
************************************************************************/
typedef enum
{
  WLZ_SIMPLE_FILL,
  WLZ_EVEN_ODD_FILL,
  WLZ_VERTEX_FILL
} WlzPolyFillMode;

/************************************************************************
* Standard 3D views.							*
************************************************************************/
typedef enum
{
  WLZ_X_Y_VIEW,
  WLZ_Y_Z_VIEW,
  WLZ_Z_X_VIEW,
  WLZ_ARBITRARY_VIEW
} WlzThreeDStdViews;

/************************************************************************
* 3D viewing modes.							*
************************************************************************/
typedef enum
{
  WLZ_STATUE_MODE,
  WLZ_UP_IS_UP_MODE,
  WLZ_FIXED_LINE_MODE,
  WLZ_ZERO_ZETA_MODE,
  WLZ_ZETA_MODE
} WlzThreeDViewMode;

/************************************************************************
* Window functions.							* 
************************************************************************/
typedef enum
{
  WLZ_WINDOWFN_NONE,
  WLZ_WINDOWFN_BLACKMAN,
  WLZ_WINDOWFN_HAMMING,
  WLZ_WINDOWFN_HANNING,
  WLZ_WINDOWFN_PARZEN,
  WLZ_WINDOWFN_RECTANGLE,
  WLZ_WINDOWFN_WELCH,
  /**********************************************************************
  * WLZ_WINDOWFN_UNSPECIFIED is not a window function, it is an error	*
  * Keep it the last in the enumeration.				*
  **********************************************************************/
  WLZ_WINDOWFN_UNSPECIFIED
} WlzWindowFnType;

/************************************************************************
* Sampling functions.							* 
************************************************************************/
typedef enum
{
  WLZ_SAMPLEFN_NONE, 				     /* No sampling function */
  WLZ_SAMPLEFN_POINT,				 	   /* Point sampling */
  WLZ_SAMPLEFN_MEAN,				/* Mean value sample of data */
  WLZ_SAMPLEFN_GAUSS,	                 /* Gaussian weighted sample of data */
  WLZ_SAMPLEFN_MIN,				   /* Minimum value sampling */
  WLZ_SAMPLEFN_MAX,				   /* Maximum value sampling */
  WLZ_SAMPLEFN_MEDIAN				    /* Median value sampling */
} WlzSampleFn;

/************************************************************************
* Vertices 								*
************************************************************************/
typedef struct
{
  int   	vtY;
  int   	vtX;
} WlzIVertex2;

typedef struct
{
  float 	vtY;
  float 	vtX;
} WlzFVertex2;

typedef struct
{
  double 	vtY;
  double 	vtX;
} WlzDVertex2;


typedef struct
{
  int		vtX;
  int		vtY;
  int		vtZ;
} WlzIVertex3;

typedef struct
{
  float		vtX;
  float		vtY;
  float		vtZ;
} WlzFVertex3;

typedef struct
{
  double	vtX;
  double	vtY;
  double	vtZ;
} WlzDVertex3;

/************************************************************************
* Bounding boxes.							*
************************************************************************/
typedef struct
{
  int		xMin;
  int		yMin;
  int		xMax;
  int		yMax;
} WlzIBox2;

typedef struct
{
  double	xMin;
  double	yMin;
  double	xMax;
  double	yMax;
} WlzDBox2;

typedef struct
{
  int		xMin;
  int		yMin;
  int		zMin;
  int		xMax;
  int		yMax;
  int		zMax;
} WlzIBox3;

typedef struct
{
  double	xMin;
  double	yMin;
  double	zMin;
  double	xMax;
  double	yMax;
  double	zMax;
} WlzDBox3;

/************************************************************************
* Woolz values union.							* 
************************************************************************/
typedef union
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

/************************************************************************
* Woolz domain union.							* 
************************************************************************/
typedef union
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
  struct _WlzWarpTrans	     *wt;
} WlzDomain;

/************************************************************************
* Simple property .							*
************************************************************************/
typedef struct
{
  WlzObjectType		type;
  int 			linkcount;
  void 			*freeptr;
  unsigned long		size;
  void 			*prop;
} WlzSimpleProperty;

/************************************************************************
* WlzObject: The Woolz object.						*
************************************************************************/
typedef struct _WlzObject
{
  WlzObjectType      type;
  int                linkcount;
  WlzDomain          domain;
  WlzValues          values;
  WlzSimpleProperty  *plist;
  struct _WlzObject *assoc;
} WlzObject;

/************************************************************************
* WlzCoreDomain: All Woolz domains must have all the fields of the core	*
* domain, in the same order and before any others.			*
************************************************************************/
typedef struct _WlzCoreDomain
{
  WlzObjectType   type;
  int             linkcount;
  void            *freeptr;
} WlzCoreDomain;

/************************************************************************
* WlzPlaneDomain: A 3D domain.						*
************************************************************************/

typedef struct _WlzPlaneDomain
{
  WlzObjectType     type;					     /* CORE */
  int               linkcount;					     /* CORE */
  void              *freeptr;					     /* CORE */
  int               plane1;     		   /* First plane coordinate */
  int               lastpl;     		    /* Last plane coordinate */
  int               line1;      		    /* First line coordinate */
  int               lastln;			     /* Last line coordinate */
  int               kol1;       	    /* First column  line coordinate */
  int               lastkl;     	     /* Last column  line coordinate */
  WlzDomain 	    *domains;       /* Array of pointers to interval domains */
  float 	    voxel_size[3];      /* Array of nominal voxel dimensions */
} WlzPlaneDomain;


/************************************************************************
* WlzIntervalDomain: 2D interval domain.				*
* The intervals in a line must be contiguous.				*
************************************************************************/
typedef struct _WlzIntervalDomain
{
  WlzObjectType     type;					     /* CORE */
  int               linkcount;					     /* CORE */
  void              *freeptr;					     /* CORE */
  int               line1;			    /* First line coordinate */
  int               lastln;			     /* Last line coordinate */
  int               kol1;			  /* First column coordinate */
  int               lastkl;			   /* Last column coordinate */
  struct _WlzIntervalLine *intvlines;   /* Array of interval line structures */
} WlzIntervalDomain;

/************************************************************************
* A line of intervals.							* 
************************************************************************/
typedef struct _WlzIntervalLine
{
  int                  nintvs;
  struct _WlzInterval *intvs;
} WlzIntervalLine;

/************************************************************************
*  A single interval.							*
************************************************************************/
typedef struct _WlzInterval
{
  int ileft;
  int iright;
} WlzInterval;

/************************************************************************
* A union of pointers to grey values.					*
************************************************************************/
typedef union
{
  long *lnp;
  int  *inp;
  short *shp;
  UBYTE *ubp;
  float *flp;
  double *dbp;
} WlzGreyP;

/************************************************************************
* A union of grey values.						*
************************************************************************/
typedef union
{
  long lnv;
  int inv;
  short shv;
  UBYTE ubv;
  float flv;
  double dbv;
} WlzGreyV;

/************************************************************************
* WlzPixelV, WlzPixelP: Pixel Structures for single typed values 	*
* eg background.							*
************************************************************************/
typedef struct
{
  WlzGreyType	type;
  WlzGreyV	v;
} WlzPixelV;

typedef struct
{
  WlzGreyType	type;
  WlzGreyP	p;
} WlzPixelP;

/************************************************************************
* Structure of grey values in a line .					*
************************************************************************/
typedef struct
{
  int      vkol1;					/* Relative left end */
  int      vlastkl;				       /* Relative right end */
  WlzGreyP values;					  /* Array of values */
} WlzValueLine;

/************************************************************************
* WlzCoreValues: All Woolz value tables must have all the fields of 	*
* the core values, in the same order and before any others.		*
************************************************************************/
typedef struct _WlzCoreValues
{
  WlzObjectType type;	
  int       linkcount;	
} WlzCoreValues;

/************************************************************************
* Ragged rectangle value table. The type encodes both the type of	*
* value table and the type of grey value.				*
************************************************************************/
typedef struct _WlzRagRValues
{
  WlzObjectType type;						     /* CORE */
  int       linkcount;						     /* CORE */
  void      *freeptr;
  WlzValues original_table; 		   /* If != NULL, valuetable which owns
				              the raw values we are using */
  int       line1;					       /* First line */
  int       lastln;					        /* Last line */
  int       kol1;					     /* First column */
  int       width;						    /* Width */
  WlzPixelV bckgrnd;		/* Background value for points not in object */
  WlzValueLine *vtblines;	      /* Array of grey table line structures */
} WlzRagRValues;

/************************************************************************
* Rectangular values table. The type encodes both the type of		*
* value table and the type of grey value.				*
************************************************************************/
typedef struct _WlzRectValues
{
  WlzObjectType type;						     /* CORE */
  int 		linkcount;					     /* CORE */
  void 		*freeptr;
  WlzValues 	original_table; 	   /* If != NULL, valuetable which owns
  					      the raw values we are using */
  int 		line1;					       /* First line */
  int 		lastln;						/* Last line */
  int 		kol1;					     /* First column */
  int 		width;						    /* Width */
  WlzPixelV 	bckgrnd; 	/* Background value for points not in object */
  WlzGreyP 	values;			       /* Contiguous array of values */
} WlzRectValues;

/************************************************************************
* One line's worth of grey-intervals 					*
************************************************************************/
typedef struct
{
  int           nintvs;				 /* Number of grey-intervals */
  WlzValueLine *vtbint;				/* Pointer to grey intervals */
} WlzValueIntervalLine;

/************************************************************************
* Interval-structured value table. The type encodes both the type of	*
* value table and the type of grey value.				*
************************************************************************/
typedef struct _WlzIntervalValues
{
  WlzObjectType type;						     /* CORE */
  int      	linkcount;      				     /* CORE */
  void    	*freeptr;
  WlzValues 	original_table;		   /* If != NULL, valuetable which owns
  					      the raw values we are using */
  int   	line1;	      				       /* First line */
  int    	lastln;	      					/* Last line */
  int     	kol1;	      				     /* First column */
  int      	width;	      					    /* Width */
  WlzPixelV 	bckgrnd;        /* Background value for points not in object */
  WlzValueIntervalLine *vil;   /* Pointers to structures of grey table lines */
} WlzIntervalValues;

/************************************************************************
* Voxel value table.							*
************************************************************************/
typedef struct _WlzVoxelValues
{
  WlzObjectType	type;						     /* CORE */
  int           linkcount;					     /* CORE */
  void*         freeptr;
  WlzValues 	original_table;		   /* If != NULL, valuetable which owns
  					      the raw values we are using */
  int           plane1;					      /* First plane */
  int           lastpl;					       /* Last plane */
  WlzPixelV     bckgrnd;	/* Background value for points not in object */
  WlzValues     *values;   		/* Array of pointers to value tables */
} WlzVoxelValues;

/************************************************************************
* 2D Polygon domain.							*
************************************************************************/
typedef struct _WlzPolygonDomain
{
  WlzObjectType type;						     /* CORE */
  int linkcount;						     /* CORE */
  void *freeptr;						     /* CORE */
  int nvertices;				      /* Number of verticies */
  int maxvertices;	   /* Maximum number of verticies (space allocated). */
  WlzIVertex2 *vtx; /* Array of verticies, may need casting according to type */
} WlzPolygonDomain;

/************************************************************************
* 3D Polygon domain.							*
************************************************************************/
typedef struct
{
  WlzObjectType type;						     /* CORE */
  int 	linkcount;						     /* CORE */
  void *freeptr;						     /* CORE */
  int nvertices;				      /* Number of verticies */
  int maxvertices;	   /* Maximum number of verticies (space allocated). */
  WlzIVertex2 *vtx; /* Array of verticies, may need casting according to type */
} WlzPolygonDomain3;

/************************************************************************
* Boundary list.							*
************************************************************************/
typedef struct _WlzBoundList
{
  WlzObjectType type;						     /* CORE */
  int linkcount;						     /* CORE */
  void *freeptr;						     /* CORE */
  struct _WlzBoundList *up;		 /* The containing hole or piece, NULL
				   	    if the universal hole (very top) */
  struct _WlzBoundList *next;		/* Next hole or piece at same level and
				   	   lying within same piece or hole,
					   NULL if no more at this level */
  struct _WlzBoundList *down;		/* First enclosed structure, NULL if
  				           none */
  int wrap;				/* Wrap number - number of points of
				   	   boundary included both at start and
					   end of polygon representation */
  WlzPolygonDomain *poly;	       /* Polygon representation of boundary */
} WlzBoundList;


/************************************************************************
* Chord and convex hull parameters.					*
* The chord equation is: ccon = (acon * x) - (bcon * y).		*
************************************************************************/
typedef struct _WlzChord 
{
  int sig;			     /* Non-zero if judged to be significant */
  int acon;			 		  /* Chord equation paramter */
  int bcon;			    		  /* Chord equation paramter */
  int ccon;					  /* Chord equation paramter */
  double cl;					         /* Chord length, x8 */
  int bl;			   /* Line number of bay bottom or bulge top */
  int bk;	                 /* Column number of bay bottom or bulge top */
  int barea;					    /* Bay or bulge area, *8 */
  int bd;			/* Bay maximum depth or bulge max height, *8 */
} WlzChord;

typedef struct _WlzConvHullValues
{
  WlzObjectType type;						     /* CORE */
  int           linkcount;					     /* CORE */
  void		*freeptr;					     /* CORE */
  WlzValues original_table; 	           /* If != NULL, valuetable which owns
				   	      the raw values we are using */
  int           nchords;				 /* Number of chords */
  int           nsigchords;		     /* Number of significant chords */
  int           mdlin;		  /* Mid-line of enclosed originating object */
  int           mdkol;		/* Mid-column of enclosed originating object */
  WlzChord      *ch;
} WlzConvHullValues;

/************************************************************************
* Histogram domain.							*
***********************************************************************/
typedef struct _WlzHistogramDomain
{
  WlzObjectType	type; 			/* WLZ_HISTOGRAMDOMAIN_(INT)|(FLOAT) */
  int		linkcount; 	      /* Number of times data structure used */
  void		*freeptr; /* Points to (void *)values, used by WlzFreeDomain */
  int		maxBins;	       /* Number of histogram bins allocated */
  int		nBins;			    /* Number of histogram bins used */
  double	origin; 	 /* Lowest grey value of first histogram bin */
  double	binSize; 	       /* Grey value range for histogram bin */
  WlzGreyP	binValues; /* Histogram values: int for WLZ_HISTOGRAMDOMAIN_INT
  			      or double for WLZ_HISTOGRAMDOMAIN_FLOAT */
} WlzHistogramDomain;

/************************************************************************
* Rectangle domains.							*
* Side from (l[0],k[0]) to (l[1],k[1]) is a long side. 			*
* The vertices are cyclic.						*
************************************************************************/
typedef struct  _WlzRect
{
  WlzObjectType type;						     /* CORE */
  int linkcount;						     /* CORE */
  void *freeptr;						     /* CORE */
  int irk[4];					/* Column vertex coordinates */
  int irl[4];					  /* Line vertex coordinates */
  float rangle;			 /* Angle of long side to vertical (radians) */
} WlzIRect;

typedef struct _WlzFRect
{
  WlzObjectType type;						     /* CORE */
  int linkcount;						     /* CORE */
  void *freeptr;						     /* CORE */
  float frk[4];					/* Column vertex coordinates */
  float frl[4];					  /* Line vertex coordinates */
  float rangle;			 /* Angle of long side to vertical (radians) */
} WlzFRect;

/************************************************************************
* Vectors.								*
************************************************************************/
typedef struct 
{
  WlzObjectType type;
  int linkcount;
  int k1;
  int l1;
  int k2;
  int l2;
  int style;
} WlzIVector;

typedef struct
{
  WlzObjectType type;
  int linkcount;
  float k1;
  float l1;
  float k2;
  float l2;
  int style;
} WlzFVector;

/************************************************************************
* Points.								*
************************************************************************/
typedef struct
{
  WlzObjectType type;
  int linkcount;
  int k;
  int l;
  int style;
} WlzIPoint;

typedef struct
{
  WlzObjectType type;
  int linkcount;
  float k;
  float l;
  int style;
} WlzFPoint;


/************************************************************************
* Convolution mask.							*
* To reduce computational cost at the expense of data storage the	*
* complete convolution must be specified even if highly symmetrical.	*
************************************************************************/
typedef struct
{
  WlzObjectType type;
  int linkcount;
  int xsize, ysize;			       /* Size of mask (must be odd) */
  int *cv;			      /* Size*size convolution mask elements */
  int divscale;   	     /* Divide by this after forming the convolution */
  int offset;     				 /* ... then add this offset */
  int modflag;		    /* ... and take the modulus if this is non-zero. */
} WlzConvolution;
 
/************************************************************************
* Affine tranform: WLZ_TRANSFORM_2D_AFFINE or WLZ_TRANSFORM_3D_AFFINE.	*
************************************************************************/
typedef struct _WlzAffineTransform
{
  WlzTransformType type;        		  		     /* CORE */
  int           linkcount;      				     /* CORE */
  void *freeptr;						     /* CORE */
  double        tx,             			    /* X translation */
		ty,             			    /* Y translation */
		tz,             			    /* Z translation */
		scale,          		     /* Scale transformation */
		theta,          		    /* Rotation about z-axis */
		phi,            	            /* Rotation about y-axis */
		alpha,          			   /* Shear strength */
		psi,            		 /* Shear angle in x-y plane */
		xsi;					   /* 3D shear angle */
  int           invert;               /* Non-zero if reflection about y-axis */
  double        **mat;
} WlzAffineTransform;

/************************************************************************
* Basis function transform: WLZ_TRANSFORM_2D_BASISFN,			*
* WLZ_TRANSFORM_2D5_BASISFN or WLZ_TRANSFORM_3D_BASISFN.		*
* The delta is used by the MQ and Gauss basis functions:  For the MQ	*
* basis fn delta = R^2, and for the Gaussian basis fn delta = 1/s^2.	*
************************************************************************/
typedef struct _WlzBasisFnTransform
{
  WlzTransformType type;					     /* CORE */
  int           linkcount;					     /* CORE */
  void          *freeptr;       				     /* CORE */
  WlzBasisFnType basisFn;       	     /* The transform basis function */
  int           nPoly;          	             /* Polynomial order + 1 */
  int           nBasis;             /* Number of basis function coefficients */
  int           nVtx;                   /* Number of control point verticies */
  double	delta;		  /* Used by the MQ and Gauss basis functions*/
  WlzDVertex2    *poly;          		  /* Polynomial coefficients */
  WlzDVertex2    *basis;         	      /* Basis function coefficients */
  WlzDVertex2    *verticies;     		  /* Control point verticies */
} WlzBasisFnTransform;

/************************************************************************
* Mesh transforms: WLZ_TRANSFORM_2D_MESH, WLZ_TRANSFORM_2D5_MESH or	*
* WLZ_TRANSFORM_3D_MESH.						*
* The nodes and neighbours are indexed such that:			*
* Neighbour 0 shares nodes 1 and 2, neighbour 1 shares nodes 2 and 0	*
* and neighbour 2 shares nodes 0 and 1. All the nodes are stored in	*
* counter clockwise (CCW) order.					*
************************************************************************/
typedef struct _WlzMeshNode
{
  unsigned int	flags;					  /* Mesh node flags */
  WlzDVertex2	position;				    /* Node position */
  WlzDVertex2	displacement;				/* Node displacement */
} WlzMeshNode;

typedef struct _WlzMeshElem
{
  WlzMeshElemType type;         		     /* Type of mesh element */
  int           idx;            		    /* Index of this element */
  unsigned int  flags;          		       /* Mesh element flags */
  int           nodes[3];       		/* Node indicies (CCW order) */
  int           neighbours[3];          /* Indicies of neighbouring elements */
  double        strainU[3];          /* Constants of strain energy function: */
  double        strainA[3];
} WlzMeshElem;

typedef struct _WlzMeshTransform
{
  WlzTransformType type;					     /* CORE */
  int           linkcount;      				     /* CORE */
  void          *freeptr;       				     /* CORE */
  int           nElem;          		       /* Number of elements */
  int           nNodes;         		   /* Number of vertex nodes */
  int           maxElem;        	     /* Space allocated for elements */
  int           maxNodes;       	 /* Space allocated for vertex nodes */
  WlzMeshElem   *elements;      			    /* Mesh elements */
  WlzMeshNode	*nodes;					       /* Mesh nodes */
} WlzMeshTransform;

/************************************************************************
* Sequential/local transformation workspace structure.			*
************************************************************************/
typedef struct 
{
  int **adrptr;
  int kdelta;
  int ldelta;
  int brdrsz;
} WlzSeqParWSpace;

typedef struct
{
  int mask_size;
  int *mask_values;
  int norm_factor;
} Wlz1DConvMask;

typedef struct
{
  WlzPixelP	inbuf;
  WlzPixelP	outbuf;
  int		len;
  WlzPixelV	bckgrnd;
} WlzSepTransWSpace;

/************************************************************************
* Standard workspace structure for interval objects.			*
************************************************************************/
typedef struct
{
  WlzObject *objaddr;				       /* The current object */
  int dmntype;						      /* Domain type */
  int lineraster;			    /* Line scan direction as follows:
					        1 increasing rows
					       -1 decreasing rows */
  int colraster;			  /* Column scan direction as follows:
					      1 increasing columns
					     -1 decreasing columns */
  WlzIntervalDomain *intdmn;		    /* Pointer to interval structure */
  WlzIntervalLine *intvln;	     /* Pointer to current line of intervals */
  WlzInterval	*intpos;	     /* Pointer to current interval - in the
  				        case of WLZ_INTERVALDOMAIN_RECT this
					is set up to point to the column bounds
					in the interval domain structure */
  int colpos;						  /* Column position */
  int colrmn;						/* Columns remaining */
  int linbot;						       /* First line */
  int linpos;						    /* Line position */
  int linrmn;						 /* Lines rnemaining */
  int intrmn;				      /* Intervals remaining in line */
  int lftpos;					     /* Left end of interval */
  int rgtpos;					    /* Right end of interval */
  int nwlpos;	  /* Non-zero if new line, counts line increment since the last
  		     interval */
  struct _WlzGreyWSpace *gryptr;	  /* Pointer to grey table workspace */
} WlzIntervalWSpace;

/************************************************************************
* Standard workspace for grey-table manipulations 			*
************************************************************************/
typedef struct _WlzGreyWSpace
{
  int gvio;				 /* Grey value i/o switch :
					    0 = input to object only
					    1 = output from object only
					    Only relevant if tranpl set, as all
					    grey-tables are unpacked. */
  int tranpl;				 /* If non-zero, transplant values to a
					    buffer whose address is u_grintptr.
					    Direction of transplant in gvio */
  WlzGreyType pixeltype;					/* Grey type */
  WlzObjectType gdomaintype;				/* Value table type. */
  WlzValues gtable;					 /* Grey value table */
  WlzValueLine *gline;	       /* Pointer to current grey table line pointer */
  WlzIntervalWSpace *intptr;	      /* Pointer to interval table workspace */
  WlzGreyP u_grintptr;	    /* Pointer to interval grey table. Always points to
  			       lowest order column, whatever the value of
			       raster */
} WlzGreyWSpace;

/********************************************************************
 * Workspace for random access grey value manipulations
 ********************************************************************/
typedef struct _WlzGreyValueWSpace
{
  WlzObjectType objType;          /* Type of object: Either WLZ_2D_DOMAINOBJ or
                                     WLZ_3D_DOMAINOBJ */
  WlzDomain     domain;         		      /* The object's domain */
  WlzValues     values;         		      /* The object's values */
  WlzObjectType	*gTabTypes3D;	/* efficiency hack while value types are
				   computed */
  WlzAffineTransform *invTrans;        /* If the object is a WLZ_TRANS_OBJ then
                                   	  the inverse transform is set */
  WlzIntervalDomain *iDom2D;       /* Current/last plane or 2D object domain */
  WlzValues     values2D;          /* Current/last plane or 2D object values */
  int           plane;                        /* Current/last plane position */
  WlzGreyType   gType;          			        /* Grey type */
  WlzObjectType gTabType2D;           /* Current/last plane or 2D grey table */
  WlzGreyV      gBkd;              /* Background grey value, always of gType */
  WlzGreyP      gPtr[8];        	 /* One, four or eight grey pointers */
  WlzGreyV      gVal[8];        	   /* One, four or eight grey values */
} WlzGreyValueWSpace;

/************************************************************************
* Compound objects: implemented as either ARRAYS or LINKED LISTS of	*
* other objects.  There is a distinction between an compound of the 	*
* same type (e.g. resulting from a labelling) and a compound of		*
* different types (e.g. resulting from a range of image processes from	*
* a single original object).						*
************************************************************************/
typedef struct 
{
  WlzObjectType type;
  int           linkcount;
  WlzObjectType otype;		       /* The permitted type if constrained. */
  int           n;		 /* Maximum number of objects (array length) */
  WlzObject     **o;			/* The list of woolz object pointers */
  WlzSimpleProperty *p;			     /* A non-specific property list */
  WlzObject     *assoc;
} WlzCompoundArray;


/************************************************************************
* Finite Element and Warping structures.				*
************************************************************************/

typedef enum
{
    WLZ_LINEAR_RECT 		= 11,
    WLZ_INCOMPRESSIBLE_RECT,
    WLZ_COMPRESSIBLE_RECT,
    WLZ_LINEAR_TRI 		= 21,
    WLZ_INCOMPRESSIBLE_TRI,
    WLZ_COMPRESSIBLE_TRI
} WlzElementType;


#define WLZ_LINEAR 1
#define WLZ_INCOMPRESSIBLE 2
#define WLZ_COMPRESSIBLE 3

#define WLZ_RECTANGULAR 1
#define WLZ_TRIANGULAR 2

typedef struct
{
    WlzElementType type;	/* type of element (linear, compressible,
				   incompressible) */
    int n;			/* global element number */
    int nodes[3];		/* global node numbers - in anti-clockwise
    				   order! */
    float u[3];			/* E = u[0] if type linear */
    float a[3];			/* mu = a[0] if type linear */
} WlzTElement;

typedef struct
{
    WlzElementType type;	/* type of element (linear, compressible,
				   incompressible) */
    int n;			/* global element number */
    int nodes[4];		/* global node numbers - in anti-clockwise
    				   order! */
    float u[3];			/* E = u[0] if type linear */
    float a[3];			/* mu = a[0] if type linear */
} WlzRElement;

typedef struct _WlzWarpTrans
{
    int type;			/* WLZ_WARP_TRANS */
    int linkcount;
    int nelts;			/* number of elements */
    int nodes;			/* number of nodes */
    WlzDVertex2 *ncoords; 	/* array of nodal coordinates */
    WlzTElement *eltlist;	/* list of elements */
    WlzDVertex2 *displacements;	/* array of nodal displacements */
    float imdisp;		/* max displacement in warped image */
    float iterdisp;		/* max displacement during last iteration */
} WlzWarpTrans;

/************************************************************************
* Feature matching structures.						*
************************************************************************/
typedef enum
{
    WLZ_DISCARD_POINT 	= -1,
    WLZ_NODE_ATTACH   	= 0,
    WLZ_ELEMENT_ATTACH 	= 1
} WlzMatchType;

#define WLZ_MAX_NODAL_DEGREE 20

typedef struct
{
    int direction;
    float magnitude;
    float mean1;
    float mean2;
    float std1;
    float std2;
} WlzFeatureVector;


typedef struct
{
    int	 vkol1;				/* left end */
    int	 vlastkl;			/* right end */
    WlzFeatureVector *values;		/* array of values */
} WlzFeatValueLine;


typedef struct _WlzFeatValues	/* ragged rectangular featurevalue table */
{
  WlzObjectType 	type;	/* should be 50 */
  int		 	linkcount;
  void 			*freeptr;
  WlzValues 	original_table; 	/* If != NULL, valuetable which owns
				   	   the raw values we are using */
  int		       	line1;
  int			lastln;
  int			kol1;
  int			width;
  WlzFeatureVector 	backgrnd;
  WlzFeatValueLine 	*vtblines;
} WlzFeatValues;
    

typedef struct _WlzRectFeatValues	/* rectangular */
{
  WlzObjectType	type;			/* should be 60 */
  int			linkcount;
  void		*freeptr;
  WlzValues original_table; 		/* If != NULL, valuetable which owns
				   	   the raw values we are using */
  int			line1;
  int			lastln;
  int			kol1;
  int			width;
  WlzFeatureVector 	backgrnd;
  WlzFeatureVector 	*values;
} WlzRectFeatValues;
    

typedef struct
{
    WlzMatchType type;
    int node;			/* node or element to which point attached */
    WlzFVertex2 ptcoords;	/* coordinate of interesting point */
    int elements[WLZ_MAX_NODAL_DEGREE];   /* list of elements in which to
					      search for point */
    WlzFeatureVector *features;	/* pointer to features of point */
} WlzFMatchPoint;

typedef struct
{
    WlzObjectType type;	/* WLZ_FMATCHOBJ */
    int linkcount;
    int nopts;			/* number of interesting points */
    WlzFMatchPoint *matchpts;	/* list of interesting points */
} WlzFMatchObj;

typedef struct
{
  WlzObjectType		type;			/* WLZ_3D_WARP_TRANS */
  int			linkcount;
  WlzPlaneDomain 	*pdom;
  WlzFMatchObj 		**intptdoms; /* array of pointers to interesting point
					 lists */
  int 			iteration;
  int 			currentplane;
  float 		maxdisp;
  WlzSimpleProperty 	*plist;
  WlzObject 		*assoc;
} Wlz3DWarpTrans;

/************************************************************************
* 3D section structure.							*
************************************************************************/
typedef struct {
  WlzObjectType	type;
  int		linkcount;
  void		*freeptr;
  int		initialised;
  WlzDVertex3	fixed;
  double	theta;
  double	phi;
  double	zeta;
  double	dist;
  double	scale;
  WlzThreeDViewMode	view_mode;
  WlzDVertex3	up;
  WlzDVertex3	fixed_2;
  double	fixed_line_angle;
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
