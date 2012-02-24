#ifndef WLZ_EXTFFTYPE_H
#define WLZ_EXTFFTYPE_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFType_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFType.h
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
* \brief	Header file with data types for external data file format
*		support for the MRC Human Genetics Unit Woolz library.
* \ingroup	WlzExtFF
*/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
* Can be used to initialise any WlzExtFF enumaerated type to an invalid
* value.
*/
#define WLZEFF_INVALID -1

/*!
* \enum		_WlzEffFormat
* \ingroup	WlzExtFF
* \brief	enumeration of the file formats supported.
*		Typedef: ::WlzEffFormat
*/
typedef enum _WlzEffFormat
{
  WLZEFF_FORMAT_NONE = 0,	/*!< None, used to indicating no match (equal
  				     to zero).  */
  WLZEFF_FORMAT_BMP,		/*!< Microscoft bitmap. */
  WLZEFF_FORMAT_DEN,		/*!< Stanford density. */
  WLZEFF_FORMAT_ICS,		/*!< International cytometry standard.  */
  WLZEFF_FORMAT_PNM,		/*!< Portable any map. */
  WLZEFF_FORMAT_PIC,		/*!< Biorad confocal pic format. */
  WLZEFF_FORMAT_SLC,		/*!< SLC volume files. */
  WLZEFF_FORMAT_VFF,		/*!< Sunvision volumes. */
  WLZEFF_FORMAT_VTK,		/*!< Visualization Toolkit. */
  WLZEFF_FORMAT_WLZ,		/*!< Woolz. */
  WLZEFF_FORMAT_IPL,		/*!< IP Lab. */
  WLZEFF_FORMAT_TIFF,		/*!< Tiff. */
  WLZEFF_FORMAT_RAW,		/*!< Raw data. */
  WLZEFF_FORMAT_AM,		/*!< Amira. */
  WLZEFF_FORMAT_JPEG,		/*!< Jpeg. */
  WLZEFF_FORMAT_ANL,		/*!< Analyze 7.5. */
  WLZEFF_FORMAT_GIF,		/*!< Graphics Interchange Format. */
  WLZEFF_FORMAT_MESH,           /*!< Pascal Frey's tetrahedral mesh format. */
  WLZEFF_FORMAT_NODEELE,        /*!< Jonathan Shewchuk's mesh format. */
  WLZEFF_FORMAT_VMESH,          /*!< GRUMMP tetrahedral mesh format. */
  WLZEFF_FORMAT_PLY2,           /*!< Riken PLY2 mesh format. */
  WLZEFF_FORMAT_OBJ,            /*!< Wavefront geometry format. */
  WLZEFF_FORMAT_TXT,		/*!< Simple ASCII text listing, csv format. */
  WLZEFF_FORMAT_NIFTI,		/*!< Neuroimaging Informatics Technology
  				     Initiative format based on Analyze 7.5. */
  WLZEFF_FORMAT_SMESH,          /*!< GRUMMP surface mesh format. */
  WLZEFF_FORMAT_EMT,            /*!< Netgen neutral mesh format. */
  WLZEFF_FORMAT_STL,            /*!< Stereolithography file format. */
  WLZEFF_FORMAT_COUNT 		/*!< Keep last: Number of formats (including
  				     WLZEFF_FORMAT_NONE). */
} WlzEffFormat;

#ifndef WLZ_EXT_BIND

#define WLZEFF_BMP_MAGIC_0		('B')
#define WLZEFF_BMP_MAGIC_1		('M')

#define WLZEFF_BMP_WIN_NEW		(40)

#define WLZEFF_BMP_CMP_RGB		(0)
#define WLZEFF_BMP_CMP_RLE8		(1)
#define WLZEFF_BMP_CMP_ESC		(0)
#define WLZEFF_BMP_CMP_EOL		(0)
#define WLZEFF_BMP_CMP_EOB		(1)
#define WLZEFF_BMP_CMP_DELTA		(2)

typedef unsigned char	WLZEFF_BMP_BYTE;

typedef short		WLZEFF_BMP_WORD;

typedef	unsigned short	WLZEFF_BMP_UINT;

typedef int		WLZEFF_BMP_DWORD;

typedef int		WLZEFF_BMP_LONG;

typedef struct _WlzEffBmpFileHead
{
  WLZEFF_BMP_UINT	bfType;
  WLZEFF_BMP_DWORD 	bfSize;
  WLZEFF_BMP_UINT	bfReserved1;
  WLZEFF_BMP_UINT	bfReserved2;
  WLZEFF_BMP_DWORD	bfOffBits;
} WlzEffBmpFileHead;

typedef struct _WlzEffBmpInfoHead
{
  WLZEFF_BMP_DWORD	biSize;
  WLZEFF_BMP_LONG	biWidth;
  WLZEFF_BMP_LONG	biHeight;
  WLZEFF_BMP_WORD	biPlanes;
  WLZEFF_BMP_WORD	biBitCount;
  WLZEFF_BMP_DWORD	biCompression;
  WLZEFF_BMP_DWORD	biSizeImage;
  WLZEFF_BMP_LONG	biXPelsPerMeter;
  WLZEFF_BMP_LONG	biYPelsPerMeter;
  WLZEFF_BMP_DWORD	biClrUsed;
  WLZEFF_BMP_DWORD	biClrImportant;
} WlzEffBmpInfoHead;

typedef struct _WlzEffBmpRGBQuad
{
  WLZEFF_BMP_BYTE	rgbBlue;
  WLZEFF_BMP_BYTE	rgbGreen;
  WLZEFF_BMP_BYTE	rgbRed;
  WLZEFF_BMP_BYTE	rgbReserved;
} WlzEffBmpRGBQuad;


#define WLZEFF_DEN_VERSION		(1)

typedef struct _WlzEffDenHeader
{
  short		version;
  short		orgMin[3];		/* Dimensions of original data file */
  short		orgMax[3];
  short		orgLen[3];
  short		extrMin[3];	        /* Extracted portion of orig file */
  short		extrMax[3];
  short		extrLen[3];
  short		mapMin[3];	  	/* Dimensions of this map */
  short		mapMax[3];
  short		mapLen[3];
  short		mapWarps;	 	/* Number of warps since extraction */
  unsigned int	mapLength;		/* Total number of densities in map */
} WlzEffDenHeader;

#define WLZEFF_ICS_REC_LEN_MAX		(256)
#define WLZEFF_ICS_VERSION_MAJOR	(1)
#define WLZEFF_ICS_VERSION_MINOR	(0)
#define WLZEFF_ICS_PARAMETERS_MAX	(8)

typedef enum _WlzEffIcsToken
{
  WLZEFF_ICS_TKN_NONE =			(0),
  WLZEFF_ICS_TKN_BITS,
  WLZEFF_ICS_TKN_COORDS,
  WLZEFF_ICS_TKN_FILENAME,
  WLZEFF_ICS_TKN_FLOAT,
  WLZEFF_ICS_TKN_FORMAT,
  WLZEFF_ICS_TKN_G3D,
  WLZEFF_ICS_TKN_INT,
  WLZEFF_ICS_TKN_LAYOUT,
  WLZEFF_ICS_TKN_L3D,
  WLZEFF_ICS_TKN_ORDER,
  WLZEFF_ICS_TKN_PARAMETERS,
  WLZEFF_ICS_TKN_REPRESENTATION,
  WLZEFF_ICS_TKN_SCIL,
  WLZEFF_ICS_TKN_SIGBITS,
  WLZEFF_ICS_TKN_SIGN,
  WLZEFF_ICS_TKN_SIGNED,
  WLZEFF_ICS_TKN_SIZES,
  WLZEFF_ICS_TKN_UNSIGNED,
  WLZEFF_ICS_TKN_VERSION,
  WLZEFF_ICS_TKN_VIDEO,
  WLZEFF_ICS_TKN_X,
  WLZEFF_ICS_TKN_Y,
  WLZEFF_ICS_TKN_Z
} WlzEffIcsToken;

typedef struct _WlzEffIcsHeader
{
  int		versionMajor;
  int		versionMinor;
  int		parameters;
  int		sizes[WLZEFF_ICS_PARAMETERS_MAX];
  int		sigBits;
  WlzEffIcsToken coords;
  WlzEffIcsToken order[WLZEFF_ICS_PARAMETERS_MAX];
  WlzEffIcsToken format;
  WlzEffIcsToken sign;
  WlzEffIcsToken scil;
} WlzEffIcsHeader;

#define WLZEFF_PIC_MAGIC		(12345)  /* BioRad .pic magic number */
#define WLZEFF_PIC_HEADBYTES		(76)

#define WLZEFF_PIC_OFF_NX		(0)
#define WLZEFF_PIC_OFF_NY		(2)
#define WLZEFF_PIC_OFF_NPIC		(4)
#define WLZEFF_PIC_OFF_BLACKVAL		(6)
#define WLZEFF_PIC_OFF_WHITEVAL		(8)
#define WLZEFF_PIC_OFF_BYTEPIX		(14)
#define WLZEFF_PIC_OFF_MERGED		(50)
#define WLZEFF_PIC_OFF_RGB		(52)
#define WLZEFF_PIC_OFF_MAGIC		(54)
#define WLZEFF_PIC_OFF_BLACKVALMERGE 	(56)
#define WLZEFF_PIC_OFF_WHITEVALMERGE 	(58)
#define WLZEFF_PIC_OFF_RGBMERGED	(60)
#define WLZEFF_PIC_OFF_LENSPOWER	(64)
#define WLZEFF_PIC_OFF_MAGFACTOR	(66)

#define WLZEFF_PIC_HEAD_WORD_SET(H,O,V) \
{ \
  *((unsigned char *)H + O) = (unsigned char )(V & 0xff); \
  *((unsigned char *)H + O + 1) = (unsigned char )((V >> 8) & 0xff); \
}

#define WLZEFF_PIC_HEAD_WORD_GET(V,O,H) \
{ \
  V = *((unsigned char *)H + O + 1); \
  V <<= 8; \
  V |= *((unsigned char *)H + O); \
}

typedef struct _WlzEffPicHeader
{
  short		nX;
  short		nY;
  short		nPic;
  short		blackVal;
  short		whiteVal;
  short		nBytes;		  /* 1 for 8 bit pixels, 0 for 16 bit pixels */
  short		merged;
  short		colSelect;
  short		magic;
  short		blackValMerged;
  short		whiteValMerged;
  short		colSelectMerged;
  short		lens;
  float		magFactor;
} WlzEffPicHeader;

#define WLZEFF_PGM_MAGIC		"P5"

typedef enum _WlzEffPnmType
{
  WLZEFF_PNM_TYPE_NONE =		(0),
  WLZEFF_PNM_TYPE_PBM_ASC,
  WLZEFF_PNM_TYPE_PBM_BIN,
  WLZEFF_PNM_TYPE_PGM_ASC,
  WLZEFF_PNM_TYPE_PGM_BIN,
  WLZEFF_PNM_TYPE_PPM_ASC,
  WLZEFF_PNM_TYPE_PPM_BIN
} WlzEffPnmType;

typedef struct _WlzEffStackCtrHeader
{
  WlzIVertex3	volOrigin;
  WlzIVertex3	volSize;
  WlzFVertex3	voxSize;
} WlzEffStackCtrHeader;

#define WLZEFF_STACK_NAMEDIGITS		(8)

#define WLZEFF_STACK_CTR_IDENT		"ident"
#define WLZEFF_STACK_CTR_IDENTSTR	"WLZSTACKCTR"
#define WLZEFF_STACK_CTR_IDENTVERSION	"1.0"
#define WLZEFF_STACK_CTR_VOLORIGIN	"volume origin"
#define WLZEFF_STACK_CTR_VOLSIZE	"volume size"
#define WLZEFF_STACK_CTR_VOXSIZE	"voxel size"
#define WLZEFF_STACK_CTR_FILES		"files"
#define WLZEFF_STACK_CTR_COMMENT	"#"
#define WLZEFF_STACK_CTR_FIELDSEP	":"
#define WLZEFF_STACK_CTR_RECORDMAX	(1024)

#define WLZEFF_SLC_MAGIC		(11111)

typedef enum _WlzEffSlcDataUnits
{
  WLZEFF_SLC_DATAUNITS_METER,
  WLZEFF_SLC_DATAUNITS_MILLIMETER,
  WLZEFF_SLC_DATAUNITS_MICROMETER,
  WLZEFF_SLC_DATAUNITS_FOOT
} WlzEffSlcDataUnits;

typedef enum _WlzEffSlcDataSource
{
  WLZEFF_SLC_DATASRC_BIORAD,
  WLZEFF_SLC_DATASRC_MRI,
  WLZEFF_SLC_DATASRC_CT,
  WLZEFF_SLC_DATASRC_SIM,
  WLZEFF_SLC_DATASRC_BINVOX,
  WLZEFF_SLC_DATASRC_FUZVOX,
  WLZEFF_SLC_DATASRC_OTHER
} WlzEffSlcDataSource;

typedef enum _WlzEffSlcDataMod
{
  WLZEFF_SLC_DATAMOD_ORIGINAL,
  WLZEFF_SLC_DATAMOD_RESAMPLED,
  WLZEFF_SLC_DATAMOD_FILTERED,
  WLZEFF_SLC_DATAMOD_RESANDFILT,
  WLZEFF_SLC_DATAMOD_OTHER
} WlzEffSlcDataMod;

typedef struct _WlzEffSlcHeader
{
  int		magic;
  WlzIVertex3	size;
  int		bits;
  WlzFVertex3	spacing;
  int		units;
  int		source;
  int		modification;
  int		compression;
  WlzIVertex2	iconSize;
  unsigned char	*icon; 		/* Three components R, G, B each of iconSize */
} WlzEffSlcHeader;

#define WLZEFF_VFF_REC_LEN_MAX		(256)

typedef enum _WlzEffVffRecord
{
  WLZEFF_VFF_REC_NONE =			(0),
  WLZEFF_VFF_REC_NCAA,
  WLZEFF_VFF_REC_TYPE,
  WLZEFF_VFF_REC_FORMAT,
  WLZEFF_VFF_REC_RANK,
  WLZEFF_VFF_REC_BANDS,
  WLZEFF_VFF_REC_BITS,
  WLZEFF_VFF_REC_RAWSIZE,
  WLZEFF_VFF_REC_SIZE,
  WLZEFF_VFF_REC_ORIGIN,
  WLZEFF_VFF_REC_EXTENT,
  WLZEFF_VFF_REC_ASPECT
} WlzEffVffRecord;

typedef enum _WlzEffVffFormat
{
  WLZEFF_VFF_FORMAT_NONE =		(0),
  WLZEFF_VFF_FORMAT_BASE,
  WLZEFF_VFF_FORMAT_SLICE
} WlzEffVffFormat;

typedef enum _WlzEffVffType
{
  WLZEFF_VFF_TYPE_NONE =		(0),
  WLZEFF_VFF_TYPE_CONNECTIVITY,
  WLZEFF_VFF_TYPE_INCLUDE,
  WLZEFF_VFF_TYPE_NURBPATCH,
  WLZEFF_VFF_TYPE_RASTER,
  WLZEFF_VFF_TYPE_VERTICIES
} WlzEffVffType;

typedef struct _WlzEffVffHeader
{
  int		ncaa;
  WlzEffVffType   type;
  WlzEffVffFormat format;
  int		rank;
  int		bands;
  int		bits;
  int		rawsize;
  WlzIVertex3	size;
  WlzDVertex3	origin;
  WlzDVertex3	extent;
  WlzDVertex3	aspect;
} WlzEffVffHeader;

#define	WLZEFF_VTK_VERSION_MAJOR	(1)
#define	WLZEFF_VTK_VERSION_MINOR	(0)

typedef enum _WlzEffVtkDataType
{
  WLZEFF_VTK_DATATYPE_ASCII,
  WLZEFF_VTK_DATATYPE_BINARY
} WlzEffVtkDataType;

typedef enum _WlzEffVtkType
{
  WLZEFF_VTK_TYPE_STRUCTURED_POINTS,
  WLZEFF_VTK_TYPE_STRUCTURED_GRID,
  WLZEFF_VTK_TYPE_UNSTRUCTURED_GRID,
  WLZEFF_VTK_TYPE_POLYDATA,
  WLZEFF_VTK_TYPE_RECTILNEAR_GRID
} WlzEffVtkType;

typedef enum _WlzEffVtkPolyDataType
{
  WLZEFF_VTK_POLYDATATYPE_POINTS,
  WLZEFF_VTK_POLYDATATYPE_VERTICIES,
  WLZEFF_VTK_POLYDATATYPE_LINES,
  WLZEFF_VTK_POLYDATATYPE_POLYGONS,
  WLZEFF_VTK_POLYDATATYPE_TRIANGLE_STRIPS
} WlzEffVtkPolyDataType;

typedef enum _WlzEffVtkUnstructuredGridType
{
  WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_CELLS,
  WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_CELL_TYPES,
  WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_LOOKUP_TABLE,
  WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_POINT_DATA,
  WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_POINTS,
  WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_SCALARS,
  WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_VECTORS
} WlzEffVtkUnstructuredGridType;

typedef struct _WlzEffVtkHeader
{
  int		versionMajor;
  int		versionMinor;
  char		title[256];
  WlzEffVtkDataType dataType;
  WlzEffVtkType type;
} WlzEffVtkHeader;

/* IPLab format */

typedef enum _WlzEffIPLType
{
  WLZEFF_IPL_TYPE_UBYTE		= 0,
  WLZEFF_IPL_TYPE_SHORT		= 1,
  WLZEFF_IPL_TYPE_INT		= 2,
  WLZEFF_IPL_TYPE_FLOAT		= 3,
  WLZEFF_IPL_TYPE_COL_16	= 4,
  WLZEFF_IPL_TYPE_COL_24	= 5,
  WLZEFF_IPL_TYPE_U_16		= 6,
  WLZEFF_IPL_TYPE_LAST
} WlzEffIPLType;

typedef void * WlzIPLCSpecArray;

typedef struct _WlzEffIPLHeader
{
  char		version[5];
  unsigned char	format;
  WlzEffIPLType	dType;
  int		nWidth;
  int		nHeight;
  int		nFrames;
  char		fileClutID;
  int		overlayInFile;
  char		viewMode;
  double	delta;
  char		units[11];
  char		normType;
  char		normSource;
  int		numRegNarks;
  double	normMin;
  double	normMax;
  WlzIPLCSpecArray	colorTable;
}WlzEffIPLHeader;

/*!
* \enum		_WlzEffAmToken
* \ingroup	WlzExtFF
* \brief	Tokens for parsing the headers of Amira lattice files.
*		Typedef: ::WlzEffAmToken
*/
typedef enum _WlzEffAmToken
{
  WLZEFF_AM_TOKEN_NONE		= (0),
  WLZEFF_AM_TOKEN_BOUNDINGBOX,
  WLZEFF_AM_TOKEN_CLOSE,
  WLZEFF_AM_TOKEN_COLOR,
  WLZEFF_AM_TOKEN_COLORMAP,
  WLZEFF_AM_TOKEN_CONTENT,
  WLZEFF_AM_TOKEN_COORDTYPE,
  WLZEFF_AM_TOKEN_DEFINE,
  WLZEFF_AM_TOKEN_EXPRESSION,
  WLZEFF_AM_TOKEN_HASH,
  WLZEFF_AM_TOKEN_ID,
  WLZEFF_AM_TOKEN_IMAGEDATA,
  WLZEFF_AM_TOKEN_LATTICE,
  WLZEFF_AM_TOKEN_MATERIALS,
  WLZEFF_AM_TOKEN_NAME,
  WLZEFF_AM_TOKEN_OPEN,
  WLZEFF_AM_TOKEN_PARAMETERS,
  WLZEFF_AM_TOKEN_SEEDS,
  WLZEFF_AM_TOKEN_LIMITS,
  WLZEFF_AM_TOKEN_TIFF,
  WLZEFF_AM_TOKEN_TRANSFORM
} WlzEffAmToken;

/*!
* \enum		_WlzEffAmDim
* \ingroup	WlzExtFF
* \brief	Dimension of the data.
*		Typedef: ::WlzEffAmDim
*/
typedef	enum _WlzEffAmDim
{
  WLZEFF_AM_DIM_NONE 		= (0),	/*!< Dimension unknown. */
  WLZEFF_AM_DIM_2		= (2),	/*!< 2D. */
  WLZEFF_AM_DIM_3		= (3)	/*!< 3D. */
} WlzEffAmDim;

/*!
* \enum		_WlzEffAmEndian
* \ingroup	WlzExtFF
* \brief	Big or little endian binary data. The values endian enum
* 		must be kept distinct from the dimension enum values
* 		(apart from none).
*		Typedef: ::WlzEffAmFormat
*/
typedef enum _WlzEffAmEndian
{
  WLZEFF_AM_ENDIAN_NONE		= (0),	/*!< Unknown endianness. */
  WLZEFF_AM_ENDIAN_BIG		= (10),	/*!< Big endian. */
  WLZEFF_AM_ENDIAN_LITTLE	= (11)	/*!< Little endian. */
} WlzEffAmEndian;

/*!
* \enum		_WlzEffAmFormat
* \ingroup	WlzExtFF
* \brief	ASCII or binary data.
*		Typedef: ::WlzEffAmFormat
*/
typedef enum _WlzEffAmFormat
{
  WLZEFF_AM_FMT_NONE		= (0),	/*!< Unknown data format. */
  WLZEFF_AM_FMT_BINARY		= (1),	/*!< Header ascii, Binary data. */
  WLZEFF_AM_FMT_ASCII		= (2)	/*!< Ascii header and data. */
} WlzEffAmFormat;

/*!
* \enum		_WlzEffAmDatType
* \ingroup	WlzExtFF
* \brief	Type of data: byte, ...
*		Typedef: ::WlzEffAmDatType
*/
typedef enum _WlzEffAmDatType
{
  WLZEFF_AM_DATTYPE_NONE	= (0),	/*!< Unknown data type. */
  WLZEFF_AM_DATTYPE_BYTE	= (1),	/*!< Byte (8 bit) data. */
  WLZEFF_AM_DATTYPE_SHORT	= (2)	/*!< Short (16 bit) data. */
} WlzEffAmDatType;

/*!
* \enum		_WlzEffAmCoordType
* \ingroup      WlzExtFF
* \brief        Type of coordinate system.
*		Typedef: ::WlzEffAmCoordType
*/
typedef enum _WlzEffAmCoordType
{
  WLZEFF_AM_COORDTYPE_NONE	= (0),	/*!< Unknown coordinate type. */
  WLZEFF_AM_COORDTYPE_UNITFORM	= (1)	/*!< Uniform coordinates. */
} WlzEffAmCoordType;

/*!
* \enum		_WlzEffAmLatComp
* \ingroup      WlzExtFF
* \brief        Type of compression used.
*		Typedef: ::WlzEffAmLatComp
*/
typedef enum _WlzEffAmLatComp
{
  WLZEFF_AM_LATCOMP_NONE	= (0),  /*!< No compression. */
  WLZEFF_AM_LATCOMP_HXBYTERLE	= (1)	/*!< Run length encoded bytes. */
} WlzEffAmLatComp;

/*!
* \enum		_WlzEffAmLatType
* \ingroup	WlzExtFF
* \brief	Type of lattice: uniform, ...
*		Typedef: ::WlzEffAmLatType
*/
typedef enum	_WlzEffAmLatType
{
  WLZEFF_AM_LATTYPE_NONE	= (0),  /*!< Unknown lattice type. */
  WLZEFF_AM_LATTYPE_DATA	= (1),	/*!< Voxel lattice data. */
  WLZEFF_AM_LATTYPE_LABELS	= (2)	/*!< Domain lattice data. */
} WlzEffAmLatType;

/*!
* \struct	_WlzEffAmMaterial
* \ingroup	WlzExtFF
* \brief	Item in an Amira material list.
*		Typedef: ::WlzEffAmMaterial
*/
typedef struct _WlzEffAmMaterial
{
  int			id;		/*!< Index in lattice labels. */
  double		color[3];	/*!< RGB colour components. */
  char			*name;		/*!< Material name. Should be free'd
  					    using AlcFree(). */
  struct _WlzEffAmMaterial *next;	/*!< Next material in list. */
  struct _WlzEffAmMaterial *prev;	/*!< Previous material in list */
} WlzEffAmMaterial;

/*!
* \struct	_WlzEffAmHead
* \ingroup	WlzExtFF
* \brief	Head of a list of Amira materials.
*		Typedef: ::WlzEffAmHead
*/
typedef struct _WlzEffAmHead
{
  int			versionMajor;
  int			versionMinor;
  WlzEffAmDim 		dim; 		/*!< Dimension of the data. */
  WlzEffAmFormat	fmt;		/*!< Data format. */
  WlzEffAmEndian	endian;		/*!< Whether the binary data in a file
  					     is big or little endian? */
  WlzEffAmLatType	latType;	/*!< Lattice type. */
  WlzEffAmDatType	latDatType;	/*!< Lattice data type. */
  WlzEffAmCoordType	coordType;	/*!< Coordinate type. */
  WlzDBox3		bBox;		/*!< Real world coordinates of the
  					     bounding box. */
  WlzIVertex3		latSize; 	/*!< Lattice size. */
  int			latBytes;	/*!< Number of bytes to read. */
  int			latLabel;	/*!< Label for lattice. */
  WlzEffAmLatComp	latComp;	/*!< Lattice compression. */
  int			matCount;	/*!< Number of materials. */
  WlzEffAmMaterial	*materials;	/*!< Linked list of materials, with
  					     the first item in the list having
					     a NULL 'prev' entry and the last
					     having a NULL 'next' entry. */
  char			*imageData;	/*!< Associated image file. */
} WlzEffAmHead;

/*!
* \struct	_WlzEffBibWarpInputThresholdParamsStruct
* \brief	Bibfile threshold parameters record.
*		Typedef: ::WlzEffBibWarpInputThresholdParamsStruct
*/
typedef struct _WlzEffBibWarpInputThresholdParamsStruct
{
  WlzRGBAThresholdType	thresholdType;	/*!< Threshold type - single, multi etc. */
  WlzRGBAColorSpace	threshRGBSpace;	/*!< Colour space */
  WlzRGBAColorChannel	threshColorChannel; /*!< Colour channel */
  int		threshRangeLow;		/*!< Single channel low threshold value */
  int		threshRangeHigh;	/*!< Single channel high threshold value */
  int		threshRangeRGBLow[3];	/*!< Multi channel low threshold value */
  int		threshRangeRGBHigh[3];	/*!< Multi channel high threshold value */
  WlzUInt	threshRGBCombination;	/*!< Colour combination logic mask */
  WlzPixelV	lowRGBPoint;	/*!< Low-point for slice/box/ball threshold */
  WlzPixelV	highRGBPoint;	/*!< High-point for slice/box/ball threshold */
  double	colorEllipseEcc;	/*!< Ball eccentricity */
  int		globalThreshFlg;	/*!< Global thresholding flag */
  WlzIVertex2	globalThreshVtx;	/*!< Global threshold vertex */
  int		incrThreshFlg;		/*!< Incremental threshold flag */
  int		pickThreshFlg; 		/*!< Pick mode (endpoint values) flag */
  int		distanceThreshFlg;	/*!< Distance mode flag */

} WlzEffBibWarpInputThresholdParamsStruct;

/* Analyze 7.5 file format types. */

/*!
* \enum		_WlzEffAnlDataType
* \ingroup	WlzExtFF
* \brief	Valid Analyze data type values.
*		Typedef: ::WlzEffAnlDataType
*/
typedef enum _WlzEffAnlDataType
{
  WLZEFF_ANL_DT_NONE = 		0,
  WLZEFF_ANL_DT_UNKOWN =	0,
  WLZEFF_ANL_DT_BINARY =	1,
  WLZEFF_ANL_DT_UNSIGNED_CHAR =	2,
  WLZEFF_ANL_DT_SIGNED_SHORT =	4,
  WLZEFF_ANL_DT_SIGNED_INT =	8,
  WLZEFF_ANL_DT_FLOAT =		16,
  WLZEFF_ANL_DT_COMPLEX =	32,
  WLZEFF_ANL_DT_DOUBLE =	64,
  WLZEFF_ANL_DT_RGB =		128,
  WLZEFF_ANL_DT_ALL =		255
} WlzEffAnlDataType;

/*!
* \enum		_WlzEffAnlDataType
* \ingroup	WlzExtFF
* \brief	Valid Analyze data type values.
*		Typedef: ::WlzEffAnlOrient
*/
typedef enum _WlzEffAnlOrient
{
  WLZEFF_ANL_ORIENT_TU = 	0, 	/*!< Transverse unflipped. */
  WLZEFF_ANL_ORIENT_CU = 	1, 	/*!< Coronal unflipped. */
  WLZEFF_ANL_ORIENT_SU = 	2, 	/*!< Sagital unflipped. */
  WLZEFF_ANL_ORIENT_TF = 	3, 	/*!< Transverse flipped. */
  WLZEFF_ANL_ORIENT_CF = 	4, 	/*!< Coronal flipped. */
  WLZEFF_ANL_ORIENT_SF = 	5 	/*!< Sagital flipped. */
} WlzEffAnlOrient;

/*!
* \struct	_WlzEffAnlComplex
* \ingroup	WlzExtFF
* \brief	Analyze 7.5 complex number representation.
*		Typedef: ::WlzEffAnlComplex
*/
typedef struct _WlzEffAnlComplex
{
  float		real;
  float		imag;
} WlzEffAnlComplex;

/*!
* \struct	_WlzEffAnlHeaderKey
* \ingroup	WlzExtFF
* \brief	Analyze 7.5 file header key.
*		Typedef: ::WlzEffAnlFileKey
*/
typedef struct _WlzEffAnlHeaderKey
{
  int		hdrSz;			/*!< Size of header in bytes: 4
  					    bytes. */
  char		dataType[10];		/*!< Data type: 10 bytes. */
  char		dbName[18];		/*!< 18 bytes. */
  int		extents;		/*!< Should be 16384: 4 bytes. */
  short		sessionErr;		/*!< Session error: 2 bytes. */
  char		regular;		/*!< Must be 'r' to indicate that the
  					    images are the same size: 1
					    byte. */
  char		hKeyUn0;		/*!< 1 byte.
					    Total = 40 bytes. */
} WlzEffAnlFileKey;

/*!
* \struct	_WlzEffAnlImageDim
* \ingroup	WlzExtFF
* \brief	Analyze 7.5 file header structure for the image dimensions.
*		Typedef: ::WlzEffAnlImageDim
*/
typedef struct _WlzEffAnlImageDim
{
  short		dim[8];			/*!< Array of the image dimensions.
					    The number of dimensions is
					    usually 4 with;
					    [0] = number of dimensions,
                                            [1] = number of columns (X),
					    [2] = number of lines (Y),
					    [3] = number of planes (Z),
					    [4] = number of time points.
					    16 bytes. */
  short		unused8;		/*!< Unused. 2 bytes. */
  short		unused9;		/*!< Unused. 2 bytes. */
  short		unused10;		/*!< Unused. 2 bytes. */
  short		unused11;		/*!< Unused. 2 bytes. */
  short		unused12;		/*!< Unused. 2 bytes. */
  short		unused13;		/*!< Unused. 2 bytes. */
  short		unused14;		/*!< Unused. 2 bytes. */
  short		dataType;		/*!< The data type for the image, which
  					    must be a member of
  					    enum::_WlzEffAnlDataType.
					    2 bytes. */
  short		bitPix;			/*!< Number of bits per pixel, which
  					    must be one of; 1, 8, 16, 32 or
					    64. 2 bytes. */
  short		dimUn0;			/*!< 2 bytes. */
  float		pixdim[8];		/*!< Real world pixel dimensions in mm
  					    and ms. The number of dimensions is
					    usually 4 with;
					    [0] = number of dimensions,
                                            [1] = x size (pixel width),
					    [2] = y size (pixel height),
					    [3] = z size (voxel depth),
					    [4] = time interval.
					    32 bytes.*/
  float		voxOffset;		/*!< The byte offset in the ".img" file
  					    at which the pixels start.
					    If negative then the absolute value
					    applies to every image in the file.
					    4 bytes. */
  float		fUnused1;		/*!< 4 bytes. */
  float		fUnused2;		/*!< 4 bytes. */
  float		fUnused3;		/*!< 4 bytes. */
  float		calMax;			/*!< Maximum calibratin value.
             				    4 bytes. */
  float		calMin;			/*!< Minimum calibratin value.
  					    4 bytes. */
  float		compressed;		/*!< 4 bytes. */
  float		verified;		/*!< 4 bytes. */
  int		glMax;			/*!< Maximum pixel value. 4 bytes. */
  int		glMin;			/*!< Minimum pixel value. 4 bytes.
  					    Total = 108 bytes.*/
} WlzEffAnlImageDim;


/*!
* \struct	_WlzEffAnlDataHistory
* \ingroup	WlzExtFF
* \brief	Analyze 7.5 file header structure for the image history.
*		Typedef: ::WlzEffAnlDataHistory
*/
typedef struct _WlzEffAnlDataHistory
{
  char		descrip[80];		/*!< Description: 80 bytes. */
  char		auxFile[24];		/*!< 24 bytes. */
  char		orient;			/*!< Slice orientation, which must
  					    be a member of
					    enum::_WlzEffAnlOrient.
					    This may be used by Analyze to
					    determine whether images should be
					    flipped before being displayed.
					    1 byte. */
  char		originator[10];		/*!< 10 bytes. */
  char		generated[10];		/*!< 10 bytes. */
  char		scanNum[10];		/*!< 10 bytes. */
  char		patientId[10];		/*!< 10 bytes. */
  char		expDate[10];		/*!< 10 bytes. */
  char		expTime[10];		/*!< 10 bytes. */
  char		hisUn0[3];		/*!< 3 bytes. */
  int		views;			/*!< 4 bytes. */
  int		volsAdded;		/*!< 4 bytes. */
  int		startField;		/*!< 4 bytes. */
  int		fieldSkip;		/*!< 4 bytes. */
  int		oMax;			/*!< 4 bytes. */
  int		oMin;			/*!< 4 bytes. */
  int		sMax;			/*!< 4 bytes. */
  int		sMin;			/*!< 4 bytes. */
} WlzEffAnlDataHistory;

/*!
* \struct	_WlzEffAnlDsr
* \ingroup	WlzExtFF
* \brief	Analyze 7.5 file header.
*		Typedef: ::WlzEffAnlDsr
*/
typedef struct _WlzEffAnlDsr
{
  struct _WlzEffAnlHeaderKey	hk;	/*!< 40 bytes. */
  struct _WlzEffAnlImageDim	dim;	/*!< 108 bytes. */
  struct _WlzEffAnlDataHistory	hist;	/*!< 200 bytes. */
} WlzEffAnlDsr;

#endif /* WLZ_EXT_BIND */

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* ! WLZ_EXTFFTYPE_H */
