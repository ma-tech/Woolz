#ifndef WLZ_EXTFFTYPE_H
#define WLZ_EXTFFTYPE_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFFType.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Header file with data types for external data file
*		format support for the MRC Human Genetics Unit Woolz
*		library.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef enum
{
  WLZEFF_FORMAT_NONE,
  WLZEFF_FORMAT_BMP,
  WLZEFF_FORMAT_DEN,
  WLZEFF_FORMAT_ICS,
  WLZEFF_FORMAT_PNM,
  WLZEFF_FORMAT_PIC,
  WLZEFF_FORMAT_SLC,
  WLZEFF_FORMAT_VFF,
  WLZEFF_FORMAT_VTK,
  WLZEFF_FORMAT_WLZ,
  WLZEFF_FORMAT_COUNT 			     /* Keep last: Number of formats */
} WlzEffFormat;

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

typedef struct
{
  WLZEFF_BMP_UINT	bfType;
  WLZEFF_BMP_DWORD 	bfSize;
  WLZEFF_BMP_UINT	bfReserved1;
  WLZEFF_BMP_UINT	bfReserved2;
  WLZEFF_BMP_DWORD	bfOffBits;
} WlzEffBmpFileHead;

typedef struct
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

typedef struct
{
  WLZEFF_BMP_BYTE	rgbBlue;
  WLZEFF_BMP_BYTE	rgbGreen;
  WLZEFF_BMP_BYTE	rgbRed;
  WLZEFF_BMP_BYTE	rgbReserved;
} WlzEffBmpRGBQuad;


#define WLZEFF_DEN_VERSION		(1)

typedef struct
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

typedef enum
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

typedef struct
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

typedef struct
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

typedef enum
{
  WLZEFF_PNM_TYPE_NONE =		(0),
  WLZEFF_PNM_TYPE_PBM_ASC,
  WLZEFF_PNM_TYPE_PBM_BIN,
  WLZEFF_PNM_TYPE_PGM_ASC,
  WLZEFF_PNM_TYPE_PGM_BIN,
  WLZEFF_PNM_TYPE_PPM_ASC,
  WLZEFF_PNM_TYPE_PPM_BIN
} WlzEffPnmType;

typedef struct
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

typedef enum
{
  WLZEFF_SLC_DATAUNITS_METER,
  WLZEFF_SLC_DATAUNITS_MILLIMETER,
  WLZEFF_SLC_DATAUNITS_MICROMETER,
  WLZEFF_SLC_DATAUNITS_FOOT
} WlzEffSlcDataUnits;

typedef enum
{
  WLZEFF_SLC_DATASRC_BIORAD,
  WLZEFF_SLC_DATASRC_MRI,
  WLZEFF_SLC_DATASRC_CT,
  WLZEFF_SLC_DATASRC_SIM,
  WLZEFF_SLC_DATASRC_BINVOX,
  WLZEFF_SLC_DATASRC_FUZVOX,
  WLZEFF_SLC_DATASRC_OTHER
} WlzEffSlcDataSource;

typedef enum
{
  WLZEFF_SLC_DATAMOD_ORIGINAL,
  WLZEFF_SLC_DATAMOD_RESAMPLED,
  WLZEFF_SLC_DATAMOD_FILTERED,
  WLZEFF_SLC_DATAMOD_RESANDFILT,
  WLZEFF_SLC_DATAMOD_OTHER
} WlzEffSlcDataMod;

typedef struct
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

typedef enum
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

typedef enum
{
  WLZEFF_VFF_FORMAT_NONE =		(0),
  WLZEFF_VFF_FORMAT_BASE,
  WLZEFF_VFF_FORMAT_SLICE
} WlzEffVffFormat;

typedef enum
{
  WLZEFF_VFF_TYPE_NONE =		(0),
  WLZEFF_VFF_TYPE_CONNECTIVITY,
  WLZEFF_VFF_TYPE_INCLUDE,
  WLZEFF_VFF_TYPE_NURBPATCH,
  WLZEFF_VFF_TYPE_RASTER,
  WLZEFF_VFF_TYPE_VERTICIES
} WlzEffVffType;

typedef struct
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

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* ! WLZ_EXTFFTYPE_H */
