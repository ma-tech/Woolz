#ifndef RECONSTRUCTTYPE_H
#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructType.h				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Header file with type definitions for the MRC Human
*		Genetics Unit reconstruction library.		
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#define RECONSTRUCTTYPE_H

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef enum
{
  REC_ERR_NONE		= (0),
  REC_ERR_USAGE,			       /* Command line usage invalid */
  REC_ERR_ARGS,				   /* Command line arguments invalid */
  REC_ERR_READ,						     /* Read failure */
  REC_ERR_WRITE,					    /* Write failure */
  REC_ERR_MALLOC,				/* Memory allocation failure */
  REC_ERR_UNIMPL,				    /* Unimplemented feature */
  REC_ERR_WLZ,						      /* Woolz error */
  REC_ERR_SYNTAX,					/* File syntax error */
  REC_ERR_FUNC,				      /* Function parameters invalid */
  REC_ERR_CANCEL,	   /* Not an error as such, more of a cancel request */
  REC_ERR_LIST,					       /* Section list error */
  REC_ERR_MAX
} RecError;

typedef enum
{
  REC_FILE_NONE		= (0),
  REC_FILE_COMPND,
  REC_FILE_SECLIST,
  REC_FILE_3D,
  REC_FILE_MAX
} RecFileType;

typedef enum
{
  REC_MTHD_NONE		= (0),
  REC_MTHD_PRINC	= (1),       /* Principle axes (with centre of mass) */
  REC_MTHD_TRANS	= (1<<1),			      /* Translation */
  REC_MTHD_ROTATE	= (1<<2),				 /* Rotation */
  REC_MTHD_IDENTITY	= (1<<3)	/* Transform initialised to identity */
} RecMethod;

typedef enum
{
  REC_PP_NONE	 	= (0),
  REC_PP_BACKGROUND	= (1),			   /* Force background value */
  REC_PP_SOBEL	 	= (1<<1),	             /* Sobel edge detection */
  REC_PP_ERODE		= (1<<2),			   /* Domain erosion */
  REC_PP_INVERT		= (1<<3),		     /* Gray value inversion */
  REC_PP_LAPLAC		= (1<<4),	       /* Laplacian edge enhancement */
  REC_PP_NOISE		= (1<<5),	/* Noise used to fill outside domain */ 
  REC_PP_THRESH	 	= (1<<6),   /* Threshold (for above) mean gray value */
  REC_PP_WINDOW	 	= (1<<7),	   /* Spatial domain window function */
  REC_PP_SAMPLE	 	= (1<<8) 	     /* Sample prior to registration */
} RecPPMethod;

typedef enum
{
  REC_DBG_NONE  	= (0),
  REC_DBG_LVL_1  	= (1),
  REC_DBG_LVL_2  	= (1<<1),
  REC_DBG_LVL_3  	= (1<<2),
  REC_DBG_LVL_FN  	= (1<<3),
  REC_DBG_3D		= (1<<4),
  REC_DBG_AUTO		= (1<<5),
  REC_DBG_CROSS		= (1<<6),
  REC_DBG_DBL2		= (1<<7),
  REC_DBG_FILE		= (1<<8),
  REC_DBG_FOUR		= (1<<9),
  REC_DBG_MISC		= (1<<10),
  REC_DBG_WLZ		= (1<<11),
  REC_DBG_PPROC		= (1<<12),
  REC_DBG_PRINC		= (1<<13),
  REC_DBG_REG		= (1<<14),
  REC_DBG_ROT		= (1<<15),
  REC_DBG_TRAN		= (1<<16),
  REC_DBG_SEC		= (1<<17)
} RecDbgMask;

typedef enum
{
  REC_SECMSK_NONE 		= 0,
  REC_SECMSK_INDEX		= (1),
  REC_SECMSK_ITERATIONS 	= (1<<1),
  REC_SECMSK_CORREL		= (1<<2),
  REC_SECMSK_IMAGEFILE		= (1<<3),
  REC_SECMSK_TRANSF_TX		= (1<<4),
  REC_SECMSK_TRANSF_TY		= (1<<5),
  REC_SECMSK_TRANSF_TZ		= (1<<6),
  REC_SECMSK_TRANSF_SCALE	= (1<<7),
  REC_SECMSK_TRANSF_THETA	= (1<<8),
  REC_SECMSK_TRANSF_PHI		= (1<<9),
  REC_SECMSK_TRANSF_ALPHA	= (1<<10),
  REC_SECMSK_TRANSF_PSI		= (1<<11),
  REC_SECMSK_TRANSF_XSI		= (1<<12),
  REC_SECMSK_TRANSF_INVERT	= (1<<13),
  REC_SECMSK_ALL		= 0x0003fff   /* Mask for all section fields */
} RecSecMask;

typedef enum
{
  REC_CCFLAG_NONE 	= (0),
  REC_CCFLAG_DATA0VALID	= (1),
  REC_CCFLAG_DATA1VALID	= (1<<1)
} RecCcFlag;

typedef enum
{
  REC_TRMODE_REL,
  REC_TRMODE_ABS
} RecTransformMode;

/* Constraints on image size for Fourier transforms */
#define REC_FOUR_DIMP2_MIN	(3)
#define REC_FOUR_DIMP2_MAX	(12)
#define REC_FOUR_DIM_MIN	(1<<(REC_FOUR_DIMP2_MIN))
#define REC_FOUR_DIM_MAX	(1<<(REC_FOUR_DIMP2_MAX))

#ifdef REC_THREADS_USED
#define REC_THREADS_MAX	(8)
#else /* ! REC_THREADS_USED */
#define REC_THREADS_MAX	(1)
#endif /* REC_THREADS_USED */

#define REC_MAX_WINSZ	(200)
#define REC_MIN_WINSZ	(10)
#define	REC_MAX_ITLIM	(1000)
#define	REC_MIN_XLIM	(4)
#define REC_MAX_XLIM	(512)
#define	REC_MIN_YLIM	(4)
#define REC_MAX_YLIM	(512)
#define	REC_MIN_RLIM	(1.0)
#define REC_MAX_RLIM	(180.0)
#define REC_MAX_RECORD	(256)
#define REC_MAX_LIST	(1000)

#define REC_TRANS_TX_MAX (1000.0)
#define REC_TRANS_TX_MIN (-(REC_TRANS_TX_MAX))
#define REC_TRANS_TY_MAX (1000.0)
#define REC_TRANS_TY_MIN (-(REC_TRANS_TY_MAX))
#define REC_TRANS_TZ_MAX (0.0)
#define REC_TRANS_TZ_MIN (-(REC_TRANS_TZ_MAX))
#define REC_TRANS_SCALE_MAX (100.0)
#define REC_TRANS_SCALE_MIN (1/(REC_TRANS_SCALE_MAX))
#define REC_TRANS_THETA_MAX (WLZ_M_2_PI)
#define REC_TRANS_THETA_MIN (-(REC_TRANS_THETA_MAX))
#define REC_TRANS_PHI_MAX (WLZ_M_2_PI)
#define REC_TRANS_PHI_MIN (-(REC_TRANS_PHI_MAX))
#define REC_TRANS_ALPHA_MAX (WLZ_M_2_PI)
#define REC_TRANS_ALPHA_MIN (-(REC_TRANS_ALPHA_MAX))
#define REC_TRANS_PSI_MAX (WLZ_M_2_PI)
#define REC_TRANS_PSI_MIN (-(REC_TRANS_PSI_MAX))
#define REC_TRANS_XSI_MAX (WLZ_M_2_PI)
#define REC_TRANS_XSI_MIN (-(REC_TRANS_XSI_MAX))
#define REC_TRANS_INVERT_MAX (1)
#define REC_TRANS_INVERT_MIN (0)

#define REC_DEF_PP	((RecPPMethod )(REC_PP_BACKGROUND|REC_PP_INVERT| \
					REC_PP_THRESH))
#define REC_DEF_METHOD	((RecMethod )(REC_MTHD_PRINC|REC_MTHD_TRANS| \
				      REC_MTHD_ROTATE))
#define REC_DEF_WINFN	(REC_WINFN_WELCH)
#define REC_DEF_WINSZ	(100)
#define REC_DEF_ITLIM	(10)
#define REC_DEF_XLIM	(30)
#define REC_DEF_YLIM	(30)
#define REC_DEF_RLIM	(60.0)
#define REC_DEF_DBG	REC_DBG_NONE

typedef struct
{
  double		weight;	          /* Weighting of the tie point pair */
  WlzDVertex2		first;		       /* Tie point in first section */
  WlzDVertex2		second;		      /* Tie point in second section */
} RecTiePointPair;

typedef struct
{
  int		linkcount;	      /* Number of times the section is used */
  int		index;		    /* Section index used for searching, etc */
  int		iterations;   /* Number of automatic registration iterations */
  double 	correl;				        /* Cross-correlation */
  char		*imageFile;		  /* Path for the section image file */
  WlzAffineTransform *transform;      /* Relative transform between sections */
  WlzObject	*obj;     /* Section image, maybe NULL if not read from file */
  WlzObject	*transObj;          /* Transformed section image, maybe NULL */
  WlzAffineTransform *cumTransform;      /* Cumulative transform, maybe NULL */
  WlzObject	*cumTransObj;  /* Transformed (cumulative) image, maybe NULL */
} RecSection;

typedef struct
{
  WlzWindowFnType function;
  WlzIVertex2	size;
  WlzIVertex2	offset;
} RecPPWindow;

typedef struct
{
  WlzSampleFn	function;
  int		factor;
} RecPPSample;

typedef struct
{
  int		approach;
  int		iteration;
  RecMethod	lastMethod;
  WlzAffineTransform *transform;
  double	correl;
  RecError	errFlag;
} RecState;

typedef struct
{
  RecTransformMode trMode;
  HGUDlpListItem *currentItem;
} RecSectionListAtrb;

typedef struct
{
  WlzObject	*obj;
  char 		*fileName;
  WlzEffFormat	fileFormat;
  int		gaussFlt;
  int		fastSam;
  int		intScale;
  int		greedy;
  WlzDVertex3	scale;
  int		matchHst;
  int		matchSec;
  int		clipSrc;
  WlzDBox3	srcBox;
  int		clipDst;
  WlzDBox3	dstBox;
} RecReconstruction;

typedef struct
{
  HGUDlpList	*list;
  RecSectionListAtrb attributes;
  RecReconstruction reconstruction;
} RecSectionList;

typedef struct
{
  RecMethod	method;
  double	xLim;
  double	yLim;
  double	rLim;
  int		itLim;
  int		firstIdx;
  int		lastIdx;
} RecControl;

typedef struct
{
  RecPPMethod	method;
  RecPPWindow	window;
  RecPPSample	sample;
  int		erode;
} RecPPControl;

typedef void    	(*RecSecUpdateFunction)(RecSection *, void *);
typedef void    	(*RecWorkFunction)(RecState *, void *);

typedef RecError	(*RecDbgFn)(char *, ...);
typedef RecError	(*RecDbgWlzFn)(WlzObject *, int);

extern RecDbgFn	recDbgOutFn;
#define REC_DBG_FN	(*recDbgOutFn)
#define REC_DBG(F,M)	((((F)&(recDbgMask))==(F))?REC_DBG_FN M:REC_ERR_NONE)

extern RecDbgWlzFn recDbgOutWlzFn;
#define REC_DBGW_FN	(*recDbgOutWlzFn)
#define REC_DBGW(F,O,X) \
		((((F)&(recDbgWlzMask))==(F))?REC_DBGW_FN((O),(X)):REC_ERR_NONE)


#ifdef  __cplusplus 
}
#endif /* __cplusplus */

#endif /* RECONSTRUCTTYPE_H */
