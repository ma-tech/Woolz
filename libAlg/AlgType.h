#ifndef ALGTYPE_H
#define ALGTYPE_H
#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas
* Title:        AlgType.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Header file with type definitions for the MRC Human
*		Genetics Unit numerical algorithm library.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif

/* Standard min, max, absolute value and nearest integer macros */
#define	ALG_MAX(X,Y)	(((X)>=(Y))?(X):(Y))
#define	ALG_MIN(X,Y)	(((X)<=(Y))?(X):(Y))
#define	ALG_ABS(X)	(((X)>0)?(X):(-(X)))
#define	ALG_NINT(X)	((int)(((X)<0)?((X)-(0.5)):((X)+(0.5))))


/* Standard math constants */
#define	ALG_M_E		(2.7182818284590452354)
#define	ALG_M_LOG2E	(1.4426950408889634074)
#define	ALG_M_LOG10E	(0.43429448190325182765)
#define	ALG_M_LN2	(0.69314718055994530942)
#define	ALG_M_LN10	(2.30258509299404568402)
#define ALG_M_PI	(3.1415926535897932384626433832795028841972)
#define ALG_M_PI_2	(1.5707963267948966192313216916397514420986)
#define ALG_M_PI_4	(0.7853981633974483096156608458198757210493)
#define	ALG_M_1_PI	(0.31830988618379067154)
#define	ALG_M_2_PI	(0.63661977236758134308)
#define	ALG_M_2_SQRTPI	(1.12837916709551257390)
#define ALG_M_SQRT2	(1.4142135623730950488016887242096980785696)
#define	ALG_M_SQRT1_2	(0.70710678118654752440)

/* Error codes */
typedef enum
{
  ALG_ERR_NONE		= (0),
  ALG_ERR_FUNC,				      /* Function parameters invalid */
  ALG_ERR_MALLOC,				/* Memory allocation failure */
  ALG_ERR_SINGULAR,					  /* Singular matrix */
  ALG_ERR_HOMOGENEOUS,				       /* Homogeneous matrix */
  ALG_ERR_MAX
} AlgError;

/* Debug mask values */
typedef enum
{
  ALG_DBG_NONE          = (0),
  ALG_DBG_LVL_1         = (1),
  ALG_DBG_LVL_2         = (1<<1),
  ALG_DBG_LVL_3         = (1<<2),
  ALG_DBG_LVL_FN        = (1<<3)
} AlgDbgMask;

typedef AlgError        (*AlgDbgFn)(char *, ...);
 
extern AlgDbgFn		algDbgOutFn;
#define ALG_DBG_FN      (*algDbgOutFn)
#define ALG_DBG(F,M)    ((((F)&(algDbgMask))==(F))?ALG_DBG_FN M:ALG_ERR_NONE)
 


#ifdef  __cplusplus 
}
#endif /* __cplusplus */

#endif /* ! ALGTYPE_H */
