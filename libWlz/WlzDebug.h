#ifndef WLZ_DEBUG_H
#define WLZ_DEBUG_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzDebug.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Defines the Woolz debug masks and function prototypes.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************************************************************
* Woolz debugging type definitions.					*
************************************************************************/

typedef enum
{
  WLZ_DBG_NONE		= (0),				  /* No debug output */
   WLZ_DBG_LVL_1	= (1),			       /* Least debug output */
   WLZ_DBG_LVL_2	= (1<<1),		/* Intermediate debug output */
   WLZ_DBG_LVL_3	= (1<<2),			/* Most debug output */
   WLZ_DBG_LVL_FN	= (1<<3),		/* Function entry and return */
   WLZ_DBG_ALLOC	= (1<<4) 		   /* Allocation and freeing */
} WlzDbgMask;

typedef WlzErrorNum	(*WlzDbgFn)(char *, ...);
typedef WlzErrorNum	(*WlzDbgObjFn)(WlzObject *, int);

/************************************************************************
* Woolz debugging prototypes.						*
************************************************************************/
extern WlzDbgMask	wlzDbgMask;
extern WlzDbgMask	wlzDbgObjMask;
extern void		*wlzDbgData;
extern void		*wlzDbgObjData;
extern WlzDbgFn	wlzDbgOutFn;
extern WlzDbgObjFn	wlzDbgOutObjFn;

extern WlzErrorNum	WlzDbgWrite(char *, ...);
extern WlzErrorNum	WlzDbgObjWrite(WlzObject *, int);

/************************************************************************
* Woolz debugging macros.						*
************************************************************************/
#define WLZ_DBG(F,M) \
		      ((((F)&(wlzDbgMask))==(F))?(*wlzDbgOutFn) M:WLZ_ERR_NONE)
#define WLZ_DBGOBJ(F,O,X) \
	 ((((F)&(wlzDbgObjMask))==(F))?(*wlzDbgOutObjFn)((O),(X)):WLZ_ERR_NONE)


#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif	/* !WLZ_DEBUG_H Don't put anything after this line */
