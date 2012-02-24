#ifndef WLZ_DEBUG_H
#define WLZ_DEBUG_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDebug_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzDebug.h
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
* \brief	Defines the Woolz debug masks and function prototypes.
* \ingroup	WlzDebug
*/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
* \enum		_WlzDbgMask
* \ingroup      WlzDebug
* \brief	Woolz debug bit masks.
*		Typedef: ::WlzDbgMask
*/
typedef enum _WlzDbgMask
{
  WLZ_DBG_NONE		= (0),	  	/*!< No debug output */
   WLZ_DBG_LVL_1	= (1),		/*!< Least debug output */
   WLZ_DBG_LVL_2	= (1<<1),	/*!< Intermediate debug output */
   WLZ_DBG_LVL_3	= (1<<2),	/*!< Most debug output */
   WLZ_DBG_LVL_FN	= (1<<3),	/*!< Function entry and return */
   WLZ_DBG_ALLOC	= (1<<4) 	/*!< Allocation and freeing */
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
