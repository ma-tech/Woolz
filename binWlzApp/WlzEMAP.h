#ifndef	WLZ_EMAP_H
#define WLZ_EMAP_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzEMAP_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzEMAP.h
* \author       Richard Baldock
* \date         December 2006
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
* \ingroup      BinWlzApp
* \brief        Special defines for the EMAP extension to the Woolz
		library.
*/

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#include <sys/types.h>

typedef struct _WLZ_EMAP_WarpTransformStruct{
  /* source and destination models */
  char	*srcModel;
  char	*dstModel;

  /* last modified time of current data */
  time_t	mtime;

  /* source domains and projection */
  WlzObject	*srcDoms[3];
  WlzThreeDViewStruct	*srcProj;

  /* destination domains and projection */
  WlzObject	*dstDoms[3];
  WlzThreeDViewStruct	*dstProj;

  /* mesh transforms */
  WlzObject	*meshObj[3];
} WLZ_EMAP_WarpTransformStruct;

/* globals */
extern char *EMAP_WarpTransformsDir; /* defined in WlzEMAPDomainTransform.c */

extern int WlzEMAPIsMapping(char	*srcModel,
			    char	*dstModel,
			    char	*transformDir,
			    WlzErrorNum	*dstErr);
extern WlzErrorNum WlzEMAPFreeMapping(
  WLZ_EMAP_WarpTransformStruct	*mapping);

extern WLZ_EMAP_WarpTransformStruct *WlzEMAPGetMapping(
  char		*srcModel,
  char		*dstModel,
  char		*transformDir,
  WlzErrorNum	*dstErr);

extern WlzObject *WlzEMAPDomainTransform(char		*srcModel,
					 char		*dstModel,
					 WlzObject	*obj,
					 WlzErrorNum	*dstErr);

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
}
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#endif	/* !WLZ_EMAP_H Don't put anything after this line */
