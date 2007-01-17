#ifndef WLZ_EMAP_H
#define WLZ_EMAP_H
#if defined(__GNUC__)
#ident "MRC HGU $Id:"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id:"
#else static char _WlzEMAP.h[] = "MRC HGU $Id:";
#endif
#endif
/*!
* \file         WlzEMAP.h
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Mon Dec 18 13:27:27 2006
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par Copyright:
* Copyright (C) 2005 Medical research Council, UK.
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
* \ingroup      EMAP
* \brief        Special defines for the EMAP extension to the woolz library.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

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

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif	/* !WLZ_EMAP_H Don't put anything after this line */
