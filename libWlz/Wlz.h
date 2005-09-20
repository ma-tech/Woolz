#ifndef WLZ_H
#define WLZ_H
#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/Wlz.h
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
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
* \brief	Main (top-level) Woolz header file which includes all
* 		other header files required by Woolz.
* \ingroup	Wlz
* \todo         -
* \bug          None known.
*/

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <math.h>

#include <Alc.h>

#include <Alg.h>

#include <WlzError.h>
#include <WlzType.h>
#include <WlzDebug.h>
#include <WlzProto.h>
#include <WlzMacro.h>

#ifdef __cplusplus
}
#endif

#endif	/* !WLZ_H Don't put anything after this line */
