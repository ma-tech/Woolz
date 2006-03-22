#ifndef ALG_H
#define ALG_H
#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _Alg_h[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/Alg.h
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
* \brief        Main (top-level) header file for the Woolz numerical
* 		algorithms library.
* \ingroup     	Alg
* \todo         -
* \bug          None known.
*/

#ifndef __EXTENSIONS__
#define __EXTENSIONS__
#endif

#ifdef __cplusplus
using namespace std;
#else
#include <stdlib.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <Alc.h>
#include <AlgType.h>
#include <AlgProto.h>

#ifdef __cplusplus
}
#endif

#endif /* ! ALG_H */
