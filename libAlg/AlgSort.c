#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgSort_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         AlgSort.c
* \author       Bill Hill
* \date         June 2018
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2018],
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
* \brief	Basic sorting functions.
* \ingroup	AlgSort
*/

#include <stddef.h>
#include <Alg.h>

/*!
* \ingroup      AlgSort
* \brief	A wrapper for qsort()/mergesort() which aims to be consistent
*		across platforms. The default sort is qsort but on some
*		platforms, eg Apple OSX this is unreliable and mergesort
*		is used. The parameters are for qsort/mergesort.
* \param        base			Start of array to be sorted.
* \param        nElm			Number of elements in the array.
* \param        elmSz			Size of an array element.
* \param        compar			Comparison  function.
*/
void            AlgSort(void *base, size_t nElm, size_t elmSz,
                        int (*compar)(const void *, const void *))
{
#if defined(__APPLE__) && defined(__MACH__)
  mergesort(base, nElm, elmSz, compar);
#else
  qsort(base, nElm, elmSz, compar);
#endif
}

