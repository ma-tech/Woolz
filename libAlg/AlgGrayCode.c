#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgGrayCode_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgGrayCode.c
* \author       Bill Hill
* \date         August 2011
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
* \brief	Provides functions for computing Gray codes and
* 		their inverse. See Christopher H. Hamilton.
* 		"Range Searching Data Structures with Cache Locality"
* 		PhD Thesis, Dalhousie University, March 20011.
* 		The code within this file is derived from this thesis
* 		and the software it refers to.
* \ingroup	AlgBits
*/
#include <Alg.h>

/*!
* \return	Gray code of the given integer.
* \ingroup	AlgBits
* \brief	Computes the Gray code of the given value.
* \param	g			Given integer.
*/
unsigned int	AlgGrayCode(unsigned int g)
{

  g ^= g >> 1;
  return(g);
}

/*!
* \return	Gray code of the given integer.
* \ingroup	AlgBits
* \brief	Computes the Gray code of the given value.
* \param	g			Given long integer.
*/
unsigned long long AlgGrayCodeLL(unsigned long long g)
{

  g ^= g >> 1;
  return(g);
}

/*!
* \return	Integer corresponding to the given Gray code.
* \ingroup	AlgBits
* \brief	Computes the inverse Gray code, ie the value
* 		corresponding to the given Gray code.
* \param	g			Given Gray code.
*/
unsigned int	AlgGrayCodeInv(unsigned int g)
{
  g ^= (g >>  1);
  g ^= (g >>  2);
  g ^= (g >>  4);
  g ^= (g >>  8);
  g ^= (g >> 16);
  return(g);
}

/*!
* \return	Integer corresponding to the given Gray code.
* \ingroup	AlgBits
* \brief	Computes the inverse Gray code, ie the value
* 		corresponding to the given Gray code.
* \param	g			Given Gray code.
*/
unsigned long long AlgGrayCodeInvLL(unsigned long long g)
{
  g ^= (g >>  1);
  g ^= (g >>  2);
  g ^= (g >>  4);
  g ^= (g >>  8);
  g ^= (g >> 16);
  g ^= (g >> 32);
  return(g);
}
