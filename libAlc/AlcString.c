#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcString_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlc/AlcString.c
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
* \brief        Provides functions for string duplication.
* \ingroup	AlcString
*/
#include <stdio.h>
#include <string.h>
#include <Alc.h>

/*!
* \return	Duplicated string or NULL on error.
* \ingroup	AlcString
* \brief	Allocates space for and duplicates the given NULL
*		terminated character string.
* \param	srcStr			 Given string.
*/
char		*AlcStrDup(const char *srcStr)
{
  char		*dstStr = NULL;

  if(srcStr &&
     ((dstStr = AlcMalloc((strlen(srcStr) + 1) * sizeof(char))) != NULL))
  {
    (void )strcpy(dstStr, srcStr);
  }
  return(dstStr);
}

/*!
* \return	New concatonated string that should be freed using AlcFree()
* 		or NULL on error.
* \ingroup	AlcString
* \brief	Concatonates the three given strings in order into a new
* 		allocated string buffer. If all strings are null or zero
* 		length then NULL is returned, otherwise if a string is
* 		null it is omitted.
* \param	s0			First string.
* \param	s1			Second string.
* \param	s2			Third string.
*/
char 		*AlcStrCat3(const char *s0, const char *s1, const char *s2)
{
  int		totLen;
  int		len[3];
  char 		*dstStr = NULL;
  
   len[0] = ((s0)? strlen(s0): 0);
   len[1] = ((s1)? strlen(s1): 0);
   len[2] = ((s2)? strlen(s2): 0);
   if((totLen = len[0] + len[1] + len[2]) > 0)
   {
     dstStr = AlcMalloc(sizeof(char) * (totLen + 1));
   }
   if(dstStr)
   {
     if(len[0])
     {
       strcpy(dstStr, s0);
     }
     if(len[1])
     {
       strcpy(dstStr + len[0], s1);
     }
     if(len[2])
     {
       strcpy(dstStr + len[0] + len[1], s2);
     }
     dstStr[totLen] = '\0';
   }
   return(dstStr);
}

/*!
* \return	Numeric hash value for the given sting.
* \ingroup	AlcString
* \brief	A hash function for strings based on Paul Hsieh's
* 		SuperFastHash function which is covered by the LGPL 2.1
* 		license amongst others. The returned value will always
* 		be zero if the given string is null of has zero length.
* \param	str			Given string.
*/
unsigned int	AlcStrSFHash(const char *str)
{
  unsigned int	l,
		r,
		t,
  		h = 0;
  const char	*s;

  if(str && *str)
  {
    l = strlen(str);
    r = l & 3;
    h = l;
    l >>= 2;
    s = str;
    while(l > 0)
    {
      h += (s[1] << 8) | s[0];
      t  = (((s[3] << 8) | s[2]) << 11) ^ h;
      h  = (h << 16) ^ t;
      h += h >> 11;
      s += 4;
      --l;
    }
    switch(r)
    {
      case 3:
	h += (s[1] << 8) | s[0];
	h ^= h << 16;
	h ^= s[2] << 18;
	h += h >> 11;
	break;
      case 2:
	h += (s[1] << 8) | s[0];
	h ^= h << 11;
	h += h >> 17;
	break;
      case 1:
	h += s[0];
	h ^= h << 10;
	h += h >> 1;
      default:
	break;
    }
    h ^= h << 3;
    h += h >> 5;
    h ^= h << 4;
    h += h >> 17;
    h ^= h << 25;
    h += h >> 6;
  }
  return(h);
}
