#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _bibFileParse_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libbibfile/bibFileParse.c
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
* \brief	Functions to parse the bibtex based file syntax.
* $Revision$
* Maintenance:	Log changes 
* \ingroup	bibfile
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <bibFile.h>

int		BibFileFieldParseFmt(BibFileField *topField,
				     void *value, char *fmt, char *name, ...);

/*!
* \return	Count of fields parsed, zero on error.
* \ingroup	bibfile
* \brief	Given the top most field from which to parse fields and
*		a NULL terminated varargs list of field value pointers,
*		format strings (for sscanf) and field name strings. All
*		fields are parsed for all name strings. Parsed values
*		are allowed to overwrite those previously parsed.
*		Note that string values are parsed using a char ** and
*		are AlcStrDup()'d.						
* \param	topField		Top most field in hierarchy of
* 					fields to be parsed.
* \param	value			First value pointer for result,
* 					passed directly to sscanf.
* \param	fmt			First format string for sscanf.
* \param	name			First field name to be matched.
* \param	...			NULL terminated varargs list of
* 					value, format and name triples.
*/
/*VARARGS*/
int		BibFileFieldParseFmt(BibFileField *topField,
				     void *value, char *fmt, char *name, ...)
{
  va_list	ap;
  BibFileField	*field;
  int		ok = 1,
		count = 0;
  const char	strFmt[] = "%s"; /* Strings are special case, see note above */

  if(topField && value && fmt && name)
  {
    va_start(ap, name);
    while(ok && value && fmt && name)
    {
      field = topField;
      while(field && ok)
      {
	if(field->name && field->value)
	{
	  if(strcmp(field->name, name) == 0)
	  {
	    if(strcmp(fmt, strFmt) == 0)
	      *((char **)value) = AlcStrDup(field->value);
	    else
	    {
	      if(sscanf(field->value, fmt, value) != 1)
		ok = 0;
	    }
	    if(ok)
	      ++count;
	  }
	  field = field->next;
	}
	else
	  ok = 0;
      }
      if(((value = va_arg(ap, void **)) != NULL) &&
	 ((fmt = va_arg(ap, char *)) != NULL))
	name = va_arg(ap, char *);
    }
    va_end(ap);
    if(ok == 0)
      count = 0;
  }
  return(count);
}
