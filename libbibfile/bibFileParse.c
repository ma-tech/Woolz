#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        bibFileParse.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions to parse the bibtex based file syntax.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <bibFile.h>

int		BibFileFieldParseFmt(BibFileField *topField,
				     void *value, char *fmt, char *name, ...);

/************************************************************************
* Function:	BibFileFieldParseFmt					*
* Returns:	int:			Count of fields parsed, zero on	*
*					error.				*
* Purpose:	Given the top most field from which to parse fields and	*
*		a NULL terminated varargs list of field value pointers,	*
*		format strings (for sscanf) and field name strings. All	*
*		fields are parsed for all name strings. Parsed values 	*
*		are allowed to overwrite those previously parsed.	*
*		Note that string values are parsed using a char ** and	*
*		are strdup()'d.						*
* Global refs:	-							*
* Parameters:	BibFileField *topField:	Top most field in hierarchy of	*
*					fields to be parsed.		*
*		void *value:		First value pointer for result,	*
*					passed directly to sscanf.	*
*		char *fmt:		First format string for sscanf.	*
*		char *name:		First field name to be matched.	*
*		...:			NULL terminated varargs list of	*
*					value, format and name triples.	*
************************************************************************/
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
