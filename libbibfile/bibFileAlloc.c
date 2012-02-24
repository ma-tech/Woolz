#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _bibFileAlloc_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libbibfile/bibFileAlloc.c
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
* \brief	Functions for allocation and freeing of the bibtex
*		based record and field data structures.
* \ingroup	bibfile
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <bibFile.h>

BibFileRecord			*BibFileRecordMake(
				  char *name,
				  char *id,
				  BibFileField *field);
BibFileField			*BibFileFieldMake(
				  char *name,
				  char *value,
				  BibFileField *next);
BibFileField			*BibFileFieldMakeVa(
				  char *name,
				  char *value,
				  ...);
BibFileField			*BibFileFieldJoin(
				  BibFileField *field0,
				  BibFileField *field1,
				  ...);
void				BibFileRecordFree(
				  BibFileRecord **record);
void				BibFileFieldFree(
				  BibFileField **field);

/*!
* \return	void
* \ingroup	bibfile
* \brief	Free's the given record and sets it to NULL.
* \param	record			Record ptr to free.
*/
void		BibFileRecordFree(BibFileRecord **record)
{
  BibFileRecord	*tmpPtr;

  if(record && ((tmpPtr = *record) != NULL))
  {
    if(tmpPtr->name)
      (void )AlcFree(tmpPtr->name);
    if(tmpPtr->id)
      (void )AlcFree(tmpPtr->id);
    BibFileFieldFree(&(tmpPtr->field));
    (void )AlcFree(tmpPtr);
    *record = NULL;
  }
}

/*!
* \return	void
* \ingroup	bibfile
* \brief	Recursively free's the given field and sets it to NULL.	
* \param	field			Field ptr to free.
*/
void		BibFileFieldFree(BibFileField **field)
{
  BibFileField	*tmpPtr;

  if(field && ((tmpPtr = *field) != NULL))
  {
    if(tmpPtr->name)
      (void )AlcFree(tmpPtr->name);
    if(tmpPtr->value)
      (void )AlcFree(tmpPtr->value);
    tmpPtr = tmpPtr->next;
    (void )AlcFree(*field);
    *field = NULL;
    BibFileFieldFree(&tmpPtr);
  }
}

/*!
* \return	New record, NULL on error.
* \ingroup	bibfile
* \brief	Given record name and id strings and a field pointer
*		a new record is created.				
* \param	name			Record name string.
* \param	id			Record id string.
* \param	field			Top field of record.
*/
BibFileRecord	*BibFileRecordMake(char *name, char *id, BibFileField *field)
{
  int		ok = 1;
  char		*newName = NULL,
		*newId = NULL;
  BibFileRecord	*newRecord = NULL;

  if((newRecord = (BibFileRecord *)AlcCalloc(1,
  					     sizeof(BibFileRecord))) == NULL)
    ok = 0;
  if(ok && name && ((newName = AlcStrDup(name)) == NULL))
    ok = 0;
  if(ok && id && ((newId = AlcStrDup(id)) == NULL))
    ok = 0;
  if(ok)
  {
    newRecord->name = newName;
    newRecord->id = newId;
    newRecord->field = field;
  }
  else
  {
    if(newName)
      (void )AlcFree(newName);
    if(newId)
      (void )AlcFree(newId);
    if(newRecord)
    {
      (void )AlcFree(newRecord);
      newRecord = NULL;
    }
  }
  return(newRecord);
}

/*!
* \return	New field, NULL on error.
* \ingroup	bibfile
* \brief	Given field name and value strings and a field pointer
*		for the next field, a new field is created.		
* \param	name			Field name string.
* \param	value			Field value string.
* \param	next			Next field below this one.
*/
BibFileField	*BibFileFieldMake(char *name, char *value, BibFileField *next)
{
  int		ok = 1;
  char		*newName = NULL,
		*newValue = NULL;
  BibFileField	*newField = NULL;

  if((newField = (BibFileField *)AlcCalloc(1,
  					   sizeof(BibFileField))) == NULL) 
    ok = 0;
  if(ok && name && ((newName = AlcStrDup(name)) == NULL))
    ok = 0;
  if(ok && value && ((newValue = AlcStrDup(value)) == NULL))
    ok = 0;
  if(ok)
  {
    newField->name = newName;
    newField->value = newValue;
    newField->next = next;
  }
  else
  {
    if(newName)
      (void )AlcFree(newName);
    if(newValue)
      (void )AlcFree(newValue);
    if(newField)
    {
      (void )AlcFree(newField);
      newField = NULL;
    }
  }
  return(newField);
}

/*!
* \return	Top of new fields, maybe NULL on error.
* \ingroup	bibfile
* \brief	Given a list of field name and value string pairs which
* 		is terminated by a NULL. A hierarchy of new fields is
* 		created.
* \param	name			First field name string.
* \param	value			First field value string.
* \param	 ...			NULL terminated varargs list of
* 					name value pairs.
*/
/*VARARGS*/
BibFileField	*BibFileFieldMakeVa(char *name, char *value, ...)
{
  va_list	ap;
  int		ok = 0;
  BibFileField	*fieldTop = NULL,
		*field1,
		*field2;

  va_start(ap, value);
  if((fieldTop = BibFileFieldMake(name, value, NULL)) != NULL)
  {
    ok = 1;
    field1 = fieldTop;
  }
  while(ok && field1 && ((name = va_arg(ap, char *)) != NULL))
  {
    if(((value = va_arg(ap, char *)) == NULL) ||
       ((field2 = BibFileFieldMake(name, value, NULL)) == NULL))
      ok = 0;
    else
    {
      field1->next = field2;
      field1 = field2;
    }
  }
  va_end(ap);
  if(ok == 0)
    BibFileFieldFree(&fieldTop);
  return(fieldTop);
}

/*!
* \return	Top field of the join, NULL on error.
* \ingroup	bibfile
* \brief	Given a NULL terminated varargs list of fields, these
* 		are joined by filling in the appropriate next fields.
* \param	field0			Top most field of join.
* \param	field1			Field to be added to top field.
* \param	 ...			NULL terminated varargs list of
* 					fields to be joined.
*/
/*VARARGS*/
BibFileField	*BibFileFieldJoin(BibFileField *field0,
				  BibFileField *field1, ...)
{
  va_list	ap;
  BibFileField	*topField = NULL,
		*tmpField;

  if(field1 && field0)
  {
    topField = field0;
    va_start(ap, field1);
    while(field1 && field0)
    {
      while((tmpField = field0->next) != NULL)
	field0 = tmpField;
      field0->next = field1;
      field0 = field1;
      field1 = va_arg(ap, BibFileField *);
    }
    va_end(ap);
  }
  return(topField);
}
