#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        bibFileAlloc.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for allocation and freeing of the bibtex
*		based record and field data structures.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <bibFile.h>

BibFileRecord	*BibFileRecordMake(char *name, char *id, BibFileField *field);
BibFileField	*BibFileFieldMake(char *name, char *value, BibFileField *next),
		*BibFileFieldMakeVa(char *name, char *value, ...),
		*BibFileFieldJoin(BibFileField *field0,
				  BibFileField *field1, ...);
void		BibFileRecordFree(BibFileRecord **record),
		BibFileFieldFree(BibFileField **field);

/************************************************************************
* Function:	BibFileRecordFree					*
* Returns:	void							*
* Purpose:	Free's the given record and sets it to NULL.		*
* Global refs:	-							*
* Parameters:	BibFileRecord **record:	Record ptr to be free'd.	*
************************************************************************/
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

/************************************************************************
* Function:	BibFileFieldFree					*
* Returns:	void							*
* Purpose:	Recursively free's the given field and sets it to NULL.	*
* Global refs:	-							*
* Parameters:	BibFileField **field:	Field ptr to be free'd.		*
************************************************************************/
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

/************************************************************************
* Function:	BibFileRecordMake					*
* Returns:	BibFileRecord *:	New record, NULL on error.	*
* Purpose:	Given record name and id strings and a field pointer	*
*		a new record is created.				*
* Global refs:	-							*
* Parameters:	char *name:		Record name string.		*
*		char *id:		Record id string.		*
*		BibFileField *field:	Top field of record.		*
************************************************************************/
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

/************************************************************************
* Function:	BibFileFieldMake					*
* Returns:	BibFileField *:		New field, NULL on error.	*
* Purpose:	Given field name and value strings and a field pointer	*
*		for the next field, a new field is created.		*
* Global refs:	-							*
* Parameters:	char *name:		Field name string.		*
*		char *value:		Field value string.		*
*		BibFileField *next:	Next field below this one.	*
************************************************************************/
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

/************************************************************************
* Function:	BibFileFieldMakeVa					*
* Returns:	BibFileField *:		Top of new fields, maybe NULL	*
*					on error.			*
* Purpose:	Given a list of field name and value string pairs which	*
*		is terminated by a NULL. A hierarchy of new fields is 	*
*		created.						*
* Global refs:	-							*
* Parameters:	char *name:		First field name string.	*
*		char *value:		First field value string.	*
*		...			NULL terminated varargs list of	*
*					name value pairs.		*
************************************************************************/
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

/************************************************************************
* Function:	BibFileFieldJoin					*
* Returns:	BibFileField *:		Top field of the join, NULL	*
*					on error.			*
* Purpose:	Given a NULL terminated varargs list of fields, these	*
*		are joined by filling in the appropriate next fields.	*
* Global refs:	-							*
* Parameters:	BibFileField *field0:	Top most field of join.		*
*		BibFileField *field1:	Field to be added to top field.	*
*		...			NULL terminated varargs list of	*
*					fields to be joined.		*
************************************************************************/
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
