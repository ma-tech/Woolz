#ifndef BIBFILE_H
#define BIBFILE_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _bibFile_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libbibfile/bibFile.h
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
* \brief	Types and constants for the bibtex based file syntax
*		used for serial section data, ....
* \ingroup	bibfile
*/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <Alc.h>

typedef enum
{
  BIBFILE_ER_NONE		= 0,
  BIBFILE_ER_MALLOC,
  BIBFILE_ER_SYNTAX,
  BIBFILE_ER_WRITE,
  BIBFILE_ER_READ,
  BIBFILE_ER_EOF
} BibFileError;


typedef struct _BibFileField
{
  char		*name;
  char		*value;
  struct _BibFileField *next;
} BibFileField;

typedef struct
{
  char		*name;
  char		*id;
  BibFileField	*field;
} BibFileRecord;

/* From	bibFileAlloc.c */
extern BibFileRecord 		*BibFileRecordMake(
				  char *name,
				  char *id,
				  BibFileField *field);
extern BibFileField 		*BibFileFieldMake(
				  char *name,
				  char *value,
				  BibFileField *next);
extern BibFileField		*BibFileFieldMakeVa(
				  char *name,
				  char *value,
				  ...);
extern BibFileField         	*BibFileFieldJoin(
				  BibFileField *field0,
                                  BibFileField *field1,
				  ...);
extern void     		BibFileRecordFree(
				  BibFileRecord **record);
extern void               	BibFileFieldFree(
				  BibFileField **field);
extern char			*BibFileStrDup(
				  const char *s1);		

/* From bibFileIO.c */
extern BibFileError 		BibFileRecordRead(
				  BibFileRecord **record,
				  char **eMsg,
		                  FILE *fP);
extern BibFileError		BibFileRecordWrite(
				  FILE *fP,
				  char **eMsg,
		                  BibFileRecord *record);
extern BibFileError		BibFileFieldRead(
				  BibFileField **field,
				  char **eMsg,
		                  int *endFlag,
				  FILE *fP);
extern BibFileError		BibFileFieldWrite(
				  FILE *fP,
				  char **eMsg,
		                  BibFileField *field);
extern BibFileError      	BibFileEscapeRestrictedChar(
				  char *pString,
				  char **outString);
extern BibFileError  		BibFileUnEscapeRestrictedChar(
				  char *pString,
				  char **outString);
/* From bibFileParse.c */
int				BibFileFieldParseFmt(
				  BibFileField *topField,
                                  void *value,
				  char *fmt,
				  char *name,
				  ...);


#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* BIBFILE_H */
