#ifndef BIBFILE_H
#define BIBFILE_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        bibFile.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Types and constants for the bibtex based file syntax
*		used for serial section data, ....
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

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
extern BibFileRecord *BibFileRecordMake(char *name, char *id, BibFileField *field);
extern BibFileField *BibFileFieldMake(char *name, char *value, BibFileField *next),
		*BibFileFieldMakeVa(char *name, char *value, ...),
                *BibFileFieldJoin(BibFileField *field0,
                                  BibFileField *field1, ...);
extern void     BibFileRecordFree(BibFileRecord **record),
                BibFileFieldFree(BibFileField **field);
extern char	*BibFileStrDup(const char *s1);		

/* From bibFileIO.c */
extern BibFileError BibFileRecordRead(BibFileRecord **record, char **eMsg,
		                      FILE *fP),
		BibFileRecordWrite(FILE *fP, char **eMsg,
		                   BibFileRecord *record),
		BibFileFieldRead(BibFileField **field, char **eMsg,
		                 int *endFlag, FILE *fP),
		BibFileFieldWrite(FILE *fP, char **eMsg,
		                  BibFileField *field),
                BibFileEscapeRestrictedChar(char *pString,
					    char **outString),
                BibFileUnEscapeRestrictedChar(char *pString,
					      char **outString);
/* From bibFileParse.c */
int		BibFileFieldParseFmt(BibFileField *topField,
                                     void *value, char *fmt, char *name, ...);


#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* BIBFILE_H */
