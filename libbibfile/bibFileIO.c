#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _bibFileIO_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libbibfile/bibFileIO.c
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
* \brief	File I/O functions to parse the bibtex based file
*		syntax used for serial section data, ....
*		
*		The syntax is:
\verbatim
                  file :=       (comment | record)
                  comment :=    COMMENT text EOL
                  record :=     NEW name OPEN id field (SEP field)*
                                CLOSE
                  field :=      name EQUAL OPEN value CLOSE
\endverbatim
*               Example:
\verbatim
                  % this is a comment which is ignored
                  @Comment
                  { 0,
                    Text = {This is a comment too, but it is parsed}
                  }
                  @Fred{1, FFred={First bit},
                    SFred={Second line of text} % comment here too
                  }
                  % another comment
                  @Bert{2,FBert={bert's 1st value has a newline
                  and an escaped newline \
                  too},
                  SBert =
                  {bert's 2nd includes a \} character}}
\endverbatim
*		All the tokens (eg OPEN, CLOSE, ...) are parsed using
*		regular expressions and character lists held in tables.
*		The parsing could be changed/extended by modifying the
*		tables, eg changing the definition of OPEN/CLOSE to
*		allow the double quote character.
*		Using regular expressions has made coding easier and
*		the parsing more flexible at the cost of speed. On a
*		Sparc 10/40 with -O4 optimisation about 20% of the
*		execution time is spent matching regular expressions.
* \ingroup	bibfile
* \todo         -
* \bug          None known.
*/

#define	BIBFILE_REGEXPR_CACHE 			/* Cache regular expressions */
#undef	BIBFILE_REGEXPR_REGEX		       /* Use regex(3) and regcmp(3) */
#define	BIBFILE_REGEXPR_REGEXEC		    /* Use regexec(3) and regcomp(3) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef BIBFILE_REGEXPR_REGEX
#include <libgen.h>
#endif /* BIBFILE_REGEXPR_REGEX */
#ifdef BIBFILE_REGEXPR_REGEXEC
#include <sys/types.h>
#include <regex.h>
#endif /* BIBFILE_REGEXPR_REGEXEC */
#include <bibFile.h>

#define BIBFILE_EOL		'\n'
#define BIBFILE_RESTRICTED_NUM	4
#define BIBFILE_STRBUF_START	256
#define BIBFILE_STRBUF_INC	1024

typedef enum 	           /* Warning: Keep consistent with the tables below */
{
  BIBFILE_TOK_MIN = 0,
  BIBFILE_TOK_NEW = 0,
  BIBFILE_TOK_SEP,
  BIBFILE_TOK_OPEN,
  BIBFILE_TOK_CLOSE,
  BIBFILE_TOK_EQUAL,
  BIBFILE_TOK_ESC,
  BIBFILE_TOK_COMMENT,
  BIBFILE_TOK_RECNAME,
  BIBFILE_TOK_RECID,
  BIBFILE_TOK_FLDNAME,
  BIBFILE_TOK_FLDVAL,
  BIBFILE_TOK_MAX
} BibFileTok;

BibFileError	BibFileRecordRead(BibFileRecord **record, char **eMsg,
				  FILE *fP),
		BibFileRecordWrite(FILE *fP, char **eMsg,
				   BibFileRecord *record),
		BibFileFieldRead(BibFileField **field, char **eMsg,
				 int *endFlag, FILE *fP),
		BibFileFieldWrite(FILE *fP, char **eMsg,
				  BibFileField *field),
                BibFileEscapeRestrictedChar(char *inString,
					    char **outString),
                BibFileUnEscapeRestrictedChar(char *inString,
					      char **outString);

static BibFileError BibFileCharRead(char *newC, FILE *fP, int skipSp,
			            int escape),
		BibFileCharUnread(char newC, FILE *fP),
		BibFileCharRegEx(char given, BibFileTok tok),
		BibFileFieldError(char **eMsg, BibFileField *field),
		BibFileRecordError(char **eMsg, BibFileRecord *record),
    		BibFileSkipSpaces(FILE *fP),
		BibFileStrRead(char **str, BibFileTok tok, int skipSp,
			       FILE *fP),
		BibFileCharMatch(BibFileTok tok, int skipFlag, int skipSp,
				 FILE *fP);

static char	*bibFileTokIn[BIBFILE_TOK_MAX] = /* Warning: Keep consistent */
{
  "@",							  /* BIBFILE_TOK_NEW */
  ",",							  /* BIBFILE_TOK_SEP */
  "\\{",						 /* BIBFILE_TOK_OPEN */
  "\\}",						/* BIBFILE_TOK_CLOSE */
  "=",							/* BIBFILE_TOK_EQUAL */
  "\\",							  /* BIBFILE_TOK_ESC */
  "%",						      /* BIBFILE_TOK_COMMENT */
  "[A-Za-z0-9_]",				      /* BIBFILE_TOK_RECNAME */
  "[-+A-Za-z0-9_.]",			        /* BIBFILE_TOK_RECID */
  "[A-Za-z0-9_]",				      /* BIBFILE_TOK_FLDNAME */
  "^[^{}]$",					       /* BIBFILE_TOK_FLDVAL */
};
static char    *bibFileTokOut[BIBFILE_TOK_MAX] = /* Warning: Keep consistent */
{
  "@",							  /* BIBFILE_TOK_NEW */
  ",",							  /* BIBFILE_TOK_SEP */
  "{",						 	 /* BIBFILE_TOK_OPEN */
  "}",							/* BIBFILE_TOK_CLOSE */
  "=",							/* BIBFILE_TOK_EQUAL */
  "\\",							  /* BIBFILE_TOK_ESC */
  "%",						      /* BIBFILE_TOK_COMMENT */
  "?",				      		      /* BIBFILE_TOK_RECNAME */
  "?",						        /* BIBFILE_TOK_RECID */
  "?",						      /* BIBFILE_TOK_FLDNAME */
  "?",						       /* BIBFILE_TOK_FLDVAL */
};
static char	*bibFileRestrictedChars[BIBFILE_RESTRICTED_NUM] =
{
    "%",
    "}",
    "{",
    "\\"
};

/*!
* \return	Non zero if read fails.
* \ingroup	bibfile
* \brief	Read a bibtex style record from the given stream and
*		allocate storage as required.				
* \param	record			Destination ptr for the record.
* \param	eMsg			Ptr for error message strings.
* \param	fP			Input file stream.
*/
BibFileError	BibFileRecordRead(BibFileRecord **record, char **eMsg,
				  FILE *fP)
{
  int		endFlag = 0;
  BibFileError	errFlag = BIBFILE_ER_NONE;
  BibFileField	**field;

  errFlag = BibFileCharMatch(BIBFILE_TOK_NEW, 1, 1, fP);
  if(errFlag != BIBFILE_ER_EOF)
  {
    if((errFlag == BIBFILE_ER_NONE) &&
       ((*record = (BibFileRecord *)AlcCalloc(1,
       					      sizeof(BibFileRecord))) == NULL))
    {
      errFlag = BIBFILE_ER_MALLOC;
    }
    if(errFlag == BIBFILE_ER_NONE)
    {
      errFlag = BibFileStrRead(&((*record)->name), BIBFILE_TOK_RECNAME,
			       0, fP);
    }
    if(errFlag == BIBFILE_ER_NONE)
    {
      errFlag = BibFileCharMatch(BIBFILE_TOK_OPEN, 1, 1, fP);
    }
    if(errFlag == BIBFILE_ER_NONE)
    {
      errFlag = BibFileStrRead(&((*record)->id), BIBFILE_TOK_RECID, 1, fP);
    }
    if(errFlag == BIBFILE_ER_NONE)
    {
      errFlag = BibFileCharMatch(BIBFILE_TOK_SEP, 1, 1, fP);
    }
    if(errFlag == BIBFILE_ER_NONE)
    {
      field = &((*record)->field);
      errFlag = BibFileFieldRead(field, eMsg, &endFlag, fP); /* read field 1 */
    }
    while((errFlag == BIBFILE_ER_NONE) && (endFlag == 0))
    {
      field = &((*field)->next);
      errFlag = BibFileFieldRead(field, eMsg, &endFlag,
				 fP);		/* read any remaining fields */
    }
    if(errFlag == BIBFILE_ER_NONE)
    {
      errFlag = BibFileCharMatch(BIBFILE_TOK_CLOSE, 1, 1, fP);
    }
    if(errFlag != BIBFILE_ER_NONE)
    {
      (void )BibFileRecordError(eMsg, *record);
    }
  }
  if((errFlag != BIBFILE_ER_NONE) && *record)
  {
    BibFileRecordFree(record);
  }
  return(errFlag);
}

/*!
* \return	Non zero if write fails.
* \ingroup	bibfile
* \brief	Write a bibtex style record to the given stream.
* \param	fP			Output file stream.
* \param	eMsg			Ptr for error message strings.
* \param	record			Record ptr.
*/
BibFileError	BibFileRecordWrite(FILE *fP, char **eMsg,
				   BibFileRecord *record)
{
  BibFileError	errFlag = BIBFILE_ER_NONE;
  BibFileField	*field;

  if(!(record && record->name && record->id && record->field))
  {
    errFlag = BIBFILE_ER_SYNTAX;
  }
  if((errFlag == BIBFILE_ER_NONE) &&
     (fprintf(fP, "%s%s%c%s%c  %s%s%c", bibFileTokOut[BIBFILE_TOK_NEW],
	      record->name, BIBFILE_EOL,
	      bibFileTokOut[BIBFILE_TOK_OPEN], 
	      BIBFILE_EOL, record->id,
	      bibFileTokOut[BIBFILE_TOK_SEP],
	      BIBFILE_EOL) == EOF))
  {
      errFlag = BIBFILE_ER_WRITE;
  }
  if((errFlag == BIBFILE_ER_NONE) && ((field = record->field) == NULL))
  {
    errFlag = BIBFILE_ER_SYNTAX;
  }
  while((errFlag == BIBFILE_ER_NONE) && field)
  {
    (void )fprintf(fP, "  ");
    errFlag = BibFileFieldWrite(fP, eMsg, field);
    field = field->next;
    if(field)
    {
      (void )fprintf(fP, "%s", bibFileTokOut[BIBFILE_TOK_SEP]);
    }
    (void )fprintf(fP, "%c", BIBFILE_EOL);
  }
  if((errFlag == BIBFILE_ER_NONE) &&
     (fprintf(fP, "%s%c",
	      bibFileTokOut[BIBFILE_TOK_CLOSE], BIBFILE_EOL) == EOF))
  {
    errFlag = BIBFILE_ER_WRITE;
  }
  if(errFlag != BIBFILE_ER_NONE)
  {
    (void )BibFileRecordError(eMsg, record);
  }
  return(errFlag);
}

/*!
* \return	Non zero if fails to build an error string!
* \ingroup	bibfile
* \brief	Build an error string showing with some record info.
* \param	eMsg			Ptr for error message strings.
* \param	record			Record ptr.
*/
static BibFileError BibFileRecordError(char **eMsg, BibFileRecord *record)
{
  char		*tmpPtr;
  BibFileError	errFlag = BIBFILE_ER_NONE;
  char		errBuf[256];

  if(record && eMsg)
  {
    tmpPtr = *eMsg;
    (void )sprintf(errBuf, "%s%srecord %s%s%s%s%s ...",
		   tmpPtr?tmpPtr:"",
		   tmpPtr?", ":"",
		   bibFileTokOut[BIBFILE_TOK_NEW],
		   ((record)->name)?((record)->name):"?",
		   bibFileTokOut[BIBFILE_TOK_OPEN],
		   ((record)->id)?((record)->id):"?",
		   bibFileTokOut[BIBFILE_TOK_SEP]);
    if(tmpPtr)
    {
      (void )AlcFree(tmpPtr);
    }
    if((*eMsg = AlcStrDup(errBuf))== NULL)
    {
      errFlag = BIBFILE_ER_MALLOC;
    }
  }
  return(errFlag);
}

/*!
* \return	Non zero if read fails.
* \ingroup	bibfile
* \brief	Read a bibtex style field from the given stream and
*		allocate storage as required.				
* \param	field			Destination ptr for the field.
* \param	eMsg			Ptr for error message strings.
* \param	endFlag			Set to non zero if last field.
* \param	fP			Input file stream.
*/
BibFileError	BibFileFieldRead(BibFileField **field, char **eMsg, 
				 int *endFlag, FILE *fP)
{
  BibFileError	errFlag = BIBFILE_ER_NONE;
  BibFileField	*newField;

  if(endFlag)
  {
    *endFlag = 0;
  }
  if((newField = (BibFileField *)AlcCalloc(1,
  					   sizeof(BibFileField))) == NULL)
  {
    errFlag = BIBFILE_ER_MALLOC;
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    errFlag = BibFileStrRead(&(newField->name), BIBFILE_TOK_FLDNAME, 1, fP);
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    errFlag = BibFileCharMatch(BIBFILE_TOK_EQUAL, 1, 1, fP);
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    errFlag = BibFileCharMatch(BIBFILE_TOK_OPEN, 1, 1, fP);
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    errFlag = BibFileStrRead(&(newField->value), BIBFILE_TOK_FLDVAL,
			     0, fP);
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    errFlag = BibFileCharMatch(BIBFILE_TOK_CLOSE, 1, 1, fP);
    if(errFlag == BIBFILE_ER_NONE)
    {
      errFlag = BibFileCharMatch(BIBFILE_TOK_SEP, 0, 1, fP);
    }
    if(errFlag == BIBFILE_ER_NONE)
    {
      (void )BibFileCharMatch(BIBFILE_TOK_SEP, 1, 1, fP);
    }
    else
    {
      errFlag = BibFileCharMatch(BIBFILE_TOK_CLOSE, 0, 1, fP);
      if((errFlag == BIBFILE_ER_NONE) && endFlag)
      {
        *endFlag = 1;
      }
    }
  }
  if(errFlag == BIBFILE_ER_EOF)
  {
    errFlag = BIBFILE_ER_SYNTAX;
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    *field = newField;
  }
  else
  {
    (void )BibFileFieldError(eMsg, newField);
    *field = NULL;
    if(newField)
    {
      (void )AlcFree(newField);
    }
  }
  return(errFlag);
}

/*!
* \return	Non zero if write fails.
* \ingroup	bibfile
* \brief	Write a bibtex style field to the given stream.
* \param	fP			Output file stream.
* \param	eMsg			Ptr for error message strings.
* \param	field			Field ptr.
*/
BibFileError	BibFileFieldWrite(FILE *fP, char **eMsg, BibFileField *field)
{
  char		*tmpStr;
  BibFileError	errFlag = BIBFILE_ER_NONE;

  if(!(field && field->name && field->value))
  {
    errFlag = BIBFILE_ER_SYNTAX;
  }
  if((errFlag == BIBFILE_ER_NONE) &&
     (fprintf(fP, "%s %s %s",
	      field->name,
	      bibFileTokOut[BIBFILE_TOK_EQUAL],
	      bibFileTokOut[BIBFILE_TOK_OPEN]) == EOF))
  {
    errFlag = BIBFILE_ER_WRITE;
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    tmpStr = field->value;
    while(*tmpStr)
    {
      if((*tmpStr == *bibFileTokOut[BIBFILE_TOK_OPEN]) ||
	 (*tmpStr == *bibFileTokOut[BIBFILE_TOK_CLOSE]) ||
	 (*tmpStr == *bibFileTokOut[BIBFILE_TOK_COMMENT]))
      {
	(void )putc((int )*bibFileTokOut[BIBFILE_TOK_ESC], fP);
      }
      (void )putc((int )*tmpStr, fP);
      ++tmpStr;
    }
  }
  if((errFlag == BIBFILE_ER_NONE) &&
     (fprintf(fP, "%s", bibFileTokOut[BIBFILE_TOK_CLOSE]) == EOF))
  {
    errFlag = BIBFILE_ER_WRITE;
  }
  if(errFlag != BIBFILE_ER_NONE)
  {
    (void )BibFileFieldError(eMsg, field);
  }
  return(errFlag);
}

/*!
* \return	Non zero if fails to build an error string!
* \ingroup	bibfile
* \brief	Build an error string showing with some field info.
* \param	eMsg			Ptr for error message strings.
* \param	field			Field ptr.
*/
static BibFileError BibFileFieldError(char **eMsg, BibFileField *field)
{
  char		*tmpPtr;
  BibFileError	errFlag = BIBFILE_ER_NONE;
  char		errBuf[256];

  if(field && eMsg)
  {
    tmpPtr = *eMsg;
    (void )sprintf(errBuf, "%s%sfield %s %s %s...",
		   tmpPtr?tmpPtr:"",
		   tmpPtr?", ":"",
		   (field->name)?(field->name):"?",
		   bibFileTokOut[BIBFILE_TOK_EQUAL],
		   bibFileTokOut[BIBFILE_TOK_OPEN]);
    if(tmpPtr)
    {
      (void )AlcFree(tmpPtr);
    }
    if((*eMsg = AlcStrDup(errBuf))== NULL)
    {
      errFlag = BIBFILE_ER_MALLOC;
    }
  }
  return(errFlag);
}

/*!
* \return	Non zero if read fails.
* \ingroup	bibfile
* \brief	Read a character string from the given stream until a
*		character not matched by the given charcter list is
*		read, this character is then put back so that it can be
*		read again. If the skip spaces flag is non zero then
*		leading white space characters will be skipped.
*		Characters in the string may be escaped by the escape
*		character. String length is only limited by the memory
*		available to AlcMalloc.					
* \param	str			Ptr for the string to be read.
* \param	tok			Regular expression for valid
* 					characters within the string.
* \param	skipSp			Skip leading white space
* 					characters if non zero.
* \param	fP			Given input stream.
*/
static BibFileError BibFileStrRead(char **str, BibFileTok tok, int skipSp,
				   FILE *fP)
{
  int		escape,
		endOfStr = 0,
		totIdx = BIBFILE_STRBUF_START,
		newIdx = 0;
  char		old = '\0',
		newC;
  BibFileError	errFlag = BIBFILE_ER_NONE;
  char		*newBuf;

  if((newBuf = (char *)AlcMalloc(BIBFILE_STRBUF_START * sizeof(char))) == NULL)
  {
    errFlag = BIBFILE_ER_MALLOC;
  }
  while((errFlag == BIBFILE_ER_NONE) && (endOfStr == 0))
  {
    escape = (old == *bibFileTokIn[BIBFILE_TOK_ESC]);
    errFlag = BibFileCharRead(&newC, fP, skipSp, escape);
    if(errFlag == BIBFILE_ER_NONE)
    {
      skipSp = 0;
      if(escape)	      /* escaped character replaces escape character */
      {
	--newIdx;
      }
      else
      {
	errFlag = BibFileCharRegEx(newC, tok);
	if(errFlag == BIBFILE_ER_SYNTAX)
	{
	  endOfStr = 1;
	  errFlag = BibFileCharUnread(newC, fP);
	}
      }
    }
    if((errFlag == BIBFILE_ER_NONE) && (endOfStr == 0))
    {
      if((old != *bibFileTokIn[BIBFILE_TOK_ESC]) ||
	 (newC != BIBFILE_EOL)) 		 /* ignore backslash'd eol's */
      {
	*(newBuf + newIdx) = newC;
	++newIdx;
      }
      old = newC;
      if((totIdx - newIdx) < 2)
      {
	totIdx += BIBFILE_STRBUF_INC;
	if((newBuf = (char *)AlcRealloc(newBuf, totIdx)) == NULL)
	{
	  errFlag = BIBFILE_ER_MALLOC;
	}
      }
    }
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    *(newBuf + newIdx) = '\0';
    if((*str = AlcStrDup(newBuf)) == NULL)
    {
      errFlag = BIBFILE_ER_MALLOC;
    }
  }
  if(newBuf)
  {
    (void )AlcFree(newBuf);
  }
  return(errFlag);
}

/*!
* \return	Non zero if match fails.
* \ingroup	bibfile
* \brief	Check for a character that matches the given charcter
*		list. If the skip flag is zero this character is put
*		back so that it can be read again. If the skip spaces
*		flag is non zero then leading white space characters
*		will be skipped.					
* \param	tok			Regular expression token.
* \param	skipFlag		Put the character back into the
* 					input stream if non zero.
* \param	skipSp			Skip leading white space characters
* 					if non zero.
* \param	fP			Given input stream.
*/
static BibFileError BibFileCharMatch(BibFileTok tok, int skipFlag, int skipSp,
				     FILE *fP)
{
  char		newC;
  BibFileError	errFlag = BIBFILE_ER_NONE;
  
  errFlag = BibFileCharRead(&newC, fP, skipSp, 0);
  if(errFlag == BIBFILE_ER_NONE)
  {
    errFlag = BibFileCharRegEx(newC, tok);
  }
  if((skipFlag == 0) &&
     ((errFlag == BIBFILE_ER_NONE) || (errFlag == BIBFILE_ER_SYNTAX)))
  {
    if(errFlag == BIBFILE_ER_NONE)
    {
      errFlag = BibFileCharUnread(newC, fP);
    }
    else if (errFlag == BIBFILE_ER_SYNTAX) 
    {
      (void )BibFileCharUnread(newC, fP);
    }
  }
  return(errFlag);
}

/*!
* \return	Non zero if read fails.
* \ingroup	bibfile
* \brief	Read character from the given stream. If the skip
*		spaces flag is non zero then leading white space
*		characters will be skipped. If escape is non zero then
*		comments will be skipped, this is the only effect of
*		the escape flag.					
* \param	newC			Ptr for the character.
* \param	fP			Given input stream.
* \param	skipSp			Skip leading white space
* 					characters if non zero.
* \param	escape			Escape next character if non zero.
*/
static BibFileError BibFileCharRead(char *newC, FILE *fP, int skipSp,
				    int escape)
{
  int 		tmp;
  BibFileError 	errFlag = BIBFILE_ER_NONE;

  if(skipSp)
  {
    errFlag = BibFileSkipSpaces(fP);
  }
  if((tmp = fgetc(fP)) == EOF)
  {
    errFlag = BIBFILE_ER_EOF;
  }
  else
  {
    while((tmp == *bibFileTokIn[BIBFILE_TOK_COMMENT]) && !escape)
    {
      escape = 0;
      while(((tmp = fgetc(fP)) != EOF) && (tmp != BIBFILE_EOL))
      {
	;
      }
      if(tmp == EOF)
      {
        errFlag = BIBFILE_ER_EOF;
      }
      else if(skipSp)
      {
	errFlag = BibFileSkipSpaces(fP);
      }
      if((errFlag == BIBFILE_ER_NONE) && ((tmp = fgetc(fP)) == EOF))
      {
	errFlag = BIBFILE_ER_SYNTAX;
      }
    }
  }
  *newC = tmp;
  return(errFlag);
}

/*!
* \return	Non zero if cant put it back.
* \ingroup	bibfile
* \brief	Put a character back into the stream from which it was
*		read.							
* \param	newC			The character to put back.
* \param	fP			Given input stream.
*/
static BibFileError BibFileCharUnread(char newC, FILE *fP)
{
  BibFileError	errFlag = BIBFILE_ER_NONE;

  if(ungetc((int )newC, fP) == EOF)
  {
    errFlag = BIBFILE_ER_READ;
  }
  return(errFlag);
}

/*!
* \return	Non zero if match read fails.
* \ingroup	bibfile
* \brief	Match the given character to the given regular
*		expression. The regular expression may either be
*		compilled or a charcater string. 			
* \param	given			Character to match.
* \param	tok			Regular expression for valid
* 					characters to match.
*/
static BibFileError BibFileCharRegEx(char given, BibFileTok tok)
{
  BibFileError	errFlag = BIBFILE_ER_NONE;
  char		tmpStr[2];
#ifdef BIBFILE_REGEXPR_CACHE
  static int	regExpCacheFlags[BIBFILE_TOK_MAX];    /* Initialized to zero */
#ifdef BIBFILE_REGEXPR_REGEX
  static char	*regExpCache[BIBFILE_TOK_MAX];
#endif /* BIBFILE_REGEXPR_REGEX */
#ifdef BIBFILE_REGEXPR_REGEXEC
  static regex_t regExpCache[BIBFILE_TOK_MAX];
#endif /* BIBFILE_REGEXPR_REGEXEC */
#else /* ! BIBFILE_REGEXPR_CACHE */
#ifdef BIBFILE_REGEXPR_REGEX
  char		*regExp = NULL;
#endif /* BIBFILE_REGEXPR_REGEX */
#ifdef BIBFILE_REGEXPR_REGEXEC
  regex_t	regExp;
#endif /* BIBFILE_REGEXPR_REGEXEC */
#endif /* BIBFILE_REGEXPR_CACHE */

  if((tok < BIBFILE_TOK_MIN) || (tok >= BIBFILE_TOK_MAX))
  {
    errFlag = BIBFILE_ER_MALLOC;  /* Not really, but it should never happen! */
  }
#ifdef BIBFILE_REGEXPR_CACHE
  if((errFlag == BIBFILE_ER_NONE) && (regExpCacheFlags[tok] == 0))
#else /* ! BIBFILE_REGEXPR_CACHE */
  if(errFlag == BIBFILE_ER_NONE)
#endif /* BIBFILE_REGEXPR_CACHE */
  {
#ifdef BIBFILE_REGEXPR_REGEX
#ifdef BIBFILE_REGEXPR_CACHE
    if((regExpCache[tok] = regcmp(bibFileTokIn[tok], NULL)) == NULL)
#else /* ! BIBFILE_REGEXPR_CACHE */
    if((regExp = regcmp(bibFileTokIn[tok], NULL)) == NULL)
#endif /* BIBFILE_REGEXPR_CACHE */
#endif /* BIBFILE_REGEXPR_REGEX */
#ifdef BIBFILE_REGEXPR_REGEXEC
#ifdef BIBFILE_REGEXPR_CACHE
    if(regcomp(&regExpCache[tok], bibFileTokIn[tok], REG_EXTENDED|REG_NOSUB))
#else /* ! BIBFILE_REGEXPR_CACHE */
    if(regcomp(&regExp, bibFileTokIn[tok], REG_EXTENDED|REG_NOSUB))
#endif /* BIBFILE_REGEXPR_CACHE */
#endif /* BIBFILE_REGEXPR_REGEXEC */
    {
        errFlag = BIBFILE_ER_MALLOC;
    }
#ifdef BIBFILE_REGEXPR_CACHE
    regExpCacheFlags[tok] = 1;
#endif /* BIBFILE_REGEXPR_CACHE */
  }
  if(errFlag == BIBFILE_ER_NONE)
  {
    tmpStr[0] = given;
    tmpStr[1] = '\0';
#ifdef BIBFILE_REGEXPR_REGEX
#ifdef BIBFILE_REGEXPR_CACHE
    if(regex(regExpCache[tok], tmpStr) == NULL)
#else /* ! BIBFILE_REGEXPR_CACHE */
    if(regex(regExp, tmpStr) == NULL)
#endif /* BIBFILE_REGEXPR_CACHE */
#endif /* BIBFILE_REGEXPR_REGEX */
#ifdef BIBFILE_REGEXPR_REGEXEC
    if(regexec(&regExpCache[tok], tmpStr, (size_t)0, NULL, 0))
#ifdef BIBFILE_REGEXPR_CACHE
    if(regexec(&regExpCache[tok], tmpStr, (size_t)0, NULL, 0))
#else /* ! BIBFILE_REGEXPR_CACHE */
    if(regexec(&regExp, tmpStr, (size_t)0, NULL, 0))
#endif /* BIBFILE_REGEXPR_CACHE */
#endif /* BIBFILE_REGEXPR_REGEXEC */
    {
      errFlag = BIBFILE_ER_SYNTAX;
    }
#ifndef BIBFILE_REGEXPR_CACHE 
#ifdef BIBFILE_REGEXPR_REGEX
    AlcFree(regExp);
#endif /* BIBFILE_REGEXPR_REGEX */
#ifdef BIBFILE_REGEXPR_REGEXEC
    regfree(&regExp);
#endif /* BIBFILE_REGEXPR_REGEXEC */
#endif /* ! BIBFILE_REGEXPR_CACHE */
  }
  return(errFlag);
}

/*!
* \return	Non zero if read fails.
* \ingroup	bibfile
* \brief	Read and skip white space characters.
* \param	fP			Given input stream.
*/
static BibFileError BibFileSkipSpaces(FILE *fP)
{
  int 		tmp;
  BibFileError	errFlag = BIBFILE_ER_NONE;

  while(((tmp = fgetc(fP)) != EOF) &&
	((tmp == ' ') || (tmp == '\t') || (tmp == '\n')))
  {
    ;
  }
  if((tmp == EOF) || (ungetc(tmp, fP) == EOF))
  {
    errFlag = BIBFILE_ER_EOF;
  }
  return(errFlag);
}

/*!
* \return	Non zero if read fails.
* \ingroup	bibfile
* \brief	Replace any special character with ESC+spacial char.
* \param	inString		Given string.
* \param	outString		Destination pointer for the
* 					escaped string.
*/
BibFileError BibFileEscapeRestrictedChar(char *inString, char **outString)
{

  int 		stringLen,
		newIdx = 0,
		oldIdx = 0,
		needsExpansion = 0,
		currentNewLen = 0,
		endOfStr = 0,
		tokNum = 0;
  char		newC,
  		tokChar;
  char 		*newString = NULL;
  BibFileError 	errorF = BIBFILE_ER_NONE;

  if(inString)
  {
    stringLen   = strlen(inString);
    if((newString =(char *)AlcMalloc((stringLen + BIBFILE_STRBUF_START) *
				     sizeof(char))) != NULL)
    {
      currentNewLen = stringLen + BIBFILE_STRBUF_START;
      while((errorF == BIBFILE_ER_NONE) && (endOfStr == 0))
      {
	newC = inString[oldIdx];
	if(newC != '\0')
	{
	  tokNum =0;
	  needsExpansion =0;
	  while((errorF == BIBFILE_ER_NONE) && (endOfStr ==0) &&
	        (tokNum < BIBFILE_RESTRICTED_NUM))
	  {
	    tokChar = *bibFileRestrictedChars[tokNum];
	    if(newC == tokChar)
	    {
	      needsExpansion = 1;
	    }
	    ++tokNum;
	  }
	  if(needsExpansion)
	  {
	    newString[newIdx] = '\\';
	    ++newIdx;
	    newString[newIdx] = newC;
	    ++newIdx;
	  }
	  else
	  {
	    newString[newIdx] = newC;
	    ++newIdx;
	  }
	  if(newIdx >= currentNewLen)
	  {
	    currentNewLen += BIBFILE_STRBUF_INC;
	    if((newString = (char *)AlcRealloc(newString,
					       currentNewLen)) == NULL)
	    {
	      errorF = BIBFILE_ER_MALLOC;
	    }
	  }
	}
	else
	{
	  endOfStr = 1;
	}
	++oldIdx;
      }
    }
    if(errorF == BIBFILE_ER_NONE)
    {
      newString[newIdx] = '\0';
      if((*outString = AlcStrDup(newString)) == NULL)
      {
	errorF = BIBFILE_ER_MALLOC;
      }
    }
    if(newString)
    {
      (void )AlcFree(newString);
    }
  }
  return(errorF);
}

/*!
* \return	Non zero if read fails.
* \ingroup	bibfile
* \brief	Replace any special caracter with ESC+spacial char.
* \param	inString		Original string.
* \param	outString		The escaped string.
*/
BibFileError BibFileUnEscapeRestrictedChar(char *inString, char **outString)
{
  char 		*newString = NULL;
  char 		newChar;
  int		stringLen,
 		newIdx = 0,
		oldIdx = 0;
  BibFileError	errorF = BIBFILE_ER_NONE;

  if(inString)
  {
    stringLen   = strlen(inString);  
    if((newString = (char *)AlcMalloc((stringLen) * sizeof(char))) != NULL) 
    {
      while((newChar = inString[oldIdx]) != '\0')
      {
	if(newChar != '\\')
	{
	  newString[newIdx] = newChar;
	  ++newIdx;
	}
	++oldIdx;
      }
      newString[newIdx] = '\0';
      if((*outString = AlcStrDup(newString)) == NULL)
      {
	errorF = BIBFILE_ER_MALLOC;
      }
    }
    if(newString)
    {
      (void )AlcFree(newString);
    }
  }
  return(errorF);
}
