#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlzExtFF/WlzExtFFAm.c
* \author       Bill Hill
* \date         October 2002
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
* \brief	Functions for reading and writting Woolz objects to
* 		the Amira lattice file format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <Wlz.h>
#include <WlzExtFF.h>

#define WLZ_EFF_AM_REC_MAXC	(1024)
#define WLZ_EFF_AM_LABEL_MAXC	(32)

static WlzEffAmHead 		*WlzEffAmNewHead(
				  void);
static void			WlzEffAmFreeHead(
				  WlzEffAmHead *head);
static WlzErrorNum		WlzEffAmReadAToken(
				  FILE *fP,
				  const char *sep,
				  char *buf,
				  const int len);
static WlzErrorNum 		WlzEffAmSkip(
				  FILE *fP,
				  const char *skip);
static WlzErrorNum 		WlzEffAmSkipEOL(
				  FILE *fP);
static WlzErrorNum 		WlzEffAmReadAValI(
				  FILE *fP,
				  const char *sep,
				  char *buf,
				  const int bufLen,
				  int *dVal);
static WlzErrorNum 		WlzEffAmReadAValD(
				  FILE *fP,
				  const char *sep,
				  char *buf,
				  const int bufLen,
				  double *dVal);
static WlzErrorNum 		WlzEffAmReadAndCheckAToken(
				  FILE *fP,
				  const char *sep,
				  char *buf,
				  const int bufLen,
				  char *gWord);
static WlzErrorNum 		WlzEffAmReadHead(
				  FILE *fP,
				  WlzEffAmHead *head);
static WlzErrorNum 		WlzEffAmSeekLabel(
				  FILE *fP,
				  int gLabel);
static WlzErrorNum		WlzEffAmReadBBox(
				  FILE *fP,
				  char *buf,
				  const int bufLen,
				  WlzDBox3 *bBox);
static WlzErrorNum		WlzEffAmCheckContent(
				  FILE *fP,
				  char *buf,
				  const int bufLen);
static WlzErrorNum 		WlzEffAmReadCoordType(
				  FILE *fP,
				  char *buf,
				  const int bufLen,
				  WlzEffAmCoordType *type);
static WlzErrorNum 		WlzEffAmReadLatticeDatType(
				  FILE *fP,
				  char *buf,
				  const int bufLen,
				  WlzEffAmDatType *type);
static WlzErrorNum 		WlzEffAmReadLatticeType(
				  FILE *fP,
				  char *buf,
				  const int bufLen,
				  WlzEffAmLatType *type,
				  int *unknown);
static WlzErrorNum 		WlzEffAmReadLatticeLabel(
				  FILE *fP,
				  char *buf,
				  const int bufLen,
				  int skip,
				  WlzEffAmHead *head);
static WlzErrorNum		WlzEffAmReadMaterials(
				  FILE *fP,
				  char *buf,
				  const int bufLen,
				  WlzEffAmHead *head);
static WlzErrorNum 		WlzEffAmReadAndSkipValues(
				  FILE *fP,
				  char *buf,
				  const int bufLen);
static WlzErrorNum 		WlzEffAmReadArray3D(
				  FILE *fP,
				  void ***data,
				  WlzEffAmHead *head,
				  WlzGreyType gType);
static void			WlzEffAmBufDecodeHXByteRLEUByte(
				  UBYTE *dst,
				  UBYTE *src,
				  int dstCnt,
				  int srcCnt);
static void			WlzEffAmBufDecodeHXByteRLEShort(
				  short *dst,
				  short *src,
				  int dstCnt,
				  int srcCnt);
static WlzObject 		*WlzEffAmReadArray3DToObj(
				  FILE *fP,
				  WlzEffAmHead *head,
				  WlzErrorNum *dstErr);
static WlzCompoundArray 	*WlzEffAmSplitLabelObj(
				  WlzObject *gObj,
				  WlzEffAmHead *head,
				  WlzErrorNum *dstErr);
static WlzPropertyList 		*WlzEffAmMakeMaterialPropList(
				  WlzEffAmMaterial *mat,
				  WlzErrorNum *dstErr);


/*!
* \return	New Woolz object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given file using the
*		Amira lattice (.am file) format.
* \param	fP			Input file pointer.
* \param	split			Split label lattice files into a
*					compound object of domains if non
*					zero.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzObject	*WlzEffReadObjAm(FILE *fP, int split, WlzErrorNum *dstErr)
{
  WlzIVertex3	org;
  WlzDVertex3	voxSz;
  WlzEffAmHead	*head;
  WlzGreyType	gType;
  WlzObject	*tObj,
  		*obj = NULL;
  void		***data = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((head = WlzEffAmNewHead()) == NULL)
  {
    errNum =  WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadHead(fP, head);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check type of lattice file can be read and call appropriate function. */
    switch(head->latType)
    {
      case WLZEFF_AM_LATTYPE_DATA:
	obj = WlzEffAmReadArray3DToObj(fP, head, &errNum);
	break;
      case WLZEFF_AM_LATTYPE_LABELS:
	obj = WlzEffAmReadArray3DToObj(fP, head, &errNum);
	if(split && (errNum == WLZ_ERR_NONE))
	{
	  tObj = (WlzObject *)WlzEffAmSplitLabelObj(WlzAssignObject(obj, NULL),
	  				head, &errNum);
	  WlzFreeObj(obj);
	  obj = tObj;
	}
        break;
      default:
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
    }
  }
  WlzEffAmFreeHead(head);
  if(data)
  {
    if(obj)
    {
      **data = NULL;
    }
    Alc3Free(data);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(obj)
    {
      (void )WlzFreeObj(obj);
      obj = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file using the
*		Amira lattice (.am file) format.
* \param	fP			Output file pointer
* \param	obj			Object to be written.
*/
WlzErrorNum	WlzEffWriteObjAm(FILE *fP, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Parses the Amira file header filling in the fields
*		of the header structure. 
* \param	fP			Input file ready to read first
*					datum.
* \param	head			Pointer to an Amira file header
*					data structure.
*/
static WlzErrorNum WlzEffAmReadHead(FILE *fP, WlzEffAmHead *head)
{
  int		skip,
  		labelId,
  		done = 0;
  long		fRecPos;
  WlzEffAmDatType latDatType;
  int		tok[8];
  char		tokBuf[WLZ_EFF_AM_REC_MAXC];
  const int	tokBufMax = WLZ_EFF_AM_REC_MAXC;
  const char   	tokSepWS[] = " \t\n\r\f\v",
  	   	tokSepWSC[] = " \t\n,\r\f\v",
		tokSepWSQ[] = " \t\n\"\r\f\v";
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Read identifier, data format and version numbers. */
  errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
  if((errNum == WLZ_ERR_NONE) && strcmp(tokBuf, "#"))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
  }
  if((errNum == WLZ_ERR_NONE) && strcmp(tokBuf, "AmiraMesh"))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tok[0] = WLZEFF_AM_DIM_NONE;
    (void )WlzStringMatchValue((int *)&(tok[0]), tokBuf,
			       "2D", WLZEFF_AM_DIM_2,
			       "3D", WLZEFF_AM_DIM_3,
			       NULL);
    if(tok[0] == WLZEFF_AM_DIM_NONE)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      head->dim = tok[0];
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tok[0] = WLZEFF_AM_FMT_NONE;
    (void )WlzStringMatchValue((int *)&(tok[0]), tokBuf,
			       "ASCII", WLZEFF_AM_FMT_ASCII,
			       "BINARY", WLZEFF_AM_FMT_BINARY,
			       NULL);
    if(tok[0] == WLZEFF_AM_TOKEN_NONE)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      head->fmt = tok[0];
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(sscanf(tokBuf, "%d.%d",
              &(head->versionMajor), &(head->versionMinor)) != 2)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Read rest of header. */
  while((errNum == WLZ_ERR_NONE) && (done == 0))
  {
    errNum = WlzEffAmSkip(fP, tokSepWS);
    fRecPos = ftell(fP);
    tok[0] = WLZEFF_AM_TOKEN_NONE;
    errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzStringMatchValue((int *)&(tok[0]), tokBuf,
				 "#", WLZEFF_AM_TOKEN_HASH,
				 "define", WLZEFF_AM_TOKEN_DEFINE,
				 "Parameters", WLZEFF_AM_TOKEN_PARAMETERS,
				 "Lattice", WLZEFF_AM_TOKEN_LATTICE,
				 "Materials", WLZEFF_AM_TOKEN_MATERIALS,
				 NULL);
      switch(tok[0])
      {
	case WLZEFF_AM_TOKEN_HASH:
	  /* Comment. */
	  errNum = WlzEffAmSkipEOL(fP);
	  break;
	case WLZEFF_AM_TOKEN_DEFINE:
	  errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tok[1] = WLZEFF_AM_TOKEN_NONE;
	    (void )WlzStringMatchValue((int *)&(tok[1]), tokBuf,
				       "Lattice", WLZEFF_AM_TOKEN_LATTICE,
				       NULL);
	    if(tok[1] == WLZEFF_AM_TOKEN_NONE)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sscanf(tokBuf, "%d", &(head->latSize.vtX)) != 1)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sscanf(tokBuf, "%d", &(head->latSize.vtY)) != 1)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf, tokBufMax);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sscanf(tokBuf, "%d", &(head->latSize.vtZ)) != 1)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	    }
	  }
	  break;
	case WLZEFF_AM_TOKEN_PARAMETERS:
	  errNum = WlzEffAmReadAndCheckAToken(fP, tokSepWS, tokBuf, 
					      tokBufMax, "{");
	  /* Read parameters. */
	  do
	  {
	    errNum = WlzEffAmReadAToken(fP, tokSepWS, tokBuf,
	    			        tokBufMax);
	    if((errNum == WLZ_ERR_NONE) && strcmp(tokBuf, "}"))
	    {
              tok[1] = WLZEFF_AM_TOKEN_NONE;
	      (void )WlzStringMatchValue((int *)&(tok[1]), tokBuf,
				   "BoundingBox", WLZEFF_AM_TOKEN_BOUNDINGBOX,
				   "Content", WLZEFF_AM_TOKEN_CONTENT,
				   "CoordType", WLZEFF_AM_TOKEN_COORDTYPE,
				   "ImageData", WLZEFF_AM_TOKEN_IMAGEDATA,
				   "Materials", WLZEFF_AM_TOKEN_MATERIALS,
				   "Seeds", WLZEFF_AM_TOKEN_SEEDS,
				   "Limits", WLZEFF_AM_TOKEN_LIMITS,
				   NULL);
	      switch(tok[1])
	      {
	        case WLZEFF_AM_TOKEN_BOUNDINGBOX:
		  errNum = WlzEffAmReadBBox(fP, tokBuf, tokBufMax,
					    &(head->bBox));
		  *tokBuf = '\0';
		  break;
	        case WLZEFF_AM_TOKEN_CONTENT:
		  errNum = WlzEffAmCheckContent(fP, tokBuf, tokBufMax);
		  *tokBuf = '\0';
		  break;
	        case WLZEFF_AM_TOKEN_COORDTYPE:
		  errNum = WlzEffAmReadCoordType(fP, tokBuf, tokBufMax,
		  				 &(head->coordType));
		  *tokBuf = '\0';
		  break;
	        case WLZEFF_AM_TOKEN_MATERIALS:
		  errNum = WlzEffAmReadMaterials(fP, tokBuf, tokBufMax, head);
		  *tokBuf = '\0';
		  break;
		case WLZEFF_AM_TOKEN_IMAGEDATA:
		  errNum = WlzEffAmReadAToken(fP, tokSepWSQ, tokBuf,
		  			      tokBufMax);
		  if(errNum == WLZ_ERR_NONE)
		  {
		    errNum = WlzEffAmSkip(fP, tokSepWSQ);
		  }
		  if(errNum == WLZ_ERR_NONE)
		  {
		    if((head->imageData = AlcStrDup(tokBuf)) == NULL)
		    {
		      errNum = WLZ_ERR_MEM_ALLOC;
		    }
		  }
		  *tokBuf = '\0';
		  break;
		case WLZEFF_AM_TOKEN_SEEDS:
		  errNum = WlzEffAmReadAndSkipValues(fP, tokBuf, tokBufMax);
		  *tokBuf = '\0';
		  break;
		case WLZEFF_AM_TOKEN_LIMITS:
		  errNum = WlzEffAmReadAndSkipValues(fP, tokBuf, tokBufMax);
		  *tokBuf = '\0';
		  break;
		default:
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		  break;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzEffAmSkip(fP, tokSepWSC);
	    }
	  } while((errNum == WLZ_ERR_NONE) && strcmp(tokBuf, "}"));
	  break;
	case WLZEFF_AM_TOKEN_LATTICE:
	  errNum = WlzEffAmReadAndCheckAToken(fP, tokSepWS, tokBuf,
					      tokBufMax, "{");
	  /* Read lattice spec. */
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzEffAmReadLatticeDatType(fP, tokBuf, tokBufMax,
	    					&latDatType);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzEffAmReadLatticeType(fP, tokBuf, tokBufMax,
					     &(head->latType), &skip);
	    if(skip == 0)
	    {
	      head->latDatType = latDatType;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzEffAmReadAndCheckAToken(fP, tokSepWS, tokBuf,
						tokBufMax, "}");
	  }
	  /* Read lattice label. */
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzEffAmReadLatticeLabel(fP, tokBuf, tokBufMax, skip,
	    				      head);
	  }
	  break;
	case WLZEFF_AM_TOKEN_MATERIALS:
	  errNum = WlzEffAmReadMaterials(fP, tokBuf, tokBufMax, head);
	  break;
	default:
	  /* Check for label. */
	  if(sscanf(tokBuf, "@%d", &labelId) == 1)
	  {
	    /* Rewind back before label. */
	    if(fseek(fP, fRecPos, SEEK_SET) == -1)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      done = 1;
	    }
	  }
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	New Amira file header data structure or NULLon error.
* \ingroup	WlzExtFF
* \brief	Allocates a new Amira file header data structure with
*		all fields cleared.
* \param        <void>
*/
static WlzEffAmHead *WlzEffAmNewHead(void)
{
  WlzEffAmHead	*head;

  head = AlcCalloc(1, sizeof(WlzEffAmHead));
  return(head);
}

/*!
* \return	<void>
* \ingroup	WlzExtFF
* \brief	Free's entries in the Amira file header data structure
*		before freeing the header data structure itself.
* \param	head			The given file header data structure.
*/
static void	WlzEffAmFreeHead(WlzEffAmHead *head)
{
  WlzEffAmMaterial *e0;

  AlcFree(head->imageData);
  if(head)
  {
    e0 = head->materials;
    while(e0)
    {
      AlcFree(e0->name);
      e0 = e0->next;
    }
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads the next ascii token from the given file, using the
*		token seperators of the given string.
* \param	fP			Given file.
* \param	sep			String of token seperators.
* \param	buf			Buffer for workspace, an error will
*					be returned if buffer not big enough
*					for word and terminating NULL.
* \param	len			Length of buffer which must be greater
* 					than one.
*/
static WlzErrorNum WlzEffAmReadAToken(FILE *fP, const char *sep,
				      char *buf, const int len)
{
  int 		c,
  		idx,
		lim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  lim = len - 1;
  do
  {
    c = getc(fP);
  } while((c != EOF) && strchr(sep, c));
  if(c == EOF)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    idx = 0;
    do
    {
      *(buf + idx++) = c;
      c = getc(fP);
      if((c == EOF) || (idx >= lim))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    } while((errNum == WLZ_ERR_NONE) && (strchr(sep, c) == NULL));
    *(buf + idx) = '\0';
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads a double value encoded as an ascii token from
*		the given file.
* \param	fP			Given file.
* \param	sep			Token seperator.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	dVal			Destination pointer for the double
*					value, must not be NULL.
*/
static WlzErrorNum WlzEffAmReadAValD(FILE *fP, const char *sep,
				     char *buf, const int bufLen,
				     double *dVal)
{
  WlzErrorNum	errNum;

  errNum = WlzEffAmReadAToken(fP, sep, buf, bufLen);
  if(errNum == WLZ_ERR_NONE)
  {
    if(sscanf(buf, "%lg", dVal) != 1)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads an int value encoded as an ascii token from
*		the given file.
* \param	fP			Given file.
* \param	sep			Token seperator.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	dVal			Destination pointer for the double
*					value, must not be NULL.
*/
static WlzErrorNum WlzEffAmReadAValI(FILE *fP, const char *sep,
				     char *buf, const int bufLen,
				     int *dVal)
{
  WlzErrorNum	errNum;

  errNum = WlzEffAmReadAToken(fP, sep, buf, bufLen);
  if(errNum == WLZ_ERR_NONE)
  {
    if(sscanf(buf, "%d", dVal) != 1)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads an ascii token from the given file and then checks
*		that it is the same as the given word. If it's not an
*		error is returned.
* \param	fP			Given file.
* \param	sep			Token seperator.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	gWord			Given word to check for.
*/
static WlzErrorNum WlzEffAmReadAndCheckAToken(FILE *fP, const char *sep,
					char *buf, const int bufLen,
					char *gWord)
{
  WlzErrorNum	errNum;

  if((errNum = WlzEffAmReadAToken(fP, sep, buf, bufLen)) == WLZ_ERR_NONE)
  {
    if(strcmp(buf, gWord))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Seek forward from the current file position until the
*		given label has been read. Label in file is terminated
*		by a white space character.
* \param	fP			Given file.
* \param	gLabel			Given label to seek.
*/
static WlzErrorNum WlzEffAmSeekLabel(FILE *fP, int gLabel)
{
  int		c,
		idx,
		lim,
  		found = 0,
		fLabel;
  char		buf[WLZ_EFF_AM_LABEL_MAXC];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	bufMax = WLZ_EFF_AM_LABEL_MAXC;
  const char   	tokSepWS[] = " \t\n\r\f\v";

  do
  {
    do
    {
      c = getc(fP);
    } while((c != EOF) && (c != '@'));
    if(c == EOF)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      idx = 0;
      lim = bufMax - 1;
      do
      {
	*(buf + idx++) = c;
	c = getc(fP);
	if((c == EOF) || (idx >= lim))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      } while((errNum == WLZ_ERR_NONE) && isdigit(c));
      *(buf + idx) = '\0';
      if(errNum == WLZ_ERR_NONE)
      {
	if(sscanf(buf, "@%d", &fLabel) != 1)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else if (fLabel = gLabel)
	{
	  found = 1;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      while(c != '\n')
      {
        c = getc(fP);
      }
    }
  } while((found == 0) && (errNum == WLZ_ERR_NONE));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Read a double precision bounding box from the file as ascii,
*		with the order being xmin, xmax, ymin, ymax, zmin and zmax.
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	bBox			Destination pointer for the bounding
* 					box, must not be NULL.
*/
static WlzErrorNum WlzEffAmReadBBox(FILE *fP, char *buf, const int bufLen,
				    WlzDBox3 *bBox)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char   tokSepWS[] = " \t\n,\r\f\v";

  errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen, &(bBox->xMin));
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen, &(bBox->xMax));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen, &(bBox->yMin));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen, &(bBox->yMax));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen, &(bBox->zMin));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen, &(bBox->zMax));
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Skip past given characters.
* \param	fP			Given file.
* \param	skip			Characters to skip.
*/
static WlzErrorNum WlzEffAmSkip(FILE *fP, const char *skip)
{
  int		c;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  do
  {
    c = getc(fP);
  } while((c != EOF) && strchr(skip, c));
  if(c == EOF)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    (void )(ungetc(c, fP));
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Skip the current line.
* \param	fP			Given file.
* \param	skip			Characters to skip.
*/
static WlzErrorNum WlzEffAmSkipEOL(FILE *fP)
{
  int		c;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  do
  {
    c = getc(fP);
  } while((c != '\n') && (c != EOF));
  if(c == EOF)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads the lattice paremeters content record and then ignores
*		it.
* 		A valid contense record looks like:
*		  "361x391x154 byte, uniform coordinates"
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
*/
static WlzErrorNum WlzEffAmCheckContent(FILE *fP, char *buf, const int bufLen)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	tokSepWS[] = " \t\n\r\f\v",
  		tokSepDQ[] = "\"";

  errNum = WlzEffAmSkip(fP, tokSepWS);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAToken(fP, tokSepDQ, buf, bufLen);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads the lattice paremeters coordinate type record and
*		then checks taht it is of type uniform.
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
*/
static WlzErrorNum WlzEffAmReadCoordType(FILE *fP, char *buf, const int bufLen,
					 WlzEffAmCoordType *type)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	tokSepWS[] = " \t\n\r\f\v",
  		tokSepDQ[] = "\"";

  errNum = WlzEffAmSkip(fP, tokSepWS);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAToken(fP, tokSepDQ, buf, bufLen);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *type = WLZEFF_AM_COORDTYPE_NONE;
    (void )WlzStringMatchValue((int *)type, buf,
			       "uniform", WLZEFF_AM_COORDTYPE_UNITFORM,
			       NULL);
    if(*type == WLZEFF_AM_COORDTYPE_NONE)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads and parses the lattice data type. This is either
*		"byte" -> WLZEFF_AM_DATTYPE_BYTE or
*		"short" -> WLZEFF_AM_DATTYPE_SHORT
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	type			Destination pointer for type.
*/
static WlzErrorNum WlzEffAmReadLatticeDatType(FILE *fP, char *buf,
					const int bufLen,
					WlzEffAmDatType *type)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	tokSepWS[] = " \t\n\r\f\v";

  errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
  if(errNum == WLZ_ERR_NONE)
  {
    *type = WLZEFF_AM_DATTYPE_NONE;
    (void )WlzStringMatchValue((int *)type, buf,
			       "byte", WLZEFF_AM_DATTYPE_BYTE,
			       "short", WLZEFF_AM_DATTYPE_SHORT,
			       NULL);
    if(*type == WLZEFF_AM_DATTYPE_NONE)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads and parses the lattice data type. This is either
*		"Data" -> WLZEFF_AM_LATTYPE_DATA or
*		"Labels" -> WLZEFF_AM_LATTYPE_LABELS.
*		any other lattice types cause the unknown flag to be set.
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	dstType			Destination pointer for type.
* \param	unknown			Destination for unknown flag.
*/
static WlzErrorNum WlzEffAmReadLatticeType(FILE *fP, char *buf, 
					const int bufLen,
					WlzEffAmLatType *dstType,
					int *unknown)
{
  int		type;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	tokSepWS[] = " \t\n\r\f\v";

  *unknown = 0;
  errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
  if(errNum == WLZ_ERR_NONE)
  {
    type = (int )WLZEFF_AM_LATTYPE_NONE;
    (void )WlzStringMatchValue(&type, buf,
			       "Data", WLZEFF_AM_LATTYPE_DATA,
			       "Labels", WLZEFF_AM_LATTYPE_LABELS,
			       NULL);
    if(type == (int )WLZEFF_AM_LATTYPE_NONE)
    {
      *unknown = 1;
    }
    else
    {
      *dstType = type;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads and parses the lattice label which may include an
*		encoding method and an encoded lattice data byte count.
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	skip			Skip this lattice if non zero,
*					parsing it but not putting data into
*					lattice header data structure.
* \param	head			The Amira lattice header data
*					structure.
*/
static WlzErrorNum WlzEffAmReadLatticeLabel(FILE *fP, char *buf, 
					const int bufLen, int skip,
					WlzEffAmHead *head)
{
  WlzErrorNum	errNum;
  const char	tokSepWS[] = " \t\n\r\f\v";

  errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
  if(skip == 0)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      if(sscanf(buf, "@%d(HxByteRLE,%d)", &(head->latLabel),
	    &(head->latBytes)) == 2)
      {
	head->latComp = WLZEFF_AM_LATCOMP_HXBYTERLE;
      }
      else if(sscanf(buf, "@%d", &(head->latLabel)) == 1)
      {
	head->latComp = WLZEFF_AM_LATCOMP_NONE;
	head->latBytes = head->latSize.vtX * head->latSize.vtY *
			 head->latSize.vtZ;
      }
    }
    else
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads the materials list from the file into the file header
*		data structure.
*		There are two formats for materials lists:
*		\verbatim
		  {
		    { Id 0, Name "Exterior" )
		    { Id 4, Name "Somites" }
		  }
		\endverbatim
*		and		
*		\verbatim
		  {
		    Exterior {
		    }
		    Inside {
		    }
		    Somites {
		      Id 4,
		      Color 0.2 0.8 0.4
		    }
		    Ventricles
		  }
		\endverbatim
*		Dispite id's being given in the file they are ignored and
*		the id's are given values from 0, incrementing by 1 with
*		each material parsed.
*		The default colour is 0.0 0.0 0.0 for the first material
*		and 1.0 1.0 1.0 for all others.
*		This function has been written to parse either of these
*		materials formats, to generate a materials list in the
*		header data structure.
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
* \param	head			The Amira lattice header data
*					structure.
*/
static WlzErrorNum WlzEffAmReadMaterials(FILE *fP, char *buf, const int bufLen,
				  WlzEffAmHead *head)
{

  int		tmpI,
		matId = 0;
  WlzEffAmToken	tok;
  WlzEffAmMaterial *lastMat = NULL,
  		*newMat = NULL;
  WlzErrorNum   errNum;
  const char	tokSepDQ[] = "\"",
  		tokSepWS[] = " \t\n\r\f\v",
            	tokSepWSC[] = " \t\n,\r\f\v",
            	tokSepWSCQ[] = " \t\n,\"\r\f\v",
            	tokSepWSQ[] = " \t\n\"\r\f\v";

  errNum =  WlzEffAmReadAndCheckAToken(fP, tokSepWS, buf, bufLen, "{");
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
  }
  do
  {
    if((errNum == WLZ_ERR_NONE) && (strcmp(buf, "}")))
    {
      /* Create a new material and set the defaults. */
      if((newMat = (WlzEffAmMaterial *)
      		   AlcCalloc(1, sizeof(WlzEffAmMaterial))) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else if(head->matCount > 0)
      {
        newMat->color[0] = newMat->color[1] = newMat->color[2] = 1.0;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(strcmp(buf, "{"))
      {
	if((newMat->name = AlcStrDup(buf)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
	}
      }
    }
    if((errNum == WLZ_ERR_NONE) && strcmp(buf, "{") && strcmp(buf, "}"))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    if((errNum == WLZ_ERR_NONE) && strcmp(buf, "}"))
    {
      do
      {
	errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
	if((errNum == WLZ_ERR_NONE) && (strcmp(buf, "}")))
	{
	  tok = WLZEFF_AM_TOKEN_NONE;
	  (void )WlzStringMatchValue((int *)&tok, buf,
				     "Color", WLZEFF_AM_TOKEN_COLOR,
				     "Id", WLZEFF_AM_TOKEN_ID,
				     "Name", WLZEFF_AM_TOKEN_NAME,
				     NULL);
	  switch(tok)
	  {
	    case WLZEFF_AM_TOKEN_COLOR:

	      errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen,
	      				 &(newMat->color[0]));
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen,
					   &(newMat->color[1]));
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = WlzEffAmReadAValD(fP, tokSepWS, buf, bufLen,
					   &(newMat->color[2]));
	      }
	      break;
	    case WLZEFF_AM_TOKEN_ID:
	      errNum = WlzEffAmReadAValI(fP, tokSepWSC, buf, bufLen, &tmpI);
	      if(errNum == WLZ_ERR_NONE)
	      {
                errNum = WlzEffAmSkip(fP, tokSepWSC);
	      }
	      break;
	    case WLZEFF_AM_TOKEN_NAME:
	      errNum = WlzEffAmSkip(fP, tokSepWS);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        errNum = WlzEffAmReadAToken(fP, tokSepDQ, buf, bufLen);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        errNum = WlzEffAmSkip(fP, tokSepWSCQ);
	      }
	      if((errNum == WLZ_ERR_NONE) &&
	      	 ((newMat->name = AlcStrDup(buf)) == NULL))
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	  }
	}
      }while((errNum == WLZ_ERR_NONE) && strcmp(buf, "}"));
      *buf = '\0';
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Insert new material into list keeping the order of the materials
       * in which they were read from the file. */
      newMat->id = (head->matCount)++;
      if(lastMat == NULL)
      {
	head->materials = newMat;
      }
      else
      {
        lastMat->next = newMat;
	newMat->prev = lastMat;
      }
      lastMat = newMat;
      newMat = NULL;
    }
    if((errNum == WLZ_ERR_NONE) && strcmp(buf, "}"))
    {
      errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
    }
  } while((errNum == WLZ_ERR_NONE) && strcmp(buf, "}"));
  if(newMat)
  {
    AlcFree(newMat->name);
    AlcFree(newMat);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads and skips everything between the '{' '}' next pair.
* \param	fP			Given file.
* \param	buf			Buffer for workspace.
* \param	bufLen			Buffer length.
*/
static WlzErrorNum WlzEffAmReadAndSkipValues(FILE *fP,
					     char *buf, const int bufLen)
{

  int		idN;
  WlzEffAmToken	tok;
  WlzErrorNum   errNum;
  const char	tokSepWS[] = " \t\n\r\f\v";

  idN = 1;
  errNum =  WlzEffAmReadAndCheckAToken(fP, tokSepWS, buf, bufLen, "{");
  while((errNum == WLZ_ERR_NONE) && (idN > 0))
  {
    errNum = WlzEffAmReadAToken(fP, tokSepWS, buf, bufLen);
    if(strcmp(buf, "{") == 0)
    {
      ++idN;
    }
    else if(strcmp(buf, "}") == 0)
    {
      --idN;
    }
  }
  return(errNum);
}

/*!
* \return	New Woolz object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given file using the
*		previously parsed Amira lattice file header.
* \param	fP			Input file pointer.
* \param	head			Amira lattice file header.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
static WlzObject *WlzEffAmReadArray3DToObj(FILE *fP, WlzEffAmHead *head,
				 	WlzErrorNum *dstErr)
{
  WlzIVertex3	org;
  WlzDVertex3	voxSz;
  WlzGreyType	gType;
  WlzObject	*obj = NULL;
  void		***data = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((head->dim != WLZEFF_AM_DIM_3) ||
     (head->fmt != WLZEFF_AM_FMT_BINARY) ||
     (head->latSize.vtX < 1) ||
     (head->latSize.vtY < 1) ||
     (head->latSize.vtZ < 1) ||
     ((voxSz.vtX = head->bBox.xMax - head->bBox.xMin) < DBL_EPSILON) ||
     ((voxSz.vtY = head->bBox.yMax - head->bBox.yMin) < DBL_EPSILON) ||
     ((voxSz.vtZ = head->bBox.zMax - head->bBox.zMin) < DBL_EPSILON))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(head->latComp)
    {
      case WLZEFF_AM_LATCOMP_NONE: /* FALLTHROUGH */
      case WLZEFF_AM_LATCOMP_HXBYTERLE:
        break;
      default:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmSeekLabel(fP, head->latLabel);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(head->latDatType)
    {
      case WLZEFF_AM_DATTYPE_BYTE:
	gType = WLZ_GREY_UBYTE;
	if(AlcChar3Malloc((char ****)&data,
			  head->latSize.vtZ, head->latSize.vtY,
			  head->latSize.vtX) != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZEFF_AM_DATTYPE_SHORT:
	gType = WLZ_GREY_SHORT;
	if(AlcShort3Malloc((short ****)&data,
			   head->latSize.vtZ, head->latSize.vtY,
			   head->latSize.vtX) != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffAmReadArray3D(fP, data, head, gType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute voxel size and origin.
     * The Amira bounding box is for voxel centres and NOT
     * voxel boundaries. */
    if(head->latSize.vtX > 1)
    {
      voxSz.vtX /= (double )(head->latSize.vtX - 1);
    }
    if(head->latSize.vtY > 1)
    {
      voxSz.vtY /= (double )(head->latSize.vtY - 1);
    }
    if(head->latSize.vtZ > 1)
    {
      voxSz.vtZ /= (double )(head->latSize.vtZ - 1);
    }
    org.vtX = head->bBox.xMin / voxSz.vtX;
    org.vtY = head->bBox.yMin / voxSz.vtY;
    org.vtZ = head->bBox.zMin / voxSz.vtZ;
    obj = WlzFromArray3D(data, head->latSize, org,
			 gType, gType, 0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->domain.p->voxel_size[0] = voxSz.vtX;
    obj->domain.p->voxel_size[1] = voxSz.vtY;
    obj->domain.p->voxel_size[2] = voxSz.vtZ;
  }
  if(data)
  {
    if(obj)
    {
      **data = NULL;
    }
    Alc3Free(data);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(obj)
    {
      (void )WlzFreeObj(obj);
      obj = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Read the block of data at the current file position into the
*		given Alc style array which has already been allocated.
* \param	fP			File pointer at ready to read 
*					first datum.
* \param	data			Alc style 3D array.
* \param	head			Amira lattice file header.
* \param	gType			Grey type to be read.
*/
static WlzErrorNum WlzEffAmReadArray3D(FILE *fP, void ***data,
				WlzEffAmHead *head, WlzGreyType gType)
{
  int		nDst;
  size_t	datumSz;
  void		*buf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_UBYTE:
      datumSz = sizeof(UBYTE);
      break;
    case WLZ_GREY_SHORT:
      datumSz = sizeof(short);
      break;
    default:
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(head->latComp)
    {
      case WLZEFF_AM_LATCOMP_NONE:
	if(fread(**data, datumSz, head->latBytes, fP) != head->latBytes)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	break;
      case WLZEFF_AM_LATCOMP_HXBYTERLE:
	if((buf = AlcMalloc(datumSz * head->latBytes)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else if(fread(buf, datumSz, head->latBytes, fP) != head->latBytes)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  nDst = head->latSize.vtX * head->latSize.vtY * head->latSize.vtZ;
	  switch(gType)
	  {
	    case WLZ_GREY_UBYTE:
	      WlzEffAmBufDecodeHXByteRLEUByte(**(UBYTE ***)data, (UBYTE *)buf,
	      				      nDst, head->latBytes);
	      break;
	    case WLZ_GREY_SHORT:
	      WlzEffAmBufDecodeHXByteRLEShort(**(short ***)data, (short *)buf,
	      				      nDst, head->latBytes);
	      break;
	  }
	}
	AlcFree(buf);
	break;
      default:
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzExtFF
* \brief	Decodes the run length encoded source buffer into the
*		destination buffer, where the buffers are both of unsigned
*		bytes.
* \param	dst			Destination buffer.
* \param	src			Source buffer.
* \param	dstCnt			Destination buffer size.
* \param	srcCnt			Source buffer size.
*/
static void	WlzEffAmBufDecodeHXByteRLEUByte(UBYTE *dst, UBYTE *src,
						int dstCnt, int srcCnt)
{
  int		vCnt;
  UBYTE 	vVal;

  vVal = *src++;
  while((vCnt = vVal & 0x7f) != 0)
  {
    if(vVal & 0x80)
    {
      while(vCnt--)
      {
	*dst++ = *src++;
      }
    }
    else
    {
      vVal = *src++;
      while(vCnt--)
      {
	*dst++ = vVal;
      }
    }
    vVal = *src++;
  } 
}

/*!
* \return	<void>
* \ingroup	WlzExtFF
* \brief	Decodes the run length encoded source buffer into the
*		destination buffer, where the buffers are both of short
*		values and the sizes are the number of shorts.
* \param	dst			Destination buffer.
* \param	src			Source buffer.
* \param	dstCnt			Destination buffer size.
* \param	srcCnt			Source buffer size.
*/
static void	WlzEffAmBufDecodeHXByteRLEShort(short *dst, short *src,
						int dstCnt, int srcCnt)
{
  int		vCnt;
  short 	vVal;

  vVal = *src++;
  while((vCnt = vVal & 0x7f) != 0)
  {
    if(vVal & 0x80)
    {
      while(vCnt--)
      {
	*dst++ = *src++;
      }
    }
    else
    {
      vVal = *src++;
      while(vCnt--)
      {
	*dst++ = vVal;
      }
    }
    vVal = *src++;
  } 
}

/*!
* \return	Woolz compound array object.
* \ingroup	WlzExtFF
* \brief	Given an object in which distinct voxel values represent
*		different domains and each of these has a corresponding
*		material, this function creates a new compound array
*		object. Each of the materials has it's own domain and
*		the material properties are encoded in a simple ascii
*		property list of the domain.
*		The materials are known to have indices which increase
*		from 0.
* \param	gObj			Given index object.
* \param	head			Amira header data structure with the
*					list of materials.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
static WlzCompoundArray *WlzEffAmSplitLabelObj(WlzObject *gObj,
					WlzEffAmHead *head,
					WlzErrorNum *dstErr)
{
  int		idx,
  		empty = 0;
  WlzPixelV	thrV;
  WlzObject	*tObj0 = NULL,
  		*tObj1 = NULL,
		*tObj2,
		*tObj3;
  WlzEffAmMaterial *mat;
  WlzPropertyList *pList = NULL;
  WlzCompoundArray *aObj = NULL;
  WlzDomain	nullDom;
  WlzValues	dValues,
  		nullVal;
  WlzErrorNum errNum = WLZ_ERR_NONE;
  const int	allEmptyAfterFirst = 0;

  idx = 0;
  nullDom.core = NULL;
  nullVal.core = NULL;
  dValues.core = NULL;
  mat = head->materials;
  thrV.type = WLZ_GREY_INT;
  /* Create compound array object. */
  aObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, head->matCount, NULL,
  			      WLZ_3D_DOMAINOBJ, &errNum);
  /* For each material create a new domain corresponding to the index. */
  while((errNum == WLZ_ERR_NONE) && (idx < head->matCount))
  {
    thrV.v.inv = mat->id + 1;
    pList = WlzEffAmMakeMaterialPropList(mat, &errNum);
    if(empty && allEmptyAfterFirst)
    {
      tObj3 = WlzMakeMain(WLZ_EMPTY_OBJ, nullDom, nullVal, pList,
      			  NULL, &errNum);
      
    }
    else
    {
      if(errNum == WLZ_ERR_NONE)
      {
	if(idx == 0)
	{
	  tObj1 = WlzAssignObject(
		  WlzThreshold(gObj, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tObj2 = WlzDiffDomain(gObj, tObj1, &errNum);
	  }
	}
	else
	{
	  tObj1 = WlzAssignObject(
		  WlzThreshold(tObj0, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tObj2 = WlzDiffDomain(tObj0, tObj1, &errNum);
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if(WlzIsEmpty(tObj1, &errNum) || WlzIsEmpty(tObj2, &errNum))
	{
	  empty = 1;
	  tObj3 = WlzMakeMain(WLZ_EMPTY_OBJ, nullDom, nullVal, pList,
	                      NULL, &errNum);
	}
	else
	{
	  tObj3 = WlzMakeMain(WLZ_3D_DOMAINOBJ, tObj2->domain,
			      dValues, pList, NULL, &errNum);
	}
      }
      (void )WlzFreeObj(tObj2); tObj2 = NULL;
      (void )WlzFreeObj(tObj0);
      tObj0 = tObj1;
      tObj1 = NULL;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      aObj->o[idx] = WlzAssignObject(tObj3, NULL);
    }
    ++idx;
    mat = mat->next;
  }
  (void )WlzFreeObj(tObj0);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(aObj);
}

/*!
* \return	New property list.
* \ingroup	WlzExtFF
* \brief	Creates a new property list from an Amira material
*		specification.
*		The WlzNameProperty is used for the material name and
*		WlzGreyProperty is used for the color.
* \param	mat			Material specification.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
static WlzPropertyList *WlzEffAmMakeMaterialPropList(WlzEffAmMaterial *mat,
					WlzErrorNum *dstErr)
{
  int		idx;
  unsigned int	tI0;
  WlzProperty	prop;
  WlzPixelV	pix;
  WlzPropertyList *pList = NULL;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  pix.v.rgbv = 0;
  for(idx = 0; idx < 3; ++idx)
  {
    if(mat->color[idx] < (0.0 + DBL_EPSILON))
    {
      tI0 = 0;
    }
    else if(mat->color[idx] > (1.0 - DBL_EPSILON))
    {
      tI0 = 255;
    }
    else
    {
      tI0 = (int )floor((mat->color[idx] * 255.0) + 0.5); 
    }
    pix.v.rgbv |= tI0 << idx;
  }
  if((pList = WlzMakePropertyList(NULL)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    prop.name = WlzMakeNameProperty(mat->name, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
			     WlzFreePropertyListEntry) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    pix.type = WLZ_GREY_RGBA;
    prop.greyV = WlzMakeGreyProperty("Colour", pix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
			     WlzFreePropertyListEntry) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum != WLZ_ERR_NONE) && pList)
  {
    (void )WlzFreePropertyList(pList);
    pList = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pList);
}
