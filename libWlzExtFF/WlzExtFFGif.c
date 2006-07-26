#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFGif_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzExtFFGif.c
* \author       Bill Hill
* \date         June 2006
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* The GIF decoding software in this file is derived from the pbmplus
* utilities GIF decoding, which in turn is based on GIFDECOD by
* David Koblas. The most significant differences between the pbmplus
* software and the software below are the inclusion of Woolz and
* a significant rewrite of the code to make it re-enterant. See the
* notice below:
* \verbatim
  +-------------------------------------------------------------------+
  | Copyright 1990, David Koblas.                                     |
  |   Permission to use, copy, modify, and distribute this software   |
  |   and its documentation for any purpose and without fee is hereby |
  |   granted, provided that the above copyright notice appear in all |
  |   copies and that both that copyright notice and this permission  |
  |   notice appear in supporting documentation.  This software is    |
  |   provided "as is" without express or implied warranty.           |
  +-------------------------------------------------------------------+
* \endverbatim
* \par
* Copyright (C) 2006 Medical research Council, UK.
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
* \brief	Functions for reading Woolz objects from Graphics
*		Interchange Format (GIF) images.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <string.h>
#include <Wlz.h>


#define	WLZEFF_GIF_MAXCMAPSZ	256
#define	WLZEFF_GIF_BUFSZ	256
#define WLZEFF_GIF_CM_RED	0
#define WLZEFF_GIF_CM_GREEN	1
#define WLZEFF_GIF_CM_BLUE	2
#define	WLZEFF_GIF_MAXLWZBITS	12
#define WLZEFF_GIF_INTERLACE	0x40
#define WLZEFF_GIF_LOCALCMAP	0x80

#define BitSet(byte, bit)	(((byte) & (bit)) == (bit))
#define LM_to_uint(a,b)		(((b)<<8)|(a))

typedef struct _WlzEffGifScreen
{
  unsigned int	Width;
  unsigned int	Height;
  unsigned int	BitPixel;
  unsigned int	ColorResolution;
  unsigned int	Background;
  unsigned int	AspectRatio;
  unsigned char	ColorMap[3][WLZEFF_GIF_MAXCMAPSZ];
} WlzEffGifScreen;

typedef struct _WlzEffGif89
{
  int		transparent;
  int		delayTime;
  int		inputFlag;
  int		disposal;
} WlzEffGif89;

typedef struct _WlzEffGifLWZWSp
{
  int		fresh;
  int		codeSz;
  int		setCodeSz;
  int		maxCode;
  int		maxCodeSz;
  int		fstCode;
  int		oldCode;
  int		clrCode;
  int		endCode;
  int		*sp;
  int		table[2][(1<< WLZEFF_GIF_MAXLWZBITS)];
  int		stack[(1<<(WLZEFF_GIF_MAXLWZBITS))*2];
} WlzEffGifLWZWSp;


static WlzErrorNum 		WlzEffGifColorMapRead(
				  FILE *fP,
				  int nC,
				  WlzUByte map[3][WLZEFF_GIF_MAXCMAPSZ]);
static WlzErrorNum 		WlzEffGifExtensionRead(
				  FILE *fP,
				  int label,
				  WlzEffGif89 *gif89,
				  char *comment,
				  int *zdb);
static int 			WlzEffGifDataBlkRead(
				  FILE *fP,
				  WlzUByte *buf,
				  int *zdb);
static int			WlzEffGifCodeRead(
				  FILE *fP,
				  int code_size,
				  int flag,
				  int *zdb);
static WlzEffGifLWZWSp		*WlzEffGifMakeLZWWSp(void);
static int 			WlzEffGifLWZReadByte(
				  FILE *fP,
				  int flag,
                                  int inCodeSz,
				  int *zdb,
				  WlzEffGifLWZWSp *wSp);
static WlzErrorNum		WlzEffGifImageRead(
				  FILE *fP,
                      		  WlzIVertex2 imgSz,
				  unsigned **img,
                      		  WlzUByte cmap[3][WLZEFF_GIF_MAXCMAPSZ],
				  int interlace,
		      		  int *zdb);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the file pointer using the '.gif'
*               file format.
*		This function currently has the following limitations:
*		<ul>
*		  <li>the aspect ratio of the pixels. is ignored.
*                 <li>only the first frame of an animated GIF is read.
*                 <li>transparency is ignored.
*		</ul>
* \param	fP			File pointer.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzEffReadObjGif(FILE *fP, WlzErrorNum *dstErr)
{
  WlzUByte	c;
  int		imageNumber = 1;
  int		useGlobalColormap;
  int		bitPixel;
  int		zdb = 0;
  int		imgCnt = 0;
  unsigned	**img = NULL;
  WlzIVertex2	imgSz,
  		imgOrg;
  char		version[4];
  char		comment[WLZEFF_GIF_BUFSZ];
  WlzUByte	buf[16];
  WlzUByte	localColorMap[3][WLZEFF_GIF_MAXCMAPSZ];
  WlzEffGif89	gif89 = { -1, -1, -1, 0 };
  WlzEffGifScreen gifScreen;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  imgOrg.vtX = imgOrg.vtY = 0.0;
  /* Check for valid GIF identifier. */
  if((fread(buf, 6, 1, fP) != 1) || (strncmp((char *)buf,"GIF",3) != 0))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    (void )strncpy(version, (char *)buf + 3, 3);
    version[3] = '\0';
    if((strcmp(version, "87a") != 0) && (strcmp(version, "89a") != 0))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Read screen descriptor. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fread(buf, 7, 1, fP) != 1)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gifScreen.Width           = LM_to_uint(buf[0],buf[1]);
    gifScreen.Height          = LM_to_uint(buf[2],buf[3]);
    gifScreen.BitPixel        = 2<<(buf[4]&0x07);
    gifScreen.ColorResolution = (((buf[4]&0x70)>>3)+1);
    gifScreen.Background      = buf[5];
    gifScreen.AspectRatio     = buf[6];
    /* Check for global colourmap */
    if(BitSet(buf[4], WLZEFF_GIF_LOCALCMAP))
    {
      /* Reading global colourmap. */
      errNum = WlzEffGifColorMapRead(fP,gifScreen.BitPixel,gifScreen.ColorMap);
    }
  }
  while((errNum == WLZ_ERR_NONE) && (imgCnt < imageNumber))
  {
    if(fread(&c, 1, 1, fP) != 1)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      if(c == ';') /* GIF terminator */
      {
	if(imgCnt < imageNumber)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  continue;
        }
      }
      if(c == '!') /* GIF extension */
      {
	if(fread(&c, 1, 1, fP) != 1)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  continue;
	}
	errNum = WlzEffGifExtensionRead(fP, c, &gif89, comment, &zdb);
	continue;
      }
      if(c != ',') /* Not a valid start character */
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	continue;
      }
      ++imgCnt;
      /* Read left/top/width/height. */
      if(fread(buf, 9, 1, fP) != 1)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      useGlobalColormap = ! BitSet(buf[8], WLZEFF_GIF_LOCALCMAP);
      bitPixel = 1<<((buf[8]&0x07)+1);
      imgSz.vtY = LM_to_uint(buf[6],buf[7]);
      imgSz.vtX = LM_to_uint(buf[4],buf[5]);
      if(AlcInt2Malloc((int ***)&img, imgSz.vtY, imgSz.vtX) != ALC_ER_NONE)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	if(!useGlobalColormap)
	{
	  if(WlzEffGifColorMapRead(fP, bitPixel,
	                           localColorMap) != WLZ_ERR_NONE)
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
	  }
	  errNum = WlzEffGifImageRead(fP, imgSz, img, localColorMap,
			     BitSet(buf[8], WLZEFF_GIF_INTERLACE), &zdb);
	}
	else
	{
	  errNum = WlzEffGifImageRead(fP, imgSz, img, gifScreen.ColorMap,
			     BitSet(buf[8], WLZEFF_GIF_INTERLACE), &zdb);
	}
      }
    }
  }
  obj = WlzFromArray2D((void **)img, imgSz, imgOrg,
  		       WLZ_GREY_RGBA, WLZ_GREY_RGBA, 0.0, 1.0, 0, 0, &errNum);
  (void )Alc2Free((void **)img);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads a GIF colour index map.
* \param	fP			File pointer.
* \param	nC			Number of colour entries.
* \param	buf			The map to be filled in.
*/
static WlzErrorNum WlzEffGifColorMapRead(FILE *fP, int nC,
				       WlzUByte map[3][WLZEFF_GIF_MAXCMAPSZ])
{
  int		idx;
  WlzUByte	rgb[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  for(idx = 0; idx < nC; ++idx)
  {
    if(fread(rgb, sizeof(rgb), 1, fP) != 1)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
    }
    map[WLZEFF_GIF_CM_RED][idx] = rgb[0] ;
    map[WLZEFF_GIF_CM_GREEN][idx] = rgb[1] ;
    map[WLZEFF_GIF_CM_BLUE][idx] = rgb[2] ;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads FIF89 extentions.
* \param	fP			File pointer.
* \param	label			Extention label.
* \param	gif89			GIF89 data structure to fill.
* \param	comment			Buffer for comments (256 bytes).
* \param	zdb			Zero data block flag.
*/
static WlzErrorNum WlzEffGifExtensionRead(FILE *fP, int label,
					  WlzEffGif89 *gif89, char *comment,
					  int *zdb)
{
  int		cInc,
		cFirst = 1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char	buf[WLZEFF_GIF_BUFSZ];

  switch(label)
  {
    case 0x01: /* Plain Text Extension */
      break;
    case 0xff: /* Application Extension */
      break;
    case 0xfe: /* Comment Extension */
      *comment = '\0';
      while((cInc = WlzEffGifDataBlkRead(fP, (WlzUByte *)buf, zdb)) > 0)
      {
        if(cFirst)
	{
	  cFirst = 0;
	  (void )strcpy(comment, buf);
	}
      }
      if(cInc < 0)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      break;
    case 0xf9: /* Graphic Control Extension */
      (void) WlzEffGifDataBlkRead(fP, (unsigned char*) buf, zdb);
      gif89->disposal    = (buf[0] >> 2) & 0x7;
      gif89->inputFlag   = (buf[0] >> 1) & 0x1;
      gif89->delayTime   = LM_to_uint(buf[1],buf[2]);
      if((buf[0] & 0x1) != 0)
      {
	gif89->transparent = buf[3];
      }
      while(WlzEffGifDataBlkRead(fP, (unsigned char*) buf, zdb) != 0)
	;
      break;
    default:
      /* Try to ignore any unknown extensions! */
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads the image data.
* \param	fP			File pointer.
* \param	imgSz			Image size.
* \param	img			Allocated image data.
* \param	cmap			Colour map.
* \param	interlace		Interlace flag.
* \param	zdb			Zero data block flag.
*/
static WlzErrorNum WlzEffGifImageRead(FILE *fP,
		      WlzIVertex2 imgSz, unsigned **img,
                      WlzUByte cmap[3][WLZEFF_GIF_MAXCMAPSZ],
                      int interlace, int *zdb)
{
  WlzUByte	alpha,
  		codeSz,
		red,
		green,
		blue;	
  int		v,
     		pass = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIVertex2	imgPos;
  WlzEffGifLWZWSp *wSp = NULL;

  /* Create workspace and initialize the Compression routines. */
  if((wSp = WlzEffGifMakeLZWWSp()) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    goto done;
  }
  if(fread(&codeSz, 1, 1, fP) != 1)
  {
    /* EOF / read error on image data. */
    errNum = WLZ_ERR_READ_INCOMPLETE;
    goto done;
  }
  if(WlzEffGifLWZReadByte(fP, 1, codeSz, zdb, wSp) < 0)
  {
    /* Error reading image. */
    errNum = WLZ_ERR_READ_INCOMPLETE;
    goto done;
  }
  pass = 0;
  imgPos.vtX = 0;
  imgPos.vtY = 0;
  while(((v = WlzEffGifLWZReadByte(fP, 0, codeSz, zdb, wSp)) >= 0) &&
         (imgPos.vtY < imgSz.vtY))
  {
    red = cmap[WLZEFF_GIF_CM_RED][v];
    green = cmap[WLZEFF_GIF_CM_GREEN][v];
    blue = cmap[WLZEFF_GIF_CM_BLUE][v];
    alpha = 0;
    WLZ_RGBA_RGBA_SET(img[imgPos.vtY][imgPos.vtX], red, green, blue, alpha);
    ++(imgPos.vtX);
    if(imgPos.vtX == imgSz.vtX)
    {
      imgPos.vtX = 0;
      if(interlace)
      {
	switch(pass)
	{
	  case 0: /* FALLTHROUGH */
	  case 1:
	    imgPos.vtY += 8;
	    break;
	  case 2:
	    imgPos.vtY += 4;
	    break;
	  case 3:
	    imgPos.vtY += 2;
	    break;
	}
	if(imgPos.vtY >= imgSz.vtY)
	{
	  switch(++pass)
	  {
	    case 1:
	      imgPos.vtY = 4;
	      break;
	    case 2:
	      imgPos.vtY = 2;
	      break;
	    case 3:
	      imgPos.vtY = 1;
	      break;
	    default:
	      goto done;
	  }
	}
      }
      else
      {
	++(imgPos.vtY);
      }
    }
  }
done:
  AlcFree(wSp);
  return(errNum);
}

/*!
* \return	Number of bytes or -ve on error.
* \ingroup	WlzExtFF
* \brief	Reads a GIF data block.
* \param	fP			File pointer.
* \param	buf			Buffer for data (256 bytes).
* \param	zdb			Zero data block flag.
*/
static int 	WlzEffGifDataBlkRead(FILE *fP, WlzUByte *buf, int *zdb)
{
  int		cnt = 0;

  if(fread(buf, 1, 1, fP) != 1)
  {
    /* Error in getting DataBlock size. */
    cnt = -1;
  }
  else
  {
    cnt = *buf;
    *zdb = (cnt == 0);
    if(cnt && (fread(buf, cnt, 1, fP) != 1))
    {
      /* Error in reading DataBlock. */
      cnt = -1;
    }
  }
  return(cnt);
}

/*!
* \return	New GIF LZW workspace.
* \ingroup	WlzExtFF
* \brief	Creates a new initialized GIF LZW workspace for use in
*		decoding the LZW encoded data.
*/
static WlzEffGifLWZWSp *WlzEffGifMakeLZWWSp(void)
{
  WlzEffGifLWZWSp *wSp;

  wSp = (WlzEffGifLWZWSp *)AlcCalloc(1, sizeof(WlzEffGifLWZWSp));
  return(wSp);
}

/*!
* \return	Byte from file or -ve on error.
* \ingroup	WlzExtFF
* \brief	Reads a byte from the LZW encoded data.
* \param	fP			File pointer.
* \param	flag
* \param	inCodeSz		Code size read from file.
* \param	zdb			Zero data block flag.
* \param	wSp			GIF LZW workspace.
*/
static int 	WlzEffGifLWZReadByte(FILE *fP, int flag, int inCodeSz,
				int *zdb, WlzEffGifLWZWSp *wSp)
{
  int		idx,
		cnt,
  		code,
		inCode;
  WlzUByte	buf[260];


  if(flag)
  {
    wSp->setCodeSz = inCodeSz;
    wSp->codeSz = wSp->setCodeSz + 1;
    wSp->clrCode = 1 << wSp->setCodeSz ;
    wSp->endCode = wSp->clrCode + 1;
    wSp->maxCodeSz = 2 * wSp->clrCode;
    wSp->maxCode = wSp->clrCode + 2;
    WlzEffGifCodeRead(fP, 0, 1, zdb);
    wSp->fresh = 1;
    for(idx = 0; idx < wSp->clrCode; ++idx)
    {
      wSp->table[0][idx] = 0;
      wSp->table[1][idx] = idx;
    }
    while(idx < (1 << WLZEFF_GIF_MAXLWZBITS))
    {
      wSp->table[0][idx++] = wSp->table[1][0] = 0;
    }
    wSp->sp = wSp->stack;
    code = 0;
    goto done;
  }
  else if(wSp->fresh)
  {
    wSp->fresh = 0;
    do
    {
      wSp->fstCode = wSp->oldCode = WlzEffGifCodeRead(fP, wSp->codeSz, 0, zdb);
    } while (wSp->fstCode == wSp->clrCode);
    code = wSp->fstCode;
    goto done;
  }
  if(wSp->sp > wSp->stack)
  {
    code = *--(wSp->sp);
    goto done;
  }
  while((code = WlzEffGifCodeRead(fP, wSp->codeSz, 0, zdb)) >= 0)
  {
    if(code == wSp->clrCode)
    {
      for(idx = 0; idx < wSp->clrCode; ++idx)
      {
	wSp->table[0][idx] = 0;
	wSp->table[1][idx] = idx;
      }
      for(; idx < (1 << WLZEFF_GIF_MAXLWZBITS); ++idx)
      {
	wSp->table[0][idx] = wSp->table[1][idx] = 0;
      }
      wSp->codeSz = wSp->setCodeSz + 1;
      wSp->maxCodeSz = 2 * wSp->clrCode;
      wSp->maxCode = wSp->clrCode + 2;
      wSp->sp = wSp->stack;
      wSp->fstCode = wSp->oldCode = WlzEffGifCodeRead(fP, wSp->codeSz, 0, zdb);
      code = wSp->fstCode;
      goto done;
    }
    else if(code == wSp->endCode)
    {
      if(*zdb)
      {
	code = -2;
	goto done;
      }
      while((cnt = WlzEffGifDataBlkRead(fP, buf, zdb)) > 0);
      /* Error: Missing EOD in data stream (cnt != 0) is a common occurence. */
      code = -2;
      goto done;
    }
    inCode = code;
    if(code >= wSp->maxCode)
    {
      *wSp->sp++ = wSp->fstCode;
      code = wSp->oldCode;
    }
    while(code >= wSp->clrCode)
    {
      *wSp->sp++ = wSp->table[1][code];
      if(code == wSp->table[0][code])
      {
	/* Error: Circular table entry! */
	code = -3;
	goto done;
      }
      code = wSp->table[0][code];
    }
    *wSp->sp++ = wSp->fstCode = wSp->table[1][code];
    if((code = wSp->maxCode) <(1 << WLZEFF_GIF_MAXLWZBITS))
    {
      wSp->table[0][code] = wSp->oldCode;
      wSp->table[1][code] = wSp->fstCode;
      ++wSp->maxCode;
      if((wSp->maxCode >= wSp->maxCodeSz) &&
	 (wSp->maxCodeSz < (1 << WLZEFF_GIF_MAXLWZBITS)))
      {
	wSp->maxCodeSz *= 2;
	++wSp->codeSz;
      }
    }
    wSp->oldCode = inCode;
    if(wSp->sp > wSp->stack)
    {
      code = *--(wSp->sp);
      goto done;
    }
  }
done:
  return(code);
}

/*!
* \return	New code.
* \ingroup	WlzExtFF
* \brief	reads code from file.
* \param	fP			File pointer.
* \param	codeSz			Code size read from file.
* \param	flag
* \param	zdb			Zero data block flag.
*/
static int WlzEffGifCodeRead(FILE *fP, int codeSz, int flag, int *zdb)
{
  static int		curbit, lastbit, done, last_byte;
  int			i, j, ret;
  unsigned char		cnt;
  static unsigned char	buf[280];

  if(flag)
  {
    curbit = 0;
    lastbit = 0;
    done = 0;
    ret = 0;
  }
  else
  {
    if((curbit + codeSz) >= lastbit)
    {
      if(done)
      {
	if(curbit >= lastbit)
	{
	  /* Ran off the end of the bits. */
	  ret = -1;
	  goto done;
	}
      }
      buf[0] = buf[last_byte-2];
      buf[1] = buf[last_byte-1];
      done = (cnt = WlzEffGifDataBlkRead(fP, &buf[2], zdb)) == 0;
      last_byte = 2 + cnt;
      curbit = (curbit - lastbit) + 16;
      lastbit = (2+cnt)*8 ;
    }
    ret = 0;
    for (i = curbit, j = 0; j < codeSz; ++i, ++j)
    {
      ret |= ((buf[ i / 8 ] & (1 << (i % 8))) != 0) << j;
    }
    curbit += codeSz;
  }
done:
  return(ret);
}
