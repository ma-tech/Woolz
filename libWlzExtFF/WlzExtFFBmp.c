#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFBmp.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2005 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for reading and writting Woolz objects to
*		and from the portable anymap '.pnm' data format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/
#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static int			WlzExtFFBmpEncode8(
				  unsigned char **rleData,
				  unsigned char *imgData,
				  WlzIVertex2 imgSz,
				  WlzErrorNum *dstErr);
static void			WlzEffBmpSwap2(
				  void *data);
static void			WlzEffBmpSwap4(
				  void *data);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given file(s) using the '.pnm'
* 		file format.
* \param	gvnFileName		Given file name.
* \param	dstErr			Dst error pointer, may be NULL.
*/
WlzObject	*WlzEffReadObjBmp(const char *gvnFileName, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;

  obj = WlzEffReadObjStack(gvnFileName, WLZEFF_FORMAT_BMP, &errNum);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file(s)
*		using the '.pnm' file format.
* \param	gvnFileName		Given file name.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjBmp(const char *gvnFileName, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzEffWriteObjStack(gvnFileName, WLZEFF_FORMAT_BMP, obj);
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given 2D Woolz object to the given file
*		using the '.pnm' file format.
* \param	fNameStr		Given file name.
* \param	obj			Given woolz object.
* \param	imgSz			Required image size.
* \param	imgOrg			Required image origin.
* \param	data			Buffer of imgSz bytes.
* \param	bgd			Background value.
*/
WlzErrorNum 	WlzEffWriteObjBmp2D(const char *fNameStr, WlzObject *obj,
				    WlzIVertex2 imgSz, WlzIVertex2 imgOrg,
				    unsigned char *data,
				    unsigned char bgd)
{
  int		tI0,
		idN,
		rleSz;
  unsigned char *rleData;
  FILE		*fP = NULL;
  WlzObject	*cutObj = NULL;
  WlzIBox2	cutBox;
  WlzEffBmpRGBQuad mapC;
  WlzEffBmpFileHead fHead;
  WlzEffBmpInfoHead iHead;
  const unsigned char padBuf[] = {0x00, 0x00, 0x00, 0x00};
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fP = fopen(fNameStr, "w")) == NULL)
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else
  {

  #ifdef _WIN32
	if (fP != NULL){
		if(_setmode(_fileno(fP), 0x8000) == -1)
		{
			errNum = WLZ_ERR_READ_EOF;
		}
	}
	#endif

	if((obj == NULL) || (obj->values.core == NULL))
    {
      (void )memset(data, (unsigned int )bgd, imgSz.vtX * imgSz.vtY);
    }
    else
    {
      cutBox.xMin = imgOrg.vtX;
      cutBox.yMin = imgOrg.vtY;
      cutBox.xMax = imgOrg.vtX + imgSz.vtX - 1;
      cutBox.yMax = imgOrg.vtY + imgSz.vtY - 1;
      cutObj = WlzCutObjToValBox2D(obj, cutBox, WLZ_GREY_UBYTE, data,
      			        0, 0.0, 0.0, &errNum);
      if(cutObj)
      {
	if(cutObj->type == WLZ_2D_DOMAINOBJ)
	{
	  /* (void )WlzPopFreePtr(cutObj->values.r->freeptr, NULL, NULL); */
	  AlcFree(cutObj->values.r->freeptr);
	  cutObj->values.r->freeptr = NULL;
	  cutObj->values.r->values.inp = NULL;
	}
	WlzFreeObj(cutObj);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Run length encode the data. */
    rleSz = WlzExtFFBmpEncode8(&rleData, data, imgSz, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute and write the headers. */
    iHead.biSize = WLZEFF_BMP_WIN_NEW;
    iHead.biWidth = imgSz.vtX;
    iHead.biHeight = imgSz.vtY;
    iHead.biPlanes = 1;
    iHead.biBitCount = 8;
    iHead.biCompression = WLZEFF_BMP_CMP_RLE8;
    iHead.biSizeImage = rleSz;
    iHead.biXPelsPerMeter = 75 * 39; 		    /* 75dpi * 29" per metre */
    iHead.biYPelsPerMeter = 75 * 39; 		    /* 75dpi * 29" per metre */
    iHead.biClrUsed = 256;
    iHead.biClrImportant = 256;
    fHead.bfType = (WLZEFF_BMP_MAGIC_1) | ((WLZEFF_BMP_MAGIC_0) << 8);
    tI0 = sizeof(WLZEFF_BMP_UINT) + sizeof(WLZEFF_BMP_DWORD) +
	  sizeof(WLZEFF_BMP_UINT) + sizeof(WLZEFF_BMP_UINT) +
	  sizeof(WLZEFF_BMP_DWORD) + WLZEFF_BMP_WIN_NEW +
	  (iHead.biClrUsed * 4);
    fHead.bfSize = tI0 + rleSz;
    fHead.bfOffBits = tI0;
	
#ifndef _WIN32
    WlzEffBmpSwap4((void *)&(fHead.bfSize));
    WlzEffBmpSwap4((void *)&(fHead.bfOffBits));
    WlzEffBmpSwap4((void *)&(iHead.biSize));
    WlzEffBmpSwap4((void *)&(iHead.biWidth));
    WlzEffBmpSwap4((void *)&(iHead.biHeight));
    WlzEffBmpSwap2((void *)&(iHead.biPlanes));
    WlzEffBmpSwap2((void *)&(iHead.biBitCount));
    WlzEffBmpSwap4((void *)&(iHead.biCompression));
    WlzEffBmpSwap4((void *)&(iHead.biSizeImage));
    WlzEffBmpSwap4((void *)&(iHead.biXPelsPerMeter));
    WlzEffBmpSwap4((void *)&(iHead.biYPelsPerMeter));
    WlzEffBmpSwap4((void *)&(iHead.biClrUsed));
    WlzEffBmpSwap4((void *)&(iHead.biClrImportant));
#else
	WlzEffBmpSwap2((void *)&(fHead.bfType));
	WlzEffBmpSwap2((void *)&(fHead.bfReserved1));
	WlzEffBmpSwap2((void *)&(fHead.bfReserved2));
#endif 
	
    if((fwrite(&(fHead.bfType), sizeof(WLZEFF_BMP_UINT),
    	       1, fP) != 1) ||
       (fwrite(&(fHead.bfSize), sizeof(WLZEFF_BMP_DWORD),
       	       1, fP) != 1) ||
       (fwrite(&(fHead.bfReserved1), sizeof(WLZEFF_BMP_UINT),
       	       1, fP) != 1) ||
       (fwrite(&(fHead.bfReserved2), sizeof(WLZEFF_BMP_UINT),
       	       1, fP) != 1) ||
       (fwrite(&(fHead.bfOffBits), sizeof(WLZEFF_BMP_DWORD),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biSize), sizeof(WLZEFF_BMP_DWORD),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biWidth), sizeof(WLZEFF_BMP_LONG),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biHeight), sizeof(WLZEFF_BMP_LONG),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biPlanes), sizeof(WLZEFF_BMP_WORD),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biBitCount), sizeof(WLZEFF_BMP_WORD),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biCompression), sizeof(WLZEFF_BMP_DWORD),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biSizeImage), sizeof(WLZEFF_BMP_DWORD),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biXPelsPerMeter), sizeof(WLZEFF_BMP_LONG),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biYPelsPerMeter), sizeof(WLZEFF_BMP_LONG),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biClrUsed), sizeof(WLZEFF_BMP_DWORD),
       	       1, fP) != 1) ||
       (fwrite(&(iHead.biClrImportant), sizeof(WLZEFF_BMP_DWORD),
       	       1, fP) != 1))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Write the colormap. */
    idN = 0;
    mapC.rgbReserved = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < 256))
    {
      mapC.rgbBlue = mapC.rgbGreen = mapC.rgbRed = idN;
      if((fwrite(&(mapC.rgbBlue), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1) ||
	 (fwrite(&(mapC.rgbGreen), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1) ||
	 (fwrite(&(mapC.rgbRed), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1) ||
	 (fwrite(&(mapC.rgbReserved), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1))
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      ++idN;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Write the runlength encoded data. */
    if(fwrite(rleData, sizeof(unsigned char), rleSz, fP) != rleSz)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(rleData)
  {
    AlcFree(rleData);
  }
  if(fP)
  {
    fclose(fP);
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Reads the  data from a pgm file into the data array
*		given.
* \param	fP			Given file stream.
* \param	gvnImgSz		Dst ptr for image size which is
*					assumed valid if height and width are
*					non zero.
* \param	 data			Ptr to 2D array for data.
*/
WlzErrorNum 	WlzEffReadObjBmpData2D(FILE *fP, WlzIVertex2 *gvnImgSz,
				       unsigned char ***data)
{
  int		tI0,
  		tI1,
		idN,
  		idX,
		idY,
		endBM,
		imgSz,
		padByte,
		rlCnt,
  		swapFlg = 0;
  unsigned char	*dataPtr;
  unsigned char	**dataBuf;
  WlzEffBmpRGBQuad *tMap,
  		*rgbMap = NULL;
  WLZEFF_BMP_BYTE *imgPtr,
  		*imgBuf = NULL;
  WlzEffBmpFileHead fHead;
  WlzEffBmpInfoHead iHead;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Read the bitmap file header. */
  if((fread(&(fHead.bfType), sizeof(WLZEFF_BMP_UINT), 1, fP) != 1) ||
     (fread(&(fHead.bfSize), sizeof(WLZEFF_BMP_DWORD), 1, fP) != 1) ||
     (fread(&(fHead.bfReserved1), sizeof(WLZEFF_BMP_UINT), 1, fP) != 1) ||
     (fread(&(fHead.bfReserved2), sizeof(WLZEFF_BMP_UINT), 1, fP) != 1) ||
     (fread(&(fHead.bfOffBits), sizeof(WLZEFF_BMP_DWORD), 1, fP) != 1))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    /* Check magic number and determine byte order. */
    if(((fHead.bfType & 0xff) == WLZEFF_BMP_MAGIC_1) &&
       ((fHead.bfType >> 8) == WLZEFF_BMP_MAGIC_0))
    {
      swapFlg = 1;
    }
    else if(((fHead.bfType & 0xff) != WLZEFF_BMP_MAGIC_0) ||
            ((fHead.bfType >> 8) != WLZEFF_BMP_MAGIC_1))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(swapFlg)
    {
      WlzEffBmpSwap2((void *)&(fHead.bfType));
      WlzEffBmpSwap4((void *)&(fHead.bfSize));
      WlzEffBmpSwap4((void *)&(fHead.bfOffBits));
    }
    /* Read the bitmap info header size to determine version. */
    if(fread(&(iHead.biSize), sizeof(WLZEFF_BMP_DWORD), 1, fP) != 1)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      if(swapFlg)
      {
        WlzEffBmpSwap4((void *)&(iHead.biSize));
      }
      /* Only bother to read versions created with MS Windows 3.0 or above! */
      if(iHead.biSize != WLZEFF_BMP_WIN_NEW)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Read the rest of bitmap info header. */
    if((fread(&(iHead.biWidth), sizeof(WLZEFF_BMP_LONG), 1,
    	      fP) != 1) ||
       (fread(&(iHead.biHeight), sizeof(WLZEFF_BMP_LONG), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biPlanes), sizeof(WLZEFF_BMP_WORD), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biBitCount), sizeof(WLZEFF_BMP_WORD), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biCompression), sizeof(WLZEFF_BMP_DWORD), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biSizeImage), sizeof(WLZEFF_BMP_DWORD), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biXPelsPerMeter), sizeof(WLZEFF_BMP_LONG), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biYPelsPerMeter), sizeof(WLZEFF_BMP_LONG), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biClrUsed), sizeof(WLZEFF_BMP_DWORD), 1,
       	      fP) != 1) ||
       (fread(&(iHead.biClrImportant), sizeof(WLZEFF_BMP_DWORD), 1,
       	      fP) != 1))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      if(swapFlg)
      {
	WlzEffBmpSwap4((void *)&(iHead.biWidth));
	WlzEffBmpSwap4((void *)&(iHead.biHeight));
	WlzEffBmpSwap2((void *)&(iHead.biPlanes));
	WlzEffBmpSwap2((void *)&(iHead.biBitCount));
	WlzEffBmpSwap4((void *)&(iHead.biCompression));
	WlzEffBmpSwap4((void *)&(iHead.biSizeImage));
	WlzEffBmpSwap4((void *)&(iHead.biXPelsPerMeter));
	WlzEffBmpSwap4((void *)&(iHead.biYPelsPerMeter));
	WlzEffBmpSwap4((void *)&(iHead.biClrUsed));
	WlzEffBmpSwap4((void *)&(iHead.biClrImportant));
      }
      /* Only accept 8 or 24 bit images. */
      if((iHead.biBitCount != 8) && (iHead.biBitCount != 24))
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Read colour map if it's an 8 bit image. */
    if(iHead.biBitCount == 8)
    {
      if((rgbMap = (WlzEffBmpRGBQuad *)
      		   AlcCalloc(sizeof(WlzEffBmpRGBQuad), 256)) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	if(iHead.biClrUsed == 0)
	{
	  iHead.biClrUsed = 1 << iHead.biBitCount;
	}
	idN = 0;
        while((errNum == WLZ_ERR_NONE) && (idN < iHead.biClrUsed))
	{
	  tMap = rgbMap + idN;
	  if((fread(&(tMap->rgbBlue), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1) ||
	     (fread(&(tMap->rgbGreen), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1) ||
	     (fread(&(tMap->rgbRed), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1) ||
	     (fread(&(tMap->rgbReserved), sizeof(WLZEFF_BMP_BYTE), 1, fP) != 1))
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  ++idN;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Read image data. */
    /* Skip to the start of the image data. */
    tI0 = fHead.bfOffBits - ftell(fP);
    while(tI0-- > 0)
    {
      (void )fgetc(fP);
    }
    /* Allocate image buffer. */
    imgSz = (((iHead.biWidth * iHead.biBitCount) + 31) /32) * 4 *
    	    iHead.biHeight;
    if((imgBuf = (WLZEFF_BMP_BYTE *)
    		 AlcMalloc(imgSz * sizeof(WLZEFF_BMP_BYTE))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      switch(iHead.biCompression)
      {
	case WLZEFF_BMP_CMP_RGB:			/* Uncompressed data */
	  if(fread(imgBuf,  1, imgSz, fP) != imgSz)
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  break;
	case WLZEFF_BMP_CMP_RLE8:			  /* Compressed data */
	  idX = 0;
	  idY = 0;
	  endBM = 0;
	  imgPtr = imgBuf;
	  while((idY < iHead.biHeight) && (endBM == 0) &&
		((rlCnt = fgetc(fP)) != EOF))
	  {
	    if(rlCnt == WLZEFF_BMP_CMP_ESC)		      /* Escape mode */
	    {
	      switch(rlCnt = fgetc(fP))
	      {
		case WLZEFF_BMP_CMP_EOL:		     /* End of line  */
		  idX = 0;
		  if(++idY >= iHead.biHeight)
		  {
		    endBM = 1;
		  }
		  else
		  {
		    imgPtr = imgBuf + (idY * iHead.biWidth);
		  }
		  break;
		case WLZEFF_BMP_CMP_EOB:		    /* End of bitmap */
		  endBM = 1;
		  break;
		case WLZEFF_BMP_CMP_DELTA:	  /* Delta mode, read offset */
		  /* Delta mode.  */
		  if(((tI0 = fgetc(fP)) == EOF) ||
		     ((idX += tI0) >= iHead.biWidth) ||
		     ((tI1 = fgetc(fP)) == EOF) ||
		     ((idY += tI1) >= iHead.biHeight))
		  {
		    endBM = 1;
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		  else
		  {
		    imgPtr = imgBuf + (idY * iHead.biWidth) + idX;
		  }
		  break;
		case EOF:
		  endBM = 1;
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		  break;
		default:				    /* Absolute mode */
		  padByte = rlCnt & 1;
		  while((rlCnt-- > 0) && ((tI0 = fgetc(fP)) != EOF))
		  {
		    *imgPtr++ = tI0;
		    ++idX;
		  }
		  if(tI0 == EOF)
		  {
		    endBM = 1;
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		  else
		  {
		    if(padByte)				    /* Read pad byte */
		    {
		      (void )fgetc(fP);
		    }
		  }
		  break;
	      }
	    }
	    else					    /* Absolute mode */
	    {
	      if((tI0 = fgetc(fP)) == EOF)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		while(rlCnt-- > 0)
		{
		  *imgPtr++ = tI0;
		  ++idX;
		}
	      }
	    }
	  }
	  break;
        default:
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check image size and allocate the destination data array if required. */
    if((gvnImgSz->vtX && (iHead.biWidth != gvnImgSz->vtX)) ||
       (gvnImgSz->vtY && (iHead.biHeight!= gvnImgSz->vtY)))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      if(*data == NULL)
      {
	if((AlcUnchar2Malloc(data, iHead.biHeight,
			     iHead.biWidth) != ALC_ER_NONE) ||
	   (*data == NULL))
        {
	  WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        dataBuf = *data;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Set image size and copy the BMP data into the array. */
    gvnImgSz->vtX = iHead.biWidth;
    gvnImgSz->vtY = iHead.biHeight;
    imgPtr = imgBuf;
    for(idY = 0; idY < iHead.biHeight; ++idY)
    {
      dataPtr = *(dataBuf + (iHead.biHeight - (idY + 1)));
      for(idX = 0; idX < iHead.biWidth; ++idX)
      {
        if(iHead.biBitCount == 8)        /* 8 bit image data with colour map */
	{
	  idN = (*imgPtr > iHead.biClrUsed)? iHead.biClrUsed - 1: *imgPtr;
	  tMap = rgbMap + idN;
	  *dataPtr++  = (tMap->rgbBlue + tMap->rgbGreen + tMap->rgbRed) / 3;
	  ++imgPtr;
	}
	else 				 /* 24 bit image data, no colour map */
	{
	  *dataPtr++ = (*imgPtr++ + *imgPtr++ + *imgPtr++) / 3;
	}
      }
    }
  }
  if(rgbMap)
  {
    AlcFree(rgbMap);
  }
  if(imgBuf)
  {
    AlcFree(imgBuf);
  }
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzExtFF
* \brief	Byte swaps a datum of two bytes.
* \param	data			Pointer to the datum to be byte
*					swapped.
*/
static void	WlzEffBmpSwap2(void *data)
{
  unsigned short inS,
  		outS;

  inS = *(unsigned short *)data;
  outS = (inS >> 8) | (inS << 8);
  *(unsigned short *)data = outS;
}

/*!
* \return	void
* \ingroup	WlzExtFF
* \brief	Byte swaps a datum of four bytes.
* \param	data			Pointer to the datum to be byte
*					swapped.
*/
static void	WlzEffBmpSwap4(void *data)
{
  unsigned char	tB0;
  unsigned char	 *dataP;

  dataP = (unsigned char *)data;
  tB0 = *(dataP + 0);
  *(dataP + 0) = *(dataP + 3);
  *(dataP + 3) = tB0;
  tB0 = *(dataP + 1);
  *(dataP + 1) = *(dataP + 2);
  *(dataP + 2) = tB0;
}

/*!
* \return	Size of run-length encoded data.
* \ingroup	WlzExtFF
* \brief	Run-length encodes the image data for a eight bit
*		Microsoft bitmap file.
* \param	rleData			Destination ptr for the run- length
*					encoded data.
* \param	imgData			Image data to be encoded.
* \param	imgSz			Image size.
* \param	dstErr			Destination error ptr.
*/
static int	WlzExtFFBmpEncode8(unsigned char **rleData,
				   unsigned char *imgData, WlzIVertex2 imgSz,
				   WlzErrorNum *dstErr)
{
  int		idR,
  		idX,
		idY,
  		rleBufSz,
  		rleSz = 0;
  unsigned char *imgPtr0,
		*imgPtr1,
    		*rlePtr0,
		*rlePtr1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rleBufSz = (imgSz.vtX * (imgSz.vtY + 3) + 3) * 2;
  if((rlePtr0 = AlcMalloc(sizeof(unsigned char) * rleBufSz)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {

    idR = 0;
    idY = 0;
    imgPtr0 = imgData + ((imgSz.vtY - 1) * imgSz.vtX);
    rlePtr1 = rlePtr0;
    while(idY < imgSz.vtY)
    {
      idX = 0;
      imgPtr1 = imgPtr0;
      while(idX < imgSz.vtX)
      {
	/* Find runlength. */
	idR = 1;
	while(((idX + idR) < imgSz.vtX) &&
	      (*(imgPtr1 + idR) == *imgPtr1 ) && (idR < 255))
	{
	  ++idR;
	}
	*rlePtr1++ = idR;
	*rlePtr1++ = *imgPtr1;
	imgPtr1 += idR;
        idX +=idR;
      }
      imgPtr0 -= imgSz.vtX;
      *rlePtr1++ = WLZEFF_BMP_CMP_ESC;
      *rlePtr1++ = WLZEFF_BMP_CMP_EOL;
      ++idY;
    }
    *rlePtr1++ = WLZEFF_BMP_CMP_ESC;
    *rlePtr1++ = WLZEFF_BMP_CMP_EOB;
    *rleData = rlePtr0;
    rleSz = rlePtr1 - rlePtr0;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rleSz);
}
