#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFTiff.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Tue Dec 18 10:50:29 2001
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        I/O functions for tiff format images
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <Wlz.h>
#include <WlzExtFF.h>
#include <tiffio.h>

#define	CVT(x)		(((x) * 255) / ((1L<<16)-1))

/************************************************************************
* Function:	WlzEffReadObjTiff					*
* Returns:	WlzObject *:		Object read from file.		*
* Purpose:	Reads a Woolz object from the given file using	 	*
*		the TIFF format. 					*
* Global refs:	-							*
* Parameters:	const char *tiffFileName: Given file name.		*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
************************************************************************/
WlzObject	*WlzEffReadObjTiff(
  const char 	*tiffFileName,
  WlzErrorNum 	*dstErr)
{
  WlzObject	*obj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  TIFF 		*tif=NULL;
  short		bitspersample;
  short		samplesperpixel;
  short		photometric;
  unsigned short *Map=NULL;
  unsigned short *redcolormap, *bluecolormap, *greencolormap;
  unsigned char red[256], green[256], blue[256];
  int		width, height, depth, numcolors;
  int		wlzDepth;
  WlzGreyType	newpixtype;
  WlzPixelV	bckgrnd;
  int		i, row, col, offset;
  float		xPosition, yPosition, xResolution, yResolution;
  int		colMin, rowMin;
  unsigned char	*buf=NULL;
  unsigned char *inp;
  WlzGreyP	wlzData;

  if((tiffFileName == NULL) || (*tiffFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if( (tif = TIFFOpen(tiffFileName, "r")) == NULL )
  {
    errNum = WLZ_ERR_READ_EOF;
  }

  if(errNum == WLZ_ERR_NONE)
  {
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
    switch( samplesperpixel ){
    case 1:
      if (bitspersample == 1){
	depth = 1;
	wlzDepth = sizeof(char);
	newpixtype = WLZ_GREY_UBYTE;
	bckgrnd.v.ubv = 0;
      }
      else if (bitspersample <= 8){
	depth = 8;
	wlzDepth = sizeof(char);
	newpixtype = WLZ_GREY_UBYTE;
	bckgrnd.v.ubv = 0;
      }
      else {
	depth = 16;
	wlzDepth = sizeof(short);
	newpixtype = WLZ_GREY_SHORT;
	bckgrnd.v.shv = 0;
      }
      break;

    case 3:
      if( bitspersample > 8 ){
	errNum = WLZ_ERR_IMAGE_TYPE;
      }
      else {
	depth = 24;
	wlzDepth = sizeof(int);
	newpixtype = WLZ_GREY_RGBA;
	bckgrnd.v.rgbv = 0;
      }
      break;

    case 4:
      if( bitspersample > 8 ){
	errNum = WLZ_ERR_IMAGE_TYPE;
      }
      else {
	wlzDepth = sizeof(int);
	newpixtype = WLZ_GREY_RGBA;
	depth = 32;
	bckgrnd.v.inv = 0;
      }
      break;

    default:
      errNum = WLZ_ERR_IMAGE_TYPE;
      break;
    }
    bckgrnd.type = newpixtype;
  }

  if(errNum == WLZ_ERR_NONE)
  {
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    numcolors = (1 << bitspersample);
    TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photometric);
    if( numcolors > 2 ){
      switch (photometric) {
      case PHOTOMETRIC_MINISBLACK:
	Map = (unsigned short *) AlcMalloc(numcolors * sizeof(unsigned short));
	for (i = 0; i < numcolors; i++)
	  Map[i] = ((numcolors-1) * i) / numcolors;
	break;

      case PHOTOMETRIC_MINISWHITE:
	Map = (unsigned short *) AlcMalloc(numcolors * sizeof(unsigned short));
	for (i = 0; i < numcolors; i++)
	  Map[i] = (numcolors-1) - (((numcolors-1) * i) / numcolors);
	break;

      case PHOTOMETRIC_RGB:
	break;

      case PHOTOMETRIC_PALETTE:
	memset(red, 0, sizeof(red));
	memset(green, 0, sizeof(green));
	memset(blue, 0, sizeof(blue));
	TIFFGetField(tif, TIFFTAG_COLORMAP,
		     &redcolormap, &greencolormap, &bluecolormap);
	for (i = 0; i < numcolors; i++) {
	  red[i] = (unsigned char) CVT(redcolormap[i]);
	  green[i] = (unsigned char) CVT(greencolormap[i]);
	  blue[i] = (unsigned char) CVT(bluecolormap[i]);
	}
	break;

      case PHOTOMETRIC_MASK:
      default:
	errNum = WLZ_ERR_IMAGE_TYPE;
	break;
      }
    }
  }

  if(errNum == WLZ_ERR_NONE)
  {
    if( (buf = (unsigned char *) AlcMalloc(TIFFScanlineSize(tif))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if(errNum == WLZ_ERR_NONE)
  {
    if( (wlzData.ubp = (UBYTE *) AlcCalloc(width*height, wlzDepth)) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    for (row = 0, offset=0; row < height; row++) {
      if (TIFFReadScanline(tif, buf, row, 0) < 0){
	errNum = WLZ_ERR_FILE_FORMAT;
	break;
      }
      inp = buf;
      switch (photometric) {
      case PHOTOMETRIC_RGB:
	if (samplesperpixel == 4){
	  for (col = 0; col < width; col++, offset++) {
	    wlzData.rgbp[offset] = 
	      (inp[3]<<24)|(inp[0]<<16)|(inp[1]<<8)|(inp[2]);
	    inp += 4;	/* skip to next values */
	  }
	}
	else {
	  for (col = 0; col < width; col++, offset++) {
	    wlzData.rgbp[offset] = 
	      (0<<24)|(inp[0]<<16)|(inp[1]<<8)|(inp[2]);
	    inp += 4;	/* skip to next values */
	  }
	}
	break;

      case PHOTOMETRIC_MINISWHITE:
      case PHOTOMETRIC_MINISBLACK:
	switch (bitspersample) {
	case 1:
	  for (col = 0; col < ((width + 7) / 8); col++, offset++){
	    wlzData.ubp[offset] = *inp++;
	  }
	  break;
	case 2:
	  for (col = 0; col < ((width + 3) / 4); col++) {
	    wlzData.ubp[offset++] = (*inp >> 6) & 3;
	    wlzData.ubp[offset++] = (*inp >> 4) & 3;
	    wlzData.ubp[offset++] = (*inp >> 2) & 3;
	    wlzData.ubp[offset++] = *inp++ & 3;
	  }
	  break;
	case 4:
	  for (col = 0; col < width / 2; col++) {
	    wlzData.ubp[offset++] = *inp >> 4;
	    wlzData.ubp[offset++] = *inp++ & 0xf;
	  }
	  break;
	case 8:
	  for (col = 0; col < width; col++, offset++){
	    wlzData.ubp[offset] = *inp++;
	  }
	  break;
	case 16:
	  for (col = 0; col < width; col++, offset++){
	    wlzData.shp[offset] = (inp[0]<<8) | (inp[1]);
	    inp += 2;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_FILE_FORMAT;
	  break;
	}
	break;

      case PHOTOMETRIC_PALETTE:
	break;

      default:
	errNum = WLZ_ERR_FILE_FORMAT;
	break;
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    if( TIFFGetField(tif, TIFFTAG_XPOSITION, &xPosition) != 1 ){
      xPosition = 0.0;
    }
    if( TIFFGetField(tif, TIFFTAG_XRESOLUTION, &xResolution) != 1 ){
      xResolution = 1.0;
    }
    if( TIFFGetField(tif, TIFFTAG_YPOSITION, &yPosition) != 1 ){
      yPosition = 0.0;
    }
    if( TIFFGetField(tif, TIFFTAG_YRESOLUTION, &yResolution) != 1 ){
      yResolution = 1.0;
    }
    colMin = WLZ_NINT(xPosition*xResolution);
    rowMin = WLZ_NINT(yPosition*yResolution);
    obj = WlzMakeRect(rowMin, rowMin+height-1, colMin, colMin+width-1,
		      newpixtype, wlzData.inp, bckgrnd,
		      NULL, NULL, &errNum);
  }  

  if( buf ){
    AlcFree(buf);
  }
  if( Map ){
    AlcFree( Map );
  }
  if( tif ){
    TIFFClose( tif );
  }
  
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return obj;
}

/************************************************************************
* Function:	WlzEffWriteObjTiff					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given Woolz object to the given file 	*
*		using the TIFF format. 					*
* Global refs:	-							*
* Parameters:	const char *tiffFileName: Given file name with .tif,	*
*		WlzObject *obj:		Given woolz object.		*
************************************************************************/
WlzErrorNum WlzEffWriteObjTiff(
  const char *tiffFileName,
  WlzObject *obj)
{
  TIFF 		*out;
  WlzObject	*rectObj=NULL;
  int		width, height;
  float		xPosition, yPosition, xResolution, yResolution;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  unsigned char	*buf=NULL;
  int		i, j;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((tiffFileName == NULL) || (*tiffFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if( obj->type != WLZ_2D_DOMAINOBJ )
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else
  {
    if( (out = TIFFOpen(tiffFileName, "w")) == NULL ){
      errNum = WLZ_ERR_FILE_OPEN;
    }
  }

  if(errNum == WLZ_ERR_NONE)
  {
    WlzIBox2	cutBox;
    WlzGreyType gType;

    cutBox.xMin = obj->domain.i->kol1;
    cutBox.yMin = obj->domain.i->line1;
    cutBox.xMax = obj->domain.i->lastkl;
    cutBox.yMax = obj->domain.i->lastln;
    gType = WlzGreyTypeFromObj(obj, &errNum);
    if((errNum == WLZ_ERR_NONE) &&
       (rectObj = WlzCutObjToBox2D(obj, cutBox, gType, 0, 0.0, 0.0, &errNum)) ){

      /* output image size parameters */
      width = rectObj->domain.i->lastkl - rectObj->domain.i->kol1 + 1;
      height = rectObj->domain.i->lastln - rectObj->domain.i->line1 + 1;
      xPosition = rectObj->domain.i->kol1;
      yPosition = rectObj->domain.i->line1;
      xResolution = 1.0;
      yResolution = 1.0;
      TIFFSetField(out, TIFFTAG_IMAGEWIDTH, (uint32) width);
      TIFFSetField(out, TIFFTAG_IMAGELENGTH, (uint32) height);
      TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
      TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
      TIFFSetField(out, TIFFTAG_XPOSITION, xPosition);
      TIFFSetField(out, TIFFTAG_YPOSITION, yPosition);
      TIFFSetField(out, TIFFTAG_XRESOLUTION, xResolution);
      TIFFSetField(out, TIFFTAG_YRESOLUTION, yResolution);

      /* out put pixel type, no compression for now */
      switch( gType ){
      case WLZ_GREY_SHORT:
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 16);
	break;

      case WLZ_GREY_UBYTE:
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
	break;

      case WLZ_GREY_RGBA:
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 4);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
	break;

      case WLZ_GREY_LONG:
      case WLZ_GREY_INT:
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
      case WLZ_GREY_BIT:
      default:
	errNum = WLZ_ERR_FILE_FORMAT;
	break;
      }
      TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

      /* write the data line by line */
      if( (buf = (unsigned char *)AlcMalloc(TIFFScanlineSize(out))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else if((errNum = WlzInitGreyScan(rectObj, &iwsp, &gwsp)) == WLZ_ERR_NONE){
	while((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE){
	  switch( gType ){
	  case WLZ_GREY_SHORT:
	    for(i=0, j=0; i < width; i++){
	      buf[j++] = ((gwsp.u_grintptr.shp[i]) & 0xff00) >> 8;
	      buf[j++] = (gwsp.u_grintptr.shp[i]) & 0xff;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    for(i=0, j=0; i < width; i++){
	      buf[j++] = gwsp.u_grintptr.ubp[i];
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    for(i=0, j=0; i < width; i++){
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0xff0000)>>16;
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0xff00)>>8;
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0xff);
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0xff000000)>>24;
	    }
	    break;
	  }

	  if( TIFFWriteScanline(out, buf, iwsp.linpos, 0) < 0 ){
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
	if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	}
	AlcFree(buf);
      }

      WlzFreeObj(rectObj);
    }
  }

  if( out ){
    TIFFClose( out );
  }
  return errNum;
}

