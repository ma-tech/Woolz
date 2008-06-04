#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFTiff_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFTiff.c
* \author       Richard Baldock
* \date         December 2001
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
* \brief        I/O functions for tiff format images.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <Wlz.h>
#include <WlzExtFF.h>
#include <tiffio.h>

#define	CVT(x)		(((x) * 255) / ((1L<<16)-1))

static WlzErrorNum setPixelProperties(
  short		bitspersample,
  short 	samplesperpixel,
  short 	sampleformat,
  int		*wlzDepthRtn,
  WlzGreyType	*newpixtypeRtn,
  WlzPixelV	*bckgrndRtn)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		wlzDepth;
  WlzGreyType	newpixtype;
  WlzPixelV	bckgrnd;

  switch( samplesperpixel ){
  case 1:
    if (bitspersample == 1){
      wlzDepth = sizeof(char);
      newpixtype = WLZ_GREY_UBYTE;
      bckgrnd.v.ubv = 0;
    }
    else if (bitspersample <= 8){
      wlzDepth = sizeof(char);
      newpixtype = WLZ_GREY_UBYTE;
      bckgrnd.v.ubv = 0;
    }
    else if (bitspersample <= 16){
      switch( sampleformat ){
      case SAMPLEFORMAT_UINT:
	wlzDepth = sizeof(int);
	newpixtype = WLZ_GREY_INT;
	bckgrnd.v.inv = 0;
	break;

      case SAMPLEFORMAT_INT:
	wlzDepth = sizeof(short);
	newpixtype = WLZ_GREY_SHORT;
	bckgrnd.v.shv = 0;
	break;

      case SAMPLEFORMAT_IEEEFP:
	wlzDepth = sizeof(float);
	newpixtype = WLZ_GREY_FLOAT;
	bckgrnd.v.flv = 0;
	break;

      default:
	errNum = WLZ_ERR_IMAGE_TYPE;
	break;
      }
    }
    else {
      switch( sampleformat ){
      case SAMPLEFORMAT_UINT:
      case SAMPLEFORMAT_INT:
	wlzDepth = sizeof(int);
	newpixtype = WLZ_GREY_INT;
	bckgrnd.v.inv = 0;
	break;

      case SAMPLEFORMAT_IEEEFP:
	wlzDepth = sizeof(float);
	newpixtype = WLZ_GREY_FLOAT;
	bckgrnd.v.flv = 0;
	break;

      default:
	errNum = WLZ_ERR_IMAGE_TYPE;
	break;
      }
    }
    break;

  case 3:
    if( bitspersample > 8 ){
      errNum = WLZ_ERR_IMAGE_TYPE;
    }
    else {
      wlzDepth = sizeof(int);
      newpixtype = WLZ_GREY_RGBA;
      bckgrnd.v.rgbv = 0x0;
    }
    break;

  case 4:
    if( bitspersample > 8 ){
      errNum = WLZ_ERR_IMAGE_TYPE;
    }
    else {
      wlzDepth = sizeof(int);
      newpixtype = WLZ_GREY_RGBA;
      bckgrnd.v.inv = 0;
    }
    break;

  default:
    errNum = WLZ_ERR_IMAGE_TYPE;
    break;
  }
  bckgrnd.type = newpixtype;

  /* set return values */
  *wlzDepthRtn = wlzDepth;
  *newpixtypeRtn = newpixtype;
  *bckgrndRtn = bckgrnd;
 
  return errNum;
}

static unsigned char *WlzEFFTiffToWlzRowData(
  unsigned char	*inp,
  unsigned char	*dstPtr,
  int		width,
  short		photometric,
  short		samplesperpixel,
  short		bitspersample,
  WlzGreyType	newpixtype,
  unsigned char	red[],
  unsigned char	green[],
  unsigned char	blue[],
  WlzErrorNum	*dstErr)
{
  WlzGreyP	wlzData;
  int		col, offset;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( inp == NULL ){
    errNum = WLZ_ERR_VALUES_NULL;
  }
  if( dstPtr ){
    wlzData.ubp = dstPtr;
  }
  else {
    errNum = WLZ_ERR_VALUES_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) ){
    offset = 0;
    switch (photometric) {
    case PHOTOMETRIC_RGB:
      if (samplesperpixel == 4){
	for (col = 0; col < width; col++, offset++) {
	  WLZ_RGBA_RGBA_SET(wlzData.rgbp[offset],
			    inp[0], inp[1], inp[2], inp[3]);
	  inp += 4;	/* skip to next values */
	}
      }
      else {
	for (col = 0; col < width; col++, offset++) {
	  WLZ_RGBA_RGBA_SET(wlzData.rgbp[offset],
			    inp[0], inp[1], inp[2], 0xff);
	  inp += 3;	/* skip to next values */
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
	switch( newpixtype ){
	case WLZ_GREY_SHORT:
	  for (col = 0; col < width; col++, offset++){
#if defined (__sparc) || defined (__mips) || defined (__ppc)
	    wlzData.shp[offset] = (inp[0]<<8) | (inp[1]);
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
	    wlzData.shp[offset] = (inp[0]) | (inp[1]<<8);
#endif /* __x86 || __alpha */
	    inp += 2;
	  }
	  break;

	case WLZ_GREY_INT:
	  for (col = 0; col < width; col++, offset++){
#if defined (__sparc) || defined (__mips) || defined (__ppc)
	    wlzData.inp[offset] = (inp[0]<<8) | (inp[1]);
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
	    wlzData.inp[offset] = (inp[0]) | (inp[1]<<8);
#endif /* __x86 || __alpha */
	    inp += 2;
	  }
	  break;

	case WLZ_GREY_FLOAT: /* FALLTHROUGH */
	default:
	  errNum = WLZ_ERR_FILE_FORMAT;
	  break;
	}
	break;
      default:
	errNum = WLZ_ERR_FILE_FORMAT;
	break;
      }
      break;

    case PHOTOMETRIC_PALETTE:
      switch (bitspersample) {
      case 4:
	for (col = 0; col < width / 2; col++) {
	  wlzData.ubp[offset++] = *inp >> 4;
	  wlzData.ubp[offset++] = *inp++ & 0xf;
	}
	break;
      case 8:
	switch( newpixtype ){
	case WLZ_GREY_UBYTE:
	  for (col = 0; col < width; col++, offset++){
	    wlzData.ubp[offset] = *inp++;
	  }
	  break;
	case WLZ_GREY_RGBA:
	  for (col = 0; col < width; col++, offset++){
	    WLZ_RGBA_RGBA_SET(wlzData.rgbp[offset],
			      red[*inp], green[*inp], blue[*inp],
			      255);
	    inp += 1;	/* skip to next value */
	  }
	  break;
	default:
	  errNum = WLZ_ERR_FILE_FORMAT;
	}
	break;

      default:
	errNum = WLZ_ERR_FILE_FORMAT;
	break;
      }
      break;

    default:
      errNum = WLZ_ERR_FILE_FORMAT;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }

  /* set the pointer for next data */
  if( dstPtr ){
    switch( newpixtype ){
    case WLZ_GREY_UBYTE:
      dstPtr = (unsigned char *) &(wlzData.ubp[offset]);
	break;

    case WLZ_GREY_SHORT:
	dstPtr = (unsigned char *) &(wlzData.shp[offset]);
	break;
	
    default:
    case WLZ_GREY_INT:
	dstPtr = (unsigned char *) &(wlzData.inp[offset]);
	break;

    case WLZ_GREY_RGBA:
	dstPtr = (unsigned char *) &(wlzData.rgbp[offset]);
	break;

    }
  }

  return dstPtr;
}
  

static WlzObject *WlzExtFFReadTiffDirObj(
  TIFF 		*tif, 
  int		dir,
  int		split,
  WlzErrorNum	*dstErr)
{
  WlzObject	*obj=NULL, **objs;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  short		bitspersample;
  short		samplesperpixel;
  short		sampleformat;
  short		photometric;
  unsigned short *Map=NULL;
  unsigned short *redcolormap, *bluecolormap, *greencolormap;
  unsigned char red[256], green[256], blue[256];
  int		width, height, numcolors;
  int		len, tileWidth, tileHeight;
  int		wlzDepth;
  WlzGreyType	newpixtype;
  WlzPixelV	bckgrnd;
  int		i, y, row, col, tileIndx;
  float		xPosition, yPosition;
  int		colMin, rowMin;
  unsigned char	*buf=NULL;
  unsigned char *dstPtr, *srcPtr;
  WlzGreyP	wlzData;

  /* don't need to check tif file pointer since already checked
     but do need to set the directory */
  TIFFSetDirectory(tif, dir);

  /* determine depth and pixel type */
  if(errNum == WLZ_ERR_NONE)
  {
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
    if( TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sampleformat) == 0 ){
      sampleformat = SAMPLEFORMAT_UINT;
    }
    errNum = setPixelProperties(bitspersample, samplesperpixel, sampleformat,
				&wlzDepth, &newpixtype, &bckgrnd);
  }

  /* determine color mapping */
  if(errNum == WLZ_ERR_NONE)
  {
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
	  if( (red[i] != green[i]) || (red[i] != blue[i]) ){
	    wlzDepth = sizeof(int);
	    newpixtype = WLZ_GREY_RGBA;
	  }
	}
	break;

      case PHOTOMETRIC_MASK:
      default:
	errNum = WLZ_ERR_IMAGE_TYPE;
	break;
      }
    }
  }

  /* establish data size */
  if(errNum == WLZ_ERR_NONE){
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
  }

  /* check for tiles */
  if(errNum == WLZ_ERR_NONE){
    tileWidth = -1;
    tileHeight = -1;
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileHeight);
  }

  /* allocate space for the woolz data */
  if(errNum == WLZ_ERR_NONE)
  {
    if((wlzData.ubp = (WlzUByte *)AlcCalloc(width*height, wlzDepth)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  /* read data */
  if( errNum == WLZ_ERR_NONE )
  {
    if( tileWidth > 0 ){
      /* read in tiles */
      if( (buf = (unsigned char *) AlcMalloc(TIFFTileSize(tif))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if( split ){
	objs = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * TIFFNumberOfTiles(tif));
	dstPtr = wlzData.ubp;
      }

      if( errNum == WLZ_ERR_NONE ){
	tileIndx = 0;
	for (row = 0; row < height; row += tileHeight){
	  for (col = 0; col < width; col += tileWidth){
	    if( TIFFReadTile(tif, buf, col, row, 0, 0) < 0 ){
	      errNum = WLZ_ERR_FILE_FORMAT;
	      break;
	    }
	    tileIndx++;

	    /* convert and fill wlz data buffer */
	    if( split ){
	      len = ((col + tileWidth) > width) ? width - col : tileWidth;
	      obj = WlzMakeRect(0, ((row+tileHeight > height)?height-row:tileHeight)-1,
				0, len-1, newpixtype, (int *) dstPtr, bckgrnd,
				NULL, NULL, &errNum);
	      objs[tileIndx-1] = WlzAssignObject(obj, NULL);
	      if( tileIndx == 1 ){
		obj->values.r->freeptr = 
		  AlcFreeStackPush(obj->values.r->freeptr,
				   (void *) wlzData.ubp, NULL);
	      }
	      else {
		obj->assoc = WlzAssignObject(objs[0], &errNum);
	      }

	      for(y = row; (y < (row + tileHeight)) && (y < height); y++){
		srcPtr = buf + TIFFTileRowSize(tif) * (y - row);
		dstPtr =  WlzEFFTiffToWlzRowData(srcPtr, dstPtr, len,
						 photometric, samplesperpixel,
						 bitspersample, newpixtype,
						 red, green, blue, &errNum);
	      }

	    }
	    else {
	      len = ((col + tileWidth) > width) ? width - col : tileWidth;
	      for(y = row; (y < (row + tileHeight)) && (y < height); y++){
		dstPtr = &(wlzData.ubp[(y*width + col)*wlzDepth]);
		srcPtr = buf + TIFFTileRowSize(tif) * (y - row);
		(void) WlzEFFTiffToWlzRowData(srcPtr, dstPtr, len,
					      photometric, samplesperpixel,
					      bitspersample, newpixtype,
					      red, green, blue, &errNum);
	      }
	    }
	  }
	}
      }
    }
    else {
      /* use scanline read */
      if( (buf = (unsigned char *) AlcMalloc(TIFFScanlineSize(tif))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
      }

      if( errNum == WLZ_ERR_NONE ){
	dstPtr = wlzData.ubp;
	for (row = 0; row < height; row++) {

	  /* read a scanline */
	  if (TIFFReadScanline(tif, buf, row, 0) < 0){
	    errNum = WLZ_ERR_FILE_FORMAT;
	    break;
	  }

	  /* convert a row of data */
	  dstPtr = WlzEFFTiffToWlzRowData(buf, dstPtr, width, photometric, samplesperpixel,
					  bitspersample, newpixtype, red, green, blue, &errNum);
	}
      }
    }
  }

  /* I don't think we should multiply with the resolution here - RAB */
  if( errNum == WLZ_ERR_NONE ){
    if( split ){
      WlzCompoundArray	*cobj;
      cobj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 2, tileIndx, objs,
  			      WLZ_2D_DOMAINOBJ, &errNum);
      obj = (WlzObject *) cobj;
    }
    else {
      if( TIFFGetField(tif, TIFFTAG_XPOSITION, &xPosition) != 1 ){
	xPosition = 0.0;
      }
      if( TIFFGetField(tif, TIFFTAG_YPOSITION, &yPosition) != 1 ){
	yPosition = 0.0;
      }
      colMin = WLZ_NINT(xPosition);
      rowMin = WLZ_NINT(yPosition);
      if((obj = WlzMakeRect(rowMin, rowMin+height-1, colMin, colMin+width-1,
			    newpixtype, wlzData.inp, bckgrnd,
			    NULL, NULL, &errNum)) != NULL){
	AlcErrno	errAlcNum;
	obj->values.r->freeptr = 
	  AlcFreeStackPush(obj->values.r->freeptr,
			   (void *) wlzData.ubp, &errAlcNum);
	if( errAlcNum != ALC_ER_NONE ){
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }  

  if( buf ){
    AlcFree(buf);
  }
  if( Map ){
    AlcFree( Map );
  }
  
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return obj;
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given file using the TIFF format.
* \param	tiffFileName		Given file name.
* \param	split			If the image is tiled return individual
*					tiles within a compound object.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjTiff(
  const char 	*tiffFileName,
  int		split,
  WlzErrorNum 	*dstErr)
{
  WlzObject	*obj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  TIFF 		*tif=NULL;
  int		p, numPlanes;

  if((tiffFileName == NULL) || (*tiffFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if( (tif = TIFFOpen(tiffFileName, "rb")) == NULL )
  {
    errNum = WLZ_ERR_READ_EOF;
  }

  /* check for a tiff-stack.
     Not quite sure what a TIFF directory is - one day read the manual,
     but this is the way tiffsplit works and seems to assume that each
     directory is a seperate image which we take as a seperate plane.
     This could also of course be a patched object so maybe we should
     allow a compound object return option.
  */
  numPlanes = TIFFNumberOfDirectories(tif);
  if( numPlanes > 1 ){
    if( split ){
      /* create a compound object here */
      WlzCompoundArray	*cobj;

      cobj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_2, 1, numPlanes, NULL,
				  WLZ_2D_DOMAINOBJ, &errNum);
      for(p=0; p < numPlanes; p++){
	obj = WlzExtFFReadTiffDirObj(tif, p, split, &errNum);
	cobj->o[p] = WlzAssignObject(obj, &errNum);
      }
      obj = (WlzObject *) cobj;
    }
    else {
      WlzObject	*tmpObj;
      WlzDomain	domain, *domains;
      WlzValues	values, *valuess;
      WlzPixelV 	bckgrnd;

      /* use width, height and position of first object for the plane-domain
	 standardise later */
      if((tmpObj = WlzExtFFReadTiffDirObj(tif, 0, split, &errNum)) != NULL){

	/* build the planedomain and voxelvaluetable */
	if((domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					  0, numPlanes - 1,
					  tmpObj->domain.i->line1,
					  tmpObj->domain.i->lastln,
					  tmpObj->domain.i->kol1,
					  tmpObj->domain.i->lastkl,
					  &errNum)) != NULL){
	  bckgrnd = WlzGetBackground(tmpObj, &errNum);
	  if((values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					       0, numPlanes - 1,
					       bckgrnd, NULL,
					       &errNum)) != NULL){
	    domains = domain.p->domains;
	    domains[0] = WlzAssignDomain(tmpObj->domain, NULL);
	    valuess = values.vox->values;
	    valuess[0] = WlzAssignValues(tmpObj->values, NULL);
	    if((obj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
				  NULL, NULL, &errNum)) == NULL){
	      WlzFreeValues(values);
	      WlzFreeDomain(domain);
	    }
	  }
	  else{
	    WlzFreeDomain(domain);
	  }
	}
	WlzFreeObj(tmpObj);
      }

      /* put in voxel size  - have to assume plane separation of 1.0 */
      if( errNum == WLZ_ERR_NONE ){
	TIFFGetField(tif, TIFFTAG_XRESOLUTION, &(domain.p->voxel_size[0]));
	TIFFGetField(tif, TIFFTAG_YRESOLUTION, &(domain.p->voxel_size[1]));
      }

      /* now put in remaining planes */
      if( errNum == WLZ_ERR_NONE ){
	for(p=1; p < numPlanes; p++){
	  if((tmpObj = WlzExtFFReadTiffDirObj(tif, p, split, &errNum)) != NULL){
	    domains[p] = WlzAssignDomain(tmpObj->domain, NULL);
	    valuess[p] = WlzAssignValues(tmpObj->values, NULL);
	    WlzFreeObj(tmpObj);
	  }
	  else {
	    /* if it is an image-type error then it is probably some
	       proprietory information. If there is only one plane
	       then convert to 2D but still return incomplete read */
	    if((errNum == WLZ_ERR_IMAGE_TYPE) && (p == 1)){
	      tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[0],
				   valuess[0], NULL, NULL, NULL);
	      WlzFreeObj(obj);
	      obj = tmpObj;
	    }
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
	  }
	}
      }

      /* should standardise here (if 3D) */
    }
  }
  else {
    obj = WlzExtFFReadTiffDirObj(tif, 0, split, &errNum);
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

static WlzErrorNum WlzEffWriteTiffDirObj(
  TIFF		*out,
  WlzObject	*obj,
  int		dir)
{
  WlzObject	*rectObj=NULL;
  int		width, height;
  float		xPosition, yPosition, xResolution, yResolution;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  unsigned char	*buf=NULL;
  int		i, j;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIBox2	cutBox;
  WlzPixelV	minP, maxP;
  WlzGreyType gType;

  /* no checks since called only from WlzEffWriteObjTiff*/
  cutBox.xMin = obj->domain.i->kol1;
  cutBox.yMin = obj->domain.i->line1;
  cutBox.xMax = obj->domain.i->lastkl;
  cutBox.yMax = obj->domain.i->lastln;
  gType = WlzGreyTypeFromObj(obj, &errNum);

  if((errNum == WLZ_ERR_NONE) &&
     (rectObj = WlzCutObjToBox2D(obj, cutBox, gType,
				 0, 0.0, 0.0, &errNum)) ){

    /* output image size parameters */
    width = rectObj->domain.i->lastkl - rectObj->domain.i->kol1 + 1;
    height = rectObj->domain.i->lastln - rectObj->domain.i->line1 + 1;
    xPosition = WLZ_MAX(rectObj->domain.i->kol1, 0);
    yPosition = WLZ_MAX(rectObj->domain.i->line1, 0);
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
    case WLZ_GREY_INT:
      /* check if it can fit into a short and fall through
	 should really be done externally */
      if( WlzGreyRange(obj, &minP, &maxP) == WLZ_ERR_NONE ){
	if((minP.v.inv < 0) || (maxP.v.inv > (1<<15))){
	  errNum = WLZ_ERR_FILE_FORMAT;
	  break;
	}
      }
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
    case WLZ_GREY_FLOAT:
    case WLZ_GREY_DOUBLE:
    case WLZ_GREY_BIT:
    default:
      errNum = WLZ_ERR_FILE_FORMAT;
      break;
    }
    TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

    /* write the data line by line */
    if( errNum == WLZ_ERR_NONE ){
      if( (buf = (unsigned char *)AlcMalloc(TIFFScanlineSize(out))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else if((errNum = WlzInitGreyScan(rectObj, &iwsp, &gwsp))
	      == WLZ_ERR_NONE){
	while((errNum == WLZ_ERR_NONE) &&
	      (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE){
	  switch( gType ){
	  case WLZ_GREY_INT:
	    for(i=0, j=0; i < width; i++){
#if defined (__sparc) || defined (__mips) || defined (__ppc)
	      buf[j++] = ((gwsp.u_grintptr.inp[i]) & 0xff00) >> 8;
	      buf[j++] = (gwsp.u_grintptr.inp[i]) & 0xff;
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
	      buf[j++] = (gwsp.u_grintptr.inp[i]) & 0xff;
	      buf[j++] = ((gwsp.u_grintptr.inp[i]) & 0xff00) >> 8;
#endif /* __x86 || __alpha */
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    for(i=0, j=0; i < width; i++){
#if defined (__sparc) || defined (__mips) || defined (__ppc)
	      buf[j++] = ((gwsp.u_grintptr.shp[i]) & 0xff00) >> 8;
	      buf[j++] = (gwsp.u_grintptr.shp[i]) & 0xff;
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
	      buf[j++] = (gwsp.u_grintptr.shp[i]) & 0xff;
	      buf[j++] = ((gwsp.u_grintptr.shp[i]) & 0xff00) >> 8;
#endif /* __x86 || __alpha */
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    for(i=0, j=0; i < width; i++){
	      buf[j++] = gwsp.u_grintptr.ubp[i];
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    for(i=0, j=0; i < width; i++){
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0x000000ff);
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0x0000ff00)>>8;
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0x00ff0000)>>16;
	      buf[j++] = (gwsp.u_grintptr.rgbp[i]&0xff000000)>>24;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	  }

	  if(errNum == WLZ_ERR_NONE)
	  {
	    if( TIFFWriteScanline(out, buf, iwsp.linpos-iwsp.linbot, 0) < 0 ){
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      break;
	    }
	  }
	}
	if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	}
	TIFFWriteDirectory( out );
	AlcFree(buf);
      }
    }

    WlzFreeObj(rectObj);
  }

  return errNum;
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file using the TIFF
* 		format.
* \param	tiffFileName		Given file name with .tif.
* \param	obj			Given woolz object.
*/
WlzErrorNum WlzEffWriteObjTiff(
  const char *tiffFileName,
  WlzObject *obj)
{
  TIFF 		*out;
  WlzObject	*tmpObj;
  WlzDomain	*domains;
  WlzValues	*valuess;
  int		p;
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

  if(errNum == WLZ_ERR_NONE)
  {
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( (out = TIFFOpen(tiffFileName, "wb")) == NULL ){
	errNum = WLZ_ERR_FILE_OPEN;
      }
      else {
	errNum = WlzEffWriteTiffDirObj(out, obj, 0);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check data types */
      switch( obj->domain.p->type ){
      case WLZ_PLANEDOMAIN_DOMAIN:
	switch( obj->values.vox->type ){
	case WLZ_VOXELVALUETABLE_GREY:
	  /* OK 3D grey-level voxel image */
	  if( (out = TIFFOpen(tiffFileName, "w")) == NULL ){
	    errNum = WLZ_ERR_FILE_OPEN;
	  }
	  else {
	    /* should shift so that all planes are positive or zero
	       offset */
	    domains = obj->domain.p->domains;
	    valuess = obj->values.vox->values;
	    tmpObj = NULL;
	    for(p=0; p < (obj->domain.p->lastpl - obj->domain.p->plane1 + 1) &&
		  (errNum == WLZ_ERR_NONE);
		p++){
	      if( domains[p].core && valuess[p].core ){
		if( tmpObj ){
		  WlzFreeObj(tmpObj);
		  tmpObj = NULL;
		}
		if((tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
					 domains[p], valuess[p],
					 NULL, NULL, &errNum)) != NULL){
		  errNum = WlzEffWriteTiffDirObj(out, tmpObj, p);
		}
	      }
	      else {
		if( tmpObj ){
		  WlzGreySetValue(tmpObj, WlzGetBackground(tmpObj, NULL));
		  errNum = WlzEffWriteTiffDirObj(out, tmpObj, p);
		}
	      }
	    }
	  }
	  break;

	default:
	  errNum = WLZ_ERR_VALUES_TYPE;
	  break;
	}
	break;

      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      /* should be able to do these as tiff stacks */
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( out ){
    TIFFClose( out );
  }
  return errNum;
}
