#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFJpeg.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Mon Aug  4 16:40:12 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzIO
* \brief        Read and write procedures for Joint Expert Group image format JPEG.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>
#include <WlzExtFF.h>

#include <string.h>

#include <jpeglib.h>
#include <setjmp.h>

/* copied from example.c in the jpeglib distribution */
struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  /*(*cinfo->err->output_message) (cinfo);*/

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}


WlzObject *WlzEffReadObjJpeg(
  FILE 		*fP,
  WlzErrorNum 	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  struct 	jpeg_decompress_struct cinfo;
  struct 	my_error_mgr jerr;
  JSAMPARRAY 	buffer;		/* Output row buffer */
  int 		row_stride;	/* physical row width in output buffer */
  int		width, height, depth;
  int		wlzDepth;
  WlzGreyType	newpixtype;
  WlzPixelV	bckgrnd;
  int		i, rOff;
  WlzGreyP	wlzData;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check input */
  if( fP == NULL ){
    errNum = WLZ_ERR_PARAM_NULL;
  }

  /* We set up the normal JPEG error routines, then override error_exit. */
  if( errNum == WLZ_ERR_NONE ){
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    /* Establish the setjmp return context for my_error_exit to use. */
    if (setjmp(jerr.setjmp_buffer)) {
      /* If we get here, the JPEG code has signaled an error.
       * We need to clean up the JPEG object, close the input file, and return.
       */
      jpeg_destroy_decompress(&cinfo);
      if( rtnObj ){
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    /* initialize the JPEG decompression object. */
    jpeg_create_decompress(&cinfo);

    /* specify data source (eg, a file) */
    jpeg_stdio_src(&cinfo, fP);

    /* read file parameters with jpeg_read_header() */
    (void) jpeg_read_header(&cinfo, TRUE);
    /* We can ignore the return value from jpeg_read_header since
     *   (a) suspension is not possible with the stdio data source, and
     *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
     * See libjpeg.doc for more info.
     */

    /* default set parameters for decompression */
    (void) jpeg_start_decompress(&cinfo);

    /* JSAMPLEs per row in output buffer */
    row_stride = cinfo.output_width * cinfo.output_components;
    /* Make a one-row-high sample array that will go away when done with image */
    buffer = (*cinfo.mem->alloc_sarray)
      ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

    /* create the appropriate woolz object, read scanlines and copy data */
    switch( cinfo.jpeg_color_space ){
    default:
    case JCS_UNKNOWN:
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;

    case JCS_GRAYSCALE:
      bckgrnd.type = WLZ_GREY_UBYTE;
#if BITS_IN_JSAMPLE == 8
      if( cinfo.data_precision != 8 ){
	errNum = WLZ_ERR_FILE_FORMAT;
      }
      else {
	wlzDepth = sizeof(char);
	newpixtype = WLZ_GREY_UBYTE;
	depth = cinfo.num_components * 8;
	bckgrnd.v.ubv = 0;
      }
#endif /* BITS_IN_SAMPLE == 8 */
#if BITS_IN_JSAMPLE == 12
      if( cinfo.data_precision != 12 ){
	errNum = WLZ_ERR_FILE_FORMAT;
      }
      else {
	wlzDepth = sizeof(short);
	newpixtype = WLZ_GREY_SHORT;
	depth = cinfo.num_components * 12;
	bckgrnd.v.shv = 0;
      }
#endif /* BITS_IN_SAMPLE == 12 */
      break;

    case JCS_YCbCr:
    case JCS_CMYK:
    case JCS_YCCK:
    case JCS_RGB:
      bckgrnd.type = WLZ_GREY_RGBA;
      cinfo.out_color_space = JCS_RGB;
#if BITS_IN_JSAMPLE == 8
      if( cinfo.data_precision != 8 ){
	errNum = WLZ_ERR_FILE_FORMAT;
      }
      else {
	wlzDepth = sizeof(int);
	newpixtype = WLZ_GREY_RGBA;
	depth = cinfo.num_components * 8;
	bckgrnd.v.inv = 0;
      }
#endif /* BITS_IN_SAMPLE == 8 */
#if BITS_IN_JSAMPLE == 12
      errNum = WLZ_ERR_UNIMPLEMENTED;
#endif /* BITS_IN_SAMPLE == 8 */
      break;

    }

    if( errNum == WLZ_ERR_NONE ){
      /* make the woolz object */
      width = cinfo.image_width;
      height = cinfo.image_height;
      if( wlzData.ubp = (UBYTE *) AlcCalloc(width*height, wlzDepth) ){
	if( rtnObj = WlzMakeRect(0, height-1, 0, width-1, newpixtype,
				  wlzData.inp, bckgrnd,
				  NULL, NULL, &errNum) ){
	  AlcErrno	errAlcNum;
	  rtnObj->values.r->freeptr = 
	    AlcFreeStackPush(rtnObj->values.r->freeptr,
			     (void *) wlzData.ubp, &errAlcNum);
	  if( errAlcNum != ALC_ER_NONE ){
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	else{
	  AlcFree((void *) wlzData.ubp);
	}
      }
      else {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   * scan object copying in scanline values */
  if( errNum == WLZ_ERR_NONE ){
    while (cinfo.output_scanline < cinfo.output_height) {
      /* jpeg_read_scanlines expects an array of pointers to scanlines.
       * Here the array is only one element long, but you could ask for
       * more than one scanline at a time if that's more convenient.
       */
      (void) jpeg_read_scanlines(&cinfo, buffer, 1);
      /* copy to the woolz object buffer */
      switch( newpixtype ){
      case WLZ_GREY_UBYTE:
	memcpy((void *) wlzData.ubp, (const void *) buffer[0],
	       sizeof(char) * width);
	wlzData.ubp += width;
	break;

      case WLZ_GREY_SHORT:
	memcpy((void *) wlzData.shp, (const void *) buffer[0],
	       sizeof(short) * width);
	wlzData.shp += width;
	break;

      case WLZ_GREY_RGBA:
	/* this we need to decode */
	rOff = 0;
	for(i=0; i < width; i++, wlzData.rgbp++, rOff += 3){
	  WLZ_RGBA_RGBA_SET(*wlzData.rgbp,
			    buffer[0][rOff],
			    buffer[0][rOff+1],
			    buffer[0][rOff+2],
			    255);
	}
	break;
      }
    }

    /* Finish decompression */
    (void) jpeg_finish_decompress(&cinfo);

    /* Release JPEG decompression object */
    jpeg_destroy_decompress(&cinfo);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

WlzErrorNum WlzEffWriteObjJpeg(
  FILE 		*fP,
  WlzObject 	*obj,
  char		*params)
{
  WlzErrorNum			errNum=WLZ_ERR_NONE;
  struct jpeg_compress_struct cinfo;
  struct my_error_mgr		jerr;
  JSAMPARRAY 	buffer;		/* Output row buffer */
  int 		row_stride;	/* physical row width in input buffer */
  int		width, height;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  int		i, j;
  WlzObject	*rectObj=NULL;
  WlzGreyType 	gType;

  /* check input */
  if( fP == NULL ){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core ){
	if( obj->values.core == NULL ){
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else {
	  WlzIBox2	cutBox;
	  
	  cutBox.xMin = obj->domain.i->kol1;
	  cutBox.yMin = obj->domain.i->line1;
	  cutBox.xMax = obj->domain.i->lastkl;
	  cutBox.yMax = obj->domain.i->lastln;
	  gType = WlzGreyTypeFromObj(obj, &errNum);
	  if( rectObj = WlzCutObjToBox2D(obj, cutBox, gType,
					 0, 0.0, 0.0, &errNum) ){
	    width = rectObj->domain.i->lastkl - rectObj->domain.i->kol1 + 1;
	    height = rectObj->domain.i->lastln - rectObj->domain.i->line1 + 1;
	  }
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_TRANS_OBJ:
      errNum = WlzEffWriteObjJpeg(fP, obj->values.obj, params);
      break;

    case WLZ_EMPTY_OBJ:
      return WLZ_ERR_NONE;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* We set up the normal JPEG error routines, then override error_exit. */
  if( errNum == WLZ_ERR_NONE ){
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    /* Establish the setjmp return context for my_error_exit to use. */
    if (setjmp(jerr.setjmp_buffer)) {
      /* If we get here, the JPEG code has signaled an error.
       * We need to clean up the JPEG object, close the input file, and return.
       */
      jpeg_destroy_compress(&cinfo);
      if( rectObj ){
	WlzFreeObj(rectObj);
      }
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, fP);

    /* set the image parameters */
    width = obj->domain.i->lastkl - obj->domain.i->kol1 + 1;
    height = obj->domain.i->lastln - obj->domain.i->line1 + 1;
    cinfo.image_width = width;
    cinfo.image_height = height;
    switch( gType ){
    case WLZ_GREY_INT:
    case WLZ_GREY_SHORT:
    case WLZ_GREY_FLOAT:
    case WLZ_GREY_DOUBLE:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;

    case WLZ_GREY_UBYTE:
      cinfo.input_components = 1;
      cinfo.in_color_space = JCS_GRAYSCALE;
      break;

    case WLZ_GREY_RGBA:
      cinfo.input_components = 3;
      cinfo.in_color_space = JCS_RGB;
      break;
    }

    if( errNum == WLZ_ERR_NONE ){
      /* set the default parameters */
      jpeg_set_defaults(&cinfo);

      /* check for other parameters */
      if( params ){
	int	quality;
	sscanf(params, "%d", &quality);
	jpeg_set_quality(&cinfo, quality, TRUE);
      }

      /* JSAMPLEs per row in output buffer */
      row_stride = cinfo.image_width * cinfo.input_components;
      /* Make a one-row-high sample array that will go away when done with image */
      buffer = (*cinfo.mem->alloc_sarray)
	((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

    }
    else {
      jpeg_destroy_compress(&cinfo);
      if( rectObj ){
	WlzFreeObj(rectObj);
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    jpeg_start_compress(&cinfo, TRUE);

    /* loop through the woolz object setting the output buffer
       note this is the "input" to the compressor */
    if((errNum = WlzInitGreyScan(rectObj, &iwsp, &gwsp)) == WLZ_ERR_NONE){
      while((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE){
	UINT val;
	switch( gType ){
	case WLZ_GREY_UBYTE:
	  for(i=0, j=0; i < width; i++){
	    buffer[0][j++] = gwsp.u_grintptr.ubp[i];
	  }
	  break;
	case WLZ_GREY_RGBA:
	  for(i=0, j=0; i < width; i++){
	    val = gwsp.u_grintptr.rgbp[i];
	    buffer[0][j++] = WLZ_RGBA_RED_GET(val);
	    buffer[0][j++] = WLZ_RGBA_GREEN_GET(val);
	    buffer[0][j++] = WLZ_RGBA_BLUE_GET(val);
	  }
	  break;
	}
	(void) jpeg_write_scanlines(&cinfo, buffer, 1);
      }

      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }

    }
    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
    WlzFreeObj(rectObj);
  }

  return errNum;
}
