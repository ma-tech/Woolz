#pragma ident "MRC HGU $Id$"
/*!
* \file         tif2Wlz.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Dec 14 14:39:56 2001
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
* \brief        convert a tiff image file to woolz format
*              
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Wlz.h>
#include <tiffio.h>

#define True (1)
#define False (0)
#define	CVT(x)		(((x) * 255) / ((1L<<16)-1))

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s -h -v"
	  "<input file>\n"
	  "\tRead the tiff file and convert to woolz\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str);
  return;
}


main(argc, argv)
  int         argc;
  char       *argv[];
{
  char       *inf = NULL;
  char       *outf = NULL;
  long        width,
    height;
  int         depth,
    numcolors;
  TIFF 	*tif;
  unsigned char	*inp;
  unsigned char	*outp;
  int 	col,
    i;
  long 	row;
  short     *Map = NULL;
  unsigned char     *buf;
  short	bitspersample;
  short	samplesperpixel;
  short	photometric;
  unsigned short    *redcolormap,
    *bluecolormap,
    *greencolormap;

  unsigned char      red[256],
    green[256],
    blue[256];
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzGreyP	wlzData;
  int		wlzDepth;
  int		offset;
  WlzGreyType	newpixtype;
  WlzObject	*obj;
  WlzPixelV	bckgrnd;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
      usage(argv[0]);
      return WLZ_ERR_NONE;

    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  if( optind < argc ){
    inf = *(argv+optind);
  }
  else {
    fprintf(stderr, "%s: can't read input file from a stream.\n", argv[0]);
    return 1;
  }

  if (verboseFlg)
    fprintf(stderr, "Reading %s...\n", inf);

  tif = TIFFOpen(inf, "r");

  if (tif == NULL)
    fprintf(stderr, "%s: error opening TIFF file %s\n", argv[0], inf);

  if (verboseFlg)
    TIFFPrintDirectory(tif, stderr, TIFFPRINT_NONE);
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);

  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
  switch (samplesperpixel) {
  case 1:
    if (bitspersample == 1){
      depth = 1;
      wlzDepth = sizeof(char);
      newpixtype = WLZ_GREY_UBYTE;
    }
    else if (bitspersample <= 8){
      depth = 8;
      wlzDepth = sizeof(char);
      newpixtype = WLZ_GREY_UBYTE;
    }
    else {
      depth = 16;
      wlzDepth = sizeof(short);
      newpixtype = WLZ_GREY_SHORT;
    }
    break;
  case 3:
    if( bitspersample > 8 ){
      fprintf(stderr, "%s: can't handle > 8bit colour\n", argv[0]);
      return 1;
    }
    depth = 24;
    wlzDepth = sizeof(int);
    newpixtype = WLZ_GREY_RGBA;
    break;
  case 4:
    if( bitspersample > 8 ){
      fprintf(stderr, "%s: can't handle > 8bit colour\n", argv[0]);
      return 1;
    }
    wlzDepth = sizeof(int);
    newpixtype = WLZ_GREY_RGBA;
    depth = 32;
    break;
  default:
    fprintf(stderr, "%s: only handle 1-channel gray scale or 3-channel color\n", argv[0]);
    return 1;
  }

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

  if (verboseFlg)
    fprintf(stderr, "%dx%dx%d image, \n", width, height, depth);
  if (verboseFlg)
    fprintf(stderr, "%d bits/sample, %d samples/pixel,\n",
	    bitspersample, samplesperpixel);

  numcolors = (1 << bitspersample);

  TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photometric);
  if (numcolors == 2) {
    if (verboseFlg)
      fprintf(stderr, "monochrome \n");
  } else {
    switch (photometric) {
    case PHOTOMETRIC_MINISBLACK:
      if (verboseFlg)
	fprintf(stderr, "%d graylevels (min=black), \n", numcolors);
      Map = (short *) AlcMalloc(numcolors * sizeof(short));
      for (i = 0; i < numcolors; i++)
	Map[i] = ((numcolors-1) * i) / numcolors;
      break;
    case PHOTOMETRIC_MINISWHITE:
      if (verboseFlg)
	fprintf(stderr, "%d graylevels (min=white), \n", numcolors);
      Map = (short *) AlcMalloc(numcolors * sizeof(short));
      for (i = 0; i < numcolors; i++)
	Map[i] = (numcolors-1) - (((numcolors-1) * i) / numcolors);
      break;
    case PHOTOMETRIC_RGB:
      if (verboseFlg)
	fprintf(stderr, "truecolor \n");
      break;
    case PHOTOMETRIC_PALETTE:
      if (verboseFlg)
	fprintf(stderr, "colormapped \n");
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
      fprintf(stderr, "%s: Don't know how to handle PHOTOMETRIC_MASK\n", argv[0]);
      return 1;
    default:
      fprintf(stderr, "%s: unknown photometric (cmap): %d\n", argv[0], photometric);
    }
  }

  if( verboseFlg ){
    fprintf(stderr, "%s: ScanlineSize = %d\n", argv[0], TIFFScanlineSize(tif));
  }
  buf = (unsigned char *) AlcMalloc(TIFFScanlineSize(tif));
  if (buf == NULL)
    fprintf(stderr, "%s: can't allocate memory for scanline buffer...\n", argv[0]);

  /* allocate space for the data */
  if( (wlzData.ubp = (UBYTE *) AlcCalloc(width*height, wlzDepth)) == NULL ){
    fprintf(stderr,
	    "%s: not enough memory for the image data please"
	    "check size parameters\n",
	    argv[0]);
    return 1;
  }

  for (row = 0, offset=0; row < height; row++) {
    if (TIFFReadScanline(tif, buf, row, 0) < 0)
      fprintf(stderr, "%s: bad data read on line: %d\n", argv[0], row);
    inp = buf;
    switch (photometric) {
    case PHOTOMETRIC_RGB:
      if (samplesperpixel == 4)
	for (col = 0; col < width; col++, offset++) {
	  wlzData.inp[offset] = 
	    (inp[3]<<24)|(inp[0]<<16)|(inp[1]<<8)|(inp[2]);
	  inp += 4;	/* skip to next values */
	}
      else
	for (col = 0; col < width; col++, offset++) {
	  wlzData.inp[offset] = 
	    (0<<24)|(inp[0]<<16)|(inp[1]<<8)|(inp[2]);
	  inp += 4;	/* skip to next values */
	}
      break;
    case PHOTOMETRIC_MINISWHITE:
    case PHOTOMETRIC_MINISBLACK:
      switch (bitspersample) {
      case 1:
	for (col = 0; col < ((width + 7) / 8); col++, offset++)
	  wlzData.ubp[offset] = *inp++;
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
	for (col = 0; col < width; col++, offset++)
	  wlzData.ubp[offset] = *inp++;
	break;
      case 16:
	for (col = 0; col < width; col++, offset++){
	  wlzData.shp[offset] = (inp[0]<<8) | (inp[1]);
	  inp += 2;
	}
	break;
      default:
	fprintf(stderr, "%s: bad bits/sample: %d\n", argv[0], bitspersample);
	return 1;
      }
      break;
    case PHOTOMETRIC_PALETTE:
      break;
    default:
      fprintf(stderr, "%s: unknown photometric (write): %d\n", argv[0],
	      photometric);
    }
  }

  free((char *) buf);

  /* create woolz object and write to stdout */
  bckgrnd.type = WLZ_GREY_INT;
  bckgrnd.v.inv = 0;
  
  if( obj = WlzMakeRect(0, height-1, 0, width-1,
			newpixtype, wlzData.inp, bckgrnd,
			NULL, NULL, &errNum) ){
    WlzWriteObj(stdout, obj);
  }

  if (verboseFlg)
    fprintf(stderr, "done.\n");

  exit(0);
}
