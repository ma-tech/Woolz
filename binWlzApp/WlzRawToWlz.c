#pragma ident "MRC HGU $Id$"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   Woolz Library					*
*   File       :   WlzRawToWlz.c					*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Fri Jan  5 10:05:47 2001				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-b] [-l] [-t#] [-u] "
	  "[-h] [-v] <width> <height> <type> [<raw-data file>]\n"
	  "\tConvert the given raw data to a woolz image.\n"
	  "\tInput can be from standard input, the width,\n"
	  "\theight and type are mandatory parameters.\n"
	  "\tWidth and height must be integer > 0\n"
	  "\tType must one of:\n"
	  "\t      1: integer (32-bit signed)\n"
	  "\t      2: short (16-bit signed)\n"
	  "\t      3: unsigned byte\n"
	  "\t      4: float (\n"
	  "\t      5: double\n"
	  "\t      6: unsigned 32-bit integer\n"
	  "\t      7: unsigned 16-bit integer\n"
	  "\t      8: signed byte\n"
	  "\tThe output type will be a suittable woolz grey-type\n"
	  "\tselected to hold the data if possible - unsigned int\n"
	  "\tmay result in loss of data\n"
	  "\tOptions are:\n"
	  "\t  -b                big-endian byte ordering (default)\n"
	  "\t  -l                little-endian byte ordering\n"
	  "\t  -h                Help - prints this usage message\n"
	  "\t  -v                verbose operation\n"
	  "",
	  proc_str);
  return;
}

typedef union _inputDataU {
  int			*i;
  short			*s;
  unsigned char		*uc;
  float			*f;
  double		*d;
  unsigned int		*ui;
  unsigned short	*us;
  char			*c;
  void			*v;
} InputDataU;
    
int main(int	argc,
	 char	**argv)
{
  FILE		*inFile;
  char 		optList[] = "blt:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*errMsg;
  int		byteOrderFlg=0;
  WlzGreyType	newpixtype;
  int		verboseFlg=0;
  int		width, height, type, depth, wlzDepth;
  InputDataU	data;
  WlzGreyP	wlzData;
  WlzObject	*obj;
  WlzPixelV	bckgrnd;
  int		i;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      byteOrderFlg = 0;
      break;

    case 'l':
      byteOrderFlg = 1;
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* check for width, height and type - all must be present */
  if( (argc - optind) <= 2 ){
    fprintf(stderr, "%s: width, height and type parameters are mandatory\n",
	    argv[0]);
    usage(argv[0]);
    return 1;
  }
  else {
    if( (width = atoi(*(argv+optind))) <= 0 ){
      fprintf(stderr, "%s: width must be > zero\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
    optind++;
    if( (height = atoi(*(argv+optind))) <= 0 ){
      fprintf(stderr, "%s: height must be > zero\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
    optind++;
    switch( type = atoi(*(argv+optind)) ){
    case 1:
      depth = sizeof(int);
      newpixtype = WLZ_GREY_INT;
      wlzDepth = sizeof(int);
      break;

    case 2:
      depth = sizeof(short);
      newpixtype = WLZ_GREY_SHORT;
      wlzDepth = sizeof(short);
      break;

    case 3:
      depth = sizeof(char);
      newpixtype = WLZ_GREY_UBYTE;
      wlzDepth = sizeof(char);
      break;

    case 4:
      depth = sizeof(float);
      newpixtype = WLZ_GREY_FLOAT;
      wlzDepth = sizeof(float);
      break;

    case 5:
      depth = sizeof(double);
      newpixtype = WLZ_GREY_DOUBLE;
      wlzDepth = sizeof(double);
      break;

    case 6:
      depth = sizeof(int);
      newpixtype = WLZ_GREY_INT;
      wlzDepth = sizeof(int);
      break;

    case 7:
      depth = sizeof(short);
      newpixtype = WLZ_GREY_INT;
      wlzDepth = sizeof(int);
      break;

    case 8:
      depth = sizeof(char);
      newpixtype = WLZ_GREY_SHORT;
      wlzDepth = sizeof(short);
      break;

    default:
      fprintf(stderr, "%s: %d is an invalid data type\n", argv[0], type);
      usage(argv[0]);
      return 1;
    }
    optind++;
  }

  /* allocate space for data */
  if( (data.v = AlcCalloc(width*height, depth)) == NULL ){
    fprintf(stderr,
	    "%s: not enough memory for the image data please"
	    "check size parameters\n",
	    argv[0]);
    return 1;
  }
  if( (wlzData.ubp = (UBYTE *) AlcCalloc(width*height, wlzDepth)) == NULL ){
    fprintf(stderr,
	    "%s: not enough memory for the image data please"
	    "check size parameters\n",
	    argv[0]);
    return 1;
  }

  /* open the raw-data file and read the data */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }
  if( fread(data.v, depth, width*height, inFile) < (width*height) ){
    fprintf(stderr,
	    "%s: not enough data in the input file please"
	    " check size parameters\n"
	    "An object will be created with the available data anyway\n",
	    argv[0]);
  }
  fclose( inFile );

  /* convert to woolz type - worry about byte ordering later */
  if( type < 7 ){
    memcpy((void *) wlzData.ubp, data.v, width*height*depth);
  }
  else if( type == 7 ){
    int		*dataPtr=wlzData.inp;
    for(i=0; i < width*height; i++, data.us++, dataPtr++){
      *dataPtr = *data.us;
    }
  }
  else {
    short	*dataPtr=wlzData.shp;
    for(i=0; i < width*height; i++, data.c++, dataPtr++){
      *dataPtr = *data.c;
    }
  }

  /* make the woolz object */
  bckgrnd.type = WLZ_GREY_INT;
  bckgrnd.v.inv = 0;
  
  if( obj = WlzMakeRect(0, height-1, 0, width-1,
			newpixtype, wlzData.inp, bckgrnd,
			NULL, NULL, &errNum) ){
    WlzWriteObj(stdout, obj);
  }

  return WLZ_ERR_NONE;
}
