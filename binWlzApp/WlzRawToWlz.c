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
*   Date        :  Fri Jan  5 10:05:47 2001    				*
*   $Revision$			       				*
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
	  "Usage:\t%s [-b] [-d#] [-l] [-o#,#,#] [-u] "
	  "[-h] [-v] <width> <height> [planes] <type> [<raw-data file>]\n"
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
	  "\t  -d<dimensions>    image dimensions (2 or 3)(default 2)\n"
	  "\t                    If d=3 then planes must be specified\n"
	  "\t  -l                little-endian byte ordering\n"
	  "\t  -oxoff,yoff,zoff  x, y, z offsets (default 0,0,0)\n"
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
  char 		optList[] = "bd:lo:t:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*errMsg;
  int		byteOrderFlg=0;
  WlzGreyType	newpixtype;
  int		verboseFlg=0;
  int		width, height, type, depth, wlzDepth;
  InputDataU	data;
  void		***data3D;
  WlzGreyP	wlzData;
  WlzObject	*obj;
  WlzPixelV	bckgrnd;
  int		i, numDims=2, planes=1;
  int		offX=0, offY=0, offZ=0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      byteOrderFlg = 0;
      break;

    case 'd':
      if( sscanf(optarg, "%d", &numDims) == 1 ){
	switch( numDims ){
	case 2:
	case 3:
	  break;
	default:
	  fprintf(stderr, "%s:number of dimensions must be 2 or 3\n",
		  argv[0]);
	  usage(argv[0]);
	  return 1;
	}
      }
      else {
	fprintf(stderr, "%s: failed to read number of dimensions\n",
	    argv[0]);
	usage(argv[0]);
	return 1;
      }
      break;

    case 'l':
      byteOrderFlg = 1;
      break;

    case 'o':
      sscanf(optarg, "%d,%d,%d", &offX, &offY, &offZ);
      break;

    case 't':
      fprintf(stderr,
	      "%s: sorry, this option ignored, type must be provided"
	      "as the last command line argument\n",
	      argv[0]);
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
  else if( (numDims == 3) && ((argc - optind) <= 3) ){
    fprintf(stderr, "%s: width, height, planes and type parameters are mandatory\n"
	    "for images of dimension 3\n",
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
    if( numDims == 3 ){
      if( (planes = atoi(*(argv+optind))) <= 0 ){
	fprintf(stderr, "%s: planes must be > zero\n", argv[0]);
	usage(argv[0]);
	return 1;
      }
      optind++;
    }      
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

  /* allocate space for input data data */
  if( (data.v = AlcCalloc(width*height*planes, depth)) == NULL ){
    fprintf(stderr,
	    "%s: not enough memory for the image data, please"
	    "check size parameters\n",
	    argv[0]);
    return 1;
  }

  /* open the raw-data file and read the data */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }
  if( fread(data.v, depth, width*height*planes, inFile) < (width*height*planes) ){
    fprintf(stderr,
	    "%s: not enough data in the input file please"
	    " check size parameters\n"
	    "An object will be created with the available data anyway\n",
	    argv[0]);
  }
  fclose( inFile );

  /* check for byte swapping */
#if defined (__sparc) || defined (__mips)
  if( byteOrderFlg ){
    InputDataU tmpBuf;
    int i, j, k;

    tmpBuf.v = AlcCalloc(width*height*planes, depth);
    switch( depth ){
    case 2:
      swab(data.v, tmpBuf.v, width*height*planes*depth);
      break;

    case 4:
      for(i=0; i < width*height*planes; i++){
	for(j=0, k=3; j < depth; j++, k--){
	  ((char *) (tmpBuf.i + i))[j] = ((char *) (data.i + i))[k];
	}
      }
      break;

    case 8:
      for(i=0; i < width*height*planes; i++){
	for(j=0, k=7; j < depth; j++, k--){
	  ((char *) (tmpBuf.d + i))[j] = ((char *) (data.d + i))[k];
	}
      }
      break;

    case 1:
    default:
      break;
    }
    memcpy(data.v, tmpBuf.v, width*height*planes*depth);
    AlcFree(tmpBuf.v);
  }
#endif /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  if( !byteOrderFlg ){
    InputDataU tmpBuf;
    int i, j, k;

    tmpBuf.v = AlcCalloc(width*height*planes, depth);
    switch( depth ){
    case 2:
      swab(data.v, tmpBuf.v, width*height*planes*depth);
      break;

    case 4:
      for(i=0; i < width*height*planes; i++){
	for(j=0, k=3; j < depth; j++, k--){
	  ((char *) (tmpBuf.i + i))[j] = ((char *) (data.i + i))[k];
	}
      }
      break;

    case 8:
      for(i=0; i < width*height*planes; i++){
	for(j=0, k=7; j < depth; j++, k--){
	  ((char *) (tmpBuf.d + i))[j] = ((char *) (data.d + i))[k];
	}
      }
      break;

    case 1:
    default:
      break;
    }
    memcpy(data.v, tmpBuf.v, width*height*planes*depth);
    AlcFree(tmpBuf.v);
  }
#endif /* __x86 || __alpha */

  /* space for woolz image data */
  switch( wlzDepth ){
  case 1:
    if((AlcUnchar3Malloc((unsigned char ****)&data3D, planes, height, width)
	!= ALC_ER_NONE) || (data3D == NULL) ){
      fprintf(stderr,
	      "%s: not enough memory for the image data, please"
	      "check size parameters\n",
	      argv[0]);
      return 1;
    }
    break;

  case 2:
    if((AlcShort3Malloc((short ****)&data3D, planes, height, width)
	!= ALC_ER_NONE) || (data3D == NULL) ){
      fprintf(stderr,
	      "%s: not enough memory for the image data, please"
	      "check size parameters\n",
	      argv[0]);
      return 1;
    }
    break;

  case 4:
    if((AlcInt3Malloc((int ****)&data3D, planes, height, width)
	!= ALC_ER_NONE) || (data3D == NULL) ){
      fprintf(stderr,
	      "%s: not enough memory for the image data, please"
	      "check size parameters\n",
	      argv[0]);
      return 1;
    }
    break;

  case 8:
    if((AlcDouble3Malloc((double ****)&data3D, planes, height, width)
	!= ALC_ER_NONE) || (data3D == NULL) ){
      fprintf(stderr,
	      "%s: not enough memory for the image data, please"
	      "check size parameters\n",
	      argv[0]);
      return 1;
    }
    break;
  }
  wlzData.ubp = (UBYTE *) **data3D;

  /* copy to woolz type */
  if( type < 7 ){
    memcpy((void *) wlzData.ubp, data.v, width*height*planes*depth);
  }
  else if( type == 7 ){
    int		*dataPtr=wlzData.inp;
    for(i=0; i < width*height*planes; i++, data.us++, dataPtr++){
      *dataPtr = *data.us;
    }
  }
  else {
    short	*dataPtr=wlzData.shp;
    for(i=0; i < width*height*planes; i++, data.c++, dataPtr++){
      *dataPtr = *data.c;
    }
  }

  /* make the woolz object */
  bckgrnd.type = WLZ_GREY_INT;
  bckgrnd.v.inv = 0;
  
  if( numDims == 2 ){
    if( obj = WlzMakeRect(offY, offY+height-1, offX, offX+width-1,
			  newpixtype, wlzData.inp, bckgrnd,
			  NULL, NULL, &errNum) ){
      WlzWriteObj(stdout, obj);
    }
  }
  else if( numDims == 3 ){
    WlzIVertex3	origin, size;
    size.vtX = width;
    size.vtY = height;
    size.vtZ = planes;
    origin.vtX = offX;
    origin.vtY = offY;
    origin.vtZ = offZ;
    obj = WlzFromArray3D((void ***)data3D, size, origin, newpixtype,
			 newpixtype, 0.0, 1.0, 0, 1, &errNum);
    obj->domain.p->voxel_size[0] = 1.0;
    obj->domain.p->voxel_size[1] = 1.0;
    obj->domain.p->voxel_size[2] = 1.0; 
    WlzWriteObj(stdout, obj);
  }

  return WLZ_ERR_NONE;
}
