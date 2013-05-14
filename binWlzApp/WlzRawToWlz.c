#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRawToWlz_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzRawToWlz.c
* \author	Richard Baldock
* \date         January 2001
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \brief	Converts raw data to a Woolz object.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzrawtowlz "WlzRawToWlz"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzrawtowlz WlzRawToWlz
\par Name
WlzRawToWlz - converta raw data to a Woolz object.
\par Synopsis
\verbatim
WlzRawToWlz  [-h] [-v] [-b] [-d#] [-l] [-o#,#,#] [-s#]
             <width> <height> [planes] <type> [<raw-data file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose output.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Big-endian byte ordering, default.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Image dimensions (2 or 3), default 2.</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Little-endian byte ordering.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>The x,y,z offsets, default 0,0,0.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Number of bytes to skip at the start of the file.</td>
  </tr>
  <tr> 
    <td><b>-X</b></td>
    <td>XXX.</td>
  </tr>
</table>
\par Description
WlzRawToWlz converts raw data to a woolz object.
Input can be from standard input, the width,
height and type are mandatory parameters.
Width and height must be integer and > 0.
Type must one of:
<table width="500" border="0">
  <tr><td> </td> </tr>
  <tr><td><b>Parameter Value</b></td> <td><b>Data type</b></td></tr>
  <tr><td>1</td> <td>integer (32-bit signed)</td></tr>
  <tr><td>2</td> <td>short (16-bit signed)</td></tr>
  <tr><td>3</td> <td>unsigned byte</td></tr>
  <tr><td>4</td> <td>float (</td></tr>
  <tr><td>5</td> <td>double</td></tr>
  <tr><td>6</td> <td>unsigned 32-bit integer</td></tr>
  <tr><td>7</td> <td>unsigned 16-bit integer</td></tr>
  <tr><td>8</td> <td>signed byte</td></tr>
</table>
\par Examples
\verbatim
\endverbatim
\par File
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref wlzextffconvert "WlzExtFFConvert(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
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
  (void )fprintf(stderr,
	  "Usage:\t%s [-b] [-d#] [-l] [-o#,#,#] [-s#] [-h] [-v]\n"
	  "\t\t<width> <height> [planes] <type> [<raw-data file>]\n"
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
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -b                little-endian byte ordering (default big-endian)\n"
	  "\t  -d<dimensions>    image dimensions (2 or 3)(default 2)\n"
	  "\t                    If d=3 then planes must be specified\n"
	  "\t  -l                little-endian byte ordering\n"
	  "\t  -oxoff,yoff,zoff  x, y, z offsets (default 0,0,0)\n"
	  "\t  -scount           Skip count bytes at start of file.\n"
	  "\t  -h                Help - prints this usage message\n"
	  "\t  -v                verbose operation\n",
	  proc_str,
	  WlzVersion());
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
  char 		optList[] = "bd:lo:s:t:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		byteOrderFlg=1;
  WlzGreyType	newpixtype;
  /* int	verboseFlg=0; */
  int		width, height, type, depth, wlzDepth, skip = 0;
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
      sscanf(optarg, "%d", &skip);
      break;

    case 's':
      sscanf(optarg, "%d,%d,%d", &offX, &offY, &offZ);
      break;

    case 't':
      fprintf(stderr,
	      "%s: sorry, this option ignored, type must be provided"
	      "as the last command line argument\n",
	      argv[0]);
      break;

    case 'v':
      /* verboseFlg = 1; */
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
  if(skip > 0) {
    (void )fseek(inFile, skip, SEEK_SET);
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
  wlzData.ubp = (WlzUByte *) **data3D;

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
    if((obj = WlzMakeRect(offY, offY+height-1, offX, offX+width-1,
			  newpixtype, wlzData.inp, bckgrnd,
			  NULL, NULL, &errNum)) != NULL){
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
