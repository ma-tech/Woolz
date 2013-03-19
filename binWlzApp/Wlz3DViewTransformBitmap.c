#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _Wlz3DViewTransformBitmap_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/Wlz3DViewTransformBitmap.c
* \author       Richard Baldock
* \date         September 2001
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
* \brief	Transform a bitmap on a section view to a 3D object.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlz3dviewtransformbitmap "Wlz3DViewTransformBitmap"
*/

/*!
\ingroup BinWlzApp
\defgroup wlz3dviewtransformbitmap Wlz3DViewTransformBitmap
\par Name
Wlz3DViewTransformBitmap - transform a bitmap on a section view to a 3D object.
\par Synopsis
\verbatim
Wlz3DViewTransformBitmap  [-h] [-v]
                 [-a <pitch>,<yaw[>,<roll]>] [-f <fx>,<fy>,<fz>]
                 [-d <dist>] [-b <bibfile>] [-m <mode>]
		 -s <width>,<height> -o <x_off>,<y_off>
		 [<bitmap data file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Viewing angles in degrees, default (0.0, 0.0, 0.0).</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Bibfile defining the view parameters e.g. from MAPaint
        or warp input I/O.
	Override all other parameter input.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Fixed point position, default (0,0,0).</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Distance parameter, default 0.0.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Viewing mode, possible values:
    <table width="500" border="0">
      <tr> <td><b>Parameter value</b></td> <td><b>Viewing mode</b></td> </tr>
      <tr> <td>"up-is-up"</td> <td>up-is-up, default</td> </tr>
      <tr> <td>"statue"</td> <td>statue</td> </tr>
      <tr> <td>"absolute"</td> <td>absolute</td> </tr>
    </table>
    </td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Size of the bitmap.
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>The x and y offsets of the bitmap.</td>
  </tr>
</table>
\par Description
Transforms a bitmap on a section view to a 3D object writing the 3D
object to standard output.
\par Examples
\verbatim
\endverbatim
\par File
\ref Wlz3DViewTransformBitmap.c "Wlz3DViewTransformBitmap.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref Wlz3DViewTransformBitmap "Wlz3DViewTransformBitmap(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Wlz.h>
#include <bibFile.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  (void )fprintf(stderr,
	  "Usage:\t%s [-a<pitch>,<yaw[>,<roll]>] [-b<bibfile>]\n"
	  "\t[-f<fx>,<fy>,<fz>] [-d<dist>] -o <x_off>,<y_off>\n"
	  "\t[-h] [-m<mode>] -s<width>,<height> <bitmap data file>\n"
	  "\tTransform a bitmap on a section view to a 3D object\n"
	  "\twriting the 3D object to standard output\n"
	  "Version: %s\n"
	  "Required are:\n"
	  "\t  -s<width,height>   size of bitmap image\n"
	  "\t  -o<x_off,y_off>    x and y offsets of the bitmap image\n"
	  "Options:\n"
	  "\t  -a<pitch,yaw[,roll]> viewing angles in degrees - default 0.0\n"
	  "\t  -b<view-bibfile>   input parameters from the view\n"
	  "\t                     bibfile - e.g. saved from MAPaint\n"
	  "\t  -f<fx,fy,fz>       fixed point position - default zero\n"
	  "\t  -d<dist>           distance parameter - default zero\n"
	  "\t  -m<mode>           viewing mode, one of: up-is-up, statue, absolute\n"
	  "\t  -u<ux,uy,uz>       up vector - default (0.0, 0.0, 1.0)\n"
	  "\t  -h                 Help - prints this usage message\n"
	  "\t  -v                 Verbose operation\n",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject		*obj;
  FILE			*inFile;
  char 			optList[] = "a:b:d:f:hm:o:s:v";
  int			option;
  double		dist=0.0, theta=0.0, phi=0.0, zeta=0.0;
  WlzDVertex3		fixed={0.0,0.0,0.0};
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  int			verboseFlg=0;
  WlzThreeDViewMode	viewMode=WLZ_UP_IS_UP_MODE;
  double		scale=1.0;
  WlzDVertex3		upVectorVtx={0.0,0.0,1.0};
  int			width, height, xOffset, yOffset, numBytes;
  int			offsetFlg=1, sizeFlg=1;
  WlzUByte		*bitData=NULL;
  char			*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'a':
      if( sscanf(optarg, "%lg,%lg,%lg", &phi, &theta, &zeta) < 2 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'b':
      if((inFile = fopen(optarg, "r")) != NULL){
	char		viewModeStr[64];
	int		numParsedFields=0;
	BibFileRecord	*bibfileRecord;
	BibFileField	*bibfileField;
	BibFileError	bibFileErr;

	/* read the bibfile - get the first section view entry */
	bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
	while((bibFileErr == BIBFILE_ER_NONE) &&
	      (strncmp(bibfileRecord->name, "Wlz3DSectionViewParams", 22))){
	  BibFileRecordFree(&bibfileRecord);
	  bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
	}
	(void) fclose(inFile);
	if( bibFileErr != BIBFILE_ER_NONE ){
	  (void )fprintf(stderr,
	                 "%s: can't read bibfile: %s\n", argv[0], optarg);
	  usage(argv[0]);
	  return 1;
	}

	/* parse the record */
	numParsedFields = BibFileFieldParseFmt
	  (bibfileRecord->field,
	   (void *) &(fixed.vtX), "%lg ,%*lg ,%*lg", "FixedPoint",
	   (void *) &(fixed.vtY), "%*lg ,%lg ,%*lg", "FixedPoint",
	   (void *) &(fixed.vtZ), "%*lg ,%*lg ,%lg", "FixedPoint",
	   (void *) &dist, "%lg", "Distance",
	   (void *) &phi, "%lg", "Pitch",
	   (void *) &theta, "%lg", "Yaw",
	   (void *) &zeta, "%lg", "Roll",
	   (void *) &scale, "%lg", "Scale",
	   (void *) &(upVectorVtx.vtX), "%lg ,%*lg ,%*lg", "UpVector",
	   (void *) &(upVectorVtx.vtY), "%*lg ,%lg ,%*lg", "UpVector",
	   (void *) &(upVectorVtx.vtZ), "%*lg ,%*lg ,%lg", "UpVector",
	   (void *) &(viewModeStr[0]), "%s", "ViewMode",
	   NULL);
	if(numParsedFields < 1) {
	  (void )fprintf(stderr,
	                 "%s: can't read section transform from bibfile: %s\n",
	                 argv[0], optarg);
	  usage(argv[0]);
	  return 1;
	}
	/* doesn't read the view mode correctly - ask Bill */
	bibfileField = bibfileRecord->field;
	while( bibfileField ){
	  if( strncmp(bibfileField->name, "ViewMode", 8) == 0 ){
	    strcpy(viewModeStr, bibfileField->value);
	    break;
	  }
	  bibfileField = bibfileField->next;
	}
	BibFileRecordFree(&bibfileRecord);

	/* convert angles to radians and get mode */
	phi *= (WLZ_M_PI/180.0);
	theta *= (WLZ_M_PI/180.0);
	zeta *= (WLZ_M_PI/180.0);
	if( !strncmp(viewModeStr, "up-is-up", 8) ){
	  viewMode = WLZ_UP_IS_UP_MODE;
	}
	else if( !strncmp(viewModeStr, "statue", 8) ){
	  viewMode = WLZ_STATUE_MODE;
	}
	else {
	  viewMode = WLZ_ZETA_MODE;
	}


      }
      else {
	(void )fprintf(stderr,
	               "%s: can't open bibfile: %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
      break;

    case 'd':
      if( sscanf(optarg, "%lg", &dist) < 1 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'f':
      if( sscanf(optarg, "%lg,%lg,%lg", &(fixed.vtX), &(fixed.vtY),
		 &(fixed.vtZ)) < 3 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'm':
      if( strcmp(optarg, "up-is-up") == 0 ){
	viewMode = WLZ_UP_IS_UP_MODE;
      }
      else if( strcmp(optarg, "statue") == 0 ){
	viewMode = WLZ_STATUE_MODE;
      }
      else if( strcmp(optarg, "absolute") == 0 ){
	viewMode = WLZ_ZETA_MODE;
      }
      else {
	(void )fprintf(stderr, "%s: invalid view mode: %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
      break;

    case 'o':
      if( sscanf(optarg, "%d,%d", &xOffset, &yOffset) < 2 ){
	usage(argv[0]);
	return 1;
      }
      offsetFlg = 0;
      break;

    case 's':
      if( sscanf(optarg, "%d,%d", &width, &height) < 2 ){
	usage(argv[0]);
	return 1;
      }
      sizeFlg = 0;
      break;

    case 'u':
      if( sscanf(optarg, "%lg,%lg,%lg", &(upVectorVtx.vtX), &(upVectorVtx.vtY),
		 &(upVectorVtx.vtZ)) < 3 ){
	usage(argv[0]);
	return 1;
      }
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

  /* check size and offset have been set */
  if( offsetFlg || sizeFlg ){
    (void )fprintf(stderr,
                   "%s: bitmap size and offset must be set\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* get the bitmap data */
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      (void )fprintf(stderr,
                     "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    numBytes = (int) (((double) width*height)/8.0 + 1.0);
    if((bitData = AlcMalloc(sizeof(char) * (numBytes + 1) )) != NULL){
      if( fread(bitData, sizeof(char), numBytes, inFile) < numBytes ){
	(void )fprintf(stderr,
	               "%s: not enough data\n", argv[0]);
      }
    }
    else {
      (void )fprintf(stderr,
                     "%s: can't alloc memory for bit data\n", argv[0]);
      return 1;
    }
  }
  else {
    (void )fprintf(stderr, "%s: please supply a bitmap file\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  if(verboseFlg) {
    (void )fprintf(stderr,
		   "%s: Bitmap and section parameters are\n"
                   "  width = %d\n"
		   "  height = %d\n"
		   "  xOffset = %d\n"
		   "  yOffset = %d\n"
		   "  fixed = %lg,%lg,%lg\n"
		   "  theta = %lg\n"
		   "  phi = %lg\n"
		   "  dist = %lg\n"
		   "  viewMode = %s\n",
		   argv[0],
		   width,
		   height,
		   xOffset,
		   yOffset,
		   fixed.vtX, fixed.vtY, fixed.vtZ,
		   theta,
		   phi,
		   dist,
		   WlzStringFromThreeDViewMode(viewMode, NULL));
  }

  /* read bitmap data and section if possible */
  if((obj = Wlz3DViewTransformBitmap((width / 8) * height, bitData,
  				     width, height,
				     xOffset, yOffset,
				     fixed.vtX, fixed.vtY, fixed.vtZ,
				     theta, phi, dist, &errNum)) != NULL){
    WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
  }
  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
