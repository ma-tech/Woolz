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
*   File       :   Wlz3DViewTransformBitmap.c				*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Thu Sep 20 08:36:53 2001				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

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
  fprintf(stderr,
	  "Usage:\t%s [-t<pitch,yaw[,roll]>] [-b<bibfile>]"
	  "[-f<fx,fy,fz>] [-d<dist>]"
	  " [-h] [-m<mode>] -s<width>,<height> <bitmap data file>\n"
	  "\tTransform a bitmap on a section view to a 3D object\n"
	  "\twriting the 3D object to standard output\n"
	  "\tRequired are:\n"
	  "\t  -s<width,height>   size of bitmap image\n"
	  "\t  -o<x_off,y_off>    x and y offsets of the bitmap image\n"
	  "\tOptions are:\n"
	  "\t  -a<pitch,yaw[,roll]> viewing angles in degrees - default 0.0\n"
	  "\t  -b<view-bibfile>   input parameters from the view\n"
	  "\t                     bibfile - e.g. saved from MAPaint\n"
	  "\t  -f<fx,fy,fz>       fixed point position - default zero\n"
	  "\t  -d<dist>           distance parameter - default zero\n"
	  "\t  -m<mode>           viewing mode, one of: up-is-up, statue, absolute\n"
	  "\t  -u<ux,uy,uz>       up vector - default (0.0, 0.0, 1.0)\n"
	  "\t  -h                 Help - prints this usage message\n"
	  "\t  -v                 Verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject		*obj, *nobj;
  FILE			*inFile;
  char 			optList[] = "a:b:d:f:hm:o:s:v";
  int			option;
  double		dist=0.0, theta=0.0, phi=0.0, zeta=0.0;
  WlzDVertex3		fixed={0.0,0.0,0.0};
  WlzThreeDViewStruct	*viewStr=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  int			verboseFlg=0;
  WlzThreeDViewMode	viewMode=WLZ_UP_IS_UP_MODE;
  double		scale=1.0;
  WlzDVertex3		upVectorVtx={0.0,0.0,1.0};
  int			width, height, xOffset, yOffset, numBytes;
  int			offsetFlg=1, sizeFlg=1;
  UBYTE			*bitData=NULL;
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
      if( inFile = fopen(optarg, "r") ){
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
	  fprintf(stderr, "%s: can't read bibfile: %s\n", argv[0], optarg);
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
	fprintf(stderr, "%s: can't open bibfile: %s\n", argv[0], optarg);
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
	fprintf(stderr, "%s: invalid view mode: %s\n", argv[0], optarg);
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
    fprintf(stderr, "%s: bitmap size and offset must be set\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* get the bitmap data */
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    numBytes = (int) (((double) width*height)/8.0 + 1.0);
    if( bitData = AlcMalloc(sizeof(char) * (numBytes + 1) ) ){
      if( fread(bitData, sizeof(char), numBytes, inFile) < numBytes ){
	fprintf(stderr, "%s: not enough data\n", argv[0]);
      }
    }
    else {
      fprintf(stderr, "%s: can't alloc memory for bit data\n", argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: please supply a bitmap file\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* read bitmap data and section if possible */
  if( obj = Wlz3DViewTransformBitmap(bitData, width, height,
				     xOffset, yOffset,
				     fixed.vtX, fixed.vtY, fixed.vtZ,
				     theta, phi, dist, &errNum) ){
    WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
  }
  else {
  }

  return WLZ_ERR_NONE;
}
