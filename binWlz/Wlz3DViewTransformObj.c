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
*   File       :   Wlz3DViewTransformObj.c				*
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Tue Nov  2 15:45:35 1999				*
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
	  " [-h] [-m<mode>] [<2D object input file>]\n"
	  "\tTransform an section view to a 3D object\n"
	  "\twriting the 3D object to standard output\n"
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
  char 			optList[] = "a:b:d:f:hm:v";
  int			option;
  double		dist=0.0, theta=0.0, phi=0.0, zeta=0.0;
  WlzDVertex3		fixed={0.0,0.0,0.0};
  WlzThreeDViewStruct	*viewStr=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  int			verboseFlg=0;
  WlzThreeDViewMode	viewMode=WLZ_UP_IS_UP_MODE;
  double		scale=1.0;
  WlzDVertex3		upVectorVtx={0.0,0.0,1.0};
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

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* create the view structure */
  if( errNum == WLZ_ERR_NONE ){
    if( viewStr = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum) ){
      viewStr->fixed 	= fixed;
      viewStr->theta 	= theta;
      viewStr->phi 	= phi;
      viewStr->zeta 	= zeta;
      viewStr->dist 	= dist;
      viewStr->scale 	= 1.0;
      viewStr->view_mode = viewMode;
      viewStr->up 	= upVectorVtx;
    }
    else {
      fprintf(stderr, "%s: can't create view structure\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects and section if possible */
  while((obj = WlzReadObj(inFile, &errNum)) != NULL) 
  {
    obj = WlzAssignObject(obj, &errNum);
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
	nobj = NULL;
	/* initialise the view structure and transform the section */
	WlzInit3DViewStruct(viewStr, obj);
	nobj = Wlz3DViewTransformObj(obj, viewStr, &errNum);

	if( nobj != NULL){
	  WlzWriteObj(stdout, nobj);
	}
	else {
	  return errNum;
	}
	WlzFreeObj(nobj);
	break;

      default:
	WlzWriteObj(stdout, obj);
	break;
    }

    WlzFreeObj(obj);
  }
  WlzFree3DViewStruct(viewStr);

  return WLZ_ERR_NONE;
}
