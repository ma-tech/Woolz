#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        Wlz3DGetSection.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Gets an arbitrary slice from a 3D object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>
#include <WlzExtFF.h>
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
	  "Usage:\t%s [-a <pitch,yaw[,roll]>] [-f <fx,fy,fz>] [-d <dist>]"
	  " [-b <parameter bibfile>] [-m <mode>] [-s <scale>]"
	  " [-o <output file>] [-u<ux,uy,uz>]"
	  " [<3D object input file>]\n"
	  "\tGet an arbitrary slice from a 3D object\n"
	  "\twriting the 2D object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -A                 Output all sections in the transformed view.\n"
	  "\t                     This will output each section with the given\n"
	  "\t                     view orientation into files using the output\n"
	  "\t                     filename as a stub. If no output file is\n"
	  "\t                     defined plane_????.wlz will be used\n"
	  "\t  -a<pitch,yaw[,roll]> viewing angles in degrees. If roll\n"
	  "\t                       is defined then the mode is \"absolute\"\n"
	  "\t  -b<bibfile>        bibfile defining the view parameters e.g.\n"
	  "\t                     from MAPaint view or warp input I/O\n"
	  "\t                     Override all other parameter input\n"
	  "\t  -f<fx,fy,fz>       fixed point position, default - (0,0,0)\n"
	  "\t  -d<dist>           distance parameter, default - 0.0\n"
	  "\t  -m<mode>           viewing mode, possible values:\n"
	  "\t                       mode = 0 - up-is-up (default)\n"
	  "\t                       mode = 1 - statue\n"
	  "\t                       mode = 2 - absolute\n"
	  "\t  -o<output file>    Output filename, default to stdout\n"
	  "\t  -s<scale>          Scale factor, default - 1.0\n"
	  "\t  -u<ux,uy,uz>       Up vector for up-is-up mode.\n"
	  "\t			  Default: (0,0,-1)\n"
	  "\t  -h                 Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFP, *outFP, *bibFP;
  char		*outFile, *bibFile;
  char 		optList[] = "Aa:b:d:f:m:o:s:u:h";
  int		option;
  int		iVal;
  int		allFlg=0;
  int		i, j;
  double	dist=0.0, pitch=0.0, yaw=0.0, roll=0.0;
  double	scale=1.0;
  WlzDVertex3	fixed={0.0,0.0,0.0};
  WlzDVertex3	up={0.0,0.0,-1.0};
  WlzThreeDViewStruct	*viewStr=NULL;
  WlzThreeDViewMode mode=WLZ_UP_IS_UP_MODE;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  char		*errMsg;
    
  /* additional defaults */
  outFile = "-";
  bibFile = NULL;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'A':
      allFlg = 1;
      break;

    case 'a':
      switch( sscanf(optarg, "%lg,%lg", &pitch, &yaw, &roll) ){
      default:
	usage(argv[0]);
	return 1;
      case 2:
	break;
      case 3:
	mode = WLZ_ZETA_MODE;
	break;
      }
      break;

    case 'b':
      bibFile = optarg;
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
      if( sscanf(optarg, "%d", &iVal) < 1 ){
	usage(argv[0]);
	return 1;
      }
      else if( mode != WLZ_ZETA_MODE ){
	switch( iVal ){
	default:
	  usage(argv[0]);
	  return 1;
	case 0:
	  mode = WLZ_UP_IS_UP_MODE;
	  break;
	case 1:
	  mode = WLZ_STATUE_MODE;
	  break;
	case 2:
	  mode = WLZ_ZETA_MODE;
	  break;
	}
      }
      break;

    case 'o':
      outFile = optarg;
      break;

    case 's':
      if( sscanf(optarg, "%lg", &scale) < 1 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'u':
      if( sscanf(optarg, "%lg,%lg,%lg", &(up.vtX), &(up.vtY),
		 &(up.vtZ)) < 3 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 0;

    }
  }

  /* check input file/stream */
  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* check output file/stream */
  if( allFlg ){
    if(strcmp(outFile, "-"))
    {
      /* strip any file extension */
      for(i=0, j=strlen(outFile); i < strlen(outFile); i++){
	if( outFile[i] == '.' ){
	  j = i;
	}
      }
      outFile[j] = '\0';
  }
    else {
      outFile = "plane";
    }
  }
  else {
    if(strcmp(outFile, "-"))
    {
      if((outFP = fopen(outFile, "w")) == NULL)
      {
	errNum = WLZ_ERR_WRITE_EOF;
      }
    }
    else
    {
      outFP = stdout;
    }
  }

  /* create view structure */
  if( viewStr = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum) ){
    viewStr->theta = yaw * WLZ_M_PI / 180.0;
    viewStr->phi = pitch * WLZ_M_PI / 180.0;
    viewStr->zeta = roll * WLZ_M_PI / 180.0;
    viewStr->dist = dist;
    viewStr->fixed = fixed;
    viewStr->up = up;
    viewStr->view_mode = mode;
    viewStr->scale = scale;
  }

  /* check bibfile - select first section parameters in the file */
  if((errNum == WLZ_ERR_NONE) && (bibFile != NULL)){
    if( bibFP = fopen(bibFile, "r") ){
      bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, bibFP);
      while((bibFileErr == BIBFILE_ER_NONE) &&
	    (strncmp(bibfileRecord->name, "Wlz3DSectionViewParams", 22))){
	BibFileRecordFree(&bibfileRecord);
	bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, bibFP);
      }
      fclose( bibFP );
      if( bibFileErr != BIBFILE_ER_NONE ){
	fprintf(stderr, "%s: error reading bibfile: %s\n", argv[0], errMsg);
	AlcFree((void *) errMsg);
	return 1;
      }
      WlzEffBibParse3DSectionViewParamsRecord(bibfileRecord, viewStr);
      BibFileRecordFree(&bibfileRecord);
    }
    else {
      fprintf(stderr, "%s: can't open parameter bibfile %s\n", argv[0],
	      bibFile);
      return 1;
    }
  }
  

  /* read objects and section if possible */
  while((errNum == WLZ_ERR_NONE) &&
        ((obj = WlzReadObj(inFP, &errNum)) != NULL))
  {
    switch( obj->type )
    {
      case WLZ_CONTOUR:
      case WLZ_3D_DOMAINOBJ:
	WlzInit3DViewStruct(viewStr, obj);
	if( allFlg ){
	  /* loop through all possible planes */
	  for(i=WLZ_NINT(viewStr->minvals.vtZ), j=0;
	      i <= WLZ_NINT(viewStr->maxvals.vtZ); i++, j++){
	    viewStr->dist = i;
	    WlzInit3DViewStruct(viewStr, obj);
	    nobj = WlzGetSectionFromObject(obj, viewStr, &errNum);
	    if( nobj != NULL){
	      char	fileBuf[256];
	      sprintf(fileBuf, "%s%04.4d.wlz", outFile, j);
	      if( outFP = fopen(fileBuf, "w") ){
		 WlzWriteObj(outFP, nobj);
		 fclose(outFP);
	      }
	    }
	    else {
	      return errNum;
	    }
	    WlzFreeObj(nobj);
	  }
	}
	else {
	  nobj = WlzGetSectionFromObject(obj, viewStr, &errNum);
	  if( nobj != NULL){
	    WlzWriteObj(outFP, nobj);
	  }
	  else {
	    return errNum;
	  }
	  WlzFreeObj(nobj);
	}
	break;

      default:
	WlzWriteObj(outFP, obj);
	break;
    }

    WlzFreeObj(obj);
  }
  if(errNum == WLZ_ERR_READ_EOF)
  {
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
}
