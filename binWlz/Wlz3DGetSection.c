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

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-t<theta,phi>] [-f<fx,fy,fz>] [-d<dist>]"
	  " [<3D object input file>]\n"
	  "\tGet an arbitrary slice from a 3D object\n"
	  "\twriting the 2D object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -a<theta,phi>      viewing angles in degrees\n"
	  "\t  -f<fx,fy,fz>       fixed point position\n"
	  "\t  -d<dist>           distance parameter\n"
	  "\t  -h                 Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFP,
  		*outFP;
  char		*outFile;
  char 		optList[] = "a:d:f:o:h";
  int		option;
  double	dist=0.0, theta=0.0, phi=0.0;
  WlzDVertex3	fixed={0.0,0.0,0.0};
  WlzThreeDViewStruct	*viewStr=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
    
  outFile = "-";
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'a':
      if( sscanf(optarg, "%lg,%lg", &theta, &phi) < 2 ){
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

    case 'o':
      outFile = optarg;
      break;
    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }
  if(strcmp(outFile, "_"))
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

  /* read objects and section if possible */
  while((errNum == WLZ_ERR_NONE) &&
        ((obj = WlzReadObj(inFP, &errNum)) != NULL))
  {
    switch( obj->type )
    {
      case WLZ_CONTOUR:
      case WLZ_3D_DOMAINOBJ:
	if( viewStr == NULL ){
	  viewStr = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum);
	  viewStr->theta = theta * WLZ_M_PI / 180.0;
	  viewStr->phi = phi * WLZ_M_PI / 180.0;
	  viewStr->dist = dist;
	  viewStr->fixed = fixed;
	  WlzInit3DViewStruct(viewStr, obj);
	}
	nobj = WlzGetSectionFromObject(obj, viewStr, &errNum);
	if( nobj != NULL){
	  WlzWriteObj(outFP, nobj);
	}
	else {
	  return errNum;
	}
	WlzFreeObj(nobj);
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
