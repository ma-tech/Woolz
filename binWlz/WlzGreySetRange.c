#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreySetRange.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Resets the grey-range of a grey-level object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
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
	  "Usage:\t%s [-u#] [-l#] [-U#] [-L#] [-h] [-v] [<input file>]\n"
	  "\tReset the grey-range of a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tThe original grey-values are linearly transformed\n"
	  "\tfrom their original range. The transformation is \n"
	  "\tnew_grey = ((U - L)/(u - l))*old_grey + L where\n"
	  "\t U and L are the new upper and lower values and \n"
	  "\t u and l are the old upper and lower values. If u and l\n"
	  "\t are set it is up to the user to ensure that the actual\n"
	  "\timage values are within the range [u,l]. Use WlzGreyRange\n"
	  "\tto check."
	  "\tOptions are:\n"
	  "\t  -l#       low grey value in source image, default min value"
	  "in source\n"
	  "\t  -L#       low grey value in dest image, default 0"
	  "\t  -u#       upper grey value in source image, default max value"
	  "in source\n"
	  "\t  -U#       upper grey value in dest image, default 0"
	  "in source\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "l:L:u:U:hv";
  int		option;
  WlzPixelV	max, min, Max, Min;
  WlzPixelV	gmin, gmax;
  int		getminFlg=1, getmaxFlg=1, verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  int		objCount=0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  Max.type = Min.type = max.type = min.type = WLZ_GREY_DOUBLE;
  Max.v.dbv = 255.0;
  Min.v.dbv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'l':
      min.v.dbv = atof(optarg);
      getminFlg = 0;
      break;

    case 'L':
      Min.v.dbv = atof(optarg);
      break;

    case 'u':
      max.v.dbv = atof(optarg);
      getmaxFlg = 0;
      break;

    case 'U':
      Max.v.dbv = atof(optarg);
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

  /* read objects and threshold if possible */
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    objCount++;
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      /* get the existing min and max grey values */
      if( (errNum = WlzGreyRange(obj, &gmin, &gmax)) == WLZ_ERR_NONE ){
	if( getminFlg ){
	  WlzValueConvertPixel(&min, gmin, WLZ_GREY_DOUBLE);
	}
	if( getmaxFlg ){
	  WlzValueConvertPixel(&max, gmax, WLZ_GREY_DOUBLE);
	}

	if( verboseFlg ){
	  fprintf(stderr,
		  "%s:\nconverting object %d with parameters:\n"
		  "\tSource (l, u) = (%f, %f)\n"
		  "\tDestination (L, U) = (%f, %f)\n",
		  argv[0], objCount, min.v.dbv, max.v.dbv,
		  Min.v.dbv, Max.v.dbv);
	}

	errNum = WlzGreySetRange(obj, min, max, Min, Max);
	if( errNum == WLZ_ERR_NONE ){
	  errNum = WlzWriteObj(stdout, obj);
	}
      }
      if( errNum != WLZ_ERR_NONE ){
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(errNum);
      }
      break;

    default:
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(1);
      }
      break;
    }
    WlzFreeObj(obj);
  }

  return 0;
}
