#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyInvert.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Inverts the grey-range of a grey-level object.
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
	  "Usage:\t%s [-u#] [-l#] [-h] [-v] [<input file>]\n"
	  "\tInvert the grey-range of a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tThe original grey-values are transformed\n"
	  "\tfrom their original range. The transformation is \n"
	  "\tnew_grey = u + l - old_grey. where u and l are"
	  "\tthe upper and lower values of the inversion range."
	  "\tOptions are:\n"
	  "\t  -l#       low grey value, default min value in source\n"
	  "\t  -u#       upper grey value, default max value in source\n"
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
  char 		optList[] = "l:u:hv";
  int		option;
  WlzPixelV	max, min;
  WlzPixelV	gmin, gmax;
  int		getminFlg=1, getmaxFlg=1, verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  int		objCount=0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  max.type = min.type = WLZ_GREY_DOUBLE;
  max.v.dbv = 255.0;
  min.v.dbv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'l':
      min.v.dbv = atof(optarg);
      getminFlg = 0;
      break;

    case 'u':
      max.v.dbv = atof(optarg);
      getmaxFlg = 0;
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
      if( getminFlg || getmaxFlg ){
	if( (errNum = WlzGreyRange(obj, &gmin, &gmax)) == WLZ_ERR_NONE ){
	  if( getminFlg ){
	    WlzValueConvertPixel(&min, gmin, WLZ_GREY_DOUBLE);
	  }
	  if( getmaxFlg ){
	    WlzValueConvertPixel(&max, gmax, WLZ_GREY_DOUBLE);
	  }
	}
      }

      if( errNum == WLZ_ERR_NONE ){
	if( verboseFlg ){
	  fprintf(stderr,
		  "%s:\nconverting object %d with parameters:\n"
		  "\t(l, u) = (%f, %f)\n",
		  argv[0], objCount, min.v.dbv, max.v.dbv);
	}

	errNum = WlzGreyInvertMinMax(obj, min, max);
	if( errNum == WLZ_ERR_NONE ){
	  if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
			   "%s: failed to write object (%s).\n",
			   argv[0], errMsg);
	    return(1);
	  }
	}
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
