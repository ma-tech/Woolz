#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzThreshold.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Thresholds the input grey-level objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
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
	  "Usage:\t%s [-t#] [-v#] [-H] [-L] [-h] [<input file>]\n"
	  "\tThreshold a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -H        Threshold high, keep pixels above threshold value "
	  "(default).\n"
	  "\t  -L        Threshold low, keep pixels below threshold value.\n"
	  "\t  -t#       Threshold pixel type:\n"
	  "\t            # = %d: integer (default)\n"
	  "\t                %d: short\n"
	  "\t                %d: unsigned byte\n"
	  "\t                %d: float\n"
	  "\t                %d: double\n"
	  "\t            Note -t option must precede -v\n"
	  "\t  -v#       threshold value  - integer unless -t used\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str, WLZ_GREY_INT, WLZ_GREY_SHORT, WLZ_GREY_UBYTE,
	  WLZ_GREY_FLOAT, WLZ_GREY_DOUBLE);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "HLht:v:";
  int		option;
  WlzThresholdType highLow = WLZ_THRESH_HIGH;
  WlzGreyType	threshpixtype = WLZ_GREY_INT;
  WlzPixelV	thresh;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  thresh.type = threshpixtype;
  thresh.v.inv = 170;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 't':
      switch( threshpixtype = (WlzGreyType) atoi(optarg) ){

      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_UBYTE:
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
	break;

      default:
        fprintf(stderr, "%s: grey type = %d is invalid\n",
		argv[0], threshpixtype );
        usage(argv[0]);
        return 1;

      }
      break;

    case 'v':
      switch( threshpixtype ){

      case WLZ_GREY_INT:
	thresh.v.inv = atoi(optarg);
	break;

      case WLZ_GREY_SHORT:
	thresh.v.shv = atoi(optarg);
	break;

      case WLZ_GREY_UBYTE:
	thresh.v.ubv = atoi(optarg);
	break;

      case WLZ_GREY_FLOAT:
	thresh.v.flv = atof(optarg);
	break;

      case WLZ_GREY_DOUBLE:
	thresh.v.dbv = atof(optarg);
	break;

      }
      break;

    case 'H':
      highLow = WLZ_THRESH_HIGH;
      break;

    case 'L':
      highLow = WLZ_THRESH_LOW;
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
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
	if( (nobj = WlzThreshold(obj, thresh,
			         highLow, &errNum)) == NULL ){
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr, "%s: failed to threshold object (%s).\n",
	  		 argv[0], errMsg);
	  return(1);
	}
	else {
	  if((errNum = WlzWriteObj(stdout, nobj)) != WLZ_ERR_NONE) {
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
	    		   "%s: failed to write thresholded object (%s).\n",
			   argv[0], errMsg);
	    return(1);
	  }
	  WlzFreeObj(nobj);
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
