#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzLabel.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Labels (segments) the input objects.
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
	  "Usage:\t%s [i#] [-v] [-h] [<input file>]\n"
	  "\tLabel (segment) the input objects and write the result\n"
	  "\tto stdout. Non-domain objects are ignored, the number\n"
	  "\tof segments found is written to stderr\n"
	  "\tOptions are:\n"
	  "\t  -i#       Ignore objects with number of lines < #\n"
	  "\t  -v        Verbose flag\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
#define MAXOBJS 2000    /* Was 10000 */

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  WlzObject	**objlist = NULL;
  FILE		*inFile;
  char 		optList[] = "i:vh";
  int		option;
  int		count, numobj, i, verbose = 0;
  int		ignw = -1;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'i':
      ignw = atoi(optarg);
      if( ignw < 0 ){
        fprintf(stderr, "%s: ignw = %d is invalid\n", argv[0], ignw);
        usage(argv[0]);
        return( 1 );
      }
      break;
 
    case 'v':
      verbose = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      (void )fprintf(stderr,
      		     "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  /* read objects and segment if possible */
  ignw = (ignw < 0) ? 1 : ignw;
  count = 0;
  while((obj = WlzReadObj(inFile, NULL)) != NULL) {
    count++;
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      errNum = WlzLabel(obj, &numobj, &objlist, MAXOBJS, ignw, WLZ_8_CONNECTED);
      if(errNum == WLZ_ERR_DOMAIN_TYPE) {
	errNum = WlzWriteObj(stdout, obj);
      }
      else {
	if(verbose) {
	  fprintf(stderr,"%s: writing %d objects from input object %d\n",
		  argv[0], numobj, count);
	}
	for(i=0; (i < numobj) && (errNum == WLZ_ERR_NONE); i++){
	  errNum = WlzWriteObj(stdout, *(objlist + i));
	  WlzFreeObj(*(objlist + i));
	}
      }
      break;

    default:
      errNum = WlzWriteObj(stdout, obj);
      break;

    }

    WlzFreeObj(obj);

  }
  if(errNum != WLZ_ERR_NONE) {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,"%s: failed to label object (%s).\n",
    		   argv[0],
		   errMsg);
    return( 1 );
  }

  return( 0 );
}
