#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzVerifyObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Checks each input object, fixes it if possible and
*		write to stdout.
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
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tCheck each input object, fix if possible and write to stdout\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "h";
  int		option;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char    *errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case '~': /* dummy to avoid compiler warning */
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
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  /* loop reading objects correcting as required */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL )
  {
    if((errNum = WlzVerifyObject(obj, 1)) == WLZ_ERR_NONE) {
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: error when writing verified object (%s).\n",
		       argv[0], errMsg);
      }
    }
    else {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: error when verifying object (%s).\n",
      		     argv[0], errMsg);
    }
    WlzFreeObj(obj);
  }

  return( 0 );
}
