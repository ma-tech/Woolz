#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzFacts.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Print out facts about Woolz objects.
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
	  "Usage:\t%s [-h] [-f] [-m] [<input file>]\n"
	  "\tPrint facts about the input woolz objects\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -f        Few facts (default)\n"
	  "\t  -m        Many facts - print out more detailed information\n"
	  "",
	  proc_str);
  return;
}
 
main(
     int		argc,
     char**	argv)
{
  WlzObject	*obj;
  FILE	*inFile;
  char 	optList[] = "fmh";
  int		option;
  int		manyFactsFlg=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  inFile = stdin;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'f':
      manyFactsFlg=0;
      break;

    case 'm':
      manyFactsFlg=1;
      break;

    case 'h':
      usage(argv[0]);
      return( WLZ_ERR_NONE );

    default:
      usage(argv[0]);
      return( WLZ_ERR_PARAM_TYPE );

    }
  }
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( WLZ_ERR_PARAM_DATA );
    }
  }

  /* read objects until EOF or an error */
  while( obj = WlzReadObj(inFile, &errNum) ){
    (void ) WlzObjectFacts(obj, stderr, NULL, manyFactsFlg);
    WlzFreeObj(obj);
    printf("\n");
  }


  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if(errNum == WLZ_ERR_READ_EOF){
    errNum = WLZ_ERR_NONE;
  }
  return errNum;
}
