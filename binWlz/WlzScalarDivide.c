#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzScalarDivide.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Divides the pixel values of a grey-level object by a
*		scalar value.
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
	  "Usage:\t%s [-v#] [-h] [<input file>]\n"
	  "\tDivide the pixel values of a grey-level woolz object\n"
	  "\tby a scalar value\n"
	  "\tOptions are:\n"
	  "\t  -v#       divisor value (0 is an error), default 1.0:\n"
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
  char 		optList[] = "v:";
  int		option;
  double	div;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  div = 1.0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      div = atof(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_EXIT_FAILURE;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_EXIT_FAILURE;
    }
  }

  /* read objects and threshold if possible */
  while((obj = WlzReadObj(inFile)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
      if( WlzScalarDivide(obj, div, obj) == WLZ_ERR_NONE ){
	WlzWriteObj(stdout, obj);
      }
      break;

    default:
      WlzWriteObj(stdout, obj);
      break;
    }

    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
