#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGaussNoise.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Adds Gaussian noise to a grey value'd object.
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
	  "Usage:\t%s [-s#] [-h] [-v] [<input mask> [<input obj>]]\n"
	  "\tAdd gaussian noise to the grey value of an object.\n"
	  "\tOptions are:\n"
	  "\t  -s#       standard deviation of the noise - default 5\n"
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
  char 		optList[] = "hs:v";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzPixelV	val;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  val.type = WLZ_GREY_FLOAT;
  val.v.flv = 5.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 's':
      val.v.flv = atof(optarg);
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

  /* check for read from file */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects and convert if possible */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    if( obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ ){
      if( (errNum = WlzGaussNoise(obj, val)) == WLZ_ERR_NONE ){
	(void )WlzWriteObj(stdout, obj);
      }
      else {
	return errNum;
      }
    }
    else {
      (void )WlzWriteObj(stdout, obj);
    }

    WlzFreeObj(obj);
  }

  return 0;
}
