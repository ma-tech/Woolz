#pragma ident "MRC HGU $Id$"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   Woolz Library					*
*   File       :   WlzHasIntersect.c					*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Mon Dec  4 10:03:47 2000				*
*   $Revision$								*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
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
	  "Usage:\t%s [-n] [-h] [-v] [<input file>]\n"
	  "\tTest if the input objects have a non-zero intersection\n"
	  "\tOptions are:\n"
	  "\t  -n        numerical output: 0 - no intersect, 1 - intersect\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj1, *obj2;
  FILE		*inFile;
  char 		optList[] = "hnv";
  int		option;
  int		verboseFlg=0;
  int		numOutputFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'n':
      numOutputFlg = 1;
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

  /* read objects */
  if( obj1 = WlzAssignObject(WlzReadObj(inFile, NULL), NULL) ){
    if( obj2 = WlzAssignObject(WlzReadObj(inFile, NULL), NULL) ){

      if( obj1->type != obj2->type ){
	fprintf(stderr, "%s: objects must be of the same type\n", argv[0]);
	return 1;
      }

      switch( obj1->type ){
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
	if( WlzHasIntersection(obj1, obj2, &errNum) ){
	  if( numOutputFlg ){
	    fprintf(stdout, "1");
	  }
	  else {
	    fprintf(stdout, "Objects intersect\n");
	  }
	}
	else {
	  if( errNum == WLZ_ERR_NONE ){
	    if( numOutputFlg ){
	      fprintf(stdout, "0");
	    }
	    else {
	      fprintf(stdout, "Objects do not intersect\n");
	    }
	  }
	  else {
	    fprintf(stderr, "%s: some sort of error\n", argv[0]);
	  }
	}
	break;

      default:
	fprintf(stderr, "%s: wrong object types\n", argv[0]);
	break;

      }
    }
    else {
      fprintf(stderr, "%s: second object not found\n", argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: first object not found\n", argv[0]);
    return 1;
  }

  return 0;
}
