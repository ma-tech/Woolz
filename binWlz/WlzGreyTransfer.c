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
*   File       :   WlzGreyTransfer.c					*
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Wed Jan 16 13:39:56 2002				*
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
	  "Usage:\t%s [-h] [-v] [<source image> [<destination obj>]]\n"
	  "\tCopy grey values from the source object to the destination\n"
	  "\tobject within the domain of intersection. The output object\n"
	  "\thas the domain and grey-values of the destination objecr \n"
	  "\totherwise. If both objects are read from stdin then the\n"
	  "\tfirst object is the destination and the second the source.\n"
	  "\tOptions are:\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *tmpl, *newobj;
  FILE		*inFile;
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* check for read from file, note if only one file is defined then
     it is the template and the object should be on stdin */
  if( (argc - optind) >= 2 ){
    /* read the template */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (tmpl = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read template from file %s\n", argv[0],
	      *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
    optind++;

    /* read the object */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read obj from file %s\n", argv[0],
	      *(argv+optind));

      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
  }
  else if( (argc - optind) == 1 ){
    /* read the template */
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    if( (tmpl = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read template from file %s\n", argv[0],
	      *(argv+optind));
      usage(argv[0]);
      return 1;
    }
    fclose( inFile );
    inFile = stdin;
    if( (obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read obj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    /* else read objects from stdin */
    inFile = stdin;
    if( (obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read obj from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
    if( (tmpl = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) == NULL ){
      fprintf(stderr, "%s: can't read template from stdin\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }
    
  /* apply template and write resultant object */
  if( newobj = WlzGreyTransfer(obj, tmpl, &errNum) ){
    if((errNum = WlzWriteObj(stdout, newobj)) != WLZ_ERR_NONE) {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write object (%s).\n",
		     argv[0], errMsg);
      return(1);
    }
    WlzFreeObj(newobj);
  }
  WlzFreeObj(tmpl);
  WlzFreeObj(obj);

  return 0;
}
