#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyTemplate.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Applies a template to the object values.
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
	  "Usage:\t%s [-t#] [-h] [-v] [<input template> [<input obj>]]\n"
	  "\tApply the template to the given object setting all pixels\n"
	  "\twithin the template to the object values. The object must be\n"
	  "\t a woolz grey-level object, the template can be a domain \n"
	  "\tobject, polyline or boundary list which will be filled to make\n"
	  "\t a domain. The template is always the second object if read\n"
	  "\tfrom the standard input. Pixels outside the domain of the\n"
	  "\tobject are set to the template value.\n"
	  "\tOptions are:\n"
	  "\t  -m#       the template value - default 0\n"
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
  char 		optList[] = "t:hv";
  int		option;
  WlzPixelV	tmplVal;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  tmplVal.type = WLZ_GREY_FLOAT;
  tmplVal.v.flv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 't':
      tmplVal.v.flv = atof(optarg);
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
  if( newobj = WlzGreyTemplate(obj, tmpl, tmplVal, &errNum) ){
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
