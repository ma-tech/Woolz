#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzTransposeObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Transposes a 2D object.
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
	  "\tTranspose a 2D woolz object\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "h";
  int		option;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  if( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

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
  while((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    case WLZ_2D_POLYGON:
    case WLZ_BOUNDLIST:
    case WLZ_EMPTY_OBJ:
      if( (nobj = WlzTransposeObj(obj, &errNum)) == NULL ){
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr, "%s: failed to transpose object (%s).\n",
		       argv[0], errMsg);
      }
      else {
	if((errNum = WlzWriteObj(stdout, nobj)) != WLZ_ERR_NONE) {
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr, "%s: failed to write object (%s).\n",
			 argv[0], errMsg);
	  return(1);
	}
	WlzFreeObj(nobj);
      }
      break;

    default:
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
        (void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr, "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(1);
      }
      break;
    }

    WlzFreeObj(obj);
  }

  return 0;
}
