#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzDiffDomain.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Calculates the difference of the input Woolz objects.
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
	  "\tCalculate the difference of the input woolz objects\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *obj1, *obj2;
  WlzObjectType	type = (WlzObjectType) -1;
  int 		n;
  FILE		*inFile;
  char 		optList[] = "h";
  int		option;
  WlzErrorNum	errNum;
  const char	*errMsg;
    
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

  /* read objects accepting the first two compatible types */
  n = 0;
  while( ((obj = WlzReadObj(inFile, NULL)) != NULL) && (n < 2) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type != type ){
      WlzFreeObj( obj );
      continue;
    }

    n++;
    if( n == 1 ){
      obj1 = obj;
    }
    else {
      obj2 = obj;
    }

  }

  if( type == (WlzObjectType) -1 ){
    fprintf(stderr, "%s: no domain objects input\n", argv[0]);
    return( WLZ_ERR_NONE );
  }

  if( n < 2 ){
    fprintf(stderr, "%s: insufficient objects of type %d\n", argv[0], type);
    return( WLZ_ERR_NONE );
  }
    
  obj = WlzDiffDomain(obj1, obj2, &errNum);
  errNum = WlzWriteObj(stdout, obj);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    fprintf(stderr, "%s: failed to write output object (%s)\n",
    	    argv[0],
	    errMsg);

  }

  /* freespace so purify can check for leaks */
  WlzFreeObj(obj1);
  WlzFreeObj(obj2);
  WlzFreeObj(obj);

  return( 0 );
}
