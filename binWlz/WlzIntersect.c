#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzIntersect.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Calculates the intersection of the input objects.
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
	  "\tCalculate the intersection of the input woolz objects\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -n#       Maximum number of objects -default=100\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj1, *obj, **objlist;
  WlzObjectType	type = (WlzObjectType) -1;
  int 		n, nmax;
  FILE		*inFile;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		optList[] = "n:h";
  int		option;
  int		uvt = 1;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  nmax = 100;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'n':
      nmax = atoi(optarg);
      if( nmax < 1 ){
	fprintf(stderr, "%s: nmax = %d is invalid\n", argv[0], nmax);
	usage(argv[0]);
	return( 1 );
      }
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

  /* allocate space for the object pointers */
  if( (objlist = (WlzObject **)
       AlcMalloc(sizeof(WlzObject *) * nmax)) == NULL ){
    (void )fprintf(stderr, "%s: Memory allocation failed - nmax too high?\n",
    		   argv[0]);
    return( 1 );
  }

  /* read objects accumulating compatible types */
  n = 0;
  while( ((obj = WlzReadObj(inFile, NULL)) != NULL) && (n < nmax) ) {

    if( type == -1 &&
	(obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ) ){
      type = obj->type;
    }

    if( obj->type == type ){
      objlist[n++] = obj;
      if( obj->values.core == NULL ){
	uvt = 0;
      }
    } else {
      WlzFreeObj( obj );
    }
  }

  if( type == (WlzObjectType) -1 ){
    return( WLZ_ERR_NONE );
  }

  if(((obj1 = WlzIntersectN(n, objlist, uvt, &errNum)) != NULL) &&
     (errNum == WLZ_ERR_NONE)) {
    errNum = WlzWriteObj(stdout, obj1);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr, "%s: failed to find intersection (%s).\n",
    		   argv[0], errMsg);
  }
  /* freespace so purify can check for leaks */
  WlzFreeObj(obj1);
  while( n-- ){
    WlzFreeObj(objlist[n]);
  }
  AlcFree((void *) objlist);

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
