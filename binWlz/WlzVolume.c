#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzVolume.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Calculates the volume of the input 3D objects.
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
	  "\tCalculate the volume of the input 3D woolz objects\n"
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
  int		count, vol;
  const char    *errMsg;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
    
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

  count = 0;
  while( (obj = WlzReadObj(inFile, NULL)) != NULL )
  {
    count++;

    switch( obj->type )
    {
     case WLZ_3D_DOMAINOBJ:
       if(((vol = WlzVolume(obj , &errNum)) < 0) ||
          (errNum != WLZ_ERR_NONE)){
	 (void )WlzStringFromErrorNum(errNum, &errMsg);
	 fprintf(stderr,
	 	 "%s: Object %d: error in calculating the volume (%s).\n",
		 *argv, count, errMsg);
         return 1 ;
       }
       else {
	 fprintf(stderr, "Object %d: number of voxels = %d\n", count, vol);
       }
       break;

    case WLZ_EMPTY_OBJ:
      fprintf(stderr, "Object %d: number of voxels = %d\n", count, 0);
      break;

     default:
       fprintf(stderr, "Object %d: not 3D object type\n", count);
       break;

    }

    WlzFreeObj( obj );
  }

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
