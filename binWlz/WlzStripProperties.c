#pragma ident "MRC HGU $Id$"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   MRC HGU Image Processing Utilities			*
*   File       :   WlzStripProperties.c					*
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Thu Aug  1 11:53:58 2002				*
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
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tStrip the property list from a woolz object.\n"
	  "\tThis could be useful if you need to read an object\n"
	  "\tfrom a newer version of woolz with old code e.g. if\n"
	  "\tyou have binaries for a different architecture than\n"
	  "\tsupplied. Object written to stdout\n"
	  "",
	  proc_str);
  return;
}

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'h':
      usage(argv[0]);
      return( 0 );

    case '-v':
      verboseFlg = 1;
      break;

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
      return 1;
    }
  }

  /* read objects and strip properties */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    case WLZ_3D_WARP_TRANS:
    case WLZ_PROPERTY_OBJ:
      if( obj->plist ){
	WlzFreePropertyList(obj->plist);
	obj->plist = NULL;
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      if( ((WlzCompoundArray *) obj)->p ){
	WlzFreePropertyList(((WlzCompoundArray *) obj)->p);
	((WlzCompoundArray *) obj)->p = NULL;
      }
      break;

    default:
      break;
    }

    (void )WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
		

