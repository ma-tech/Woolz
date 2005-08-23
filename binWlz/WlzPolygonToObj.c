#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzPolygonToObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Converts a 2D or 3D polygon object to the 
*		corresponding domain * object.
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
	  "\tConvert a 2D or 3D polygon woolz object\n"
	  "\tto the corresponding domain object\n"
	  "\tOptions are:\n"
	  "\t  -t#       Polygon fill-mode:\n"
	  "\t            # = %d: simple fill (default)\n"
	  "\t                %d: even-odd fill\n"
	  "\t                %d: vertex-fill\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str, WLZ_SIMPLE_FILL, WLZ_EVEN_ODD_FILL,
	  WLZ_VERTEX_FILL);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "f:hv";
  int		option;
  int		verboseFlg=0;
  WlzPolyFillMode	fillMode=WLZ_SIMPLE_FILL;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  if( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'f':
      switch( fillMode = (WlzPolyFillMode) atoi(optarg) ){
      case WLZ_SIMPLE_FILL:
      case WLZ_EVEN_ODD_FILL:
      case WLZ_VERTEX_FILL:
	break;

      default:
	fprintf(stderr, "%s: grey typefill-mode = %d is invalid\n",
		argv[0], fillMode);
        usage(argv[0]);
        return 1;
      }
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

  /* read objects and threshold if possible */
  while((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_2D_POLYGON:
      if( (nobj = WlzPolygonToObj(obj, WLZ_SIMPLE_FILL,
      				NULL)) != NULL ){
	(void )WlzWriteObj(stdout, nobj);
	WlzFreeObj(nobj);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.p->type == WLZ_PLANEDOMAIN_POLYGON ){
	if( (nobj = WlzPolygonToObj(obj, WLZ_SIMPLE_FILL,
				    NULL)) != NULL ){
	  (void )WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
      }
      else {
	(void )WlzWriteObj(stdout, obj);
      }
      break;
	
    default:
      (void )WlzWriteObj(stdout, obj);
      break;
    }

    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
