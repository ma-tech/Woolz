#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzBoundaryToObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Filter to convert a boundary list to the corresponding
*		domain object.
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
	  "\tConvert a boundary list woolz object\n"
	  "\tto the corresponding domain object\n"
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
    case WLZ_BOUNDLIST:
      if( (nobj = WlzBoundaryToObj(obj, WLZ_EVEN_ODD_FILL,
      				NULL)) != NULL ){
	(void )WlzWriteObj(stdout, nobj);
	WlzFreeObj(nobj);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.p->type == WLZ_PLANEDOMAIN_BOUNDLIST ){
	if( (nobj = WlzBoundaryToObj(obj, WLZ_EVEN_ODD_FILL,
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
