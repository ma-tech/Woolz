#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyRange.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Prints grey-range of a grey-level object.
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
	  "Usage:\t%s [-h] [-v] [<input file>]\n"
	  "\tPrint grey-range of a grey-level woolz object\n"
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

  WlzObject	*obj;
  WlzPixelV	gmin, gmax;
  FILE		*inFile;
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		objCount=0;
    
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

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects and threshold if possible */
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    objCount++;
    fprintf(stdout, "Grey ranges for input objects:\n");
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      /* get the existing min and max grey values */
      if( (errNum = WlzGreyRange(obj, &gmin, &gmax)) == WLZ_ERR_NONE ){
	switch( gmin.type ){
	case WLZ_GREY_INT:
	  fprintf(stdout,
		  "Obj %d: type: %s, gtype: WLZ_GREY_INT, "
		  "range: (%d, %d)\n", objCount,
		  WlzStringFromObjType(obj, NULL),
		  gmin.v.inv, gmax.v.inv);
	  break;
	case WLZ_GREY_SHORT:
	  fprintf(stdout,
		  "Obj %d: type: %s, gtype: WLZ_GREY_SHORT, "
		  "range: (%d, %d)\n", objCount,
		  WlzStringFromObjType(obj, NULL),
		  gmin.v.shv, gmax.v.shv);
	  break;
	case WLZ_GREY_UBYTE:
	  fprintf(stdout,
		  "Obj %d: type: %s, gtype: WLZ_GREY_UBYTE, "
		  "range: (%d, %d)\n", objCount,
		  WlzStringFromObjType(obj, NULL),
		  gmin.v.ubv, gmax.v.ubv);
	  break;
	case WLZ_GREY_FLOAT:
	  fprintf(stdout,
		  "Obj %d: type: %s, gtype: WLZ_GREY_FLOAT, "
		  "range: (%f, %f)\n", objCount,
		  WlzStringFromObjType(obj, NULL),
		  gmin.v.flv, gmax.v.flv);
	  break;
	case WLZ_GREY_DOUBLE:
	  fprintf(stdout,
		  "Obj %d: type: %s, gtype: WLZ_GREY_DOUBLE, "
		  "range: (%f, %f)\n", objCount,
		  WlzStringFromObjType(obj, NULL),
		  gmin.v.dbv, gmax.v.dbv);
	  break;
	}
      }
      break;

    default:
      fprintf(stdout,
	      "Obj %d: type: %s\n", objCount,
	      WlzStringFromObjType(obj, NULL));
      break;
    }
    WlzFreeObj(obj);
  }

  return 0;
}
