#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzConvertPix.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Converts the pixel type of a grey-level woolz object.
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
	  "Usage:\t%s [-t#] [-h] [<input file>]\n"
	  "\tConvert the pixel type of a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -t#       Output pixel type:\n"
	  "\t            # = %d: integer\n"
	  "\t                %d: short\n"
	  "\t                %d: unsigned byte (default)\n"
	  "\t                %d: float\n"
	  "\t                %d: double\n"
	  "\t                %d: rgba\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str, WLZ_GREY_INT, WLZ_GREY_SHORT, WLZ_GREY_UBYTE,
	  WLZ_GREY_FLOAT, WLZ_GREY_DOUBLE, WLZ_GREY_RGBA);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *newobj;
  FILE		*inFile;
  char 		optList[] = "ht:";
  int		option;
  WlzGreyType	newpixtype = WLZ_GREY_UBYTE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 't':
      switch( newpixtype = (WlzGreyType) atoi(optarg) ){

      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_UBYTE:
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
      case WLZ_GREY_RGBA:
	break;

      default:
        fprintf(stderr, "%s: grey type = %d is invalid\n",
		argv[0],newpixtype );
        usage(argv[0]);
        return 1;

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
      return 1;
    }
  }

  /* read objects and convert if possible */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    if( obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ ){
      if( (newobj = WlzConvertPix(obj, newpixtype, NULL)) != NULL ){
	(void )WlzWriteObj(stdout, newobj);
	WlzFreeObj(newobj);
      }
    }
    else {
      (void )WlzWriteObj(stdout, obj);
    }

    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
		

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
