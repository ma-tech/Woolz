#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzDilation.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Erodes a domain woolz object.
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
	  "Usage:\t%s [-c#] [-h] [<input file>]\n"
	  "\tErode a domain woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -c#       Dilation connectivity:\n"
	  "\t            # =  4: 4-connected (2D)\n"
	  "\t                 8: 8-connected (2D) - default\n"
	  "\t                 6: 6-connected (3D)\n"
	  "\t                18: 18-connected (3D)\n"
	  "\t                26: 26-connected (3D)\n"
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
  char 		optList[] = "hc:";
  int		option;
  WlzConnectType	connectivity;
  int		conn;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* read the argument list and check for an input file */
  opterr = 0;
  connectivity = WLZ_8_CONNECTED;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'c':
      switch( conn = atoi(optarg) ){

      case 8:
	connectivity = WLZ_8_CONNECTED;
	break;
	
      case 4:
	connectivity = WLZ_4_CONNECTED;
	break;
	
      case 6:
	connectivity = WLZ_6_CONNECTED;
	break;
	
      case 18:
	connectivity = WLZ_18_CONNECTED;
	break;
	
      case 26:
	connectivity = WLZ_26_CONNECTED;
	break;
	
      default:
        fprintf(stderr, "%s: connectivity = %d is invalid\n",
		argv[0], conn );
        usage(argv[0]);
        return WLZ_ERR_PARAM_DATA;

      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  /* read objects and threshold if possible */
  while((obj = WlzAssignObject(WlzReadObj(inFile, &errNum), NULL)) != NULL) 
  {
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
	if( (nobj = WlzDilation(obj, connectivity, &errNum)) != NULL ){
	  errNum = WlzWriteObj(stdout, nobj);
	  (void) WlzFreeObj(nobj);
	}
	break;

      default:
	errNum = WlzWriteObj(stdout, obj);
	break;
    }

    WlzFreeObj(obj);
  }

  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if( errNum == WLZ_ERR_READ_EOF ){
    errNum = WLZ_ERR_NONE;
  }
  return errNum;
}
