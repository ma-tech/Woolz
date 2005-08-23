#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGauss.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Applies a Gaussian filter to a grey level object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

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
	  "Usage:\t%s [-w#[#]] [-x#] [-y#] [-h] [<input file>]\n"
	  "\tApply a Gaussian filter to a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -w#[,#]   x_width[y_width] gaussian widths defined as\n"
	  "\t            full width half maximum in pixels\n"
	  "\t            default value 3.0, if y_width omitted\n"
	  "\t            then it is set equal to x_width\n"
	  "\t  -x#       x derivative - possible values 0,1,2, default - 0\n"
	  "\t  -y#       y derivative - possible values 0,1,2, default - 0\n"
	  "\t  -h        Help - prints this usage message\n",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "hw:x:y:";
  int		option;
  double	x_width, y_width;
  int		x_deriv, y_deriv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* set defaults, read the argument list and check for an input file */
  opterr = 0;
  x_width = 3.0;
  y_width = 3.0;
  x_deriv = 0;
  y_deriv = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'w':
      switch( sscanf(optarg, "%lg,%lg", &x_width, &y_width) ){

      default:
      case 0:
	fprintf(stderr, "%s: no gaussian width set\n", argv[0]);
        usage(argv[0]);
        return 1;

      case 1:
	y_width = x_width;
	break;

      case 2:
	break;

      }
      break;

    case 'x':
      x_deriv = atoi(optarg);
      break;

    case 'y':
      y_deriv = atoi(optarg);
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
      case WLZ_2D_DOMAINOBJ:
	if( (nobj = WlzGauss2(obj, x_width, y_width,
			      x_deriv, y_deriv, &errNum)) != NULL ){
	  errNum = WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
	break;

      case WLZ_3D_DOMAINOBJ:
      default:
	errNum = WlzWriteObj(stdout, obj);
	break;
    }

    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
