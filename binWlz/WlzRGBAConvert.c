#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzRGBAConvert.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Thu Jun  3 08:39:00 2004
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      binwlz
* \brief        Access to conversion functions for RGBA data.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

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
	  "Usage:\t%s [-c] [-m] [-h] [<input file>]\n"
	  "\tConvert the RGBA woolz object to a compound object\n"
	  "\tor to the modulus of the rgb values,\n"
	  "\twriting the new object to standard output.\n"
	  "\tNote input object MUST have grey-value type RGBA\n."
	  "\tOptions are:\n"
	  "\t  -c       convert to compound (default)\n"
	  "\t  -m       convert to modulus\n"
	  "\t  -s#	select colour space:\n"
	  "\t           # = %d: RGB (default)\n"
	  "\t               %d: HSB\n"
	  "\t               %d: CMY\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str, WLZ_RGBA_SPACE_RGB, WLZ_RGBA_SPACE_HSB,
	  WLZ_RGBA_SPACE_CMY);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *newobj;
  FILE		*inFile;
  char 		optList[] = "cmh";
  int		option;
  int		compoundFlg=1;
  WlzRGBAColorSpace	colSpc=WLZ_RGBA_SPACE_RGB;
  WlzErrorNum	errNum;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'c':
      compoundFlg = 1;
      break;

    case 'm':
      compoundFlg = 0;
      break;

    case 's':
      switch( colSpc = (WlzRGBAColorSpace) atoi(optarg) ){

      case WLZ_RGBA_SPACE_RGB:
      case WLZ_RGBA_SPACE_HSB:
      case WLZ_RGBA_SPACE_CMY:
	break;

      default:
        fprintf(stderr, "%s: colour space = %d is invalid\n",
		argv[0], colSpc );
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
      if( compoundFlg ){
	if( (newobj = (WlzObject *) 
	     WlzRGBAToCompound(obj, colSpc, &errNum)) != NULL ){
	  (void )WlzWriteObj(stdout, newobj);
	  WlzFreeObj(newobj);
	}
	else {
	  if( errNum == WLZ_ERR_VALUES_TYPE ){
	    fprintf(stderr, "%s: wrong values type, object not converted\n",
		    argv[0]);
	    (void )WlzWriteObj(stdout, obj);
	  }
	  else {
	    fprintf(stderr, "%s: something wrong, time to quit.\n",
		    argv[0]);
	    usage(argv[0]);
	    return 1;
	  }
	}
      }
      else {
	if( (newobj = WlzRGBAToModulus(obj, &errNum)) != NULL ){
	  (void )WlzWriteObj(stdout, newobj);
	  WlzFreeObj(newobj);
	}
	else {
	  if( errNum == WLZ_ERR_VALUES_TYPE ){
	    fprintf(stderr, "%s: wrong values type, object not converted\n",
		    argv[0]);
	    (void )WlzWriteObj(stdout, obj);
	  }
	  else {
	    fprintf(stderr, "%s: something wrong, time to quit.\n",
		    argv[0]);
	    usage(argv[0]);
	    return 1;
	  }
	}
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
