#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzShadeCorrect.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Jun 16 17:29:29 2004
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
* \ingroup      binWlz
* \brief        Apply shade correction to input objects.
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
	  "Usage:\t%s -b <bright-field image> -d <dark-field image>"
	  "-h -v [<input file>]\n"
	  "\tApply a shade correction to the input objects. The shade\n"
	  "\tcorrection using the bright-field and dark-field images\n"
	  "\tassumes intensity images and calculates a normalised\n"
	  "\ttransmission coefficient.\n"
	  "\tOptions are:\n" 
	  "\t  -b <file> Bright-field image - REQUIRED\n"
	  "\t  -d <file> Dark-field image - optional\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str);
  return;
}

int main(int	argc,
	 char	**argv)
{
  char 		optList[] = "b:d:hv";
  int		option;
  int		verboseFlg=0;
  FILE		*fp;
  WlzObject	*inObj=NULL, *outObj, *bfObj=NULL, *dfObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      if( (fp = fopen(optarg, "rb")) == NULL ){
	fprintf(stderr, "%s: can't open bright-field file %s\n",
		argv[0], argv[optind]);
	usage(argv[0]);
	return WLZ_ERR_UNSPECIFIED;
      }
      else {
	bfObj = WlzReadObj(fp, &errNum);
	fclose(fp);
      }
      break;

    case 'd':
      if( (fp = fopen(optarg, "rb")) == NULL ){
	fprintf(stderr, "%s: can't open dark-field file %s\n",
		argv[0], argv[optind]);
	usage(argv[0]);
	return WLZ_ERR_UNSPECIFIED;
      }
      else {
	dfObj = WlzReadObj(fp, &errNum);
	fclose(fp);
	break;
      }
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  /* bright field must have been defined */
  if( bfObj == NULL ){
    fprintf(stderr, "%s: bright field image required\n",
	    argv[0]);
    usage(argv[0]);
    return WLZ_ERR_UNSPECIFIED;
  }

  /* get the input stream */
  if( optind < argc ){
    if( (fp = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open input file %s\n",
	      argv[0], argv[optind]);
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }
  else {
    fp = stdin;
  }

  /* now read the objects and shade correct */
  while( inObj = WlzReadObj(fp, &errNum) ){
    switch( inObj->type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      if( outObj = WlzShadeCorrectBFDF(inObj, bfObj, dfObj,
				       255.0, 0, &errNum) ){
	WlzWriteObj(stdout, outObj);
	WlzFreeObj(outObj);
      }
      else {
	fprintf(stderr, "%s: shade correction failed\n", argv[0]);
	return errNum;
      }
      break;

    default:
      WlzWriteObj(stdout, inObj);
      break;
    }
    WlzFreeObj(inObj);
  }

  if( bfObj ){
    WlzFreeObj(bfObj);
  }
  if( dfObj ){
    WlzFreeObj(dfObj);
  }

  return 0;
}
