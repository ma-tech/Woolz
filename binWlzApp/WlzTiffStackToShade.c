#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzTiffStackToShade.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Thu Jun 10 15:11:22 2004
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
* \ingroup      binWlzApp
* \brief        Convert a tiff stack to a shade image.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdio.h>
#include <stdlib.h>

#include <Wlz.h>
#include <WlzExtFF.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s <input file>\n"
	  "\tConvert a tiff stack to a shade image\n"
	  "\tby finding the maximal pixel values in each\n"
	  "\tchannel. The tiff is read from the given file,\n"
	  "\toutput to stdout.\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str);
  return;
}

int main(int	argc,
	 char	**argv)
{
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*tiffFile=NULL;
  WlzObject	*inObj, *outObj, *obj, **objVec;
  int		i, objVecCount;
  WlzGreyType	gType;

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
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  if( optind < argc ){
    tiffFile = *(argv+optind);
  }
  else {
    usage(argv[0]);
    return WLZ_ERR_UNSPECIFIED;
  }

  /* read the tiff file */
  if( (inObj = WlzAssignObject(WlzEffReadObj(NULL, tiffFile, WLZEFF_FORMAT_TIFF, 0,
					     &errNum), NULL)) == NULL ){
    usage(argv[0]);
    return errNum;
  }

  /* if 2D then that is the shade file, else take max of 3D stack */
  switch( inObj->type ){
  case WLZ_2D_DOMAINOBJ:
    outObj = inObj;
    break;

  case WLZ_3D_DOMAINOBJ:
    if( (errNum = WlzExplode3D(&objVecCount, &objVec, inObj)) != WLZ_ERR_NONE ){
      usage(argv[0]);
      return errNum;
    }
    gType = WlzGreyTypeFromObj(objVec[0], &errNum);
    outObj = WlzAssignObject(objVec[0], NULL);
    for(i=1; i < objVecCount; i++){
      if( gType == WLZ_GREY_RGBA ){
	obj = WlzRGBAImageArithmetic(outObj, objVec[i], WLZ_BO_MAX, 0, &errNum);
      }
      else {
	obj = WlzImageArithmetic(outObj, objVec[i], WLZ_BO_MAX, 0, &errNum);
      }
      if( obj ){
	WlzFreeObj(outObj);
	outObj = WlzAssignObject(obj, &errNum);
      }
      else {
	break;
      }
    }
    
    WlzFreeObj(inObj);
    break;

  default:
    WlzFreeObj(inObj);
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  /* write shade object */
  if( errNum == WLZ_ERR_NONE ){
    WlzWriteObj(stdout, outObj);
    WlzFreeObj(outObj);
  }
  else {
    usage(argv[0]);
  }
  return errNum;
}
