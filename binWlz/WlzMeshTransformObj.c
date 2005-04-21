#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMeshTransformObj.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Mar 11 13:51:31 2005
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
* \ingroup      WlzTransform
* \brief        Read in a mesh transform and warp given woolz object.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
	  "Usage:\t%s -m <mesh transform file>"
	  " [-o <output file>] [-h] [-v]"
	  " [<2D object input file>]\n"
	  "\tApply a mesh transform to given input objects\n"
	  "\twriting the warped objects to standard output.\n"
	  "\tOptions are:\n"
	  "\t  -L                 Use linear interpolation instead of nearest-neighbour\n"
	  "\t  -m<meshfile>       Mesh transform object\n"
	  "\t  -o<output file>    Output filename, default to stdout\n"
	  "\t  -h                 Help - this message\n"
	  "\t  -v                 verbose operation\n"
	  "",
	  proc_str);
  return;
}

int main(int	argc,
	 char	**argv)
{
  WlzObject	*inObj, *meshObj, *outObj;
  FILE		*inFP, *outFP, *meshFP;
  char		*outFile, *meshFile;
  char 		optList[] = "Lm:o:hv";
  int		option;
  int		verboseFlg = 0;
  const char    *errMsg;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* additional defaults */
  outFile = "-";
  meshFile = NULL;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'L':
      interp = WLZ_INTERPOLATION_LINEAR;
      break;

    case 'm':
      meshFile = optarg;
      break;

    case 'o':
      outFile = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 0;

    case 'v':
      verboseFlg = 1;
      break;

    }
  }

  /* check input file/stream */
  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* check output file/stream */
  if(strcmp(outFile, "-"))
  {
    if((outFP = fopen(outFile, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    outFP = stdout;
  }

  /* check meshfile */
  if((errNum == WLZ_ERR_NONE) && (meshFile != NULL)){
    if( meshFP = fopen(meshFile, "r") ){
      if( meshObj = WlzReadObj(meshFP, &errNum) ){
	if( meshObj->type != WLZ_MESH_TRANS ){
	  fprintf(stderr,
		  "%s: mesh file does not have a mesh transform object\n",
		  argv[0]);
	  return 1;
	}
      }
      else {
	fprintf(stderr, "%s: failed to read mesh file\n",
	    argv[0]);
	return 1;
      }
    }
    else {
      fprintf(stderr, "%s: failed to open mesh file\n",
	    argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: mesh input file required\n",
	    argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* transform any suitable input objects */
  while((errNum == WLZ_ERR_NONE) &&
        ((inObj = WlzReadObj(inFP, &errNum)) != NULL))
  {
    switch( inObj->type )
    {
    case WLZ_2D_DOMAINOBJ:
      if( outObj = WlzMeshTransformObj(inObj, meshObj->domain.mt,
				       interp, &errNum) ){
	WlzWriteObj(outFP, outObj);
	WlzFreeObj(outObj);
      }
      else {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to transform object (%s).\n",
		       *argv, errMsg);
	return 1;
      }
      break;

    default:
      WlzWriteObj(outFP, inObj);
      break;
    }

    WlzFreeObj(inObj);
  }

  if(errNum == WLZ_ERR_READ_EOF)
  {
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
}
