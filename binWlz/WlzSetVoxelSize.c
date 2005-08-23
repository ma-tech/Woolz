#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzSetVoxelSize.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Resets the voxel sizes of the input 3D object.
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
	  "Usage:\t%s [-x#] [-y#] [-z#] [-h] [-v] [<input file>]\n"
	  "\tReset the voxel sizes of the input 3D object\n"
	  "\twriting the new object to standard output.\n"
	  "\tThis is required until the 3D objects are converted to\n"
	  "\tWLZ_TRANS_OBJ type. If a voxel size is not defined then\n"
	  "\tthe original size is retained\n"
	  "\tOptions are:\n"
	  "\t  -x#       x voxel size\n"
	  "\t  -y#       y voxel size\n"
	  "\t  -z#       z voxel size\n"
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
  FILE		*inFile;
  char 		optList[] = "x:y:z:hv";
  int		option;
  int		xFlg=0, yFlg=0, zFlg=0;
  float		x_size=1.0, y_size=1.0, z_size=1.0;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'x':
      x_size = atof(optarg);
      xFlg = 1;
      break;

    case 'y':
      y_size = atof(optarg);
      yFlg = 1;
      break;

    case 'z':
      z_size = atof(optarg);
      zFlg = 1;
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

  if( (x_size <= 0.0) || (y_size <= 0.0) || (z_size <= 0.0) ){
    fprintf(stderr, "%s: voxel sizes must be non-zero and positive\n",
	    argv[0]);
    return 1;
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
    switch( obj->type )
    {
    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.p && (obj->domain.p->type == WLZ_2D_DOMAINOBJ) ){
	if( xFlg ){
	  obj->domain.p->voxel_size[0] = x_size;
	}
	if( yFlg ){
	  obj->domain.p->voxel_size[1] = y_size;
	}
	if( zFlg ){
	  obj->domain.p->voxel_size[2] = z_size;
	}
      }
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(1);
      }
      break;

    default:
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(1);
      }
      break;
    }
    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
