#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzRegisterICP.c
* Date:         December 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Attempts to register two objects using an itterative
*		closest point algorithm.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		idx,
		option,
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzTransformType trType = WLZ_TRANSFORM_2D_REG;
  WlzDomain	outDom;
  WlzValues	nullVal;
  WlzObject	*outObj = NULL;
  WlzObject	*inObj[2];
  FILE		*fP = NULL;
  char 		*outObjFileStr;
  char  	*inObjFileStr[2];
  const char	*errMsg;
  static char	optList[] = "o:hrt",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outDom.core = NULL;
  nullVal.core = NULL;
  inObj[0] = NULL;
  inObj[1] = NULL;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr[0] = inObjFileStrDef;
  inObjFileStr[1] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'r':
        trType = WLZ_TRANSFORM_2D_REG;
	break;
      case 't':
        trType = WLZ_TRANSFORM_2D_TRANS;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    idx = 0;
    while((idx < 2) && (optind < argc))
    {
      inObjFileStr[idx] = *(argv + optind);
      ++optind;
      ++idx;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  if(ok)
  {
    /* Read objects. */
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx]= WlzAssignObject(WlzReadObj(fP,
	  					   &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to read object %d from file %s (%s)\n",
		       *argv, idx, inObjFileStr[idx], errMsg);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
      }
      ++idx;
    }
  }
  if(ok)
  {
    /* Check object types. */
    if(inObj[0]->type != inObj[1]->type)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((inObj[0]->domain.core == NULL) ||
            (inObj[1]->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      switch(inObj[0]->type)
      {
        case WLZ_2D_POLYGON: /* FALLTHROUGH */
        case WLZ_BOUNDLIST:
	  break;
	case WLZ_CONTOUR:
	  if((inObj[0]->domain.ctr->model == NULL) ||
	     (inObj[1]->domain.ctr->model == NULL))
	  {
	    errNum = WLZ_ERR_DOMAIN_NULL;
	  }
	  else if(inObj[0]->domain.ctr->model->type !=
	          inObj[1]->domain.ctr->model->type)
	  {
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	  }
	  else
	  {
	    switch(inObj[0]->domain.ctr->model->type)
	    {
	      case WLZ_GMMOD_2I:
	      case WLZ_GMMOD_2D:
	        break;
	      case WLZ_GMMOD_3I:
	      case WLZ_GMMOD_3D:
		/* Make sure the type of affine transform is appropriate. */
		switch(trType)
		{
		  case WLZ_TRANSFORM_2D_REG:
		    trType = WLZ_TRANSFORM_3D_REG;
		    break;
		  case WLZ_TRANSFORM_2D_TRANS:
		    trType = WLZ_TRANSFORM_3D_TRANS;
		    break;
		}
	        break;
	      default:
		errNum = WLZ_ERR_DOMAIN_TYPE;
	        break;
	    }
	  }
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: input object(s) not appropriate\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    outDom.t = WlzRegICPObjs(inObj[0], inObj[1], trType,
			     NULL, NULL, 1000, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to register objects (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_AFFINE_TRANS, outDom, nullVal, NULL, NULL,
    			 &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to make affine transform object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?  fopen(outObjFileStr, "w"):
	      				    stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(inObj[0])
  {
    (void )WlzFreeObj(inObj[0]);
  }
  if(inObj[1])
  {
    (void )WlzFreeObj(inObj[1]);
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(outDom.core)
  {
    (void )WlzFreeAffineTransform(outDom.t);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s%s",
    *argv,
    " [-o<out obj>] [-t] [-r] [<in obj 0>] [<in obj 1>]\n"
    "Options:\n"
    "  -o  Output file name for affine transform.\n"
    "  -t  Find the translation only transform.\n"
    "  -r  Find the rigid body (aka registration) transform, default.\n"
    "  -h  Help, prints this usage message.\n"
    "Attempts to register two objects using an itterative closest point\n"
    "(ICP) algorithm.  The two objects must be contours, boundary lists\n"
    "or polygons.\n"
    "The input objects are read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    "cat in0.wlz | \\\n             ",
    *argv,
    " -o out.wlz -a - in1.wlz\n"
    "An affine transform is found by registering in1.wlz to in0.wlz using\n"
    "an ICP algorithm. The affine transform is then written to out.wlz.\n"
    "obj2.wlz.\n");
  }
  return(!ok);
}
