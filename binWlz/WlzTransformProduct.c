#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzTransformProduct.c
* \author       Bill Hill
* \date         August 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Computes the product of a pair of transforms.
* \todo         -
* \bug          None known.
*/
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
  WlzAffineTransform *t[3];
  WlzValues	trVal;
  WlzDomain	trDom;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outObjFileStr = NULL;
  char		*inObjFileStr[2];
  static char	optList[] = "ho:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  t[0] = t[1] = t[2] = NULL;
  outObjFileStr = outFileStrDef;
  inObjFileStr[0] = inObjFileStr[1] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if(optind < argc)
    {
      if((optind + 1) == argc)
      {
	inObjFileStr[0] = *(argv + optind);
      }
      else if((optind + 2) == argc)
      {
	inObjFileStr[0] = *(argv + optind);
	inObjFileStr[1] = *(argv + optind + 1);
      }
      else
      {
	usage = 1;
	ok = 0;
      }
    }
  }
  if(ok)
  {
    idx = 0;
    while(ok && (idx < 2))
    {
      if((inObjFileStr[idx] == NULL) ||
	 (*inObjFileStr[idx] == '\0') ||
	 ((fP = (strcmp(inObjFileStr[idx], "-")?
		fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	 ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	 (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: Failed to read object from file %s.\n",
		       *argv, inObjFileStr[idx]);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
      }
      if(ok)
      {
	if((obj->type != WLZ_AFFINE_TRANS) ||
	   (obj->domain.core == NULL) ||
	   ((obj->domain.core->type != WLZ_TRANSFORM_2D_AFFINE) &&
	    (obj->domain.core->type != WLZ_TRANSFORM_3D_AFFINE)))
	{
	  ok = 0;
	  (void )fprintf(stderr,
			 "%s: Invalid object read from file %s.\n",
			 *argv, inObjFileStr[idx]);
	}
      }
      if(ok)
      {
	t[idx] = WlzAssignAffineTransform(obj->domain.t, NULL);
      }
      if(obj)
      {
	(void )WlzFreeObj(obj);
	obj = NULL;
      }
      ++idx;
    }
  }
  if(ok)
  {
    t[3] = WlzAffineTransformProduct(t[0], t[1], &errNum);
    if((t[3] == NULL) || (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Failed to compute product.\n", *argv);
    }
  }
  if(ok)
  {
    trDom.t = WlzAssignAffineTransform(t[3], NULL);
    trVal.core = NULL;
    obj = WlzMakeMain(WLZ_AFFINE_TRANS, trDom, trVal, NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to make transform object.\n",
		     *argv);
    }
    else
    {
      t[3] = NULL;
      if(((fP = (strcmp(outObjFileStr, "-")?
		fopen(outObjFileStr, "w"):
		stdout)) == NULL) ||
	  (WlzWriteObj(fP, obj) != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write output transform object\n",
		       *argv);
      }
      if(fP && strcmp(outObjFileStr, "-"))
      {
	fclose(fP);
      }
      (void )WlzFreeObj(obj);
      obj = NULL;
    }
  }
  for(idx = 0; idx < 3; ++idx)
  {
    (void )WlzFreeAffineTransform(t[idx]);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-o<output object>] [transform 0] [transform 1]\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "Computes the product of a pair of transforms.\n",
    *argv,
    " -o t2.wlz t0.wlz t1.wlz\n"
    "Computes the product of t0 and t1 such that t2 = t0 t1\n");
  }
  return(!ok);
}
