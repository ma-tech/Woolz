#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzCrossCorValue.c
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
* \brief	Computes the cross correlation value (with no shift) of
*		a pair of 2D spatial domain objects with grey values.
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
		usage = 0,
  		unionFlg = 0;
  double 	cCor = 0.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzTransformType trType = WLZ_TRANSFORM_2D_REG;
  WlzDomain	outDom;
  WlzValues	nullVal;
  WlzObject	*inObj[2];
  FILE		*fP = NULL;
  char 		*outFileStr;
  char  	*inObjFileStr[2];
  const char	*errMsg;
  static char	optList[] = "ho:u",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outDom.core = NULL;
  nullVal.core = NULL;
  inObj[0] = NULL;
  inObj[1] = NULL;
  outFileStr = outFileStrDef;
  inObjFileStr[0] = inObjFileStrDef;
  inObjFileStr[1] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'u':
        unionFlg = 1;
	break;
      case 'h':
      default:
        usage = 1;
	break;
    }
  }
  ok = !usage;
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outFileStr == NULL) || (*outFileStr == '\0'))
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
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP,
	  					    &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to read object %d from file %s (%s)\n",
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
    if((inObj[0]->type != WLZ_2D_DOMAINOBJ) ||
       (inObj[0]->type != inObj[1]->type))
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((inObj[0]->domain.core == NULL) ||
            (inObj[1]->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else if((inObj[0]->values.core == NULL) ||
            (inObj[1]->values.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: Input object(s) not appropriate\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    cCor = WlzCCorS2D(inObj[0], inObj[1], unionFlg, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Compute cross correlation of objects (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    else
    {
      (void )fprintf(fP, "%g\n", cCor);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-o<out obj>] [-u] [<in obj 0>] [<in obj 1>]\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file name for affine transform.\n"
    "  -u  Compute value in union of objects domains. The default is to\n"
    "      compute the cross correlation value within the intersection of\n"
    "      the objects domains.\n"
    "Computes the cross correlation value (with zero shift) of two\n"
    "2D spatial domain objects with grey values.\n"
    "The input objects are read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    *argv,
    " -o out.num in0.wlz in1.wlz\n"
    "The cross correlation value computed from the objects read from\n"
    "in0.wlz and in1.wlz is written to the file out.num.\n");
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
