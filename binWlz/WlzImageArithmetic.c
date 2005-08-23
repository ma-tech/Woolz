#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzImageArithmetic.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Performs binary (ie two images) image arithmetic on a
*		pair of Woolz domain objects to produce a new object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

typedef enum
{
  WLZ_IMARNORM_NONE,
  WLZ_IMARNORM_INPUT,
  WLZ_IMARNORM_256
} WlzImArNorm;

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
  WlzImArNorm	norm = WLZ_IMARNORM_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*outObj = NULL;
  WlzObject	*inObj[2];
  WlzBinaryOperatorType operator = WLZ_BO_ADD;
  char 		*outObjFileStr;
  char  	*inObjFileStr[2];
  WlzPixelV	gMin[3],
  		gMax[3];
  const char	*errMsg;
  static char	optList[] = "o:adglmsnNh",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
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
      case 'a':
        operator = WLZ_BO_ADD;
	break;
      case 'd':
        operator = WLZ_BO_DIVIDE;
	break;
      case 'g':
      	operator = WLZ_BO_MAGNITUDE;
	break;
      case 'l':
        operator = WLZ_BO_MODULUS;
	break;
      case 'm':
        operator = WLZ_BO_MULTIPLY;
	break;
      case 's':
        operator = WLZ_BO_SUBTRACT;
	break;
      case 'n':
        norm = WLZ_IMARNORM_INPUT;
	break;
      case 'N':
        norm = WLZ_IMARNORM_256;
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
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      if((inObj[idx]->type != WLZ_2D_DOMAINOBJ) &&
         (inObj[idx]->type != WLZ_3D_DOMAINOBJ))
      {
        ok = 0;
	(void )fprintf(stderr, "%s: input object %d is not a domain object\n",
		       *argv, idx);
      }
      ++idx;
    }
  }
  if(ok && (norm == WLZ_IMARNORM_INPUT))
  {
    if(((errNum = WlzGreyRange(inObj[0], gMin, gMax)) != WLZ_ERR_NONE) ||
       ((errNum = WlzGreyRange(inObj[1], gMin + 1, gMax + 1)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to find input objects grey range (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if((outObj = WlzImageArithmetic(inObj[0], inObj[1], operator, 1,
				    &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute new object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok && (norm != WLZ_IMARNORM_NONE))
  {
    if((errNum = WlzGreyRange(outObj, gMin+ 2, gMax + 2)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to find output object's grey range (%s)\n",
		     *argv, errMsg);
    }
    if(ok)
    {
      if(norm == WLZ_IMARNORM_256)
      {
        gMin[0].type = WLZ_GREY_UBYTE;
	gMin[0].v.ubv = 0;
        gMax[0].type = WLZ_GREY_UBYTE;
	gMax[0].v.ubv = 255;
	(void )WlzValueConvertPixel(gMin, gMin[0], gMin[2].type);
	(void )WlzValueConvertPixel(gMax, gMax[0], gMax[2].type);
      }
      else
      {
	(void )WlzValueConvertPixel(gMin, gMin[0], gMin[2].type);
	(void )WlzValueConvertPixel(gMax, gMax[0], gMax[2].type);
	(void )WlzValueConvertPixel(gMin + 1, gMin[1], gMin[2].type);
	(void )WlzValueConvertPixel(gMax + 1, gMax[1], gMax[2].type);
	switch(gMin[2].type)
	{
	  case WLZ_GREY_INT:
	    if(gMin[1].v.inv < gMin[0].v.inv)
	    {
	      gMin[0].v.inv = gMin[1].v.inv;
	    }
	    if(gMax[1].v.inv > gMax[0].v.inv)
	    {
	      gMax[0].v.inv = gMax[1].v.inv;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    if(gMin[1].v.shv < gMin[0].v.shv)
	    {
	      gMin[0].v.shv = gMin[1].v.shv;
	    }
	    if(gMax[1].v.shv > gMax[0].v.shv)
	    {
	      gMax[0].v.shv = gMax[1].v.shv;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if(gMin[1].v.ubv < gMin[0].v.ubv)
	    {
	      gMin[0].v.ubv = gMin[1].v.ubv;
	    }
	    if(gMax[1].v.ubv > gMax[0].v.ubv)
	    {
	      gMax[0].v.ubv = gMax[1].v.ubv;
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    if(gMin[1].v.flv < gMin[0].v.flv)
	    {
	      gMin[0].v.flv = gMin[1].v.flv;
	    }
	    if(gMax[1].v.flv > gMax[0].v.flv)
	    {
	      gMax[0].v.flv = gMax[1].v.flv;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    if(gMin[1].v.dbv < gMin[0].v.dbv)
	    {
	      gMin[0].v.dbv = gMin[1].v.dbv;
	    }
	    if(gMax[1].v.dbv > gMax[0].v.dbv)
	    {
	      gMax[0].v.dbv = gMax[1].v.dbv;
	    }
	    break;
	}
      }
    }
    if((errNum = WlzGreySetRange(outObj, gMin[2], gMax[2],
    				 gMin[0], gMax[0])) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to set output object's grey range (%s)\n",
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
    WlzFreeObj(inObj[0]);
  }
  if(inObj[1])
  {
    WlzFreeObj(inObj[1]);
  }
  if(outObj)
  {
    WlzFreeObj(outObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out file>] [-a] [-d] [-l] [-m] [-s] [-h] "
    "[<in object 0>] [<in object 1>]\n"
    "Options:\n"
    "  -o  Output file name.\n"
    "  -a  Add the object's grey values.\n"
    "  -d  Divide the grey values of the 1st object by those of the 2nd.\n"
    "  -g  Vector magnitude of horizontal and vertical component objects\n"
    "  -l  Compute the modulus of the grey values of the 1st object wrt\n"
    "      those of the 2nd.\n"
    "  -m  Multiply the object's grey values.\n"
    "  -s  Subtract the grey values of the 2nd object from those of the 1st.\n"
    "  -n  Normalises the output object to the range of the imput objects.\n"
    "  -N  Normalises the output objects to the range [0-255].\n"
    "  -h  Help, prints this usage message.\n"
    "Computes an arithmetic binary (two objects) operation on two domain\n"
    "objects. The default operator is add.\n"
    "The input objects are read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    "cat obj1.wlz | %s",
    " -o obj3.wlz -a - obj2.wlz\n"
    "A new object 'obj3.wlz is formed by adding the grey values of obj1 and\n"
    "obj2.wlz.\n", *argv);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
