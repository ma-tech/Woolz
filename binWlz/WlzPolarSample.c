#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzPolarSample.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which applies a rectangular to polar
*		resampling filter to the given object.
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
  int		option,
		angleSamples,
		outFlag = 1,
		ok = 1,
		usage = 0;
  double	angleInc;
  WlzIVertex2	center,
  		objSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "o:ih",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'i':
        outFlag = 0;
	break;
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
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if(inObj->type == WLZ_EMPTY_OBJ)
    {
      center.vtX = 0;
      center.vtY = 0;
      angleInc = 1.0;
      angleSamples = 1;
    }
    else if(inObj->type == WLZ_2D_DOMAINOBJ)
    {
      if(inObj->domain.core == NULL)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: input object read from %s has null domain\n",
		       *argv, inObjFileStr);
      }
      else if((inObj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
      	      (inObj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: input object read from %s has invalid domain\n",
		       *argv, inObjFileStr);
      }
      else
      {
        center.vtX = (inObj->domain.i->kol1 + inObj->domain.i->lastkl) / 2;
        center.vtY = (inObj->domain.i->line1 + inObj->domain.i->lastln) / 2;
	objSz.vtX = inObj->domain.i->lastkl - inObj->domain.i->kol1 + 1;
	objSz.vtY = inObj->domain.i->lastln - inObj->domain.i->line1 + 1;
	if(outFlag)
	{
	  angleSamples = WLZ_MAX(objSz.vtX, objSz.vtY) / 2;
	}
	else
	{
	  angleSamples = WLZ_MIN(objSz.vtX, objSz.vtY) / 2;
	}
	if(angleSamples > 0)
	{
	  angleInc = (2.0 * WLZ_M_PI) / angleSamples;
	}
	else
	{
	  angleInc = 1.0;
	}
      }
    }
    else
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: input object read from %s has invalid type\n",
		     *argv, inObjFileStr);
    }
  }
  if(ok)
  {
    if((outObj = WlzPolarSample(inObj, center, angleInc, 1,
    				angleSamples, outFlag, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to polar resample object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output object (%s)\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out object>] [-i] [-h] [<in object>]\n"
    "Options:\n"
    "  -o  Output object file name.\n"
    "  -i  Resampling circle inside object, default is to use the outer\n"
    "      bounding circle.\n"
    "  -h  Help, prints this usage message.\n"
    "Polar resamples the input Woolz object.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o polar.wlz rect.wlz\n"
    "The input Woolz object is read from rect.wlz, resampled and written\n"
    "to polar.wlz.\n");
  }
  return(!ok);
}
