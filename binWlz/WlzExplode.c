#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExplode.c
* Date:         March 1999
* Author:       Bill Hill, Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter to explode 3D and compound objects into
*		2D domain objects.
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
  int		objPlane,
  		objVecIdx,
  		objVecCount,
  		option,
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  WlzObject	**objVec = NULL;
  WlzCompoundArray	*cmpObj;
  char 		*outObjFileBodyStr = NULL,
		*outObjFileExtStr = NULL,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "b:e:h",
		outObjFileStr[FILENAME_MAX],
		outObjFileExtStrDef[] = "wlz",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outObjFileExtStr = outObjFileExtStrDef;
  strcpy(outObjFileStr, outObjFileStrDef);
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        outObjFileBodyStr = optarg;
	break;
      case 'e':
        outObjFileExtStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outObjFileExtStr && (outObjFileBodyStr == NULL)))
  {
    ok = 0;
    usage = 1;
  }
  if(outObjFileBodyStr && outObjFileExtStr &&
     (strlen(outObjFileExtStr) +
      strlen(outObjFileExtStr)) > (FILENAME_MAX - 24))
  {
    ok = 0;
    (void )fprintf(stderr,
    		   "%s: output filename length exceeded!\n",
		   *argv);
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
       ((inObj= WlzReadObj(fP, &errNum)) == NULL))
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
    switch( inObj->type )
    {
    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      cmpObj = (WlzCompoundArray *) inObj;
      objVec = cmpObj->o;
      objVecCount = cmpObj->n;
      break;

    default:
      if((errNum = WlzExplode3D(&objVecCount, &objVec, inObj)) !=
	 WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to explode object (%s)\n",
		       *argv, errMsg);
      }
      break;
    }
  }
  if(ok)
  {
    objVecIdx = 0;
    objPlane = 0;
    if( inObj->type == WLZ_3D_DOMAINOBJ )
    {
      objPlane = inObj->domain.p->plane1;
    }
    while((objVecIdx < objVecCount) && ok)
    {
      if(outObjFileBodyStr)
      {
        sprintf(outObjFileStr,
		"%s%06d.%s", outObjFileBodyStr, objPlane, outObjFileExtStr);
        if((fP = fopen(outObjFileStr, "w")) == NULL)
	{
	  ok = 0;
	  (void )fprintf(stderr,
	  		 "%s: failed to open output file %s.\n",
			 *argv, outObjFileStr);
	}
      }
      else
      {
        fP = stdout;
      }
      if(ok)
      {
        if((errNum = WlzWriteObj(fP, *(objVec + objVecIdx))) != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	  		 "%s: failed to write output object (%s).\n",
			 *argv, errMsg);
	}
      }
      if(*(objVec + objVecIdx))
      {
        WlzFreeObj(*(objVec + objVecIdx));
      }
      if(outObjFileBodyStr && fP)
      {
        fclose(fP);
      }
      ++objPlane;
      ++objVecIdx;
    }
    AlcFree(objVec);
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-b<file body>] [-e <file extension>] [-h] [<in object>]\n"
    "Options:\n"
    "  -b  Output object file body.\n"
    "  -e  Output object file extension (default wlz).\n"
    "  -h  Help, prints this usage message.\n"
    "Explodes the input 3D Woolz domain object into 2D domain objects or\n"
    "explodes an input compound object into its constituent objects.\n"
    "An output file extension can only be given if an output file body is\n"
    "also given.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -b plane myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz (minimum plane coordinate\n"
    "30), exploded and written to plane000030.wlz, plane000031.wlz, etc.\n");
  }
  return(!ok);
}
