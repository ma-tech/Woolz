#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzShiftObj.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Shifts a Woolz object by applying an integral
*		translation.
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
		ok = 1,
		usage = 0,
		sftToOrg = 0,
		trX = 0,
		trY = 0,
		trZ = 0;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzIBox3	bBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inObjFileStr;
  static char	optList[] = "hgo:x:y:z:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'g':
        sftToOrg = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'x':
	if(sscanf(optarg, "%d", &trX) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'y':
	if(sscanf(optarg, "%d", &trY) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'z':
	if(sscanf(optarg, "%d", &trZ) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
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
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
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
      if((inObjFileStr == NULL) ||
	 (*inObjFileStr == '\0') ||
	 ((fP = (strcmp(inObjFileStr, "-")?
		fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	 ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	 (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read object from file %s\n",
		       *argv, inObjFileStr);
      }
      if(fP && strcmp(inObjFileStr, "-"))
      {
	fclose(fP);
      }
    }
  }
  if(ok)
  {
    if(sftToOrg)
    {
      bBox = WlzBoundingBox3I(inObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        trX = -(bBox.xMin);
        trY = -(bBox.yMin);
        trZ = -(bBox.zMin);
      }
      else
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to get object origin\n",
		       *argv);
      }
    }
  }
  if(ok)
  {
    outObj = WlzAssignObject(WlzShiftObject(inObj, trX, trY, trZ,
					    &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to shift object\n",
		     *argv);
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"):
	      stdout)) == NULL) ||
       (WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to write output object\n",
		     *argv);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
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
    " [-o<output object>] [-h] [-g] [-x#] [-y#] [-z#] [<input object>]\n"
    "Options:\n"
    "  -g  Shift object to the origin\n"
    "  -h  Prints this usage information\n"
    "  -o  Output object file name.\n"
    "  -x  Column (x) translation.\n"
    "  -y  Row (y) translation.\n"
    "  -z  Plane (z) translation.\n"
    "Shifts a Woolz object by applying an integer translation.\n"
    "The input object is read from stdin and the shifted object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -x100 -y200 -o shifted.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, shifted 100 columns\n"
    "and 200 lines and then written to shifted.wlz\n");
  }
  return(!ok);
}
