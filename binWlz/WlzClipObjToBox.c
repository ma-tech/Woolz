#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzClipObjToBox.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter for clipping either a 2D or 3D domain
*		objects domain so that it lies within a given box.
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
  WlzIBox3	clipBox,
		clipSet;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  int		clipVal[2];
  char		*clipStr[2];
  static char	optList[] = "o:x:y:z:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  clipSet.xMin = 0;
  clipSet.xMax = 0;
  clipSet.yMin = 0;
  clipSet.yMax = 0;
  clipSet.zMin = 0;
  clipSet.zMax = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'x':
      case 'y':
      case 'z':
	if(optarg)
	{
	  while(*optarg && isspace(*optarg))
	  {
	    ++optarg;
	  }
	  if(*optarg == ',')
	  {
	    clipStr[0] = NULL;
	    clipStr[1] = strtok(optarg, ",");
	  }
	  else
	  {
	    clipStr[0] = strtok(optarg, ",");
	    clipStr[1] = strtok(NULL, ",");
	  }
	  if((clipStr[0] == NULL) && (clipStr[1] == NULL))
	  {
	    usage = 1;
	    ok = 0;
	  }
	  else
	  {
	    idx = 0;
	    while(ok && (idx < 2))
	    {
	      if(clipStr[idx] && (sscanf(clipStr[idx], "%d",
					 clipVal + idx) != 1))
	      {
		usage = 1;
		ok = 0;
	      }
	      ++idx;
	    }
	  }
	}
	if(ok)
	{
	  switch(option)
	  {
	    case 'x':
	      if(clipStr[0])
	      {
		clipSet.xMin = 1;
		clipBox.xMin = clipVal[0];
	      }
	      if(clipStr[1])
	      {
		clipSet.xMax = 1;
		clipBox.xMax = clipVal[1];
	      }
	      break;
	    case 'y':
	      if(clipStr[0])
	      {
		clipSet.yMin = 1;
		clipBox.yMin = clipVal[0];
	      }
	      if(clipStr[1])
	      {
		clipSet.yMax = 1;
		clipBox.yMax = clipVal[1];
	      }
	      break;
	    case 'z':
	      if(clipStr[0])
	      {
		clipSet.zMin = 1;
		clipBox.zMin = clipVal[0];
	      }
	      if(clipStr[1])
	      {
		clipSet.zMax = 1;
		clipBox.zMax = clipVal[1];
	      }
	      break;
	  }
	}
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
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzReadObj(fP, &errNum)) == NULL) ||
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
  if(ok)
  {
    switch(inObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	if(clipSet.xMin == 0)
	{
          clipBox.xMin = inObj->domain.i->kol1;
	}
	if(clipSet.xMax == 0)
	{
          clipBox.xMax = inObj->domain.i->lastkl;
	}
	if(clipSet.yMin == 0)
	{
          clipBox.yMin = inObj->domain.i->line1;
	}
	if(clipSet.yMax == 0)
	{
          clipBox.yMax = inObj->domain.i->lastln;
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	if(clipSet.xMin == 0)
	{
          clipBox.xMin = inObj->domain.p->kol1;
	}
	if(clipSet.xMax == 0)
	{
          clipBox.xMax = inObj->domain.p->lastkl;
	}
	if(clipSet.yMin == 0)
	{
          clipBox.yMin = inObj->domain.p->line1;
	}
	if(clipSet.yMax == 0)
	{
          clipBox.yMax = inObj->domain.p->lastln;
	}
	if(clipSet.zMin == 0)
	{
          clipBox.zMin = inObj->domain.p->plane1;
	}
	if(clipSet.zMax == 0)
	{
          clipBox.zMax = inObj->domain.p->lastpl;
	}
        break;
      default:
        ok = 0;
	(void )fprintf(stderr,
		       "%s: invalid object type %d\n", 
		       *argv, inObj->type);
    }
  }
  if(ok)
  {
    if(inObj->domain.core == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: invalid object, null domain\n",
		     *argv);
    }
  }
  if(ok)
  {
    if(((outObj = WlzClipObjToBox3D(inObj, clipBox, &errNum)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to clip object\n",
		     *argv);
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to write clipped object\n",
		     *argv);
    }
    if(fP && strcmp(outObjFileStr, "-"))
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
    " [-h] [-o<out object>]\n"
    "       [-x<x min>,<x max>] [-y<y min>,<y max>] [-z<z min>,<z max>]\n"
    "       [<in object>]\n"
    "Clips the input Woolz object using specified clip box.\n"
    "Options:\n"
    "  -o  Output object file name.\n"
    "  -x  Column clip box limits.\n"
    "  -y  Line clip box limits.\n"
    "  -z  Plane clip box limits.\n"
    "  -h  Help, prints this usage message.\n"
    "All clipping limits default to the input objects domain limits.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    "-o clipped.wlz -x 100,800 -y ,800 myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, clipped and written\n"
    "to clipped.wlz.\n");
  }
  return(!ok);
}
