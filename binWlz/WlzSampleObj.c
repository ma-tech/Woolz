#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzSampleObj.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which applies a sub-sampling filter to
*		the given object.
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
		usage = 0;
  WlzIVertex2	samFac;
  WlzSampleFn	samFn;
  char		*outObjFileStr,
  		*inObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  const char	*errMsg;
  static char	optList[] = "o:x:y:aegipmh",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  samFac.vtX = 1;
  samFac.vtY = 1;
  samFn = WLZ_SAMPLEFN_POINT;
  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outObjFileStr = outObjFileStrDef;
  while(ok && (usage == 0) &&
  	((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
	samFn = WLZ_SAMPLEFN_MAX;
	break;
      case 'e':
	samFn = WLZ_SAMPLEFN_MEDIAN;
	break;
      case 'i':
	samFn = WLZ_SAMPLEFN_MIN;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'g':
	samFn = WLZ_SAMPLEFN_GAUSS;
	break;
      case 'm':
	samFn = WLZ_SAMPLEFN_MEAN;
	break;
      case 'p':
	samFn = WLZ_SAMPLEFN_POINT;
	break;
      case 'x':
        if(sscanf(optarg, "%d", &(samFac.vtX)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'y':
        if(sscanf(optarg, "%d", &(samFac.vtY)) != 1)
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
    if((outObj = WlzSampleObj(inObj, samFac, samFn, &errNum)) == NULL)
    {
    ok = 0;
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
    		   "%s: failed to sub-sample object (%s).\n",
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
                     "%s: failed to write output object (%s).\n",
                     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    fprintf(stderr,
	    "Usage: "
	    "%s [-x#] [-y#] [-k#] [-l#] [-g] [-m] [-p] [-h] [<in object>]\n"
	    "Sub-samples a woolz interval domain object with grey values.\n"
	    "Options are:\n"
	    "  -x#   Sampling factor for rows (set to %d).\n"
	    "  -y#   Sampling factor for columns (set to %d).\n"
	    "  -o#   Output object file name.\n"
	    "  -a    Use max sampling kernel (%sset).\n"
	    "  -e    Use median sampling kernel (%sset).\n"
	    "  -i    Use min sampling kernel (%sset).\n"
	    "  -g    Use gaussian sampling kernel (%sset).\n"
	    "  -m    Use mean sampling kernel (%sset).\n"
	    "  -p    Use point sampling (%sset).\n"
	    "  -h    Display this usage information.\n",
	    *argv,
	    samFac.vtX, samFac.vtY,
	    (samFn == WLZ_SAMPLEFN_MAX) ? "" : "not ",
	    (samFn == WLZ_SAMPLEFN_MEDIAN) ? "" : "not ",
	    (samFn == WLZ_SAMPLEFN_MIN) ? "" : "not ",
	    (samFn == WLZ_SAMPLEFN_GAUSS) ? "" : "not ",
	    (samFn == WLZ_SAMPLEFN_MEAN) ? "" : "not ",
	    (samFn == WLZ_SAMPLEFN_POINT) ? "" : "not ");
  }
  return(!ok);
}
