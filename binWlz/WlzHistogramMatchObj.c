#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzHistogramMatchObj.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which modifies the grey values of the
*		given Woolz domain object so that the histogram of
*		the modified object's grey values matches the given
*		histogram object.
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
		dither = 0,
		smoothingSrcHist = 0,
		smoothingTargetHist = 0,
		independentFlag = 0,
		ok = 1,
		usage = 0;
  double	distVal[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL,
		*matchObj = NULL;
  char 		*outObjFileStr,
		*matchObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  char		*distStr[2];
  static char	optList[] = "Dm:s:t:o:id:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  matchObjFileStr = NULL;
  outObjFileStr = outObjFileStrDef;
  distVal[0] = 0.0;
  distVal[1] = 1.0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'D':
        dither = 1;
	break;
      case 'm':
        matchObjFileStr = optarg;
	break;
      case 's':
        if(sscanf(optarg, "%d", &smoothingSrcHist) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 't':
        if(sscanf(optarg, "%d", &smoothingTargetHist) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'i':
      	independentFlag = 1;
      case 'd':
        if(optarg)
	{
	  while(*optarg && isspace(*optarg))
	  {
	    ++optarg;
	  }
	  if(*optarg == ',')
	  {
	    distStr[0] = NULL;
	    distStr[1] = strtok(optarg, ",");
	  }
	  else
	  {
	    distStr[0] = strtok(optarg, ",");
	    distStr[1] = strtok(NULL, ",");
	  }
	  if((distStr[0] == NULL) && (distStr[1] == NULL))
	  {
	    usage = 1;
	    ok = 0;
	  }
	  else
	  {
	    idx = 0;
	    while(ok && (idx < 2))
	    {
	      if(distStr[idx] && (sscanf(distStr[idx], "%lg",
	      				 distVal + idx) != 1))
	      {
	        usage = 1;
		ok = 0;
	      }
	      ++idx;
	    }
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
     (matchObjFileStr == NULL) || (*matchObjFileStr == '\0') ||
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
    if(((fP = (strcmp(matchObjFileStr, "-")?
        fopen(matchObjFileStr, "r"): stdin)) == NULL) ||
       ((matchObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		    "%s: failed to read histogram object from file %s (%s).\n",
		    *argv, matchObjFileStr, errMsg);
    }
    if(fP && strcmp(matchObjFileStr, "-"))
    {
      fclose(fP);
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
    if((inObj->type != WLZ_2D_DOMAINOBJ) && (inObj->type != WLZ_3D_DOMAINOBJ))
    {
      ok = 0;
      (void )fprintf(stderr,
		 "%s: input object read from file %s is not a domain object\n",
		     *argv, inObjFileStr);
    }
  }
  if(ok)
  {
    if(matchObj->type != WLZ_HISTOGRAM)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: input object read from file %s is not a histogram\n",
		     *argv, inObjFileStr);
    }
  }
  if(ok)
  {
    if((errNum = WlzHistogramSmooth(matchObj,
    				    smoothingTargetHist)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to smooth target histogram (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    outObj = WlzAssignObject(inObj, NULL);
    if((errNum = WlzHistogramMatchObj(outObj, matchObj,
    				      independentFlag, smoothingSrcHist,
    				      distVal[0], distVal[1],
				      dither)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to match object to histogram %s).\n",
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
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(matchObj)
  {
    WlzFreeObj(matchObj);
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
    " -m <histogram object> [-h] [-o<output object>]\n"
    "                              [-D] [-i] [-s#] [-t#] [-d#,#]\n"
    "                              [<input object>]\n"
    "Options:\n"
    "  -D  Dither pixel values when mapping.\n"
    "  -m  Histogram object to match domain objects values to.\n"
    "  -s  Input object histogram smoothing.\n"
    "  -t  Target histogram smoothing.\n"
    "  -o  Output object/data file name.\n"
    "  -i  If the domain object is a 3D objects then the planes are matched\n"
    "      independently of each other.\n"
    "  -d  Histogram distance range for matching.\n"
    "  -h  Help, prints this usage message.\n"
    "Modifies the grey values of the Woolz domain object so that the\n"
    "histogram of the modified object's grey values matches the given\n"
    "histogram object.\n"
    "If a histogram distance range is given then only domain objects/planes\n"
    "within the range are matched. The range is specified as <min>,<max>\n"
    "where both the minimum and maximum are optional. A valid histogram\n"
    "range is contained by the closed interval [0.0-1.0].\n"
    "Smoothing parameters are gaussian convolution half height full widths\n"
    "specified in histogram bins.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " -m myhist.wlz -i -o matched.wlz mydom.wlz\n"
    "The input Woolz domain object is read from mydom.wlz, matched\n"
    "to the histogram read from the file myhist.wlz and written to the\n",
    "file matched.wlz. If the domain object is a 3D object then it's\n",
    "planes are matched independently.\n");
  }
  return(!ok);
}
