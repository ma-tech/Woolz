#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzHistogramObj.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which computes a Woolz histogram object
*		from the given Woolz domain object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 1999-09-19 bill	Error reporting could cause segfault, fixed bug.
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
		asciiFlag = 0,
		nBins = 0,
		ok = 1,
		usage = 0;
  double	binOrigin = 0.0,
  		binSize = 1.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzHistogramDomain *histDom;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  const char	*errMsg;
  char 		*outObjFileStr,
  		*inObjFileStr;
  static char	optList[] = "o:n:s:w:ah",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        asciiFlag = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'n':
	if(sscanf(optarg, "%d", &nBins) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
	if(sscanf(optarg, "%lg", &binOrigin) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'w':
	if(sscanf(optarg, "%lg", &binSize) != 1)
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
		     "%s: failed to read object from file %s (%s)\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if(((outObj = WlzHistogramObj(inObj, nBins, binOrigin,
    				  binSize, &errNum)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute histogram object (%s)\n",
		     *argv,
		     WlzStringFromErrorNum(errNum, NULL));
    }
  }
  if(ok)
  {
    if(asciiFlag)
    {
      if((fP = (strcmp(outObjFileStr, "-")?
               fopen(outObjFileStr, "w"): stdout)) == NULL)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write output data\n",
		       *argv);
      }
      else
      {
	histDom = outObj->domain.hist;
	switch(histDom->type)
	{
	  case WLZ_HISTOGRAMDOMAIN_INT:
	    for(idx = 0; idx < histDom->nBins; ++idx)
	    {
	      fprintf(fP, "%8g %8d\n",
		      histDom->origin + (idx * histDom->binSize),
		      *(histDom->binValues.inp + idx));
	    }
	    break;
	  case WLZ_HISTOGRAMDOMAIN_FLOAT:
	    for(idx = 0; idx < histDom->nBins; ++idx)
	    {
	      fprintf(fP, "%8g %8g\n",
		      histDom->origin + (idx * histDom->binSize),
		      *(histDom->binValues.dbp + idx));
	    }
	    break;
	}
      }
    }
    else
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
    "Usage: %s%sExample 1: %s%sExample 2: %s%s",
    *argv,
    " [-h] [-a] [-o<out file>] [-n#] [-s#] [-w#] [<in object>]\n"
    "Options:\n"
    "  -a  Output the histogram as ascii data not a Woolz histogram object.\n"
    "  -o  Output object/data file name.\n"
    "  -n  Number of histogram bins.\n"
    "  -s  Start value, lower limit of histogram.\n"
    "  -w  Histogram bin width.\n"
    "  -h  Help, prints this usage message.\n"
    "Computes a Woolz histogram object from the input Woolz object.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " -a myobj.wlz | xgraph\n"
    "The input Woolz object is read from myobj.wlz, a histogram is\n"
    "computed and written as ascii data for display by the program xgraph.\n",
    *argv,
    " -o myHist.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, a histogram is computed\n"
    "and written as a Woolz histogram object to the file myHist.wlz\n");
  }
  return(!ok);
}
