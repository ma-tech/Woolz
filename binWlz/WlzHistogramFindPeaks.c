#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzHistogramGauss.c
* Date:         February 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which finds histogram peak positions.
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
		pkCnt = 0,
		ok = 1,
		usage = 0;
  double	sigma = 1.0,
  		thresh = 1.0;
  int		*pkPos = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzHistogramDomain *histDom;
  WlzObject	*inObj = NULL;
  char 		*outDatFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "o:s:t:h",
		outDatFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outDatFileStr = outDatFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outDatFileStr = optarg;
	break;
      case 's':
	if(sscanf(optarg, "%lg", &sigma) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 't':
	if(sscanf(optarg, "%lg", &thresh) != 1)
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
     (outDatFileStr == NULL) || (*outDatFileStr == '\0'))
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
    if(inObj->type != WLZ_HISTOGRAM)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: input object read from file %s not a histogram\n",
		     *argv, inObjFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzHistogramFindPeaks(inObj, sigma, thresh,
    				   &pkCnt, &pkPos);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to find histogram peaks (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if(pkCnt > 0)
    {
      if((fP = (strcmp(outDatFileStr, "-")?
               fopen(outDatFileStr, "w"): stdout)) == NULL)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write output data\n",
		       *argv);
      }
      else
      {
	histDom = inObj->domain.hist;
	for(idx = 0; idx < pkCnt; ++idx)
	{
	  fprintf(fP, "%8g\n",
		  histDom->origin + (*(pkPos + idx) * histDom->binSize));
	}
      }
    }
    if(fP && strcmp(outDatFileStr, "-"))
    {
      fclose(fP);
    }
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
    " [-h] [-o<out file>] [-s#] [-t#] [<in object>]\n"
    "Options:\n"
    "  -o  Output data file name.\n"
    "  -s  Gaussian sigma (standard deviation).\n"
    "  -t  Threshold value for peak height.\n",
    "  -h  Help, prints this usage message.\n"
    "Finds the positions of histogram peaks.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " -s2 -t 10 myhist.wlz\n"
    "The input Woolz histogram object is read from myhist.wlz and peak\n"
    "positions are written th the standard output.\n");
  }
  return(!ok);
}
