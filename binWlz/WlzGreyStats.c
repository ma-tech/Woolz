#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyStats.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes simple statistics about a 2D or 3D Woolz
*		domain object's grey values.
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
  int		area,
  		option,
		ok = 1,
		verbose = 0,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  double	min,
  		max,
		sum,
		sumSq,
		mean,
		stdDev;
  WlzGreyType	gType;
  char 		*outFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "o:vh",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
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
		     *argv, inObjFileStr,
		     errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    area = WlzGreyStats(inObj, &gType, &min, &max, &sum, &sumSq, &mean,
    		        &stdDev, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute object's grey value statistics"
		     " (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"):
	      stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to write output file\n",
		     *argv);
    }
    else
    {
      if(verbose)
      {
        (void )fprintf(fP,
		       "area        % 14d\n"
		       "grey type   % 14s\n"
		       "min         % 14g\n"
		       "max         % 14g\n"
		       "sum         % 14g\n"
		       "sum sq      % 14g\n"
		       "mean        % 14g\n"
		       "std dev     % 14g\n",
		       area, WlzStringFromGreyType(gType, NULL),
		       min, max, sum, sumSq, mean, stdDev);
      }
      else
      {
      (void )fprintf(fP, "%d %s %g %g %g %g %g\n",
      		     area, WlzStringFromGreyType(gType, NULL),
		     min, max, sum, sumSq, mean, stdDev);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out file>] [-v] [-h] [<in object>]\n"
    "Options:\n"
    "  -o  Output file name.\n"
    "  -v  Verbose output flag.\n"
    "  -h  Help, prints this usage message.\n"
    "Calculates the area, minimum, maximum, sum, sum of squares, mean\n"
    "and the standard deviation of the input 2D or 3D Woolz object's\n"
    "grey values.\n"
    "If the verbose output flag is not set then the following are written\n"
    "to the output file in the following order:\n"
    "  <area> <grey type> <min> <max> <sum> <sum of sq> <mean> <std dev>.\n"
    "The input object is read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    *argv,
    " -o stats.txt myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz. The statistics are\n"
    "calculated and  written to out.txt.\n");
  }
  return(!ok);
}
