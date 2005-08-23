#pragma ident "MRC HGU $Id$"
/*!
\ingroup      BinWlz
\defgroup     wlzhistogramsmooth WlzHistogramSmooth
\par Name
WlzHistogramSmooth -  Low-pass filters a Woolz histogram object.
\par Synopsis
\verbatim
WlzHistogramSmooth [-o<output object file>] [-a] [-h]
           [-s<smoothing>] [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>output the histogram as ascii data not a Woolz histogram object.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>smoothing factor. </td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>
By  default  the input histogram object is read from the standard input
and the output histogram data/object  is written to the  standard  out-
put.
\par
The  smoothing  factor is the gaussian kernel half-height full-width in
bins (default 1).

\par Description
Smooths a Woolz histogram object by convolving the histogram bins  with
a gaussian.
\par
Output  ascii  data  are written as a set of records one per line, with
each record composed of two fields separated by white space characters:
\par
    \<grey value\>   \<histogram bin occupancy\>

\par Examples
\verbatim
# An example which uses WlzHistogramSmooth and xgraph to plot a
# smoothed histogram.

WlzHistogramSmooth -a -s5 myHist.wlz | xgraph

\endverbatim

\par See Also
WlzHistoGram(3)
\par Bugs
None known
\author       richard <Richard.Baldock@hgu.mrc.ac.uk>
\date         Fri Jul 29 11:30:35 2005
\version      MRC HGU $Id$
              $Revision$
              $Name$
\par Copyright:
             1994-2003 Medical Research Council, UK.
              All rights reserved.
\par Address:
              MRC Human Genetics Unit,
              Western General Hospital,
              Edinburgh, EH4 2XU, UK.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/***********************************************************************
* Project:      Woolz
* Title:        WlzHistogramSmooth.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which smooths a Woolz histogram object.
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
		asciiFlag = 0,
		smoothing = 1,
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzHistogramDomain *histDom;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "ao:s:h",
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
      case 's':
	if(sscanf(optarg, "%d", &smoothing) != 1)
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
    outObj = WlzAssignObject(inObj, NULL);
    if((errNum = WlzHistogramSmooth(outObj, smoothing)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to smooth histogram (%s).\n",
		     *argv, errMsg);
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
		       "%s: failed to write output object (%s).\n",
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
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-a] [-o<out file>] [-s#] [<in object>]\n"
    "Options:\n"
    "  -a  Output the histogram as ascii data not a Woolz histogram object.\n"
    "  -o  Output object/data file name.\n"
    "  -s  Smoothing factor (bins for gaussian half height full width).\n"
    "  -h  Help, prints this usage message.\n"
    "Smooths a Woolz histogram object.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " -s3 -a myhist.wlz | xgraph\n"
    "The input Woolz histogram object is read from myhist.wlz, smoothed\n"
    "by convolving the histogram bins with a gaussian which has a half\n",
    "height half width of 3 histogram bins, and written as ascii to the\n",
    "program xgraph for display.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
