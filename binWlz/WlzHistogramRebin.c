#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzHistogramRebin_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzHistogramRebin.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	Re-bins a histogram.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzhistogramrebin "WlzHistogramRebin"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzhistogramrebin WlzHistogramRebin
\par Name
WlzHistogramRebin - Re-bins a Woolz histogram object.
\par Synopsis
\verbatim
WlzHistogramRebin [-o<output object file>] [-a] [-h]
           [-i] [-f] [-n<bins>] [-s<start>] [-w<width>]
           [<input object file>]


\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>output the histogram as ascii data not a Woolz histogram object.</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>create integer bin values.</td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>create (double precision) floating point bin values.</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>number of histogram bins.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>start value, lower limit of histogram.</td>
  </tr>
  <tr>
    <td><b>-w</b></td>
    <td>width of the histogram bins.</td>
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
put.  The default number of bins, origin and bin width values are those
of the input histogram object.
\par Description
Re-bins a Woolz histogram object as required.
\par
Output ascii data are written as a set of records one  per  line,  with
each record composed of two fields separated by white space characters:
\par
\<grey value\>   \<histogram bin occupancy\>

\par Examples
\verbatim
# An example which uses WlzHistogramRebin and xgraph to plot a
# histogram with 256 bins, starting at 0 and with a bin width
# of 1.

WlzHistogramRebin -a -s0 -n256 -w1 myHist.wlz | xgraph
\endverbatim

\par File
\ref WlzHistogramRebin.c "WlzHistogramRebin.c"
\par See Also
\ref wlzhistogramobj "WlzHistogramObj(1)"
\ref WlzHistogramRebin "WlzHistogramRebin(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
		nBins,
		asciiFlag = 0,
		typeFlag = 0,
		nBinsFlag = 0,
		originFlag = 0,
		widthFlag = 0,
		ok = 1,
		usage = 0;
  WlzObjectType	binType;
  double	binOrigin = 0.0,
  		binWidth = 0.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzHistogramDomain *histDom;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "aifo:n:s:w:h",
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
      case 'i':
        typeFlag = 1;
	binType = WLZ_HISTOGRAMDOMAIN_INT;
	break;
      case 'f':
        typeFlag = 1;
	binType = WLZ_HISTOGRAMDOMAIN_FLOAT;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'n':
	nBinsFlag = 1;
	if(sscanf(optarg, "%d", &nBins) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
	originFlag = 1;
	if(sscanf(optarg, "%lg", &binOrigin) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'w':
	widthFlag = 1;
	if(sscanf(optarg, "%lg", &binWidth) != 1)
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

    histDom = inObj->domain.hist;
    binType = (typeFlag)? binType: histDom->type; 
    nBins = (nBinsFlag)? nBins: histDom->nBins; 
    binOrigin = (originFlag)? binOrigin: histDom->origin; 
    binWidth = (widthFlag)? binWidth: histDom->binSize; 
    if(((outObj = WlzAssignObject(WlzHistogramRebin(inObj, binType,
    						    nBins, nBins,
						    binOrigin, binWidth,
						   &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to rebin histogram (%s)\n",
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
	  default:
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
	 ((errNum = WlzWriteObj(fP, outObj)) == 0))
      {
	ok = 0;
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
    " [-h] [-a] [-o<out file>] [-i] [-f] [-n#] [-s#] [-w#] [<in object>]\n"
    "Options:\n"
    "  -a  Output the histogram as ascii data not a Woolz histogram object.\n"
    "  -o  Output object/data file name.\n"
    "  -i  Integer bin values.\n"
    "  -f  Floating point (double precision) bin values.\n"
    "  -n  Number of histogram bins.\n"
    "  -s  Start value, lower limit of histogram.\n"
    "  -w  Histogram bin width.\n"
    "  -h  Help, prints this usage message.\n"
    "Re-bins a Woolz histogram object.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " -s 0.0 -w 10.0 -n 30 -a myhist.wlz | xgraph\n"
    "The input Woolz histogram object is read from myhist.wlz, re-bined\n"
    "so that it starts at 0.0, has 30 bins of width 10.0, and is written\n"
    "as ascii to the program xgraph for display.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
