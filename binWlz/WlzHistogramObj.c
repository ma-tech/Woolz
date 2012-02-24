#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzHistogramObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzHistogramObj.c
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
* \brief	Computes a histogram from a 2D or 3D domain object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzhistogramobj "WlzHistogramObj"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzhistogramobj WlzHistogramObj
\par Name
WlzHistogramObj -  computes a Woolz histogram object from a 2D or 3D
       domain object.

\par Synopsis
\verbatim
WlzHistogramObj [-o<output object file>] [-a] [-h]
           [-n<bins>] [-s<start>] [-w<width>]
           [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>output the histogram as ascii data not a Woolz histogram object.</td>
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
By default the input object is read from the  standard  input  and  the
output  histogram  data/object   is written to the standard output.  If
the number of histogram bins is not specified then the number of  bins,
the origin and the bin width are computed as follows:
\verbatim
    number of bins     =     ceil(max(g) - min(g) + 1.0)
    start value        =     min(g)
    bin width          =     1.0
\endverbatim
where  max(g) and min(g) are the maximum and minimum grey values in the
domain object.
\par Description
Computes a Woolz histogram object from the input Woolz 2D or 3D  domain
object.
\par
Output  ascii  data  are written as a set of records one per line, with
each record composed of two fields separated by white space characters:
\par
    \<grey value\>   \<histogram bin occupancy\>
\par Examples
\verbatim
# An example which uses WlzHistogramObj and xgraph to plot the
# histogram of a domain object.

WlzHistogramObj -a myObj.wlz | xgraph

\endverbatim

\par File
\ref WlzHistogramObj.c "WlzHistogramObj.c"
\par See Also
\ref WlzHistogramObj "WlzHistogram(3)"
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
