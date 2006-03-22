#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzHistogramData_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzHistogramData.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Reads a histogram object and writes it as ascii data.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzhistogramdata "WlzHistogramData"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzhistogramdata WlzHistogramData
\par Name
WlzHistogramData - Reads  a  Woolz  histogram object and writes it as
       ascii data.

\par Synopsis
\verbatim
WlzHistogramData [-o<output file>] [-h] [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
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
By default the input histogram object is read from the  standard  input
       and the output data  is written to the standard output.

\par Description
Reads a Woolz histogram object and writes it as ascii data.
The output data are written as a set of records one per line, with each
record composed of two fields separated by white space characters:
\par
\verbatim
    <grey value> <histogram bin occupancy>
\endverbatim

\par Examples
\verbatim
# An example which uses WlzHistogramData and xgraph to plot a
# histogram.

WlzHistogramData myHist.wlz | xgraph

\endverbatim

\par File
\ref WlzHistogramData.c "WlzHistogramData.c"
\par See Also
\ref wlzhistogramobj "WlzHistogramObj(1)"
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
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzHistogramDomain *histDom;
  WlzObject	*inObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "o:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
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
      histDom = inObj->domain.hist;
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
    if(fP && strcmp(outObjFileStr, "-"))
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
    " [-h] [-o<out file>] [<in object>]\n"
    "Options:\n"
    "  -o  Output data file name.\n"
    "  -h  Help, prints this usage message.\n"
    "Reads a Woolz histogram object and writes it as ascii data.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " myHist.wlz | xgraph\n"
    "The input Woolz histogram object is read from myHist.wlz, and\n"
    "written as ascii data for display by the program xgraph.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
