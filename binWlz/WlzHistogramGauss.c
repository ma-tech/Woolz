#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzHistogramGauss.c
* \author       Bill Hill
* \date         August 2005
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
* \brief	Convolves the bin values of a histogram object with
*               a Gaussian's 0th, 1st or nd derivative.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzhistogramgauss "WlzHistogramGauss"
*/

/*!
\ingroup BinWlz
\defgroup wlzhistogramgauss WlzHistogramGauss
\par Name
WlzHistogramGauss - convolves a histogram with a Gaussian.
\par Synopsis
\verbatim
WlzHistogramGauss [-h] [-a] [-r] [-o<out file>] [-d#] [-s#] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Use a recursive IIR filter.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object/data file name.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Derivative, range [0-2].</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Gaussian sigma (standard deviation).</td>
  </tr>
</table>
\par Description
Convolves the bin values of a Woolz histogram object with a Gaussian
or it's 1st or 2nd derivtive.
Objects/data are read from stdin and written to stdout unless the
filenames are given.
\par Examples
\verbatim
WlzHistogramGauss -s3 -a myhist.wlz | xgraph
\endverbatim
The input Woolz histogram object is read from myhist.wlz, filtered
by convolving the histogram bins with a Gaussian (sigma = 3.0) and
written as ascii to the program xgraph for display.
\par File
\ref WlzHistogramGauss.c "WlzHistogramGauss.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzhistogramobj "WlzHistogramObj(1)"
\ref wlzhistogramsmooth "WlzHistogramSmooth(1)"
\ref WlzHistogramCnvGauss "WlzHistogramCnvGauss(3)"
\ref WlzHistogramRsvGauss "WlzHistogramRsvGauss(3)"
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
		deriv = 0,
		asciiFlag = 0,
		rsvFlag = 0,
		ok = 1,
		usage = 0;
  double	sigma = 1.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzHistogramDomain *histDom;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "aro:d:s:h",
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
      case 'r':
        rsvFlag = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'd':
	if((sscanf(optarg, "%d", &deriv) != 1) ||
	   (deriv < 0) || (deriv > 2))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
	if(sscanf(optarg, "%lg", &sigma) != 1)
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
    if(rsvFlag)
    {
      errNum = WlzHistogramRsvGauss(outObj, sigma, deriv);
    }
    else
    {
      errNum = WlzHistogramCnvGauss(outObj, sigma, deriv);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to filter histogram (%s).\n",
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
    " [-h] [-a] [-r] [-o<out file>] [-d#] [-s#] [<in object>]\n"
    "Options:\n"
    "  -a  Output the histogram as ascii data not a Woolz histogram object.\n"
    "  -r  Use a recursive IIR filter.\n"
    "  -o  Output object/data file name.\n"
    "  -d  Derivative, range [0-2].\n"
    "  -s  Gaussian sigma (standard deviation).\n"
    "  -h  Help, prints this usage message.\n"
    "Convolves the bin values of a Woolz histogram object with a Gaussian\n"
    "or it's 1st or 2nd derivtive.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " -s3 -a myhist.wlz | xgraph\n"
    "The input Woolz histogram object is read from myhist.wlz, filtered\n"
    "by convolving the histogram bins with a Gaussian (sigma = 3.0) and\n"
    "written as ascii to the program xgraph for display.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
