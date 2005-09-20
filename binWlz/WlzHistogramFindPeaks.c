#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzHistogramFindPeaks.c
* \author       Bill Hill
* \date         February 2000
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
* \brief	Finds histogram peak positions.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzhistogramfindpeaks "WlzHistogramFindPeaks"
*/

/*!
\ingroup BinWlz
\defgroup wlzhistogramfindpeaks WlzHistogramFindPeaks
\par Name
WlzHistogramFindPeaks - finds histogram peak positions.
\par Synopsis
\verbatim
WlzHistogramFindPeaks [-h] [-o<out file>] [-f #] [-s#] [-t#]
                      [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output data file name.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Features to find, either p (peak), t (trough) or p,t (both).</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Gaussian sigma (standard deviation).</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Threshold value for peak height.</td>
  </tr>
</table>
\par Description
Finds the positions of histogram peaks.
Objects/data are read from stdin and written to stdout unless the
filenames are given.
\par Examples
\verbatim
WlzHistogramFindPeaks -s2 -t 10 myhist.wlz
The input Woolz histogram object is read from myhist.wlz and peak
positions are written to the standard output.
\endverbatim
\par File
\ref WlzHistogramFindPeaks.c "WlzHistogramFindPeaks.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzhistogramobj "WlzHistogramObj(1)"
\ref wlzcompthresh "WlzCompThresh(1)"
\ref WlzHistogramFindPeaks "WlzHistogramFindPeaks(3)"
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
		pkCnt = 0,
		ok = 1,
		usage = 0;
  double	sigma = 1.0,
  		thresh = 1.0;
  int		*pkPos = NULL;
  WlzHistFeature feat = WLZ_HIST_FEATURE_PEAK;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzHistogramDomain *histDom;
  WlzObject	*inObj = NULL;
  char 		*argStr0,
  		*outDatFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "f:o:s:t:h",
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
      case 'f':
	feat = WLZ_HIST_FEATURE_NONE;
	if((argStr0 = strtok(optarg, ", \t")) == NULL)
	{
	  usage = 1;
	  ok = 0;
	}
	else
	{
	  while(ok && argStr0)
	  {
	    if(strcmp(argStr0, "p") == 0)
	    {
	      feat |= WLZ_HIST_FEATURE_PEAK;
	    }
	    else if(strcmp(argStr0, "t") == 0)
	    {
	      feat |= WLZ_HIST_FEATURE_TROUGH;
	    }
	    else
	    {
	      usage = 1;
	      ok = 0;
	    }
	    if(ok)
	    {
	      argStr0 = strtok(NULL, ", \t");
	    }
	  }
	}
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
    				   &pkCnt, &pkPos, feat);
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
    " [-h] [-o<out file>] [-f #] [-s#] [-t#] [<in object>]\n"
    "Options:\n"
    "  -o  Output data file name.\n"
    "  -f  Features to find, either p (peak), t (trough) or p,t (both).\n"
    "  -s  Gaussian sigma (standard deviation).\n"
    "  -t  Threshold value for peak height.\n"
    "  -h  Help, prints this usage message.\n"
    "Finds the positions of histogram peaks.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " -s2 -t 10 myhist.wlz\n"
    "The input Woolz histogram object is read from myhist.wlz and peak\n"
    "positions are written to the standard output.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
