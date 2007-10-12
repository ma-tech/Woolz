#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzScalarFeatures_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzScalarFeatures.c
* \author       Bill Hill
* \date         November 2002
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
* \brief	Extracts scalar features from Woolz domain objects.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzscalarfeatures "WlzScalarFeatures"
*/

/*!
\ingroup BinWlz
\defgroup wlzscalarfeatures WlzScalarFeatures
\par Name
WlzScalarFeatures - extracts scalar features from Woolz domain objects.
\par Synopsis
\verbatim
WlzScalarFeatures [-h] [-o<out file>] [-H] [-L] [-a] [-v#] [-d#]
                  [-t<feature type>] [-f#] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name, default stdout.</td>
  </tr>
  <tr> 
    <td><b>-H</b></td>
    <td>Feature values will be at or above the threshold value, default.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Feature values will be below threshold value, not default.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Compute a threshold value automatically, not default.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Specified threshold value, default 128.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Minimum distance between features, default 20.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Type of feature, valid types are "value" and "grad",
        default "value".</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Filter value for feature extraction, default 1.</td>
  </tr>
</table>
\par Description
Extracts scalar features from Woolz domain objects. These features
are separated by a minimum distance and are either maximal or minimal
within a region (which approximates to the Vorinoi cell) containing
the feature.
\par Examples
\verbatim
WlzScalarFeatures -t grad -d 25 -H -v 200 obj.wlz >out.num
\endverbatim
A 2D domain object with values is read from obj.wlz and a list of
2D cooordinates is then output to the file out.num.
These coordinates are positions within obj.wlz which have high image gradients
and are seperated by at least 25 pixels.
\par File
\ref WlzScalarFeatures.c "WlzScalarFeatures.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzScalarFeatures2D "WlzScalarFeatures2D(3)"
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
  int		tI0,
		idx,
  		autoThr = 0,
		nFeat = 0,
  		option,
  		ok = 1,
		usage = 0;
  WlzScalarFeatureType fType = WLZ_SCALARFEATURE_VALUE;
  WlzThresholdType thrHL = WLZ_THRESH_LOW;
  WlzPixelV	thrV;
  WlzIVertex2	*feat = NULL;
  double	filterV = 1.0,
  		minDist = 20.0;
  char		*tS0,
  		*inObjFileStr,
  		*outFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  const char	*errMsg;
  static char	optList[] = "ahHLd:f:o:t:v:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  thrV.type = WLZ_GREY_DOUBLE;
  thrV.v.dbv = 128;
  inObjFileStr = inObjFileStrDef;
  outFileStr = outFileStrDef;
  while(ok && (usage == 0) &&
  	((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
	autoThr = 1;
	break;
      case 'H':
	thrHL = WLZ_THRESH_HIGH;
	break;
      case 'L':
	thrHL = WLZ_THRESH_LOW;
	break;
      case 'd':
        if((sscanf(optarg, "%lg", &minDist) != 1) || (minDist < 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'f':
        if(sscanf(optarg, "%lg", &filterV) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 't':
	if(WlzStringMatchValue(&tI0, optarg,
			"value", WLZ_SCALARFEATURE_VALUE,
			"grad", WLZ_SCALARFEATURE_GRADIENT,
			NULL))
	{
	  fType = tI0;
	}
	else
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'v':
	if(sscanf(optarg, "%lg", &(thrV.v.dbv)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h': /* FALLTHROUGH */
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
    errNum = WlzScalarFeatures2D(inObj, &nFeat, &feat, fType, thrHL, thrV,
    			 	 filterV, minDist);
    if(errNum != WLZ_ERR_NONE)
    {
    ok = 0;
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
    		   "%s: failed to sub-sample object (%s).\n",
		   *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: failed to write output object (%s).\n",
                     *argv, errMsg);
    }
    else
    {
      for(idx = 0; idx < nFeat; ++idx)
      {
        (void )fprintf(fP, "%d %d\n", (feat + idx)->vtX, (feat + idx)->vtY);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(feat)
  {
    AlcFree(feat);
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: "
    "%s [-o<out file>] [-h] [-H] [-L] [-a] [-v#] [-d#]\n"
    "                   [-t<feature type>] [-f#] [<in object>]\n"
    "Extracts scalar features from Woolz domain objects. These features\n"
    "are separated by a minimum distance and are either maximal or minimal\n"
    "within a region (which approximates to the Vorinoi cell) containing\n"
    "the feature.\n"
    "Options are:\n"
    "  -o Output file name, set to %s.\n"
    "  -h Display this usage information.\n"
    "  -H Feature values will be at or above the threshold value, %sset.\n"
    "  -L Feature values will be below threshold value, %sset.\n"
    "  -a Compute a threshold value automatically, %sset.\n"
    "  -v Specified threshold value, set to %g.\n"
    "  -d Minimum distance between features, set to %g.\n"
    "  -t Type of feature, valid types are \"value\" and \"grad\", set\n"
    "     to %s.\n"
    "  -f Filter value for feature extraction, set to %g\n",
    *argv,
    strcmp(outFileStr, "-")? outFileStr: "<stdout>",
    (thrHL == WLZ_THRESH_HIGH)? "": "not ",
    (thrHL == WLZ_THRESH_LOW)? "": "not ",
    (autoThr)? "": "not ",
    thrV.v.dbv,
    minDist,
    (tS0 = (char *)WlzStringFromScalarFeatureType(fType, NULL))?
           tS0: "Invalid",
    filterV);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
