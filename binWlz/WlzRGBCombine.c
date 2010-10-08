#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRGBCombine_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzRGBCombine.c
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
* \brief	Combines upto three domain objects with values to make a
*		single RGBA domain object in which the input objects are
*		used for the RGB components.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzrgbcombine "WlzRGBCombine"
*/

/*!
\ingroup BinWlz
\defgroup wlzrgbcombine WlzRGBCombine
\par Name
WlzRGBCombine - Combines domain objects to make a single RGBA domain object.
\par Synopsis
\verbatim
WlzRGBCombine [-h] [-i] [-o<out file>] [-n] [-m <match object>]
              [<red object 0>] [<green object 1>] [<blue object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Invert grey values.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Normalise the grey value range of the objects prior to combining them
        instead of the default which is to clamp them.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Match the histograms of the input objects to that of the given
        object.</td>
  </tr>
</table>
\par Description
Combines the input domain objects with values to form an RGBA object
in which the red, green and blue components are the input objects.
Histogram matching and grey value normalisation may be used to allow
better visual comparison between the components.
By default the input grey values are clamped to the range 0 - 255,
but they will be normalised to this range if the normalise flag is set.
The optional object for histogram matching may either be a domain object
with values or a histogram.
By default the objects are read from the standard input (in the order
given) and written to the standard output. Null objects may be specified
using the word null inplace of a object filename.
\par Example
\verbatim
WlzRGBCombine -m hist.wlz -o rgb.wlz red.wlz null blue.wlz
\endverbatim
Creates a new object rgb.wlz in which the red and blue channels are
created from the grey value objects read from red.wlz and blue.wlz.
The green channel is left empty because null is given on the command
line. The read and blue channels have their values matched to the
histogram of hist.wlz.
\par File
\ref WlzRGBCombine.c "WlzRGBCombine.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCompoundToRGBA "WlzCompoundToRGBA(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>


/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int             main(int argc, char **argv)
{
  int		idx,
		option,
		invFlg = 0,
		nrmFlg = 0,
		ok = 1,
		usage = 0;
  FILE		*fP = NULL;
  WlzPixelV	gMin,
  		gMax;
  WlzObject	*tmpObj = NULL,
  		*histObj = NULL,
  		*outObj = NULL;
  WlzObject	*inObj[4];
  WlzCompoundArray *cpdObj = NULL;
  WlzObjectType	lastObjType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*histObjFileStr = NULL,
  		*outObjFileStr = NULL;
  char  	*inObjFileStr[3];
  const char	*errMsg;
  static char	optList[] = "hinsm:o:",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObj[0] = NULL;
  inObj[1] = NULL;
  inObj[2] = NULL;
  inObj[3] = NULL;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr[0] = inObjFileStrDef;
  inObjFileStr[1] = inObjFileStrDef;
  inObjFileStr[2] = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
    case 'i':
      invFlg = 1;
      break;
    case 'n':
      nrmFlg = 1;
      break;
    case 'm':
      histObjFileStr = optarg;
      break;
    case 'o':
      outObjFileStr = optarg;
      break;
    case 'h': /* FALLTHROUGH */
    default:
      usage = 1;
      ok = 0;
      break;
    }
  }
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    idx = 0;
    while((idx < 3) && (optind < argc))
    {
      inObjFileStr[idx] = argv[optind];
      ++optind;
      ++idx;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  /* Read histogram match object if required. */
  if(ok && (histObjFileStr != NULL))
  {
    errNum = WLZ_ERR_READ_EOF;
    if((*histObjFileStr == '\0') ||
       ((fP = (strcmp(histObjFileStr, "-")?
		fopen(histObjFileStr, "r"): stdin)) == NULL) ||
       ((histObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to read histogram object from file %s (%s)\n",
		     argv[0], inObjFileStr[idx], errMsg);
    }
    if(fP && strcmp(histObjFileStr, "-"))
    {
      fclose(fP);
      fP = NULL;
    }
    /* Check object type and create histogram if required. */
    if(ok)
    {
      switch(histObj->type)
      {
	case WLZ_HISTOGRAM:
	  break;
        case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
	case WLZ_3D_DOMAINOBJ:
	  tmpObj = WlzHistogramObj(histObj, 256, 0.0, 1.0, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    (void )WlzFreeObj(histObj);
	    histObj = tmpObj;
	  }
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to create histogram object (%s)\n",
		     argv[0], errMsg);
    }
  }
  /* Read input objects. */
  if(ok)
  {
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 3))
    {
      if(strcmp(inObjFileStr[idx], "null") != 0)
      {
	errNum = WLZ_ERR_READ_EOF;
	if((*inObjFileStr[idx] == '\0') ||
	   ((fP = (strcmp(inObjFileStr[idx], "-")?
		    fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	   ((inObj[idx]= WlzAssignObject(WlzReadObj(fP,
						     &errNum), NULL)) == NULL))
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: Failed to read object %d from file %s (%s)\n",
			 argv[0], idx, inObjFileStr[idx], errMsg);
	}
	if(fP && strcmp(inObjFileStr[idx], "-"))
	{
	  fclose(fP);
	  fP = NULL;
	}
      }
      ++idx;
    }
  }
  /* Check input objects. */
  if(ok)
  {
    idx = 0;
    lastObjType = WLZ_NULL;
    while((errNum == WLZ_ERR_NONE) && (idx < 4))
    {
      if(inObj[idx] == NULL)
      {
	inObj[idx] = WlzAssignObject(WlzMakeEmpty(&errNum), NULL);
      }
      else
      {
	if((inObj[idx]->type != WLZ_2D_DOMAINOBJ) &&
	   (inObj[idx]->type != WLZ_3D_DOMAINOBJ))
	{
	  ok = 0;
	  (void )fprintf(stderr, "%s: Input object %d not a domain object.\n",
			 argv[0], idx);
	}
	else if(inObj[idx]->domain.core == NULL)
	{
	  ok = 0;
	  (void )fprintf(stderr, "%s: Input object %d has null domain.\n",
			 argv[0], idx);
	}
	else if(inObj[idx]->values.core == NULL)
	{
	  ok = 0;
	  (void )fprintf(stderr, "%s: Input object %d has null values.\n",
			 argv[0], idx);
	}
	else if((lastObjType != WLZ_NULL) &&
	        (inObj[idx]->type != lastObjType))
	{
	  ok = 0;
	  (void )fprintf(stderr, "%s: Input objects have different types.\n",
			 argv[0]);
	}
	else
	{
	  lastObjType = inObj[idx]->type;
	}
      }
      ++idx;
    }
    if(lastObjType == WLZ_NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: No input objects given.\n", argv[0]);
    }
  }
  /* Process input objects as required. */
  if(ok && (invFlg || nrmFlg || (histObj != NULL)))
  {
    idx = 0;
    lastObjType = WLZ_NULL;
    while((errNum == WLZ_ERR_NONE) && (idx < 3))
    {
      if(inObj[idx]->type != WLZ_EMPTY_OBJ)
      {
	if(invFlg)
	{
	  errNum = WlzGreyRange(inObj[idx], &gMin, &gMax);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzGreyInvertMinMax(inObj[idx], gMin, gMax);
	  }
	}
	if(histObj)
	{
	  errNum = WlzHistogramMatchObj(inObj[idx], histObj, 0, 0, 0.0, 1.0,
	                                1);
	}
	if(nrmFlg)
	{
	  errNum = WlzGreyNormalise(inObj[idx], 1);
	}
      }
      ++idx;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to process object %d (%s)\n",
		     argv[0], idx, errMsg);
    }
  }
  /* Create a compound object. */
  if(ok)
  {
    cpdObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_2, 3, 4, inObj, lastObjType,
                                  &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	      "%s: Failed to create compound array from input objects (%s).\n",
	      argv[0], errMsg);
    }
  }
  /* Create RGBA object. */
  if(ok)
  {
    outObj = WlzCompoundToRGBA(cpdObj, WLZ_RGBA_SPACE_RGB, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to create RGBA object (%s).\n",
		     argv[0], errMsg);
    }
    (void )WlzFreeObj(tmpObj);
  }
  /* Output the RGBA object. */
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?  fopen(outObjFileStr, "w"):
	      				    stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write output object (%s).\n",
		     argv[0], errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  /* Free objects. */
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  (void )WlzFreeObj(inObj[2]);
  (void )WlzFreeObj(inObj[3]);
  (void )WlzFreeObj(outObj);
  (void )WlzFreeObj(histObj);
  (void )WlzFreeObj((WlzObject *)cpdObj);
  /* Report usage if required. */
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-i] [-o<out file>] [-n] [-m <match object>]\n"
    "       [<red object 0>] [<green object 1>] [<blue object>]\n"
    "Options:\n"
    "  -h        Help, prints this usage message.\n"
    "  -i        Invert grey values.\n"
    "  -o        Output file name.\n"
    "  -n        Normalise the grey value range of the objects prior to\n"
    "            combining them instead of the default which is to clamp\n"
    "            them.\n"
    "  -m        Match the histograms of the input objects to that of the\n"
    "            given object.\n"
    "Combines the input domain objects with values to form an RGBA object\n"
    "in which the red, green and blue components are the input objects. \n"
    "Histogram matching and grey value normalisation may be used to allow\n"
    "better visual comparison between the components.  By default the input\n"
    "grey values are clamped to the range 0 - 255, but they will be\n"
    "normalised to this range if the normalise flag is set. The optional\n"
    "object for histogram matching may either be a domain object with\n"
    "values or a histogram.\n"
    "By default the objects are read from the standard input (in the\n"
    "order given) and written to the standard output. Null objects may\n"
    "be specified using the word null inplace of a object filename.\n"
    "Example:\n"
    "  %s -m hist.wlz -o rgb.wlz red.wlz null blue.wlz\n"
    "Creates a new object rgb.wlz in which the red and blue channels\n"
    "are created from the grey value objects read from red.wlz and\n"
    "blue.wlz. The green channel is left empty because null is given\n"
    "on the command line. The read and blue channels have their values\n"
    "matched to the histogram of hist.wlz.\n",
    argv[0], argv[0]);
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
