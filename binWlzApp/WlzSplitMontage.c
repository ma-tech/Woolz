#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzSplitImage_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzApp/WlzSplitMontage.c
* \author       Bill Hill
* \date         October 2004
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
* \brief	Splits a montage object - one that is composed of several
*		components separated by a gap.
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzSplitMontage "WlzSplitMontage"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzSplitMontage WlzSplitMontage
\par Name
WlzSplitMontage - splits a montage image, one that is composed of several
                axis aligned rectangular component images separated a gap.
\par Synopsis
\verbatim
WlzSplitMontage [-h] [-v] [-a <min area>] [-b <base>] [-g <gap value>]
                [-f <fmt>] [-F <fmt>] [-t <tol>]
		[-w <border width>] <path to image file>
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Be verbose, probably only useful for debugging.</td>
  </tr>
  <tr> 
    <td><b>-A</b></td>
    <td>Use letters rather than numbers to distinguish component files.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Minimum area for the each of the extracted components.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Base for output file name (the default is the input file name
        base).</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Input image format, default is to use the input file extension.</td>
  </tr>
  <tr> 
    <td><b>-F</b></td>
    <td>Output image format, default is same as input.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Border value which must either be a single value of a
        comma separated triple of values (red,green,blue).</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Tolerance (percentage) for the gap value.</td>
  </tr>
  <tr> 
    <td><b>-w</b></td>
    <td>Width of border added to extracted components.</td>
  </tr>
</table>
\par Description
Splits a montage image, one that is composed of several axis aligned
rectangular component images separated by a gap, with the gap having
either a  single value or a  narrow range of values.  The component
images are each output to a separate file in raster scan order from
top left to bottom right (approximately).
\par Examples
\verbatim
WlzSplitMontage -g 255 img.jpg
\endverbatim
Splits the montage file img.jpg,
which has pure white gaps separating the component images.
The component images are
written to the files img000001.jpg, img000002.jpg, etc..
\par File
\ref WlzSplitMontage.c "WlzSplitMontage.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzLabel "WlzLabel(1)"
\ref WlzThreshold "WlzThreshold(1)"
\ref WlzSplitImage "WlzSplitImage(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <Wlz.h>
#include <WlzExtFF.h>


extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

static int 			WlzSplitImageSortBBFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
static void			WlzSplitImageIdxToAlpha(
 				  char *idxBuf,
				  int nComp,
				  int idx);

int             main(int argc, char *argv[])
{
  int		idC,
		len,
		nComp,
  		option,
		useAlpha = 0,
		bWidth = 2,
		minArea = 2500,
  		ok = 1,
		usage = 0,
		verbose = 0;
  double	tol = 0.05;
  int		*orderTb = NULL;
  char		*ext,
  		*dot,
  		*sep,
		*str,
		*pathDir,
		*base = NULL,
  		*gapStr = NULL,
		*inFile = NULL,
  		*inPath = NULL;
  int		tI[4];
  WlzEffFormat	inFmt = WLZEFF_FORMAT_NONE,
  		outFmt = WLZEFF_FORMAT_NONE;
  WlzObject	*inObj = NULL;
  WlzObject	**compObj = NULL;
  WlzPixelV	gapVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		idxBuf[64],
  		pathBuf[FILENAME_MAX];
  const char	*errMsg;
  const int	maxComp = 1024;
  static char	optList[] = "Adhva:b:f:F:g:t:w:",
  		defPathDir[] = ".";

  gapVal.type = WLZ_GREY_INT;
  gapVal.v.inv = 255;
  gapStr = AlcStrDup("255");
  while((usage == 0) && ok &&
        ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        if(sscanf(optarg, "%d", &minArea) != 1)
	{
	  usage = 1;
	}
        break;
      case 'A':
        useAlpha = 1;
	break;
      case 'b':
        base = optarg;
	break;
      case 'h':
        usage = 1;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'f':
        if((inFmt = WlzEffStringExtToFormat(optarg)) == WLZEFF_FORMAT_NONE)
	{
	  usage = 1;
	}
	break;
      case 'F':
        if((outFmt = WlzEffStringExtToFormat(optarg)) == WLZEFF_FORMAT_NONE)
	{
	  usage = 1;
	}
	break;
      case 'g':
	gapStr = AlcStrDup(optarg);
	if(strchr(optarg, ',') == NULL)
	{
	  gapVal.type = WLZ_GREY_INT;
	  if(sscanf(optarg, "%d", &(gapVal.v.inv)) != 1)
	  {
	    usage = 1;
	  }
	}
	else
	{
	  if((str = strtok(optarg, ",")) != NULL)
	  {
	    if(sscanf(str, "%d", tI + 0) != 1)
	    {
	      str = NULL;
	    }
	    else
	    {
	      if((str = strtok(NULL, ",")) != NULL)
	      {
		if(sscanf(str, "%d", tI + 1) != 1)
		{
		  str = NULL;
		}
		else
		{
		  if((str = strtok(NULL, ",")) != NULL)
		  {
		    if(sscanf(str, "%d", tI + 2) != 1)
		    {
		      str = NULL;
		    }
		  }
		}
	      }
	    }
	  }
	  if(str == NULL)
	  {
	    usage = 1;
	  }
	  else
	  {
	    gapVal.type = WLZ_GREY_RGBA;
	    WLZ_RGBA_RGBA_SET(gapVal.v.rgbv, tI[0], tI[1], tI[2], 255);
	  }
	}
	break;
      case 't':
        if(sscanf(optarg, "%lg", &tol) == 1)
	{
	  tol = fabs(tol) * 0.01; /* Convert from percentage to fraction. */
	}
	else
	{
	  usage = 1;
	}
        break;
      case 'w':
        if(sscanf(optarg, "%d", &bWidth) != 1)
	{
	  usage = 1;
	}
        break;
      default:
        usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok)
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      inPath = *(argv + optind);
    }
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "inPath = %s\n", inPath);
  }
  /* Parse input file path into path + name + ext. */
  if(ok)
  {
    if(((len = strlen(inPath)) < 3) || /* Minimum of 3 chars. */
        (*(inPath + len - 1) == '/'))  /* Directory not plain file. */
    {
      usage = 1;
      ok = 0;
    }
  }
  if(ok)
  {
    sep = strrchr(inPath, '/');
    dot = strrchr(inPath, '.');
    if(sep == NULL)
    {
      pathDir = defPathDir;
      inFile = inPath;
    }
    else
    {
      *sep = '\0';
      pathDir = inPath;
      inFile = sep + 1;
    }
    if(dot)
    {
      *dot = '\0';
      ext = dot + 1;
    }
    else
    {
      ext = NULL;
    }
    if(base == NULL)
    {
      base = inFile;
    }
  }
  if(ok)
  {
    if(inFmt == WLZEFF_FORMAT_NONE)
    {
      inFmt = WlzEffStringExtToFormat(ext);
      if(inFmt == WLZEFF_FORMAT_NONE)
      {
        usage = 1;
	ok = 0;
      }
    }
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "pathDir = %s\n", pathDir);
    (void )fprintf(stderr, "inFile = %s\n", inFile);
    (void )fprintf(stderr, "ext = %s\n", (ext)? ext: "(null)");
    (void )fprintf(stderr, "inFmt = %s\n",
    		   WlzEffStringFromFormat(inFmt, NULL));
  }
  /* Read image. */
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if(ext)
    {
      sprintf(pathBuf, "%s/%s.%s", pathDir, inFile, ext);
    }
    else
    {
      sprintf(pathBuf, "%s/%s", pathDir, inFile);
    }
    if((inObj= WlzAssignObject(WlzEffReadObj(NULL, pathBuf, inFmt, 0,
    				             &errNum), NULL)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to read object from file %s (%s)\n",
		     *argv, pathBuf, errMsg);
    }
    else
    {
      if(outFmt == WLZEFF_FORMAT_NONE)
      {
	outFmt = inFmt;
      }
      (void )WlzEffStringFromFormat(outFmt, (const char **)&ext);
    }
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "read input image %s ok.\n", pathBuf);
  }
  /* Split image. */
  if(ok)
  {
    errNum = WlzSplitMontageObj(inObj, gapVal, tol, bWidth, minArea,
    			        maxComp, &nComp, &compObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to split image (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "found %d component images.\n", nComp);
  }
  /* Sort the objects into left -> right then top to bottom order. */
  if(ok)
  {
    if((orderTb = (int *)AlcMalloc(sizeof(int) * nComp)) == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to sort components (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    for(idC = 0; idC < nComp; ++idC)
    {
      orderTb[idC] = idC;
    }
    (void )AlgHeapSortIdx(compObj, orderTb, nComp, WlzSplitImageSortBBFn);
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "order table is:\n");
    for(idC = 0; idC < nComp; ++idC)
    {
      (void )fprintf(stderr, " % 6d", orderTb[idC]);
      if(idC % 10 == 9)
      {
        (void )fprintf(stderr, "\n");
      }
    }
    if(idC % 10 != 0)
    {
      (void )fprintf(stderr, "\n");
    }
  }
  /* Write the components to files. */
  if(ok)
  {
    idC = 0;
    while((errNum == WLZ_ERR_NONE) && (idC < nComp))
    {
      /* Compute output file path.. */
      if(useAlpha)
      {
	WlzSplitImageIdxToAlpha(idxBuf, nComp, idC + 1);
      }
      else
      {
	sprintf(idxBuf, "06d", idC + 1);
      }
      if(ext)
      {
	(void )sprintf(pathBuf, "%s%s.%s",
		       base, idxBuf, ext);
      }
      else
      {
	(void )sprintf(pathBuf, "%s%s",
		       base, idxBuf);
      }
      /* Write image component to file. */
      if((errNum = WlzEffWriteObj(NULL, pathBuf,
      			          compObj[orderTb[idC]],
				  outFmt)) != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: failed to write object to file(s) %s (%s).\n",
		       argv[0], pathBuf, errMsg);
      }
      ++idC;
    }
  }
  /* Tidy up. */
  AlcFree(orderTb);
  (void )WlzFreeObj(inObj);
  if(compObj)
  {
    for(idC = 0; idC < nComp; ++idC)
    {
      (void )WlzFreeObj(*(compObj + idC));
    }
    AlcFree(compObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-v] [-a <min area>] [-A] [-b <base>] [-g <gap\n"
    "                       value> [-f <fmt>] [-F <fmt>] [-t <tol>]\n"
    "                       [-w <border width>] <path to image file>\n"
    "Splits a montage image, one that is composed of several axis aligned\n"
    "rectangular component images separated by a gap, with the gap having\n"
    "either a  single value or a  narrow range of values.  The component\n"
    "images are each output to a separate file in raster scan order from\n"
    "top left to bottom right (approximately).\n"
    "Command line options are:\n"
    "  -h  Help - prints this usage massage.\n"
    "  -v  Be verbose, probably only useful for debugging.\n"
    "  -a  Minimum area for each of the extracted components, value = %d.\n"
    "  -A  Use letters rather than numbers to distinguish component files\n"
    "      (default %s).\n"
    "  -b  Base for output file name, default is the input file name base.\n"
    "  -g  Gap value which must either be a single value or a comma\n"
    "      separated triple of values (red,green,blue), value = %s.\n"
    "  -f  Input image format, default is to use the input file extension.\n"
    "  -F  Output image format, default is same as input.\n"
    "  -t  Tolerance (percentage) for the gap value, value %g.\n"
    "  -w  Width of border added to extracted components, value = %d.\n"
    "Example:\n"
    "  %s -g 255 img.jpg\n"
    "Splits the montage file img.jpg which has a pure white gaps,\n"
    "separating it into component images. The component images are\n"
    "written to files img000001.jpg, img000002.jpg, etc..\n",
    argv[0], minArea, (useAlpha)? "set": "not set", gapStr,
    tol * 100.0, bWidth, argv[0]);
  }
  AlcFree(gapStr);
  return(!ok);
}

/*!
* \return       Difference between column origin of object's bounding boxes.
* \ingroup	BinWlzApp
* \brief        Comparison function for AlgHeapSortIdx().
*               Sorted data will approximately have top left entry first and
*		bottom right last. A small allowance is made for slightly
*		missaligned images.
* \param        data                    Data array.
* \param        idx                     Index array.
* \param        id0                     First index.
* \param        id1                     Second index.
*/
static int      WlzSplitImageSortBBFn(void *data, int *idx, int id0, int id1)
{
  WlzIBox3      bBox0,
                bBox1;
  WlzObject     *obj0,
                *obj1;
  int           cmp0,
  		cmp1,
		cmp = 0;

  obj0 = *((WlzObject **)data + *(idx + id0));
  obj1 = *((WlzObject **)data + *(idx + id1));
  bBox0 = WlzBoundingBox3I(obj0, NULL);
  bBox1 = WlzBoundingBox3I(obj1, NULL);
  cmp0 = bBox0.xMin - bBox1.xMin;
  cmp1 = bBox0.yMin - bBox1.yMin;
  if(abs(cmp0) > 10 * abs(cmp1))
  {
    cmp = cmp0;
  }
  else if(cmp1)
  {
    cmp = cmp1;
  }
  else
  {
    cmp = cmp0;
  }
  return(cmp);
}

/*!
* \return	String encoding the given index.
* \ingroup	BinWlzApp
* \brief	Encodes the given index in a string using A-Z and numbers.
* \param	idxBuf			Buffer for string.
* \param	nComp			Maximum value of the index.
* \param	idx			Given index (always >= 1).
*/
static void	WlzSplitImageIdxToAlpha(char *idxBuf, int maxIdx, int idx)
{
  int		n;
  char		c;

  --idx;
  if(maxIdx <= 26)
  {
    *(idxBuf + 0) = 'A' + idx;
    *(idxBuf + 1) = '\0';
  }
  else
  {
    c = 'A' + idx % 26;
    n = (idx / 26) + 1;
    sprintf(idxBuf, "%c%06d", c, n);
  }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
