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
* \file         binWlzApp/WlzSplitImage.c
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
* \brief	Splits an object that is composed of several components
* 		seperated by background into seperate objects.
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzSplitImage "WlzSplitImage"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzSplitImage WlzSplitImage
\par Name
WlzSplitImage - splits an object that is composed of several components
                seperated by background into seperate objects.
\par Synopsis
\verbatim
WlzSplitImage [-h] [-v] [-b#] [-d] [-f#] [-F#] [-n#] [-s#]
              <path to image file>

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
    <td><b>-b</b></td>
    <td>Border width for each component image object.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Write the component image objects to separate directories
        using MA Editorial office conventions.</td>
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
    <td><b>-n</b></td>
    <td>Number of component images to extract, value 2.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Histogram smoothing parameter, value 5.</td>
  </tr>
</table>
\par Description
WlzSplitImage splits the given image into component images, with a border region
around each one and the components being numbered 1, 2, .... The
image is split using the image values or their modulus of the image
is colour. At least 20% of the image should be background and there
should be a good background/foreground brightness contrast. The output
images are sorted by column origin of the bounding boxes, ie the
component numbrs increase from left to right.
origins of their bounding boxes increase.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzSplitImage.c "WlzSplitImage.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzLabel "WlzLabel(1)"
\ref WlzThreshold "WlzThreshold(1)"
\ref WlzSplitMontage "WlzSplitMontage(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
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

int             main(int argc, char *argv[])
{
  int		idC,
  		len,
  		option,
		fileNum,
		bWidth = 50,
		nComp = 2,
		useDirs = 0,
  		ok = 1,
		usage = 0,
		verbose = 0;
  double	sigma = 5.0;
  char		*dot,
  		*ext,
  		*sep,
		*pathDir,
		*inFile,
		*dirPFx,
  		*inPath = NULL;
  int		*orderTb = NULL;
  WlzEffFormat	inFmt = WLZEFF_FORMAT_NONE,
  		outFmt = WLZEFF_FORMAT_NONE;
  WlzObject	*ppObj = NULL,
  		*inObj = NULL;
  WlzObject	**compObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  struct stat	statBuf;
  char		pathBuf[FILENAME_MAX];
  const char	*errMsg;
  static char	optList[] = "dhvb:f:F:n:s:",
  		defPathDir[] = ".",
		defDirPFx[] = "";
  const double	bgdFrac = 0.20;

  dirPFx = defDirPFx;
  while((usage == 0) && ok &&
        ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
        useDirs = 1;
	break;
      case 'h':
        usage = 1;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'b':
        if(sscanf(optarg, "%d", &bWidth) != 1)
	{
	  usage = 1;
	}
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
      case 'n':
        if(sscanf(optarg, "%d", &nComp) != 1)
	{
	  usage = 1;
	}
	break;
      case 's':
        if(sscanf(optarg, "%lg", &sigma) != 1)
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
    if(sscanf(inFile, "%d%s", &fileNum, dirPFx) != 2)
    {
      dirPFx = defDirPFx;
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
    (void )fprintf(stderr, "read input image ok.\n");
  }
  /* Convert to grey and normalise 0 - 255. */
  if(ok)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      if(WlzGreyTypeFromObj(inObj, NULL) == WLZ_GREY_RGBA)
      {
        ppObj = WlzAssignObject(
	        WlzRGBAToModulus(inObj, &errNum), NULL);
      }
      else
      {
        ppObj = WlzAssignObject(
		WlzCopyObject(inObj, &errNum), NULL);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzGreyNormalise(ppObj);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to normalise object (%s)\n",
		     *argv, errMsg);
    }
  }
  /* Split image. */
  if(ok)
  {
    errNum = WlzSplitObj(inObj, ppObj, bWidth, bgdFrac,
    			 sigma, WLZ_COMPTHRESH_GRADIENT, nComp,
    			 &nComp, &compObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to split image (%s)\n",
		     *argv, errMsg);
    }
  }
  /* Sort the objects into left -> right order. */
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
  /* Create directory structure and write the components to files. */
  if(ok)
  {
    idC = 0;
    while((errNum == WLZ_ERR_NONE) && (idC < nComp))
    {
      if(useDirs)
      {
	/* Create directory name. */
        (void )sprintf(pathBuf, "%s/%s%d", pathDir, dirPFx, idC + 1);
	/* Check if diectory exists and create it if it doesn't. */
	if(stat(pathBuf, &statBuf) == 0)
	{
#ifdef __S_IFDIR
          if((statBuf.st_mode &
	      (__S_IFDIR|__S_IREAD|__S_IWRITE|__S_IEXEC)) == 0)
#else
          if((statBuf.st_mode & (S_IFDIR | S_IRWXU)) == 0)
#endif
	  {
	    errNum = WLZ_ERR_WRITE_EOF;
	  }
	}
	else
	{
	  if(mkdir(pathBuf, S_IRWXU) != 0)
	  {
	    errNum = WLZ_ERR_WRITE_EOF;
	  }
	}
	/* Compute output file path. */
	if(ext)
	{
	  (void )sprintf(pathBuf, "%s/%s%d/%s%d.%s",
	  		 pathDir, dirPFx, idC + 1, inFile, idC + 1, ext);
	}
	else
	{
	  (void )sprintf(pathBuf, "%s/%s%d/%s%d",
	                 pathDir, dirPFx, idC + 1, inFile, idC + 1);
	}
      }
      else
      {
        /* Compute output file path.. */
	if(ext)
	{
	  (void )sprintf(pathBuf, "%s/%s%d.%s",
	  		 pathDir, inFile, idC + 1, ext);
	}
	else
	{
	  (void )sprintf(pathBuf, "%s/%s%d",
	  		 pathDir, inFile, idC + 1);
	}
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
  (void )WlzFreeObj(ppObj);
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
    "Usage: %s [-b#] [-d] [-h] [-f#] [-F#] [-n#] [-s#] [-v]\n"
    "       <path to image file>\n"
    "Splits the given image into component images, with a border region\n"
    "around each one and the components being numbered 1, 2, .... The\n"
    "image is split using the image values or their modulus of the image\n"
    "is colour. At least 20%% of the image should be background and there\n"
    "should be a good background/foreground brightness contrast. The output\n"
    "images are sorted by column origin of the bounding boxes, ie the\n"
    "component numbrs increase from left to right.\n"
    "origins of their bounding boxes increase.\n"
    "Command line options are:\n"
    "  -b  Border width for each component image.\n"
    "  -d  Write component images to separate directories using MA Editorial\n"
    "      office conventions.\n"
    "  -h  Help - prints this usage masseage.\n"
    "  -f  Input image format, default is to use the input file extension.\n"
    "  -F  Output image format, default is same as input.\n"
    "  -n  Number of component images to extract, value %d.\n"
    "  -s  Histogram smoothing parameter, value %g.\n"
    "  -v  Be verbose, probably only useful for debugging.\n",
    argv[0], nComp, sigma);
  }
  return(!ok);
}

/*!
* \return	Difference between column origin of object's bounding boxes.
* \brief	Comparison function for AlgHeapSortIdx().
*		Sorted data will have left-most entry first and right-most
*		last.
* \param	data			Data array.
* \param	idx			Index array.
* \param	id0			First index.
* \param	id1			Second index.
*/
static int	WlzSplitImageSortBBFn(void *data, int *idx, int id0, int id1)
{
  WlzIBox3	bBox0,
  		bBox1;
  WlzObject	*obj0,
  		*obj1;
  int		cmp = 0;

  obj0 = *((WlzObject **)data + *(idx + id0));
  obj1 = *((WlzObject **)data + *(idx + id1));
  bBox0 = WlzBoundingBox3I(obj0, NULL);
  bBox1 = WlzBoundingBox3I(obj1, NULL);
  cmp = bBox0.xMin - bBox1.xMin;
  return(cmp);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
