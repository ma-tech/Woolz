#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFConstruct3D_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzExtFF/WlzExtFFConstruct3D.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Constructs a 3D image from 2D images, assuming perfect
* 		alignment.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzconstruct3d "WlzConstruct3D"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzextffconstruct3d WlzExtFFConstruct3D
\par Name
WlzExtFFConstruct3D - Constructs a 3D image from 2D images, assuming perfect
alignment.
\par Synopsis
\verbatim
WlzExtFFConstruct3D [-h] [-f] [-o<output file>] [-p #] [-s #,#,#]
                    [<file>|<input file list>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-o</b></td>
    <td>Output file name, default to standard out.</td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>Input file is a list of image files</td>
  </tr>
  <tr>
    <td><b>-p</b></td>
    <td>Coordinate of the first plane</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Voxel size (x, y, z).</td>
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
Read a single 2D image from each image file and writes a 3D image.

\par Description
Given a set of 2D image files,
these are read in turn and used to build a 3D image.
Alignment of the 2D images is assumed to be perfect and is not addressed.
The 2D image files can either be given on the command line or in a file.
If a file is used to give the 2D image file names then there should
be one image file name per line,
although lines which start with a '#' character will be treated as comments
and ignored. 
File formats are specified using the file 'extension'
(.bmp, .jpg, .wlz, etc) in the file name.
The dash character ('-') can be used to specify the standard input and
standard output,
in which case the Woolz object format is used.
The word 'NULL' may be used to specify an empty section.

\par Examples
\verbatim
WlzExtFFConstruct3D img-00.bmp img-01.bmp img-02.bmp img-03.bmp >out.wlz
\endverbatim
Constructs a 3D Woolz object from the 2D images in img-0X.bmp.

\verbatim
WlzExtFFConstruct3D -f -o out.am  img-list.txt
\endverbatim
Constructs a 3D Amira lattice file from the 2D images in given (one per line)
in the file img-list.txt.

\par File
\ref WlzExtFFConstruct3D.c "WlzExtFFConstruct3D.c"
\par See Also
\ref wlzexplode "WlzExplode(1)"
\ref wlzmakeempty "WlzMakeEmpty(1)"
\ref WlzConstruct3DObjFromFile "WlzConstruct3DObjFromFile(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

#define WLZEFFCON_FILE_STEP	1024
#define WLZEFFCON_BUF_LEN	256

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static WlzErrorNum 		WlzEFFConAddImg(
				  WlzObject **objP,
				  const char *fStr);

int		main(int argc, char *argv[])
{
  int		idx,
  		ok = 1,
  		option,
  		usage = 0,
		fileList = 0,
		endOfFileList,
		nObj,
		nObjMax,
		nFiles,
		plane1 = 0;
  FILE		*fP;
  char		*fStr,
  		*fLstStr = NULL,
  		*outFileStr = NULL;
  const char	*errMsgStr;
  WlzEffFormat	outFmt = WLZEFF_FORMAT_WLZ;
  WlzDVertex3	voxSz;
  WlzObject	*obj = NULL;
  WlzObject	**objs = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		fStrBuf[WLZEFFCON_BUF_LEN];
  static char   optList[] = "fho:p:s:";
  const char    outFileStrDef[] = "-";

  opterr = 0;
  fStr = "";
  voxSz.vtX = voxSz.vtY = voxSz.vtZ = 1.0;
  outFileStr = (char *)outFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'f':
        fileList = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'p':
	if(sscanf(optarg, "%d", &plane1) != 1)
	{
	  usage = 1;
	}
        break;
      case 's':
	if(sscanf(optarg, "%lg,%lg,%lg",
	          &(voxSz.vtX), &(voxSz.vtY), &(voxSz.vtZ)) != 3)
	{
	  usage = 1;
	}
        break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok)
  {
    if((outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
  }
  if(ok)
  {
    if(((nFiles = argc - optind) <= 0) ||
       ((fileList != 0) && (nFiles != 1)))
    {
      ok = 0;
      usage = 1;
    }
  }
  if(ok)
  {
    nObj = 0;
    nObjMax = (fileList)? WLZEFFCON_FILE_STEP: nFiles;
    if((objs = (WlzObject **)AlcMalloc(nObjMax * sizeof(WlzObject *))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to allocate object list.\n",
		     *argv);
    }
  }
  /* Read the 2D image files. */
  if(ok)
  {
    if(fileList == 0)
    {
      idx = 0;
      while((errNum == WLZ_ERR_NONE) && (idx < nFiles))
      {
	fStr = *(argv + optind + idx);
	errNum = WlzEFFConAddImg(objs + idx, fStr);
	++idx;
      }
      nObj = idx;
    }
    else /* fileList != 0 */
    {
      fLstStr = *(argv + optind);
      if((fP = (strcmp(fLstStr, "-")?
	        fopen(fLstStr, "r"): stdin)) == NULL)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to open file list in file %s.\n",
		       *argv, fStr);
      }
      else
      {
        idx = 0;
	nObj = 0;
	endOfFileList = 0;
	while((errNum == WLZ_ERR_NONE) &&
	      (fgets(fStrBuf, WLZEFFCON_BUF_LEN, fP) != NULL))
	{
	  if((idx + 1) > nObjMax)
	  {
	    nObjMax += WLZEFFCON_FILE_STEP;
	    if((objs = (WlzObject **)
	               AlcRealloc(objs,
		                  nObjMax * sizeof(WlzObject *))) == NULL)
            {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(((fStr = strtok(fStrBuf, " \t\n")) != NULL) && (*fStr != '#'))
	    {
	      errNum = WlzEFFConAddImg(objs + idx, fStr);
	      ++idx;
	    }
	  }
	}
	nObj = idx;
      }
      if(fP && (strcmp(fLstStr, "-")))
      {
        (void )fclose(fP);
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read images at image %d\n",
		     *argv, idx);
    }
  }
  /* Construct the 3D object. */
  if(ok)
  {
    obj = WlzConstruct3DObjFromObj(nObj, objs, plane1,
                                   voxSz.vtX, voxSz.vtY, voxSz.vtZ,
				   &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to construct 3D object (%s).\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  /* Write out the 3D object. */
  if(ok)
  {
    if(strcmp(outFileStr, "-") == 0)
    {
      fP = stdout;
      outFileStr = NULL;
      outFmt = WLZEFF_FORMAT_WLZ;
    }
    else
    {
      fP = NULL;
      outFmt = WlzEffStringFormatFromFileName(outFileStr);
    }
    if(outFmt == WLZEFF_FORMAT_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Unrecognized file format for output file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzEffWriteObj(NULL, outFileStr, obj, outFmt);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  if(obj)
  {
    (void )WlzFreeObj(obj);
  }
  if(objs)
  {
    for(idx = 0; idx < nObj; ++idx)
    {
      (void )WlzFreeObj(*(objs + idx));
    }
    AlcFree(objs);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-f] [-o<output file>] [-p #] [-s #,#,#]\n"
    "                           [<file>|<input file list>]\n"
    "  -h  Output this usage message.\n"
    "  -f  Input file is a list of image files.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -p  Coordinate of the first plane.\n"
    "  -s  Voxel size (x,y,z).\n"
    "Given a set of 2D image files, these are read in turn and used\n"
    "to build a 3D image. The 2D image files can either be given on\n"
    "the command line or in a file.\n"
    "Alignment of the 2D images is assumed to be perfect and is not\n"
    "addressed.\n"
    "If a file is used to give the 2D image file names then there\n"
    "should be one image file name per line, although lines which\n"
    "start with a '#' character will be treated as comments and\n"
    "ignored.\n"
    "File formats are specified using the file 'extension' (.bmp,\n"
    ".jpg, .wlz, etc) in the file name.  The dash character ('-')\n"
    "can be used to specify the standard input and standard output,\n"
    "in which case the Woolz object format is used."
    "The word 'NULL' may be used to specify an empty section.\n"
    "Example:\n"
    "%s img00.bmp img01.bmp NULL img03.bmp img04.bmp >out.wlz\n"
    "Constructs a 3D woolz object from the 2D BMP image files\n"
    "img0X.bmp, with an empty 3rd section.\n",
    argv[0], argv[0]);
  }
  return(!ok);
}

static WlzErrorNum WlzEFFConAddImg(WlzObject **objP, const char *fStr)
{
  FILE		*fP;
  WlzEffFormat	inFmt = WLZEFF_FORMAT_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(strcmp(fStr, "NULL") == 0)
  {
    *objP = NULL;
  }
  else
  {
    if(strcmp(fStr, "-") == 0)
    {
      fP = stdin;
      fStr = NULL;
      inFmt = WLZEFF_FORMAT_WLZ;
    }
    else
    {
      fP = NULL;
      inFmt = WlzEffStringFormatFromFileName(fStr);
    }
    if(inFmt == WLZEFF_FORMAT_NONE)
    {
      errNum = WLZ_ERR_FILE_OPEN;
    }
    else
    {
      *objP = WlzAssignObject(
              WlzEffReadObj(fP, fStr, inFmt, 0, &errNum), NULL);
    }
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
