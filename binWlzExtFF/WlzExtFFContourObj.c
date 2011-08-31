#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFContourObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzExtFF/WlzExtFFContourObj.c
* \author       Bill Hill
* \date         September 2001
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
* \brief	A program to extract boundary, maximal gradient and isosurface
* 		contours from 2 or 3D Woolz and saving them in VTK format
*		files.
* \ingroup	BinWlzExtFF
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzextffcontourobj "WlzExtFFContourObj"
*/

/*!
\ingroup BinWlzExtFF
\defgroup wlzextffcontourobj WlzExtFFContourObj
\par Name
WlzExtFFContourObj - computes a VTK polydata file from a Woolz object.
\par Synopsis
\verbatim
WlzExtFFContourObj  [-h] [-o<output object>] [-f#] [-F#]
                    [-b] [-g] [-i] [-l] [-L] [-N] [-U]
		    [-p#] [-s#] [-n#] [-v#] [-w#] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output VTK polydata file name.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Compute object boundary contours.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Input file format.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Compute maximal gradient contours.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Compute iso-value contours.</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Flip orientation (normals will be reversed).</td>
  </tr>
  <tr> 
    <td><b>-F</b></td>
    <td>Output file format.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Use geometry filter.</td>
  </tr>
  <tr> 
    <td><b>-N</b></td>
    <td>Allow non manifold vertices to be filtered.</td>
  </tr>
  <tr> 
    <td><b>-U</b></td>
    <td>Use unit voxel size.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Generate normals (if possible).</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Geometry filter itterations.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Geometry filter low band value.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Geometry filter stop band value.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Contour iso-value or minimum gradient value.</td>
  </tr>
  <tr> 
    <td><b>-w</b></td>
    <td>Contour (Deriche) gradient operator width.</td>
  </tr>
</table>
\par Description
WlzExtFFContourObj computes a contour object from the given Woolz object and saves it
using the VTK ascii polydata format.
See the documentation for WlzContourGeomFilter(1) for an explaination
of the filter parameters.
The input object is read from stdin and output data are written
to stdout unless filenames are given.
\par Examples
\verbatim
WlzExtFFContourObj -o out.vtk -i -v 100.0 in.wlz
\endverbatim
The input Woolz object is read from in.wlz, and the iso-value
(iso-value = 100.0) contour is written in VTK ascii polydata format
to out.vtk.
\par File
\ref WlzFacts.c "WlzFacts.c"
\par See Also
\ref BinWlzExtFF "WlzIntro(1)"
\ref wlzextffconvert "WlzExtFFConvert(1)"
\ref wlzcontourobj "WlzContourObj(1)"
\ref wlzcontourgeomfilter "WlzContourGeomFilter(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int           option,
  		ok = 1,
		usage = 0,
		flip = 0,
		nrm = 0,
		nItr = 10,
		nonMan = 0,
		unitVoxelSz = 0,
		filterGeom = 0;
  double	lambda,
  		mu,
  		ctrVal = 100,
  		ctrWth = 1.0,
		filterPB = 0.1,
		filterSB = 1.1;
  FILE		*fP = NULL;
  char		*fStr,
		*inObjFileStr,
  		*outObjFileStr;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  WlzDomain	ctrDom;
  WlzValues	dumVal;
  WlzEffFormat  inFmt = WLZEFF_FORMAT_NONE,
                outFmt = WLZEFF_FORMAT_NONE;
  WlzContourMethod ctrMtd = WLZ_CONTOUR_MTD_ISO;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const double	filterDPB = 0.25,
  		filterDSB = 0.10;
  const char	*errMsgStr;
  static char	optList[] = "bfghilmFLNUo:p:s:n:v:w:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  ctrDom.core = NULL;
  dumVal.core = NULL;
  outObjFileStr = (char *)outFileStrDef;
  inObjFileStr = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        ctrMtd = WLZ_CONTOUR_MTD_BND;
	break;

      case 'g':
        ctrMtd = WLZ_CONTOUR_MTD_GRD;
	break;
      case 'f':
        if((inFmt = WlzEffStringExtToFormat(optarg)) == 0)
        {
          usage = 1;
          ok = 0;
        }
        break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'i':
        ctrMtd = WLZ_CONTOUR_MTD_ISO;
	break;
      case 'l':
        flip = 1;
	break;
      case 'm':
        nrm = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'v':
        if(sscanf(optarg, "%lg", &ctrVal) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'n':
        if(sscanf(optarg, "%d", &nItr) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'p':
        if(sscanf(optarg, "%lg", &filterPB) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
        if(sscanf(optarg, "%lg", &filterSB) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'F':
        if((outFmt = WlzEffStringExtToFormat(optarg)) == 0)
        {
          usage = 1;
          ok = 0;
        }
        break;
      case 'L':
        filterGeom = 1;
	break;
      case 'N':
        nonMan = 1;
	break;
      case 'U':
        unitVoxelSz = 1;
	break;
      case 'w':
        if(sscanf(optarg, "%lg", &ctrWth) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
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
  }
  if(ok)
  {
    /* Try and determine file formats from extensions if not known already. */
    if(inFmt == WLZEFF_FORMAT_NONE)
    {
      if((inFmt = WlzEffStringFormatFromFileName(inObjFileStr)) == 0)
      {
        usage = 1;
        ok = 0;
      }
    }
    if(outFmt == WLZEFF_FORMAT_NONE)
    {
      if((outFmt = WlzEffStringFormatFromFileName(outObjFileStr)) == 0)
      {
        usage = 1;
        ok = 0;
      }
    }
  }
  if(ok)
  {
    /* If the output format is TIFF then a file name must be
       given */
    if((outFmt == WLZEFF_FORMAT_TIFF) && !strcmp(outObjFileStr, "-"))
    {
      usage = 1;
      ok = 0;
    }
  }
  if(ok)
  {
    if(strcmp(inObjFileStr, "-") == 0)
    {
      fP = stdin;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = inObjFileStr;
    }
    errNum = WLZ_ERR_READ_EOF;
    if((inObj = WlzAssignObject(WlzEffReadObj(fP, fStr, inFmt, 0,
                                              &errNum), NULL)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok && unitVoxelSz)
  {
    if(inObj && (inObj->type == WLZ_3D_DOMAINOBJ) &&
       inObj->domain.core && (inObj->domain.core->type = WLZ_2D_DOMAINOBJ))
    {
      inObj->domain.p->voxel_size[0] = 1.0;
      inObj->domain.p->voxel_size[1] = 1.0;
      inObj->domain.p->voxel_size[2] = 1.0;
    } 
  }
  if(ok)
  {
    ctrDom.ctr = WlzContourObj(inObj, ctrMtd, ctrVal, ctrWth, nrm, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to compute contour (%s).\n",
      	     	     argv[0], errMsgStr);
    }
  }
  if(ok && filterGeom)
  {
    errNum = WlzGMFilterGeomLPParam(&lambda, &mu, &nItr,
    				    filterPB, filterSB, filterDPB, filterDSB);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzGMFilterGeomLPLM(ctrDom.ctr->model, lambda, mu,
      				   nItr, nonMan);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      "%s: Failed to filter output woolz object's geometry (%s).\n",
      		     argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_CONTOUR, ctrDom, dumVal, NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to create output woolz object (%s).\n",
		     argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    if(flip && ctrDom.core && ctrDom.ctr->model)
    {
      if((errNum = WlzGMFilterFlipOrient(ctrDom.ctr->model)) != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		       "%s: Failed to flip orientation (%s).\n",
		       argv[0], errMsgStr);

      }
    }
  }
  if(ok)
  {
    if(strcmp(outObjFileStr, "-") == 0)
    {
      fP = stdout;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = outObjFileStr;
    }
  }
  if(ok)
  {
    errNum = WlzEffWriteObj(fP, fStr, outObj, outFmt);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to write output object (%s).\n",
      		     argv[0], errMsgStr);
    }
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(ctrDom.ctr)
  {
    (void )WlzFreeContour(ctrDom.ctr);
  }
  if(fP && strcmp(outObjFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
    char *fmtStr = NULL;

    fmtStr = WlzEffFormatTable(2, 50, 10, NULL);

    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s%s%s",
    *argv,
    " [-h] [-o<output object>] [-b] [-g] [-i] [-l]\n"
    "        [-L] [-N] [-U] [-p#] [-s#] [-n#] [-v#] [-w#]\n"
    "        [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output VTK polydata file name.\n"
    "  -b  Compute object boundary contours.\n"
    "  -f  Input file format.\n"
    "  -g  Compute maximal gradient contours.\n"
    "  -i  Compute iso-value contours.\n"
    "  -l  Flip orientation (normals will be reversed).\n"
    "  -F  Ouput file format.\n"
    "  -L  Use geometry filter.\n"
    "  -N  Allow non manifold vertices to be filtered.\n"
    "  -U  Use unit voxel size.\n"
    "  -m  Generate normals (if possible).\n"
    "  -n  Geometry filter itterations.\n"
    "  -p  Geometry filter low band value.\n"
    "  -s  Geometry filter stop band value.\n"
    "  -v  Contour iso-value or minimum gradient.\n"
    "  -w  Contour (Deriche) gradient operator width.\n",
    "The known file formats are:\n"
    "  Description                                       Extension\n"
    "  ***********                                       *********\n",
    fmtStr,
    "Computes a contour object from the given Woolz object and saves it\n"
    "using the VTK ascii polydata format.\n"
    "The input object is read from stdin and output data are written\n"
    "to stdout unless filenames are given.\n",
    *argv,
    " -o out.vtk -i -v 100.0 in.wlz\n"
    "The input Woolz object is read from in.wlz, and the iso-value\n"
    "(iso-value = 100.0) contour is written in VTK ascii polydata format\n"
    "to out.vtk.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
