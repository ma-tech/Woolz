#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshTrStrainTensor_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzCMeshTrStrainTensor.c
* \author       Bill Hill
* \date         January 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Creates a conforming mesh with a strain tensor attached
* 		to each element from the given conforming mesh transform.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzcmeshtrstraintensor "WlzCMeshTrStrainTensor"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshtrstraintensor WlzCMeshTrStrainTensor
\par Name
WlzCMeshTrStrainTensor - computes strain tensors for the given conforming
                         mesh transform.
\par Synopsis
\verbatim
WlzCMeshTrStrainTensor [-d] [h] [-i] [-p] [-o<out file>] [-s<x>,<y>,<z>]
                       [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-d</b></td>
    <td>If computing the tensor at points, dither the point locations.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Invert the transform, by default the tensors are computed for the
        transform from the mesh to the displaced mesh.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Make the output object a points object rather than a conforming
        mesh.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Sampling interval. Giving a sampling interval will also result
        in a points output object.</td>
  </tr>
</table>
\par Description
Computes strain tensors for the given conforming mesh transform.
These then output as either a conforming mesh with the tensors
attached to the elements or a points object.
By default all files are read from the standard input and are written
to the standard output.
\par Examples
\verbatim
WlzCMeshTrStrainTensor -o meshten.wlz meshtr.wlz
\endverbatim
Creates a new conforming mesh object (written to meshten.wlz) with the
same mesh domain as the input conforming mesh transform (read from
meshtr.wlz) but with tensor indexed values attached to the output mesh's
elements.
\par File
\ref WlzCMeshTrStrainTensor.c "WlzCMeshTrStrainTensor.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshStrainTensor "WlzCMeshStrainTensor(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>


/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		option,
		dither = 0,
		invert = 0,
  		ok = 1,
  		usage = 0;
  FILE		*fP = NULL;
  char		*inFileStr = NULL,
  		*outFileStr = NULL;
  WlzObjectType	outType = WLZ_CMESH_3D;
  WlzDVertex3   sampleDist;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char   optList[] = "dhipo:s:";
  const char    fileStrDef[] = "-";

  opterr = 0;
  inFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  WLZ_VTX_3_SET(sampleDist, 1.0, 1.0, 1.0);
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
        dither = 1;
	break;
      case 'i':
        invert = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'p':
        outType = WLZ_POINTS;
        break;
      case 's':
        outType = WLZ_POINTS;
        if(sscanf(optarg,
                  "%lg,%lg,%lg",
                  &(sampleDist.vtX),
                  &(sampleDist.vtY),
                  &(sampleDist.vtZ)) < 3)
        {
          usage = 1;
        }
        break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if(optind < argc)
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inFileStr = *(argv + optind);
      }
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if(((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    switch(outType)
    {
      case WLZ_CMESH_3D:
        outObj = WlzCMeshStrainTensor(inObj, invert, &errNum);
	break;
      case WLZ_POINTS:
        outObj = WlzCMeshStrainTensorAtPts(inObj, invert, sampleDist,
					   dither, &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		"%s: Failed to compute tensor mesh object (%s)\n",
		argv[0],
		errMsgStr);
    }
  }
  WlzFreeObj(inObj);
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to write object to file %s\n",
                     *argv, outFileStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
      "Usage: %s [-d] [-h] [-i]  [-p] [-o<out file>] [-s<x>,<y>,<z>]\n"
      "\t\t[<input file>]\n"
      "Computes strain tensors for the given conforming mesh transform.\n"
      "If the output object type is a conforming mesh these tensor values\n"
      "are attached to the output object's elements.\n"
      "By default all files are read from the standard input and are written\n"
      "to the standard output.\n"
      "Version: %s\n"
      "Options are:\n"
      "  -d  If computing the tensor at points, dither the point locations.\n"
      "  -h  Help, prints this usage message.\n"
      "  -i  Invert the transform, by default the tensors are computed for\n"
      "      the transform from the mesh to the displaced mesh.\n"
      "  -p  Make the output object a points object rather than a conforming\n"
      "      mesh.\n"
      "  -s  Sampling interval. Giving a sampling interval will also result\n"
      "      in a points output object.\n"
      "  -o  Output object.\n"
      "%s -o meshten.wlz meshtr.wlz\n"
      "Creates a new conforming mesh object (written to meshten.wlz) with\n"
      "the same mesh domain as the input conforming mesh transform (read\n"
      "from meshtr.wlz) but with tensor indexed values attached to the\n"
      "output mesh's elements.\n",
      argv[0],
      WlzVersion(),
      argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
