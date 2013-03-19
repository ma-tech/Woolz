#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshSurfaceMap_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshSurfaceMap.c
* \author       Bill Hill
* \date         May 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
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
* \brief	Computes a conforming mesh transform which maps a surface
* 		in 3D to a plane.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshsurfacemap "WlzCMeshSurfaceMap"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshsurfacemap WlzCMeshSurfaceMap
\par Name
WlzCMeshSurfaceMap - computes a conforming mesh transform which maps a
                     surface in 3D to a plane.
\par Synopsis
\verbatim
WlzCMeshSurfaceMap [-2] [-h] [-i#] [-n #] [-o<output object>] [-p<points file>]
                   [-r#] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-2</b></td>
    <td>Output a 2D mesh rather than a 3D mesh with displacements.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>Maximum number of itterations, 0 implies one itteration and
        mapping will preserve angles.</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Weight for preserving angles, range 0.0 - 1.0.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-p</b></td>
    <td>Displacements file.</td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Weight for preserving areas, range 0.0 - 1.0.</td>
  </tr>
</table>
\par Description
Computes a conforming mesh transform which maps the surface defined by
the 2D5 conforming mesh onto another surface defined by the mesh with
displacements. The given displacements are used to compute the least
squares best conformal (minimum angular distortion) mapping. All files
are read from the standard input and written to the standard output
unless filenames are given.
The displacements must be specified as a list of lines with each line
having six double precission numbers that are white space seperated.
In each of these lines the values must have the format:
\verbatim
  Sx Sy Sz Dx Dy Dz
\endverbatim
where S and D signify source and displacement.
\par Examples
\verbatim
WlzCMeshSurfaceMap -o out.wlz -p points.tie in.wlz
\endverbatim
The displacements in the file points.tie are used to compute a least
squares conformal mapping which is output as a mesh with the same
domain as the input (read from in.wlz), but with the computed
displacements. The output mesh is written to the file out.wlz
\par File
\ref WlzCMeshSurfaceMap.c "WlzCMeshSurfaceMap.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshFromContour "WlzCMeshFromContour(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

#define IN_RECORD_MAX   (1024)

static WlzErrorNum		WlzCMSMGetVertices3D(
				  int *dstNV,
				  WlzDVertex3 **dstSrcV,
				  WlzDVertex3 **dstDstV,
				  FILE *fP);

extern int      getopt(int argc, char * const *argv, const char *optstring);
extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int		option,
		maxItr = 0,
		flg2D = 0,
		nV = 0,
  		ok = 1,
		usage = 0;
  double	angleW = 1.0,
  		areaW = 0.0;
  WlzDVertex3	*srcV = NULL,
  		*dstV = NULL;
  WlzObject     *inObj = NULL,
		*mapObj = NULL,
  		*outObj = NULL;
  FILE		*fP = NULL;
  char		*ptsFileStr,
  		*inObjFileStr,
  		*outObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "2hi:n:o:p:r:",
  		fileStrDef[] = "-";

  opterr = 0;
  ptsFileStr = fileStrDef;
  inObjFileStr = fileStrDef;
  outObjFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
        flg2D = 1;
	break;
      case 'i':
        if(sscanf(optarg, "%d", &maxItr) != 1)
	{
	  usage = 1;
	}
	break;
      case 'n':
        if(sscanf(optarg, "%lg", &angleW) != 1)
	{
	  usage = 1;
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'p':
        ptsFileStr = optarg;
	break;
      case 'r':
        if(sscanf(optarg, "%lg", &areaW) != 1)
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
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inObjFileStr = *(argv + optind);
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    if((inObjFileStr == NULL) ||
	(*inObjFileStr == '\0') ||
	((fP = (strcmp(inObjFileStr, "-")?
	       fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
    }
    if(fP)
    {
      if(strcmp(inObjFileStr, "-"))
      {
	(void )fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    if(inObj->type != WLZ_CMESH_2D5)
    {
      ok = 0;
      errNum = WLZ_ERR_OBJECT_TYPE;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Input object must be of type WLZ_CMESH_2D5 (%s).\n",
		     *argv, errMsg);
    }
  }
  /* Read displacements. */
  if(ok)
  {
    if((fP = (strcmp(ptsFileStr, "-")?
             fopen(ptsFileStr, "r"): stdin)) == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_READ_EOF;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s: Failed to open points file %s (%s)\n",
      *argv, ptsFileStr, errMsg);
    }
    else
    {
      if((errNum = WlzCMSMGetVertices3D(&nV, &srcV, &dstV,
                                        fP)) != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	"%s: Failed to read points from the file  %s (%s)\n",
	*argv, ptsFileStr, errMsg);
      }
    }
    if(fP && strcmp(ptsFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  /* Compute mapping transform. */
  if(ok)
  {
    mapObj = WlzCMeshCompSurfMap(inObj, nV, dstV, nV, srcV, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to compute surface mapping (%s)\n",
		     *argv, errMsg);

    }
  }
  AlcFree(srcV);
  AlcFree(dstV);
  if(ok)
  {
    if(flg2D)
    {
      outObj = WlzCMeshExtract2D(mapObj, 1, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to make 2D mesh object (%s)\n",
		       *argv, errMsg);
      }
    }
    else
    {
      outObj = mapObj;
      mapObj = NULL;
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(mapObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-2] [-h] [-i#] [-n #] [-o<output object>]\n"
    "                          [-p<points file>] [-r#] [<input object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -2  Output a 2D mesh rather than a 3D mesh with displacements.\n"
    "  -h  Help, prints usage message.\n"
    "  -i  Maximum number of itterations, 0 implies one itteration and\n"
    "      mapping will preserve angles.\n"
    "  -n  Weight for preserving angles, range 0.0 - 1.0.\n"
    "  -p  Displacements file.\n"
    "  -o  Output object file.\n"
    "  -r  Weight for preserving areas, range 0.0 - 1.0.\n"
    "Computes a conforming mesh transform which maps the surface defined by\n"
    "the 2D5 conforming mesh onto another surface defined by the mesh with\n"
    "displacements. The given displacements are used to compute the least\n"
    "squares best conformal (minimum angular distortion) mapping. All files\n"
    "are read from the standard input and written to the standard output\n"
    "unless filenames are given.\n"
    "The displacements must be specified as a list of lines with each line\n"
    "having six double precission numbers that are white space seperated.\n"
    "In each of these lines the values must have the format:\n"
    "  Sx Sy Sz Dx Dy Dz\n"
    "where S and D signify source and displacement.\n",
    *argv,
    " -o out.wlz -p points.tie in.wlz\n"
    "The displacements in the file points.tie are used to compute a least\n"
    "squares conformal mapping which is output as a mesh with the same\n"
    "domain as the input (read from in.wlz), but with the computed\n"
    "displacements. The output mesh is written to the file out.wlz\n");
  }
  return(!ok);
}

static WlzErrorNum WlzCMSMGetVertices3D(int *dstNV,
				WlzDVertex3 **dstSrcV, WlzDVertex3 **dstDstV,
				FILE *fP)
{
  size_t	sX,
  		sY;
  double	**ary = NULL;
  WlzDVertex3	*srcV = NULL,
  		*dstV = NULL;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  alcErr = AlcDouble2ReadAsci(fP, &ary, &sY, &sX);
  switch(alcErr)
  {
    case ALC_ER_NONE:
      errNum = WLZ_ERR_NONE;
      break;
    case ALC_ER_ALLOC:
      errNum = WLZ_ERR_MEM_ALLOC;
      break;
    default:
      errNum = WLZ_ERR_READ_INCOMPLETE;
      break;
  }
  if((errNum == WLZ_ERR_NONE) && ((sX != 6) || (sY < 1)))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((srcV = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * sY)) == NULL) ||
       ((dstV = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * sY)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx;

    for(idx = 0; idx < sY; ++idx)
    {
      srcV[idx].vtX = ary[idx][0];
      srcV[idx].vtY = ary[idx][1];
      srcV[idx].vtZ = ary[idx][2];
      dstV[idx].vtX = ary[idx][3] + ary[idx][0];
      dstV[idx].vtY = ary[idx][4] + ary[idx][1];
      dstV[idx].vtZ = ary[idx][5] + ary[idx][2];
    }
  }
  AlcDouble2Free(ary);
  if(errNum == WLZ_ERR_NONE)
  {
    *dstNV = sY;
    *dstSrcV = srcV;
    *dstDstV = dstV;
  }
  else
  {
    AlcFree(srcV);
    AlcFree(dstV);
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
