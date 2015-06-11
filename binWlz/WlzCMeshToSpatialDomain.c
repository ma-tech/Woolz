#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshToSpatialDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshToSpatialDomain.c
* \author       Bill Hill
* \date         January 2009
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
* \brief        Constructs a 2D or 3D spatial domain without any
* 		values, but which corresponds to the given conforming
* 		mesh.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzcmeshtospatialdomain "WlzCMeshToSpatialDomain"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshtospatialdomain WlzCMeshToSpatialDomain
\par Name
WlzCMeshToSpatialDomain - computes a new spatial domain corresponding to the
			  conforming mesh.
\par Synopsis
\verbatim
WlzCMeshToSpatialDomain [-h] [-i] [-L] [-o<out obj file>] [-s#]
                        [<input mesh file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Set values using first indexed value.</td>
  </tr>
  <tr>
    <td><b>-L</b></td>
    <td>Use kriging rather than barycentric interpolation for values./td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Scale factor from mesh to spatial domain
        (can not be used if setting values).</td>
  </tr>
</table>
\par Description
Constructs a 2D or 3D spatial domain object that corresponds
to the given conforming mesh. The domain of the output object covers the input
mesh.
\par Examples
\verbatim
WlzCMeshToSpatialDomain -o out.wlz mesh.wlz
\endverbatim
Creates a new spatial domain object which covers the given mesh read from
the file mesh.wlz. The new domain is written to the file out.wlz.
\par File
\ref WlzCMeshToSpatialDomain.c "WlzCMeshToSpatialDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshDistance2D "WlzCMeshDistance2D(3)"
\ref WlzCMeshDistance3D "WlzCMeshDistance3D(3)"
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
  		ok = 1,
		interp = 0,
  		usage = 0,
		setVal = 0;
  double	scale = 1.0;
  FILE		*fP = NULL;
  char		*inObjFileStr = NULL,
  		*outObjFileStr = NULL;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  static char   optList[] = "hiLo:s:";
  const double	eps = 0.000001;
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'i':
        setVal = 1;
	break;
      case 'L':
        interp = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 's':
        if((sscanf(optarg, "%lg", &scale) != 1) || (fabs(scale) < eps))
	{
	  usage = 1;
	}
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if((fabs(scale - 1.0) > eps) && (setVal != 0))
    {
      usage = 1;
    }
  }
  ok = usage == 0;
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
    if(((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read mesh object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    switch(inObj->type)
    {
      case WLZ_CMESH_2D: /* FALLTROUGH */
      case WLZ_CMESH_3D:
	break;
      default:
        ok = 0;
	errNum = WLZ_ERR_OBJECT_TYPE;
	(void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		 "%s: Invalid mesh object, must be WLZ_CMESH_[23]D (%s),\n",
		 argv[0],
		 errMsgStr);
	break;
    }
  }
  if(ok)
  {
    outObj = WlzAssignObject(
             WlzCMeshToDomObj(inObj, 0, scale, &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		"%s: Failed to create spatial domain object (%s)\n",
		argv[0],
		errMsgStr);
    }
  }
  if(ok && setVal)
  {
    WlzObject *tObj;

    tObj = WlzAssignObject(
           WlzCMeshToDomObjValues(outObj, inObj,
	                          (interp)? WLZ_INTERPOLATION_LINEAR:
				            WLZ_INTERPOLATION_NEAREST,
			          &errNum), NULL);
    WlzFreeObj(outObj);
    outObj = tObj;
  }
  WlzFreeObj(inObj);
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to write output object to file %s\n",
                     *argv, outObjFileStr);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
      "Usage: %s [-h] [-o<out obj file>] [-i] [-s #]\n"
      "                               [<input mesh file>]\n"
      "Constructs a 2D or 3D spatial domain object which covers the given\n"
      "conforming mesh.\n"
      "Version: %s\n"
      "Options are:\n"
      "  -h  Help, prints this usage message.\n"
      "  -o  Output object.\n"
      "  -i  Set values using first indexed element value.\n"
      "  -L  Use kriging rather than barycentric interpolation for values.\n"
      "  -s  Additional scale factor from the mesh to the spatial domain\n"
      "      (can not be used if setting values).\n",
      argv[0],
      WlzVersion());

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
