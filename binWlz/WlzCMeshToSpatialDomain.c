#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzMeshGen_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzCMeshToSpatialDomain.c
* \author       Bill Hill
* \date         January 2009
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
* \brief        Constructs a 2D or 3D spatial domain without any
* 		values, but which corresponds to the given conforming
* 		mesh.
* \ingroup	Wlz
* \todo         -
* \bug          None known.
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
WlzCMeshToSpatialDomain [-h] [-o<out obj file>] [<input mesh file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
</table>
\par Description
Constructs a 2D or 3D spatial domain object (without values) that corresponds
to the given conforming mesh. The domain of the output object covers the input
mesh.
\par Examples
\verbatim
WlzCMeshToSpatialDomain -o out.wlz mesh.wlz
\endverbatim
Creates a new spatial domain object which covers the given mesh read from
the file mesh.wlz. The new domain is written to the file out.wlz.
\par File
\ref WlzCMeshToSpatialDomain.c "WlzDistanceTransform.c"
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
  		usage = 0;
  FILE		*fP = NULL;
  char		*inObjFileStr = NULL,
  		*outObjFileStr = NULL;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzCMeshTransform *mTr;
  static char   optList[] = "ho:";
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
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
    mTr = WlzMakeCMeshTransform(WLZ_TRANSFORM_2D_CMESH, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	       "%s: Failed to allocate mesh transform (%s),\n",
	       argv[0],
	       errMsgStr);
    }
  }
  if(ok)
  {
    switch(inObj->type)
    {
      case WLZ_CMESH_2D:
	mTr->mesh.m2 = inObj->domain.cm2;
	break;
      case WLZ_CMESH_3D:
	mTr->mesh.m3 = inObj->domain.cm3;
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
    outObj = WlzCMeshToDomObj(mTr, 0, &errNum);
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
  WlzFreeObj(outObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    fprintf(stderr,
      "Usage: %s [-b] [-h] [-o<out obj file>] [-r<ref obj file>]\n"
      "                        [-s #,#,#] [<input mesh file>]\n"
      "Constructs a 2D or 3D spatial domain object without values which\n"
      "cvovers the given conforming mesh.\n"
      "Options are:\n"
      "  -h  Help, prints this usage message.\n"
      "  -o  Output object.\n",
      argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
