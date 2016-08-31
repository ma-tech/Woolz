#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshValueTransfer_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshValueTransfer.c
* \author       Bill Hill
* \date         August 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Transfers the values of one CMesh object to another
* 		within the intersection of their (mesh) domains.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshvaluetransfer "WlzCMeshValueTransfer"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshvaluetransfer WlzCMeshValueTransfer
\par Name
WlzCMeshValueTransfer - creates a new object  with the mesh domain of
			the target object and (within the domain intersection)
			the values of the source object.
\par Synopsis
\verbatim
WlzCMeshValueTransfer [-h] [-o<output object>] [-L] [-x]
                      [<source object>] [<target object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object.</td>
  </tr>
  <tr>
    <td><b>-L</b></td>
    <td>Use linear (barycentric) interpolation instead of nearest
        neighbour.</td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>Value for indexed values external to the source mesh
        (default is 0).</td>
  </tr>
</table>
\par Description
Given a pair of 2, 2.5 or 3D conforming meshes,
with both of the same type,
values are transfered from the source to the target mesh
within their intersection.
Outside of the intersection values are set to the value
for external values.
By default files are read from the standard input and written to the
standard output. If both the source and target files are read from the
standard input, the source is read before the target object.
\par Examples
\verbatim
WlzCMeshValueTransfer -x NaN -o out.wlz mesh.wlz supermesh.wlz
\endverbatim
Reads source and target meshes from the files mesh.wlz and
supermesh.wlz respectively then creates a new mesh setting all
values within the intersection to the nearest value in source.wlz.
Values outside the intersection are set to the given external value
NaN. The resulting object is then written to the file out.wlz.
\par File
\ref WlzCMeshValueTransfer.c "binWlz/WlzCMeshValueTransfer.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshValueTransfer "WlzCMeshValueTransfer(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int		option,
  		ok = 1,
		usage = 0;
  WlzObject     *srcObj = NULL,
  		*tgtObj = NULL,
  		*outObj = NULL;
  WlzPixelV	extVal;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  char		*srcObjFileStr,
  		*tgtObjFileStr,
  		*outObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "hLo:x:",
  		fileStrDef[] = "-";

  opterr = 0;
  extVal.v.dbv = 0.0;
  extVal.type = WLZ_GREY_DOUBLE;
  srcObjFileStr = fileStrDef;
  tgtObjFileStr = fileStrDef;
  outObjFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'x':
	if(strcasecmp(optarg, "NaN") == 0)
	{
	  extVal.v.dbv = NAN;
	}
	else if(sscanf(optarg, "%lg", &(extVal.v.dbv)) != 1)
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
    int		cnt;
    
    cnt = argc - optind;
    switch(cnt)
    {
      case 0:
        break;
      case 2:
        srcObjFileStr = *(argv + optind);
        tgtObjFileStr = *(argv + optind + 1);
	break;
      default:
        usage = 1;
	break;
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    int		idx;
    const char	*fileStr[2];
    WlzObject   **obj[2];

    obj[0] = &srcObj;
    obj[1] = &tgtObj;
    fileStr[0] = srcObjFileStr;
    fileStr[1] = tgtObjFileStr;
    for(idx = 0; idx < 2; ++idx)
    {
      FILE	*fP = NULL;

      errNum = WLZ_ERR_FILE_OPEN;
      if((fileStr[idx] == NULL) ||
	  (*(fileStr[idx]) == '\0') ||
	  ((fP = (strcmp(fileStr[idx], "-")?
		 fopen(fileStr[idx], "r"): stdin)) == NULL) ||
	  ((*(obj[idx]) = WlzAssignObject(
		          WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
      }
      if(fP && strcmp(fileStr[idx], "-"))
      {
	(void )fclose(fP);
      }
      if(!ok)
      {
        break;
      }
    }
    if(!ok)
    {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	     "%s: Failed to read %s mesh object from file %s (%s).\n",
	     *argv, (idx)? "target": "source", fileStr[idx], errMsg);
    }
  }
  if(ok)
  {
    outObj = WlzAssignObject(
             WlzCMeshValueTransfer(srcObj, tgtObj, extVal, interp,
	                           &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	     "%s: Failed to transfer values (%s).\n",
	     *argv, errMsg);

    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

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
  (void )WlzFreeObj(srcObj);
  (void )WlzFreeObj(tgtObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-h] [-o<output object>] [-L] [-x<value>]\n"
    "\t\t{<source object>] [<target object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints usage message.\n"
    "  -o  Output object file.\n"
    "  -L  Use linear (barycentric) interpolation instead of nearest\n"
    "      neighbour.\n"
    "  -x  Value for indexed values external to the source mesh (default\n"
    "      is 0).\n"
    "Creates a new object with the mesh domain of the target object and\n"
    "(within the domain intersection) the values of the source object.\n"
    "Given a pair of 2, 2.5 or 3D conforming meshes, with both of the same\n"
    "type, a new mesh object is created using the domain of the target\n"
    "and with the value type of the source object. Values are transfered\n"
    "from the source to the target mesh within their intersection. Outside\n"
    "of the intersection values are set to the value for external values.\n"
    "By default files are read from the standard input and written to the\n"
    "standard output. If both the source and target files are read from the\n"
    "standard input, the source is read before the target object.\n",
    *argv,
    " -x NaN -o out.wlz mesh.wlz supermesh.wlz\n"
    "Reads source and target meshes from the files mesh.wlz and\n"
    "supermesh.wlz respectively then creates a new mesh setting all\n"
    "values within the intersection to the nearest value in source.wlz.\n"
    "Values outside the intersection are set to the given external value\n"
    "NaN. The resulting object is then written to the file out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
