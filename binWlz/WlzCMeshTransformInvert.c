#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshTransformInvert_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzCMeshTransformInvert.c
* \author       Bill Hill
* \date         April 2011
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2011 Medical research Council, UK.
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
* \brief	Inverts a constrained mesh transform.
* \ingroup	BinWlz
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshtransforminvert WlzCMeshTransformInvert
\par Name
WlzCMeshTransformInvert - inverts a constrained mesh transform.
\par Synopsis
\verbatim
WlzCMeshTransformInvert [-h] [-o<output transform>] [<input transform>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-0</b></td>
    <td>Output transform file.</td>
  </tr>
</table>
\par Description
Reads a constrained mesh transform object, inverts it and then
writes out the inverted transform.
Unless given input and output files are read from the standard input
and written to the standard output.
\par Example
\verbatim
WlzCMeshTransformInvert -o out.wlz in.wlz
\endverbatim
Reads a constrained mesh transform read from in.wlz and writes the
inverted transform to out.wlz
\par File
\ref WlzCMeshTransformInvert.c "WlzCMeshTransformInvert.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */


int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0;
  FILE		*fP = NULL;
  char          *inFileStr = NULL,
                *outFileStr = NULL;
  const char	*errMsgStr;
  WlzObject	*iObj = NULL,
  		*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "ho:",
  		inFileStrDef[] = "-",
                outFileStrDef[] = "-";

  opterr = 0;
  inFileStr = (char *)inFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
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
      inFileStr = argv[optind];
    }
  }
  ok = usage == 0;
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if(((fP = (strcmp(inFileStr, "-")?
	      fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((iObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(iObj->type)
      {
        case WLZ_CMESH_2D:  /* FALLTHROUGH */
        case WLZ_CMESH_2D5: /* FALLTHROUGH */
	case WLZ_CMESH_3D:
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if((fP != NULL) && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    oObj = WlzCMeshTransformInvert(iObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to invert transform (%s).\n",
                     argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"):
	      stdout)) == NULL) ||
	((errNum = WlzWriteObj(fP, oObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	  "%s: Failed to write CMesh object (%s).\n",
	  *argv, errMsgStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }

  }
  (void )WlzFreeObj(iObj);
  (void )WlzFreeObj(oObj);
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s [-h] [-o<output transform>] [<input transform>]\n"
      "Reads a constrained mesh transform object, inverts it and then\n"
      "writes out the inverted transform.\n"
      "Unless given input and output files are read from the standard input\n"
      "and written to the standard output.\n"
      "Example:\n"
      "  %s -o out.wlz in.wlz\n"
      "Reads a constrained mesh transform read from in.wlz and writes the\n"
      "inverted transform to out.wlz\n",
      argv[0],argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
