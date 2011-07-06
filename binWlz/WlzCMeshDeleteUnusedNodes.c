#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshDeleteUnusedNodess_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzCMeshDeleteUnusedNodes.c
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
* \brief	Deletes unused nodes, ie those with no edge use
* 		from the given mesh.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzcmeshdeleteunusednodes "WlzCMeshDeleteUnusedNodes"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshdeleteunusednodes WlzCMeshDeleteUnusedNodes
\par Name
WlzCMeshDeleteUnusedNodes - deletes unused nodes from the given mesh.
\par Synopsis
\verbatim
WlzCMeshDeleteUnusedNodes [-h] [-o<output file>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name, default standard output.</td>
  </tr>
</table>
\par Description
Reads a conforming mesh object and deletes any unused nodes from it,
ie those with no edge uses. By default the input object is read from
the standard input and the output object is written to the standard
output.
\par File
\ref WlzCMeshDeleteUnusedNodes.c "WlzCMeshDeleteUnusedNodes.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshDelUnusedNodes2D	WlzCMeshDelUnusedNodes2D(3)
\ref WlzCMeshDelUnusedNodes2D5	WlzCMeshDelUnusedNodes2D5(3)
\ref WlzCMeshDelUnusedNodes3D	WlzCMeshDelUnusedNodes3D(3)
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <Wlz.h>

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0;
  FILE		*fP = NULL;
  char		*iFile,
  		*oFile;
  const char	*errMsg;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "ho:";
  const char    defFile[] = "-";

  opterr = 0;
  iFile = (char *)defFile;
  oFile = (char *)defFile;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        oFile = optarg;
	break;
      case 'h': /* FALLTROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if((oFile == NULL) || (*oFile == '\0') ||
       (iFile == NULL) || (*iFile == '\0'))
    {
      usage = 1;
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
      iFile = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    if((fP = (strcmp(iFile, "-")? fopen(iFile, "r"): stdin)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open input file %s.\n",
		     argv[0], iFile);
    }
  }
  if(ok)
  {
    obj = WlzReadObj(fP, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read input mesh object, %s.\n",
		     argv[0], errMsg);
    }
  }
  if(fP && strcmp(iFile, "-"))
  {
    (void )fclose(fP);
  }
  if(ok)
  {
    switch(obj->type)
    {
      case WLZ_CMESH_2D:
        WlzCMeshDelUnusedNodes2D(obj->domain.cm2);
        break;
      case WLZ_CMESH_2D5:
        WlzCMeshDelUnusedNodes2D5(obj->domain.cm2d5);
        break;
      case WLZ_CMESH_3D:
        WlzCMeshDelUnusedNodes3D(obj->domain.cm3);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to remove unused nodes (%s).\n",
      		     argv[0],
		     errMsg);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(oFile, "-")? fopen(oFile, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], oFile);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsg);
    }
  }
  if(fP && strcmp(oFile, "-"))
  {
    (void )fclose(fP);
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output>] [<input>]\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "Reads a conforming mesh object and deletes and unused nodes from it,\n"
    "ie those with no edge uses. By default the input object is read from\n"
    "the standard input and the output object is written to the standard\n"
    "output.\n",
    argv[0]);
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
