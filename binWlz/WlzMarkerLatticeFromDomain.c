#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMarkerLatticeFromDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzMarkerLatticeFromDomain.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Covers a spatial domain with a lattice of markers.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzmarkerlatticefromdomain "WlzMarkerLatticeFromDomain"
*/

/*!
\ingroup BinWlz
\defgroup wlzmarkerlatticefromdomain WlzMarkerLatticeFromDomain
\par Name
WlzMarkerLatticeFromDomain - covers a spatial domain with a lattice of markers.
\par Synopsis
\verbatim
WlzMarkerLatticeFromDomain [-h] [-o<output file>] [-s #] [-t <type>] [-z #]
		   [<input file>]
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
  <tr> 
    <td><b>-2</b></td>
    <td>marker separation.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Marker type: Only valid type is 'sphere'.</td>
  </tr>
  <tr> 
    <td><b>-2</b></td>
    <td>marker size.</td>
  </tr>
</table>
\par Description
Reads either a 2D or 3D spatial domain object from a given file or
the standard input (default). A new domain is then constructed which
covers the given domain with a lattice of markers.
\par File
\ref WlzMarkerLatticeFromDomain.c "WlzMarkerLatticeFromDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzMarkerLattice WlzMarkerLattice(3)
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
  		usage = 0,
		markerSep = 16,
		markerSz = 4;
  FILE		*fP = NULL;
  char		*iFile,
  		*oFile;
  const char	*errMsg;
  WlzObject	*iObj = NULL,
  		*oObj = NULL;
  WlzMarkerType	markerType = WLZ_MARKER_SPHERE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "ho:s:t:z:";
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
      case 's':
	if(sscanf(optarg, "%d", &markerSep) != 1)
	{
	  usage = 1;
	}
        break;
      case 't':
	markerType = WlzStringToMarkerType(optarg, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  usage = 1;
	}
        break;
      case 'z':
	if(sscanf(optarg, "%d", &markerSz) != 1)
	{
	  usage = 1;
	}
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
    iObj = WlzReadObj(fP, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read input domain object, %s.\n",
		     argv[0], errMsg);
    }
  }
  if(fP && strcmp(iFile, "-"))
  {
    (void )fclose(fP);
  }
  if(ok)
  {
    oObj = WlzMarkerLattice(iObj, markerType, markerSz, markerSep, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to create domain from vertices (%s).\n",
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
    errNum = WlzWriteObj(fP, oObj);
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
  (void )WlzFreeObj(iObj);
  (void )WlzFreeObj(oObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output>] [-s #] [-t <type>] [-z #] [<input>]\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -s  Marker separation.\n"
    "  -t  Marker type: Only valid type is 'sphere'.\n"
    "  -z  Marker size.\n"
    "Reads either a 2D or 3D spatial domain object from a given file or\n"
    "the standard input (default). A new domainis then constructed which\n"
    "covers the given domain with a lattice of markers.\n",
    argv[0]);
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
