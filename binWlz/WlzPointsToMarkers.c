#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPointsToMarkers_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzPointsToMarkers.c
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
* \brief	Creates a domain with a marker located at the position
* 		of each point vertex using a points object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzpointstomarkers "WlzPointsToMarkers"
*/

/*!
\ingroup BinWlz
\defgroup wlzmarkerstodomain WlzPointsToMarkers
\par Name
WlzPointsToMarkers - creates a domain with a marker located at the position
		     of the vertex of each point in a points object.
\par Synopsis
\verbatim
WlzPointsToMarkers [-h] [-T] [-o<output file>] [-s #] [-t <type>]
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
    <td><b>-s</b></td>
    <td>marker size.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Marker type: Only valid type is 'sphere'.</td>
  </tr>
</table>
\par Description
Reads a 2D or 3D points object from either a given file or the
standard input (default). The vertices of the points are then used
to create 2D or 3D domain with a markers at the positions of the
point vertices.
\par Examples
\verbatim
WlzPointsToMarkers -o out.wlz -s 3 in.num
\endverbatim
Creates a new domain with spheres or radius 3 at the coordinates of
the points in the points object read from the file in.wlz.
\par File
\ref WlzPointsToMarkers.c "WlzPointsToMarkers.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzPointsFromDomain "WlzPointsFromDomain(1)"
\ref WlzPointsToMarkers  "WlzPointsToMarkers(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <Wlz.h>

#define WLZ_CFP_READLN_LEN	(1024)

static WlzVertexP 		WlzMTDReadVtxArray(
				  FILE *fP,
				  int dim,
				  int *dstNVtx, 
				  WlzVertexType *dstVType,
				  WlzErrorNum *dstErr);

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0,
		markerSz = 1,
		timeFlg = 0;
  FILE		*fP = NULL;
  char		*ptsFile,
  		*outFile;
  const char	*errMsg;
  struct timeval times[3];
  WlzObject	*ptsObj = NULL,
  		*outObj = NULL;
  WlzMarkerType	markerType = WLZ_MARKER_SPHERE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "hTo:s:t:";
  const char    fileDef[] = "-";

  opterr = 0;
  ptsFile = (char *)fileDef;
  outFile = (char *)fileDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outFile = optarg;
	break;
      case 's':
	if(sscanf(optarg, "%d", &markerSz) != 1)
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
      case 'T':
        timeFlg = 1;
	break;
      case 'h': /* FALLTROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    if((outFile == NULL) || (*outFile == '\0') ||
       (ptsFile == NULL) || (*ptsFile == '\0'))
    {
      usage = 1;
    }
  }
  if(!usage && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      ptsFile = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    if((ptsFile == NULL) ||
       (*ptsFile == '\0') ||
       ((fP = (strcmp(ptsFile, "-")?
              fopen(ptsFile, "r"): stdin)) == NULL) ||
       ((ptsObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_READ_EOF;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to read input file %s (%s).\n",
		     *argv, ptsFile, errMsg);
    }
  }
  if(fP)
  {
    if(strcmp(ptsFile, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(ok)
  {
    if(ptsObj->type != WLZ_POINTS)
    {
      ok = 0;
      errNum = WLZ_ERR_OBJECT_TYPE;
      (void )fprintf(stderr,
	             "%s: Object read from input file %s is not a points "
		     "object.\n",
		     *argv, ptsFile);
    }
  }
  if(ok)
  {
    if(timeFlg)
    {
      gettimeofday(times + 0, NULL);
    }
    outObj = WlzPointsToMarkers(ptsObj->domain.pts, markerType, markerSz,
    			        &errNum);
    if(timeFlg)
    {
      gettimeofday(times + 1, NULL);
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: Elapsed time for WlzPointsToMarkers() %gus\n",
		     argv[0],
		     (1000000.0 * times[2].tv_sec) + times[2].tv_usec);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to create marker domain from points (%s).\n",
      		     argv[0],
		     errMsg);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outFile, "-")? fopen(outFile, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outFile);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, outObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsg);
    }
  }
  if(fP && strcmp(outFile, "-"))
  {
    (void )fclose(fP);
  }
  (void )WlzFreeObj(ptsObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-T] [-o<output file>] [-s #] [-t <type>]\n"
    "\t\t[<input file>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -s  Marker size.\n"
    "  -t  Marker type: Only valid type is 'sphere'.\n"
    "  -T  Report elapsed time.\n"
    "Reads a 2D or 3D points object from either a given file or the\n"
    "standard input (default). The vertices of the points are then used\n"
    "to create 2D or 3D domain with a markers at the positions of the\n"
    "point vertices.\n",
    argv[0],
    WlzVersion());
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
