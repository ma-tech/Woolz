#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstCMeshCellStats_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstCMeshCellStats.c
* \author       Bill Hill
* \date         May 2009
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
* \brief	Outputs constrained mesh cell occupancy statistics.
* \ingroup	BinWlzTst
*/


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
	        minNodPerCell = 0,
		maxNodPerCell = 0,
	        minElmPerCell = 0,
	        maxElmPerCell = 0,
  		ok = 1,
		verbose = 0,
  		usage = 0;
  double	meanNodPerCell = 0.0,
  		meanElmPerCell = 0.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL;
  WlzCMeshP 	mesh;
  static char   optList[] = "hvo:";
  const char    inObjFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  opterr = 0;
  mesh.v = NULL;
  inObjFileStr = (char *)inObjFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'v':
        verbose = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      usage = 1;
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
  }
  ok = usage == 0;
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
      (void )fprintf(stderr,
                     "%s: Failed to read object from file (%s)\n",
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
      case WLZ_CMESH_2D:
        mesh.m2 = inObj->domain.cm2;
	break;
      case WLZ_CMESH_3D:
        mesh.m3 = inObj->domain.cm3;
	break;
      default:
        ok = 0;
	errNum = WLZ_ERR_OBJECT_TYPE;
        (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		   "%s: Invalid object type, must be WLZ_CMESH_[23]D (%s),\n",
		       argv[0],
		       errMsgStr);
        break;
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    WlzCMeshGetCellStats(mesh,
    			 &minNodPerCell, &maxNodPerCell, &meanNodPerCell,
			 &minElmPerCell, &maxElmPerCell, &meanElmPerCell);
    if(verbose != 0)
    {
      (void )fprintf(fP,
                     "Minimum nodes per cell    %d\n"
		     "Maximum nodes per cell    %d\n"
		     "Mean nodes per cell       %lg\n"
		     "Minimum elements per cell %d\n"
		     "Maximum elements per cell %d\n"
		     "Mean elements per cell    %lg\n",
		     minNodPerCell, maxNodPerCell, meanNodPerCell,
		     minElmPerCell, maxElmPerCell, meanElmPerCell);
    }
    else
    {
      (void )fprintf(fP,
                     "% 8d % 8d % 8lg % 8d % 8d % 8lg\n",
		     minNodPerCell, maxNodPerCell, meanNodPerCell,
		     minElmPerCell, maxElmPerCell, meanElmPerCell);
    }
  }
  (void )WlzFreeObj(inObj);
  if((fP != NULL) && (fP != stdout))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-v] [-o<output file>] [<input cmesh object>]\n"
    "Computes and outputs constrained mesh cell occupancy statistics.\n"
    "Options are:\n"
    "  -h  Help, prints this usage message.\n"
    "  -v  Verbose output.\n"
    "  -o  Output file.\n",
    argv[0]);

  }
  return(!ok);
}
