#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzSnapFitObjs.c
* \author       Bill Hill
* \date         July 2004
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
* \brief	Finds a set of correspondences between two objects based
* 		on closest points.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzsnapfitobjs "WlzSnapFitObjs"
*/

/*!
\ingroup BinWlz
\defgroup wlzsnapfitobjs WlzSnapFitObjs
\par Name
WlzSnapFitObjs - finds a set of correspondences between two objects based
                 on closest points.
\par Synopsis
\verbatim
WlzSnapFitObjs [-h] [-o<output file>]
               [-A] [-i <transform>] [-d #] [-t #] [-s #]
	       <target object> <source object>

\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file for correspondences.</td>
  </tr>
  <tr> 
    <td><b>-A</b></td>
    <td>Correspondences output in absolute format.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Initial affine transform to be applied to the source object.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Maximum distance between any target source correspondences.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Minimum distance between target correspondence points.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Minimum distance between source correspondence points.</td>
  </tr>
</table>
\par Description
Computes a set of correspondences (tie points) from the given
source object to the given target object based on closest points.
The correspondences are output as ascii vertices either in relative
format:
\verbatim
  <sx> <sy>[ <sz> ]<tx - sx> <ty - sy>[ <tz - sz>]
\endverbatim
or absolute format:
\verbatim
    <tx> <ty>[ <tz> ]<sx> <sy>[ <sz>],
\endverbatim
where sx, sy, sz, tx, ty and tz are the source and target coordinates.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzSnapFitObjs.c "WlzSnapFitObjs.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzSnapFit "WlzSnapFit(3)"
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
  int		idN,
  		nVtx,
  		absTPMode = 0,
  		option,
		ok = 1,
		usage = 0;
  double	maxCDist = DBL_MAX,
  		minTDist = 10.0,
		minSDist = 10.0;
  WlzVertexType vType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzVertexP	tVtxP,
  		sVtxP;
  WlzObject	*tObj = NULL,
  		*sObj = NULL,
		*trObj = NULL;
  FILE		*fP = NULL;
  char 		*tObjFileStr,
  		*sObjFileStr,
		*trObjFileStr,
		*outFileStr;
  static char	optList[] = "hAi:d:s:t:o:";

  opterr = 0;
  trObjFileStr = NULL;
  tVtxP.v = sVtxP.v = NULL;
  trObjFileStr = NULL;
  sObjFileStr = tObjFileStr = outFileStr = "-";
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'A':
        absTPMode = 1; 
	break;
      case 'i':
        trObjFileStr = optarg;
	break;
      case 'd':
        if(sscanf(optarg, "%lg", &maxCDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 's':
        if(sscanf(optarg, "%lg", &minSDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 't':
        if(sscanf(optarg, "%lg", &minTDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default: /* FALLTHROUGH */
        usage = 1;
	break;
    }
  }
  if(ok && (optind < argc))
  {
    if((optind + 2) != argc)
    {
      usage = 1;
    }
    else
    {
      tObjFileStr = *(argv + optind);
      sObjFileStr = *(argv + optind + 1);
    }
  }
  ok = !usage;
  /* Read the target object. */
  if(ok)
  {
    if((tObjFileStr == NULL) ||
       (*tObjFileStr == '\0') ||
       ((fP = (strcmp(tObjFileStr, "-")?
	      fopen(tObjFileStr, "r"): stdin)) == NULL) ||
       ((tObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read target object from file %s.\n",
		     *argv, tObjFileStr);
    }
    if(fP && strcmp(tObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the source object. */
  if(ok)
  {
    if((sObjFileStr == NULL) ||
       (*sObjFileStr == '\0') ||
       ((fP = (strcmp(sObjFileStr, "-")?
	      fopen(sObjFileStr, "r"): stdin)) == NULL) ||
       ((sObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read source object from file %s.\n",
		     *argv, sObjFileStr);
    }
    if(fP && strcmp(sObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the optional initial transform object. */
  if(ok && trObjFileStr)
  {
    if((*trObjFileStr == '\0') ||
       ((fP = (strcmp(trObjFileStr, "-")?
	      fopen(trObjFileStr, "r"): stdin)) == NULL) ||
       ((trObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (trObj == NULL) || (trObj->type != WLZ_AFFINE_TRANS) ||
       (trObj->domain.core == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		 "%s: Failed to read affine transform object from file %s.\n",
		     *argv, trObjFileStr);
    }
    if(fP && strcmp(trObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Compute correspondences. */
  if(ok)
  {
    errNum = WlzSnapFit(tObj, sObj, (trObj)? trObj->domain.t: NULL,
    			&vType, &nVtx, &tVtxP, &sVtxP,
			maxCDist, minTDist, minSDist);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to compute correspondences.\n",
		     *argv);
    }
  }
  /* Output the correspondences. */
  if(ok)
  {
    if((outFileStr == NULL) ||
       (*outFileStr == '\0') ||
       ((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to write correspondences to file %s\n",
		     *argv, outFileStr);
    }
  }
  if(ok)
  {
    if(vType == WLZ_VERTEX_D2)
    {
      if(absTPMode)
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g\n", 
	  		 (tVtxP.d2 + idN)->vtX,
			 (tVtxP.d2 + idN)->vtY,
	  		 (sVtxP.d2 + idN)->vtX,
			 (sVtxP.d2 + idN)->vtY);
	}
      }
      else
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g\n", 
	  		 (sVtxP.d2 + idN)->vtX,
			 (sVtxP.d2 + idN)->vtY,
	  		 (tVtxP.d2 + idN)->vtX - (sVtxP.d2 + idN)->vtX,
	  		 (tVtxP.d2 + idN)->vtY - (sVtxP.d2 + idN)->vtY);
	}
      }
    }
    else /* vType == WLZ_VERTEX_D3 */
    {
      if(absTPMode)
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g %g %g\n", 
	  		 (tVtxP.d3 + idN)->vtX,
			 (tVtxP.d3 + idN)->vtY,
			 (tVtxP.d3 + idN)->vtZ,
	  		 (sVtxP.d3 + idN)->vtX,
			 (sVtxP.d3 + idN)->vtY,
			 (sVtxP.d3 + idN)->vtZ);
	}
      }
      else
      {
	for(idN = 0; idN < nVtx; ++idN)
	{
	  (void )fprintf(fP, "%g %g %g %g %g %g\n", 
	  		 (sVtxP.d3 + idN)->vtX,
			 (sVtxP.d3 + idN)->vtY,
			 (sVtxP.d3 + idN)->vtZ,
	  		 (tVtxP.d3 + idN)->vtX - (sVtxP.d3 + idN)->vtX,
	  		 (tVtxP.d3 + idN)->vtY - (sVtxP.d3 + idN)->vtY,
	  		 (tVtxP.d3 + idN)->vtZ - (sVtxP.d3 + idN)->vtZ);
	}
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP); fP = NULL;
    }
    
  }
  (void )WlzFreeObj(tObj);
  (void )WlzFreeObj(sObj);
  (void )WlzFreeObj(trObj);
  (void )AlcFree(tVtxP.v);
  (void )AlcFree(sVtxP.v);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-A] [-i <transform>] [-d #] [-t #] [-s #]\n"
    "        [-o<output file>] <target> <source>\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -A  Correspondences output in absolute format.\n"
    "  -i  Initial affine transform to be applied to the source object.\n"
    "  -d  Maximum distance between any target source correspondences.\n"
    "  -t  Minimum distance between target correspondence points.\n"
    "  -s  Minimum distance between source correspondence points.\n"
    "  -o  Output file for correspondences.\n"
    "Computes a set of correspondences (tie points) from the given\n"
    "source object to the given target object based on closest points.\n"
    "The correspondences are output as ascii vertices either in relative\n"
    "format:\n"
    "  <sx> <sy>[ <sz> ]<tx - sx> <ty - sy>[ <tz - sz>]\n"
    "or absolute format:\n"
    "  <tx> <ty>[ <tz> ]<sx> <sy>[ <sz>],\n"
    "where sx, sy, sz, tx, ty and tz are the source and target coordinates.\n",
    *argv,
    " -3 -d 10 -i tr.wlz -o out.num -t 20 -s 20 trg.wlz src.wlz\n"
    "Reads an initial source affine transform from tr.wlz and computes\n"
    "a set of corresponding closest points such that the minimum distance\n"
    "between and target and source pair is less than 10, no target points\n"
    "are closer than 20 and no target points have corresponding points\n"
    "closer than 30 in the source object. all distances are in pixels/\n"
    "voxels.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
