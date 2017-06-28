#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzThinObjToPoints_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzThinObjToPoints.c
* \author       Bill Hill
* \date         April 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Outputs a point object given a 2 or 3D spatial domain
* 		with the points computed by thinning the spatial domain
* 		object..
* \ingroup	wlzthinobjtopoints "WlzThinObjToPoints"
*/


/*!
\ingroup BinWlz
\defgroup wlzthinobjtopoints WlzThinObjToPoints
\par Name
WlzThinObjToPoints - computes points by thinning spatial domain object.
\par Synopsis
\verbatim
WlzThinObjToPoints [-G] [-g[<start>][,<inc>]] [-o<output object>]
                   [-h] [-T] [-D #,#[#]] [-x] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-D</b></td>
    <td>Dither the points by applying a random offset. Supplied values
        are the maximum dither displacement.</td>
  </tr>
  <tr> 
    <td><b>-G</b></td>
    <td>Use given object grey values if the object has them.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>The initial and increment grey values.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Report elapsed time.</td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>Use voxel size scaling.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Computes points from a given spatial domain object.
If using grey values
the object is erroded by successive thresholding
while collecting all disconnected regions
or if the object has no values the object is simply labeled
to give a collection of disconnected regions.
The disconnected regions are then eroded (morphologicaly)
again collecting disconnected regions.
Finally the centroids of the disconnected regions are used to compute
the points object.
The default initial grey value is 0 and the default increment is 1.
By default all files are read from the standard input and written to
the standard output.
\par Examples
\verbatim
WlzThinObjToPoints -o points.wlz -g 250,-1 in.wlz
\endverbatim
Creates a points object containing locations within the given
spatial domain object (read from the file in.wlz) in which the
points at located at the centroids of the thinned object's regions,
where grey value thinning is used with an initial threshold of 250
reducing by 1 at each threshold.
The resulting points object is written to the file out.wlz
\par File
\ref WlzThinObjToPoints.c "WlzThinObjToPoints.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzThinToPoints "WlzThinToPoints(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <Wlz.h>

#if !defined(HAVE_STRSEP)
#define strsep(B,S) strtok((B),(S))
#endif

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
		usage = 0,
		dither = 0,
		timeFlg = 0,
		gStart = 0,
		gInc = 1,
		useGrey = 0,
		voxelScaling = 0;
  char		*inFileStr,
		*outFileStr;
  FILE		*fP = NULL;
  WlzDVertex3	ditherVal = {0.0, 0.0, 0.0};
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  struct timeval times[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hGTxD:g:o:";
  const char    fileStrDef[] = "-";

  /* Parse the argument list and check for input files. */
  opterr = 0;
  inFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  while(ok && !usage && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'D':
        dither = 1;
	if(sscanf(optarg, "%lg,%lg,%lg",
	          &(ditherVal.vtX), &(ditherVal.vtY), &(ditherVal.vtZ)) < 2)
        {
	  usage = 1;
	}
	break;
      case 'G':
	useGrey = 1;
	break;
      case 'g':
	{
	  char  *tok,
	  	*buf,
		*savbuf;

	  buf = savbuf = AlcStrDup(optarg);
	  tok = strsep(&buf, ",");
	  if(tok && *tok)
	  {
	    if(sscanf(tok, "%d", &gStart) != 1)
	    {
	      usage = 1;
	    }
	  }
	  if(usage == 0)
	  {
	    tok = strsep(&buf, ",");
	    if(tok && *tok)
	    {
	      if(sscanf(tok, "%d", &gInc) != 1)
	      {
	        usage = 1;
	      }
	    }
	  }
	  AlcFree(savbuf);
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'T':
        timeFlg = 1;
	break;
      case 'x':
        voxelScaling = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(ok)
  {
    ok = !usage;
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
      inFileStr = argv[optind];
    }
  }
  /* Read the input object. */
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read input object from file %s.\n",
                     argv[0], inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Compute the points object. */
  if(ok)
  {
    if(timeFlg)
    {
      gettimeofday(times + 0, NULL);
    }
    outObj = WlzThinToPoints(inObj, useGrey, gStart, gInc, &errNum);
    if(timeFlg)
    {
      gettimeofday(times + 1, NULL);
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: Elapsed time for WlzPointsFromDomObj()  %gus\n",
		     argv[0],
		     (1000000.0 * times[2].tv_sec) + times[2].tv_usec);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Perform voxel scaling if required. */
      if(voxelScaling &&
	 (inObj->type == WLZ_3D_DOMAINOBJ) &&
         (outObj->domain.pts->type == WLZ_POINTS_3I))
      {
	int	   nPts;
	WlzDomain  sclDom;
	WlzVertexP nullP;

	nullP.v = NULL;
	nPts = outObj->domain.pts->nPoints;
	sclDom.pts = WlzMakePoints(WLZ_POINTS_3D, 0, nullP,
	                           outObj->domain.pts->nPoints, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
	  int	i;
	  WlzIVertex3 *p;
	  WlzDVertex3 *q;
	  WlzDVertex3  s;
	  
	  sclDom.pts->nPoints = nPts;
	  q = sclDom.pts->points.d3;
	  p = outObj->domain.pts->points.i3;
	  s.vtX = inObj->domain.p->voxel_size[0];
	  s.vtY = inObj->domain.p->voxel_size[1];
	  s.vtZ = inObj->domain.p->voxel_size[2];
	  for(i = 0; i < nPts; ++i)
	  {
	    q[i].vtX = p[i].vtX * s.vtX;
	    q[i].vtY = p[i].vtY * s.vtY;
	    q[i].vtZ = p[i].vtZ * s.vtZ;
	  }
	  (void )WlzFreeDomain(outObj->domain);
	  outObj->domain = WlzAssignDomain(sclDom, NULL);
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(dither)
      {
	WlzDomain dtrDom;

	dtrDom.pts = WlzPointsDither(outObj->domain.pts, ditherVal, 
	                             inObj, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  (void )WlzFreeDomain(outObj->domain);
	  outObj->domain = WlzAssignDomain(dtrDom, NULL);
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to compute points object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  /* Output the points object. */
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
    errNum = WlzWriteObj(fP, outObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
                     argv[0], errMsgStr);
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-G] [-g[<start>][,<inc>]] [[-D#,#[,#]]\n"
    "\t\t[-o<output file>] [-T] [-x] [-h] [<input object file>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -G  Use given object grey values if the object has them.\n"
    "  -D  Dither the points by applying a random offset. Supplied values\n"
    "      are the maximum dither displacement.\n"
    "  -g  The initial and increment grey values.\n"
    "  -o  Output file.\n"
    "  -T  Report elapsed time.\n"
    "  -x  Use voxel size scaling.\n"
    "  -h  Help, prints usage message.\n"
    "Computes points from a given spatial domain object. If using grey\n"
    "values the object is erroded by successive thresholding while\n"
    "collecting all disconnected regions or if the object has no values\n"
    "the object is simply labeled to give a collection of disconnected\n"
    "regions. The disconnected regions are then eroded (morphologicaly)\n"
    "again collecting disconnected regions. Finally the centroids of the\n"
    "disconnected regions are used to compute the points object.\n"
    "The default initial grey value is 0 and the default increment is 1.\n"
    "By default all files are read from the standard input and written to\n"
    "the standard output.\n"
    "Example:\n"
    "  %s -o points.wlz -g 250,-1 in.wlz\n"
    "creates a points object containing locations within the given\n"
    "spatial domain object (read from the file in.wlz) in which the\n"
    "points at located at the centroids of the thinned object's regions,\n"
    "where grey value thinning is used with an initial threshold of 250\n"
    "reducing by 1 at each threshold.\n"
    "The resulting points object is written to the file out.wlz\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
