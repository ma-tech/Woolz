#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPointsFromDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzPointsFromDomain.c
* \author       Bill Hill
* \date         October 2014
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
* \brief	Outputs a point object given a 2 or 3D spatial domain.
* \ingroup	wlzpointsfromdomain "WlzPointsFromDomain"
*/


/*!
\ingroup BinWlz
\defgroup wlzpointsfromdomain WlzPointsFromDomain
\par Name
WlzPointsFromDomain - computes points from a given spatial domain object.
\par Synopsis
\verbatim
WlzPointsFromDomain [-d#] [-D #,#[#]] [-g[<min>][,<max>][,<gam>]] 
                [-o<output object>] [-h] [-T] [-x] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-d</b></td>
    <td>Minmum distance between points.</td>
  </tr>
  <tr> 
    <td><b>-D</b></td>
    <td>Dither the points by applying a random offset. Supplied values
        are the maximum dither displacement.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Use given object grey values to determine point density.</td>
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
If the given objects grey values are used, then the probability of
a point being placed is proportional to 
\f[
  \left\{ \begin{array}{ll}
  0.0 & g_i < g_{min} \\
  frac{(g_i - g_{min})^\gamma}{g_{max} - g_{min}} & g_i > g_{min}
  \end{array} \right.
\f]
By default all files are read from the standard input and written to
the standard output.
\par Examples
\verbatim
WlzPointsFromDomain -o out.wlz -d 33 in.wlz
\endverbatim
Creates a points object containing locations within the given
spatial domain object (read from the file in.wlz) which are at least
33 pixels or voxels distant from any other points. This object is
written to the file out.wlz
\par File
\ref WlzPointsFromDomain.c "WlzPointsFromDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzPointsFromDomObj "WlzPointsFromDomObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <Wlz.h>

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
		useGrey = 0,
		voxelScaling = 0;
  char		*inFileStr,
		*outFileStr;
  double	dMin = 32.0,
                gMin = 0.0,
		gMax = 255.0,
		gGam = 1.0;
  FILE		*fP = NULL;
  WlzDVertex3	ditherSz = {0.0};
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  struct timeval times[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hTxd:D:g:o:";
  const char    fileStrDef[] = "-";

  /* Parse the argument list and check for input files. */
  opterr = 0;
  inFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  while(ok && !usage && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	if((sscanf(optarg, "%lg", &dMin) != 1) || (dMin < 1.0))
	{
	  usage = 1;
	}
        break;
      case 'D':
	dither = 1;
	if(sscanf(optarg, "%lg,%lg,%lg",
	          &(ditherSz.vtX), &(ditherSz.vtY), &(ditherSz.vtZ)) < 2)
	{
	  usage = 1;
	}
	break;
      case 'g':
	useGrey = 1;
	{
	  char  *tok,
	  	*buf;

	  buf = AlcStrDup(optarg);
	  tok = strsep(&buf, ",");
	  if(tok && *tok)
	  {
	    if(sscanf(tok, "%lg", &gMin) != 1)
	    {
	      usage = 1;
	    }
	  }
	  if(usage == 0)
	  {
	    tok = strsep(&buf, ",");
	    if(tok && *tok)
	    {
	      if(sscanf(tok, "%lg", &gMax) != 1)
	      {
	        usage = 1;
	      }
	    }
	  }
	  if(usage == 0)
	  {
	    tok = strsep(&buf, ",");
	    if(tok && *tok)
	    {
	      if(sscanf(tok, "%lg", &gGam) != 1)
	      {
	        usage = 1;
	      }
	    }
	  }
	  AlcFree(buf);
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
    WlzDomain outDom;
    WlzValues nullVal;

    outDom.core = NULL;
    nullVal.core = NULL;
    if(timeFlg)
    {
      gettimeofday(times + 0, NULL);
    }
    outDom.pts = WlzPointsFromDomObj(inObj, dMin, voxelScaling,
                                     useGrey, gMin, gMax, gGam, &errNum);
    if((errNum == WLZ_ERR_NONE) && dither)
    {
      WlzPoints	*pts;

      pts = WlzPointsDither(outDom.pts, ditherSz, &errNum);
      (void )WlzFreeDomain(outDom);
      outDom.pts = pts;
    }
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
      outObj = WlzAssignObject(
      	       WlzMakeMain(WLZ_POINTS, outDom, nullVal, NULL, NULL, &errNum),
	       NULL);
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
    "Usage: %s [-d#] [-D#,#[,#]] [-g[<min>][,<max>][,<gam>]]\n"
    "\t\t[-o<output file>] [-T] [-x] [-h] [<Reference object file>]\n"
    "Computes points from a given spatial domain object.\n"
    "If the given objects grey values are used, then the probability of\n"
    "a point being placed is proportional to:\n"
    " 0.0 if g_i < g_min else ((g_i - g_min)^gam)/(g_max - g_min).\n"
    "Version: %s\n"
    "Options:\n"
    "  -d  Minimum distance between points.\n"
    "  -D  Dither the points by applying a random offset. Supplied values\n"
    "      are the maximum dither displacement.\n"
    "  -g  Use given object grey values to determine point density.\n"
    "  -o  Output object file.\n"
    "  -T  Report elapsed time.\n"
    "  -x  Use voxel size scaling.\n"
    "  -h  Help - prints this usage message\n"
    "By default all files are read from the standard input and written to\n"
    "the standard output.\n"
    "Example:\n"
    "  %s -o out.wlz -d 33 in.wlz\n"
    "which creates a points object containing locations within the given\n"
    "spatial domain object (read from the file in.wlz) which are at least\n"
    "33 pixels or voxels distant from any other points. This object is\n"
    "written to the file out.wlz\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
