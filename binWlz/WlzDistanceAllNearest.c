#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDistanceAllNearest_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzDistanceAllNearest.c
* \author       Bill Hill
* \date         April 2021
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
* \brief	Computes location of closest point in reference domain
* 		for all points in foreground domain.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzdistanceallnearest "WlzDistanceAllNearest"
*/

//HACK TODO UPDATE doxygen usage

/*!
\ingroup BinWlz
\defgroup wlzdistanceallnearest WlzDistanceAllNearest
\par Name
WlzDistanceAllNearest - computes the location of the closest point in a
			reference domain for all points in a foreground domain.
\par Synopsis
\verbatim
WlzDistanceAllNearest [-d#] [-f<foreground object>] [-o<output>] [-h]
                     [<reference object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Distance function:
      <table width="500" border="0">
      <tr> <td>1</td> <td>octagonal (2D and 3D) - default</td></tr>
      <tr> <td>4</td> <td>4-connected (2D)</td></tr>
      <tr> <td>8</td> <td>8-connected (2D)</td></tr>
      <tr> <td>6</td> <td>6-connected (3D)</td></tr>
      <tr> <td>18</td> <td>18-connected (3D)</td></tr>
      <tr> <td>26</td> <td>26-connected (3D)</td></tr>
    </td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>The foreground object file.</td>
  </tr>
</table>
\par Description
Computes the location of the closest point in a reference domain for
all points within the foreground domain.
 The output object is a compound object in which the first component
objects (2 for 2D and 3 for 3D) are the coordinates of the closest point
in the reference domain and the next (3rd for 2D and 4th for 3D) is the
distance.
\par Examples
\verbatim
WlzDistanceAllNearest -f fish4.wlz fish4eye.wlz > out.wlz
\endverbatim
Reads the given 2D domain objects and then computes the closest point
to and distance from the domain of fish4eye.wlz within the domain of
fish4.wlz. The output compound array will have 3 component objects,
with the first two encoding the integer column (x) and row (y) followed
by the third object with float distances from the closest point.
\par File
\ref WlzDistanceAllNearest.c "WlzDistanceAllNearest.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzDistAllNearest "WlzDistAllNearest(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok,
  		con,
  		option,
		usage = 0;
  const char	*forFileStr,
  		*refFileStr,
		*outFileStr;
  FILE		*fP = NULL;
  WlzObject	*disObj = NULL,
  		*forObj = NULL,
		*refObj = NULL;
  WlzCompoundArray *nrpObj = NULL,
                   *outObj = NULL;
  WlzDistanceType dFn;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hd:f:o:";
  const char    fileStrDef[] = "-";
  const WlzConnectType defDFn = WLZ_OCTAGONAL_DISTANCE;

  /* Parse the argument list and check for input files. */
  opterr = 0;
  refFileStr = fileStrDef;
  forFileStr = fileStrDef;
  outFileStr = fileStrDef;
  dFn = defDFn;
  while((option = getopt(argc, argv, optList)) != EOF)
  {
    switch(option)
    {
      case 'd':
	if(sscanf(optarg, "%d", &con) != 1)
	{
	  usage = 1;
	}
	else
	{
	  switch(con)
	  {
	    case 1:
	      dFn = WLZ_OCTAGONAL_DISTANCE;
	      break;
	    case 4:
	      dFn = WLZ_4_DISTANCE;
	      break;
	    case 6:
	      dFn = WLZ_6_DISTANCE;
	      break;
	    case 8:
	      dFn = WLZ_8_DISTANCE;
	      break;
	    case 18:
	      dFn = WLZ_18_DISTANCE;
	      break;
	    case 26:
	      dFn = WLZ_26_DISTANCE;
	      break;
	    default:
	      usage = 1;
	      break;
	  }
	}
	break;
      case 'f':
	forFileStr = optarg;
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
  ok = !usage;
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      refFileStr = argv[optind];
    }
  }
  /* Read the reference object. */
  if(ok)
  {
    if((refFileStr == NULL) ||
       (*refFileStr == '\0') ||
       ((fP = (strcmp(refFileStr, "-")?
              fopen(refFileStr, "r"): stdin)) == NULL) ||
       ((refObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read reference object from file %s.\n",
                     argv[0], refFileStr);
    }
    if(fP && strcmp(refFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the foreground object. */
  if(ok)
  {
    if((forFileStr == NULL) ||
       (*forFileStr == '\0') ||
       ((fP = (strcmp(forFileStr, "-")?
              fopen(forFileStr, "r"): stdin)) == NULL) ||
       ((forObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read foreground object from file %s.\n",
                     argv[0], forFileStr);
    }
    if(fP && strcmp(forFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Compute the nearest point and distance compound object. */
  if(ok)
  {
    nrpObj = WlzDistAllNearest(forObj, refObj, dFn, &disObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to compute nearest points object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  /* Make new compound object which includes distances. */
  if(ok)
  {
    outObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, nrpObj->n + 1, NULL,
                                  disObj->type, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to create output compound object, %s.\n",
		     argv[0], errMsgStr);
    }
    else
    {
      for(int i = 0; i < nrpObj->n; ++i)
      {
	outObj->o[i] = WlzAssignObject(nrpObj->o[i], NULL);
      }
      outObj->o[nrpObj->n] = WlzAssignObject(disObj, NULL);
      disObj = NULL;
    }
  }
  /* Output the compound object. */
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
    errNum = WlzWriteObj(fP, (WlzObject *)outObj);
    (void )fclose(fP);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
                     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(disObj);
  (void )WlzFreeObj(forObj);
  (void )WlzFreeObj(refObj);
  (void )WlzFreeCompoundArray(nrpObj);
  (void )WlzFreeCompoundArray(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-d#] [-f<foreground object>] [-o<output>] [-h]\n"
    "\t\t[<reference object>]\n"
    "Computes the location of the closest point in a reference domain for\n"
    "all points within the foreground domain.\n"
    "The output object is a compound object in which the first component\n"
    "objects (2 for 2D and 3 for 3D) are the coordinates of the closest point\n"
    "in the reference domain and the next (3rd for 2D and 4th for 3D) is the\n"
    "distance. Distances are correct for the given distance function and are\n"
    "not normalised to remove over estimates with respect to constrained\n"
    "Euclidean distance.\n"
    "Version: %s\n"
    "Options:\n"
    "  -d  Distance function:\n"
    "              1: octagonal (2D and 3D) - default\n"
    "              4: 4-connected (2D)\n"
    "              8: 8-connected (2D)\n"
    "              6: 6-connected (3D)\n"
    "             18: 18-connected (3D)\n"
    "             26: 26-connected (3D)\n"
    "  -f  The foreground object file.\n" 
    "  -o  Output object file.\n"
    "  -h  Help - prints this usage message\n"
    "Example:\n"
    "\t%s -f fish4.wlz fish4eye.wlz > out.wlz\n"
    "Reads the given 2D domain objects and then computes the closest point\n"
    "to and distance from the domain of fish4eye.wlz within the domain of\n"
    "fish4.wlz. The output compound array will have 3 component objects,\n"
    "with the first two encoding the integer column (x) and row (y) followed\n"
    "by the third object with float distances from the closest point.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
