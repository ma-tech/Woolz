#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzOffsetDistance_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzOffsetDistance.c
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
* \brief	Computes an equidistant object between a pair of doimain
* 		objects.
* \ingroup	wlzoffsetdistance "WlzOffsetDistance"
*/

/*!
\ingroup BinWlz
\defgroup wlzoffsetdistance wlzOffsetDistance
\par Name
WlzOffsetDistance - computes equidistant object between a pair of domain
		    objects.
\par Synopsis
\verbatim
WlzOffsetDistance [-m#] [-S] [-o<out file>] [-h] [<object 1>] [<object 2>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-m</b></td>
    <td>Maximum distance to consider (default 0, limit is convex hull).</td>
  </tr>
  <tr> 
    <td><b>-S</b></td>
    <td>Output statistics rather than distance object.
        Output statistics are printed on a single line with the order:
	volume, min, max, sum, sum of squares, mean, standard deviation.
    </td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Computes an object with domain and values that are the
set of minimum distances between the two given objects.
The equidistant boundary is computed between the domains
of the two given objects, within the maximum distance limit
and within the convex hull of the union of the two given object's
domains.
By default all files are read from the standard input and written to
the standard output.
\par Examples
\verbatim
WlzOffsetDistance -o out.wlz obj1.wlz obj2.wlz
\endverbatim
Creates a new domain object corresponding to the minimum equidistant
voxel set between the given objects, with maximum distance constrained
only be the convex hull of the two given objects.
The resulting object is written to the file out.wlz.
\par File
\ref WlzOffsetDistance.c "WlzOffsetDistance.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzOffsetDist "WlzOffsetDist(3)"
\ref wlzclassifyrcc "WlzClassifyRCC(1)"
\ref wlzgreystats "WlzGreyStats(1)"
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
  int		ok = 1,
  		option,
		usage = 0,
		stats = 0;
  char		*inFileStr[2];
  char		*outFileStr;
  int		maxDist = 0;
  WlzObject	*outObj = NULL;
  WlzObject	*inObj[2] = {NULL};
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hSm:o:";
  const char    fileStrDef[] = "-";

  /* Parse the argument list and check for input files. */
  opterr = 0;
  outFileStr = (char *)fileStrDef;
  inFileStr[0] = inFileStr[1] = (char *)fileStrDef;
  while(ok && !usage && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'm':
	if((sscanf(optarg, "%d", &maxDist) != 1) || (maxDist < 0))
	{
	  usage = 1;
	}
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
  if(ok)
  {
    ok = !usage;
  }
  if(ok && (optind < argc))
  {
    if((optind + 2) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      inFileStr[0] = argv[optind + 0];
      inFileStr[1] = argv[optind + 1];
    }
  }
  /* Read the reference object. */
  if(ok)
  {
    int		idx;

    for(idx = 0; idx < 2; ++idx)
    {
      FILE	*fP = NULL;

      if((inFileStr[idx] == NULL) ||
	 (*inFileStr[idx] == '\0') ||
	 ((fP = (strcmp(inFileStr[idx], "-")?
		fopen(inFileStr[idx], "r"): stdin)) == NULL) ||
	 ((inObj[idx] = WlzAssignObject(
	                WlzReadObj(fP, &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: Failed to read object from file %s.\n",
		       argv[0], inFileStr[idx]);
      }
      if(fP && strcmp(inFileStr[idx], "-"))
      {
	(void )fclose(fP); fP = NULL;
      }
    }
  }
  /* Compute the equi-distant distance object. */
  if(ok)
  {
    outObj = WlzAssignObject(
    	     WlzOffsetDist(inObj[0], inObj[1], maxDist, &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		 "%s: Failed to compute equi-distant distance object, %s.\n",
		 argv[0], errMsgStr);
    }
  }
  /* Output the nearby domain object. */
  if(ok)
  {
    FILE	*fP = NULL;

    if((fP = (strcmp(outFileStr, "-")?
	    fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
	  "%s: Failed to open output file %s.\n",
	  argv[0], outFileStr);
    }
    if(ok)
    {
      if(stats == 0)
      {
	/* Output equi-distant distance object. */
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
      }
      else
      {
	int	n;
	double	dMin,
		  dMax,
		  dSum,
		  dSSq,
		  dMean,
		  dSDev;
	n = WlzGreyStats(outObj,
			 NULL, &dMin, &dMax, &dSum, &dSSq, &dMean, &dSDev,
			 &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	  (void )fprintf(stderr,
	      "%s: Failed to compute equi-distant distance statistics, %s.\n",
	      argv[0], errMsgStr);
	}
	else
	{
	  (void )fprintf(fP, "%d %lg %lg %lg %lg %lg %lg\n",
			 n, dMin, dMax, dSum, dSSq, dMean, dSDev);
	}
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-m#] [-S] [-o<output file>] [-h]\n"
    "\t\t[<object 1>] [<object 2>]\n"
    "Version: %s\n"
    "Computes an object with domain and values that are the set of minimum\n"
    "distances between the two given objects. The equidistant boundary is\n"
    "computed between the domains of the two given objects, within the\n"
    " distance limit and within the convex hull of the union of the two\n"
    "given object's domains.\n"
    "By default all files are read from the standard input and written\n"
    "to the standard output.\n"
    "Options:\n"
    "  -m  Maximum distance to consider (default 0, limit is convex hull).\n"
    "  -S  Output statistics rather than distance object. Output statistics\n"
    "      are printed on a single line with the order:\n"
    "      volume, min, max, sum, sum of squares, mean, standard deviation.\n"
    "  -o  Output object file.\n"
    "  -h  Help - prints this usage message\n"
    "By default all files are read from the standard input and written to\n"
    "the standard output.\n"
    "Example:\n"
    "  %s -o out.wlz -o out.wlz obj1.wlz obj2.wlz\n"
    "which creates a domain object corresponding to the minimum equidistant\n"
    "voxel set between the given objects, with maximum distance constrained\n"
    "only be the convex hull of the two given objects. The resulting object\n"
    "is written to the file out.wlz.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
