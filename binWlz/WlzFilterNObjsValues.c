#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFilterNObjsValues_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzFilterNObjsValues.c
* \author       Bill Hill
* \date         September 2007
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
* \brief	Filters the values of multiple domain objects to construct a
* 		new domainn object in which each value is formed by
*		filtering the values of all input objects at the same
*		position as the values of the resulting object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzfilternobjsvalues "WlzFilterNObjsValues"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzfilternobjsvalues  WlzFilterNObjsValues
\par Name
WlzFilterNObjsValues - Filters the values of multiple domain objects to
                     construct a construct a new domainn object.

\par Synopsis
\verbatim
WlzFilterNObjsValues [-h] [-o<output file>] [-f <ref object>] [-m] [-r #]
                   [-i] [-u] [<input objects>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-o</b></td>
    <td>Output file name, default to standard out.</td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>Reference object.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Use a mean filter.</td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Use a rank filter. The parameter (range [0.0-1.0]) specifies the
        rank with: 0.0 minimum, 0.5 median and 1.0 maximum - intermediate
	values specify intermediate ranks.</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>The computed reference domain is the intersection of all input object
        domains.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>The computed reference domain is the union of all input object
        domains.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
</table>

\par Description
Each object is read from the given files in turn:
First from the files specified on the
commandline and then from the standard input.

The output object is a domain object of the same type
and with values of the same type of as the input objects.
All input objects must be the same type and have the same grey value type.
The output object has the same domain as the reference object,
which may either be computed or given on the command line.
The computed
domain is either the intersection or union of all input objects depending
on the commandline flags used.
The available filter functions are mean and rank,
where rank is a generalisation of minimum, median, maximum, etc.

For an output object \f$O\f$ having values at coordinate \f$w\f$
with \f$w \in \Omega\f$ where \f$\Omega\f$ is the domain of \f$O\f$,
the value at \f$w\f$ is
\f$O_w = f(I_{0,w}, I_{1,w}, I_{2,w}, \cdots)\f$,
where \f$\{I_0, I_1, I_2, \cdots)\f$ are the input objects
and \f$f\f$ is the filter function.

\par Examples
\verbatim
cat obj2.wlz | WlzFilterNObjsValues -i -r 0.5 obj0.wlz obj1.wlz >out.wlz
\endverbatim
Creates a new objectct, the domain of which is the intersection of the
domains of the input objects (obj0.wlz obj1.wlz obj2.wlz) and the values
are the median of the input object values.

\par File
\ref WlzFilterNObjValues.c "WlzFilterNObjValues.c"
\par See Also
\ref WlzFilterNObjValues "WlzFilterNObjValues(3)"
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
  int		idx,
		nObjs = 0,
		nFiles = 0,
		filterFn = 0, 			       /* 0 = rank, 1 = mean */
		refFn = 1, /* 0 = intersection, 1 = union, 2 = user supplied */
		objsMax = 0,
  		ok = 1,
  		option,
  		usage = 0;
  double	rank = 0.5;
  FILE		*fP = NULL;
  char		*fStr,
  		*outObjFileStr = NULL,
		*refObjFileStr = NULL;
  const char	*errMsgStr;
  WlzObject	*refObj = NULL,
  		*outObj = NULL;
  WlzObject	**objs = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "himuf:o:r:";
  const char    outObjFileStrDef[] = "-";
  const int	objsStep = 1024;

  opterr = 0;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'f':
	refFn = 2;
        refObjFileStr = optarg;
	break;
      case 'm':
        filterFn = 1;
	break;
      case 'r':
        filterFn = 0;
        break;
      case 'i':
	refFn = 0;
	break;
      case 'u':
	refFn = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok)
  {
    if((outObjFileStr == NULL) || (*outObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
  }
  if((objs = AlcMalloc(objsStep * sizeof(WlzObject *))) == NULL)
  {
    ok = 0;
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  /* Read objects from the command line. */
  if(ok) 
  {
    idx = 0;
    objsMax = objsStep;
    nFiles = argc - optind;
    while((errNum == WLZ_ERR_NONE) && (idx < nFiles))
    {
      if(nObjs >= objsMax)
      {
        objsMax += objsStep;
	if((objs = AlcRealloc(objs, objsStep * sizeof(WlzObject *))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(ok)
      {
        fStr = *(argv + optind + idx);
	if((fP = (strcmp(fStr, "-")?  fopen(fStr, "rb"): stdin)) == NULL)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
	else
	{
	  objs[nObjs] = WlzReadObj(fP, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    ++nObjs;
	  }
	  if(strcmp(fStr, "-"))
	  {
	    (void )fclose(fP);
	  }
	}
      }
      ++idx;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	       "%s: failed to read input objects from command line (%s).\n",
		     *argv, errMsgStr);
    }
  }
  /* Read objects from the standard input. */
  if(ok) 
  {
    while(errNum == WLZ_ERR_NONE)
    {
      if((nObjs + 1) >= objsMax)
      {
        objsMax += objsStep;
	if((objs = AlcRealloc(objs, objsStep * sizeof(WlzObject *))) == NULL)
	{
	  ok = 0;
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      objs[nObjs] = WlzReadObj(fP, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        ++nObjs;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      	     "%s: failed to read input objects from standard input (%s).\n",
		     *argv, errMsgStr);
    }
  }
  /* Compute or read the reference domain object. */
  if(ok)
  {
    switch(refFn)
    {
      case 0: /* intersection */
	refObj = WlzIntersectN(nObjs, objs, 0, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	  (void )fprintf(stderr,
		   "%s: failed to compute intersection (%s).\n",
			 *argv, errMsgStr);
	}
        break;
      case 1: /* union */
	refObj = WlzUnionN(nObjs, objs, 0, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	  (void )fprintf(stderr,
		   "%s: failed to compute union (%s).\n",
			 *argv, errMsgStr);
	}
        break;
      case 2: /* user supplied */
	if((fP = (strcmp(refObjFileStr, "-")?
	         fopen(refObjFileStr, "rb"): stdin)) == NULL)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
	else
	{
	  refObj = WlzReadObj(fP, &errNum);
	  if(strcmp(refObjFileStr, "-"))
	  {
	    (void )fclose(fP);
	  }
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	  (void )fprintf(stderr,
		   "%s: failed to read reference object from file %s (%s).\n",
			 *argv, refObjFileStr, errMsgStr);
	}
        break;
      default:
        break;
    }
  }
  if(ok)
  {
    outObj = WlzFilterNObjValues(refObj, nObjs, objs, filterFn, rank, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to filter objects (%s).\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outObjFileStr, "-")?
	     fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outObjFileStr);
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
  /* Clear up. */
  if(objs)
  {
    for(idx = 0; idx < nObjs; ++idx)
    {
      (void )WlzFreeObj(objs[idx]);
    }
    AlcFree(objs);
  }
  (void )WlzFreeObj(refObj);
  (void )WlzFreeObj(outObj);
  /* Report usage. */
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output file>] [-f <ref object>]\n"
    "                          [-m] [-r #] [-i] [-u]\n"
    "                          [<input objects>]\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -f  Reference object.\n"
    "  -m  Use mean filter.\n"
    "  -r  Use rank filter. The parameter (range [0.0-1.0]) specifies the\n"
    "      rank with: 0.0 minimum, 0.5 median and 1.0 maximum - intermediate\n"
    "      values specify intermediate ranks.\n"
    "  -i  The computed reference domain is the intersection of all input\n"
    "      object domains.\n"
    "  -u  The computed reference domain is the union of all input\n"
    "      object domains.\n"
    "\n"
    "      Each object is read from the given files in turn: First from\n"
    "      the files specified on the commandline and then from the\n"
    "      standard input.\n"
    "      The output object is a domain object of the same type and\n"
    "      with values of the same type of as the input objects.  All\n"
    "      input objects must be the same type and have the same grey\n"
    "      value type.  The output object has the same domain as the\n"
    "      reference object, which may either be computed or given on\n"
    "      the command line. The computed domain is either the intersection\n"
    "      or union of all input objects depending on the commandline\n"
    "      flags used.  The available filter functions are mean and rank,\n"
    "      where rank is a generalisation of minimum, median, maximum,\n"
    "      etc.\n"
    "      For an output object O having values at coordinate w with w\n"
    "      in Omega where Omega is the domain of O, the value at w is\n"
    "        O_w = f(I_{0,w}, I_{1,w}, I_{2,w}, ...)\n"
    "      where {I_0, I_1, I_2, ...} are the input objects and f is\n"
    "      the filter function.\n"
    "Example:\n"
    "  cat obj2.wlz | %s -i -r 0.5 obj0.wlz obj1.wlz >out.wlz\n"
    "Reads objects from obj0.wlz, obj1.wlz and obj2.wlz and filters them\n"
    "to create a new objectct, the domain of which is the intersection of\n"
    "the domains of the input objects and the values of which values are \n"
    "the median of the input object values.\n",
    argv[0], argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
