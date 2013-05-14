#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzNObjsGreyStats_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzNObjsGreyStats.c
* \author       Bill Hill
* \date         March 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Computes basic grey value statistics across multiple objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlznobjsgreystats "WlzNObjsGreyStats"
*/

/*!
\ingroup      BinWlz
\defgroup     wlznobjsgreystats  WlzNObjsGreyStats
\par Name
WlzNObjsGreyStats - computes basic grey value statistics across multiple
		    objects.

\par Synopsis
\verbatim
WlzNObjsGreyStats [-h] [-m] [-s] [-o<output file>] [<input objects>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-m</b></td>
    <td>Output mean not sum of object values.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Output standard deviation not sum of squares of object values.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file name, by default the output objecta are written to
        the standard out.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
</table>

\par Description
Reads objects (including compound objects) from the given files in turn:
First from the files specified on the
commandline and then from the standard input.

The output objects are written to a single compound object. Each object
of the compound object will have the same domain (which is the intersection
of the input object domains) and values for the simple statistics. The
statistics generated are min value, max value, sum of values and
sum of squares of values.

\par Examples
\verbatim
cat obj1.wlz obj2.wlz obj3.wlz | WlzNObjsGreyStats -s -m obj4.wlz >out.wlz
\endverbatim
Creates a compound object with four component objects, the domains of which
are all the intersection of the domains of the input objects (obj[1-4].wlz)
and the values are the minimum, maximum, mean and standard deviation of
the input object values across the objects within the intersection domain.

\par File
\ref WlzNObjsGreyStats.c "WlzNObjsGreyStats.c"
\par See Also
\ref WlzNObjGreyStats "WlzNObjGreyStats(3)"
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
  int		idC,
  		idO,
		idP,
  		option,
		ok = 1,
  		usage = 0,
		mean = 0,
		stddev = 0,
		nObjs = 0,
		maxObjs = 0;
  FILE		*fP = NULL;
  char		*fStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzObject	**objs = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "hmso:";
  const char    fileStrDef[] = "-";
  const int	stepObjs = 1024;

  opterr = 0;
  outObjFileStr = (char *)fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'm':
	mean = 1;
	break;
      case 's':
        stddev = 1;
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
  if(ok)
  {
    maxObjs += stepObjs;
    if((objs = AlcMalloc(maxObjs * sizeof(WlzObject *))) == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Read objects from the command line. */
  if(ok) 
  {
    int		nFiles;

    idO = 0;
    nFiles = argc - optind;
    while((errNum == WLZ_ERR_NONE) && (idO < nFiles))
    {
      if(nObjs >= maxObjs)
      {
        maxObjs += stepObjs;
	if((objs = AlcRealloc(objs, stepObjs * sizeof(WlzObject *))) == NULL)
	{
	  ok = 0;
	  errNum = WLZ_ERR_MEM_ALLOC;
	  break;
	}
      }
      if(ok)
      {
        fStr = *(argv + optind + idO);
	if((fP = (strcmp(fStr, "-")?  fopen(fStr, "rb"): stdin)) == NULL)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
	else
	{
	  objs[nObjs] = WlzAssignObject(
	  		WlzReadObj(fP, &errNum), NULL);
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
      ++idO;
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
  /* Create a single compound object from the input objects. */
  if(ok)
  {
    int		oCnt = 0,
    		oType = 0;

    for(idO = 0; idO < nObjs; ++idO)
    {
      if(objs[idO])
      {
        switch(objs[idO]->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	    ++oCnt;
	    if(oType == 0)
	    {
	      oType = 2;
	    }
	    else if(oType != 2)
	    {
	      errNum = WLZ_ERR_OBJECT_TYPE;
	    }
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    ++oCnt;
	    if(oType == 0)
	    {
	      oType = 3;
	    }
	    else if(oType != 3)
	    {
	      errNum = WLZ_ERR_OBJECT_TYPE;
	    }
	    break;
	  case WLZ_COMPOUND_ARR_1: /* FALLTHOUGH */
	  case WLZ_COMPOUND_ARR_2:
	    oCnt += ((WlzCompoundArray *)(objs[idO]))->n;
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      inObj = (WlzObject *)WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1,
      				1, oCnt, NULL, WLZ_NULL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      idC = 0;
      for(idO = 0; idO < nObjs; ++idO)
      {
	if(objs[idO])
	{
	  switch(objs[idO]->type)
	  {
	    case WLZ_2D_DOMAINOBJ: /* FALLTHOUGH */
	    case WLZ_3D_DOMAINOBJ:
	      ((WlzCompoundArray *)inObj)->o[idC++] =
	        WlzAssignObject(objs[idO], NULL);
	      break;
	    case WLZ_COMPOUND_ARR_1: /* FALLTHOUGH */
	    case WLZ_COMPOUND_ARR_2:
	      for(idP = 0; idP < ((WlzCompoundArray *)(objs[idO]))->n; ++idP)
	      {
		((WlzCompoundArray *)inObj)->o[idC++] =
		  WlzAssignObject(((WlzCompoundArray *)(objs[idO]))->o[idP],
		                  NULL);
	      }
	      break;
	    default:
	      break;
	  }
	}
      }
      ((WlzCompoundArray *)(inObj))->n = idC;
      if((((WlzCompoundArray *)(inObj))->n < 1) ||
         (((WlzCompoundArray *)(inObj))->o[0] == NULL))
      {
        errNum = WLZ_ERR_OBJECT_NULL;
      }
      else
      {
        ((WlzCompoundArray *)(inObj))->otype = ((WlzCompoundArray *)
	                                        (inObj))->o[0]->type;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	       "%s: Failed to create compound object of input object (%s).\n",
		     *argv, errMsgStr);
    }
  }
  if(objs)
  {
    for(idO = 0; idO < nObjs; ++idO)
    {
      (void )WlzFreeObj(objs[idO]);
    }
    AlcFree(objs);
    objs = NULL;
  }
  /* Create compound object for collecting thye statistics. */
  if(ok)
  {
    outObj = (WlzObject *)WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1,
      				1, 4, NULL, 
				((WlzCompoundArray *)(inObj))->otype,
				&errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	   "%s: Failed to allocate compound object for output object (%s).\n",
		     *argv, errMsgStr);
    }
  }
  /* Compute the cross object statistics. */
  if(ok)
  {
    errNum = WlzNObjGreyStats(inObj, mean, stddev, NULL,
                              &(((WlzCompoundArray *)outObj)->o[0]),
                              &(((WlzCompoundArray *)outObj)->o[1]),
                              &(((WlzCompoundArray *)outObj)->o[2]),
                              &(((WlzCompoundArray *)outObj)->o[3]));
    if(errNum == WLZ_ERR_NONE)
    {
      for(idO = 0; idO < 4; ++idO)
      {
        (void )WlzAssignObject(((WlzCompoundArray *)outObj)->o[idO], NULL);
      }
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      "%s: Failed to compute cross object statistics (%s).\n",
                     *argv, errMsgStr);
    }
  }
  /* Output cross object statistics compound object. */
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
    if(strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  /* Report usage. */
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-m] [-s] [-o<output file>] [<input objects>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -m  Output mean instead of sum of values.\n"
    "  -s  Output standard deviation instead of sum of squares of values.\n"
    "Reads objects (including compound objects) from the given files in\n"
    "turn: First from the files specified on the commandline and then from\n"
    "the standard input. The output objects are written to a single compound\n"
    "object. Each object of the compound object will have the same domain\n"
    "(which is the intersection of the input object domains) and values for\n"
    "the simple statistics. The statistics generated are min value, max\n"
    "value, sum of values and sum of squares of values.\n"
    "Example:\n"
    "  cat obj1.wlz obj2.wlz obj3.wlz | %s -s -m obj4.wlz >out.wlz\n"
    "Creates a compound object with four component objects, the domains of\n"
    "which are all the intersection of the domains of the input objects\n"
    "(obj[1-4].wlz) and the values are the minimum, maximum, mean and\n"
    "standard deviation of the input object values across the objects within\n"
    "the intersection domain.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
