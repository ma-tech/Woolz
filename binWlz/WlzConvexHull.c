#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConvexHull_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzConvexHull.c
* \author       Bill Hill
* \date         May 2014
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
* \brief	Computes convex hulls.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzconvexhull "WlzConvexHull"
*/

/*!
\ingroup BinWlz
\defgroup wlzconvexhull WlzConvexHull
\par Name
WlzConvexHull - computes the convex hull of the given input object.
\par Synopsis
\verbatim
WlzConvexHull [-o<output object>] [-h] [-t<output object type>] [-u]
              [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Reports usage information.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Output object type, which must be one of:
    </td>
  </tr>
  <tr> 
    <td><b>-u</b></td>
    <td>Ignore voxel size and use unit voxel size.</td>
  </tr>
</table>
\par Description
Computes the convex hull of the given input object.
Degenerate convex hulls (on a single plane for 3D or on a
single line for 2D are only allowed when creating pixel/voxel
spatial domain objects.
The input object is read from stdin and output data are written
to stdout unless filenames are given.
\par Examples
\verbatim
WlzConvexHull -o out.wlz -t t in.wlz
\endverbatim
The input Woolz object is read from in.wlz, and the contour model
corresponding to the convex hull of the input object is written to
the file out.wlz.
\par File
\ref WlzConvexHull.c "WlzConnvexHull.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzObjToConvexHull "WlzObjToConvexHull(3)"
\ref WlzConvexHullToObj "WlzConvexHullToObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int           option,
  		ok = 1,
		usage = 0,
		unitVoxelSz = 0;
  FILE		*fP = NULL;
  char		*inFileStr,
  		*outFileStr;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  WlzObjectType	outObjType = WLZ_CONV_HULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "huo:t:";
  const char	outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  outFileStr = (char *)outFileStrDef;
  inFileStr = (char *)inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 't':
	if(optarg && (strlen(optarg) == 1))
	{
	  switch(*optarg)
	  {
	    case 't':
	      outObjType = WLZ_CONTOUR;
	      break;
	    case 'v':
	      outObjType = WLZ_CONV_HULL;
	      break;
	    case 'x':
	      outObjType = WLZ_2D_DOMAINOBJ;
	      break;
	    default:
	      usage = 1;
	      break;
	  }
	}
	else
	{
	  usage = 1;
	}
      case 'u':
        unitVoxelSz = 1;
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inFileStr == NULL) || (*inFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
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
        inFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
	      fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s\n",
		     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      fclose(fP);
    }

  }
  if(ok && unitVoxelSz)
  {
    if(inObj &&
       (inObj->type == WLZ_3D_DOMAINOBJ) &&
       inObj->domain.core &&
       (inObj->domain.core->type == WLZ_PLANEDOMAIN_DOMAIN))
    {
      inObj->domain.p->voxel_size[0] = 1.0;
      inObj->domain.p->voxel_size[1] = 1.0;
      inObj->domain.p->voxel_size[2] = 1.0;
    }
  }
  if(ok)
  {
    outObj = WlzAssignObject(
    	     WlzObjToConvexHull(inObj, &errNum), NULL);
    if((outObj != NULL) &&
       (errNum == WLZ_ERR_DEGENERATE) && (outObjType == WLZ_2D_DOMAINOBJ))
    {
      errNum = WLZ_ERR_NONE;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to compute convex hull (%s).\n",
      	     	     argv[0], errMsgStr);
    }
  }
  if(ok && (outObjType != WLZ_CONV_HULL))
  {
    WlzObject	*tObj;

    if(outObj &&
       (outObj->type == WLZ_CONV_HULL) &&
       outObj->domain.core &&
       (outObj->domain.core->type == WLZ_CONVHULL_DOMAIN_3D) &&
       (outObjType == WLZ_2D_DOMAINOBJ))
    {
      outObjType = WLZ_3D_DOMAINOBJ;
    }
    tObj = WlzAssignObject(
    	   WlzConvexHullToObj(outObj, outObjType, &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	         "%s: Failed convert convex hull to the required type (%s).\n",
		 argv[0], errMsgStr);
    }
    if(tObj)
    {
      (void )WlzFreeObj(outObj);
      outObj = tObj;
    }
  }
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
    errNum = WlzWriteObj(fP, outObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to write output object (%s).\n",
      		     argv[0], errMsgStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-o<output object>] [-h] [-o] [-t <type>] [-u] [<input object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output object file name.\n"
    "  -t  Output object type given as a single character:\n"
    "        t - contour.\n"
    "        v - convex hull (default)\n"
    "        x - (pixel/voxel) spatial domain\n"
    "  -u  Use unit voxel size.\n"
    "Computes the convex hull of the given input object.\n"
    "Degenerate convex hulls (on a single plane for 3D or on a\n"
    "single line for 2D are only allowed when creating pixel/voxel\n"
    "spatial domain objects.\n"
    "The input object is read from stdin and output data are written\n"
    "to stdout unless filenames are given.\n",
    *argv,
    " -o out.wlz -t t in.wlz\n"
    "The input Woolz object is read from in.wlz, and the contour\n"
    "model corresponding to the convex hull of the input object is\n"
    "written to the file out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
