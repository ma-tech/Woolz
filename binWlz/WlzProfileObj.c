#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzProfileObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzProfileObj.c
* \author       Bill Hill
* \date         September 2018
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
* \brief	Extracts an arbitrary line of values from the given object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzprofileobj "WlzProfileObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzprofileobj WlzProfileObj
\par Name
WlzProfileObj - Extracts an arbitrary line of values from the given object.
\par Synopsis
\verbatim
WlzProfileObj [-h] [-o<output object file>]
           [-s<x>,<y>,<z>] [-e<x>,<y>,<z>] [<input object file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file</td>
  </tr>
  <tr> 
    <td><b>-e</b></td>
    <td>End of line segment.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Start of line segment.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation (not enabled).</td>
  </tr>
</table>
By  default  the  input  object  is read from the standard
input and the output object is  written  to  the  standard output.
\par Description
Extracts a single line 2D domain object (with values if they exist in the
given object) from the given object. The line segment runs from the given
start to the given end point. 
The components values of the start and end points of the line segment
have defaults of the minimum and maximum bounding box values
respectively.
\par Examples
Extracts a line of values from (0,100) to (200,300) through the given
2D object read from in.wlz.
The output object is written to the file out.wlz.
\verbatim
WlzProfileObj -o out.wlz -s 0,100 -e 200,300 in.wlz
\endverbatim

\par File
\ref WlzProfileObj.c "WlzProfileObj.c"
\par See Also
\ref WlzProfileLine "WlzProfileLine(3)"
\ref WlzProfileLineIDom "WlzProfileLineIDom(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0;
  WlzIVertex3	startPos,
  		endPos,
		startSet = {0},
		endSet = {0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "e:o:s:h",
		fileStrDef[] = "-";

  outObjFileStr = fileStrDef;
  inObjFileStr = fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'e': /* FALLTHROUGH */
      case 's':
	if(optarg)
	{
	  char	*sep[2];
	  WlzIVertex3 set = {0},
	  	      pos = {0};

	  (void )WlzStringWhiteSpSkip(optarg);
	  sep[0] = index(optarg, ',');
	  if(sep[0])
	  {
	    *(sep[0]) = '\0';
	  }
	  set.vtX = (sscanf(optarg, "%d", &(pos.vtX)) == 1);
	  if(sep[0])
	  {
	    sep[1] = index(sep[0] + 1, ',');
	    if(sep[1])
	    {
	      *(sep[1]) = '\0';
	    }
	    set.vtY = (sscanf(sep[0] + 1, "%d", &(pos.vtY)) == 1);
	    if(sep[1])
	    {
	      set.vtZ = (sscanf(sep[1] + 1, "%d", &(pos.vtZ)) == 1);
	    }
	  }
	  if(option == 'e')
	  {
	    endPos = pos;
	    endSet = set;
	  }
	  else /* option == 's' */
	  {
	    startPos = pos;
	    startSet = set;
	  }
	}
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
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
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzReadObj(fP, &errNum)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if(!(startSet.vtX || startSet.vtY || startSet.vtZ ||
         endSet.vtX || endSet.vtY || endSet.vtZ))
    {
      WlzIBox3 box;
      
      box = WlzBoundingBox3I(inObj, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0; 
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: failed to find bounding box of input object(%s).\n",
		       *argv, errMsg);
      }
      else
      {
        startPos.vtX = (startSet.vtX)? startPos.vtX: box.xMin;
        startPos.vtY = (startSet.vtY)? startPos.vtY: box.yMin;
        startPos.vtZ = (startSet.vtZ)? startPos.vtZ: box.zMin;
        endPos.vtX = (endSet.vtX)? endPos.vtX: box.xMax;
        endPos.vtY = (endSet.vtY)? endPos.vtY: box.yMax;
        endPos.vtZ = (endSet.vtZ)? endPos.vtZ: box.zMax;
      }
    }
  }
  if(ok)        
  {
    if(((outObj = WlzProfileLine(inObj, startPos, endPos, &errNum)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed extract profile object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write profile object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-h] [-o<out object>]\n"
    "       [-e<x>,<y>,<z>] [-s<x>,<y>,<z>] [<input object file>]\n"
    "Extracts a single line 2D domain object (with values if they exist in\n"
    "the given object) from the given object. The line segment runs from the\n"
    "given start to the given end point.\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -e  End of line segment.\n"
    "  -o  Output object file name.\n"
    "  -s  Start of line segment.\n"
    "The components values of the start end and end points of the line\n"
    "segment have defaults of the minimum and maximum bounding box\n"
    "values respectively.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o out.wlz -s 0,100 -e 200,300 in.wlz\n"
    "Extracts a line of values from (0,100) to (200,300) through the given\n"
    "2D object read from in.wlz. The output object is written to the file\n"
    "out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
