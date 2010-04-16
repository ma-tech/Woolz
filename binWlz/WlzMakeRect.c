#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzMakeRect_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzMakeRect.c
* \author       Bill Hill
* \date         April 2010
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
* \brief	Makes a rectangular or cuboid domain.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref WlzMakeRect "WlzMakeRect"
*/

/*!
\ingroup BinWlz
\defgroup wlzmakerect WlzMakeRect
\par Name
WlzMakeRect - Makes a either a 2D or 3D domain object with either
a rectanglar domain or a plane domain with rectangular domains.
\par Synopsis
\verbatim
WlzMakeRect [-23] [-o<output file>] [-x<x min>,<x max>] [-y<y min>,<y max>]
            [-z<z min>,<z max>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-2</b></td>
    <td>Create a 2D rectanglar domain object (default).</td>
  </tr>
  <tr> 
    <td><b>-3</b></td>
    <td>Create a 3D cuboid domain object.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file</td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>Column domain limits (default 0,0).</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Line domain limits (default 0,0).</td>
  </tr>
  <tr> 
    <td><b>-z</b></td>
    <td>Plane domain limits (default 0,0 and makes the default dimension 3).</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints message.</td>
  </tr>
</table>
\par Description
Makes a either a 2D or 3D domain object with either
a rectanglar domain or a plane domain with rectangular domains.
By default an object is created for a single 2D pixel at the origin
and the output object is written to the standard output.
\par Examples
A Single voxel object with origin (0,0,0) is written to the
standard output:
\verbatim
WlzMakeRect -3
\endverbatim
A single pixel object with origin (1,2) is written to the
standard output.
\verbatim
WlzMakeRect -x 1 -y 2
\endverbatim
A 3D domain object with planes from 3 to 6, lines from 4 to 7 and
column from 5 to 8 is written to the file dom.wlz:
\verbatim
WlzMakeRect -x5,8, -y4,7 -z3,6 -o dom.wlz
\endverbatim  

\par File
\ref WlzMakeRect.c "WlzMakeRect.c"
\par See Also
\ref WlzMakeRect "WlzMakeRect(3)"
\ref WlzMakeCuboid "WlzMakeCuboid(3)"
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
  int		idx,
  		option,
		dim = 2,
		ok = 1,
		usage = 0;
  WlzIBox3	domBox;
  WlzPixelV	dumV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*obj = NULL;
  char 		*outFileStr;
  int		domVal[2];
  char		*domStr[2];
  const char	*errMsg;
  static char	optList[] = "23ho:x:y:z:",
		outFileStrDef[] = "-";

  opterr = 0;
  domBox.xMin = 0;
  domBox.xMax = 0;
  domBox.yMin = 0;
  domBox.yMax = 0;
  domBox.zMin = 0;
  domBox.zMax = 0;
  dumV.type = WLZ_GREY_UBYTE;
  dumV.v.ubv = 0;
  outFileStr = outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'x':
      case 'y':
      case 'z':
	if(optarg)
	{
	  while(*optarg && isspace(*optarg))
	  {
	    ++optarg;
	  }
	  if(*optarg == ',')
	  {
	    domStr[0] = NULL;
	    domStr[1] = strtok(optarg, ",");
	  }
	  else
	  {
	    domStr[0] = strtok(optarg, ",");
	    domStr[1] = strtok(NULL, ",");
	  }
	  if((domStr[0] == NULL) && (domStr[1] == NULL))
	  {
	    usage = 1;
	    usage = 1;
	  }
	  else
	  {
	    idx = 0;
	    while((usage == 0) && (idx < 2))
	    {
	      if(domStr[idx] && (sscanf(domStr[idx], "%d",
					 domVal + idx) != 1))
	      {
		usage = 1;
		usage = 1;
	      }
	      ++idx;
	    }
	  }
	}
	if(usage == 0)
	{
	  switch(option)
	  {
	    case 'x':
	      if(domStr[0])
	      {
		domBox.xMin = domVal[0];
	      }
	      if(domStr[1])
	      {
		domBox.xMax = domVal[1];
	      }
	      break;
	    case 'y':
	      if(domStr[0])
	      {
		domBox.yMin = domVal[0];
	      }
	      if(domStr[1])
	      {
		domBox.yMax = domVal[1];
	      }
	      break;
	    case 'z':
	      dim = 3;
	      if(domStr[0])
	      {
		domBox.zMin = domVal[0];
	      }
	      if(domStr[1])
	      {
		domBox.zMax = domVal[1];
	      }
	      break;
	  }
	}
	break;
      case 'h':
      default:
        usage = 1;
	break;
    }
  }
  if((usage == 0) &&
     ((outFileStr == NULL) || (*outFileStr == '\0')))
  {
    usage = 1;
  }
  if(ok && (optind != argc))
  {
    usage = 1;
  }
  ok = !usage;
  if(ok)
  {
    switch(dim)
    {
      case 2:
	obj = WlzMakeRect(domBox.yMin, domBox.yMax,
	                  domBox.xMin, domBox.xMax,
			  WLZ_GREY_ERROR, NULL, dumV, NULL, NULL,
			  &errNum);
        break;
      case 3:
	obj = WlzMakeCuboid(domBox.zMin, domBox.zMax,
			    domBox.yMin, domBox.yMax,
	                    domBox.xMin, domBox.xMax,
			    WLZ_GREY_ERROR, dumV, NULL, NULL,
			    &errNum);
        break;
      default:
        break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to make %s domain object (%s).\n",
		     *argv, (dim == 2)? "rectangular": "cuboid", errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, obj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write domain object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(obj)
  {
    WlzFreeObj(obj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExamples:\n%s%s%s%s%s%s",
    *argv,
    " [-23h] [-o<out object>]\n"
         "       [-x<x min>,<x max>] [-y<y min>,<y max>] [-z<z min>,<z max>]\n"
    "Makes a either a 2D or 3D domain object with either a rectanglar domain\n"
    "or a plane domain with rectangular domains. By default an object is\n"
    "created for a single 2D pixel at the origin and the output object is\n"
    "written to the standard output.\n"
    "Options:\n"
    "  -2  Create a 2D rectangular domain object (default).\n"
    "  -3  Create a 3D cuboid domain object.\n"
    "  -x  Column domain limits (default 0,0).\n"
    "  -y  Line domain limits (default 0,0).\n"
    "  -z  Plane domain limits (default 0,0 and makes the default dimension\n"
    "      3.\n",
    *argv,
    " -3\n"
    "A Single voxel object with origin (0,0,0) is written to the standard\n"
    "output.\n",
    *argv,
    " -x 1 -y 2\n"
    "A single pixel object with origin (1,2) is written to the standard\n"
    "output.\n",
    *argv,
    " -x5,8, -y4,7 -z3,6 -o dom.wlz\n"
    "A 3D domain object with planes from 3 to 6, lines from 4 to 7 and\n"
    "column from 5 to 8 is written to the file dom.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
