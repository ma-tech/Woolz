#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzInsideDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzInsideDomain.c
* \author       Bill Hill
* \date         May 2004
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
* \brief	Determines whether a vertex is within an object's domain.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzinsidedomain "WlzInsideDomain"
*/

/*!
\ingroup BinWlz
\defgroup wlzinsidedomain WlzInsideDomain
\par Name
WlzInsideDomain - determines whether a vertex is within an object's domain.
\par Synopsis
\verbatim
WlznsideDomain [-x#] [-y#] [-z#] [-h] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>Column position, default 0.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Line position, default 0.</td>
  </tr>
  <tr> 
    <td><b>-z</b></td>
    <td>Plane position, default 0.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
</table>
\par Description
Establishes whether the given vertex is inside or outside the domain
of the given object's domain. If inside 1 is output otherwise 0, with
either of these digits being followed by a new line character.
\par Examples
\verbatim
WlzInsideDomain -x 100 -y 100 toucan.wlz
\endverbatim
Outputs 1 if the vertex with coordiantes (100,100) is inside the domain
of the object read from the.wlz, otherwise outputs 0.
\par File
\ref WlzInsideDomain.c "WlzInsideDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreyvalue "WlzGreyValue(1)"
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
  		inside = 0,
  		ok = 1,
		usage = 0;
  WlzDVertex3	pos;
  char		*outFileStr,
  		*inObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  const char	*errMsg;
  static char	optList[] = "o:x:y:z:h",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  pos.vtX = 0;
  pos.vtY = 0;
  pos.vtZ = 0;
  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outFileStr = outFileStrDef;
  while(ok && (usage == 0) &&
  	((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'x':
        if(sscanf(optarg, "%lg", &(pos.vtX)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'y':
        if(sscanf(optarg, "%lg", &(pos.vtY)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'z':
        if(sscanf(optarg, "%lg", &(pos.vtZ)) != 1)
	{
	  usage = 1;
	  ok = 0;
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
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
                     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    inside = WlzInsideDomain(inObj, pos.vtZ, pos.vtY, pos.vtX, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to establish inside or outside (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"):
              stdout)) == NULL) ||
       (fprintf(fP, "%d\n", inside != 0) != 2))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output file (%s).\n",
                     *argv, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: "
    "%s [-x#] [-y#] [-z#] [-h] [<in object>]\n"
    "Establishes whether the given vertex is inside or outside the domain\n"
    "of the given object's domain. If inside 1 is output otherwise 0, with\n"
    "either of these digits being followed by a new line character\n"
    "Options are:\n"
    "  -x#   Column position (set to %g).\n"
    "  -y#   Line positio (set to %g).\n"
    "  -z#   Plane position (set to %g).\n"
    "  -o#   Output file name.\n"
    "  -h    Display this usage information.\n",
    *argv,
    pos.vtX, pos.vtY, pos.vtZ);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
