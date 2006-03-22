#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzContourFlipOrient_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzContourFlipOrient.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Flips orientation of edges and faces in contour models.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzcontourfliporient "WlzContourFlipOrient"
*/


/*!
\ingroup      BinWlz
\defgroup     wlzcontourfliporient WlzContourFlipOrient
\par Name
WlzContourFlipOrient - flips orientation of edges and faces in contour models.
\par Synopsis
\verbatim
WlzContourFlipOrient [-h] [-o<output file>] [<input file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>
Reads from standard input and writes to standard output by default.

\par Description
Flips the orientation of edges in 2D contour models and
faces in 3D contour models.

\par Examples
\verbatim
WlzContourFlipOrient -o out.wlz lobster.wlz
\endverbatim
Reverses the orientation of the faces in the 3D contour model lobster.wlz
and writes the output to out.wlz.

\par File
\ref WlzContourFlipOrient.c "WlzContourFlipOrient.c"
\par See Also
\ref wlzcontourobj "WlzContourObj(1)"
\ref WlzGMFilterFlipOrient "WlzGMFilterFlipOrient(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  static char   optList[] = "ho:";
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
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
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if(obj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    if(obj->type != WLZ_CONTOUR)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if(obj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      errNum = WlzGMFilterFlipOrient(obj->domain.ctr->model);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to flip orientation (%s).\n",
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
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output file>] [<input file>]\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "Flips the orientation of a contour's edges if 2D or faces if 3D\n",
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
