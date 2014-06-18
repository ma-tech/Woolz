#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRasterObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzRasterObj.c
* \author       Bill Hill
* \date         March 2001
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
* \brief	Rasterizes geometric objects into domain objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzrasterobj "WlzRasterObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzrasterobj WlzRasterObj
\par Name
WlzRasterObj  -  Scan converts a geometric object into a domain object.
\par Synopsis
\verbatim
WlzRasterObj [-o<output object>] [-h] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
</table>
\par Description
Scan converts the given geometric object into a domain object.
The input object is read from stdin and the output output is written
to stdout unless filenames are given.
\par Examples
\verbatim
WlzRasterObj in.wlz
\endverbatim
The input geometric object is read from in.wlz, and the scan converted
object is written to stdout.
\par File
\ref WlzRasterObj.c "WlzRasterObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzpolygontoobj "WlzPolyToObj(1)"
\ref WlzRasterObj "WlzRasterObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <Wlz.h>
#include <limits.h>
#include <string.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           option,
  		ok = 1,
		usage = 0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  static char	optList[] = "ho:";
  const char	*errMsgStr;
  const char	outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  outObjFileStr = (char *)outObjFileStrDef;
  inObjFileStr = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
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
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
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
    fP = NULL;
  }
  if(ok)
  {
    outObj = WlzRasterObj(inObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed to scan convert object (%s).\n",
      	     	     argv[0], errMsgStr);
    }
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(ok)
  {
    if((fP = (strcmp(outObjFileStr, "-")?
	     fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Failed to open output file %s.\n",
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
                     "%s: Failed to write output object (%s).\n",
                     argv[0], errMsgStr);
    }
  }
  if(fP && strcmp(outObjFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s%s%sExample: %s%s",
      *argv,
      " [-o<output object>] [-h] [<input object>]\n"
      "Version: ",
      WlzVersion(),
      "\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -o  Output object file name.\n"
      "Scan converts the given geometric object into a domain object.\n"
      "The input object is read from stdin and the output output is written\n"
      "to stdout unless filenames are given.\n",
      *argv,
      " in.wlz\n"
      "The input Woolz object is read from in.wlz, and the scan converted\n"
      "object is written to stdout.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
