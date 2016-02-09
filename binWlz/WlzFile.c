#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFile_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzFile.c
* \author       Bill Hill
* \date         December 2008
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
* \brief	Identifies the type of object in a file.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzfile "WlzFile"
*/

/*!
\ingroup BinWlz
\defgroup wlzfile WlzFile
\par Name
WlzFile  -  identifies the type of the first object in a file.
\par Synopsis
\verbatim
WlzFile  [-h] [<input object file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
WlzFile reads the type of a woolz object from the  standard  input  or
       the  given  input  file without reading the entire object and
       writes the type of the first object in the file to the stadard
       output.
\par Examples
\verbatim
WlzFile toucan.wlz
\endverbatim
Prints the type od the object read in toucan.wlz (WLZ_2D_DOMAINOBJ)
to the standard error output.
\par File
\ref WlzFile.c "WlzFile.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzfacts "WlzFacts(1)"
\ref WlzReadObjType "WlzReadObjType(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int           option,
                ok = 1,
                usage = 0;
  WlzObjectType	objType = WLZ_NULL;
  char          *objFileStr;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  FILE          *fP = NULL;
  const char    *errMsg,
  		*objTypeStr;
  static char   optList[] = "h",
                objFileStrDef[] = "-";

  opterr = 0;
  objFileStr = objFileStrDef;
  while(ok && (usage == 0) &&
        ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'h':
      default:
        usage = 1;
        break;
    }
  }
  if((objFileStr == NULL) || (*objFileStr == '\0'))
  {
    usage = 1;
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      objFileStr = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if(strcmp(objFileStr, "-") == 0)
    {
      objType = WlzReadObjType(stdin, &errNum);
    }
    else
    {
      if((fP = fopen(objFileStr, "r")) != NULL)
      {
        objType = WlzReadObjType(fP, &errNum);
	(void )fclose(fP);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      objTypeStr = WlzStringFromObjTypeValue(objType, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object type from file %s (%s).\n",
                     *argv, objFileStr, errMsg);
    }
  }
  if(ok)
  {
    (void )printf("%s\n", objTypeStr);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: "
    "%s [-h] [<in object>]\n"
    "Reads the type of the first Woolz object in the given file without\n"
    "reading the entire object and writes the corresponding object type\n"
    "string to the standard output.\n"
    "By default the input file is the standard input.\n"
    "Version %s\n"
    "Options:\n"
    "  -h    Display this usage information.\n",
    *argv,
    WlzVersion());
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
