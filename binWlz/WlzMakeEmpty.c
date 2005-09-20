#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzMakeEmpty.c
* \author       Bill Hill
* \date         September 2005
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
* \brief	Makes an empty woolz object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzmakeempty "WlzMakeEmpty"
*/

/*!
\ingroup BinWlz
\defgroup wlzmakeempty WlzMakeEmpty
\par Name
WlzMakeEmpty - makes an empty woolz object.
\par Synopsis
\verbatim
WlzIsEmpty [-h] [-o <output file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
</table>
\par Description
Makes an empty woolz object.
\par Examples
\verbatim
WlzMakeEmpty >empty.wlz
\endverbatim
An empty object is written to the file empty.wlz.
\par File
\ref WlzMakeEmpty.c "WlzMakeEmpty.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzisempty "WlzIsEmpty(1)"
\ref WlzMakeEmpty() "WlzMakeEmpty(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
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
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  FILE		*fP = NULL;
  char		*outFileStr;
  const char	*errMsg;
  static char	optList[] = "ho:";

  opterr = 0;
  outFileStr = "-";
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default: /* FALLTHROUGH */
        usage = 1;
	break;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
  }
  ok = !usage;
  /* Create the empty object. */
  if(ok)
  {
    obj = WlzMakeEmpty(&errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	             "%s: Failed to create empty object (%s)\n",
		     *argv, errMsg);

    }
  }
  /* Open the output file. */
  if(ok)
  {
    if((outFileStr == NULL) ||
       (*outFileStr == '\0') ||
       ((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     *argv, outFileStr);
    }
  }
  /* Write out the object. */
  if(ok)
  {
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write emptty object to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-o<output file>]\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file.\n"
    "Makes an empty object. By default the object is written to the standard\n"
    "output.\n",
    *argv,
    " >empty.wlz\n"
    "Writes an empty object to the file empty.wlz.\n");
  }
  return((ok)? 0: 1);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
