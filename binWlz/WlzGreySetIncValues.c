#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreySetIncValues_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGreySetIncValues.c
* \author       Bill Hill
* \date         February 2008
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
* \brief	Sets incrementing (and hence unique) values in scan
* 		order throughout a domain object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzgreysetincvalues "WlzGreySetIncValues"
*/

/*!
\ingroup BinWlz
\defgroup wlzgreysetincvalues  WlzGreySetIncValues
\par Name
WlzGreySetIncValues  -  sets incrementing (and hence unique) values in scan
                        order throughout a domain object.
\par Synopsis
\verbatim
WlzGreySetIncValues [-bh] [-o<output file>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-b</b></td>
    <td>Use Hilbert ranking.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file (default is the standard output).</td>
  </tr>
</table>
\par Description
Sets incrementing (and hence unique) values in scan order throughout a domain
object. The values are of type WLZ_GREY_INT.
\par Examples
\verbatim
WlzGreySetIncValues -o out.wlz in.wlz
\endverbatim
Creates a new object which is written to the file out.wlz.
This output object has grey values from 1 upwards incrementing throuhout it in
scan line order,
\par File
\ref WlzGreySetIncValues.c "WlzGreySetIncValues.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>


/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		hilbert = 0,
  		ok = 1,
  		option,
  		usage = 0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  static char   optList[] = "bho:";
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'b':
        hilbert = 1;
	break;
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
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
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
    if(hilbert)
    {
      outObj = WlzGreyNewHilbertRankValues(inObj, &errNum);
    }
    else
    {
      outObj = WlzGreyNewIncValues(inObj, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed to set object values, %s.\n",
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
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
            "Usage: %s [-b] [-h] [-o<output file>] [<input file>]\n"
            "Sets incrementing (and hence unique) values in scan order\n"
	    "throughout the domain of the input object.\n"
	    "Version: %s\n"
	    "Options:\n"
	    "  -b  Use Hilbert ranking.\n"
	    "  -h  Output this usage message.\n"
	    "  -o  Output file.\n",
	    argv[0],
	    WlzVersion());

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
