#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzAutoCorrelate_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzAutoCorrelate.c
* \author       Bill Hill
* \date         November 2007
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2007 Medical research Council, UK.
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
* \brief	Computes the autocorrelation of the given object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*/

/*!
\ingroup BinWlz
\defgroup wlzautocorrelate WlzAutoCorrelate
\par Name
WlzAutoCorrelate - computes the autocorrelation of the given object.
\par Synopsis
\verbatim
WlzAutoCorrelate [-h] [-o<output file>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
</table>
\par Description
Computes the autocorrelation of the given object.
\par Examples
\verbatim
WlzAutoCorrelate -o out.wlz in.wlz
\endverbatim
Creates a new object which is written to the file out.wlz.
This object is the autocorrelation of the object read from the file in.wlz.
\par File
\ref WlzAutoCorrelate.c "WlzAutoCorrelate.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgauss "WlzGauss(1)"
\ref wlzrsvfilterobj "WlzRsvFilterObj(1)"
\ref WlzRankFilter "WlzRankFilter(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>


/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

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
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
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
    outObj = WlzAutoCor(inObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed objects autocorrelation, %s.\n",
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
    fprintf(stderr,
            "Usage: %s [-h] [-o<output file>] [<input object>]\n"
            "Computes the autocorrelation of the given object.\n"
	    "Options are:\n"
	    "  -h  Output this usage message.\n"
	    "  -o  Output file.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
