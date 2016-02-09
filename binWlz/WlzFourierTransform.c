#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFourierTransform_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzFourierTransform.c
* \author       Bill Hill
* \date         December 2012
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
* \brief	Computes the Fourier transform (or the inverse transform)
* 		of an object.
* \ingroup	BinWlz
* \par Binary
* \ref		wlzfouriertransform "WlzFourierTransform"
*/

/*!
\ingroup BinWlz
\defgroup wlzfouriertransform WlzFourierTransform
\par Name
WlzFourierTransform - computes the Fourier transform (or it's inverse)
                      of an object.
\par Synopsis
\verbatim
WlzFourierTransform [-h] [-i] [-o<output object>] [input object]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Computes the inverse transform.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
</table>
\par Description
Computes the Fourier transform (or it's inverse) of an spatial domain
object with scalar values.
\par Examples
\verbatim
WlzFourierTransform in.wlz | WlzFourierTransform -i -o out.wlz -
\endverbatim
Computes the Fourier transform of the object in in.wlz then through
a pipe computes the inverse transform writing the output to out.wlz.
\par File
\ref WlzFourierTransform.c "WlzFourierTransform.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzFourierTransformObj "WlzFourierTransformObj(3)"
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
		invFlg = 0,
		ok = 1,
		usage = 0;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*inObjFileStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  static char	optList[] = "hio:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outObjFileStr = outFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'i':
        invFlg = 1;
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
    if(optind < argc)
    {
      if((optind + 1) == argc)
      {
	inObjFileStr = *(argv + optind);
      }
      else
      {
	usage = 1;
	ok = 0;
      }
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
    if(!ok)
    {
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsgStr);
    }
  }
  if(ok)
  {
    outObj = WlzFourierTransformObj(inObj, !invFlg, &errNum);
    if((outObj == NULL) || (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s: Failed to compute Fourier transform (%s).\n",
                     *argv, errMsgStr);
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outObjFileStr, "-")?
	      fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
	(WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		     "%s: Failed to write output object (%s).\n",
		     *argv, errMsgStr);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s%s%s",
    *argv,
    " [-h] [-i] [-o<output object>] [input object]\n" 
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -i  Computes the inverse transform.\n"
    "  -o  Output file name.\n"
    "Computes the Fourier transform (or it's inverse) of an spatial domain\n"
    "object with scalar values. If the input object is a compound object\n"
    "then the domains must be identical and the values are treated as\n"
    "complex values (with the first component being real).\n",
    *argv,
    " in.wlz | ",
    *argv,
    " -i -o out.wlz -\n"
    "Computes the Fourier transform of the object in in.wlz then through\n"
    "a pipe computes the inverse transform writing the output to out.wlz\n.");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
