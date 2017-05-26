#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTensorComponent_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzTensorComponent.c
* \author       Bill Hill
* \date         April 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Extract tensor component values.
* \ingroup	wlztensorcomponent "WlzTensorComponent"
*/


/*!
\ingroup BinWlz
\defgroup wlztensorcomponent WlzTensorComponent
\par Name
WlzTensorComponent - extracts tensor component values.
\par Synopsis
\verbatim
WlzTensorComponent [-c <cpt>] [-o<output object>] [-h] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-c</b></td>
    <td>Component index which must be in the range \f$[0-(1-V_n)]\f$
        where \f$V_n\f$ is the number of component values per element
	in the object, eg for a scalar image 0 and for a rank 2 3x3 tensor
	this would be 9.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Extracts tensor component values from a tensor image object.
Because the resulting object is always a non-tiled object this
filter may also be used to convert a tiled to non-tiled object.
It is safe to supply a non-tiled object as the input (provided that
the component index is zero).
By default all files are read from the standard input and written to
the standard output.
\par Examples
\verbatim
WlzTensorComponent -o c_0.wlz -c 0 in.wlz
\endverbatim
Creates a scalar object containing just the first component of
the given tensor image. Alternatively this would also convert
a tiled object (in.wlz) to a non-tiled object (retaining only
the zeroth component of it's values).
The object is written to the file c_0.wlz
\par File
\ref WlzTensorComponent.c "WlzTensorComponent.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzTensorGetComponent "WlzTensorGetComponent(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
		usage = 0,
		cpt = 0;
  char		*inFileStr,
		*outFileStr;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hc:o:";
  const char    fileStrDef[] = "-";

  /* Parse the argument list and check for input files. */
  opterr = 0;
  inFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  while(ok && !usage && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'c':
	if(sscanf(optarg, "%d", &cpt) != 1)
        {
	  usage = 1;
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(ok)
  {
    ok = !usage;
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
      inFileStr = argv[optind];
    }
  }
  /* Read the input object. */
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read input object from file %s.\n",
                     argv[0], inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Extract component to a new object. */
  if(ok)
  {
    outObj = WlzTensorGetComponent(inObj, cpt, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to extract component with index %d, %s.\n",
                     argv[0], cpt, errMsgStr);
    }
  }
  /* Output the component object. */
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outFileStr);
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
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-c<cpt>] [-o<output object>] [-h]\n"
    "\t\t[<input object>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -c  Component index which must be in the range [0-(1-V_n)] where\n"
    "      V_n is the number of component values per element in the object,\n"
    "      eg for a scalar image 0 and for a rank 2 3x3 tensor this would\n"
    "      be 9.\n"
    "  -o  Output file.\n"
    "  -h  Help, prints usage message.\n"
    "Extracts tensor component values from a tensor image object.\n"
    "Because the resulting object is always a non-tiled object this\n"
    "filter may also be used to convert a tiled to non-tiled object.\n"
    "It is safe to supply a non-tiled object as the input (provided that\n"
    "the component index is zero). By default all files are read from the\n"
    "standard input and written to the standard output.\n"
    "Example:\n"
    "  %s -o c_0.wlz -c 0 in.wlz\n"
    "creates a scalar object containing just the first component of\n"
    "the given tensor image. Alternatively this would also convert\n"
    "a tiled object (in.wlz) to a non-tiled object (retaining only\n"
    "the zeroth component of it's values). The object is written to the\n"
    "file c_0.wlz\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
