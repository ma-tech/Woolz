#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzScalarScale_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzScalarScale.c
* \author       Bill Hill
* \date         April 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
* \brief	Creates an object with a tiled value table from an
* 		object with a valid spatial domain.
* \ingroup	BinWlz
*/

/*!
\ingroup BinWlz
\defgroup WlzScalarScale WlzScalarScale
\par Name
WlzScalarScale  - scales the grey values of an object.
\par Synopsis
\verbatim
WlzScalarScale  [-a #] [-m #] [-g #] [-h] [-o<output file>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>Value for addition.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Value for product.</td>
  </tr>
  <tr>
    <td><b>-g</b></td>
    <td>Grey type specified using one of the characters:
        i, s, u, f, d, r for int, short, unsigned byte, float, double or
	red-green-blue-alpha.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
WlzScalarScale scales the grey values of an object using a linear scaling
with \f$v_{out} = m v_{in} + a.\f$
\par Examples
\verbatim
WlzScalarScale -a 256 -s 2 -g s -o out.wlz in.wlz
\endverbatim
Creates a new object (out.wlz) with the same doamin as the input object
(in.wlz) but with short grey values which vale been scaled so that the
output grey values are given by: g_out = g_in * 2 + 256.
\par File
\ref WlzScalarScale.c "WlzScalarScale.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzScalarMulAdd "WlzScalarMulAdd(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <sys/time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		option,
  		ok = 1,
		usage = 0;
  WlzPixelV	a,
  		m;
  WlzGreyType	gType = WLZ_GREY_UBYTE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  char		*inFileStr,
		*outFileStr;
  const char	*errMsg;
  static char	optList[] = "a:g:hm:o:",
  		inFileStrDef[] = "-",
		outFileStrDef[] = "-";

  opterr = 0;
  a.type = WLZ_GREY_DOUBLE;
  a.v.dbv = 0.0;
  m.type = WLZ_GREY_DOUBLE;
  m.v.dbv = 1.0;
  inFileStr = inFileStrDef;
  outFileStr = outFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
	if(sscanf(optarg, "%lg", &(a.v.dbv)) != 1)
	{
	  usage = 1;
	}
	break;
      case 'g':
        switch(*optarg)
	{
	  case 'i':
	    gType = WLZ_GREY_INT;
	    break;
	  case 's':
	    gType = WLZ_GREY_SHORT;
	    break;
	  case 'u':
	    gType = WLZ_GREY_UBYTE;
	    break;
	  case 'f':
	    gType = WLZ_GREY_FLOAT;
	    break;
	  case 'd':
	    gType = WLZ_GREY_DOUBLE;
	    break;
	  case 'r':
	    gType = WLZ_GREY_RGBA;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'm':
	if(sscanf(optarg, "%lg", &(m.v.dbv)) != 1)
	{
	  usage = 1;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inFileStr = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    fP = NULL;
    errNum = WLZ_ERR_READ_EOF;
    if(((fP = (strcmp(inFileStr, "-")? fopen(inFileStr, "r"):
                                       stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    outObj = WlzScalarMulAdd(inObj, m, a, gType, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to create scaled object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok && (outFileStr != NULL))
  {
    if(((fP = (strcmp(outFileStr, "-")? fopen(outFileStr, "w"):
                                        stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write scaled object to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    if(fP!= NULL)
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(outObj);
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-a#] [-g#] [-h] [-o<output object>] [-m #] [<input object>]\n"
    "Scales the grey values of an object using a linear scaling with\n"
    "v_out = m v_{in} + a.\n"
    "Options:\n"
    "  -a  Value for addition.\n"
    "  -g  Grey type specified using one of the characters:\n"
    "      i, s, u, f, d, r for int, short, unsigned byte, float, double\n"
    "      or RGBA.\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output tiled object.\n"
    "  -m Value for product.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
