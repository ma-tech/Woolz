#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSepFilterObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzSepFilterObj.c
* \author       Bill Hill
* \date         November 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Seperable spatial convolution kernel base filters
* 		for 2 and 3D objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzsepfilterobj "WlzSepFilterObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzsepfilterobj WlzSepFilterObj
\par Name
WlzSepFilterObj - seperable spatial convolution filter for domain objects
		  with grey values.
\par Synopsis
\verbatim
WlzSepFilterObj - [-h] [-v] [-o<output object>] [-m#,#,#] [-n#,#,#]
                  [-g <t>] [-t <flags>] [-P <p>] [-p#] [-G] [-S a|c|s]
		  [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose output with parameter values and execution time.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Filter parameters; for a Gaussian these are the sigma
        values for each orthogonal direction (x, y and z). Default
	values are 1,1,1.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Order parameters; for a Gaussian these are the order of the
        Gaussian's derivatives (eg 0 for no derivative, 1 for 1st
        derivative), default 0,0,0.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Required grey type, default is the same as the input object's,
        with 'i', 's', 'u', 'f' and 'd' used to request int, short,
	unsigned byte, float or double grey values.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Filter directions, eg x filter along lines and x, y filter
        along lines then through columns.</td>
  </tr>
  <tr> 
    <td><b>-P</b></td>
    <td>Type of padding to use, default is to use the background value,
        with 'b', 'e', 'v' and 'z' used to select background, end, given
        value and zero padding.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Used to supply given padding value, default 0.0.</td>
  </tr>
  <tr> 
    <td><b>-G</b></td>
    <td>Use Gaussian filter, default.</td>
  </tr>
  <tr>
    <td><b>-S</b></td>
    <td>Use the input object for each of the directional filters rather
        than the output of the previous directional filter. The seperate
	filter outputs are then combined to form a compound object (<b>c</b>),
	the sum (<b>a</b>) or the square root of the sum of the squared
	values of the filter passes (<b>s</b>).</td>
  </tr>
</table>
\par Description
Applies a seperable filter to a Woolz domain object with grey values.
The input object is read from stdin and the filtered object is written
to stdout unless the filenames are given.
\par Examples
\verbatim
WlzSepFilterObj -m ,3.0,2.0 -t x,y,z -o smooth.wlz in.wlz
\endverbatim
The input Woolz object is read from in.wlz, and filtered using a
Gaussian filter with sigma values of 1.0, 3.0 and 2.0 in the x, y and z
directions. The filtered object is then written to the file smoothed.wlz.
\par File
\ref WlzSepFilterObj.c "WlzSepFilterObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgauss "WlzGauss(1)"
\ref WlzRsvFilterObj "WlzRsvFilterObj(1)"
\ref WlzGaussFilter "WlzGaussFilter(3)"
\ref WlzSepFilter "WlzSepFilter(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sys/time.h>
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
		sep = 0,
		usage = 0,
		verbose = 0;
  double	padVal = 0.0;
  WlzDVertex3	param[] = {{1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}};
  WlzIVertex3	order = {0, 0, 0};
  AlgPadType	pad = ALG_PAD_NONE;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzIVertex3	action = {1, 1, 1};
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  struct timeval times[3];
  const char    *errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*outFileStr,
  		*inObjFileStr;
  static char	optList[] = "Ghvg:o:m:n:P:p:t:S:",
		defFile[] = "-";

  opterr = 0;
  outFileStr = defFile;
  inObjFileStr = defFile;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
        usage = 1;
	break;
      case 'v':
        verbose = 1;
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
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'm':
	if((optarg == NULL) ||
	   (sscanf(optarg, "%lg,%lg,%lg",
	           &(param[0].vtX), &(param[0].vtY), &(param[0].vtZ)) < 1))
	{
	  usage = 1;
	}
        break;
      case 'n':
	if((optarg == NULL) ||
	   (sscanf(optarg, "%d,%d,%d",
	           &(order.vtX), &(order.vtY), &(order.vtZ)) < 1))
	{
	  usage = 1;
	}
        break;
      case 'P':
        switch(*optarg)
	{
	  case 'b':
	    pad = ALG_PAD_NONE;
	    break;
	  case 'e':
	    pad = ALG_PAD_END;
	    break;
	  case 'v':
	    pad = ALG_PAD_VALUE;
	    break;
	  case 'z':
	    pad = ALG_PAD_ZERO;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'p':
        if((optarg == NULL) || (sscanf(optarg, "%lg", &padVal) != 1))
	{
	  usage = 1;
	}
	break;
      case 'G':
	/* Only have Gaussian implemented. */
        break;
      case 't':
	if(optarg)
	{
	  char	*actStr;

	  action.vtX = 0;
	  action.vtY = 0;
	  action.vtZ = 0;
	  actStr = strtok(optarg, ",");
	  while(actStr && ok)
	  {
	    if(strcmp(actStr, "x") == 0)
	    {
	      action.vtX = 1;
	    }
	    else if(strcmp(actStr, "y") == 0)
	    {
	      action.vtY = 1;
	    }
	    else if(strcmp(actStr, "z") == 0)
	    {
	      action.vtZ = 1;
	    }
	    else
	    {
	      usage = 1;
	    }
	    actStr = strtok(NULL, ",");
	  }
	}
	break;
      case 'S':
	if(!optarg || (strlen(optarg) > 1))
	{
	  usage = 1;
	}
	else
	{
	  switch(*optarg)
	  {
	    case 'c':
	      sep = 1;
	      break;
	    case 'a':
	      sep = 2;
	      break;
	    case 's':
	      sep = 3;
	      break;
	    default:
	      usage = 1;
	      break;
	  }
	}
        break;
      default:
        usage = 1;
	break;
    }
  }
  if(!usage)
  {
    if(optind < argc)
    {
      if((optind + 1) != argc)
      {
	usage = 1;
      }
      else
      {
	inObjFileStr = *(argv + optind);
      }
    }
  }
  ok = !usage;
  if(ok)
  {
    FILE	*fP = NULL;

    if(verbose)
    {
      (void )fprintf(stderr,
                     "%s: Reading object from %s\n",
		     *argv, strcmp(inObjFileStr, "-")? inObjFileStr: "stdin");
    }
    if(((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok && (pad == ALG_PAD_NONE))
  {
    WlzPixelV	bgd;

    bgd = WlzGetBackground(inObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValueConvertPixel(&bgd, bgd, WLZ_GREY_DOUBLE);
      padVal = bgd.v.dbv;
      pad = ALG_PAD_VALUE;
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to get object background value (%s).\n",
		     *argv, errMsgStr);
    }
  }
  if(ok)
  {
    if(verbose)
    {
      (void )fprintf(stderr,
                     "%s: Filtering object with:\n"
		     "  param[0] = (%g, %g, %g)\n"
		     "  param[1] = (%g, %g, %g)\n"
		     "  order    = (%d, %d, %d)\n"
		     "  action   = (%d, %d, %d)\n"
		     "  gType    = %d\n"
		     "  pad      = %d\n"
		     "  padVal   = %lg\n"
		     "  sep      = %d\n",
		     *argv,
		     param[0].vtX, param[0].vtY, param[0].vtZ,
		     param[1].vtX, param[1].vtY, param[1].vtZ,
		     order.vtX, order.vtY, order.vtZ,
		     action.vtX, action.vtY, action.vtZ,
		     gType, pad, padVal, sep);
      gettimeofday(times + 0, NULL);
    }
    outObj = WlzAssignObject(
	     WlzGaussFilter(inObj, param[0], order, action, gType,
	     		    pad, padVal, sep, &errNum),
	     NULL);
    if(verbose)
    {
      gettimeofday(times + 1, NULL);
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: Elapsed time for WlzGaussFilter() %gus\n",
		     argv[0], (1.0e06 * times[2].tv_sec) + times[2].tv_usec);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to filter object (%s).\n",
		     *argv, errMsgStr);
    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

    if(verbose)
    {
      (void )fprintf(stderr,
                     "%s: Writing object to %s\n",
		     *argv, strcmp(outFileStr, "-")? outFileStr: "stdout");
    }
    errNum = WLZ_ERR_FILE_OPEN;
    if(((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		     "%s: Failed to write output object (%s).\n",
		     *argv,
		     errMsgStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-h] [-o<output object>] [-m#,#,#] [-n#,#,#]\n"
    "\t[-g <t>] [-t <flags>] [-P <p>] [-p#] [-G]\n"
    "\t[-S a|c|s] [-v] [<input object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Prints this usage information\n"
    "  -o  Output object file name.\n"
    "  -m  Filter parameters; for a Gaussian these are the sigma values\n"
    "      for each orthogonal direction (x, y and z). Default values are\n"
    "      1,1,1.\n"
    "  -n  Order parameters; for a Gaussian these are the order of the\n"
    "      Gaussian's derivatives (eg 0 for no derivative, 1 for 1st\n"
    "      derivative), default 0,0,0.\n"
    "  -g  Required grey type, default is the same as the input object's,\n"
    "      with 'i', 's', 'u', 'f' and 'd' used to request int, short,\n"
    "      unsigned byte, float or double grey values.\n"
    "  -P  Type of padding to use, default is to use the background value,\n"
    "      with 'b', 'e', 'v' and 'z' used to select background, end, given\n"
    "      value and zero padding.\n"
    "  -p  Used to supply given padding value, default 0.0.\n"
    "  -G  Use a Gaussain filter (default).\n"
    "  -S  Use the input object for each of the directional filters rather\n"
    "      than the output of the previous directional filter. The seperate\n"
    "      filter outputs are then combined to form a compound object (c),\n"
    "      the sum (a) or the square root of the sum of the squared values\n"
    "      of the filter passes (s).\n"
    "  -t  Filter directions, eg x filter along lines and x, y filter\n"
    "      along lines then through columns.\n"
    "  -v  Verbose output with parameter values and execution time.\n"
    "Applies a seperable filter to a Woolz domain object with grey values.\n"
    "The input object is read from stdin and the filtered object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -m ,3.0,2.0 -t x,y,z -o smooth.wlz in.wlz\n"
    "The input Woolz object is read from in.wlz, and filtered using a\n"
    "Gaussian filter with sigma values of 1.0, 3.0 and 2.0 in the x, y\n"
    "and z directions. The filtered object is then written to the file\n"
    "smoothed.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
