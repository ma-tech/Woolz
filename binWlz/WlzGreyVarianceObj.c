#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreyVarianceObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzGreyVarianceObj.c
* \author       Bill Hill
* \date         February 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2008 Medical research Council, UK.
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
* \brief	Variance filters domain objects with values.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*/

/*!
\ingroup BinWlz
\defgroup wlzgreyvarianceobj WlzGreyVarianceObj
\par Name
WlzGreyVarianceObj  -  variance filters domain objects with values.
\par Synopsis
\verbatim
WlzGreyVarianceObj [-h] [-o<output file>] [-k#] [-s#] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-k</b></td>
    <td>Kernel size, must be 3, 5 or 7 (default 3).</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Variance scale factor which is applied to the computed variance
        before converting it to an integral value.</td>
  </tr>
</table>
\par Description
A grey value variance filter which computes the variance at each pixel
within and object by considering the values in a small kernel surrounding
the pixel.
\par Examples
\verbatim
WlzGreyVarianceObj -o out.wlz -k 5 in.wlz
\endverbatim
Creates a new object which is written to the file out.wlz.
This object is computed by applying a variance filter,
with kernel size 5x5 to
the object read from the file in.wlz.
\par File
\ref WlzGreyVarianceObj.c "WlzGreyVarianceObj.c"
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
  int		ok = 1,
  		option,
  		usage = 0,
		kSz = 3;
  double 	scale = 1.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  static char   optList[] = "ho:k:k:";
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
      case 'k':
	if(sscanf(optarg, "%d", &kSz) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
	if(sscanf(optarg, "%lg", &scale) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
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
    (void )WlzGreyVariance(obj, 0, kSz, scale, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed to filter object, %s.\n",
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
    fprintf(stderr,
            "Usage: %s [-h] [-o<output file>] [-r#] [-s#] [<input file>]\n"
    	    "Variance filters the grey values of a Woolz domain object.\n"
            "This computes the variance at each pixel within and object by\n"
	    "considering the values in a small kernel surrounding the pixel.\n"
	    "Options are:\n"
	    "  -h  Output this usage message.\n"
	    "  -o  Output file.\n"
	    "  -k  Size of filter kernel size, must be 3, 5 or 7.\n"
	    "      Default 3 for 3x3 region.\n"
	    "  -s  Variance scale factor which is applied to the computed\n"
	    "      variance before converting it to an integral value.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
