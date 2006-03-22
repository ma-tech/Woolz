#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRankObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzRankObj.c
* \author       Bill Hill
* \date         March 2002
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
* \brief	Rank filters domain objects with values.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzrankobj "WlzRankObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzrankobj WlzRankObj
\par Name
WlzRankObj  -  rank filters domain objects with values.
\par Synopsis
\verbatim
WlzRankObj [-h] [-o<output file>] [-r#] [-s#] [<input file>]
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
  <tr> 
    <td><b>-r</b></td>
    <td>Required rank. Range [0.0-1.0] with 0.0 minimum, 0.5
        median and 1.0 maximum value. Default 0.5.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Size of filter region, must be greater than zero.
        Default 3 for 3x3 region.</td>
  </tr>
</table>
\par Description
Rank filters the grey values of a domain object.
By selecting the appropriate rank, the filter can be used as a
maximum, minimum or median filter.
\par Examples
\verbatim
WlzRankObj -o out.wlz -r 0.5 in.wlz
\endverbatim
Creates a new object which is written to the file out.wlz.
This object is computed by applying a median filter to
the object read from the file in.wlz.
\par File
\ref WlzRankObj.c "WlzRankObj.c"
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


extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0,
		rSz = 3;
  double 	rank = 0.5;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  static char   optList[] = "ho:r:s:";
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
      case 'r':
	if((sscanf(optarg, "%lg", &rank) != 1) || (rank < 0.0) || (rank > 1.0))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
	if((sscanf(optarg, "%d", &rSz) != 1) || (rSz < 0))
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
    errNum = WlzRankFilter(obj, rSz, rank);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed to rank filter object, %s.\n",
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
    	    "Rank filters the grey values of a Woolz domain object.\n"
	    "Options are:\n"
	    "  -h  Output this usage message.\n"
	    "  -o  Output file.\n"
	    "  -r  Required rank. Range [0.0-1.0] with 0.0 minimum, 0.5\n"
	    "      median and 1.0 maximum value. Default 0.5.\n"
	    "  -s  Size of filter region, must be greater than zero.\n"
	    "      Default 3 for 3x3 region.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
