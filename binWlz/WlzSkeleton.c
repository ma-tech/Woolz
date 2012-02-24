#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSkeleton_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzSkeleton.c
* \author       Bill Hill
* \date         November 1998
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
* \brief	Computes the skeleton of the a domain object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzskeleton "WlzSkeleton"
*/

/*!
\ingroup BinWlz
\defgroup wlzskeleton WlzSkeleton
\par Name
WlzSkeleton - computes the skeleton of the a domain object.
\par Synopsis
\verbatim
WlzSkeleton [-h] [-o<out object>] [-c #] [-s #] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Minimum connectivity (4, 6, 8, 18, 26), default 8.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Smoothing passes, default 1XXX.</td>
  </tr>
</table>
\par Description
Computes the skeleton of a 2D domain object.
Objects are read from stdin and written to stdout unless the filenames
are given.
\par Examples
\verbatim
WlzSkeleton -o skel.wlz myobj.wlz
\endverbatim
The input Woolz object is read from myobj.wlz, filtered and written
to skel.wlz.
\par File
\ref WlzSkeleton.c "WlzSkeleton.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzerosion "WlzErosion(1)"
\ref wlzdilation "WlzDilation(1)"
\ref WlzSkeleton "WlzSkeleton(3)"
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
		cValI,
		smPass,
		ok = 1,
		usage = 0;
  WlzConnectType con;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	defCValI = 8,
  		defSmPass = 1;
  const char    *errMsg;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  static char	optList[] = "o:c:s:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  cValI = defCValI;
  smPass = defSmPass;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'c':
	if(sscanf(optarg, "%d", &cValI) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 's':
	if(sscanf(optarg, "%d", &smPass) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
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
  if(ok)
  {
    switch(cValI)
    {
      case 4:
	con = WLZ_4_CONNECTED;
	break;
      case 6:
	con = WLZ_6_CONNECTED;
	break;
      case 8:
	con = WLZ_8_CONNECTED;
	break;
      case 18:
	con = WLZ_18_CONNECTED;
	break;
      case 26:
	con = WLZ_26_CONNECTED;
	break;
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if((outObj = WlzSkeleton(inObj, smPass, con, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute skeleton (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%d%s%d%sExample: %s%s",
    *argv,
    " [-o<out object>] [-c #] [-s #] [-h] [<in object>]\n"
    "Options:\n"
    "  -o  Output object file name.\n"
    "  -c  Minimum connectivity (4, 6, 8, 18, 26), default ",
    defCValI,
    ".\n"
    "  -s  Smoothing passes, default ",
    defSmPass,
    ".\n"
    "  -h  Help, prints this usage message.\n"
    "Computes the skeleton of a 2D domain object.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o skel.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, filtered and written\n"
    "to skel.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
