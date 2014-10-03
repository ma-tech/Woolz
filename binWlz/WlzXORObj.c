#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzXORObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzXORObj.c
* \author       Bill Hill
* \date         March 2014
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
* \brief	Morphological set exclusive or operation on a pair of Woolz
* 		objects..
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzxorobj "WlzXORObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzimagearithmetic WlzImageArithmetic
\par Name
WlzXORObj - computes the set exclusive or of the two objects.
\par Synopsis
\verbatim
WlzXORObj [-o<out file>] [-h] [<in file 0>] [<in file 1>]
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
Computes the set exclusive or of the two given objects.
The input objects are read from stdin and values are written
to stdout unless the filenames are given.
\par Examples
\verbatim
WlzXORObj -o out.wlz in0.wlz in1.wlz
\endverbatim
The exclusive of the the domains objects read from in0.wlz and in1.wlz is
computed and written to out.wlz.
\par File
\ref WlzXORObj.c "WlzXORObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzxordom "WlzXORDom(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <Wlz.h>

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*oObj = NULL;
  WlzObject	*iObj[2];
  char 		*oObjFileStr;
  char  	*iObjFileStr[2];
  const char	*errMsg;
  static char	optList[] = "o:h",
		fileStrDef[] = "-";

  opterr = 0;
  iObj[0] = NULL;
  iObj[1] = NULL;
  oObjFileStr = fileStrDef;
  iObjFileStr[0] = fileStrDef;
  iObjFileStr[1] = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
	oObjFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if((iObjFileStr[0] == NULL) || (*iObjFileStr[0] == '\0') ||
     (iObjFileStr[1] == NULL) || (*iObjFileStr[1] == '\0') ||
     (oObjFileStr == NULL) || (*oObjFileStr == '\0'))
  {
    usage = 1;
  }
  if((usage == 0) && (optind < argc))
  {
    int 	i;

    i = 0;
    while((i < 2) && (optind < argc))
    {
      iObjFileStr[i] = *(argv + optind);
      ++optind;
      ++i;
    }
  }
  if((usage == 0) && (optind != argc))
  {
    usage = 1;
  }
  ok = !usage;
  if(ok)
  {
    int		i;

    i = 0;
    while((errNum == WLZ_ERR_NONE) && (i < 2))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((iObjFileStr[i] == NULL) ||
	  (*iObjFileStr[i] == '\0') ||
	  ((fP = (strcmp(iObjFileStr[i], "-")?
		  fopen(iObjFileStr[i], "r"): stdin)) == NULL) ||
	  ((iObj[i]= WlzAssignObject(WlzReadObj(fP,
	  					   &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to read object %d from file %s (%s).\n",
		       *argv, i, iObjFileStr[i], errMsg);
      }
      if(fP && strcmp(iObjFileStr[i], "-"))
      {
	fclose(fP);
      }
      ++i;
    }
  }
  if(ok)
  {
    oObj = WlzXORDom(iObj[0], iObj[1], &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to compute set xor object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(oObjFileStr, "-")? fopen(oObjFileStr, "w"):
	      				 stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, oObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(oObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  WlzFreeObj(iObj[0]);
  WlzFreeObj(iObj[1]);
  WlzFreeObj(oObj);
  if(usage)
  {

    (void )fprintf(stderr,
	"Usage: %s"
	" [-o <out file>] [-h] [<in file 0>] [<in file 1>]\n"
	"Version: %s\n"
	"Options:\n"
	"  -o        Output file name.\n"
	"  -h        Help, prints this usage message.\n"
	"Computes the set exclusive or of the two given objects.\n"
	"The input objects are read from stdin and values are written\n"
	"to stdout unless the filenames are given.\n"
	"Example:\n%s -o out.wlz in0.wlz in1.wlz\n"
	"The exclusive of the the domains objects read from in0.wlz\n"
	"and in1.wlz is computed and written to out.wlz\n",
	*argv,
	WlzVersion(),
	*argv);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
