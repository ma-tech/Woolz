#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzIndexObjToCompound_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzIndexObjToCompound.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Creates an grey value index object from the given compound
* 		array object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzindexobjtocompound "WlzIndexObjToCompound"
*/

/*!
\ingroup BinWlz
\defgroup wlzindexobjtocompound WlzIndexObjToCompound
\par Name
WlzIndexObjToCompound - creates a grey value index object from a compound
		          array object..
\par Synopsis
\verbatim
WlzIndexObjToCompound [-h] [-o<output compound object>]
                      [<input index object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name, default standard output.</td>
  </tr>
</table>
\par Description
Reads a grey values index object and uses it to create a a new
compound array object in which each object of the array is either
an empty object or a spatial domain without grey values. The
spatial domain objects will be of the same type as the given
object and their domains will correspond those of the grey
values of the given object, where the i'th object of the compound
object's array will have the domain for which all grey values have
the value i.  This function can be considered the inverse of
This binary can be considered the inverse of WlzIndexObjFromCompound.
\par File
\ref WlzIndexObjToCompound.c "WlzIndexObjToCompound.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzcompound WlzCompound(1)
\ref wlzexplode WlzExplode(1)
\ref wlzindexobjfromcompound WlzIndexObjFromCompound(1)
\ref WlzIndexObjToCompound WlzIndexObjToCompound(3)
\ref WlzIndexObjFromCompound WlzIndexObjFromCompound(3)
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0;
  FILE		*fP = NULL;
  char		*iFile,
  		*oFile;
  const char	*errMsg;
  WlzObject	*iObj = NULL,
  		*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "ho:";
  const char    fileDef[] = "-";

  opterr = 0;
  iFile = (char *)fileDef;
  oFile = (char *)fileDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        oFile = optarg;
	break;
      case 'h': /* FALLTROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    if((oFile == NULL) || (*oFile == '\0') ||
       (iFile == NULL) || (*iFile == '\0'))
    {
      usage = 1;
    }
  }
  if(!usage && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      iFile = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    if((*iFile == '\0') ||
       ((fP = (strcmp(iFile, "-")?  fopen(iFile, "r"): stdin)) == NULL) ||
       ((iObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
    }
    if(fP)
    {
      if(strcmp(iFile, "-"))
      {
        (void )fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    if(strcmp(iFile, "-"))
    {
      if((fP = fopen(iFile, "r")) == NULL)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open input file %s (%s).\n",
		       *argv, iFile, errMsg);
      }
    }
    else
    {
      fP = stdin;
    }
  }
  if(ok)
  {
    oObj = (WlzObject *)WlzIndexObjToCompound(iObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to create compound object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(oFile, "-")?  fopen(oFile, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, oObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object to file %s (%s).\n",
                     *argv, oFile, errMsg);
    }
    if(fP && strcmp(oFile, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(iObj);
  (void )WlzFreeObj(oObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output index object>] [<input compound object>]\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "Reads a grey values index object and uses it to create a a new\n"
    "compound array object in which each object of the array is either\n"
    "an empty object or a spatial domain without grey values. The\n"
    "spatial domain objects will be of the same type as the given\n"
    "object and their domains will correspond those of the grey\n"
    "values of the given object, where the i'th object of the compound\n"
    "object's array will have the domain for which all grey values have\n"
    "the value i.  This function can be considered the inverse of\n"
    "This binary can be considered the inverse of WlzIndexObjFromCompound.\n",
    argv[0]);
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
