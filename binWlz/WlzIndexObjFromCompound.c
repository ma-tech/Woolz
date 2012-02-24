#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzIndexObjFromCompound_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzIndexObjFromCompound.c
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
* \ref wlzindexobjfromcompound "WlzIndexObjFromCompound"
*/

/*!
\ingroup BinWlz
\defgroup wlzindexobjfromcompound WlzIndexObjFromCompound
\par Name
WlzIndexObjFromCompound - creates a grey value index object from a compound
		          array object..
\par Synopsis
\verbatim
WlzIndexObjFromCompound [-h] [-o<output index object>]
                        [<input compound object>]
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
Reads a compound array object, which must contain only 2D or 3D
domain objects or empty objects and uses the domain objects to
create a new spatial domain object with grey values; the grey
values being indices of the domains in the given compound array
object.

The domain of the new object is the union of the domains of the
objects in the given compound array object. The values are set so
that all values in the domain of the i'th compound array object are
set to i, with i incremented from 0 through to one less than the
number of objects in the compound array object. If the domains
overlap then higher index objects will overwrite lower index objects
within their intersection.
This binary can be considered the inverse of WlzIndexObjToCompound.
\par File
\ref WlzIndexObjFromCompound.c "WlzIndexObjFromCompound.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzcompound WlzCompound(1)
\ref wlzindexobjtocompound WlzIndexObjToCompound(1)
\ref WlzIndexObjFromCompound WlzIndexObjFromCompound(3)
\ref WlzIndexObjToCompound WlzIndexObjToCompound(3)
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
    oObj = WlzIndexObjFromCompound((WlzCompoundArray *)iObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to create index object (%s).\n",
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
    "Reads a compound array object, which must contain only 2D or 3D\n"
    "domain objects or empty objects and uses the domain objects to\n"
    "create a new spatial domain object with grey values; the grey\n"
    "values being indices of the domains in the given compound array\n"
    "object.\n"
    "The domain of the new object is the union of the domains of the\n"
    "objects in the given compound array object. The values are set so\n"
    "that all values in the domain of the i'th compound array object are\n"
    "set to i, with i incremented from 0 through to one less than the\n"
    "number of objects in the compound array object. If the domains\n"
    "overlap then higher index objects will overwrite lower index objects\n"
    "within their intersection.\n"
    "This binary can be considered the inverse of WlzIndexObjToCompound.\n",
    argv[0]);
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
