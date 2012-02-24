#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCompDispIncGrey_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCompDispIncGrey.c
* \author       Bill Hill
* \date         February 2008
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
* \brief	Computes a displacement field as a compound array object
* 		from two given objects. The displacements are between
* 		pixels/voxels with the same grey value.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcompdispincgrey "WlzCompDispIncGrey"
*/

/*!
\ingroup BinWlz
\defgroup wlzcompdispincgrey WlzCompDispIncGrey
\par Name
WlzCompDispIncGrey  -  computes the displacements between the pixels/voxels of
		       that have the same value from a pair of objects.
\par Synopsis
\verbatim
WlzCompDispIncGrey [-h] [-o<output file>] [<input file 1>] [<input file 2>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file (default is the standard output).</td>
  </tr>
</table>
\par Description
Computes a displacement field as a compound array object from two given
objects. The displacements are between pixels/voxels with the same grey
value. If the values of one of the objects were set using WlzGreySetIncValues
and the other object is the result of a spatial transformation of that object
then the displacements field will capture the transformation. The dispplacement
values map object 1 to object 2.
\par Examples
\verbatim
WlzCompDispIncGrey -o out.wlz obj.wlz obj_t.wlz
\endverbatim
Creates a new object which is written to the file out.wlz.
The output object is a compound array containing the displacements from
obj.wlz to obj1.wlz.
\par File
\ref WlzCompDispIncGrey.c "WlzCompDispIncGrey.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref BinWlz "WlzGreySetIncValues(1)"
\ref BinWlz "WlzCompoundArrayToScalar(1)"
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
  int		idN,
  		ok = 1,
  		option,
  		usage = 0;
  FILE		*fP = NULL;
  char 		*outObjFileStr;
  char		*inObjFileStr[2];
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj[2];
  WlzObject	*outObj = NULL;
  static char   optList[] = "ho:";
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  opterr = 0;
  inObj[0] = inObj[1] = NULL;
  inObjFileStr[0] = inObjFileStr[1] = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
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
    if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
       (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
       (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 2) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inObjFileStr[0] = *(argv + optind + 0);
        inObjFileStr[1] = *(argv + optind + 1);
      }
    }
  }
  for(idN = 0; (ok != 0) && (idN < 2); ++idN)
  {
    if((inObjFileStr[idN] == NULL) ||
       (*inObjFileStr[idN] == '\0') ||
       ((fP = (strcmp(inObjFileStr[idN], "-")?
	      fopen(inObjFileStr[idN], "r"): stdin)) == NULL) ||
       ((inObj[idN] = WlzAssignObject(WlzReadObj(fP, &errNum),
				      NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr[idN]);
    }
    if(fP && strcmp(inObjFileStr[idN], "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    outObj = (WlzObject *)WlzCompDispIncGrey(inObj[0], inObj[1], &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed to compute displacement values, %s.\n",
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
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<output object>] [<object 1>] [<object 2>]\n"
	    "Computes a displacement field as a compound array object from\n"
	    "two given objects. The displacements are between pixels/voxels\n"
	    "with the same grey value. If the values of one of the objects\n"
	    "were set using WlzGreySetIncValues and the other object is the\n"
	    "result of a spatial transformation of that object then the\n"
	    "displacements field will capture the transformation.\n"
	    "The dispplacement values map object 1 to object 2.\n"
	    "Options are:\n"
	    "  -h  Output this usage message.\n"
	    "  -o  Output file.\n"
	    "Example:\n"
            "  %s -o out.wlz obj.wlz obj_t.wlz\n"
	    "Creates a new object which is written to the file out.wlz. The\n"
	    "output object is a compound array containing the displacements\n"
	    "from obj.wlz to obj1.wlz.\n",
	    argv[0], argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
