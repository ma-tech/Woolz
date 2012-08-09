#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCopyObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCopyObj.c
* \author       Bill Hill
* \date         December 2009
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
* \brief	Reads an object and writes it out again. Useful for moving
* 		to a new file format.
* \ingroup	BinWlz
*
* \par Binary
* \ref 		wlzcopyobj "WlzCopyObj"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcopyobj WlzCopyObj
\par Name
WlzCopyObj - reads an object and writes it out again which can be useful for
             moving to a new file format.
\par Synopsis
\verbatim
WlzCopyObj [-h] [-o<output file>] [<input file> [... <input file>]]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
</table>
\par Description
Reads objects and writes them out again which can be useful for
moving to a new file format.
\par Examples
\verbatim
WlzCopyObj -o new.wlz old.wlz
\endverbatim
\verbatim
WlzCopyObj <old.wlz >old.wlz
\endverbatim
Both these examples copy the object(s) from the file old.wlz to the file
out.wlz.
\par File
\ref WlzCopyObj.c "WlzCopyObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static WlzErrorNum 		WlzCopyObj(
				  FILE *outFP,
				  const char *inFile);

int		main(int argc, char *argv[])
{
  int		idx,
  		ok = 1,
  		option,
  		usage = 0;
  FILE		*fP = NULL;
  char		*inFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "ho:";
  const char    inFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  opterr = 0;
  inFileStr = (char *)inFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = usage == 0;
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
    for(idx = 0; (errNum == WLZ_ERR_NONE) && (optind + idx < argc); ++idx)
    {
      inFileStr = *(argv + optind + idx);
      errNum = WlzCopyObj(fP, inFileStr);
    }
    if((errNum == WLZ_ERR_NONE) && (idx == 0))
    {
      errNum = WlzCopyObj(fP, inFileStr);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to copy object from file %s, %s.\n",
      		     argv[0], inFileStr, errMsgStr);
    }
  }
  if((fP != NULL) && (strcmp(outFileStr, "-")))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<out file>] [<in file> [... <in file>]]\n"
            "Reads objects and writes them out again which can be useful for\n"
            "moving to a new file format.\n"
	    "Version: %s\n"
	    "Options:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -o  Output file.\n"
            "Examples:\n"
	    "  %s -o new.wlz old.wlz\n"
	    "  %s <old.wlz >new.wlz\n"
            "Both these examples copy the object(s) from the file old.wlz to\n"
	    "the file out.wlz.\n",
	    argv[0], WlzVersion(), argv[0], argv[0]);

  }
  return(!ok);
}

/*!
* \return	Woolz error code.
* \ingroup	BinWlz
* \brief	Reads objects from the file with the given name and writes
* 		them to the already opened output file.
* \param	outFP			Output file pointer.
* \param	inFile			Input file string, which may use "-"
* 					to specify the standard input.
*/
static WlzErrorNum WlzCopyObj(FILE *outFP, const char *inFile)
{
  WlzObject	*obj = NULL;
  FILE		*inFP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

    if((inFile == NULL) ||
       (*inFile == '\0') ||
       ((inFP = (strcmp(inFile, "-")?  fopen(inFile, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(inFP, &errNum), NULL)) == NULL))
    {
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
    while((errNum == WLZ_ERR_NONE) && (obj != NULL))
    {
      errNum = WlzWriteObj(outFP, obj);
      (void )WlzFreeObj(obj); obj = NULL;
      if(errNum == WLZ_ERR_NONE)
      {
        obj = WlzReadObj(inFP, &errNum);
	if(errNum == WLZ_ERR_READ_EOF)
	{
	  obj = NULL;
	  errNum = WLZ_ERR_NONE;
	}
      }
    }
    if(inFP && strcmp(inFile, "-"))
    {
      (void )fclose(inFP);
    }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
