#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _lzCompoundArrayToScalar_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzCompoundArrayToScalar.c
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
* \brief	Reduces a compound array to a scalar valued object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*/

/*!
\ingroup BinWlz
\defgroup wlzcompoundarraytoscalar WlzCompoundArrayToScalar
\par Name
WlzCompoundArrayToScalar  -  reduces a compound array to a scalar valued
                             object.
\par Synopsis
\verbatim
WlzCompoundArrayToScalar [-h] [-o<output object>] [-l[ [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Computes the modulus of the compound array values (default).</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file (default is the standard output).</td>
  </tr>
</table>
\par Description
Reduces a compound array to an object with scalar values.
\par Examples
\verbatim
WlzCompoundArrayToScalar -o out.wlz ca.wlz
\endverbatim
Creates a new object which is written to the file out.wlz.
The output object is a domain object in which the values are the modulus
of the values in the compound array read from ca.wlz.
\par File
\ref WlzCompoundArrayToScalar.c "WlzCompoundArrayToScalar.c"
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
  		usage = 0;
  FILE		*fP = NULL;
  char 		*outObjFileStr;
  char		*inObjFileStr;
  WlzBinaryOperatorType bop = WLZ_BO_MODULUS;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj;
  WlzObject	*outObj = NULL;
  static char   optList[] = "lho:";
  const char    inObjFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  opterr = 0;
  inObj = NULL;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'l':
        bop = WLZ_BO_MODULUS;
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
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum),
				      NULL)) == NULL) ||
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
    outObj = (WlzObject *)WlzCompoundArrayToScalar((WlzCompoundArray *)inObj,
                                          bop, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr, "%s Failed to compute scalar object, %s.\n",
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
    else
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
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<output object>] [-l] [<in object>]\n"
            "Reduces a compound array to an object with scalar values.\n"
	    "Options are:\n"
	    "  -h  Output this usage message.\n"
	    "  -l  Modulus operator to convert to scalar (default).\n"
	    "  -o  Output file.\n"
	    "Example:\n"
	    "  %s -o out.wlz ca.wlz\n"
	    "Creates a new object which is written to the file out.wlz. The\n"
	    "output object is a domain object in which the values are the\n"
	    "modulus of the values in the compound array read from ca.wlz.\n",
	    argv[0], argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
