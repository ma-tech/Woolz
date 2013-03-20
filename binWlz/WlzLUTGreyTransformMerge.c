#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLUTGreyTransformMerge_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzLUTGreyTransformMerge.c
* \author       Bill Hill
* \date         November 2011
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
* \brief	Merges integer look up table grey transforms
* 		to form a RGB\f$\alpha\f$ grey transform.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzlutgreytransformmerge "WlzLUTGreyTransformMerge"
*/
 
/*!
\ingroup      BinWlz
\defgroup     wlzlutgreytransformmerge WlzLUTGreyTransformMerge
\par Name
WlzLUTGreyTransformMerge -
Merges integer look up table grey transforms to form a RGB\f$\alpha\f$ grey
transform.
\par Synopsis
\verbatim
WlzLUTGreyTransformMerge [-o<output file>] [-h]
           [<red LUT>] [<green LUT>] [<blue LUT>] [<alpha LUT>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints a usage message.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>output object file name.</td>
  </tr>
</table>

\par Description
Reads integer look up table objects given on the command line and
merges them to form a red,green,blue,\f$\alpha\f$ look up table
object. Any of the look up table object file names can be replaced
by the string \verbatim NULL \endverbatim. Any of the input files
may be read from the standard input. If no files are given on the
command line all files will be read from the standard input, otherwise
a \verbatim - \endverbatim in place of a file name may be used to specify
standard input. When reading from the standard input an empty object
may be used as well as integer look up table objects, with these read
in the order red, green, blue, alpha (the same as the command line order).
When a NULL or empty object is used it implies an identity look up table.
By default the merged look up table is written to the standard output.

\par Examples
\verbatim
  # Merges three look up tables red.wlz, green.wlz and alpha.wlz
  # together with and identity look up table for the blue channel
  # to form a RGBA look up table object that is written to rgba.wlz.
  cat green.wlz | WlzLUTGreyTransformMerge -o rgba.wlz red.wlz - NULL alpha.wlz
\endverbatim.
\par File
\ref WlzLUTGreyTransformMerge.c "WlzLUTGreyTransformMerge.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzfacts "WlzFacts(1)"
\ref WlzLUTGreyTransformNew "WlzLUTGreyTransformNew(3)"
\ref WlzLUTMergeToRGBA "WlzLUTMergeToRGBA(3)"
\ref WlzLUTTransformObj "WlzLUTTransformObj(1)"
\ref WlzMakeEmpty "WlzMakeEmpty(1)"
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
  int		idx,	
		nFiles,
		option,
		ok = 1,
		usage = 0;
  WlzObject	*outObj = NULL;
  WlzObject	*inObj[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  const char	*errMsg;
  char 		*outFileStr;
  char		*inFileStr[4];
  static char	optList[] = "ho:",
		outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObj[0] = inObj[1] = inObj[2] = inObj[3] = NULL;
  inFileStr[0] = inFileStr[1] = inFileStr[2] = inFileStr[3] = inFileStrDef;
  /* Parse the command line. */
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
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
  if(usage == 0)
  {
    if((outFileStr == NULL) || (*outFileStr == '\0'))
    {
      usage = 1;
    }
    if(usage == 0)
    {
      if((nFiles = argc - optind) < 0)
      {
	usage = 1;
      }
      else
      {
	for(idx = 0; idx < nFiles; ++idx)
	{
	  char *fStr;
	  fStr = *(argv + optind + idx);
	  if(fStr)
	  {
	    inFileStr[idx] = fStr;
	  }
	}
      }
    }
  }
  ok = (usage == 0);
  /* Read if input LUT files and check that they are valid. */
  if(ok)
  {
    for(idx = 0; idx < 4; ++idx)
    {
      if(strcmp(inFileStr[idx], "NULL"))
      {
	errNum = WLZ_ERR_READ_EOF;
	if(((fP = (strcmp(inFileStr[idx], "-")?
	          fopen(inFileStr[idx], "r"): stdin)) == NULL) ||
		  ((inObj[idx] = WlzAssignObject(
		                 WlzReadObj(fP, &errNum), NULL)) == NULL))
        {
	  ok = 0;        
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	                 "%s: failed to read object from file %s (%s)\n",
			 *argv, inFileStr[idx], errMsg);
	}
	if(fP && strcmp(outFileStr, "-"))
	{
	  fclose(fP);
	  fP = NULL;
	}
	if(!ok)
	{
	  break;
	}
      }
    }
  }
  if(ok)
  {
    for(idx = 0; idx < 4; ++idx)
    {
      if(inObj[idx])
      {
        if((inObj[idx]->type != WLZ_LUT) ||
	   (inObj[idx]->domain.core == NULL) ||
	   (inObj[idx]->values.core == NULL) ||
	   (inObj[idx]->values.lut->vType != WLZ_GREY_INT))
        {
	  ok = 0;
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	                 "%s: invalid object type read from file %s (%s)\n",
			 *argv, inFileStr[idx], errMsg);
	  break;
	}
      }
    }
  }
  /* Merge the integer look up tables to form a RGBA look up table. */
  if(ok)
  {
    outObj = WlzAssignObject(
             WlzLUTMergeToRGBA(inObj[0], inObj[1], inObj[2], inObj[3],
    			       &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to merge integer look up tables (%s)\n",
		     *argv, errMsg);
    }
  }
  /* Output the merged look up table object. */
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL) ||
	(WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write LUT object to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
      fP = NULL;
    }
  }
  for(idx = 0; idx < 4; ++idx)
  {
    WlzFreeObj(inObj[idx]);
  }
  WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%s",
    *argv,
    " [-o<output object>] [-h]\n"
    "         [<red lut>[ [<green lut>[ [<blue lut>] [<alpha lut>]\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "Reads integer look up table objects given on the command line and\n"
    "merges them to form a red,green,blue,alpha look up table object.\n"
    "Any of the look up table object file names can be replaced by the\n"
    "string \"NULL\". Any of the input files may be read from the\n"
    "standard input. If no files are given on the command line all files\n"
    "will be read from the standard input, otherwise a \"-\" in place of\n"
    "a file name may be used to specify standard input. When reading from\n"
    "the standard input an empty object may be used as well as integer\n"
    "look up table objects, with these read in the order red, green, blue,\n"
    "alpha (the same as the command line order). When a NULL or empty\n"
    "object is used it implies an identity look up table.\n"
    "By default the merged look up table is written to the standard output.\n"
    "Example:\n"
    "  cat green.wlz | ",
    *argv,
    "-o rgba.wlz red.wlz - NULL alpha.wlz\n"
    "Merges three look up tables red.wlz, green.wlz and alpha.wlz together\n"
    "with and identity look up table for the blue channel to form a RGBA\n"
    "look up table object that is written to rgba.wlz.\n");
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
