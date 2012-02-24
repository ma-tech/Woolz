#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRMSDiff_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzRMSDiff.c
* \author       Bill Hill
* \date         October 2003
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
* \brief	Computes the RMS difference between the grey values of objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzrmsdiff "WlzRMSDiff"
*/

/*!
\ingroup BinWlz
\defgroup wlzrmsdiff WlzRMSDiff
\par Name
WlzRMSDiff - computes the RMS difference between the grey values of objects.
\par Synopsis
\verbatim
WlzRMSDiff [-o<out>] [<in obj 0>] [<in obj 1>]
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
Computes the RMS difference in grey values within the intersection
of the two given domain objects.
The input objects are read from stdin and values are written to stdout
unless the filenames are given.
\par Examples
\verbatim
WlzRMSDiff -o out.num  in0.wlz in1.wlz
\endverbatim
An RMS difference in grey values between in0.wlz and in1.wlz is written
to the file out.num.
\par File
\ref WlzRMSDiff.c "WlzRMSDiff.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		idx,
  		vol,
		option,
		ok = 1,
		usage = 0;
  double	rms;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*difObj = NULL;
  WlzObject	*inObj[2];
  FILE		*fP = NULL;
  char 		*outFileStr;
  char  	*inObjFileStr[2];
  const char	*errMsg;
  static char	optList[] = "o:h",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObj[0] = NULL;
  inObj[1] = NULL;
  outFileStr = outFileStrDef;
  inObjFileStr[0] = inObjFileStrDef;
  inObjFileStr[1] = inObjFileStrDef;
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
	ok = 0;
	break;
    }
  }
  if((inObjFileStr[0] == NULL) || (*inObjFileStr[0] == '\0') ||
     (inObjFileStr[1] == NULL) || (*inObjFileStr[1] == '\0') ||
     (outFileStr == NULL) || (*outFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    idx = 0;
    while((idx < 2) && (optind < argc))
    {
      inObjFileStr[idx] = *(argv + optind);
      ++optind;
      ++idx;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  if(ok)
  {
    /* Read objects. */
    idx = 0;
    while((errNum == WLZ_ERR_NONE) && (idx < 2))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((inObjFileStr[idx] == NULL) ||
	  (*inObjFileStr[idx] == '\0') ||
	  ((fP = (strcmp(inObjFileStr[idx], "-")?
		  fopen(inObjFileStr[idx], "r"): stdin)) == NULL) ||
	  ((inObj[idx] = WlzAssignObject(WlzReadObj(fP,
	  					    &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to read object %d from file %s (%s)\n",
		       *argv, idx, inObjFileStr[idx], errMsg);
      }
      if(fP && strcmp(inObjFileStr[idx], "-"))
      {
	fclose(fP);
      }
      ++idx;
    }
  }
  if(ok)
  {
    /* Check object types. */
    if(inObj[0]->type != inObj[1]->type)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((inObj[0]->domain.core == NULL) ||
            (inObj[1]->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      switch(inObj[0]->type)
      {
        case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
        case WLZ_3D_DOMAINOBJ:
	  break;
        default:
          errNum = WLZ_ERR_OBJECT_TYPE;
          ok = 0;
          break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: input object(s) not appropriate (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    difObj = WlzImageArithmetic(inObj[0], inObj[1], WLZ_BO_SUBTRACT,
    				0, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      vol = WlzGreyStats(difObj, NULL, NULL, NULL, NULL, &rms, NULL, NULL,
                         &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rms = (rms > DBL_EPSILON)? sqrt(rms / vol): 0.0;
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute RMS (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(outFileStr, "-")? fopen(outFileStr, "w"):
	      			       stdout)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to open output file (%s).\n",
		     *argv, errMsg);
    }
    else
    {
      (void )fprintf(fP, "%g\n", rms);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  (void )WlzFreeObj(difObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out>] [<in obj 0>] [<in obj 1>]\n"
    "Options:\n"
    "  -o  Output file name.\n"
    "  -h  Help, prints this usage message.\n"
    "Computes the RMS difference in grey values within the intersection\n"
    "of the two given domain objects.\n"
    "The input objects are read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    *argv,
    " -o out.num  in0.wlz in1.wlz\n"
    "An RMS difference in grey values between in0.wlz and in1.wlz is written\n"
    "to the file out.num.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
