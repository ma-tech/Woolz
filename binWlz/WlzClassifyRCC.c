#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzClassifyRCC_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzClassifyRCC.c
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
* \brief	Classifies an ordered pair of domain objects using a
* 		region connected calculus.
* \ingroup	BinWlz
* \par Binary
* \ref wlzclassifyrcc   "WlzClassifyRCC"
*/

/*!
\ingroup BinWlz
\defgroup wlzimagearithmetic WlzImageArithmetic
\par Name
WlzImageArithmetic - binary image arithmetic on a pair of domain objects.
\par Synopsis
\verbatim
WlzImageArithmetic [-o<out file>] [-a] [-d] [-l] [-g] [-i] [-m] [-s]
                   [-n] [-N] [-h] [-O#] [<in object 0>] [<in object 1>]
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
  <tr> 
    <td><b>-a</b></td>
    <td>Add the object's grey values.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Divide the grey values of the 1st object by those of the 2nd.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Vector magnitude of horizontal and vertical component objects.</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>Dither values when setting range.</td>
  </tr>
  <tr>
    <td><b>-l</b></td>
    <td>Compute the modulus of the grey values of the 1st object wrt
        those of the 2nd.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Multiply the object's grey values.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Subtract the grey values of the 2nd object from those of the 1st.</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Normalises the output object to the range of the input objects.</td>
  </tr>
  <tr>
    <td><b>-N</b></td>
    <td>Normalises the output objects to the range [0-255].</td>
  </tr>
  <tr>
    <td><b>-O</b></td>
    <td>Overwrite option (only useful for debugging).</td>
  </tr>
</table>
\par Description
Computes an arithmetic binary (two objects) operation on two domain
objects. The default operator is add.
The input objects are read from stdin and values are written to stdout
unless the filenames are given.
\par Examples
\verbatim
cat obj1.wlz | WlzImageArithmetic -o obj3.wlz -a  obj2.wlz
\endverbatim
A new object obj3.wlz is formed by adding the grey values of obj1 and
obj2.wlz.
\par File
\ref WlzImageArithmetic.c "WlzImageArithmetic.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzImageArithmetic "WlzImageArithmetic(3)"
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
		noEnc = 0,
		ok = 1,
		usage = 0;
  double	nrmVol = 0.0;
  WlzRCCClass 	cls = WLZ_RCC_EMPTY;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr = NULL;
  WlzObject	*iObj[2] = {NULL};
  char  	*iObjFileStr[2] = {NULL};
  const char	*errMsg;
  static char	optList[] = "o:eh",
		fileStrDef[] = "-";

  opterr = 0;
  outFileStr = fileStrDef;
  iObjFileStr[0] = fileStrDef;
  iObjFileStr[1] = fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'e':
        noEnc = 1;
	break;
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
  if((iObjFileStr[0] == NULL) || (*iObjFileStr[0] == '\0') ||
     (iObjFileStr[1] == NULL) || (*iObjFileStr[1] == '\0') ||
     (outFileStr      == NULL) || (*outFileStr      == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    int		i = 0;

    while((i < 2) && (optind < argc))
    {
      iObjFileStr[i] = *(argv + optind);
      ++optind;
      ++i;
    }
  }
  if(ok && (optind != argc))
  {
    usage = 1;
    ok = 0;
  }
  if(ok)
  {
    int		i = 0;

    while((errNum == WLZ_ERR_NONE) && (i <= 1))
    {
      errNum = WLZ_ERR_READ_EOF;
      if((iObjFileStr[i] == NULL) ||
	  (*iObjFileStr[i] == '\0') ||
	  ((fP = (strcmp(iObjFileStr[i], "-")?
		  fopen(iObjFileStr[i], "r"): stdin)) == NULL) ||
	  ((iObj[i] = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to read object %d from file %s (%s)\n",
		       *argv, i, iObjFileStr[i], errMsg);
      }
      if(fP && strcmp(iObjFileStr[i], "-"))
      {
	(void )fclose(fP);
      }
      ++i;
    }
  }
  if(ok)
  {
    cls = WlzRegConCalcRCC(iObj[0], iObj[1], noEnc, &nrmVol, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute classification (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?  fopen(outFileStr, "w"):
	      				 stdout)) != NULL))
    {
      int	first = 1;
      unsigned int m;

      errNum = WLZ_ERR_NONE;
      for(m = 1; m <= WLZ_RCC_MSK; m <<= 1)
      {
        if(m & cls)
	{
	  const char *str;

	  str = WlzStringFromRCC((WlzRCCClass )m, NULL);
	  if(str == NULL)
	  {
	    errNum = WLZ_ERR_PARAM_TYPE;
	  }
	  else
	  {
	    if(fprintf(fP, "%s%s", (first)? "": "|", str) < 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	    first = 0;
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
      if(first)
      {
        if(fprintf(fP, "WLZ_RCC_EMPTY") < 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if(fprintf(fP, "\t%g\n", nrmVol) < 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  WlzFreeObj(iObj[0]);
  WlzFreeObj(iObj[1]);
  if(usage)
  {

    fprintf(stderr,
	    "Usage: %s"
	    "  [-o<out file>] [-e] [-h] [<in object 0>] [<in object 1>]\n"
	    "Version: %s\n"
	    "Options:\n"
	    "  -o        Output file name.\n"
	    "  -e        Don't include enclosure classifications.\n"
	    "  -h        Help, prints this usage message.\n"
	    "Classifies the given domain objects using the region connected\n"
	    "calculus defined in Woolz:\n"
	    "  WLZ_RCC_EMPTY - One or both of the domains is EMPTY.\n"
	    "  WLZ_RCC_DC    - DisConnected.\n"
	    "  WLZ_RCC_EC    - Edge Connected.\n"
	    "  WLZ_RCC_EQ    - EQual.\n"
	    "  WLZ_RCC_PO    - Partial Overlap.\n"
	    "  WLZ_RCC_TPP   - Tangential Proper Part.\n"
	    "  WLZ_RCC_NTPP  - Non-Tangential Proper Part.\n"
	    "  WLZ_RCC_TPPI  - Non-Tangential Proper Part Inverse.\n"
	    "  WLZ_RCC_SUR   - SURrounds.\n"
	    "  WLZ_RCC_SURI  - SURrounds Inverse.\n"
	    "  WLZ_RCC_ENC   - ENCloses.\n"
	    "  WLZ_RCC_ENCI  - ENCloses Inverse.\n"
	    "See WlzRegConCalcRCC(3) for more information.\n"
	    "Example:\n%s obj0.wlz obj1.wlz\n"
	    "Classifies the object read from obj0.wlz with respect to the\n"
	    "object read from obj1.wlz\n",
	    *argv,
	    WlzVersion(),
	    *argv);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
