#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzAffineTransformDiff_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzAffineTransformDiff.c
* \author       Bill Hill
* \date         July 2004
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
* \brief	Compares two given affine transforms.
* \ingroup 	BinWlz
*
* \par Binary:
* \ref wlzaffinetransformdiff "WlzAffineTransformDiff"
*/

/*!
\ingroup BinWlz
\defgroup wlzaffinetransformdiff WlzAffineTransformDiff
\par Name
WlzAffineTransformDiff - compares two given affine transforms.

\par Synopsis
\verbatim
WlzAffineTransformDiff [-h] [-e #] [-o<output file>] <transform 0> <transform 1>
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-e</b></td>
    <td>Tolerance value.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file for differences.</td>
  </tr>
</table>
\par Description

Compares two given affine transforms for significant
differences using the given tollerance value.

The exit status is 0 if the two transforms do not significantly
differ, otherwise it is 1.
\par Examples
\verbatim
WlzAffineTransformDiff -e 0.00001 -o diff.out tr0.wlz tr1.wlz
\endverbatim
Reads a pair of affine transforms from tr0.wlz and tr1.wlz, then
writes any differences to diff.out.
\par File
\ref WlzAffineTransformDiff.c "WlzAffineTransformDiff.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzAffineTransformObj "WlzAffineTransformObj(3)"
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
  int		idM,
  		idN,
  		option,
		ok = 1,
		usage = 0;
  double	dif,
  		tol = 1e-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		dim[3];
  WlzObject	*obj[2];
  FILE		*fP = NULL;
  char 		*objFileStr[2];
  char		*outFileStr;
  static char	optList[] = "he:o:";

  opterr = 0;
  obj[0] = obj[1] = NULL;
  objFileStr[0] = objFileStr[1] = "-";
  outFileStr = "-";
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'e':
        if(sscanf(optarg, "%lg", &tol) != 1)
	{
	  usage = 1;
	}
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default: /* FALLTHROUGH */
        usage = 1;
	break;
    }
  }
  if(ok && (optind < argc))
  {
    if((optind + 2) != argc)
    {
      usage = 1;
    }
    else
    {
      objFileStr[0] = *(argv + optind);
      objFileStr[1] = *(argv + optind + 1);
    }
  }
  ok = !usage;
  /* Read the affine transforms objects. */
  if(ok)
  {
    idN = 0;
    while(ok && (idN < 2))
    {
      if((objFileStr[idN] == NULL) ||
	 (*(objFileStr[idN]) == '\0') ||
	 ((fP = (strcmp(objFileStr[idN], "-")?
		fopen(objFileStr[idN], "r"): stdin)) == NULL) ||
	 ((obj[idN] = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: Failed to read object from file %s.\n",
		       *argv, objFileStr[idN]);
      }
      if(fP && strcmp(objFileStr[idN], "-"))
      {
	(void )fclose(fP); fP = NULL;
      }
      ++idN;
    }
  }
  /* Check object types. */
  if(ok)
  {
    idN = 0;
    while(ok && (idN < 2))
    {
      if(obj[idN]->type != WLZ_AFFINE_TRANS)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: Object read from %s is not an affine transform.\n",
		       *argv, objFileStr[idN]);
      }
      else if(obj[idN]->domain.core == NULL)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: Object read from %s has a NULL domain.\n",
		       *argv, objFileStr[idN]);
      }
      else
      {
        switch(obj[idN]->domain.t->type)
	{
          case WLZ_TRANSFORM_2D_AFFINE:  /* FALLTHROUGH */
  	  case WLZ_TRANSFORM_2D_REG:     /* FALLTHROUGH */
  	  case WLZ_TRANSFORM_2D_TRANS:   /* FALLTHROUGH */
  	  case WLZ_TRANSFORM_2D_NOSHEAR:
	    dim[idN] = 2;
	    break;
  	  case WLZ_TRANSFORM_3D_AFFINE:  /* FALLTHROUGH */
  	  case WLZ_TRANSFORM_3D_REG:     /* FALLTHROUGH */
  	  case WLZ_TRANSFORM_3D_TRANS:   /* FALLTHROUGH */
          case WLZ_TRANSFORM_3D_NOSHEAR:
	    dim[idN] = 3;
	    break;
	  default:
	    ok = 0;
	    (void )fprintf(stderr,
		       "%s: Object read from %s is not an affine transform.\n",
			   *argv, objFileStr[idN]);
	    break;
	}
      }
      ++idN;
    }
  }
  /* Open the output file. */
  if(ok)
  {
    if((outFileStr == NULL) ||
       (*outFileStr == '\0') ||
       ((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     *argv, outFileStr);
    }
  }
  /* Compare the transform dimensions. */
  if(ok)
  {
    if(dim[0] != dim[1])
    {
      ok = 0;
      (void )fprintf(fP, "dimension\n< %d\n---\n> %d\n", dim[0], dim[1]);
    }
    dim[2] = WLZ_MAX(dim[0], dim[1]);
    for(idM = 0; idM < dim[2]; ++idM)
    {
      for(idN = 0; idN < dim[2]; ++idN)
      {
        dif = obj[0]->domain.t->mat[idM][idN] -
	      obj[1]->domain.t->mat[idM][idN];
        dif = fabs(dif);
	if(dif > tol)
	{
	  ok = 0;
	  (void )fprintf(fP, "[%d][%d]\n< %1.8e\n---\n> %1.8e\n",
			 idM, idN,
	  		 obj[0]->domain.t->mat[idM][idN],
			 obj[1]->domain.t->mat[idM][idN]);
	}
      }
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(obj[0]);
  (void )WlzFreeObj(obj[1]);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-e #] [-o<output file>] <transform 0> <transform 1>\n" 
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -e  Tolerance value.\n"
    "  -o  Output file for differences.\n"
    "Compares two given affine transforms for significant differences\n"
    "using the given tollerance value.\n"
    "The exit status is 0 if the two transforms do not significantly\n"
    "differ, otherwise it is 1.\n"
    "The output file (if given) will be touched even if there are no\n"
    "differences found.\n",
    *argv,
    " -e 0.00001 -o diff.out tr0.wlz tr1.wlz\n"
    "Reads a pair of affine transforms from tr0.wlz and tr1.wlz, then\n"
    "writes any differences to diff.out.\n");
  }
  return((ok)? 0: 1);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
