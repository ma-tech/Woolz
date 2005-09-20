#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzGreyGradient.c
* \author       Bill Hill
* \date         July 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Computes a grey gradient object from the input object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzgreygradient "WlzGreyGradient"
*/

/*!
\ingroup BinWlz
\defgroup wlzgreygradient WlzGreyGradient
\par Name
WlzGreyGradient - computes a grey gradient object from the input object.
\par Synopsis
\verbatim
WlzGreyGradient [-o<output object>] [-h] [-o] [-g] [-d] [-m#] [-p#]
                [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Use approximate Gaussian filter.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Use Deriche filter.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Filter parameter for Gaussian (sigma) and Deriche (alpha) filters.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Extra multiplication factor for gradient grey levels.</td>
  </tr>
</table>
\par Description
Computes a grey gradient object from the input object.
The input object is read from stdin and the gradient object is
written to stdout unless the filenames are given.
\par Examples
\verbatim
WlzGreyGradient -d -m 2.0 -o grad.wlz in.wlz
\endverbatim
The input Woolz object is read from in.wlz, and the grey gradient
image is computed using a Deriche filter and an additional multiplying
factor of 2.0. The gradient object is written to grad.wlz.
Objects with unsigned byte grey values will be promoted to short
grey valued objects.
\par File
\ref WlzGreyGradient.c "WlzGreyGradient.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzcannyderiche "WlzCannyDeriche(1)"
\ref wlzgreycrossing "WlzGreyCrossing(1)"
\ref WlzGreyGradient "WlzGreyGradient(3)"
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
		ok = 1,
		usage = 0,
		printFtrPrm = 0;
  unsigned int	actMsk = WLZ_RSVFILTER_ACTION_X | WLZ_RSVFILTER_ACTION_Y |
  			 WLZ_RSVFILTER_ACTION_Z;
  double	ftrMlt = 1.0,
  		ftrPrm = 1.0;
  WlzRsvFilterName ftrName = WLZ_RSVFILTER_NAME_DERICHE_1;
  WlzRsvFilter *ftr = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inObjFileStr;
  static char	optList[] = "ho:gdm:p:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'g':
        ftrName = WLZ_RSVFILTER_NAME_GAUSS_1;
	break;
      case 'd':
        ftrName = WLZ_RSVFILTER_NAME_DERICHE_1;
	break;
      case 'm':
	if(sscanf(optarg, "%lg", &ftrMlt) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'p':
	if(sscanf(optarg, "%lg", &ftrPrm) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
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
    if(ok)
    {
      if(((ftr = WlzRsvFilterMakeFilter(ftrName, ftrPrm, &errNum)) == NULL) ||
         (errNum != WLZ_ERR_NONE))
      {
        ok = 0;
	(void )fprintf(stderr, "%s Failed to create filter (%s)\n",
		       *argv, WlzStringFromErrorNum(errNum, NULL));
      }
      else
      {
	ftr->c *= ftrMlt;
      }
    }
    if(ok)
    {
      if((inObjFileStr == NULL) ||
	 (*inObjFileStr == '\0') ||
	 ((fP = (strcmp(inObjFileStr, "-")?
		fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	 ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
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
      if(ok)
      {
        if(((outObj = WlzAssignObject(
		      WlzGreyGradient(NULL, NULL, NULL,
		      		      inObj, ftr, &errNum), NULL)) == NULL) ||
	   (errNum != WLZ_ERR_NONE))

        {
	  ok = 0;
	  (void )fprintf(stderr, "%s Failed to compute gradient object (%s)\n",
	  		 *argv, WlzStringFromErrorNum(errNum, NULL));
	}
      }
      if(ok)
      {
	if(((fP = (strcmp(outFileStr, "-")?
		  fopen(outFileStr, "w"):
		  stdout)) == NULL) ||
	   (WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
	{
	  ok = 0;
	  (void )fprintf(stderr,
			 "%s: failed to write output object\n",
			 *argv);
	}
	if(fP && strcmp(outFileStr, "-"))
	{
	  fclose(fP);
	}
      }
    }
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(outObj)
  {
    WlzFreeObj(outObj);
  }
  if(ftr)
  {
    WlzRsvFilterFreeFilter(ftr);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<output object>] [-h] [-o] [-g] [-d] [-m#] [-p#]\n"
    "        [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information\n"
    "  -o  Output object file name.\n"
    "  -g  Use approximate Gaussian filter.\n"
    "  -d  Use Deriche filter.\n"
    "  -p  Filter parameter for Gaussian (sigma) and Deriche (alpha) filters\n"
    "  -m  Extra multiplication factor for gradient grey levels.\n"
    "Computes a grey gradient object from the input object.\n"
    "The input object is read from stdin and the gradient object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -d -m 2.0 -o grad.wlz in.wlz\n"
    "The input Woolz object is read from in.wlz, and the grey gradient\n"
    "image is computed using a Deriche filter and an additional multiplying\n"
    "factor of 2.0. The gradient object is written to grad.wlz.\n"
    "Objects with unsigned byte grey values will be promoted to short\n"
    "grey valued objects.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
