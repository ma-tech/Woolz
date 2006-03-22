#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRsvFilterObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzRsvFilterObj.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Recursive filter for domain objects with grey values.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzrsvfilterobj "WlzRsvFilterObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzrsvfilterobj WlzRsvFilterObj
\par Name
WlzRsvFilterObj - recursive filter for domain objects with grey values.
\par Synopsis
\verbatim
WlzRsvFilterObj - [-h] [-o<output object>]
                  [-a#,#,#,#] [-b#,#] [-c#]
		  [-P] [-n] [-g] [-d] [-p#] [-r#] [-t#]
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
    <td><b>-a</b></td>
    <td>Filter feed forward parameters.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Filter feed back parameters.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Filter normalization parameter.</td>
  </tr>
  <tr> 
    <td><b>-P</b></td>
    <td>Print filter parameters to the standard output.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Dont filter object.</td>
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
    <td>Filter parameter for Gaussian (sigma) and Deriche (alpha).</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Order of derivative for Gaussian and Deriche filters [0-2].</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Filter directions, eg x filter along lines and x,y filter
        along lines then through columns.</td>
  </tr>
</table>
\par Description
Applies a recursive filter to a Woolz domain object with grey values.
The input object is read from stdin and the filtered object is
written to stdout unless the filenames are given.
Objects with unsigned byte grey values will be promoted to short
grey valued objects.
\par Examples
\verbatim
WlzRsvFilterObj -g -p 2.0 -r 0 -t y -o smoothed.wlz in.wlz
\endverbatim
The input Woolz object is read from in.wlz, and filtered using a
Gaussian filter with varience = 2.0. The is only applied through
the objects columns. The filtered object is then written to
smoothed.wlz.
\par File
\ref WlzRsvFilterObj.c "WlzRsvFilterObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgauss "WlzGauss(1)"
\ref WlzRsvFilterObj "WlzRsvFilterObj(3)"
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
		order = 0,
		useA = 0,
		useB = 0,
		useC = 0,
		useFtr = 1,
		printFtrPrm = 0;
  unsigned int	actMsk = WLZ_RSVFILTER_ACTION_X | WLZ_RSVFILTER_ACTION_Y |
  			 WLZ_RSVFILTER_ACTION_Z;
  double	ftrPrm = 1.0;
  double 	c;
  double	a[4],
  		b[2];
  WlzRsvFilterName ftrName = WLZ_RSVFILTER_NAME_GAUSS_0;
  WlzRsvFilter *ftr = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inObjFileStr,
		*actStr;
  static char	optList[] = "ho:a:b:c:Pngdp:r:t:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  c = 0.0;
  a[0] = a[1] = a[2] = a[3] = 0.0;;
  b[0] = b[1] = 0.0;

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
      case 'a':
	useA = 1;
	if((optarg == NULL) ||
	   (sscanf(optarg, "%lg,%lg,%lg,%lg",
	           &(a[0]), &(a[1]), &(a[2]), &(a[3])) != 4))
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'b':
	useB = 1;
	if((optarg == NULL) ||
	   (sscanf(optarg, "%lg,%lg", &(b[0]), &(b[1])) != 2))
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'c':
        useC = 1;
	if((optarg == NULL) || (sscanf(optarg, "%lg", &c) != 1))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'g':
        ftrName = WLZ_RSVFILTER_NAME_GAUSS_0;
	break;
      case 'd':
        ftrName = WLZ_RSVFILTER_NAME_DERICHE_0;
	break;
      case 'P':
        printFtrPrm = 1;
	break;
      case 'n':
        useFtr = 0;
        break;
      case 'p':
	if(sscanf(optarg, "%lg", &ftrPrm) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'r':
	if(sscanf(optarg, "%d", &order) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 't':
	actMsk = WLZ_RSVFILTER_ACTION_NONE;
	if(optarg)
	{
	  actStr = strtok(optarg, ",");
	  while(actStr && ok)
	  {
	    if(strcmp(actStr, "x") == 0)
	    {
	      actMsk |= WLZ_RSVFILTER_ACTION_X;
	    }
	    else if(strcmp(actStr, "y") == 0)
	    {
	      actMsk |= WLZ_RSVFILTER_ACTION_Y;
	    }
	    else if(strcmp(actStr, "z") == 0)
	    {
	      actMsk |= WLZ_RSVFILTER_ACTION_Z;
	    }
	    else
	    {
	      usage = 1;
	      ok = 0;
	    }
	    actStr = strtok(NULL, ",");
	  }
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
      switch(order)
      {
        case 0:
	  ftrName = (ftrName == WLZ_RSVFILTER_NAME_DERICHE_0)?
	  	    WLZ_RSVFILTER_NAME_DERICHE_0:
		    WLZ_RSVFILTER_NAME_GAUSS_0;
          break;
        case 1:
	  ftrName = (ftrName == WLZ_RSVFILTER_NAME_DERICHE_0)?
	  	    WLZ_RSVFILTER_NAME_DERICHE_1:
		    WLZ_RSVFILTER_NAME_GAUSS_1;
          break;
        case 2:
	  ftrName = (ftrName == WLZ_RSVFILTER_NAME_DERICHE_0)?
	  	    WLZ_RSVFILTER_NAME_DERICHE_2:
		    WLZ_RSVFILTER_NAME_GAUSS_2;
          break;
        default:
	  ok = 0;
	  (void )fprintf(stderr,
	  		 "%s: derivative order must be in range [0-2]\n",
			 *argv);
	  break;
      }
      if(((ftr = WlzRsvFilterMakeFilter(ftrName, ftrPrm, &errNum)) == NULL) ||
         (errNum != WLZ_ERR_NONE))
      {
        ok = 0;
	(void )fprintf(stderr, "%s Failed to create filter (%s)\n",
		       *argv, WlzStringFromErrorNum(errNum, NULL));
      }
      else
      {
	if(useA || useB || useC)
	{
	  ftr->name = WLZ_RSVFILTER_NAME_NONE;
	  if(useA)
	  {
	    ftr->a[0] = a[0];
	    ftr->a[1] = a[1];
	    ftr->a[2] = a[2];
	    ftr->a[3] = a[3];
	  }
	  if(useB)
	  {
	    ftr->b[0] = b[0];
	    ftr->b[1] = b[1];
	  }
	  if(useC)
	  {
	    ftr->c = c;
	  }
	}
      }
    }
    if(ok)
    {
      if(printFtrPrm)
      {
        switch(ftr->name)
	{
	  case WLZ_RSVFILTER_NAME_NONE:
	    (void )printf("name WLZ_RSVFILTER_NAME_NONE\n");
	    break;
	  case WLZ_RSVFILTER_NAME_DERICHE_0:
	    (void )printf("name WLZ_RSVFILTER_NAME_DERICHE_0\n");
	    break;
	  case WLZ_RSVFILTER_NAME_DERICHE_1:
	    (void )printf("name WLZ_RSVFILTER_NAME_DERICHE_1\n");
	    break;
	  case WLZ_RSVFILTER_NAME_DERICHE_2:
	    (void )printf("name WLZ_RSVFILTER_NAME_DERICHE_2\n");
	    break;
	  case WLZ_RSVFILTER_NAME_GAUSS_0:
	    (void )printf("name WLZ_RSVFILTER_NAME_GAUSS_0\n");
	    break;
	  case WLZ_RSVFILTER_NAME_GAUSS_1:
	    (void )printf("name WLZ_RSVFILTER_NAME_GAUSS_1\n");
	    break;
	  case WLZ_RSVFILTER_NAME_GAUSS_2:
	    (void )printf("name WLZ_RSVFILTER_NAME_GAUSS_2\n");
	    break;
	  default:
	    (void )fprintf(stderr, "%s Unknown filter name(%d)\n",
	    		   *argv, ftr->name);
	    ok = 0;
	    break;
	}
	if(ok)
	{
	  (void )printf("a %g %g %g %g\n",
	  		ftr->a[0], ftr->a[1], ftr->a[2], ftr->a[3]);
	  (void )printf("b %g %g\n",
	  		ftr->b[0], ftr->b[1]);
	  (void )printf("c %g\n",
	  		ftr->c);
	}
      }
    }
    if(ok && useFtr)
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
        if(((outObj = WlzAssignObject(WlzRsvFilterObj(inObj, ftr,
						  actMsk,
						  &errNum), NULL)) == NULL) ||
	   (errNum != WLZ_ERR_NONE))

        {
	  ok = 0;
	  (void )fprintf(stderr, "%s Failed to filter object (%s)\n",
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
    " [-h] [-o<output object>] [-a#,#,#,#] [-b#,#] [-c#]\n"
    "        [-P] [-n] [-g] [-d] [-p#] [-r#] [-t#]\n"
    "        [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information\n"
    "  -o  Output object file name.\n"
    "  -a  Filter feed forward parameters.\n"
    "  -b  Filter feed back parameters.\n"
    "  -c  Filter normalization parameter.\n"
    "  -P  Print filter parameters to the standard output.\n"
    "  -n  Dont filter object.\n"
    "  -g  Use approximate Gaussian filter.\n"
    "  -d  Use Deriche filter.\n"
    "  -p  Filter parameter for Gaussian (sigma) and Deriche (alpha) filters\n"
    "  -r  Order of derivative for Gaussian and Deriche filters [0-2].\n"
    "  -t  Filter directions, eg x filter along lines and x,y filter\n"
    "      along lines then through columns.\n"
    "Applies a recursive filter to a Woolz domain object with grey values.\n"
    "The input object is read from stdin and the filtered object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -g -p 2.0 -r 0 -t y -o smoothed.wlz in.wlz\n"
    "The input Woolz object is read from in.wlz, and filtered using a\n"
    "Gaussian filter with varience = 2.0. The is only applied through\n"
    "the objects columns. The filtered object is then written to\n"
    "smoothed.wlz.\n"
    "Objects with unsigned byte grey values will be promoted to short\n"
    "grey valued objects.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
