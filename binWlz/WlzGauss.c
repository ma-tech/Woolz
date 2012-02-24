#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGauss_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGauss.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Applies a Gaussian filter to an objects grey values.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzgauss "WlzGauss"
*/

/*!
\ingroup BinWlz
\defgroup wlzgauss WlzGauss
\par Name
WlzGauss - applies a Gaussian filter to an objects grey values.
\par Synopsis
\verbatim
WlzGauss [-w #[ #]] [-x#] [-y#] [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-w</b></td>
    <td>
    Gaussian widths in the x and (optionaly) y directions,
    specified as full width half maximum in pixels,
    with a default value of 3.0. If a single width parameter is given
    the y width is set equal to the x width.
    </td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>
    Order of the x derivative with possible values 0,1,2 and
    a default value of 0.
    </td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>
    Order of the x derivative with possible values 0,1,2 and
    a default value of 0.
    </td>
  </tr>
</table>
\par Description
Applies a Gaussian filter to the grey values of a Woolz object.
\par Examples
\verbatim
WlzGauss -w 5 in.wlz >smooth.wlz
\endverbatim
Smooths the grey values of the object read from in.wlz using a Gaussian filter
with a full width half maximum of 5 pixels and both derivatives equal to zero.
The resulting object is written to smooth.wlz.
\par File
\ref WlzGauss.c "WlzGauss.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzrankobj "WlzRankObj(1)"
\ref wlzrsvfilterobj "WlzRsvFilterObj(1)"
\ref WlzGauss2 "WlzGauss2(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-w#[#]] [-x#] [-y#] [-h] [<input file>]\n"
	  "\tApply a Gaussian filter to a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -w#[,#]   x_width[y_width] gaussian widths defined as\n"
	  "\t            full width half maximum in pixels\n"
	  "\t            default value 3.0, if y_width omitted\n"
	  "\t            then it is set equal to x_width\n"
	  "\t  -x#       x derivative - possible values 0,1,2, default - 0\n"
	  "\t  -y#       y derivative - possible values 0,1,2, default - 0\n"
	  "\t  -h        Help - prints this usage message\n",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "hw:x:y:";
  int		option;
  double	x_width, y_width;
  int		x_deriv, y_deriv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* set defaults, read the argument list and check for an input file */
  opterr = 0;
  x_width = 3.0;
  y_width = 3.0;
  x_deriv = 0;
  y_deriv = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'w':
      switch( sscanf(optarg, "%lg,%lg", &x_width, &y_width) ){

      default:
      case 0:
	fprintf(stderr, "%s: no gaussian width set\n", argv[0]);
        usage(argv[0]);
        return 1;

      case 1:
	y_width = x_width;
	break;

      case 2:
	break;

      }
      break;

    case 'x':
      x_deriv = atoi(optarg);
      break;

    case 'y':
      y_deriv = atoi(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects and threshold if possible */
  while((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) 
  {
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
	if( (nobj = WlzGauss2(obj, x_width, y_width,
			      x_deriv, y_deriv, &errNum)) != NULL ){
	  errNum = WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
	break;

      case WLZ_3D_DOMAINOBJ:
      default:
	errNum = WlzWriteObj(stdout, obj);
	break;
    }

    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
