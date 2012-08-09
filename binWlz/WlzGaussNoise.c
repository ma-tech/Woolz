#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGaussNoise_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGaussNoise.c
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
* \brief	Adds Gaussian noise to the grey values of an object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzgaussnoise "WlzGaussNoise"
*/

/*!
\ingroup BinWlz
\defgroup wlzgaussnoise WlzGaussNoise
\par Name
WlzGaussNoise - add Gaussian noise to the grey values of an object.
\par Synopsis
\verbatim
WlzGaussNoise [-s#] [-h] [-v] [<input mask> [<input obj>]]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Standard deviation of the noise with a default of 5.</td>
  </tr>
</table>
\par Description
Add gaussian noise to the grey value of an object.
\par Examples
\verbatim
WlzGaussNoise -s 10 in.wlz >noisy.wlz
\endverbatim
Adds Gaussian noise to the object read from in.wlz and writes the noisy
object to noisy.wlz.
\par File
\ref WlzGaussNoise.c "WlzGaussNoise.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzGaussNoise "WlzGaussNoise(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
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
	  "Usage:\t%s [-s#] [-h] [-v] [<input mask> [<input obj>]]\n"
	  "\tAdd gaussian noise to the grey value of an object.\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -s#       standard deviation of the noise - default 5\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n",
	  proc_str,
	  WlzVersion());
  return;
}


int main(int	argc,
	 char	**argv)
{
  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "hs:v";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzPixelV	val;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  val.type = WLZ_GREY_FLOAT;
  val.v.flv = 5.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 's':
      val.v.flv = atof(optarg);
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* check for read from file */
  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects and convert if possible */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    if( obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ ){
      if( (errNum = WlzGaussNoise(obj, val)) == WLZ_ERR_NONE ){
	(void )WlzWriteObj(stdout, obj);
      }
      else {
	return errNum;
      }
    }
    else {
      (void )WlzWriteObj(stdout, obj);
    }

    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
