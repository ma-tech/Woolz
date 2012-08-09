#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzVolume_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzVolume.c
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
* \brief	Computes the volume of 3D objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzvolume "WlzVolume"
*/

/*!
\ingroup BinWlz
\defgroup wlzvolume WlzVolume
\par Name
WlzVolume - computes the volume of 3D objects.
\par Synopsis
\verbatim
WlzVolume [-h] [-n] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-n</b></td>
    <td>Numeric output only - number of voxels, zero on error</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Computes the volume of 3D objects.
\par Examples
\verbatim
WlzThreshold -v 200 -L rec.wlz | WlzLabel | WlzVolume
\endverbatim
Reads a 3D object from rec.wlz,
thresholds the objectt to remove background,
labels the thresholded object to compute seperate objects for the
disconnected regions
and then outputs the volume of each of these regions as ascii text.
\par File
\ref WlzVolume.c "WlzVolume.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzarea "WlzArea(1)"
\ref wlzthreshold "WlzThreshold(1)"
\ref wlzlabel "WlzLabel(1)"
\ref WlzVolume "WlzVolume(3)"
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
	  "Usage:\t%s [-n] [-h] [<input file>]\n"
	  "\tCalculate the volume of the input 3D woolz objects\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -n        Numeric output only - number of voxels, zero on error\n"
	  "\t  -h        Help - prints this usage message\n",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "nh";
  int		option;
  int		count, vol;
  int		numericFlg=0;
  const char    *errMsg;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case '~': /* dummy to avoid compiler warning */
      break;

    case 'n':
      numericFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      if( numericFlg ){
	fprintf(stdout, "0\n");
      }
      else {
	fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
	usage(argv[0]);
      }
	return( 1 );
    }
  }

  count = 0;
  while( (obj = WlzReadObj(inFile, NULL)) != NULL )
  {
    count++;

    switch( obj->type )
    {
     case WLZ_3D_DOMAINOBJ:
       if(((vol = WlzVolume(obj , &errNum)) < 0) ||
          (errNum != WLZ_ERR_NONE)){
	 if( numericFlg ){
	   fprintf(stdout, "0\n");
	 }
	 else {
	   (void )WlzStringFromErrorNum(errNum, &errMsg);
	   fprintf(stderr,
		   "%s: Object %d: error in calculating the volume (%s).\n",
		   *argv, count, errMsg);
	 }
         return 1 ;
       }
       else {
	 if( numericFlg ){
	   fprintf(stdout, "%d\n", vol);
	 }
	 else {
	   fprintf(stderr, "Object %d: number of voxels = %d\n", count, vol);
	 }
       }
       break;

    case WLZ_EMPTY_OBJ:
      if( numericFlg ){
	fprintf(stdout, "0\n");
      }
      else {
	fprintf(stderr, "Object %d: number of voxels = %d\n", count, 0);
      }
      break;

     default:
       fprintf(stderr, "Object %d: not 3D object type\n", count);
       break;

    }

    WlzFreeObj( obj );
  }

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
