#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzEuclidDistanceMap_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzEuclidDistanceMap.c
* \author       Richard Baldock
* \date         January 2006
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
* \ingroup      BinWlzApp
* \brief        Calculate the Euclidean distance from the input vertex to
*               all parts of the domain.
*               
* \par Binary
* \ref wlzeucliddistancemap "WlzEuclidDistanceMap"
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzeucliddistancemap WlzEuclidDistanceMap
\par Name
WlzEuclidDistanceMap - calculate the Euclidean distance from a vertex.
\par Synopsis
\verbatim
WlzEuclidDistanceMap -p x,y,z  -h -v 

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-p x,y,z</b></td>
    <td>Input vertex, default (0,0,0).</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>

\par Description
Calculate the Euclidean distance from an input vertex to all pixel locations
 in a given domain. The input object values are converted to float
 (WLZ_GREY_FLOAT) type and overwritten with the distance value. Currently
 the object must have a grey-table and of type WLZ_2D_DOMAINOBJ.

\par Examples
\verbatim
\endverbatim

\par See Also
\par Bugs
None known
\todo
Add a grey-table if required and extend to 3D.
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(
  char	*str)
{
  fprintf(stderr,
	  "Usage:\n"
	  "%s -p x,y,z  -h -v \n"
	  "\tRead in domain object from stdin and calculate the Euclidean\n"
	  "\tdistance from the input vertex to each point in the domain,\n"
	  "\twriting to stdout.\n"
	  "Arguments:\n"
	  "\t-p#,#,#    input vertex, default (0,0,0)\n"
	  "\t-h         print this message\n"
	  "\t-v         verbose operation\n"
	  "\n",
	  str);

  return;
}

int main(
  int   argc,
  char  **argv)
{
  FILE		*inFile=stdin;
  char 		optList[] = "p:hv";
  int		option;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		verboseFlg=0;
  int		type;
  WlzObject	*obj, *obj1;
  int		i;
  double	x, y, z, val;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyP		gptr;

  /* read the argument list and check for an input file */
  opterr = 0;
  x = y = z = 0.0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'p':
      if( sscanf(optarg, "%lf,%lf,%lf", &x, &y, &z) != 3 ){
	fprintf(stderr, "%s: not enough values for vertex\n", argv[0]);
	usage(argv[0]);
	return 1;
      }
      break;

    case 't':
      type = atoi(optarg);
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

  /* get objects from stdin */
  if((obj = WlzReadObj(inFile, &errNum)) != NULL){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      break;

    case WLZ_3D_DOMAINOBJ: /* to be done */
    default:
      fprintf(stderr, "%s: invalid object type: %d\n", argv[0],
	      obj->type);
      usage(argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: can't read object\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* scan through object setting distance values
     use integer values */
  if((obj1 = WlzConvertPix(obj, WLZ_GREY_FLOAT, &errNum)) != NULL){
    WlzInitGreyScan(obj1, &iwsp, &gwsp);
    while( WlzNextGreyInterval(&iwsp) == 0 ){

      gptr = gwsp.u_grintptr;
      switch (gwsp.pixeltype) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp.colrmn; i++, gptr.inp++){
	  val = (iwsp.colpos + i - x)*(iwsp.colpos + i - x) \
	    + (iwsp.linpos - y) * (iwsp.linpos - y);
	  val = sqrt(val);
	  *gptr.flp = WLZ_NINT(val);
	}
	break;

      default:
	/* something gone badly wrong - quit */
	fprintf(stderr, "%s: something wrong, oops\n", argv[0]);
	return 1;
      }
    }
    
    WlzWriteObj(stdout, obj1);
    WlzFreeObj(obj1);
  }

  WlzFreeObj(obj);
  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
