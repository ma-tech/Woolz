#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzErosion_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzErosion.c
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
* \brief	Morphological erosion of a woolz domain object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzerosion "WlzErosion"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzerosion WlzErosion
\par Name
WlzErosion - Morphological erosion of a woolz domain object
\par Synopsis
\verbatim
WlzErosion [-c<n>] [-r<radius>] [-h] [input file]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-c4</b></td>
    <td>4-connected erosion (2D)</td>
  </tr>
  <tr>
    <td><b>-c8</b></td>
    <td>8-connected erosion (2D)</td>
  </tr>
  <tr>
    <td><b>-c6</b></td>
    <td>6-connected erosion (3D)</td>
  </tr>
  <tr>
    <td><b>-c18</b></td>
    <td>18-connected erosion (3D)</td>
  </tr>
  <tr>
    <td><b>-c26</b></td>
    <td>26-connected erosion (3D)</td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Structuring element radius (default 1)</td>
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
By default the input object is read from the standard input
and the output is written to the standard output.

\par Description
WlzErosion erodes the input object using a structuring element defined
by the given connectivity. If a 2D connectivity is applied to a 3D
object then each plane is eroded independently.

\par Examples  
\verbatim
WlzErosion -r8 -c26 obj.wlz >eroded.wlz
\endverbatim
Erodes the 3D domain of object read from obj.wlz using a 26-connected
sphere of radius 8 and then writes the eroded object to eroded.wlz.
\par File
\ref WlzErosion.c "WlzErosion.c"
\par See Also
\ref wlzdilation "WlzDilation(1)"
\ref wlzdomainfill "WlzDomainFill(1)"
\ref wlzstructdilation WlzStructDilation(1)"
\ref wlzstructerosion WlzStructErosion(1)"
\ref WlzDilation "WlzDilation(3)"
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
	  "Usage:\t%s [-c#] [-h] [<input file>]\n"
	  "\tErode a domain woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -c#       Erosion connectivity:\n"
	  "\t            # =  4: 4-connected (2D)\n"
	  "\t                 8: 8-connected (2D) - default\n"
	  "\t                 6: 6-connected (3D)\n"
	  "\t                18: 18-connected (3D)\n"
	  "\t                26: 26-connected (3D)\n"
	  "\t  -r#       Structuring element radius (default 1)\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "hc:r:";
  int		option;
  WlzConnectType	connectivity;
  int		conn;
  int		radius = 1;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* read the argument list and check for an input file */
  opterr = 0;
  connectivity = WLZ_8_CONNECTED;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'c':
      switch( conn = atoi(optarg) ){

      case 8:
	connectivity = WLZ_8_CONNECTED;
	break;
	
      case 4:
	connectivity = WLZ_4_CONNECTED;
	break;
	
      case 6:
	connectivity = WLZ_6_CONNECTED;
	break;
	
      case 18:
	connectivity = WLZ_18_CONNECTED;
	break;
	
      case 26:
	connectivity = WLZ_26_CONNECTED;
	break;
	
      default:
        fprintf(stderr, "%s: connectivity = %d is invalid\n",
		argv[0], conn );
        usage(argv[0]);
        return WLZ_ERR_PARAM_DATA;

      }
      break;

    case 'r':
      radius = atoi(optarg);
      if( radius < 1 ){
	radius = 1;
	fprintf(stderr, "%s: radius reset to %d\n",
		argv[0], radius );
      }
      else if( radius > 100 ){
	radius = 100;
	fprintf(stderr, "%s: radius reset to %d\n",
		argv[0], radius );
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  /* read objects and threshold if possible */
  while((obj = WlzAssignObject(WlzReadObj(inFile, &errNum), NULL)) != NULL) 
  {
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
	if( radius > 1 ){
	  WlzObject	*structElem;

	  if((obj->type == WLZ_2D_DOMAINOBJ) ||
	     (connectivity == WLZ_8_CONNECTED) ||
	     (connectivity == WLZ_4_CONNECTED)){
	    
	    structElem = WlzMakeSphereObject(WLZ_2D_DOMAINOBJ, radius,
					     0.0, 0.0, 0.0, &errNum);
	  }
	  else {
	    structElem = WlzMakeSphereObject(WLZ_3D_DOMAINOBJ, radius,
					     0.0, 0.0, 0.0, &errNum);
	  }

	  if(structElem && (errNum == WLZ_ERR_NONE) &&
	     (nobj = WlzStructErosion(obj, structElem, &errNum)) ){
	    errNum = WlzWriteObj(stdout, nobj);
	    (void) WlzFreeObj(nobj);
	    (void) WlzFreeObj(structElem);
	  }

	}
	else {
	  if( (nobj = WlzErosion(obj, connectivity, &errNum)) != NULL ){
	    errNum = WlzWriteObj(stdout, nobj);
	    (void) WlzFreeObj(nobj);
	  }
	}
	break;

      default:
	errNum = WlzWriteObj(stdout, obj);
	break;
    }

    WlzFreeObj(obj);
  }

  /* trap the WLZ_ERR_READ_EOF since this is a legal way of indicating
     the end of objects in a file */
  if( errNum == WLZ_ERR_READ_EOF ){
    errNum = WLZ_ERR_NONE;
  }
  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
