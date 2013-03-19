#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzVtxToObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzVtxToObj.c
* \author	Richard Baldock
* \date         November 2000
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
* \brief	Builds a Woolz domain object from a vertex.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzvtxtoobj "WlzVtxToObj"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzvtxtoobj WlzVtxToObj
\par Name
WlzVtxToObj - builds a Woolz domain object from a vertex.
\par Synopsis
\verbatim
WlzVtxToObj  [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-X</b></td>
    <td>XXX.</td>
  </tr>
</table>
\par Description
WlzVtxToObj builds Woolz domain objects from vertices represented in text
format.
The vertices are read from the standard input and for each vertex a
corresponding domain object is written to the standard output.
\par Examples
\verbatim
echo '10 20 30' | WlzVtxToObj >out.wlz
\endverbatim
Writes a single 3D domain object to the file out.wlz
corresponding to the single voxel at (10,20,30).
\par File
\ref WlzVtxToObj.c "WlzVtxToObj.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
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
	  "Usage:\t%s [-d#] [-h] [-v]\n"
	  "\tBuild a Woolz domain object from a vertex reading the vertex\n"
	  "\tvalues from standard input\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t-d#     input dimension, can be 2 or 3 (default 2)\n"
	  "\t-h      help message\n"
	  "\t-v     verbose operation\n",
	  str,
	  WlzVersion());

  return;
}

WlzObject *WlzCoordsToObject(
  WlzObjectType	type,
  double	x,
  double	y,
  double	z,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL, *tmpObj;
  WlzDomain	domain;
  WlzValues	values;
  int		k, l, p;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check type */
  switch( type ){

  case WLZ_2D_DOMAINOBJ:
    k = x;
    l = y;
    if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					 l, l, k, k, &errNum)) != NULL){
      values.core = NULL;
      rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			   NULL, NULL, &errNum);
    }
    break;

  case WLZ_3D_DOMAINOBJ:
    if((tmpObj = WlzCoordsToObject(WLZ_2D_DOMAINOBJ,
				   x, y, z, &errNum)) != NULL){
      k = x;
      l = y;
      p = z;
      if((domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
					p, p, l, l, k, k, &errNum)) != NULL){
	values.core = NULL;
	domain.p->domains[0] = WlzAssignDomain(tmpObj->domain, NULL);
	rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			     NULL, NULL, &errNum);
      }
      WlzFreeObj(tmpObj);
    }
    break;

  default:
    errNum = WLZ_ERR_OBJECT_TYPE;
    break;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

int main(
  int   argc,
  char  **argv)
{
  
  WlzObject     	*obj;
  WlzDVertex3		vtx;
  char 	optList[] = "d:h";
  int		option;
  int		dimension=2;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* parse the command line */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'd':
      dimension = atoi(optarg);
      if( (dimension != 2) && (dimension != 3) ){
	fprintf(stderr, "%s: must have dimension equal to 2 or 3\n",
		argv[0]);
	usage(argv[0]);
	return 1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  /* read vertex from stdin and write the objects to stdout */
  while( scanf("%lg %lg %lg", &vtx.vtX, &vtx.vtY, &vtx.vtZ) != EOF ){
    if( dimension == 2 ){
      obj = WlzCoordsToObject(WLZ_2D_DOMAINOBJ, vtx.vtX, vtx.vtY, vtx.vtZ,
			      &errNum);
    }
    else {
      obj = WlzCoordsToObject(WLZ_3D_DOMAINOBJ, vtx.vtX, vtx.vtY, vtx.vtZ,
			      &errNum);
    }
    WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
