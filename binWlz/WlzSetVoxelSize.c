#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSetVoxelSize_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzSetVoxelSize.c
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
* \brief	Sets the voxel size of a 3D domain object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzsetvoxelsize "WlzSetVoxelSize"
*/

/*!
\ingroup BinWlz
\defgroup wlzsetvoxelsize WlzSetVoxelSize
\par Name
WlzSetVoxelSize - sets the voxel size of a 3D domain object.
\par Synopsis
\verbatim
WlzSetVoxelSize [-h] [-v] [-x#] [-y#] [-z#] [<input file>]
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
    <td><b>-x</b></td>
    <td>Voxel column (x) size.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Voxel line (y) size.</td>
  </tr>
  <tr> 
    <td><b>-z</b></td>
    <td>Voxel plane (z) thickness.</td>
  </tr>
</table>
\par Description
Sets the voxel size of the input 3D object writing the new object
to standard output.
If a voxel size is not defined then the original size is retained.
\par Examples
\verbatim
WlzSetVoxelSize -x 1 -y 1 -z 4 in.wlz >out.wlz
\endverbatim
Reads a 3D domain object from in.wlz,
sets the objects voxel size to be 1.0, 1.0 and 4.0 units in the
x, y and z dimensions.
The resulting object is then written to the out.wlz.
\par File
\ref WlzSetVoxelSize.c "WlzSetVoxelSize.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzfacts "WlzFacts(1)"
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
	  "Usage:\t%s [-x#] [-y#] [-z#] [-h] [-v] [<input file>]\n"
	  "\tReset the voxel sizes of the input 3D object\n"
	  "\twriting the new object to standard output.\n"
	  "\tThis is required until the 3D objects are converted to\n"
	  "\tWLZ_TRANS_OBJ type. If a voxel size is not defined then\n"
	  "\tthe original size is retained\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -x#       x voxel size\n"
	  "\t  -y#       y voxel size\n"
	  "\t  -z#       z voxel size\n"
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
  char 		optList[] = "x:y:z:hv";
  int		option;
  int		xFlg=0, yFlg=0, zFlg=0;
  float		x_size=1.0, y_size=1.0, z_size=1.0;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'x':
      x_size = atof(optarg);
      xFlg = 1;
      break;

    case 'y':
      y_size = atof(optarg);
      yFlg = 1;
      break;

    case 'z':
      z_size = atof(optarg);
      zFlg = 1;
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

  if( (x_size <= 0.0) || (y_size <= 0.0) || (z_size <= 0.0) ){
    fprintf(stderr, "%s: voxel sizes must be non-zero and positive\n",
	    argv[0]);
    return 1;
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
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    switch( obj->type )
    {
    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.p && (obj->domain.p->type == WLZ_2D_DOMAINOBJ) ){
	if( xFlg ){
	  obj->domain.p->voxel_size[0] = x_size;
	}
	if( yFlg ){
	  obj->domain.p->voxel_size[1] = y_size;
	}
	if( zFlg ){
	  obj->domain.p->voxel_size[2] = z_size;
	}
      }
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(1);
      }
      break;

    default:
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(1);
      }
      break;
    }
    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
