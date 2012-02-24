#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTransposePlanes_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzTransposePlanes.c
* \author       Richard Baldock
* \date         July 2000
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
* \brief	Transposes the planes of 3D objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlztransposeplanes "WlzTransposePlanes"
*/

/*!
\ingroup      BinWlz
\defgroup     wlztransposeplanes WlzTransposePlanes
\par Name
WlzTransposePlanes - Transposes the planes of 3D Woolz objects.
\par Synopsis
\verbatim
WlzTransposePlanes [-o<offset>] [-h] [-v] [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-o</b></td>
    <td>offset for transpose. </td>
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
By  default  the  input  object is read from the standard input and the
output object is written to the standard output.

\par Description
Transpose the planes of a 3D woolz object. The  transpose  relation  is
new-plane  = offset - old-plane. This is equivalent to reflection about
the X-Y plane through offset/2. By default the offset is the sum of the
first and last planes of the input object.

\par Examples
\verbatim
# The following command transposes the planes of the object leaving
# the first and last planes unchanged.

WlzTransposePlanes  infile.wlz > outfile.wlz


# The following command transposes the planes of the object relecting
# the object about the origin

WlzTransposePlanes  -o0 infile.wlz > outfile.wlz
\endverbatim

\par File
\ref WlzTransposePlanes.c "WlzTransposePlanes.c"
\par See Also
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
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

static WlzObject *WlzTransposePlanes(WlzObject	*obj,
				     int	offset,
				     WlzErrorNum	*dstErr);

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-o#] [-h] [-v] [<input file>]\n"
	  "\tTranspose planes in the 3D object with respect\n"
	  "\tto the offset so that p' = offset - p\n"
	  "\tThis is equivalent to reflection about offset/2\n"
	  "\tOptions are:\n"
	  "\t  -o#       offset(default plane1+lastpl)\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  int		option;
  char 		optList[] = "ho:v";
  int		offset;
  int		offsetFlg;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		verboseFlg=0;
    
  /* set the defaults */
  offsetFlg = 0;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'o':
      offset = atoi(optarg);
      offsetFlg = 1;
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_PARAM_TYPE;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return WLZ_ERR_PARAM_TYPE;
    }
  }

  /* read objects and select planes if possible */
  while((obj = WlzReadObj(inFile, &errNum)) != NULL) 
  {
    if( errNum == WLZ_ERR_NONE ){
      switch( obj->type )
      {
      case WLZ_3D_DOMAINOBJ:
	if( obj->domain.p->type == WLZ_EMPTY_OBJ ){
	  errNum = WlzWriteObj(stdout, obj);
	  break;
	}
	/* check the transpose parameters */
	if( !offsetFlg ){
	  offset = obj->domain.p->plane1 + obj->domain.p->lastpl;
	}

	if((nobj = WlzTransposePlanes(obj, offset, &errNum)) != NULL){
	  errNum = WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
	break;

      default:
	errNum = WlzWriteObj(stdout, obj);
	break;
      }
    }

    WlzFreeObj(obj);
  }

  return errNum;
}

WlzObject *WlzTransposePlanes(
  WlzObject	*obj,
  int		offset,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  WlzDomain	domain;
  WlzValues	values;
  int		indx, indx1, p;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	 errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( obj->domain.p->type ){
	case WLZ_EMPTY_OBJ:
	case WLZ_EMPTY_DOMAIN:
	  return WlzMakeEmpty(dstErr);

	default:
	  domain.p = WlzMakePlaneDomain(obj->domain.p->type,
					offset-obj->domain.p->lastpl,
					offset-obj->domain.p->plane1,
					obj->domain.p->line1,
					obj->domain.p->lastln,
					obj->domain.p->kol1,
					obj->domain.p->lastkl,
					&errNum);
	  if( errNum == WLZ_ERR_NONE ){
	    domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
	    domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
	    domain.p->voxel_size[2] = obj->domain.p->voxel_size[2];
	    if( obj->values.core ){
	      values.vox =
		WlzMakeVoxelValueTb(obj->values.vox->type,
				    offset-obj->domain.p->lastpl,
				    offset-obj->domain.p->plane1,
				    obj->values.vox->bckgrnd,
				    NULL, &errNum);
	    }
	    else {
	      values.core = NULL;
	    }
	  }
	  break;
	}
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* make links */
  if( errNum == WLZ_ERR_NONE ){
    indx = obj->domain.p->lastpl - obj->domain.p->plane1;
    indx1 = 0;
    for(p=domain.p->plane1; p <= domain.p->lastpl; p++, indx--, indx1++){
      if( obj->domain.p->domains[indx].core ){
	domain.p->domains[indx1] =
	  WlzAssignDomain(obj->domain.p->domains[indx], &errNum);
      }
      else {
	domain.p->domains[indx1].core = NULL;
      }

      if( values.core ){
	if( obj->values.vox->values[indx].core ){
	  values.vox->values[indx1] =
	    WlzAssignValues(obj->values.vox->values[indx], &errNum);
	}
	else {
	  values.vox->values[indx1].core = NULL;
	}
      }
    }
  }

  /* create object */
  if( errNum == WLZ_ERR_NONE ){
    rtnObj = WlzMakeMain(obj->type, domain, values,
			 NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
