#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzMeshTransformObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzMeshTransformObj.c
* \author       Richard Baldock
* \date         March 2005
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
* \brief	Applies a mesh transform to an object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzmeshtransformobj "WlzMeshTransformObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzmeshtransformobj WlzMeshTransformObj
\par Name
WlzMeshTransformObj - applies a mesh transform to an object.
\par Synopsis
\verbatim
WlzMeshTransformObj -m <mesh transform file> [-o <output file>] [-h] [-v]
                    <2D object input file>
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
    <td><b>-L</b></td>
    <td>Use linear interpolation, default is nearest-neighbour.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Mesh transform object.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output filename, default is stdout.</td>
  </tr>
</table>
\par Description
Applies a mesh transform to an object.
\par Examples
\verbatim
WlzMeshTransformObj -m mesh.wlz -o out.wlz obj.wlz
\endverbatim
Uses the mesh transform read from mesh.wlz to transform the 2D domain object
with read from obj.wlz. The resulting warped object is written to out.wlz.
\par File
\ref WlzMeshTransformObj.c "WlzMeshTransformObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzbasisfntransformobj "WlzBasisFnTransformObj(1)"
\ref WlzMeshTransformObj "WlzMeshTransformObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
	  "Usage:\t%s -m <mesh transform file>"
	  " [-o <output file>] [-h] [-v]"
	  " [<2D object input file>]\n"
	  "\tApply a mesh transform to given input objects\n"
	  "\twriting the warped objects to standard output.\n"
	  "\tOptions are:\n"
	  "\t  -L                 Use linear interpolation instead of nearest-neighbour\n"
	  "\t  -m<meshfile>       Mesh transform object\n"
	  "\t  -o<output file>    Output filename, default to stdout\n"
	  "\t  -h                 Help - this message\n"
	  "\t  -v                 verbose operation\n"
	  "",
	  proc_str);
  return;
}

int main(int	argc,
	 char	**argv)
{
  WlzObject	*inObj, *meshObj, *outObj;
  FILE		*inFP, *outFP, *meshFP;
  char		*outFile, *meshFile;
  char 		optList[] = "Lm:o:hv";
  int		option;
  int		verboseFlg = 0;
  const char    *errMsg;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* additional defaults */
  outFile = "-";
  meshFile = NULL;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'L':
      interp = WLZ_INTERPOLATION_LINEAR;
      break;

    case 'm':
      meshFile = optarg;
      break;

    case 'o':
      outFile = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 0;

    case 'v':
      verboseFlg = 1;
      break;

    }
  }

  /* check input file/stream */
  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* check output file/stream */
  if(strcmp(outFile, "-"))
  {
    if((outFP = fopen(outFile, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    outFP = stdout;
  }

  /* check meshfile */
  if((errNum == WLZ_ERR_NONE) && (meshFile != NULL)){
    if( meshFP = fopen(meshFile, "r") ){
      if( meshObj = WlzReadObj(meshFP, &errNum) ){
	if( meshObj->type != WLZ_MESH_TRANS ){
	  fprintf(stderr,
		  "%s: mesh file does not have a mesh transform object\n",
		  argv[0]);
	  return 1;
	}
      }
      else {
	fprintf(stderr, "%s: failed to read mesh file\n",
	    argv[0]);
	return 1;
      }
    }
    else {
      fprintf(stderr, "%s: failed to open mesh file\n",
	    argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: mesh input file required\n",
	    argv[0]);
    usage(argv[0]);
    return 1;
  }

  /* transform any suitable input objects */
  while((errNum == WLZ_ERR_NONE) &&
        ((inObj = WlzReadObj(inFP, &errNum)) != NULL))
  {
    switch( inObj->type )
    {
    case WLZ_2D_DOMAINOBJ:
      if( outObj = WlzMeshTransformObj(inObj, meshObj->domain.mt,
				       interp, &errNum) ){
	WlzWriteObj(outFP, outObj);
	WlzFreeObj(outObj);
      }
      else {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to transform object (%s).\n",
		       *argv, errMsg);
	return 1;
      }
      break;

    default:
      WlzWriteObj(outFP, inObj);
      break;
    }

    WlzFreeObj(inObj);
  }

  if(errNum == WLZ_ERR_READ_EOF)
  {
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
