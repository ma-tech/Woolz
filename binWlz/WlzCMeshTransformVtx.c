#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshTransformVtx_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzCMeshTransformVtx.c
* \author       Bill Hill
* \date         December 2009
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2009 Medical research Council, UK.
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
* \brief	Transforms vertices using a constrained mesh transform.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \ref 		wlzcmeshtransformvtx "WlzCMeshTransformVtx"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcopyobj WlzCMeshTransformVtx
\par Name
WlzCMeshTransformVtx - transforms vertices using a constrained mesh transform.
\par Synopsis
\verbatim
WlzCMeshTransformVtx [-h] [-i] [-o<output file>] [-x<input file>] [<transform>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Invert the transform after reading.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-x</b></td>
    <td>Input file.</td>
  </tr>
</table>
\par Description
Reads a constrained mesh transform and input vertices objects and outputs
transformed vertices.
\par Example
\verbatim
WlzCMeshTransformVtx -i -o new.num -x in.num tx.wlz
\endverbatim
Reads a constrained mesh transform from the file tx.wlz and inverts it as
instructed by the -i flag. Vertices are then read from the file in.num,
transformed and written to the file out.num.

The input vertices must be in the format:
\verbatim
<x> <y> <z>
\endverbatim
That is three double precission numbers per line seperated by whitespace.
If the transform is in 2D then the last number (which must be present)
can have any value as it will not be used.
The output vertices will be in the format:
\verbatim
<e> <x> <y> <z>
\endverbatim
Here the first number is an integer corresponding to the Woolz error code
for the transformation of the vertex (value is zero for no error).
The remaining three numbers are the double precission transformed vertex.
For 2D the z coordinate will still be present, but will have value zero.
\par File
\ref WlzCMeshTransformVtx.c "WlzCMeshTransformVtx.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		idx,
		dim = 0,
  		inv = 0,
  		ok = 1,
  		option,
  		usage = 0;
  FILE		*inFP = NULL,
  		*outFP = NULL;
  char		*txFileStr,
		*inFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzObject	*tx = NULL;
  WlzDVertex2	v2;
  WlzDVertex3	v3;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "iho:x:";
  const char    txFileStrDef[] = "-",
  		inFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  opterr = 0;
  txFileStr = (char *)txFileStrDef;
  inFileStr = (char *)inFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'i':
        inv = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'x':
        inFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      txFileStr = argv[optind];
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if(((inFP = (strcmp(txFileStr, "-")?
	        fopen(txFileStr, "r"): stdin)) == NULL) ||
       ((tx = WlzAssignObject(WlzReadObj(inFP, &errNum), NULL)) == NULL))
    {
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(tx->type)
      {
        case WLZ_CMESH_2D:
	  dim = 2;
	  break;
	case WLZ_CMESH_3D:
	  dim = 3;
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
	  "%s: Failed to read constrained mesh transform object from\n"
	  "file %s (%s)\n",
	  *argv, txFileStr, errMsgStr);
    }
    if((inFP != NULL) && strcmp(txFileStr, "-"))
    {
      (void )fclose(inFP); inFP = NULL;
    }
  }
  if(ok)
  {
    if((inFP = (strcmp(inFileStr, "-")?
               fopen(inFileStr, "r"): stdin)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open input vertex file %s.\n",
		     argv[0], inFileStr);
    }
  }
  if(ok)
  {
    if((outFP = (strcmp(outFileStr, "-")?
                fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output vertex file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok && inv)
  {
    if((errNum = WlzCMeshTransformInvert(tx)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to invert transform (%s).\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    idx = 0;
    do
    {
      if(fscanf(inFP, "%lg %lg %lg", &(v3.vtX), &(v3.vtY), &(v3.vtZ)) != 3)
      {
	if(feof(inFP))
	{
	  errNum = WLZ_ERR_EOO;
	}
	else
	{
          errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else
      {
        switch(dim)
	{
	  case 2:
	    v2.vtX = v3.vtX;
	    v2.vtY = v3.vtY;
	    errNum = WlzCMeshTransformVtxAry2D(tx, 1, &v2);
	    v3.vtX = v2.vtX;
	    v3.vtY = v2.vtY;
	    v3.vtZ = 0.0;
	    break;
	  case 3:
	    errNum = WlzCMeshTransformVtxAry3D(tx, 1, &v3);
	    break;
	  default:
	    break;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if(fprintf(outFP, "%d %lg %lg %lg\n",
	           (int )errNum, v3.vtX, v3.vtY, v3.vtZ) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      ++idx;
    } while(errNum == WLZ_ERR_NONE);
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to transform %d'th vertex (%s).\n",
		     argv[0],
		     idx,
		     errMsgStr);
    }
  }
  (void )WlzFreeObj(tx);
  if(inFP && strcmp(inFileStr, "-"))
  {
    (void )fclose(inFP);
  }
  if(outFP && strcmp(outFileStr, "-"))
  {
    (void )fclose(outFP);
  }
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-i] [-o<output vtx file>] [-x<input vtx file>]\n"
	    "                   [<constrained mesh transform object>]\n"
	    "Reads a constrained mesh transform object and input vertices\n"
	    "and outputs transformed vertices.\n"
	    "Options are:\n"
            "  -h  Help, prints this usage message.\n"
            "  -i  Invert the transform after reading.\n"
	    "  -o  Output vertex file.\n"
	    "  -x  Input vertex file.\n"
	    "The input vertices must be in the format:\n"
	    "  <x> <y> <z>\n"
	    "That is three double precission numbers per line seperated by\n"
	    "whitespace.  If the transform is in 2D then the last number\n"
	    "(which must be present) can have any value as it will not be\n"
	    "used.  The output vertices will be in the format:\n"
	    "  <e> <x> <y> <z>\n"
	    "Here the first number is an integer corresponding to the Woolz\n"
	    "error code for the transformation of the vertex (value is zero\n"
	    "for no error). The remaining three numbers are the double\n"
	    "precission transformed vertex.  For 2D the z coordinate will\n"
	    "still be present, but will have value zero.\n"
	    "By default all files are read from and written to the standard\n"
	    "input and standard output. The transform object is read before\n"
	    "any vertices are read.\n"
	    "Example:\n"
	    "  %s -i -o new.num -x in.num tx.wlz\n"
	    "Reads a constrained mesh transform from the file tx.wlz and\n"
	    "inverts it as instructed by the -i flag. Vertices are then read\n"
	    "from the file in.num, transformed and written to the file\n"
	    "out.num.\n",
	    argv[0], argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */






