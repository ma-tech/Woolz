#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshTransformVtx_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshTransformVtx.c
* \author       Bill Hill
* \date         December 2009
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
* \brief	Transforms vertices using a constrained mesh transform.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshtransformvtx "WlzCMeshTransformVtx"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshtransformvtx WlzCMeshTransformVtx
\par Name
WlzCMeshTransformVtx - transforms vertices using a constrained mesh transform.
\par Synopsis
\verbatim
WlzCMeshTransformVtx [-h] [-i] [-o<output file>] [-f<input file>] [<transform>]
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
    <td><b>-f</b></td>
    <td>Input file.</td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>Interpolation value, with 0 <= ivalue <= 1.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Number of intermediate interpolations.</td>
  </tr>
  <tr>
    <td><b>-b</b></td>
    <td>Output object file body</td>
  </tr>
  <tr>
    <td><b>-e</b></td>
    <td>Output object file extension (default wlz)</td>
  </tr>
</table>
\par Description
Reads a constrained mesh transform and input vertices objects and outputs
transformed vertices.
\par Example
\verbatim
WlzCMeshTransformVtx -i -o new.num -f in.num tx.wlz
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
Partly transformed vertices can be generated if interpolation
values between 0 and 1 are given(0 means no transformation)
Alternatively a set of interpolations covering the full range
of transformation is computed if the number of interpolations
is given. For this, the output base filename and its extension
must be supplied.
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

static WlzErrorNum 		WlzCMeshTransPrintEncAndVtx(
				  char *outFileStr,
				  int dim,
				  int *elm,
				  WlzVertexP v,
				  int nV);
static WlzErrorNum 		WlzCMeshTransVtxEncAndTrans(
				  WlzObject *trObj,
				  int dim,
				  int *elm,
				  WlzVertexP v,
				  int nV);

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		idx,
		nV = 0,
		dim = 0,
  		inv = 0,
		maxV = 0,
		nStep = 1,
		useStep = 0,
  		ok = 1,
  		option,
  		usage = 0;
  double	transition = 1.0;
  int		*elm = NULL;
  FILE		*fP = NULL;
  char		*txFileStr,
		*inFileStr,
		*outFileBodyStr = NULL,
		*outFileExtStr = NULL;
  char		outFileStr[FILENAME_MAX];
  const char	*errMsgStr;
  WlzVertexP	v;
  WlzObject	*trObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "b:e:f:iho:s:x:";
  const char    txFileStrDef[] = "-",
  		inFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  v.v = NULL;
  opterr = 0;
  txFileStr = (char *)txFileStrDef;
  inFileStr = (char *)inFileStrDef;
  (void )strcpy(outFileStr, outFileStrDef);
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'b':
        outFileBodyStr = optarg;
        break;
      case 'e':
        outFileExtStr = optarg;
        break;
      case 'i':
        inv = 1;
	break;
      case 'o':
        (void )strncpy(outFileStr, optarg, FILENAME_MAX - 1);
	outFileStr[FILENAME_MAX - 1] = '\0';
	break;
      case 'f':
        inFileStr = optarg;
	break;
      case 's':
        useStep = 1;
        if(sscanf(optarg, "%d", &nStep) != 1)
        {
          usage = 1;
        } 
        break;
      case 'x':
        useStep = 0;
        if(sscanf(optarg, "%lg", &transition) != 1)
        {
          usage = 1;
        }
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
  /* Read the transform. */
  if(ok)
  {
    if(((fP = (strcmp(txFileStr, "-")?
	        fopen(txFileStr, "r"): stdin)) == NULL) ||
       ((trObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      if(errNum != WLZ_ERR_NONE)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(trObj->type)
      {
        case WLZ_CMESH_2D:
	  dim = 2;
	  break;
	case WLZ_CMESH_2D5: /* FALLTHROUGH */
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
    if((fP != NULL) && strcmp(txFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Invert the transform if required. */
  if(ok && inv && trObj->values.core)
  {
    WlzObject *newTr = NULL;

    newTr = WlzCMeshTransformInvert(trObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzFreeObj(trObj);
      trObj = WlzAssignObject(newTr, NULL);
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to invert transform (%s).\n",
		     argv[0],
		     errMsgStr);
    }
  }
  /* Read the vertices. */
  if(ok)
  {
    if((fP = (strcmp(inFileStr, "-")?
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
    idx = 0;
    if(dim == 2)
    {
      do
      {
	if(idx >= maxV)
	{
	  maxV = (maxV <= 0)? 1024: 2 * maxV;
	  if((v.v = AlcRealloc(v.v, sizeof(WlzDVertex2) * maxV)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  WlzDVertex2 *u;

	  u = v.d2 + idx;
	  if(fscanf(fP, "%lg %lg", &(u->vtX), &(u->vtY)) != 2)
	  {
	    if(feof(fP))
	    {
	      errNum = WLZ_ERR_EOO;
	    }
	    else
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	  }
	}
        ++idx;
      } while(errNum == WLZ_ERR_NONE);
    }
    else /* dim == 3 */
    {
      do
      {
	if(idx >= maxV)
	{
	  maxV = (maxV <= 0)? 1024: 2 * maxV;
	  if((v.v = AlcRealloc(v.v, sizeof(WlzDVertex3) * maxV)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  WlzDVertex3 *u;

	  u = v.d3 + idx;
	  if(fscanf(fP, "%lg %lg %lg", &(u->vtX), &(u->vtY), &(u->vtZ)) != 3)
	  {
	    if(feof(fP))
	    {
	      errNum = WLZ_ERR_EOO;
	    }
	    else
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	  }
	}
        ++idx;
      } while(errNum == WLZ_ERR_NONE);
      nV = idx;
    }
    if(errNum == WLZ_ERR_EOO)
    {
      --nV;
      errNum = WLZ_ERR_NONE;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if((elm = (int *)AlcMalloc(sizeof(int) * nV)) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to read %d'th vertex (%s).\n",
		     argv[0],
		     idx,
		     errMsgStr);
    }
  }
  if(fP && strcmp(inFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(ok)
  {
    if(useStep)
    {
      int	step;
      WlzVertexP w;

      w.v = (dim == 2)? AlcMalloc(sizeof(WlzDVertex2) * maxV)
                      : AlcMalloc(sizeof(WlzDVertex3) * maxV);
      if(w.v == NULL)
      {
	  errNum = WLZ_ERR_MEM_ALLOC;
      }
      step = 0;
      while((errNum == WLZ_ERR_NONE) && (step < nStep))
      {
	WlzObject *stepTrObj = NULL;

        (void )sprintf(outFileStr, "%s%08d.%s",
	               outFileBodyStr, step, outFileExtStr);
        stepTrObj = WlzAssignObject(
	            WlzCopyScaleCMeshValue(
		      (double )step / (double )(nStep - 1),
		      trObj, &errNum), NULL);
        if(errNum == WLZ_ERR_NONE)
	{
	  if(dim == 2)
	  {
	    WlzValueCopyDVertexToDVertex(w.d2, v.d2, nV);
	  }
	  else
	  {
	    WlzValueCopyDVertexToDVertex3(w.d3, v.d3, nV);
	  }
	  errNum = WlzCMeshTransVtxEncAndTrans(stepTrObj, dim, elm, w, nV);
	}
	if(errNum == WLZ_ERR_NONE)
	{
          errNum =  WlzCMeshTransPrintEncAndVtx(outFileStr, dim, elm, w, nV);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    ok = 0;
	    (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	    (void )fprintf(stderr,
			   "%s: Failed to write vertices to file %s (%s).\n",
			   argv[0], outFileStr, errMsgStr);
	  }
	}
	WlzFreeObj(stepTrObj);
      }
      AlcFree(w.v);
    }
    else
    {
      if(fabs(transition - 1.0) > WLZ_MESH_TOLERANCE)
      {
        errNum = WlzScaleCMeshValue(transition, trObj);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzCMeshTransVtxEncAndTrans(trObj, dim, elm, v, nV);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	errNum =  WlzCMeshTransPrintEncAndVtx(outFileStr, dim, elm, v, nV);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
          (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	  (void )fprintf(stderr,
			 "%s: Failed to write vertices to file %s (%s).\n",
			 argv[0], outFileStr, errMsgStr);
	}
      }
    }
  }
  AlcFree(v.v);
  AlcFree(elm);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-i] [-o<output vtx file>]\n"
	    "                   [-b <body>] [-e<extension>]\n"
	    "                   [-f<input vtx file>] [-x #] [-s #]\n"
	    "                   [<constrained mesh transform object>]\n"
	    "Reads a constrained mesh transform object and input vertices\n"
	    "and outputs transformed vertices.\n"
	    "Version: %s\n"
	    "Options are:\n"
            "  -h  Help, prints this usage message.\n"
            "  -i  Invert the transform after reading.\n"
	    "  -o  Output vertex file.\n"
	    "  -b  Output object file body\n"
	    "  -e  Output object file extension (default wlz)\n"
	    "  -f  Input vertex file.\n"
	    "  -s  Number of intermediate interpolations.\n"
	    "  -x  Interpolation value, with 0 <= ivalue <= 1.\n"
	    "The input vertices must be in the format:\n"
	    "  <x> <y>\n"
	    "for 2D transforms and\n"
	    "  <x> <y> <z>\n"
	    "for 3D transforms.\n"
	    "The output vertices will be in the format:\n"
	    "  <e> <x> <y>\n"
	    "for 2D transforms and\n"
	    "  <e> <x> <y> <z>\n"
	    "For 3D transforms. In both output cases the first number is an\n"
	    "integer corresponding to the index of the mesh element\n"
	    "enclosing the vertex.\n"
	    "Partly transformed vertices can be generated if interpolation\n"
	    "values between 0 and 1 are given(0 means no transformation)\n"
	    "Alternatively a set of interpolations covering the full range\n"
	    "of transformation is computed if the number of interpolations\n"
	    "is given. For this, the output base filename and its extension\n"
	    "must be supplied.\n"
	    "By default all files are read from and written to the standard\n"
	    "input and standard output. The transform object is read before\n"
	    "any vertices are read.\n"
	    "Example:\n"
	    "  %s -i -o new.num -x in.num tx.wlz\n"
	    "Reads a constrained mesh transform from the file tx.wlz and\n"
	    "inverts it as instructed by the -i flag. Vertices are then read\n"
	    "from the file in.num, transformed and written to the file\n"
	    "out.num.\n",
	    argv[0],
	    WlzVersion(),
	    argv[0]);

  }
  return(!ok);
}

/*!
* \return	Woolz error code.
* \ingroup	BinWlz
* \brief	Finds the element enclosing each vertex and transforms
* 		the vertex using the given mesh transform.
* \param	trObj			Conforming mesh.
* \param	elm			Array for element indices, on
* 					return set to the those enclosing the
* 					vertices.
* \param	v			Given vertices to be transformed.
* \param	nV			Number of vertices (and element
*					indices).
*/
static WlzErrorNum WlzCMeshTransVtxEncAndTrans(WlzObject *trObj, int dim,
				int *elm, WlzVertexP v, int nV)
{
  int		idx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(dim == 2)
  {
    for(idx = 0; idx < nV; ++idx)
    {
      WlzDVertex2 *u;

      u = v.d2 + idx;             
      elm[idx] = WlzCMeshElmEnclosingPos2D(trObj->domain.cm2, -1,
					   u->vtX, u->vtY, 0, NULL);
    }
    if(trObj->values.core)
    {
      errNum = WlzCMeshTransformVtxAry2D(trObj, nV, v.d2);
    }
  }
  else /* dim == 3 */
  {
    for(idx = 0; idx < nV; ++idx)
    {
      WlzDVertex3 *u;

      u = v.d3 + idx;             
      if(trObj->domain.core->type == WLZ_CMESH_2D5)
      {
	elm[idx] = WlzCMeshElmEnclosingPos2D5(trObj->domain.cm2d5, -1,
					      u->vtX, u->vtY, u->vtZ, 0, NULL);
      }
      else /* WLZ_CMESH_3D */
      {
	elm[idx] = WlzCMeshElmEnclosingPos3D(trObj->domain.cm3, -1,
					     u->vtX, u->vtY, u->vtZ, 0, NULL);
      }
    }
    if(trObj->values.core)
    {
      if(trObj->domain.core->type == WLZ_CMESH_2D5)
      {
        errNum = WlzCMeshTransformVtxAry2D5(trObj, nV, v.d3);
      }
      else /* WLZ_CMESH_3D */
      {
        errNum = WlzCMeshTransformVtxAry3D(trObj, nV, v.d3);
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	BinWlz
* \brief	Outputs the enclosing element index and vertex coordinates
* 		to the named file.
* \param	outFileStr		Output file name, "-" implies the
* 					standard output.
* \param	dim			Dimension of the vertices.
* \param	elm			Element indices.
* \param	v			Vertices.
* \param	nV			Number of element indices and
* 					vertices.
*/
static WlzErrorNum WlzCMeshTransPrintEncAndVtx(char *outFileStr, int dim,
				int *elm, WlzVertexP v, int nV)
{
  int		idx;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"): stdout)) == NULL)
  {
    errNum = WLZ_ERR_FILE_OPEN;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dim == 2)
    {
      WlzDVertex2 *u;

      for(idx = 0; idx < nV; ++idx)
      {
	u = v.d2 + idx;
	if(fprintf(fP, "%d %lg %lg\n",
		   elm[idx], u->vtX, u->vtY) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
    else /* dim == 3 */
    {
      WlzDVertex3 *u;

      for(idx = 0; idx < nV; ++idx)
      {
	u = v.d3 + idx;
	if(fprintf(fP, "%d %lg %lg %lg\n",
		   elm[idx], u->vtX, u->vtY, u->vtZ) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP);
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
