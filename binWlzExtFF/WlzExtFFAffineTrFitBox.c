#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFAffineTrFitBox_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzExtFFAffineTrFitBox.c
* \author       Bill Hill
* \date         July 2015
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2015],
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
* \brief	Computes the Woolz affine transform which makes the
* 		bounding box of the first object equal to that of the
* 		second.
* \ingroup	BinWlzExtFF
*
* \par Binary
* \ref wlzextffaffinetransformfitbox "WlzExtFFAffineTrFitBox"
*/

/*!
\ingroup BinWlzExtFF
\defgroup wlzextffaffinetransformfitbox WlzExtFFAffineTrFitBox
\par Name
WlzExtFFAffineTrFitBox - computes an affine transform which makes
			        bounding boxes equal.
\par Synopsis
\verbatim
WlzExtFFAffineTrFitBox [-h] [-T] [-o<output file>]
                              [-s<source format>]
			      [-t<target format>]
			      [-u<output format>]
			      <source object> <target object>
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
  <td><b>-h</b></md>
    <td>Prints usage information.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Source file format.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Target file format.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>Output file format.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Transform the source object and write it to the output
        file rather than the affine transform.</td>
  </tr>
</table>
Computes the Woolz affine transform which makes the bounding box of the
source object equal to that of the target object. Valid file formats are
the same as for WlzExtFFConvert.

\par Examples
\verbatim
WlzExtFFAffineTrFitBox -o out.wlz small.vtk big.stl
\endverbatim
Reads the source object from the file small.vtk and the target object
from big.stl then computes the affine transform which scales the source
objects bounding box to fit that of the target object. The output transform
is written to the file out.wlz.
\par File
WlzExtFFAffineTrFitBox.c "WlzExtFFAffineTrFitBox.c"
\par See Also
\ref BinWlzExtFF "WlzIntro(1)"
\ref WlzExtFFConvert "WlzExtFFConvert(1)"
\ref WlzAffineTransformLSq "WlzAffineTransformLSq"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Wlz.h>
#include <WlzExtFF.h>

/*!
* \return	Dimension of the given object:
* \ingroup	WlzFeatures
* \brief	Returns the dimension of the given object.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzObjectDimension(WlzObject *obj, WlzErrorNum *dstErr)
{
  int		dim = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(obj->type)
    {
      /* Objects which don't have a dimension. */
      case WLZ_CONVOLVE_FLOAT: /* FALLTHROUGH */
      case WLZ_CONVOLVE_INT:   /* FALLTHROUGH */
      case WLZ_EMPTY_OBJ:      /* FALLTHROUGH */
      case WLZ_HISTOGRAM:      /* FALLTHROUGH */
      case WLZ_LUT:            /* FALLTHROUGH */
      case WLZ_PROPERTY_OBJ:   /* FALLTHROUGH */
      case WLZ_TEXT:
        dim = 0;
	break;
      /* Objects which are clearly 2D. */
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_2D_POLYGON:   /* FALLTHROUGH */
      case WLZ_BOUNDLIST:    /* FALLTHROUGH */
      case WLZ_CMESH_2D:     /* FALLTHROUGH */
      case WLZ_RECTANGLE:
        dim = 2;
	break;
      /* Objects which are clearly 3D. */
      case WLZ_3D_DOMAINOBJ:  /* FALLTHROUGH */
      case WLZ_3D_POLYGON:    /* FALLTHROUGH */
      case WLZ_3D_WARP_TRANS: /* FALLTHROUGH */
      case WLZ_CMESH_2D5:     /* FALLTHROUGH */
      case WLZ_CMESH_3D:
        dim = 3;
	break;
      /* Objects which need to have their domain checked. */
      case WLZ_TRANS_OBJ:
	if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  dim = WlzAffineTransformDimension(obj->domain.t, &errNum);
	}
        break;
      case WLZ_CONV_HULL:
	if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(obj->domain.core->type)
	  {
	    case WLZ_CONVHULL_DOMAIN_2D:
	      dim = 2;
	      break;
	    case WLZ_CONVHULL_DOMAIN_3D:
	      dim = 3;
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	}
        break;
      case WLZ_CONTOUR:
        if((obj->domain.core == NULL) || (obj->domain.ctr->model == NULL))
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  dim = WlzGMModelGetDimension(obj->domain.ctr->model, &errNum);
	}
	break;
      case WLZ_POINTS:
        if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(obj->domain.pts->type)
	  {
	    case WLZ_POINTS_2I: /* FALLTHROUGH */
	    case WLZ_POINTS_2D:
	      dim = 2;
	      break;
	    case WLZ_POINTS_3I: /* FALLTHROUGH */
	    case WLZ_POINTS_3D:
	      dim = 3;
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;


	  }
	}
	break;
      case WLZ_AFFINE_TRANS:
	if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  dim = WlzAffineTransformDimension(obj->domain.t, &errNum);
	}
        break;
      case WLZ_COMPOUND_ARR_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_ARR_2:
	{
	  int		idx;
	  WlzCompoundArray *ca;

	  ca = (WlzCompoundArray *)obj;
	  if(ca->n < 0)
	  {
	    errNum = WLZ_ERR_OBJECT_NULL;
	  }
	  else
	  {
	    dim = WlzObjectDimension(ca->o[0], &errNum);
	    for(idx = 1; (errNum == WLZ_ERR_NONE) && (idx < ca->n); ++idx)
	    {
	      int	d;

	      d = WlzObjectDimension(ca->o[idx], &errNum);
	      if((errNum == WLZ_ERR_NONE) && (d != dim))
	      {
		errNum = WLZ_ERR_OBJECT_TYPE;
	      }
	    }
	  }
	}
	break;
      case WLZ_MESH_TRANS:
	if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(obj->domain.core->type)
	  {
	    case WLZ_TRANSFORM_2D_MESH:
	      dim = 2;
	      break;
	    case WLZ_TRANSFORM_2D5_MESH: /* FALLTHROUGH */
	    case WLZ_TRANSFORM_3D_MESH:
	      dim = 3;
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	}
        break;
      case WLZ_CMESH_TRANS:
	if(obj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(obj->domain.core->type)
	  {
	    case WLZ_TRANSFORM_2D_CMESH:
	      dim = 2;
	      break;
	    case WLZ_TRANSFORM_2D5_CMESH: /* FALLTHROUGH */
	    case WLZ_TRANSFORM_3D_CMESH:
	      dim = 3;
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	}
        break;
	/* Objects which aren't supported. */
      case WLZ_COMPOUND_LIST_1: /* FALLTHROUGH */
      case WLZ_COMPOUND_LIST_2: /* FALLTHROUGH */
      case WLZ_FMATCHOBJ:       /* FALLTHROUGH */
      case WLZ_WARP_TRANS:      /* FALLTHROUGH */
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dim);
}

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int           option,
                ok = 1,
                usage = 0,
		dim = 0,
		transform = 0;
  WlzIBox3	bBox[2];
  WlzEffFormat  outFmt = WLZEFF_FORMAT_NONE;
  WlzEffFormat  inFmt[2] = {WLZEFF_FORMAT_NONE};
  WlzObject     *outObj = NULL;
  WlzObject	*inObj[2] = {NULL};
  char          *outFileStr;
  char		*inFileStr[2] = {NULL};
  WlzAffineTransform *tr = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  static char   optList[] = "s:t:o:u:hT",
                fileStrDef[] = "-";

  opterr = 0;
  inFileStr[0] = fileStrDef;
  inFileStr[1] = fileStrDef;
  outFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 's':
	if((inFmt[0] = WlzEffStringExtToFormat(optarg)) == 0)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 't':
	if((inFmt[1] = WlzEffStringExtToFormat(optarg)) == 0)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'u':
	if((outFmt = WlzEffStringExtToFormat(optarg)) == 0)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'T':
	transform = 1;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    usage = (optind + 2 != argc);
  }
  if(usage == 0)
  {
    inFileStr[0] = *(argv + optind);
    inFileStr[1] = *(argv + optind + 1);
  }
  if(usage == 0)
  {
    if(inFmt[0] == WLZEFF_FORMAT_NONE)
    {
      inFmt[0] = WlzEffStringFormatFromFileName(inFileStr[0]);
    }
    if(inFmt[1] == WLZEFF_FORMAT_NONE)
    {
      inFmt[1] = WlzEffStringFormatFromFileName(inFileStr[1]);
    }
    if(outFmt == WLZEFF_FORMAT_NONE)
    {
      outFmt = WlzEffStringFormatFromFileName(outFileStr);
    }
    if((inFmt[0] == WLZEFF_FORMAT_NONE) || (inFmt[1] == WLZEFF_FORMAT_NONE) ||
       (outFmt == WLZEFF_FORMAT_NONE))
    {
      usage = 1;
    }
  }
  ok = (usage == 0);
  /* Read source and target objects. */
  if(ok)
  {
    int		idx;

    for(idx = 0; ok && (idx < 2); ++idx)
    {
      char *fStr;
      FILE *fP;

      if(strcmp(inFileStr[idx], "-") == 0)
      {
        fP = stdin;
	fStr = NULL;
      }
      else
      {
        fP = NULL;
	fStr = inFileStr[idx];
      }
      if((inObj[idx] = WlzAssignObject(WlzEffReadObj(fP, fStr, inFmt[idx],
                                       0, 0, 0, &errNum), NULL)) == NULL)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	    "%s: Failed to read object from file(s) %s (%s)\n",
	    *argv, inFileStr[idx], errMsg);
      }
    }
  }
  /* Check object dimension */
  if(ok)
  {
    int		idx;
    int		d[2];

    d[0] = 0;
    for(idx = 0; ok && (idx < 2); ++idx)
    {
      d[idx] = WlzObjectDimension(inObj[idx], &errNum);
      ok = (errNum == WLZ_ERR_NONE);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if((d[0] == 0) || (d[0] != d[1]))
      {
	ok = 0;
        errNum = WLZ_ERR_OBJECT_TYPE;
      }
    }
    dim = d[0];
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to establish object dimensions (%s)\n",
	  *argv, errMsg);
    }
  }
  /* Compute the bounding box of the source and target objects. */
  if(ok)
  {
    int		idx;

    for(idx = 0; ok && (idx < 2); ++idx)
    {
      bBox[idx] = WlzBoundingBox3I(inObj[idx], &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	    "%s: Failed to compute bounding box of object read from file\n"
	    "%s (%s).\n",
	    *argv, inFileStr[idx], errMsg);
      }
    }
  }
  /* Compute affine transform which will register the source to target
   * bounding box. */
  if(ok)
  {
    int		idx,
    		nVtx;
    WlzVertexType vType;
    WlzTransformType tType;
    WlzVertexP vP[2];
    WlzDVertex3	wSp[48];

    vP[0].d3 = wSp;
    vP[1].d3 = wSp + 24;
    if(dim == 2)
    {
      nVtx = 4;
      vType = WLZ_VERTEX_D2;
      tType = WLZ_TRANSFORM_2D_AFFINE;
      for(idx = 0; idx < 2; ++idx)
      {
        (vP[idx].d2)[0].vtX = bBox[idx].xMin;
	(vP[idx].d2)[0].vtY = bBox[idx].yMin;
        (vP[idx].d2)[1].vtX = bBox[idx].xMax;
	(vP[idx].d2)[1].vtY = bBox[idx].yMin;
        (vP[idx].d2)[2].vtX = bBox[idx].xMax;
	(vP[idx].d2)[2].vtY = bBox[idx].yMax;
        (vP[idx].d2)[3].vtX = bBox[idx].xMin;
	(vP[idx].d2)[3].vtY = bBox[idx].yMax;
      }
    }
    else /* dim == 3 */
    {
      nVtx = 8;
      vType = WLZ_VERTEX_D3;
      tType = WLZ_TRANSFORM_3D_AFFINE;
      for(idx = 0; idx < 2; ++idx)
      {
        (vP[idx].d3)[0].vtX = bBox[idx].xMin;
	(vP[idx].d3)[0].vtY = bBox[idx].yMin;
	(vP[idx].d3)[0].vtZ = bBox[idx].zMin;
        (vP[idx].d3)[1].vtX = bBox[idx].xMax;
	(vP[idx].d3)[1].vtY = bBox[idx].yMin;
	(vP[idx].d3)[1].vtZ = bBox[idx].zMin;
        (vP[idx].d3)[2].vtX = bBox[idx].xMax;
	(vP[idx].d3)[2].vtY = bBox[idx].yMax;
	(vP[idx].d3)[2].vtZ = bBox[idx].zMin;
        (vP[idx].d3)[3].vtX = bBox[idx].xMin;
	(vP[idx].d3)[3].vtY = bBox[idx].yMax;
	(vP[idx].d3)[3].vtZ = bBox[idx].zMin;
        (vP[idx].d3)[4].vtX = bBox[idx].xMin;
	(vP[idx].d3)[4].vtY = bBox[idx].yMin;
	(vP[idx].d3)[4].vtZ = bBox[idx].zMax;
        (vP[idx].d3)[5].vtX = bBox[idx].xMax;
	(vP[idx].d3)[5].vtY = bBox[idx].yMin;
	(vP[idx].d3)[5].vtZ = bBox[idx].zMax;
        (vP[idx].d3)[6].vtX = bBox[idx].xMax;
	(vP[idx].d3)[6].vtY = bBox[idx].yMax;
	(vP[idx].d3)[6].vtZ = bBox[idx].zMax;
        (vP[idx].d3)[7].vtX = bBox[idx].xMin;
	(vP[idx].d3)[7].vtY = bBox[idx].yMax;
	(vP[idx].d3)[7].vtZ = bBox[idx].zMax;
      }
    }
    tr = WlzAffineTransformLSq(vType, nVtx, vP[1], nVtx, vP[0], 0, NULL,
                               tType, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	  "%s: Failed to compute affine transform (%s).\n",
	  *argv, errMsg);
    }
  }
  /* Create affine transform object or affine transform the source object
   * as required. */
  if(ok)
  {
    int		idx,
    		idy;
    const double eps = 1.0e-12;

    /* Tidy up the transform, |t_{ij}| < eps => t_{ij} = 0.0. */
    for(idy = 0; idy < 4; ++idy)
    {
      for(idx = 0; idx < 4; ++idx)
      {
        if(fabs(tr->mat[idy][idx]) < eps)
	{
	  tr->mat[idy][idx] = 0.0;
	}
      }
    }
    if(transform)
    {
      outObj = WlzAffineTransformObj(inObj[0], tr, WLZ_INTERPOLATION_NEAREST,
                                     &errNum);
      if(outObj)
      {
        tr = NULL;
      }
    }
    else
    {
      WlzDomain dom;
      WlzValues val;

      dom.t = tr;
      val.core = NULL;
      outObj = WlzMakeMain(WLZ_AFFINE_TRANS, dom, val, NULL, NULL, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to %saffine transform object (%s).\n",
		     *argv, (transform)? "create ": "", errMsg);
    }
  }
  /* Write the output object. */
  if(ok)
  {
    char *fStr;
    FILE *fP;

    if(strcmp(outFileStr, "-") == 0)
    {
      fP = stdout;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = outFileStr;
    }
    errNum = WlzEffWriteObj(fP, fStr, outObj, outFmt);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	  "%s: Failed to write object to file %s (%s)\n",
	  *argv, outFileStr, errMsg);
    }
  }
  (void )WlzFreeObj(inObj[0]);
  (void )WlzFreeObj(inObj[1]);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    char *fmtStr = NULL;

    fmtStr = WlzEffFormatTable(2, 50, 10, NULL);
    (void )fprintf(
        stderr,
	"Usage: %s%s%s%s%s%s%s%s\n",
	*argv,
	" [-h] [-T]\n"
	"\t\t[-o<out file>] [-s<src fmt>] [-t<tgt fmt>] [-u<out fmt>]\n"
	"\t\t<source object> <target object>\n"
        "Computes the Woolz affine transform which makes the bounding box of\n"
	"the source object equal to that of the target object.\n"
	"Version: ",
	WlzVersion(),
	"\n"
	"  -h Help, prints this usage information.\n"
	"  -o Output file.\n"
	"  -s Source file format.\n"
	"  -t Target file format.\n"
	"  -u Output file format.\n"
	"  -T Transform the source object and write it to the output file\n"
	"     rather than the affine transform.\n"
        "The known file formats are:\n"
	"  Description                                       Extension\n"
	"  ***********                                       *********\n",
	fmtStr,
	"Simple example:\n  ",
	*argv,
	"-o out.wlz small.vtk big.stl\n"
	"Reads the source object from the file small.vtk and the target\n"
	"object from big.stl then computes the affine transform which scales\n"
	"the source objects bounding box to fit that of the target object.\n"
	"The output transform is written to the file out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
