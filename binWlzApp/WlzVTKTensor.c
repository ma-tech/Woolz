#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzVTKTensor_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzVTKTensor.c
* \author       Bill Hill
* \date         January 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief 	Converts a Woolz conforming mesh object with tensor values
* 		to a VTK polydata file with tensor dataset attributes.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzvtktensor "WlzVTKTensor"
*/

/*!
 \ingroup BinWlzApp
 \defgroup wlzvtktensor WlzVTKTensor
 \par Name
 WlzVTKTensor - Converts a Woolz object with tensor values to a VTK polydata
                file with tensor dataset attributes.
\par Synopsis
\verbatim
WlzVTKTensor [-h] [-o<out vtk file>] [<in woolz file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output VTK file name.</td>
  </tr>
</table>
\par Description
Reads either a 3D Woolz conforming mesh object with tensor data attached
to it's elements or a points object with 3D points and tensor values,
then writes the tensor data to a VTK polydata file with points at either
the conforming mesh centroids or at given points.
By default files are read from the standard input and written to the
standard output.
\par Example
\verbatim
WlzVTKTensor -o out.vtk in.wlz
\endverbatim
Creates a VTK polydata file with tensor dataset attributes from the
given Woolz file.
\par File
\ref WlzVTKTensor.c "WlzVTKTensor.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshTrStrainTensor "WlzCMeshTrStrainTensor(1)"
*/

#ifndef  DOXYGEN_SHOULD_SKIP_THIS

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
extern char     *optarg;
extern int      optind,
                opterr,
		optopt;

static WlzErrorNum		WlzVTKTensorFromCMesh(
				  char *prog,
				  char *inFileStr,
				  FILE *fP,
				  WlzObject *inObj);
static WlzErrorNum		WlzVTKTensorFromPoints(
				  char *prog,
				  char *inFileStr,
				  FILE *fP,
				  WlzObject *inObj);

int		main(int argc, char *argv[])
{
  int		ok,
  		option,
  		usage = 0;
  char		*inFileStr,
  		*outFileStr;
  WlzObject	*inObj = NULL;
  FILE		*fP = NULL;
  const char	*errMsg;
  static char	optList[] = "o:h",
  		fileStrDef[] = "-";
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  opterr = 0;
  inFileStr = outFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
	outFileStr = optarg;
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
      inFileStr = *(argv + optind);
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if(((fP = (strcmp(inFileStr, "-")?
        fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s)\n",
		     argv[0], inFileStr, errMsg);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
  }
  if(ok)
  {
    /* Check object is valid for tensor export. */
    if(inObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(inObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else if(inObj->values.core == NULL)
    {
      errNum = WLZ_ERR_VALUES_NULL;
    }
    else
    {
      switch(inObj->type)
      {
        case WLZ_CMESH_3D:
	  if(inObj->domain.core->type  != WLZ_CMESH_3D)
	  {
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	  }
	  if(inObj->values.core->type != WLZ_INDEXED_VALUES)
	  {
	    errNum = WLZ_ERR_VALUES_TYPE;
	  }
	  else if((inObj->values.x->attach != WLZ_VALUE_ATTACH_ELM) ||
		  (inObj->values.x->vType != WLZ_GREY_DOUBLE) ||
		  (inObj->values.x->rank != 2) ||
		  (inObj->values.x->dim[0] != 3) ||
		  (inObj->values.x->dim[1] != 3))
	  {
	    errNum = WLZ_ERR_VALUES_DATA;
	  }
	  break;
        case WLZ_POINTS:
	  if(inObj->domain.core->type  != WLZ_POINTS_3D)
	  {
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	  }
	  else if(inObj->values.core->type != WLZ_POINT_VALUES)
	  {
	    errNum = WLZ_ERR_VALUES_TYPE;
	  }
	  else if((inObj->values.pts->vType != WLZ_GREY_DOUBLE) ||
	          (inObj->values.pts->rank != 2) ||
		  (inObj->values.pts->dim[0] != 3) ||
		  (inObj->values.pts->dim[1] != 3))
	  {
	    errNum = WLZ_ERR_VALUES_DATA;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Input object read from %s is inappropriate.\n"
		     "Only 3D conforming mesh objects with tensor values\n"
		     "attached to elements or points objects with tensor\n"
		     "values are aupported. (%s)\n",
		     argv[0], inFileStr, errMsg);
    }
  }
  if(ok)   
  {
    if((fP = (strcmp(outFileStr, "-")?
        fopen(outFileStr, "w"): stdin)) == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_WRITE_EOF;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s (%s)\n",
		     argv[0], outFileStr, errMsg);
    }
  }
  if(ok)   
  {
    switch(inObj->type)
    {
      case WLZ_CMESH_3D:
        errNum = WlzVTKTensorFromCMesh(argv[0], inFileStr, fP, inObj);
	break;
      case WLZ_POINTS:
        errNum = WlzVTKTensorFromPoints(argv[0], inFileStr, fP, inObj);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;   
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write tesnsor data to file %s (%s).\n",
		     argv[0], outFileStr, errMsg);
      
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP);
  }
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<out vtk file>] [-s<x>,<y>[,<z>]]\n"
    "\t\t[<in woolz file>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h        Help, prints this usage message.\n"
    "  -o        Output vtk file name.\n"
    "  -s	 Sampling interval.\n"
    "Reads either a 3D Woolz conforming mesh object with tensor data\n"
    "attached to it's elements or a points object with 3D points and\n"
    "tensor values, then writes the tensor data to a VTK polydata file\n"
    "with points at either the conforming mesh centroids or at given\n"
    "points.\n"
    "By default files are read from the standard input and written to\n"
    "the standard output.\n"
    "Example:\n"
    "  %s -o out.vtk in.wlz\n"
    "Creates a VTK polydata file with tensor dataset attributes from the\n"
    "given Woolz file.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}

/*!
* \return	Woolz error code.
* \ingroup	BinWlzApp
* \brief	Outputs a VTK tensor file from a conforming mesh object.
* \param	prog:			Program name.
* \param	inFileStr		Name of input file object was read
* 					from.
* \param	fP			Opened file.
* \param	obj			Given conforming mesh object.
*/
static WlzErrorNum WlzVTKTensorFromCMesh(char *prog, char *inFileStr,
					 FILE *fP, WlzObject *obj)
{
  int		idE,
		maxElm;
  WlzCMesh3D	*mesh;
  AlcVector	*elmVec;
  WlzIndexedValues *ixv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  mesh = obj->domain.cm3;
  ixv = obj->values.x;
  elmVec = mesh->res.elm.vec;
  maxElm = mesh->res.elm.maxEnt;
  if(fprintf(fP,
	     "# vtk DataFile Version 1.0\n"
	     "Output from Woolz WLZ_CMESH_3D file %s via %s\n"
	     "ASCII\n"
	     "DATASET POLYDATA\n"
	     "POINTS %d float\n",
	     inFileStr, prog,  maxElm) <= 0)
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < maxElm; ++idE)
    {
      WlzDVertex3 c;
      WlzCMeshElm3D *elm;
      WlzCMeshNod3D *nod[4];
      
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(elmVec, idE);
      if(elm->idx >= 0)
      {
	nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
	nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
	nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
	nod[3] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
	c.vtX = 0.25 * (nod[0]->pos.vtX + nod[1]->pos.vtX + 
			nod[2]->pos.vtX + nod[3]->pos.vtX);
	c.vtY = 0.25 * (nod[0]->pos.vtY + nod[1]->pos.vtY + 
			nod[2]->pos.vtY + nod[3]->pos.vtY);
	c.vtZ = 0.25 * (nod[0]->pos.vtZ + nod[1]->pos.vtZ + 
			nod[2]->pos.vtZ + nod[3]->pos.vtZ);
      }
      else
      {
	WLZ_VTX_3_ZERO(c);
      }
      if(fprintf(fP,
		 "%f %f %f\n", c.vtX, c.vtY, c.vtZ) <= 0)
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
	break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
	       "POINT_DATA %d\n"
	       "TENSORS tensors float\n",
	       maxElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < maxElm; ++idE)
    {
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(elmVec, idE);
      if(elm->idx >= 0)
      {
	double *t;

	t = (double *)WlzIndexedValueGet(ixv, elm->idx);
	if(fprintf(fP,
		   "%f %f %f "
		   "%f %f %f "
		   "%f %f %f \n",
		   t[0], t[1], t[2],
		   t[3], t[4], t[5],
		   t[6], t[7], t[8]) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
      else
      {
	if(fprintf(fP,
		   "0 0 0\n"
		   "0 0 0\n"
		   "0 0 0\n\n") <= 0)
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
	break;
      }
    }
  }
  return(errNum);
}


/*!
* \return	Woolz error code.
* \ingroup	BinWlzApp
* \brief	Outputs a VTK tensor file from a points object.
* \param	prog:			Program name.
* \param	inFileStr		Name of input file object was read
* 					from.
* \param	fP			Opened file.
* \param	obj			Given points object.
*/
static WlzErrorNum	WlzVTKTensorFromPoints(char *prog, char *inFileStr,
					       FILE *fP, WlzObject *obj)
{
  int		idN,
		nP;
  WlzPoints	*pD;
  WlzPointValues *pV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  pD = obj->domain.pts;
  pV = obj->values.pts;
  nP = pD->nPoints;
  if(fprintf(fP,
	     "# vtk DataFile Version 1.0\n"
	     "Output from Woolz WLZ_POINTS file %s via %s\n"
	     "ASCII\n"
	     "DATASET POLYDATA\n"
	     "POINTS %d float\n",
	     inFileStr, prog,  nP) <= 0)
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDVertex3 *p;

    p = pD->points.d3;
    for(idN = 0; idN < nP; ++idN)
    {
      if(fprintf(fP,
		 "%f %f %f\n", p->vtX, p->vtY, p->vtZ) <= 0)
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
	break;
      }
      ++p;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
	       "POINT_DATA %d\n"
	       "TENSORS tensors float\n",
	       nP) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < nP; ++idN)
    {
      double	*t;

      t = WlzPointValueGet(pV, idN);
      if(fprintf(fP,
		 "%f %f %f "
		 "%f %f %f "
		 "%f %f %f \n",
		 t[0], t[1], t[2],
		 t[3], t[4], t[5],
		 t[6], t[7], t[8]) <= 0)
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
	break;
      }
    }
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
