#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFDefGrdExportVTK_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzExtFF/WlzExtFFDefGrdExportVTK.c
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
* \brief	Exports tensor deformation features as computed by
* 		WlzDGTensorFeatures(3) to a legacy ASCII VTK format file.
* \ingroup	BinWlzExtFF
*
* \par Binary
* \ref wlzextffdefgrdexportvtk "WlzExtFFDefGrdExportVTK"
*/

/*!
\ingroup BinWlzExtFF
\defgroup wlzextffdefgrdexportvtk WlzExtFFDefGrdExportVTK
\par Name
WlzExtFFDefGrdExportVTK - exports tensor deformation features to an
                          extended VTK format file.
\par Synopsis
\verbatim
WlzExtFFDefGrdExportVTK [-h] [-m #] [-o<output file>] [input file]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
  <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr> <tr>
    <td><b>-m</b></td>
    <td>Value multiplication factor.</td>
  </tr> <tr>
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
</table>
Reads a tensor deformation gradient features compound object file
as produced by WlzDefGrdTensorFeatures(1) or WlzDGTensorFeatures(3) and
exports the features to an extended legacy ASCII VTK format file.
If given all output values are multiplied by the given factor.
The feature to export is identified by a single chaacter and must
be one of the following:
<table width="500" border="0">
  <tr>
  <td><b>ID Character</b></td><td><b>Feature</b></td>
  <td><b>j</b></td>           <td>Determinant of the Jacobian matrix
                                  (default)</td>
  <td><b>e</b></td>           <td>eigen system</td>
</table>
Here jacobian refers to the determinant of the Jacobian matrix and
eigen system to the three eigen values followed by the first two
eigen vectors (where the eigen vectors are sorted by increasing
eigen value).
For the Jacobians the input object must contain Jacobians and the
VTK POINT_DATA output is of the form:
\verbatim
  POINT_DATA <n>
  SCALARS jacobian float
  LOOKUP_TABLE default
  <j0>
  ...
\endverbatim
For the eigen system input object must contain both eigen values and
eigen vectors; the POINT_DATA output is of the form:"
\verbatim
  POINT_DATA <n>
  FIELD eigen_system float 1 6
  <l00> <l01> <l02> <e001> ... <e012>
  ...
\endverbatim
Here:
<table width="500" border="0">
  <tr><td><b>n</b></td>   <td>number of points</td></tr>
  <tr><td><b>j0</b></td>  <td>Jacobian of at the first point</td></tr>
  <tr><td><b>l0i</b></td> <td>i'th eigen value at the first point</td></tr>
  <tr><td><b>e0ij</b></td><td>j'th component of the i'th eigen vector
                              at the first point</td></tr>
</table>
Only the first two eigen vectors are output, the remaining eigen
vector, is redundant since it is orthogonal to the first two.
By default files are read from the standard input and written to the
standard output.
\par Examples
\verbatim
WlzExtFFDefGrdExportVTK -o features.vtk features.wlz
\endverbatim
Reads tensor deformation gradient features from features.wlz and
writes the Jacobians to an ASCII VTK file features.vtk.
\par File
WlzExtFFDefGrdExportVTK.c "WlzExtFFDefGrdExportVTK.c"
\par See Also
\ref BinWlzExtFF "WlzIntro(1)"
\ref wlzextffconvert "WlzExtFFConvert(1)"
\ref wlzdefgrdtensorfeatures "WlzDefGrdTensorFeatures(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Wlz.h>
#include <WlzExtFF.h>

extern int      		getopt(
				  int argc,
				  char * const *argv,
				  const char *optstring);

static WlzErrorNum		WlzExtFFDefGrdJacScale(
				  WlzObject *jObj,
				  double mf);
static WlzErrorNum		WlzExtFFDefGrdWriteEigen(
				  FILE *fP,
				  WlzObject *lObj,
				  WlzObject *eObj,
				  double mf);
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             		main(
				  int argc,
				  char **argv)
{
  int           option,
                ok = 1,
                usage = 0,
      		feature = 'j';
  double	mulFac = 1.0;
  WlzCompoundArray *inCpd = NULL;
  char		*inFileStr,
  		*outFileStr;
  FILE		*fP = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  static char   optList[] = "hf:m:o:",
                fileStrDef[] = "-";

  opterr = 0;
  inFileStr = fileStrDef;
  outFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'f':
        switch(*optarg)
	{
	  case 'e': /* FALLTHROUGH */
	  case 'j':
	    feature = *optarg;
	    break;
	  default:
	    usage = 1;
	    break;
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'm':
        if(sscanf(optarg, "%lg", &mulFac) != 1)
	{
	  usage = 1;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    usage = (optind + 1 > argc);
  }
  if(!usage && (argc > optind))
  {
    inFileStr = *(argv + optind);
  }
  ok = !usage;
  /* Read input label image. */
  if(ok)
  {
    WlzObject *inObj = NULL;

    errNum = WLZ_ERR_FILE_OPEN;
    if(((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) != NULL) &&
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) != NULL) &&
       (errNum == WLZ_ERR_NONE))
    {
      inCpd = (WlzCompoundArray *)inObj;
      inObj = NULL;
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      /* Check object. */
      if(inCpd == NULL)
      {
        errNum = WLZ_ERR_OBJECT_NULL;
      }
      else if(((inCpd->type != WLZ_COMPOUND_ARR_1) &&
               (inCpd->type != WLZ_COMPOUND_ARR_2)) ||
	      (inCpd->n < 1))
      {
        errNum = WLZ_ERR_OBJECT_TYPE;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file(s) %s (%s)\n",
	             *argv, inFileStr, errMsg);

    }
  }
  if(ok)
  {
    /* Open the output file. */
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      errNum = WLZ_ERR_FILE_OPEN;
    }
  }
  if(ok)
  {
    /* Write VTK file header and point locations to VTK file. */
    int		idN;
    WlzObject	*obj = NULL;

    for(idN = 0; idN < inCpd->n; ++idN)
    {
      if(((obj = inCpd->o[idN]) != NULL) &&
         (obj->type == WLZ_POINTS))
      {
        break;
      }
      obj = NULL;
    }
    if(obj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else
    {
      errNum = WlzEffWritePointsVtk(fP, obj, 1);
    }
  }
  if(ok)
  {
    /* Find appropriate object and write out the values. */
    int		idN;
    int		featIndex[3]  = {-1, -1, -1};
    char	*featName[3] = {
    			      "jacobian",
			      "eigen vectors",
			      "eigen values"
			    };

    for(idN = 0; idN < 3; ++idN)
    {
      int		idO;
      WlzObject		*obj = NULL;

      for(idO = 0; idO < inCpd->n; ++idO)
      {
	
	if(((obj = inCpd->o[idO]) != NULL) &&
	   (obj->type == WLZ_POINTS) &&
	   (WlzPropertyListContainsName(obj->plist, featName[idN])) != NULL)
	{
	  featIndex[idN] = idO;
	  break;
	}
      }
    }
    errNum = WLZ_ERR_VALUES_DATA;
    switch(feature)
    {
      case 'e':
	if((featIndex[1] >= 0) && (featIndex[2] >= 0))
	{
	  errNum = WlzExtFFDefGrdWriteEigen(fP,
	                               inCpd->o[featIndex[2]],
				       inCpd->o[featIndex[1]],
				       mulFac);
	}
        break;
      case 'j':
	if(featIndex[0] >= 0)
	{
	  const	double eps = 1.0e-06;

	  errNum = WLZ_ERR_NONE;
	  if(fabs(mulFac - 1.0) > eps)
	  {
	    errNum = WlzExtFFDefGrdJacScale(inCpd->o[featIndex[0]], mulFac);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzEffWritePointsVtkScalarValues(fP,
						      inCpd->o[featIndex[0]]);
	  }
	}
        break;
      default:
        break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write vlues to file %s (%s)\n",
	             *argv, outFileStr, errMsg);
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP);
  }
  (void )WlzFreeObj((WlzObject *)inCpd);
  if(usage)
  {
    (void )fprintf(
        stderr,
	"Usage: %s [-h] [-f <feat>] [-m #] [-o<output file>]\n"
	"\t\t[input file]\n"
	"Version: %s\n"
	"Options:\n"
	"  -f  Feature to export.\n"
	"  -o  Output file.\n"
	"  -m  Value multiplication factor.\n"
	"  -h  Help, prints usage message.\n"
        "Exports tensor deformation features to an extended VTK format file.\n"
	"Reads a tensor deformation gradient features compound object file\n"
	"as produced by WlzDefGrdTensorFeatures(1) or WlzDGTensorFeatures(3)\n"
	"and exports the features to an extended legacy ASCII VTK format\n"
	"file.\n"
	"If given all output values are multiplied by the given factor.\n"
	"The feature to export is identified by a single chaacter and must\n"
	"be one of the following:\n"
	"  *id char* *feature*\n"
        "   j         Determinant of the Jacobian matrix (default)\n"
	"   e         Eigen system\n"
	"Here jacobian refers to the determinant of the Jacobian matrix and\n"
	"eigen system to the three eigen values followed by the first two\n"
	"eigen vectors (where the eigen vectors are sorted by increasing\n"
	"eigen value).\n"
	"For the Jacobians the input object must contain Jacobians and the\n"
	"VTK POINT_DATA output is of the form:\n"
	"  POINT_DATA <n>\n"
	"  SCALARS jacobian float\n"
	"  LOOKUP_TABLE default\n"
	"  <j0>\n"
	"  ...\n"
	"For the eigen system the eigen system input object must contain\n"
	"both eigen values and eigen vectors; the POINT_DATA output is of\n"
	"the form:\n"
	"  POINT_DATA <n>\n"
	"  FIELD eigen_system float 1 6\n"
	"  <l00> <l01> <l02> <e001> ... <e012>\n"
	"...\n"
	"Here:\n"
	"  n    number of points\n"
	"  j0   Jacobian of at the first point\n"
	"  l0i  i'th eigen value at the first point\n"
	"  e0ij j'th component of the i'th eigen vector at the first point\n"
	"Only the first two eigen vectors are output, the remaining eigen\n"
	"vector, is redundant since it is orthogonal to the first two.\n"
        "By default files are read from the standard input and written to\n"
	"the standard output.\n"
	"Example:\n"
	"  %s -o features.vtk features.wlz\n"
	"Reads tensor deformation gradient features from features.wlz and\n"
	"writes the Jacobians to an ASCII VTK file features.vtk.\n",
	argv[0],
	WlzVersion(),
	argv[0]);
  }
  return(!ok);
}

/*!
* \return	Woolz error code.
* \brief	Scales the given Jacobian values in place.
* \param	jObj		Jacobian object.
* \param	mf		Value multiplication factor.
*/
static WlzErrorNum		WlzExtFFDefGrdJacScale(
				  WlzObject *jObj,
				  double mf)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(jObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(jObj->type != WLZ_POINTS)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(jObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(jObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  {
    size_t	idP,
  		nPts;
    WlzPointValues *jV;

    jV = jObj->values.pts;
    nPts = jObj->domain.pts->nPoints;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idP = 0; idP < nPts; ++idP)
    {
      double	*j;

      j = (double *)WlzPointValueGet(jV, idP);
      j[0] = j[0] * mf;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Writes the "eigen system" values to the given file.
* \param	fP		File pointer opened for write.
* \param	lObj		Eigen value object.
* \param	eObj		Eigen vector object.
* \param	mf		Value multiplication factor.
*/
static WlzErrorNum		WlzExtFFDefGrdWriteEigen(
				  FILE *fP,
				  WlzObject *lObj,
				  WlzObject *eObj,
				  double mf)
{
  size_t	nPts;
  WlzValues 	nVal;
  WlzObject	*nObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nVal.core = NULL;
  if((lObj == NULL) || (eObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((lObj->type != WLZ_POINTS) || (lObj->type != eObj->type))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((lObj->domain.core == NULL) || (eObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((lObj->values.core == NULL) || (eObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((nPts = lObj->domain.pts->nPoints) != eObj->domain.pts->nPoints)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if((lObj->values.pts->rank != 1) || (lObj->values.pts->dim[0] != 3) ||
	  (eObj->values.pts->rank != 2) || (eObj->values.pts->dim[0] != 3) ||
	  (eObj->values.pts->dim[1] != 3))
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    int		dim[1] = {9};

    nVal.pts = WlzMakePointValues(nPts, 1, dim, WLZ_GREY_DOUBLE, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nObj = WlzMakeMain(WLZ_POINTS, lObj->domain, nVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzProperty	prop = {0};
    WlzPropertyList *pList = NULL;

    if((pList = WlzMakePropertyList(NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      prop.name = WlzMakeNameProperty("eigen system", &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
                               WlzFreePropertyListEntry) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      nObj->plist = WlzAssignPropertyList(pList, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(pList->list == NULL)
      {
        WlzFreeProperty(prop);
      }
      else
      {
        (void )WlzFreeProperty(prop);
	(void )WlzFreePropertyList(pList);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    size_t	idP;
    WlzPointValues *lV,
                   *eV,
		   *nV;

    lV = lObj->values.pts;
    eV = eObj->values.pts;
    nV = nObj->values.pts;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idP = 0; idP < nPts; ++idP)
    {
      double	*l,
      		*e,
		*n;

      l = (double *)WlzPointValueGet(lV, idP);
      e = (double *)WlzPointValueGet(eV, idP);
      n = (double *)WlzPointValueGet(nV, idP);
      n[0] = l[0] * mf; n[1] = l[1] * mf; n[2] = l[2] * mf; 
      n[3] = e[0] * mf; n[4] = e[1] * mf; n[5] = e[2] * mf; 
      n[6] = e[3] * mf; n[7] = e[4] * mf; n[8] = e[5] * mf; 
    }
    errNum = WlzEffWritePointsVtkFieldValues(fP, nObj);
  }
  if(nObj)
  {
    (void )WlzFreeObj(nObj);
    nObj = NULL;
  }
  else if(nVal.core)
  {
    (void )WlzFreePointValues(nVal.pts);
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
