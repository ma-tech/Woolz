#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFDefGrdExportVTK_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzExtFFDefGrdExportVTK.c
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
WlzExtFFDefGrdExportVTK - exports tensor deformation features to a VTK format
file.
\par Synopsis
\verbatim
WlzExtFFDefGrdExportVTK [-h] [-o<output file>] [input file]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
  <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
</table>
Reads a tensor deformation gradient features compound object file
as produced by WlzDefGrdTensorFeatures(1) or WlzDGTensorFeatures(3) and
exports the features to a legacy ASCII VTK format file.
By default files are read from the standard input and written to the
standard output.
\par Examples
\verbatim
WlzExtFFDefGrdExportVTK -o features.vtk features.wlz
\endverbatim
Reads a tensor deformation gradient features from features.wlz and
writes them to an ASCII VTK file features.vtk.
\par File
WlzExtFFDefGrdExportVTK.c "WlzExtFFDefGrdExportVTK.c"
\par See Also
\ref BinWlzExtFF "WlzIntro(1)"
\ref wlzextffconvert "WlzExtFFConvert(1)"
\ref wlzdefgrdtensorfeatures "WlzDefGrdTensorFeatures(1)"
\ref wlzdgtensorfeatures "WlzDGTensorFeatures(3)"
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
                usage = 0;
  WlzCompoundArray *inCpd = NULL;
  char		*inFileStr,
  		*outFileStr;
  FILE		*fP = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  static char   optList[] = "ho:",
                fileStrDef[] = "-";

  opterr = 0;
  inFileStr = fileStrDef;
  outFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
	outFileStr = optarg;
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
    char	*names[3] = {
    			      "jacobian",
			      "direction vectors",
			      "stretch values"
			    };

    for(idN = 0; idN < 3; ++idN)
    {
      int		idO;
      WlzObject		*obj = NULL;

      for(idO = 0; idO < inCpd->n; ++idO)
      {
	if(((obj = inCpd->o[idO]) != NULL) &&
	   (obj->type == WLZ_POINTS) &&
	   WlzPropertyListContainsName(obj->plist, names[idN]))
	{
	  break;
	}
	obj = NULL;
      }
      if(obj)
      {
	switch(idN)
	{
	  case 0:
	    errNum = WlzEffWritePointsVtkScalarValues(fP, obj);
	    break;
	  case 1: /* FALLTHROUGH */
	  case 2:
	    errNum = WlzEffWritePointsVtkFieldValues(fP, obj);
	    break;
	  default:
	    break;
	}
      }
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
	"Usage: %s [-h] [-o<output file>] [input file]\n"
	"Version: %s\n"
	"Options:\n"
	"  -o  Output file.\n"
	"  -h  Help, prints usage message.\n"
        "Exports tensor deformation features to a VTK format file.\n"
	"Reads a tensor deformation gradient features compound object file\n"
	"as produced by WlzDefGrdTensorFeatures(1) or WlzDGTensorFeatures(3)\n"
	"and exports the features to a legacy ASCII VTK format file.\n"
        "By default files are read from the standard input and written to\n"
	"the standard output.\n"
	"Example:\n"
	"  %s -o features.vtk features.wlz\n"
	"Reads a tensor deformation gradient features from features.wlz and\n"
	"writes them to an ASCII VTK file features.vtk.\n",
	argv[0],
	WlzVersion(),
	argv[0]);
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
