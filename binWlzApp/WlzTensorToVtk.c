#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTensorToVtk_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTensorToVtk.c
* \author       Angus Murray
* \date         September 2005
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
* \brief
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <Wlz.h>
 
#define MAXY 712

/* externals required by getopt  - not in ANSI C standard */

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;

static WlzErrorNum WlzTensorToVtk(FILE *, WlzObject *, WlzObject *, 
				  WlzObject *, int);

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-h] [-t] [<tensor file>] [-v] [<vtk file>] [-n]\n"
	  "\tConverts Woolz tensor object to VTK format\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -t        tensor file\n"
	  "\t  -v        VTK file\n"
          "\t  -n        normalise tensors\n"
	  "",
	  proc_str);
  return;
}

int		main(int argc, char **argv)
{
  WlzCompoundArray *inObj;
  WlzObject        *objT11,
                   *objT12,
                   *objT22;
  WlzErrorNum      errNum = WLZ_ERR_NONE;
  FILE	           *inFile = NULL, 
                   *outFile = NULL;
  static char 	   optList[] = "t:v:hn";
  const char       *inFileStr = NULL,
                   *outFileStr = NULL,
                   *errMsg;
  int	           option,
    normalise = 0;
    
  /* read the argument list and check for an input and output file */
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option )
    {
      case 't':
        inFileStr = optarg;
	break;
      case 'v':
        outFileStr = optarg;
	break;
      case 'n':
        normalise = 1;
	break;
      case 'h':
        usage(argv[0]);
        return( WLZ_ERR_NONE );

      default:
        usage(argv[0]);
        return( WLZ_ERR_PARAM_TYPE );
    }
  }
    
  errNum = WLZ_ERR_READ_EOF;
  if((inFileStr == NULL) || (*inFileStr == '\0') ||
       ((inFile = fopen(inFileStr, "r")) == NULL) ||
       ((inObj= (WlzCompoundArray *)WlzReadObj(inFile, &errNum)) == NULL))
  {
      errNum = WLZ_ERR_PARAM_DATA ;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
  }
  if(inFile)
  {
    fclose(inFile);
  }
  if((outFileStr == NULL) || (*outFileStr == '\0') ||
       ((outFile = fopen(outFileStr, "w")) == NULL))
  {
      errNum = WLZ_ERR_PARAM_DATA;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to open file %s (%s).\n",
		     *argv, outFileStr, errMsg);
  }
  /* read  X and Y objects */
  if(errNum == WLZ_ERR_NONE)
  {
    switch(inObj->type)
    {
      case WLZ_COMPOUND_ARR_1:
	if (inObj->n == 3)
	{
	  if ((objT11 = inObj->o[0]) != NULL)
	  {
	    if ((objT22 = inObj->o[1]) != NULL)
	    {
	      if ((objT12 = inObj->o[2]) == NULL)
	      {
		objT11 = NULL;
		objT22 = NULL;
		errNum = WLZ_ERR_OBJECT_NULL;
	      }
	    }
	    else
	    {
	      objT11 = NULL;
	      errNum = WLZ_ERR_OBJECT_NULL;
	    }
	  }
	  else
	    errNum = WLZ_ERR_OBJECT_NULL;
	}
	else
	  errNum = WLZ_ERR_OBJECT_DATA;
        break;

      default:
	errNum = WLZ_ERR_OBJECT_DATA;
	break;
    }
  }
  /* extract tensor components and store in a VTK file format */
  if (errNum == WLZ_ERR_NONE)
  {
    errNum = WlzTensorToVtk(outFile, objT11, objT22, objT12, normalise);
  }

  if (outFile)
  {
    fclose(outFile);
  }
  
  return errNum;
}

/*!
* \return	WlzErrorNum.
* \ingroup	WlzTransform
* \brief	Outputs tensors to a file in VTK format.
* \param	outFile			VTK format file.
* \param        objT11                  tensor component t11 object.
* \param        objT22                  tensor component t22 object.
* \param        objT12                  tensor component t12 object.
*/
static WlzErrorNum WlzTensorToVtk(FILE *outFile, WlzObject *objT11, 
				  WlzObject *objT22, WlzObject *objT12,
				  int normalise)
{
  WlzIntervalWSpace     iWspT11,
                        iWspT22,
                        iWspT12;
  WlzGreyWSpace	        gWspT11,
                        gWspT22,
                        gWspT12;
  WlzGreyP	        gPixT11,
                        gPixT22,
                        gPixT12;
  WlzErrorNum           errNum = WLZ_ERR_NONE;
  int                   numpoints,
                        k,
                        l;
  AlgMatrix		aM;
  double                **aA,
                        vM[2];
  float                 norm,
                        eccentricity;
 
  aM.core = NULL;

  /* initialise workspaces */
  if ((errNum = WlzInitGreyScan(objT11, &iWspT11, &gWspT11)) == WLZ_ERR_NONE)
  {
    if ((errNum = WlzInitGreyScan(objT22, &iWspT22, &gWspT22)) == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(objT12, &iWspT12, &gWspT12);
    }
  }

  /* Calculate number of points */
  if (errNum == WLZ_ERR_NONE)
  {
    numpoints = 0;
    while((errNum = WlzNextGreyInterval(&iWspT11)) == WLZ_ERR_NONE)
    {
      for (k = iWspT11.lftpos; k <= iWspT11.rgtpos; k++)
      {
        numpoints++;
      }
    }
    if(errNum == WLZ_ERR_EOO)        
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if (errNum == WLZ_ERR_NONE)
  {
    fprintf(outFile, "# vtk DataFile Version 4.0\n");
    fprintf(outFile, "Tensors of an MAPaint mapping\n");
    fprintf(outFile, "ASCII\n");
    fprintf(outFile, "DATASET POLYDATA\n");
    fprintf(outFile, "%s %d %s\n","POINTS",numpoints,"  float");
  }
  /* now output pixel locations */
  if (errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(objT11, &iWspT11, &gWspT11);
  }
  if (errNum == WLZ_ERR_NONE)
  {
    while((errNum = WlzNextGreyInterval(&iWspT11)) == WLZ_ERR_NONE)
    {
      l = iWspT11.linpos;
      for (k = iWspT11.lftpos; k <= iWspT11.rgtpos; k++)
      {
        fprintf(outFile, "%lf  %lf  0.0\n", (float)k, (float)l);
      }
    }
    if(errNum == WLZ_ERR_EOO)        
    {
      errNum = WLZ_ERR_NONE;
    }   
  }
  /* output tensors vectors */
  if (errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(objT11, &iWspT11, &gWspT11);
  }
  if (errNum == WLZ_ERR_NONE)
  {
    fprintf(outFile,  "%s %d\n","POINT_DATA", numpoints);
    fprintf(outFile, "%s %s %s\n","TENSORS","tensors", "float" );
    while(((errNum = WlzNextGreyInterval(&iWspT11)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iWspT22)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iWspT12)) == WLZ_ERR_NONE))
    {
      gPixT11 = gWspT11.u_grintptr;
      gPixT22 = gWspT22.u_grintptr;
      gPixT12 = gWspT12.u_grintptr;
      l = iWspT11.linpos;
      for (k = iWspT11.lftpos; k <= iWspT11.rgtpos; k++)
      {
	if (normalise)
	{
	  norm = *(gPixT11.flp) * (*(gPixT11.flp));
	  norm += *(gPixT22.flp) * (*(gPixT22.flp));
	  norm += *(gPixT12.flp) * (*(gPixT12.flp));
	  norm += *(gPixT12.flp) * (*(gPixT12.flp));
	  norm = sqrt(norm);
	}
	else
	  norm = 1.0;
        fprintf(outFile, "%lf  %lf 0.0\n", *(gPixT11.flp)/norm, *(gPixT12.flp)/norm);
        fprintf(outFile, "%lf  %lf 0.0\n", *(gPixT12.flp)/norm, *(gPixT22.flp)/norm);
	fprintf(outFile, "0.0 0.0 0.0\n");
	++(gPixT11.flp);
	++(gPixT22.flp);
	++(gPixT12.flp);
      }
    }
    if(errNum == WLZ_ERR_EOO)        
    {
      errNum = WLZ_ERR_NONE;
    }   
  }
  /* output scalars */
  if (errNum == WLZ_ERR_NONE)
  {
    if ((errNum = WlzInitGreyScan(objT11, &iWspT11, &gWspT11)) == WLZ_ERR_NONE)
    {
      if ((errNum = WlzInitGreyScan(objT22, &iWspT22, &gWspT22)) == 
	  WLZ_ERR_NONE)
      {
	errNum = WlzInitGreyScan(objT12, &iWspT12, &gWspT12);
      }
    }
  }
  if (errNum == WLZ_ERR_NONE)
  {
    if ((aM.rect = AlgMatrixRectNew(2, 2, NULL)) != NULL)
    {
      aA = aM.rect->array;
      fprintf(outFile, "%s %s %s\n","SCALARS","scalars", "float" );
      fprintf(outFile, "%s %s\n","LOOKUP_TABLE","default");
      while(((errNum = WlzNextGreyInterval(&iWspT11)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&iWspT22)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&iWspT12)) == WLZ_ERR_NONE))
      {
	gPixT11 = gWspT11.u_grintptr;
	gPixT22 = gWspT22.u_grintptr;
	gPixT12 = gWspT12.u_grintptr;
	l = iWspT11.linpos;
	for (k = iWspT11.lftpos; k <= iWspT11.rgtpos; k++)
	{
	  if (normalise)
	  {
	    norm = *(gPixT11.flp) * (*(gPixT11.flp));
	    norm += *(gPixT22.flp) * (*(gPixT22.flp));
	    norm += *(gPixT12.flp) * (*(gPixT12.flp));
	    norm += *(gPixT12.flp) * (*(gPixT12.flp));
	    norm = sqrt(norm);
	  }
	  else
	    norm = 1.0;
	  aA[0][0] = *(gPixT11.flp) / norm;
	  aA[0][1] = *(gPixT12.flp) / norm;
	  aA[1][0] = *(gPixT12.flp) / norm;

	  AlgMatrixRSEigen(aM, vM, 1);
	  if (vM[0] < 0)
	    vM[0] = -vM[0];
	  if (vM[1] < 0)
	    vM[1] = -vM[1];
	  if (vM[0] > vM[1])
	    eccentricity = sqrt(1.0 - vM[1]/vM[0]);
	  else
	    eccentricity = sqrt(1.0 - vM[0]/vM[1]);
	  fprintf(outFile, "%lf\n ",eccentricity);
	  
	  ++(gPixT11.flp);
	  ++(gPixT22.flp);
	  ++(gPixT12.flp);
	}
      }
      AlgMatrixFree(aM);
    }
  }

  if(errNum == WLZ_ERR_EOO)        
  {
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
} 

