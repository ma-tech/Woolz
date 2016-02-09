#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDisplacementToVtk_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzDisplacementToVtk.c
* \author       Angus Murray
* \date         September 2005
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
* \ingroup	BinWlzApp
* \brief	Computes VTK visualisation of displacements.
*
* \par Binary
* \ref wlzdisplacementtovtk "WlzDisplacementToVtk"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzdisplacementtovtk WlzDisplacementToVtk
\par Name
WlzDisplacementToVtk
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <Wlz.h>
 
#define MAXY 712

/* externals required by getopt  - not in ANSI C standard */

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;

static WlzErrorNum WlzDisplacementsToVtk(FILE *, WlzObject *, WlzObject *);

static void usage(char *proc_str)
{
  (void )fprintf(stderr,
	  "Usage:\t%s [-h] [-c] [<compound array file>] [-v] [<vtk file>]\n"
	  "\tConverts Compound Array object to VTK format\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -c        Compound Array file\n"
	  "\t  -v        VTK file\n",
	  proc_str,
	  WlzVersion());
  return;
}

int		main(int argc, char **argv)
{
  WlzCompoundArray *inObj;
  WlzObject        *objX,
                   *objY;
  WlzErrorNum      errNum = WLZ_ERR_NONE;
  FILE	           *inFile = NULL, 
                   *outFile = NULL;
  static char 	   optList[] = "c:v:h";
  const char       *inFileStr = NULL,
                   *outFileStr = NULL,
                   *errMsg;
  int	           option;
    
  /* read the argument list and check for an input and output file */
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option )
    {
      case 'c':
        inFileStr = optarg;
	break;
      case 'v':
        outFileStr = optarg;
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
	if (inObj->n == 2)
	{
	  if ((objX = inObj->o[0]) != 0)
	  {
	    if ((objY = inObj->o[1]) == 0)
	    {
	      objX = NULL;
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
  /* combine x and y displacements into a VTK file format */
  if (errNum == WLZ_ERR_NONE)
  {
    errNum = WlzDisplacementsToVtk(outFile, objX, objY);
  }

  if (outFile)
  {
    fclose(outFile);
  }
  
  return(errNum);
}

/*!
* \return	WlzErrorNum.
* \ingroup	WlzTransform
* \brief	Outputs displacements to a file in VTK format.
* \param	outFile			VTK format file.
* \param        objX                    X displacement object.
* \param	objY			Y displacement object.
*/
static WlzErrorNum WlzDisplacementsToVtk(FILE *outFile, 
                                         WlzObject *objX, WlzObject *objY)
{
  WlzIntervalWSpace     iWspX,
                        iWspY;
  WlzGreyWSpace	        gWspX,
                        gWspY;
  WlzGreyP	        gPixX,
                        gPixY;
  WlzErrorNum           errNum = WLZ_ERR_NONE;
  int                   numpoints,
                        k,
                        l;
  float                 m;

  /* initialise workspaces */
  if ((errNum = WlzInitGreyScan(objX, &iWspX, &gWspX)) == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(objY, &iWspY, &gWspY);
  }

  /* Calculate number of points */
  if (errNum == WLZ_ERR_NONE)
  {
    numpoints = 0;
    while((errNum = WlzNextGreyInterval(&iWspX)) == WLZ_ERR_NONE)
    {
      for (k = iWspX.lftpos; k <= iWspX.rgtpos; k++)
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
    fprintf(outFile, "Displacements of an MAPaint mapping\n");
    fprintf(outFile, "ASCII\n");
    fprintf(outFile, "DATASET POLYDATA\n");
    fprintf(outFile, "%s %d %s\n","POINTS",numpoints,"  float");
  }
  /* now output pixel locations */
  if (errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(objX, &iWspX, &gWspX);
  }
  if (errNum == WLZ_ERR_NONE)
  {
    while((errNum = WlzNextGreyInterval(&iWspX)) == WLZ_ERR_NONE)
    {
      l = iWspX.linpos;
      for (k = iWspX.lftpos; k <= iWspX.rgtpos; k++)
      {
        fprintf(outFile, "%lf  %lf  0.0\n", (float)k, (float)l);
      }
    }
    if(errNum == WLZ_ERR_EOO)        
    {
      errNum = WLZ_ERR_NONE;
    }   
  }
  /* output displacement vectors */
  if (errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(objX, &iWspX, &gWspX);
  }
  if (errNum == WLZ_ERR_NONE)
  {
    fprintf(outFile,  "%s %d\n","POINT_DATA", numpoints);
    fprintf(outFile, "%s %s %s\n","VECTORS","vectors", "float" );
    while(((errNum = WlzNextGreyInterval(&iWspX)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iWspY)) == WLZ_ERR_NONE))
    {
      gPixX = gWspX.u_grintptr;
      gPixY = gWspY.u_grintptr;
      l = iWspX.linpos;
      for (k = iWspX.lftpos; k <= iWspX.rgtpos; k++)
      {
        fprintf(outFile, "%lf  %lf  0.0\n", *(gPixX.flp), *(gPixY.flp));
	++(gPixX.flp);
	++(gPixY.flp);
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
    if ((errNum = WlzInitGreyScan(objX, &iWspX, &gWspX)) == WLZ_ERR_NONE)
      errNum = WlzInitGreyScan(objY, &iWspY, &gWspY);
  }
  if (errNum == WLZ_ERR_NONE)
  {
    fprintf(outFile, "%s %s %s\n","SCALARS","scalars", "float" );
    fprintf(outFile, "%s %s\n","LOOKUP_TABLE","default");
    while(((errNum = WlzNextGreyInterval(&iWspX)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iWspY)) == WLZ_ERR_NONE))
    {
      gPixX = gWspX.u_grintptr;
      gPixY = gWspY.u_grintptr;
      l = iWspX.linpos;
      for (k = iWspX.lftpos; k <= iWspX.rgtpos; k++)
      {
	m = sqrt(pow(*(gPixX.flp),2) + pow(*(gPixY.flp),2));
        fprintf(outFile, "%lf\n", m);
	++(gPixX.flp);
	++(gPixY.flp);
      }
    }
  }

  if(errNum == WLZ_ERR_EOO)        
  {
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
} 
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
