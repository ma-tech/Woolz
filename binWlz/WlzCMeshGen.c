#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzMeshGen_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzCMeshGen.c
* \author       Bill Hill
* \date         January 2009
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
* \brief	Constructs a 2D or 3D conforming simplical mesh from
* 		a domain object.
*
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \ref 		wlzcmeshgen "WlzCMeshGen"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshgen WlzCMeshGen
\par Name
WlzCMeshGen - computes a conforming mesh that covers the given domain object.
\par Synopsis
\verbatim
WlzCMeshGen [-h] [-o<output file>]
            [-L#] [-a#] [-W#] [-l#] [-u#] [-B] [-C]
	    [-m#] [-M#] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Number of Laplacian smoothing iterations.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Laplacian alpha parameter.</td>
  </tr>
  <tr> 
    <td><b>-W</b></td>
    <td>Number of low pass filter smoothing iterations.</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Low pass filter \f$\lambda\f$ value.</td>
  </tr>
  <tr> 
    <td><b>-u</b></td>
    <td>Low pass filter \f$\mu\f$ value.</td>
  </tr>
  <tr> 
    <td><b>-B</b></td>
    <td>Smooth boundary (requires a smoothing method to
        be selected).</td>
  </tr>
  <tr> 
    <td><b>-C</b></td>
    <td>Don't make the mesh conform to the object's domain.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Initial minimum mesh element size.</td>
  </tr>
  <tr> 
    <td><b>-M</b></td>
    <td>Maximum mesh element size.</td>
  </tr>
</table>
\par Description
Computes a conforming mesh that covers the given domain object.
\par Examples
\verbatim
WlzCMeshGen -m 10 -M 30 -o mesh.wlz dom.wlz
\endverbatim
Creates a new conforming mesh object which covers the input domain
read from dom.wlz. The mesh will have maximum edge length 30 and
minimum edge length 10 before it is made to conform to the objects
domain.
\par File
\ref WlzCMeshGen.c "WlzCMeshGen.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshFromObj "WlzCMeshFromObj(3)"
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
  int		ok = 1,
  		option,
  		usage = 0,
		conform = 1,
		laplacianItr = 0,
		lowPassItr = 0,
		smoothBnd = 0;
  double	laplacianAlpha = 0.1,
		lowPassLambda = 0.33,
		lowPassMu = -0.34,
  		minElmSz = 25.0,
  		maxElmSz = 100.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*obj = NULL,
  		*outObj = NULL;
  WlzCMeshP 	mesh;
  WlzObjectType	cMType;
  static char   optList[] = "a:BChl:L:m:M:o:u:W";
  const char    inObjFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  val.core = NULL;
  mesh.v = NULL;
  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'a':
	if(sscanf(optarg, "%lg", &laplacianAlpha) != 1)
	{
	  usage = 1;
	}
        break;
      case 'B':
        smoothBnd = 1;
	break;
      case 'C':
        conform = 0;
	break;
      case 'l':
	if(sscanf(optarg, "%lg", &lowPassLambda) != 1)
	{
	  usage = 1;
	}
        break;
      case 'L':
	if(sscanf(optarg, "%d", &laplacianItr) != 1)
	{
	  usage = 1;
	}
        break;
      case 'm':
        if(sscanf(optarg, "%lg", &minElmSz) != 1)
	{
	  usage = 1;
	}
	break;
      case 'M':
        if(sscanf(optarg, "%lg", &maxElmSz) != 1)
	{
	  usage = 1;
	}
        break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'u':
	if(sscanf(optarg, "%lg", &lowPassMu) != 1)
	{
	  usage = 1;
	}
        break;
      case 'W':
	if(sscanf(optarg, "%d", &lowPassItr) != 1)
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
  ok = usage == 0;
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    (void )WlzAssignObject(obj, NULL);
    mesh = WlzCMeshFromObj(obj, minElmSz, maxElmSz, NULL, conform, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to create conforming mesh, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    WlzCMeshSetBoundNodFlags(mesh);
  }
  if(ok && (laplacianItr > 0))
  {
    errNum = WlzCMeshLaplacianSmooth(mesh, laplacianItr, laplacianAlpha,
    				     smoothBnd, 1);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to Laplacian smooth mesh, %s.\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok && (lowPassItr > 0))
  {
    errNum = WlzCMeshLPFilterLM(mesh, lowPassLambda, lowPassMu,
    				lowPassItr, smoothBnd, 1);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to low pass filter mesh, %s.\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
	     fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_2D:
	cMType = WLZ_CMESH_2D;
	dom.cm2 = mesh.m2;
	break;
      case WLZ_CMESH_3D:
	cMType = WLZ_CMESH_3D;
	dom.cm3 = mesh.m3;
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      outObj = WlzMakeMain(cMType, dom, val, NULL, NULL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzWriteObj(fP, outObj);
    }
    (void )WlzFreeObj(outObj);
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<output file>]\n"
	    "       [-L#] [-a#] [-W#] [-l#] [-u#] [-B] [-C]\n"
	    "       [-m#] [-M#] [<input object>]\n"
    	    "Computes a conforming mesh that covers the given domain object.\n"
	    "Options are:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -o  Output file.\n"
	    "  -a  Laplacian alpha parameter.\n"
	    "  -W  Number of low pass filter smoothing iterations.\n"
	    "  -l  Low pass filter lambda value.\n"
	    "  -u  Low pass filter mu value.\n"
	    "  -B  Smooth boundary (requires a smoothing method to be\n"
	    "      selected).\n"
	    "  -C  Don't make the mesh conform to the object's domain.\n"
	    "  -m  Minimum mesh element size.\n"
	    "  -L  Number of Laplacian smoothing iterations.\n"
	    "  -M  Maximum mesh element size.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
