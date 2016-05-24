#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFitPlane_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzFitPlane.c
* \author       Bill Hill
* \date         May 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Computes a best fit plane to given vertices.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzfitplane "WlzFitPlane"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzfitplane WlzFitPlane
\par Name
WlzFitPlane - Calculates the best (least squares) plane through the given
              input vertices.
\par Synopsis
\verbatim
WlzFitPlane [-o <output file>] [-h]
            [-a] [-A <alg>]  [-T] [-u<x>,<y>,<z>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>Output is an affine transform.</td>
  </tr>
  <tr>
    <td><b>-A</b></td>
    <td>Force the algorithm selection, valid algorithms are:
      <ul>
        <li>SVD - Singular Value Decomposition based least squares
	          (the default).</li>
      </ul></td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Output is 3D view struct (default).</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Output is a text description of the plane.</td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Use test data, probably only useful for debugging.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>Up vector, the default is (0, 0, 0) which implies that
        the up vector is not explicitly set.</td>
  </tr>
  <tr>
    <td><b>-w</b></td>
    <td>Input is a Woolz object rather than ascii vertices.</td>
  </tr>
</table>
  
\par Description
WlzFitPlane computes a least squares best fit plane through the given
input vertices or object. 
Text output is of the form:
\verbatim
  <nx> <ny> <nz> <fx> <fy> <fz> <pitch> <yaw>
\endverbatim
where these are the normal components, centroid location in the plane
and the Euler angles (in degrees).
The input vertices are read from an ascii file with the format:
\verbatim
  <vx> <vy> <vz>
\endverbatim
The input data are read from stdin and the output object is written
to stdout unless the filenames are given.

\par Examples
\verbatim
cat tst.num
181 62 39  
105 192 97  
281 185 282  
342 81 235  
118 139 54  
221 146 173  
312 147 274  
225 69 94  
190 148 142  
141 208 154 
WlzFitPlane tst.num -t
-0.593865 -0.592398 0.544417 211.6 137.7 154.4 57.0152 224.929
\endverbatim
Note that although the output plane is a least squares solution is is not
unique (eg pitch=123.0, yaw=45.0 are also valid solutions).

\par File
\ref WlzFitPlane.c "WlzFitPlane.c"
\par See Also
\ref wlzfacts "WlzFacts(1)"
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
\ref WlzFitPlaneSVD "WlzFitPlaneSVD(3)"
\ref Wlz3DViewStructFromNormal "Wlz3DViewStructFromNormal(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

/*!
* \enum		_WlzFitPlaneAlg
* \brief 	Encodes the selected algorithm.
*/
typedef enum _WlzFitPlaneAlg
{
  WLZ_FITPLANE_ALG_SVD	 		/*!< Use WlzFitPlaneSVD(). */
} WlzFitPlaneAlg;

/*!
* \enum		_WlzFitPlaneOut
* \brief 	Encodes the selected output object type.
*/
typedef enum _WlzFitPlaneOut
{
  WLZ_FITPLANE_OUT_AFFINE,		/*!< Output an affine transform. */
  WLZ_FITPLANE_OUT_SECTION,		/*!< Output a 3D view struct. */
  WLZ_FITPLANE_OUT_TEXT			/*!< Output text. */
} WlzFitPlaneOut;

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0,
      		nVx = 0,
		inputWlz = 0,
		testFlg = 0;
  WlzVertexP	vxp;
  WlzDVertex3	up,
  		nrm,
		pip;
  WlzObject	*outObj = NULL;
  WlzVertexType vtxType = WLZ_VERTEX_D3;
  WlzObjectType	outObjType = WLZ_NULL;
  WlzDomain	outDomain;
  WlzValues	outValues;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzFitPlaneAlg alg = WLZ_FITPLANE_ALG_SVD;
  WlzFitPlaneOut out = WLZ_FITPLANE_OUT_SECTION;
  char		*inFileStr = NULL,
		*outFileStr = NULL;
  const char	*errMsg;
  static char	optList[] = "ahstTwA:o:u:",
		fileStrDef[] = "-";

  /* These vertices correspond to a plane in EMA27 with
   * fp=0,0,0, dist=123, pitch=123, yaw=45 which should give
   * the following text output:
   * -0.593865 -0.592398 0.544417 211.6 137.7 154.4 57.0152 224.929
   */
  const int	     testVxCount = 10;
  static WlzDVertex3 testVx[10] =
  {
    {181,62,39},
    {105,192,97},
    {281,185,282},
    {342,81,235},
    {118,139,54},
    {221,146,173},
    {312,147,274},
    {225,69,94},
    {190,148,142},
    {141,208,154}
  };

  opterr = 0;
  vxp.v = NULL;
  up.vtX = 0.0;
  up.vtY = 0.0;
  up.vtZ = 0.0;
  outDomain.core = NULL;
  outValues.core = NULL;
  inFileStr = fileStrDef;
  outFileStr = fileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'A':
	if(WlzStringMatchValue((int *)&alg, optarg,
			       "SVD",  WLZ_FITPLANE_ALG_SVD,
			       NULL) == 0)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'a':
        out = WLZ_FITPLANE_OUT_AFFINE;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 's':
        out = WLZ_FITPLANE_OUT_SECTION;
	break;
      case 't':
        out = WLZ_FITPLANE_OUT_TEXT;
	break;
      case 'T':
        testFlg = 1;
	break;
      case 'u':
        if(sscanf(optarg, "%lg,%lg,%lg", &(up.vtX), &(up.vtY), &(up.vtZ)) < 3)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'w':
        inputWlz = 1;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
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
      inFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    if(testFlg)
    {
      nVx = testVxCount;
      vxp.d3 = testVx;
    }
    else
    {
      FILE	*fP = NULL;

      if((*inFileStr == '\0') ||
	 ((fP = (strcmp(inFileStr, "-")?
		fopen(inFileStr, "r"): stdin)) == NULL))
      {
	ok = 0;
	errNum = WLZ_ERR_READ_EOF;
      }
      else
      {
        if(inputWlz)
	{
	  WlzObject *inObj = NULL;

	  if((inObj = WlzReadObj(fP, &errNum)) != NULL)
	  {
	    WlzVertexP inVxP;
	    WlzVertexType inVxType;

	    inVxP = WlzVerticesFromObj(inObj, NULL, &nVx, &inVxType, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      switch(inVxType)
	      {
	        case WLZ_VERTEX_I3: /* FALLTHROUGH */
		case WLZ_VERTEX_F3: /* FALLTHROUGH */
		case WLZ_VERTEX_D3:
		  vtxType = inVxType;
		  vxp.v = inVxP.v;
		  inVxP.v = NULL;
		  break;
	        default:
		  errNum = WLZ_ERR_OBJECT_TYPE;
		  break;
	      }
	      AlcFree(inVxP.v);
	    }
	    (void )WlzFreeObj(inObj);
	  }
	}
	else
	{
	  size_t nM = 0,
		 nN = 0;
	  double **inData = NULL;

	  if((AlcDouble2ReadAsci(fP, &inData, &nM, &nN) != ALC_ER_NONE) ||
	     (nM < 1) || (nN != 3))
	  {
	    ok = 0;
	    errNum = WLZ_ERR_PARAM_DATA;
	  }
	  (void )fclose(fP);
	  if(ok)
	  {
	    nVx = nM;
	    if((vxp.d3 = (WlzDVertex3 *)
			 AlcMalloc(nVx * sizeof(WlzDVertex3))) == NULL)
	    {
	      ok = 0;
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(ok)
	  {
	    int	idx;

	    for(idx = 0; idx < nVx; ++idx)
	    {
	      vxp.d3[idx].vtX = inData[idx][0];
	      vxp.d3[idx].vtY = inData[idx][1];
	      vxp.d3[idx].vtZ = inData[idx][2];
	    }
	  }
	  (void )AlcDouble2Free(inData);
	}
	if(strcmp(inFileStr, "-"))
	{
	  (void )fclose(fP);
	}
      }
      if(!ok)
      {
        (void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: failed to input vertices from file %s (%s).\n",
		       *argv, inFileStr, errMsg);
      }
    }
  }
  if(ok)
  {
    switch(alg)
    {
      case WLZ_FITPLANE_ALG_SVD:
	errNum = WlzFitPlaneSVD(vtxType, nVx, vxp, &pip, &nrm);
	if(errNum == WLZ_ERR_NONE)
	{
	  outDomain.vs3d = Wlz3DViewStructFromNormal(nrm, pip, up, &errNum);
	}
	break;
      default:
	errNum = WLZ_ERR_PARAM_TYPE;
        break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(out)
      {
        case WLZ_FITPLANE_OUT_AFFINE:
	  outObjType = WLZ_AFFINE_TRANS;
	  {
	    WlzAffineTransform *tr = NULL;

	    errNum = WlzInit3DViewStructAffineTransform(outDomain.vs3d);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      tr = WlzAffineTransformCopy(outDomain.vs3d->trans, &errNum);
	      (void )WlzFree3DViewStruct(outDomain.vs3d);
	      outDomain.t = tr;
	    }
	  }
	  break;
        case WLZ_FITPLANE_OUT_SECTION:
	  outObjType = WLZ_3D_VIEW_STRUCT;
	  break;
	case WLZ_FITPLANE_OUT_TEXT:
	  outObjType = WLZ_NULL;
	  break;
        default:
	  errNum = WLZ_ERR_PARAM_TYPE;
	  break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute best fit plane (%s).\n",
		     *argv, errMsg);
    }
  }
  if((testFlg == 0) && (vxp.v != NULL))
  {
    AlcFree(vxp.v);
  }
  if(ok && (outObjType != WLZ_NULL))
  {
    outObj = WlzMakeMain(outObjType, outDomain, outValues,
    		        NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to make output object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) != NULL)
    {
      if(outObjType == WLZ_NULL)
      {
	double	theta,
		phi;

	theta = outDomain.vs3d->theta * 180.0 / ALG_M_PI;
	phi =   outDomain.vs3d->phi   * 180.0 / ALG_M_PI;
	while(theta < 0.0)
	{
	  theta += 360.0;
	}
	while(phi < 0.0)
	{
	  phi += 360.0;
	}

        if(fprintf(fP, "%g %g %g %g %g %g %g %g\n",
	           nrm.vtX, nrm.vtY, nrm.vtZ,
		   pip.vtX, pip.vtY, pip.vtZ,
		   phi,
		   theta) > 16)
        {
	  errNum = WLZ_ERR_NONE;
	}
      }
      else
      {
        errNum = WlzWriteObj(fP, outObj);
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output object "
		     "to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(outObj)
  {
    WlzFreeObj(outObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%s",
    *argv,
    " [-a] [-A <alg>] [-o <output file>] [-h]\n"
    "\t\t[-s] [-t] [-T] [-u<x>,<y>,<z>] [<input file>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -A  Force the algorithm selection, valid algorithms are:\n"
    "        SVD  Singular Value Decomposition based least squares\n"
    "             (the default).\n"
    "  -a  Output is an affine transform.\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file name.\n"
    "  -s  Output is a 3D view struct (default).\n"
    "  -t  Output is a text description of the plane.\n"
    "  -T  Use test data, probably only useful for debugging.\n"
    "  -u  Up vector, the default is (0, 0, 0) which implies that\n"
    "      the up vector is not explicitly set.\n"
    "  -w  Input is a Woolz object rather than ascii vertices.\n"
    "Calculates the best (least squares) plane through the given input\n"
    "vertices or object.\n"
    "Text output is of the form:\n"
    "  <nx> <ny> <nz> <cx> <cy> <cz> <pitch> <yaw>\n"
    "where these are the normal components, centroid location in the plane\n"
    "and the Euler angles (in degrees).\n"
    "The input vertices are read from an ascii file with the format:\n"
    "  <vtx x> <vtx y> <vtx z>\n"
    "The input data are read from stdin and the output object is written\n"
    "to stdout unless the filenames are given.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
