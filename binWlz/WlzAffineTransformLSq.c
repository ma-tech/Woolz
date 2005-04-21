#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzAffineTransformLSq.c
* \author       Bill Hill
* \date         June 2004
* \version      $Id$
* \note
*               Copyright
*               2004 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Computes an affine transform from a list of vertices
*		and vertex displacements in the 2D format
*		  <vtx x> <vtx y> <delta x> <delta y>
*		or the 3D format
*		  <vtx x> <vtx y> <vtx z> <delta x> <delta y> <delta z>
*		Comment lines start with a '#' character.
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

/*!
* \enum		_WlzAffineTransLSqAlg
* \brief 	Encodes the selected algorithm.
*/
typedef enum _WlzAffineTransLSqAlg
{
  WLZ_AFFINETRANSLSQ_ALG_DEFAULT,	/*!< Algorithm default. */
  WLZ_AFFINETRANSLSQ_ALG_WLZ,	 	/*!< Woolz algorithm for 2D */
  WLZ_AFFINETRANSLSQ_ALG_DQ     	/*!< Walker's dual quaternion
  					 *   algorithm 2D or 3D rigid body
					 *   with optional normal vectors
					 *   and weights */
} WlzAffineTransLSqAlg;

#define IN_RECORD_MAX	(1024)

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		tI0,
  		option,
		ok = 1,
		usage = 0,
      		inR,
		inC,
		iVtx = 0,
      		nV = 0,
      		nN = 0,
		vtxSz,
		testFlg = 0,
		testVtxCount;
  double	*vWgt = NULL,
  		*nWgt = NULL;
  double	**inData = NULL;
  WlzVertexType	vtxType = WLZ_VERTEX_D2;
  WlzObject	*trObj = NULL;
  WlzDomain	trDomain;
  WlzValues	trValues;
  WlzVertexP	vS,
  		vT,
		nS,
		nT;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzTransformType trType = WLZ_TRANSFORM_2D_AFFINE;
  WlzAffineTransLSqAlg alg = WLZ_AFFINETRANSLSQ_ALG_DEFAULT;
  FILE		*fP = NULL;
  char 		*outObjFileStr = NULL,
  		*inFileStr = NULL,
		*wgtFileStr = NULL;
  const char	*errMsg;
  static char	optList[] = "A:o:w:23Tanrth",
		outObjFileStrDef[] = "-",
  		inFileStrDef[] = "-";
  const WlzDVertex2 testVx2D[4] =
  {
    {-100.0, -100.0},
    {-100.0,  100.0},
    { 100.0, -100.0},
    { 100.0,  100.0}
  };
  const WlzDVertex3 testVx3D[8] =
  {
    {-100.0, -100.0, -100.0},
    {-100.0, -100.0,  100.0},
    {-100.0,  100.0, -100.0},
    {-100.0,  100.0,  100.0},
    { 100.0, -100.0, -100.0},
    { 100.0, -100.0,  100.0},
    { 100.0,  100.0, -100.0},
    { 100.0,  100.0,  100.0}
  };

  opterr = 0;
  vS.v = vT.v = NULL;
  nS.v = nT.v = NULL;
  trDomain.core = NULL;
  trValues.core = NULL;
  outObjFileStr = outObjFileStrDef;
  inFileStr = inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
        vtxType = WLZ_VERTEX_D2;
	break;
      case '3':
        vtxType = WLZ_VERTEX_D3;
	break;
      case 'A':
	if(WlzStringMatchValue((int *)&alg, optarg,
			       "WLZ", WLZ_AFFINETRANSLSQ_ALG_WLZ,
			       "DQ",  WLZ_AFFINETRANSLSQ_ALG_DQ,
			       NULL) == 0)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'T':
        testFlg = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'a':
        trType = WLZ_TRANSFORM_2D_AFFINE;
	break;
      case 'n':
        trType = WLZ_TRANSFORM_2D_NOSHEAR;
	break;
      case 'r':
        trType = WLZ_TRANSFORM_2D_REG;
	break;
      case 't':
        trType = WLZ_TRANSFORM_2D_TRANS;
	break;
      case 'w':
        wgtFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inFileStr == NULL) || (*inFileStr == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
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
      inFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    if(vtxType == WLZ_VERTEX_D2)
    {
      testVtxCount = 4;
      vtxSz = sizeof(WlzDVertex2);
    }
    else /* vtxType == WLZ_VERTEX_D3 */
    {
      testVtxCount = 8;
      vtxSz = sizeof(WlzDVertex3);
      switch(trType)
      {
	case WLZ_TRANSFORM_2D_AFFINE:
	  trType = WLZ_TRANSFORM_3D_AFFINE;
	  break;
	case WLZ_TRANSFORM_2D_NOSHEAR:
	  trType = WLZ_TRANSFORM_3D_NOSHEAR;
	  break;
	case WLZ_TRANSFORM_2D_REG:
	  trType = WLZ_TRANSFORM_3D_REG;
	  break;
	case WLZ_TRANSFORM_2D_TRANS:
	  trType = WLZ_TRANSFORM_3D_TRANS;
	  break;
      }
    }
  }
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
	      fopen(inFileStr, "r"): stdin)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_READ_EOF;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to open input file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
  }
  if(ok)
  {
    if(testFlg)
    {
      nV = testVtxCount;
      if(((vS.v = AlcMalloc(nV * vtxSz)) == NULL) ||
	 ((vT.v = AlcMalloc(nV * vtxSz)) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if((trObj = WlzAssignObject(
			WlzReadObj(fP, &errNum), NULL)) == NULL)
        {
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else if(trObj->type != WLZ_AFFINE_TRANS)
	{
	  errNum = WLZ_ERR_OBJECT_TYPE;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(vtxType == WLZ_VERTEX_D2)
	{
	  (void )memcpy(vS.v, testVx2D, nV * vtxSz);
	  iVtx = 0;
	  while((errNum == WLZ_ERR_NONE) && (iVtx < nV))
	  {
	    *(vT.d2 + iVtx) = WlzAffineTransformVertexD2(
	    					trObj->domain.t,
						*(vS.d2 + iVtx),
						&errNum);
	    ++iVtx;
	  }
	}
	else /* vtxType == WLZ_VERTEX_D2 */
	{
	  (void )memcpy(vS.v, testVx3D, nV * vtxSz);
	  iVtx = 0;
	  while((errNum == WLZ_ERR_NONE) && (iVtx < nV))
	  {
	    *(vT.d3 + iVtx) = WlzAffineTransformVertexD3(
	    					trObj->domain.t,
						*(vS.d3 + iVtx),
						&errNum);
	    ++iVtx;
	  }
	}
      }
      if(trObj)
      {
        WlzFreeObj(trObj);
	trObj = NULL;
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to create test data (%s).\n",
		       *argv, inFileStr, errMsg);
      }
    }
    else
    {
      if(AlcDouble2ReadAsci(fP, &inData, (size_t *) &inR, (size_t *) &inC) != ALC_ER_NONE)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read input file %s.\n",
		       *argv, inFileStr);
      }
      if(fP && strcmp(inFileStr, "-"))
      {
	fclose(fP);
      }
      if(ok)
      {
        if(vtxType == WLZ_VERTEX_D2)
	{
	  if(inC != 4)
	  {
	    ok = 0;
	    (void )fprintf(stderr,
		   "%s: invalid format in input file %s, the required\n"
		   "format is:\n"
		   "  <vtx x> <vtx y> <disp x> <disp y>\n",
		   *argv, inFileStr);
	  }
	}
	else
	{
	  if(inC != 6)
	  {
	    ok = 0;
	    (void )fprintf(stderr,
		   "%s: invalid format in input file %s, the required\n"
		   "format is:\n"
		   "  <vtx x> <vtx y> <vtx z> <disp x> <disp y> <disp z>\n",
		   *argv, inFileStr);
	  }
	}
      }
      if(ok)
      {
	nV = inR;
	if(((vS.v = AlcMalloc(nV * vtxSz)) == NULL) ||
	   ((vT.v = AlcMalloc(nV * vtxSz)) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	if(vtxType == WLZ_VERTEX_D2)
	{
	  for(iVtx = 0; iVtx < nV; ++iVtx)
	  {
	    (vS.d2 + iVtx)->vtX = inData[iVtx][0];
	    (vS.d2 + iVtx)->vtY = inData[iVtx][1];
	    (vT.d2 + iVtx)->vtX = inData[iVtx][0] + inData[iVtx][2];
	    (vT.d2 + iVtx)->vtY = inData[iVtx][1] + inData[iVtx][3];
	  }
	}
	else
	{
	  for(iVtx = 0; iVtx < nV; ++iVtx)
	  {
	    (vS.d3 + iVtx)->vtX = inData[iVtx][0];
	    (vS.d3 + iVtx)->vtY = inData[iVtx][1];
	    (vS.d3 + iVtx)->vtZ = inData[iVtx][2];
	    (vT.d3 + iVtx)->vtX = inData[iVtx][0] + inData[iVtx][3];
	    (vT.d3 + iVtx)->vtY = inData[iVtx][1] + inData[iVtx][4];
	    (vT.d3 + iVtx)->vtZ = inData[iVtx][2] + inData[iVtx][5];
	  }
	}
      }
      AlcDouble2Free(inData);
      if(wgtFileStr)
      {
	fP = NULL;
        if(ok)
	{
          if((fP = (strcmp(wgtFileStr, "-")?
	      fopen(wgtFileStr, "r"): stdin)) == NULL)
	  {
	    ok = 0;
	    errNum = WLZ_ERR_READ_EOF;
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
			   "%s: failed to open input weights file %s (%s).\n",
			   *argv, wgtFileStr, errMsg);
	  }
	}
	if(ok)
	{
	  if((AlcDouble1ReadAsci(fP, &vWgt, (size_t *) &tI0) != ALC_ER_NONE) ||
	     (tI0 < nV))
	  {
	    ok = 0;
	    (void )fprintf(stderr,
			   "%s: failed to read weights file %s.\n",
			   *argv, wgtFileStr);
	  }
	}
	if(fP && strcmp(wgtFileStr, "-"))
	{
	  fclose(fP);
	}
      }
    }
  }
  if(ok)
  {
    switch(alg)
    {
      case WLZ_AFFINETRANSLSQ_ALG_DEFAULT:
	trDomain.t = WlzAffineTransformLSq(vtxType,
					   nV, vT, nV, vS, nV, vWgt,
					   trType, &errNum);
	break;
      case WLZ_AFFINETRANSLSQ_ALG_WLZ:
	trDomain.t = WlzAffineTransformLSqRegWlz2D(vT.d2, vS.d2, nV,
					           &errNum);
        break;
      case WLZ_AFFINETRANSLSQ_ALG_DQ:
	switch(vtxType)
	{
	  case WLZ_VERTEX_D2:
	    trDomain.t = WlzAffineTransformLSqDQ2D(
	    				 nV, vWgt, vT.d2, vS.d2,
					 nN, nWgt, nT.d2, nS.d2,
					 &errNum);
	    break;
	  case WLZ_VERTEX_D3:
	    trDomain.t = WlzAffineTransformLSqDQ3D(
	    				 nV, vWgt, vT.d3, vS.d3,
					 nN, nWgt, nT.d3, nS.d3,
					 &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_TYPE;
	    break;
	}
        break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute transform (%s).\n",
		     *argv, errMsg);
    }
  }
  if(vS.v)
  {
    AlcFree(vS.v);
  }
  if(vT.v)
  {
    AlcFree(vT.v);
  }
  if(vWgt)
  {
    AlcFree(vWgt);
  }
  if(ok)
  {
    trObj = WlzMakeMain(WLZ_AFFINE_TRANS, trDomain, trValues,
    		        NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to make affine transform object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, trObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output affine transform object "
		     "to file %s (%s).\n",
		     *argv, outObjFileStr, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(trObj)
  {
    WlzFreeObj(trObj);
  }
  else if(trDomain.core)
  {
    WlzFreeDomain(trDomain);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<out obj>] [-A<algorithm>]\n"
    "                             [-2] [-3] [-T] [-a] [-r] [-t] [-h]\n"
    "                             [-w<weights>] [<in data>]\n"
    "Options:\n"
    "  -2  Compute a 2D transform, requires verticies in 2D format.\n"
    "  -3  Compute a 3D transform, requires verticies in 3D format.\n"
    "  -A  Force the algorithm selection, valid algorithms are:\n"
    "        WLZ  Woolz algorithm for 2D rigid body.\n"
    "        DQ   Walker's dual quaternion algorithm for 2D and 3D\n"
    "             rigid body.\n"
    "  -T  Test, probably only useful for debugging.\n"
    "  -o  Output transform object file name.\n"
    "  -a  Compute affine transform.\n"
    "  -n  Compute no-shear transform, translate, rotate, scale only.\n"
    "  -r  Compute registration transform, affine but no scale or shear.\n"
    "  -t  Compute translation transform.\n"
    "  -w  Weights file with a weight value for each vertex displacement\n"
    "      pair.\n"
    "  -h  Help, prints this usage message.\n"
    "Calculates the best (least squares) affine transform from the given\n"
    "input vertex/vertex displacement list.\n"
    "The input verticies/displacements are read from an ascii file with\n"
    "the 2D format:\n"
    "  <vtx x> <vtx y> <disp x> <disp y>\n"
    "or the 3D format:\n"
    "  <vtx x> <vtx y> <vtx z> <disp x> <disp y> <disp z>\n"
    "The input data are read from stdin and the transform object is written\n"
    "to stdout unless the filenames are given.\n",
    *argv);
  }
  return(!ok);
}
