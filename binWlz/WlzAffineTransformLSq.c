#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzAffineTransformLSq.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes an affine transform from a list of verticies
*		and vertex displacements in the 2D format
*		  <vtx x> <vtx y> <delta x> <delta y>
*		or the 3D format
*		  <vtx x> <vtx y> <vtx z> <delta x> <delta y> <delta z>
*		Comment lines start with a '#' character.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 13-12-00 bill Change members of WlzVertex and WlzVertexP.
* 30-11-00 bill Make changes for 3D least squares affine transforms
*		and add test code.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

/*!
* \enum		_WlzAffineTransLSqAlg
* \brief 	Encodes the selected algorithm.
*/
enum _WlzAffineTransLSqAlg
{
  WLZ_AFFINETRANSLSQ_ALG_WLZ,	/*!< Woolz algorithms for 2D */
  WLZ_AFFINETRANSLSQ_ALG_SVD,	/*!< Arun's singular value decomposition
  				 *   algorithm for 3D verticies */
  WLZ_AFFINETRANSLSQ_ALG_DQ     /*!< Walker's dual quaternion algorithm for
  				 *   3D verticies with optional normal
				 *   vectors and weights */
};
typedef enum _WlzAffineTransLSqAlg WlzAffineTransLSqAlg;

#define IN_RECORD_MAX	(1024)

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
      		inR,
		inC,
		iVtx = 0,
      		nVtx = 0,
      		nNrm = 0,
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
  WlzVertexP	vtx0,
  		vtx1,
		nrm0,
		nrm1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzTransformType trType = WLZ_TRANSFORM_2D_AFFINE;
  WlzAffineTransLSqAlg alg = WLZ_AFFINETRANSLSQ_ALG_WLZ;
  FILE		*fP = NULL;
  char 		*outObjFileStr,
  		*inFileStr;
  const char	*errMsg;
  static char	optList[] = "A:o:23Tanrth",
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
  vtx0.v = vtx1.v = NULL;
  nrm0.v = nrm1.v = NULL;
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
			       "SVD", WLZ_AFFINETRANSLSQ_ALG_SVD,
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
      if(alg == WLZ_AFFINETRANSLSQ_ALG_WLZ)
      {
        alg = WLZ_AFFINETRANSLSQ_ALG_SVD;
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
      nVtx = testVtxCount;
      if(((vtx0.v = AlcMalloc(nVtx * vtxSz)) == NULL) ||
	 ((vtx1.v = AlcMalloc(nVtx * vtxSz)) == NULL))
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
	  (void )memcpy(vtx0.v, testVx2D, nVtx * vtxSz);
	  iVtx = 0;
	  while((errNum == WLZ_ERR_NONE) && (iVtx < nVtx))
	  {
	    *(vtx1.d2 + iVtx) = WlzAffineTransformVertexD2(
	    					trObj->domain.t,
						*(vtx0.d2 + iVtx),
						&errNum);
	    ++iVtx;
	  }
	}
	else /* vtxType == WLZ_VERTEX_D2 */
	{
	  (void )memcpy(vtx0.v, testVx3D, nVtx * vtxSz);
	  iVtx = 0;
	  while((errNum == WLZ_ERR_NONE) && (iVtx < nVtx))
	  {
	    *(vtx1.d3 + iVtx) = WlzAffineTransformVertexD3(
	    					trObj->domain.t,
						*(vtx0.d3 + iVtx),
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
      if(AlcDouble2ReadAsci(fP, &inData, &inR, &inC) != ALC_ER_NONE)
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
	nVtx = inR;
	if(((vtx0.v = AlcMalloc(nVtx * vtxSz)) == NULL) ||
	   ((vtx1.v = AlcMalloc(nVtx * vtxSz)) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	if(vtxType == WLZ_VERTEX_D2)
	{
	  for(iVtx = 0; iVtx < nVtx; ++iVtx)
	  {
	    (vtx0.d2 + iVtx)->vtX = inData[iVtx][0];
	    (vtx0.d2 + iVtx)->vtY = inData[iVtx][1];
	    (vtx1.d2 + iVtx)->vtX = inData[iVtx][0] + inData[iVtx][2];
	    (vtx1.d2 + iVtx)->vtY = inData[iVtx][1] + inData[iVtx][3];
	  }
	}
	else
	{
	  for(iVtx = 0; iVtx < nVtx; ++iVtx)
	  {
	    (vtx0.d3 + iVtx)->vtX = inData[iVtx][0];
	    (vtx0.d3 + iVtx)->vtY = inData[iVtx][1];
	    (vtx0.d3 + iVtx)->vtZ = inData[iVtx][2];
	    (vtx1.d3 + iVtx)->vtX = inData[iVtx][0] + inData[iVtx][3];
	    (vtx1.d3 + iVtx)->vtY = inData[iVtx][1] + inData[iVtx][4];
	    (vtx1.d3 + iVtx)->vtZ = inData[iVtx][2] + inData[iVtx][5];
	  }
	}
      }
      AlcDouble2Free(inData);
    }
  }
  if(ok)
  {
    switch(alg)
    {
      case WLZ_AFFINETRANSLSQ_ALG_WLZ: /* FALLTHROUGH */
      case WLZ_AFFINETRANSLSQ_ALG_SVD:
	trDomain.t = WlzAffineTransformLSq(vtxType,  nVtx, vtx0,
					   nVtx, vtx1, trType, &errNum);
        break;
      case WLZ_AFFINETRANSLSQ_ALG_DQ:
	trDomain.t = WlzAffineTransformLSq2(vtxType,
					   nVtx, vWgt, vtx0, vtx1,
					   nNrm, nWgt, nrm0, nrm1,
					   trType, &errNum);
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
  if(vtx0.v)
  {
    AlcFree(vtx0.v);
  }
  if(vtx1.v)
  {
    AlcFree(vtx1.v);
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
    "                             [<in data>]\n"
    "Options:\n"
    "  -2  Compute a 2D transform, requires verticies in 2D format.\n"
    "  -3  Compute a 3D transform, requires verticies in 3D format.\n"
    "  -A  Force the algorithm selection, valid algorithms are:\n"
    "        WLZ  Woolz algorithms for 2D.\n"
    "        SVD  Arun's singular value decomposition algorithm for 3D.\n"
    "        DQ   Walker's dual quaternion algorithm for 3D.\n"
    "  -T  Test, probably only useful for debugging.\n"
    "  -o  Output transform object file name.\n"
    "  -a  Compute affine transform.\n"
    "  -n  Compute no-shear transform, translate, rotate, scale only.\n"
    "  -r  Compute registration transform, affine but no scale or shear.\n"
    "  -t  Compute translation transform.\n"
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
