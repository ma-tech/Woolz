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
      		vtxCount = 0,
      		vtxLimit = 0,
		vtxSz,
		testFlg = 0,
		testVtxCount;
  WlzVertexType	vtxType = WLZ_VERTEX_D2;
  WlzObject	*trObj = NULL;
  WlzDomain	trDomain;
  WlzValues	trValues;
  WlzVertexP	vtx0,
  		vtx1,
		vtxVec0,
  		vtxVec1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzTransformType trType = WLZ_TRANSFORM_2D_AFFINE;
  char 		*rec,
  		*outObjFileStr,
  		*inFileStr;
  char		inRecord[IN_RECORD_MAX];
  const char	*errMsg;
  static char	optList[] = "o:23Tanrth",
		outObjFileStrDef[] = "-",
  		inFileStrDef[] = "-";
  const WlzDVertex2 testVx2D[4] =
  {
    {1.0, 1.0},
    {1.0, 2.0},
    {2.0, 1.0},
    {2.0, 2.0}
  };
  const WlzDVertex3 testVx3D[8] =
  {
    {1.0, 1.0, 1.0},
    {1.0, 1.0, 2.0},
    {1.0, 2.0, 1.0},
    {1.0, 2.0, 2.0},
    {2.0, 1.0, 1.0},
    {2.0, 1.0, 2.0},
    {2.0, 2.0, 1.0},
    {2.0, 2.0, 2.0}
  };

  opterr = 0;
  vtxVec0.v = NULL;
  vtxVec1.v = NULL;
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
      vtxLimit = testVtxCount;
      if(((vtxVec0.v = AlcRealloc(vtxVec0.v,
				  vtxLimit * vtxSz)) == NULL) ||
	 ((vtxVec1.v = AlcRealloc(vtxVec1.v,
				  vtxLimit * vtxSz)) == NULL))
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
	  (void )memcpy(vtxVec0.v, testVx2D, vtxLimit * vtxSz);
	  vtxCount = 0;
	  while((errNum == WLZ_ERR_NONE) && (vtxCount < vtxLimit))
	  {
	    *(vtxVec1.d2 + vtxCount) = WlzAffineTransformVertexD2(
	    					trObj->domain.t,
						*(vtxVec0.d2 + vtxCount),
						&errNum);
	    ++vtxCount;
	  }
	}
	else /* vtxType == WLZ_VERTEX_D2 */
	{
	  (void )memcpy(vtxVec0.v, testVx3D, vtxLimit * vtxSz);
	  vtxCount = 0;
	  while((errNum == WLZ_ERR_NONE) && (vtxCount < vtxLimit))
	  {
	    *(vtxVec1.d3 + vtxCount) = WlzAffineTransformVertexD3(
	    					trObj->domain.t,
						*(vtxVec0.d3 + vtxCount),
						&errNum);
	    ++vtxCount;
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
      while((errNum == WLZ_ERR_NONE) &&
	    (fgets(inRecord, IN_RECORD_MAX - 1, fP) != NULL))
      {
	inRecord[IN_RECORD_MAX - 1] = '\0';
	rec = inRecord;
	while(*rec && isspace(*rec))
	{
	  ++rec;
	}
	if(*rec && (*rec != '#'))
	{
	  if(vtxCount >= vtxLimit)
	  {
	    vtxLimit = (vtxLimit + 1024) * 2;
	    if(((vtxVec0.v = AlcRealloc(vtxVec0.v,
					vtxLimit * vtxSz)) == NULL) ||
	       ((vtxVec1.v = AlcRealloc(vtxVec1.v,
					vtxLimit * vtxSz)) == NULL))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(vtxType == WLZ_VERTEX_D2)
	      {
		vtx0.d2 = vtxVec0.d2 + vtxCount;
		vtx1.d2 = vtxVec1.d2 + vtxCount;
	      }
	      else /* vtxType == WLZ_VERTEX_D3 */
	      {
		vtx0.d3 = vtxVec0.d3 + vtxCount;
		vtx1.d3 = vtxVec1.d3 + vtxCount;
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(vtxType == WLZ_VERTEX_D2)
	    {
	      if(sscanf(rec, "%lg %lg %lg %lg",
			&(vtx0.d2->vtX), &(vtx0.d2->vtY),
			&(vtx1.d2->vtX), &(vtx1.d2->vtY)) != 4)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		++(vtx0.d2);
		++(vtx1.d2);
		++vtxCount;
	      }
	    }
	    else /* vtxType == WLZ_VERTEX_D3 */
	    {
	      if(sscanf(rec, "%lg %lg %lg %lg %lg %lg ",
			&(vtx0.d3->vtX), &(vtx0.d3->vtY),
			&(vtx0.d3->vtZ), &(vtx1.d3->vtX),
			&(vtx1.d3->vtY), &(vtx1.d3->vtZ)) != 6)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		++(vtx0.d3);
		++(vtx1.d3);
		++vtxCount;
	      }
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  vtxLimit = vtxCount;
	  if(vtxType == WLZ_VERTEX_D2)
	  {
	    vtx0.d2 = vtxVec0.d2;
	    vtx1.d2 = vtxVec1.d2;
	    while(vtxCount-- > 0)
	    {
	      vtx1.d2->vtX += vtx0.d2->vtX;
	      vtx1.d2->vtY += vtx0.d2->vtY;
	      ++(vtx0.d2);
	      ++(vtx1.d2);
	    }
	  }
	  else /* vtxType == WLZ_VERTEX_D3 */
	  {
	    vtx0.d3 = vtxVec0.d3;
	    vtx1.d3 = vtxVec1.d3;
	    while(vtxCount-- > 0)
	    {
	      vtx1.d3->vtX += vtx0.d3->vtX;
	      vtx1.d3->vtY += vtx0.d3->vtY;
	      vtx1.d3->vtZ += vtx0.d3->vtZ;
	      ++(vtx0.d3);
	      ++(vtx1.d3);
	    }
	  }
	}
	else /* errNum != WLZ_ERR_NONE */
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: failed to read input file %s (%s).\n",
			 *argv, inFileStr, errMsg);
	}
      }
      if(fP && strcmp(inFileStr, "-"))
      {
	fclose(fP);
      }
    }
  }
  if(ok)
  {
    trDomain.t = WlzAffineTransformLSq(vtxType,  vtxLimit, vtxVec0,
				       vtxLimit, vtxVec1,
				       trType, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute transform (%s).\n",
		     *argv, errMsg);
    }
  }
  if(vtxVec0.v)
  {
    AlcFree(vtxVec0.v);
  }
  if(vtxVec1.v)
  {
    AlcFree(vtxVec1.v);
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
  if(vtxVec0.v)
  {
    AlcFree(vtxVec0.v);
  }
  if(vtxVec1.v)
  {
    AlcFree(vtxVec1.v);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<out obj>] [-2] [-3] [-a] [-r] [-t]\n"
    "                             [-h] [<in data>]\n"
    "Options:\n"
    "  -2  Compute a 2D transform, requires verticies in 2D format.\n"
    "  -3  Compute a 3D transform, requires verticies in 3D format.\n"
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
