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
*		and vertex displacements in the format:
*		  <vertex x> <vertex y> <delta x> <delta y>
*		Comment lines start with a '#' character.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
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
      		vtxLimit = 0;
  WlzObject	*trObj = NULL;
  WlzDomain	trDomain;
  WlzValues	trValues;
  WlzDVertex2	*vtx0,
  		*vtx1,
		*vtxVec0 = NULL,
  		*vtxVec1 = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzTransformType trType = WLZ_TRANSFORM_2D_AFFINE;
  char 		*rec,
  		*outObjFileStr,
  		*inFileStr;
  char		inRecord[IN_RECORD_MAX];
  const char	*errMsg;
  static char	optList[] = "o:anrth",
		outObjFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  trDomain.core = NULL;
  trValues.core = NULL;
  outObjFileStr = outObjFileStrDef;
  inFileStr = inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
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
    if(errNum == WLZ_ERR_NONE)
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
	    if(((vtxVec0 = (WlzDVertex2 *)AlcRealloc(vtxVec0,
	    			   vtxLimit * sizeof(WlzDVertex2))) == NULL) ||
	       ((vtxVec1 = (WlzDVertex2 *)AlcRealloc(vtxVec1,
	    			   vtxLimit * sizeof(WlzDVertex2))) == NULL))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    else
	    {
	      vtx0 = vtxVec0 + vtxCount;
	      vtx1 = vtxVec1 + vtxCount;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(sscanf(rec, "%lg %lg %lg %lg", &(vtx0->vtX), &(vtx0->vtY),
		      &(vtx1->vtX), &(vtx1->vtY)) != 4)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      ++vtx0;
	      ++vtx1;
	      ++vtxCount;
	    }
	  }
	}
      }
      if(errNum != WLZ_ERR_NONE)
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
  if(ok)
  {
    vtxLimit = vtxCount;
    vtx0 = vtxVec0;
    vtx1 = vtxVec1;
    while(vtxCount-- > 0)
    {
      vtx1->vtX += vtx0->vtX;
      vtx1->vtY += vtx0->vtY;
      ++vtx0;
      ++vtx1;
    }
  }
  if(ok)
  {
    trDomain.t = WlzAffineTransformLSq(vtxLimit, vtxVec0,
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
  if(vtxVec0)
  {
    AlcFree(vtxVec0);
  }
  if(vtxVec1)
  {
    AlcFree(vtxVec1);
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
  if(vtxVec0)
  {
    AlcFree(vtxVec0);
  }
  if(vtxVec1)
  {
    AlcFree(vtxVec1);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<out obj>] [-a] [-r] [-t] [-h] [<in data>]\n"
    "Options:\n"
    "  -o  Output transform object file name.\n"
    "  -a  Compute 2D affine transform.\n"
    "  -n  Compute 2D no-shear transform, translate, rotate, scale only.\n"
    "  -r  Compute 2D registration transform, affine but no scale or shear.\n"
    "  -t  Compute 2D translation transform.\n"
    "  -h  Help, prints this usage message.\n"
    "Calculates the best (least squares) affine transform from the given\n"
    "input vertex/vertex displacement list.\n"
    "The input verticies/displacements are read from an ascii file with\n"
    "the format:\n"
    "  <vertex x> <vertex y> <displacement x> <displacement y>\n"
    "The input data are read from stdin and the transform object is written\n"
    "to stdout unless the filenames are given.\n",
    *argv);
  }
  return(!ok);
}
