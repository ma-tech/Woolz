#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:	Woolz
* Title:	WlzBuildContour.c
* Date: 	August 2001
* Author:	Bill Hill
* Copyright:	2001 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:	Builds a contour from input simplicies.
*		These are in one of the two the formats,
*		2D:
*                 x0 y0 x1 y1
*		or 3D
*                 x0 y0 z0 x1 y1 z1 x2 y2 z3
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int           idR,
  		nRec,
  		nFld,
		option,
  		ok = 1,
		usage = 0;
  FILE		*fP = NULL;
  char		*inFileStr,
  		*outObjFileStr;
  double	**bufD = NULL;
  WlzGMModelType modelType = WLZ_GMMOD_2D;
  WlzDVertex2	bufD2[2];
  WlzDVertex3	bufD3[3];
  WlzObject     *outObj = NULL;
  WlzDomain	ctrDom;
  WlzValues	dumVal;
  WlzGMModel	*model = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "23diho:";
  const char	outObjFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  ctrDom.core = NULL;
  dumVal.core = NULL;
  outObjFileStr = (char *)outObjFileStrDef;
  inFileStr = (char *)inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
	switch(modelType)
	{
	  case WLZ_GMMOD_2I: /* FALLTHROUGH */
	  case WLZ_GMMOD_2D:
	    break;
	  case WLZ_GMMOD_3I:
	    modelType = WLZ_GMMOD_2I;
	    break;
	  case WLZ_GMMOD_3D:
	    modelType = WLZ_GMMOD_2D;
	    break;
	}
	break;
      case '3':
	switch(modelType)
	{
	  case WLZ_GMMOD_2I:
	    modelType = WLZ_GMMOD_3I;
	    break;
	  case WLZ_GMMOD_2D:
	    modelType = WLZ_GMMOD_3D;
	    break;
	  case WLZ_GMMOD_3I: /* FALLTHROUGH */
	  case WLZ_GMMOD_3D:
	    break;
	}
	break;
      case 'd':
	switch(modelType)
	{
	  case WLZ_GMMOD_2I:
	    modelType = WLZ_GMMOD_2D;
	  case WLZ_GMMOD_2D:
	    break;
	  case WLZ_GMMOD_3I:
	    modelType = WLZ_GMMOD_3D;
	    break;
	  case WLZ_GMMOD_3D:
	    break;
	}
	break;
      case 'i':
	switch(modelType)
	{
	  case WLZ_GMMOD_2I:
	    break;
	  case WLZ_GMMOD_2D:
	    modelType = WLZ_GMMOD_2I;
	    break;
	  case WLZ_GMMOD_3I:
	    break;
	  case WLZ_GMMOD_3D:
	    modelType = WLZ_GMMOD_3I;
	    break;
	}
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
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
  }
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
	      fopen(inFileStr, "r"): stdin)) == NULL) ||
       (AlcDouble2ReadAsci(fP, &bufD, (size_t *) &nRec, (size_t *) &nFld) != ALC_ER_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read from file %s\n",
		     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    ctrDom.ctr = WlzMakeContour(&errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      ctrDom.ctr->model = model = WlzAssignGMModel(
      		   WlzGMModelNew(modelType, 0, 0, &errNum), NULL);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to make contour (%s).\n",
      	     	     argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    switch(modelType)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	if(nFld != 4)
	{
	  ok = 0;
	}
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	if(nFld != 9)
	{
	  ok = 0;
	}
	break;
    }
    if(ok == 0)
    {
      (void )fprintf(stderr,
      		     "%s: Failed to read simplicies (syntax error).\n",
		     argv[0]);
    }
  }
  if(ok)
  {
    idR = 0;
    while((errNum == WLZ_ERR_NONE) && (idR < nRec))
    {
      switch(modelType)
      {
        case WLZ_GMMOD_2I: /* FALLTHROUGH */
        case WLZ_GMMOD_2D:
          bufD2[0].vtX = *(*(bufD + idR) + 0);
          bufD2[0].vtY = *(*(bufD + idR) + 1);
          bufD2[1].vtX = *(*(bufD + idR) + 2);
          bufD2[1].vtY = *(*(bufD + idR) + 3);
	  errNum = WlzGMModelConstructSimplex2D(model, bufD2);
	  break;
        case WLZ_GMMOD_3I: /* FALLTHROUGH */
        case WLZ_GMMOD_3D:
          bufD3[0].vtX = *(*(bufD + idR) + 0);
          bufD3[0].vtY = *(*(bufD + idR) + 1);
          bufD3[0].vtZ = *(*(bufD + idR) + 2);
          bufD3[1].vtX = *(*(bufD + idR) + 3);
          bufD3[1].vtY = *(*(bufD + idR) + 4);
          bufD3[1].vtZ = *(*(bufD + idR) + 5);
          bufD3[2].vtX = *(*(bufD + idR) + 6);
          bufD3[2].vtY = *(*(bufD + idR) + 7);
          bufD3[2].vtZ = *(*(bufD + idR) + 8);
	  errNum = WlzGMModelConstructSimplex3D(model, bufD3);
	  break;
      }
      ++idR;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to build model (%s).\n",
		     argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outObjFileStr, "-")?
	     fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Failed to open output file %s.\n",
      	      	     argv[0], outObjFileStr);
    }
  }
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_CONTOUR, ctrDom, dumVal, NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to create output woolz object (%s).\n",
		     argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, outObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to write output object (%s).\n",
      		     argv[0], errMsgStr);
    }
  }
  if(bufD)
  {
    (void )AlcDouble2Free(bufD);
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(ctrDom.ctr)
  {
    (void )WlzFreeContour(ctrDom.ctr);
  }
  if(fP && strcmp(outObjFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%sExample: %s%s",
      *argv,
      " [-o<output object>] [-h] [-2] [-3] [-d] [-i] [-o#] \n"
      "        [<input data>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -o  Output object file name.\n"
      "  -2  Build a 2D contour.\n"
      "  -3  Build a 3D contour.\n"
      "  -d  Build a contour with double floating point geometry.\n"
      "  -i  Build a contour with integer geometry.\n"
      "Builds a contour from the given input simplicies.\n"
      "The input data are read from stdin and the output object is written\n"
      "to stdout unless filenames are given.\n",
      *argv,
      " -d -2 -o tri.wlz\n"
      "0.0 0.0 1.0 0.0\n"
      "1.0 0.0 0.5 0.866025\n"
      "0.5 0.866025 0.0 0.0\n"
      "<eof>\n"
      "A 2D triangle is built from the three simplicies read from the\n"
      "standard input and the contour is written to tri.wlz.\n");
  }
  return(!ok);
}
