#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzRsvFilterObj.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Canny/Deriche edge detection filter for Woolz domain
*		objects with grey values.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

typedef enum
{
  WLZ_CANNYDERICHE_FLTACT_DIR,
  WLZ_CANNYDERICHE_FLTACT_EDG,
  WLZ_CANNYDERICHE_FLTACT_GRD
} WlzCannyDericheFilterAction;

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0;
  WlzObject	*inObj = NULL,
  		*outObj = NULL,
		*grdObj = NULL,
		*edgObj = NULL;
  double	alpha = 1.0,
  		mult = 1.0;
  WlzPixelV	pThrV,
  		sThrV;
  WlzCannyDericheFilterAction action = WLZ_CANNYDERICHE_FLTACT_EDG;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inObjFileStr;
  static char	optList[] = "ho:a:m:l:u:deg",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  pThrV.type = sThrV.type = WLZ_GREY_INT;
  pThrV.v.inv = sThrV.v.inv = 1.0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'a':
	if((optarg == NULL) || (sscanf(optarg, "%lg", &alpha) != 1))
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'd':
        action = WLZ_CANNYDERICHE_FLTACT_DIR;
	break;
      case 'e':
        action = WLZ_CANNYDERICHE_FLTACT_EDG;
	break;
      case 'g':
        action = WLZ_CANNYDERICHE_FLTACT_GRD;
        break;
      case 'm':
        if((optarg == NULL) || (sscanf(optarg, "%lg", &mult) != 1))
	{
	  usage = 1;
	  ok = 0;
	}
      case 'l':
        if((optarg == NULL) || (sscanf(optarg, "%d", &(sThrV.v.inv)) != 1))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'u':
        if((optarg == NULL) || (sscanf(optarg, "%d", &(pThrV.v.inv)) != 1))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
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
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    switch(action)
    {
      case WLZ_CANNYDERICHE_FLTACT_DIR:
	outObj = WlzCannyDeriche(NULL, inObj, alpha, mult,
				 pThrV, sThrV, &errNum);
        break;
      case WLZ_CANNYDERICHE_FLTACT_EDG:
	edgObj = WlzCannyDeriche(&grdObj, inObj, alpha, mult,
				 pThrV, sThrV, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzMakeMain(edgObj->type, edgObj->domain, grdObj->values,
			       NULL, NULL, &errNum);
	}
	WlzFreeObj(edgObj);
	WlzFreeObj(grdObj);
        break;
      case WLZ_CANNYDERICHE_FLTACT_GRD:
	edgObj = WlzCannyDeriche(&outObj, inObj, alpha, mult,
				 pThrV, sThrV, &errNum);
	WlzFreeObj(edgObj);
        break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
       (void )fprintf(stderr, "%s Failed to filter object (%s)\n",
       	              *argv, WlzStringFromErrorNum(errNum, NULL));
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
	      fopen(outFileStr, "w"):
	      stdout)) == NULL) ||
       (WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to write output object\n",
		     *argv);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(outObj)
  {
    WlzFreeObj(outObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<output object>] [-h] [-o] [-a#] [-d] [-e] [-g]\n"
    "        [-m#] [-l#] [-u#] [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information\n"
    "  -o  Output object file name.\n"
    "  -a  Deriche's alpha parameter.\n"
    "  -d  Output object with edge direction encoding.\n"
    "  -e  Output object with edge gradients.\n"
    "  -g  Output object with all gradients.\n"
    "  -m  Filter multiplier parameter.\n"
    "  -l  Lower hysteresis threshold value.\n"
    "  -u  Upper hysteresis threshold value.\n"
    "Applies a Canny/Deriche edge detection filter to a Woolz domain object\n"
    "with grey values.\n"
    "The input object is read from stdin and the filtered object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -e -a 1.0 -l 15 -u 40 -m 10 -o edges.wlz in.wlz\n"
    "The input Woolz object is read from in.wlz, and filtered using a\n"
    "Deriche edge operator with alpha=1.0, gradients are multiplied by\n"
    "10.0 before non-maximal suppression or hysteresis thresholding. The\n"
    "hysteresis threshold values of 15 and 40 are used to threshold the\n"
    "edge image which is written to edge.wlz.\n");
  }
  return(!ok);
}
