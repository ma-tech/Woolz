#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFContourObj.c
* \author       Bill Hill
* \date         September 2001
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	A program to extract boundary, maximal gradient and isosurface
* 		contours from 2 or 3D Woolz and saving them in VTK format
*		files.
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int           option,
  		ok = 1,
		usage = 0;
  double	ctrVal = 100,
  		ctrWth = 1.0;
  FILE		*fP = NULL;
  char		*fStr,
		*inObjFileStr,
  		*outFileStr;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  WlzDomain	ctrDom;
  WlzValues	dumVal;
  WlzContourMethod ctrMtd = WLZ_CONTOUR_MTD_ISO;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "bghio:v:w:";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  ctrDom.core = NULL;
  dumVal.core = NULL;
  outFileStr = (char *)outFileStrDef;
  inObjFileStr = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        ctrMtd = WLZ_CONTOUR_MTD_BND;
	break;
      case 'g':
        ctrMtd = WLZ_CONTOUR_MTD_GRD;
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'i':
        ctrMtd = WLZ_CONTOUR_MTD_ISO;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'v':
        if(sscanf(optarg, "%lg", &ctrVal) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'w':
        if(sscanf(optarg, "%lg", &ctrWth) != 1)
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
    ctrDom.ctr = WlzContourObj(inObj, ctrMtd, ctrVal, ctrWth, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to compute contour (%s).\n",
      	     	     argv[0], errMsgStr);
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
    if(strcmp(outFileStr, "-") == 0)
    {
      fP = stdout;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = outFileStr;
    }
  }
  if(ok)
  {
    errNum = WlzEffWriteObj(fP, fStr, outObj, WLZEFF_FORMAT_VTK);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to write output object (%s).\n",
      		     argv[0], errMsgStr);
    }
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(outObj)
  {
    (void )WlzFreeObj(outObj);
  }
  else if(ctrDom.ctr)
  {
    (void )WlzFreeContour(ctrDom.ctr);
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%sExample: %s%s",
      *argv,
      " [-o<output object>] [-h] [-o] [-g] [-i] [-o#] [-v#] [-w#]\n"
      "        [<input object>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -o  Output object file name.\n"
      "  -b  Compute object boundary contours.\n"
      "  -g  Compute maximal gradient contours.\n"
      "  -i  Compute iso-value contours.\n"
      "  -v  Contour iso-value or minimum gradient.\n"
      "  -w  Contour (Deriche) gradient operator width.\n"
      "Computes a contour object from the given Woolz object and saves it\n"
      "using the VTK ascii polydata format.\n"
      "The input object is read from stdin and output data are written\n"
      "to stdout unless filenames are given.\n",
      *argv,
      " -o out.vtk -i -v 100.0 in.wlz\n"
      "The input Woolz object is read from in.wlz, and the iso-value\n"
      "(iso-value = 100.0) contour is written in VTK ascii polydata format\n"
      "to out.vtk.\n");
  }
  return(!ok);
}
