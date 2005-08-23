#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzRsvFilterObj.c
* Date:         July 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes a gradient image.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

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
		printFtrPrm = 0;
  unsigned int	actMsk = WLZ_RSVFILTER_ACTION_X | WLZ_RSVFILTER_ACTION_Y |
  			 WLZ_RSVFILTER_ACTION_Z;
  double	ftrMlt = 1.0,
  		ftrPrm = 1.0;
  WlzRsvFilterName ftrName = WLZ_RSVFILTER_NAME_DERICHE_1;
  WlzRsvFilter *ftr = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inObjFileStr;
  static char	optList[] = "ho:gdm:p:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
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
      case 'g':
        ftrName = WLZ_RSVFILTER_NAME_GAUSS_1;
	break;
      case 'd':
        ftrName = WLZ_RSVFILTER_NAME_DERICHE_1;
	break;
      case 'm':
	if(sscanf(optarg, "%lg", &ftrMlt) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'p':
	if(sscanf(optarg, "%lg", &ftrPrm) != 1)
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
    if(ok)
    {
      if(((ftr = WlzRsvFilterMakeFilter(ftrName, ftrPrm, &errNum)) == NULL) ||
         (errNum != WLZ_ERR_NONE))
      {
        ok = 0;
	(void )fprintf(stderr, "%s Failed to create filter (%s)\n",
		       *argv, WlzStringFromErrorNum(errNum, NULL));
      }
      else
      {
	ftr->c *= ftrMlt;
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
      if(ok)
      {
        if(((outObj = WlzAssignObject(
		      WlzGreyGradient(NULL, NULL, NULL,
		      		      inObj, ftr, &errNum), NULL)) == NULL) ||
	   (errNum != WLZ_ERR_NONE))

        {
	  ok = 0;
	  (void )fprintf(stderr, "%s Failed to compute gradient object (%s)\n",
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
  if(ftr)
  {
    WlzRsvFilterFreeFilter(ftr);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<output object>] [-h] [-o] [-g] [-d] [-m#] [-p#]\n"
    "        [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information\n"
    "  -o  Output object file name.\n"
    "  -g  Use approximate Gaussian filter.\n"
    "  -d  Use Deriche filter.\n"
    "  -p  Filter parameter for Gaussian (sigma) and Deriche (alpha) filters\n"
    "  -m  Extra multiplication factor for gradient grey levels.\n"
    "Computes a grey gradient object from the input object.\n"
    "The input object is read from stdin and the gradient object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -d -m 2.0 -o grad.wlz in.wlz\n"
    "The input Woolz object is read from in.wlz, and the grey gradient\n"
    "image is computed using a Deriche filter and an additional multiplying\n"
    "factor of 2.0. The gradient object is written to grad.wlz.\n"
    "Objects with unsigned byte grey values will be promoted to short\n"
    "grey valued objects.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
