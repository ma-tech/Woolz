#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzGreyValue.c
* \author       Bill Hill
* \date         May 2004
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Woolz binary to determine a single grey value in an object.
* \ingroup
* \todo         -
* \bug          None known.
*/
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
		setVal = 0,
  		ok = 1,
		usage = 0;
  WlzDVertex3	pos;
  WlzPixelV	pix0,
  		pix1;
  WlzGreyValueWSpace *gVWSp = NULL;
  char		*outFileStr,
  		*inObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  const char	*errMsg;
  static char	optList[] = "o:v:x:y:z:h",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  pos.vtX = 0;
  pos.vtY = 0;
  pos.vtZ = 0;
  opterr = 0;
  pix0.type = WLZ_GREY_DOUBLE;
  pix0.v.dbv = 0.0;
  inObjFileStr = inObjFileStrDef;
  outFileStr = outFileStrDef;
  while(ok && (usage == 0) &&
  	((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
	break;
      case 'v':
	setVal = 1;
        if(sscanf(optarg, "%lg", &(pix0.v.dbv)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'x':
        if(sscanf(optarg, "%lg", &(pos.vtX)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'y':
        if(sscanf(optarg, "%lg", &(pos.vtY)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'z':
        if(sscanf(optarg, "%lg", &(pos.vtZ)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
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
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
                     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    gVWSp = WlzGreyValueMakeWSp(inObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to make workspace (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    WlzGreyValueGet(gVWSp, pos.vtZ, pos.vtY, pos.vtX);
    if(setVal)
    {
      errNum = WlzValueConvertPixel(&pix1, pix0, gVWSp->gType);
      if(errNum == WLZ_ERR_NONE)
      {
	switch(gVWSp->gType)
	{
	  case WLZ_GREY_LONG:
	    *(gVWSp->gPtr[0].lnp) = pix1.v.lnv;
	    break;
	  case WLZ_GREY_INT:
	    *(gVWSp->gPtr[0].inp) = pix1.v.inv;
	    break;
	  case WLZ_GREY_SHORT:
	    *(gVWSp->gPtr[0].shp) = pix1.v.shv;
	    break;
	  case WLZ_GREY_UBYTE:
	    *(gVWSp->gPtr[0].ubp) = pix1.v.ubv;
	    break;
	  case WLZ_GREY_FLOAT:
	    *(gVWSp->gPtr[0].flp) = pix1.v.flv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    *(gVWSp->gPtr[0].dbp) = pix1.v.dbv;
	    break;
	  case WLZ_GREY_RGBA:
	    *(gVWSp->gPtr[0].rgbp) = pix1.v.rgbv;
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
    }
    else
    {
      pix1.type = gVWSp->gType;
      pix1.v = gVWSp->gVal[0];
      errNum = WlzValueConvertPixel(&pix0, pix1, WLZ_GREY_DOUBLE);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(setVal == 0)
    {
      if(((fP = (strcmp(outFileStr, "-")?
		fopen(outFileStr, "w"):
		stdout)) == NULL) ||
	 (fprintf(fP, "%g\n", pix0.v.dbv) < 2))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to write output file (%s).\n",
		       *argv, errMsg);
      }
    }
    else
    {
      if(((fP = (strcmp(outFileStr, "-")?
		fopen(outFileStr, "w"):
		stdout)) == NULL) ||
	 ((errNum = WlzWriteObj(fP, inObj)) != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to write output object (%s).\n",
		       *argv, errMsg);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(gVWSp)
  {
    WlzGreyValueFreeWSp(gVWSp);
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: "
    "%s [-v#] [-x#] [-y#] [-z#] [-h] [<in object>]\n"
    "Either sets or gets a grey value within th given objects values.\n"
    "If the -v flag is given the objects value is set and the modified\n"
    "object is written to the output file, but if the flag is not set\n"
    "then the grey value is written out in ascii form with a new line\n"
    "character.\n"
    "Options are:\n"
    "  -v#   Grey value to be set (set to %g).\n"
    "  -x#   Vector column position (set to %g).\n"
    "  -y#   Vector line positio (set to %g).\n"
    "  -z#   Vector plane position (set to %g).\n"
    "  -o#   Output file name.\n"
    "  -h    Display this usage information.\n",
    *argv,
    pix0.v.dbv,
    pos.vtX, pos.vtY, pos.vtZ);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
