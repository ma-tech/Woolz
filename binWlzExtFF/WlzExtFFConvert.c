#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFFConvert.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz filter which uses the external file formats
*		extension to convert object data file formats.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Wlz.h>
#include <WlzExtFF.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

static WlzEffFormat WlzEffStringExtFromFileName(const char *fNameStr)
{
  WlzEffFormat	fFmt = WLZEFF_FORMAT_NONE;
  int		fNameLen;

  if(fNameStr && ((fNameLen = strlen(fNameStr)) > 4) &&
     (*(fNameStr + fNameLen - 4) == '.') &&
     (*(fNameStr + fNameLen - 1) != '/'))
  {
    fFmt = WlzEffStringExtToFormat(fNameStr + fNameLen - 3);
  }
  return(fFmt);
}

int             main(int argc, char **argv)
{
  int		bgdFlag = 0,
  		option,
		ok = 1,
		usage = 0;
  WlzEffFormat  inFmt = WLZEFF_FORMAT_NONE,
  		outFmt = WLZEFF_FORMAT_NONE;
  WlzFVertex3	voxelSz;
  WlzIVertex3	voxelSzFlags;
  FILE		*fP;
  WlzObject     *inObj = NULL;
  char		*fStr,
  		*outObjFileStr,
  		*inObjFileStr;
  WlzPixelV	tV0,
  		bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  static char	optList[] = "b:f:F:o:x:y:z:h",
		outObjFileStrDef[] = "-",
		inObjFileStrDef[] = "-";
 
  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  voxelSzFlags.vtX = voxelSzFlags.vtY = voxelSzFlags.vtZ = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        if(sscanf(optarg, "%lg", &(bgdV.v.dbv)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	else
	{
	  bgdV.type = WLZ_GREY_DOUBLE;
	  bgdFlag = 1;
	}
        break;
      case 'f':
        if((inFmt = WlzEffStringExtToFormat(optarg)) == NULL)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'F':
        if((outFmt = WlzEffStringExtToFormat(optarg)) == NULL)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'x':
	voxelSzFlags.vtX = 1;
	if(sscanf(optarg, "%g", &(voxelSz.vtX)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'y':
	voxelSzFlags.vtY = 1;
	if(sscanf(optarg, "%g", &(voxelSz.vtY)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'z':
	voxelSzFlags.vtZ = 1;
	if(sscanf(optarg, "%g", &(voxelSz.vtZ)) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
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
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    /* Try and determine file formats from extensions if not known already. */
    if(inFmt == WLZEFF_FORMAT_NONE)
    {
      if((inFmt = WlzEffStringExtFromFileName(inObjFileStr)) == NULL)
      {
        usage = 1;
	ok = 0;
      }
    }
    if(outFmt == WLZEFF_FORMAT_NONE)
    {
      if((outFmt = WlzEffStringExtFromFileName(outObjFileStr)) == NULL)
      {
        usage = 1;
	ok = 0;
      }
    }
  }
  if(ok)
  {
    if(strcmp(inObjFileStr, "-") == 0)
    {
      fP = stdin;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = inObjFileStr;
    }
    errNum = WLZ_ERR_READ_EOF;
    if((inObj= WlzAssignObject(WlzEffReadObj(fP, fStr,
       					     inFmt, &errNum), NULL)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: failed to read object from file(s) %s (%s)\n",
                     *argv, inObjFileStr, errMsg);
    }
  }
  if(ok && bgdFlag)
  {
    /* Set background if requested on command line */
    tV0 = WlzGetBackground(inObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzValueConvertPixel(&bgdV, bgdV, tV0.type);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum =  WlzSetBackground(inObj, bgdV);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
               "%s: failed to set background of object from file(s) %s (%s)\n",
		     *argv, inObjFileStr, errMsg);
    }
  }
  if(ok && (inObj->type == WLZ_3D_DOMAINOBJ) && inObj->domain.core &&
     (inObj->domain.core->type == WLZ_PLANEDOMAIN_DOMAIN))
  {
    if(voxelSzFlags.vtX)
    {
      inObj->domain.p->voxel_size[0] = voxelSz.vtX;
    }
    if(voxelSzFlags.vtY)
    {
      inObj->domain.p->voxel_size[1] = voxelSz.vtY;
    }
    if(voxelSzFlags.vtZ)
    {
      inObj->domain.p->voxel_size[2] = voxelSz.vtZ;
    }
  }
  if(ok)
  {
    if(strcmp(outObjFileStr, "-") == 0)
    {
      fP = stdout;
      fStr = NULL;
    }
    else
    {
      fP = NULL;
      fStr = outObjFileStr;
    }
    if((errNum = WlzEffWriteObj(fP, fStr,
    				inObj, outFmt)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: failed to write object to file(s) %s (%s).\n",
                     *argv, outObjFileStr, errMsg);
    }
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(!ok)
  {
    if(usage)
    {
      (void )fprintf(stderr,
      "Usage: %s%s%s%s%s%s%s%s",
      *argv,
      " [-h] [-b<background>]\n"
      "        [-f<input format>] [-F<output format>]\n"
      "        [-x<x size>] [-y<y size>] [-z<z size>]\n"
      "        [-o<output file>] [<input file>)]\n"
      "Converts image objects between one file format and another,\n"
      "neither of which need be the Woolz data file format.\n"
      "Options are:\n"
      "  -h    Help, prints this usage information.\n"
      "  -b#   Set background to vale,\n"
      "  -f#   Input file format.\n"
      "  -F#   Ouput file format.\n"
      "  -o#   Output file name.\n"
      "  -x#   X voxel/pixel size.\n"
      "  -y#   Y voxel/pixel size.\n"
      "  -z#   Z voxel/pixel size.\n"
      "Valid formats are:\n"
      "  Description                 Fmt     Extension\n"
      "  ***********                 ***     *********\n"
      "  Microsoft Bitmap            bmp     .bmp\n"
      "  Stanford Density            den     .den\n"
      "  ICS                         ics     .ics/.ids\n"
      "  PNM                         pnm     .pgm\n"
      "  BioRad Confocal             pic     .pic\n"
      "  SLC                         slc     .slc\n"
      "  Sunvision VFF               vff     .vff\n"
      "  Visualization Toolkit VTK   vtk     .vtk\n"
      "  MRC HGU Woolz               wlz     .wlz\n",
      "Simple example:\n  ",
      *argv,
      " -f wlz -F slc <in.wlz >out.slc\n"
      "  Converts the Woolz object in.wlz to an SLC data file out.slx\n",
      "More complex example:\n  ",
      *argv,
      " -f den -F pnm -o out.pgm in.den\n"
      "  Converts the Stanford density file in.den to a series of PGM files\n"
      "  each with a name of the form out000001.pgm where the number\n"
      "  encodes the image plane and a control file which specifies the\n"
      "  volume origin, size and voxel dimensions.\n"
      "By default objects are read from the standard input and written to\n"
      "the standard output.\n"
      "File formats which use more than one file can not be read or written\n"
      "using the standard input or standard output.\n");
    }
  }
  return(!ok);
}
