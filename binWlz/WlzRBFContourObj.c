#ifndef DOXYGEN_SHOULD_SKIP_THIS
#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzRBFContourObj.c
* \author       Bill Hill
* \date         April 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Computes a contour from a Woolz domain by using radial
*		basis function to approximate the signed distance to
*		the boundary of the domain and then extracting the zero
*		level set through the volume as a contour.
* \todo         -
* \bug          None known.
*/

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
  int           option,
		bErosion = 4,
		bDilation = 4,
		sDilation = 10,
  		sFac = 10,
  		oFac = 10,
  		ok = 1,
		usage = 0,
		unitVoxelSz = 0;
  double	sAlpha = 0.1,
		oAlpha = 2.0,
		delta = 0.1,
		tau = 0.1,
		samFac = 1.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  WlzDomain	ctrDom;
  WlzValues	dumVal;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "o:s:S:Z:f:F:a:A:d:t:M:hU";
  const char	outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  ctrDom.core = NULL;
  dumVal.core = NULL;
  outFileStr = (char *)outFileStrDef;
  inObjFileStr = (char *)inObjFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'h':
        usage = 1;
	break;
      case 'U':
        unitVoxelSz = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 's':
        usage = sscanf(optarg, "%d", &bErosion) != 1;
	break;
      case 'S':
        usage = sscanf(optarg, "%d", &bDilation) != 1;
	break;
      case 'Z':
        usage = sscanf(optarg, "%d", &sDilation) != 1;
	break;
      case 'f':
        usage = sscanf(optarg, "%d", &sFac) != 1;
	break;
      case 'F':
        usage = sscanf(optarg, "%d", &oFac) != 1;
	break;
      case 'a':
        usage = sscanf(optarg, "%lg", &sAlpha) != 1;
	break;
      case 'A':
        usage = sscanf(optarg, "%lg", &oAlpha) != 1;
	break;
      case 'd':
        usage = sscanf(optarg, "%lg", &delta) != 1;
	break;
      case 't':
        usage = sscanf(optarg, "%lg", &tau) != 1;
	break;
      case 'M':
        usage = sscanf(optarg, "%lg", &samFac) != 1;
	break;
      default:
        usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  ok = !usage;
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
  if(ok && unitVoxelSz)
  {
    if((inObj->type == WLZ_3D_DOMAINOBJ) &&
       (inObj->domain.core) &&
       (inObj->domain.core->type == WLZ_PLANEDOMAIN_DOMAIN))
    {
      inObj->domain.p->voxel_size[0] = 1.0;
      inObj->domain.p->voxel_size[1] = 1.0;
      inObj->domain.p->voxel_size[2] = 1.0;
    }
  }
  if(ok)
  {
    ctrDom.ctr = WlzContourRBFBndObj3D(inObj, bErosion, bDilation, sDilation,
    				       sFac, oFac, sAlpha, oAlpha,
				       delta, tau, samFac, &errNum);
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
  if((fP = (strcmp(outFileStr, "-")?
	   fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Failed to open output file %s.\n",
      	      	     argv[0], outFileStr);
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
    " [-o<output object>] [-h] [-o] [-s#] [-S#] [-Z#] [-f#] [-F#]\n"
    "        [-a#] [-A#] [-d#] [-t#] [-M#] [-U] [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output object file name.\n"
    "  -s  Distance for points inside domain surface.\n"
    "  -S  Distance for points outside domain surface.\n"
    "  -Z  Distance from domain surface within which to evaluate the\n"
    "      radial basis function.\n"
    "  -f  Sampling factor for surface points.\n"
    "  -F  Sampling factor for interior and extrior points.\n"
    "  -a  Surface alpha value, must be greater than zero.\n"
    "  -A  Interior and extrior alpha value, must be greater than zero.\n"
    "  -d  Multiorder spline delta value.\n"
    "  -t  Multiorder spline tau value.\n"
    "  -M  Distance object sampling factor (can be used to reduce the\n"
    "      number of surface facets.)\n"
    "  -U  Use unit voxel size.\n"
    "Computes a contour from a Woolz domain by using radial basis\n"
    "function to approximate the signed distance to the boundary of\n"
    "the domain and then extracting the zero level set through the\n"
    "volume as a contour.\n"
    "The input object is read from stdin and output data are written\n"
    "to stdout unless filenames are given.\n",
    *argv,
    " -o out.wlz -Z20 in.wlz\n"
    "The input Woolz object is read from in.wlz; a signed distance\n"
    "distance function is evaluated within the volume formed by dilating\n"
    "the surface of the input domain by 20 and the contour is written\n"
    "the file out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
