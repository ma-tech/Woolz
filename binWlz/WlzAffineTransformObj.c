#pragma ident "MRC HGU $Id$"
#define _WlzAffineTransformObj_c
/************************************************************************
* Project:      Woolz							*
* Title:        WlzAffineTransformObj.c			                *
* Date:         March 1999	                                    	*
* Author:       Bill Hill 				    		*
* Copyright:	1999 Medical Research Council, UK.			*
*		All rights reserved.					*
* Address:	MRC Human Genetics Unit,				*
*		Western General Hospital,				*
*		Edinburgh, EH4 2XU, UK.					*
* Purpose:      Applies an affine transform to a Woolz object.		*
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.	*
* 04-12-00 bill Add affine transform output.
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
		trInvert = 0,
		matrixValuesFlag = 0,
		noTransformationFlag = 0,
		primitivesFlag = 0,
		radiansFlag = 0,
                inverseTransformFlag = 0;
  double	trX = 0.0,
  		trY = 0.0,
		trZ = 0.0,
		trScale = 1.0,
  		trTheta = 0.0,
  		trPhi = 0.0,
		trAlpha = 0.0,
		trPsi = 0.0,
		trXsi = 0.0;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzTransformType trType = WLZ_TRANSFORM_2D_AFFINE;
  WlzAffineTransform *trans = NULL;
  WlzValues	trVal;
  WlzDomain	trDom;
  WlzObject	*inObj = NULL,
  		*outObj = NULL,
		*trObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inObjFileStr,
		*inTrObjFileStr = NULL,
		*outTrObjFileStr = NULL;
  const char	*trTypeStr;
  WlzAffineTransformPrim prim;
  static char	optList[] = "3LMNPRhiIo:a:b:s:t:T:u:v:w:x:y:z:",
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
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      case 'M':
        matrixValuesFlag = 1;
	break;
      case 'N':
        noTransformationFlag = 1;
	break;
      case 'P':
        primitivesFlag = 1;
	break;
      case 'R':
        radiansFlag = 1;
	break;
      case 'i':
        trInvert = 1;
	break;
      case 'I':
        inverseTransformFlag = 1;
	break;
      case '3':
	trType = WLZ_TRANSFORM_3D_AFFINE;
        break;
      case 'a':
	if(sscanf(optarg, "%lg", &trTheta) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'b':
	if(sscanf(optarg, "%lg", &trPhi) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 's':
	if(sscanf(optarg, "%lg", &trScale) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 't':
        inTrObjFileStr = optarg;
	break;
      case 'T':
        outTrObjFileStr = optarg;
	break;
      case 'u':
	if(sscanf(optarg, "%lg", &trAlpha) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'v':
	if(sscanf(optarg, "%lg", &trPsi) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'w':
	if(sscanf(optarg, "%lg", &trXsi) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'x':
	if(sscanf(optarg, "%lg", &trX) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'y':
	if(sscanf(optarg, "%lg", &trY) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'z':
	if(sscanf(optarg, "%lg", &trZ) != 1)
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
  if(ok && (radiansFlag == 0))
  {
    trTheta *= WLZ_M_PI / 180;
    trPhi *= WLZ_M_PI / 180;
    trPsi *= WLZ_M_PI / 180;
  }
  if(ok && (noTransformationFlag == 0))
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
  }

  /* get the transform - either from file or from input primitives */
  if(ok)
  {
    if(inTrObjFileStr)
    {
      if((*inTrObjFileStr == '\0') ||
	 ((fP = (strcmp(inTrObjFileStr, "-")?
		fopen(inTrObjFileStr, "r"): stdin)) == NULL) ||
	 ((trObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read object from file %s\n",
		       *argv, inObjFileStr);
      }
      if(fP && strcmp(inTrObjFileStr, "-"))
      {
	fclose(fP);
      }
      if(trObj)
      {
        if((trObj->type != WLZ_AFFINE_TRANS) ||
	   (trObj->domain.core == NULL))
	{
	  ok = 0;
	  (void )fprintf(stderr,
	  		 "%s: invalid transform object read from file %s\n",
			 *argv, inTrObjFileStr);
	}
	else
	{
	  trans = WlzAssignAffineTransform(trObj->domain.t, NULL);
	}
      }
    }
    else
    {
      trans = WlzAffineTransformFromPrimVal(trType, trX, trY, trZ,
					    trScale, trTheta, trPhi,
					    trAlpha, trPsi, trXsi,
					    trInvert, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )fprintf(stderr,
	            "%s: failed to make an affine transform from primitives\n",
		       *argv);
      }
    }
  }
  /* check for inverse transform */
  if( ok && inverseTransformFlag )
  {
    WlzAffineTransform	*tmpTrans;

    if((tmpTrans = WlzAffineTransformInverse(trans, &errNum)) == NULL)
    {
      ok = 1;
    }
    else 
    {
      (void) WlzFreeAffineTransform(trans);
      trans = tmpTrans;
    }
  }
  /* check for output of primitives or matrix */
  if(ok && primitivesFlag)
  {
    switch(WlzAffineTransformDimension(trans, NULL))
    {
      case 2:
      case 3:
	errNum = WlzAffineTransformPrimGet(trans, &prim);
	break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to get transform primitives\n",
		     *argv);

    }
    if(ok)
    {
      trTypeStr = WlzStringFromTransformType(trans->type, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to get transform type string\n",
		       *argv);
      }
    }
    if(ok)
    {
      (void )fprintf(stderr,
		     "type   %s\n"
		     "tx     %-10g\nty     %-10g\ntz     %-10g\nscale  %-10g\n"
		     "theta  %-10g\n"
		     "phi    %-10g\n"
		     "alpha  %-10g\npsi    %-10g\nxsi    %-10g\n"
		     "invert %-10d\n",
		     trTypeStr,
		     prim.tx, prim.ty, prim.tz, prim.scale,
		     radiansFlag? prim.theta: prim.theta * 180 / WLZ_M_PI,
		     radiansFlag? prim.phi: prim.phi * 180 / WLZ_M_PI,
		     prim.alpha, prim.psi, prim.xsi,
		     prim.invert);
    }
  }
  if(ok && matrixValuesFlag)
  {
    switch(WlzAffineTransformDimension(trans, NULL))
    {
      case 2:
	(void )fprintf(stderr,
		       "matrix %-10g %-10g %-10g\n"
		       "       %-10g %-10g %-10g\n"
		       "       %-10g %-10g %-10g\n",
		       trans->mat[0][0], trans->mat[0][1],
		       trans->mat[0][2],
		       trans->mat[1][0], trans->mat[1][1],
		       trans->mat[1][2],
		       trans->mat[2][0], trans->mat[2][1],
		       trans->mat[2][2]);
        break;
      case 3:
	(void )fprintf(stderr,
		       "matrix %-10g %-10g %-10g %-10g\n"
		       "       %-10g %-10g %-10g %-10g\n"
		       "       %-10g %-10g %-10g %-10g\n"
		       "       %-10g %-10g %-10g %-10g\n",
		       trans->mat[0][0], trans->mat[0][1],
		       trans->mat[0][2], trans->mat[0][3],
		       trans->mat[1][0], trans->mat[1][1],
		       trans->mat[1][2], trans->mat[1][3],
		       trans->mat[2][0], trans->mat[2][1],
		       trans->mat[2][2], trans->mat[2][3],
		       trans->mat[3][0], trans->mat[3][1],
		       trans->mat[3][2], trans->mat[3][3]);
        break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: invalid transform type\n",
		     *argv);

    }
  }
  if(ok && (outTrObjFileStr != NULL))
  {
    if(trObj == NULL)
    {
      trDom.t = trans;
      trVal.core = NULL;
      trObj = WlzMakeMain(WLZ_AFFINE_TRANS, trDom, trVal, NULL, NULL,
      			  &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to make transform object.\n",
		       *argv);
      }
    }
    if(ok)
    {
      if(((fP = (strcmp(outTrObjFileStr, "-")?
		fopen(outTrObjFileStr, "w"):
		stdout)) == NULL) ||
	  (WlzWriteObj(fP, trObj) != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write output transform object\n",
		       *argv);
      }
      if(fP && strcmp(outTrObjFileStr, "-"))
      {
	fclose(fP);
      }
    }
  }
  if(ok && (noTransformationFlag == 0))
  {
    outObj = WlzAssignObject(WlzAffineTransformObj(inObj, trans, interp,
    						   &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to transform object\n",
		     *argv);
    }
  }
  if(ok && (noTransformationFlag == 0))
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
  if(trObj)
  {
    WlzFreeObj(trObj);
  }
  if(trans)
  {
    WlzFreeAffineTransform(trans);
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
    " [-o<output object>] [-L] [-M] [-N] [-P] [-R]\n"
    "        [-h] [-i] [-3] [-x#] [-y#] [-z#] [-s#]"
    "[-t#] [-T#] [-a#] [-b#] [-u#] [-v#] [-w#]\n"
    "        [<input object>]\n" 
    "Options:\n"
    "  -3  3D transform instead of 2D.\n"
    "  -L  Use linear interpolation instead of nearest neighbour.\n"
    "  -M  Print matix values.\n"
    "  -N  No transformation.\n"
    "  -P  Print transform primatives.\n"
    "  -R  Use radians for angles instead of degrees.\n"
    "  -a  Rotation about the z-axis.\n"
    "  -b  Rotation about the y-axis.\n"
    "  -h  Help, prints this usage message.\n"
    "  -i  Invert: reflect about the y-axis.\n"
    "  -I  Inverse: use the inverse of the input transform.\n"
    "  -o  Output object file name.\n"
    "  -s  Scale factor.\n" 
    "  -t  Input affine transform object.\n" 
    "  -T  Output affine transform object.\n" 
    "  -u  Shear strength.\n"
    "  -v  Shear angle in x-y plane.\n"
    "  -w  3D shear angle.\n"
    "  -x  Column (x) translation.\n"
    "  -y  Row (y) translation.\n"
    "  -z  Plane (z) translation.\n"
    "Applies an affine transform to a Woolz object.\n"
    "A composite transform is applied to the object with the order of\n"
    "composition being scale (applied first), shear, rotation and then\n"
    "translation (applied last).\n"
    "If a transform object is specified on the command line then none\n"
    "of the command line transform primatives are used.\n"
    "The input object is read from stdin and the transformed object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -x100 -y200 -o shifted.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, shifted 100 columns\n"
    "and 200 lines and then written to shifted.wlz\n");
  }
  return(!ok);
}
