#pragma ident "MRC HGU $Id$"
/*!
 * \file	WlzMatchICPPlane.c
 * \author      Bill Hill
 * \date        June 2002
 * \version	$Id$
 * \note
 *             	Copyright:
 *             	2001 Medical Research Council, UK.
 *             	All rights reserved.
 * \par  Address:
 *             	MRC Human Genetics Unit,
 *             	Western General Hospital,
 *             	Edinburgh, EH4 2XU, UK.
 * \brief      	Compute tie points for a computed plane of section in a
 *		reference object and a 2D object, where these are read from
 *		an 'MAPaint 2D warp input parameters' bibfile.
 * \todo
 * \bug
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <Wlz.h>
#include <bibFile.h>
#include <WlzExtFF.h>
#include <string.h>

/* TODO remove this #define makeing this the only code available. */
#define WLZ_MATCHICPPLANE_OWNCB

static FILE			*WlzMatchICPPlaneSecParFile(
				  FILE *lFP,
				  int multi,
				  int index,
				  char *secParFile);
static WlzErrorNum 		WlzMatchICPPlaneReadSecParam(
				  FILE *fP,
				  WlzThreeDViewStruct **dstView,
				  char **dstRefFileStr,
				  WlzEffFormat *dstRefFileType,
				  WlzObject **dstRefObj,
				  char **dstSrcFileStr,
				  WlzEffFormat *dstSrcFileType,
				  WlzObject **dstSrcObj);
static WlzErrorNum		WlzMatchICPPlaneWriteSecParam(
				  FILE *fP,
			          WlzThreeDViewStruct *view,
				  char *refObjFileStr,
				  WlzEffFormat refObjFileType,
				  char *srcObjFileStr,
				  WlzEffFormat srcObjFileType,
			      	  WlzAffineTransform *srcTr,
				  int nMatch,
				  WlzDVertex2 *tieRP,
				  WlzDVertex2 *tieSP,
				  WlzDVertex2 *refObj2DOrg);
static WlzObject 		*WlzMatchICPPlaneCreateContourObj(
				  WlzObject *gObj,
				  int medianSz,
				  double smooth,
				  WlzAffineTransform *tr,
				  double cThr,
				  int minSpx,
				  int debug,
				  char *objDbgFileName,
				  WlzErrorNum *dstErr);

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int		option,
  		ok = 1,
		debug = 0,
		usage = 0,
		verbose = 0,
		verboseObj = 0,
		index,
  		nMatch = 0,
		useCOfM = 0,
		maxItr = 200,
		minSpx = 15,
		minSegSpx = 10,
		matchImpNN = 7,
		noMatching = 0,
		decompLimit = INT_MAX,
		removeRefOrg = 1,
		refMedianSz = 0,
		srcMedianSz = 0,
		multipleFiles = 0;
  double	maxAng = 30 * (ALG_M_PI / 180.0),
  		maxDeform = 0.5,
  		maxDisp = 20.0,
		refCThr = 20.0,
		srcCThr = 20.0,
		matchImpThr = 2.5,
		refSmooth = 3.0,
		srcSmooth = 3.0;
  double	parseDbl[2];
  FILE		*vFP,
  		*lFP = NULL,
  		*fP = NULL;
  char		*inFileStr = NULL,
		*inTrFileStr = NULL,
  		*outFileBaseStr = NULL,
		*ctrFileBaseStr = NULL;
  char		secParFile[256];
  char		*parseStr[2];
  char		*refObjFileStr = NULL,
  		*srcObjFileStr = NULL;
  WlzEffFormat	refObjFileType = WLZEFF_FORMAT_WLZ,
  		srcObjFileType = WLZEFF_FORMAT_WLZ;
  WlzObject	*tObj0 = NULL,
		*tObj1 = NULL,
  		*refObj3D = NULL,
		*refObj2D = NULL,
  		*srcObj2D = NULL,
		*refCObj2D = NULL,
  		*srcCObj2D = NULL;
  WlzDVertex2	tDV0,
  		refCOfM,
  		srcCOfM,
		refObj2DOrg;
  WlzAffineTransform *tTr0,
		*cOfMTr = NULL,
  		*inTr = NULL,
  		*invTr = NULL;
  WlzThreeDViewStruct *view = NULL;
  WlzVertexP	matchRP,
  		matchSP;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzMatchICPWeightCbData cbData;
  WlzAffineTransformPrim inTrPrim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		fileNameBuf[FILENAME_MAX];
  const char	*errMsg;
  static char	optList[] = "dhvVo:r:Yt:x:y:a:s:efg:k:u:i:p:A:S:F:P:m:n:c:NL";
  const int	nScatter = 5;
  const char	nullStr[] = "<NULL>",
  		inFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  matchRP.v = matchSP.v = NULL;
  inFileStr = (char *)inFileStrDef;
  outFileBaseStr = (char *)outFileStrDef;
  (void )memset(&inTrPrim, 0, sizeof(WlzAffineTransformPrim));
  while((usage == 0) && ok &&
        ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
        debug = 1;
	break;
      case 'h':
        usage = 1;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'V':
        verboseObj = 1;
	break;
      case 'o':
        outFileBaseStr = optarg;
	break;
      case 'r':
        if((refObjFileStr = AlcStrDup(optarg)) == NULL)
	{
	  ok = 0;
	  (void )fprintf(stderr, "%s Failed to allocate enough memory.\n",
	     	         *argv);
	}
        break;
      case 'Y':
        multipleFiles = 1;
	break;
      case 't':
        inTrFileStr = optarg;
	break;
      case 'x':
        if(sscanf(optarg, "%lg", &(inTrPrim.tx)) != 1)
	{
	  usage = 1;
	}
	break;
      case 'y':
        if(sscanf(optarg, "%lg", &(inTrPrim.ty)) != 1)
	{
	  usage = 1;
	}
	break;
      case 'a':
        if(sscanf(optarg, "%lg", &(inTrPrim.theta)) != 1)
	{
	  usage = 1;
	}
	else
	{
	  inTrPrim.theta *= ALG_M_PI / 180.0;
	}
	break;
      case 's':
        if(sscanf(optarg, "%lg", &(inTrPrim.scale)) != 1)
	{
	  usage = 1;
	}
	break;
      case 'e':
        useCOfM = 1;
	break;
      case 'f':
        removeRefOrg = 0;
	break;
      case 'g': /* FALLTHROUGH */
      case 'k': /* FALLTHROUGH */
      case 'u':
        /* Parse <ref value>,<src value> with both optional. */
	if(optarg)
	{
          while(*optarg && isspace(*optarg))
          {
            ++optarg;
          }
          if(*optarg == ',')
          {
            parseStr[0] = NULL;
            parseStr[1] = strtok(optarg, ",");
          }
          else
          {
            parseStr[0] = strtok(optarg, ",");
            parseStr[1] = strtok(NULL, ",");
          }
          if((parseStr[0] == NULL) && (parseStr[1] == NULL))
          {
            usage = 1;
          }
          else
          {
	    if((parseStr[0] &&
	        (sscanf(parseStr[0], "%lg", parseDbl + 0) != 1)) ||
	       (parseStr[1] &&
	        (sscanf(parseStr[1], "%lg", parseDbl + 1) != 1)))
	    {
	      usage = 1;
	    }
          }
	  if(usage == 0)
	  {
	    switch(option)
	    {
	      case 'g':
	        refCThr = parseDbl[0];
	        srcCThr = parseDbl[1];
		break;
	      case 'k':
	        refMedianSz = WLZ_NINT(parseDbl[0]);
	        srcMedianSz = WLZ_NINT(parseDbl[1]);
		break;
	      case 'u':
	        refSmooth = parseDbl[0];
	        srcSmooth = parseDbl[1];
		break;
	    }
	  }
        }
        break;
      case 'i':
        if((sscanf(optarg, "%d", &maxItr) != 1) || (maxItr <= 0))
	{
	  usage = 1;
	}
	break;
      case 'p':
        if((sscanf(optarg, "%d", &minSpx) != 1) || (minSpx <= 10))
	{
	  usage = 1;
	}
	break;
      case 'P':
        if(sscanf(optarg, "%d", &minSegSpx) != 1)
	{
	  usage = 1;
	}
	break;
      case 'A':
        if((sscanf(optarg, "%lg", &maxAng) != 1) ||
	   (maxAng < 0.0) || (maxAng > 180.0))
	{
	  usage = 1;
	}
	else
	{
	  maxAng *= ALG_M_PI / 180;
	}
	break;
      case 'F':
        if((sscanf(optarg, "%lg", &maxDeform) != 1) || (maxDeform < 0.0))
	{
	  usage = 1;
	}
	break;
      case 'S':
        if((sscanf(optarg, "%lg", &maxDisp) != 1) || (maxDisp <= 0.0))
	{
	  usage = 1;
	}
	break;
      case 'm':
        if((sscanf(optarg, "%lg", &matchImpThr) != 1) || (matchImpThr < 0.0))
	{
	  usage = 1;
	}
	break;
      case 'n':
        if((sscanf(optarg, "%d", &matchImpNN) != 1) || (matchImpNN < 1))
	{
	  usage = 1;
	}
	break;
      case 'c':
        ctrFileBaseStr = optarg;
        if(strlen(ctrFileBaseStr) > (FILENAME_MAX - 16))
	{
	  usage = 1;
	}
	break;
      case 'N':
        noMatching = 1;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      default:
        usage = 1;
	break;
    }
  }
  if(usage)
  {
    ok = 0;
  }
  if(ok)
  {
    if((inFileStr == NULL) || (*inFileStr == '\0') ||
       (outFileBaseStr == NULL) || (*outFileBaseStr == '\0'))
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
  if(verbose)
  {
    (void )fprintf(stderr, "Parameter and other internal variable values:\n");
    (void )fprintf(stderr, "  ok = %d\n", ok);
    (void )fprintf(stderr, "  debug = %d\n", debug);
    (void )fprintf(stderr, "  verbose = %d\n", verbose);
    (void )fprintf(stderr, "  verboseObj = %d\n", verboseObj);
    (void )fprintf(stderr, "  outFileBaseStr = %s\n",
		   outFileBaseStr? outFileBaseStr: nullStr);
    (void )fprintf(stderr, "  refObjFileStr = %s\n",
		   refObjFileStr? refObjFileStr: nullStr);
    (void )fprintf(stderr, "  multipleFiles = %d\n", multipleFiles);
    (void )fprintf(stderr, "  inTrFileStr = %s\n",
                   inTrFileStr? inTrFileStr: nullStr);
    (void )fprintf(stderr, "  inTrPrim.tx = %g\n", inTrPrim.tx);
    (void )fprintf(stderr, "  inTrPrim.ty = %g\n", inTrPrim.ty);
    (void )fprintf(stderr, "  inTrPrim.theta = %g (Radians)\n",
    		   inTrPrim.theta);
    (void )fprintf(stderr, "  inTrPrim.scale = %g\n", inTrPrim.scale);
    (void )fprintf(stderr, "  useCOfM = %d\n", useCOfM);
    (void )fprintf(stderr, "  removeRefOrg = %d\n", removeRefOrg);
    (void )fprintf(stderr, "  refCThr = %g\n", refCThr);
    (void )fprintf(stderr, "  srcCThr = %g\n", srcCThr);
    (void )fprintf(stderr, "  refMedianSz = %g\n", refMedianSz);
    (void )fprintf(stderr, "  srcMedianSz = %g\n", srcMedianSz);
    (void )fprintf(stderr, "  refSmooth = %g\n", refSmooth);
    (void )fprintf(stderr, "  srcSmooth = %g\n", srcSmooth);
    (void )fprintf(stderr, "  maxItr = %d\n", maxItr);
    (void )fprintf(stderr, "  minSpx = %d\n", minSpx);
    (void )fprintf(stderr, "  maxAng = %g (Radians)\n", maxAng);
    (void )fprintf(stderr, "  maxDeform = %g\n", maxDeform);
    (void )fprintf(stderr, "  maxDisp = %g\n", maxDisp);
    (void )fprintf(stderr, "  matchImpThr = %g\n", matchImpThr);
    (void )fprintf(stderr, "  matchImpNN = %d\n", matchImpNN);
    (void )fprintf(stderr, "  ctrFileBaseStr = %s\n",
		   ctrFileBaseStr? ctrFileBaseStr: nullStr);
    (void )fprintf(stderr, "  noMatching = %d\n", noMatching);
    (void )fprintf(stderr, "  usage = %d\n", usage);
  }
  /* Create the initial affine transform and it's inverse. */
  if(ok)
  {
    if(inTrFileStr)
    {
      tObj0 = NULL;
      if(((fP = fopen(inTrFileStr, "r")) == NULL) ||
         ((tObj0 = WlzReadObj(fP, &errNum)) == NULL) ||
	 (tObj0->type != WLZ_AFFINE_TRANS) ||
	 (tObj0->domain.core == NULL))
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
        inTr = WlzAffineTransformCopy(tObj0->domain.t, &errNum);
      }
      if(fP)
      {
        (void )fclose(fP);
	fP = NULL;
      }
      if(tObj0)
      {
        WlzFreeObj(tObj0);
	tObj0 = NULL;
      }
    }
    else
    {
      inTr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzAffineTransformPrimSet(inTr, inTrPrim);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(verbose)
      {
        (void )fprintf(stderr, "Affine transform inTr = \n");
	(void )AlcDouble2WriteAsci(stderr, inTr->mat, 3, 3);
      }
      invTr = WlzAffineTransformInverse(inTr, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to set initial affine transform (%s).\n",
		     *argv, errMsg);
    }
    else
    {
      if(verbose)
      {
        (void )fprintf(stderr, "Affine transform invTr = \n");
	(void )AlcDouble2WriteAsci(stderr, invTr->mat, 3, 3);
      }
    }
  }
  /* If the reference file name was given on the command line then it use it
   * in place of the reference file name in the section parameters bibfile. */
  if(ok && refObjFileStr)
  {
    if(((refObj3D = WlzAssignObject(
		    WlzEffReadObj(NULL, refObjFileStr, refObjFileType,
				  0, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s Failed to read the reference object from file %s (%s).\n",
      *argv, refObjFileStr, errMsg);
    }
  }
  if(ok)
  {
    if((lFP = (strcmp(inFileStr, "-")?
             fopen(inFileStr, "r"): stdin)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		   "%s Failed to open the input section parameters file %s.\n",
		    *argv, inFileStr);
    }
  }
  /* Get each input section parameters file and process it. */
  index = 0;
  while(ok &&
        ((fP = WlzMatchICPPlaneSecParFile(lFP, multipleFiles, index,
					  secParFile)) != NULL))
  {
    errNum = WlzMatchICPPlaneReadSecParam(fP, &view,
	refObjFileStr? NULL: &refObjFileStr, &refObjFileType,
	refObj3D? NULL: &refObj3D,
	&srcObjFileStr, &srcObjFileType, &srcObj2D);
    if(errNum == WLZ_ERR_NONE)
    {
      if(verboseObj)
      {
	if(verbose)
	{
	  (void )fprintf(stderr,
			 "Writting srcObj2DOrg to dbg-srcObj2DOrg.wlz.\n");
	}
	if((vFP = fopen("dbg-srcObj2DOrg.wlz", "w")) != NULL)
	{
	  (void )WlzWriteObj(vFP, srcObj2D);
	  (void )fclose(vFP);
	}
      }
      errNum = WlzInit3DViewStruct(view, refObj3D);
    }
    if(fP)
    {
      fclose(fP);
      fP = NULL;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
       "%s Failed to read input section parameters from file %s (%s).\n",
       argv[0], secParFile, errMsg);
    }
    /* Create the section object. */
    if(ok)
    {
      refObj2D = WlzGetSectionFromObject(refObj3D, view, interp, &errNum);
      if((errNum == WLZ_ERR_NONE) && (refObj2D != NULL) &&
	  (refObj2D->type = WLZ_2D_DOMAINOBJ) && (refObj2D->domain.core))
      {
	refObj2DOrg.vtX = refObj2D->domain.i->kol1;
	refObj2DOrg.vtY = refObj2D->domain.i->line1;
	if(verboseObj)
	{
	  if(verbose)
	  {
	    (void )fprintf(stderr,
			   "Writting refObj2DOrg to dbg-refObj2DOrg.wlz.\n");
	  }
	  if((vFP = fopen("dbg-refObj2DOrg.wlz", "w")) != NULL)
	  {
	    (void )WlzWriteObj(vFP, refObj2D);
	    (void )fclose(vFP);
	  }
	}
      }
      else
      {
	ok = 0;
        (void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	"%s Failed to get section from reference object (%s).\n",
	argv[0], errMsg);
      }
    }
    /* Create the contour objects. */
    if(ok)
    {
      refCObj2D = WlzMatchICPPlaneCreateContourObj(refObj2D, refMedianSz,
	  refSmooth, NULL, refCThr, minSpx, debug,
	  verboseObj? "dbg-refObj2D.wlz": NULL,
	  &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	srcCObj2D = WlzMatchICPPlaneCreateContourObj(srcObj2D, srcMedianSz,
	    srcSmooth, inTr, srcCThr, minSpx, debug,
	    verboseObj? "dbg-srcObj2D.wlz": NULL,
	    &errNum);
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
        (void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr, "%s Failed to compute contours (%s).\n",
		       argv[0], errMsg);
      }
    }
    /* Compute translation using centre of mass if required. Then modify
     * the initial and inverse transforms. */
    if(ok && useCOfM)
    {
      if(verbose)
      {
        (void )fprintf(stderr,
        "Using centre of mass to refine initial affine transform.\n");
      }
      refCOfM = WlzCentreOfMass2D(refCObj2D, 0, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        srcCOfM = WlzCentreOfMass2D(srcCObj2D, 0, NULL, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(verbose)
	{
	  (void )fprintf(stderr,
	  "centres of mass are: refCOfM = {%g,%g}, srcCOfM = {%g,%g}.\n",
	  refCOfM.vtX, refCOfM.vtY, srcCOfM.vtX, srcCOfM.vtY);
	}
	WLZ_VTX_2_SUB(tDV0, refCOfM, srcCOfM);
	cOfMTr = WlzAffineTransformFromTranslation(WLZ_TRANSFORM_2D_AFFINE,
						   tDV0.vtX, tDV0.vtY, 0.0,
						   &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(verbose)
	{
	  (void )fprintf(stderr, "Affine transform cOfMTr = \n");
	  (void )AlcDouble2WriteAsci(stderr, cOfMTr->mat, 3, 3);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* Transform the source model using the centre of mass offset. */
        (void )WlzAffineTransformGMModel(srcCObj2D->domain.ctr->model,
					 cOfMTr, 0, &errNum);
      }
      /* Now modify the inverse transform. */
      if(errNum == WLZ_ERR_NONE)
      {
        /* Invert the centre of mass offset transform. No need to call
	 * a function just invert the translation matrix elements. */
	cOfMTr->mat[0][2] = -(cOfMTr->mat[0][2]);
	cOfMTr->mat[1][2] = -(cOfMTr->mat[1][2]);
	/* ... and modify the inverse transform. */
	tTr0 = WlzAffineTransformProduct(cOfMTr, invTr, &errNum);
	(void )WlzFreeAffineTransform(invTr);
	invTr = tTr0;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(verbose)
	{
	  (void )fprintf(stderr, "Affine transform invTr = \n");
	  (void )AlcDouble2WriteAsci(stderr, invTr->mat, 3, 3);
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	"%s Failed to modify transforms using centre of mass (%s).\n",
		       argv[0], errMsg);
      }
    }
    if(ok && ctrFileBaseStr)
    {
      if(multipleFiles)
      {
        sprintf(fileNameBuf, "%s_%06d_ctr_ref.wlz", ctrFileBaseStr, index);
      }
      else
      {
        sprintf(fileNameBuf, "%s_ctr_ref.wlz", ctrFileBaseStr);
      }
      if((fP = fopen(fileNameBuf, "w")) == NULL)
      {
	ok = 0;
	(void )fprintf(stderr, "%s Failed to open contour file %s\n",
		       argv[0], fileNameBuf);
      }
      else
      {
	if((errNum = WlzWriteObj(fP, refCObj2D)) != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	  		 "%s Failed to write contour to file %s (%s)\n",
			 argv[0], fileNameBuf, errMsg);
	}
	(void )fclose(fP);
	fP = NULL;
      }
      if(ok)
      {
	if(multipleFiles)
	{
	  sprintf(fileNameBuf, "%s_%06d_ctr_src.wlz", ctrFileBaseStr, index);
	}
	else
	{
	  sprintf(fileNameBuf, "%s_ctr_src.wlz", ctrFileBaseStr);
	}
	if((fP = fopen(fileNameBuf, "w")) == NULL)
	{
	  ok = 0;
	  (void )fprintf(stderr, "%s Failed to open contour file %s\n",
			 argv[0], fileNameBuf);
	}
	else
	{
	  if((errNum = WlzWriteObj(fP, srcCObj2D)) != WLZ_ERR_NONE)
	  {
	    ok = 0;
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
	    		   "%s Failed to write contour to file %s (%s)\n",
			   argv[0], fileNameBuf, errMsg);
	  }
	  (void )fclose(fP);
	  fP = NULL;
	}
      }
    }
    if(ok && (noMatching == 0))
    {
#ifdef WLZ_MATCHICPPLANE_OWNCB
      /* Set up weighting function callback data. */
      cbData.tGM = refCObj2D->domain.ctr->model;
      cbData.sGM = srcCObj2D->domain.ctr->model;
      cbData.maxDisp = maxDisp;
      cbData.nScatter = nScatter;
      errNum = WlzMatchICPCtr(refCObj2D->domain.ctr, srcCObj2D->domain.ctr,
			      NULL, maxItr, minSpx, minSegSpx,
			      &nMatch, &matchRP, &matchSP, decompLimit,
			      maxDisp, maxAng, maxDeform,
			      matchImpNN, matchImpThr,
			      WlzMatchICPWeightMatches, &cbData);
#else
      errNum = WlzMatchICPObjs(refCObj2D, srcCObj2D, NULL,
			       &nMatch, &matchRP, &matchSP,
			       maxItr, minSpx, minSegSpx, decompLimit,
			       maxDisp, maxAng, maxDeform,
			       matchImpNN, matchImpThr);
#endif
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s Failed to compute tie-points from contours (%s).\n",
		       argv[0], errMsg);
      }
    }
    if(ok && ctrFileBaseStr)
    {
      if(multipleFiles)
      {
	sprintf(fileNameBuf, "%s_%06d_ctr_dcp.wlz", ctrFileBaseStr, index);
      }
      else
      {
	sprintf(fileNameBuf, "%s_ctr_dcp.wlz", ctrFileBaseStr);
      }
      if((fP = fopen(fileNameBuf, "w")) == NULL)
      {
	ok = 0;
	(void )fprintf(stderr, "%s Failed to open contour file %s\n",
		       argv[0], fileNameBuf);
      }
      else
      {
	if((errNum = WlzWriteObj(fP, srcCObj2D)) != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	  		 "%s Failed to write contour to file %s (%s)\n",
			 argv[0], fileNameBuf, errMsg);
	}
	(void )fclose(fP);
	fP = NULL;
      }
    }
    /* Write the new section parameters file together with the computed
     * tie-points. */
    if(ok)
    {
      if(multipleFiles)
      {
        sprintf(fileNameBuf, "%s_%06d.bib", outFileBaseStr, index);
      }
      else
      {
        sprintf(fileNameBuf, "%s.bib", outFileBaseStr);
      }
      if((fP = fopen(fileNameBuf, "w")) == NULL)
      {
	ok = 0;
	(void )fprintf(stderr,
		 "%s Failed to open the output section parameters file %s.\n",
		 *argv, fileNameBuf);
      }
    }
    if(ok)
    {
      errNum = WlzMatchICPPlaneWriteSecParam(fP, view,
	  refObjFileStr, refObjFileType,
	  srcObjFileStr, srcObjFileType,
	  invTr,
	  nMatch, matchRP.d2, matchSP.d2,
	  removeRefOrg? &refObj2DOrg: NULL);
      if(fP && strcmp(outFileBaseStr, "-"))
      {
	fclose(fP);
	fP = NULL;
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	"%s Failed to write output section parameters to file %s (%s).\n",
	argv[0], fileNameBuf, errMsg);
      }
    }
    ++index;
  }
  (void )WlzFree3DViewStruct(view);
  AlcFree(refObjFileStr);
  AlcFree(srcObjFileStr);
  AlcFree(matchRP.v);
  AlcFree(matchSP.v);
  (void )WlzFreeAffineTransform(inTr);
  (void )WlzFreeAffineTransform(invTr);
  (void )WlzFreeAffineTransform(cOfMTr);
  (void )WlzFreeObj(refObj3D);
  (void )WlzFreeObj(refObj2D);
  (void )WlzFreeObj(srcObj2D);
  (void )WlzFreeObj(refCObj2D);
  (void )WlzFreeObj(srcCObj2D);
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s",
      *argv,
      " [-h] [-d] [-v] [-V]\n"
      "          [-o<output file base>] [-r <reference file>] [-Y]\n" 
      "          [-t#] [-x#] [-y#] [-a#] [-s#] [-e]\n"
      "          [-g#,#] [-k#,#] [-u#,#]\n"
      "          [-f] [-i#] [-s#] [-A#] [-S#] [-F#] [-m#] [-n#]\n"
      "          [-c<contour file base>] [-N] [-L]\n"
      "          [<section parameters file>]\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -d  Perform extra tests to aid debuging.\n"
      "  -v  Be verbose, lots of text output to the standard error output!\n"
      "  -V  Be verbose with Woolz objects, objects written file name names\n"
      "      prefixed by dbg-.\n"
      "  -o  Output file base.\n"
      "  -r  Reference file.\n"
      "  -Y  The section parameters file is a list of section parameters\n"
      "      files, one per line.\n"
      "  -e  Use centres of mass of geometric models to compute translation.\n"
      "  -t  Initial affine transform (if given all initial affine\n"
      "      transform primitives are ignored).\n"
      "  -x  Initial horizontal translation.\n"
      "  -y  Initial vertical translation.\n"
      "  -a  Initial angle of rotation (degrees).\n"
      "  -s  Initial scale factor.\n"
      "  -g  Maximal gradient contour thresholds, with the format:\n"
      "      <ref threshold>,<src threshold>, either may be omitted.\n"
      "  -k  Median filter size, with the format:\n"
      "      <ref size>, <src size>, either may be omitted.\n"
      "  -u  Gaussian smoothing factors, with the format:\n"
      "      <ref smooth>,<src smooth>, either may be omitted.\n"
      "  -f  Keep the reference offset in the reference tie-points (MAPaint\n"
      "      doesn't do this.\n"
      "  -i  Maximum number of iterations.\n"
      "  -p  Minimum number of simplices per shell.\n"
      "  -P  Minimum number of simplices per matched shell segment, with a\n"
      "      pair of correspondence points possibly being generated per\n"
      "      matched shell segment.\n"
      "  -A  Maximum angle (degrees) from a global transformation.\n"
      "  -S  Maximum displacement from a global transformed position.\n"
      "  -F  Maximum deformation from a global transformation.\n"
      "  -m  Implausibility threshold for rejecting implausible\n"
      "      correspondence points which should be greater than zero,\n"
      "      although the useful range is probably [0.5-5.0]. Higher\n"
      "      values allow more implausible matches to be returned.\n"
      "  -n  Number of match points in neighbourhood when checking the\n"
      "	     plausibility of the correspondence points.\n"
      "  -c  Outputs the computed and decomposed geometric models using\n"
      "	     the given file base.\n"
      "  -N  Don't compute the tie-points.\n"
      "  -L  Use linear interpolation (instead of nearest neighbour) when\n"
      "      cuting sections.\n"
      "  Reads a 3D reference object and a MAPaint section parameters file.\n"
      "Computes tie-points and then writes a new MAPaint section parameters\n"
      "file which includes the tie-points.\n"
      "  An initial affine transform is computed from the (optional) initial\n"
      "translation, rotation and scale parameters. A centre of mass can be\n"
      "computed to improve the initial translation estimates. If a centre of\n"
      "mass computation is used then the images must have background with\n"
      "high values and foreground with low values. This initial affine\n"
      "transform is appiled to the source image before computing the\n"
      "tie-points. Once computed, the tie-points are transformed using the\n"
      "inverse of the intial transform, before output.\n"
      "  The tie-points are computed using an ICP based matching algorithm\n"
      "in which geometric models built from the maximal gradient edges\n"
      "extracted from the computed section of the reference object and\n"
      "the source object (refered to in the section parameters file).\n"
      "  To aid rejection of poor tie-points, the tie-points are ranked by\n"
      "plausibility, with the most plausible first\n");
  }
  return(!ok);
}

/*!
* \return	Woolz error code.
* \brief	Read the input section parameters file to find the reference
* 		and source file object file names, read the reference and
* 		source objects from these files and parse the section view
* 		parameters.
* \param	fP			Input section parameters file.
* \param	dstView			Destination pointer for the 3D view
*					transform.
* \param	dstRefFileStr		Destination pointer for the reference
* 					object file name. If the reference file
* 					string is NULL and the reference object
* 					is not NULL then the reference object
* 					is assumed valid.
* \param	dstRefFileType		Destination pointer for the reference
*					object file type.
* \param	dstRefObj		Destination pointer for the reference
*					object.
* \param	dstSrcFileStr		Destination pointer for the source
* 					object file name. If the source file
*					string is NULL and the source object
*					is not NULL then the source object is
*					assumed valid.
* \param	dstSrcFileType		Destination pointer for the source
*					object file type.
* \param	dstsrcObj		Destination pointer for the source
*					object.
*/
static WlzErrorNum WlzMatchICPPlaneReadSecParam(FILE *fP,
		      WlzThreeDViewStruct **dstView,
		      char **dstRefFileStr, WlzEffFormat *dstRefFileType,
		      WlzObject **dstRefObj, 
		      char **dstSrcFileStr, WlzEffFormat *dstSrcFileType,
		      WlzObject **dstSrcObj)
{
  int		idx;
  char		*fileStr;
  WlzEffFormat	fileType;
  FILE		*wFP;
  BibFileRecord	*bibRec;
  BibFileError  bibErr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  *dstView = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    bibErr = BibFileRecordRead(&bibRec, NULL, fP);
  }
  while((bibErr == BIBFILE_ER_NONE) && (errNum == WLZ_ERR_NONE))
  {
    if(!strcmp(bibRec->name, "Wlz3DSectionViewParams"))
    {
      errNum = WlzEffBibParse3DSectionViewParamsRecord(bibRec, *dstView);
    }
    if(!strcmp(bibRec->name, "MAPaintWarpInputSourceFile"))
    {
      errNum = WlzEffBibParseFileRecord(bibRec, &idx, &fileStr, &fileType);
      if(dstSrcFileStr)
      {
        *dstSrcFileStr = fileStr;
	if(dstSrcFileType)
	{
	  *dstSrcFileType = fileType;
	}
	if(dstSrcObj)
	{
	  *dstSrcObj = WlzAssignObject(
		       WlzEffReadObj(NULL, fileStr, fileType,
		       		     0, &errNum), NULL);
	}
      }
    }
    if(!strcmp(bibRec->name, "MAPaintWarpInputReferenceFile"))
    {
      errNum = WlzEffBibParseFileRecord(bibRec, &idx, &fileStr, &fileType);
      if(dstRefFileStr)
      {
        *dstRefFileStr = fileStr;
	if(dstRefFileType)
	{
	  *dstRefFileType = fileType;
	}
	if(dstRefObj)
	{
	  *dstRefObj = WlzAssignObject(
		       WlzEffReadObj(NULL, fileStr, fileType, 0,
		       		     &errNum), NULL);
	}
      }
    }
    BibFileRecordFree(&bibRec);
    bibErr = BibFileRecordRead(&bibRec, NULL, fP);
    if((errNum == WLZ_ERR_NONE) && (bibErr != BIBFILE_ER_NONE))
    {
      if(bibErr != BIBFILE_ER_EOF)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Outputs a section parameters file with the reference and source
* 		file object file names, the 3D view transform and the
*		tie-points.
* \param	fP			Output section parameters file.
* \param	view			3D view transform.
* \param	refObjFileStr		Reference object file name.
* \param	refObjFileType		Reference object file type.
* \param	srcObjFileStr		Source object file name.
* \param	srcObjFileType		Source object file type.
* \param	srcTr			Affine transform with which to
* 					transform the source tie-points
*					before output, NULL is equivalent to
*					the identity transform.
* \param	nMatch			Number of tie-points.
* \param	tieRP			Reference object tie-points.
* \param	tieSP			Source object tie-points.
* \param	refObj2DOrg		Origin of the 2D reference object to
*					be subtracted if non NULL.
*/
static WlzErrorNum WlzMatchICPPlaneWriteSecParam(FILE *fP,
			      WlzThreeDViewStruct *view,
			      char *refObjFileStr, WlzEffFormat refObjFileType,
			      char *srcObjFileStr, WlzEffFormat srcObjFileType,
			      WlzAffineTransform *srcTr, int nMatch,
			      WlzDVertex2 *tieRP, WlzDVertex2 *tieSP,
			      WlzDVertex2 *refObj2DOrg)
{
  int		idx;
  char		*tmpS,
		*dateS = NULL,
		*hostS = NULL,
		*userS = NULL,
		*refFileS = NULL,
		*srcFileS = NULL,
		*sgnlFileS = NULL;
  time_t	tmpTime;
  BibFileRecord	*bibRec;
  WlzDVertex2	tVx;
  WlzDVertex3	refVx,
  		srcVx;
  char		tmpBufS[256];
  static char	unknownS[] = "unknown";
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((bibRec = BibFileRecordMake("Ident", "0",
  			     BibFileFieldMakeVa("Text",
					"MAPaint 2D warp input parameters",
					"Version", "1",
					NULL))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  /* Comment with user, machine, date etc. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(BibFileRecordWrite(fP, NULL, bibRec) != BIBFILE_ER_NONE)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
    BibFileRecordFree(&bibRec);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tmpS = getenv("USER");
    (void )sprintf(tmpBufS, "User: %s", tmpS? tmpS: unknownS);
    if((userS = AlcStrDup(tmpBufS)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tmpTime = time(NULL);
    tmpS = ctime(&tmpTime);
    *(tmpS + strlen(tmpS) - 1) = '\0';
    (void )sprintf(tmpBufS, "Date: %s", tmpS? tmpS: unknownS);
    if((dateS = AlcStrDup(tmpBufS)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tmpS = getenv("HOST");
    (void )sprintf(tmpBufS, "Host: %s", tmpS? tmpS: unknownS);
    hostS = AlcStrDup(tmpBufS);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )sprintf(tmpBufS, "RefFile: %s", refObjFileStr);
    refFileS = AlcStrDup(tmpBufS);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )sprintf(tmpBufS, "SrcFile: %s", srcObjFileStr);
    srcFileS = AlcStrDup(tmpBufS);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )sprintf(tmpBufS, "SignalFile: %s", unknownS);
    sgnlFileS = AlcStrDup(tmpBufS);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((bibRec = BibFileRecordMake("Comment", "0",
			     BibFileFieldMakeVa("Text", userS,
					        "Text", dateS,
					        "Text", hostS,
					        "Text", refFileS,
					        "Text", srcFileS,
					        "Text", sgnlFileS,
					        NULL))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  AlcFree(dateS);
  AlcFree(hostS);
  AlcFree(userS);
  AlcFree(refFileS);
  AlcFree(srcFileS);
  AlcFree(sgnlFileS);
  if(errNum == WLZ_ERR_NONE)
  {
    if(BibFileRecordWrite(fP, NULL, bibRec) != BIBFILE_ER_NONE)
    {
      WLZ_ERR_WRITE_INCOMPLETE;
    }
    BibFileRecordFree(&bibRec);
  }
  /* Source file string. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffBibWriteFileRecord(fP, "MAPaintWarpInputSourceFile",
				      srcObjFileStr, srcObjFileType);
  }
  /* Reference file string. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzEffBibWriteFileRecord(fP, "MAPaintWarpInputReferenceFile",
				      refObjFileStr, refObjFileType);
  }
  /* View parameters. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzEffBibWrite3DSectionViewParamsRecord(fP, "Wlz3DSectionViewParams",
					    view);
    
  }
  /* Tie points. */
  idx = 0;
  refVx.vtZ = 0.0;
  srcVx.vtZ = 0.0;
  while((errNum == WLZ_ERR_NONE) && (idx < nMatch))
  {
    refVx.vtX = (tieRP + idx)->vtX;
    refVx.vtY = (tieRP + idx)->vtY;
    if(refObj2DOrg)
    {
      refVx.vtX -= refObj2DOrg->vtX;
      refVx.vtY -= refObj2DOrg->vtY;
    }
    tVx = WlzAffineTransformVertexD2(srcTr, *(tieSP + idx), &errNum);
    srcVx.vtX = tVx.vtX;
    srcVx.vtY = tVx.vtY;
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzEffBibWriteTiePointVtxsRecord(fP, "WlzTiePointVtxs", idx,
	  refVx, srcVx, 0);
    }
    ++idx;
  }
  return(errNum);
}

/*!
* \return	New contour object.
* \brief	Create a contour object from a 2D domain object with values,
*		together with a set of parameters.
* \param	gObj			Given 2D domain object.
* \param	medianSz		Median filter size if > 0.
* \param	smooth			Gaussian smoothing value.
* \param	tr			Initial affine transform.
* \param	cThr			Contour threshold value.
* \param	minSpx			Minimum number of simplicies per shell.
* \param	objDbgFileName		If non-null used as the name of a
*					file for the image object just prior
*					to computing the geometric model.
* \param	dstErr			Destination ptr for error, may be NULL.
*/
static WlzObject *WlzMatchICPPlaneCreateContourObj(WlzObject *gObj,
					int medianSz, double smooth,
					WlzAffineTransform *tr,
					double cThr, int minSpx,
					int debug, char *objDbgFileName,
					WlzErrorNum *dstErr)
{
  WlzObject	*tObj0 = NULL,
  		*tObj1 = NULL,
		*cObj = NULL;
  WlzDomain	tDom;
  WlzValues	tVal;
  FILE		*dFP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const		nrmFlg = 1;

  tVal.core = NULL;
  if(medianSz > 0)
  {
    errNum = WlzRankFilter(gObj, medianSz, 0.5);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(tr)
    {
      tObj0 = WlzAssignObject(
      	      WlzAffineTransformObj(gObj, tr, WLZ_INTERPOLATION_NEAREST,
	                            &errNum), NULL);
    }
    else
    {
      tObj0 = WlzAssignObject(gObj, NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj1 = WlzAssignObject(
	    WlzGauss2(tObj0, smooth, smooth, 0, 0, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(objDbgFileName)
    {
      if((dFP = fopen(objDbgFileName, "w")) != NULL)
      {
	(void )WlzWriteObj(dFP, tObj1);
	(void )fclose(dFP);
      }
    }
    tDom.ctr = WlzContourObj(tObj1, WLZ_CONTOUR_MTD_GRD, cThr, 1.0, nrmFlg,
			     &errNum);
  }
  WlzFreeObj(tObj0);
  WlzFreeObj(tObj1);
  if(errNum == WLZ_ERR_NONE)
  {
    cObj = WlzMakeMain(WLZ_CONTOUR, tDom, tVal, NULL, NULL, &errNum);
  }
  /* There's a bug somewhere in the deletion of small shells. Delete small
   * shells and then copy the contours. */
  if(debug && (errNum == WLZ_ERR_NONE))
  {
    errNum = WlzGMVerifyModel(cObj->domain.ctr->model, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMFilterRmSmShells(cObj->domain.ctr->model, minSpx * 3);
  }
  if(debug && (errNum == WLZ_ERR_NONE))
  {
    errNum = WlzGMVerifyModel(cObj->domain.ctr->model, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj0 = WlzAssignObject(WlzCopyObject(cObj, &errNum), NULL);
    WlzFreeObj(cObj);
    if(errNum == WLZ_ERR_NONE)
    {
      cObj = tObj0;
    }
    tObj0 = NULL;
  }
  if(debug && (errNum == WLZ_ERR_NONE))
  {
    errNum = WlzGMVerifyModel(cObj->domain.ctr->model, NULL);
  }
  return(cObj);
}

/*!
* \return	File pointer for section parameters or NULL.
* \brief	Gets the section parameters file pointer or NULL if no
*		section parameters files are left.
* \param	lFP			Pointer to open file.
* \param	multi			Multiple files if non-zero.
* \param	index			File index, starts with zero.
* \param	secParFile		Pointer for section parameters file
*					name.
*/
static FILE	*WlzMatchICPPlaneSecParFile(FILE *lFP, int multi, int index,
					    char *secParFile)
{
  int		len;
  FILE		*fP = NULL;

  if(multi)
  {
    if(fgets(secParFile, 255, lFP) != NULL)
    {
      *(secParFile + 255) = '\0';
      if((len = strlen(secParFile)) > 0)
      {
	if(*(secParFile + len - 1) == '\n')
	{
	  *(secParFile + len - 1) = '\0';
	}
	fP = fopen(secParFile, "r");
      }
    }
  }
  else
  {
    fP = index? NULL: lFP;
  }
  return(fP);
}
