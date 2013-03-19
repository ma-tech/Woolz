#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMatchICPPlane_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzMatchICPPlane.c
* \author       Bill Hill
* \date         June 2002
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	Compute tie points for a computed plane of section in a
*               reference object and a 2D object, where these are read from
*               an 'MAPaint 2D warp input parameters' bibfile.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzmatchicpplane "WlzMatchICPPlane"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzmatchicpplane WlzMatchICPPlane
\par Name
WlzMatchICPPlane - matches objects and computes tie points using the
                   warp parameters bibfile.
\par Synopsis
\verbatim
WlzMatchICPPlane [-h] [-d] [-v] [-V]
                 [-o<output file base>] [-r <reference file>] [-Y]
		 [-t#] [-x#] [-y#] [-a#] [-s#] [-e]
		 [-g#,#] [-k#,#] [-u#,#] [-b #] [-B #]
		 [-f] [-i#] [-s#] [-A#] [-S#] [-F#] [-m#] [-n#]
		 [-c<contour file base>] [-N] [-L]
		 [<section parameters file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Perform extra tests to aid debuging.</td>
  </tr>
  <tr>
    <td>-E</td>
    <td>Tolerance in mean registration metric value.</td>
  </tr>
  <tr>
    <td>-v</td>
    <td>Be verbose, lots of text output to the standard error output!</td>
  </tr>
  <tr>
  <td>V</td> <td>Be verbose with Woolz objects and data, with file name names
      prefixed by dbg-.</td>
  </tr>
  <tr>
  <td>o</td> <td>Output file base.</td>
  </tr>
  <tr>
  <td>r</td> <td>Reference file.</td>
  </tr>
  <tr>
  <td>Y</td> <td>The section parameters file is a list of section parameters
      files, one per line.</td>
  </tr>
  <tr>
  <td>e</td> <td>Use centres of mass of geometric models to compute translation,
      with this translaton being applied after the optional initial
      affine transform.</td>
  </tr>
  <tr>
  <td>t</td> <td>Initial affine transform (if given all initial affine
      transform primitives are ignored).</td>
  </tr>
  <tr>
  <td>x</td> <td>Initial horizontal translation.</td>
  </tr>
  <tr>
  <td>y</td> <td>Initial vertical translation.</td>
  </tr>
  <tr>
  <td>a</td> <td>Initial angle of rotation (degrees).</td>
  </tr>
  <tr>
  <td>s</td> <td>Initial scale factor.</td>
  </tr>
  <tr>
  <td>g</td> <td>Maximal gradient contour thresholds, with the format:
\verbatim
      <ref threshold>,<src threshold>
\endverbatim
      either may be omitted.</td>
  </tr>
  <tr>
  <td>k</td> <td>Median filter size, with the format:
\verbatim
  <ref size>,<src size>
\endverbatim
  either may be omitted.</td>
  </tr>
  <tr>
  <td>u</td> <td>Gaussian smoothing factors, with the format:
\verbatim
  <ref smooth>,<src smooth>
\endverbatim
  either may be omitted.</td>
  </tr>
  <tr>
  <td>b</td> <td>Low threshold values used to produce binary mask images for
      matching instead of the grey valued images, with objects
      having values at or below the given thresholds.</td>
  </tr>
  <tr>
  <td>B</td> <td>High threshold values used to produce binary mask images for
      matching instead of the grey valued images, with objects
      having values above the given thresholds.</td>
  </tr>
  <tr>
  <td>f</td> <td>Keep the reference offset in the reference tie-points (MAPaint
      doesn't do this).</td>
  </tr>
  <tr>
  <td>i</td> <td>Maximum number of iterations.</td>
  </tr>
  <tr>
  <td>p</td> <td>Minimum number of simplices per shell.</td>
  </tr>
  <tr>
  <td>P</td> <td>Minimum number of simplices per matched shell segment, with a
      pair of correspondence points possibly being generated per
      matched shell segment.</td>
  </tr>
  <tr>
  <td>A</td> <td>Maximum angle (degrees) from a global transformation.</td>
  </tr>
  <tr>
  <td>S</td> <td>Maximum displacement from a global transformed position.</td>
  <tr>
  <td>F</td> <td>Maximum deformation from a global transformation.</td>
  </tr>
  <tr>
  <td>m</td> <td>Implausibility threshold for rejecting implausible
      correspondence points which should be greater than zero,
      although the useful range is probably [0.5-5.0]. Higher
      values allow more implausible matches to be returned.</td>
  </tr>
  <tr>
  <td>n</td> <td>Number of match points in neighbourhood when checking the
             plausibility of the correspondence points.</td>
  </tr>
  <tr>
  <td>c</td> <td>Outputs the computed and decomposed geometric models using
             the given file base.</td>
  </tr>
  <tr>
  <td>N</td> <td>Don't compute the tie-points.</td>
  </tr>
  <tr>
  <td>L</td> <td>Use linear interpolation (instead of nearest neighbour) when
      cuting sections.</td>
  </tr>
</table>
\par Description
WlzMatchICPPlane reads either a 2D or 3D reference object and an MAPaint section
parameters file. Computes tie-points and then writes a new MAPaint
section parameters file which includes the tie-points.

An initial affine transform is computed from the (optional) initial
translation, rotation and scale parameters. A centre of mass can be
computed to improve the initial translation estimates. If a centre of
mass computation is used then the images must have background with
high values and foreground with low values. This initial affine
transform is appiled to the source image before computing the
tie-points. Once computed, the tie-points are transformed using the
inverse of the intial transform, before output.

The tie-points are computed using an ICP based matching algorithm
in which geometric models built from the maximal gradient edges
extracted from the computed section of the reference object and
the source object (refered to in the section parameters file).

To aid rejection of poor tie-points, the tie-points are ranked by
plausibility, with the most plausible first.
\par Examples
\verbatim
\endverbatim

\par File
\ref WlzMatchICPPlane.c "WlzMatchICPPlane.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzMatchICPCtr "WlzMatchICPCtr(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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

static int			WlzMatchICPPlaneParseDPair(
				  char *arg,
				  double *val0,
				  double *val1);
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
				  int nMatch,
				  WlzDVertex2 *tieRP,
				  WlzDVertex2 *tieSP,
				  WlzDVertex2 *refObj2DOrg);
static WlzObject 		*WlzMatchICPPlaneCreateContourObj(
				  WlzObject *gObj,
				  int binFlg,
				  WlzThresholdType thrType,
				  double thrVal,
				  int medianSz,
				  double smooth,
				  double cThr,
				  int minSpx,
				  int debug,
				  char *objDbgFileName,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzMatchICPPlaneGetSection(
				  WlzObject *gRefObj,
				  WlzThreeDViewStruct *view,
				  WlzInterpolationType interp,
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
		verboseDat = 0,
		idx0,
		index,
		binFlg = 0,
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
  double	tD0,
  		tD1,
		delta = 0.01,	
  		maxAng = 30 * (ALG_M_PI / 180.0),
  		maxDeform = 0.5,
  		maxDisp = 20.0,
		refCThr = 20.0,
		srcCThr = 20.0,
		refThr = 254.0,
		srcThr = 254.0,
		matchImpThr = 2.5,
		refSmooth = 3.0,
		srcSmooth = 3.0;
  FILE		*vFP,
  		*lFP = NULL,
  		*fP = NULL;
  char		*inFileStr = NULL,
		*inTrFileStr = NULL,
  		*outFileBaseStr = NULL,
		*ctrFileBaseStr = NULL;
  char		secParFile[256];
  char		*refObjFileStr = NULL,
  		*srcObjFileStr = NULL;
  WlzThresholdType thrType = WLZ_THRESH_LOW;
  WlzEffFormat	refObjFileType = WLZEFF_FORMAT_WLZ,
  		srcObjFileType = WLZEFF_FORMAT_WLZ;
  WlzObject	*tObj0 = NULL,
  		*refObj3D = NULL,
		*refObj2D = NULL,
  		*srcObj2D = NULL,
		*refCObj2D = NULL,
  		*srcCObj2D = NULL;
  WlzDVertex2	tDV0,
  		refCOfM,
  		srcCOfM,
		refObj2DOrg;
  WlzAffineTransform *inTr = NULL,
		*inPTr = NULL;
  WlzThreeDViewStruct *view = NULL;
  WlzVertexP	matchRP,
  		matchSP;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  WlzMatchICPWeightCbData cbData;
  WlzAffineTransformPrim inTrPrim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		fileNameBuf[FILENAME_MAX];
  const char	*errMsg;
  static char	optList[] =
                "dhvVo:r:Yt:x:y:a:s:efg:k:u:b:B:i:p:A:E:S:F:P:m:n:c:NL";
  const int	nScatter = 5;
  const char	nullStr[] = "<NULL>",
  		inFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  matchRP.v = matchSP.v = NULL;
  srcCOfM.vtX = srcCOfM.vtY = 0.0;
  inFileStr = (char *)inFileStrDef;
  outFileBaseStr = (char *)outFileStrDef;
  (void )memset(&inTrPrim, 0, sizeof(WlzAffineTransformPrim));
  inTrPrim.scale = 1.0;
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
        verboseDat = 1;
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
      case 'b':
	binFlg = 1;
        thrType = WLZ_THRESH_LOW;
        usage = WlzMatchICPPlaneParseDPair(optarg, &refThr, &srcThr) != 3;
	break;
      case 'B':
	binFlg = 1;
        thrType = WLZ_THRESH_HIGH;
        usage = WlzMatchICPPlaneParseDPair(optarg, &refThr, &srcThr) != 3;
	break;
      case 'g':
        usage = WlzMatchICPPlaneParseDPair(optarg, &refCThr, &srcCThr) != 3;
	break;
      case 'k':
        usage = WlzMatchICPPlaneParseDPair(optarg, &tD0, &tD1) != 3;
	refMedianSz = WLZ_NINT(tD0);
	srcMedianSz = WLZ_NINT(tD1);
	break;
      case 'u':
        usage = WlzMatchICPPlaneParseDPair(optarg,
					   &refSmooth, &srcSmooth) != 3;
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
      case 'E':
        if((sscanf(optarg, "%lg", &delta) != 1) || (delta < 0.0))
	{
	  usage = 1;
	}
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
    (void )fprintf(stderr, "  verboseDat = %d\n", verboseDat);
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
    (void )fprintf(stderr, "  binFlg = %d\n", binFlg);
    (void )fprintf(stderr, "  thrType = %s\n",
	    (thrType == WLZ_THRESH_LOW)? "WLZ_THRESH_LOW": "WLZ_THRESH_HIGH");
    (void )fprintf(stderr, "  refThr = %g\n", refThr);
    (void )fprintf(stderr, "  srcThr = %g\n", srcThr);
    (void )fprintf(stderr, "  refCThr = %g\n", refCThr);
    (void )fprintf(stderr, "  srcCThr = %g\n", srcCThr);
    (void )fprintf(stderr, "  refMedianSz = %d\n", refMedianSz);
    (void )fprintf(stderr, "  srcMedianSz = %d\n", srcMedianSz);
    (void )fprintf(stderr, "  refSmooth = %g\n", refSmooth);
    (void )fprintf(stderr, "  srcSmooth = %g\n", srcSmooth);
    (void )fprintf(stderr, "  delta = %g\n", delta);
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
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to set initial affine transform (%s).\n",
		     *argv, errMsg);
    }
  }
  /* If the reference file name was given on the command line then it use it
   * in place of the reference file name in the section parameters bibfile. */
  if(ok && refObjFileStr)
  {
    if(((refObj3D = WlzAssignObject(
		    WlzEffReadObj(NULL, refObjFileStr, refObjFileType,
				  0, 0, 0, &errNum), NULL)) == NULL) ||
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
    /* Note: srcObj2D assigned in WlzMatchICPPlaneReadSecParam(). */
    if(errNum == WLZ_ERR_NONE)
    {
      if(verboseDat)
      {
	if(verbose)
	{
	  (void )fprintf(stderr,
			 "Writing srcObj2D to dbg-srcObj2D.wlz.\n");
	}
	if((vFP = fopen("dbg-srcObj2D.wlz", "w")) != NULL)
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
      refObj2D = WlzAssignObject(
      		 WlzMatchICPPlaneGetSection(refObj3D, view, interp,
		 		            &errNum), NULL);
      if((errNum == WLZ_ERR_NONE) && (refObj2D != NULL) &&
	  (refObj2D->type = WLZ_2D_DOMAINOBJ) && (refObj2D->domain.core))
      {
	refObj2DOrg.vtX = refObj2D->domain.i->kol1;
	refObj2DOrg.vtY = refObj2D->domain.i->line1;
	if(verboseDat)
	{
	  if(verbose)
	  {
	    (void )fprintf(stderr,
			   "Writing refObj2D to dbg-refObj2D.wlz.\n");
	  }
	  if((vFP = fopen("dbg-refObj2D.wlz", "w")) != NULL)
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
      refCObj2D = WlzMatchICPPlaneCreateContourObj(refObj2D,
				binFlg, thrType, refThr,
      				refMedianSz, refSmooth,
				refCThr, minSpx, debug,
				verboseDat? "dbg-refObj2D.wlz": NULL,
				&errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	srcCObj2D = WlzMatchICPPlaneCreateContourObj(srcObj2D,
				binFlg, thrType, srcThr,
				srcMedianSz, srcSmooth,
				srcCThr, minSpx, debug,
				verboseDat? "dbg-srcObj2D.wlz": NULL,
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
    /* Compute translation using centre of mass if required. Then create
     * modified initial and inverse transforms for this plane. */
    if(ok)
    {
      if(useCOfM)
      {
	if(verbose)
	{
	  (void )fprintf(stderr,
	  "Using centre of mass to refine initial affine transform.\n");
	}
	refCOfM = WlzCentreOfMass2D(refCObj2D, 0, NULL, &errNum);
	tObj0 = NULL;
	if(errNum == WLZ_ERR_NONE)
	{
	  tObj0 = WlzAssignObject(
		  WlzAffineTransformObj(srcCObj2D, inTr,
					WLZ_INTERPOLATION_NEAREST,
					&errNum), NULL);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  srcCOfM = WlzCentreOfMass2D(tObj0, 0, NULL, &errNum);
	}
	(void )WlzFreeObj(tObj0);
	if(errNum == WLZ_ERR_NONE)
	{
	  if(verbose)
	  {
	    (void )fprintf(stderr,
	    "centres of mass are: refCOfM = {%g,%g}, srcCOfM = {%g,%g}.\n",
	    refCOfM.vtX, refCOfM.vtY, srcCOfM.vtX, srcCOfM.vtY);
	  }
	  inPTr = WlzAffineTransformCopy(inTr, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  WLZ_VTX_2_SUB(tDV0, refCOfM, srcCOfM);
	  inPTr->mat[0][2] += tDV0.vtX;
	  inPTr->mat[1][2] += tDV0.vtY;
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	  "%s Failed to modify transform using centre of mass (%s).\n",
			 argv[0], errMsg);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(verbose)
	{
	  (void )fprintf(stderr, "Affine transform inPTr = ");
	  if(inPTr)
	  {
	    (void )fprintf(stderr, "\n");
	    (void )AlcDouble2WriteAsci(stderr, inPTr->mat, 3, 3);
	  }
	  else
	  {
	    (void )fprintf(stderr, "Identity\n");
	  }
	}
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
      errNum = WLZ_ERR_FILE_OPEN;
      ok = (fP = fopen(fileNameBuf, "w")) != NULL;
      if(ok)
      {
	ok = (errNum = WlzWriteObj(fP, refCObj2D)) == WLZ_ERR_NONE;
      }
      if(fP)
      {
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
	ok = (fP = fopen(fileNameBuf, "w")) != NULL;
      }
      if(ok)
      {
	ok = (errNum = WlzWriteObj(fP, srcCObj2D)) == WLZ_ERR_NONE;
      }
      if(fP)
      {
	(void )fclose(fP);
	fP = NULL;
      }
      tObj0 = NULL;
      if(ok)
      {
        tObj0 = WlzAffineTransformObj(srcCObj2D, inPTr,
				      WLZ_INTERPOLATION_NEAREST, &errNum);
        ok = errNum == WLZ_ERR_NONE;
      }
      if(ok)
      {
        if(multipleFiles)
	{
	  sprintf(fileNameBuf, "%s_%06d_ctr_itrsrc.wlz", ctrFileBaseStr,
	  	  index);
	}
	else
	{
	  sprintf(fileNameBuf, "%s_ctr_itrsrc.wlz", ctrFileBaseStr);
	}
	ok = (fP = fopen(fileNameBuf, "w")) != NULL;
      }
      if(ok)
      {
	ok = (errNum = WlzWriteObj(fP, tObj0)) == WLZ_ERR_NONE;
      }
      (void )WlzFreeObj(tObj0);
      if(fP)
      {
        (void )fclose(fP);
	fP = NULL;
      }
      if(!ok)
      {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s Failed to write contour to file %s (%s)\n",
		       argv[0], fileNameBuf, errMsg);
      }
    }
    if(ok && (noMatching == 0))
    {
      /* Set up weighting function callback data. */
      cbData.tGM = refCObj2D->domain.ctr->model;
      cbData.sGM = srcCObj2D->domain.ctr->model;
      cbData.maxDisp = maxDisp;
      cbData.nScatter = nScatter;
      errNum = WlzMatchICPCtr(refCObj2D->domain.ctr, srcCObj2D->domain.ctr,
			      inPTr, maxItr, minSpx, minSegSpx,
			      &nMatch, &matchRP, &matchSP, decompLimit,
			      maxDisp, maxAng, maxDeform,
			      matchImpNN, matchImpThr,
			      WlzMatchICPWeightMatches, &cbData,
			      delta);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s Failed to compute tie-points from contours (%s).\n",
		       argv[0], errMsg);
      }
    }
    if(ok && verboseDat)
    {
	if(verbose)
	{
	  (void )fprintf(stderr,
			 "Writing tie points to dbg-tiepoints.num.\n");
	}
	if((vFP = fopen("dbg-tiepoints.num", "w")) != NULL)
	{
	  for(idx0 = 0; idx0 < nMatch; ++ idx0)
	  {
	    (void )fprintf(vFP, "%g %g %g %g\n",
		       (matchSP.d2 + idx0)->vtX,
		       (matchSP.d2 + idx0)->vtY,
		       (matchRP.d2 + idx0)->vtX - (matchSP.d2 + idx0)->vtX,
		       (matchRP.d2 + idx0)->vtY - (matchSP.d2 + idx0)->vtY);
	  }
	  (void )fclose(vFP);
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
      if(verbose)
      {
	(void )fprintf(stderr, "Writing section parameters.\n");
      }
      errNum = WlzMatchICPPlaneWriteSecParam(fP, view,
	  refObjFileStr, refObjFileType,
	  srcObjFileStr, srcObjFileType,
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
    WlzFreeAffineTransform(inPTr); inPTr = NULL;
    ++index;
  }
  (void )WlzFree3DViewStruct(view);
  AlcFree(refObjFileStr);
  AlcFree(srcObjFileStr);
  AlcFree(matchRP.v);
  AlcFree(matchSP.v);
  (void )WlzFreeAffineTransform(inTr);
  (void )WlzFreeObj(refObj3D);
  (void )WlzFreeObj(refObj2D);
  (void )WlzFreeObj(srcObj2D);
  (void )WlzFreeObj(refCObj2D);
  (void )WlzFreeObj(srcCObj2D);
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s%s%s%s",
      *argv,
      " [-h] [-d] [-v] [-V]\n"
      "                        [-o<output file base>] [-r <reference file>]\n"
      "                        [-Y] [-t#] [-x#] [-y#] [-a#] [-s#] [-e]\n"
      "                        [-g#,#] [-k#,#] [-u#,#] [-b #] [-B #]\n"
      "                        [-f] [-i#] [-s#] [-A#] [-S#] [-F#] [-m#] [-n#]\n"
      "                        [-c<contour file base>] [-N] [-L]\n"
      "                        [<section parameters file>]\n"
      "Version: ",
      WlzVersion(),
      "\n"
      "Options:\n"
      "  -h  Prints this usage information.\n"
      "  -d  Perform extra tests to aid debuging.\n"
      "  -E  Tolerance in mean registration metric value.\n"
      "  -v  Be verbose, lots of text output to the standard error output!\n"
      "  -V  Be verbose with Woolz objects and data, with file name names\n"
      "      prefixed by dbg-.\n"
      "  -o  Output file base.\n"
      "  -r  Reference file.\n"
      "  -Y  The section parameters file is a list of section parameters\n"
      "      files, one per line.\n"
      "  -e  Use centres of mass of geometric models to compute translation,\n"
      "      with this translaton being applied after the optional initial\n"
      "      affine transform.\n"
      "  -t  Initial affine transform (if given all initial affine\n"
      "      transform primitives are ignored).\n"
      "  -x  Initial horizontal translation.\n"
      "  -y  Initial vertical translation.\n"
      "  -a  Initial angle of rotation (degrees).\n"
      "  -s  Initial scale factor.\n"
      "  -g  Maximal gradient contour thresholds, with the format:\n"
      "      <ref threshold>,<src threshold>, either may be omitted.\n"
      "  -k  Median filter size, with the format:\n"
      "      <ref size>,<src size>, either may be omitted.\n"
      "  -u  Gaussian smoothing factors, with the format:\n"
      "      <ref smooth>,<src smooth>, either may be omitted.\n"
      "  -b  Low threshold values used to produce binary mask images for\n"
      "      matching instead of the grey valued images, with objects\n"
      "      having values at or below the given thresholds.\n"
      "  -B  High threshold values used to produce binary mask images for\n"
      "      matching instead of the grey valued images, with objects\n"
      "      having values above the given thresholds.\n"
      "  -f  Keep the reference offset in the reference tie-points (MAPaint\n"
      "      doesn't do this).\n"
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
      "  Reads either a 2D or 3D reference object and an MAPaint section\n"
      "parameters file. Computes tie-points and then writes a new MAPaint\n"
      "section parameters file which includes the tie-points.\n"
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
      "plausibility, with the most plausible first.\n");
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
  BibFileRecord	*bibRec;
  BibFileError  bibErr = BIBFILE_ER_NONE;
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
		       		     0, 0, 0, &errNum), NULL);
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
		       WlzEffReadObj(NULL, fileStr, fileType,
		       		     0, 0, 0, &errNum), NULL);
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
			      int nMatch,
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
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
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
    refVx.vtX = (tieRP + idx)->vtX - view->fixed.vtX;
    refVx.vtY = (tieRP + idx)->vtY - view->fixed.vtY;
    if(refObj2DOrg)
    {
      refVx.vtX -= refObj2DOrg->vtX;
      refVx.vtY -= refObj2DOrg->vtY;
    }
    srcVx.vtX = (tieSP + idx)->vtX;
    srcVx.vtY = (tieSP + idx)->vtY;
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
* \param	binFlg			Generate contours from binary
*					(thresholded) image.
* \param	thrType			Threshold type.
* \param	thrVal			Threshold value.
* \param	medianSz		Median filter size if > 0.
* \param	smooth			Gaussian smoothing value.
* \param	cThr			Contour threshold value.
* \param	minSpx			Minimum number of simplicies per shell.
* \param	objDbgFileName		If non-null used as the name of a
*					file for the image object just prior
*					to computing the geometric model.
* \param	dstErr			Destination ptr for error, may be NULL.
*/
static WlzObject *WlzMatchICPPlaneCreateContourObj(WlzObject *gObj,
					int binFlg, WlzThresholdType thrType,
					double thrVal,
					int medianSz, double smooth,
					double cThr, int minSpx,
					int debug, char *objDbgFileName,
					WlzErrorNum *dstErr)
{
  WlzObject	*tObj0 = NULL,
		*cObj = NULL;
  WlzDomain	tDom;
  WlzValues	tVal;
  FILE		*dFP = NULL;
  WlzPixelV	thrV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	nrmFlg = 1;

  tVal.core = NULL;
  if(binFlg)
  {
    thrV.type = WLZ_GREY_DOUBLE;
    thrV.v.dbv = thrVal;
    tObj0 = WlzThreshold(gObj, thrV, thrType, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      thrV.v.dbv = 0.0;
      errNum = WlzGreySetValue(gObj, thrV);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      thrV.v.dbv = 255.0;
      errNum = WlzGreySetValue(tObj0, thrV);
    }
    (void )WlzFreeObj(tObj0);
    tObj0 = NULL;
  }
  if(medianSz > 0)
  {
    errNum = WlzRankFilter(gObj, medianSz, 0.5);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(smooth > DBL_EPSILON)
    {
      tObj0 = WlzAssignObject(
	      WlzGauss2(gObj, smooth, smooth, 0, 0, &errNum), NULL);
    }
    else
    {
      tObj0 = WlzAssignObject(gObj, NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(objDbgFileName)
    {
      if((dFP = fopen(objDbgFileName, "w")) != NULL)
      {
	(void )WlzWriteObj(dFP, tObj0);
	(void )fclose(dFP);
      }
    }
    tDom.ctr = WlzContourObj(tObj0, WLZ_CONTOUR_MTD_GRD, cThr, 1.0, nrmFlg,
			     &errNum);
  }
  WlzFreeObj(tObj0);
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

/*!
* \return	2D reference section cut from the reference object.
* \brief	Cuts a 2D reference section from the given reference object
*		which may either be a 2D or 3D object. If the given object
*		is 2D then it i shifted to respect the given 3D view structures
*		fixed point.
* \param	gRefObj			Given reference object.
* \param	view			Given 3D view data structure.
* \param	interp			Required interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzMatchICPPlaneGetSection(WlzObject *gRefObj,
				WlzThreeDViewStruct *view,
				WlzInterpolationType interp,
				WlzErrorNum *dstErr)
{
  WlzObject	*refObj2D = NULL;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if(gRefObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gRefObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	tr = WlzAffineTransformFromTranslation(WLZ_TRANSFORM_2D_AFFINE,
				-(view->fixed.vtX), -(view->fixed.vtY), 0.0,
				&errNum);
	if(errNum == WLZ_ERR_NONE)
	{
          refObj2D = WlzAffineTransformObj(gRefObj, tr, interp, &errNum);
	}
	(void )WlzFreeAffineTransform(tr);
	break;
      case WLZ_3D_DOMAINOBJ:
        refObj2D = WlzGetSectionFromObject(gRefObj, view, interp, &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(refObj2D);
}

/*!
* \return	Status of parsing indicated by:
*		  -1 if a parsing error occured,
*		   0 if neither value set,
*		   1 if just the first value set,
*		   2 if just  the second value set
*		   and 3 if both values set.
* \brief	Parses the given string for a comma seperated pair of
*		double values, each of which may be missing.
* \param	arg			Given string to parse.
* \param	dstVal0			First destination pointer.
* \param	dstVal1			Second destination pointer.
*/
static int	WlzMatchICPPlaneParseDPair(char *arg,
				           double *val0, double *val1)
{
  int		stat = 0;
  char		*parseStr[2];

  if(arg)
  {
    while(*arg && isspace(*arg))
    {
      ++arg;
    }
    if(*arg == ',')
    {
      parseStr[0] = NULL;
      parseStr[1] = strtok(arg, ",");
    }
    else
    {
      parseStr[0] = strtok(arg, ",");
      parseStr[1] = strtok(NULL, ",");
    }
    if(parseStr[0])
    {
      if(sscanf(parseStr[0], "%lg", val0) == 1)
      {
        stat = 1;
      }
      else
      {
        stat = -1;
	parseStr[1] = NULL;
      }
    }
    if(parseStr[1])
    {
      if(sscanf(parseStr[1], "%lg", val1) == 1)
      {
        stat += 2;
      }
      else
      {
        stat = -1;
      }
    }
  }
  return(stat);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
