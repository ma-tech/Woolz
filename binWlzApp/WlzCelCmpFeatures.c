#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCelCmpFeatures_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzCelCmpFeatures.c
* \author       Bill Hill
* \date         April 2008
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
* \brief	Extracts domains for clumps and cells within clumps
* 		both from within the given image. The cells and
* 		clumps are used to compute numerical features
* 		and these are output together with "colorised"
* 		images showing the segmentation.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzcelcmpfeatures "WlzCelCmpFeatures"
*/

/*!
\ingroup BinWlz
\defgroup wlzcelcmpfeatures WlzCelCmpFeatures
\par Name
WlzCelCmpFeatures - computes the features of cells within clumps.
\par Synopsis
\verbatim
WlzCelCmpFeatures [-b] [-c<col 0>[,<col 1>[,col 2]]] [-C] [-d<debug mask>]
                  [-f<input image format>] [-h]
		  [-l<col 0>[,<col 1>[,<col 2>]]] [-L] [-n#]
		  [-o<out file>] [-s<seg file base>]
		  [-S<seg file format>] [-v] [<input image>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-b</b></td>
    <td>Treats the segmented cells and clumps as binary objects, ignoring
        their image values.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Colour space for segmenting the cells. The colour is specified
        by a single letter which is one of - <b>r</b>ed, <b>g</b>reen,
	<b>b</b>lue, <b>y</b>ellow, <b>m</b>agenta, <b>c</b>yan,
	<b>h</b>ue, <b>s</b>aturation, <b>v</b>alue and gr<b>e</b>y.
	All processing is done on single (grey) valued images. Two
	methods may be used to get the grey valued image from the
	input image: If one colour is specified then that colour
	space is used to construct a grey valued image, 
	if a colour pair is given then the ratio of the first to second
	is used,
	alternatively if a colour triple is given the first two
	are used for the ratio and the third as a multiplier.
	See WlzRGBChannelRatio(1).</td>
  </tr>
  <tr> 
    td><b>-C</b></td>
    <td>Invert cell image values.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Probably only useful for debugging, this is a bit mask which
        controls the debug output.
	Mask bit 1 set for text output and mask bit 2 set for object
	output.
	The default is no debug output.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Input image format.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help - prints a usage message.</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Colour space for segmenting the clumps. The colour is specified
        by a single letter as for the cells segmentation.
  </tr>
  <tr> 
    td><b>-L</b></td>
    <td>Invert clump image values.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Number of angular increments to use (default 360).</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file for the feature values.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>The segmentation output image base is used to provide a colorised
        image showing the segmentation for each segmented clump. Unless a
	file base is given no segmentation images are output.</td>
  </tr>
  <tr> 
    <td><b>-S</b></td>
    <td>Output colorised image format.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose output.</td>
  </tr>
</table>
\par Description
Extracts domains for clumps and cells within clumps, both from within
the given image. The cells and clumps are used to compute numerical
features which are output to the given file. Colorised images showing
the segmentation may also be output if required.
The features computed are:
  The ratio of cell / clump areas.
  The ratio of clump minimum to maximum diameter.
  The centrality of the cells with respect to the clump, with the centrality
  \f[
  c = \frac{\sum_{i,j}{m_{i,j}(R_i - r_{i,j})}}
           {\sum_{i,j}{m_{i,j}R}}
  \f]
The image formats recognised are: jpg, tif and wlz.
By default the image is read from the standard input and the features are
output to the stanard output.

\par Example
\verbatim
WlzCelCmpCentrality -b -c b,r -l r,g -s image-seg-.jpg -v image.jpg
\endverbatim
Reads an image from the file image.jpg.
Segments clumps using the ratio of red to green and cells using the ratio
of blue to red.
For each clump an colorised image showing the segmentation is output
to files image-seg-000000.jpg, image-seg-000001.jpg, ....
The numerical features are preceded by column labels and output to the
standard output.
\par File
\ref WlzCelCmpFeatures.c "WlzCelCmpFeatures.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCentrality "WlzCentrality(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static void			WlzCelCmpDbgFileParse(
				  char *prog,
				  char *base,
				  char *fStr,
				  char *fPath,
				  char *fBody,
				  char *fExt,
				  WlzEffFormat fFmt);
static int			WlzCelCmpDbgObjOut(
				  WlzObject *obj,
				  char *fName,
				  char *prog,
				  char *id);
static int			WlzCelCmpParseColours(
				  WlzRGBAColorChannel *dstCNum,
				  WlzRGBAColorChannel *dstCDen,
				  WlzRGBAColorChannel *dstCMul,
				  int *dstCCnt,
				  char *gStr);
static int			WlzCelCmpParseFilename(
				  char **dstFPath,
				  char **dstFBody,
				  char **dstFExt,
				  WlzEffFormat *dstFFmt,
				  char *fStr,
				  WlzEffFormat fFmt,
				  WlzEffFormat defFFmt);
static double			**WlzCelCmpCompFeat(
				  WlzCompoundArray *cmpCAObj,
				  WlzCompoundArray *celCAObj,
				  char ***dstLabels,
				  WlzIVertex2 *dstFeatSz,
				  int nRay, int binFlg,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCelCmpExtractDomain(
				  WlzObject *iObj,
				  WlzRGBAColorChannel numC,
				  WlzRGBAColorChannel denC,
				  WlzRGBAColorChannel mulC,
				  int cRatio,
				  int invFlg,
				  int smFlg,
				  double minHSm,
				  double maxHSm,
				  double dRad,
				  double eRad,
				  double tExtra,
				  int dbgMsk,
				  char *prog,
				  char *dbgId,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCelCmpColorise(
				  WlzObject *gObj,
				  WlzObject *cmpObj,
				  WlzObject *celObj,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCelCmpNormalise(
				  WlzObject *gObj,
				  int dither,
				  WlzErrorNum *dstErr);
static WlzCompoundArray 	*WlzCelCmpSplitCmp(
				  WlzObject *cmpObj,
				  int minCmpSz,
				  WlzErrorNum *dstErr);
static WlzCompoundArray 	*WlzCelCmpSplitCel(
				  WlzObject *celObj,
				  WlzCompoundArray *cmpCAObj,
				  int minCelSz,
				  WlzErrorNum *dstErr);
/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		idN,
  		idX,
		idY,
		len,
		binFlg = 1,
		dbgMsk = 0,
  		nAng = 360,
		celCCnt = 3,
		celInv = 0,
		cmpCCnt = 3,
		cmpInv = 0,
  		option,
		ok = 1,
		usage = 0,
		vrbFlg = 0;
  WlzIVertex2	featSz;
  WlzRGBAColorChannel celCNum = WLZ_RGBA_CHANNEL_RED,
  		celCDen = WLZ_RGBA_CHANNEL_BLUE,
		celCMul = WLZ_RGBA_CHANNEL_GREY,
		cmpCNum = WLZ_RGBA_CHANNEL_GREEN,
		cmpCDen = WLZ_RGBA_CHANNEL_RED,
		cmpCMul = WLZ_RGBA_CHANNEL_GREY;
  WlzEffFormat	iFileFmt = WLZEFF_FORMAT_NONE,
  		sFileFmt = WLZEFF_FORMAT_NONE;
  WlzObject	*celObj = NULL,
  		*cmpObj = NULL,
		*gObj = NULL,
		*iObj = NULL,
		*nObj = NULL,
  		*sObj = NULL,
		*tObj0 = NULL,
		*tObj1 = NULL;
  WlzCompoundArray *celCAObj = NULL,
  		*cmpCAObj = NULL;
  WlzPixelV	cmpC,
  		celC;
  char		*iFileStr,
  		*iFilePath = NULL,
		*iFileBody = NULL,
		*iFileExt = NULL,
		*oFileStr,
		*sFileStr = NULL,
  		*sFilePath = NULL,
		*sFileBody = NULL,
		*sFileExt = NULL;
  char		**labels = NULL;
  double	**feat = NULL;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char	buffer[256],
  		optList[] = "bcC:d:f:hl:Ln:o:s:S:v",
		fileDef[] = "-";
  const char	fSep = '/';
  const char	*errMsg;
  const int	outFldWth = 16,
		minCmpSz = 10000,
		minCelSz = 100;
  const double	celDRad = 4.0,
		cmpDRad = 16.0,
  		celERad = 2.0,
  		cmpERad = 4.0,
		minHSm = 1.0,
		maxHSm = 8.0,
  		tExtra = 0.00;
  const	WlzEffFormat defFileFmt = WLZEFF_FORMAT_WLZ;

  cmpC.type = celC.type = WLZ_GREY_RGBA;
  WLZ_RGBA_RGBA_SET(cmpC.v.rgbv, 255, 0, 0, 255);
  WLZ_RGBA_RGBA_SET(celC.v.rgbv, 0, 0, 255, 255);
  /* Parse the argument list and check for input files. */
  opterr = 0;
  iFileStr = fileDef;
  oFileStr = fileDef;
  while((option = getopt(argc, argv, optList)) != EOF)
  {
    switch(option)
    {
      case 'b':
        binFlg = 1;
	break;
      case 'c':
        usage = WlzCelCmpParseColours(&celCNum, &celCDen, &celCMul, &celCCnt,
				optarg);
        break;
      case 'C':
	celInv = 1;
        break;
      case 'd':
        if(sscanf(optarg, "%d", &dbgMsk) != 1)
	{
	  usage = 1;
	}
        break;
      case 'f':
        if((iFileFmt = WlzEffStringExtToFormat(optarg)) == WLZEFF_FORMAT_NONE)
	{
	  usage = 1;
	}
        break;
      case 'l':
        usage = WlzCelCmpParseColours(&cmpCNum, &cmpCDen, &cmpCMul, &cmpCCnt,
				optarg);
        break;
      case 'L':
	cmpInv = 1;
        break;
      case 'n':
	if((sscanf(optarg, "%d", &nAng) != 1) || (nAng <= 0))
	{
	  usage = 1;
	}
	break;
      case 'o':
	oFileStr = optarg;
	break;
      case 's':
        sFileStr = optarg;
        break;
      case 'S':
        if((sFileFmt = WlzEffStringExtToFormat(optarg)) == WLZEFF_FORMAT_NONE)
	{
	  usage = 1;
	}
	break;
      case 'v':
        vrbFlg = 1;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if((dbgMsk & 1) != 0)
  {
    (void )fprintf(stderr, "%s DBG binFlg = %d.\n",
                   argv[0], binFlg);
    (void )fprintf(stderr,
                   "%s DBG celCNum = %d, celCDen = %d, celCMul = %d.\n",
		   argv[0], celCNum, celCDen, celCMul);
    (void )fprintf(stderr, "%s DBG celInv = %d.\n",
                   argv[0], celInv);
    (void )fprintf(stderr, "%s DBG dbgMsk = 0x%x.\n",
                   argv[0], dbgMsk);
    (void )fprintf(stderr, "%s DBG iFileFmt = %d.\n",
                   argv[0], iFileFmt);
    (void )fprintf(stderr,
                   "%s DBG cmpCNum = %d, cmpCDen = %d cmpCMul = %d\n",
		   argv[0], cmpCNum, cmpCDen, cmpCMul);
    (void )fprintf(stderr, "%s DBG nAng = %d.\n",
                   argv[0], nAng);
    (void )fprintf(stderr, "%s DBG oFileStr = \"%s\".\n",
                   argv[0], (oFileStr)? oFileStr: "(null)");
    (void )fprintf(stderr, "%s DBG sFileStr = \"%s\".\n",
                   argv[0], (sFileStr)? sFileStr: "(null)");
    (void )fprintf(stderr, "%s DBG sFileFmt = %d.\n",
                   argv[0], sFileFmt);
    (void )fprintf(stderr, "%s DBG vrbFlg = %d.\n",
                   argv[0], vrbFlg);
  }
  if(usage == 0)
  {
    if(sFileStr != NULL)
    {
      usage = WlzCelCmpParseFilename(&sFilePath, &sFileBody, &sFileExt,
				      &sFileFmt, sFileStr, sFileFmt,
				      defFileFmt);
      if((dbgMsk & 1) != 0)
      {
	WlzCelCmpDbgFileParse(*argv, "s", sFileStr,
			      sFilePath, sFileBody, sFileExt,
			      sFileFmt);
      }
    }
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) == argc)
    {
      iFileStr = argv[optind];
      usage = WlzCelCmpParseFilename(&iFilePath, &iFileBody, &iFileExt,
					&iFileFmt, iFileStr, iFileFmt,
					defFileFmt);
      if((dbgMsk & 1) != 0)
      {
	WlzCelCmpDbgFileParse(*argv, "i", iFileStr,
	                      iFilePath, iFileBody, iFileExt,
			      iFileFmt);
      }
    }
    else
    {
      usage = 1;
    }
  }
  ok = !usage;
  /* Read the input object. */
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if(iFileExt)
    {
      sprintf(buffer, "%s%c%s.%s", iFilePath, fSep, iFileBody, iFileExt);
    }
    else
    {
      sprintf(buffer, "%s%c%s", iFilePath, fSep, iFileBody);
    }
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr, "%s DBG Reading input object from file %s\n",
                     argv[0], buffer);
    }
    if((iObj = WlzAssignObject(WlzEffReadObj(NULL, buffer, iFileFmt, 0,
                                              &errNum), NULL)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read image from file %s (%s).\n",
                     *argv, buffer, errMsg);
    }
  }
  if(ok)
  {
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr, "%s DBG Normalising the input object.\n",
                     argv[0]);
    }
    nObj = WlzAssignObject(WlzCelCmpNormalise(iObj, 1, &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to normalise input object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok && ((dbgMsk & 2) != 0))
  {
    ok = WlzCelCmpDbgObjOut(nObj, "nObj", *argv, "");
  }
  /* Create a grey image to use when composing the segmentation images
   * if it will be required. */
  if(ok && (sFileBody != NULL))
  {
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr,
                     "%s DBG Creating grey image to use when composing the\n"
		     "segmentation images.\n",
                     argv[0]);
    }
    gObj = WlzAssignObject(WlzRGBAToModulus(nObj, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzGreyNormalise(gObj, 1);
    }
  }
  if(ok && ((dbgMsk & 2) != 0))
  {
    ok = WlzCelCmpDbgObjOut(gObj, "gObj", *argv, "");
  }
  /* Compute clump and cell domains. */
  if(ok)
  {
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr,
                     "%s DBG Extracting clump domain.\n",
                     argv[0]);
    }
    cmpObj = WlzAssignObject(
             WlzCelCmpExtractDomain(nObj, cmpCNum, cmpCDen, cmpCMul,
    				    cmpCCnt, cmpInv, 1,
				    minHSm, maxHSm,
				    cmpDRad, cmpERad, tExtra,
				    dbgMsk, *argv, "-cmp", &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to extract clump domain (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr,
                     "%s DBG Extracting cell domain.\n",
                     argv[0]);
    }
    /* I can't see why this intersection should be needed as the domains
     * should be the same, but they're not! */
    tObj1 = NULL;
    tObj0 = WlzIntersect2(cmpObj, nObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      tObj1 = WlzAssignObject(
	      WlzMakeMain(nObj->type, tObj0->domain, nObj->values, NULL,
			  NULL, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      celObj = WlzAssignObject(
      	       WlzCelCmpExtractDomain(tObj1, celCNum, celCDen, celCMul,
      				      celCCnt, celInv, 1,
				      minHSm, maxHSm,
				      celDRad, celERad,
				      tExtra, dbgMsk, *argv, "-cel",
				      &errNum), NULL);
    }
    (void )WlzFreeObj(tObj0);
    (void )WlzFreeObj(tObj1);
    tObj0 = tObj1 = NULL;
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to extract cell domain (%s).\n",
		     *argv, errMsg);
    }
  }
  /* Get split clump and cell domains. */
  if(ok)
  {
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr,
                     "%s DBG Splitting clump and cell domains.\n",
                     argv[0]);
    }
    cmpCAObj = WlzCelCmpSplitCmp(cmpObj, minCmpSz, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to split clump domain (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    celCAObj = WlzCelCmpSplitCel(celObj, cmpCAObj, minCelSz, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to split cell domain (%s).\n",
		     *argv, errMsg);
    }
  }
  /* Compute and output colour washed segmented images if required. */ 
  if(ok && (sFileBody != NULL))
  {
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr,
                     "%s DBG Computing %d colorised segmentation images.\n",
                     argv[0], cmpCAObj->n);
    }
    for(idN = 0; idN < cmpCAObj->n; ++idN)
    {
      sObj = WlzAssignObject(
	     WlzCelCmpColorise(gObj,
	                       cmpCAObj->o[idN],
	                       celCAObj->o[idN],
	                       &errNum), NULL);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to compute colorised object (%s).\n",
		       *argv, errMsg);
        break;
      }
      else
      {
	(void )sprintf(buffer, "%s%c%s%06d.%s",
	               sFilePath, fSep, sFileBody, idN, sFileExt);
	if((dbgMsk & 1) != 0)
	{
	  (void )fprintf(stderr,
			 "%s DBG Writing colorised segmentation image\n"
			 "to file \"%s\".\n",
			 argv[0], buffer);
	}
        errNum = WlzEffWriteObj(NULL, buffer, sObj, sFileFmt);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
		     "%s: Failed to write colorised object to file %s (%s).\n",
			 *argv, buffer, errMsg);
	  break;
	}
      }
      (void )WlzFreeObj(sObj);
    }
  }
  /* Compute features for split clumps and cells. */
  if(ok)
  {
    if((dbgMsk & 1) != 0)
    {
      (void )fprintf(stderr,
		     "%s DBG Computing features.\n",
		     argv[0]);
    }
    feat = WlzCelCmpCompFeat(cmpCAObj, celCAObj, &labels, &featSz,
    			     nAng, binFlg, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to compute features (%s).\n",
		     *argv, errMsg);
    }
  }
  /* Output computed features. */
  if(ok)
  {
    if((fP = (strcmp(oFileStr, "-")?
             fopen(oFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
		     argv[0], oFileStr);
    }
    /* Output labels. */
    if(vrbFlg)
    {
      for(idX = 0; idX < featSz.vtX; ++idX)
      {
	len = outFldWth - strlen(labels[idX]);
	(void )fprintf(fP, "%s", labels[idX]);
        for(idN = 0; idN < len; ++idN)
	{
	  (void )fprintf(fP, " ");
	}
      }
      (void )fprintf(fP, "\n");
    }
    /* Output features. */
    for(idY = 0; idY < featSz.vtY; ++idY)
    {
      for(idX = 0; idX < featSz.vtX; ++idX)
      {
        (void )sprintf(buffer, "%g", feat[idY][idX]);
	len = outFldWth - strlen(buffer);
        for(idN = 0; idN < len; ++idN)
	{
	  (void )strcat(buffer, " ");
	}
	(void )fprintf(fP, "%s", buffer);
      }
      (void )fprintf(fP, "\n");
    }
    if(fP && strcmp(oFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  /* Free storage. */
  if((dbgMsk & 1) != 0)
  {
    (void )fprintf(stderr,
		   "%s DBG Freeing storage.\n",
		   argv[0]);
  }
  AlcFree(iFilePath);
  AlcFree(iFileBody);
  AlcFree(iFileExt);
  AlcFree(sFilePath);
  AlcFree(sFileBody);
  AlcFree(sFileExt);
  (void )WlzFreeObj(cmpObj);
  (void )WlzFreeObj(celObj);
  (void )WlzFreeObj(gObj);
  (void )WlzFreeObj(iObj);
  (void )WlzFreeObj(nObj);
  (void )WlzFreeObj((WlzObject *)celCAObj);
  (void )WlzFreeObj((WlzObject *)cmpCAObj);
  if(labels)
  {
    for(idX = 0; idX < featSz.vtX; ++idX)
    {
      AlcFree(labels[idX]);
    }
    AlcFree(labels);
  }
  (void )AlcDouble2Free(feat);
  if(usage)
  {
    (void )fprintf(stderr, "%s help TODO\n", argv[0]);
  }
  return(!ok);
}


/*!
* \ingroup	WlzCelCmpFeatures
* \brief	Outputs information for debuging the parsing of given
* 		file string.
* \param	prog			Program name.
* \param	base			Base with which to prefix output.
* \param	fStr			File string prior to parsing.
* \param	fPath			File path.
* \param	fBody			File body.
* \param	fExt			File extension.
* \param	fFmt			File format.
*/
static void	WlzCelCmpDbgFileParse(char *prog, char *base, char *fStr,
	                              char *fPath, char *fBody, char *fExt,
			              WlzEffFormat fFmt)
{
  static char	nullStr[] = "(null)";

  (void )fprintf(stderr, "%s DBG %sFileStr  = %s.\n",
		 prog, base, (fStr)? fStr: nullStr);
  (void )fprintf(stderr, "%s DBG %sFilePath = %s.\n",
		 prog, base, (fPath)? fPath: nullStr);
  (void )fprintf(stderr, "%s DBG %sFileBody = %s.\n",
		 prog, base, (fBody)? fBody: nullStr);
  (void )fprintf(stderr, "%s DBG %sFileExt  = %s.\n",
		 prog, base, (fExt)? fExt: nullStr);
  (void )fprintf(stderr, "%s DBG %sFileFmt = %s.\n",
		 prog, base, WlzEffStringFromFormat(fFmt, NULL));
}

/*!
* \return	Non-zero if no error.
* \ingroup	WlzCelCmpFeatures
* \brief	Outputs an object for debuging.
* \param	obj			Given object.
* \param	fName			Given file name.
* \param	prog			Given program name.
* \param	id			Identifier string.
*/
static int	WlzCelCmpDbgObjOut(WlzObject *obj, char *fName,
				   char *prog, char *id)
{
  int		ok = 1;
  FILE		*fP;
  WlzErrorNum	errNum = WLZ_ERR_FILE_OPEN;
  char		buffer[256];
  const char	*errMsg;

  (void )sprintf(buffer, "%s%s.wlz", fName, id);
  if((fP = fopen(buffer, "w")) == NULL)
  {
    ok = 0;
  }
  else
  {
    if((errNum = WlzWriteObj(fP, obj)) != WLZ_ERR_NONE)
    {
      ok =  0;
    }
    (void )fclose(fP);
  }
  if(ok == 0)
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
		   "%s: Failed to write debug object %s (%s).\n",
		   prog, fName, errMsg);
  }
  return(ok);
}
/*!
* \return	Non-zero on error.
* \ingroup	WlzCelCmpFeatures
* \brief	Parses the given colour channel specification string
* 		for either a single colour channel character or a
* 		pair of comma separated colour channel characters.
* 		All destination pointers must be valid.
* \param	dstCNum			Destination pointer for numerator.
* \param	dstCDen			Destination pointer for denominator.
* \param	dstCMul			Destination pointer for multiplier.
* \param	dstCCnt		        Destination pointer for parsed
* 					colour channel count.
* 					set to non-zero value if pair parsed.
* \param	gStr			Given colour channel specification
* 					string.
*/
static int	WlzCelCmpParseColours(WlzRGBAColorChannel *dstCNum,
				WlzRGBAColorChannel *dstCDen,
				WlzRGBAColorChannel *dstCMul,
				int *dstCCnt, char *gStr)
{
  int		idx,
  		cnt = 0,
  		err = 0;
  char		*tStr;
  WlzRGBAColorChannel cC;
  int		cCC[3];

  *dstCNum = WLZ_RGBA_CHANNEL_DUMMY;
  *dstCDen = WLZ_RGBA_CHANNEL_DUMMY;
  *dstCMul = WLZ_RGBA_CHANNEL_DUMMY;
  if(gStr)
  {
    while(*gStr && isspace(*gStr))
    {
      ++gStr;
    }
    if((tStr = strtok(gStr, ",")) != NULL)
    {
      cCC[cnt++] = toupper(*tStr);
      if((tStr = strtok(NULL, ",")) != NULL)
      {
        cCC[cnt++] = toupper(*tStr);
      }
    }
    if((tStr = strtok(gStr, ",")) != NULL)
    {
      cCC[cnt++] = toupper(*tStr);
      if((tStr = strtok(NULL, ",")) != NULL)
      {
        cCC[cnt++] = toupper(*tStr);
      }
    }
    idx = 0;
    while((idx < cnt) && (err == 0))
    {
      switch(cCC[idx])
      {
	case 'E':
	  cC = WLZ_RGBA_CHANNEL_GREY;
	  break;
	case 'R':
	  cC = WLZ_RGBA_CHANNEL_RED;
	  break;
	case 'G':
	  cC = WLZ_RGBA_CHANNEL_GREEN;
	  break;
	case 'B':
	  cC = WLZ_RGBA_CHANNEL_BLUE;
	  break;
	case 'C':
	  cC = WLZ_RGBA_CHANNEL_CYAN;
	  break;
	case 'M':
	  cC = WLZ_RGBA_CHANNEL_MAGENTA;
	  break;
	case 'Y':
	  cC = WLZ_RGBA_CHANNEL_YELLOW;
	  break;
	case 'H':
	  cC = WLZ_RGBA_CHANNEL_HUE;
	  break;
	case 'S':
	  cC = WLZ_RGBA_CHANNEL_SATURATION;
	  break;
	case 'V':
	  cC = WLZ_RGBA_CHANNEL_BRIGHTNESS;
	  break;
	default:
	  err = 1;
	  break;
      }
      if(err == 0)
      {
	switch(idx)
	{
	  case 0:
	    *dstCNum = cC;
	    break;
	  case 1:
	    *dstCDen = cC;
	    break;
	  case 2:
	    *dstCMul = cC;
	    break;
	  default:
	    err = 1;
	    break;
	}
        ++idx;
      }
    }
    *dstCCnt = cnt;
  }
  else
  {
    err = 1;
  }
  return(err);
}

/*!
* \return	Non-zero on error.
* \ingroup	WlzCelCmpFeatures
* \brief	Parses the given file string and file format for the
* 		file path, file body, file extension and file format.
* 		All destination pointers must be valid.
* 		All strings need to be freed using AlcFree().
* \param	dstFPath		Destination pointer for file path.
* \param	dstFBody		Destination pointer for file body.
* \param	dstFExt			Destination pointer for file extension.
* \param	dstFFmt			Destination pointer for
* \param	dstFStr			Given file string.
* \param	fFmt			Given file format which overrides
* 					the computed file format.
* \param	defFFmt			default file format if not set by
* 					given file format or extension.
*/
static int	WlzCelCmpParseFilename(char **dstFPath, char **dstFBody,
				char **dstFExt, WlzEffFormat *dstFFmt,
				char *fStr, WlzEffFormat fFmt,
				WlzEffFormat defFFmt)
{
  int		len,
  		err = 0;
  char		*sep,
  		*dot,
		*fPath = NULL,
		*fBody = NULL,
		*fExt = NULL;
  static char	defPathDir[] = ".";
  
  *dstFPath = NULL;
  *dstFBody = NULL;
  *dstFExt = NULL;
  if((fStr == NULL) ||
     ((len = strlen(fStr)) < 3) || /* Minimum of 3 chars. */
     (*(fStr + len - 1) == '/'))  /* Directory not plain file. */
  {
    err = 1;
  }
  if(err == 0)
  {
    sep = strrchr(fStr, '/');
    dot = strrchr(fStr, '.');
    if(sep == NULL)
    {
      fPath = defPathDir;
      fBody = fStr;
    }
    else
    {
      *sep = '\0';
      fPath = fStr;
      fBody = sep + 1;
    }
    if(dot)
    {
      *dot = '\0';
      fExt = dot + 1;
    }
  }
  if(err == 0)
  {
    if((fFmt == WLZEFF_FORMAT_NONE) && (fExt != NULL) && (strlen(fExt) > 0))
    {
      fFmt = WlzEffStringExtToFormat(fExt);
    }
    if(fFmt == WLZEFF_FORMAT_NONE)
    {
      fFmt = defFFmt;
    }
  }
  if(err == 0)
  {
    *dstFFmt = fFmt;
    if(((*dstFPath = AlcStrDup(fPath)) == NULL) ||
       ((*dstFBody = AlcStrDup(fBody)) == NULL))
    {
      err = 1;
    }
  }
  if(err == 0)
  {
    if(fExt == NULL)
    {
      *dstFExt = NULL;
    }
    else if((*dstFExt = AlcStrDup(fExt)) == NULL)
    {
      err = 1;
    }
  }
  if(err != 0)
  {
    *dstFFmt = WLZEFF_FORMAT_NONE;
    AlcFree(*dstFPath); *dstFPath = NULL;
    AlcFree(*dstFBody); *dstFBody = NULL;
    AlcFree(*dstFExt); *dstFExt = NULL;
  }
  return(err);
}

/*!
* \return	Object with forground domain.
* \ingroup	WlzCelCmpFeatures
* \brief	Given a set of parameters segments the foreground from
* 		the input onbject using colour channels, thresholding
* 		and morphological operators.
* \param	iObj			Input object.
* \param	numC			Numerator (or only) clour channel.
* \param	denC			Denominator colour channel (not
* 					used if cRatio == 0).
* \param	cCnt			Number of  colour channels that are
* 					valid.
* \param	invFlg			Invert the grey data if non-zero.
* 					Without this flag set the threshold
* 					sellects high (bright) image values.
* \param	smFlg			Smooth object if non-zero.
* \param	minHSm			Minimum histogram smoothing for
* 					WlzHistogramSmooth().
* \param	minHSm			Minimum histogram smoothing for
* \param	maxHSm			Maximum histogram smoothing for
* 					WlzHistogramSmooth().
* \param	dRad			Radius of dilation structuring object
* 					used to fill small gaps.
* \param	eRad			Radius of erosion structuring object
* 					used to get rid of grot.
* \param	tExtra			Extra fraction added to the  computed
* 					threshold value.
* \param	dbgMsk			Debug output mask.
* \param	prog			Program name for debug output.
* \param	dbgId			Identifier string for debug output.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCelCmpExtractDomain(WlzObject *iObj,
				WlzRGBAColorChannel numC,
				WlzRGBAColorChannel denC,
				WlzRGBAColorChannel mulC,
				int cCnt, int invFlg,
				int smFlg,
				double minHSm, double maxHSm,
				double dRad, double eRad, double tExtra,
				int dbgMsk, char *prog, char *dbgId,
				WlzErrorNum *dstErr) 
{
  WlzGreyType	inGType;
  WlzObject	*dSEObj = NULL,
		*eObj = NULL,
		*dObj = NULL,
  		*eSEObj = NULL,
	        *hObj = NULL,
		*nObj = NULL,
  		*rObj = NULL,
		*tObj = NULL;
  WlzPixelV	lV,
  		uV,
		thrV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzCompThreshType thrType = WLZ_COMPTHRESH_SMOOTHSPLIT;
  
  thrV.type = WLZ_GREY_DOUBLE;
  lV.type = uV.type = WLZ_GREY_INT;
  lV.v.inv = 0;
  uV.v.inv = 255;
  /* Preprocess input object to get normalised smoothed grey object. */
  inGType = WlzGreyTypeFromObj(iObj, &errNum);
  if(inGType == WLZ_GREY_RGBA)
  {
    switch(cCnt)
    {
      case 1:
	nObj = WlzAssignObject(
	       WlzRGBAToChannel(iObj, numC, &errNum), NULL);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzGreyNormalise(nObj, 1);
	}
	break;
      case 2:
	nObj = WlzAssignObject(
	       WlzRGBChanRatio(iObj, numC, denC, WLZ_RGBA_CHANNEL_DUMMY,
			       0, 1, &errNum), NULL);
	break;
      case 3:
	nObj = WlzAssignObject(
	       WlzRGBChanRatio(iObj, numC, denC, mulC,
			       1, 1, &errNum), NULL);
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  else
  {
    nObj = WlzAssignObject(WlzCopyObject(iObj, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzGreyNormalise(nObj, 1);
    }
  }
  if((errNum == WLZ_ERR_NONE) && (smFlg != 0))
  {
    tObj = WlzAssignObject(
           WlzGauss2(nObj, 3, 3, 0, 0, &errNum), NULL);
    WlzFreeObj(nObj);
    nObj = tObj;
    tObj = NULL;
  }
  if((errNum == WLZ_ERR_NONE) && (invFlg != 0))
  {
    errNum = WlzGreyInvertMinMax(nObj, lV, uV);
  }
  if((errNum == WLZ_ERR_NONE) && ((dbgMsk & 2) != 0))
  {
    (void )WlzCelCmpDbgObjOut(nObj, "nObj", prog, dbgId);
  }
  /* Compute threshold value. */
  if(errNum == WLZ_ERR_NONE)
  {
    hObj = WlzAssignObject(
           WlzHistogramObj(nObj, 256, 0, 1, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCompThresholdVT(hObj, thrType, minHSm, maxHSm, tExtra,
    				&thrV, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzValueConvertPixel(&thrV, thrV, WLZ_GREY_DOUBLE);
  }
  if((errNum == WLZ_ERR_NONE) && ((dbgMsk & 1) != 0))
  {
    (void )fprintf(stderr, "%s DBG thrV.v.dbv = %g.\n",
		   prog, thrV.v.dbv);
  }
  /* Threshold the object to get domain, erode to get rid of grot,
   * dilate and fill holes. */
  if(errNum == WLZ_ERR_NONE)
  {
    tObj = WlzAssignObject(
           WlzThreshold(nObj, thrV, WLZ_THRESH_LOW, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dSEObj = WlzAssignObject(
	     WlzMakeSphereObject(WLZ_2D_DOMAINOBJ, dRad,
				   0.0, 0.0, 0.0, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    eSEObj = WlzAssignObject(
	     WlzMakeSphereObject(WLZ_2D_DOMAINOBJ, eRad,
				   0.0, 0.0, 0.0, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    eObj = WlzAssignObject(
    	   WlzStructErosion(tObj, eSEObj, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dObj = WlzAssignObject(
    	   WlzStructDilation(eObj, dSEObj, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzDomainFill(dObj, &errNum);
  }
  (void )WlzFreeObj(dSEObj);
  (void )WlzFreeObj(eObj);
  (void )WlzFreeObj(dObj);
  (void )WlzFreeObj(eSEObj);
  (void )WlzFreeObj(hObj);
  (void )WlzFreeObj(nObj);
  (void )WlzFreeObj(tObj);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New compound array object.
* \ingroup	WlzCelCmpFeatures
* \brief	Splits the given clump object using connectivity and returns
* 		a compound array of objects which are larger than the
* 		given threshold value.
* \param	cmpObj			Given clump object to split.
* \param	minCmpSz		Minimum clump size threshold.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray *WlzCelCmpSplitCmp(WlzObject *cmpObj, int minCmpSz,
				WlzErrorNum *dstErr)
{
  int		idN,
  		idM,
		mass,
		nObj = 0;
  WlzObject	**objs = NULL;
  WlzCompoundArray *cmpCAObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minLn = 4,
  		maxObj = 1024;

  /* Split into seperate objects using connectivity. */
  errNum = WlzLabel(cmpObj, &nObj, &objs, maxObj, minLn, WLZ_8_CONNECTED);
  /* Remove objects which are smaller than the given threshold. */
  if(errNum == WLZ_ERR_NONE)
  {
    idN = idM = 0;
    while(idN < nObj)
    {
      switch(cmpObj->type)
      {
	case WLZ_2D_DOMAINOBJ:
          mass = WlzArea(objs[idN], &errNum);
	  break;
	case WLZ_3D_DOMAINOBJ:
          mass = WlzVolume(objs[idN], &errNum);
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(mass < minCmpSz)
	{
	  (void )WlzFreeObj(objs[idN]);
	  objs[idN] = NULL;
	}
	else
	{
	  objs[idM] = objs[idN];
	  ++idM;
	}
      }
      ++idN;
    }
    nObj = idM;
  }
  /* Create a compoind array of the objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(nObj == 0)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else
    {
      cmpCAObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, nObj, objs,
                                      objs[0]->type, &errNum);
    }
  }
  if(objs != NULL)
  {
    for(idN = 0; idN < nObj; ++idN)
    {
      (void )WlzFreeObj(objs[idN]);
    }
    AlcFree(objs);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cmpCAObj);
}

/*!
* \return	New compound array object.
* \ingroup	WlzCelCmpFeatures
* \brief	Creates a new compoind object of cells objects, each
* 		of which is paired with a clump object. Cell objects
* 		may be NULL if a clump object does not intersect the
* 		cells object.
* \param	celObj			Given cells object.
* \param	cmpCAObj		Given clump compound array.
* \param	minCelSz		Minimum threshold for cell size.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCompoundArray *WlzCelCmpSplitCel(WlzObject *celObj,
				WlzCompoundArray *cmpCAObj,
				int minCelSz,
				WlzErrorNum *dstErr)
{
  int		idN,
		mass;
  WlzObject	**objs = NULL;
  WlzCompoundArray *celCAObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((objs = AlcCalloc(cmpCAObj->n, sizeof(WlzObject *))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idN = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < cmpCAObj->n))
    {
      mass = 0;
      objs[idN] = WlzIntersect2(celObj, cmpCAObj->o[idN], &errNum);
      if((errNum == WLZ_ERR_NONE) && (objs[idN] != NULL))
      {
	switch(objs[idN]->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	    mass = WlzArea(objs[idN], &errNum);
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    mass = WlzVolume(objs[idN], &errNum);
	    break;
	  case WLZ_EMPTY_OBJ:
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
      if((errNum == WLZ_ERR_NONE) && (mass <= 0))
      {
        (void )WlzFreeObj(objs[idN]);
        objs[idN] = NULL;
      }
      ++idN;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    celCAObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, cmpCAObj->n,
    				    objs, cmpCAObj->o[0]->type, &errNum);
    
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if((objs != NULL) && (cmpCAObj != NULL))
    {
      for(idN = 0; idN < cmpCAObj->n; ++idN)
      {
	(void )WlzFreeObj(objs[idN]);
      }
    }
    else
    {
      (void )WlzFreeObj((WlzObject *)celCAObj);
    }
  }
  AlcFree(objs);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(celCAObj);
}

/*!
* \return	2D double array of features.
* \ingroup	WlzCelCmpFeatures
* \brief	Computes numerical features of the segmented clump and
* 		cell objects. The ojects are paired and may be NULL.
* 		The features computed are clump mass, cell mass,
* 		cell centrality wrt the cells.
* 		The 2d array should be freed using AlcDouble2Free()
* 		and the label strings should each be freed before
* 		freeing the label array using AlcFree().
* \param	cmpCAObj		Clump objects.
* \param	celCAObj		Cell objects.
* \param	dstLabels		Destination pointer for feature
* 					label strings. May be NULL.
* \param	dstFeatSz		Destination pointer for the size
* 					of the features array. Must not
* 					be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static double	**WlzCelCmpCompFeat(WlzCompoundArray *cmpCAObj,
				WlzCompoundArray *celCAObj,
				char ***dstLabels,
				WlzIVertex2 *dstFeatSz,
				int nRay, int binFlg,
				WlzErrorNum *dstErr)
{
  int		idX,
  		idY;
  double	maxR = 0.0,
  		tD;
  double	**feat = NULL;
  WlzIVertex2	featSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		**labels = NULL;

  featSz.vtX = 4;
  featSz.vtY = cmpCAObj->n;
  if(dstLabels)
  {
    if(((labels = AlcCalloc(featSz.vtX, sizeof(char *))) == NULL) ||
       ((labels[0] = AlcStrDup("clump mass")) == NULL) ||
       ((labels[1] = AlcStrDup("cell mass")) == NULL) ||
       ((labels[2] = AlcStrDup("centrality")) == NULL) ||
       ((labels[3] = AlcStrDup("circularity")) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) &&
     (AlcDouble2Malloc(&feat, featSz.vtY, featSz.vtX) != ALC_ER_NONE))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  idY = 0;
  while((errNum == WLZ_ERR_NONE) && (idY < featSz.vtY))
  {
    for(idX = 0; idX < featSz.vtX; ++idX)
    {
      feat[idY][idX] = 0.0;
    }
    if(cmpCAObj->o[idY] != NULL)
    {
      switch(cmpCAObj->o[idY]->type)
      {
	case WLZ_2D_DOMAINOBJ:
	  feat[idY][0] = WlzArea(cmpCAObj->o[idY], &errNum);
	  break;
	case WLZ_3D_DOMAINOBJ:
	  feat[idY][0] = WlzVolume(cmpCAObj->o[idY], &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(celCAObj->o[idY] != NULL)
      {
	switch(celCAObj->o[idY]->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	    feat[idY][1] = WlzArea(celCAObj->o[idY], &errNum);
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    feat[idY][1] = WlzVolume(celCAObj->o[idY], &errNum);
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if((feat[idY][0] > 0) && (feat[idY][1] > 0))
      {
        feat[idY][2] = WlzCentrality(celCAObj->o[idY], cmpCAObj->o[idY],
				     nRay, binFlg, &maxR, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
	  tD = ALG_M_PI * maxR * maxR;
	  feat[idY][3] = (tD > DBL_EPSILON)? feat[idY][0] / tD: 0.0;
	}
      }
    }
    ++idY;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstLabels = labels;
    *dstFeatSz = featSz;
  }
  else
  {
    if(labels)
    {
      for(idX = 0; idX < featSz.vtX; ++idX)
      {
        (void )AlcFree(labels[idX]);
      }
    }
    (void )AlcDouble2Free(feat);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(feat);
}

/*!
* \return	New Woolz object.
* \ingroup	WlzCelCmpFeatures
* \brief	Creates and returns a new Woolz object in which the given
* 		object is combined with washes for the clump and mask
* 		domains.
* 		Returned object has given object green, clump red and
* 		cells blue.
* \param	gObj			Object with grey valued image.
* \param	cmpObj			Object with clump domain.
* \param	celObj			Object with cell domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCelCmpColorise(WlzObject *gObj,
				  WlzObject *cmpObj,
				  WlzObject *celObj,
				  WlzErrorNum *dstErr)
{
  int		idN;
  WlzObjectType	gTType;
  WlzObject	*cObj = NULL,
  		*tObj0 = NULL,
		*tObj1 = NULL;
  WlzObject	*objA[4];
  WlzValues	tVal;
  WlzCompoundArray *cpd = NULL;
  WlzPixelV	pV0,
  		pV255;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  tVal.core = NULL;
  pV0.type = pV255.type = WLZ_GREY_UBYTE;
  pV0.v.ubv = 0;
  pV255.v.ubv = 255;
  objA[0] = objA[1] = objA[2] = objA[3] = NULL;
  /* Make sure given object is 2D, 3D not supported yet. */
  if(gObj->type != WLZ_2D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  /* Set given image array object. */
  if(errNum == WLZ_ERR_NONE)
  {
    /* Green */
    objA[1] = WlzAssignObject(gObj, NULL);
    gTType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_UBYTE, &errNum);
  }
  /* Create clump array object. */
  if(errNum == WLZ_ERR_NONE)
  {
    /* Red */
    objA[0] = WlzAssignObject(
              WlzCopyObject(gObj, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((cmpObj != NULL) && (cmpObj->type != WLZ_EMPTY_OBJ))
    {
      if((celObj == NULL) || (celObj->type == WLZ_EMPTY_OBJ))
      {
        tObj0 = WlzAssignObject(cmpObj, NULL);
      }
      else
      {
	tObj0 = WlzAssignObject(
		WlzDiffDomain(cmpObj, celObj, &errNum), NULL);
      }
      if((errNum == WLZ_ERR_NONE) &&
	 (tObj0 != NULL) && (tObj0->type != WLZ_EMPTY_OBJ))
      {
	tObj1 = WlzAssignObject(
		WlzMakeMain(gObj->type, tObj0->domain, objA[0]->values,
			     NULL, NULL, &errNum), NULL);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzGreySetValue(tObj1, pV255);
	}
	(void )WlzFreeObj(tObj1); tObj1 = NULL;
      }
      (void )WlzFreeObj(tObj0); tObj0 = NULL;
    }
  }
  /* Create cell mask object. */
  if(errNum == WLZ_ERR_NONE)
  {
    /* Blue */
    objA[2] = WlzAssignObject(
              WlzCopyObject(gObj, &errNum), NULL);
  }
  if((errNum == WLZ_ERR_NONE) &&
     (celObj != NULL) && (celObj->type != WLZ_EMPTY_OBJ))
  {
    tObj0 = WlzMakeMain(gObj->type, celObj->domain, objA[2]->values,
		         NULL, NULL, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzGreySetValue(tObj0, pV255);
    }
    (void )WlzFreeObj(tObj0); tObj0 = NULL;
  }
  /* Create alpha mask object. */
  if(errNum == WLZ_ERR_NONE)
  {
    tVal.v = WlzNewValueTb(gObj, gTType, pV0, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Alpha */
    objA[3] = WlzAssignObject(
              WlzMakeMain(gObj->type, gObj->domain, tVal,
                          NULL, NULL, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetValue(objA[3], pV255);
  }
  /* Make compound object. */
  if(errNum == WLZ_ERR_NONE)
  {
    cpd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_2, 3, 4, objA,
                               gObj->type, &errNum);
  }
  for(idN = 0; idN < 4; ++idN)
  {
    (void )WlzFreeObj(objA[idN]);
  }
  /* Create RGB object. */
  if(errNum == WLZ_ERR_NONE)
  {
    cObj = WlzCompoundToRGBA(cpd, WLZ_RGBA_SPACE_RGB, &errNum);
  }
  (void )WlzFreeObj((WlzObject *)cpd);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cObj);
}

/*!
* \return	New Normalised object.
* \ingroup	WlzCelCmpFeatures
* \brief	Normalises the given object's grey values. If the object
* 		has RGBA grey values then the red, green and blue channels
* 		are independently normalised (the alpha channel is not
* 		changed).
* \param	gObj			Given Object.
* \param	dither			Non-zero for grey value dithering.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCelCmpNormalise(WlzObject *gObj, int dither,
				     WlzErrorNum *dstErr)
{
  int		idx;
  WlzGreyType 	gType;
  WlzObject	*nObj = NULL;
  WlzCompoundArray *cpdObj = NULL;;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	gType = WlzGreyTypeFromObj(gObj, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  if(gType == WLZ_GREY_RGBA)
	  {
	    idx = 0;
	    cpdObj = WlzRGBAToCompound(gObj, WLZ_RGBA_SPACE_RGB, &errNum);
	    while((errNum == WLZ_ERR_NONE) && (idx < 3))
	    {
	      errNum = WlzGreyNormalise(cpdObj->o[idx], dither);
	      ++idx;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      nObj = WlzCompoundToRGBA(cpdObj, WLZ_RGBA_SPACE_RGB, &errNum);
	    }
	    (void )WlzFreeObj((WlzObject *)cpdObj);
	  }
	  else
	  {
	    nObj = WlzCopyObject(gObj, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzGreyNormalise(nObj, dither);
	    }
	  }
	}
        break;
      default:
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nObj);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
