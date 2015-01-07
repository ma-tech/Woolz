#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRadialDistribution_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzRadialDistribution.c
* \author       Bill Hill
* \date         January 2015
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2015],
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
* \brief	Outputs the radial distribution of components of
* 		the input image.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref WlzRadialDistribution "WlzRadialDistribution"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzradialdistribution WlzRadialDistribution
\par Name
WlzRadialDistribution - outputs the radial distribution of components of
		        the input image.
\par Synopsis
\verbatim
WlzRadialDistribution [-h] [-v] [-A] [-D] [-G] [-H] [-E] [-L] [-R]
                      [-c #,#] [-d <debug image>] [-n #]  [-o <out file>]
		      [-t #] [<input image>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Be verbose output.</td>
  </tr>
  <tr>
    <td><b>-A</td>
    <td>Sort output by area (default).</td>
  </tr>
  <tr>
    <td><b>-D</td>
    <td>Sort output by distance from boundary.</td>
  </tr>
  <tr>
    <td><b>-G</td>
    <td>Sort output by angle.</td>
  </tr>
  <tr>
    <td><b>-R</td>
    <td>Sort output by radial distance from centre.</td>
  </tr>
  <tr>
    <td><b>-H</td>
    <td>Threshold high, use pixels at or above threshold (default).</td>
  </tr>
  <tr>
    <td><b>-E</td>
    <td>Threshold equal, use pixels at threshold.</td>
  </tr>
  <tr>
    <td><b>-L</td>
    <td>Threshold low, use pixels below threshold.</td>
  </tr>
  <tr>
    <td><b>-c</td>
    <td>Centre (default is image centre).</td>
  </tr>
  <tr>
    <td><b>-d</td>
    <td>Debug image.</td>
  </tr>
  <tr>
    <td><b>-n</td>
    <td>Minimum area (default %g).</td>
    </tr>
  <tr>
    <td><b>-t</td>
    <td>Threshold value (default is to compute using Otsu's method).</td>
  </tr>
</table>
\par Description
Segments the given object using a threshold value and outputs the
radial distribution of the thresholded components.
By default the input image object is read from the standard input and
the radial distribution is written to the standard output.
The image formats understood include wlz, jpg and tif.
The output format is:
\verbatim
  <angle> <dist from centre> <area> <x pos>,<y pos> <dist form boundary>
\endverbatim
\par Examples
\verbatim
WlzRadialDistribution -o out.txt -d debug.jpg in.tif
\endverbatim
The input image is read from in.tif, a debug image showing the
segmented regions is written to debug.jpg and the radial distribution
statistics are written to the file out.txt. With the output in
out.txt, the following R code would plot the data as a set of circles
with radius proportional to the square root of the component area:
\verbatim
  data <- read.table("out.txt")
  attach(data)
  symbols(x=data$V1, y=data$V2, circles=sqrt(data$V3))
\endverbatim
\par File
\ref WlzRadialDistribution.c "WlzRadialDistribution.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzLabel "WlzLabel(1)"
\ref WlzThreshold "WlzThreshold(1)"
\ref WlzExtFFConvert "WlzExtFFConvert(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

typedef struct _WlzRadDistRec
{
  double	angle;
  double	radius;
  double	area;
  WlzDVertex2	pos;
  double	dist;
} WlzRadDistRec;

typedef enum _WlzRadDistVal
{
  WLZ_RADDISTVAL_ANGLE,
  WLZ_RADDISTVAL_RADIUS,
  WLZ_RADDISTVAL_AREA,
  WLZ_RADDISTVAL_DIST
} WlzRadDistVal;

static int      		WlzRadDistRecSortAngle(
				  const void *p0,
				  const void *p1);
static int      		WlzRadDistRecSortArea(
                                  const void *p0,
				  const void *p1);
static int      		WlzRadDistRecSortRadius(
                                  const void *p0,
				  const void *p1);
static int      		WlzRadDistRecSortDist(
                                  const void *p0,
				  const void *p1);
static int 			WlzRadDistParsePath(
				  char *path,
				  char **dstDir,
				  char **dstFile,
				  char **dstExt,
				  WlzEffFormat *dstFmt);

int             main(int argc, char *argv[])
{
  int		option,
		nReg = 0,
		tNReg = 0,
  		ok = 1,
		usage = 0,
		verbose = 0,
		threshSet = 0,
		centreSet = 0;
  double	minArea = 2;
  char		*inExt,
		*dbgExt,
		*inDir,
		*dbgDir,
		*inFile,
		*dbgFile,
  		*inPath = NULL,
		*dbgPath = NULL,
		*outFile = NULL;
  WlzRadDistVal distSort = WLZ_RADDISTVAL_AREA;
  WlzRadDistRec	*distData = NULL;
  WlzPixelV	thrVal;
  WlzDVertex2	centre;
  WlzCompThreshType thrMtd = WLZ_COMPTHRESH_OTSU;
  WlzThresholdType thrMod = WLZ_THRESH_HIGH;
  WlzEffFormat	inFmt = WLZEFF_FORMAT_NONE,
  		dbgFmt = WLZEFF_FORMAT_NONE;
  WlzObject	*inObj = NULL,
		*disObj = NULL,
  		*segObj = NULL;
  WlzGreyValueWSpace *disGVWSp = NULL;
  WlzObject	**regObjs = NULL;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	maxObj = 1000000;
  char		pathBuf[FILENAME_MAX];
  const double	eps = 1.0e-06;
  const char	*errMsg;
  static char	optList[] = "hvAGDHELR:c:d:n:o:t:",
		defFile[] = "-";

  thrVal.type = WLZ_GREY_DOUBLE;
  thrVal.v.dbv = 0.0;
  outFile = defFile;
  while((usage == 0) && ok &&
        ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'A':
        distSort = WLZ_RADDISTVAL_AREA;
	break;
      case 'D':
        distSort = WLZ_RADDISTVAL_DIST;
	break;
      case 'G':
        distSort = WLZ_RADDISTVAL_ANGLE;
	break;
      case 'H':
        thrMod = WLZ_THRESH_HIGH;
	break;
      case 'E':
        thrMod = WLZ_THRESH_EQUAL;
	break;
      case 'L':
        thrMod = WLZ_THRESH_LOW;
	break;
      case 'R':
        distSort = WLZ_RADDISTVAL_RADIUS;
	break;
      case 'h':
        usage = 1;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'c':
	centreSet = 1;
        if(sscanf(optarg, "%lg,%lg", &(centre.vtX), &(centre.vtY)) != 2)
	{
	  usage = 1;
	}
        break;
      case 'd':
        dbgPath = optarg;
	break;
      case 'o':
        outFile = optarg;
	break;
      case 'n':
        if(sscanf(optarg, "%lg", &minArea) != 1)
	{
	  usage = 1;
	}
	break;
      case 't':
	threshSet = 1;
        if(sscanf(optarg, "%lg", &(thrVal.v.dbv)) != 1)
	{
	  usage = 1;
	}
	break;
      default:
        usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok)
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      inPath = *(argv + optind);
    }
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "inPath = %s\n", inPath);
  }
  /* Parse input file path into path + name + ext. */
  if(ok)
  {
    ok = (usage = WlzRadDistParsePath(inPath, &inDir, &inFile, &inExt,
                                      &inFmt)) == 0;
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "inDir = %s\n", inDir);
    (void )fprintf(stderr, "inFile = %s\n", inFile);
    (void )fprintf(stderr, "inExt = %s\n", (inExt)? inExt: "(null)");
    (void )fprintf(stderr, "inFmt = %s\n",
    		   WlzEffStringFromFormat(inFmt, NULL));
  }
  /* Read image. */
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if(inExt)
    {
      (void )sprintf(pathBuf, "%s/%s.%s", inDir, inFile, inExt);
    }
    else
    {
      (void )sprintf(pathBuf, "%s/%s", inDir, inFile);
    }
    if(((inObj = WlzAssignObject(WlzEffReadObj(NULL, pathBuf, inFmt,
    					       0, 0, 0,
					       &errNum), NULL)) == NULL) ||
       (inObj->type != WLZ_2D_DOMAINOBJ))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: Failed to read 2D image object from file %s (%s)\n",
		     *argv, pathBuf, errMsg);
    }
  }
  if(ok && verbose)
  {
    (void )fprintf(stderr, "read input image ok.\n");
  }
  /* Convert to grey if needed, normalise 0 - 255 if needed and compute
   * threshold value unless already known. */
  if(ok)
  {
    if(WlzGreyTypeFromObj(inObj, NULL) == WLZ_GREY_RGBA)
    {
      WlzObject *ppObj;

      ppObj = WlzAssignObject(
	      WlzRGBAToModulus(inObj, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzFreeObj(inObj);
	inObj = ppObj;
      }
    }
    if(threshSet == 0)
    {
      WlzObject *hObj = NULL;

      errNum = WlzGreyNormalise(inObj, 1);
      if(errNum == WLZ_ERR_NONE)
      {
        hObj = WlzHistogramObj(inObj, 256, 0.0, 1.0, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	threshSet = 1;
        errNum = WlzCompThreshold(&thrVal.v.dbv, hObj, thrMtd, 0);
      }
      (void )WlzFreeObj(hObj);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to normalise object (%s)\n",
		     *argv, errMsg);
    }
  }
  /* Segment the object. */
  if(ok)
  {
    if(inObj->values.core == NULL)
    {
      segObj = WlzAssignObject(inObj, NULL);
    }
    else
    {
      segObj = WlzAssignObject(
               WlzThreshold(inObj, thrVal, thrMod, &errNum), NULL);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr, "%s: failed to segment image (%s)\n",
		       *argv, errMsg);
      }
    }
  }
  /* Compute object with the same domain as the input object but in which
   * the values are the minimum distance from an edge. */
  if(ok)
  {
    WlzObject	*bObj = NULL;

    bObj = WlzBoundaryDomain(inObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      disObj = WlzAssignObject(       
               WlzDistanceTransform(inObj, bObj, WLZ_OCTAGONAL_DISTANCE,
	       			    0.0, 0.0, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      disGVWSp = WlzGreyValueMakeWSp(disObj, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to compute distance object (%s)\n",
		     *argv, errMsg);
    }
    (void )WlzFreeObj(bObj);
  }
  /* Output the debug image if required. */
  if(ok && dbgPath)
  {
    WlzObject	*dbgObj;

    dbgObj = WlzAssignObject(WlzCopyObject(inObj, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzPixelV	iMin,
		iMax,
		oMin,
		oMax;

      if(dbgObj->values.core == NULL)
      {
        WlzValues tmpVal;

	oMax.type = WLZ_GREY_UBYTE;
	oMax.v.ubv = 255;
	tmpVal.v = WlzNewValueTb(dbgObj,
				 WlzGreyTableType(WLZ_GREY_TAB_RAGR,
				                  WLZ_GREY_UBYTE, NULL),
	                         oMax, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
	  dbgObj->values = WlzAssignValues(tmpVal, NULL);
	}
      }
      else
      {
        WlzObject *tmpObj = NULL;

	oMin.type = WLZ_GREY_UBYTE;
	oMin.v.ubv = 0;
	oMax.type = WLZ_GREY_UBYTE;
	oMax.v.ubv = 200;
	errNum = WlzGreyRange(dbgObj, &iMin, &iMax);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzGreySetRange(dbgObj, iMin, iMax, oMin, oMax, 0);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  tmpObj = WlzMakeMain(inObj->type, segObj->domain, dbgObj->values,
	                       NULL, NULL, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  oMax.v.ubv = 255;
	  errNum = WlzGreySetValue(tmpObj, oMax);
	}
	(void )WlzFreeObj(tmpObj);
	if(errNum == WLZ_ERR_NONE)
	{
	  tmpObj = WlzConvertPix(dbgObj, WLZ_GREY_UBYTE, &errNum);
	  (void )WlzFreeObj(dbgObj);
	  dbgObj = WlzAssignObject(tmpObj, NULL);
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzRadDistParsePath(dbgPath, &dbgDir, &dbgFile, &dbgExt,
      			         &dbgFmt);
      if(dbgExt)
      {
	(void )sprintf(pathBuf, "%s/%s.%s", dbgDir, dbgFile, dbgExt);
      }
      else
      {
	(void )sprintf(pathBuf, "%s/%s", dbgDir, dbgFile);
      }
      errNum = WlzEffWriteObj(NULL, pathBuf, dbgObj, dbgFmt);
    }
    (void )WlzFreeObj(dbgObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to output the debug image (%s)\n",
		     *argv, errMsg);
    }
  }
  /* Label the segmented object. */
  if(ok)
  {
    errNum = WlzLabel(segObj, &nReg, &regObjs, maxObj, 0, WLZ_8_CONNECTED);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to split into components (%s)\n",
		     *argv, errMsg);
    }
    if(ok && verbose)
    {
      (void )fprintf(stderr, "nReg = %d\n", nReg);
    }
  }
  /* Compute centre of mass if not known. */
  if(ok)
  {
    if(centreSet == 0)                          
    {
      centre = WlzCentreOfMass2D(inObj, 1, NULL, &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr, "%s: failed to compute centre of mass (%s)\n",
		       *argv, errMsg);
      }
    }
    if(ok && verbose)
    {
      (void )fprintf(stderr, "centre = %lg,%lg\n", centre.vtX, centre.vtY);
    }
  }
  /* Allocate a radial distribution table. */
  if(ok)
  {
    if((distData = (WlzRadDistRec *)
                   AlcCalloc(nReg, sizeof(WlzRadDistRec))) == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to allocate result lable (%s)\n",
		     *argv, errMsg);
    }
    
  }
  /* Compute the redial distribution data. */
  if(ok)
  {
    int		idR = 0,
    		idS = 0;

    while((errNum == WLZ_ERR_NONE) && (idR < nReg))
    {
      double	mass;
      WlzDVertex2 com;

      com = WlzCentreOfMass2D(regObjs[idR], 1, &mass, NULL);
      if(mass > minArea - eps)
      {
	WlzGreyValueGet(disGVWSp, 0.0, com.vtY, com.vtX);
	distData[idS].pos = com;
	distData[idS].area = mass;
	WLZ_VTX_2_SUB(com, centre, com);
	distData[idS].radius = WLZ_VTX_2_LENGTH(com);
	distData[idS].angle = ALG_M_PI + atan2(com.vtY, com.vtX);
	switch(disGVWSp->gType)
	{
	  case WLZ_GREY_LONG:
	    distData[idS].dist = *(disGVWSp->gPtr[0].lnp);
	    break;
	  case WLZ_GREY_INT:
	    distData[idS].dist = *(disGVWSp->gPtr[0].inp);
	    break;
	  case WLZ_GREY_SHORT:
	    distData[idS].dist = *(disGVWSp->gPtr[0].shp);
	    break;
	  case WLZ_GREY_UBYTE:
	    distData[idS].dist = *(disGVWSp->gPtr[0].ubp);
	    break;
	  case WLZ_GREY_FLOAT:
	    distData[idS].dist = *(disGVWSp->gPtr[0].flp);
	    break;
	  case WLZ_GREY_DOUBLE:
	    distData[idS].dist = *(disGVWSp->gPtr[0].dbp);
	    break;
	  default:
	    distData[idS].dist = 0.0;
	    break;
	}
	++idS;
      }
      ++idR;
    }
    tNReg = idS;
    switch(distSort)
    {
      case WLZ_RADDISTVAL_AREA:
        (void )qsort(distData, tNReg, sizeof(WlzRadDistRec),
		     WlzRadDistRecSortArea);
	break;
      case WLZ_RADDISTVAL_ANGLE:
        (void )qsort(distData, tNReg, sizeof(WlzRadDistRec), 
		     WlzRadDistRecSortAngle);
	break;
      case WLZ_RADDISTVAL_RADIUS:
        (void )qsort(distData, tNReg, sizeof(WlzRadDistRec),
		     WlzRadDistRecSortRadius);
	break;
      case WLZ_RADDISTVAL_DIST:
        (void )qsort(distData, tNReg, sizeof(WlzRadDistRec),
		     WlzRadDistRecSortDist);
	break;
    }
  }
  /* Output the sorted radial distribution table. */
  if(ok)
  {
    if(((fP = strcmp(outFile, "-")?
              fopen(outFile, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: failed to open output file %s\n",
                     *argv, outFile);
    }
  }
  if(ok)
  {
    int		idR;

    for(idR = 0; idR < tNReg; ++idR)
    {
      double a;

      a = (distData[idR].angle > 0.0)?
	  0   + (180 * distData[idR].angle / ALG_M_PI):
          360 + (180 * distData[idR].angle / ALG_M_PI);
      (void )fprintf(fP, "%g %g %g %g,%g %g\n",
		     a,
                     distData[idR].radius,
		     distData[idR].area,
		     distData[idR].pos.vtX,
		     distData[idR].pos.vtY,
		     distData[idR].dist);
    }
  }
  if(strcmp(outFile, "-"))
  {
    (void )fclose(fP);
  }
  /* Tidy up. */
  AlcFree(distData);
  WlzGreyValueFreeWSp(disGVWSp);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(disObj);
  (void )WlzFreeObj(segObj);
  if(regObjs)
  {
    int		idR;

    for(idR = 0; idR < nReg; ++idR)
    {
      (void )WlzFreeObj(regObjs[idR]);
    }
    AlcFree(regObjs);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-v] [-A] [-D] [-G] [-H] [-E] [-L] [-R]\n"
    "\t\t[-c #,#] [-d <debug image>] [-n #]  [-o <out file>]\n"
    "\t\t[-t #] [<input image>]\n"
    "Segments the given object using a threshold value and outputs the \n"
    "radial distribution of the thresholded components.\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Help - prints this usage masseage.\n"
    "  -v  Verbose output.\n"
    "  -A  Sort output by area (default).\n"
    "  -D  Sort output by distance from boundary.\n"
    "  -G  Sort output by angle.\n"
    "  -H  Threshold high, use pixels at or above threshold (default).\n"
    "  -E  Threshold equal, use pixels at threshold.\n"
    "  -L  Threshold low, use pixels below threshold.\n"
    "  -R  Sort output by radial distance from centre.\n"
    "  -c  Centre (default is image centre).\n"
    "  -d  Debug image.\n"
    "  -n  Minimum area (default %g).\n"
    "  -t  Threshold value (default is to compute using Otsu's method).\n"
    "By default the input image object is read from the standard input and\n"
    "the radial distribution is written to the standard output.\n"
    "The image formats understood include wlz, jpg and tif.\n"
    "The output format is:\n"
    "  <angle> <dist from centre> <area> <x pos>,<y pos> <dist form boundary>\n"
    "Example:\n"
    "  %s -o out.txt -d debug.jpg in.tif\n"
    "The input image is read from in.tif, a debug image showing the\n"
    "segmented regions is written to debug.jpg and the radial distribution\n"
    "statistics are written to the file out.txt. With the output in\n"
    "out.txt, the following R code would plot the data as a set of circles\n"
    "with radius proportional to the square root of the component area:\n"
    "  data <- read.table(\"out.txt\")\n"
    "  attach(data)\n"
    "  symbols(x=data$V1, y=data$V2, circles=sqrt(data$V3))\n",
    argv[0],
    WlzVersion(),
    minArea,
    argv[0]);
  }
  return(!ok);
}

/*!
* \return	Difference between angle of records.
* \brief	Comparison function for qsort() to sort radial distribution
* 		records.
* \param	p0			Pointer to first record.
* \param	p1			Pointer to second record.
*/
static int      		WlzRadDistRecSortAngle(
				  const void *p0,
				  const void *p1)
{
  int		cmp;
  WlzRadDistRec *r0,
  		*r1;
  
  r0 = (WlzRadDistRec *)p0;
  r1 = (WlzRadDistRec *)p1;
  cmp = (r1->angle > r0->angle);
  return(cmp);
}

/*!
* \return	Difference between angle of records.
* \brief	Comparison function for qsort() to sort radial distribution
* 		records.
* \param	p0			Pointer to first record.
* \param	p1			Pointer to second record.
*/
static int      		WlzRadDistRecSortArea(
				  const void *p0,
				  const void *p1)
{
  int		cmp;
  WlzRadDistRec *r0,
  		*r1;
  
  r0 = (WlzRadDistRec *)p0;
  r1 = (WlzRadDistRec *)p1;
  cmp = (r1->area > r0->area);
  return(cmp);
}

/*!
* \return	Difference between angle of records.
* \brief	Comparison function for qsort() to sort radial distribution
* 		records.
* \param	p0			Pointer to first record.
* \param	p1			Pointer to second record.
*/
static int      		WlzRadDistRecSortRadius(
				  const void *p0,
				  const void *p1)
{
  int		cmp;
  WlzRadDistRec *r0,
  		*r1;
  
  r0 = (WlzRadDistRec *)p0;
  r1 = (WlzRadDistRec *)p1;
  cmp = (r1->radius > r0->radius);
  return(cmp);
}

/*!
* \return	Difference between distance of records.
* \brief	Comparison function for qsort() to sort radial distribution
* 		records.
* \param	p0			Pointer to first record.
* \param	p1			Pointer to second record.
*/
static int      		WlzRadDistRecSortDist(
				  const void *p0,
				  const void *p1)
{
  int		cmp;
  WlzRadDistRec *r0,
  		*r1;
  
  r0 = (WlzRadDistRec *)p0;
  r1 = (WlzRadDistRec *)p1;
  cmp = (r1->dist > r0->dist);
  return(cmp);
}

static int 			WlzRadDistParsePath(
				  char *path,
				  char **dstDir,
				  char **dstFile,
				  char **dstExt,
				  WlzEffFormat *dstFmt)
{
  int		len,
  		usage = 0;
  char		*dir,
  		*file,
		*ext;
  WlzEffFormat  fmt = WLZEFF_FORMAT_NONE;

  if(*dstFmt)
  {
    fmt = *dstFmt;
  }
  if(((len = strlen(path)) < 3) || /* Minimum of 3 chars. */
      (*(path + len - 1) == '/'))  /* Directory not plain file. */
  {
    usage = 1;
  }
  if(!usage)
  {
    char	*sep,
    		*dot;

    sep = strrchr(path, '/');
    dot = strrchr(path, '.');
    if(sep == NULL)
    {
      dir = ".";
      file = path;
    }
    else
    {
      *sep = '\0';
      dir = path;
      file = sep + 1;
    }
    if(dot)
    {
      *dot = '\0';
      ext = dot + 1;
    }
    else
    {
      ext = NULL;
    }
  }
  if(!usage)
  {
    if(fmt == WLZEFF_FORMAT_NONE)
    {
      char	extBuf[16];

      strncpy(extBuf, ext, 15); extBuf[15] = '\0';
      (void )WlzStringToLower(extBuf);
      fmt = WlzEffStringExtToFormat(extBuf);
      if(fmt == WLZEFF_FORMAT_NONE)
      {
        usage = 1;
      }
    }
  }
  if(!usage)
  {
    if(dstDir)
    {
      *dstDir = dir;
    }
    if(dstFile)
    {
      *dstFile = file;
    }
    if(dstExt)
    {
      *dstExt = ext;
    }
    if(dstFmt)
    {
      *dstFmt = fmt;
    }
  }
  return(usage);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
