#pragma ident "MRC HGU $Id$"
/*!
\ingroup      BinWlz
\defgroup     wlzcontourfrompoints WlzContourFromPoints
\par Name
WlzContourFromPoints - Compute contour from point sets.
\par Synopsis
\verbatim
WlzContourFromPoints [-h] [-2] [-3] [-o<out object>]
                     [-S<surface points file>]
                     [-I<inside points file>] [-O<outside points file>]
                     [-a#] [-b#] [-c#] [-d#] [-e#] [-D#] [-M#] [-T#]
                     [-X#]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-2</b></td>
    <td> Points are 2D and contour is a curve in 2D. </td>
  </tr>
  <tr>
    <td><b>-3</b></td>
    <td> Points are 3D and contour is a surface in 3D. </td>
  </tr>
  <tr>
    <td><b>-S</b></td>
    <td> On surface points file. </td>
  </tr>
  <tr>
    <td><b>-I</b></td>
    <td>Inside points file. </td>
  </tr>
  <tr>
    <td><b>-O</b></td>
    <td>Outside points file. </td>
  </tr>
  <tr>
    <td><b>-a</b></td>
    <td>Surface points alpha value. </td>
  </tr>
  <tr>
    <td><b>-b</b></td>
    <td>Inside points alpha value. </td>
  </tr>
  <tr>
    <td><b>-c</b></td>
    <td>Outside points alpha value. </td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Inside points distance value. </td>
  </tr>
  <tr>
    <td><b>-e</b></td>
    <td>Outside points distance value. </td>
  </tr>
  <tr>
    <td><b>-D</b></td>
    <td>Delta multiorder spline soothing value. </td>
  </tr>
  <tr>
    <td><b>-M</b></td>
    <td>Distance object sampling factor (can be used to reduce the
      number of surface facets.) </td>
  </tr>
  <tr>
    <td><b>-T</b></td>
    <td>Tau multiorder spline soothing value. </td>
  </tr>
  <tr>
    <td><b>-X</b></td>
    <td>Plane for 2D contour from 3D point sets. </td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>
Unless specified, output to standard out.

\par Description
Computes a contour from three point sets: inside; outside and on a
surface.
All points are given using the following ascii format:
\par
  \<vertex x\> \<vertex y\> [\<vertex z\>]
\par Examples
\verbatim
\endverbatim

\par See Also
\par Bugs
None known
\author       richard <Richard.Baldock@hgu.mrc.ac.uk>
\date         Mon Aug  1 08:39:03 2005
\version      MRC HGU $Id$
              $Revision$
              $Name$
\par Copyright:
             1994-2003 Medical Research Council, UK.
              All rights reserved.
\par Address:
              MRC Human Genetics Unit,
              Western General Hospital,
              Edinburgh, EH4 2XU, UK.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*!
* \file         binWlz/WlzContourFromPoints.c
* \author       Bill Hill
* \date         February 2003
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Woolz filter which computes a contour from three
*		point sets: inside; outside and on a surface.
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

#define WLZ_CFP_READLN_LEN	1024

static WlzVertexP 		WlzCFPReadVtxArray(
				  FILE *fP,
				  int dim,
				  int *nIPts,
				  WlzErrorNum *dstErr);

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		dim = 3,
  		nSPts,
  		nIPts,
		nOPts,
		option,
		planeFlg = 0,
		planeIdx = 0,
		ok = 1,
		usage = 0;
  double	delta,
  		tau,
		sAlpha,
		iAlpha,
		oAlpha,
		iDist,
		oDist,
		samFac = 1.0;
  WlzObject 	*dObj = NULL;
  WlzBox	bBox0,
  		bBox1;
  WlzIVertex2	tDASz;
  WlzVertexP	sPts,
  		iPts,
		oPts;
  WlzPixelV	bgdV;
  WlzDomain	outDom;
  WlzValues	dumVal;
  WlzObject	*outObj = NULL;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*sFileStr,
  		*iFileStr,
  		*oFileStr,
		*outObjFileStr;
  const char    *errMsg;
  static char	optList[] = "h23o:S:I:O:a:b:c:d:e:D:M:T:X:",
		defFileStr[] = "-";

  bgdV.type = WLZ_GREY_DOUBLE;
  outDom.core = NULL;
  dumVal.core = NULL;
  sPts.v = iPts.v = oPts.v = NULL;
  /* Parse the command line. */
  opterr = 0;
  sFileStr = iFileStr = oFileStr = outObjFileStr = defFileStr;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'S':
        sFileStr = optarg;
	break;
      case 'I':
        iFileStr = optarg;
	break;
      case 'O':
        oFileStr = optarg;
	break;
      case 'a':
        if(sscanf(optarg, "%lg", &sAlpha) != 1)
	{
	  usage = 1;
	}
	break;
      case 'b':
        if(sscanf(optarg, "%lg", &iAlpha) != 1)
	{
	  usage = 1;
	}
	break;
      case 'c':
        if(sscanf(optarg, "%lg", &oAlpha) != 1)
	{
	  usage = 1;
	}
	break;
      case 'd':
        if(sscanf(optarg, "%lg", &iDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 'e':
        if(sscanf(optarg, "%lg", &oDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 'D':
        if(sscanf(optarg, "%lg", &delta) != 1)
	{
	  usage = 1;
	}
	break;
      case 'M':
        if(sscanf(optarg, "%lg", &samFac) != 1)
	{
	  usage = 1;
	}
	break;
      case 'T':
        if(sscanf(optarg, "%lg", &tau) != 1)
	{
	  usage = 1;
	}
	break;
      case 'X':
        if(sscanf(optarg, "%d", &planeIdx) != 1)
	{
	  usage = 1;
	}
	else
	{
	  planeFlg = 1;
	}
	break;
      case 'h':
      default:
        usage = 1;
	break;
    }
  }
  if(ok && (optind < argc))
  {
    usage = 1;
  }
  if((outObjFileStr == NULL) || (*outObjFileStr == '\0') ||
     (sFileStr == NULL) || (*sFileStr == '\0') ||
     (iFileStr == NULL) || (*iFileStr == '\0') ||
     (oFileStr == NULL) || (*oFileStr == '\0'))
  {
    usage = 1;
  }
  ok = ok && !usage;
  /* Read the surface, inside and outside points. */
  if(ok)
  {
    if(strcmp(sFileStr, "-"))
    {
      if((fP = fopen(sFileStr, "r")) == NULL)
      {
	ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open surface points file %s (%s).\n",
		       *argv, sFileStr, errMsg);
      }
    }
    else
    {
      fP = stdin;
    }
  }
  if(ok)
  {
    sPts = WlzCFPReadVtxArray(fP, dim, &nSPts, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	             "%s: failed to read points from file %s (%s).\n",
		     *argv, sFileStr, errMsg);
    }
  }
  if(fP)
  {
    if(strcmp(sFileStr, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(ok)
  {
    if(strcmp(sFileStr, "-"))
    {
      if((fP = fopen(iFileStr, "r")) == NULL)
      {
	ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open inside points file %s (%s).\n",
		       *argv, iFileStr, errMsg);
      }
    }
    else
    {
      fP = stdin;
    }
  }
  if(ok)
  {
    iPts = WlzCFPReadVtxArray(fP, dim, &nIPts, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	             "%s: failed to read points from file %s (%s).\n",
		     *argv, sFileStr, errMsg);
    }
  }
  if(fP)
  {
    if(strcmp(iFileStr, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(ok)
  {
    if(strcmp(sFileStr, "-"))
    {
      if((fP = fopen(oFileStr, "r")) == NULL)
      {
	ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open outside points file %s (%s).\n",
		       *argv, oFileStr, errMsg);
      }
    }
    else
    {
      fP = stdin;
    }
  }
  if(ok)
  {
    oPts = WlzCFPReadVtxArray(fP, dim, &nOPts, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	             "%s: failed to read points from file %s (%s).\n",
		     *argv, oFileStr, errMsg);
    }
  }
  if(fP)
  {
    if(strcmp(oFileStr, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(ok)
  {
    bgdV.v.dbv = oDist;
    bgdV.type = WLZ_GREY_DOUBLE;
    switch(dim)
    {
      case 2:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      case 3:
	/* Compute ounding box of points. */
	bBox0.d3 = WlzBoundingBoxVtx3D(nSPts, sPts.d3, NULL);
	bBox1.d3 = WlzBoundingBoxVtx3D(nIPts, iPts.d3, NULL);
	bBox0.d3 = WlzBoundingBoxUnion3D(bBox0.d3, bBox1.d3);
	bBox1.d3 = WlzBoundingBoxVtx3D(nOPts, oPts.d3, NULL);
	bBox0.d3 = WlzBoundingBoxUnion3D(bBox0.d3, bBox1.d3);
	/* Create distance object that will be used to accumulate distance
	   values. If only a 2D contour is required then create a distance
	   object with a single plane. */
	dObj = WlzMakeCuboid((planeFlg)?
	                     planeIdx: (int )floor(bBox0.d3.zMin) - 2,
			     (planeFlg)?
			     planeIdx: (int )ceil(bBox0.d3.zMax) + 2,
			     (int )floor(bBox0.d3.yMin) - 2,
			     (int )ceil(bBox0.d3.yMax) + 2,
			     (int )floor(bBox0.d3.xMin) - 2,
			     (int )ceil(bBox0.d3.xMax) + 2,
			     WLZ_GREY_DOUBLE, bgdV,
			     NULL, NULL, &errNum);
	if((dObj == NULL) || (errNum != WLZ_ERR_NONE))
	{
	  ok = 0;
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	  		 "%s: failed to allocate distance object (%s).\n",
			 *argv, errMsg);

	}
	/* Compute contour using given distance object. */
	outDom.ctr = WlzContourFromPoints(dObj, WLZ_VERTEX_D3,
					  nSPts, sPts, sAlpha,
					  nIPts, iPts, iDist, iAlpha,
					  nOPts, oPts, oDist, oAlpha,
					  delta, tau, samFac,
					  &errNum);
	WlzFreeObj(dObj);
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	             "%s: failed to compute contour from points (%s).\n",
		     *argv, errMsg);
    }
  }
  /* Compute the contour from the point sets. */
  if(ok)
  {
    outObj = WlzMakeMain(WLZ_CONTOUR, outDom, dumVal, NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to create output woolz object (%s).\n",
                     argv[0], errMsg);
    }
    if(outObj)
    {
      outDom.core = NULL;
    }
  }
  /* Write out the contour object. */
  if(ok)
  {
    if(strcmp(outObjFileStr, "-"))
    {
      if((fP = fopen(outObjFileStr, "w")) == NULL)
      {
	ok = 0;
	errNum = WLZ_ERR_WRITE_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to open output object file %s (%s).\n",
		       *argv, outObjFileStr, errMsg);
      }
    }
    else
    {
      fP = stdout;
    }
    if(ok)
    {
      if((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write output object to file %s (%s).\n",
		       *argv, outObjFileStr, errMsg);
      }
    }
    if(strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  AlcFree(sPts.v);
  AlcFree(iPts.v);
  AlcFree(oPts.v);
  (void )WlzFreeDomain(outDom);
  (void )WlzFreeObj(outObj);
  /* Tidy up and exit. */
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    "WlzContourFromPoints [-h] [-2] [-3] [-o<out object>]\n"
    "                     [-S<surface points file>]\n"
    "                     [-I<inside points file>] [-O<outside points file>]\n"
    "                     [-a#] [-b#] [-c#] [-d#] [-e#] [-D#] [-M#] [-T#]\n"
    "                     [-X#]\n"
    "Options:\n"
    "  -h  Print this usage message.\n"
    "  -2  Points are 2D and contour is a curve in 2D.\n"
    "  -3  Points are 3D and contour is a surface in 3D.\n"
    "  -o  Output object file name.\n"
    "  -S  On surface points file.\n"
    "  -I  Inside points file.\n"
    "  -O  Outside points file.\n"
    "  -a  Surface points alpha value.\n"
    "  -b  Inside points alpha value.\n"
    "  -c  Outside points alpha value.\n"
    "  -d  Inside points distance value.\n"
    "  -e  Outside points distance value.\n"
    "  -D  Delta multiorder spline soothing value.\n"
    "  -M  Distance object sampling factor (can be used to reduce the\n"
    "      number of surface facets.)\n"
    "  -T  Tau multiorder spline soothing value.\n"
    "  -X  Plane for 2D contour from 3D point sets.\n"
    "Computes a contour from three point sets: inside; outside and on a\n"
    "surface.\n"
    "All points are given using the following ascii format:\n"
    "  <vertex x> <vertex y> [<vertex z>]\n");
  }
  return(!ok);
}

/*!
* \return	New vertex array.
* \brief	Reads an array of vertices (2D or 3D) from the given file.
* \param	fP			File pointer.
* \param	dim			Dimension which MUST be either 2 or 3.
* \param	dstNVtx			Destination pointer for the number of
*					vertices read, MUST not be NULL.
* \param	dstErr			Used to return Woolz error code,
*					MUST not be NULL.
*/
static WlzVertexP 	WlzCFPReadVtxArray(FILE *fP, int dim, int *dstNVtx, 
					   WlzErrorNum *dstErr)
{
  int		datCnt = 0,
  		datMax = 0;
  char		lnBuf[WLZ_CFP_READLN_LEN];
  WlzVertex	datV;
  WlzVertexP	datVP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  datVP.v = NULL;
  while((errNum == WLZ_ERR_NONE) &&
	(fgets(lnBuf, WLZ_CFP_READLN_LEN, fP) != NULL))
  {
    lnBuf[WLZ_CFP_READLN_LEN - 1] = '\0';
    switch(dim)
    {
      case 2:
	errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      case 3:
        if(sscanf(lnBuf, "%lg %lg %lg",
	          &(datV.d3.vtX), &(datV.d3.vtY), &(datV.d3.vtZ)) != 3)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  if(datCnt >= datMax)
	  {
	    datMax = (datMax + 1024) * 2;
	    if((datVP.v = AlcRealloc(datVP.v,
	    			     sizeof(WlzDVertex3) * datMax)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *(datVP.d3 + datCnt++) = datV.d3;
	}
	break;
    }
  }
  *dstErr = errNum;
  *dstNVtx = datCnt;
  return(datVP);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
