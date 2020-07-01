#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstFitBSpline_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzTstFitBSpline.c
* \author       Bill Hill
* \date         June 2020
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2020],
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
* \brief	Tests for the Woolz B-spline fitting code.
* \ingroup	BinWlzTst
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>

static void			WlzTstPrintPoints(
				  FILE *fP,
				  WlzObject *obj);

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		option,
		order = 3,
		usage = 0,
		closed = 0,
		nEval = -1,
		verbose = 0,
		showFacts = 0,
		printPoints = 0,
		useTestData = 0;
  double	smooth = 0.01;
  WlzBSpline	*bs = NULL;
  WlzObject	*iObj = NULL,
  		*oObj = NULL,
		*sObj = NULL;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		*iFile = NULL,
  		*oFile = NULL,
		*pFile = NULL;
  const char 	*ots = NULL;
  static char	optList[] = "cfFhpPtTve:k:o:O:s:",
		defFile[] = "-";
  const int	testDataSz = 20;
  const double  testData[3][20] =
                {
                  {
                    1.5 , 1.36747888, 1.017313 , 0.57262637, 0.18300587,
                    -0.03447032, -0.0394885, 0.12006954, 0.3337372, 0.47972797,
                    0.47972797, 0.3337372, 0.12006954, -0.0394885, -0.03447032,
                    0.18300587, 0.57262637, 1.017313, 1.36747888, 1.5
                  },
                  {
                    0.0, .469456091, .791806489, .876469902,
                    .722673829, .415994951, .090024708, -.130430291,
                    -.180609543, -.08005244, .08005244, .180609543,
                    .130430291, -.090024708, -.41599495, -.722673829,
                    -.876469902, -.791806489, -.469456091, 0.0
                  },
                  {
                    0.0 , 0.33069396, 0.66138793, 0.99208189, 1.32277585,
                    1.65346982, 1.98416378, 2.31485774, 2.64555171, 2.97624567,
                    3.30693964, 3.6376336 , 3.96832756, 4.29902153, 4.62971549,
                    4.96040945, 5.29110342, 5.62179738, 5.95249134, 6.28318531
                  }
                };

  iFile = oFile = pFile = defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'c':
        closed = 1;
	break;
      case 'f':
        showFacts = 1;
	break;
      case 'F':
        showFacts = 2;
	break;
      case 'p':
        printPoints |= 1;
	break;
      case 'P':
        printPoints |= 2;
	break;
      case 't':
        useTestData = 2;
	break;
      case 'T':
        useTestData = 3;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'e':
        usage = (sscanf(optarg, "%d", &nEval) != 1);
	break;
      case 'k':
         usage = (sscanf(optarg, "%d", &order) != 1) ||
	         (order < WLZ_BSPLINE_ORDER_MIN) ||
		 (order > WLZ_BSPLINE_ORDER_MAX);
        break;
      case 'o':
        oFile = optarg;
	break;
      case 'O':
        pFile = optarg;
	break;
      case 's':
        usage = (sscanf(optarg, "%lg", &smooth) != 1);
	break;
      case 'h':   /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(!usage && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      iFile = argv[optind];
    }
  }
  if(usage)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(useTestData)
    {
      case 2: /* FALLTHROUGH */
      case 3:
	{
	  WlzDomain tDom = {0};
	  WlzVertexP nullP = {0};

	  tDom.pts = WlzMakePoints(
	      (useTestData == 2)?WLZ_POINTS_2D: WLZ_POINTS_3D, 0, nullP,
	      testDataSz, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzValues nullVal = {0};

	    iObj = WlzMakeMain(WLZ_POINTS, tDom, nullVal, NULL, NULL, &errNum);
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    if(iObj)
	    {
	      (void )WlzFreeObj(iObj);
	    }
	    else
	    {
	      (void )WlzFreeDomain(tDom);
	    }
	  }
	  else
          {
	    if(useTestData == 2)
	    {
	      int	i;
	      WlzDVertex2 *p;

	      p = tDom.pts->points.d2;
	      for(i = 0; i < testDataSz; ++i)
	      {
	        p[i].vtX = testData[2][i];
		p[i].vtY = testData[1][i];
	      }
	    }
	    else
	    {
	      int	i;
	      WlzDVertex3 *p;

	      p = tDom.pts->points.d3;
	      for(i = 0; i < testDataSz; ++i)
	      {
		p[i].vtX = testData[0][i];
		p[i].vtY = testData[1][i];
	        p[i].vtZ = testData[2][i];
	      }
	    }
	    tDom.pts->nPoints = testDataSz;
	  }
	  ots = WlzStringFromErrorNum(errNum, NULL);
	  if(verbose)
	  {
	    (void )fprintf(stderr,
		"%s: Input object created (%dD) %s\n",
		*argv, useTestData, ots);
	  }
	}
        break;
      default:
	errNum = WLZ_ERR_READ_EOF;
	if((fP = (strcmp(iFile, "-")?  fopen(iFile, "r"): stdin)) != NULL)
	{
	  iObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL);
	}
	if(fP && strcmp(iFile, "-"))
	{
	  (void )fclose(fP);
	  fP = NULL;
	}
        ots = WlzStringFromErrorNum(errNum, NULL);
	if(verbose)
	{
	  (void )fprintf(stderr,
	      "%s: Input object read from >%s< %s\n", *argv, iFile, ots);
	}
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(verbose)
    {
      ots = WlzStringFromObjType(iObj, &errNum);
      (void )fprintf(stderr,
          "%s: Input object type = %s\n", *argv, (ots)? ots: "NULL");
      if(showFacts)
      {
        (void )WlzObjectFacts(iObj, stderr, NULL, showFacts > 1);
      }
    }
    if(((printPoints & 1) != 0) && (iObj->type == WLZ_POINTS))
    {
      if((fP = (strcmp(pFile, "-")?  fopen(pFile, "w"): stdout)) != NULL)
      {
	(void )fprintf(fP, "{\n\"indata\": {\n\"points\": [\n");
        WlzTstPrintPoints(fP, iObj);
	(void )fprintf(fP, "]}");
	if(printPoints & 2)
	{
	  (void )fprintf(fP, ",\n");
	}
	else
	{
	  (void )fprintf(fP, "}\n");
	}
      }
      if(fP && strcmp(pFile, "-"))
      {
        (void )fclose(fP);
	fP = NULL;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(verbose)
    {
      (void )fprintf(stderr,
          "%s: Calling WlzBSplineFromObj(iObj, %d, %d, %g, &errNum)\n",
	  *argv, order, closed, smooth);

    }
    if(iObj->type != WLZ_SPLINE)
    {
      bs = WlzBSplineFromObj(iObj, order, closed, smooth, &errNum);
    }
    if(verbose)
    {
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr, "%s: bs = %p, errNum = %s\n", *argv, bs, ots);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain dom;
    WlzValues nullVal = {0};

    dom.bs = bs;
    sObj = WlzMakeMain(WLZ_SPLINE, dom, nullVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(verbose)
    {
      ots = WlzStringFromObjType(sObj, &errNum);
      (void )fprintf(stderr,
          "%s: Output object type = %s\n", *argv, (ots)? ots: "NULL");
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr, "%s: errNum = %s\n", *argv, ots);
      if((errNum == WLZ_ERR_NONE) && showFacts)
      {
        (void )WlzObjectFacts(sObj, stderr, NULL, showFacts > 1);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nEval < 0)
    {
      oObj = sObj;
      sObj = NULL;
    }
    else
    {
      int	n;
      WlzDomain dom = {0};
      WlzVertexP nullP = {0};

      n = (nEval == 0)? bs->nKnots: nEval;
      dom.pts = WlzMakePoints(
          (bs->type == WLZ_BSPLINE_C2D)?WLZ_POINTS_2D: WLZ_POINTS_3D,
	  0, nullP, nEval, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	dom.pts->nPoints = n;
        errNum = WlzBSplineEval(bs, n, dom.pts->points);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        WlzValues nullVal = {0};

        oObj = WlzMakeMain(WLZ_POINTS, dom, nullVal, NULL, NULL, &errNum);
      }
      if(errNum != WLZ_ERR_NONE)
      {
        (void )WlzFreeDomain(dom);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((printPoints & 2) != 0) && (oObj->type == WLZ_POINTS))
    {
      if((fP = (strcmp(pFile, "-")?  fopen(pFile, "a"): stdout)) != NULL)
      {
	if((printPoints & 1) == 0)
	{
	  (void )fprintf(fP, "{\n");
	}
	(void )fprintf(fP,
	    "\"outdata\": {\n\"knots\": %d,\n\"smooth\": %g,\n"
	    "\"order\": %d,\n\"dim\": %d,\n\"points\": [\n",
	    bs->nKnots, smooth, order, (bs->type == WLZ_BSPLINE_C2D)? 2: 3);
        WlzTstPrintPoints(fP, oObj);
	(void )fprintf(fP, "]}\n}\n");
      }
      if(fP && strcmp(pFile, "-"))
      {
        (void )fclose(fP);
	fP = NULL;
      }
    }
    if(verbose)
    {
      (void )fprintf(stderr, "%s: Opening output file >%s<\n", *argv, oFile);
    }
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(oFile, "-")?  fopen(oFile, "w"): stdout)) != NULL)
    {
      if(verbose)
      {
	(void )fprintf(stderr,
	    "%s: Calling WlzWriteObj(%p, %p)\n", *argv, fP, oObj);
      }
      errNum = WlzWriteObj(fP, oObj);
    }
    if(fP && strcmp(iFile, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
    if(verbose)
    {
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr, "%s: errNum = %s\n", *argv, ots);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(verbose)
    {
      (void )fprintf(stderr,
          "%s: Output object written to file  >%s<\n", *argv, oFile);
    }
  }
  (void )WlzFreeObj(iObj);
  (void )WlzFreeObj(oObj);
  (void )WlzFreeObj(sObj);
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s%s%s%s%s",
    *argv,
    " [-h] [-c] [-f] [-F] [-p] [-P] [-t] [-T] [-v]\n"
    "\t\t[-e<eval>] [-k<order>] [-o<out object>] [-O<out print>]\n"
    "\t\t[-s<smoothing>] [<in object>]\n",
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -c  Fit a closed (periodic) B-spline curve.\n"
    "  -f  Show object facts if verbose output.\n"
    "  -F  Show object 'many' facts if verbose output.\n"
    "  -p  Print input points that the spline is fitted to (only for\n"
    "      input object of type WLZ_POINTS).\n"
    "  -P  Print output points from evaluating the spline (only for\n"
    "      input (and consequently the) output object of type WLZ_POINTS).\n"
    "  -t  Ignore input file (if given) and use in built 2D test data.\n"
    "  -T  Ignore input file (if given) and use in built 3D test data.\n"
    "  -v  Show verbose output.\n"
    "  -e  Evaluate the spline and output it's evaluation instead of\n"
    "      outputting the spline.\n"
    "  -k  Spline order.\n"
    "  -o  Output object file name.\n"
    "  -O  Output print data file name.\n"
    "  -s  Smoothing parameter (0.0 implies no smoothing).\n"
    "Test binary to read / fit / write B-spline objects.\n");
  }
  exit(errNum);
}

static void			WlzTstPrintPoints(
				  FILE *fP,
				  WlzObject *obj)
{
  if(obj && (obj->type == WLZ_POINTS) && (obj->domain.core != NULL))
  {
    int		i,
		n,
    		n1;
    WlzPoints	*pts;

    pts = obj->domain.pts;
    n = pts->nPoints;
    n1 = n - 1;
    switch(pts->type)
    {
      case WLZ_POINTS_2I:
	for(i = 0; i < n; ++i)
	{
	  WlzIVertex2 *p;

	  p = pts->points.i2 + i;
	  (void )fprintf(fP, "[%d,%d]", p->vtX, p->vtY);
	  if(i < n1)
	  {
	    (void )fprintf(fP, ",");
	  }
	  (void )fprintf(fP, "\n");
	}
        break;
      case WLZ_POINTS_2D:
	for(i = 0; i < n; ++i)
	{
	  WlzDVertex2 *p;

	  p = pts->points.d2 + i;
	  (void )fprintf(fP, "[%g,%g]", p->vtX, p->vtY);
	  if(i < n1)
	  {
	    (void )fprintf(fP, ",");
	  }
	  (void )fprintf(fP, "\n");
	}
        break;
      case WLZ_POINTS_3I:
	for(i = 0; i < n; ++i)
	{
	  WlzIVertex3 *p;

	  p = pts->points.i3 + i;
	  (void )fprintf(fP, "[%d,%d,%d]", p->vtX, p->vtY, p->vtZ);
	  if(i < n1)
	  {
	    (void )fprintf(fP, ",");
	  }
	  (void )fprintf(fP, "\n");
	}
        break;
      case WLZ_POINTS_3D:
	for(i = 0; i < n; ++i)
	{
	  WlzDVertex3 *p;

	  p = pts->points.d3 + i;
	  (void )fprintf(fP, "[%g,%g,%g]", p->vtX, p->vtY, p->vtZ);
	  if(i < n1)
	  {
	    (void )fprintf(fP, ",");
	  }
	  (void )fprintf(fP, "\n");
	}
        break;
      default:
        break;
    }
  }
}
