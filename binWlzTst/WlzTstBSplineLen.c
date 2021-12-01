#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstBSplineLen_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzTstBSplineLen.c
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
* \brief	Tests for the Woolz B-spline length code.
* \ingroup	BinWlzTst
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		option,
		n = 100,
		dim = 2,
		order = 3,
		usage = 0,
		verbose = 0;
  double	tB = 0.0,
  		tE = 1.0,
		bLen = 0.0,
		cLen = 0.0,
  		rad = 10.0,
		smooth = 0.0;
  WlzVertexP	eval = {0},
  		circle = {0};
  WlzVertexType	vType = WLZ_VERTEX_ERROR;
  WlzDomain	bsDom = {0};
  double	*tbuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char 	*ots = NULL;
  static char	optList[] = "hvb:d:e:n:s:";

  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'v':
        verbose = 1;
	break;
      case 'b':
        usage = (sscanf(optarg, "%lg", &tB) != 1) || (tB < 0.0) || (tB > 1.0);
	break;
      case 'd':
        usage = (sscanf(optarg, "%d", &dim) != 1) || (dim < 2) || (dim > 3);
	break;
      case 'e':
        usage = (sscanf(optarg, "%lg", &tE) != 1) || (tE < 0.0) || (tE > 1.0);
	break;
      case 'n':
        usage = (sscanf(optarg, "%d", &n) != 1) || (n < 1);
	break;
      case 's':
        usage = (sscanf(optarg, "%lg", &smooth) != 1) || (smooth < 0.0);
	break;
      case 'h':   /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    double	d;

    if(dim == 2)
    {
      vType = WLZ_VERTEX_D2;
      if(((tbuf = (double *)
                  AlcMalloc(n * sizeof(double))) == NULL) ||
         ((circle.d2 = (WlzDVertex2 *)
                       AlcMalloc(2 * n * sizeof(WlzDVertex2))) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	double	d1;
	WlzDVertex2 *v;

	v = circle.d2;
	eval.d2 = v + n;
	d1 = 1.0 / (n - 1);
	d = 2.0 * ALG_M_PI * d1;
	for(i = 0; i < n; ++i)
	{
	  tbuf[i] = i * d1;
	  v[i].vtX = rad * cos(d * i);
	  v[i].vtY = rad * sin(d * i);
	}
	if(verbose)
        {
	  (void )fprintf(stderr, "%s: i tbuf circle.d2.vtX circle.d2.vtY\n",
	      *argv);
	}
	for(i = 0; i < n; ++i)
	{
	  if(i > 1)
	  {
	    WlzDVertex2 t;

	    WLZ_VTX_2_SUB(t, v[i - 1], v[i]);
	    cLen += WLZ_VTX_2_LENGTH(t);
	  }
	  if(verbose)
          {
	    (void )fprintf(stderr, "%d %g %g %g\n",
	        i, tbuf[i], v[i].vtX, v[i].vtY);
	  }
	}
      }
    }
    else
    {
      vType = WLZ_VERTEX_D3;
      if(((tbuf = (double *)
                  AlcMalloc(n * sizeof(double))) == NULL) ||
         ((circle.d3 = (WlzDVertex3 *)
                       AlcMalloc(2 * n * sizeof(WlzDVertex3))) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	double	d1;
	WlzDVertex3 *v;

	v = circle.d3;
	eval.d3 = v + n;
	d1 = 1.0 / (n - 1);
	d = 2.0 * ALG_M_PI * d1;
	for(i = 0; i < n; ++i)
	{
	  tbuf[i] = i * d1;
	  v[i].vtX = rad * cos(d * i);
	  v[i].vtY = rad * sin(d * i);
	}
	if(verbose)
        {
	  (void )fprintf(stderr,
	      "%s: i tbuf circle.d3.vtX circle.d3.vtY circle.d3.vtZ\n",
	      *argv);
	}
	for(i = 0; i < n; ++i)
	{
	  if(i > 1)
	  {
	    WlzDVertex3 t;

	    WLZ_VTX_3_SUB(t, v[i - 1], v[i]);
	    cLen += WLZ_VTX_3_LENGTH(t);
	  }
	  if(verbose)
          {
	    (void )fprintf(stderr, "%d %g %g %g %g\n",
	        i, tbuf[i], v[i].vtX, v[i].vtY, v[i].vtZ);
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(verbose)
    {
      (void )fprintf(stderr,
          "%s: Calling "
	  "WlzBSplineFromVertices(%s, %d, circle, %d, 0, %g, &errNum)\n",
	  *argv, WlzStringFromVertexType(vType, NULL), n, order, smooth);
    }
    bsDom.bs = WlzBSplineFromVertices(vType, n, circle, order, 0, smooth,
	&errNum);
    if(verbose)
    {
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr, "%s: errNum = %s\n", *argv, ots);
    }
  }
  if((errNum == WLZ_ERR_NONE) && verbose)
  {
    int		i;

    (void )fprintf(stderr,
        "%s: Calling WlzBSplineEval(bsDom.bs, %d, tbuf, 0, eval)\n",
	*argv, n);
    errNum = WlzBSplineEval(bsDom.bs, n, tbuf, 0, eval);
    ots = WlzStringFromErrorNum(errNum, NULL);
    (void )fprintf(stderr, "%s: errNum = %s\n", *argv, ots);
    if(dim == 2)
    {
      WlzDVertex2 *v;

      v = eval.d2;
      (void )fprintf(stderr, "%s: i, tbuf, eval.d2.vtX, eval.d2.vtY\n", *argv);
      for(i = 0; i < n; ++i)
      {
	(void )fprintf(stderr, "%d %g %g %g\n",
	    i, tbuf[i], v[i].vtX, v[i].vtY);
      }
    }
    else
    {
      WlzDVertex3 *v;

      v = eval.d3;
      (void )fprintf(stderr,
          "%s: i, tbuf, eval.d3.vtX, eval.d3.vtY, eval.d3.vtZ\n", *argv);
      for(i = 0; i < n; ++i)
      {
	(void )fprintf(stderr, "%d %g %g %g %g\n",
	    i, tbuf[i], v[i].vtX, v[i].vtY, v[i].vtZ);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    double	circ;

    circ = 2 * rad * ALG_M_PI;
    if(verbose)
    {
      (void )fprintf(stderr,
          "%s: Calling \n"
	  "\t\tWlzBSplineLength(bsDom.bs, 0.0, 1.0, 1.0, 1.0, 1.0, &errNum)\n",
	  *argv);

    }
    bLen = WlzBSplineLength(bsDom.bs, tB, tE, 1.0, 1.0, 1.0, &errNum);
    if(verbose)
    {
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
          "%s: circ = %lg, bLen = %lg, cLen = %lg, errNum = %s\n",
          *argv, circ, bLen, cLen, ots);
    }
    (void )printf("%lg, %lg, %lg, %lg\n",
        circ, cLen, circ * (tE - tB), bLen);
  }
  AlcFree(tbuf);
  AlcFree(circle.v);
  (void )WlzFreeDomain(bsDom);
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s%s%s%s",
    *argv,
    " [-h] [-v] [-b#] [-d (2|3)] [-e#] [-n#] [-s#]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -n  Number of vertices.\n"
    "  -v  Show verbose output.\n"
    "  -b  Begining of length computation in parametric coordinates.\n"
    "  -d  Dimension of test object.\n"
    "  -e  End of length computation in parametric coordinates.\n"
    "  -s  Smoothing value.\n");
  }
  exit(errNum);
}

