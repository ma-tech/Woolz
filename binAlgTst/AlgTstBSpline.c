#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstBSpline_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstBSpline.c
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
* \brief	Simple test for B-Splie fitting code from AlgBSpline.c
* 		fitting and evaluating periodic, non-periodic and three
* 		dimensional splines.
* \ingroup	binAlgTst
*/
#include <float.h>
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

#define DATA_SZ	(20)

int				main(int argc, char *argv[])
{

  int		i,
		nc,
		nd,
		nest,
		lwrk,
		k = 3,
		n = 0,
		dim = 2,
		iopt = 0,
		deriv = -1,
		m = DATA_SZ,
		usePData = 0,
		option,
		periodic = 0,
		usage = 0,
		verbose = 0;
  double	fp = 0.0,
  		s = 0.01;
  int		*iwrk = NULL;
  char		*fn = "unknown";
  char  	optList[] = "23chpqvd:";
  double	*c = NULL,
  		*t = NULL,
		*u = NULL,
  		*w = NULL,
		*tmp = NULL,
		*wrk = NULL;
  double	**out = NULL,
  		**din = NULL;
  AlgError	errNum = ALG_ERR_NONE;
  const double pData[3][20] =
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

  opterr = 0;
  while((option = getopt(argc, argv, optList)) != EOF )
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'c':
        usePData = 1;
	break;
      case 'd':
        usage = (sscanf(optarg, "%d", &deriv) != 1) || (deriv < 0);
	break;
      case 'p':
        periodic = 1;
	break;
      case 'q':
        periodic = 0;
	break;
      case 'v':
        verbose = 1;
	break;
      default:
	usage = 1;
        break;
    }
  }
  if(verbose)
  {
    (void )fprintf(stderr,
	"%s: dim = %d, periodic = %d, usePData = %d verbose = %d\n",
	*argv, dim, periodic, usePData, verbose);
  }
  /* Allocate storage and generate test data. */
  nd = (dim == 3)?3: 1;
  nest = m + 2 * (k + 1);
  nc = nd * nest;
  lwrk = m * (k + 1) + nest * (7 + nd + 5 * k);
  iwrk = (int *)AlcCalloc(nest, sizeof(int));
  fn = "malloc";
  if(((c = (double *)AlcCalloc(nc, sizeof(double))) == NULL) ||
     ((t = (double *)AlcCalloc(nest, sizeof(double))) == NULL) ||
     ((u = (double *)AlcCalloc(m, sizeof(double))) == NULL) ||
     ((w = (double *)AlcCalloc(m, sizeof(double))) == NULL) ||
     ((tmp = (double *)AlcCalloc(m * 3, sizeof(double))) == NULL) ||
     ((wrk = (double *)AlcCalloc(lwrk, sizeof(double))) == NULL) ||
     (AlcDouble2Malloc(&out, 3, m) != ALC_ER_NONE) ||
     (AlcDouble2Malloc(&din, 3, m) != ALC_ER_NONE))
  {
    errNum = ALG_ERR_MALLOC;
  }
  else
  {
    for(i = 0; i < m; ++i)
    {
      double	phi,
		  r;

      u[i] = 0.0;
      w[i] = 1.0;
      if(usePData && (nd == 3) && (m == 20))
      {
	din[0][i] = pData[0][i];
	din[1][i] = pData[1][i];
	din[2][i] = pData[2][i];
      }
      else
      {
	phi = (i * 2.0 * ALG_M_PI ) / (m - 1);
	r = 0.5 + cos(phi);
	din[0][i] = r * cos(phi);
	din[1][i] = r * sin(phi);
	din[2][i] = phi + 10.0;
      }
      out[0][i] = out[1][i] = out[2][i] = 0.0;
      (void )fprintf(stdout, "%1.6f %1.6f %1.6f\n",
	  din[0][i], din[1][i], din[2][i]);
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    if(nd == 1)
    {
      if(periodic)
      {
        fn = "AlgBSplinePerFit";
	errNum = AlgBSplinePerFit(iopt, m, din[2], din[1], w, k, s, nest, &n,
	    t, c, &fp, wrk, iwrk);
	(void )fprintf(stdout, "AlgBSplinePerFit() n = %d, errNum = %d\n",
	    n, errNum);
      }
      else
      {
        fn = "AlgBSplineFit";
	errNum = AlgBSplineFit(iopt, m, din[2], din[1], w, 
	    din[2][0], din[2][m - 1], k, s, nest, &n, t, c,
	    &fp, wrk, iwrk);
	(void )fprintf(stdout, "AlgBSplineFit() n = %d, errNum = %d\n",
	    n, errNum);
      }
    }
    else /* nd == 3 */
    {
      int	ipar = 0,
		mx = nd * m;
      double	ub = 0.0,
		ue = 0.0;

      if(periodic)
      {
	fn = "unknown";
        errNum = ALG_ERR_UNIMPLEMENTED;
      }
      else
      {
	int	ij = 0;

	for(i = 0; i < m; ++i)
	{
	  int	j;

	  for(j = 0; j < 3; ++j)
	  {
	    tmp[ij++] = din[j][i];
	  }
	}
	fn = "AlgBSplineNDFit";
	errNum = AlgBSplineNDFit(iopt, ipar, nd, m, u, mx, tmp, w,
	    ub, ue, k, s, nest, &n, t, &nc, c, &fp, wrk, iwrk);
        (void )fprintf(stdout,
	    "AlgBSplineNDFit() n = %d, nc = %d, errNum = %d\n", n, nc, errNum);
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    (void )fprintf(stdout, "n = %d\n", n);
    (void )fprintf(stdout, "t = \n");
    for(i = 0; i < DATA_SZ; ++i)
    {
      (void )fprintf(stdout, "%1.6f\n", t[i]);
    }
    if(nd != 1)
    {
      (void )fprintf(stdout, "u =\n");
      for(i = 0; i < DATA_SZ; ++i)
      {
	(void )fprintf(stdout, "%1.6f\n", u[i]);
      }
    }
    (void )fprintf(stdout, "c =\n");
    for(i = 0; i < nc; ++i)
    {
      (void )fprintf(stdout, "%1.6f\n", c[i]);
    }
    if(nd == 1)
    {
      if(deriv >= 0)
      {
        fn = "AlgBSplineDeriv";
	errNum = AlgBSplineDer(t, n, c, k, deriv, din[2], out[0], m, wrk);
      }
      else
      {
	fn = "AlgBSplineEval";
	errNum = AlgBSplineEval(t, n, c, k, (nd == 1)? din[2]: u, out[0], m);
      }
    }
    else
    {
      if(deriv >= 0)
      {
        fn = "AlgBSplineDeriv";
	for(i = 0; i < 3; ++i)
	{
	  errNum = AlgBSplineDer(t, n, &(c[i * n]), k, deriv,
	                         u, out[i], m, wrk);
	  if(errNum != ALG_ERR_NONE)
	  {
	    break;
	  }
	}
      }
      else
      {
	fn = "AlgBSplineEval";
	for(i = 0; i < 3; ++i)
	{
	  errNum = AlgBSplineEval(t, n, &(c[i * n]), k, u, out[i], m);
	  if(errNum != ALG_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    for(i = 0; i < DATA_SZ; ++i)
    {
      if(nd == 1)
      {
	(void )fprintf(stdout, "%1.6f %1.6f\n", din[0][i], out[0][i]);
      }
      else /* nd == 3 */
      {
	(void )fprintf(stdout, "%1.6f %1.6f %1.6f\n",
	    out[0][i], out[1][i], out[2][i]);
      }
    }
  }
  AlcFree(c);
  AlcFree(t);
  AlcFree(u);
  AlcFree(w);
  AlcFree(tmp);
  AlcFree(wrk);
  Alc2Free((void **)out);
  Alc2Free((void **)din);
  if(usage)
  {
    (void )fprintf(stderr,
		   "Test AlgBSpline code.\n"
		   "Usage: %s: [-2|3] [-c] [-h] [-p|q]\n"
		   "Options:\n"
		   "  -2  Test in 2D\n"
		   "  -3  Test in 3D\n"
		   "  -c  Use identical data to python code for 3D\n"
		   "  -h  Help prints this usage message.\n"
		   "  -p  Use fit periodic spline.\n"
		   "  -q  Use fit non-periodic spline.\n",
		   *argv);
  }
  if(errNum != ALG_ERR_NONE)
  {
    (void )fprintf(stdout, "%s: Error exit, error (%d) in function %s\n",
        *argv, errNum, fn);
  }
  exit(errNum);
}

/* For the 3D case the ouput can be compared to the following Python code
 * since Python uses Dierckx FITPACK too:
 *

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

phi = np.linspace(0, 2.0 * np.pi, 100)
r = 0.5 + np.cos(phi)
x, y, z = r * np.cos(phi), r * np.sin(phi), phi

# interpolate
tck,u = interpolate.splprep([x,y,z],k=3,s=0.01) # splprep <- FITPACK parcur
out = interpolate.splev(u, tck) # splev

print("# knots = " + str(len(tck[0])))

fig = plt.figure()
ax3d = fig.add_subplot(111, projection='3d')
ax3d.plot(x, y, z, 'k--',
          label='Input data', marker='.',markerfacecolor='k')
ax3d.plot(out[0], out[1], out[2],
          'b',linewidth=2.0,label='Interpolated B-spline')
plt.legend(['Input data', 'Interpolated B-spline', 'True'],
           loc='best')
plt.title('B-Spline interpolation')
plt.show()

 *
 */
