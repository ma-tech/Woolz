#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstFourier_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         AlgTstFourier.c
* \author       Bill Hill
* \date         December 2012
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
* \brief	Tests for the libAlg fourier transform code.
* \ingroup	binAlgTst
*/

#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#include <Alc.h>
#include <Alg.h>

static void    			AlgTstFourOutput(
				  int d,
				  int n,
				  int realFlg,
				  double z,
				  double s,
                                  double *a10,
				  double *a11,
				  double **a20,
				  double **a21,
				  double ***a30,
				  double ***a31);
static void			AlgTstFourOutput1(
				  char *str,
				  double z,
				  double s,
				  double *a,
				  int n);
static void			AlgTstFourOutput2(
				  char *str,
				  double z,
				  double s,
				  double **a,
				  int n);
static void			AlgTstFourOutput3(
				  char *str,
				  double z,
				  double s,
				  double ***a,
				  int n);

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int           d = 1,
  		n = 8,
		asyFlg = 0,
		outFlg = 1,
		rndFlg = 0,
		realFlg = 0,
		scaleFlg = 0,
		timeFlg = 0,
		useBuf = 0,
		option,
  		ok = 1,
		usage = 0;
  double	scale = 1.0,
  		z = 0.0;
  double	*a1[2];
  double	**a2[2];
  double	***a3[2];
  struct timeval times[3];
  const	double	asy = ((1.0 / 7.0) + ALG_M_PI) * ALG_M_E / 13;
  static char	optList[] = "chrsuyNRTd:n:z:";

  a1[0] = a1[1] = NULL;
  a2[0] = a2[1] = NULL;
  a3[0] = a3[1] = NULL;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'c':
	realFlg = 0;
	break;
      case 'r':
	realFlg = 1;
	break;
      case 's':
	scaleFlg = 1;
	break;
      case 'u':
	useBuf = 1;
	break;
      case 'y':
	asyFlg = 1;
	break;
      case 'N':
	outFlg = 0;
	break;
      case 'R':
	asyFlg = 0;
	rndFlg = 1;
	break;
      case 'T':
	timeFlg = 1;
	break;
      case 'd':
	if((sscanf(optarg, "%d", &d) != 1) || (d < 1) || (d > 3))
	{
	  usage = 1;
	}
	break;
      case 'n':
	if((sscanf(optarg, "%d", &n) != 1) || (n < 4) ||
	   (AlgBitIsPowerOfTwo(n) == 0))
	{
	  usage = 1;
	}
	break;
      case 'z':
	if(sscanf(optarg, "%lg", &z) != 1)
	{
	  usage = 1;
	}
	else
	{
	  z = fabs(z);
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  /* Allocate array(s) and set values. */
  if(ok)
  {
    int	n2,
    	s1,
	s2;

    if(rndFlg)
    {
      s1 = 0;
      s2 = n;
    }
    else
    {
      s1 = n / 4;
      s2 = n - s1;
    }
    n2 = (realFlg)? n: 2 * n;
    if(rndFlg)
    {
      AlgRandSeed(0L);
    }
    switch(d)
    {
      case 1:
        if((a1[0] = AlcCalloc(n2, sizeof(double))) == NULL)
	{
	  ok = 0;
	}
	else
	{
	  int	i;

	  for(i = s1; i < s2; ++i)
	  {
	    a1[0][i] = (rndFlg)? AlgRandUniform(): 1.0;
	  }
	  if(asyFlg)
	  {
	    a1[0][s1] = asy;
	  }
	  if(realFlg == 0)
	  {
	    a1[1] = a1[0] + n;
	  }
	}
	break;
      case 2:
        if(AlcDouble2Calloc(&(a2[0]), n2, n) != ALC_ER_NONE)
	{
	  ok = 0;
	}
	else
	{
	  int	i,
	  	j;

	  for(j = s1; j < s2; ++j)
	  {
	    for(i = s1; i < s2; ++i)
	    {
	      a2[0][j][i] = (rndFlg)? AlgRandUniform(): 1.0;
	    }
	  }
	  if(asyFlg)
	  {
	    a2[0][s1][s1] = asy;
	  }
	  if(realFlg == 0)
	  {
	    a2[1] = a2[0] + n;
	  }
	}
	break;
      case 3:
        if(AlcDouble3Calloc(&(a3[0]), n2, n, n) != ALC_ER_NONE)
	{
	  ok = 0;
	}
	else
	{
	  int	i,
	  	j,
		k;

	  for(k = s1; k < s2; ++k)
	  {
	    for(j = s1; j < s2; ++j)
	    {
	      for(i = s1; i < s2; ++i)
	      {
		a3[0][k][j][i] = (rndFlg)? AlgRandUniform(): 1.0;
	      }
	    }
	  }
	  if(asyFlg)
	  {
	    a3[0][s1][s1][s1] = asy;
	  }
	  if(realFlg == 0)
	  {
	    a3[1] = a3[0] + n;
	  }
	}
	break;
      default:
        break;
    }
    if(ok == 0)
    {
      (void )fprintf(stderr,
		     "%s: Failed to allocate array(s).\n",
		     *argv);
    }
  }
  if(ok)
  {
    /* Output array values. */
    if(outFlg)
    {
      AlgTstFourOutput(d, n, realFlg, z, 1.0,
		       a1[0], a1[1], a2[0], a2[1], a3[0], a3[1]);
    }
    /* Do the forward transform. */
    if(timeFlg)
    {
      gettimeofday(times + 0, NULL);
    }
    switch(d)
    {
      case 1:
	scale /= (scaleFlg)? sqrt(n): 1.0;
        if(realFlg == 0)
	{
	  AlgFour1D(a1[0], a1[1], n, 1);
	}
	else
	{
	  AlgFourReal1D(a1[0], n, 1);
	}
	break;
      case 2:
	scale /= (scaleFlg)? n: 1.0;
        if(realFlg == 0)
	{
	  ok = (AlgFour2D(a2[0], a2[1], useBuf, n, n) == ALG_ERR_NONE);
	}
	else
	{
	  ok = (AlgFourReal2D(a2[0], useBuf, n, n) == ALG_ERR_NONE);
	}
        break;
      case 3:
	scale /= (scaleFlg)? n * sqrt(n): 1.0;
        if(realFlg == 0)
	{
	  ok = (AlgFour3D(a3[0], a3[1], useBuf, n, n, n) == ALG_ERR_NONE);
	}
	else
	{
	  ok = (AlgFourReal3D(a3[0], useBuf, n, n, n) == ALG_ERR_NONE);
	}
        break;
        break;
      default:
        break;
    }
    if(timeFlg)
    {
      gettimeofday(times + 1, NULL);
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: elapsed time for forward transform %gus\n",
		     *argv,
		     (1000000.0 * times[2].tv_sec) + times[2].tv_usec);
    }
  }
  if(ok)
  {
    /* Output array values. */
    if(outFlg)
    {
      AlgTstFourOutput(d, n, realFlg, z, scale,
		       a1[0], a1[1], a2[0], a2[1], a3[0], a3[1]);
    }
    /* Do the inverse transform. */
    if(timeFlg)
    {
      gettimeofday(times + 0, NULL);
    }
    switch(d)
    {
      case 1:
	scale /= (scaleFlg)? sqrt(n): 1.0;
        if(realFlg == 0)
	{
	  AlgFourInv1D(a1[0], a1[1], n, 1);
	}
	else
	{
	  AlgFourRealInv1D(a1[0], n, 1);
	}
	break;
      case 2:
	scale /= (scaleFlg)?  n: 1.0;
        if(realFlg == 0)
	{
	  ok = (AlgFourInv2D(a2[0], a2[1], useBuf, n, n) == ALG_ERR_NONE);
	}
	else
	{
	  ok = (AlgFourRealInv2D(a2[0], useBuf, n, n) == ALG_ERR_NONE);
	}
        break;
      case 3:
	scale /= (scaleFlg)? n * sqrt(n): 1.0;
        if(realFlg == 0)
	{
	  ok = (AlgFourInv3D(a3[0], a3[1], useBuf, n, n, n) == ALG_ERR_NONE);
	}
	else
	{
	  ok = (AlgFourRealInv3D(a3[0], useBuf, n, n, n) == ALG_ERR_NONE);
	}
        break;
      default:
        break;
    }
    if(timeFlg)
    {
      gettimeofday(times + 1, NULL);
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: elapsed time for inverse transform %gus\n",
		     *argv,
		     (1000000.0 * times[2].tv_sec) + times[2].tv_usec);
    }
  }
  if(ok)
  {
    /* Output array values. */
    if(outFlg)
    {
      AlgTstFourOutput(d, n, realFlg, z, scale,
		       a1[0], a1[1], a2[0], a2[1], a3[0], a3[1]);
    }
  }
  AlcFree(a1[0]);
  AlcDouble2Free(a2[0]);
  AlcDouble3Free(a3[0]);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-c] [-r] [-u] [-y] [-N] [-R] [-T]\n"
    "\t\t[-d #] [-n #] [-z #]\n"
    "Tests for the libAlg Fourier transform code.\n"
    "Options are:\n"
    "  -c  Complex transforms, as opposed to real (value %s).\n"
    "  -r  Real transforms, as opposed to complex (value %s).\n"
    "  -s  Rescale output (value %s).\n"
    "  -u  Use buffers (value %s).\n"
    "  -N  No array values output (value %s).\n"
    "  -R  Used uniform random values (value %s).\n"
    "  -T  Print execution times (value %s).\n"
    "  -Y  Use slightly asymetric data values (value %s).\n"
    "  -d  Number of dimensions (value %d).\n"
    "  -n  Size of array, must be power or two (value %d)\n"
    "  -z  Value any value less than the absolute value of this is\n"
    "      considered zero in output (value %lg)\n",
    *argv,
    (realFlg)? "false": "true",
    (realFlg)? "true": "false",
    (scaleFlg)? "true": "false",
    (useBuf)? "true": "false",
    (outFlg)? "false": "true",
    (rndFlg)? "false": "true",
    (timeFlg)? "true": "false",
    (asyFlg)? "true": "false",
    d, n, z);
    ok = 0;
  }
  return(!ok);
}

/*!
* \ingroup	AlgTst
* \brief	Outputs the array values.
* \param	d			Array dimension.
* \param	n			Array size.
* \param	realFlg			Real data if non zero.
* \param	z			Absolute minimum value for non zero
* 					values.
* \param	s			Value scale factor.
* \param	a10			Array for 1D real / complex conjugate.
* \param	a11			Array for 1D imaginary.
* \param	a20			Array for 2D real / complex conjugate.
* \param	a21			Array for 2D imaginary.
* \param	a30			Array for 3D real / complex conjugate.
* \param	a31			Array for 3D imaginary.
*/
static void    AlgTstFourOutput(int d, int n, int realFlg,
                                double z, double s,
                                double *a10, double *a11,
				double **a20, double **a21,
				double ***a30, double ***a31)
{
  switch(d)
  {
    case 1:
      AlgTstFourOutput1("a0", z, s, a10, n);
      if(realFlg == 0)
      {
	AlgTstFourOutput1("a1", z, s, a11, n);
      }
      break;
    case 2:
      AlgTstFourOutput2("a0", z, s, a20, n);
      if(realFlg == 0)
      {
	AlgTstFourOutput2("a1", z, s, a21, n);
      }
      break;
    case 3:
      AlgTstFourOutput3("a0", z, s, a30, n);
      if(realFlg == 0)
      {
	AlgTstFourOutput3("a1", z, s, a31, n);
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup	AlgTst
* \brief	Outputs a 1D array of values.
* \param	str			Array name string.
* \param	realFlg			Real data if non zero.
* 					values.
* \param	s			Value scale factor.
* \param	a			Array of values.
* \param	n			Array size.
*/
static void	AlgTstFourOutput1(char *str, double z, double s, double *a,
				  int n)
{
  int		i;

  (void )printf("%s (%d)\n", str, n);
  for(i = 0; i < n; ++i)
  {
    double v;

    v = a[i] * s;
    if(fabs(v) < z)
    {
      v = 0.0;
    }
    (void )printf("%lg ", v);
  }
  (void )printf("\n");
}

/*!
* \ingroup	AlgTst
* \brief	Outputs a 2D array of values.
* \param	str			Array name string.
* \param	realFlg			Real data if non zero.
* 					values.
* \param	s			Value scale factor.
* \param	a			Array of values.
* \param	n			Array size.
*/
static void	AlgTstFourOutput2(char *str, double z, double s, double **a,
				  int n)
{
  int		i,
  		j;

  (void )printf("%s (%d x %d)\n", str, n, n);
  for(j = 0; j < n; ++j)
  {
    for(i = 0; i < n; ++i)
    {
        double v;

	v = a[j][i] * s;
	if(fabs(v) < z)
	{
	  v = 0.0;
	}
      (void )printf("%lg ", v);
    }
    (void )printf("\n");
  }
}

/*!
* \ingroup	AlgTst
* \brief	Outputs a 3D array of values.
* \param	str			Array name string.
* \param	realFlg			Real data if non zero.
* 					values.
* \param	s			Value scale factor.
* \param	a			Array of values.
* \param	n			Array size.
*/
static void	AlgTstFourOutput3(char *str, double z, double s, double ***a,
				  int n)
{
  int		i,
  		j,
		k;

  (void )printf("%s (%d x %d)\n", str, n, n);
  for(k = 0; k < n; ++k)
  {
    for(j = 0; j < n; ++j)
    {
      for(i = 0; i < n; ++i)
      {
        double v;

	v = a[k][j][i] * s;
	if(fabs(v) < z)
	{
	  v = 0.0;
	}
	(void )printf("%lg ", v);
      }
      (void )printf("\n");
    }
    (void )printf("\n");
  }
}
