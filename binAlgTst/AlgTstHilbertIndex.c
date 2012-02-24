#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstHilbertIndex_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstHilbertIndex.c
* \author       Bill Hill
* \date         September 2011
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
* \brief	Simple test for AlgHilbertIndex() and AlgHilbertIndexInv().
* \ingroup	binAlgTst
*/
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

static int			AlcTstGetNUInt(
				  unsigned int *v,
				  char *s,
				  int n);
static int			AlgTstSort2(
				  const void *v0,
				  const void *v1);
static int			AlgTstSort3(
				  const void *v0,
				  const void *v1);

int             main(int argc, char **argv)
{
  int           i,
  		n = 2,
		o = 32,
		option,
  		ok = 1,
		usage = 0,
		flgInv = 0,
		flgLin = 0,
		flgMsh = 0;
  unsigned int  *v,
  		*h = NULL,
  		*p = NULL;
  const int	bufMax = 4096;
  static char	buf[4096];
  const char	*optList = "hilmn:o:";

  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'i':
	flgInv = 1;
	break;
      case 'l':
	flgLin = 1;
	break;
      case 'm':
	flgMsh = 1;
	break;
      case 'n':
	if((sscanf(optarg, "%d", &n) != 1) || (n < 1) || (n > 8))
	{
	  usage = 1;
	}
	break;
      case 'o':
	if((sscanf(optarg, "%d", &o) != 1) || (o < 1) || (o > 32))
	{
	  usage = 1;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-i] [-l] [-m] [-n #] [-o #]\n%s\n",
    *argv,
    "A test for AlgHilbertIndex() and AlgHilbertIndexInv().\n"
    "Options are:\n"
    "  i  Compute inverse (ie use AlgHilbertIndexInv()).\n"
    "  l  Simple linear ordering (no Hilbert index).\n"
    "  m  Create a VTK unstructured mesh with a Hilbert curve.\n"
    "  n  Number of dimensions (default 2).\n"
    "  o  Maximum number of bits for coordinate (default 32).\n"
    "Values are read from the stadard input and written to the standard\n"
    "output.\n");
  }
  else
  {
    if(((h = (unsigned int *)AlcMalloc(sizeof(unsigned int) * n)) == NULL) ||
       ((p = (unsigned int *)AlcMalloc(sizeof(unsigned int) * n)) == NULL))
    {
      (void )fprintf(stderr, "%s: Failed to allocate storage.\n", *argv);
      ok = 0;
    }
  }
  if(ok && (usage == 0))
  {
    if(flgMsh)
    {
      int	   np;
      unsigned int m;
      unsigned int *tbl = NULL;

      m = 1 << (o - 1);
      np = (n == 3)? m * m * m: m * m;
      if((flgInv != 0) || (n < 2) || (n > 3) || (o < 2) || (o > 10) ||
	 (m > 1024)|| (np > INT_MAX))
      {
	(void )fprintf(stderr,
		       "%s: Unable to write mesh for these settings.\n",
		       *argv);
	ok = 0;
      }
      else if((tbl = (unsigned int *)
		     AlcMalloc(sizeof(unsigned int) * (n + 1) * np)) == NULL)
      {
	(void )fprintf(stderr,
		       "%s: Unable to allocate table.\n",
		       *argv);
	ok = 0;
      }
      else
      {
	int	   i,
		     x,
		     y,
		     z,
		     zs;
	unsigned int *t;

	(void )printf("# vtk DataFile Version 1.0\n"
		      "%s -n %d -o %d\n"
		      "ASCII\n"
		      "DATASET UNSTRUCTURED_GRID\n"
		      "POINTS %d float\n",
		      *argv, n, o, np);
	i = 0;
	t = tbl;
	zs = (n == 3)? m: 1;
	for(z = 0; z < zs; ++z)
	{
	  for(y = 0; y < m; ++y)
	  {
	    for(x = 0; x < m; ++x)
	    {
	      (void )printf("%d %d %d\n", x, y, z);
	      p[0] = x;
	      p[1] = y;
	      if(n == 3)
	      {
		p[2] = z;
	      }
	      *t++ = i++;
	      if(flgLin == 0)
	      {
	        AlgHilbertIndex(h, p, n, o);
		*t++ = h[0];
		*t++ = h[1];
		if(n == 3)
		{
		  *t++ = h[2];
		}
	      }
	      else
	      {
	        *t++ = p[0];
		*t++ = p[1];
		if(n == 3)
		{
		  *t++ = p[2];
		}
	      }
	    }
	  }
	}
	if(n == 3)
	{
	  qsort(tbl, np, sizeof(unsigned int) * 4, AlgTstSort3);
	}
	else
	{
	  qsort(tbl, np, sizeof(unsigned int) * 3, AlgTstSort2);
	}
	(void )printf("CELLS %d %d\n", np - 1, 4 * (np - 1));
	t = tbl;
	for(i = 1; i < np; ++i)
	{
	  (void )printf("3 %d %d %d\n",
			t[0], t[n + 1], t[0]);
	  t += n + 1;
	}
	(void )printf("CELL_TYPES %d\n", np - 1);
	for(i = 1; i < np; ++i)
	{
	  (void )printf("5\n");
	}
      }
      AlcFree(tbl);
    }
    else
    {
      while(fgets(buf, bufMax, stdin) != NULL)
      {
	buf[bufMax - 1] = '\0';
	if(flgLin == 0)
	{
	  if(flgInv == 0)
	  {
	    if((ok = AlcTstGetNUInt(p, buf, n)) != 0)
	    {
	      v = h;
	      AlgHilbertIndex(h, p, n, o);
	    }
	  }
	  else
	  {
	    if((ok = AlcTstGetNUInt(h, buf, n)) != 0)
	    {
	      v = p;
	      AlgHilbertIndexInv(h, p, n, o);
	    }
	  }
	}
	if(ok)
	{
	  for(i = 0; i < n; ++i)
	  {
	    (void )printf("% 8d ", v[i]);
	  }
	  (void )printf("\n");
	}
	else
	{
	  (void )fprintf(stderr,
			 "%s: Failed to parse %d unsigned ints from input.\n",
			 *argv, n);
	}
      }
    }
  }
  return(!ok);
}

static int	AlcTstGetNUInt(unsigned int *v, char *s, int n)
{
  int		i = 0,
  		r = 0;
  char		*t;

  t = strtok(s, ",  \t");
  while((t != NULL) && (i < n))
  {
    int		u;

    if((sscanf(t, "%d", &u) == 1) || (u >= 0))
    {
      v[i++] = (unsigned int)u;
    }
    t = strtok(NULL, ",  \t");
  }
  r = (i == n);
  return(r);
}

static int	AlgTstSort2(const void *v0, const void *v1)
{
  int 		cmp;
  unsigned int	*p0,
  		*p1;

  p0 = (unsigned int *)v0;
  p1 = (unsigned int *)v1;
  if((cmp = p0[2] - p1[2]) == 0)
  {
    cmp = p0[1] - p1[1];
  }
  return(cmp);
}

static int	AlgTstSort3(const void *v0, const void *v1)
{
  int 		cmp;
  unsigned int	*p0,
  		*p1;

  p0 = (unsigned int *)v0;
  p1 = (unsigned int *)v1;
  if((cmp = p0[3] - p1[3]) == 0)
  {
    if((cmp = p0[2] - p1[2]) == 0)
    {
      cmp = p0[1] - p1[1];
    }
  }
  return(cmp);
}
