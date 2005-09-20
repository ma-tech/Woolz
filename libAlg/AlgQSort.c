#pragma ident "MRC HGU $Id$"
/*!
* \file         libAlg/AlgQSort.c
* \author       Bill Hill
* \date         August 2004
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Specialized implementation of quick sort based on
*		"Engineering a Sort Function" J.L. Bentley and M.D. McIlroy,
*		Software Practice and Experience 23 (1993) 1249-1265.
* \ingroup	AlgSort
* \todo         -
* \bug          None known.
*/
#include <stddef.h>
#include <Alg.h>


#ifndef DOXYGEN_SKIP_THIS

#define ALG_QSORT_SWAPCODE(TYPE, PI, PJ, N)\
{  \
  size_t	st; \
  TYPE 		t, \
  		*pi = (TYPE *)(PI), \
  		*pj = (TYPE *)(PJ); \
  \
  st = sizeof(TYPE); \
  do \
  { \
    t = *pi; \
    *pi++ = *pj; \
    *pj++ = t; \
  } while((N -= st) > 0); \
}

#define ALG_QSORT_SWAP(A, B) \
  if(swaptype == 0) \
  { \
    t = *(long *)(A); \
    *(long *)(A) = *(long *)(B); \
    *(long *)(B) = t; \
  } \
  else \
  { \
    AlgQSortSwapFn(A, B, es, swaptype); \
  }


#define ALG_QSORT_PVINIT(PV, PM) \
  if(swaptype != 0) \
  { \
    PV = a; \
    ALG_QSORT_SWAP(PV, PM); \
  } \
  else \
  { \
    PV = (char *)&v; \
    *(long *)PV = *(long *)PM; \
  }

#define ALG_QSORT_VECSWAP(A, B, N) \
  if(N > 0) \
  { \
    AlgQSortSwapFn(A, B, N, swaptype); \
  }

#endif /* DOXYGEN_SKIP_THIS */

static void			AlgQSortPrv(
				  char *a,
				  size_t n,
				  size_t es,
				  void *cData,
				  int (*cmpFn)(
				    const void *,
				    const void *,
				    const void *));
static void			AlgQSortSwapFn(
				  char *a,
				  char *b,
				  size_t n,
				  int swaptype);
static char			*AlgQSortMed3(
				  char *a,
				  char *b,
				  char *c,
				  void * d,
				  int (*cmpFn)(
				    const void *,
				    const void *,
				    const void *));

/*!
* \return	void
* \ingroup	AlgSort
* \brief	A qsort implementation which allows client data to be passed
*		to the sort function.
* \param	base			Array of elements to sort.
* \param	nElm			Number of elements in array.
* \param	elmSz			Element size.
* \param	cData			Client data passed to the client
*					comparison function.
* \param	(*cmpFn)		Comparison function which is given
					the client data followed by the
					two data pointers for comparison.
*/
void		AlgQSort(void *base, size_t nElm, size_t elmSz,
			 void *cData,
		       int (*cmpFn)(const void *, const void *, const void *))
{
  if(base && (nElm > 0) && cmpFn)
  {
    AlgQSortPrv((char *)base, nElm, elmSz, cData, cmpFn);
  }
}

/*!
* \return	void
* \ingroup	AlgSort
* \brief	Private function which implements qsort() without need for
*		parameter checking.
* \param	a			Array of elements to sort.
* \param	n			Number of elements in array.
* \param	es			Element size.
* \param	cData			Client data passed to the client
* 					comparison function.
* \param	cmpFn			Comparison function.
*/
static void	AlgQSortPrv(char *a, size_t n, size_t es, void *cData,
			int (*cmpFn)(const void *, const void *, const void *))
{
  char *pa, *pb, *pc, *pd, *pl, *pm, *pn, *pv;
  int r, swaptype;
  long t, v;
  size_t s;

  if((a - (char *)0 | es) % sizeof(long))
  {
    swaptype = 2;
  }
  else
  {
    swaptype = es > sizeof(long);
  }
  if(n < 7)
  {
    /* Do an insertion sort for small n. */
    for(pm = a + es; pm < a + n*es; pm += es)
    {
      for(pl = pm; (pl > a) && (cmpFn(cData, pl - es, pl) > 0); (pl -= es))
      {
	ALG_QSORT_SWAP(pl, pl - es)
      }
    }
  }
  else
  {
    pm = a + (n / 2) * es;    		/* Small arrays, middle element. */
    if(n > 7)
    {
      pl = a;
      pn = a + (n - 1) * es;
      if(n > 40)
      {
        /* Big arrays, pseudomedian of 9. */
	s = (n / 8) * es;
	pl = AlgQSortMed3(pl, pl + s, pl + (2 * s), cData, cmpFn);
	pm = AlgQSortMed3(pm - s, pm, pm + s, cData, cmpFn);
	pn = AlgQSortMed3(pn - (2 * s), pn - s, pn, cData, cmpFn);
      }
      pm = AlgQSortMed3(pl, pm, pn, cData, cmpFn); /* Mid-size, med of 3 */
    }
    ALG_QSORT_PVINIT(pv, pm);       	/* pv points to partition value */
    pa = pb = a;
    pc = pd = a + (n - 1) * es;
    for (;;)
    {
      while(pb <= pc && (r = cmpFn(cData, pb, pv)) <= 0)
      {
	if(r == 0)
	{
	  ALG_QSORT_SWAP(pa, pb);
	  pa += es;
	}
	pb += es;
      }
      while((pb <= pc) && ((r = cmpFn(cData, pc, pv)) >= 0))
      {
	if(r == 0)
	{
	  ALG_QSORT_SWAP(pc, pd);
	  pd -= es;
	}
	pc -= es;
      }
      if(pb <= pc)
      {
	ALG_QSORT_SWAP(pb, pc);
	pb += es;
	pc -= es;
      }
      else
      {
        break;
      }
    }
    pn = a + (n * es);
    s = ALG_MIN(pa-a, pb - pa);
    ALG_QSORT_VECSWAP(a,  pb - s, s);
    s = ALG_MIN(pd-pc, pn - pd - es);
    ALG_QSORT_VECSWAP(pb, pn - s, s);
    if((s = pb - pa) > es)
    {
      AlgQSortPrv(a, s / es, es, cData, cmpFn);
    }
    if((s = pd - pc) > es)
    {
      AlgQSortPrv(pn - s, s / es, es, cData, cmpFn);
    }
  }
}

/*!
* \return	void
* \ingroup	AlgSort
* \brief	Element swaping function.
* \param	a			First element param.
* \param	b			Second element pointer.
* \param	n			Number of elements.
* \param	swaptype		Control swaping by bytes or longs.
*/
static void	AlgQSortSwapFn(char *a, char *b, size_t n, int swaptype)

{ if(swaptype <= 1)
  {
    ALG_QSORT_SWAPCODE(long, a, b, n)
  }
  else
  {
    ALG_QSORT_SWAPCODE(char, a, b, n)
  }
}

/*!
* \return	Pointer to median value.
* \ingroup	AlgSort
* \brief	Finds the median of three values
* \param	a			First value.
* \param	b			Second value.
* \param	c			Third value.
* \param	d			Client data for comparison function.
* \param	(*cmpFn)		Comparison function.
*/
static char	*AlgQSortMed3(char *a, char *b, char *c, void * d,
			int (*cmpFn)(const void *, const void *, const void *))

{
  char		*cmp;

  if(cmpFn(d, a, b) < 0)
  {
    if(cmpFn(d, b, c) < 0)
    {
      cmp = b;
    }
    else
    {
      if(cmpFn(d, a, c) < 0)
      {
        cmp = c;
      }
      else
      {
        cmp = a;
      }
    }
  }
  else
  {
    if(cmpFn(d, b, c) > 0)
    {
      cmp = b;
    }
    else
    {
      if(cmpFn(d, a, c) > 0)
      {
        cmp = c;
      }
      else
      {
        cmp = a;
      }
    }
  }
  return(cmp);
}

#ifndef DOXYGEN_SKIP_THIS

#ifdef ALG_QSORT_TEST

static int			AlgQSortTestSortFn(
				  const void *cData,
				  const void *a,
				  const void *b);
extern char	*optarg;
extern int	optind,
		optind,
		opterr;

int		main(int argc, char *argv[])
{
  int		idN,
  		option,
  		nElm = 20,
  		elmSz,
		ok = 1,
		reverse = 0,
		usage = 0,
		verbose = 0;
  int		*ary = NULL;
  static char	optList[] = "hrvn:";

  opterr = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
	if(sscanf(optarg, "%d", &nElm) != 1)
	{
	  usage = 1;
	}
	break;
      case 'r':
	reverse = 1;
	break;
      case 'v':
	verbose = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok)
  {
    elmSz = sizeof(int);
    if((ary = AlcMalloc(elmSz * nElm)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
      	             "%s: Filed to allocate data array.\n",
		     *argv);
    }
  }
  if(ok)
  {
    for(idN = 0; idN < nElm; ++idN)
    {
      ary[idN] = (int )(10.0 * nElm * AlgRandUniform());
    }
    if(verbose)
    {
      for(idN = 0; idN < nElm; ++idN)
      {
        (void )printf("% 8d % 8d\n", idN, ary[idN]);
      }
    }
    (void )printf("\n");
    AlgQSort(ary, nElm, elmSz, (void *)reverse, AlgQSortTestSortFn);
    if(verbose)
    {
      for(idN = 0; idN < nElm; ++idN)
      {
        (void )printf("% 8d % 8d\n", idN, ary[idN]);
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    		   "Usage: %s [-ch [-v] [-n#]\n%s",
		   "Options:\n"
		   "  -h  Show this usage message.\n"
		   "  -r  Reverse order.\n"
		   "  -v  Be verbose.\n"
		   "  -n  Number of data to sort.\n"
		   "Test for AlgQSort().\n",
		   *argv);
  }
  return(!ok);
}

static int	AlgQSortTestSortFn(const void *cData,
				   const void *a, const void *b)
{
  int		cmp;

  if(cData)
  {
    cmp = *(int *)b - *(int *)a;
  }
  else
  {
    cmp = *(int *)a - *(int *)b;
  }
  return(cmp);
}

#endif /* ALG_QSORT_TEST */

#endif /* DOXYGEN_SKIP_THIS */
