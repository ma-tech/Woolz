#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgHilbertIndex_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgHilbertIndex.c
* \author       Bill Hill
* \date         August 2011
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
* \brief	Provides functions for Hilbert indices and their
* 		inverse. See J. K. Lawder "Calculation of Mappings
* 		Between One an n-dimensional Values Using the
* 		Hilbert Space-filling Curve", Birkbeck, University
* 		of London Research Report BBKCS-00-01, 2000.
* 		The code in this technical report has been modifiied
* 		to support n-dimensions and variable bit orders
* 		without recompilation.
* \ingroup	AlgBits
*/
#include <Alg.h>

static unsigned int 		AlgHilbertCalcP(
				  unsigned int *h,
				  int i,
				  int n,
				  int o);
static unsigned int 		AlgHilbertCalcP2(
				  unsigned int S,
				  int n);
static unsigned int 		AlgHilbertCalcJ(
				  unsigned int P,
				  int n);
static unsigned int 		AlgHilbertCalcT(
				  unsigned int P);
static unsigned int 		AlgHilbertCalcTSTT(
				  unsigned int xJ,
				  unsigned int v,
				  int n);

/*!
* \ingroup	AlgBits
* \brief	Computes the Hilbert index of a point in n dimensions.
* \param	h			Hilbert index (set on return).
* \param	p			Position of the n dimensional point.
* \param	n			Number of dimensions.
* \param	o			Order (number of bits for coordinate).
*/
void		AlgHilbertIndex(unsigned int *h, unsigned int *p, int n, int o)
{
  int 		i,
  		j,
		k,
		l;
  unsigned int	S,
		tS,
		T,
		tT,
		J,
		xJ,
		m,
		A = 0U,
		P = 0U,
  		W = 0U;

  i = (o - 1) * n;
  m = 1U << (o  - 1);
  for(j = 0; j < n; ++j)
  {
    h[j] = 0;
    if((p[j] & m) != 0U)
    {
      A |= 1U << (n - j - 1);
    }
  }
  S = tS = A;
  P = AlgHilbertCalcP2(S, n);
  k = i / o;
  l = i % o;
  if(l > (o - n))
  {
    h[k] |= P << l;
    h[k + 1] |= P >> (o - l);
  }
  else
  {
    h[k] |= P << (i - (k * o));
  }
  J = AlgHilbertCalcJ(P, n);
  xJ = J - 1;
  T = AlgHilbertCalcT(P);
  tT = T;
  for(i -= n; i >= 0; i -= n)
  {
    A = 0U;
    m >>= 1;
    for(j = 0; j < n; ++j)
    {
      if((p[j] & m) != 0U)
      {
        A |= 1 << (n - j - 1);
      }
    }
    W ^= tT;
    tS = A ^ W;
    S = AlgHilbertCalcTSTT(xJ, tS, n);
    P = AlgHilbertCalcP2(S, n);
    k = i / o;
    l = i % o;
    if(l > (o - n))
    {
      h[k] |= P << l;
      h[k + 1] |= P >> (o - l);
    }
    else
    {
      h[k] |= P << (i - (k * o));
    }
    if(i > 0)
    {
      T = AlgHilbertCalcT(P);
      tT = AlgHilbertCalcTSTT(xJ, T, n);
      J = AlgHilbertCalcJ(P, n);
      xJ += J - 1;
    }
  }
}

/*!
* \ingroup	AlgBits
* \brief	Computes the coordinates of the point in n dimensions
* 		corresponding to the given Hilbert index.
* \param	h			The given Hilbert index.
* \param	p			Decoded point coordinates (set on
* 					return).
* \param	n			Number of dimensions.
* \param	o			Order (number of bits for coordinate).
*/
void		AlgHilbertIndexInv(unsigned int *h, unsigned int *p,
				   int n, int o)
{
  int		i,
  		j;
  unsigned int	
  		m,
		A,
		S,
		tS,
		T,
		tT,
		J,
		xJ,
		P = 0U,
		W = 0U;

  m = 1U << (o - 1);
  i = (o  - 1) * n;
  P = AlgHilbertCalcP(h, i, n, o);
  J = AlgHilbertCalcJ(P, n);
  xJ = J - 1;
  A = S = tS = P ^ (P >> 1);
  T = AlgHilbertCalcT(P);
  tT = T;
  for(j = n - 1; P > 0; --j)
  {
    p[j] = 0;
    if (P & 1)
    {
      p[j] |= m;
    }
    P >>= 1;
  }
  for(i -= n; i >= 0; i -= n)
  {
    m >>= 1;
    P = AlgHilbertCalcP(h, i, n, o);
    S = P ^ (P >> 1);
    tS = AlgHilbertCalcTSTT(xJ, S, n);
    W ^= tT;
    A = W ^ tS;
    for (j = n - 1; A > 0; --j)
    {
      if(A & 1)
      {
	p[j] |= m;
      }
      A >>= 1;
    }
    if(i > 0)
    {
      T = AlgHilbertCalcT(P);
      tT = AlgHilbertCalcTSTT(xJ, T, n);
      J = AlgHilbertCalcJ(P, n);
      xJ += J - 1;
    }
  }
}


/*!
* \ingroup	AlgBits
* \brief	Computes the coordinates of the point in n dimensions
* 		corresponding to the given Hilbert index.
* \param	h			The given Hilbert index.
* \param	i			
* \param	n			Number of dimensions.
* \param	o			Order (number of bits for coordinate).
*/
static unsigned int AlgHilbertCalcP(unsigned int *h, int i, int n, int o)
{
  int 		k,
  		l;
  unsigned int 	P,
  		t;

  k = i / o;
  l = i % o;
  P = h[k];
  if(l > (o - n))
  {
    t = h[k + 1];
    P >>= l;
    t <<= o - l;
    P |= t;
  }
  else
  {
    P >>= l;
  }
  if(n < o)
  {
    P &= (1 << n) -1;
  }
  return(P);
}

static unsigned int AlgHilbertCalcP2(unsigned int S, int n)
{
  int 		i;
  unsigned int	g,
  		P;

  P = S & (1U << (n - 1));
  for(i = 1; i < n; ++i)
  {
    g = 1U << (n - i - 1);
    if(((S & g) ^ ((P >> 1) & g)) != 0)
    {
      P |= g;
    }
  }
  return(P);
}

static unsigned int AlgHilbertCalcJ(unsigned int P, int n)
{
  int 		i;
  unsigned int	J;

  J = n;
  for(i = 1; i < n; ++i)
  {
    if((P >> i & 1) == (P & 1))
    {
      continue;
    }
    else
    {
      break;
    }
  }
  if(i != n)
  {
    J -= i;
  }
  return(J);
}

static unsigned int AlgHilbertCalcT(unsigned int P)
{
  unsigned int	T = 0U;

  if(P >= 3U)
  {
    if((P & 1U) != 0)
    {
      P = P - 1U;
    }
    else
    {
      P = P - 2U;
    }
    T = P ^ (P >> 1);
  }
  return(T);
}

static unsigned int AlgHilbertCalcTSTT(unsigned int xJ, unsigned int v, int n)
{
  unsigned int	r,
		t;

  r = v;
  t = xJ % n;
  if(t != 0)
  {
    r = ((v >> t) | (v << (n - t))) & ((1U << n) - 1);
  }
  return(r);
}
