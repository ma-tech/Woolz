#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgVectorMath_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgVectorMath.c
* \author       Bill Hill
* \date         March 2003
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
* \brief	Basic vector arithmatic
* \ingroup      AlgVector
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

/*!
* \return	Norm of the vector.
* \ingroup	AlgVector
* \brief	Computes the norm of the given vector \f$\mathbf{a}\f$
*		with n elements:
*		\f[
		  \|\mathbf{a}\| = \sqrt{\mathbf{a} \cdot \mathbf{a}}.
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Vector size is limited only by address space.
* \param	aV			Given vector \f$\mathbf{a}\f$,
* \param	n			Number of elements in \f$\mathbf{a}\f$.
*/
double		AlgVectorNorm(double *aV, size_t n)
{
  double 	norm;

  norm = AlgVectorDot(aV, aV, n);
  if(norm > DBL_EPSILON)
  {
    norm = sqrt(norm);
  }
  return(norm);
}

/*!
* \return	Scalar (dot) product of the vectors.
* \ingroup	AlgVector
* \brief	Computes the scalar (dot) product of the two vectors
*		\f$\mathbf{a}\f$ and \f$\mathbf{b}\f$ each with n elements:
*		\f[
		  \mathbf{a} \cdot \mathbf{b} = \sum_{i = 0}^{n - 1}{a_i b_i}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Vector size is limited only by address space.
* \param	aV			Vector \f$\mathbf{a}\f$.
* \param	bV			Vector \f$\mathbf{b}\f$.
* \param	n			Number of elements in each of
*					the vectors.
*/
double		AlgVectorDot(double *aV, double *bV, size_t n)
{
  size_t	id0;
  double	dot = 0.0;

  for(id0 = 0; id0 < n; ++id0)
  {
    dot += *aV++ * *bV++;
  }
  return(dot);
}

/*!
* \return	void
* \ingroup	AlgVector
* \brief	Adds vector \f$\mathbf{b}\f$ to vector \f$\mathbf{c}\f$.
*		Computes \f$a_i = b_i + c_i, \forall i \in [0 \ldots n - 1]\f$.
* \note		For efficiency the given parameters are not checked.
* \note		Vector size is limited only by address space.
* \param	aV			Vector \f$\mathbf{a}\f$.
* \param	bV			Vector \f$\mathbf{b}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	n			Number of elements in each of
* 					the vectors.
*/
void		AlgVectorAdd(double *aV, double *bV, double *cV, size_t n)
{
  size_t	id0;

  for(id0 = 0; id0 < n; ++id0)
  {
    *aV++ = *bV++ + *cV++;
  }
}

/*!
* \return	void
* \ingroup	AlgVector
* \brief	Subtracts vector \f$\mathbf{c}\f$ from vector \f$\mathbf{b}\f$.
*		Computes \f$a_i = b_i - c_i, \forall i \in [0 \ldots n - 1]\f$.
* \note		For efficiency the given parameters are not checked.
* \note		Vector size is limited only by address space.
* \param	aV			Vector \f$\mathbf{a}\f$.
* \param	bV			Vector \f$\mathbf{b}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	n			Number of elements in each of
* 					the vectors.
*/
void		AlgVectorSub(double *aV, double *bV, double *cV, size_t n)
{
  size_t	id0;

  for(id0 = 0; id0 < n; ++id0)
  {
    *aV++ = *bV++ - *cV++;
  }
}

/*!
* \return	void
* \ingroup	AlgVector
* \brief	Copies one vector \f$\mathbf{b}\f$ to vector \f$\mathbf{a}\f$.
*		\f$a_i = b_i, \forall i \in [0 \ldots n - 1]\f$.
* \note		For efficiency the given parameters are not checked.
* \note		Vector size is limited only by address space.
* \param	aV			Vector \f$\mathbf{a}\f$.
* \param	bV			Vector \f$\mathbf{b}\f$.
* \param	n			Number of elements in each of
* 					the vectors.
*/
void		AlgVectorCopy(double *aV, double *bV, size_t n)
{
  memcpy(aV, bV, sizeof(double) * n);
}

/*!
* \return
* \brief	Scales a vector \f$\mathbf{b}\f$ and then adds vector
* 		\f$\mathbf{c}\f$.
*		Computes \f$a_i = b_i s + c_i, \forall i \in [0\ldots n - 1]\f$.
* \note		For efficiency the given parameters are not checked.
* \note		Vector size is limited only by address space.
* \param	aV			Vector \f$\mathbf{a}\f$.
* \param	bV			Vector \f$\mathbf{b}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	s			Scalar scale  \f$s\f$.
* \param	n			Number of elements in each of
*					the vectors.
*/
void		AlgVectorScaleAdd(double *aV, double *bV, double *cV,
				  double s, size_t n)
{
#ifdef _OPENMP
  int		id0,
  		oN;
#else
  size_t	id0;
#endif

#ifdef _OPENMP
  oN = n;
  #pragma omp parallel for default(shared)
  for(id0 = 0; id0 < oN; ++id0)
#else
  for(id0 = 0; id0 < n; ++id0)
#endif
  {
    *(aV + id0) = (*(bV + id0) * s) + *(cV + id0);
  }
}
