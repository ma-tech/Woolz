#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixMath_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixMath.c
* \author       Bill Hill
* \date         June 2001
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
* \brief        Functions for basic arithmatic with matricies.
* \ingroup      AlgMatrix
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Computes the sum of two matricies and returns
*		the result in a third supplied matrix:
*		\f[
		  \mathbf{A} = \mathbf{B} + \mathbf{C}
		\f]
*		The dimensions of the all three matricies must be
*		nR, nC. It is safe to supply the same matrix as any
*		combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \note		All matrices must be of the same type (not checked for).
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			First matrix in the sum,
*					\f$\mathbf{B}\f$.
* \param        cM			Second matrix in the sum,
*					\f$\mathbf{C}\f$.
*/
void		AlgMatrixAdd(AlgMatrix aM, AlgMatrix bM, AlgMatrix cM)
{
  size_t  	id0,
  		nR,
		nC;

  nR = aM.core->nR;
  nC = aM.core->nC;
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t  id1;
	  double  *aRow,
		  *bRow,
		  *cRow;

	  aRow = aM.rect->array[id0];
	  bRow = bM.rect->array[id0];
	  cRow = cM.rect->array[id0];
	  for(id1 = 0; id1 < nC; ++id1)
	  {
	    *aRow++ = *bRow++ + *cRow++;
	  }
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t  id1;
	  double  **aA,
		  **bA,
		  **cA;

	  aA = aM.sym->array;
	  bA = bM.sym->array;
	  cA = cM.sym->array;
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    aA[id0][id1] = bA[id0][id1] + cA[id0][id1];
	  }
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	/* Check to see if A is an alias of B or C and if not zero it's
	 * entries. */
	if((aM.core != bM.core)  && (aM.core != cM.core))
	{
	  AlgMatrixLLRZero(aM.llr);
	}
	for(id0 = 0; id0 < nR; ++id0)
	{
	  AlgMatrixLLRE aE;
	  AlgMatrixLLRE *bE,
	  		*cE;

	  bE = bM.llr->tbl[id0];
	  cE = cM.llr->tbl[id0];
	  while((bE != NULL) || (cE != NULL))
	  {
	    if(bE == NULL)
	    {
	      aE.col = cE->col; aE.val = cE->val;
	      cE = cE->nxt;
	    }
	    else if(cE == NULL)
	    {
	      aE.col = bE->col; aE.val = bE->val;
	      bE = bE->nxt;
	    }
	    else
	    {
	      if(bE->col < cE->col)
	      {
	        aE.col = bE->col; aE.val = bE->val;
		bE = bE->nxt;
	      }
	      else if(bE->col > cE->col)
	      {
	        aE.col = cE->col; aE.val = cE->val;
		cE = cE->nxt;
	      }
	      else
	      {
	        aE.col = bE->col; aE.val = bE->val + cE->val;
		bE = bE->nxt;
		cE = cE->nxt;
	      }
	    }
	    (void )AlgMatrixLLRSet(aM.llr, id0, aE.col, aE.val);
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Subtracts on matrix from another and returns the
*		result in a third supplied matrix:
*		\f[
		  \mathbf{A} = \mathbf{B} - \mathbf{C}
		\f]
*		The dimensions of the all three
*		matricies must be nR, nC. It is safe to supply
*		the same matrix as any combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			First matrix in the subtraction,
					\f$\mathbf{B}\f$.
* \param        cM			Second matrix in the subtraction,
					\f$\mathbf{C}\f$.
*/
void		AlgMatrixSub(AlgMatrix aM, AlgMatrix bM, AlgMatrix cM)
{
  size_t  	id0,
  		nR,
		nC;

  nR = aM.core->nR;
  nC = aM.core->nC;
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t  id1;
	  double  *aRow,
		  *bRow,
		  *cRow;

	  aRow = aM.rect->array[id0];
	  bRow = bM.rect->array[id0];
	  cRow = cM.rect->array[id0];
	  for(id1 = 0; id1 < nC; ++id1)
	  {
	    *aRow++ = *bRow++ - *cRow++;
	  }
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t  id1;
	  double  **aA,
		  **bA,
		  **cA;

	  aA = aM.sym->array;
	  bA = bM.sym->array;
	  cA = cM.sym->array;
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    aA[id0][id1] = bA[id0][id1] - cA[id0][id1];
	  }
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	/* Check to see if A is an alias of B or C and if not zero it's
	 * entries. */
	if((aM.core != bM.core)  && (aM.core != cM.core))
	{
	  AlgMatrixLLRZero(aM.llr);
	}
	for(id0 = 0; id0 < nR; ++id0)
	{
	  AlgMatrixLLRE aE;
	  AlgMatrixLLRE *bE,
	  		*cE;

	  bE = bM.llr->tbl[id0];
	  cE = cM.llr->tbl[id0];
	  while((bE != NULL) || (cE != NULL))
	  {
	    if(bE == NULL)
	    {
	      aE.col = cE->col; aE.val = -cE->val;
	      cE = cE->nxt;
	    }
	    else if(cE == NULL)
	    {
	      aE.col = bE->col; aE.val = bE->val;
	      bE = bE->nxt;
	    }
	    else
	    {
	      if(bE->col < cE->col)
	      {
	        aE.col = bE->col; aE.val = bE->val;
		bE = bE->nxt;
	      }
	      else if(bE->col > cE->col)
	      {
	        aE.col = cE->col; aE.val = -cE->val;
		cE = cE->nxt;
	      }
	      else
	      {
	        aE.col = bE->col; aE.val = bE->val - cE->val;
		bE = bE->nxt;
		cE = cE->nxt;
	      }
	    }
	    (void )AlgMatrixLLRSet(aM.llr, id0, aE.col, aE.val);
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Computes the product of two matricies and returns
*		the result in a third supplied matrix:
*               \f[
                  \mathbf{A} = \mathbf{B} \mathbf{C}
		\f]
*		The dimensions of the result matrix (aM) must be cR, bC.
* 		All the matrices must be valid and of the same type except
* 		in the case of multiplying tow symmetric matrices when
* 		the result matrix must be rectangular (because in general
* 		the product of two symmetric matrices is not symmetric).
* \note		For efficiency the given parameters are not fully checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$
* \param        bM 			First matrix in the product,
*					\f$\mathbf{B}\f$
* \param        cM			Second matrix in the product,
*					\f$\mathbf{C}\f$
*/
void		AlgMatrixMul(AlgMatrix aM, AlgMatrix bM, AlgMatrix cM)
{
  size_t  	id0,
  		nBR,
		nBC,
		nCC;
  AlgError	errNum = ALG_ERR_NONE;

  nBR = bM.core->nR;
  nBC = bM.core->nC;
  nCC = cM.core->nC;
  switch(bM.core->type)
  {
    case ALG_MATRIX_RECT:
      if((aM.core->type != bM.core->type) || (cM.core->type != bM.core->type))
      {
        errNum = ALG_ERR_MATRIX_TYPE;
      }
      else
      {
        double	**aA,
		**bA,
		**cA;

        aA = aM.rect->array;
        bA = bM.rect->array;
        cA = cM.rect->array;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nBR; ++id0)
	{
	  size_t  id1;
	  double  *bRow;
	  
	  bRow = bA[id0];
	  for(id1 = 0; id1 < nCC; ++id1)
	  {
	    size_t id2;
	    double v = 0.0;

	    for(id2 = 0; id2 < nBC; ++id2)
	    {
	      v += bRow[id2] * cA[id2][id1];
	    }
	    aA[id0][id1] = v;
	  }
	}
      }
      break;
    case ALG_MATRIX_SYM:
      if((aM.core->type != ALG_MATRIX_RECT) || (cM.core->type != bM.core->type))
      {
        errNum = ALG_ERR_MATRIX_TYPE;
      }
      else
      {
        double	**aA,
		**bA,
		**cA;

        aA = aM.rect->array;
        bA = bM.rect->array;
        cA = cM.rect->array;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nBR; ++id0)
	{
	  size_t  id1;
	  
	  for(id1 = 0; id1 < nCC; ++id1)
	  {
	    size_t id2;
	    double v = 0.0;

	    id2 = 0;
	    while(id2 <= id0)
	    {
	      if(id1 <= id2)
	      {
	        v += bA[id0][id2] * cA[id2][id1];
	      }
	      else
	      {
	        v += bA[id0][id2] * cA[id1][id2];
	      }
	      ++id2;
	    }
	    while(id2 < nBC)
	    {
	      if(id1 <= id2)
	      {
	        v += bA[id2][id0] * cA[id2][id1];
	      }
	      else
	      {
	        v += bA[id2][id0] * cA[id1][id2];
	      }
	      ++id2;
	    }
	    aA[id0][id1] = v;
	  }
	}
      }
      break;
    case ALG_MATRIX_LLR:
      if((aM.core->type != bM.core->type) || (cM.core->type != bM.core->type))
      {
        errNum = ALG_ERR_MATRIX_TYPE;
      }
      else
      {
	size_t id0;

	for(id0 = 0; id0 < nBR; ++id0)
	{
	  if(errNum == ALG_ERR_NONE)
	  {
	    size_t id1;
	    AlgMatrixLLRE *bRow;

	    bRow = bM.llr->tbl[id0];
	    for(id1 = 0; id1 < nCC; ++id1)
	    {
	      double v = 0.0;
	      AlgMatrixLLRE *bEnt = bRow;

	      while(bEnt != NULL)
	      {
		v += bEnt->val * AlgMatrixLLRValue(cM.llr, bEnt->col, id1);
		bEnt = bEnt->nxt;
	      }
	      errNum = AlgMatrixLLRSet(aM.llr, id0, id1, v);
	      if(errNum != ALG_ERR_NONE)
	      {
		break;
	      }
	    }
	  }
	}
      }
      break;
    default:
      errNum = ALG_ERR_MATRIX_TYPE;
      break;
  }
}

/*!
* \return       The trace of the matrix.
* \ingroup      AlgMatrix
* \brief        Computes the trace of the given matrix.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix.
*/
double		AlgMatrixTrace(AlgMatrix aM)
{
  size_t	id0,
  		nN;
  double	trace = 0.0;

  nN = ALG_MIN(aM.core->nR, aM.core->nC);
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
      for(id0 = 0; id0 < nN; ++id0)
      {
	trace += aM.rect->array[id0][id0];
      }
      break;
    case ALG_MATRIX_SYM:
      for(id0 = 0; id0 < nN; ++id0)
      {
	trace += aM.sym->array[id0][id0];
      }
      break;
    case ALG_MATRIX_LLR:
      for(id0 = 0; id0 < nN; ++id0)
      {
	trace += AlgMatrixLLRValue(aM.llr, id0, id0);
      }
      break;
    default:
      break;
  }
  return(trace);
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Computes the transpose of the given matrix:
*		\f[
		\mathbf{A} = \mathbf{B^T}
		\f]
*		aM = transpose(bM). The dimensions of the result
*		matrix (aM) must be aR == bC, aC == bR.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Matrix to transpose,
*					\f$\mathbf{B}\f$.
*/
void            AlgMatrixTranspose(AlgMatrix aM, AlgMatrix bM)
{
  size_t	nR,
  		nC;

  nR = bM.core->nR;
  nC = bM.core->nC;
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
	size_t	id0,
	        id1;
        double	**aA,
		**bA;

        aA = aM.rect->array;
	bA = bM.rect->array;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0, id1)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  for(id1 = 0; id1 < nC; ++id1)
	  {
	    aA[id0][id1] = bA[id1][id0];
	  }
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
	size_t	id0,
	        id1;
        double	**aA,
		**bA;

        aA = aM.sym->array;
	bA = bM.sym->array;
	/* Transpose is just a copy of the values for a symetric matrix. */
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0, id1)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    aA[id0][id1] = bA[id0][id1];
	  }
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	size_t	id0;

	AlgMatrixLLRZero(aM.llr);
	for(id0 = 0; id0 < nR; ++id0)
	{
	  AlgMatrixLLRE *p;

	  p = bM.llr->tbl[id0];
	  while(p != NULL)
	  {
	    AlgMatrixLLRSet(aM.llr, id0, p->col, p->val);
	    p = p->nxt;
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Copies the values of the matrix bM to the result
*		matric aM:
*		\f[
		\mathbf{A} = \mathbf{B}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Matrix to copy, \f$\mathbf{B}\f$.
*/
void            AlgMatrixCopy(AlgMatrix aM, AlgMatrix bM)
{
  size_t	nR,
  		nC;

  nR = bM.core->nR;
  nC = bM.core->nC;
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
	size_t	id0;

#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
          AlgVectorCopy(aM.rect->array[id0], bM.rect->array[id0], nC);
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
	size_t	id0;

#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
          AlgVectorCopy(aM.sym->array[id0], bM.sym->array[id0], id0);
	}
      }
      break;
    case ALG_MATRIX_LLR:
      (void )AlgMatrixLLRCopyInPlace(aM.llr, bM.llr);
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief        Multiplies the given matrix by the given scalar:
*		\f[
		\mathbf{A} = s \mathbf{B}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Given matrix to scale,
*					\f$\mathbf{B}\f$.
* \param	sv			Scalar value, \f$s\f$.
*/
void		AlgMatrixScale(AlgMatrix aM, AlgMatrix bM, double sv)
{
  size_t	id0,
  		nR,
		nC;

  nR = aM.core->nR;
  nC = aM.core->nC;
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
      for(id0 = 0; id0 < nR; ++id0)
      {
	AlgVectorScale(aM.rect->array[id0], bM.rect->array[id0], sv, nC);
      }
      break;
    case ALG_MATRIX_SYM:
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
      for(id0 = 0; id0 < nR; ++id0)
      {
	AlgVectorScale(aM.sym->array[id0], bM.sym->array[id0], sv, id0 + 1);
      }
      break;
    case ALG_MATRIX_LLR:
      if(aM.llr == bM.llr)
      {
	if(fabs(sv) < 1.0)
	{
	  for(id0 = 0; id0 < nR; ++id0)
	  {
	    AlgMatrixLLRE *p,
	    		  *q;

	    p = aM.llr->tbl[id0];
	    while(p != NULL)
	    {
	      q = p;
	      p = p->nxt;
	      q->val *= sv;
	      if(q->val < aM.llr->tol)
	      {
	        AlgMatrixLLRERemove(aM.llr, id0, q->col);
	      }
	    }
	  }
	}
	else
	{
	  for(id0 = 0; id0 < aM.llr->nR; ++id0)
	  {
	    AlgMatrixLLRE *q;

	    q = aM.llr->tbl[id0];
	    while(q != NULL)
	    {
	      q->val *= sv;
	      q = q->nxt;
	    }
	  }
	}
      }
      else
      {
	for(id0 = 0; id0 < aM.llr->nR; ++id0)
	{
	  AlgMatrixLLRE *q;

	  q = bM.llr->tbl[id0];
	  while(q != NULL)
	  {
	    AlgMatrixLLRSet(aM.llr, id0, q->col, q->val * sv);
	    q = q->nxt;
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Multiplies the a matrix by a scalar and then adds
*		another matrix:
*		\f[
		\mathbf{A} = \mathbf{B} + s \mathbf{C}.
		\f]
*		All the matrices must have the same dimensions.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Given matrix to scale,
*					\f$\mathbf{B}\f$.
* \param	cM			Matrix too add, \f$\mathbf{C}\f$.
* \param	sv			Scalar value, \f$s\f$.
*/
void		AlgMatrixScaleAdd(AlgMatrix aM, AlgMatrix bM, AlgMatrix cM,
				  double sv)
{
  size_t  	id0,
  		nR,
		nC;

  nR = aM.core->nR;
  nC = aM.core->nC;
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t  id1;
	  double  *aRow,
		  *bRow,
		  *cRow;

	  aRow = aM.rect->array[id0];
	  bRow = bM.rect->array[id0];
	  cRow = cM.rect->array[id0];
	  for(id1 = 0; id1 < nC; ++id1)
	  {
	    *aRow++ = *bRow++ + (sv * *cRow++);
	  }
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t  id1;
	  double  **aA,
		  **bA,
		  **cA;

	  aA = aM.sym->array;
	  bA = bM.sym->array;
	  cA = cM.sym->array;
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    aA[id0][id1] = bA[id0][id1] + (sv * cA[id0][id1]);
	  }
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	/* Check to see if A is an alias of B or C and if not zero it's
	 * entries. */
	if((aM.core != bM.core)  && (aM.core != cM.core))
	{
	  AlgMatrixLLRZero(aM.llr);
	}
	for(id0 = 0; id0 < nR; ++id0)
	{
	  AlgMatrixLLRE aE;
	  AlgMatrixLLRE *bE,
	  		*cE;

	  bE = bM.llr->tbl[id0];
	  cE = cM.llr->tbl[id0];
	  while((bE != NULL) || (cE != NULL))
	  {
	    if(bE == NULL)
	    {
	      aE.col = cE->col; aE.val = sv * cE->val;
	      cE = cE->nxt;
	    }
	    else if(cE == NULL)
	    {
	      aE.col = bE->col; aE.val = bE->val;
	      bE = bE->nxt;
	    }
	    else
	    {
	      if(bE->col < cE->col)
	      {
	        aE.col = bE->col; aE.val = bE->val;
		bE = bE->nxt;
	      }
	      else if(bE->col > cE->col)
	      {
	        aE.col = cE->col; aE.val = sv * cE->val;
		cE = cE->nxt;
	      }
	      else
	      {
	        aE.col = bE->col; aE.val = bE->val + (sv * cE->val);
		bE = bE->nxt;
		cE = cE->nxt;
	      }
	    }
	    (void )AlgMatrixLLRSet(aM.llr, id0, aE.col, aE.val);
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief        Sets the elements of the given square matrix so that it
*		is a scalar matrix:
*		\f[
                \mathbf{A} = s \mathbf{I}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		This function assumes that the matrix has been allocated
*		by AlcDouble2Malloc().
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param	sv			Scalar value, \f$s\f$.
*/
void		AlgMatrixScalar(AlgMatrix aM, double sv)
{
  size_t	nN;

  nN = ALG_MIN(aM.core->nR, aM.core->nC);
  switch(aM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
	size_t id0;

	AlgMatrixRectZero(aM.rect);
	for(id0 = 0; id0 < nN; ++id0)
	{
	  aM.rect->array[id0][id0] = sv;
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
	size_t id0;

	AlgMatrixSymZero(aM.sym);
	for(id0 = 0; id0 < nN; ++id0)
	{
	  aM.rect->array[id0][id0] = sv;
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	size_t id0;

	AlgMatrixLLRZero(aM.llr);
	if(AlgMatrixLLRExpand(aM.llr, nN) == ALG_ERR_NONE)
	{
	  for(id0 = 0; id0 < nN; ++id0)
	  {
	    AlgMatrixLLRE *p;

	    p = AlgMatrixLLRENew(aM.llr);
	    p->col = id0;
	    p->val = sv;
	    p->nxt = NULL;
	    aM.llr->tbl[id0] = p;
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the matrix \f$\mathbf{B}\f$ by the vector
*		\f$\mathbf{c}\f$:
*		\f[
		\mathbf{a} = \mathbf{B} \mathbf{c}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixVectorMul(double *aV, AlgMatrix bM, double *cV)
{
  switch(bM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
	size_t id0;

#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.rect->nR; ++id0)
	{
	  size_t id1;
	  double v;
	  double *bRow;

	  v = 0.0;
	  bRow = bM.rect->array[id0];
	  for(id1 = 0; id1 < bM.rect->nC; ++id1)
	  {
	    v += bRow[id1] * cV[id1];
	  }
	  aV[id0] = v;
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
	size_t id0;

#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.sym->nR; ++id0)
	{
	  size_t id1;
	  double v;

	  v = 0.0;
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    v += bM.sym->array[id0][id1] * cV[id1];
	  }
	  for( ; id1 < bM.sym->nR; ++id1)
	  {
	    v += bM.sym->array[id1][id0] * cV[id1];
	  }
	  aV[id0] = v;
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	size_t id0;

#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.llr->nR; ++id0)
	{
	  double	v;
	  AlgMatrixLLRE *p;

	  v = 0.0;
	  p = bM.llr->tbl[id0];
	  while(p != NULL)
	  {
	    v += p->val * cV[p->col];
	    p = p->nxt;
	  }
	  aV[id0] = v;
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the matrix \f$\mathbf{B}\f$ by the vector
*		\f$\mathbf{c}\f$ and adds the vector \f$\mathbf{d}\f$:
*		\f[
		\mathbf{a} = \mathbf{B} \mathbf{c} + \mathbf{d}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	dV			Vector \f$\mathbf{d}\f$.
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixVectorMulAdd(double *aV, AlgMatrix bM,
				   double *cV, double *dV)
{
  switch(bM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
	size_t id0,
		nC,
		nR;
        double	**bA;

        nR = bM.rect->nR;
        nC = bM.rect->nC;
	bA = bM.rect->array;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t id1;
	  double v;
	  double *bRow;

	  v = dV[id0];
	  bRow = bA[id0];
	  for(id1 = 0; id1 < nC; ++id1)
	  {
	    v += bRow[id1] * cV[id1];
	  }
	  aV[id0] = v;
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
	size_t id0,
		nR;
        double	**bA;

	nR = bM.sym->nR;
	bA = bM.sym->array;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  size_t id1;
	  double v;

	  v = dV[id0];
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    v += bA[id0][id1] * cV[id1];
	  }
	  for(id1 = id0 + 1; id1 < nR; ++id1)
	  {
	    v += bA[id1][id0] * cV[id1];
	  }
	  aV[id0] = v;
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	size_t id0,
		nR;

        nR = bM.llr->nR;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < nR; ++id0)
	{
	  double	v;
	  AlgMatrixLLRE *p;

	  v = dV[id0];
	  p = bM.llr->tbl[id0];
	  while(p != NULL)
	  {
	    v += p->val * cV[p->col];
	    p = p->nxt;
	  }
	  aV[id0] = v;
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the matrix \f$\mathbf{B}\f$ by the vector
*		\f$\mathbf{c}\f$ and adds the vector \f$\mathbf{d}\f$
*		using the given weights \f$s\f$ and \f$t\f$:
*		\f[
		\mathbf{a} = s \mathbf{B} \mathbf{c} + t \mathbf{d}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	dV			Vector \f$\mathbf{d}\f$.
* \param	s			First weighting scalar \f$s\f$.
* \param	t			Second weighting scalar \f$t\f$.
*/
void 		AlgMatrixVectorMulWAdd(double *aV, AlgMatrix bM,
				   double *cV, double *dV, double s, double t)
{

  size_t id0;

  switch(bM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.rect->nR; ++id0)
	{
          size_t id1;
	  double v;
	  double *bRow;

	  v = 0.0;
	  bRow = bM.rect->array[id0];
	  for(id1 = 0; id1 < bM.rect->nC; ++id1)
	  {
	    v += bRow[id1] * cV[id1];
	  }
	  aV[id0] = (s * v) + (t * dV[id0]);
	}
	break;
      }
    case ALG_MATRIX_SYM:
      {
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.sym->nR; ++id0)
	{
          size_t id1;
	  double v;

	  v = 0.0;
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    v += bM.sym->array[id0][id1] * cV[id1];
	  }
	  for(id1 = id0 + 1; id1 < bM.sym->nR; ++id1)
	  {
	    v += bM.sym->array[id1][id0] * cV[id1];
	  }
	  aV[id0] = (s * v) + (t * dV[id0]);
	}
	break;
      }
    case ALG_MATRIX_LLR:
      {
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.llr->nR; ++id0)
	{
	  double	v;
	  AlgMatrixLLRE *p;

	  v = 0.0;
	  p = bM.llr->tbl[id0];
	  while(p != NULL)
	  {
	    v += p->val * cV[p->col];
	    p = p->nxt;
	  }
	  aV[id0] = (s * v) + (t * dV[id0]);
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the transpose of matrix \f$\mathbf{B}\f$ by the
* 		vector \f$\mathbf{c}\f$:
*		\f[
		\mathbf{a} = \mathbf{B}^T \mathbf{c}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixTVectorMul(double *aV, AlgMatrix bM, double *cV)
{
  switch(bM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
	size_t id0;

#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.rect->nR; ++id0)
	{
	  size_t id1;

	  aV[id0] = 0.0;
	  for(id1 = 0; id1 < bM.rect->nC; ++id1)
	  {
	    aV[id0] += bM.rect->array[id1][id0] * cV[id1];
	  }
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
	size_t id0;

#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.sym->nR; ++id0)
	{
	  size_t id1;

	  aV[id0] = bM.sym->array[id0][0] * cV[0];
	  for(id1 = 1; id1 <= id0; ++id1)
	  {
	    aV[id0] += bM.sym->array[id0][id1] * cV[id1];
	  }
	  for(id1 = id0 + 1; id1 < bM.sym->nR; ++id1)
	  {
	    aV[id0] += bM.sym->array[id1][id0] * cV[id1];
	  }
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	size_t id0;

	AlgVectorZero(aV, bM.llr->nC);
	for(id0 = 0; id0 < bM.llr->nR; ++id0)
	{
	  double	c;
	  AlgMatrixLLRE *p;

	  c = cV[id0];
	  p = bM.llr->tbl[id0];
	  while(p != NULL)
	  {
	    aV[p->col] += p->val * c;
	    p = p->nxt;
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the transpose of matrix \f$\mathbf{B}\f$ by the
* 		vector \f$\mathbf{c}\f$:
*		\f[
		\mathbf{a} = \mathbf{B}^T \mathbf{c} + \mathbf{d}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	dV			Vector \f$\mathbf{d}\f$.
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixTVectorMulAdd(double *aV, AlgMatrix bM,
				       double *cV, double *dV)
{
  switch(bM.core->type)
  {
    case ALG_MATRIX_RECT:
      {
        size_t id0;

#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.rect->nC; ++id0)
	{
          size_t id1;

	  aV[id0] = dV[id0];
	  for(id1 = 0; id1 < bM.rect->nR; ++id1)
	  {
	    aV[id0] += bM.sym->array[id1][id0] * cV[id1];
	  }
	}
      }
      break;
    case ALG_MATRIX_SYM:
      {
        size_t id0;

#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(id0)
#endif
	for(id0 = 0; id0 < bM.sym->nR; ++id0)
	{
          size_t id1;

	  aV[id0] = dV[id0];
	  for(id1 = 0; id1 <= id0; ++id1)
	  {
	    aV[id0] += bM.sym->array[id0][id1] * cV[id1];
	  }
	  for(id1 = id0 + 1; id1 < bM.sym->nR; ++id1)
	  {
	    aV[id0] += bM.sym->array[id1][id0] * cV[id1];
	  }
	}
      }
      break;
    case ALG_MATRIX_LLR:
      {
	size_t id0;

	AlgVectorCopy(aV, dV, bM.llr->nC);
	for(id0 = 0; id0 < bM.llr->nR; ++id0)
	{
	  double	c;
	  AlgMatrixLLRE *p;

	  c = cV[id0];
	  p = bM.llr->tbl[id0];
	  while(p != NULL)
	  {
	    aV[p->col] += p->val * c;
	    p = p->nxt;
	  }
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \return	Alg error code, either ALG_ERR_NONE or ALG_ERR_MATRIX_SINGULAR.
* \ingroup	AlgMatrix
* \brief	Inverts a raw 2 x 2 matrix. All matrix values are passed
* 		using pointers and their values are set to those of the
* 		inverse matrix on return. If the given matrix is singular
* 		the ALG_ERR_MATRIX_SINGULAR error code is returned and no
* 		values are changed.
* \param	a00			Matrix element row 0, column 0.
* \param	a01			Matrix element row 0, column 1.
* \param	a10			Matrix element row 1, column 0.
* \param	a11			Matrix element row 1, column 1.
*/
AlgError	AlgMatrixRawInv2x2(double *a00, double *a01,
             	                   double *a10, double *a11)
{
  double 	d;
  double	b[4];
  const double	eps = 0.000001;
  AlgError	errNum = ALG_ERR_MATRIX_SINGULAR;

  d = (*a00 * *a11) - (*a01 * *a10);
  if(fabs(d) > eps)
  {
    errNum = ALG_ERR_NONE;
    d = 1.0 / d;
    b[0] = d * *a00;
    b[1] = d * *a01;
    b[2] = d * *a10;
    b[3] = d * *a11;
    *a00 =  b[3];
    *a01 = -b[1];
    *a10 = -b[2];
    *a11 =  b[0];
  }
  return(errNum);
}

/*!
* \return	Alg error code, either ALG_ERR_NONE or ALG_ERR_MATRIX_SINGULAR.
* \ingroup	AlgMatrix
* \brief	Inverts a raw 2 x 2 matrix. All matrix values are passed
* 		using pointers and their values are set to those of the
* 		inverse matrix on return. If the given matrix is singular
* 		the ALG_ERR_MATRIX_SINGULAR error code is returned and no
* 		values are changed.
* \param	a00			Matrix element row 0, column 0.
* \param	a01			Matrix element row 0, column 1.
* \param	a02			Matrix element row 0, column 2.
* \param	a10			Matrix element row 1, column 0.
* \param	a11			Matrix element row 1, column 1.
* \param	a12			Matrix element row 1, column 2.
* \param	a20			Matrix element row 2, column 0.
* \param	a21			Matrix element row 2, column 1.
* \param	a22			Matrix element row 2, column 2.
*/
AlgError	AlgMatrixRawInv3x3(double *a00, double *a01, double *a02,
				   double *a10, double *a11, double *a12,
				   double *a20, double *a21, double *a22)
{
  double	d;
  double	b[9];
  const double	eps = 0.000001;
  AlgError	errNum = ALG_ERR_MATRIX_SINGULAR;
  
  b[0] = (*a11 * *a22) - (*a12 * *a21);
  b[1] = (*a12 * *a20) - (*a10 * *a22);
  b[2] = (*a10 * *a21) - (*a11 * *a20);
  d = (*a00 * b[0]) + (*a01 * b[1]) + (*a02 * b[2]); 
  if(fabs(d) > eps)
  {
    errNum = ALG_ERR_NONE;
    d = 1.0 / d;
    b[0] *= d;
    b[1] *= d;
    b[2] *= d;
    b[3] = d * ((*a02 * *a21) - (*a01 * *a22));
    b[4] = d * ((*a00 * *a22) - (*a02 * *a20));
    b[5] = d * ((*a01 * *a20) - (*a00 * *a21));
    b[6] = d * ((*a01 * *a12) - (*a02 * *a11));
    b[7] = d * ((*a02 * *a10) - (*a00 * *a12));
    b[8] = d * ((*a00 * *a11) - (*a01 * *a10));
    *a00 = b[0]; *a01 = b[3]; *a02 = b[6];
    *a10 = b[1]; *a11 = b[4]; *a12 = b[7];
    *a20 = b[2]; *a21 = b[5]; *a22 = b[8];
  }
  return(errNum);
}
