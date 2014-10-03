#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixCG_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixCG.c
* \author       Bill Hill
* \date         March 2003
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
* \brief	Conjugate Gradient iterative method with preconditioning
*		for the solution of linear systems with the form
*		\f$\mathbf{A} \mathbf{x} = \mathbf{b}\f$.
*		A must be a symmetric postive definite matrix,
*		i.e. \f${\mathbf{x}}^T \mathbf{A} \mathbf{x} < 0\f$,
*		\f$\forall \mathbf{x} \not= \mathbf{0}\f$,
*               \f$\mathbf{x} \in R^n\f$.
*		
*		This code is adapted from CG.f which is described in
*		"Templates for the Solution of Linear Systems: Building Blocks
*		for Iterative Methods", Barrett, Berry, Chan, Demmel, Donato,
*		Dongarra, Eijkhout, Pozo, Romine, and van der Vorst, SIAM
*		Publications, 1993. Also downloadable from
*		http://www.netlib.org/templates/index.html.
* \ingroup	AlgMatrix
*/
#include <Alg.h>
#include <float.h>

#ifdef ALG_MATRIXCG_DEBUG
static void			AlgMatrixCGDebug(
				  FILE *fP,
				  const char *name,
				  double *dat,
				  int nDat);
#endif /* ALG_MATRIXCG_DEBUG */

/*!
* \return	Error code.
* \ingroup	AlgMatrix
* \brief	Conjugate Gradient iterative method with preconditioning
*		for the solution of linear systems with the form
*		\f$\mathbf{A} \mathbf{x} = \mathbf{b}\f$.
*		\f$\mathbf{A}\f$ must be a symmetric postive definite matrix,
*		i.e. \f${\mathbf{x}}^T \mathbf{A} \mathbf{x} < 0\f$,
*		\f$\forall \mathbf{x} \not= \mathbf{0}\f$,
*               \f$\mathbf{x} \in R^n\f$.
*  		Convergence is tested using:
*		\f$ \frac{\| \mathbf{b} - mathbf{A} \mathbf{x} \|}
                         {\|\mathbf{b}\|} < \delta\f$.
*		If the preconditioning function pFn is non NULL then
*		it is called, passing the preconditioning data pDat, as:
*		(*pFn)(void *pDat, double **aM, double *r, double *z, int n)
*		to solve \f$\mathbf{A} \mathbf{z} = \mathbf{r}\f$ for
*		\f$\mathbf{z}\f$ with the solution overwriting the initial
*		contents of z.
* \param	aM			Matrix \f$\mathbf{A}\f$.
* \param	xV			Matrix \f$\mathbf{x}\f$ which
*					should contain an initial estimate
*					although this may be \f$\mathbf{0}\f$.
* \param	bV			Matrix \f$\mathbf{b}\f$.
* \param	wM			Matrix with dimensions [4, n],
*					this must be a rectangular matrix.
* \param	pFn			Preconditioning function.
* \param	pDat			Data to be passed to preconditioning
* 					function.
* \param	tol			Tolerance required, \f$\delta\f$.
* \param	maxItr			The maximum number of itterations.
* \param	dstTol			Destination pointer for the residual
*					after the final iteration, may be NULL.
* \param	dstItr			Destination pointer for the actual
*					number of itterations performed,
*					may be NULL.

*/
AlgError	AlgMatrixCGSolve(
			         AlgMatrix aM,
				 double *xV, double *bV,
				 AlgMatrix wM,
		     		 void (*pFn)(void *, AlgMatrix,
				 	     double *, double *),
				 void *pDat, double tol, int maxItr,
				 double *dstTol, int *dstItr)
{
  int		itr = 0,
  		conv = 0;
  size_t	nN;
  double	alpha,
		beta,
  		resid = DBL_MAX,
  		nrmB;
  double	*r,
  		*p,
		*z,
		*q;
  double	rho[2];
  AlgError	errCode = ALG_ERR_NONE;

  rho[0] = rho[1] = 0.0;
  if((aM.core == NULL) || (aM.core->nR < 1) || (aM.core->nR != aM.core->nC) ||
     (wM.core == NULL) || (wM.core->type != ALG_MATRIX_RECT) ||
     (wM.core->nR != 4) || (wM.core->nC != aM.core->nC) ||
     (xV == NULL) || (bV == NULL) || (tol < 0.0) || (maxItr < 0))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    switch(aM.core->type)
    {
      case ALG_MATRIX_LLR:  /* FALLTHROUGH */
      case ALG_MATRIX_RECT: /* FALLTHROUGH */
      case ALG_MATRIX_SYM:
	break;
      default:
        errCode = ALG_ERR_FUNC;
	break;
    }
  }
  if(errCode == ALG_ERR_NONE)
  {
    nN = aM.core->nR;
    AlgMatrixZero(wM);
    r = wM.rect->array[0];
    p = wM.rect->array[1];
    z = wM.rect->array[2];
    q = wM.rect->array[3];
    nrmB = AlgVectorNorm(bV, aM.core->nR);
    /* r = b - A x */
    AlgMatrixVectorMul(r, aM, xV);
    AlgVectorSub(r, bV, r, aM.core->nR);
    if(nrmB < DBL_EPSILON) 
    {
      nrmB = 1.0;
    }
    if((resid = AlgVectorNorm(r, nN) / nrmB) <= tol)
    {
      conv = 1;
    }
    while(!conv && (itr < maxItr))
    {
#ifdef ALG_MATRIXCG_DEBUG
      (void )fprintf(stderr, "rho = {%g, %g}\n", rho[0], rho[1]);
      AlgMatrixCGDebug(stderr, "x", xV, nN);
      AlgMatrixCGDebug(stderr, "r", r, nN);
#endif /* ALG_MATRIXCG_DEBUG */
      /* Pre-conditioning: Solve aM z = r for z. */
      if(pFn)
      {
	(*pFn)(pDat, aM, r, z);
      }
      else
      {
	AlgVectorCopy(z, r, nN);
      }
      rho[0] = AlgVectorDot(r, z, nN);
      if(itr == 0)
      {
	AlgVectorCopy(p, z, nN);
      }
      else
      {
	beta = rho[0] / rho[1];
	/* p = z + beta p */
	AlgVectorScaleAdd(p, p, z, beta, nN);
      }
      /* q = A p */
      AlgMatrixVectorMul(q, aM, p);
      alpha = rho[0] / AlgVectorDot(p, q, nN);
      /* x = x + alpha p */
      AlgVectorScaleAdd(xV, p, xV, alpha, nN);
      /* r = r - alpha * q */
      AlgVectorScaleAdd(r, q, r, -alpha, nN);
      if((resid = AlgVectorNorm(r, nN) / nrmB) <= tol)
      {
	tol = resid;
	conv = 1;
      }
      else
      {
	rho[1] = rho[0];
	++itr;
      }
    }
    if(!conv)
    {
      errCode = ALG_ERR_CONVERGENCE;
    }
  }
  if(dstTol)
  {
    *dstTol = resid;
  }
  if(dstItr)
  {
    *dstItr = itr;
  }
  return(errCode);
}

#ifdef ALG_MATRIXCG_DEBUG
static void	AlgMatrixCGDebug(FILE *fP, const char *name,
				 double *dat, int nDat)
{
  const char	*str0 = ", ",
  		*str1 = "}\n";

  (void )fprintf(stderr, "%s = {", name);
  while(nDat-- > 0)
  {
    (void )fprintf(stderr, "%g%s", *dat++, (nDat > 0)? str0: str1);

  }
}
#endif /* ALG_MATRIXCG_DEBUG */
