#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixLSQR_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixLSQR.c
* \author       Bill Hill
* \date         May 2010
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
* \brief	Provides functions for solving matrix equations using LSQR.
* 		This software is based on lsqr.c, a C version of LSQR derived
* 		by James W. Howse <jhowse@lanl.gov> from the Fortran 77
* 		implementation of C. C. Paige and M. A. Saunders.
*		In most cases the extensive comments from Howse's lsqr.c
*		have been preserved with little change.
* \ingroup	AlgMatrix
*/

#include <Alg.h>
#include <float.h>


static double 			AlgMatrixLSQRNorm2(
				  double a,
				  double b);
static double 			AlgMatrixLSQRNorm4(
				  double a,
				  double b,
				  double c,
				  double d);
/*!
* \return
* \ingroup	AlgMatrix
* \brief	This function finds a solution x to the following problems:
*
*     1. Unsymmetric equations --    solve  A*x = b
*
*     2. Linear least squares  --    solve  A*x = b
*                                    in the least-squares sense
*
*     3. Damped least squares  --    solve  (   A    )*x = ( b )
*                                           ( damp*I )     ( 0 )
*                                    in the least-squares sense
*
*     where 'A' is a matrix with 'm' rows and 'n' columns, 'b' is an
*     'm'-vector, and 'damp' is a scalar.  (All quantities are real.)
*     The matrix 'A' is intended to be large and sparse.  
*
*     Notation
*     --------
*     The following quantities are used in discussing the subroutine
*     parameters:
*
*     'Abar'   =  (   A    ),          'bbar'  =  ( b )
*                 ( damp*I )                      ( 0 )
*
*     'r'      =  b  -  A*x,           'rbar'  =  bbar  -  Abar*x
*
*     'rnorm'  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
*              =  norm( rbar )
*
*     'rel_prec'  =  the relative precision of floating-point arithmetic
*                    on the machine being used.  Typically 2.22e-16
*                    with 64-bit arithmetic.
*
*     LSQR  minimizes the function 'rnorm' with respect to 'x'.
*
*     References
*     ----------
*     C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
*          linear equations and sparse least squares,
*          ACM Transactions on Mathematical Software 8, 1 (March 1982),
*          pp. 43-71.
*     C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
*          linear equations and least-squares problems,
*          ACM Transactions on Mathematical Software 8, 2 (June 1982),
*          pp. 195-209.
*     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
*          Basic linear algebra subprograms for Fortran usage,
*          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
*          pp. 308-323 and 324-325.
*
* \param	aType
* \param	aM			Matrix A with nR rows and nC columns.
* 					The values of A are not modified by
* 					this function.
* \param	nR			Number of rows in matrix A.
* \param	nC             		Number of columns in matrix A.
* \param	bV			Vector b with nR entries which are
* 					modified by this function.
* \param	xV			Vector x for with initial guess at
* 					solution and for the return of the
* 					solution with nC entries.
* \param	damping			Damping parameter, set to 0.0 for no
* 					damping.
* \param	relErrA			An estimate of the relative error in
* 					matrix A. If A is accurate to about
* 					6 digits, set to 1.0e-06.
* \param	relErrB			An estimate of the relative error in
* 					vector b. Set as for relErrA.
* \param	maxItr			An upper limit on the number of
* 					iterations. If zero set to nC.
* \param	condLim			Upper limit on the condition number
* 					of matrix 'Abar'. Iterations will be
* 					terminated if a computed estimate of
* 					exceeds this condition number.
* \param	dstTerm			Destination pointer for termination
* 					code, which may be NULL. The values
* 					used are:
*   - 0       x = x0  is the exact solution. No iterations were performed.
*   - 1      The equations A*x = b are probably compatible. Norm(A*x - b)
*             is sufficiently small, given the values of relErrA and relErrB.
*   - 2       The system A*x = b is probably not compatible.  A least-squares
*             solution has been obtained that is sufficiently accurate, given
*             the value of relErrA.
*   - 3       An estimate of cond('Abar') has exceeded condLim.  The system
*             A*x = b appears to be ill-conditioned.
*   - 4       The equations A*x = b are probably compatible.  Norm(A*x - b)
*             is as small as seems reasonable on this machine.
*   - 5       The system A*x = b is probably not compatible.  A least-squares
*   	      solution has been obtained that is as accurate as seems
*   	      reasonable on this machine.
*   - 6       Cond('Abar') seems to be so large that there is no point in
*   	      doing further iterations, given the precision of this machine.
*   - 7       The iteration limit maxItr was reached.
* \param	dstItr			Destination pointer for the number of
* 				 	iterations performed, which may be
* 				 	NULL.
*
* \param	dstFNorm		Destination pointer for an estimate of
* 					the Frobenius norm of 'Abar', may be
* 					NULL.
* 					This is the square-root of the sum of
* 					squares of the elements of 'Abar'. If
* 					'damping' is small and if the columns
* 					of 'A' have all been scaled to have
* 					length 1.0, then the Frobenius norm
* 					should increase to roughly sqrt('nC').
* \param	dstCondN		Destination pointer for an estimate of
* 					the condition number of 'Abar', may be
* 					NULL.
* \param	dstResNorm		Destination pointer for an estimate of
* 					the final value of norm('rbar'), the
* 					function being minimized. This will be
* 					small if A*x = b has a solution. May be
* 					NULL.
* \param	dstResNormA		Destination pointer for an estimate of
* 					the final value of the residual for the
* 					usual normal equations. This should be
* 					small in all cases. May be NULL.
* \param	dstNormX		Destination pointer for an estimate of
* 					the final solution vector 'x'.
*/
AlgError 	AlgMatrixSolveLSQR(AlgMatrix aM,
				double *bV, double *xV,
				double damping, double relErrA, double relErrB,
				long maxItr, long condLim,
				int *dstTerm, long *dstItr, double *dstFNorm,
				double *dstCondN, double *dstResNorm,
				double *dstResNormA, double *dstNormX)
{
  long    	idx,
          	termItrMax,
  	  	itr = 0,
		term = 0,
          	termItr = 0;
  size_t	nR,
  		nC;
  double  	bnorm,
		condNTol,
  		cs,
		cs1,
		delta,
		gamma,
		gammabar,
		phi,
		psi,
		resNorm,
		resNormA,
		resTol,
		resTolMach,
		rho,
		rhobar1,
		sn,
		sn1,
		t0,
		t1,
		t2,
		t3,
		tau,
		temp,
		theta,
		zetabar,
  		alpha = 0.0, 
          	bbnorm = 0.0,
		beta = 0.0,
		condN = 0.0,
		cs2 = -1.0,
		ddnorm = 0.0,
		fnorm = 0.0,
		normX = 0.0,
		phibar = 0.0,
		res = 0.0,
		rhobar = 0.0,
		sn2 = 0.0,
		stopCrit1,
		stopCrit2,
		stopCrit3,
		xxnorm = 0.0,
		zeta = 0.0;
  double	*vV = NULL,
  		*wV = NULL;
  AlgError	errNum = ALG_ERR_NONE;

  nR = aM.core->nR;
  nC = aM.core->nC;
  if(maxItr <= 0)
  {
    maxItr = nC;
  }
  if(condLim > 0.0)
  {
    condNTol = 1.0 / condLim;
  }
  else
  {
    condNTol = DBL_EPSILON;
  }
#ifdef ALG_MATRIXLSQR_DEBUG
  (void )fprintf(stderr,
		 "AlgMatrixSolveLSQR()\n"
                 "  Least Squares Solution of A*x = b\n"
		 "  Matrix A is %7li rows x %7li columns\n"
		 "  Damping parameter = %10.2e\n"
		 "  ATOL = %10.2e\t\tCONDLIM = %10.2e\n"
		 "  BTOL = %10.2e\t\tITERLIM = %10li\n\n",
        	 nR, nC, damping, relErrA,
        	 condLim, relErrB, maxItr);
#endif
  /* Allocate workspace vectors. */
  if(((vV = (double *)AlcCalloc(nC, sizeof(double))) == NULL) ||
     ((wV = (double *)AlcCalloc(nC, sizeof(double))) == NULL))
  {
    errNum = ALG_ERR_MALLOC;
  }
  if(errNum == ALG_ERR_NONE)
  {
    /* Set up the initial vectors u and v for bidiagonalization.  These
     * satisfy  the relations
     * beta*u = b - A*x0 
     * alpha*v = A^T*u
     */
    /* Compute Euclidean length of u and store as beta */
    beta = AlgVectorNorm(bV, nR);
    if(beta > 0.0)
    {
      /* Scale vector u by the inverse of beta */
      AlgVectorScale(bV, bV, 1.0 / beta, nR);
      /* Compute matrix-vector product A^T*u and store it in vector v */
      AlgMatrixTVectorMulAdd(vV, aM, bV, vV);
      /* Compute Euclidean length of v and store as alpha */
      alpha = AlgVectorNorm(vV, nC);
    }
    if(alpha > 0.0)
    {
      /* Scale vector v by the inverse of alpha */
      AlgVectorScale(vV, vV, 1.0 / alpha, nC);
      /* Copy vector v to vector w */
      AlgVectorCopy(wV, vV, nC);
    }    
    resNormA = alpha * beta;
    resNorm = beta;
    bnorm = beta;
    /*
     *  If the norm || A^T r || is zero, then the initial guess is the exact
     *  solution.  Exit and report this.
     */
    if(resNormA == 0.0)
    {
      term = 0;
#ifdef ALG_MATRIXLSQR_DEBUG
      (void )fprintf(stderr,
		     "AlgMatrixSolveLSQR()\n"
		     "  ISTOP = %3li\t\t\tITER = %9li\n"
		     "  || A ||_F = %13.5e\tcond(A) = %13.5e\n"
		     "  || r ||_2 = %13.5e\t|| A^T r ||_2 = %13.5e\n"
		     "  || b ||_2 = %13.5e\t|| x - x0 ||_2 = %13.5e\n\n", 
		     term, itr, fnorm, 
		     condN, resNorm, resNormA,
		     bnorm, normX);
#endif
    }
    else
    {
      rhobar = alpha;
      phibar = beta;
#ifdef ALG_MATRIXLSQR_DEBUG
      stopCrit1 = 1.0;
      stopCrit2 = alpha / beta;
      (void )fprintf(stderr,
		     "AlgMatrixSolveLSQR()\n"
		     " Er ||  ITER     || r ||    Compatible  "
		     "||A^T r|| / ||A|| ||r||  || A ||    cond(A)\n"
		     "%6li %13.5e %10.2e \t%10.2e \t%10.2e  %10.2e\n",
		     itr, resNorm, stopCrit1, stopCrit2,
		     fnorm, condN);
#endif
    }
    /*
     *  The main iteration loop is continued as long as no stopping criteria
     *  are satisfied and the number of total iterations is less than some upper
     *  bound.
     */
    while(term == 0)
    {
      ++itr;
      /*      
       *     Perform the next step of the bidiagonalization to obtain
       *     the next vectors u and v, and the scalars alpha and beta.
       *     These satisfy the relations
       *                beta*u  =  A*v  -  alpha*u,
       *                alpha*v =  A^T*u  -  beta*v.
       */      
      /* Scale vector u by -alpha */
      AlgVectorScale(bV, bV, -alpha, nR);
      /* Compute A*v - alpha*u and store in vector u */
      AlgMatrixVectorMulAdd(bV, aM, vV, bV);
      /* Compute Euclidean length of u and store as beta */
      beta = AlgVectorNorm(bV, nR);
      /* Accumulate this quantity to estimate Frobenius norm of matrix A */
      bbnorm = AlgMatrixLSQRNorm4(alpha, beta, damping, bbnorm);
      if(beta > 0.0)
      {
	/* Scale vector u by 1.0 / beta */
	AlgVectorScale(bV, bV, 1.0 / beta, nR);
	/* Scale vector v by -beta */
	AlgVectorScale(vV, vV, -beta, nC);
	/* Compute A^T*u - beta*v and store in vector v */
	AlgMatrixTVectorMulAdd(vV, aM, bV, vV);
	/* Compute Euclidean length of v and store as alpha */
	alpha = AlgVectorNorm(vV, nC);
	if(alpha > 0.0)
	{
	  /* Scale vector v by the inverse of alpha */
	  AlgVectorScale(vV, vV, 1.0 / alpha, nC); 
	}
      }
      /* Use a plane rotation to eliminate the damping parameter.
       * This alters the diagonal (RHOBAR) of the lower-bidiagonal matrix.
       */
      rhobar1 = AlgMatrixLSQRNorm2(rhobar, damping);
      cs1     = rhobar / rhobar1;
      sn1     = damping / rhobar1;
      psi     = sn1 * phibar;
      phibar  = cs1 * phibar;
      /* Use a plane rotation to eliminate the subdiagonal element (beta)
       * of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
       */
      rho    =  AlgMatrixLSQRNorm2(rhobar1, beta);
      cs     =  rhobar1 / rho;
      sn     =  beta / rho;
      theta  =  sn * alpha;
      rhobar = -cs * alpha;
      phi    =  cs * phibar;
      phibar =  sn * phibar;
      tau    =  sn * phi;
      /* Update the solution vector x, the search direction vector w.
       */     
      t3 = 1.0 / rho;
      t1 = phi * t3;
      t2 = -theta * t3;
      for(idx = 0; idx < nC; idx++)
      {
	t0 = wV[idx];
	/* Update the solution vector x */
	xV[idx] = t0 * t1 + xV[idx];
	/* Update the search direction vector w */
	wV[idx] = t0 * t2 + vV[idx];
	/* Accumulate this quantity to estimate condition number of A */
	ddnorm += ALG_SQR(t0 * t3);
      }
      ddnorm = sqrt(ddnorm);
      /* Use a plane rotation on the right to eliminate the super-diagonal
       * element (THETA) of the upper-bidiagonal matrix.  Then use the result
       * to estimate the solution norm || x ||.
       */
      delta    = sn2 * rho;
      gammabar = -cs2 * rho;
      zetabar = (phi - delta * zeta) / gammabar;
      /* Compute an estimate of the solution norm || x || */
      normX = sqrt(xxnorm + ALG_SQR(zetabar));
      gamma = AlgMatrixLSQRNorm2(gammabar, theta);
      cs2 = gammabar / gamma;
      sn2 = theta / gamma;
      zeta = (phi - delta * zeta) / gamma;
      /* Accumulate this quantity to estimate solution norm || x || */
      xxnorm = AlgMatrixLSQRNorm2(xxnorm, zeta);
      /* Estimate the norm and condition of the matrix A, and the norms of
       * the vectors r and A^T*r.
       */
      condN = fnorm * sqrt(ddnorm);
      res = AlgMatrixLSQRNorm2(res, psi);
      resNorm = AlgMatrixLSQRNorm2(phibar, res) + 1e-30;     /* avoid == 0.0 */
      resNormA = alpha * fabs(tau);

      /* Use these norms to estimate certain quantities, some of which will be
       * small near a solution.
       */
      stopCrit1 = resNorm / bnorm;
      stopCrit2 = 0.0;
      if(resNorm > 0.0)
      {
	stopCrit2 = resNormA / (bbnorm * resNorm);
      }
      stopCrit3 = 1.0 / condN;
      resTol = relErrB + relErrA * fnorm * normX / bnorm;
      resTolMach = DBL_EPSILON + DBL_EPSILON * fnorm * normX / bnorm;
      /* Check to see if any of the stopping criteria are satisfied.
       * First compare the computed criteria to the machine precision.
       * Second compare the computed criteria to the the user specified
       * precision.
       */
      if(itr >= maxItr)
      {
        /* Iteration limit reached */
	term = 7;
	errNum = ALG_ERR_CONVERGENCE;
      }
      else if(stopCrit3 <= DBL_EPSILON)
      {
        /* Condition number greater than machine precision */
	term = 6;
	errNum = ALG_ERR_MATRIX_CONDITION;
      }
      else if(stopCrit2 <= DBL_EPSILON)
      {
        /* Least squares error less than machine precision */
	term = 5;
      }
      else if(stopCrit1 <= resTolMach)
      {
        /* Residual less than a function of machine precision */
	term = 4;
      }
      else if(stopCrit3 <= condNTol)
      {
        /* Condition number greater than CONLIM */
	term = 3;
	errNum = ALG_ERR_MATRIX_CONDITION;
      }
      else if(stopCrit2 <= relErrA)
      {
        /* Least squares error less than ATOL */
	term = 2;
      }
      else if(stopCrit1 <= resTol)
      {
        /* Residual less than a function of ATOL and BTOL */
	term = 1;
      }
#ifdef ALG_MATRIXLSQR_DEBUG
      (void )fprintf(stderr,
		     "AlgMatrixSolveLSQR()\n"
		     "%2d %6li %10.5e %10.2e \t%10.2e \t%10.2e %10.2e\n",
		     errNum,
		     itr, resNorm, stopCrit1, 
		     stopCrit2,
		     fnorm, condN);
#endif
      /* The convergence criteria are required to be met on NCONV consecutive 
       * iterations, where NCONV is set below.  Suggested values are 1, 2, or
       * 3.
       */
      if(term == 0)
      {
	termItr = -1;
      }
      termItrMax = 1;
      ++termItr;
      if((termItr < termItrMax) && (itr < maxItr))
      {
	term = 0;
      }
    }
    /* Finish computing the standard error estimates vector se.
     */
    temp = 1.0;
    if(nR > nC)
    {
      temp = (double )(nR - nC);
    }
    if(ALG_SQR(damping) > 0.0)
    {
      temp = (double )nR;
    }
    temp = resNorm / sqrt(temp);
    if(dstTerm)
    {
      *dstTerm = term;
    }
    if(dstItr)
    {
      *dstItr = itr;
    }
    if(dstFNorm)
    {
      *dstFNorm = fnorm;
    }
    if(dstCondN)
    {
      *dstCondN = condN;
    }
    if(dstResNorm)
    {
      *dstResNorm = resNorm;
    }
    if(dstResNormA)
    {
      *dstResNormA = resNormA;
    }
    if(dstNormX)
    {
      *dstNormX = normX;
    }
  }
  AlcFree(vV);
  AlcFree(wV);
#ifdef ALG_MATRIXLSQR_DEBUG
  (void )fprintf(stderr,
		 "AlgMatrixSolveLSQR()\n"
		 "errcode = %3d\n"
		 "\tISTOP = %3li\t\t\tITER = %9li\n"
		 "  || A ||_F = %13.5e\tcond( A ) = %13.5e\n"
		 "  || r ||_2 = %13.5e\t|| A^T r ||_2 = %13.5e\n"
		 "  || b ||_2 = %13.5e\t|| x - x0 ||_2 = %13.5e\n\n", 
		 errNum,
		 term, itr, fnorm, 
		 condN, resNorm, resNormA,
		 bnorm, normX );
#endif
  return(errNum);
}

/*!
* \return	Norm of the given values.
* \ingroup	AlgMatrix
* \brief	Computes the norm of the given values, i.e.
*               \f$\sqrt{a^2 + b^2}\f$ but with precautions to avoid overflow.
* \param	a			First value.
* \param	b			Second value.
*/
static double 	AlgMatrixLSQRNorm2(double a, double b)
{
  double	scale,
  		norm = 0.0;

  scale = fabs(a) + fabs(b);
  if(scale > DBL_EPSILON)
  {
    a /= scale;
    b /= scale;
    norm = scale * sqrt(ALG_SQR(a) + ALG_SQR(b));
  }
  return(norm);
}

/*!
* \return	Norm of the given values.
* \ingroup	AlgMatrix
* \brief	Computes the norm of the given values, i.e.
*               \f$\sqrt{a^2 + b^2 + c^2 + d^2}\f$ but with precautions to avoid
*               overflow.
* \param	a			First value.
* \param	b			Second value.
* \param	c			Third value.
*/
static double 	AlgMatrixLSQRNorm4(double a, double b, double c, double d)
{
  double	scale,
  		norm = 0.0;

  scale = fabs(a) + fabs(b) + fabs(c) + fabs(d);
  if(scale > DBL_EPSILON)
  {
    a /= scale;
    b /= scale;
    c /= scale;
    d /= scale;
    norm = scale * sqrt(ALG_SQR(a) + ALG_SQR(b) + ALG_SQR(c) + ALG_SQR(d));
  }
  return(norm);
}
