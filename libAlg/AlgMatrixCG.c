#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgMatrixCG.c
* \author       Bill Hill
* \date         March 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
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
* \todo         -
* \bug          None known.
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
* \param	aType			The type of matrix \f$\mathbf{A}\f$.
* \param	aM			Matrix \f$\mathbf{A}\f$.
* \param	xV			Matrix \f$\mathbf{x}\f$ which
*					should contain an initial estimate
*					although this may be \f$\mathbf{0}\f$.
* \param	bV			Matrix \f$\mathbf{b}\f$.
* \param	wM			Matrix with dimensions [4, n],
*					this must be a rectangular matrix.
* \param	n			The dimension of the matrix and
*					vectors \f$n\f$.
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
AlgError	AlgMatrixCGSolve(AlgMatrixType aType, double **aM,
				 double *xV, double *bV,
				 double **wM, size_t n,
		     		 void (*pFn)(void *, double **,
				 	     double *, double *, size_t),
				 void *pDat, double tol, int maxItr,
				 double *dstTol, int *dstItr)
{
  int		itr = 0,
  		conv = 0;
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
  if((aM == NULL) || (*aM == NULL) || (xV == NULL) || (bV == NULL) || 
     (wM == NULL) || (*wM == NULL) || (n < 1) || (tol < 0.0) || (maxItr < 0))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    switch(aType)
    {
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
    r = *(wM + 0);
    p = *(wM + 1);
    z = *(wM + 2);
    q = *(wM + 3);
    AlgMatrixZero(wM, 4, n);
    nrmB = AlgVectorNorm(bV, n);
    /* r = b - A x */
    AlgMatrixVectorMul(r, aType, aM, xV, n, n);
    AlgVectorSub(r, bV, r, n);
    if(nrmB < DBL_EPSILON) 
    {
      nrmB = 1.0;
    }
    if((resid = AlgVectorNorm(r, n) / nrmB) <= tol)
    {
      conv = 1;
    }
    while(!conv && (itr < maxItr))
    {
#ifdef ALG_MATRIXCG_DEBUG
      (void )fprintf(stderr, "rho = {%g, %g}\n", rho[0], rho[1]);
      AlgMatrixCGDebug(stderr, "x", xV, n);
      AlgMatrixCGDebug(stderr, "r", r, n);
#endif /* ALG_MATRIXCG_DEBUG */
      /* Pre-conditioning: Solve aM z = r for z. */
      if(pFn)
      {
	(*pFn)(pDat, aM, r, z, n);
      }
      else
      {
	AlgVectorCopy(z, r, n);
      }
      rho[0] = AlgVectorDot(r, z, n);
      if(itr == 0)
      {
	AlgVectorCopy(p, z, n);
      }
      else
      {
	beta = rho[0] / rho[1];
	/* p = z + beta p */
	AlgVectorScaleAdd(p, p, z, beta, n);
      }
      /* q = A p */
      AlgMatrixVectorMul(q, aType, aM, p, n, n);
      alpha = rho[0] / AlgVectorDot(p, q, n);
      /* x = x + alpha p */
      AlgVectorScaleAdd(xV, p, xV, alpha, n);
      /* r = r - alpha * q */
      AlgVectorScaleAdd(r, q, r, -alpha, n);
      if((resid = AlgVectorNorm(r, n) / nrmB) <= tol)
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

#if (ALG_MATRIXCG_TEST == 1)
int		main(int argc, char *argv[])
{
  int		id0,
  		itr = 1000;
  double	tol = 0.000001;
  double	**a,
  		**w;
  double	*b,
  		*x;
  AlgError	errCode = ALG_ERR_NONE;

  (void )AlcDouble2Malloc(&a, 3, 3);
  (void )AlcDouble2Malloc(&w, 4, 3);
  b = (double *)AlcMalloc(3 * sizeof(double));
  x = (double *)AlcMalloc(3 * sizeof(double));
  a[0][0] = 6.0; 
  a[0][1] = 1.0; 
  a[0][2] = 0.0; 
  a[1][0] = 1.0; 
  a[1][1] = 5.0; 
  a[1][2] = 1.0; 
  a[2][0] = 0.0; 
  a[2][1] = 1.0; 
  a[2][2] = 8.0; 
  b[0] = 8.0;
  b[1] = 14.0;
  b[2] = 26.0;
  x[0] = 1.0;
  x[1] = 1.0;
  x[2] = 1.0;
  errCode = AlgMatrixCGSolve(ALG_MATRIX_RECT,
  			     (double **)a, (double *)x, (double *)b,
                             (double **)w, 3, NULL, NULL, tol, itr,
			     &tol, &itr);
  (void )printf("itr = %d\n", itr);
  (void )printf("tol = %g\n", tol);
  for(id0 = 0; id0 < 3; ++id0)
  {
    printf("%g\n", x[id0]);
  }
}
#endif /* ALG_MATRIXCG_TEST == 1 */

#if (ALG_MATRIXCG_TEST == 2)
int		main(int argc, char *argv[])
{
  int		id0,
  		id1,
		rpt,
		sz = 100,
  		nRpt = 1,
  		itr = 1000;
  double 	tD0,
  		del;
  double	**a0,
  		**a1,
  		**w;
  double	*b0,
  		*b1,
  		*x0,
		*x1,
		*i0;
  AlgError	errCode = ALG_ERR_NONE;
  const double	tol = 0.000001;

  /* AlgRandSeed(0); */
  (void )AlcDouble2Malloc(&a0, sz, sz);
  (void )AlcDouble2Malloc(&a1, sz, sz);
  (void )AlcDouble2Malloc(&w, 4, sz);
  b0 = (double *)AlcMalloc(sz * sizeof(double));
  b1 = (double *)AlcMalloc(sz * sizeof(double));
  x0 = (double *)AlcMalloc(sz * sizeof(double));
  x1 = (double *)AlcMalloc(sz * sizeof(double));
  i0 = (double *)AlcMalloc(sz * sizeof(double));
  for(id1 = 0; id1 < sz; ++id1)
  {
    for(id0 = 0; id0 <= id1; ++id0)
    {
      *(*(a0 + id1) + id0) = *(*(a0 + id0) + id1) = AlgRandUniform();
    }
    *(*(a0 + id1) + id1) += 1.0;
    *(b0 + id1) = AlgRandUniform();
    *(i0 + id1) = 1.0;
  }
  for(rpt = 0; rpt < nRpt; ++rpt)
  {
    AlgVectorCopy(x0, i0, sz);
    AlgVectorCopy(x1, b0, sz);
    AlgMatrixCopy(a1, a0, sz, sz);
    del = tol;
    errCode = AlgMatrixCGSolve(ALG_MATRIX_RECT,
    			       (double **)a0, (double *)x0, (double *)b0,
			       (double **)w, sz, NULL, NULL, tol, itr,
			       &del, &itr);
    errCode = AlgMatrixSVSolve((double **)a1, sz, sz, x1, tol);
    for(id1 = 0; id1 < sz; ++id1)
    {
      tD0 = *(x0 + id1) - *(x1 + id1);
      del += tD0 * tD0;
    }
    del = (del > DBL_EPSILON)? sqrt(del) / sz: 0.0;
  }
  printf("rms %g\n", del);
}
#endif /* ALG_MATRIXCG_TEST == 2 */

#if (ALG_MATRIXCG_TEST == 3)
extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		id0,
  		id1,
		rpt,
		option,
		ok = 1,
		useCG = 1,
		sz = 100,
  		nRpt = 1,
  		itr = 1000;
  AlgMatrixType	aType = ALG_MATRIX_RECT;
  double 	tD0,
  		del;
  double	**a0,
  		**a1,
  		**w;
  double	*b0,
  		*b1,
  		*x0,
		*i0;
  AlgError	errCode = ALG_ERR_NONE;
  const double	tol = 0.000001;
  const char	*optList = "CSRYr:s:";

  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'C':
        useCG = 1;
	break;
      case 'S':
        useCG = 0;
	break;
      case 'R':
        aType = ALG_MATRIX_RECT;
	break;
      case 'Y':
        aType = ALG_MATRIX_SYM;
	break;
      case 'r':
        if((sscanf(optarg, "%d", &nRpt) != 1) || (nRpt < 0))
	{
	  ok = 0;
	}
	break;
      case 's':
        if((sscanf(optarg, "%d", &sz) != 1) || (sz < 1))
	{
	  ok = 0;
	}
	break;
      default:
        ok = 0;
	break;
    }
  }
  if(ok)
  {
    if(!useCG && (aType != ALG_MATRIX_RECT))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Non-rectangular matrices can be used with the\n"
      		      "conjugate gradient method\n");
    }
  }
  if(ok)
  {
    switch(aType)
    {
      case ALG_MATRIX_RECT:
        (void )AlcDouble2Malloc(&a0, sz, sz);
	break;
      case ALG_MATRIX_SYM:
        (void )AlcSymDouble2Malloc(&a0, sz);
	break;
    }
    if(!useCG)
    {
      (void )AlcDouble2Malloc(&a1, sz, sz);
    }
    (void )AlcDouble2Malloc(&w, 4, sz);
    b0 = (double *)AlcMalloc(sz * sizeof(double));
    b1 = (double *)AlcMalloc(sz * sizeof(double));
    x0 = (double *)AlcMalloc(sz * sizeof(double));
    i0 = (double *)AlcMalloc(sz * sizeof(double));
    for(id1 = 0; id1 < sz; ++id1)
    {
      switch(aType)
      {
        case ALG_MATRIX_RECT:
	  for(id0 = 0; id0 <= id1; ++id0)
	  {
	    a0[id1][id0] = a0[id0][id1] = ((id0 + id1) % 7) * 0.01;
	  }
	  a0[id1][id1] += 1.0;
	  break;
        case ALG_MATRIX_SYM:
	  for(id0 = 0; id0 <= id1; ++id0)
	  {
	    a0[id1][id0] = ((id0 + id1) % 7) * 0.01;
	  }
	  a0[id1][id1] += 1.0;
	  break;
      }
      *(b0 + id1) = ((id1) % 5) * 0.1;
      *(i0 + id1) = 1.0;
    }
#ifdef ALG_MATRIXCG_DEBUG
    for(id1 = 0; id1 < sz; ++id1)
    {
      for(id0 = 0; id0 <= id1; ++id0)
      {
        (void )fprintf(stderr, "%g ", a0[id1][id0]);
      }
      (void )fprintf(stderr, "\n");
    }
#endif /* ALG_MATRIXCG_DEBUG */
    for(rpt = 0; rpt < nRpt; ++rpt)
    {
      if(useCG)
      {
	del = tol;
	AlgVectorCopy(x0, i0, sz);
	errCode = AlgMatrixCGSolve(aType, (double **)a0,
				   (double *)x0, (double *)b0,
				   (double **)w, sz, NULL, NULL,
				   tol, itr, &del, &itr);
      }
      else
      {
	AlgVectorCopy(x0, b0, sz);
	AlgMatrixCopy(a1, a0, sz, sz);
	errCode = AlgMatrixSVSolve((double **)a1, sz, sz, x0, tol);
      }
    }
    for(id1 = 0; id1 < sz; ++id1)
    {
      (void )printf("%g\n", *(x0 + id1));
    }
  }
  else
  {
    (void )fprintf(stderr, "Usage: %s [-C] [-S] [-R] [-Y] [-s #] [-r #]\n%s",
		   *argv,
		   "Test for timing solution of Ax = b using the Conjugate\n"
		   "Gradient and Singular Value Decomposition algorithms.\n"
                   "  -C  Use the Conjugate Gradient algorithm.\n"
                   "  -S  Use the Singular Value Decomposition algorithm.\n"
		   "  -R  Use a rectangular matrix.\n"
		   "  -Y  Use a symetric matrix.\n"
                   "  -s  matrix size.\n"
                   "  -r  Number of repeat solutions.\n");
  }
}
#endif /* ALG_MATRIXCG_TEST == 3 */
