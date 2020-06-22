#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgBSpline_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         AlgBSpline.c
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
* Copyright (C), [2020],
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
* \brief	Provides functions for fitting and evaluating B-Splines.
* 		This software is based on netlib Dierckx Fortran subroutines.
* 		In most cases the original comments have been preserved with
* 		little change and are inculded in the Doxygen documentation
* 		verbatim.
* \ingroup	AlgFit
*/

#include <Alg.h>
#include <float.h>

void				AlgBSplineBspl(
				  double *t,
				  int k,
				  double x,
				  int l, double *h);

static void			AlgBSplineBack(
				  double *a,
				  double *z,
				  int n,
				  int k,
				  double *c,
				  int nest);
static void 			AlgBSplineBacp(
				  double *a,
				  double *b,
				  double *z,
				  int n,
				  int k,
				  double *c,
				  int nest);
static void			AlgBSplineRota(
				  double cs,
				  double sn,
				  double *a,
				  double *b);
static void			AlgBSplineKnot(
				  double *x,
				  int m,
				  double *t,
				  int *n,
				  double * fpint,
				  int *nrdata,
				  int *nrint,
				  int istart);
static void			AlgBSplineGivens(
				  double piv,
				  double *ww,
				  double *c,
				  double *s);
static void			AlgBSplineDisc(
				  double *t,
				  int n,
				  int k2,
				  double *b,
				  int nest);
static int 			AlgBSplineChec(
				  double *x,
				  int m,
				  double *t,
				  int n,
				  int k);
static int			AlgBSplineChep(
				  double *x,
				  int m,
				  double *t,
				  int n,
				  int k);
static int			AlgBSplinePeri(
				  int iopt,
				  double *x,
				  double *y,
				  double *w,
				  int m,
				  int k,
				  double s,
				  int nest,
				  double tol,
				  int maxit,
				  int k1,
				  int k2,
				  int *n,
				  double *t,
				  double *c,
				  double *fp,
				  double *fpint,
				  double *z,
				  double *a1,
				  double *a2,
				  double *b,
				  double *g1,
				  double *g2,
				  double *q,
				  int *nrdata);
static int 			AlgBSplineCurf(
				  int iopt,
				  double *x,
				  double *y,
				  double *w,
				  int m,
				  double xb,
				  double xe,
				  int k,
				  double s,
				  int nest,
				  double tol,
				  int maxit,
				  int k1,
				  int k2,
				  int *n,
				  double *t,
				  double *c,
				  double *fp,
				  double *fpint,
				  double *z,
				  double *a,
				  double *b,
				  double *g,
				  double *q,
				  int *nrdata);
static int 			AlgBSplinePara(
				  int iopt,
				  int idim,
				  int m,
				  double *u,
				  int *mx,
				  double *x,
				  double *w,
				  double *ub,
				  double *ue,
				  int k,
				  double s,
				  int nest,
				  double tol,
				  int maxit,
				  int k1,
				  int k2,
				  int *n,
				  double *t,
				  int *nc,
				  double *c,
				  double *fp,
				  double *fpint,
				  double *z,
				  double *a,
				  double *b,
				  double *g,
				  double *q,
				  int *nrdata);
static double			AlgBSplineRati(
				  double *p1,
				  double *f1,
				  double p2,
				  double f2,
				  double *p3,
				  double *f3);
static AlgError			AlgErrorFromDierckx(
				  int ier);

/*!
* \ingroup	AlgFit
* \brief	Evaluates the (k+1) non-zero B-Splines of degree k at t[l]
* 		using the de Boor Cox recurrence.
* 
* This function has been derived from the netlib Dierckx function fpbspl().
* The original fortran coments are:
*
*  \verbatim
   Subroutine fpbspl evaluates the (k+1) non-zero b-splines of
   degree k at t(l) <= x < t(l+1) using the stable recurrence
   relation of de Boor and Cox.
   \endverbatim
*
* \param	t		Array of knot positions.
* \param	k		Degree of the B-Spline.
* \param	x		Point at which the spline is to be evaluated.
* \param	l		Knot index st k at t[l] <= x < t[l+1].
* \param	h		Array with length >= k in which to evaluate
* 				the spline.
*/
void				AlgBSplineBspl(
				  double *t,
				  int k,
				  double x,
				  int l,
				  double *h)
{
  int		j;
  double	hh[5];

  h[0] = 1.0;
  for(j = 1; j <= k; ++j)
  {
    int		i;

    for(i = 0; i < j; ++i)
    {
      hh[i] = h[i];
    }
    h[0] = 0.0;
    for(i = 0; i < j; ++i)
    {
      int	li,
      		lj;
      double	f;

      li = l + i;
      lj = li - j;
      f = hh[i] / (t[li] - t[lj]);
      h[i] += f * (t[li] - x);
      h[i + 1] = f * (x - t[lj]);
    }
  }
}
/*!
* \return	Alg error code.
* \ingroup	AlgFit
* \brief	Computes a smooth approximating b-spline in multi-dimensional
*		space.
*
* This function has been derived from the netlib Dierckx function curfit().
* The original fortran coments are:
*
*  \verbatim
   Given the ordered set of m points x(i) in the idim-dimensional space
   and given also a corresponding set of strictly increasing values u(i)
   and the set of positive numbers w(i),i=1,2,...,m, subroutine parcur
   determines a smooth approximating spline curve s(u), i.e.
       x1 = s1(u)
       x2 = s2(u)       ub <= u <= ue
       .........
       xidim = sidim(u)
   with sj(u),j=1,2,...,idim spline functions of degree k with common
   knots t(j),j=1,2,...,n.
   if ipar=1 the values ub,ue and u(i),i=1,2,...,m must be supplied by
   the user. if ipar=0 these values are chosen automatically by parcur
   as  v(1) = 0
       v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
       u(i) = v(i)/v(m) ,i=1,2,...,m
       ub = u(1) = 0, ue = u(m) = 1.
   if iopt=-1 parcur calculates the weighted least-squares spline curve
   according to a given set of knots.
   if iopt>=0 the number of knots of the splines sj(u) and the position
   t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
   ness of s(u) is then achieved by minimalizing the discontinuity
   jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
   n-k-1. the amount of smoothness is determined by the condition that
   f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
   negative constant, called the smoothing factor.
   the fit s(u) is given in the b-spline representation and can be
   evaluated by means of subroutine curev.
 
   calling sequence:
      call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,
     * fp,wrk,lwrk,iwrk,ier)
 
   parameters:
    iopt  : integer flag. on entry iopt must specify whether a weighted
            least-squares spline curve (iopt=-1) or a smoothing spline
            curve (iopt=0 or 1) must be determined.if iopt=0 the routine
            will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
            i=1,2,...,k+1. if iopt=1 the routine will continue with the
            knots found at the last call of the routine.
            attention: a call with iopt=1 must always be immediately
            preceded by another call with iopt=1 or iopt=0.
            unchanged on exit.
    ipar  : integer flag. on entry ipar must specify whether (ipar=1)
            the user will supply the parameter values u(i),ub and ue
            or whether (ipar=0) these values are to be calculated by
            parcur. unchanged on exit.
    idim  : integer. on entry idim must specify the dimension of the
            curve. 0 < idim < 11.
            unchanged on exit.
    m     : integer. on entry m must specify the number of data points.
            m > k. unchanged on exit.
    u     : real array of dimension at least (m). in case ipar=1,before
            entry, u(i) must be set to the i-th value of the parameter
            variable u for i=1,2,...,m. these values must then be
            supplied in strictly ascending order and will be unchanged
            on exit. in case ipar=0, on exit,array u will contain the
            values u(i) as determined by parcur.
    mx    : integer. on entry mx must specify the actual dimension of
            the array x as declared in the calling (sub)program. mx must
            not be too small (see x). unchanged on exit.
    x     : real array of dimension at least idim*m.
            before entry, x(idim*(i-1)+j) must contain the j-th coord-
            inate of the i-th data point for i=1,2,...,m and j=1,2,...,
            idim. unchanged on exit.
    w     : real array of dimension at least (m). before entry, w(i)
            must be set to the i-th value in the set of weights. the
            w(i) must be strictly positive. unchanged on exit.
            see also further comments.
    ub,ue : real values. on entry (in case ipar=1) ub and ue must
            contain the lower and upper bound for the parameter u.
            ub <=u(1), ue>= u(m). if ipar = 0 these values will
            automatically be set to 0 and 1 by parcur.
    k     : integer. on entry k must specify the degree of the splines.
            1<=k<=5. it is recommended to use cubic splines (k=3).
            the user is strongly dissuaded from choosing k even,together
            with a small s-value. unchanged on exit.
    s     : real.on entry (in case iopt>=0) s must specify the smoothing
            factor. s >=0. unchanged on exit.
            for advice on the choice of s see further comments.
    nest  : integer. on entry nest must contain an over-estimate of the
            total number of knots of the splines returned, to indicate
            the storage space available to the routine. nest >=2*k+2.
            in most practical situation nest=m/2 will be sufficient.
            always large enough is nest=m+k+1, the number of knots
            needed for interpolation (s=0). unchanged on exit.
    n     : integer.
            unless ier = 10 (in case iopt >=0), n will contain the
            total number of knots of the smoothing spline curve returned
            if the computation mode iopt=1 is used this value of n
            should be left unchanged between subsequent calls.
            in case iopt=-1, the value of n must be specified on entry.
    t     : real array of dimension at least (nest).
            on succesful exit, this array will contain the knots of the
            spline curve,i.e. the position of the interior knots t(k+2),
            t(k+3),..,t(n-k-1) as well as the position of the additional
            t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
            the b-spline representation.
            if the computation mode iopt=1 is used, the values of t(1),
            t(2),...,t(n) should be left unchanged between subsequent
            calls. if the computation mode iopt=-1 is used, the values
            t(k+2),...,t(n-k-1) must be supplied by the user, before
            entry. see also the restrictions (ier=10).
    nc    : integer. on entry nc must specify the actual dimension of
            the array c as declared in the calling (sub)program. nc
            must not be too small (see c). unchanged on exit.
    c     : real array of dimension at least (nest*idim).
            on succesful exit, this array will contain the coefficients
            in the b-spline representation of the spline curve s(u),i.e.
            the b-spline coefficients of the spline sj(u) will be given
            in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
    fp    : real. unless ier = 10, fp contains the weighted sum of
            squared residuals of the spline curve returned.
    wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
            used as working space. if the computation mode iopt=1 is
            used, the values wrk(1),...,wrk(n) should be left unchanged
            between subsequent calls.
    lwrk  : integer. on entry,lwrk must specify the actual dimension of
            the array wrk as declared in the calling (sub)program. lwrk
            must not be too small (see wrk). unchanged on exit.
    iwrk  : integer array of dimension at least (nest).
            used as working space. if the computation mode iopt=1 is
            used,the values iwrk(1),...,iwrk(n) should be left unchanged
            between subsequent calls.
    ier   : integer. unless the routine detects an error, ier contains a
            non-positive value on exit, i.e.
     ier=0  : normal return. the curve returned has a residual sum of
              squares fp such that abs(fp-s)/s <= tol with tol a relat-
              ive tolerance set to 0.001 by the program.
     ier=-1 : normal return. the curve returned is an interpolating
              spline curve (fp=0).
     ier=-2 : normal return. the curve returned is the weighted least-
              squares polynomial curve of degree k.in this extreme case
              fp gives the upper bound fp0 for the smoothing factor s.
     ier=1  : error. the required storage space exceeds the available
              storage space, as specified by the parameter nest.
              probably causes : nest too small. if nest is already
              large (say nest > m/2), it may also indicate that s is
              too small
              the approximation returned is the least-squares spline
              curve according to the knots t(1),t(2),...,t(n). (n=nest)
              the parameter fp gives the corresponding weighted sum of
              squared residuals (fp>s).
     ier=2  : error. a theoretically impossible result was found during
              the iteration proces for finding a smoothing spline curve
              with fp = s. probably causes : s too small.
              there is an approximation returned but the corresponding
              weighted sum of squared residuals does not satisfy the
              condition abs(fp-s)/s < tol.
     ier=3  : error. the maximal number of iterations maxit (set to 20
              by the program) allowed for finding a smoothing curve
              with fp=s has been reached. probably causes : s too small
              there is an approximation returned but the corresponding
              weighted sum of squared residuals does not satisfy the
              condition abs(fp-s)/s < tol.
     ier=10 : error. on entry, the input data are controlled on validity
              the following restrictions must be satisfied.
              -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
              0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
              nc>=nest*idim
              if ipar=0: sum j=1,idim (x(idim*i+j)-x(idim*(i-1)+j))**2>0
                         i=1,2,...,m-1.
              if ipar=1: ub<=u(1)<u(2)<...<u(m)<=ue
              if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
                          ub<t(k+2)<t(k+3)<...<t(n-k-1)<ue
                             (ub=0 and ue=1 in case ipar=0)
                        the schoenberg-whitney conditions, i.e. there
                        must be a subset of data points uu(j) such that
                          t(j) < uu(j) < t(j+k+1), j=1,2,...,n-k-1
              if iopt>=0: s>=0
                          if s=0 : nest >= m+k+1
              if one of these conditions is found to be violated,control
              is immediately repassed to the calling program. in that
              case there is no approximation returned.
 
   further comments:
    by means of the parameter s, the user can control the tradeoff
    between closeness of fit and smoothness of fit of the approximation.
    if s is too large, the curve will be too smooth and signal will be
    lost ; if s is too small the curve will pick up too much noise. in
    the extreme cases the program will return an interpolating curve if
    s=0 and the least-squares polynomial curve of degree k if s is
    very large. between these extremes, a properly chosen s will result
    in a good compromise between closeness of fit and smoothness of fit.
    to decide whether an approximation, corresponding to a certain s is
    satisfactory the user is highly recommended to inspect the fits
    graphically.
    recommended values for s depend on the weights w(i). if these are
    taken as 1/d(i) with d(i) an estimate of the standard deviation of
    x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
    sqrt(2*m)). if nothing is known about the statistical error in x(i)
    each w(i) can be set equal to one and s determined by trial and
    error, taking account of the comments above. the best is then to
    start with a very large value of s ( to determine the least-squares
    polynomial curve and the upper bound fp0 for s) and then to
    progressively decrease the value of s ( say by a factor 10 in the
    beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
    approximating curve shows more detail) to obtain closer fits.
    to economize the search for a good s-value the program provides with
    different modes of computation. at the first call of the routine, or
    whenever he wants to restart with the initial set of knots the user
    must set iopt=0.
    if iopt=1 the program will continue with the set of knots found at
    the last call of the routine. this will save a lot of computation
    time if parcur is called repeatedly for different values of s.
    the number of knots of the spline returned and their location will
    depend on the value of s and on the complexity of the shape of the
    curve underlying the data. but, if the computation mode iopt=1 is
    used, the knots returned may also depend on the s-values at previous
    calls (if these were smaller). therefore, if after a number of
    trials with different s-values and iopt=1, the user can finally
    accept a fit as satisfactory, it may be worthwhile for him to call
    parcur once more with the selected value for s but now with iopt=0.
    indeed, parcur may then return an approximation of the same quality
    of fit but with fewer knots and therefore better if data reduction
    is also an important objective for the user.
 
    the form of the approximating curve can strongly be affected by
    the choice of the parameter values u(i). if there is no physical
    reason for choosing a particular parameter u, often good results
    will be obtained with the choice of parcur (in case ipar=0), i.e.
         v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
    where
         q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
    other possibilities for q(i) are
         q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
         q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
         q(i)= max j=1,idim abs(xj(i)-xj(i-1))
         q(i)= 1
 
   other subroutines required:
     fpback,fpbspl,fpchec,fppara,fpdisc,fpgivs,fpknot,fprati,fprota
 
   references:
    dierckx p. : algorithms for smoothing data with periodic and
                 parametric splines, computer graphics and image
                 processing 20 (1982) 171-184.
    dierckx p. : algorithms for smoothing data with periodic and param-
                 etric splines, report tw55, dept. computer science,
                 k.u.leuven, 1981.
    dierckx p. : curve and surface fitting with splines, monographs on
                 numerical analysis, oxford university press, 1993.
 
   author:
     p.dierckx
     dept. computer science, k.u. leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be
 
   creation date : may 1979
   latest update : march 1987
   \endverbatim
* *
* \param	iopt			Option for how spline is computed.
* \param	ipar			If parameter values u[], ub and ue
* 					are supplied the must be 1, else 0.
* \param	idim			Curve dimension (0 < idim < 11).
* \param	m			Number of data points.
* \param	u			Array of parameter variables.
* \param	mx			Dimension of the array x.
* \param	x			Independent data.
* \param	w			Weights for data points.
* \param	ub			Lower bound of parameter values
* 					if ipar == 1.
* \param	ue			Upper bound of parameter values
* 					if ipar == 1.
* \param	k			Degree of spline.
* \param	s			Smoothing parameter.
* \param	nest			Over estimate of the number of knots
* 					to be returned (eg)
* \param	n			Destination pointer for the number of
* 					knots returned.
* \param	t			Array for returned computed knots.
* \param	nc			
* \param	c			Array for returned spline coefficients.
* \param	fp			Weighted sum of squared residuals of the
* 					spline approximation returned.
* \param	wrk			Workspace with minimum size
* 					m * (k + 1) + nest * (idim + 6 + k * 3).
* \param	iwrk			Workspace with minimum size
* 					nest * idim.
*/
AlgError			AlgBSplineNDFit(
				  int iopt,
				  int ipar,
				  int idim,
				  int m,
				  double *u,
				  int mx,
				  double *x,
				  double *w,
				  double ub,
				  double ue,
				  int k,
				  double s,
				  int nest,
				  int *n,
				  double *t,
				  int *nc,
				  double *c,
				  double *fp,
				  double *wrk,
				  int *iwrk)
{
    int 	k1,
    		k2,
		ncc,
		nmin,
    		maxit = 20;
    double	tol = .001;
    AlgError	errNum = ALG_ERR_NONE;

  if((iopt < -1) || (iopt > 1) ||
     (ipar < 0) || (ipar > 1) ||
     (idim <= 0) || (idim > 10) ||
     (k <= 0) || (k > 5))
  {
    errNum = ALG_ERR_FUNC;
  }
  if(errNum == ALG_ERR_NONE)
  {
    k1 = k + 1;
    k2 = k1 + 1;
    nmin = 2 * k1;
    ncc = nest * idim;
    if((m < k1) || (nest < nmin) || (mx < m * idim))
    {
      errNum = ALG_ERR_FUNC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int		i,
    		i1,
		i2;
    double	r1,
    		dist;

    i1 = 0;
    i2 = idim;
    u[0] = 0.0;
    for(i = 1; i < m; ++i)
    {
      int	j;

      dist = 0.0;
      for(j = 0; j < idim; ++j)
      {
	r1 = x[i2] - x[i1];
	dist += r1 * r1;
	++i1;
	++i2;
      }
      u[i] = u[i - 1] + sqrt(dist);
    }
    if (u[m - 1] <= 0.0)
    {
      errNum = ALG_ERR_FUNC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int		i;

    for(i = 1; i < m; ++i)
    {
      u[i] /= u[m - 1];
    }
    ub = 0.0;
    ue = 1.0;
    u[m - 1] = ue;
    if((ub > u[0]) || (ue < u[m - 1]) || (w[0] <= 0.0))
    {
      errNum = ALG_ERR_FUNC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int		i;

    for(i = 1; i < m; ++i)
    {
      if((u[i - 1] >= u[i]) || (w[i] <= 0.0))
      {
	errNum = ALG_ERR_FUNC;
	break;
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    if(iopt >= 0)
    {
      if((s < 0.f) || ((s == 0.0) && (nest < m + k1)))
      {
        errNum = ALG_ERR_FUNC;
      }
    }
    else
    {
      if((*n < nmin) || (*n > nest))
      {
	errNum = ALG_ERR_FUNC;
      }
      else
      {
	int	i,
		j,
		i1,
		ier = 0;

	j = *n - 1;
	i1 = k1;
	for(i = 0; i < i1; ++i)
	{
	  t[i] = ub;
	  t[j] = ue;
	  --j;
	}
	ier = AlgBSplineChec(u, m, t, *n, k);
	errNum = AlgErrorFromDierckx(ier);
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int		ia,
    		ib,
		ig,
		iq,
		iz,
		ifp = 1,
		ier = 0;

    iz = ifp + nest - 1;
    ia = iz + ncc;
    ib = ia + nest * k1;
    ig = ib + nest * k2;
    iq = ig + nest * k2;
    ier = AlgBSplinePara(iopt, idim, m, u, &mx, x, w,
            &ub, &ue, k, s, nest, tol, maxit,
	    k1, k2, n, t, &ncc, c, fp,
	    &wrk[ifp], &wrk[iz], &wrk[ia], &wrk[ib],
	    &wrk[ig], &wrk[iq], &iwrk[0]);
    errNum = AlgErrorFromDierckx(ier);
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgFit
*
* This function has been derived from the netlib Dierckx function percur().
* The original fortran coments are:
*
*  \verbatim
  Given the set of data points (x(i),y(i)) and the set of positive 
  numbers w(i),i=1,2,...,m-1, subroutine percur determines a smooth 
  periodic spline approximation of degree k with period per=x(m)-x(1). 
  if iopt=-1 percur calculates the weighted least-squares periodic 
  spline according to a given set of knots. 
  if iopt>=0 the number of knots of the spline s(x) and the position 
  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth- 
  ness of s(x) is then achieved by minimalizing the discontinuity 
  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,..., 
  n-k-1. the amount of smoothness is determined by the condition that 
  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non- 
  negative constant, called the smoothing factor. 
  the fit s(x) is given in the b-spline representation (b-spline coef- 
  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of 
  subroutine splev. 

  calling sequence: 
     call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk, lwrk,iwrk,ier) 

  parameters: 
   iopt  : integer flag. on entry iopt must specify whether a weighted 
           least-squares spline (iopt=-1) or a smoothing spline (iopt= 
           0 or 1) must be determined. if iopt=0 the routine will start 
           with an initial set of knots t(i)=x(1)+(x(m)-x(1))*(i-k-1), 
           i=1,2,...,2*k+2. if iopt=1 the routine will continue with 
           the knots found at the last call of the routine. 
           attention: a call with iopt=1 must always be immediately 
           preceded by another call with iopt=1 or iopt=0. 
           unchanged on exit. 
   m     : integer. on entry m must specify the number of data points. 
           m > 1. unchanged on exit. 
   x     : real array of dimension at least (m). before entry, x(i) 
           must be set to the i-th value of the independent variable x, 
           for i=1,2,...,m. these values must be supplied in strictly 
           ascending order. x(m) only indicates the length of the 
           period of the spline, i.e per=x(m)-x(1). 
           unchanged on exit. 
   y     : real array of dimension at least (m). before entry, y(i) 
           must be set to the i-th value of the dependent variable y, 
           for i=1,2,...,m-1. the element y(m) is not used. 
           unchanged on exit. 
   w     : real array of dimension at least (m). before entry, w(i) 
           must be set to the i-th value in the set of weights. the 
           w(i) must be strictly positive. w(m) is not used. 
           see also further comments. unchanged on exit. 
   k     : integer. on entry k must specify the degree of the spline. 
           1<=k<=5. it is recommended to use cubic splines (k=3). 
           the user is strongly dissuaded from choosing k even,together 
           with a small s-value. unchanged on exit. 
   s     : real.on entry (in case iopt>=0) s must specify the smoothing 
           factor. s >=0. unchanged on exit. 
           for advice on the choice of s see further comments. 
   nest  : integer. on entry nest must contain an over-estimate of the 
           total number of knots of the spline returned, to indicate 
           the storage space available to the routine. nest >=2*k+2. 
           in most practical situation nest=m/2 will be sufficient. 
           always large enough is nest=m+2*k,the number of knots needed 
           for interpolation (s=0). unchanged on exit. 
   n     : integer. 
           unless ier = 10 (in case iopt >=0), n will contain the 
           total number of knots of the spline approximation returned. 
           if the computation mode iopt=1 is used this value of n 
           should be left unchanged between subsequent calls. 
           in case iopt=-1, the value of n must be specified on entry. 
   t     : real array of dimension at least (nest). 
           on succesful exit, this array will contain the knots of the 
           spline,i.e. the position of the interior knots t(k+2),t(k+3) 
           ...,t(n-k-1) as well as the position of the additional knots 
           t(1),t(2),...,t(k+1)=x(1) and t(n-k)=x(m),..,t(n) needed for 
           the b-spline representation. 
           if the computation mode iopt=1 is used, the values of t(1), 
           t(2),...,t(n) should be left unchanged between subsequent 
           calls. if the computation mode iopt=-1 is used, the values 
           t(k+2),...,t(n-k-1) must be supplied by the user, before 
           entry. see also the restrictions (ier=10). 
   c     : real array of dimension at least (nest). 
           on succesful exit, this array will contain the coefficients 
           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x) 
   fp    : real. unless ier = 10, fp contains the weighted sum of 
           squared residuals of the spline approximation returned. 
   wrk   : real array of dimension at least (m*(k+1)+nest*(8+5*k)). 
           used as working space. if the computation mode iopt=1 is 
           used, the values wrk(1),...,wrk(n) should be left unchanged 
           between subsequent calls. 
   lwrk  : integer. on entry,lwrk must specify the actual dimension of 
           the array wrk as declared in the calling (sub)program. lwrk 
           must not be too small (see wrk). unchanged on exit. 
   iwrk  : integer array of dimension at least (nest). 
           used as working space. if the computation mode iopt=1 is 
           used,the values iwrk(1),...,iwrk(n) should be left unchanged 
           between subsequent calls. 
   ier   : integer. unless the routine detects an error, ier contains a 
           non-positive value on exit, i.e. 
    ier=0  : normal return. the spline returned has a residual sum of 
             squares fp such that abs(fp-s)/s <= tol with tol a relat- 
             ive tolerance set to 0.001 by the program. 
    ier=-1 : normal return. the spline returned is an interpolating 
             periodic spline (fp=0). 
    ier=-2 : normal return. the spline returned is the weighted least- 
             squares constant. in this extreme case fp gives the upper 
             bound fp0 for the smoothing factor s. 
    ier=1  : error. the required storage space exceeds the available 
             storage space, as specified by the parameter nest. 
             probably causes : nest too small. if nest is already 
             large (say nest > m/2), it may also indicate that s is 
             too small 
             the approximation returned is the least-squares periodic 
             spline according to the knots t(1),t(2),...,t(n). (n=nest) 
             the parameter fp gives the corresponding weighted sum of 
             squared residuals (fp>s). 
    ier=2  : error. a theoretically impossible result was found during 
             the iteration proces for finding a smoothing spline with 
             fp = s. probably causes : s too small. 
             there is an approximation returned but the corresponding 
             weighted sum of squared residuals does not satisfy the 
             condition abs(fp-s)/s < tol. 
    ier=3  : error. the maximal number of iterations maxit (set to 20 
             by the program) allowed for finding a smoothing spline 
             with fp=s has been reached. probably causes : s too small 
             there is an approximation returned but the corresponding 
             weighted sum of squared residuals does not satisfy the 
             condition abs(fp-s)/s < tol. 
    ier=10 : error. on entry, the input data are controlled on validity 
             the following restrictions must be satisfied. 
             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,...,m-1 
             x(1)<x(2)<...<x(m), lwrk>=(k+1)*m+nest*(8+5*k) 
             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k) 
                         x(1)<t(k+2)<t(k+3)<...<t(n-k-1)<x(m) 
                       the schoenberg-whitney conditions, i.e. there 
                       must be a subset of data points xx(j) with 
                       xx(j) = x(i) or x(i)+(x(m)-x(1)) such that 
                         t(j) < xx(j) < t(j+k+1), j=k+1,...,n-k-1 
             if iopt>=0: s>=0 
                         if s=0 : nest >= m+2*k 
             if one of these conditions is found to be violated,control 
             is immediately repassed to the calling program. in that 
             case there is no approximation returned. 

  further comments: 
   by means of the parameter s, the user can control the tradeoff 
   between closeness of fit and smoothness of fit of the approximation. 
   if s is too large, the spline will be too smooth and signal will be 
   lost ; if s is too small the spline will pick up too much noise. in 
   the extreme cases the program will return an interpolating periodic 
   spline if s=0 and the weighted least-squares constant if s is very 
   large. between these extremes, a properly chosen s will result in 
   a good compromise between closeness of fit and smoothness of fit. 
   to decide whether an approximation, corresponding to a certain s is 
   satisfactory the user is highly recommended to inspect the fits 
   graphically. 
   recommended values for s depend on the weights w(i). if these are 
   taken as 1/d(i) with d(i) an estimate of the standard deviation of 
   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+ 
   sqrt(2*m)). if nothing is known about the statistical error in y(i) 
   each w(i) can be set equal to one and s determined by trial and 
   error, taking account of the comments above. the best is then to 
   start with a very large value of s ( to determine the least-squares 
   constant and the corresponding upper bound fp0 for s) and then to 
   progressively decrease the value of s ( say by a factor 10 in the 
   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the 
   approximation shows more detail) to obtain closer fits. 
   to economize the search for a good s-value the program provides with 
   different modes of computation. at the first call of the routine, or 
   whenever he wants to restart with the initial set of knots the user 
   must set iopt=0. 
   if iopt=1 the program will continue with the set of knots found at 
   the last call of the routine. this will save a lot of computation 
   time if percur is called repeatedly for different values of s. 
   the number of knots of the spline returned and their location will 
   depend on the value of s and on the complexity of the shape of the 
   function underlying the data. but, if the computation mode iopt=1 
   is used, the knots returned may also depend on the s-values at 
   previous calls (if these were smaller). therefore, if after a number 
   of trials with different s-values and iopt=1, the user can finally 
   accept a fit as satisfactory, it may be worthwhile for him to call 
   percur once more with the selected value for s but now with iopt=0. 
   indeed, percur may then return an approximation of the same quality 
   of fit but with fewer knots and therefore better if data reduction 
   is also an important objective for the user. 

  other subroutines required: 
    fpbacp,fpbspl,fpchep,fpperi,fpdisc,fpgivs,fpknot,fprati,fprota 

  references: 
   dierckx p. : algorithms for smoothing data with periodic and 
                parametric splines, computer graphics and image 
                processing 20 (1982) 171-184. 
   dierckx p. : algorithms for smoothing data with periodic and param- 
                etric splines, report tw55, dept. computer science, 
                k.u.leuven, 1981. 
   dierckx p. : curve and surface fitting with splines, monographs on 
                numerical analysis, oxford university press, 1993. 

  author: 
    p.dierckx 
    dept. computer science, k.u. leuven 
    celestijnenlaan 200a, b-3001 heverlee, belgium. 
    e-mail : Paul.Dierckx@cs.kuleuven.ac.be 

  creation date : may 1979 
  latest update : march 1987 
  \endverbatim
*
* \brief	Computes a smooth periodic B-spline approximation
* 		through the given data.
* \param	iopt		Option for how spline is computed.
* \param	m		Number of data points.
* \param	x		Independent data.
* \param	y		Dependent data.
* \param	w		Weights for data points.
* \param	k		Degree of spline.
* \param	s		Smoothing parameter.
* \param	nest		Over estimate of the number of knots
* 				to be returned (eg m + 2 * (k + 1)).
* \param	n		Destination pointer for the number of
* 				knots returned.
* \param	t		Array for returned computed knots.
* \param	c		Array for returned spline coefficients.
* \param	fp		Weighted sum of squared residuals of the
* 				spline approximation returned.
* \param	wrk		Workspace with minimum size
* 				(m * (k + 1) + nest * (8 + 5 * k)).
* \param	iwrk		Workspace with minimum size nest.
*/
AlgError			AlgBSplinePerFit(
				  int iopt,
				  int m,
				  double *x,
				  double *y,
				  double *w,
				  int k,
				  double s,
				  int nest,
				  int *n,
				  double *t, 
				  double *c,
				  double *fp,
				  double *wrk,
				  int *iwrk)
{
  int		k1,
  		k2,
      		nmin,
  		maxit = 20;
  double 	tol = 0.001;
  AlgError	errNum = ALG_ERR_NONE;

  if((iopt < -1) || (iopt > 1) || (k <= 0) || (k > 5))
  {
    errNum = ALG_ERR_FUNC;
  }
  if(errNum == ALG_ERR_NONE)
  {
    k1 = k + 1;
    k2 = k1 + 1;
    nmin = 2 * k1;
    if((m < 2) || (nest < nmin))
    {
      errNum = ALG_ERR_FUNC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int		i,
    		m1;

    m1 = m - 1;
    for(i = 0; i < m1; ++i)
    {
      if(x[i] >= x[i + 1] || w[i] <= 0.f)
      {
	errNum = ALG_ERR_FUNC;
	break;
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    if (iopt >= 0)
    {
      if((s < 0.f) || (s == 0.f && nest < m + (2 * k)))
      {
        errNum = ALG_ERR_FUNC;
      }
    }
    else
    {
      if((*n <= nmin) || (*n > nest))
      {
	errNum = ALG_ERR_FUNC;
      }
      else
      {
	int	i,
		i1,
		i2,
		j1,
		j2;
        double	per;

	per = x[m - 1] - x[0];
	j1 = k1;
	t[j1 - 1] = x[0];
	i1 = *n - k;
	t[i1 - 1] = x[m - 1];
	j2 = j1;
	i2 = i1;
	for (i = 0; i < k; ++i) {
	    ++i1;
	    --i2;
	    ++j1;
	    --j2;
	    t[j2 - 1] = t[i2 - 1] - per;
	    t[i1 - 1] = t[j1 - 1] + per;
	}
	if(AlgBSplineChep(x, m, t, *n, k))
	{
	  errNum = ALG_ERR_FUNC;
	}
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int		iz,
		ia1,
		ia2,
		ib,
		ig1,
		ig2,
		iq,
		ier = 0;

    /* Partition the working space and determine the spline approximation. */
    ier = 0;
    iz = nest;
    ia1 = iz + nest;
    ia2 = ia1 + nest * k1;
    ib = ia2 + nest * k;
    ig1 = ib + nest * k2;
    ig2 = ig1 + nest * k2;
    iq = ig2 + nest * k1;
    ier = AlgBSplinePeri(iopt, x, y, w, m, k, s, nest, tol, maxit,
                 k1, k2, n, t, c, fp,
		 wrk, &wrk[iz], &wrk[ia1], &wrk[ia2], 
	         &wrk[ib], &wrk[ig1], &wrk[ig2], &wrk[iq], iwrk);
    errNum = AlgErrorFromDierckx(ier);
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgFit
* \brief	Computes a smooth B-spline approximation through the
* 		given data.
*
* This function has been derived from the netlib Dierckx function curfit().
* The original fortran coments are:
*
*  \verbatim
   given the set of data points (x(i),y(i)) and the set of positive
   numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
   approximation of degree k on the interval xb <= x <= xe.
   if iopt=-1 curfit calculates the weighted least-squares spline
   according to a given set of knots.
   if iopt>=0 the number of knots of the spline s(x) and the position
   t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
   ness of s(x) is then achieved by minimalizing the discontinuity
   jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
   n-k-1. the amount of smoothness is determined by the condition that
   f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
   negative constant, called the smoothing factor.
   the fit s(x) is given in the b-spline representation (b-spline coef-
   ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
   subroutine splev.
 
   calling sequence:
      call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
     * lwrk,iwrk,ier)
 
   parameters:
    iopt  : integer flag. on entry iopt must specify whether a weighted
            least-squares spline (iopt=-1) or a smoothing spline (iopt=
            0 or 1) must be determined. if iopt=0 the routine will start
            with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
            k+1. if iopt=1 the routine will continue with the knots
            found at the last call of the routine.
            attention: a call with iopt=1 must always be immediately
            preceded by another call with iopt=1 or iopt=0.
            unchanged on exit.
    m     : integer. on entry m must specify the number of data points.
            m > k. unchanged on exit.
    x     : real array of dimension at least (m). before entry, x(i)
            must be set to the i-th value of the independent variable x,
            for i=1,2,...,m. these values must be supplied in strictly
            ascending order. unchanged on exit.
    y     : real array of dimension at least (m). before entry, y(i)
            must be set to the i-th value of the dependent variable y,
            for i=1,2,...,m. unchanged on exit.
    w     : real array of dimension at least (m). before entry, w(i)
            must be set to the i-th value in the set of weights. the
            w(i) must be strictly positive. unchanged on exit.
            see also further comments.
    xb,xe : real values. on entry xb and xe must specify the boundaries
            of the approximation interval. xb<=x(1), xe>=x(m).
            unchanged on exit.
    k     : integer. on entry k must specify the degree of the spline.
            1<=k<=5. it is recommended to use cubic splines (k=3).
            the user is strongly dissuaded from choosing k even,together
            with a small s-value. unchanged on exit.
    s     : real.on entry (in case iopt>=0) s must specify the smoothing
            factor. s >=0. unchanged on exit.
            for advice on the choice of s see further comments.
    nest  : integer. on entry nest must contain an over-estimate of the
            total number of knots of the spline returned, to indicate
            the storage space available to the routine. nest >=2*k+2.
            in most practical situation nest=m/2 will be sufficient.
            always large enough is  nest=m+k+1, the number of knots
            needed for interpolation (s=0). unchanged on exit.
    n     : integer.
            unless ier =10 (in case iopt >=0), n will contain the
            total number of knots of the spline approximation returned.
            if the computation mode iopt=1 is used this value of n
            should be left unchanged between subsequent calls.
            in case iopt=-1, the value of n must be specified on entry.
    t     : real array of dimension at least (nest).
            on succesful exit, this array will contain the knots of the
            spline,i.e. the position of the interior knots t(k+2),t(k+3)
            ...,t(n-k-1) as well as the position of the additional knots
            t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
            the b-spline representation.
            if the computation mode iopt=1 is used, the values of t(1),
            t(2),...,t(n) should be left unchanged between subsequent
            calls. if the computation mode iopt=-1 is used, the values
            t(k+2),...,t(n-k-1) must be supplied by the user, before
            entry. see also the restrictions (ier=10).
    c     : real array of dimension at least (nest).
            on succesful exit, this array will contain the coefficients
            c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
    fp    : real. unless ier=10, fp contains the weighted sum of
            squared residuals of the spline approximation returned.
    wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
            used as working space. if the computation mode iopt=1 is
            used, the values wrk(1),...,wrk(n) should be left unchanged
            between subsequent calls.
    lwrk  : integer. on entry,lwrk must specify the actual dimension of
            the array wrk as declared in the calling (sub)program.lwrk
            must not be too small (see wrk). unchanged on exit.
    iwrk  : integer array of dimension at least (nest).
            used as working space. if the computation mode iopt=1 is
            used,the values iwrk(1),...,iwrk(n) should be left unchanged
            between subsequent calls.
    ier   : integer. unless the routine detects an error, ier contains a
            non-positive value on exit, i.e.
     ier=0  : normal return. the spline returned has a residual sum of
              squares fp such that abs(fp-s)/s <= tol with tol a relat-
              ive tolerance set to 0.001 by the program.
     ier=-1 : normal return. the spline returned is an interpolating
              spline (fp=0).
     ier=-2 : normal return. the spline returned is the weighted least-
              squares polynomial of degree k. in this extreme case fp
              gives the upper bound fp0 for the smoothing factor s.
     ier=1  : error. the required storage space exceeds the available
              storage space, as specified by the parameter nest.
              probably causes : nest too small. if nest is already
              large (say nest > m/2), it may also indicate that s is
              too small
              the approximation returned is the weighted least-squares
              spline according to the knots t(1),t(2),...,t(n). (n=nest)
              the parameter fp gives the corresponding weighted sum of
              squared residuals (fp>s).
     ier=2  : error. a theoretically impossible result was found during
              the iteration proces for finding a smoothing spline with
              fp = s. probably causes : s too small.
              there is an approximation returned but the corresponding
              weighted sum of squared residuals does not satisfy the
              condition abs(fp-s)/s < tol.
     ier=3  : error. the maximal number of iterations maxit (set to 20
              by the program) allowed for finding a smoothing spline
              with fp=s has been reached. probably causes : s too small
              there is an approximation returned but the corresponding
              weighted sum of squared residuals does not satisfy the
              condition abs(fp-s)/s < tol.
     ier=10 : error. on entry, the input data are controlled on validity
              the following restrictions must be satisfied.
              -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
              xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
              if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
                          xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
                        the schoenberg-whitney conditions, i.e. there
                        must be a subset of data points xx(j) such that
                          t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
              if iopt>=0: s>=0
                          if s=0 : nest >= m+k+1
              if one of these conditions is found to be violated,control
              is immediately repassed to the calling program. in that
              case there is no approximation returned.
 
   further comments:
    by means of the parameter s, the user can control the tradeoff
    between closeness of fit and smoothness of fit of the approximation.
    if s is too large, the spline will be too smooth and signal will be
    lost ; if s is too small the spline will pick up too much noise. in
    the extreme cases the program will return an interpolating spline if
    s=0 and the weighted least-squares polynomial of degree k if s is
    very large. between these extremes, a properly chosen s will result
    in a good compromise between closeness of fit and smoothness of fit.
    to decide whether an approximation, corresponding to a certain s is
    satisfactory the user is highly recommended to inspect the fits
    graphically.
    recommended values for s depend on the weights w(i). if these are
    taken as 1/d(i) with d(i) an estimate of the standard deviation of
    y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
    sqrt(2*m)). if nothing is known about the statistical error in y(i)
    each w(i) can be set equal to one and s determined by trial and
    error, taking account of the comments above. the best is then to
    start with a very large value of s ( to determine the least-squares
    polynomial and the corresponding upper bound fp0 for s) and then to
    progressively decrease the value of s ( say by a factor 10 in the
    beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
    approximation shows more detail) to obtain closer fits.
    to economize the search for a good s-value the program provides with
    different modes of computation. at the first call of the routine, or
    whenever he wants to restart with the initial set of knots the user
    must set iopt=0.
    if iopt=1 the program will continue with the set of knots found at
    the last call of the routine. this will save a lot of computation
    time if curfit is called repeatedly for different values of s.
    the number of knots of the spline returned and their location will
    depend on the value of s and on the complexity of the shape of the
    function underlying the data. but, if the computation mode iopt=1
    is used, the knots returned may also depend on the s-values at
    previous calls (if these were smaller). therefore, if after a number
    of trials with different s-values and iopt=1, the user can finally
    accept a fit as satisfactory, it may be worthwhile for him to call
    curfit once more with the selected value for s but now with iopt=0.
    indeed, curfit may then return an approximation of the same quality
    of fit but with fewer knots and therefore better if data reduction
    is also an important objective for the user.
 
   other subroutines required:
     fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
 
   references:
    dierckx p. : an algorithm for smoothing, differentiation and integ-
                 ration of experimental data using spline functions,
                 j.comp.appl.maths 1 (1975) 165-184.
    dierckx p. : a fast algorithm for smoothing data on a rectangular
                 grid while using spline functions, siam j.numer.anal.
                 19 (1982) 1286-1304.
    dierckx p. : an improved algorithm for curve fitting with spline
                 functions, report tw54, dept. computer science,k.u.
                 leuven, 1981.
    dierckx p. : curve and surface fitting with splines, monographs on
                 numerical analysis, oxford university press, 1993.
 
   author:
     p.dierckx
     dept. computer science, k.u. leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be
 
   creation date : may 1979
   latest update : march 1987
   \endverbatim
*
* \param	iopt			Option for how spline is computed.
* \param	m			Number of data points.
* \param	x			Independent data.
* \param	y			Dependent data.
* \param	w			Weights for data points.
* \param	xb
* \param	xe
* \param	k			Degree of spline.
* \param	s			Smoothing parameter.
* \param	nest			Over estimate of the number of knots
* 					to be returned (eg m + k + 1).
* \param	n			Destination pointer for the number of
* 					knots returned.
* \param	t			Array for returned computed knots.
* \param	c			Array for returned spline coefficients.
* \param	fp			Weighted sum of squared residuals of the
* 					spline approximation returned.
* \param	wrk			Workspace with minimum size
* 					m * k1 + nest * (k * 3 + 7)
* \param	iwrk			Workspace with minimum size nest.
*/
AlgError			AlgBSplineFit(
				  int iopt,
				  int m,
				  double *x,
				  double *y,
				  double *w,
				  double xb,
				  double xe,
				  int k,
				  double s,
				  int nest,
				  int *n,
				  double *t,
				  double *c,
				  double *fp,
				  double *wrk,
				  int *iwrk)
{
  int		k1,
  		k2,
		nmin,
		maxit = 20;
  double 	tol = 0.001f;
  AlgError	errNum = ALG_ERR_NONE;

  if((iopt < -1) || (iopt > 1) || (k <= 0) || (k > 5))
  {
    errNum = ALG_ERR_FUNC;
  }
  else
  {
    k1 = k + 1;
    k2 = k1 + 1;
    nmin = 2 * k1;
    if((m < k1) || (nest < nmin) ||
       (xb > x[0]) || (xe < x[m - 1]) || (w[0] <= 0.0))
    {
      errNum = ALG_ERR_FUNC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int	i;

    for(i = 1; i < m; ++i)
    {
      if((x[i - 1] >= x[i]) || (w[i] <= 0.0))
      {
	errNum = ALG_ERR_FUNC;
	break;
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    if(iopt >= 0)
    {
      if((s < 0.0) || ((s == 0.0) && (nest < m + k1)))
      {
        errNum = ALG_ERR_FUNC;
      }
    }
    else
    {
      if((*n < nmin) || (*n > nest))
      {
        errNum = ALG_ERR_FUNC;
      }
      else
      {
	int	i,
	      	j,
	      	ier = 0;

	j = *n - 1;
	for(i = 0; i < k1; ++i)
	{
	  t[i] = xb;
	  t[j] = xe;
	  --j;
	}
	ier = AlgBSplineChec(x, m, t, *n, k);
	errNum = AlgErrorFromDierckx(ier);
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int 	ia,
    		ib,
		ig,
		iq,
		iz,
    		ier = 0;

    /* Partition the working space and determine the spline approximation. */
    iz = nest;
    ia = iz + nest;
    ib = ia + nest * k1;
    ig = ib + nest * k2;
    iq = ig + nest * k2;
    ier = AlgBSplineCurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol,
	      maxit, k1, k2, n, t, c, fp, wrk, &wrk[iz],
	      &wrk[ia], &wrk[ib], &wrk[ig], &wrk[iq], iwrk);
    errNum = AlgErrorFromDierckx(ier);
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgFit
* \brief	Evaluates a b-spline at a number of points.
*
* This function has been derived from the netlib Dierckx function splev().
* The original fortran coments are:
*
*  \verbatim
   Subroutine splev evaluates in a number of points x(i),i=1,2,...,m
   a spline s(x) of degree k, given in its b-spline representation.
 
   calling sequence:
      call splev(t,n,c,k,x,y,m,ier)
 
   input parameters:
     t    : array,length n, which contains the position of the knots.
     n    : integer, giving the total number of knots of s(x).
     c    : array,length n, which contains the b-spline coefficients.
     k    : integer, giving the degree of s(x).
     x    : array,length m, which contains the points where s(x) must
            be evaluated.
     m    : integer, giving the number of points where s(x) must be
            evaluated.
 
   output parameter:
     y    : array,length m, giving the value of s(x) at the different
            points.
     ier  : error flag
       ier = 0 : normal return
       ier =10 : invalid input data (see restrictions)
 
   restrictions:
     m >= 1
     t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
 
   other subroutines required: fpbspl.
 
   references :
     de boor c  : on calculating with b-splines, j. approximation theory
                  6 (1972) 50-62.
     cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
                  applics 10 (1972) 134-149.
     dierckx p. : curve and surface fitting with splines, monographs on
                  numerical analysis, oxford university press, 1993.
 
   author :
     p.dierckx
     dept. computer science, k.u.leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be
 
   latest update : march 1987
   \endverbatim
*
* \param	t			Knot positions.
* \param	n			Number of knots.
* \param	c			B-spline coefficients.
* \param	k			Degree of the B-spline.
* \param	x			Points at which the B-spline must
* 					be evaluated.
* \param	y			Array for evaluated B-spline at given
* 					points.
* \param	m			Number of points at which B-spline is
* 					to be evaluated.
*/
AlgError 			AlgBSplineEval(
				  double *t,
				  int n,
				  double *c,
				  int k,
				  double * x,
				  double *y,
				  int m)
{
  AlgError	errNum = ALG_ERR_NONE;

  if (m <= 0)
  {
    errNum = ALG_ERR_FUNC;
  }
  else if (m > 1)
  {
    int		i,
    		m1;

    m1 = m - 1;
    for(i = 2; i <= m1; ++i)
    {
      if(x[i - 1] < x[i - 2])
      {
	errNum = ALG_ERR_FUNC;
	break;
      }
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    int		i,
    		k1,
		l,
		l1,
		nk1;
    double	tb,
    		te;

    /* Fetch tb and te, the boundaries of the approximation interval. */
    k1 = k + 1;
    nk1 = n - k1;
    tb = t[k];
    te = t[nk1];
    l = k1;
    l1 = l + 1;
    /*  Main loop for the different points. */
    for(i = 0; i < m; ++i)
    {
      int	j,
      		ll;
      double	sp,
      		arg;
      double 	h[6];

      /* Fetch a new x-value arg. */
      arg = x[i];
      arg = ALG_CLAMP(arg, tb, te);
      /* Search for knot interval t(l-1) <= arg < t(l+1-1) */
      while(!(arg < t[l]) && (l != nk1))
      {
        l = l1;
        ++l1;
      }
      AlgBSplineBspl(t, k, arg, l, h);
      /* Find the value of s(x=arg). */
      sp = 0.f;
      ll = l - (k + 2);
      for (j = 0; j < k1; ++j) {
	++ll;
	sp += c[ll] * h[j];
      }
      y[i] = sp;
    }
  }
  return(errNum);
}

/*!
* \return	Alg error code.
* \ingroup	AlgFit
* \brief	Returns the closest matching Alg error code for the given
* 		Dierckx error code.
* \param	ier			Dierckx error code.
*/
static AlgError			AlgErrorFromDierckx(int ier)
{
  AlgError	errNum = ALG_ERR_NONE;

  switch(ier)
  {
    case 0:   /* FALLTHROUGH */
    case -1:  /* FALLTHROUGH */
    case -2:
      break;
    case 1:
      errNum = ALG_ERR_MALLOC;
      break;
    case 3:
      errNum = ALG_ERR_CONVERGENCE;
      break;
    case 10:
      errNum = ALG_ERR_FUNC;
      break;
    case 2:  /* FALLTHROUGH */
    default:
      errNum = ALG_ERR_UNKNOWN;
      break;
  }
  return(errNum);
}

/*!
* \ingroup	ALgFit
* \brief	Adds an additional knot to a B-Spline.
*
* This function has been derived from the netlib Dierckx function fpknot().
* The original fortran coments are:
*
*  \verbatim
   Subroutine fpknot locates an additional knot for a spline of degree
   k and adjusts the corresponding parameters,i.e.
     t     : the position of the knots.
     n     : the number of knots.
     nrint : the number of knotintervals.
     fpint : the sum of squares of residual right hand sides
             for each knot interval.
     nrdata: the number of data points inside each knot interval.
   istart indicates that the smallest data point at which the new knot
   may be added is x(istart+1)
   \endverbatim
*
* \param	x			Data points
* \param	m			Number of data points.
* \param	t			Positions of the knots.
* \param	n			Number of knots.
* \param	fpint			Sum of squares of residual right hand
* 					sides for each knot interval.
* \param	nrdata			Number of data points inside each
* 					knot interval.
* \param	nrint			Number of knot intervals.
* \param	istart			Smallest data point at which the new
* 					knot may be added.
*/
static void			AlgBSplineKnot(
				  double *x,
				  int m,
				  double *t,
				  int *n,
				  double *fpint,
				  int *nrdata,
				  int *nrint,
				  int istart)
{
  int		j,
  		k,
		jk,
		ni,
		next,
		ihalf,
		maxpt,
  		jbegin,
		maxbeg,
		number,
		jpoint,
		nrx;
  double	am,
  		an,
		fpmax;


  ni = *nrint;
  k = (*n - ni - 1) / 2;
  /*  Search for knot interval t(number+k) <= x <= t(number+k+1) where
   *  fpint(number) is maximal on the condition that nrdata(number) not
   *  equals zero. */
  fpmax = 0.0;
  jbegin = istart;
  for(j = 0; j < ni; ++j)
  {
    jpoint = nrdata[j];
    if((fpmax < fpint[j]) && (jpoint != 0))
    {
      fpmax = fpint[j];
      number = j + 1;
      maxpt = jpoint;
      maxbeg = jbegin;
    }
    jbegin = jbegin + jpoint + 1;
  }
  /* Let coincide the new knot t(number+k+1) with a data point x(nrx) inside
   * the old knot interval t(number+k) <= x <= t(number+k+1). */
  ihalf = maxpt / 2 + 1;
  nrx = maxbeg + ihalf;
  next = number + 1;
  if(next <= ni)
  {
    /* Adjust the different parameters. */
    for(j = next; j <= ni; ++j) {
      int	jj,
		jk;

      jj = next + ni - j;
      fpint[jj] = fpint[jj - 1];
      nrdata[jj] = nrdata[jj - 1];
      jk = jj + k;
      t[jk] = t[jk - 1];
    }
  }
  nrdata[number - 1] = ihalf - 1;
  nrdata[next - 1] = maxpt - ihalf;
  am = (double) maxpt;
  an = (double) nrdata[number - 1];
  fpint[number - 1] = fpmax * an / am;
  an = (double) nrdata[next - 1];
  fpint[next - 1] = fpmax * an / am;
  jk = next + k;
  t[jk - 1] = x[nrx - 1];
  ++(*n);
  ++(*nrint);
}


/*!
* \return	Zero if the B-Splines knots are valid, otherwise non-zero.
* \ingroup	AlgFit
* \brief	Verifies the number and the position of the knots of a
* 		periodic spline.
* 
* This function has been derived from the netlib Dierckx function fpchep().
* The original fortran coments are:
*
*  \verbatim
   Subroutine fpchep verifies the number and the position of the knots
   t(j),j=1,2,...,n of a periodic spline of degree k, in relation to
   the number and the position of the data points x(i),i=1,2,...,m.
   if all of the following conditions are fulfilled, ier is set
   to zero. if one of the conditions is violated ier is set to ten.
       1) k+1 <= n-k-1 <= m+k-1
       2) t(1) <= t(2) <= ... <= t(k+1)
          t(n-k) <= t(n-k+1) <= ... <= t(n)
       3) t(k+1) < t(k+2) < ... < t(n-k)
       4) t(k+1) <= x(i) <= t(n-k)
       5) the conditions specified by schoenberg and whitney must hold
          for at least one subset of data points, i.e. there must be a
          subset of data points y(j) such that
              t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1
   \endverbatim
*
* \param	x			Data points.	
* \param	m			Number of data points.
* \param	t			Knots.
* \param	n			Number of knots.
* \param	k			B-Spline order.
*/
static int			AlgBSplineChep(
				  double *x,
				  int m,
				  double *t,
				  int n,
				  int k)
{
  int		i,
  		k1,
		m1,
		nk1,
		nk2,
		ier = 0;

  k1 = k + 1;
  nk1 = n - k1;
  nk2 = nk1 + 1;
  m1 = m - 1;
  /* Check condition no 1 */
  if((nk1 < k1) || (n > m + (2 * k)))
  {
    ier = 10;
  }
  /* Check condition no 2 */
  if(ier == 0)
  {
    int		j;

    j = n - 1;
    for(i = 0; i < k; ++i)
    {
      if((t[i] > t[i + 1]) || (t[j] < t[j - 1]))
      {
        ier = 10;
	break;
      }
      --j;
    }
  }
  /* Check condition no 3 */
  if(ier == 0)
  {
    for(i = k1; i < nk2; ++i)
    {
      if(t[i] <= t[i - 1])
      {
	ier = 10;
	break;
      }
    }
  }
  /* Check condition no 4 */
  if((x[0] < t[k1 - 1]) || (x[m - 1] > t[nk1]))
  {
    ier = 10;
  }
  /* Check condition no 5 */
  if(ier == 0)
  {
    int		l,
		l1,
		l2,
    		i1;
    double	xi,
    		per;

    ier = 10;
    l1 = k1;
    l2 = 1;
    for(l = 1; l <= m; ++l)
    {
      xi = x[l - 1];
L40:
      if((xi >= t[l1 + 1 - 1]) && (l != nk1))
      {
	++l1;
	++l2;
	if(l2 > k1)
	{
	  goto L60;
	}
	goto L40;
      }
    }
    l = m;
L60:
    per = t[nk1] - t[k];
    for(i1 = 2; i1 <= l; ++i1)
    {
      int	j,
      		mm;
      double	tj,
      		tl;

      i = i1 - 1;
      mm = i + m1;
      for(j = k1; j <= nk1; ++j)
      {
	int	i2;

	tj = t[j - 1];
	tl = t[j + k];
L70:
	++i;
	if(i > mm)
	{
	  goto L120;
	}
	i2 = i - m1;
	if(i2 <= 0)
	{
	  goto L80;
	}
	else
	{
	  goto L90;
	}
L80:
	xi = x[i - 1];
	goto L100;
L90:
	xi = x[i2 - 1] + per;
L100:
	if(xi <= tj)
	{
	  goto L70;
	}
	if(xi >= tl)
	{
	  goto L120;
	}
      }
      ier = 0;
      goto L130;
L120:
      ;
    }
  }
L130:
  return(ier);
}

/*!
* \ingroup	AlgFit
* \brief	Calculates the discontinuity jumps of the kth derivative
* 		of the b-splines of degree k at the knots t[k+2]..t[n-k-1].
*
* This function has been derived from the netlib Dierckx function fpdisc().
* The original fortran coments are:
*
*  \verbatim
   subroutine fpdisc calculates the discontinuity jumps of the kth
   derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
   \endverbatim
*
* \param	t			Knot array.
* \param	n			Number of knots.
* \param	k2			Degree ofspline.
* \param	b			Array of basis functions b[nest][k2].
* \param	nest			Minimum number of knots estimate.
*/
static void			AlgBSplineDisc(
				  double *t,
				  int n,
				  int k2,
				  double *b,
				  int nest)
{
  int		k,
		l,
  		k1,
  		nk1,
		nrint;
  double 	fac;
  double	h[12];

  k1 = k2 - 1;
  k = k1 - 1;
  nk1 = n - k1;
  nrint = nk1 - k;
  fac = (double )nrint / (t[nk1] - t[k]);
  for(l = k1; l < nk1; ++l)
  {
    int		j,
    		lp,
    		lmk;

    lmk = l - k - 1;
    for(j = 0; j < k1; ++j)
    {
      int	ik,
      		lj,
		lk;

      ik = j + k1;
      lj = l + j + 1;
      lk = lj - k2;
      h[j] = t[l] - t[lk];
      h[ik] = t[l] - t[lj];
    }
    lp = lmk;
    for(j = 0; j < k2; ++j)
    {
      int	i,
      		lk,
      		jk,
      		lmkb;
      double	prod;

      jk = j;
      lmkb = lmk + j * nest;
      prod = h[j];
      for(i = 1; i <= k; ++i)
      {
	++jk;
	prod = prod * h[jk] * fac;
      }
      lk = lp + k1;
      b[lmkb] = (t[lk] - t[lp]) / prod;
      ++lp;
    }
  }
}

/*!
* \ingroup	AlgFit
* \brief	Computes the parameters of a Givens transformation.
*
* This function has been derived from the netlib Dierckx function fpgivs().
* The original fortran coments are:
*
*  \verbatim
   Subroutine fpgivs calculates the parameters of a givens transformation
   \endverbatim
*
* \param	piv
* \param	ww
* \param	c
* \param	s
*/
static void			AlgBSplineGivens(
				  double piv,
				  double *ww,
				  double *c,
				  double *s)
{
  double 	r,
  		dd,
		store;

  store = fabs(piv);
  if(store >= *ww) {
    r = *ww / piv;
    dd = store * sqrt(1.0 + r * r);
  }
  if(store < *ww) {
    r = piv / *ww;
    dd = *ww * sqrt(1.0 + r * r);
  }
  *c = *ww / dd;
  *s = piv / dd;
  *ww = dd;
}

/*!
* \return	Dierckx error code.
* \ingroup	AlgFit
* \brief	Determines a smooth periodic spline approximation of degree k
* 		without checks as in AlgBSplinePerFit().
*
* This function has been derived from the netlib Dierckx function fpperi().
* The original fortran coments are included in the code below.
* 
* \param	iopt		Option for how spline is computed.
* \param	x		Independent data.
* \param	y		Dependent data.
* \param	w		Weights for data points.
* \param	m		Number of data points.
* \param	k		Degree of spline.
* \param	s		Smoothing parameter.
* \param	nest		Over estimate of the number of knots.
* \param	tol		Tolerance.
* \param	maxit		Maximum itterations.
* \param	k1
* \param	k2
* \param	n		Destination pointer for the number of knots.
* \param	t		Array for computed knots.
* \param	c		Array for computed spline coefficients.
* \param	fp		Weighted sum of squared residuals.
* \param	fpint
* \param	z
* \param	a1
* \param	a2
* \param	b
* \param	g1
* \param	g2
* \param	q
* \param	nrdata
*/
static int			AlgBSplinePeri(
				  int iopt,
				  double *x,
				  double *y,
				  double *w,
				  int m,
				  int k,
				  double s,
				  int nest,
				  double tol,
				  int maxit,
				  int k1,
				  int k2,
				  int *n,
				  double *t,
				  double *c,
				  double *fp,
				  double *fpint,
				  double *z,
				  double *a1,
				  double *a2,
				  double *b,
				  double *g1,
				  double *g2,
				  double *q,
				  int *nrdata)
{
  int		ier = 0;

  int 		i, j, l,
		ch1, ch3,
		i1, i2, i3, i4, i5, i6, i7, it,
		j1, j2,
		l0, l1, l5,
		ij, ik, jk, 
                kk, kk1, nk1, nk2,
		m1, mm,
  		n7, n8, n10, n11,
		nplus, nrint,
		jper,
		nmin,
		iter,
		nw, nmax,
		npl1,
		a1_dim1, a1_offset,
		a2_dim1, a2_offset,
		b_dim1, b_offset,
		g1_dim1, g1_offset,
		g2_dim1, g2_offset,
		q_dim1, q_offset;
  double 	fp0, fpms, term, pinv,
		fpold, fpart,
		c1, d1, f1, f2, f3, store,
		piv,
		p, p1, p2, p3,
		wi, xi, yi, rn, r1,
		acc, cs, per, sn;
  double	h[6], h1[7], h2[6];

  /* Parameter adjustments for indexing from 1. */
  --w;
  --y;
  --x;
  --z;
  --c;
  --t;
  --fpint;
  --nrdata;
  a2_dim1 = nest;
  a2_offset = 1 + a2_dim1 * 1;
  a2 -= a2_offset;
  q_dim1 = m;
  q_offset = 1 + q_dim1 * 1;
  q -= q_offset;
  g2_dim1 = nest;
  g2_offset = 1 + g2_dim1 * 1;
  g2 -= g2_offset;
  a1_dim1 = nest;
  a1_offset = 1 + a1_dim1 * 1;
  a1 -= a1_offset;
  g1_dim1 = nest;
  g1_offset = 1 + g1_dim1 * 1;
  g1 -= g1_offset;
  b_dim1 = nest;
  b_offset = 1 + b_dim1 * 1;
  b -= b_offset;

  /* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   *  part 1: determination of the number of knots and their position     c
   *  **************************************************************      c
   *  given a set of knots we compute the least-squares periodic spline   c
   *  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c
   *  the initial choice of knots depends on the value of s and iopt.     c
   *    if s=0 we have spline interpolation; in that case the number of   c
   *    knots equals nmax = m+2*k.                                        c
   *    if s > 0 and                                                      c
   *      iopt=0 we first compute the least-squares polynomial of         c
   *      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c
   *      find that s(x) is a constant function.                          c
   *      iopt=1 we start with the set of knots found at the last         c
   *      call of the routine, except for the case that s > fp0; then     c
   *      we compute directly the least-squares periodic polynomial.      c
   * cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  m1 = m - 1;
  kk = k;
  kk1 = k1;
  nmin = k1 << 1;
  /* Determine the length of the period of s(x). */
  per = x[m] - x[1];
  if(iopt < 0)
  {
    goto L50;
  }
  /* Calculation of acc, the absolute tolerance for the root of f(p)=s. */
  acc = tol * s;
  /* Determine nmax, the number of knots for periodic spline interpolation */
  nmax = m + (k << 1);
  if(s > 0.0 || nmax == nmin)
  {
    goto L30;
  }
  /* If s=0, s(x) is an interpolating spline. */
  *n = nmax;
  /* Test whether the required storage space exceeds the available one. */
  if(*n > nest)
  {
      goto L620;
  }
  /* Find the position of the interior knots in case of interpolation. */
L5:
  if(k / 2 << 1 == k)
  {
    goto L20;
  }
  for(i = 2; i <= m1; ++i)
  {
    j = i + k;
    t[j] = x[i];
  }
  if(s > 0.0)
  {
    goto L50;
  }
  kk = k - 1;
  kk1 = k;
  if(kk > 0)
  {
    goto L50;
  }
  t[1] = t[m] - per;
  t[2] = x[1];
  t[m + 1] = x[m];
  t[m + 2] = t[3] + per;
  for(i = 1; i <= m1; ++i)
  {
    c[i] = y[i];
  }
  c[m] = c[1];
  *fp = 0.f;
  fpint[*n] = fp0;
  fpint[*n - 1] = 0.f;
  nrdata[*n] = 0;
  goto L630;
L20:
  for(i = 2; i <= m1; ++i)
  {
    j = i + k;
    t[j] = (x[i] + x[i - 1]) * 0.5;
  }
  goto L50;
  /*  if s > 0 our initial choice depends on the value of iopt.
   *  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
   *  periodic polynomial. (i.e. a constant function).
   *  if iopt=1 and fp0>s we start computing the least-squares periodic
   *  spline according the set of knots found at the last call of the
   *  routine. */
L30:
  if (iopt == 0)
  {
    goto L35;
  }
  if(*n == nmin)
  {
    goto L35;
  }
  fp0 = fpint[*n];
  fpold = fpint[*n - 1];
  nplus = nrdata[*n];
  if(fp0 > s)
  {
    goto L50;
  }
  /* The case that s(x) is a constant function is treated separetely.
   *  find the least-squares constant c1 and compute fp0 at the same time. */
L35:
  d1 = 0.0;
  c1 = 0.0;
  fp0 = 0.0;
  for(it = 1; it <= m1; ++it)
  {
    wi = w[it];
    yi = y[it] * wi;
    AlgBSplineGivens(wi, &d1, &cs, &sn);
    AlgBSplineRota(cs, sn, &yi, &c1);
    fp0 += yi * yi;
  }
  c1 /= d1;
  /* Test whether that constant function is a solution of our problem. */
  fpms = fp0 - s;
  if(fpms < acc || nmax == nmin)
  {
    goto L640;
  }
  fpold = fp0;
  /* Test whether the required storage space exceeds the available one. */
  if(nmin >= nest)
  {
    goto L620;
  }
  /* Start computing the least-squares periodic spline with one interior
   * knot. */
  nplus = 1;
  *n = nmin + 1;
  mm = (m + 1) / 2;
  t[k2] = x[mm];
  nrdata[1] = mm - 2;
  nrdata[2] = m1 - mm;
  /* Main loop for the different sets of knots. m is a save upper bound for
   * the number of trials. */
L50:
  for(iter = 1; iter <= m; ++iter)
  {
    /*  Find nrint, the number of knot intervals. */
    nrint = *n - nmin + 1;
    /*  find the position of the additional knots which are needed for
     *  the b-spline representation of s(x). if we take
     *      t(k+1) = x(1), t(n-k) = x(m)
     *      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
     *      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
     *  then s(x) is a periodic spline with period per if the b-spline
     *  coefficients satisfy the following conditions
     *      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1. */
      t[k1] = x[1];
      nk1 = *n - k1;
      nk2 = nk1 + 1;
      t[nk2] = x[m];
      for(j = 1; j <= k; ++j)
      {
	i1 = nk2 + j;
	i2 = nk2 - j;
	j1 = k1 + j;
	j2 = k1 - j;
	t[i1] = t[j1] + per;
	t[j2] = t[i2] - per;
      }
      /* compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
       * periodic spline sinf(x). the observation matrix a is built up row
       * by row while taking into account condition (**) and is reduced to
       * triangular form by givens transformations .
       * at the same time fp=f(p=inf) is computed.
       * the n7 x n7 triangularised upper matrix a has the form
       *           ! a1 '    !
       *       a = !    ' a2 !
       *           ! 0  '    !
       * with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
       * matrix of bandwith k+1 ( n10 = n7-k).
       * initialization. */
      for(i = 1; i <= nk1; ++i)
      {
	z[i] = 0.f;
	for(j = 1; j <= kk1; ++j)
	{
	  a1[i + j * a1_dim1] = 0.f;
	  }
      }
      n7 = nk1 - k;
      n10 = n7 - kk;
      jper = 0;
      *fp = 0.f;
      l = k1;
      for(it = 1; it <= m1; ++it)
      {
        /* Fetch the current data point x(it),y(it) */
	xi = x[it];
	wi = w[it];
	yi = y[it] * wi;
	/* Search for knot interval t(l) <= xi < t(l+1). */
L80:
	if(xi < t[l + 1])
	{
	  goto L85;
	}
	++l;
	goto L80;
        /* Evaluate the (k+1) non-zero b-splines at xi and store them in q. */
L85:
	  AlgBSplineBspl(&t[1], k, xi, l, h);
	  for(i = 1; i <= k1; ++i)
	  {
	    q[it + i * q_dim1] = h[i - 1];
	    h[i - 1] *= wi;
	    /* L90: */
	  }
	  l5 = l - k1;
	  /* Test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all
	   * zero at xi */
	  if(l5 < n10)
	  {
	    goto L285;
	  }
	  if(jper != 0)
	  {
	    goto L160;
	  }
          /* Initialize the matrix a2. */
	  for(i = 1; i <= n7; ++i)
	  {
	    for(j = 1; j <= kk; ++j)
	    {
	      a2[i + j * a2_dim1] = 0.0;
	    }
	  }
	  jk = n10 + 1;
	  for(i = 1; i <= kk; ++i)
	  {
	    ik = jk;
	    for(j = 1; j <= kk1; ++j)
	    {
	      if(ik <= 0)
	      {
		goto L105;
	      }
	      a2[ik + i * a2_dim1] = a1[ik + j * a1_dim1];
	      --ik;
	    }
L105:
	    ++jk;
	  }
	  jper = 1;
	  /* If one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
	   * we take account of condition (**) for setting up the new row
	   * of the observation matrix a. this row is stored in the arrays h1
	   * (the part with respect to a1) and h2 (the part with respect
	   * to a2). */
L160:
	  for(i = 1; i <= kk; ++i)
	  {
	    h1[i - 1] = 0.0;
	    h2[i - 1] = 0.0;
	  }
	  h1[kk1 - 1] = 0.0;
	  j = l5 - n10;
	  for(i = 1; i <= kk1; ++i)
	  {
	    ++j;
	    l0 = j;
L180:
	    l1 = l0 - kk;
	    if(l1 <= 0)
	    {
	      goto L200;
	    }
	    if(l1 <= n10)
	    {
	      goto L190;
	    }
	    l0 = l1 - n10;
	    goto L180;
L190:
	    h1[l1 - 1] = h[i - 1];
	    goto L210;
L200:
	    h2[l0 - 1] += h[i - 1];
L210:
	      ;
	  }
	  /* Rotate the new row of the observation matrix into triangle
	   * by givens transformations. */
	  if(n10 <= 0)
	  {
	    goto L250;
	  }
          /* Rotation with the rows 1,2,...n10 of matrix a. */
	  for(j = 1; j <= n10; ++j)
	  {
	    piv = h1[0];
	    if(piv != 0.0)
	    {
	      goto L214;
	    }
	    for(i = 1; i <= kk; ++i)
	    {
	      h1[i - 1] = h1[i];
	    }
	    h1[kk1 - 1] = 0.f;
	    goto L240;
	    /* Calculate the parameters of the givens transformation. */
L214:
	    AlgBSplineGivens(piv, &a1[j + a1_dim1], &cs, &sn);
            /* Transformation to the right hand side. */
	    AlgBSplineRota(cs, sn, &yi, &z[j]);
	    /* Transformations to the left hand side with respect to a2. */
	    for(i = 1; i <= kk; ++i)
	    {
	      AlgBSplineRota(cs, sn, &h2[i - 1], &a2[j + i * a2_dim1]);
	    }
	    if(j == n10)
	    {
	      goto L250;
	    }
            /* Computing MIN */
	    i4 = n10 - j;
	    i2 = ALG_MIN(i4, kk);
	    /* Transformations to the left hand side with respect to a1. */
	    for(i = 1; i <= i2; ++i)
	    {
	      i1 = i + 1;
	      AlgBSplineRota(cs, sn, &h1[i1 - 1], &a1[j + i1 * a1_dim1]);
	      h1[i - 1] = h1[i1 - 1];
	    }
	    h1[i1 - 1] = 0.f;
L240:
	    ;
	  }
	  /* Rotation with the rows n10+1,...n7 of matrix a. */
L250:
	  for(j = 1; j <= kk; ++j)
	  {
	    ij = n10 + j;
	    if(ij <= 0)
	    {
	      goto L270;
	    }
	    piv = h2[j - 1];
	    if(piv == 0.0)
	    {
	      goto L270;
	    }
	    /* Calculate the parameters of the givens transformation. */
	    AlgBSplineGivens(piv, &a2[ij + j * a2_dim1], &cs, &sn);
	    /* Transformations to right hand side. */
	    AlgBSplineRota(cs, sn, &yi, &z[ij]);
	    if(j == kk)
	    {
	      goto L280;
	    }
	    j1 = j + 1;
	    /* Transformations to left hand side. */
	    for(i = j1; i <= kk; ++i)
	    {
	      AlgBSplineRota(cs, sn, &h2[i - 1], &a2[ij + i * a2_dim1]);
	    }
L270:
	    ;
	  }
	  /*  Add contribution of this row to the sum of squares of residual
	   *  right hand sides. */
L280:
	  *fp += yi * yi;
	  goto L290;
	  /* Rotation of the new row of the observation matrix into triangle
	   * in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero 
	   *  at xi. */
L285:
	  j = l5;
	  for(i = 1; i <= kk1; ++i)
	  {
	    ++j;
	    piv = h[i - 1];
	    if(piv == 0.0)
	    {
	      goto L140;
	    }
	    /*  calculate the parameters of the givens transformation. */
	    AlgBSplineGivens(piv, &a1[j + a1_dim1], &cs, &sn);
	    /*  transformations to right hand side. */
	    AlgBSplineRota(cs, sn, &yi, &z[j]);
	    if(i == kk1)
	    {
	      goto L150;
	    }
	    i2 = 1;
	    i3 = i + 1;
	    /*  transformations to left hand side. */
	    for(i1 = i3; i1 <= kk1; ++i1)
	    {
	      ++i2;
	      AlgBSplineRota(cs, sn, &h[i1 - 1], &a1[j + i2 * a1_dim1]);
	    }
L140:
	    ;
	  }
	  /* Add contribution of this row to the sum of squares of residual
	   * right hand sides. */
L150:
	  *fp += yi * yi;
L290:
	  ;
      }
      fpint[*n] = fp0;
      fpint[*n - 1] = fpold;
      nrdata[*n] = nplus;
      /* Backward substitution to obtain the b-spline coefficients
       * c(j),j=1,.n */
      AlgBSplineBacp(&a1[a1_offset], &a2[a2_offset],
          &z[1], n7, kk, &c[1], nest);
      /* Calculate from condition (**) the coefficients c(j+n7),j=1,2,...k. */
      for(i = 1; i <= k; ++i)
      {
	j = i + n7;
	c[j] = c[i];
      }
      if(iopt < 0)
      {
	goto L660;
      }
      /*  test whether the approximation sinf(x) is an acceptable solution. */
      fpms = *fp - s;
      if(fabs(fpms) < acc)
      {
	goto L660;
      }
      /*  if f(p=inf) < s accept the choice of knots. */
      if(fpms < 0.0)
      {
	goto L350;
      }
      /* If n=nmax, sinf(x) is an interpolating spline. */
      if(*n == nmax)
      {
	goto L630;
      }
      /* Increase the number of knots. */
      /* If n=nest we cannot increase the number of knots because of the
       * storage capacity limitation. */
      if(*n == nest)
      {
	  goto L620;
      }
      /* Determine the number of knots nplus we are going to add. */
      npl1 = nplus << 1;
      rn = (double) nplus;
      if(fpold - *fp > acc)
      {
	npl1 = rn * fpms / (fpold - *fp);
      }
      i4 = npl1, i7 = nplus / 2, i4 = ALG_MAX(i4, i7);
      i5 = nplus << 1, i6 = ALG_MAX(i4, 1);
      nplus = ALG_MIN(i5, i6);
      fpold = *fp;
      /* Compute the sum(wi*(yi-s(xi))**2) for each knot interval
       * t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint. */
      fpart = 0.0;
      i = 1;
      l = k1;
      for(it = 1; it <= m1; ++it)
      {
	if (x[it] < t[l]) {
	  goto L300;
	}
	nw = 1;
	++l;
L300:
	term = 0.0;
	l0 = l - k2;
	for(j = 1; j <= k1; ++j)
	{
	  ++l0;
	  term += c[l0] * q[it + j * q_dim1];
	}
	r1 = w[it] * (term - y[it]);
	term = r1 * r1;
	fpart += term;
	if(nw == 0)
	{
	  goto L320;
	}
	if(l > k2)
	{
	  goto L315;
	}
	fpint[nrint] = term;
	nw = 0;
	goto L320;
L315:
	store = term * 0.5;
	fpint[i] = fpart - store;
	++i;
	fpart = store;
	nw = 0;
L320:
	;
      }
      fpint[nrint] += fpart;
      for(l = 1; l <= nplus; ++l)
      {
	/* Add a new knot */
	AlgBSplineKnot(&x[1], m, &t[1], n, &fpint[1],
	    &nrdata[1], &nrint, 1);
	/* If n=nmax we locate the knots as for interpolation. */
	if(*n == nmax)
	{
	  goto L5;
	}
	/* Test whether we cannot further increase the number of knots. */
	if(*n == nest)
	{
	  goto L340;
	}
      }
      /* Restart the computations with the new set of knots. */
L340:
      ;
  }
  /* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   *  part 2: determination of the smoothing periodic spline sp(x).       c
   *  *************************************************************       c
   *  we have determined the number of knots and their position.          c
   *  we now compute the b-spline coefficients of the smoothing spline    c
   *  sp(x). the observation matrix a is extended by the rows of matrix   c
   *  b expressing that the kth derivative discontinuities of sp(x) at    c
   *  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
   *  ponding weights of these additional rows are set to 1/sqrt(p).      c
   *  iteratively we then have to determine the value of p such that      c
   *  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c
   *  the least-squares constant function corresponds to p=0, and that    c
   *  the least-squares periodic spline corresponds to p=infinity. the    c
   *  iteration process which is proposed here, makes use of rational     c
   *  interpolation. since f(p) is a convex and strictly decreasing       c
   *  function of p, it can be approximated by a rational function        c
   *  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
   *  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
   *  to calculate the new value of p such that r(p)=s. convergence is    c
   *  guaranteed by taking f1>0 and f3<0.                                 c
   * cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*  Evaluate the discontinuity jump of the kth derivative of the
   *  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b. */
L350:
  AlgBSplineDisc(&t[1], *n, k2, &b[b_offset], nest);
  /*  initial value for p. */
  p1 = 0.0;
  f1 = fp0 - s;
  p3 = -1.0;
  f3 = fpms;
  n11 = n10 - 1;
  n8 = n7 - 1;
  p = 0.0;
  l = n7;
  for(i = 1; i <= k; ++i)
  {
    j = k + 1 - i;
    p += a2[l + j * a2_dim1];
    --l;
    if(l == 0)
    {
      goto L356;
    }
  }
  for(i = 1; i <= n10; ++i)
  {
    p += a1[i + a1_dim1];
  }
L356:
  rn = (double )n7;
  p = rn / p;
  ch1 = 0;
  ch3 = 0;
  /* Iteration process to find the root of f(p) = s. */
  for(iter = 1; iter <= maxit; ++iter)
  {
    /*  Form the matrix g  as the matrix a extended by the rows of matrix b.
     *  the rows of matrix b with weight 1/p are rotated into
     *  the triangularised observation matrix a.
     *  after triangularisation our n7 x n7 matrix g takes the form
     *            ! g1 '    !
     *        g = !    ' g2 !
     *            ! 0  '    !
     *  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
     *  matrix of bandwidth k+2. ( n11 = n7-k-1) */
    pinv = 1.0 / p;
    /* Store matrix a into g */
    for(i = 1; i <= n7; ++i)
    {
      c[i] = z[i];
      g1[i + k1 * g1_dim1] = a1[i + k1 * a1_dim1];
      g1[i + k2 * g1_dim1] = 0.f;
      g2[i + g2_dim1] = 0.f;
      for(j = 1; j <= k; ++j)
      {
	g1[i + j * g1_dim1] = a1[i + j * a1_dim1];
	g2[i + (j + 1) * g2_dim1] = a2[i + j * a2_dim1];
      }
    }
    l = n10;
    for(j = 1; j <= k1; ++j)
    {
      if(l <= 0)
      {
	goto L375;
      }
      g2[l + g2_dim1] = a1[l + j * a1_dim1];
      --l;
    }
L375:
    for(it = 1; it <= n8; ++it)
    {
      /*  fetch a new row of matrix b and store it in the arrays h1 (the part
       *  with respect to g1) and h2 (the part with respect to g2). */
      yi = 0.0;
      for(i = 1; i <= k1; ++i)
      {
	h1[i - 1] = 0.0;
	h2[i - 1] = 0.0;
      }
      h1[k2 - 1] = 0.0;
      if (it > n11)
      {
	goto L420;
      }
      l = it;
      l0 = it;
      for(j = 1; j <= k2; ++j)
      {
	if(l0 == n10)
	{
	  goto L400;
	}
	h1[j - 1] = b[it + j * b_dim1] * pinv;
	++l0;
      }
      goto L470;
L400:
      l0 = 1;
      for(l1 = j; l1 <= k2; ++l1)
      {
	h2[l0 - 1] = b[it + l1 * b_dim1] * pinv;
	++l0;
      }
      goto L470;
L420:
      l = 1;
      i = it - n10;
      for(j = 1; j <= k2; ++j)
      {
	++i;
	l0 = i;
L430:
	l1 = l0 - k1;
	if(l1 <= 0)
	{
	  goto L450;
	}
	if(l1 <= n11)
	{
	  goto L440;
	}
	l0 = l1 - n11;
	goto L430;
L440:
	h1[l1 - 1] = b[it + j * b_dim1] * pinv;
	goto L460;
L450:
	h2[l0 - 1] += b[it + j * b_dim1] * pinv;
L460:
	;
      }
      if(n11 <= 0)
      {
	goto L510;
      }
      /* Rotate this row into triangle by givens transformations without
       * square roots. */
      /* Rotation with the rows l,l+1,...n11. */
L470:
      for(j = l; j <= n11; ++j)
      {
	piv = h1[0];
	/*  Calculate the parameters of the givens transformation. */
	AlgBSplineGivens(piv, &g1[j + g1_dim1], &cs, &sn);
	/*  Transformation to right hand side. */
	AlgBSplineRota(cs, sn, &yi, &c[j]);
	/*  Transformation to the left hand side with respect to g2. */
	for(i = 1; i <= k1; ++i)
	{
	  AlgBSplineRota(cs, sn, &h2[i - 1], &g2[j + i * g2_dim1]);
	}
	if(j == n11)
	{
	  goto L510;
	}
	i = n11 - j;
	i2 = ALG_MIN(i, k1);
	/* Transformation to the left hand side with respect to g1. */
	for(i = 1; i <= i2; ++i)
	{
	  i1 = i + 1;
	  AlgBSplineRota(cs, sn, &h1[i1 - 1], &g1[j + i1 * g1_dim1]);
	  h1[i - 1] = h1[i1 - 1];
	}
	h1[i1 - 1] = 0.0;
      }
      /* Rotation with the rows n11+1,...n7 */
L510:
      for(j = 1; j <= k1; ++j)
      {
	ij = n11 + j;
	if(ij <= 0)
	{
	  goto L530;
	}
	piv = h2[j - 1];
	/*  calculate the parameters of the givens transformation */
	AlgBSplineGivens(piv, &g2[ij + j * g2_dim1], &cs, &sn);
	/*  transformation to the right hand side. */
	AlgBSplineRota(cs, sn, &yi, &c[ij]);
	if(j == k1)
	{
	  goto L540;
	}
	j1 = j + 1;
	/* Transformation to the left hand side. */
	for(i = j1; i <= k1; ++i)
	{
	  AlgBSplineRota(cs, sn, &h2[i - 1], &g2[ij + i * g2_dim1]);
	}
L530:
	;
      }
L540:
      ;
    }
    /*  Backward substitution to obtain the b-spline coefficients 
     *  c(j),j=1,2,...n7 of sp(x). */
    AlgBSplineBacp(&g1[g1_offset], &g2[g2_offset],
        &c[1], n7, k1, &c[1], nest);
    /* Calculate from condition (**) the b-spline coefficients c(n7+j),j=1,. */
    for(i = 1; i <= k; ++i)
    {
      j = i + n7;
      c[j] = c[i];
    }
    /* Computation of f(p). */
    *fp = 0.0;
    l = k1;
    for(it = 1; it <= m1; ++it)
    {
      if(x[it] < t[l])
      {
	goto L550;
      }
      ++l;
L550:
      l0 = l - k2;
      term = 0.f;
      for(j = 1; j <= k1; ++j)
      {
	++l0;
	term += c[l0] * q[it + j * q_dim1];
      }
      r1 = w[it] * (term - y[it]);
      *fp += r1 * r1;
    }
    /* Test whether the approximation sp(x) is an acceptable solution. */
    fpms = *fp - s;
    if(fabs(fpms) < acc)
    {
      goto L660;
    }
    /*  test whether the maximal number of iterations is reached. */
    if(iter == maxit)
    {
      goto L600;
    }
    /* Carry out one more step of the iteration process. */
    p2 = p;
    f2 = fpms;
    if(ch3 != 0)
    {
      goto L580;
    }
    if(f2 - f3 > acc)
    {
      goto L575;
    }
    /* Our initial choice of p is too large. */
    p3 = p2;
    f3 = f2;
    p *= 0.04;
    if(p <= p1)
    {
      p = p1 * 0.9 + p2 * 0.1;
    }
    goto L595;
L575:
    if(f2 < 0.0)
    {
      ch3 = 1;
    }
L580:
    if(ch1 != 0)
    {
      goto L590;
    }
    if(f1 - f2 > acc)
    {
      goto L585;
    }
    /*  Our initial choice of p is too small */
    p1 = p2;
    f1 = f2;
    p /= 0.04;
    if(p3 < 0.0)
    {
      goto L595;
    }
    if(p >= p3)
    {
      p = p2 * 0.1 + p3 * 0.9;
    }
    goto L595;
L585:
    if(f2 > 0.0) 
    {
      ch1 = 1;
    }
    /* Test whether the iteration process proceeds as theoretically
     * expected. */
L590:
    if(f2 >= f1 || f2 <= f3)
    {
      goto L610;
    }
    /* Find the new value for p. */
    p = AlgBSplineRati(&p1, &f1, p2, f2, &p3, &f3);
L595:
    ;
  }
  /*  error codes and messages. */
L600:
  ier = 3;
  goto L660;
L610:
  ier = 2;
  goto L660;
L620:
  ier = 1;
  goto L660;
L630:
  ier = -1;
  goto L660;
L640:
  ier = -2;
  /* The least-squares constant function c1 is a solution of our problem.
   * a constant function is a spline of degree k with all b-spline
   * coefficients equal to that constant c1. */
  for(i = 1; i <= k1; ++i)
  {
    rn = (double) (k1 - i);
    t[i] = x[1] - rn * per;
    c[i] = c1;
    j = i + k1;
    rn = (double) (i - 1);
    t[j] = x[m] + rn * per;
  }
  *n = nmin;
  *fp = fp0;
  fpint[*n] = fp0;
  fpint[*n - 1] = 0.0;
  nrdata[*n] = 0;
L660:
  return(ier);
}

/*!
* \return	Interpolated value.
* \ingroup	AlgFit
* \brief	Evaluates interpolating function.
*
* This function has been derived from the netlib Dierckx function fprati().
* The original fortran coments are:
*
*  \verbatim
   Given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
   gives the value of p such that the rational interpolating function
   of the form r(p) = (u*p+v)/(p+w) equals zero at p.
   \endverbatim
*
* \param	p1			First given point (value may be
* 					modified).
* \param	f1			First given point (value may be
* 					modified).
* \param	p2			Second given point.
* \param	f2			Second given point.
* \param	p3			Third given point (value may be
* 					modified).
* \param	f3			Third given point (value may be
* 					modified).
*/
static double			AlgBSplineRati(
				  double *p1,
				  double *f1,
				  double p2,
				  double f2,
				  double *p3,
				  double *f3)
{
  double 	p;

  if(*p3 > 0.0)
  {
    double	h1,
    		h2,
		h3;

    /*  Value of p in case p3 ^= infinity. */
    h1 = *f1 * (f2 - *f3);
    h2 = f2 * (*f3 - *f1);
    h3 = *f3 * (*f1 - f2);
    p = -(*p1 * p2 * h3 + p2 * *p3 * h1 + *p3 * *p1 * h2) /
	(*p1 * h1 + p2 * h2 + *p3 * h3);
  }
  else
  {
    /*  Value of p in case p3 = infinity. */
    p = (*p1 * (*f1 - *f3) * f2 - p2 * (f2 - *f3) * *f1) /
	((*f1 - f2) * * f3);
  }
  /*  Adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0. */
  if(f2 < 0.0)
  {
    *p3 = p2;
    *f3 = f2;
  }
  else
  {
    *p1 = p2;
    *f1 = f2;
  }
  return(p);
}

/*!
* \ingroup	AlgFit
* \brief	Applies a Givens rotation.
* 
* This function has been derived from the netlib Dierckx function fprota().
* The original fortran coments are:
*
*  \verbatim
   Subroutine fprota applies a givens rotation to a and b.
   \endverbatim
*
* \param	cs		Cosine parameter.
* \param	sn		Sine parameter.
* \param	a		First parameter for rotation.
* \param	b		Second parameter for rotation.
*/
static void			AlgBSplineRota(
				  double cs,
				  double sn,
				  double *a,
				  double *b)
{
  double 	stor1,
	      stor2;

  stor1 = *a;
  stor2 = *b;
  *b = cs * stor2 + sn * stor1;
  *a = cs * stor1 - sn * stor2;
}

/*!
* \return	Dierckx error code.
* \ingroup	AlgFit
* \brief	Computes a smooth B-spline approximation through the
* 		given data. See AlgBSplineFit().
*
* This function has been derived from the netlib Dierckx function fpcurf().
* The original fortran coments are included in the code below.
* \param	iopt
* \param	x
* \param	y
* \param	w
* \param	m
* \param	xb
* \param	xe
* \param	k
* \param	s
* \param	nest
* \param	tol
* \param	maxit
* \param	k1
* \param	k2
* \param	n
* \param	t
* \param	c
* \param	fp
* \param	fpint
* \param	z
* \param	a
* \param	b
* \param	g
* \param	q
* \param	nrdata
*/
static int 			AlgBSplineCurf(
				  int iopt,
				  double *x,
				  double *y,
				  double *w,
				  int m,
				  double xb,
				  double xe,
				  int k,
				  double s,
				  int nest,
				  double tol,
				  int maxit,
				  int k1,
				  int k2,
				  int *n,
				  double *t,
				  double *c,
				  double *fp,
				  double *fpint,
				  double *z,
				  double *a,
				  double *b,
				  double *g,
				  double *q,
				  int *nrdata)
{
  int 		i, i1, i2, i3, i4, i5, i6, i7,
      		ich1, ich3, it,
      		j,
      		k3,
      		l, l0,
		nplus, nrint, n8,
		mk1, nk1,
		nmin, iter, nmax,
		npl1,
		nw,
  		a_dim1, a_offset,
		b_dim1, b_offset,
		g_dim1, g_offset,
		q_dim1, q_offset,
		ier = 0;
  double 	acc, cs, sn,
  		fp0, fpms, term, pinv,
		p, p1, p2, p3,
		fpold, fpart,
		f1, f2, f3,
		piv, rn, r1,
		wi, xi, yi, store;
  double 	h[7];

  /* Parameter adjustments for indexing from 1. */
  --c;
  --t;
  --w;
  --y;
  --x;
  --z;
  --fpint;
  --nrdata;
  q_dim1 = m;
  q_offset = 1 + q_dim1 * 1;
  q -= q_offset;
  a_dim1 = nest;
  a_offset = 1 + a_dim1 * 1;
  a -= a_offset;
  g_dim1 = nest;
  g_offset = 1 + g_dim1 * 1;
  g -= g_offset;
  b_dim1 = nest;
  b_offset = 1 + b_dim1 * 1;
  b -= b_offset;

  /* Function Body */
  /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*  part 1: determination of the number of knots and their position     c */
  /*  **************************************************************      c */
  /*  given a set of knots we compute the least-squares spline sinf(x),   c */
  /*  and the corresponding sum of squared residuals fp=f(p=inf).         c */
  /*  if iopt=-1 sinf(x) is the requested approximation.                  c */
  /*  if iopt=0 or iopt=1 we check whether we can accept the knots:       c */
  /*    if fp <=s we will continue with the current set of knots.         c */
  /*    if fp > s we will increase the number of knots and compute the    c */
  /*       corresponding least-squares spline until finally fp<=s.        c */
  /*    the initial choice of knots depends on the value of s and iopt.   c */
  /*    if s=0 we have spline interpolation; in that case the number of   c */
  /*    knots equals nmax = m+k+1.                                        c */
  /*    if s > 0 and                                                      c */
  /*      iopt=0 we first compute the least-squares polynomial of         c */
  /*      degree k; n = nmin = 2*k+2                                      c */
  /*      iopt=1 we start with the set of knots found at the last         c */
  /*      call of the routine, except for the case that s > fp0; then     c */
  /*      we compute directly the least-squares polynomial of degree k.   c */
  /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*  determine nmin, the number of knots for polynomial approximation. */
  nmin = k1 * 2;
  if(iopt < 0)
  {
    goto L60;
  }
  /*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
  acc = tol * s;
  /*  determine nmax, the number of knots for spline interpolation. */
  nmax = m + k1;
  if(s > 0.0)
  {
    goto L45;
  }
  /*  if s=0, s(x) is an interpolating spline. */
  /*  test whether the required storage space exceeds the available one. */
  *n = nmax;
  if(nmax > nest)
  {
    goto L420;
  }
  /*  find the position of the interior knots in case of interpolation. */
L10:
  mk1 = m - k1;
  if(mk1 == 0)
  {
    goto L60;
  }
  k3 = k / 2;
  i = k2;
  j = k3 + 2;
  if(k3 * 2 == k)
  {
    goto L30;
  }
  for(l = 0; l < mk1; ++l)
  {
    t[i] = x[j];
    ++i;
    ++j;
  }
  goto L60;
L30:
  for(l = 0; l < mk1; ++l)
  {
    t[i] = (x[j] + x[j - 1]) * 0.5;
    ++i;
    ++j;
  }
  goto L60;
  /*  if s>0 our initial choice of knots depends on the value of iopt. */
  /*  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares */
  /*  polynomial of degree k which is a spline without interior knots. */
  /*  if iopt=1 and fp0>s we start computing the least squares spline */
  /*  according to the set of knots found at the last call of the routine. */
L45:
  if(iopt == 0)
  {
    goto L50;
  }
  if(*n == nmin)
  {
    goto L50;
  }
  fp0 = fpint[*n];
  fpold = fpint[*n - 1];
  nplus = nrdata[*n];
  if(fp0 > s)
  {
    goto L60;
  }
L50:
  *n = nmin;
  fpold = 0.f;
  nplus = 0;
  nrdata[1] = m - 2;
  /*  main loop for the different sets of knots. m is a save upper bound */
  /*  for the number of trials. */
L60:
  for(iter = 1; iter <= m; ++iter)
  {
    if(*n == nmin)
    {
      ier = -2;
    }
    /*  find nrint, tne number of knot intervals. */
    nrint = *n - nmin + 1;
    /*  find the position of the additional knots which are needed for */
    /*  the b-spline representation of s(x). */
    nk1 = *n - k1;
    i = *n;
    for(j = 1; j <= k1; ++j)
    {
      t[j] = xb;
      t[i] = xe;
      --i;
    }
    /*  compute the b-spline coefficients of the least-squares spline */
    /*  sinf(x). the observation matrix a is built up row by row and */
    /*  reduced to upper triangular form by givens transformations. */
    /*  at the same time fp=f(p=inf) is computed. */
    *fp = 0.0;
    /*  initialize the observation matrix a. */
    for(i = 1; i <= nk1; ++i)
    {
      z[i] = 0.0;
      for (j = 1; j <= k1; ++j)
      {
	a[i + j * a_dim1] = 0.0;
      }
    }
    l = k1;
    for(it = 1; it <= m; ++it)
    {
      /*  fetch the current data point x(it),y(it). */
      xi = x[it];
      wi = w[it];
      yi = y[it] * wi;
      /*  search for knot interval t(l) <= xi < t(l+1). */
L85:
      if((xi < t[l + 1]) || (l == nk1))
      {
	goto L90;
      }
      ++l;
      goto L85;
      /*  evaluate the (k+1) non-zero b-splines at xi and store them in q. */
L90:
      AlgBSplineBspl(&t[1], k, xi, l, h);
      for(i = 1; i <= k1; ++i)
      {
	q[it + i * q_dim1] = h[i - 1];
	h[i - 1] *= wi;
      }
      /*  rotate the new row of the observation matrix into triangle. */
      j = l - k1;
      for(i = 1; i <= k1; ++i)
      {
	++j;
	piv = h[i - 1];
	if(piv == 0.0)
	{
	  goto L110;
	}
	/*  calculate the parameters of the givens transformation. */
	AlgBSplineGivens(piv, &a[j + a_dim1], &cs, &sn);
	/*  transformations to right hand side. */
	AlgBSplineRota(cs, sn, &yi, &z[j]);
	if(i == k1)
	{
	  goto L120;
	}
	i2 = 1;
	i3 = i + 1;
	for(i1 = i3; i1 <= k1; ++i1)
	{
	  ++i2;
	  /*  transformations to left hand side. */
	  AlgBSplineRota(cs, sn, &h[i1 - 1], &a[j + i2 * a_dim1]);
	}
L110:
	;
      }
      /*  add contribution of this row to the sum of squares of residual */
      /*  right hand sides. */
L120:
      r1 = yi;
      *fp += r1 * r1;
    }
    if(ier == -2)
    {
      fp0 = *fp;
    }
    fpint[*n] = fp0;
    fpint[*n - 1] = fpold;
    nrdata[*n] = nplus;
    /*  backward substitution to obtain the b-spline coefficients. */
    AlgBSplineBack(&a[a_offset], &z[1], nk1, k1, &c[1], nest);
    /*  test whether the approximation sinf(x) is an acceptable solution. */
    if (iopt < 0) {
      goto L440;
    }
    fpms = *fp - s;
    if(fabs(fpms) < acc)
    {
      goto L440;
    }
    /*  if f(p=inf) < s accept the choice of knots. */
    if(fpms < 0.0)
    {
      goto L250;
    }
    /*  if n = nmax, sinf(x) is an interpolating spline. */
    if(*n == nmax)
    {
      goto L430;
    }
    /*  increase the number of knots. */
    /*  if n=nest we cannot increase the number of knots because of */
    /*  the storage capacity limitation. */
    if(*n == nest)
    {
      goto L420;
    }
    /*  determine the number of knots nplus we are going to add. */
    if(ier == 0)
    {
      goto L140;
    }
    nplus = 1;
    ier = 0;
    goto L150;
L140:
    npl1 = nplus * 2;
    rn = (double )nplus;
    if(fpold - *fp > acc)
    {
      npl1 = rn * fpms / (fpold - *fp);
    }
    i7 = nplus / 2;
    i6 = ALG_MAX(npl1,i7);
    i5 = nplus * 2;
    i4 = ALG_MAX(i6,1);
    nplus = ALG_MIN(i5,i4);
L150:
    fpold = *fp;
    /*  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval */
    /*  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint. */
    fpart = 0.0;
    i = 1;
    l = k2;
    nw = 0;
    for(it = 1; it <= m; ++it)
    {
      if((x[it] < t[l]) || (l > nk1))
      {
	goto L160;
      }
      nw = 1;
      ++l;
L160:
      term = 0.0;
      l0 = l - k2;
      for(j = 1; j <= k1; ++j)
      {
	++l0;
	term += c[l0] * q[it + j * q_dim1];
      }
      r1 = w[it] * (term - y[it]);
      term = r1 * r1;
      fpart += term;
      if(nw == 0)
      {
	goto L180;
      }
      store = term * 0.5;
      fpint[i] = fpart - store;
      ++i;
      fpart = store;
      nw = 0;
L180:
      ;
    }
    fpint[nrint] = fpart;
    for(l = 1; l <= nplus; ++l)
    {
      /*  add a new knot. */
      AlgBSplineKnot(&x[1], m, &t[1], n, &fpint[1], &nrdata[1], &nrint, 1);
      /*  if n=nmax we locate the knots as for interpolation. */
      if(*n == nmax)
      {
	goto L10;
      }
      /*  test whether we cannot further increase the number of knots. */
      if(*n == nest)
      {
	goto L200;
      }
    }
    /*  restart the computations with the new set of knots. */
L200:
    ;
  }
  /*  test whether the least-squares kth degree polynomial is a solution */
  /*  of our approximation problem. */
L250:
  if (ier == -2)
  {
    goto L440;
  }
  /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*  part 2: determination of the smoothing spline sp(x).                c */
  /*  ***************************************************                 c */
  /*  we have determined the number of knots and their position.          c */
  /*  we now compute the b-spline coefficients of the smoothing spline    c */
  /*  sp(x). the observation matrix a is extended by the rows of matrix   c */
  /*  b expressing that the kth derivative discontinuities of sp(x) at    c */
  /*  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c */
  /*  ponding weights of these additional rows are set to 1/p.            c */
  /*  iteratively we then have to determine the value of p such that      c */
  /*  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c */
  /*  the least-squares kth degree polynomial corresponds to p=0, and     c */
  /*  that the least-squares spline corresponds to p=infinity. the        c */
  /*  iteration process which is proposed here, makes use of rational     c */
  /*  interpolation. since f(p) is a convex and strictly decreasing       c */
  /*  function of p, it can be approximated by a rational function        c */
  /*  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c */
  /*  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c */
  /*  to calculate the new value of p such that r(p)=s. convergence is    c */
  /*  guaranteed by taking f1>0 and f3<0.                                 c */
  /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  /*  evaluate the discontinuity jump of the kth derivative of the */
  /*  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b. */
  AlgBSplineDisc(&t[1], *n, k2, &b[b_offset], nest);
  /*  initial value for p. */
  p1 = 0.0;
  f1 = fp0 - s;
  p3 = -1.0;
  f3 = fpms;
  p = 0.0;
  for(i = 1; i <= nk1; ++i)
  {
    p += a[i + a_dim1];
  }
  rn = (double )nk1;
  p = rn / p;
  ich1 = 0;
  ich3 = 0;
  n8 = *n - nmin;
  /*  iteration process to find the root of f(p) = s. */
  for(iter = 1; iter <= maxit; ++iter)
  {
    /*  the rows of matrix b with weight 1/p are rotated into the */
    /*  triangularised observation matrix a which is stored in g. */
    pinv = 1.0 / p;
    for(i = 1; i <= nk1; ++i)
    {
      c[i] = z[i];
      g[i + k2 * g_dim1] = 0.0;
      for (j = 1; j <= k1; ++j)
      {
	g[i + j * g_dim1] = a[i + j * a_dim1];
      }
    }
    for(it = 1; it <= n8; ++it)
    {
      /*  the row of matrix b is rotated into triangle by givens
       *  transformation */
      for(i = 1; i <= k2; ++i)
      {
	h[i - 1] = b[it + i * b_dim1] * pinv;
      }
      yi = 0.0;
      for(j = it; j <= nk1; ++j)
      {
	piv = h[0];
	/*  calculate the parameters of the givens transformation. */
	AlgBSplineGivens(piv, &g[j + g_dim1], &cs, &sn);
	/*  transformations to right hand side. */
	AlgBSplineRota(cs, sn, &yi, &c[j]);
	if(j == nk1)
	{
	  goto L300;
	}
	i2 = k1;
	if(j > n8)
	{
	  i2 = nk1 - j;
	}
	for(i = 1; i <= i2; ++i)
	{
	  /*  transformations to left hand side. */
	  i1 = i + 1;
	  AlgBSplineRota(cs, sn, &h[i1 - 1], &g[j + i1 * g_dim1]);
	  h[i - 1] = h[i1 - 1];
	}
	h[i2] = 0.0;
      }
L300:
      ;
    }
    /*  backward substitution to obtain the b-spline coefficients. */
    AlgBSplineBack(&g[g_offset], &c[1], nk1, k2, &c[1], nest);
    /*  computation of f(p). */
    *fp = 0.0;
    l = k2;
    for(it = 1; it <= m; ++it)
    {
      if((x[it] < t[l]) || (l > nk1))
      {
	goto L310;
      }
      ++l;
L310:
      l0 = l - k2;
      term = 0.0;
      for(j = 1; j <= k1; ++j)
      {
	++l0;
	term += c[l0] * q[it + j * q_dim1];
      }
      r1 = w[it] * (term - y[it]);
      *fp += r1 * r1;
    }
    /*  test whether the approximation sp(x) is an acceptable solution. */
    fpms = *fp - s;
    if(fabs(fpms) < acc)
    {
      goto L440;
    }
    /*  test whether the maximal number of iterations is reached. */
    if(iter == maxit)
    {
      goto L400;
    }
    /*  carry out one more step of the iteration process. */
    p2 = p;
    f2 = fpms;
    if(ich3 != 0)
    {
      goto L340;
    }
    if(f2 - f3 > acc)
    {
      goto L335;
    }
    /*  our initial choice of p is too large. */
    p3 = p2;
    f3 = f2;
    p *= 0.04;
    if(p <= p1)
    {
      p = p1 * 0.9 + p2 * 0.1;
    }
    goto L360;
L335:
    if(f2 < 0.0)
    {
      ich3 = 1;
    }
L340:
    if(ich1 != 0)
    {
      goto L350;
    }
    if(f1 - f2 > acc)
    {
      goto L345;
    }
    /*  our initial choice of p is too small */
    p1 = p2;
    f1 = f2;
    p /= 0.04;
    if(p3 < 0.0)
    {
      goto L360;
    }
    if(p >= p3)
    {
      p = p2 * 0.1 + p3 * 0.9;
    }
    goto L360;
L345:
    if(f2 > 0.0)
    {
      ich1 = 1;
    }
    /*  test whether the iteration process proceeds as theoretically */
    /*  expected. */
L350:
    if((f2 >= f1) || (f2 <= f3))
    {
      goto L410;
    }
    /*  find the new value for p. */
    p = AlgBSplineRati(&p1, &f1, p2, f2, &p3, &f3);
L360:
    ;
  }
  /*  error codes and messages. */
L400:
  ier = 3;
  goto L440;
L410:
  ier = 2;
  goto L440;
L420:
  ier = 1;
  goto L440;
L430:
  ier = -1;
L440:
  return(ier);
}

/*!
* \ingroup	AlgFit
* \brief	Solves linear system with upper triangular matrix.
* 
* This function has been derived from the netlib Dierckx function fpback().
* The original fortran coments are:
*
*  \verbatim
   Subroutine fpback calculates the solution of the system of equations
   a*c = z with a a n x n upper triangular matrix of bandwidth k.
   \endverbatim
*
* \param	a			Upper triangular matrix.
* \param	z			
* \param	n
* \param	k
* \param	c
* \param	nest
*/
static void			AlgBSplineBack(
				  double *a,
				  double *z,
				  int n,
				  int k,
				  double *c,
				  int nest)
{
  int		i,
  		k1;

  k1 = k - 1;
  c[n - 1] = z[n - 1] / a[n - 1];
  i = n - 1;
  if(i)
  {
    int		j;

    for(j = 2; j <= n; ++j)
    {
      int	i1,
		l,
		m;
      double	store;

      --i;
      store = z[i];
      i1 = k1;
      if(j <= k1)
      {
	i1 = j - 1;
      }
      m = i;
      for(l = 1; l <= i1; ++l)
      {
	++m;
	store -= c[m] * a[i + (l * nest)];
      }
      c[i] = store / a[i];
    }
  }
}

/*!
* \return	Dierckx error code.
* \ingroup	AlgFit
* \brief	Checks the number and the position of the knots.
*
* This function has been derived from the netlib Dierckx function fpchec().
* The original fortran coments are:
*
*  \verbatim
  Subroutine fpchec verifies the number and the position of the knots
  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
  and the position of the data points x(i),i=1,2,...,m. if all of the
  following conditions are fulfilled, the error parameter ier is set
  to zero. if one of the conditions is violated ier is set to ten.
      1) k+1 <= n-k-1 <= m
      2) t(1) <= t(2) <= ... <= t(k+1)
         t(n-k) <= t(n-k+1) <= ... <= t(n)
      3) t(k+1) < t(k+2) < ... < t(n-k)
      4) t(k+1) <= x(i) <= t(n-k)
      5) the conditions specified by schoenberg and whitney must hold
         for at least one subset of data points, i.e. there must be a
         subset of data points y(j) such that
             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
   \endverbatim
*
* \param	x			Data points.
* \param	m			Number of data points.
* \param	t			Knots.
* \param	n			Number of knots.
* \param	k			Order of spline.
*/
static int 			AlgBSplineChec(
				  double *x,
				  int m,
				  double *t,
				  int n,
				  int k)
{
  /* Local variables */
  int		i, j, l,
  		k1, m1,
		nk1, nk2, nk3,
		ier = 0;

  /* Function Body */
  k1 = k + 1;
  nk1 = n - k1;
  nk2 = nk1 + 1;
  /* Check condition no 1 */
  if((nk1 < k1) || (nk1 > m))
  {
    ier = 10;
    goto L80;
  }
  /* Check condition no 2 */
  j = n - 1;
  for(i = 0; i < k; ++i)
  {
    if((t[i] > t[i + 1]) || (t[j] < t[j - 1]))
    {
      ier = 10;
      goto L80;
    }
    --j;
  }
  /* Check condition no 3 */
  for(i = k1; i < nk2; ++i)
  {
    if(t[i] <= t[i - 1])
    {
      ier = 10;
      goto L80;
    }
  }
  /* Check condition nos 4 & 5 */
  m1 = m - 1;
  if((x[0] <  t[k])  || (x[m1] >  t[nk1]) || 
     (x[0] >= t[k1]) || (x[m1] <= t[nk1 - 1]))
  {
    ier = 10;
    goto L80;
  }
  i = 0;
  l = k1;
  nk3 = nk1 - 1;
  if(nk3 >= 2) {
    for(j = 1; j < nk3; ++j)
    {
      double	tj,
      		tl;

      tj = t[j];
      tl = t[l];
      ++l;
      do
      {
	++i;
	if(i >= m1)
	{
	  ier = 10;
	  goto L80;
	}
      }
      while(x[i] <= tj);
      if(x[i] >= tl)
      {
	ier = 10;
	goto L80;
      }
    }
  }
L80:
  return(ier);
}

/*!
* \return	Dierckx error code.
* \ingroup	AlgFit
* \brief	Computes a smooth approximating b-spline in multi-dimensional
*		space. See AlgBSplineNDFit().
* 
* This function has been derived from the netlib Dierckx function fppara().
* The original fortran coments are included in the code below.
*
* \param	iopt
* \param	idim
* \param	m
* \param	u
* \param	mx
* \param	x
* \param	w
* \param	ub
* \param	ue
* \param	k
* \param	s
* \param	nest
* \param	tol
* \param	maxit
* \param	k1
* \param	k2
* \param	n
* \param	t
* \param	nc
* \param	c
* \param	fp
* \param	fpint
* \param	z
* \param	a
* \param	b
* \param	g
* \param	q
* \param	nrdata
*/
static int 			AlgBSplinePara(
				  int iopt,
				  int idim,
				  int m,
				  double *u,
				  int *mx,
				  double *x,
				  double *w,
				  double *ub,
				  double *ue,
				  int k,
				  double s,
				  int nest,
				  double tol,
				  int maxit,
				  int k1,
				  int k2,
				  int *n,
				  double *t,
				  int *nc,
				  double *c,
				  double *fp,
				  double *fpint,
				  double *z,
				  double *a,
				  double *b,
				  double *g,
				  double *q,
				  int *nrdata)
{
  int		i, j, l,
		i1, i2, i3, i4, i5, i6, i7,
		ich1, ich3,
		it, iter,
		j1, j2, jj, k3, l0,
		mk1, n8, nw, nk1,
		nmax, nmin, npl1, nplus, nrint,
		a_dim1, a_offset,
		b_dim1, b_offset,
		g_dim1, g_offset,
		q_dim1, q_offset,
		ier = 0;
  double	f1, f2, f3,
  		fac, fp0, fpart, fpms, fpold,
  		p, p1, p2, p3,
  		pinv, piv,
  		acc, cs, rn, sn,
  		r1, store, term,
  		ui, wi;
  double	h[7],
  		xi[10];

  /* Parameter adjustments. */
  --c;
  --t;
  --u;
  --w;
  --x;
  --z;
  --fpint;
  --nrdata;
  q_dim1 = m; q_offset = 1 + q_dim1 * 1; q -= q_offset;
  a_dim1 = nest; a_offset = 1 + a_dim1 * 1; a -= a_offset;
  g_dim1 = nest; g_offset = 1 + g_dim1 * 1; g -= g_offset;
  b_dim1 = nest; b_offset = 1 + b_dim1 * 1; b -= b_offset;
  /* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   *  part 1: determination of the number of knots and their position     c
   *  **************************************************************      c
   *  given a set of knots we compute the least-squares curve sinf(u),    c
   *  and the corresponding sum of squared residuals fp=f(p=inf).         c
   *  if iopt=-1 sinf(u) is the requested curve.                          c
   *  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
   *    if fp <=s we will continue with the current set of knots.         c
   *    if fp > s we will increase the number of knots and compute the    c
   *       corresponding least-squares curve until finally fp<=s.         c
   *    the initial choice of knots depends on the value of s and iopt.   c
   *    if s=0 we have spline interpolation; in that case the number of   c
   *    knots equals nmax = m+k+1.                                        c
   *    if s > 0 and                                                      c
   *      iopt=0 we first compute the least-squares polynomial curve of   c
   *      degree k; n = nmin = 2*k+2                                      c
   *      iopt=1 we start with the set of knots found at the last         c
   *      call of the routine, except for the case that s > fp0; then     c
   *      we compute directly the polynomial curve of degree k.           c
   * cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  /* Determine nmin, the number of knots for polynomial approximation. */
  nmin = k1 * 2;
  if(iopt < 0)
  {
    goto L60;
  }
  /* Calculation of acc, the absolute tolerance for the root of f(p)=s. */
  acc = tol * s;
  /*  determine nmax, the number of knots for spline interpolation. */
  nmax = m + k1;
  if(s > 0.0)
  {
    goto L45;
  }
  /*  If s=0, s(u) is an interpolating curve.
   *  test whether the required storage space exceeds the available one. */
  *n = nmax;
  if(nmax > nest)
  {
    goto L420;
  }
  /* Find the position of the interior knots in case of interpolation. */
L10:
  mk1 = m - k1;
  if(mk1 == 0)
  {
    goto L60;
  }
  k3 = k / 2;
  i = k2;
  j = k3 + 2;
  if(k3 << 1 == k)
  {
    goto L30;
  }
  for(l = 1; l <= mk1; ++l)
  {
    t[i] = u[j];
    ++i;
    ++j;
  }
  goto L60;
L30:
  for(l = 1; l <= mk1; ++l)
  {
    t[i] = (u[j] + u[j - 1]) * 0.5;
    ++i;
    ++j;
  }
  goto L60;
  /*  If s>0 our initial choice of knots depends on the value of iopt.
   *  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
   *  polynomial curve which is a spline curve without interior knots.
   *  if iopt=1 and fp0>s we start computing the least squares spline curve
   *  according to the set of knots found at the last call of the routine. */
L45:
  if(iopt == 0)
  {
    goto L50;
  }
  if(*n == nmin)
  {
    goto L50;
  }
  fp0 = fpint[*n];
  fpold = fpint[*n - 1];
  nplus = nrdata[*n];
  if(fp0 > s)
  {
    goto L60;
  }
L50:
  *n = nmin;
  fpold = 0.0;
  nplus = 0;
  nrdata[1] = m - 2;
  /* Main loop for the different sets of knots. m is a save upper bound
   * for the number of trials. */
L60:
  for(iter = 1; iter <= m; ++iter)
  {
    if(*n == nmin)
    {
      ier = -2;
    }
    /* Find nrint, tne number of knot intervals. */
    nrint = *n - nmin + 1;
    /* Find the position of the additional knots which are needed for
     * the b-spline representation of s(u). */
    nk1 = *n - k1;
    i = *n;
    for(j = 1; j <= k1; ++j)
    {
      t[j] = *ub;
      t[i] = *ue;
      --i;
    }
    /* Compute the b-spline coefficients of the least-squares spline curve
     *  sinf(u). the observation matrix a is built up row by row and
     *  reduced to upper triangular form by givens transformations.
     *  at the same time fp=f(p=inf) is computed. */
    *fp = 0.0;
    /* Initialize the b-spline coefficients and the observation matrix a. */
    i4 = *nc;
    for(i = 1; i <= i4; ++i)
    {
      z[i] = 0.0;
    }
    for(i = 1; i <= nk1; ++i)
    {
      for(j = 1; j <= k1; ++j)
      {
	a[i + j * a_dim1] = 0.0;
      }
    }
    l = k1;
    jj = 0;
    for(it = 1; it <= m; ++it)
    {
      /*  fetch the current data point u(it),x(it). */
      ui = u[it];
      wi = w[it];
      for(j = 1; j <= idim; ++j)
      {
	++jj;
	xi[j - 1] = x[jj] * wi;
      }
      /* Search for knot interval t(l) <= ui < t(l+1). */
L85:
      if(ui < t[l + 1] || l == nk1)
      {
	goto L90;
      }
      ++l;
      goto L85;
      /* Evaluate the (k+1) non-zero b-splines at ui and store them in q. */
L90:
      AlgBSplineBspl(&t[1], k, ui, l, h);
      for(i = 1; i <= k1; ++i)
      {
	q[it + i * q_dim1] = h[i - 1];
	h[i - 1] *= wi;
      }
      /* Rotate the new row of the observation matrix into triangle. */
      j = l - k1;
      for(i = 1; i <= k1; ++i)
      {
	++j;
	piv = h[i - 1];
	if(piv == 0.0)
	{
	  goto L110;
	}
	/* Calculate the parameters of the givens transformation. */
	AlgBSplineGivens(piv, &a[j + a_dim1], &cs, &sn);
	/* Transformations to right hand side. */
	j1 = j;
	for(j2 = 1; j2 <= idim; ++j2)
	{
	  AlgBSplineRota(cs, sn, &xi[j2 - 1], &z[j1]);
	  j1 += *n;
	}
	if(i == k1)
	{
	  goto L120;
	}
	i2 = 1;
	i3 = i + 1;
	for(i1 = i3; i1 <= k1; ++i1)
	{
	  ++i2;
	  /* Transformations to left hand side. */
	  AlgBSplineRota(cs, sn, &h[i1 - 1], &a[j + i2 * a_dim1]);
	}
L110:
	;
      }
      /* Add contribution of this row to the sum of squares of residual
       * right hand sides. */
L120:
      for(j2 = 1; j2 <= idim; ++j2)
      {
	r1 = xi[j2 - 1];
	*fp += r1 * r1;
      }
    }
    if(ier == -2)
    {
      fp0 = *fp;
    }
    fpint[*n] = fp0;
    fpint[*n - 1] = fpold;
    nrdata[*n] = nplus;
    /* Backward substitution to obtain the b-spline coefficients. */
    j1 = 1;
    for(j2 = 1; j2 <= idim; ++j2)
    {
      AlgBSplineBack(&a[a_offset], &z[j1], nk1, k1, &c[j1], nest);
      j1 += *n;
    }
    /* Test whether the approximation sinf(u) is an acceptable solution. */
    if(iopt < 0)
    {
      goto L440;
    }
    fpms = *fp - s;
    if(fabs(fpms) < acc)
    {
      goto L440;
    }
    /* If f(p=inf) < s accept the choice of knots. */
    if(fpms < 0.0)
    {
      goto L250;
    }
    /* If n = nmax, sinf(u) is an interpolating spline curve. */
    if(*n == nmax)
    {
      goto L430;
    }
    /* Increase the number of knots.
     * If n=nest we cannot increase the number of knots because of
     * the storage capacity limitation. */
    if(*n == nest)
    {
      goto L420;
    }
    /* Determine the number of knots nplus we are going to add. */
    if(ier == 0)
    {
      goto L140;
    }
    nplus = 1;
    ier = 0;
    goto L150;
L140:
    npl1 = nplus << 1;
    rn = (double ) nplus;
    if(fpold - *fp > acc)
    {
      npl1 = rn * fpms / (fpold - *fp);
    }
    i6 = npl1, i7 = nplus / 2, i6 = ALG_MAX(i6,i7);
    i5 = nplus << 1, i4 = ALG_MAX(i6,1);
    nplus = ALG_MIN(i5,i4);
L150:
    fpold = *fp;
    /* Compute the sum of squared residuals for each knot interval
     * t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint. */
    fpart = 0.0;
    i = 1;
    l = k2;
    nw = 0;
    jj = 0;
    for(it = 1; it <= m; ++it)
    {
      if(u[it] < t[l] || l > nk1)
      {
	goto L160;
      }
      nw = 1;
      ++l;
L160:
      term = 0.0;
      l0 = l - k2;
      for(j2 = 1; j2 <= idim; ++j2)
      {
	fac = 0.0;
	j1 = l0;
	for(j = 1; j <= k1; ++j)
	{
	  ++j1;
	  fac += c[j1] * q[it + j * q_dim1];
	}
	++jj;
	r1 = w[it] * (fac - x[jj]);
	term += r1 * r1;
	l0 += *n;
      }
      fpart += term;
      if(nw == 0)
      {
	goto L180;
      }
      store = term * 0.5;
      fpint[i] = fpart - store;
      ++i;
      fpart = store;
      nw = 0;
L180:
      ;
    }
    fpint[nrint] = fpart;
    for(l = 1; l <= nplus; ++l)
    {
      /* Add a new knot. */
      AlgBSplineKnot(&u[1], m, &t[1], n, &fpint[1],
	  &nrdata[1], &nrint, 1);
      /* If n=nmax we locate the knots as for interpolation */
      if(*n == nmax)
      {
	goto L10;
      }
      /* Test whether we cannot further increase the number of knots. */
      if(*n == nest) 
      {
	goto L200;
      }
    }
    /* Restart the computations with the new set of knots. */
L200:
    ;
  }
  /* Test whether the least-squares kth degree polynomial curve is a
   * solution of our approximation problem. */
L250:
  if(ier == -2)
  {
    goto L440;
  }
  /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   * c part 2: determination of the smoothing spline curve sp(u).          c
   * c **********************************************************          c
   * c we have determined the number of knots and their position.          c
   * c we now compute the b-spline coefficients of the smoothing curve     c
   * c sp(u). the observation matrix a is extended by the rows of matrix   c
   * c b expressing that the kth derivative discontinuities of sp(u) at    c
   * c the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
   * c ponding weights of these additional rows are set to 1/p.            c
   * c iteratively we then have to determine the value of p such that f(p),c
   * c the sum of squared residuals be = s. we already know that the least c
   * c squares kth degree polynomial curve corresponds to p=0, and that    c
   * c the least-squares spline curve corresponds to p=infinity. the       c
   * c iteration process which is proposed here, makes use of rational     c
   * c interpolation. since f(p) is a convex and strictly decreasing       c
   * c function of p, it can be approximated by a rational function        c
   * c r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
   * c ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
   * c to calculate the new value of p such that r(p)=s. convergence is    c
   * c guaranteed by taking f1>0 and f3<0.                                 c
   * ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
  /* Evaluate the discontinuity jump of the kth derivative of the
   * b-splines at the knots t(l),l=k+2,...n-k-1 and store in b. */
  AlgBSplineDisc(&t[1], *n, k2, &b[b_offset], nest);
  /* Initial value for p. */
  p1 = 0.0;
  f1 = fp0 - s;
  p3 = -1.0;
  f3 = fpms;
  p = 0.0;
  for(i = 1; i <= nk1; ++i)
  {
    p += a[i + a_dim1];
  }
  rn = (double ) nk1;
  p = rn / p;
  ich1 = 0;
  ich3 = 0;
  n8 = *n - nmin;
  /* Iteration process to find the root of f(p) = s. */
  for(iter = 1; iter <= maxit; ++iter)
  {
    /* The rows of matrix b with weight 1/p are rotated into the
     * triangularised observation matrix a which is stored in g. */
    pinv = 1.0 / p;
    i5 = *nc;
    for(i = 1; i <= i5; ++i)
    {
      c[i] = z[i];
    }
    for(i = 1; i <= nk1; ++i)
    {
      g[i + k2 * g_dim1] = 0.0;
      for(j = 1; j <= k1; ++j)
      {
	g[i + j * g_dim1] = a[i + j * a_dim1];
      }
    }
    for(it = 1; it <= n8; ++it)
    {
      /* The row of matrix b is rotated into triangle by Givens
       * transformation */
      for(i = 1; i <= k2; ++i)
      {
	h[i - 1] = b[it + i * b_dim1] * pinv;
      }
      for(j = 1; j <= idim; ++j)
      {
	xi[j - 1] = 0.0;
      }
      for(j = it; j <= nk1; ++j)
      {
	piv = h[0];
	/* Calculate the parameters of the givens transformation. */
	AlgBSplineGivens(piv, &g[j + g_dim1], &cs, &sn);
	/* Transformations to right hand side. */
	j1 = j;
	for(j2 = 1; j2 <= idim; ++j2)
	{
	  AlgBSplineRota(cs, sn, &xi[j2 - 1], &c[j1]);
	  j1 += *n;
	}
	if(j == nk1)
	{
	  goto L300;
	}
	i2 = k1;
	if(j > n8)
	{
	  i2 = nk1 - j;
	}
	for(i = 1; i <= i2; ++i)
	{
	  /* Transformations to left hand side. */
	  i1 = i + 1;
	  AlgBSplineRota(cs, sn, &h[i1 - 1], &g[j + i1 * g_dim1])
	    ;
	  h[i - 1] = h[i1 - 1];
	}
	h[i2] = 0.0;
      }
L300:
      ;
    }
    /* Backward substitution to obtain the b-spline coefficients. */
    j1 = 1;
    for(j2 = 1; j2 <= idim; ++j2)
    {
      AlgBSplineBack(&g[g_offset], &c[j1], nk1, k2, &c[j1], nest);
      j1 += *n;
    }
    /* Computation of f(p). */
    *fp = 0.0;
    l = k2;
    jj = 0;
    for(it = 1; it <= m; ++it)
    {
      if(u[it] < t[l] || l > nk1)
      {
	goto L310;
      }
      ++l;
L310:
      l0 = l - k2;
      term = 0.0;
      for(j2 = 1; j2 <= idim; ++j2)
      {
	fac = 0.0;
	j1 = l0;
	for(j = 1; j <= k1; ++j)
	{
	  ++j1;
	  fac += c[j1] * q[it + j * q_dim1];
	}
	++jj;
	r1 = fac - x[jj];
	term += r1 * r1;
	l0 += *n;
      }
      r1 = w[it];
      *fp += term * (r1 * r1);
    }
    /* Test whether the approximation sp(u) is an acceptable solution. */
    fpms = *fp - s;
    if(fabs(fpms) < acc)
    {
      goto L440;
    }
    /* Test whether the maximal number of iterations is reached. */
    if(iter == maxit)
    {
      goto L400;
    }
    /* Carry out one more step of the iteration process. */
    p2 = p;
    f2 = fpms;
    if(ich3 != 0)
    {
      goto L340;
    }
    if(f2 - f3 > acc)
    {
      goto L335;
    }
    /* Our initial choice of p is too large. */
    p3 = p2;
    f3 = f2;
    p *= 0.04;
    if(p <= p1)
    {
      p = p1 * 0.9 + p2 * 0.1;
    }
    goto L360;
L335:
    if(f2 < 0.0)
    {
      ich3 = 1;
    }
L340:
    if(ich1 != 0)
    {
      goto L350;
    }
    if(f1 - f2 > acc)
    {
      goto L345;
    }
    /* Our initial choice of p is too small */
    p1 = p2;
    f1 = f2;
    p /= 0.04;
    if(p3 < 0.0)
    {
      goto L360;
    }
    if(p >= p3)
    {
      p = p2 * 0.1 + p3 * 0.9;
    }
    goto L360;
L345:
    if(f2 > 0.0)
    {
      ich1 = 1;
    }
    /* Test whether the iteration process proceeds as theoretically
     * expected. */
L350:
    if(f2 >= f1 || f2 <= f3)
    {
      goto L410;
    }
    /*  find the new value for p. */
    p = AlgBSplineRati(&p1, &f1, p2, f2, &p3, &f3);
L360:
    ;
  }
  /* Error codes and messages. */
L400:
  ier = 3;
  goto L440;
L410:
  ier = 2;
  goto L440;
L420:
  ier = 1;
  goto L440;
L430:
  ier = -1;
L440:
  return(ier);
}

/*!
* \ingroup	AlgFit
* \brief	Solves a linear system with a upper triangular matrix.
*
* This function has been derived from the netlib Dierckx function fpbacp().
* The original fortran coments are:
*
*  \verbatim
   Subroutine fpbacp calculates the solution of the system of equations
   g * c = z  with g  a n x n upper triangular matrix of the form
             ! a '   !
         g = !   ' b !
             ! 0 '   !
   with b a n x k matrix and a a (n-k) x (n-k) upper triangular
   matrix of bandwidth k1.
   \endverbatim
*
* \param	a
* \param	b
* \param	z
* \param	n
* \param	k
* \param	c
* \param	nest
*/
static void			AlgBSplineBacp(
				  double *a,
				  double *b,
				  double *z,
				  int n,
				  int k,
				  double *c,
				  int nest)
{
  int a_offset, b_offset;
  static int i, j, l;

  l = n;
  a_offset = 1 + nest; a -= a_offset;
  b_offset = 1 + nest; b -= b_offset;
  for(i = 0; i < k; ++i)
  {
    double	store;

    store = z[l - 1];
    j = k + 1 - i;
    if(i > 0)
    {
      int	l0,
      		l1;

      l0 = l - 1;
      for(l1 = j; l1 <= k; ++l1)
      {
	++l0;
	store -= c[l0] * b[l + l1 * nest];
      }
    }
    c[l - 1] = store / b[l + (j - 1) * nest];
    --l;
    if(l == 0)
    {
      break;
    }
  }
  if(l != 0)
  {
    int		n2;

    n2 = n - k;
    for(i = 0; i < n2; ++i)
    {
      double	store;

      store = z[i];
      l = n2;
      for(j = 1; j <= k; ++j)
      {
	++l;
	store -= c[l - 1] * b[i + 1 + j * nest];
      }
      c[i] = store;
    }
    i = n2;
    c[i - 1] /= a[i + nest];
    if(i != 1)
    {
      for(j = 2; j <= n2; ++j)
      {
	int	i1,
		l0;
        double	store;

	--i;
	store = c[i - 1];
	i1 = k;
	if(j <= k)
	{
	  i1 = j - 1;
	}
	l = i;
	for(l0 = 1; l0 <= i1; ++l0)
	{
	  ++l;
	  store -= c[l - 1] * a[i + (l0 + 1) * nest];
	}
	c[i - 1] = store / a[i + nest];
      }
    }
  }
}

