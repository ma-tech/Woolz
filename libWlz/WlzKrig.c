#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzKrig_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzKrig.c
* \author       Bill Hill
* \date         October 2012
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
* \brief	Functions for interpolation based on kriging.
*
* 		Kriging is an interpolation method from geo-statistics
* 		(hence some terms such as nugget), which aims to find the
* 		Best Linear Unbiased Estimate (BLUE) of a value at some
* 		point in space based on the semi-variograms of the
* 		values within some local region. The procedure is first
* 		to choose a model for the semi-variogram behavior which
* 		when plotted against the lag (distance) is assumed to have
* 		a nugget (constant or intercept variance), sill (variance
* 		at an infinite lag) and range (lag for which values have
* 		little correlation). This should be set using
* 		WlzKrigSetModelFn(). The model semi-variogram should
* 		then be computed using WlzKrigOSetModelSV2D() and
* 		then for each position to be interpolated with this
* 		model (ie all positions where the model is valid)
* 		WlzKrigOSetPosSV2D() should be called followed by
* 		WlzKrigOWeightsSolve().
*
* 		For a simple system in which there are \f$N\f$ known values
* 		\f$\{v_i\}\f$ at positions \f$\{p_i\}\f$ (\f$i \in [1-N]\f$)
* 		then we want to estimate a system of weights \f$\{w_i\}\f$
* 		such that the value \f$v_x\f$ at some position \f$p_x\f$
* 		can be computed by a simple linear system of weights
* 		\f[
		v_x = \sum_i w_i v_i 
		\f]
*		with the weights computed such that \f$ \sum_i w_i = 1\f$.
*		Setting \f$h_{ij} = \|p_i - p_j\|\f$ and having choosen
*		a variance model \f$\gamma\f$ then
*		\f[
		\mathbf{M} \mathbf{w} = \mathbf{p}
		\f]
*		where
*		\f$\mathbf{M} = \gamma(h_{ij})\f$,
*		\f$\mathbf{w} = w_j\f$ and
*		\f$\mathbf{p} = h_{i}\f$.
*		In practice a slack variable \f$\lambda\f$ is added to
*		the weights column vector and both the model matrix and
*		position vector are padded using  to 1.
*		This equation is then solved for the weights.
* \ingroup	WlzValuesUtils
*/

#include <float.h>
#include <Wlz.h>

/*!
* \ingroup	WlzValuesUtils
* \brief	Sets a kriging function data stuctures parameters and
* 		function pointer. See the kriging model function types
* 		and their implementation for the meaning of the parameters.
* \param	fn			The kriging function data stucture to
* 					be set.
* \param	type			Kriging model function type.
* \param	c0			Nugget (offset) parameter.
* \param	c1			Sill (slope) parameter.
* \param	a			Range parameter.
*/
void		WlzKrigSetModelFn(WlzKrigModelFn *fn,
				   WlzKrigModelFnType type,
				   double c0, double c1, double a)
{
  fn->type = type;
  fn->c0 = c0;
  fn->c1 = c1;
  fn->a = a;
  switch(type)
  {
    case WLZ_KRIG_MODELFN_NUGGET:
      fn->fn = WlzKrigModelFnNugget;
      break;
    case WLZ_KRIG_MODELFN_LINEAR:
      fn->fn = WlzKrigModelFnLinear;
      break;
    case WLZ_KRIG_MODELFN_SPHERICAL:
      fn->fn = WlzKrigModelFnSpherical;
      break;
    case WLZ_KRIG_MODELFN_EXPONENTIAL:
      fn->fn = WlzKrigModelFnExponential;
      break;
    case WLZ_KRIG_MODELFN_GAUSSIAN:
      fn->fn = WlzKrigModelFnGaussian;
      break;
    case WLZ_KRIG_MODELFN_QUADRATIC:
      fn->fn = WlzKrigModelFnQuadratic;
      break;
    default:
      fn->type = WLZ_KRIG_MODELFN_INVALID;
      fn->fn = NULL;
      break;
  }
}

/*!
* \return	Function output value.
* \ingroup	WlzValuesUtili
* \brief	A nugget kriging variance model.
* 		\[
		\left{
		\begin{array}
		h = 0 & \gamma(h) = 0 \\
		h > 0 & \gamma(h) = c_0
		\end{array}
		\right.
		\]
* \param	f			The kriging function data stucture.
* \param	h			Function input value.
*/
double		WlzKrigModelFnNugget(WlzKrigModelFn *f, double h)
{
  double 	g;

  g = (h > DBL_EPSILON)? f->c0: 0.0;
  return(g);
}

/*!
* \return	Function output value.
* \ingroup	WlzValuesUtils
* \brief	A linear kriging variance model.
* 		\[
		\left{
		\begin{array}
		h < a    & \gamma(h) = c_0 + \frac{(c_1 - c_0)}{a} h
		h \geq a & \gamma(h) = c_1
		\end{array}
		\right.
		\]
* \param	f			The kriging function data stucture.
* \param	h			Function input value.
*/
double		WlzKrigModelFnLinear(WlzKrigModelFn *f, double h)
{
  double 	g;

  g = (h > f->a)? f->c1: f->c0 + (f->c1 - f->c0) * h / f->a;
  return(g);
}

/*!
* \return	Function output value.
* \ingroup	WlzValuesUtils
* \brief	A spherical kriging variance model.
* 		\[
		\left{
		\begin{array}
		0 < h \leq a & \gamma(h) = c_0 + c_1 \frac{1}{2}
		                           (\frac{3h}{a} - \frac{h^3}{a^3}) \\
		h > a        & \gamma(h) = c_0 + c_1
		\end{array}
		\right.
		\]
* \param	f			The kriging function data stucture.
* \param	h			Function input value.
*/
double		WlzKrigModelFnSpherical(WlzKrigModelFn *f, double h)
{
  double 	g;

  if((0.0 < h) && (h <= f->a))
  {
    double	t;
    t = h / f->a;
    g = f->c0 + (0.5 * t * f->c1 * (3.0 - (t * t)));
  }
  else
  {
    g = f->c0 + f->c1;
  }
  return(g);
}

/*!
* \return	Function output value.
* \ingroup	WlzValuesUtils
* \brief	An exponential kriging variance model.
* 		\[
		\gamma(h) = c_0 + c_1(1 - e^{-\frac{h}{a}})
		\]
* \param	f			The kriging function data stucture.
* \param	h			Function input value.
*/
double		WlzKrigModelFnExponential(WlzKrigModelFn *f, double h)
{
  double 	g;

  g = f->c0 + (f->c1 * (1.0 - exp((-h) / f->a)));
  return(g);
}

/*!
* \return	Function output value.
* \ingroup	WlzValuesUtils
* \brief	A Gaussian kriging variance model.
* 		\[
		\gamma(h) = c_0 + c_1(1 - e^{-(\frac{h}{a})^2})
		\]
* \param	f			The kriging function data stucture.
* \param	h			Function input value.
*/
double		WlzKrigModelFnGaussian(WlzKrigModelFn *f, double h)
{
  double 	g,
  		t;

  t = h / f->a;
  g = f->c0 + (f->c1 * (1.0 - exp(-(t * t))));
  return(g);
}

/*!
* \return	Function output value.
* \ingroup	WlzValuesUtils
* \brief	A Quadratic kriging variance model.
* 		\[
		\gamma(h) = c_0 + c_1
		                  \left(
				  \frac{    (\frac{h}{a})^2}
				       {1 + (\frac{h}{a})^2}
				  \right)
		\]
* \param	f			The kriging function data stucture.
* \param	h			Function input value.
*/
double		WlzKrigModelFnQuadratic(WlzKrigModelFn *f, double h)
{
  double 	g,
  		t;

  t = h / f->a;
  t = t * t;
  g = f->c0 + (f->c1 * t / (1 + t));
  return(g);
}

/*!
* \return	Woolz error code. If all input parameters are valid then
* 		this will return WLZ_ERR_NONE.
* \ingroup	WlzValuesUtils
* \brief	Computes the ordinary kriging model semi-variogram from the
* 		given neighbourhood vertices.
* \param	modelSV			Valid \f$(n + 1)\times(n + 1)\f$
* 					matrix for the model semi-variogram.
* 					This must be a square matrix.
* \param	modelFn			The model value function.
* \param	n			Number of neighbourhood vertices.
* \param	nbr			The neighbourhood vertex positions.
* \param	wSp			Workspace with room for n values
* 					(used by AlgMatrixLUDecomp()).
*/
WlzErrorNum	WlzKrigOSetModelSV2D(AlgMatrix modelSV,
				     WlzKrigModelFn *modelFn,
				     int n, WlzDVertex2 *nbr,
				     int *wSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((modelSV.core == NULL) || (nbr == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(modelSV.core->type != ALG_MATRIX_RECT)
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  else if((modelSV.core->nR != modelSV.core->nC) ||
          (modelSV.core->nR != n + 1))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    int		i,
    		j;
    double	**a;

    /* Set values in matrix. */
    a = modelSV.rect->array;
    for(j = 0; j < n; ++j)
    {
      for(i = 0; i < j; ++i)
      {
	double d;
	WlzDVertex2 s;
        
	WLZ_VTX_2_SUB(s, nbr[j], nbr[i]);
	d = WLZ_VTX_2_LENGTH(s);
        a[j][i] = a[i][j] = modelFn->fn(modelFn, d);;
      }
      a[j][j] = 0.0;
      a[n][j] = a[j][n] = 1.0;
    }
    a[n][n] = 0.0;
    /* Perform LU decomposition. */
    errNum = WlzErrorFromAlg(AlgMatrixLUDecompRaw(a, n + 1, wSp, NULL));
  }
  return(errNum);
}

/*!
* \return	Woolz error code. If all input parameters are valid then
* 		this will return WLZ_ERR_NONE.
* \ingroup	WlzValuesUtils
* \brief	Computes the ordinary kriging position semi-variogram from
* 		the given neighbourhood vertices and the given position
* 		vertex.
* \param	posSV			Valid \f$(n + 1)\f$ column vector
* 					for the position semi-variogram.
* \param	modelFn			The model value function.
* \param	n			Number of neighbourhood vertices.
* \param	nbr			The neighbourhood vertex positions.
* \param	pos			The given position.
*/
WlzErrorNum	WlzKrigOSetPosSV2D(double *posSV, WlzKrigModelFn *modelFn,
				   int n, WlzDVertex2 *nbr, WlzDVertex2 pos)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((posSV == NULL) || (nbr == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(n < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    int		j;

    for(j = 0; j < n; ++j)
    {
      double d;
      WlzDVertex2 s;

      WLZ_VTX_2_SUB(s, pos, nbr[j]);
      d = WLZ_VTX_2_LENGTH(s);
      posSV[j] = modelFn->fn(modelFn, d);;
    }
    posSV[n] = 1.0;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	Computes the ordinary kriging weights from the model and
* 		position semi-variograms.
* \param	modelSV			Valid \f$(n + 1)\times(n + 1)\f$
* 					matrix containing the model
* 					semi-variogram.
* \param	posSV			Valid \f$(n + 1)\f$ column vector
* 					containing the position semi-variogram
* 					and on return containing the weights.
* \param	wSp			Workspace for AlgMatrixLUBackSub()
* 					with space for modelSV.core->nR
* 					values.
* \param	eps			Epsilon value used to test for
* 					position semi-variogram value zero.
*/
WlzErrorNum	WlzKrigOWeightsSolve(AlgMatrix modelSV, double *posSV,
				     int *wSp, double eps)
{
  int		n;
  AlgMatrix	mat;

  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((modelSV.core == NULL) || (wSp == NULL) || (posSV == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(modelSV.core->type != ALG_MATRIX_RECT)
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  else if((n = modelSV.core->nR) < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    int		idR,
    		idZ = -1;

    /* Checking for a position sem-variogram entry beeing zero avoids
     * the degenerate solution of all weights = 0.0 and lambda = 1.0. */
    for(idR = 0; idR < n; ++idR)
    {
      if(posSV[idR] < eps)
      {
        idZ = idR;
	break;
      }
    }
    if(idZ >= 0)
    {
      for(idR = 0; idR < n; ++idR)
      {
        posSV[idR] = 0.0;
      }
      posSV[idZ] = 1.0;
    }
    else
    {
      errNum = WlzErrorFromAlg(AlgMatrixLUBackSub(modelSV, wSp, posSV));
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValuesUtils
* \brief	A convinience function for allocating and reallocating
* 		the buffers needed for 2D kriging. All given pointers
* 		must be valid and on the first call all pointer values
* 		(eg *dstNbrPosBuf) must be NULL.
* \param	dstNbrPosBuf		Destination pointer for the
* 				        neighbourhood position buffer.
* \param	dstPosSV		Destination pointer for the
* 				        position vector.
* \param	dstWSp			Destination pointer for kriging
* 					matrix workspace.
* \param	dstModelSV		Destination pointer for kriging
* 					model matrix.
* \param	dstMaxNbrIdxBuf		Destination pointer for the
* 					maximum number of neighbours the
* 					buffers have been reallocated
* 					for.
* \param	nNbrC			Current number of neighbours.
* \param	nNbrL			Last number of neighbours.
*/
WlzErrorNum	WlzKrigReallocBuffers2D(WlzDVertex2 **dstNbrPosBuf,
                                        double **dstPosSV, int **dstWSp,
					AlgMatrix *dstModelSV,
					int *dstMaxKrigBuf,
					int nNbrC, int nNbrL)
{
  int		nC1,
		maxC1,
  		maxKrigBuf = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(nNbrC > *dstMaxKrigBuf)
  {
    maxKrigBuf = nNbrC * 2;
    maxC1 = maxKrigBuf + 1;
    if(((*dstNbrPosBuf = (WlzDVertex2 *)
			 AlcRealloc(*dstNbrPosBuf,
				    maxC1 * sizeof(WlzDVertex2))) == NULL) ||
	((*dstPosSV = (double *)
		      AlcRealloc(*dstPosSV,
		                 maxC1 * sizeof(double))) == NULL) ||
	((*dstWSp = (int *)
		    AlcRealloc(*dstWSp, maxC1 * sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && (nNbrC != nNbrL))
  {
    if(maxKrigBuf > 0)
    {
      AlgMatrixFree(*dstModelSV);
      *dstModelSV = AlgMatrixNew(ALG_MATRIX_RECT, maxC1, maxC1, 0, 0.0, NULL);
      if((*dstModelSV).core == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    nC1 = nNbrC + 1;
    (*dstModelSV).rect->nR = nC1;
    (*dstModelSV).rect->nC = nC1;
  }
  if((errNum == WLZ_ERR_NONE) && (maxKrigBuf > 0))
  {
    *dstMaxKrigBuf = maxKrigBuf;
  }
  return(errNum);
}
