#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgComplexUtils.c
* \author       Richard Baldock
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Provides basic complex number utilities.
* \ingroup     	AlgComplex
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Alg.h>

/*!
* \return	Square of modulus.
* \ingroup	AlgComplex
* \brief	Computes the square of the modulus of the given complex value.
* \param	z			Given complex value.
*/
double		AlgCModSq(ComplexD z)
{
  double result;

  result = z.re * z.re + z.im * z.im;

  return(result);
}

/*!
* \return	Modulus.
* \ingroup	AlgComplex
* \brief	Computes the modulus of the given complex value.
* \param	z			Given complex value.
*/
double		AlgCMod(ComplexD z)
{
  double result;

  result = z.re * z.re + z.im * z.im;
  result = (result > 0.0)? sqrt(result): 0.0;

  return(result);
}

/*!
* \return	Argument.
* \ingroup	AlgComplex
* \brief	Computes the argument of the given complex value.
* \param	z			Given complex value.
*/
double		AlgCArg(ComplexD z)
{
  double result;

  result = atan2(z.im, z.re);

  return(result);
}

/*!
* \return	Real component.
* \ingroup	AlgComplex
* \brief	Returns the real component of the given complex value.
* \param	z			Given complex value.
*/
double		AlgCRe(ComplexD z)
{
  return(z.re);
}

/*!
* \return	Imaginary component.
* \ingroup	AlgComplex
* \brief	Returns the imaginary component of the given complex value.
* \param	z			Given complex value.
*/
double		AlgCIm(ComplexD z)
{
  return(z.im);
}

/*!
* \return	Complex conjugate.
* \ingroup	AlgComplex
* \brief	Returns the complex conjugate of the given complex value.
* \param	z			Given complex value.
*/
ComplexD 	AlgCConj(ComplexD z)
{
  ComplexD result;

  result.re = z.re;
  result.im = -z.im;

  return(result);
}

/*!
* \return	Sum.
* \ingroup	AlgComplex
* \brief	Computes the sum of the two given complex values.
* \param	z			Given complex value.
*/
ComplexD 	AlgCAdd(ComplexD z1, ComplexD z2)
{
  ComplexD result;

  result.re = z1.re + z2.re;
  result.im = z1.im + z2.im;

  return(result);
}

/*!
* \return	Difference.
* \ingroup	AlgComplex
* \brief	Subtracts the second complex value from the first.
* \param	z1			First complex value.
* \param	z2			Second complex value.
*/
ComplexD 	AlgCSub(ComplexD z1, ComplexD z2)
{
  ComplexD result;

  result.re = z1.re - z2.re;
  result.im = z1.im - z2.im;

  return(result);
}

/*!
* \return	Product.
* \ingroup	AlgComplex
* \brief	Multiplies the second complex value with the first.
* \param	z1			First complex value.
* \param	z2			Second complex value.
*/
ComplexD 	AlgCMult(ComplexD z1, ComplexD z2)
{
  ComplexD result;

  result.re = z1.re * z2.re - z1.im * z2.im;
  result.im = z1.re * z2.im + z1.im * z2.re;

  return(result);
}

/*!
* \return	Ratio.
* \ingroup	AlgComplex
* \brief	Divides the first complex value by the second.
* \param	z1			First complex value.
* \param	z2			Second complex value.
*/
ComplexD 	AlgCDiv(ComplexD z1, ComplexD z2)
{
  double	modSq;
  ComplexD 	result;

  modSq = AlgCModSq(z2);
  result = AlgCMult(z1, AlgCConj(z2));
  result.re /= modSq;
  result.im /= modSq;

  return(result);
}

/*!
* \return	Power.
* \ingroup	AlgComplex
* \brief	Pomputes the value of the given complex value to the
*		power of the given real value.
* \param	z			Given complex value.
* \param	y			Given real value.
*/
ComplexD 	AlgCPow(ComplexD z, double y)
{
  ComplexD result;
  double a, b;

  a = AlgCArg(z);
  b = AlgCMod(z);

  a *= y;
  b = pow(b, y);

  result.re = b * cos(a);
  result.im = b * sin(a);

  return(result);
}
