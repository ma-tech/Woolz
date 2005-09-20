#pragma ident "MRC HGU $Id$"
/*!
* \file         libAlg/AlgComplexUtils.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Provides basic complex number utilities.
* \ingroup	AlgComplex
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
* \param	z1			First complex value.
* \param	z2			Second complex value.
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
