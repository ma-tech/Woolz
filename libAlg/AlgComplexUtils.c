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
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup      Alg
* \defgroup     AlgComplex
* @{
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Alg.h>

double	AlgCMod(
  ComplexD z)
{
  double result;

  result = z.re * z.re + z.im * z.im;
  if( result > 0.0 ){
    result = sqrt( result );
  }

  return result;
}

double	AlgCArg(
  ComplexD z)
{
  double result;

  result = atan2(z.im, z.re);

  return result;
}

double	AlgCRe(
  ComplexD z)
{
  return z.re;
}

double	AlgCIm(
  ComplexD z)
{
  return z.im;
}

ComplexD AlgCConj(
  ComplexD z)
{
  ComplexD result;

  result.re = z.re;
  result.im = -z.im;

  return result;
}

ComplexD AlgCAdd(
  ComplexD z1,
  ComplexD z2)
{
  ComplexD result;

  result.re = z1.re + z2.re;
  result.im = z1.im + z2.im;

  return result;
}

ComplexD AlgCSub(
  ComplexD z1,
  ComplexD z2)
{
  ComplexD result;

  result.re = z1.re - z2.re;
  result.im = z1.im - z2.im;

  return result;
}

ComplexD AlgCMult(
  ComplexD z1,
  ComplexD z2)
{
  ComplexD result;

  result.re = z1.re * z2.re - z1.im * z2.im;
  result.im = z1.re * z2.im + z1.im * z2.re;

  return result;
}

ComplexD AlgCDiv(
  ComplexD z1,
  ComplexD z2)
{
  ComplexD result;

  result = AlgCMult(z1, AlgCConj(z2));
  result.re /= AlgCMod(z2);
  result.im /= AlgCMod(z2);

  return result;
}

ComplexD AlgCPow(
  ComplexD z, 
  double y)
{
  ComplexD result;
  double a, b;

  a = AlgCArg(z);
  b = AlgCMod(z);

  a *= y;
  b = pow(b, y);

  result.re = b * cos(a);
  result.im = b * sin(a);

  return result;
}

/*!
* @}
*/
