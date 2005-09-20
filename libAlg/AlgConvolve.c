#pragma ident "MRC HGU $Id$"
/*!
* \file         libAlg/AlgConvolve.c
* \author       Richard Baldock, Bill Hill
* \date         January 2000
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
* \brief	Provides functions for convolving data arrays.
* \ingroup	AlgConvolve
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <float.h>

/*!
* \return	Error code.
* \ingroup	AlgConvolve
* \brief	Convolves 1D kernel and data arrays, cnv = krn * data.
*		The return convolution array must not be aliased to
*		either the kernel or data arrays.
* \param	sizeArrayCnv		Length of return array must be
*					>= max(len(dat),len(krn)).
* \param	arrayCnv		Return convolution array.
* \param	sizeArrayKrn		Length of kernel array, must be
*					odd.
* \param	arrayKrn		Kernel array.
* \param	sizeArrayDat		Length of data array.
* \param	arrayDat		Data array.
* \param	pad			Type of padding.
*/
AlgError	AlgConvolve(int sizeArrayCnv, double *arrayCnv,
			    int sizeArrayKrn, double *arrayKrn,
			    int sizeArrayDat, double *arrayDat,
			    AlgPadType pad)
{
  int		pCnt,
  		kCnt0,
		kCnt1;
  double	cnv,
  		dat0,
		dat1;
  double	*cP,
  		*kP,
		*dP0,
		*dP1;
  AlgError	errCode = ALG_ERR_NONE;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgConvolve FE %d 0x%lx %d 0x%lx %d 0x%lx %d\n",
	   sizeArrayCnv, (unsigned long )arrayCnv,
	   sizeArrayKrn, (unsigned long )arrayKrn,
	   sizeArrayDat, (unsigned long )arrayDat,
	   (int )pad));
  if((sizeArrayCnv <= 0) || (arrayCnv == NULL) ||
     (sizeArrayKrn <= 0) || ((sizeArrayKrn % 2) != 1) || (arrayKrn == NULL) ||
     (sizeArrayDat <= 0) || (arrayDat == NULL) ||
     (sizeArrayCnv < sizeArrayKrn) || (sizeArrayCnv < sizeArrayDat))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    switch(pad)
    {
      case ALG_PAD_NONE:
	pad = ALG_PAD_ZERO;
        break;
      case ALG_PAD_ZERO:
        break;
      case ALG_PAD_END:
	dat0 = *arrayDat;
	dat1 = *(arrayDat + sizeArrayDat - 1);
        break;
      default:
        errCode = ALG_ERR_FUNC;
	break;
    }
  }
  if(errCode == ALG_ERR_NONE)
  {
    /* Pad leading data with zeros or first data value and convolve with the
     * kernel until the whole of the kernel is within the data. */
    cP = arrayCnv;
    pCnt = sizeArrayKrn / 2;
    while(pCnt > 0)
    {
      cnv = 0.0;
      if(pad == ALG_PAD_END)
      {
        kP = arrayKrn;
	kCnt0 = pCnt;
	while(kCnt0-- > 0)
	{
	  cnv += *kP++ * dat0;
	}
      }
      else
      {
        kP = arrayKrn + pCnt; 
      }
      dP0 = arrayDat;
      kCnt0 = sizeArrayKrn - pCnt;
      cnv += *kP * *dP0;
      while(--kCnt0 > 0)
      {
        cnv += *++kP * *++dP0;
      }
      *cP++ = cnv;
      --pCnt;
    }
    /* Between leading and trailing padding regions just convolue the data
     * with the kernel. */
    dP0 = arrayDat;
    cP = arrayCnv + sizeArrayKrn / 2;
    pCnt = sizeArrayDat - (sizeArrayKrn - 1);
    while(pCnt-- > 0)
    {
      kP = arrayKrn;
      dP1 = dP0++;
      kCnt0 = sizeArrayKrn;
      cnv = *kP * *dP1;
      while(--kCnt0 > 0)
      {
	cnv += *++kP * *++dP1;
      }
      *cP++ = cnv;
    }
    /* Pad trailing data with zeros or last data value and convolve with the
     * kernel until the whole of the kernel is outside the data. */
    dP0 = arrayDat + sizeArrayDat - sizeArrayKrn + 1;
    cP = arrayCnv + sizeArrayDat - (sizeArrayKrn / 2);
    pCnt = sizeArrayKrn / 2;
    kCnt0 = sizeArrayKrn - 1;
    while(pCnt-- > 0)
    {
      kP = arrayKrn; 
      kCnt1 = kCnt0;
      dP1 = dP0++;
      cnv = *kP * *dP1;
      while(--kCnt1 > 0)
      {
        cnv += *++kP * *++dP1;
      }
      if(pad == ALG_PAD_END)
      {
	kCnt1 = sizeArrayKrn - kCnt0;
	while(kCnt1-- > 0)
	{
	  cnv += *++kP * dat1;
	}
      }
      *cP++ = cnv;
      --kCnt0;
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgConvolve FX %d\n",
	   (int )errCode));
  return(errCode);
}

#ifdef ALGCONVOLVE_TEST0
int		main(int argc, char *argv[])
{
  int		idx;
  AlgPadType 	pad  = ALG_PAD_NONE;
  AlgError	errCode = ALG_ERR_NONE;
  static double	dat[25],
  		krn[5] = {0.50, 0.87, 1.00, 0.87, 0.50},
		cnv[25];
  const int	datSz = 25,
		krnSz = 5,
		cnvSz = 25;

  for(idx = 0; idx < krnSz; ++idx)
  {
    krn[idx] /= 3.74;
  }
  for(idx = 0; idx < datSz; ++idx)
  {
    dat[idx] = 1.0;
  }
  dat[0] = 2.0;
  dat[13] = 2.0;
  dat[24] = 2.0;
  errCode = AlgConvolve(25, cnv, 5, krn, 25, dat, pad);
}
#endif /* ALGCONVOLVE_TEST0 */
