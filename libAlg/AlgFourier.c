#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgFourier_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgFourier.c
* \author       Bill Hill
* \date         March 1999
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
* \brief        Fast Fourier and Hartley transform functions.
*
* \par
*		The history of this software is (as stated by Mayer
*		with my later addition):
<table width="500" border="0">
                  <tr>
		  <td>Euler</td> <td>Probable inventor of the fourier
				     transform.</td>
		  </tr>
		  <tr>
		  <td>Gauss</td> <td>Probable inventor of the FFT.</td>
		  </tr>
		  <tr>
		  <td>Hartley</td> <td>Probable inventor of the hartley
				       transform.</td>
		  </tr>
		  <tr>
		  <td>Bracewell & Buneman</td> <td>Patent holders for FHT!</td>
		  </tr>
		  <tr>
		  <td>Ron Mayer</td> <td>Produced FFT code using the FHT.</td>
		  </tr>
		  <tr>
		  <td>Bill Hill</td> <td>Hacked Ron Mayer's code to simplify
				  multi dimensional FFT's and for
				  compatability with our existing FFT
				  routines here at MRC HGU.
				  Multithreaded and integrated into
				  libAlg for Woolz.</td>
		  </tr>
</table>
* 		The multi-dimensional transform routines may be supplied
*		with a use buffers flag, which if set will allocate
*		buffers sufficient to hold copies of the data to
*		ensure the transforms are only computed for contiguous
*		data. Using buffers can result in an order of magnitude
*		lower run times depending on whether the array fit into
*		the fastest caches of the CPUs.
* \ingroup      AlgFourier
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif


/*!
* \enum 	 _AlgFourDir
* \brief	 Transform direction: Forward or inverse Fourier transform.
*/
typedef enum _AlgFourDir
{
  ALG_FOUR_DIR_FWD = 0,
  ALG_FOUR_DIR_INV = 1
} AlgFourDir;

/*!
* \enum         _AlgFourAxis
* \brief        Axis for partial transform evaluation.
*/
typedef enum _AlgFourAxis
{
  ALG_FOUR_AXIS_X = 0,
  ALG_FOUR_AXIS_Y = 1,
  ALG_FOUR_AXIS_Z = 2
} AlgFourAxis;

static AlgError			AlgFourRepXY1D(
				  double **real,
				  double **imag,
			          AlgFourAxis axis,
				  int useBuf,
			          int numX,
				  int numY,
				  AlgFourDir dir);
static AlgError			AlgFourRepXYReal1D(
				  double **data,
				  AlgFourAxis axis,
				  int useBuf,
				  int numX,
				  int numY,
				  AlgFourDir dir);
static AlgError			AlgFourRepXYZ1D(
				  double ***real,
				  double ***imag,
			          AlgFourAxis axis,
				  int useBuf,
			          int numX,
				  int numY,
				  int numZ,
				  AlgFourDir dir);
static AlgError			AlgFourRepXYZReal1D(
				  double ***data,
				  AlgFourAxis axis,
				  int useBuf,
				  int numX,
				  int numY,
				  int numZ,
				  AlgFourDir dir);

/*!
* \return	void
* \ingroup      AlgFourier
* \brief	Computes the Hartley transform of the given one
*		dimensional data, and does it in place.
* \param	data		Given data.
* \param	num		Number of data.
* \param	step		Offset in data elements between
*				the data to be transformed.
*/
void		AlgFourHart1D(double *data, int num, int step)
{
#ifdef ALG_FOUR_LONGDBL_TRIG
  long
#endif /* ALG_FOUR_LONGDBL_TRIG */
  double	cos0,
  		cos1,
		cos2,
		sin0,
		sin1,
		sin2,
		tTrig;
  double	tD0,
		tD1,
		tD2,
		tD3,
		tD4,
		tD5,
		tD6,
		tD7,
		tD8,
		tD9,
		tD10,
		tD11;
  int		pTwo0,
		pTwo1,
		pTwo2,
		pTwo3,
		pTwo4,
		pTwoX,
		count,
		index,
		tI1,
		tI2,
		tI3,
		tI4;
  double	*tDp0,
		*tDp1,
		*tDp2,
		*tDp3,
		*tGp0,
		*tGp1,
		*tGp2,
		*tGp3;
#ifdef ALG_FOUR_LONGDBL_TRIG
  static const long double regFourCosTab[20]=
  {
    0.00000000000000000000000000000000000000000000000000L,
    0.70710678118654752440084436210484903928483593768847L,
    0.92387953251128675612818318939678828682241662586364L,
    0.98078528040323044912618223613423903697393373089333L,
    0.99518472667219688624483695310947992157547486872985L,
    0.99879545620517239271477160475910069444320361470461L,
    0.99969881869620422011576564966617219685006108125772L,
    0.99992470183914454092164649119638322435060646880221L,
    0.99998117528260114265699043772856771617391725094433L,
    0.99999529380957617151158012570011989955298763362218L,
    0.99999882345170190992902571017152601904826792288976L,
    0.99999970586288221916022821773876567711626389934930L,
    0.99999992646571785114473148070738785694820115568892L,
    0.99999998161642929380834691540290971450507605124278L,
    0.99999999540410731289097193313960614895889430318945L,
    0.99999999885102682756267330779455410840053741619428L
  },
		regFourSinTab[20]=
  {
    1.00000000000000000000000000000000000000000000000000L,
    0.70710678118654752440084436210484903928483593768846L,
    0.38268343236508977172845998403039886676134456248561L,
    0.19509032201612826784828486847702224092769161775195L,
    0.09801714032956060199419556388864184586113667316749L,
    0.04906767432741801425495497694268265831474536302574L,
    0.02454122852291228803173452945928292506546611923944L,
    0.01227153828571992607940826195100321214037231959176L,
    0.00613588464915447535964023459037258091705788631738L,
    0.00306795676296597627014536549091984251894461021344L,
    0.00153398018628476561230369715026407907995486457522L,
    0.00076699031874270452693856835794857664314091945205L,
    0.00038349518757139558907246168118138126339502603495L,
    0.00019174759731070330743990956198900093346887403385L,
    0.00009587379909597734587051721097647635118706561284L,
    0.00004793689960306688454900399049465887274686668768L
  };
#else /* ! ALG_FOUR_LONGDBL_TRIG */
  static const double regFourCosTab[20]=
  {
    0.00000000000000000000,
    0.70710678118654752440,
    0.92387953251128675612,
    0.98078528040323044912,
    0.99518472667219688624,
    0.99879545620517239271,
    0.99969881869620422011,
    0.99992470183914454092,
    0.99998117528260114265,
    0.99999529380957617151,
    0.99999882345170190992,
    0.99999970586288221916,
    0.99999992646571785114,
    0.99999998161642929380,
    0.99999999540410731289,
    0.99999999885102682756
  },
		regFourSinTab[20]=
  {
    1.00000000000000000000,
    0.70710678118654752440,
    0.38268343236508977172,
    0.19509032201612826784,
    0.09801714032956060199,
    0.04906767432741801425,
    0.02454122852291228803,
    0.01227153828571992607,
    0.00613588464915447535,
    0.00306795676296597627,
    0.00153398018628476561,
    0.00076699031874270452,
    0.00038349518757139558,
    0.00019174759731070330,
    0.00009587379909597734,
    0.00004793689960306688
  };
#endif /* ALG_FOUR_LONGDBL_TRIG */

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourHart1D FE %p %d %d\n",
	   data, num, step));
  pTwo1 = 1;
  pTwo2 = 0;
  if(step == 1)
  {
    while(pTwo1 < num)
    {
      pTwo0 = num >> 1;
      while(((pTwo2 ^= pTwo0) & pTwo0) == 0)
      {
	pTwo0 >>= 1;
      }
      if(pTwo1 > pTwo2)
      {
	tDp0 = data + pTwo1;
	tDp1 = data + pTwo2;
	tD0 = *tDp0;
	*tDp0 = *tDp1;
	*tDp1 = tD0;
      }
      ++pTwo1;
    }
  }
  else
  {
    while(pTwo1 < num)
    {
      pTwo0 = num >> 1;
      while(((pTwo2 ^= pTwo0) & pTwo0) == 0)
      {
	pTwo0 >>= 1;
      }
      if(pTwo1 > pTwo2)
      {
	tDp0 = data + (pTwo1 * step);
	tDp1 = data + (pTwo2 * step);
	tD0 = *tDp0;
	*tDp0 = *tDp1;
	*tDp1 = tD0;
      }
      ++pTwo1;
    }
  }
  pTwo0 = 0;
  while((1 << pTwo0) < num)
  {
    ++pTwo0;
  }
  pTwo0  &= 1;
  tDp0 = data;
  count = num;
  if(pTwo0 == 0)
  {
    if(step == 1)
    {
      while(count > 0)
      {
	tD0 = *tDp0;
	tD1 = tD0 - *(tDp0 + 1);
	tD0 += *(tDp0 + 1);
	tD2 = *(tDp0 + 2);
	tD3 = tD2 - *(tDp0 + 3);
	tD2 += *(tDp0 + 3);
	*(tDp0 + 2) = tD0 - tD2;	
	*tDp0 = tD0 + tD2;
	*(tDp0 + 3) = tD1 - tD3;	
	*(tDp0 + 1) = tD1 + tD3;
	tDp0 = tDp0 + 4;
	count -= 4;
      }
    }
    else /* (step != 1) */
    {
      while(count > 0)
      {
	tDp1 = tDp0 + step;
	tDp2 = tDp1 + step;
	tDp3 = tDp2 + step;
	tD0 = *tDp0;
	tD1 = tD0 - *tDp1;
	tD0 += *tDp1;
	tD2 = *tDp2;
	tD3 = tD2 - *tDp3;
	tD2 += *tDp3;
	*tDp2 = tD0 - tD2;	
	*tDp0 = tD0 + tD2;
	*tDp3 = tD1 - tD3;	
	*tDp1 = tD1 + tD3;
	tDp0 = tDp3 + step;
	count -= 4;
      }
    }
  }
  else /* (pTwo0 != 0) */
  {
    if(step == 1)
    {
      while(count > 0)
      {
	tD1 = *(tDp0 + 1);
	tD0 = *tDp0 - tD1;
	tD1 += *tDp0;
	tD3 = *(tDp0 + 3);
	tD2 = *(tDp0 + 2) - tD3;
	tD3 += *(tDp0 + 2);
	tD5 = *(tDp0 + 5);
	tD4 = *(tDp0 + 4) - tD5;
	tD5 += *(tDp0 + 4);
	tD7 = *(tDp0 + 7);
	tD6 = *(tDp0 + 6) - tD7;
	tD7 += *(tDp0 + 6);
	tD8 = (tD1 - tD3);
	tD9 = (tD1 + tD3);
	tD10 = (tD5 - tD7);
	tD11 = (tD5 + tD7);
	*(tDp0 + 4) = tD9 - tD11;
	*tDp0 = tD9 + tD11;
	*(tDp0 + 6) = tD8 - tD10;
	*(tDp0 + 2) = tD8 + tD10;
	tD8 = (tD0 - tD2);
	tD9 = (tD0 + tD2);
	tD10 = ALG_M_SQRT2 * tD6;
	tD11 = ALG_M_SQRT2 * tD4;
	*(tDp0 + 5) = tD9 - tD11;
	*(tDp0 + 1) = tD9 + tD11;
	*(tDp0 + 7) = tD8 - tD10;
	*(tDp0 + 3) = tD8 + tD10;
	tDp0 += 8;
	count -= 8;
      }
    }
    else /* (step != 1) */
    {
      tI1 = 2 * step;
      tGp0 = data + step;
      while(count > 0)
      {
	tDp1 = tDp0 + tI1;
	tGp1 = tDp1 + step;
	tDp2 = tDp1 + tI1;
	tGp2 = tDp2 + step;
	tDp3 = tDp2 + tI1;
	tGp3 = tDp3 + step;
	tD1 = *tGp0;
	tD0 = *tDp0 - tD1;
	tD1 += *tDp0;
	tD3 = *tGp1;
	tD2 = *tDp1 - tD3;
	tD3 += *tDp1;
	tD5 = *tGp2;
	tD4 = *tDp2 - tD5;
	tD5 += *tDp2;
	tD7 = *tGp3;
	tD6 = *tDp3 - tD7;
	tD7 += *tDp3;
	tD8 = (tD1 - tD3);
	tD9 = (tD1 + tD3);
	tD10 = (tD5 - tD7);
	tD11 = (tD5 + tD7);
	*tDp2 = tD9 - tD11;
	*tDp0 = tD9 + tD11;
	*tDp3 = tD8 - tD10;
	*tDp1 = tD8 + tD10;
	tD8 = (tD0 - tD2);
	tD9 = (tD0 + tD2);
	tD10 = ALG_M_SQRT2 * tD6;
	tD11 = ALG_M_SQRT2 * tD4;
	*tGp2 = tD9 - tD11;
	*tGp0 = tD9 + tD11;
	*tGp3 = tD8 - tD10;
	*tGp1 = tD8 + tD10;
	tDp0 = tDp3 + tI1;
	tGp0 = tGp3 + tI1;
	count -= 8;
      }
    }
  }
  if(num >= 16)
  {
    if(step == 1)
    {
      do
      {
	pTwo0  += 2;
	pTwo1  = 1 << pTwo0;
	pTwo2  = pTwo1 << 1;
	pTwo4  = pTwo2 << 1;
	pTwo3  = pTwo2 + pTwo1;
	pTwoX  = pTwo1 >> 1;
	count = num;
	tDp0  = data;
	tGp0  = tDp0 + pTwoX;
	do
	{
	  tD0 = *tDp0;
	  tD4 = *(tDp0 + pTwo1);
	  tD1 = tD0 - tD4;
	  tD0 += tD4;
	  tD2 = *(tDp0 + pTwo2);
	  tD4 = *(tDp0 + pTwo3);
	  tD3 = tD2 - tD4;
	  tD2 += tD4;
	  *tDp0 = tD0 + tD2;
	  *(tDp0 + pTwo1) = tD1 + tD3;
	  *(tDp0 + pTwo2) = tD0 - tD2;
	  *(tDp0 + pTwo3) = tD1 - tD3;
	  tD0 = *tGp0;
	  tD4 = *(tGp0 + pTwo1);
	  tD1 = tD0 - tD4;
	  tD0 += tD4;
	  tD3 = ALG_M_SQRT2 * *(tGp0 + pTwo3);
	  tD2 = ALG_M_SQRT2 * *(tGp0 + pTwo2);
	  *tGp0 = tD0 + tD2;
	  *(tGp0 + pTwo1) = tD1 + tD3;
	  *(tGp0 + pTwo2) = tD0 - tD2;
	  *(tGp0 + pTwo3) = tD1 - tD3;
	  tGp0 += pTwo4;
	  tDp0 += pTwo4;
	  count -= pTwo4;
	} while(count > 0);
	cos0 = regFourCosTab[pTwo0];
	sin0 = regFourSinTab[pTwo0];
	cos1 = 1.0;
	sin1 = 0.0;
	for(index = 1; index < pTwoX; ++index)
	{
	  tTrig = cos1;
	  cos1 = tTrig * cos0 - sin1 * sin0;
	  sin1 = tTrig * sin0 + sin1 * cos0;
	  cos2 = (cos1 * cos1) - (sin1 * sin1);
	  sin2 = 2.0 * (cos1 * sin1);
	  tDp0 = data + index;
	  tGp0 = data + pTwo1 - index;
	  count = num;
	  do
	  {
	    tD0 = *(tDp0 + pTwo1);
	    tD1 = *(tGp0 + pTwo1);
	    tD9 = (sin2 * tD0) - (cos2 * tD1);
	    tD8 = (cos2 * tD0) + (sin2 * tD1);
	    tD0 = *tDp0;
	    tD1 = tD0 - tD8;
	    tD0 += tD8;
	    tD4 = *tGp0;
	    tD5 = tD4 - tD9;
	    tD4 += tD9;
	    tD2 = *(tDp0 + pTwo3);
	    tD3 = *(tGp0 + pTwo3);
	    tD9 = (sin2 * tD2) - (cos2 * tD3);
	    tD8 = (cos2 * tD2) + (sin2 * tD3);
	    tD2 = *(tDp0 + pTwo2);
	    tD3 = tD2 - tD8;
	    tD2 += tD8;
	    tD6 = *(tGp0 + pTwo2);
	    tD7 = tD6 - tD9;
	    tD6 += tD9;
	    tD9 = (sin1 * tD2) - (cos1 * tD7);
	    tD8 = (cos1 * tD2) + (sin1 * tD7);
	    *(tDp0 + pTwo2) = tD0 - tD8;
	    *tDp0 = tD0 + tD8;
	    *(tGp0 + pTwo3) = tD5 - tD9;
	    *(tGp0 + pTwo1) = tD5 + tD9;
	    tD9 = (cos1 * tD6) - (sin1 * tD3);
	    tD8 = (sin1 * tD6) + (cos1 * tD3);
	    *(tGp0 + pTwo2) = tD4 - tD8;
	    *tGp0 = tD4 + tD8;
	    *(tDp0 + pTwo3) = tD1 - tD9;
	    *(tDp0 + pTwo1) = tD1 + tD9;
	    tGp0 += pTwo4;
	    tDp0 += pTwo4;
	    count -= pTwo4;
	  }while(count > 0);
	}
      }while(pTwo4 < num);
    }
    else /* (step != 1) */
    {
      do
      {
	pTwo0  += 2;
	pTwo1  = 1  << pTwo0;
	pTwo2  = pTwo1 << 1;
	pTwo4  = pTwo2 << 1;
	pTwo3  = pTwo2 + pTwo1;
	pTwoX  = pTwo1 >> 1;
	tI1 = pTwo1 * step;
	tI2 = pTwo2 * step;
	tI3 = pTwo3 * step;
	tI4 = pTwo4 * step;
	count = num;
	tDp0  = data;
	tGp0  = tDp0 + (pTwoX * step);
	do
	{
	  tDp1 = tDp0 + tI1;
	  tDp2 = tDp0 + tI2;
	  tDp3 = tDp0 + tI3;
	  tD0 = *tDp0;
	  tD1 = tD0 - *tDp1;
	  tD0 += *tDp1;
	  tD2 = *tDp2;
	  tD3 = tD2 - *tDp3;
	  tD2 += *tDp3;
	  *tDp0 = tD0 + tD2;
	  *tDp1 = tD1 + tD3;
	  *tDp2 = tD0 - tD2;
	  *tDp3 = tD1 - tD3;
	  tGp1 = tGp0 + tI1;
	  tGp2 = tGp0 + tI2;
	  tGp3 = tGp0 + tI3;
	  tD0 = *tGp0;
	  tD1 = tD0 - *tGp1;
	  tD0 += *tGp1;
	  tD3 = ALG_M_SQRT2 * *tGp3;
	  tD2 = ALG_M_SQRT2 * *tGp2;
	  *tGp2 = tD0 - tD2;
	  *tGp0 = tD0 + tD2;
	  *tGp3 = tD1 - tD3;
	  *tGp1 = tD1 + tD3;
	  tGp0 += tI4;
	  tDp0 += tI4;
	  count -= pTwo4;
	} while(count > 0);
	cos0 = regFourCosTab[pTwo0];
	sin0 = regFourSinTab[pTwo0];
	cos1 = 1.0;
	sin1 = 0.0;
	for(index = 1; index < pTwoX; ++index)
	{
	  tTrig = cos1;
	  cos1 = tTrig * cos0 - sin1 * sin0;
	  sin1 = tTrig * sin0 + sin1 * cos0;
	  cos2 = (cos1 * cos1) - (sin1 * sin1);
	  sin2 = 2.0 * (cos1 * sin1);
	  tDp0 = data + (index * step);
	  tGp0 = data + ((pTwo1 - index) * step);
	  count = num;
	  do
	  {
	    tDp1 = tDp0 + tI1;
	    tGp1 = tGp0 + tI1;
	    tDp2 = tDp0 + tI2;
	    tGp2 = tGp0 + tI2;
	    tDp3 = tDp0 + tI3;
	    tGp3 = tGp0 + tI3;
	    tD9 = (sin2 * *tDp1) - (cos2 * *tGp1);
	    tD8 = (cos2 * *tDp1) + (sin2 * *tGp1);
	    tD0 = *tDp0;
	    tD1 = tD0 - tD8;
	    tD0 += tD8;
	    tD4 = *tGp0;
	    tD5 = tD4 - tD9;
	    tD4 += tD9;
	    tD9 = (sin2 * *tDp3) - (cos2 * *tGp3);
	    tD8 = (cos2 * *tDp3) + (sin2 * *tGp3);
	    tD2 = *tDp2;
	    tD3 = tD2 - tD8;
	    tD2 += tD8;
	    tD6 = *tGp2;
	    tD7 = tD6 - tD9;
	    tD6 += tD9;
	    tD9 = (sin1 * tD2) - (cos1 * tD7);
	    tD8 = (cos1 * tD2) + (sin1 * tD7);
	    *tDp2 = tD0 - tD8;
	    *tDp0 = tD0 + tD8;
	    *tGp3 = tD5 - tD9;
	    *tGp1 = tD5 + tD9;
	    tD9 = (cos1 * tD6) - (sin1 * tD3);
	    tD8 = (sin1 * tD6) + (cos1 * tD3);
	    *tGp2 = tD4 - tD8;
	    *tGp0 = tD4 + tD8;
	    *tDp3 = tD1 - tD9;
	    *tDp1 = tD1 + tD9;
	    tGp0 += tI4;
	    tDp0 += tI4;
	    count -= pTwo4;
	  }while(count > 0);
	}
      }while(pTwo4 < num);
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourHart1D FX\n"));
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup      AlgFourier
* \brief	Computes the Hartley transform of the given two
*		dimensional data, and does it in place.
* \param	data			Given data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
*/
AlgError 	AlgFourHart2D(double **data, int useBuf, int numX, int numY)
{
  double	*buf = NULL;
  AlgError	errNum = ALG_ERR_NONE;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourHart2D FE %p %d %d %d\n",
	   data, useBuf, numX, numY));
  if(useBuf && ((buf = (double*)AlcMalloc(numY * sizeof(double))) == NULL))
  {
    errNum = ALG_ERR_MALLOC;
  }
  else
  {
    int		idX,
		idY,
		idU,
		idV,
		halfX,
     		halfY;

    halfX = numX/2;
    halfY = numY/2;
    idU = numX - 1;
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idU,idV)
#endif
    for(idX = 1; idX < halfX; ++idX)
    {
      idV = numY - 1;
      for(idY = 1; idY < halfY; ++idY)
      {
	double	t[4];
	double	*p[4];

	p[0] = data[idY] + idX;
	p[1] = data[idY] + idU;
	p[2] = data[idV] + idX;
	p[3] = data[idV] + idU;
	t[0] = *p[0];
	t[1] = *p[1];
	t[2] = *p[2];
	t[3] = *p[3];
	*p[0] = ( t[0] + t[1] + t[2] - t[3]) * 0.5;
	*p[1] = ( t[0] + t[1] - t[2] + t[3]) * 0.5;
	*p[2] = ( t[0] - t[1] + t[2] + t[3]) * 0.5;
	*p[3] = (-t[0] + t[1] + t[2] + t[3]) * 0.5;
	--idV;
      }
      --idU;
    }
#ifdef _OPENMP
#pragma omp parallel for private(idY)
#endif
    for(idY = 0; idY < numY; ++idY)
    {
      AlgFourHart1D(data[idY], numX, 1);
    }
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY)
#endif
    for(idX = 0; idX < numX; ++idX)
    {
      if(useBuf)
      {
	for(idY = 0; idY < numY; ++idY)
	{
	  buf[idY] = data[idY][idX];
	}
	AlgFourHart1D(buf, numY, 1);
	for(idY = 0; idY < numY; ++idY)
	{
	  data[idY][idX] = buf[idY];
	}
      }
      else
      {
	AlgFourHart1D(data[0] + idY, numY, numX);
      }
    }
  }
  AlcFree(buf);
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourHart2D FX\n"));
  return(errNum);
}

/*!
* \return	void
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given one
*		dimensional complex data, and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n}\f$.
* \param	real			Given real data.
* \param	imag			Given imaginary data.
* \param	num			Number of data.
* \param	step			Offset in data elements between
*					the data to be transformed.
*/
void		AlgFour1D(double *real, double *imag, int num, int step)
{
  double	tD0,
		tD1,
		tD2,
		tD3,
		tD4;
  double	*tRp0,
		*tRp1,
		*tIp0,
		*tIp1;
  int		count;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour1D FE %p %p %d %d\n",
	   real, imag, num, step));
  tRp0 = real + step;
  tRp1 = real + ((num - 1) * step);
  tIp0 = imag + step;
  tIp1 = imag + ((num - 1) * step);
  count = (num / 2) - 1;
  while(count-- > 0)
  {
    tD1 = *tRp0;
    tD0 = *tRp1;
    tD2 = tD1 - tD0;
    tD1 += tD0;
    tD3 = *tIp0;
    tD0 = *tIp1;
    tD4 = tD3 - tD0;
    tD3 += tD0;
    *tRp0 = (tD1 + tD4) * 0.5;
    tRp0 += step;
    *tRp1 = (tD1 - tD4) * 0.5;
    tRp1 -= step;
    *tIp0 = (tD3 - tD2) * 0.5;
    tIp0 += step;
    *tIp1 = (tD3 + tD2) * 0.5;
    tIp1 -= step;
  }
#ifdef _OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef _OPENMP
#pragma omp section
#endif
    {
      AlgFourHart1D(real, num, step);
    }
#ifdef _OPENMP
#pragma omp section
#endif
    {
      AlgFourHart1D(imag, num, step);
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour1D FX\n"));
}

/*!
* \return	void
* \ingroup   	AlgFourier
* \brief	Computes the inverse Fourier transform of the given
*		complex one dimensional data, and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n}\f$.
* \param	real			Given real data.
* \param	imag			Given imaginary data.
* \param	num			Number of data.
* \param	step			Offset in data elements between
*					the data to be transformed.
*/
void		AlgFourInv1D(double *real, double *imag, int num, int step)
{
  double	tD0,
		tD1,
		tD2,
		tD3,
		tD4;
  double	*tRp0,
		*tRp1,
		*tIp0,
		*tIp1;
  int		count;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv1D FE %p %p %d %d\n",
	   real, imag, num, step));
#ifdef _OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef _OPENMP
#pragma omp section
#endif
    {
      AlgFourHart1D(real, num, step);
    }
#ifdef _OPENMP
#pragma omp section
#endif
    {
      AlgFourHart1D(imag, num, step);
    }
  }
  tRp0 = real + step;
  tRp1 = real + ((num - 1) * step);
  tIp0 = imag + step;
  tIp1 = imag + ((num - 1) * step);
  count = (num / 2) - 1;
  while(count-- > 0)
  {
    tD1 = *tRp0;
    tD0 = *tRp1;
    tD2 = tD1 - tD0;
    tD1 += tD0;

    tD3 = *tIp0;
    tD0 = *tIp1;
    tD4 = tD3 - tD0;
    tD3 += tD0;
    *tRp0 = (tD1 - tD4) * 0.5;
    tRp0 += step;
    *tRp1 = (tD1 + tD4) * 0.5;
    tRp1 -= step;
    *tIp0 = (tD3 + tD2) * 0.5;
    tIp0 += step;
    *tIp1 = (tD3 - tD2) * 0.5;
    tIp1 -= step;
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv1D FX\n"));
}

/*!
* \return	void
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given one
*		dimensional real data, and does it in place.
*	 	Data are returned in the array of size \f$N\f$ as:
*		The data are returned in the array of size N with the layout
*		as shown in the table (with 2M = N, r = real, i = imaginary):
 		|----------|
 		| r0\f$    |
 		| r1\f$    |
		| ..       |
 		| r(M-1)   |
		| rM       |
		| i1       |
		| i2       |
		| ...      |
		| i(M - 1) |
*		where the real and imaginary components are indexed as in
*		the arrays computed with AlgFour1D().
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n}\f$.
* \param	real			Given real data.
* \param	num			Number of data (N).
* \param	step			Offset in data elements between
*					the data to be transformed.
*/
void		AlgFourReal1D(double *real, int num, int step)
{
  double	tD0,
		tD1;
  double	*tRp0,
		*tRp1;
  int		count;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal1D FE %p %d %d\n",
	   real, num, step));
  tRp0 = real + step;
  tRp1 = real + ((num - 1) * step);
  count = num / 2;
  AlgFourHart1D(real, num, step);
  while(--count > 0)
  {
    tD0 = *tRp0;
    tD1 = *tRp1;
    *tRp0 = (tD0 + tD1) * 0.5;
    *tRp1 = (tD0 - tD1) * 0.5;
    tRp0 += step;
    tRp1 -= step;
  }
  count = (num / 2);
  tRp0 = real + ((count + 1) * step);
  tRp1 = real + ((num - 1) * step);
  while(count > 0)
  {
    tD0 = -(*tRp0);
    tD1 = -(*tRp1);
    *tRp0 = tD1;
    *tRp1 = tD0;
    tRp0 += step;
    tRp1 -= step;
    count -= 2;
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal1D FX\n"));
}

/*!
* \return	void
* \ingroup   	AlgFourier
* \brief	Computes the inverse Fourier transform of the given one
*		one dimensional real data, and does it in place.
*		The data should layed out in the array as returned by
*		AlgFourReal1D().
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n}\f$.
* \param	real			Given real/complex data.
* \param	num			Number of data.
* \param	step			Offset in data elements between
*					the data to be transformed.
*/
void		AlgFourRealInv1D(double *real, int num, int step)
{
  double	tD0,
		tD1;
  double	*tRp0,
		*tRp1;
  int		count;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv1D FE %p %d %d\n",
	   real, num, step));
  count = (num / 2);
  tRp0 = real + ((count + 1) * step);
  tRp1 = real + ((num - 1) * step);
  while(count > 0)
  {
    tD0 = -(*tRp0);
    tD1 = -(*tRp1);
    *tRp0 = tD1;
    *tRp1 = tD0;
    tRp0 += step;
    tRp1 -= step;
    count -= 2;
  }
  tRp0 = real + step;
  tRp1 = real + ((num - 1) * step);
  count = num / 2;
  while(--count > 0)
  {
    tD0 = *tRp0;
    tD1 = *tRp1;
    *tRp0 = (tD0 + tD1);
    *tRp1 = (tD0 - tD1);
    tRp0 += step;
    tRp1 -= step;
  }
  AlgFourHart1D(real, num, step);
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv1D FX\n"));
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given two
*		dimensional complex data, and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y}\f$.
* \param	real			Given real data.
* \param	imag			Given imaginary data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
*/
AlgError	AlgFour2D(double **real, double **imag,
			  int useBuf, int numX, int numY)
{
  AlgError	errNum;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour2D FE %p %p %d %d %d\n",
	   real, imag, useBuf, numX, numY));
  errNum = AlgFourRepXY1D(real, imag, ALG_FOUR_AXIS_X, useBuf,
                          numX, numY, ALG_FOUR_DIR_FWD);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXY1D(real, imag, ALG_FOUR_AXIS_Y, useBuf,
                            numX, numY, ALG_FOUR_DIR_FWD);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour2D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the inverse Fourier transform of the given two
*		dimensional complex data, and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y}\f$.
* \param	real			Given real data.
* \param	imag			Given imaginary data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
*/
AlgError	AlgFourInv2D(double **real, double **imag,
			     int useBuf, int numX, int numY)
{
  AlgError	errNum;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv2D FE %p %p %d %d %d\n",
	   real, imag, useBuf, numX, numY));
  errNum = AlgFourRepXY1D(real, imag, ALG_FOUR_AXIS_Y, useBuf,
                          numX, numY, ALG_FOUR_DIR_INV);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXY1D(real, imag, ALG_FOUR_AXIS_X, useBuf,
                            numX, numY, ALG_FOUR_DIR_INV);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv2D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given two
*		dimensional real data, and does it in place.
*
*		The layout of the array is similar to that in AlgFourReal1D()
*		and as shown in the table (with 2M = N, r = real,
*		i = imaginary):
*
	-------|--------|---|------------|-------|--------|---|-----------
	r00    |r01     |...|r0(M-1)     |r0M    |i01     |...|i0(M-1)
	...    |...     |...|...         |...    |...     |...|...
	r(M-1)0|r(M-1)1 |...|r(M-1)(M-1) |r(M-1)M|i(M-1)1 |...|i(M-1)(M-1)
	rM0    |rM1     |...|rM(M-1)     |rMM    |iM1     |...|iM(M-1)
	i10    |r(M+1)1 |...|r(M+1)(M-1) |i1M    |i(M+1)1 |...|i(M+1)(M-1)
	...    |...     |...|...         |...    |...     |...|...
	i(M-1)0|r(2M-1)1|...|r(2M-1)(M-1)|i(M-1)M|i(2M-1)1|...|i(2M-1)(M-1)

*		Using contiguous buffers has a large effect for
*		data larger than a CPU's fastest cache and little
*		cost for smaller data arrays.
*		The times below were measured for on a Lenovo T430s
*		with a i7-2640M with 4kB cache. Data sizes are the
*		number of double values. These times include buffer
*		allocation but not array allocation or file I/O
*		and are the mean of 5 runs using elapsed (wall
*		clock) time.
*	data sz | tm (no buf, 1 thr) | tm (buf, 1 thr) | tm (buf, 2 thr)
 	--------|--------------------|-----------------|-------------
          256^2 |     3ms            |    3ms          |    2ms
          512^2 |    12ms            |    6ms          |    4ms
         1024^2 |   120ms            |   44ms          |   35ms
         2048^2 |   600ms            |  170ms          |  110ms
         4096^2 |  3300ms            |  764ms          |  460ms
         8192^2 | 17000ms            | 3600ms          | 2100ms
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y}\f$.
* \param	real			Given real data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
*/
AlgError	AlgFourReal2D(double **real,
			      int useBuf, int numX, int numY)
{
  AlgError	errNum;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal2D FE %p %d %d %d\n",
	   real, useBuf, numX, numY));
  errNum = AlgFourRepXYReal1D(real, ALG_FOUR_AXIS_X, useBuf, numX, numY,
                              ALG_FOUR_DIR_FWD);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYReal1D(real, ALG_FOUR_AXIS_Y, useBuf, numX, numY,
    				ALG_FOUR_DIR_FWD);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal2D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given two
*		dimensional data which resulted from a transform using
*		AlgFourReal2D(), and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y}\f$.
* \param	real			Given real/complex data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
*/
AlgError	AlgFourRealInv2D(double **real,
				 int useBuf, int numX, int numY)
{
  AlgError	errNum;
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv2D FE %p %d %d %d\n",
	   real, useBuf, numX, numY));
  errNum = AlgFourRepXYReal1D(real, ALG_FOUR_AXIS_Y, useBuf, numX, numY,
                              ALG_FOUR_DIR_INV);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYReal1D(real, ALG_FOUR_AXIS_X, useBuf, numX, numY,
  				ALG_FOUR_DIR_INV);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv2D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given three
*		dimensional complex data, and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y} \sqrt{n_z}\f$.
*
*		Using contiguous buffers has a large effect for
*		data larger than a CPU's fastest cache and little
*		cost for smaller data arrays.
* \param	real			Given real data.
* \param	imag			Given imaginary data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
* \param	numZ			Number of data in each plane.
*/
AlgError	AlgFour3D(double ***real, double ***imag,
			  int useBuf, int numX, int numY, int numZ)
{
  AlgError	errNum;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour3D FE %p %p %d %d %d %d\n",
	   real, imag, useBuf, numX, numY, numZ));
  errNum = AlgFourRepXYZ1D(real, imag, ALG_FOUR_AXIS_X, useBuf,
                           numX, numY, numZ, ALG_FOUR_DIR_FWD);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZ1D(real, imag, ALG_FOUR_AXIS_Y, useBuf,
                             numX, numY, numZ, ALG_FOUR_DIR_FWD);
  }
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZ1D(real, imag, ALG_FOUR_AXIS_Z, useBuf,
                             numX, numY, numZ, ALG_FOUR_DIR_FWD);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour3D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the inverse Fourier transform of the given three
*		dimensional complex data, and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y} \sqrt{n_z}\f$.
* \param	real			Given real data.
* \param	imag			Given imaginary data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
* \param	numZ			Number of data in each plane.
*/
AlgError	AlgFourInv3D(double ***real, double ***imag,
			     int useBuf, int numX, int numY, int numZ)
{
  AlgError	errNum;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv3D FE %p %p %d %d %d %d\n",
	   real, imag, useBuf, numX, numY, numZ));
  errNum = AlgFourRepXYZ1D(real, imag, ALG_FOUR_AXIS_Z, useBuf,
                          numX, numY, numZ, ALG_FOUR_DIR_INV);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZ1D(real, imag, ALG_FOUR_AXIS_Y, useBuf,
                            numX, numY, numZ, ALG_FOUR_DIR_INV);
  }
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZ1D(real, imag, ALG_FOUR_AXIS_X, useBuf,
                            numX, numY, numZ, ALG_FOUR_DIR_INV);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv3D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given three
*		dimensional real data, and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y} \sqrt{n_z}\f$.
*
*		Using contiguous buffers has a large effect for
*		data larger than a CPU's fastest cache and little
* \param	real			Given real data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
* \param	numZ			Number of data in each plane.
*/
AlgError	AlgFourReal3D(double ***real,
			      int useBuf, int numX, int numY, int numZ)
{
  AlgError	errNum;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal3D FE %p %d %d %d %d\n",
	   real, useBuf, numX, numY, numZ));
  errNum = AlgFourRepXYZReal1D(real, ALG_FOUR_AXIS_X, useBuf,
                               numX, numY, numZ, ALG_FOUR_DIR_FWD);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZReal1D(real, ALG_FOUR_AXIS_Y, useBuf,
                                 numX, numY, numZ, ALG_FOUR_DIR_FWD);
  }
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZReal1D(real, ALG_FOUR_AXIS_Z, useBuf,
    				 numX, numY, numZ, ALG_FOUR_DIR_FWD);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal3D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \ingroup   	AlgFourier
* \brief	Computes the Fourier transform of the given three
*		dimensional data which resulted from a transform using
*		AlgFourReal3D(), and does it in place.
*		The transformed values data are scaled by a factor
*		of \f$\sqrt{n_x} \sqrt{n_y} \sqrt{n_z}\f$.
* \param	real			Given real/complex data.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data in each row.
* \param	numY			Number of data in each column.
* \param	numZ			Number of data in each plane.
*/
AlgError	AlgFourRealInv3D(double ***real,
				 int useBuf, int numX, int numY, int numZ)
{
  AlgError	errNum;
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv3D FE %p %d %d %d %d\n",
	   real, useBuf, numX, numY, numZ));
  errNum = AlgFourRepXYZReal1D(real, ALG_FOUR_AXIS_Z, useBuf,
                               numX, numY, numZ, ALG_FOUR_DIR_INV);
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZReal1D(real, ALG_FOUR_AXIS_Y, useBuf,
                                 numX, numY, numZ, ALG_FOUR_DIR_INV);
  }
  if(errNum == ALG_ERR_NONE)
  {
    errNum = AlgFourRepXYZReal1D(real, ALG_FOUR_AXIS_X, useBuf,
                                 numX, numY, numZ, ALG_FOUR_DIR_INV);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv3D FX\n"));
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \brief	Computes repeated Fourier transforms of a 1D complex
*	 	array within a complex array 2D.
* \param	real			Given real 2D data array.
* \param	imag			Given imaginary 2D data array.
* \param	axis			Axis for partial evaluation.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data columns in given
*					data array.
* \param	numY			Number of data rows in given
*					data array.
* \param	dir			Forward or inverse transform.
*/
static AlgError	AlgFourRepXY1D(double **real, double **imag,
			       AlgFourAxis axis, int useBuf,
			       int numX, int numY, AlgFourDir dir)
{
  AlgError	errNum = ALG_ERR_NONE;

  if(axis == ALG_FOUR_AXIS_X)                              /* Transform rows */
  {
    if(dir == ALG_FOUR_DIR_FWD)
    {
      int	idY;

#ifdef _OPENMP
#pragma omp parallel for private(idY)
#endif
      for(idY = 0; idY < numY; ++idY)
      {
	AlgFour1D(*(real + idY), *(imag + idY), numX, 1);
      }
    }
    else
    {
      int	idY;

#ifdef _OPENMP
#pragma omp parallel for private(idY)
#endif
      for(idY = 0; idY < numY; ++idY)
      {
	AlgFourInv1D(*(real + idY), *(imag + idY), numX, 1);
      }
    }
  }
  else /* axis == ALG_FOUR_AXIS_Y */                    /* Transform columns */
  {
    if(useBuf)	                           /* Use column buffers if provided */
    {
      int	nThr = 1;
      double	*bufBase = NULL;

#ifdef _OPENMP
#pragma omp parallel
      {
#pragma omp master
        {
	  nThr = omp_get_num_threads();
	}
      }
#endif
      if((bufBase = AlcMalloc(sizeof(double) * numY * 2 * nThr)) == NULL)
      {
	errNum = ALG_ERR_MALLOC;
      }
      if(errNum == ALG_ERR_NONE)
      {
        int	idX;

#ifdef _OPENMP
#pragma omp parallel for private(idX) num_threads(nThr)
#endif
	for(idX = 0; idX < numX; ++idX)
	{
	  int	  idY,
	  	  thrId = 0;
	  double  *reBuf,
		  *imBuf;

#ifdef _OPENMP
	  thrId = omp_get_thread_num();
#endif
	  reBuf = bufBase + (2 * numY * thrId);
	  imBuf = reBuf + numY;
	  /* Copy to buffer. */
	  for(idY = 0; idY < numY; ++idY)
	  {
	    *(reBuf + idY) = *(*(real + idY) + idX);
	    *(imBuf + idY) = *(*(imag + idY) + idX);
	  }
	  /* Transform buffer. */
	  if(dir == ALG_FOUR_DIR_FWD)
	  {
	    AlgFour1D(reBuf, imBuf, numY, 1);
	  }
	  else
	  {
	    AlgFourInv1D(reBuf, imBuf, numY, 1);
	  }
	  /* Copy back from buffer. */
	  for(idY = 0; idY < numY; ++idY)
	  {
	    *(*(real + idY) + idX) = *(reBuf + idY);
	    *(*(imag + idY) + idX) = *(imBuf + idY);
	  }
	}
	AlcFree(bufBase);
      }
    }
    else              /* If no column buffers transform the columns in place */
    {
      int	idX;

      if(dir == ALG_FOUR_DIR_FWD)
      {
#ifdef _OPENMP
#pragma omp parallel for private(idX)
#endif
	for(idX = 0; idX < numX; ++idX)
	{
	  AlgFour1D(*real + idX, *imag + idX, numY, numX);
        }
      }
      else
      {
#ifdef _OPENMP
#pragma omp parallel for private(idX)
#endif
	for(idX = 0; idX < numX; ++idX)
	{
	  AlgFourInv1D(*real + idX, *imag + idX, numY, numX);
        }
      }
    }
  }
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \brief	Computes repeated the Fourier transforms of the given
*	 	one dimensional real data array.
* \brief	Computes repeated Fourier transforms of a 1D real/
*	 	complex conjugate array within a 2D real/complex conjugate
*	 	array.
* \param	data			Given 2D data array.
* \param	axis			Axis for partial evaluation.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data columns in given
*					data array.
* \param	numY			Number of data rows in given
*					data array.
* \param	dir			Forward or inverse transform.
*/
static AlgError	AlgFourRepXYReal1D(double **data, AlgFourAxis axis, int useBuf,
				   int numX, int numY, AlgFourDir dir)
{
  AlgError	errNum = ALG_ERR_NONE;

  if(axis == ALG_FOUR_AXIS_X)			 	   /* Transform rows */
  {
    if(dir == ALG_FOUR_DIR_FWD)
    {
      int 	idY;

#ifdef _OPENMP
#pragma omp parallel for private(idY)
#endif
      for(idY = 0; idY < numY; ++idY)
      {
	AlgFourReal1D(*(data + idY), numX, 1);
      }
    }
    else
    {
      int	idY;

#ifdef _OPENMP
#pragma omp parallel for private(idY)
#endif
      for(idY = 0; idY < numY; ++idY)
      {
	AlgFourRealInv1D(*(data + idY), numX, 1);
      }
    }
  }
  else /* axis == ALG_FOUR_AXIS_Y */			/* Transform columns */
  {
    int		halfData;

    halfData = numX / 2;
    if(useBuf)
    {
      int	nThr = 1;
      double	*bufBase = NULL;

#ifdef _OPENMP
#pragma omp parallel
      {
#pragma omp master
        {
          nThr = omp_get_num_threads();
	}
      }
#endif
      if((bufBase = AlcMalloc(sizeof(double) * numY * 2 * nThr)) == NULL)
      {
	errNum = ALG_ERR_MALLOC;
      }
      if(errNum == ALG_ERR_NONE)
      {
#ifdef _OPENMP
#pragma omp parallel sections
#endif
	{
#ifdef _OPENMP
#pragma omp section
#endif
	  {
	    int	    idY,
		    thrId = 0;
	    double  *reBuf;

#ifdef _OPENMP
	    thrId = omp_get_thread_num();
#endif
	    reBuf = bufBase + (numY * thrId);
	    for(idY = 0; idY < numY; ++idY)
	    {
	      *(reBuf + idY) = **(data + idY);
	    }
	    if(dir == ALG_FOUR_DIR_FWD)
	    {
	      AlgFourReal1D(reBuf, numY, 1);
	    }
	    else
	    {
	      AlgFourRealInv1D(reBuf, numY, 1);
	    }
	    for(idY = 0; idY < numY; ++idY)
	    {
	      **(data + idY) = *(reBuf + idY);
	    }
	  }
#ifdef _OPENMP
#pragma omp section
#endif
	  {
	    int	    idY,
		    thrId = 0;
	    double  *reBuf;

#ifdef _OPENMP
	    thrId = omp_get_thread_num();
#endif
	    reBuf = bufBase + (numY * thrId);
	    for(idY = 0; idY < numY; ++idY)
	    {
	      *(reBuf + idY) = *(*(data + idY) + halfData);
	    }
	    if(dir == ALG_FOUR_DIR_FWD)
	    {
	      AlgFourReal1D(reBuf, numY, 1);
	    }
	    else
	    {
	      AlgFourRealInv1D(reBuf, numY, 1);
	    }
	    for(idY = 0; idY < numY; ++idY)
	    {
	      *(*(data + idY) + halfData) = *(reBuf + idY);
	    }
	  }
	}
	{
	  int	idX;

#ifdef _OPENMP
#pragma omp parallel for private(idX)
#endif
	  for(idX = 1; idX < halfData; ++idX)
	  {
	    int	   idY,
		   thrId = 0;
	    double *reBuf,
		   *imBuf;

#ifdef _OPENMP
	    thrId = omp_get_thread_num();
#endif
	    reBuf = bufBase + (2 * numY * thrId);
	    imBuf = reBuf + numY;
	    for(idY = 0; idY < numY; ++idY)
	    {
	      *(reBuf + idY) = *(*(data + idY) + idX);
	      *(imBuf + idY) = *(*(data + idY) + halfData + idX);
	    }
	    if(dir == ALG_FOUR_DIR_FWD)
	    {
	      AlgFour1D(reBuf, imBuf, numY, 1);
	    }
	    else
	    {
	      AlgFourInv1D(reBuf, imBuf, numY, 1);
	    }
	    for(idY = 0; idY < numY; ++idY)
	    {
	      *(*(data + idY) + idX) = *(reBuf + idY);
	      *(*(data + idY) + halfData + idX) = *(imBuf + idY);
	    }
	  }
	}
      }
    }
    else
    {
      /* If no column buffers then just transform the columns in place this
       * is done by: transforming column 0, transform column numX/2 and then
       * transforming the remaining columns. */
      if(dir == ALG_FOUR_DIR_FWD)
      {
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	  int	idX;
#ifdef _OPENMP
#pragma omp sections
#endif
	  {
#ifdef _OPENMP
#pragma omp section
#endif
	    {
	      AlgFourReal1D(*data, numY, numX);
	    }
#ifdef _OPENMP
#pragma omp section
#endif
	    {
	      AlgFourReal1D(*data + halfData, numY, numX);
	    }
	  }
#ifdef _OPENMP
#pragma omp for private(idX)
#endif
	  for(idX = 1; idX < halfData; ++idX)
	  {
	    AlgFour1D(*data + idX, *data + halfData + idX, numY, numX);
	  }
	}
      }
      else
      {
        int 	idX;

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
#ifdef _OPENMP
#pragma omp sections
#endif
	  {
#ifdef _OPENMP
#pragma omp section
#endif
	    {
	      AlgFourRealInv1D(*data, numY, numX);
	    }
#ifdef _OPENMP
#pragma omp section
#endif
	    {
	      AlgFourRealInv1D(*data + halfData, numY, numX);
	    }
	  }
#ifdef _OPENMP
#pragma omp for private(idX)
#endif
	  for(idX = 1; idX < halfData; ++idX)
	  {
	    AlgFourInv1D(*data + idX, *data + halfData + idX, numY,
		numX);
	  }
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \brief	Computes repeated Fourier transforms of a 1D complex
*	 	array within a complex array 3D.
* \param	real			Given 3D real data array.
* \param	imag			Given 3D imaginary data array.
* \param	axis			Axis for partial evaluation.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data columns in given
*					data array.
* \param	numY			Number of data rows in given
*					data array.
* \param	numZ			Number of data planes in given
*					data array.
* \param	dir			Forward or inverse transform.
*/
static AlgError	AlgFourRepXYZ1D(double ***real, double ***imag,
			        AlgFourAxis axis, int useBuf,
			        int numX, int numY, int numZ, AlgFourDir dir)
{
  int		idX,
  		idY,
		idZ;
  AlgError	errNum = ALG_ERR_NONE;

  switch(axis)
  {
    case ALG_FOUR_AXIS_X:
      /* Transform rows */
#ifdef _OPENMP
#pragma omp parallel for private(idY,idZ)
#endif
      for(idZ = 0; idZ < numZ; ++idZ)
      {
	for(idY = 0; idY < numY; ++idY)
	{
	  if(dir == ALG_FOUR_DIR_FWD)
	  {
	    AlgFour1D(*(*(real + idZ) + idY), *(*(imag + idZ) + idY),
		numX, 1);
	  }
	  else
	  {
	    AlgFourInv1D(*(*(real + idZ) + idY), *(*(imag + idZ) + idY),
		numX, 1);
	  }
	}
      }
      break;
    case ALG_FOUR_AXIS_Y:
      /* Transform columns */
      if(useBuf)
      {
	int	nThr = 1;
	double	*bufBase = NULL;

#ifdef _OPENMP
#pragma omp parallel
	{
#pragma omp master
	  {
	    nThr = omp_get_num_threads();
	  }
	}
#endif
	if((bufBase = AlcMalloc(sizeof(double) * numY * 2 * nThr)) == NULL)
	{
	  errNum = ALG_ERR_MALLOC;
	}
	if(errNum == ALG_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ) num_threads(nThr)
#endif
	  for(idZ = 0; idZ < numZ; ++idZ)
	  {
	    int	    thrId = 0;
	    double  *reBuf,
		    *imBuf;

#ifdef _OPENMP
	    thrId = omp_get_thread_num();
#endif
	    reBuf = bufBase + (2 * numY * thrId);
	    imBuf = reBuf + numY;
	    for(idX = 0; idX < numX; ++idX)
	    {
	      /* Copy to buffer. */
	      for(idY = 0; idY < numY; ++idY)
	      {
		*(reBuf + idY) = *(*(*(real + idZ) + idY) + idX);
		*(imBuf + idY) = *(*(*(imag + idZ) + idY) + idX);
	      }
	      /* Transform buffer. */
	      if(dir == ALG_FOUR_DIR_FWD)
	      {
		AlgFour1D(reBuf, imBuf, numY, 1);
	      }
	      else
	      {
		AlgFourInv1D(reBuf, imBuf, numY, 1);
	      }
	      /* Copy back from buffer. */
	      for(idY = 0; idY < numY; ++idY)
	      {
		*(*(*(real + idZ) + idY) + idX) = *(reBuf + idY);
		*(*(*(imag + idZ) + idY) + idX) = *(imBuf + idY);
	      }
	    }
	  }
	  AlcFree(bufBase);
	}
      }
      else
      {
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ)
#endif
	for(idZ = 0; idZ < numZ; ++idZ)
	{
	  for(idX = 0; idX < numX; ++idX)
	  {
	    if(dir == ALG_FOUR_DIR_FWD)
	    {
	      AlgFour1D(**(real + idZ) + idX, **(imag + idZ) + idX,
		  numY, numX);
	    }
	    else
	    {
	      AlgFourInv1D(**(real + idZ) + idX, **(imag + idZ) + idX,
		  numY, numX);
	    }
	  }
	}
      }
      break;
    case ALG_FOUR_AXIS_Z:
      /* Transform columns */
      if(useBuf)
      {
	int	nThr = 1;
	double	*bufBase = NULL;

#ifdef _OPENMP
#pragma omp parallel
	{
#pragma omp master
	  {
	    nThr = omp_get_num_threads();
	  }
	}
#endif
	if((bufBase = AlcMalloc(sizeof(double) * numZ * 2 * nThr)) == NULL)
	{
	  errNum = ALG_ERR_MALLOC;
	}
	if(errNum == ALG_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ) num_threads(nThr)
#endif
	  for(idY = 0; idY < numY; ++idY)
	  {
	    int	    thrId = 0;
	    double  *reBuf,
		    *imBuf;

#ifdef _OPENMP
	    thrId = omp_get_thread_num();
#endif
	    reBuf = bufBase + (2 * numZ * thrId);
	    imBuf = reBuf + numZ;
	    for(idX = 0; idX < numX; ++idX)
	    {
	      /* Copy to buffer. */
	      for(idZ = 0; idZ < numZ; ++idZ)
	      {
		*(reBuf + idZ) = *(*(*(real + idZ) + idY) + idX);
		*(imBuf + idZ) = *(*(*(imag + idZ) + idY) + idX);
	      }
	      /* Transform buffer. */
	      if(dir == ALG_FOUR_DIR_FWD)
	      {
		AlgFour1D(reBuf, imBuf, numZ, 1);
	      }
	      else
	      {
		AlgFourInv1D(reBuf, imBuf, numZ, 1);
	      }
	      /* Copy back from buffer. */
	      for(idZ = 0; idZ < numZ; ++idZ)
	      {
		*(*(*(real + idZ) + idY) + idX) = *(reBuf + idZ);
		*(*(*(imag + idZ) + idY) + idX) = *(imBuf + idZ);
	      }
	    }
	  }
	  AlcFree(bufBase);
	}
      }
      else
      {
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY)
#endif
	for(idY = 0; idY < numY; ++idY)
	{
	  for(idX = 0; idX < numX; ++idX)
	  {
	    if(dir == ALG_FOUR_DIR_FWD)
	    {
	      AlgFour1D(*(*real + idY) + idX, *(*imag + idY) + idX,
		  numZ, numX * numY);
	    }
	    else
	    {
	      AlgFourInv1D(*(*real + idY) + idX, *(*imag + idY) + idX,
		  numZ, numX * numY);
	    }
	  }
	}
      }
      break;
    default:
      break;
  }
  return(errNum);
}

/*!
* \return	Error code, may be set if buffers can not be allocated.
* \brief	Computes repeated the Fourier transforms of the given
*	 	one dimensional real data array.
* \brief	Computes repeated Fourier transforms of a 1D real/
*	 	complex conjugate array within a 3D real/complex conjugate
*	 	array.
* \param	data			Given 2D data array.
* \param	axis			Axis for partial evaluation.
* \param	useBuf			Allocate private buffers to make
* 					columns contiguous.
* \param	numX			Number of data columns in given
*					data array.
* \param	numY			Number of data rows in given
*					data array.
* \param	numZ			Number of data planes in given
*					data array.
* \param	dir			Forward or inverse transform.
*/
static AlgError	AlgFourRepXYZReal1D(double ***data, AlgFourAxis axis,
				    int useBuf,
				    int numX, int numY, int numZ,
				    AlgFourDir dir)
{
  int		idX,
		idY,
		idZ,
		halfData;
  AlgError	errNum = ALG_ERR_NONE;

  switch(axis)
  {
    case  ALG_FOUR_AXIS_X:
      /* Transform rows */
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ)
#endif
      for(idZ = 0; idZ < numZ; ++idZ)
      {
	for(idY = 0; idY < numY; ++idY)
	{
	  if(dir == ALG_FOUR_DIR_FWD)
	  {
	    AlgFourReal1D(*(*(data + idZ) + idY), numX, 1);
	  }
	  else
	  {
	    AlgFourRealInv1D(*(*(data + idZ) + idY), numX, 1);
	  }
	}
      }
      break;
    case  ALG_FOUR_AXIS_Y:
      /* Transform columns */
      halfData = numY / 2;
      if(useBuf)
      {
	int	nThr = 1;
	double	*bufBase = NULL;

#ifdef _OPENMP
#pragma omp parallel
	{
#pragma omp master
	  {
	    nThr = omp_get_num_threads();
	  }
	}
#endif
	if((bufBase = AlcMalloc(sizeof(double) * numY * 2 * nThr)) == NULL)
	{
	  errNum = ALG_ERR_MALLOC;
	}
	if(errNum == ALG_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ)
#endif
	  for(idZ = 0; idZ < numZ; ++idZ)
	  {
	    int	thrId = 0;
	    double  *reBuf,
	    	    *imBuf;
#ifdef _OPENMP
	    thrId = omp_get_thread_num();
#endif
	    reBuf = bufBase + (2 * numY * thrId);
	    imBuf = reBuf + numY;
	    /* Copy to buffer. */
	    for(idY = 0; idY < numY; ++idY)
	    {
	      *(reBuf + idY) = **(*(data + idZ) + idY);
	      *(imBuf + idY) = *(*(*(data + idZ) + idY) + halfData);
	    }
	    if(dir == ALG_FOUR_DIR_FWD)
	    {
	      AlgFourReal1D(reBuf, numY, 1);
	      AlgFourReal1D(imBuf, numY, 1);
	    }
	    else
	    {
	      AlgFourRealInv1D(reBuf, numY, 1);
	      AlgFourRealInv1D(imBuf, numY, 1);
	    }
	    /* Copy back from buffer again. */
	    for(idY = 0; idY < numY; ++idY)
	    {
	      **(*(data + idZ) + idY) = *(reBuf + idY);
	      *(*(*(data + idZ) + idY) + halfData) = *(imBuf + idY);
	    }
	    for(idX = 1; idX < halfData; ++idX)
	    {
	      /* Copy to buffer. */
	      for(idY = 0; idY < numY; ++idY)
	      {
		*(reBuf + idY) = *(*(*(data + idZ) + idY) + idX);
		*(imBuf + idY) = *(*(*(data + idZ) + idY) + halfData + idX);
	      }
	      if(dir == ALG_FOUR_DIR_FWD)
	      {
		AlgFour1D(reBuf, imBuf, numY, 1);
	      }
	      else
	      {
		AlgFourInv1D(reBuf, imBuf, numY, 1);
	      }
	      /* Copy back. */
	      for(idY = 0; idY < numY; ++idY)
	      {
		*(*(*(data + idZ) + idY) + idX) = *(reBuf + idY);
		*(*(*(data + idZ) + idY) + halfData + idX) = *(imBuf + idY);
	      }
	    }
	  }
	  AlcFree(bufBase);
        }
      }
      else
      {
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ)
#endif
	for(idZ = 0; idZ < numZ; ++idZ)
	{
	  if(dir == ALG_FOUR_DIR_FWD)
	  {
	    AlgFourReal1D(**(data + idZ), numY, numX);
	    AlgFourReal1D(**(data + idZ) + halfData, numY, numX);
	    for(idX = 1; idX < halfData; ++idX)
	    {
	      AlgFour1D(**(data + idZ) + idX,
		  **(data + idZ) + halfData + idX, numY, numX);
	    }
	  }
	  else
	  {
	    AlgFourRealInv1D(**(data + idZ), numY, numX);
	    AlgFourRealInv1D(**(data + idZ) + halfData, numY, numX);
	    for(idX = 1; idX < halfData; ++idX)
	    {
	      AlgFourInv1D(**(data + idZ) + idX,
		  **(data + idZ) + halfData + idX, numY, numX);
	    }
	  }
	}
      }
      break;
    case  ALG_FOUR_AXIS_Z:
      /* Transform planes */
      halfData = numX / 2;
      if(useBuf)
      {
	int	nThr = 1;
	double	*bufBase = NULL;

#ifdef _OPENMP
#pragma omp parallel
	{
#pragma omp master
	  {
	    nThr = omp_get_num_threads();
	  }
	}
#endif
	if((bufBase = AlcMalloc(sizeof(double) * numZ * 2 * nThr)) == NULL)
	{
	  errNum = ALG_ERR_MALLOC;
	}
	if(errNum == ALG_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ)
#endif
	  for(idY = 0; idY < numY; ++idY)
	  {
	    int	thrId = 0;
	    double  *reBuf,
	    	    *imBuf;
#ifdef _OPENMP
	    thrId = omp_get_thread_num();
#endif
	    reBuf = bufBase + (2 * numZ * thrId);
	    imBuf = reBuf + numZ;
	    /* Copy to buffer. */
	    for(idZ = 0; idZ < numZ; ++idZ)
	    {
	      *(reBuf + idZ) = **(*(data + idZ) + idY);
	      *(imBuf + idZ) = *(*(*(data + idZ) + idY) + halfData);
	    }
	    if(dir == ALG_FOUR_DIR_FWD)
	    {
	      AlgFourReal1D(reBuf, numZ, 1);
	      AlgFourReal1D(imBuf, numZ, 1);
	    }
	    else
	    {
	      AlgFourRealInv1D(reBuf, numZ, 1);
	      AlgFourRealInv1D(imBuf, numZ, 1);
	    }
	    /* Copy back from buffer. */
	    for(idZ = 0; idZ < numZ; ++idZ)
	    {
	      **(*(data + idZ) + idY) = *(reBuf + idZ);
	      *(*(*(data + idZ) + idY) + halfData) = *(imBuf + idZ);
	    }
	    for(idX = 1; idX < halfData; ++idX)
	    {
	      /* Copy to buffer. */
	      for(idZ = 0; idZ < numZ; ++idZ)
	      {
	        *(reBuf + idZ) = *(*(*(data + idZ) + idY) + idX);
	        *(imBuf + idZ) = *(*(*(data + idZ) + idY) + halfData + idX);
	      }
	      if(dir == ALG_FOUR_DIR_FWD)
	      {
		AlgFour1D(reBuf, imBuf, numZ, 1);
	      }
	      else
	      {
		AlgFourInv1D(reBuf, imBuf, numZ, 1);
	      }
	      /* Copy back from buffer. */
	      for(idZ = 0; idZ < numZ; ++idZ)
	      {
	        *(*(*(data + idZ) + idY) + idX) = *(reBuf + idZ);
		*(*(*(data + idZ) + idY) + halfData + idX) = *(imBuf + idZ);
	      }
	    }

	  }
	  AlcFree(bufBase);
	}
      }
      else
      {
#ifdef _OPENMP
#pragma omp parallel for private(idX,idY,idZ)
#endif
	for(idY = 0; idY < numY; ++idY)
	{
	  if(dir == ALG_FOUR_DIR_FWD)
	  {
	    AlgFourReal1D(*(*data + idY), numZ, numX * numY);
	    AlgFourReal1D(*(*data + idY) + halfData, numZ, numX * numY);
	    for(idX = 1; idX < halfData; ++idX)
	    {
	      AlgFour1D(*(*data + idY) + idX,
		  *(*data + idY) + halfData + idX, numZ, numX * numY);
	    }
	  }
	  else
	  {
	    AlgFourRealInv1D(*(*data + idY), numZ, numX * numY);
	    AlgFourRealInv1D(*(*data + idY) + halfData, numZ, numX * numY);
	    for(idX = 1; idX < halfData; ++idX)
	    {
	      AlgFourInv1D(*(*data + idY) + idX,
		  *(*data + idY) + halfData + idX, numZ, numX * numY);
	    }
	  }
	}
      }
      break;
    default:
      break;
  }
  return(errNum);
}
