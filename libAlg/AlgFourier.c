#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas
* Title:        AlgFourier.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Fast Fourier and Hartley transform functions for the
*	  	MRC Human Genetics Unit numerical algorithm library.
*		The history of this software is (as stated by Mayer):
*		  Euler		- Probable inventor of the fourier
*				  transform.
*		  Gauss		- Probable inventor of the FFT.
*		  Hartley	- Probable inventor of the hartley
*				  transform.
*		  Bracewell & Buneman - Patent holders for FHT!
*		  Ron Mayer	- Produced FFT code using the FHT.
*		  Bill Hill	- Hacked Ron Mayer's code to simplify
*				  multi dimensional FFT's and for
*				  compatability with our existing FFT
*				  routines here at MRC HGU.
*				  Added multithreading for increased
*				  speed on multi-cpu system VR4 (Sun
*				  Solaris 2.X) machines.
* 		The two dimensional transform routines may be supplied
*		with buffers sufficient to hold a column of complex
*		data, or NULLS may be given. If supplied then these
*		buffers will be used to hold each column of data during
*		the transforms of the columns. When the size of the
*		data being transformed is comparable with or greater
*		than the size of the data cache of the target machine,
*		then it is far more efficient (a factor of 10 faster on
*		a Sun 10/411) to supply these buffers.
*		All code which depends on threads and can't be ignored
*		is controlled by ALC_THREADS_USED. All functions pass
*		the number of concurrent threads available as a
*		parameter, this can be 0 or 1 if multi-threading is not
*		required.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.
************************************************************************/
#include <Alg.h>
#include <float.h>

#ifdef ALG_THREADS_USED
#include <pthread.h>
#define ALG_THREADS_MAX	(64)
#endif /* ALG_THREADS_USED */

typedef enum			     /* Forward or inverse Fourier transform */
{
  ALG_FOUR_DIR_FWD = 0,
  ALG_FOUR_DIR_INV = 1
} AlgFourDirection;

#ifdef ALG_THREADS_USED
#define ALG_FOUR_THR_NUM1D 512  /* Min size of 1D complex FT, to use threads */

typedef struct    /* Used for args by AlgFourThrHart1D(), AlgFourThrReal1D() */
{						/* and AlgFourThrRealInv1D() */
  double	*data;
  int		num;
  int		step;
  int		cThr;
} AlgFourArgs1;

typedef struct        /* Used for args by AlgFourThr1D() and AlgFourThrInv1D */
{
  double	*real;
  double	*imag;
  int		num;
  int		step;
  int		cThr;
} AlgFourArgs2;

typedef struct			 /* Used for args by AlgFourThrRepXYReal1D() */
{
  double	**data;
  double	*reBuf;
  double	*imBuf;
  int		numData;
  int		stepData;
  int		repX;
  int		repY;
  AlgFourDirection dir;
  int		cThr;
} AlgFourArgs3;

typedef struct        		     /* Used for args by AlgFourThrRepXY1D() */
{
  double	**real;
  double	**imag;
  double	*reBuf;
  double	*imBuf;
  int		numData;
  int		stepData;
  int		repX;
  int		repY;
  AlgFourDirection dir;
  int		cThr;
} AlgFourArgs4;
#endif /* ALG_THREADS_USED */

static void	AlgFourRepXY1D(double **real, double **imag,
			       double *reBuf, double *imBuf,
			       int numData, int stepData,
			       int repX, int repY,
			       AlgFourDirection dir, int cThr),
		AlgFourRepXYReal1D(double **data,
				   double *reBuf, double *imBuf,
				   int numData, int stepData,
				   int repX, int repY,
				   AlgFourDirection dir, int cThr);
#ifdef ALG_THREADS_USED
static void	*AlgFourThrHart1D(AlgFourArgs1 *args),
		*AlgFourThrReal1D(AlgFourArgs1 *args),
		*AlgFourThrRealInv1D(AlgFourArgs1 *args),
		*AlgFourThr1D(AlgFourArgs2 *args),
		*AlgFourThrInv1D(AlgFourArgs2 *args),
		*AlgFourThrRepXYReal1D(AlgFourArgs3 *args),
		*AlgFourThrRepXY1D(AlgFourArgs4 *args);
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFourHart1D						*
* Returns:	void							*
* Purpose:	Computes the Hartley transform of the given one		*
*		dimensional data, and does it in place.			*
* Global refs:	-							*
* Parameters:	double *data:		Given data.			*
*		int num:		Number of data.			*
*		int step:		Offset in data elements between	*
*					the data to be transformed.	*
*		int cThr:		Concurrent threads available,	*
*					if cThr <= 1 then no threads	*
*					will be created.		*
************************************************************************/
void		AlgFourHart1D(double *data, int num, int step, int cThr)
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
	  ("AlgFourHart1D FE 0x%lx %d %d %d\n",
	   (unsigned long )data, num, step, cThr));
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
    ++pTwo0;
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

#ifdef ALG_THREADS_USED
/************************************************************************
* Function:	AlgFourThrHart1D					*
* Returns:	void *:			Always NULL.			*
* Purpose:	Simple wrapper for AlgFourHart1D(), used for thread 	*
*		creation.						*
* Global refs:	-							*
* Parameters:	AlgFourArgs1 *args:	Parameter list.			*
************************************************************************/
static void	*AlgFourThrHart1D(AlgFourArgs1 *args)
{
  AlgFourHart1D(args->data, args->num, args->step, args->cThr);
  return(NULL);
}
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFour1D						*
* Returns:	void							*
* Purpose:	Computes the Fourier transform of the given one		*
*		dimensional complex data, and does it in place.		*
* Global refs:	-							*
* Parameters:	double *real:		Given real data.		*
*		double *imag:		Given imaginary data.		*
*		int num:		Number of data.			*
*		int step:		Offset in data elements between	*
*					the data to be transformed.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFour1D(double *real, double *imag, int num, int step,
			  int cThr)
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
  int		count,
  		threadCreated = 0;
#ifdef ALG_THREADS_USED
  AlgFourArgs1	thrArgs;
  pthread_t	thrId;
#endif /* ALG_THREADS_USED */

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour1D FE 0x%lx 0x%lx %d %d %d\n",
	   (unsigned long )real, (unsigned long )imag, num, step, cThr));
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
#ifdef ALG_THREADS_USED
  if((cThr > 1) && (num >= ALG_FOUR_THR_NUM1D))
  {
    --cThr;
    thrArgs.data = real;
    thrArgs.num = num;
    thrArgs.step = step;
    thrArgs.cThr = cThr;
    if(pthread_create(&thrId, NULL, (void *(*)(void *) )AlgFourThrHart1D,
    		      (void *)&thrArgs) == 0)
    {
      threadCreated = 1;
      AlgFourHart1D(imag, num, step, cThr);
      (void )pthread_join(thrId, NULL);
      ++cThr;
    }
  }
#endif /* ALG_THREADS_USED */
  if(threadCreated == 0)
  {
    AlgFourHart1D(real, num, step, cThr);
    AlgFourHart1D(imag, num, step, cThr);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour1D FX\n"));
}

#ifdef ALG_THREADS_USED
/************************************************************************
* Function:	AlgFourThr1D						*
* Returns:	void *:			Always NULL.			*
* Purpose:	Simple wrapper for AlgFour1D(), used for thread 	*
*		creation.						*
* Global refs:	-							*
* Parameters:	AlgFourArgs2 *args:	Parameter list.			*
************************************************************************/
static void	*AlgFourThr1D(AlgFourArgs2 *args)
{
  AlgFour1D(args->real, args->imag, args->num, args->step, args->cThr);
  return(NULL);
}
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFourInv1D						*
* Returns:	void							*
* Purpose:	Computes the inverse Fourier transform of the given	*
*		complex one dimensional data, and does it in place.	*
* Global refs:	-							*
* Parameters:	double *real:		Given real data.		*
*		double *imag:		Given imaginary data.		*
*		int num:		Number of data.			*
*		int step:		Offset in data elements between	*
*					the data to be transformed.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFourInv1D(double *real, double *imag, int num, int step,
			     int cThr)
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
  int		count,
  		threadCreated = 0;
#ifdef ALG_THREADS_USED
  AlgFourArgs1  thrArgs;
  pthread_t      thrId;
#endif /* ALG_THREADS_USED */

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv1D FE 0x%lx 0x%lx %d %d %d\n",
	   (unsigned long )real, (unsigned long )imag, num, step, cThr));
#ifdef ALG_THREADS_USED
  if((cThr > 1) && (num >= ALG_FOUR_THR_NUM1D))
  {
    --cThr;
    thrArgs.data = real;
    thrArgs.num = num;
    thrArgs.step = step;
    thrArgs.cThr = cThr;
    if(pthread_create(&thrId, NULL, (void *(*)(void *) )AlgFourThrHart1D,
    		      (void *)&thrArgs) == 0)
    {
      threadCreated = 1;
      AlgFourHart1D(imag, num, step, cThr);
      (void )pthread_join(thrId, NULL);
      ++cThr;
    }
  }
#endif /* ALG_THREADS_USED */
  if(threadCreated == 0)
  {
    AlgFourHart1D(real, num, step, cThr);
    AlgFourHart1D(imag, num, step, cThr);
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

#ifdef ALG_THREADS_USED
/************************************************************************
* Function:	AlgFourThrInv1D						*
* Returns:	void *:			Always NULL.			*
* Purpose:	Simple wrapper for AlgFourInv1D(), used for thread 	*
*		creation.						*
* Global refs:	-							*
* Parameters:	AlgFourArgs2 *args:	Parameter list.			*
************************************************************************/
void		*AlgFourThrInv1D(AlgFourArgs2 *args)
{
  AlgFourInv1D(args->real, args->imag, args->num, args->step, args->cThr);
  return(NULL);
}
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFourRepXY1D						*
* Returns:	void							*
* Purpose:	Computes repeated the Fourier transforms of the given	*
*	 	one dimensional complex data sets.			*
*		These may either be done wrt the rows or columns of	*
*		the data. Buffers may be provided for the columns to	*
*		improve efficiency, if provided (non NULL) then these	*
*		should be large enough to hold a single column of data.	*
*		Multiple threads may be used (if cThr > 1) for repeated	*
*		rows, but (because of the column buffers) not for the	*
*		columns. Although AlgFour1D()/AlgFourInv1D() can make	*
*		use of two threads in themselves.			*
* Global refs:	-							*
* Parameters:	double **real:		Given real data sets.		*
*		double **imag:		Given imaginary data sets.	*
*		double *reBuf:		Given buffer(s) for real data.	*
*		double *imBuf:		Given buffer(s) for imaginary	*
*					data.				*
*		int numData:		Number of data in a row/column.	*
*		int stepData:		Offset in data elements between	*
*					the data to be transformed.	*
*		int repX:		Number of data columns to be	*
*					transformed.			*
*		int repY:		Number of data rows to be	*
*					transformed.			*
*		AlgFourDirection dir:	Forward or inverse transform.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
static void	AlgFourRepXY1D(double **real, double **imag,
			       double *reBuf, double *imBuf,
			       int numData, int stepData,
			       int repX, int repY,
			       AlgFourDirection dir, int cThr)
{
  int		idX,
		idY,
		threadCreated = 0;
  double	*tDp0,
		*tDp1,
		*tDp2,
		*tDp3;
#ifdef ALG_THREADS_USED
  int           offY,
                numThr,
                repThr;
  AlgFourArgs4  *thrArgP;
  AlgFourArgs4  thrArgs[ALG_THREADS_MAX];
  pthread_t      thrIds[ALG_THREADS_MAX];
#endif /* ALG_THREADS_USED */

  if(repY)			                           /* Transform rows */
  {
#ifdef ALG_THREADS_USED
    if(cThr > 1)  /* Calc num threads and num of repeated FT for each thread */
    {
      if(ALG_THREADS_MAX > cThr)
      {
        if(cThr < repY)
	{
	  numThr = cThr;
	}
	else
	{
	  numThr = repY;
	}
      }
      else
      {
        if(ALG_THREADS_MAX > repY)
	{
	  numThr = ALG_THREADS_MAX;
	}
	else
	{
	  numThr = repY;
	}
      }
      repThr = repY / numThr;
    }
    else
    {
      numThr = 1;
      repThr = repY;
    }
    offY = 0;                      /* Build parameter structures for threads */
    thrArgP = thrArgs;
    for(idY = 0; idY < numThr; ++idY)
    {
      thrArgP->real = real + offY;
      thrArgP->imag = imag + offY;
      thrArgP->reBuf = NULL;
      thrArgP->imBuf = NULL;
      thrArgP->numData = numData;
      thrArgP->stepData = stepData;
      thrArgP->repX = 0;
      thrArgP->repY = repThr;
      thrArgP->dir = dir;
      thrArgP->cThr = 1;
      offY += repThr;
      ++thrArgP;
    }
    if(numThr > 1)
    {
      thrArgP = thrArgs + numThr - 1;
      if((cThr - numThr) > 0)           /* Don't waste any left over threads */
      {
        thrArgP->cThr += cThr - numThr;
      }
      thrArgP->repY += repY % numThr;     /* Add on any remaining FT repeats */
      for(idY = 1; idY < numThr; ++idY)            /* Create the new threads */
      {
	if(pthread_create(thrIds + idY, NULL,
			  (void *(*)(void *))AlgFourThrRepXY1D,
			  (void *)(thrArgs + idY)) == 0)
	{
	  threadCreated = 1;
	}
      }
    }
    if(threadCreated)
    {
      thrArgP = thrArgs;      /* Use this thread and terminate any recursion */
      if(dir == ALG_FOUR_DIR_FWD)
      {
	for(idY = 0; idY < thrArgP->repY; ++idY)
	{
	  AlgFour1D(*(thrArgP->real + idY), *(thrArgP->imag + idY),
		    thrArgP->numData, thrArgP->stepData, thrArgP->cThr);
	}
      }
      else
      {
	for(idY = 0; idY < thrArgP->repY; ++idY)
	{
	  AlgFourInv1D(*(thrArgP->real + idY), *(thrArgP->imag + idY),
		       thrArgP->numData, thrArgP->stepData, thrArgP->cThr);
	}
      }
      if(numThr > 1) /* Wait for all of the threads created for repeated FTs */
      {
	for(idY = 1; idY < numThr; ++idY)
	{
	  (void )pthread_join(thrIds[idY], NULL);
        }
      }
    }
#endif /* ALG_THREADS_USED */
    if(threadCreated == 0)
    {
      if(dir == ALG_FOUR_DIR_FWD)
      {
	for(idY = 0; idY < repY; ++idY)
	{
	  AlgFour1D(*(real + idY), *(imag + idY), numData, stepData, cThr);
	}
      }
      else
      {
	for(idY = 0; idY < repY; ++idY)
	{
	  AlgFourInv1D(*(real + idY), *(imag + idY), numData, stepData, cThr);
	}
      }
    }
  }
  else if(repX)						/* Transform columns */
  {
    if(reBuf && imBuf)                     /* Use column buffers if provided */
    {
      tDp0 = reBuf;				          /* Copy 1st column */
      tDp1 = imBuf;
      for(idY = 0; idY < numData; ++idY)
      {
	*tDp0++ = **(real + idY);
	*tDp1++ = **(imag + idY);
      }
      if(dir == ALG_FOUR_DIR_FWD)		     /* Transform 1st column */
      {
	AlgFour1D(reBuf, imBuf, numData, 1, cThr);
      }
      else
      {
	AlgFourInv1D(reBuf, imBuf, numData, 1, cThr);
      }
      for(idX = 1; idX < repX; ++idX)
      {
	tDp0 = reBuf;	         /* Copy back previous, and copy next column */
	tDp1 = imBuf;
	for(idY = 0; idY < numData; ++idY)
	{
	  tDp2 = *(real + idY) + idX - 1;
	  *tDp2++ = *tDp0;
	  *tDp0++ = *tDp2;
	  tDp3 = *(imag + idY) + idX - 1;
	  *tDp3++ = *tDp1;
	  *tDp1++ = *tDp3;
	}
	if(dir == ALG_FOUR_DIR_FWD)		 /* Transform current column */
	{
	  AlgFour1D(reBuf, imBuf, numData, 1, cThr);
	}
	else
	{
	  AlgFourInv1D(reBuf, imBuf, numData, 1, cThr);
        }
      }
      tDp0 = reBuf;                                 /* Copy back last column */
      tDp1 = imBuf;
      for(idY = 0; idY < numData; ++idY)
      {
	*(*(real + idY) + repX - 1) = *tDp0++;
	*(*(imag + idY) + repX - 1) = *tDp1++;
      }
    }
    else              /* If no column buffers transform the columns in place */
    {
      if(dir == ALG_FOUR_DIR_FWD)
      {
	for(idX = 0; idX < repX; ++idX)
	{
	  AlgFour1D(*real + idX, *imag + idX, numData, stepData, cThr);
        }
      }
      else
      {
	for(idX = 0; idX < repX; ++idX)
	{
	  AlgFourInv1D(*real + idX, *imag + idX, numData, stepData, cThr);
        }
      }
    }
  }
}

#ifdef ALG_THREADS_USED
/************************************************************************
* Function:	AlgFourThrRepXY1D					*
* Returns:	void *:			Always NULL.			*
* Purpose:	Simple wrapper for AlgFourRepXY1D(), used for thread 	*
*		creation.						*
* Global refs:	-							*
* Parameters:	AlgFourArgs4 *args:	Parameter list.			*
************************************************************************/
static void	*AlgFourThrRepXY1D(AlgFourArgs4 *args)
{
  AlgFourRepXY1D(args->real, args->imag, args->reBuf, args->imBuf,
  		 args->numData, args->stepData,
		 args->repX, args->repY, args->dir, args->cThr);
  return(NULL);
}
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFourReal1D						*
* Returns:	void							*
* Purpose:	Computes the Fourier transform of the given one 	*
*		dimensional real data, and does it in place.		*
* Global refs:	-							*
* Parameters:	double *real:		Given real data.		*
*		int num:		Number of data.			*
*		int step:		Offset in data elements between	*
*					the data to be transformed.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFourReal1D(double *real, int num, int step, int cThr)
{
  double	tD0,
		tD1;
  double	*tRp0,
		*tRp1;
  int		count;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal1D FE 0x%lx %d %d %d\n",
	   (unsigned long )real, num, step, cThr));
  tRp0 = real + step;
  tRp1 = real + ((num - 1) * step);
  count = num / 2;
  AlgFourHart1D(real, num, step, cThr);
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

#ifdef ALG_THREADS_USED
/************************************************************************
* Function:	AlgFourThrReal1D					*
* Returns:	void *:			Always NULL.			*
* Purpose:	Simple wrapper for AlgFourReal1D(), used for thread 	*
*		creation.						*
* Global refs:	-							*
* Parameters:	AlgFourArgs1 *args:	Parameter list.			*
************************************************************************/
void		*AlgFourThrReal1D(AlgFourArgs1 *args)
{
  AlgFourReal1D(args->data, args->num, args->step, args->cThr);
  return(NULL);
}
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFourRealInv1D					*
* Returns:	void							*
* Purpose:	Computes the inverse Fourier transform of the given one	*
*		one dimensional real data, and does it in place.	*
* Global refs:	-							*
* Parameters:	double *real:		Given real/complex data.	*
*		int num:		Number of data.			*
*		int step:		Offset in data elements between	*
*					the data to be transformed.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFourRealInv1D(double *real, int num, int step, int cThr)
{
  double	tD0,
		tD1;
  double	*tRp0,
		*tRp1;
  int		count;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv1D FE 0x%lx %d %d %d\n",
	   (unsigned long )real, num, step, cThr));
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
  AlgFourHart1D(real, num, step, cThr);
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv1D FX\n"));
}

#ifdef ALG_THREADS_USED
/************************************************************************
* Function:	AlgFourThrRealInv1D					*
* Returns:	void *:			Always NULL.			*
* Purpose:	Simple wrapper for AlgFourRealInv1D(), used for thread 	*
*		creation.						*
* Global refs:	-							*
* Parameters:	AlgFourArgs1 *args:	Parameter list.			*
************************************************************************/
void		*AlgFourThrRealInv1D(AlgFourArgs1 *args)
{
  AlgFourRealInv1D(args->data, args->num, args->step, args->cThr);
  return(NULL);
}
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFourRepXYReal1D					*
* Returns:	void							*
* Purpose:	Computes repeated the Fourier transforms of the given	*
*	 	one dimensional real data sets.				*
*		These may either be done wrt the rows or columns of	*
*		the data. Buffers may be provided for the columns to	*
*		improve efficiency, if provided (non NULL) then these	*
*		should be large enough to hold a single column of data.	*
*		Multiple threads may be used (if cThr > 1) for repeated	*
*		rows, but (because of the column buffers) not for the	*
*		columns. Although AlgFour1D()/AlgFourInv1D() can make	*
*		use of two threads in themselves.			*
* Global refs:	-							*
* Parameters:	double **data:		Given data sets.		*
*		double *reBuf:		Given buffer(s) for real data.	*
*		double *imBuf:		Given buffer(s) for imaginary	*
*					data.				*
*		int numData:		Number of data in row/column.	*
*		int stepData:		Offset in data elements between	*
*					the data to be transformed.	*
*		int repX:		Number of data columns to be	*
*					transformed.			*
*		int repY:		Number of data rows to be	*
*					transformed.			*
*		AlgFourDirection dir:	Forward or inverse transform.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
static void	AlgFourRepXYReal1D(double **data, double *reBuf, double *imBuf,
				   int numData, int stepData,
				   int repX, int repY,
				   AlgFourDirection dir, int cThr)
{
  int		idX,
		idY,
		halfData,
		threadCreated = 0;
  double	*tDp0,
		*tDp1,
		*tDp2,
		*tDp3;
#ifdef ALG_THREADS_USED
  int		offY,
  		numThr,
  		repThr;
  AlgFourArgs3	*thrArgP;
  AlgFourArgs3	thrArgs[ALG_THREADS_MAX];
  pthread_t	thrIds[ALG_THREADS_MAX];
#endif /* ALG_THREADS_USED */

  if(repY)						   /* Transform rows */
  {
#ifdef ALG_THREADS_USED
    if(cThr > 1)  /* Calc num threads and num of repeated FT for each thread */
    {
      if(ALG_THREADS_MAX > cThr)
      {
        if(cThr < repY)
	{
	  numThr = cThr;
	}
	else
	{
	  numThr = repY;
	}
      }
      else
      {
        if(ALG_THREADS_MAX > repY)
	{
	  numThr = ALG_THREADS_MAX;
	}
	else
	{
	  numThr = repY;
	}
      }
      repThr = repY / numThr;
    }
    else
    {
      numThr = 1;
      repThr = repY;
    }
    offY = 0;			   /* Build parameter structures for threads */
    thrArgP = thrArgs;
    for(idY = 0; idY < numThr; ++idY)
    {
      thrArgP->data = data + offY;
      thrArgP->reBuf = NULL;
      thrArgP->imBuf = NULL;
      thrArgP->numData = numData;
      thrArgP->stepData = stepData;
      thrArgP->repX = 0;
      thrArgP->repY = repThr;
      thrArgP->dir = dir; 
      thrArgP->cThr = 1;
      offY += repThr;
      ++thrArgP;
    }
    if(numThr > 1)
    {
      thrArgP = thrArgs + numThr - 1;
      if((cThr - numThr) > 0)		/* Don't waste any left over threads */
      {
	thrArgP->cThr += cThr - numThr;
      }
      thrArgP->repY += repY % numThr;	  /* Add on any remaining FT repeats */
      for(idY = 1; idY < numThr; ++idY) 	   /* Create the new threads */
      {
	if(pthread_create(thrIds + idY, NULL,
			  (void *(*)(void *) )AlgFourThrRepXYReal1D,
			  (void *)(thrArgs + idY)) == 0)
	{
	  threadCreated = 1;
	}
      }
    }
    if(threadCreated)
    {
      thrArgP = thrArgs;      /* Use this thread and terminate any recursion */
      if(dir == ALG_FOUR_DIR_FWD)
      {
	for(idY = 0; idY < thrArgP->repY; ++idY)
	{
	  AlgFourReal1D(*(thrArgP->data + idY), thrArgP->numData,
			thrArgP->stepData, thrArgP->cThr);
	}
      }
      else
      {
	for(idY = 0; idY < thrArgP->repY; ++idY)
	{
	  AlgFourRealInv1D(*(thrArgP->data + idY), thrArgP->numData,
			   thrArgP->stepData, thrArgP->cThr);
	}
      }
      if(numThr > 1) /* Wait for all of the threads created for repeated FTs */
      {
	for(idY = 1; idY < numThr; ++idY) 
	{
	  (void )pthread_join(thrIds[idY], NULL);
	}
      }
    }
#endif /* ALG_THREADS_USED */
    if(threadCreated == 0)
    {
      if(dir == ALG_FOUR_DIR_FWD)
      {
	for(idY = 0; idY < repY; ++idY)
	{
	  AlgFourReal1D(*(data + idY), numData, 1, cThr);
	}
      }
      else
      {
	for(idY = 0; idY < repY; ++idY)
	{
	  AlgFourRealInv1D(*(data + idY), numData, 1, cThr);
	}
      }
    }
  }
  else if(repX)						/* Transform columns */
  {
    halfData = repX / 2;
    if(reBuf && imBuf)			   /* Use column buffers if provided */
    {
      tDp0 = reBuf;
      for(idY = 0; idY < numData; ++idY)		    /* Copy column 0 */
      {
	*tDp0++ = **(data + idY);
      }
      if(dir == ALG_FOUR_DIR_FWD)		       /* Transform column 0 */
      {
	AlgFourReal1D(reBuf, numData, 1, cThr);
      }
      else
      {
	AlgFourRealInv1D(reBuf, numData, 1, cThr);
      }
      tDp0 = reBuf;
      for(idY = 0; idY < numData; ++idY) /* Copy back col 0, copy col numX/2 */
      {
	**(data + idY) = *tDp0;
	*tDp0++ = *(*(data + idY)+ halfData);
      }
      if(dir == ALG_FOUR_DIR_FWD)		  /* Transform column numX/2 */
      {
	AlgFourReal1D(reBuf, numData, 1, cThr);
      }
      else
      {
	AlgFourRealInv1D(reBuf, numData, 1, cThr);
      }
      tDp0 = reBuf;
      for(idY = 0; idY < numData; ++idY)	  /* Copy back column numX/2 */
      {
	*(*(data + idY)+ halfData) = *tDp0++;
      }
      tDp0 = reBuf;
      tDp1 = imBuf;
      for(idY = 0; idY < numData; ++idY)    /* Copy columns 1 and numX/2 + 1 */
      {
	tDp2 = *(data + idY) + 1;
	*tDp0++ = *tDp2;
	*tDp1++ = *(tDp2 + halfData);
      }
      if(dir == ALG_FOUR_DIR_FWD)      /* Transform columns 1 and numX/2 + 1 */
      {
        AlgFour1D(reBuf, imBuf, numData, 1, cThr); 
      }
      else
      {
        AlgFourInv1D(reBuf, imBuf, numData, 1, cThr);
      }
      for(idX = 2; idX < halfData;
          ++idX)    /* Copy back previous, copy and transform current column */
      {
	tDp0 = reBuf;
	tDp1 = imBuf;
	for(idY = 0; idY < numData; ++idY)
	{
	  tDp2 = *(data + idY) + idX - 1;
	  tDp3 = tDp2 + halfData;
	  *tDp2++ = *tDp0;
	  *tDp0++ = *tDp2;
	  *tDp3++ = *tDp1;
	  *tDp1++ = *tDp3;
	}
	if(dir == ALG_FOUR_DIR_FWD)
	{
	  AlgFour1D(reBuf, imBuf, numData, 1, cThr);
	}
        else
	{
	  AlgFourInv1D(reBuf, imBuf, numData, 1, cThr);
        }
      }
      tDp0 = reBuf;
      tDp1 = imBuf;
      for(idY = 0; idY < numData; ++idY)      /* Copy back  last two columns */
      {
	*(*(data + idY) + halfData - 1) = *tDp0++;
	*(*(data + idY) + repX - 1) = *tDp1++;
      }
    }
    else    /* If no column buffers then just transform the columns in place */
    {
      if(dir == ALG_FOUR_DIR_FWD)
      {
	AlgFourReal1D(*data, numData, stepData, cThr); /* Transform column 0 */
	AlgFourReal1D(*data + halfData, numData, stepData,
		      cThr);			  /* Transform column numX/2 */
	for(idX = 1; idX < halfData; ++idX)	  /* Transform other columns */
	{
	  AlgFour1D(*data + idX, *data + halfData + idX, numData,
		    stepData, cThr);
	}
      }
      else
      {
	AlgFourRealInv1D(*data, numData, stepData, cThr); /* Transform col 0 */
	AlgFourRealInv1D(*data + halfData, numData, stepData,
		      cThr);			  /* Transform column numX/2 */
	for(idX = 1; idX < halfData; ++idX)	  /* Transform other columns */
	{
	  AlgFourInv1D(*data + idX, *data + halfData + idX, numData,
		    stepData, cThr);
	}
      }
    }
  }
}

#ifdef ALG_THREADS_USED
/************************************************************************
* Function:	AlgFourThrRepXYReal1D					*
* Returns:	void *:			Always NULL.			*
* Purpose:	Simple wrapper for AlgFourRepXYReal1D(), used for 	*
*		thread creation.					*
* Global refs:	-							*
* Parameters:	AlgFourArgs3 *args:	Parameter list.			*
************************************************************************/
static void	*AlgFourThrRepXYReal1D(AlgFourArgs3 *args)
{
  AlgFourRepXYReal1D(args->data, args->reBuf, args->imBuf,
  		     args->numData, args->stepData,
  	             args->repX, args->repY, args->dir, args->cThr);
  return(NULL);
}
#endif /* ALG_THREADS_USED */

/************************************************************************
* Function:	AlgFour2D						*
* Returns:	void							*
* Purpose:	Computes the Fourier transform of the given two		*
*		dimensional complex data, and does it in place.		*
*		If the two buffer pointers are NULL then the transform	*
*		will be done without copying the data between temporary	*
*		buffers.						*
* Global refs:	-							*
* Parameters:	double **real:		Given real data.		*
*		double **imag:		Given imaginary data.		*
*		double *reBuf:		Given buffer for real data,	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		double *imBuf:		Given buffer for imaginary data	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		int numX:		Number of data in each row.	*
*		int numX:		Number of data in each column.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFour2D(double **real, double **imag,
			  double *reBuf, double *imBuf, int numX, int numY,
			  int cThr)
{
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour2D FE 0x%lx 0x%lx 0x%lx 0x%lx %d %d %d\n",
	   (unsigned long )real, (unsigned long )imag,
	   (unsigned long )reBuf, (unsigned long)imBuf, numX, numY, cThr));
  AlgFourRepXY1D(real, imag, NULL,  NULL, numX, 1,
  		 0, numY, ALG_FOUR_DIR_FWD, cThr);	   /* Transform rows */
  AlgFourRepXY1D(real, imag, reBuf, imBuf, numY, numX,
  		 numX, 0, ALG_FOUR_DIR_FWD, cThr);	/* Transform columns */
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFour2D FX\n"));
}

/************************************************************************
* Function:	AlgFourInv2D						*
* Returns:	void							*
* Purpose:	Computes the inverse Fourier transform of the given two	*
*		dimensional complex data, and does it in place.		*
*		If the two buffer pointers are NULL then the transform	*
*		will be done without copying the data between temporary	*
*		buffers.						*
* Global refs:	-							*
* Parameters:	double **real:		Given real data.		*
*		double **imag:		Given imaginary data.		*
*		double *reBuf:		Given buffer for real data,	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		double *imBuf:		Given buffer for imaginary data	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		int numX:		Number of data in each row.	*
*		int numX:		Number of data in each column.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFourInv2D(double **real, double **imag,
			     double *reBuf, double *imBuf,
			     int numX, int numY, int cThr)
{
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv2D FE 0x%lx 0x%lx 0x%lx 0x%lx %d %d %d\n",
	   (unsigned long )real, (unsigned long )imag,
	   (unsigned long )reBuf, (unsigned long)imBuf, numX, numY, cThr));
  AlgFourRepXY1D(real, imag, NULL,  NULL, numX, 1,
  		 0, numY, ALG_FOUR_DIR_INV, cThr);	   /* Transform rows */
  AlgFourRepXY1D(real, imag, reBuf, imBuf, numY, numX,
  		 numX, 0, ALG_FOUR_DIR_INV, cThr);	/* Transform columns */
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourInv2D FX\n"));
}

/************************************************************************
* Function:	AlgFourReal2D						*
* Returns:	void							*
* Purpose:	Computes the Fourier transform of the given two		*
*		dimensional real data, and does it in place.		*
*		If the two buffer pointers are NULL then the transform	*
*		will be done without copying the data between temporary	*
*		buffers.						*
* Global refs:	-							*
* Parameters:	double **real:		Given real data.		*
*		double *reBuf:		Given buffer for real data,	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		double *imBuf:		Given buffer for imaginary data	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		int numX:		Number of data in each row.	*
*		int numX:		Number of data in each column.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFourReal2D(double **real, double *reBuf, double *imBuf,
			      int numX, int numY, int cThr)
{
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal2D FE 0x%lx 0x%lx 0x%lx %d %d %d\n",
	   (unsigned long )real, (unsigned long )reBuf, (unsigned long)imBuf,
	   numX, numY, cThr));
  AlgFourRepXYReal1D(real, NULL,  NULL, numX, 1,
  		     0, numY, ALG_FOUR_DIR_FWD, cThr);	   /* Transform rows */
  AlgFourRepXYReal1D(real, reBuf, imBuf, numY, numX,
  		     numX, 0, ALG_FOUR_DIR_FWD, cThr);	   /* Transform cols */
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourReal2D FX\n"));
}

/************************************************************************
* Function:	AlgFourRealInv2D					*
* Returns:	void							*
* Purpose:	Computes the Fourier transform of the given two		*
*		dimensional data which resulted from a transform using	*
*		AlgFourReal2D(), and does it in place.			*
*		If the two buffer pointers are NULL then the transform	*
*		will be done without copying the data between temporary	*
*		buffers.						*
* Global refs:	-							*
* Parameters:	double **real:		Given real/complex data.	*
*		double *reBuf:		Given buffer for real data,	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		double *imBuf:		Given buffer for imaginary data	*
*					may be NULL or a buffer region	*
*					suitable for a single column.	*
*		int numX:		Number of data in each row.	*
*		int numX:		Number of data in each column.	*
*		int cThr:		Concurrent threads available,   *
*					if cThr <= 1 then no threads    *
*					will be created.		*
************************************************************************/
void		AlgFourRealInv2D(double **real,double *reBuf, double *imBuf,
				 int numX, int numY, int cThr)
{
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv2D FE 0x%lx 0x%lx 0x%lx %d %d %d\n",
	   (unsigned long )real, (unsigned long )reBuf, (unsigned long)imBuf,
	   numX, numY, cThr));
  AlgFourRepXYReal1D(real, reBuf, imBuf, numY, numX,
  		     numX, 0, ALG_FOUR_DIR_INV, cThr);	   /* Transform cols */
  AlgFourRepXYReal1D(real, NULL,  NULL, numX, 1,
  		     0, numY, ALG_FOUR_DIR_INV, cThr);	   /* Transform rows */
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgFourRealInv2D FX\n"));
}
