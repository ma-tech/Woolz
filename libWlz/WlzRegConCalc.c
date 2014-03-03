#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRegConCalc_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzRegConCalc.c
* \author       Bill Hill
* \date         July 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Computes region connected calculus spatial classifications.
* \ingroup	WlzBinaryOps
*/
#include <Wlz.h>

/*!
* \return	RCC classification of the given objects, ie object 0 is a
*               returned classification of object 1.
* \ingroup	WlzBinaryOps
* \brief	The given pair of spatial domain objects are classified
*		using a RCC.
*		For an explanation of RCC8 classifications
*		see the type definition ::WlzRegConRCC and the paper:
*		D.A. Randell, etal,
*		"Discrete Mereotopology for Spatial Reasoning in
*		Automated Histological Image Analysis", PAMI 35(3) 2013.
*		The RCC8 has been extended to include encloses and
*		surrounds.
* 		The classification is performed using simple combinations
* 		of the Woolz union, intersection, difference, dilation
* 		and convex hull operators on an ordered pair of
* 		spatial domains(\f$\Omega_0\f$ and \f$\Omega_1\f$):
*               \f{eqnarray*}{
                  C_0 &\leftarrow&
		    \Omega_0   \cap \Omega_1
		    \neq \emptyset \\
                  C_1 &\leftarrow&
		    \Omega_0^+ \cap \Omega_1
		    \neq \emptyset \\
                  C_2 &\leftarrow&
		    (\Omega_0 \oplus \Omega_1)
		    \neq \emptyset \\
                  C_3 &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus \Omega_0 
		    \neq \emptyset \\
                  C_4 &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus \Omega_1 
		    \neq \emptyset \\
                  C_5 &\leftarrow&
		    (\Omega_0^+ \cup \Omega_1)   \oplus \Omega_1 
		    \neq \emptyset \\
                  C_6 &\leftarrow&
		    (\Omega_0   \cup \Omega_1^+) \oplus \Omega_0 
		    \neq \emptyset \\
		  C_7 &\leftarrow&
		    2|\Omega_0^v \cap \Omega_1|   \ge |\Omega_1| \\
		  C_8 &\leftarrow&
		    2|\Omega_0   \cap \Omega_1^v| \ge |\Omega_0|
  		\f}
*		where
*		  are the \f$\cup\f$, \f$\cap\f$ and \f$\oplus\f$
*		  are the set union (logical or), intersection (logical and)
*		  and difference (logical exclusive or) operators;
*		  \f$\Omega^+\f$ indicates the dilation of \f$\Omega\f$,
*		  \f$\Omega^v\f$ the convex hull of \f$\Omega\f$ and
*		  \f$|\Omega|\f$ the cardinality (area or volume) of
*		  \f$\Omega\f$.
* 		The decision tree for the classification excluding
* 		enclosure is:
* 		\f[
C_0?
\left\{
\begin{array}{ll}
0 & C_1?
    \left\{
    \begin{array}{ll}
    0 & DC \\
      & \\
      & \\
    1 & C_4?
	\left\{
	\begin{array}{ll}
	0 & SUR \\
	  & \\
	  & \\
	1 & C_3?
	    \left\{
	    \begin{array}{ll}
	    0 & SURI \\
	      & \\
	    1 & EC
	    \end{array}
	    \right. \\
	\end{array}
	\right. \\
    \end{array}
    \right. \\
  & \\
1 & C_2?
    \left\{
    \begin{array}{ll}
    0 & EQ \\
      & \\
      & \\
    1 & C_4?
	\left\{
	\begin{array}{ll}
	0 & C_5?
	    \left\{
	    \begin{array}{ll}
	    0 & NTPP \\
	      & \\
	    1 & TPP
	    \end{array}
	    \right. \\
	  & \\
	  & \\
	1 & C_3?
	    \left\{
	    \begin{array}{ll}
	    0 & C_6?
		\left\{
		\begin{array}{ll}
		0 & NTPPI \\
		 & \\
		1 & TPPI
		\end{array}
		\right. \\
	      & \\
	    1 & PO
	    \end{array}
	    \right. \\
	\end{array}
	\right. \\
    \end{array}
    \right. \\
\end{array}
\right.
		\f]
*		Following this classification, enclosure may be added
*		using \f$C_7\f$ and \f$C_8\f$ to test for enclosure
*		and it's inverse respectively.
*		The normalised volumes are computed for each classification
*		as below:
* 		<table width="500" border="0">
		<caption>
		  Basic Morphological Operations for RCC8 Spatial
		  Relationships
		</caption>
                <tr>
                  <td>RCC</td>
		  <td>Normalised Volume</td>
                </tr>
		<tr>
		  <td>\f$EMPTY(\Omega_0,\Omega_1)\f$</td>
		  <td>0.0</td>
		</tr>
		<tr>
		  <td>\f$DC(\Omega_0,\Omega_1)\f$</td>
		  <td>0.0</td>
		</tr>
		<tr>
		  <td>\f$EC(\Omega_0,\Omega_1)\f$</td>
		  <td>0.0</td>
		</tr>
		<tr>
		  <td>\f$EQ(\Omega_0,\Omega_1)\f$</td>
		  <td>1.0</td>
		</tr>
		<tr>
		  <td>\f$PO(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0 \cap \Omega_1|/
		         |\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$TPP(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$NTPP(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$TPPI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$NTPPI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$ENC(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0 \cup \Omega_1^v|\f$</td>
		</tr>
		<tr>
		  <td>\f$ENCI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0^v \cup \Omega_1|\f$</td>
		</tr>
		</table>
		When computing the normalised volumes enclosure has the
		higher priority.
* \param	obj0			First given spatial domain object.
* \param	obj1			Second given spatial domain object.
* \param	dstNrmVol		Destination pointer for the normalized
* 					volume (see above), may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzRegConRCC 	WlzRegConCalcRCC(WlzObject *obj0, WlzObject *obj1,
				 double *dstNrmVol, WlzErrorNum *dstErr)
{
  int		n = 1;
  double	vol = 0.0;
  WlzObject	*objI = NULL,
  		*objD = NULL,
  		*objU = NULL,
		*objU0 = NULL,
		*objU1 = NULL;
  WlzRegConRCC con = WLZ_REGCON_RCC_EMPTY;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((obj0 == NULL) || (obj1 == NULL) ||
     (WlzIsEmpty(obj0, NULL) != 0) ||
     (WlzIsEmpty(obj1, NULL) != 0))
  {
    con = WLZ_REGCON_RCC_EMPTY;
  }
  else if((obj0->domain.core == NULL) || (obj1->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj0->type != obj1->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else
  {
    objI = WlzIntersect2(obj0, obj1, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if(WlzIsEmpty(objI, NULL))
      {
        /* DC || EC */
	objD = WlzDilation(obj0, WLZ_8_CONNECTED, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  if(!WlzIsEmpty(objD, NULL))
	  {
	    con = WLZ_REGCON_RCC_EC;
	  }
	}
      }
      else
      {
	/* EQ || PO || TPP || NTPP || TPPI || NTPPI */
	objU = WlzUnion2(obj0, obj1, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  objU0 = WlzDiffDomain(objU, obj0, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  objU1 = WlzDiffDomain(objU, obj1, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  int	oU0,
		oU1,
		msk;

	  oU0 = WlzIsEmpty(objU0, NULL)? 0: 1;
	  oU1 = WlzIsEmpty(objU1, NULL)? 0: 1;
	  msk = (oU1 << 1) | oU0;
	  if((msk > 0) && (dstNrmVol != NULL))
	  {
	    double volN,
	    	   volD;

	    volD = WlzVolume(objU, &errNum);
	    if((errNum == WLZ_ERR_NONE) && (volD < 0.5))
	    {
	      errNum = WLZ_ERR_DOMAIN_DATA;
	    }
	    else
	    {
	      volN =  WlzVolume(objI, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		vol = volN / volD;
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    switch(msk)
	    {
	      case 0:
		con = WLZ_REGCON_RCC_EQ;
		vol = 1.0;
		break;
	      case 1:
		/* TPP || NTPP */
		objD = WlzDilation(obj0, WLZ_8_CONNECTED, &errNum);
		if(errNum == WLZ_ERR_NONE)
		{
		  (void )WlzFreeObj(objU); objU = NULL;
		  (void )WlzFreeObj(objU0); objU0 = NULL;
		  objU = WlzUnion2(objD, obj1, &errNum);
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  objU0 = WlzDiffDomain(objU, obj0, &errNum);
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  if(WlzIsEmpty(objU0, NULL))
		  {
		    con = WLZ_REGCON_RCC_NTPP;
		  }
		  else
		  {
		    con = WLZ_REGCON_RCC_TPP;
		  }
		}
		break;
	      case 2:
		/* TPPI || NTPPI */
		objD = WlzDilation(obj1, WLZ_8_CONNECTED, &errNum);
		if(errNum == WLZ_ERR_NONE)
		{
		  (void )WlzFreeObj(objU); objU = NULL;
		  (void )WlzFreeObj(objU1); objU1 = NULL;
		  objU = WlzUnion2(objD, obj0, &errNum);
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  objU1 = WlzDiffDomain(objU, obj1, &errNum);
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  if(WlzIsEmpty(objU1, NULL))
		  {
		    con = WLZ_REGCON_RCC_NTPP;
		  }
		  else
		  {
		    con = WLZ_REGCON_RCC_TPP;
		  }
		}
		break;
	      case 3:
		con = WLZ_REGCON_RCC_PO;
		break;
	      default:
		break;
	    }
	  }
	}
      }
    }
  }
  (void )WlzFreeObj(objI);
  (void )WlzFreeObj(objD);
  (void )WlzFreeObj(objU);
  (void )WlzFreeObj(objU0);
  (void )WlzFreeObj(objU1);
  if((errNum == WLZ_ERR_NONE) && (dstNrmVol != NULL))
  {
    *dstNrmVol = vol;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(con);
}
