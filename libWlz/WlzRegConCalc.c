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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*!
* \enum		_WlzRCCTOIdx
* \ingroup	WlzBinaryOps
* \brief	Enumeration of temporary objects used in within
* 		WlzRegConCalcRCC(), WlzRCCMakeC() and WlzRCCMakeT().
* 		Typedef::WlzRCCTOIdx.
*/
typedef	enum	_WlzRCCTOIdx
{
  WLZ_RCCTOIDX_O0O1U	= 0,		/*!<\f$ o_0   \cup \o_1  \f$*/
  WLZ_RCCTOIDX_O0O1I,			/*!<\f$ o_0   \cap \o_1  \f$*/
  WLZ_RCCTOIDX_O0D,			/*!<\f$ o_0^+            \f$*/
  WLZ_RCCTOIDX_O1D,			/*!<\f$ o_1^+            \f$*/
  WLZ_RCCTOIDX_O0DO1U,			/*!<\f$ o_0^+ \cup o_1   \f$*/
  WLZ_RCCTOIDX_O0O1DU,			/*!<\f$ o_0   \cup o_1^+ \f$*/
  WLZ_RCCTOIDX_O0VO1I,                  /*!<\f$ o_0^v \cap o_1   \f$*/
  WLZ_RCCTOIDX_O0O1VI,			/*!<\f$ o_0   \cap o_1^v \f$*/
  WLZ_RCCTOIDX_CNT			/*!< Not an index but their number. */
} WlzRCCTOIdx;
#endif

static WlzErrorNum 		WlzRCCMakeC(
                                  WlzObject **c,
				  WlzObject **o,
                                  WlzObject **t,
				  int i);
static WlzErrorNum 		WlzRCCMakeT(
				  WlzObject **t,
				  WlzObject **o,
				  WlzRCCTOIdx i);

/*!
* \return	RCC classification of the given objects, ie object 0 is a
*               returned classification of object 1.
* \ingroup	WlzBinaryOps
* \brief	The given pair of spatial domain objects are classified
*		using a RCC.
*		For an explanation of RCC8 classifications
*		see the type definition ::WlzRCCClass and the paper:
*		D.A. Randell, etal,
*		"Discrete Mereotopology for Spatial Reasoning in
*		Automated Histological Image Analysis", PAMI 35(3) 2013.
*		The RCC8 has been extended to include encloses and
*		surrounds.
* 		The classification is performed using simple combinations
* 		of the Woolz union, intersection, exclusive or, dilation
* 		and convex hull operators on an ordered pair of
* 		spatial domains(\f$\Omega_0\f$ and \f$\Omega_1\f$):
*               \f{eqnarray*}{
                  C_0    &\leftarrow&
		    \Omega_0   \cap \Omega_1
		    \neq \emptyset \\
                  C_1    &\leftarrow&
		    \Omega_0^+ \cap \Omega_1
		    \neq \emptyset \\
                  C_2    &\leftarrow&
		    (\Omega_0   \oplus \Omega_1)
		    \neq \emptyset \\
                  C_3    &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus
		    (\Omega_0^+ \cup \Omega_1)
		    \neq \emptyset \\
                  C_4    &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus
		    \Omega_1 
		    \neq \emptyset \\
                  C_5    &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus
		    (\Omega_0   \cup \Omega_1^+) 
		    \neq \emptyset \\
                  C_6    &\leftarrow&
		    (\Omega_0^+ \cup \Omega_1)   \oplus
		    \Omega_1 
		    \neq \emptyset \\
                  C_7    &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus
		    \Omega_0 
		    \neq \emptyset \\
                  C_8    &\leftarrow&
		    (\Omega_0   \cup \Omega_1^+) \oplus
		    \Omega_0 
		    \neq \emptyset \\
		  C_9    &\leftarrow&
		    2|\Omega_0^v \cap \Omega_1|   \ge |\Omega_1| \\
		  C_{10} &\leftarrow&
		    2|\Omega_0   \cap \Omega_1^v| \ge |\Omega_0|
  		\f}
*		where
*		  are the \f$\cup\f$, \f$\cap\f$ and \f$\oplus\f$
*		  are the set union (logical or), intersection (logical and)
*		  and xor (logical exclusive or) operators;
*		  \f$\Omega^+\f$ indicates the dilation of \f$\Omega\f$,
*		  \f$\Omega^v\f$ the convex hull of \f$\Omega\f$ and
*		  \f$|\Omega|\f$ the cardinality (area or volume) of
*		  \f$\Omega\f$.
* 		The decision tree for the classification excluding
* 		enclosure is:
* 		\f[
C_0
\left\{
\begin{array}{ll}
0 & C_1
    \left\{
    \begin{array}{ll}
    0 & DC \\
      & \\
      & \\
    1 & C_3
	\left\{
	\begin{array}{ll}
	0 & SUR \\
	  & \\
	  & \\
	1 & C_5
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
1 & C_2
    \left\{
    \begin{array}{ll}
    0 & EQ \\
      & \\
      & \\
    1 & C_4
	\left\{
	\begin{array}{ll}
	0 & C_6
	    \left\{
	    \begin{array}{ll}
	    0 & NTPP \\
	      & \\
	    1 & TPP
	    \end{array}
	    \right. \\
	  & \\
	  & \\
	1 & C_7
	    \left\{
	    \begin{array}{ll}
	    0 & C_8
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
		  Basic Morphological Operations for the RCC Spatial
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
		  <td>\f$SUR(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$SURI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$ENC(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0|/|\Omega_0 \cup \Omega_1^v|\f$</td>
		</tr>
		<tr>
		  <td>\f$ENCI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0^v \cup \Omega_1|\f$</td>
		</tr>
		</table>
*		Many of the objects that are computed during the classification
*		are done so using a lazy evaluation with the functions
*		WlzRCCMakeC() and WlzRCCMakeT().
*		When computing the normalised volumes enclosure has the
*		higher priority.
* \param	obj0			First given spatial domain object.
* \param	obj1			Second given spatial domain object.
* \param	noEnc			Don't include enclosure if non-zero.
* \param	dstNrmVol		Destination pointer for the normalized
* 					volume (see above), may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzRCCClass 	WlzRegConCalcRCC(WlzObject *obj0, WlzObject *obj1, int noEnc,
				 double *dstNrmVol, WlzErrorNum *dstErr)
{
  int 		i;
  double	nrmVol = 0.0;		/* Normalized volume */
  WlzLong	u01 = 0; 		/* |\Omega_0 \cup \Omega_1| */
  WlzLong	u[2] = {0},		/* |\Omega_i|, i \in 0 \cdots 1 */
  		v[2] = {0};		/* |c_9|, |c_{10}| */
  WlzObject	*c[8] = {NULL},		/* c_i, i \in 0 \cdots 8 */
		*o[2] = {NULL},		/* \Omega_i, i \in 0 \cdots 1 */
		*t[WLZ_RCCTOIDX_CNT] = {NULL}; /* Temporary object as
					in the enum WlzRCCTOIdx. */
  WlzValues	nullValues;
  WlzRCCClass	cls = WLZ_RCC_EMPTY; /* Classification mask. */
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nullValues.core = NULL;
  /* Compute classification using the decision tree. */
  if((obj0 == NULL) || (obj1 == NULL) ||
     (WlzIsEmpty(obj0, NULL) != 0) ||
     (WlzIsEmpty(obj1, NULL) != 0))
  {
    cls = WLZ_RCC_EMPTY;
  }
  else if((obj0->domain.core == NULL) || (obj1->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((obj0->type != obj1->type) ||
          ((obj0->type != WLZ_2D_DOMAINOBJ) &&
	   (obj0->type != WLZ_3D_DOMAINOBJ)))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(((o[0] = WlzMakeMain(obj0->type, obj0->domain, nullValues,
  			       NULL, NULL, &errNum)) != NULL) &&
          ((o[1] = WlzMakeMain(obj1->type, obj1->domain, nullValues,
	  		       NULL, NULL, &errNum)) != NULL))
  {
    errNum = WlzRCCMakeC(c, o, t, 0);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(WlzIsEmpty(c[0], NULL))
    {
      errNum = WlzRCCMakeC(c, o, t, 1);
      if(errNum == WLZ_ERR_NONE)
      {
	if(WlzIsEmpty(c[1], NULL))
	{
	  cls = WLZ_RCC_DC;
	}
	else
	{
	  errNum = WlzRCCMakeC(c, o, t, 3);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(WlzIsEmpty(c[3], NULL))
	    {
	      cls = WLZ_RCC_SUR;
	    }
	    else
	    {
	      errNum = WlzRCCMakeC(c, o, t, 5);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(WlzIsEmpty(c[5], NULL))
		{
		  cls = WLZ_RCC_SURI;
		}
		else
		{
		  cls = WLZ_RCC_EC;
		}
	      }
	    }
	  }
	}
      }
    }
    else
    {
      errNum = WlzRCCMakeC(c, o, t, 2);
      if(errNum == WLZ_ERR_NONE)
      {
        if(WlzIsEmpty(c[2], NULL))
	{
	  cls = WLZ_RCC_EQ;
	}
	else
	{
	  errNum = WlzRCCMakeC(c, o, t, 4);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(WlzIsEmpty(c[4], NULL))
	    {
	      errNum = WlzRCCMakeC(c, o, t, 6);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(WlzIsEmpty(c[6], NULL))
		{
		  cls = WLZ_RCC_NTPP;
		}
		else
		{
		  cls = WLZ_RCC_TPP;
		}
	      }
	    }
	    else
	    {
	      errNum = WlzRCCMakeC(c, o, t, 7);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(WlzIsEmpty(c[7], NULL))
		{
		  errNum = WlzRCCMakeC(c, o, t, 8);
		  if(errNum == WLZ_ERR_NONE)
		  {
		    if(WlzIsEmpty(c[8], NULL))
		    {
		      cls = WLZ_RCC_NTPPI;
		    }
		    else
		    {
		      cls = WLZ_RCC_TPPI;
		    }
		  }
		}
		else
		{
		  cls = WLZ_RCC_PO;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  /* If enclosure is required check for it and add to classification mask. */
  if((errNum == WLZ_ERR_NONE) && (noEnc == 0) &&
     ((cls &
       (WLZ_RCC_EQ | WLZ_RCC_SUR | WLZ_RCC_SURI)) == 0))
  {
    for(i = 0; i <= 1; ++i)
    {
      errNum = WlzRCCMakeT(t, o,
                           (i == 0)? WLZ_RCCTOIDX_O0VO1I: WLZ_RCCTOIDX_O0O1VI);
      if(errNum == WLZ_ERR_NONE)
      {
        u[i] = WlzVolume(o[i], &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        v[i] = WlzVolume(t[(i == 0)?
	                 WLZ_RCCTOIDX_O0VO1I: WLZ_RCCTOIDX_O0O1VI],
			 &errNum);
	if((errNum == WLZ_ERR_NONE) && (v[i] < 1))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if((2 * v[0]) >= u[0])
      {
        cls |= WLZ_RCC_ENC;
      }
      if((2 * v[1]) >= u[1])
      {
        cls |= WLZ_RCC_ENCI;
      }
    }
  }
  /* Compute the maximum normalized volume for the classification(s) in the
   * classification mask. */
  if((errNum == WLZ_ERR_NONE) && (dstNrmVol != NULL))
  {
    unsigned int m;

    for(m = 1; m <= WLZ_RCC_MSK; m <<= 1)
    {
      if(m & cls)
      {
	double	nV = 0.0;

        switch(cls)
	{
	  case WLZ_RCC_EQ:
	    nV = 1.0;
	    break;
	  case WLZ_RCC_PO:
	    /* |\Omega_0 \cap \Omega_1| / |\Omega_0 \cup \Omega_1|  =
	     * u_0 / u_01 */
	    if(u[0] <= 0)
	    {
	      errNum = WlzRCCMakeC(c, o, t, 0);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        u[0] = WlzVolume(c[0], &errNum);
	      }
	    }
	    if((errNum == WLZ_ERR_NONE) && (u01 <= 0))
	    {
	      errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1U);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        u01 = WlzVolume(t[WLZ_RCCTOIDX_O0O1U], &errNum);
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      nV = (double )(u[0]) / (double )u01;
            }
	    break;
	  case WLZ_RCC_SUR: /* FALLTHROUGH */
	  case WLZ_RCC_TPP: /* FALLTHROUGH */
	  case WLZ_RCC_NTPP:
	    /* |\Omega_0| / |\Omega_0 \cup \Omega_1|  =
	     * u_0 / u_01 */
	    if(u[0] <= 0)
	    {
	      u[0] = WlzVolume(o[0], &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(u01 <= 0)
	      {
	        errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1U);
		if(errNum == WLZ_ERR_NONE)
		{
		  u01 = WlzVolume(t[WLZ_RCCTOIDX_O0O1U], &errNum);
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      nV = (double )(u[0]) / (double )u01;
	    }
	    break;
	  case WLZ_RCC_SURI: /* FALLTHROUGH */
	  case WLZ_RCC_TPPI: /* FALLTHROUGH */
	  case WLZ_RCC_NTPPI:
	    /* |\Omega_1| / |\Omega_0 \cup \Omega_1|  =
	     * u_1 / u_01 */
	    if(u[1] <= 0)
	    {
	      u[1] = WlzVolume(o[1], &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(u01 <= 0)
	      {
		errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1U);
		if(errNum == WLZ_ERR_NONE)
		{
		  u01 = WlzVolume(t[WLZ_RCCTOIDX_O0O1U], &errNum);
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      nV = (double )(u[1]) / (double )u01;
	    }
	    break;
	  case WLZ_RCC_ENC:
	    /* |\Omega_0|/|\Omega_0 \cup \Omega_1^v| =
	     * u_0 / v_1 */
	    nV = (double )(u[0]) / (double )(v[1]);
	    break;
	  case WLZ_RCC_ENCI:
	    /* |\Omega_1|/|\Omega_0^v \cup \Omega_1| =
	     * u_1 / v_0 */
	    nV = (double )(u[1]) / (double )(v[0]);
	    break;
	  default:
	    break;
	}
        if((errNum == WLZ_ERR_NONE) && (nV > nrmVol))
	{
	  nrmVol = nV;
	}
      }
    }
    *dstNrmVol = nrmVol;
  }
  /* Free objects. */
  for(i = 0; i < WLZ_RCCTOIDX_CNT; ++i)
  {
    (void )WlzFreeObj(t[i]);
  }
  for(i = 0; i <= 8; ++i)
  {
    (void )WlzFreeObj(c[i]);
  }
  for(i = 0; i < 2; ++i)
  {
    (void )WlzFreeObj(o[i]);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cls);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	Test to see if the classification object \f$c_i\f$ as in
* 		WlzRegConCalcRCC() has been computed and if not do so.
* \param	c			Array c as in WlzRegConCalcRCC().
* \param	o			Array of objects, o[0] and o[1].
* \param	t			Array of temporary objects as in
* 					WlzRCCTOIdx.
* \param	i			Index for the classification object.
*/
static WlzErrorNum  WlzRCCMakeC(WlzObject **c, WlzObject **o, WlzObject **t,
				int i)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((i >= 0) && (i <= 8) && (c[i] == NULL))
  {
    switch(i)
    {
      case 0: /* \Omega_0 \cap \Omega_1 */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1I);
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(t[WLZ_RCCTOIDX_O0O1I], NULL);
	}
	break;
      case 1: /* \Omega_1^+ \cap \Omega_1 */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0D);
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(
	         WlzIntersect2(t[WLZ_RCCTOIDX_O0D], o[1], &errNum), NULL);
	}
	break;
      case 2: /* \Omega_0 \oplus \Omega_1 */
	c[i] = WlzAssignObject(
	       WlzXORDom(o[0], o[1], &errNum), NULL);
	break;
      case 3: /* (\Omega_0^+ \cup \Omega_1) \oplus (\Omega_0 \cup \Omega_1) */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0DO1U);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1U);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(
	         WlzXORDom(t[WLZ_RCCTOIDX_O0DO1U],
		           t[WLZ_RCCTOIDX_O0O1U ], &errNum), NULL);
	}
	break;
      case 4: /* (\Omega_0 \cup \Omega_1) \oplus \Omega_1 */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1U);
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(
	         WlzXORDom(t[WLZ_RCCTOIDX_O0O1U],
		           o[1], &errNum), NULL);
	}
	break;
      case 5: /* (\Omega_0 \cup \Omega_1^+) \oplus (\Omega_0 \cup \Omega_1) */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1U);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1DU);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(
	         WlzXORDom(t[WLZ_RCCTOIDX_O0O1U], t[WLZ_RCCTOIDX_O0O1DU],
		           &errNum), NULL);
	}
	break;
      case 6: /* (\Omega_0^+ \cup \Omega_1) \oplus \Omega_1 */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0DO1U);
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(
	  	 WlzXORDom(t[WLZ_RCCTOIDX_O0DO1U], o[1], &errNum), NULL);
	}
	break;
      case 7: /* (\Omega_0 \cup \Omega_1) \oplus \Omega_0 */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1U);
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(
	  	 WlzXORDom(o[0], t[WLZ_RCCTOIDX_O0O1U], &errNum), NULL);
	}
	break;
      case 8: /* (\Omega_0 \cup \Omega_1^+) \oplus \Omega_0 */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0O1DU);
	if(errNum == WLZ_ERR_NONE)
	{
	  c[i] = WlzAssignObject(
	  	 WlzXORDom(o[0], t[WLZ_RCCTOIDX_O0O1DU], &errNum), NULL);
	}
	break;
      default:
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	Make temporary objects as in WlzRCCTOIdx for
* 		WlzRegConCalcRCC().
* \param	t			Array of temporary objects.
* \param	o			Array of objects, o[0] and o[1].
* \param	i			Temporary object index.
*/
static WlzErrorNum WlzRCCMakeT(WlzObject **t, WlzObject **o, WlzRCCTOIdx i)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((i >= 0) && (i < WLZ_RCCTOIDX_CNT) && (t[i] == NULL))
  {
    switch(i)
    {
      case WLZ_RCCTOIDX_O0O1U:			/* o_0 \cup \o_1 */
        t[i] = WlzAssignObject(
	       WlzUnion2(o[0], o[1], &errNum), NULL);
        break;
      case WLZ_RCCTOIDX_O0O1I:			/* o_0 \cap \o_1 */
	t[i] = WlzAssignObject(
	       WlzIntersect2(o[0], o[1], &errNum), NULL);
        break;
      case WLZ_RCCTOIDX_O0D:			/* o_0^+ */
      t[i] = WlzAssignObject(
      	     WlzDilation(o[0],
			 (o[0]->type == WLZ_2D_DOMAINOBJ)?
			 WLZ_8_CONNECTED: WLZ_26_CONNECTED,
			 &errNum), NULL);
        break;
      case WLZ_RCCTOIDX_O1D:			/* o_1^+ */
      t[i] = WlzAssignObject(
             WlzDilation(o[1],
			 (o[1]->type == WLZ_2D_DOMAINOBJ)?
			 WLZ_8_CONNECTED: WLZ_26_CONNECTED,
			 &errNum), NULL);
        break;
      case WLZ_RCCTOIDX_O0DO1U:			/* o_0^+ \cup o_1 */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O0D);
	if(errNum == WLZ_ERR_NONE)
	{
	  t[i] = WlzAssignObject(
	         WlzUnion2(t[WLZ_RCCTOIDX_O0D], o[1], &errNum), NULL);
	}
        break;
      case WLZ_RCCTOIDX_O0O1DU:			/* o_0   \cup o_1^+ */
	errNum = WlzRCCMakeT(t, o, WLZ_RCCTOIDX_O1D);
	if(errNum == WLZ_ERR_NONE)
	{
	  t[i] = WlzAssignObject(
	         WlzUnion2(o[0], t[WLZ_RCCTOIDX_O1D], &errNum), NULL);
	}
        break;
      case WLZ_RCCTOIDX_O0VO1I:                   /* o_0^v \cap o_1   */
      case WLZ_RCCTOIDX_O0O1VI: /* FALLTHROUGH */ /* o_0   \cap o_1^v */
	{
	  int		i0,
	  		i1;
	  WlzObject	*x;

	  i0 = (i == WLZ_RCCTOIDX_O0VO1I)? 0: 1;
	  x = WlzObjToConvexHullSpDomObj(o[i0], &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    t[i] = WlzAssignObject(
	           WlzIntersect2(t, o[!i0], &errNum), NULL);
	  }
	  (void )WlzFreeObj(t);
	}
        break;
      default:
        break;
    }
  }
  return(errNum);
}
