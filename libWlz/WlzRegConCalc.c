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

#include <string.h>
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
  WLZ_RCCTOIDX_O0O1U	= 0,		/*!<\f$ o_0     \cup \o_1        \f$*/
  WLZ_RCCTOIDX_O0O1I,			/*!<\f$ o_0     \cap \o_1        \f$*/
  WLZ_RCCTOIDX_O0D,			/*!<\f$ o_0^+                    \f$*/
  WLZ_RCCTOIDX_O1D,			/*!<\f$ o_1^+                    \f$*/
  WLZ_RCCTOIDX_O0F,			/*!<\f$ o_0^{\bullet}            \f$*/
  WLZ_RCCTOIDX_O1F,			/*!<\f$ o_1^{\bullet}            \f$*/
  WLZ_RCCTOIDX_O0DO1U,			/*!<\f$ o_0^+   \cup o_1         \f$*/
  WLZ_RCCTOIDX_O0O1DU,			/*!<\f$ o_0     \cup o_1^+       \f$*/
  WLZ_RCCTOIDX_O0CO1I,                  /*!<\f$ o_0^{\circ} \cap o_1     \f$*/
  WLZ_RCCTOIDX_O0O1CI,			/*!<\f$ o_0     \cap o_1^{\circ} \f$*/
  WLZ_RCCTOIDX_O0FO1U,                  /*!<\f$ o_0^{\bullet} \cup o_1   \f$*/
  WLZ_RCCTOIDX_O0O1FU,                  /*!<\f$ o_0 \cup o_1^{\bullet}   \f$*/
  WLZ_RCCTOIDX_CNT			/*!< Not an index but their number. */
} WlzRCCTOIdx;
#endif

static WlzErrorNum 		WlzRCCMakeC(
				  WlzObject **o,
                                  WlzObject **c,
                                  WlzObject **t,
				  int i);
static WlzErrorNum 		WlzRCCMakeT(
				  WlzObject **o,
				  WlzObject **t,
				  WlzRCCTOIdx i);
static WlzErrorNum		WlzRCCOffset(
				  WlzObject **o,
				  WlzObject **t,
				  int maxDist,
				  double *dQ0,
				  double *dQ1,
				  double *dQ2);
static WlzErrorNum	WlzRCCCompDistHist2D(
			  int maxDist,
			  int *dHist,
			  WlzObject *dobj);
static WlzErrorNum	WlzRCCCompDistHist3D(
			  int maxDist,
			  int *dHist,
			  WlzObject *dobj);

/*!
* \return	RCC classification of the given objects, ie object 0 is a
*               returned classification of object 1.
* \ingroup	WlzBinaryOps
* \brief	The given pair of spatial domain objects are classified
*		using a RCC with optional enclosure and offset classifications.
*
*		For an explanation of RCC8 classifications
*		see the type definition ::WlzRCCClass and the paper:
*		D.A. Randell, etal,
*		"Discrete Mereotopology for Spatial Reasoning in
*		Automated Histological Image Analysis", PAMI 35(3) 2013.
*		The RCC8 has been extended to include both tangential and
*		non-tangential surrounds.
*
* 		The RCC classification is performed using simple combinations
* 		of the Woolz union, intersection, exclusive or, dilation,
* 		fill and convex hull operators on an ordered pair of
* 		spatial domains(\f$\Omega_0\f$ and \f$\Omega_1\f$):
*               \f{eqnarray*}{
                  C_0 &\leftarrow&
		    \Omega_0   \cap \Omega_1
		    \neq \emptyset \\
                  C_1 &\leftarrow&
		    \Omega_0^+ \cap \Omega_1
		    \neq \emptyset \\
                  C_2 &\leftarrow&
		    (\Omega_0   \oplus \Omega_1)
		    \neq \emptyset \\
                  C_3 &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus
		    \Omega_1 
		    \neq \emptyset \\
                  C_4 &\leftarrow&
		    (\Omega_0^+ \cup \Omega_1)   \oplus
		    \Omega_1 
		    \neq \emptyset \\
                  C_5 &\leftarrow&
		    (\Omega_0   \cup \Omega_1)   \oplus
		    \Omega_0 
		    \neq \emptyset \\
                  C_6 &\leftarrow&
		    (\Omega_0   \cup \Omega_1^+) \oplus
		    \Omega_0 
		    \neq \emptyset \\
		  C_7 &\leftarrow&
		    (\Omega_0^{\bullet} \cup \Omega_1) \oplus
		    \Omega_0^{\bullet}
		    \neq \emptyset \\
		  C_8 &\leftarrow&
		    (\Omega_0 \cup \Omega_1^{\bullet}) \oplus
		    \Omega_1^{\bullet}
		    \neq \emptyset \\
		  C_9 &\leftarrow&
		    2|\Omega_0 \cap \Omega_1^{\circ}|   \ge |\Omega_0| \\
		  C_{10} &\leftarrow&
		    2|\Omega_0^{\circ}   \cap \Omega_1| \ge |\Omega_1|
  		\f}
*		where
*		  are the \f$\cup\f$, \f$\cap\f$ and \f$\oplus\f$
*		  are the set union (logical or), intersection (logical and)
*		  and xor (logical exclusive or) operators;
*		  \f$\Omega^+\f$ indicates the dilation of \f$\Omega\f$,
*		  \f$\Omega^{\circ}\f$ the convex hull of \f$\Omega\f$,
*		  \f$\Omega^{\bullet}\f$ indicates \f$\Omega\f$ filled and
*		  \f$|\Omega|\f$ the cardinality (area or volume) of
*		  \f$\Omega\f$.
* 		The decision tree for the classification excluding
* 		enclosure and offset is:
* 		\f[
C_0
\left\{
\begin{array}{ll}
0 & C_1
    \left\{
    \begin{array}{ll}
    0 & C_7
	\left\{
	\begin{array}{ll}
	0 & NTSURI \\
	  & \\
	  & \\
	1 & C_8
	    \left\{
	    \begin{array}{ll}
	    0 & NTSUR \\
	      & \\
	    1 & DC
	    \end{array}
	    \right. \\
	\end{array}
	\right. \\
      & \\
      & \\
    1 & C_7
	\left\{
	\begin{array}{ll}
	0 & TSURI \\
	  & \\
	  & \\
	1 & C_8
	    \left\{
	    \begin{array}{ll}
	    0 & TSUR \\
	      & \\
	    1 & EC
	    \end{array}
	    \right. \\
	\end{array}
	\right. \\
    \end{array}
    \right. \\
  & \\
  & \\
1 & C_2
    \left\{
    \begin{array}{ll}
    0 & EQ \\
      & \\
      & \\
    1 & C_3
	\left\{
	\begin{array}{ll}
	0 & C_4
	    \left\{
	    \begin{array}{ll}
	    0 & NTPP \\
	      & \\
	    1 & TPP
	    \end{array}
	    \right. \\
	  & \\
	  & \\
	1 & C_5
	    \left\{
	    \begin{array}{ll}
	    0 & C_6
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
*		The statistics are computed for each classification
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
		  <td>\f$TSUR(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$TSURI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$NTSUR(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$NTSURI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_1|/|\Omega_0 \cup \Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$ENC(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0 \cap \Omega_1^{\circ}|/|\Omega_0|\f$</td>
		</tr>
		<tr>
		  <td>\f$ENCI(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$|\Omega_0^{\circ} \cap \Omega_1|/|\Omega_1|\f$</td>
		</tr>
		<tr>
		  <td>\f$OST(\Omega_0,\Omega_1)\f$</td>
		  <td>\f$q_1/(q_1 + q_2 - q_0)\f$</td>
		</tr>
		</table>
*
*		Many of the objects that are computed during the
*		classification are done so using a lazy evaluation with the
*		functions WlzRCCMakeC() and WlzRCCMakeT().
*
*		Enclosure and offset are somwhat more expensive to compute
*		than the other classifications, for this reason and because
*		they are not strictly part of a RCC they can be avoided by
*		setting the noEnc or noOst flags.
*
*		Enclosure will be computed if the noEnc has not been set and
*		the classification is not one of WLZ_RCC_EQ, WLZ_RCC_TSUR,
*		WLZ_RCC_TSURI, WLZ_RCC_NTSUR or WLZ_RCC_NTSURI.
*		Enclosure is computed using:
*		\f[
		  |\Omega_0 \cap \Omega_1^{\circ}|/|\Omega_0|
		\f]
*		for \f$\Omega_0\f$ to be encloded by \f$\Omega_1\f$ then
*		at least half of \f$\Omega_0\f$ must intersect the convex
*		hull of \f$\Omega_1\f$.
*
*		Offset will be computed if the noOst parameter has not
*		been set and the classification is not WLZ_RCC_EQ.
*		Offset is computed within a restricted domain in which
*		all pixels/voxels are equidistant for the domains of the
*		given objects:
		\f[
		\Omega_e = (\Omega_0 \cup \Omega_1)^\circ \cap
		           \Omega_0^{+d_{max}} \cap
			   \Omega_1^{+d_{max}} \cap
			   \Omega(D(\Omega_0) = D(\Omega_1))
		\f]
*		where \f$D(\Omega)\f$ is the distance transform of the domain
*		\f$\Omega\f$. Within \f$\Omega_e\f$ the first, second and
*		third quantiles (\f$q_0\f$, \f$q_1\f$ and \f$q_2\f$) of
*		the distances \f$D(\Omega_0)\f$ (or equivalently
*		\f$D(\Omega_1)\f$) are computed. The ratio of the median
*		to the median plus interquartile range is then computed
*		and the domains are classified as offset if this ratio
*		is greater than or equal to one half:
*		\f[
		\frac{q_1}{q_1 + q_2 - q_0} \geq 0.5
		\f]
* \param	obj0			First given spatial domain object.
* \param	obj1			Second given spatial domain object.
* \param	noEnc			Don't include enclosure if non-zero.
* \param	noOst			Don't include offset if non-zero.
* \param	maxOstDist		Maximum distance for offset, not
* 					used if noOst is non-zero.
* \param	dstStatCnt		Destination pointer for the number
* 					of elements returned in the array of
* 					statistics (see above), may be NULL.
* 					Ignored if dstStatAry is NULL.
* \param	dstStatAry		Destination pointer for an array of
* 					statistics (see above), may be NULL.
* 					If an array is returned it should be
* 					freed using AlcFree().
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzRCCClass 			WlzRegConCalcRCC(
				  WlzObject *obj0,
				  WlzObject *obj1,
				 int noEnc,
				 int noOst,
				 int maxOstDist,
				 int *dstStatCnt,
				 double **dstStatAry,
				 WlzErrorNum *dstErr)
{
  int 		i;
  WlzLong	i01 = 0, 		/* |\Omega_0 \cap \Omega_1| */
  		u01 = 0; 		/* |\Omega_0 \cup \Omega_1| */
  WlzLong	u[2] = {0},		/* |\Omega_i|, i \in 0 \cdots 1 */
  		v[2] = {0};		/* |c_9|, |c_{10}| */
  WlzObject	*c[11] = {NULL},		/* c_i, i \in 0 \cdots 10 */
		*o[2] = {NULL},		/* \Omega_i, i \in 0 \cdots 1 */
		*t[WLZ_RCCTOIDX_CNT] = {NULL}; /* Temporary object as
					in the enum WlzRCCTOIdx. */
  double	stats[WLZ_RCCIDX_CNT] = {0.0}; /* Classification statistics. */
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
  else if(((o[0] = WlzAssignObject(
                   WlzMakeMain(obj0->type, obj0->domain, nullValues,
  			       NULL, NULL, &errNum), NULL)) != NULL) &&
          ((o[1] = WlzAssignObject(
		   WlzMakeMain(obj1->type, obj1->domain, nullValues,
	  		       NULL, NULL, &errNum), NULL)) != NULL))
  {
    errNum = WlzRCCMakeC(o, c, t, 0);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(WlzIsEmpty(c[0], NULL))
    {
      errNum = WlzRCCMakeC(o, c, t, 1);
      if(errNum == WLZ_ERR_NONE)
      {
	if(WlzIsEmpty(c[1], NULL))
	{
	  errNum = WlzRCCMakeC(o, c, t, 7);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(WlzIsEmpty(c[7], NULL))
	    {
	      cls = WLZ_RCC_NTSURI;
	    }
	    else
	    {
	      errNum = WlzRCCMakeC(o, c, t, 8);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(WlzIsEmpty(c[8], NULL))
		{
		  cls = WLZ_RCC_NTSUR;
		}
		else
		{
		  cls = WLZ_RCC_DC;
		}
	      }
	    }
	  }
	}
	else
	{
	  errNum = WlzRCCMakeC(o, c, t, 7);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(WlzIsEmpty(c[7], NULL))
	    {
	      cls = WLZ_RCC_TSURI;
	    }
	    else
	    {
	      errNum = WlzRCCMakeC(o, c, t, 8);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(WlzIsEmpty(c[8], NULL))
		{
		  cls = WLZ_RCC_TSUR;
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
      errNum = WlzRCCMakeC(o, c, t, 2);
      if(errNum == WLZ_ERR_NONE)
      {
        if(WlzIsEmpty(c[2], NULL))
	{
	  cls = WLZ_RCC_EQ;
	}
	else
	{
	  errNum = WlzRCCMakeC(o, c, t, 3);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(WlzIsEmpty(c[3], NULL))
	    {
	      errNum = WlzRCCMakeC(o, c, t, 4);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(WlzIsEmpty(c[4], NULL))
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
	      errNum = WlzRCCMakeC(o, c, t, 5);
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(WlzIsEmpty(c[5], NULL))
		{
		  errNum = WlzRCCMakeC(o, c, t, 6);
		  if(errNum == WLZ_ERR_NONE)
		  {
		    if(WlzIsEmpty(c[6], NULL))
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
       (WLZ_RCC_EQ |
        WLZ_RCC_TSUR | WLZ_RCC_TSURI |
	WLZ_RCC_NTSUR | WLZ_RCC_NTSURI)) == 0))
  {
    for(i = 0; i <= 1; ++i)
    {
      errNum = WlzRCCMakeT(o, t,
                           (i == 0)? WLZ_RCCTOIDX_O0O1CI: WLZ_RCCTOIDX_O0CO1I);
      if(errNum == WLZ_ERR_NONE)
      {
        u[i] = WlzVolume(o[i], &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        v[i] = WlzVolume(t[(i == 0)?
	                 WLZ_RCCTOIDX_O0O1CI: WLZ_RCCTOIDX_O0CO1I],
			 &errNum);
	if((errNum == WLZ_ERR_NONE) && (v[i] < 0))
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
  if((errNum == WLZ_ERR_NONE) && (dstStatAry != NULL))
  {
    int 	i,
    		m;

    for(i = 0; i < WLZ_RCCIDX_CNT; ++i)
    {
      m = 1<<i;
      if(m & cls)
      {
	double	s = 0.0;

        switch(m)
	{
	  case WLZ_RCC_EQ:
	    s = 1.0;
	    break;
	  case WLZ_RCC_PO:
	    /* |\Omega_0 \cap \Omega_1| / |\Omega_0 \cup \Omega_1|  =
	     * u_0 / u_01 */
	    if(i01 <= 0)
	    {
	      errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1I);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        i01 = WlzVolume(t[WLZ_RCCTOIDX_O0O1I], &errNum);
	      }
	    }
	    if((errNum == WLZ_ERR_NONE) && (u01 <= 0))
	    {
	      errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1U);
	      if(errNum == WLZ_ERR_NONE)
	      {
	        u01 = WlzVolume(t[WLZ_RCCTOIDX_O0O1U], &errNum);
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      s = (double )(i01) / (double )u01;
            }
	    break;
	  case WLZ_RCC_TSUR: /* FALLTHROUGH */
	  case WLZ_RCC_NTSUR: /* FALLTHROUGH */
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
	        errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1U);
		if(errNum == WLZ_ERR_NONE)
		{
		  u01 = WlzVolume(t[WLZ_RCCTOIDX_O0O1U], &errNum);
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      s = (double )(u[0]) / (double )u01;
	    }
	    break;
	  case WLZ_RCC_TSURI: /* FALLTHROUGH */
	  case WLZ_RCC_NTSURI: /* FALLTHROUGH */
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
		errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1U);
		if(errNum == WLZ_ERR_NONE)
		{
		  u01 = WlzVolume(t[WLZ_RCCTOIDX_O0O1U], &errNum);
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      s = (double )(u[1]) / (double )u01;
	    }
	    break;
	  case WLZ_RCC_ENC:
	    /* |\Omega_0 \cup \Omega_1^{\circ}|/|\Omega_0| =
	     * v_0 / u_0 */
	    if(u[1] >= 0)
	    {
	      s = (double )(v[0]) / (double )(u[0]);
	    }
	    break;
	  case WLZ_RCC_ENCI:
	    /* |\Omega_0^{\circ} \cup \Omega_1|/|\Omega_1| =
	     * v_1 / u_1 */
	    if(v[1] >= 0)
	    {
	      s = (double )(v[1]) / (double )(u[1]);
	    }
	    break;
	  default:
	    break;
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  stats[i] = s;
	}
      }
    }
  }
  /* If offset is required check for it and add to both the classification
   * mask and statistics. */
  if((errNum == WLZ_ERR_NONE) && (noOst == 0) &&
     ((cls & WLZ_RCC_EQ) == 0))
  {
    double		ostQ[3];

    errNum = WlzRCCOffset(o, t,
                          maxOstDist, &(ostQ[0]), &(ostQ[1]), &(ostQ[2]));
    if(errNum == WLZ_ERR_NONE)
    {
      if((ostQ[1] > 0) && (ostQ[1] < maxOstDist) && (ostQ[2] > ostQ[0]))
      {
	const double eps = 1.0e-06;

	if(ostQ[2] > ostQ[0])
	{
	  stats[WLZ_RCCIDX_OST] = ostQ[1] / (ostQ[2] + ostQ[1] - ostQ[0]);
	}
	else
	{
	  stats[WLZ_RCCIDX_OST] = 1.0;
	}
	if(stats[WLZ_RCCIDX_OST] > (0.5 - eps))
	{
	  cls |= WLZ_RCC_OST;
	}
      }
    }
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
  if((errNum == WLZ_ERR_NONE) && (dstStatAry != NULL))
  {
    if((*dstStatAry = (double *)
                        AlcMalloc(sizeof(double) * WLZ_RCCIDX_CNT)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      (void )memcpy(*dstStatAry, stats, sizeof(double) * WLZ_RCCIDX_CNT);
      if(dstStatCnt)
      {
        *dstStatCnt = WLZ_RCCIDX_CNT;
      }
    }
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
* \param	o			Array of objects, o[0] and o[1].
* \param	c			Array c as in WlzRegConCalcRCC().
* \param	t			Array of temporary objects as in
* 					WlzRCCTOIdx.
* \param	i			Index for the classification object.
*/
static WlzErrorNum  		WlzRCCMakeC(
				  WlzObject **o,
				  WlzObject **c,
				  WlzObject **t,
				int i)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((i >= 0) && (i <= 8))
  {
    if(c[i] == NULL)
    {
      switch(i)
      {
	case 0: /* \Omega_0 \cap \Omega_1 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1I);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    c[i] = WlzAssignObject(t[WLZ_RCCTOIDX_O0O1I], NULL);
	  }
	  break;
	case 1: /* \Omega_1^+ \cap \Omega_1 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0D);
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
	case 3: /* (\Omega_0 \cup \Omega_1) \oplus \Omega_1 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1U);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    c[i] = WlzAssignObject(
		   WlzXORDom(t[WLZ_RCCTOIDX_O0O1U],
			     o[1], &errNum), NULL);
	  }
	  break;
	case 4: /* (\Omega_0^+ \cup \Omega_1) \oplus \Omega_1 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0DO1U);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    c[i] = WlzAssignObject(
		   WlzXORDom(t[WLZ_RCCTOIDX_O0DO1U], o[1], &errNum), NULL);
	  }
	  break;
	case 5: /* (\Omega_0 \cup \Omega_1) \oplus \Omega_0 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1U);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    c[i] = WlzAssignObject(
		   WlzXORDom(o[0], t[WLZ_RCCTOIDX_O0O1U], &errNum), NULL);
	  }
	  break;
	case 6: /* (\Omega_0 \cup \Omega_1^+) \oplus \Omega_0 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1DU);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    c[i] = WlzAssignObject(
		   WlzXORDom(o[0], t[WLZ_RCCTOIDX_O0O1DU], &errNum), NULL);
	  }
	  break;
	case 7:/*(\Omega_0^{\bullet} \cup \Omega_1) \oplus \Omega_0^{\bullet}*/
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0F);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0FO1U);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    c[i] = WlzAssignObject(
		   WlzXORDom(t[WLZ_RCCTOIDX_O0FO1U], t[WLZ_RCCTOIDX_O0F],
			     &errNum), NULL);
	  }
	  break;
	case 8:/*(\Omega_0 \cup \Omega_1^{\bullet}) \oplus \Omega_1^{\bullet}*/
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O1F);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1FU);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    c[i] = WlzAssignObject(
		   WlzXORDom(t[WLZ_RCCTOIDX_O0O1FU], t[WLZ_RCCTOIDX_O1F],
			     &errNum), NULL);
	  }
	  break;
	default:
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
      }
    }
  }
  else
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	Make temporary objects as in WlzRCCTOIdx for
* 		WlzRegConCalcRCC().
* \param	o			Array of objects, o[0] and o[1].
* \param	t			Array of temporary objects.
* \param	i			Temporary object index.
*/
static WlzErrorNum 		WlzRCCMakeT(
				  WlzObject **o,
				  WlzObject **t,
				  WlzRCCTOIdx i)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((i >= 0) && (i < WLZ_RCCTOIDX_CNT))
  {
    if(t[i] == NULL)
    {
      WlzConnectType con;

      con = (o[0]->type == WLZ_2D_DOMAINOBJ)?
            WLZ_8_CONNECTED: WLZ_26_CONNECTED;
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
		 WlzDilation(o[0], con, &errNum), NULL);
	  break;
	case WLZ_RCCTOIDX_O1D:			/* o_1^+ */
	  t[i] = WlzAssignObject(
		 WlzDilation(o[1], con, &errNum), NULL);
	  break;

	case WLZ_RCCTOIDX_O0F:			/* o_0^{\bullet} */
	  t[i] = WlzAssignObject(
		 WlzDomainFill(o[0], &errNum), NULL);
	  break;
	case WLZ_RCCTOIDX_O1F:			/* o_1^{\bullet} */
	  t[i] = WlzAssignObject(
		 WlzDomainFill(o[1], &errNum), NULL);
	  break;

	case WLZ_RCCTOIDX_O0DO1U:		/* o_0^+ \cup o_1 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0D);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    t[i] = WlzAssignObject(
		   WlzUnion2(t[WLZ_RCCTOIDX_O0D], o[1], &errNum), NULL);
	  }
	  break;
	case WLZ_RCCTOIDX_O0O1DU:		/* o_0   \cup o_1^+ */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O1D);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    t[i] = WlzAssignObject(
		   WlzUnion2(o[0], t[WLZ_RCCTOIDX_O1D], &errNum), NULL);
	  }
	  break;
	case WLZ_RCCTOIDX_O0CO1I:               /* o_0^{\circ} \cap o_1   */
	case WLZ_RCCTOIDX_O0O1CI:  /* FALLTHROUGH  o_0   \cap o_1^{\circ} */
	  {
	    int		i0;
	    WlzObject	*c = NULL,
			  *x = NULL;

	    i0 = (i == WLZ_RCCTOIDX_O0CO1I)? 0: 1;
	    c = WlzObjToConvexHull(o[i0], &errNum);
	    if((errNum == WLZ_ERR_NONE) || (errNum == WLZ_ERR_DEGENERATE))
	    {
	      x = WlzAssignObject(
		  WlzConvexHullToObj(c, o[i0]->type, &errNum), NULL);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      t[i] = WlzAssignObject(
		     WlzIntersect2(o[!i0], x, &errNum), NULL);
	    }
	    (void )WlzFreeObj(c);
	    (void )WlzFreeObj(x);
	  }
	  break;
	case WLZ_RCCTOIDX_O0FO1U:		/* o_0^{\bullet} \cup o_1 */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0F);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    t[i] = WlzAssignObject(
		   WlzUnion2(t[WLZ_RCCTOIDX_O0F], o[1], &errNum), NULL);
	  }
	  break;
	case WLZ_RCCTOIDX_O0O1FU:		/* o_0 \cup o_1^{\bullet} */
	  errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O1F);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    t[i] = WlzAssignObject(
		   WlzUnion2(o[0], t[WLZ_RCCTOIDX_O1F], &errNum), NULL);
	  }
	  break;
	default:
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
      }
    }
  }
  else
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	Computes a metric which quantifies the extent to which
* 		the domain of one of the given objects is offset from 
* 		the domain of the second. This is a symetric metric, ie
* 		\f$OST(\Omega_0, \Omega_1) \equiv OST(\Omega_0, \Omega_1)\f$.
*		A domain is computed which is equidistant from the domains
*		of the two given objects, is within maxDist of each object's
*		domain and is within the convex hull of the union of the
*		domains of the two given objects. Within this domain the
*		1st, 2nd and 3rd quantiles of the distance
*		(\f$q_0\f$, \f$q_1\f$ and \f$q_2\f$) are found.
*		The object's domains are classified as offset if
*		\f[
		frac{q_1}{q_1 + (q_1  - q_0) + (q_2 - q_1)} \geq 0.5
		\f]
*		ie
*		\f[
		frac{q_1}{q_1 + q_2 - q_0} \geq 0.5
		\f]
*		Small equi-distant domains with a volume less than half
*		the maximum distance do not classify the relationship as
*		an overlap.
* \param	o			Array with the two given spatial
* 					domain objects, must not be NULL
* 					and nor must the objects.
* \param	t			Array of temporary objects as in
* 					WlzRCCTOIdx.
* \param	maxDist			Maximum distance for offset. This
* 					is used to compute a distance object,
* 					large distances will significantly
* 					increase the processing time.
* \param	dQ0			Destination pointer for 1st quantile
* 					offset distance, must not be NULL.
* \param	dQ1			Destination pointer for 2nd quantile
* 					(ie median) offset distance,
* 					must not be NULL.
* \param	dQ2			Destination pointer for 3rd quantile
* 					offset distance, must not be NULL.
*/
static WlzErrorNum		WlzRCCOffset(
				  WlzObject **o,
				  WlzObject **t,
				  int maxDist,
				  double *dQ0,
				  double *dQ1,
				  double *dQ2)
{
  int		empty = 0;
  double	q[3] = {0.0};
  int		*dHist = NULL;
  WlzObject	*eObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((o[0]->type == WLZ_2D_DOMAINOBJ) &&
      (o[1]->type == WLZ_2D_DOMAINOBJ)) ||
     ((o[0]->type == WLZ_3D_DOMAINOBJ) &&
      (o[1]->type == WLZ_3D_DOMAINOBJ)))
  {
    if((o[0]->domain.core == NULL) || (o[1]->domain.core == NULL))
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
  }
  else
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  /* Compute distance transforms of the two given objects out to a given
   * maximum distance and then using these distances the equi-distant
   * domain between these two objects. The values of the eqi-distant object
   * are those of the distance between the objects.*/
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    WlzObject	*sObj = NULL,
    		*cObj = NULL;
    WlzObject	*dObj[2];

    dObj[0] = dObj[1] = NULL;
    /* Create structuring element with which to dilate the given object
     * domains(by maxDist). */
    sObj = WlzAssignObject(
	   WlzMakeSphereObject(o[0]->type, maxDist, 0, 0, 0,
	                       &errNum), NULL);
    /* Create domain for convex hull of the union of the two given object
     * domains. */
    if(errNum == WLZ_ERR_NONE)
    {
      WlzObject	*uObj = NULL,
      		*xObj = NULL;

      errNum = WlzRCCMakeT(o, t, WLZ_RCCTOIDX_O0O1U);
      if(errNum == WLZ_ERR_NONE)
      {
	uObj = WlzAssignObject(t[WLZ_RCCTOIDX_O0O1U], NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        xObj = WlzAssignObject(
	       WlzObjToConvexHull(uObj, &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        cObj = WlzAssignObject(
	       WlzConvexHullToObj(xObj, o[0]->type, &errNum), NULL);
      }
      (void )WlzFreeObj(xObj);
      (void )WlzFreeObj(uObj);
    }
    /* Dilate the two given objects and find the intersection of the
     * dilated domains with each other and the convex hull computed
     * above. Within his domain compute the distances. */
    if(errNum == WLZ_ERR_NONE)
    {
      for(i = 0; i < 2; ++i)
      {
	WlzObject *tObj = NULL,
		  *rObj = NULL;

	tObj = WlzAssignObject(
	       WlzStructDilation(o[i], sObj, &errNum), NULL);
        if(errNum == WLZ_ERR_NONE)
	{
	  rObj = WlzAssignObject(
	         WlzIntersect2(tObj, cObj, &errNum), NULL);
	}
        (void )WlzFreeObj(tObj);
        if(errNum == WLZ_ERR_NONE)
	{
	  dObj[i] = WlzAssignObject(
		    WlzDistanceTransform(rObj, o[!i],
					 WLZ_OCTAGONAL_DISTANCE,
					 0.0, maxDist, &errNum), NULL);
	}
        (void )WlzFreeObj(rObj);
        if(errNum == WLZ_ERR_NONE)
	{
	  WlzPixelV bgdV;

	  bgdV.type = WLZ_GREY_INT;
	  bgdV.v.inv = maxDist;
	  errNum = WlzSetBackground(dObj[i], bgdV);
	}
        if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
    /* Find the domain which is equi-distant from the two given objects,
     * within the xDist range and within the convex hull of the union of
     * the two given object's domains. */
    (void )WlzFreeObj(sObj); sObj = NULL;
    if(errNum == WLZ_ERR_NONE)
    {
      WlzLong	vol = 0;
      WlzObject *qObj = NULL,
      		*tObj = NULL;

      qObj = WlzAssignObject(
             WlzImageArithmetic(dObj[0], dObj[1], WLZ_BO_EQ, 0, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	WlzPixelV thrV;

	thrV.type = WLZ_GREY_INT;
	thrV.v.inv = 1;
        tObj = WlzAssignObject(
	       WlzThreshold(qObj, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
      }
      /* Check that the eqi-distant domain is of a reasonable size ie has
       * a area or volume greater than half the maximum distance. */
      if(errNum == WLZ_ERR_NONE)
      {
        vol = WlzVolume(tObj, &errNum);
	if((maxDist / 2) >= vol)
	{
	  empty = 1;
	}
      }
      if((errNum == WLZ_ERR_NONE) && !empty)
      {
	WlzObject *hObj = NULL,
		  *mObj = NULL;
	WlzPixelV tmpV;

	tmpV.type = WLZ_GREY_INT;
	tmpV.v.inv = 0;
        mObj = WlzAssignObject(
	       WlzGreyTemplate(dObj[0], tObj, tmpV, &errNum), NULL);
        if(errNum == WLZ_ERR_NONE)
	{
	  tmpV.v.inv = 1;
	  hObj = WlzAssignObject(
	  	 WlzThreshold(mObj, tmpV, WLZ_THRESH_HIGH, &errNum), NULL);
	}
        if(errNum == WLZ_ERR_NONE)
	{
	  tmpV.v.inv = maxDist - 1;
	  eObj = WlzAssignObject(
	  	WlzThreshold(hObj, tmpV, WLZ_THRESH_LOW, &errNum), NULL);
	}
	(void )WlzFreeObj(hObj);
	(void )WlzFreeObj(mObj);
      }
      (void )WlzFreeObj(tObj);
      (void )WlzFreeObj(qObj);
      if((errNum == WLZ_ERR_NONE) && !empty)
      {
	WlzLong vol;

	vol = WlzVolume(eObj, &errNum);
	if((maxDist / 2) >= vol)
	{
	  empty = 1;
	}
      }
    }
    (void )WlzFreeObj(cObj);
    (void )WlzFreeObj(sObj);
    (void )WlzFreeObj(dObj[0]);
    (void )WlzFreeObj(dObj[1]);
  }
  /* Compute a quantised distance histogram in which equi-distant distances
   * are quantized to integer values. */
  if((errNum == WLZ_ERR_NONE) && !empty)
  {
    if((dHist = (int *)AlcCalloc(maxDist + 1, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && !empty)
  {
    if(eObj->type == WLZ_2D_DOMAINOBJ)
    {
      errNum = WlzRCCCompDistHist2D(maxDist, dHist, eObj);
    }
    else
    {
      errNum = WlzRCCCompDistHist3D(maxDist, dHist, eObj);
    }
  }
  WlzFreeObj(eObj);
  if((errNum == WLZ_ERR_NONE) && !empty)
  {
    int		i,
    		j,
		n,
		nq,
		s0 = 0,
		s1 = 0;

    /* Compute the median, first and third quantile offset distances,
     * the ratio of median to the median plus inner inter-quantile range
     * using linear interpolation for the ranges. */
    n = 0;
    for(i = 0; i < maxDist; ++i)
    {
      n += dHist[i];
    }
    i = 0;
    for(j = 1; j <= 3; ++j)
    {
      nq = (n * j) / 4;
      while(s1 < nq)
      {
	s0 = s1;
	s1 += dHist[i++];
      }
      q[j - 1] = i;
      if(s1 > nq)
      {
	q[j - 1] -= (double )(s1 - nq) / (double )(s1 - s0);
      }
    }
  }
  AlcFree(dHist);
  *dQ0 = q[0];
  *dQ1 = q[1];
  *dQ2 = q[2];
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	For each value in the given 2D distance object the count in
* 		the given distance histogram (array) is incremented. The
* 		distance if clipped to the range [0-maxDist].
* \param	maxDist			Maximum distance value.
* \param	dHist			Distance histogram (array).
* \param	dObj			Given 2D distance object.
*/
static WlzErrorNum	WlzRCCCompDistHist2D(
			  int maxDist,
			  int *dHist,
			  WlzObject *dObj)
{
  WlzGreyWSpace	gWSp;
  WlzIntervalWSpace iWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  errNum = WlzInitGreyScan(dObj, &iWSp, &gWSp);
  if(errNum == WLZ_ERR_NONE)
  {
    while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
    {
      int	i;
      int	*p;

      p = gWSp.u_grintptr.inp;
      for(i = iWSp.lftpos; i <= iWSp.rgtpos; ++i)
      {
	int	d;

	d = *p++;
	d = WLZ_CLAMP(d, 0, maxDist);
	++(dHist[d]);
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	For each value in the given 3D distance object the count in
* 		the given distance histogram (array) is incremented. This
* 		is done by forming 2D objects for each plane and calling
* 		WlzRCCCompDistHist2D().
* \param	maxDist			Maximum distance value.
* \param	dHist			Distance histogram (array).
* \param	dObj			Given 3D distance object.
*/
static WlzErrorNum	WlzRCCCompDistHist3D(
			  int maxDist,
			  int *dHist,
			  WlzObject *dObj)
{
  int		n,
  		p;
  WlzVoxelValues *pVal;
  WlzPlaneDomain *pDom;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  pVal = dObj->values.vox;
  pDom = dObj->domain.p;
  n = pDom->lastpl - pDom->plane1 + 1;
  for(p = 0; p < n; ++p)
  {
    if(pDom->domains[p].core != NULL)
    {
      WlzObject	*dObj2;

      dObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, pDom->domains[p], pVal->values[p],
			  NULL, NULL, &errNum);
      errNum = WlzRCCCompDistHist2D(maxDist, dHist, dObj2);
      (void )WlzFreeObj(dObj2);
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  return(errNum);
}
