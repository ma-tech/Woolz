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
* \return	RCC8 classification of the given objects.
* \ingroup	WlzBinaryOps
* \brief	The given pair of spatial domain objects are classified
*		using the RCC8.
*		For an explanation of RCC8 classifications
*		see the type definition ::WlzRegConRCC8 and the paper:
*		D.A. Randell, etal,
*		"Discrete Mereotopology for Spatial Reasoning in
*		Automated Histological Image Analysis", PAMI 35(3) 2013.
* 		The classification is performed using simple combinations
* 		of the Woolz union, intersection and difference
* 		morphological operators:
*               \f{eqnarray*}{
                  C_0 &=& O_0   \cap O_1 \\
                  C_1 &=& O_0^+ \cap O_1 \\
                  C_2 &=& (O_0   \cup O_1)   \oplus o_0 \\
                  C_3 &=& (O_0   \cup O_1)   \oplus o_1 \\
                  C_4 &=& (O_0^+ \cup O_1)   \oplus o_0 \\
                  C_5 &=& (O_0   \cup O_1^+) \oplus o_1
  		\f}
*		where
*		  \f$O^+\f$ indicates the dilation of \f$O\f$.
* 		<table width="500" border="0">
		<caption>
		  Basic Morphological Operations for RCC8 Spatial
		  Relationships
		</caption>
                <tr>
                  <td>RCC8</td>
                  <td>\f$C_0\f$</td>
                  <td>\f$C_1\f$</td>
                  <td>\f$C_2\f$</td>
                  <td>\f$C_3\f$</td>
                  <td>\f$C_4\f$</td>
                  <td>\f$C_5\f$</td>
		  <td>Normalised Volume</td>
                </tr>
		<tr>
		  <td>\f$DC(O_0,O_1)\f$</td>
		  <td>0</td>
		  <td>0</td>
		  <td>-</td>
		  <td>-</td>
		  <td>-</td>
		  <td>-</td>
		  <td>0.0</td>
		</tr>
		<tr>
		  <td>\f$EC(O_0,O_1)\f$</td>
		  <td>0</td>
		  <td>1</td>
		  <td>-</td>
		  <td>-</td>
		  <td>-</td>
		  <td>-</td>
		  <td>0.0</td>
		</tr>
		<tr>
		  <td>\f$EQ(O_0,O_1)\f$</td>
		  <td>1</td>
		  <td>-</td>
		  <td>0</td>
		  <td>0</td>
		  <td>-</td>
		  <td>-</td>
		  <td>1.0</td>
		</tr>
		<tr>
		  <td>\f$PO(O_0,O_1)\f$</td>
		  <td>1</td>
		  <td>-</td>
		  <td>1</td>
		  <td>1</td>
		  <td>-</td>
		  <td>-</td>
		  <td>\f$V(O_0 \cap O_1)/V(O_0 \cup O_1)\f$</td>
		</tr>
		<tr>
		  <td>\f$TPP(O_0,O_1)\f$</td>
		  <td>1</td>
		  <td>-</td>
		  <td>1</td>
		  <td>0</td>
		  <td>-</td>
		  <td>1</td>
		  <td>\f$V(O_0)/V(O_0 \cup O_1)\f$</td>
		</tr>
		<tr>
		  <td>\f$NTPP(O_0,O_1)\f$</td>
		  <td>1</td>
		  <td>-</td>
		  <td>1</td>
		  <td>0</td>
		  <td>-</td>
		  <td>0</td>
		  <td>\f$V(O_0)/V(O_0 \cup O_1)\f$</td>
		</tr>
		<tr>
		  <td>\f$TPPI(O_0,O_1)\f$</td>
		  <td>1</td>
		  <td>-</td>
		  <td>0</td>
		  <td>1</td>
		  <td>1</td>
		  <td>-</td>
		  <td>\f$V(O_1)/V(O_0 \cup O_1)\f$</td>
		</tr>
		<tr>
		  <td>\f$NTPPI(O_0,O_1)\f$</td>
		  <td>1</td>
		  <td>-</td>
		  <td>0</td>
		  <td>1</td>
		  <td>0</td>
		  <td>-</td>
		  <td>\f$V(O_1)/V(O_0 \cup O_1)\f$</td>
		</tr>
		</table>
* \param	obj0			First given spatial domain object.
* \param	obj1			Second given spatial domain object.
* \param	dstNrmVol		Destination pointer for the normalized
* 					volume (see above), may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzRegConRCC8 	WlzRegConCalcRCC8(WlzObject *obj0, WlzObject *obj1,
				  double *dstNrmVol, WlzErrorNum *dstErr)
{
  double	vol = 0.0;
  WlzObject	*objI = NULL,
  		*objD = NULL,
  		*objU = NULL,
		*objU0 = NULL,
		*objU1 = NULL;
  WlzRegConRCC8 con = WLZ_REGCON_RCC8_DC;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((obj0 == NULL) || (obj1 == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((obj0->domain.core == NULL) || (obj1->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((obj0->type == WLZ_EMPTY_OBJ) || (obj1->type == WLZ_EMPTY_OBJ))
  {
    con = WLZ_REGCON_RCC8_DC;
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
	    con = WLZ_REGCON_RCC8_EC;
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
		con = WLZ_REGCON_RCC8_EQ;
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
		    con = WLZ_REGCON_RCC8_NTPP;
		  }
		  else
		  {
		    con = WLZ_REGCON_RCC8_TPP;
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
		    con = WLZ_REGCON_RCC8_NTPP;
		  }
		  else
		  {
		    con = WLZ_REGCON_RCC8_TPP;
		  }
		}
		break;
	      case 3:
		con = WLZ_REGCON_RCC8_PO;
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
