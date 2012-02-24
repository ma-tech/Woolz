#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCompoundArrayToScalar_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCompoundArrayToScalar.c
* \author       Bill Hill
* \date         November 2009
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
* \brief	Functions to convert vector values held in a compound
*               array to scalar values in a domain object.
* \ingroup	WlzValueUtils
*/
#include <Wlz.h>


static WlzObject 		*WlzCompoundArrayToScalar2D(
				  WlzCompoundArray *cpd,
				  WlzObject *iObj,
				  WlzBinaryOperatorType bop,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCompoundArrayToScalar3D(
				  WlzCompoundArray *cpd,
				  WlzObject *iObj,
				  WlzBinaryOperatorType bop,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzCompoundArrayToScalarAdd1(
				  WlzObject *dst,
				  WlzObject *src,
				  WlzBinaryOperatorType bop,
				  int srcObjIdx);
static WlzErrorNum 		WlzCompoundArrayToScalarAdd2(
				  WlzObject *obj,
				  WlzBinaryOperatorType bop);
/*!
* \return	New object with double scalar values but the same domain
* 		as the given object.
* \ingroup	WlzValueUtils
* \brief	Computes an object with scalar values from the given
* 		compound array object using the given binary operator.
* \param	cpd			Given compound array in which all
* 					objects must be of the same type
* 					and have the same domain.
* \param	bop			Binary operator. Currently only
* 					WLZ_BO_MODULUS is implimented.
* \param	dstErr			Destiantion error pointer, may be NULL.
*/
WlzObject	*WlzCompoundArrayToScalar(WlzCompoundArray *cpd,
					  WlzBinaryOperatorType bop,
				          WlzErrorNum *dstErr)
{
  int		idN;
  WlzObject	*iObj = NULL,
  		*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cpd == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cpd->type != WLZ_COMPOUND_ARR_1)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(cpd->n < 1)
  {
    errNum = WLZ_ERR_OBJECT_DATA;
  }
  else
  {
    for(idN = 0; idN < cpd->n; ++idN)
    {
      if(cpd->o[idN] == NULL)
      {
        errNum = WLZ_ERR_OBJECT_NULL;
      }
      else if(cpd->o[idN]->domain.core == NULL)
      {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if(cpd->o[idN]->values.core == NULL)
      {
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if(idN > 0)
      {
        if(cpd->o[idN]->type != cpd->o[idN - 1]->type)
	{
	  errNum = WLZ_ERR_OBJECT_TYPE;
	}
	else if(cpd->o[idN]->domain.core->type !=
	        cpd->o[idN - 1]->domain.core->type)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute intersection of the compound object domains. */
    iObj = WlzIntersectN(cpd->n, cpd->o, 0, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(cpd->o[0]->type)
    {
      case WLZ_2D_DOMAINOBJ:
	rObj = WlzCompoundArrayToScalar2D(cpd, iObj, bop, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	rObj = WlzCompoundArrayToScalar3D(cpd, iObj, bop, &errNum);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  (void )WlzFreeObj(iObj);
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New object with double scalar values but the same domain
* 		as the given object. This function is for 2D domain objects.
* \ingroup	WlzValueUtils
* \brief	Computes an object with scalar values from the given
* 		compound array object using the given binary operator.
* \param	cpd			Given compound array in which all
* 					objects must be of the same type
* 					and have the same domain.
* \param	iObj			Object, the domain of which is the
* 					intersection of the compound object's
* 					object domains.
* \param	bop			Binary operator. Currently only
* 					WLZ_BO_MODULUS is implimented.
* \param	dstErr			Destiantion error pointer, may be NULL.
* */
static WlzObject *WlzCompoundArrayToScalar2D(WlzCompoundArray *cpd,
					     WlzObject *iObj,
					     WlzBinaryOperatorType bop,
				             WlzErrorNum *dstErr)
{
  int		idN;
  WlzObject	*tObj,
  		*rObj = NULL;
  WlzObjectType	gTT;
  WlzPixelV	bgd;
  WlzValues	values;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  bgd.type = WLZ_GREY_DOUBLE;
  bgd.v.dbv = 0.0;
  values.core = NULL;
  gTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    values.v = WlzNewValueTb(iObj, gTT, bgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, iObj->domain, values,
                      NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < cpd->n; ++idN)
    {
      tObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, iObj->domain, cpd->o[idN]->values,
      			 NULL, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzCompoundArrayToScalarAdd1(rObj, tObj, bop, idN);
      }
      (void )WlzFreeObj(tObj);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCompoundArrayToScalarAdd2(rObj, bop);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj != NULL)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else if(values.core != NULL)
    {
      (void )WlzFreeValues(values);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New object with double scalar values but the same domain
* 		as the given object. This function is for 3D domain objects.
* \ingroup	WlzValueUtils
* \brief	Computes an object with scalar values from the given
* 		compound array object using the given binary operator.
* \param	cpd			Given compound array in which all
* 					objects must be of the same type
* 					and have the same domain.
* \param	iObj			Object, the domain of which is the
* 					intersection of the compound object's
* 					object domains.
* \param	bop			Binary operator. Currently only
* 					WLZ_BO_MODULUS is implimented.
* \param	dstErr			Destiantion error pointer, may be NULL.
* */
static WlzObject *WlzCompoundArrayToScalar3D(WlzCompoundArray *cpd,
					     WlzObject *iObj,
					     WlzBinaryOperatorType bop,
				             WlzErrorNum *dstErr)
{
  int		idI,
  		idN,
		idP,
		idR;
  WlzObject	*in2D,
  		*out2D,
		*rObj = NULL;
  WlzObjectType	gTT;
  WlzPixelV	bgd;
  WlzValues	values;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  bgd.type = WLZ_GREY_DOUBLE;
  bgd.v.dbv = 0.0;
  values.core = NULL;
  gTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    values.vox = WlzNewValuesVox(iObj, gTT, bgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, iObj->domain, values,
                      NULL, NULL, &errNum);
  }
  for(idN = 0; (errNum == WLZ_ERR_NONE) && (idN < cpd->n); ++idN)
  {
    for(idP = rObj->domain.p->plane1; idP <= rObj->domain.p->lastpl; ++idP)
    {
      in2D = out2D = NULL;
      idR = idP - rObj->domain.p->plane1;
      idI = idP - cpd->o[idN]->domain.p->plane1;
      in2D = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			 *(iObj->domain.p->domains + idR),
			 *(cpd->o[idN]->values.vox->values + idI),
			 NULL, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	out2D = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			    *(rObj->domain.p->domains + idR),
			    *(rObj->values.vox->values + idR),
			    NULL, NULL, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzCompoundArrayToScalarAdd1(out2D, in2D, bop, idN);
      }
      (void )WlzFreeObj(in2D);
      (void )WlzFreeObj(out2D);
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idP = rObj->domain.p->plane1; idP <= rObj->domain.p->lastpl; ++idP)
    {
      out2D = NULL;
      idR = idP - rObj->domain.p->plane1;
      out2D = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			  *(rObj->domain.p->domains + idR),
			  *(rObj->values.vox->values + idR),
			  NULL, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzCompoundArrayToScalarAdd2(out2D, bop);
      }
      (void )WlzFreeObj(out2D);
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValueUtils
* \brief	Adds values of a 2D source object to the 2D destination
* 		object according to the manner of the binary operator.
* \param	dst			Destination object.
* \param	src			Source object.
* \param	bop			The binary operator.
* \param	srcObjIdx		Source object index [0-...].
*/
static WlzErrorNum WlzCompoundArrayToScalarAdd1(WlzObject *dst,
				WlzObject *src, WlzBinaryOperatorType bop,
				int srcObjIdx)
{
  int		idK,
  		width;
  double	*bP,
  		*dP;
  WlzGreyP	bufP;
  WlzIntervalWSpace dIWSp,
  		sIWSp;
  WlzGreyWSpace dGWSp,
  		sGWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  width = dst->domain.i->lastkl - dst->domain.i->kol1 + 1;
  if((bufP.v = AlcMalloc(sizeof(double) * width)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(dst, &dIWSp, &dGWSp);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(src, &sIWSp, &sGWSp);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      while(((errNum = WlzNextGreyInterval(&dIWSp)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&sIWSp)) == WLZ_ERR_NONE))
      {
	width = sIWSp.rgtpos - sIWSp.lftpos + 1;
	WlzValueCopyGreyToGrey(bufP, 0, WLZ_GREY_DOUBLE,
			       sGWSp.u_grintptr, 0, sGWSp.pixeltype,
			       width);
	bP = bufP.dbp;
	dP = dGWSp.u_grintptr.dbp;
	if(srcObjIdx == 0)
	{
	  switch(bop)
	  {
	    case WLZ_BO_MODULUS:
	      for(idK = 0; idK < width; ++idK)
	      {
		*dP++ = *bP * *bP;
		++bP;
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_PARAM_TYPE;
	      break;
	  }
	}
	else
	{
	  switch(bop)
	  {
	    case WLZ_BO_MODULUS:
	      for(idK = 0; idK < width; ++idK)
	      {
		*dP++ += *bP * *bP;
		++bP;
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_PARAM_TYPE;
	      break;
	  }
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }
  AlcFree(bufP.v);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValueUtils
* \brief	Having already added the values to the 2D destination
* 		object this function performs a second pass though it's
* 		values as required by the binary operator.
* \param	dst			Destination object.
* \param	bop			The binary operator.
*/
static WlzErrorNum WlzCompoundArrayToScalarAdd2(WlzObject *obj,
				WlzBinaryOperatorType bop)
{
  int		idK,
  		iWidth;
  double	*dP;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace gWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if(bop == WLZ_BO_MODULUS)
  {
    errNum = WlzInitGreyScan(obj, &iWSp, &gWSp);
    while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
    {
      dP = gWSp.u_grintptr.dbp;
      iWidth = iWSp.rgtpos - iWSp.lftpos + 1;
      for(idK = 0; idK < iWidth; ++idK)
      {
	*dP = sqrt(*dP);
	++dP;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  return(errNum);
}
