#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzScalarArithmeticOp_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzScalarArithmeticOp.c
* \author       Richard Baldock, Bill Hill
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
* \brief	Functions which apply scalar arithmetic operations
* 		to domain objects.
* \ingroup	WlzArithmetic
*/

#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <Wlz.h>

static WlzObject 		*WlzScalarMulAdd2D(
				  WlzObject *iObj,
				  WlzPixelV m,
				  WlzPixelV a,
				  WlzGreyType rGType,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzScalarMulAdd3D(
				  WlzObject *iObj,
				  WlzPixelV m,
				  WlzPixelV a,
				  WlzGreyType rGType,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzScalarMulAddSet2D(
				  WlzObject *rObj,
				  WlzObject *iObj,
				  double m,
				  double a);
static WlzErrorNum 		WlzGreyIncValuesInDomain2D(
				  WlzObject *gObj,
				  WlzObject *dObj);
static WlzErrorNum 		WlzGreyIncValuesInDomain3D(
				  WlzObject *gObj,
				  WlzObject *dObj);

/*!
* \return	Woolz error code.
* \ingroup	WlzArithmetic
* \brief	Increments all valus of the firstobjct which are within
* 		the domain of the second object. The domain of the first
* 		object must cover that of the second.
* \param	gObj		First object.
* \param	dObj		Second object.
*/
WlzErrorNum	WlzGreyIncValuesInDomain(WlzObject *gObj, WlzObject *dObj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gObj == NULL) || (dObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((gObj->domain.core == NULL) || (dObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->type != dObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	errNum = WlzGreyIncValuesInDomain2D(gObj, dObj);
        break;
      case WLZ_3D_DOMAINOBJ:
	errNum = WlzGreyIncValuesInDomain3D(gObj, dObj);
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArithmetic
* \brief	Increments all valus of the firstobjct which are within
* 		the domain of the second object. The domain of the first
* 		object must cover that of the second.
*		Because this is a static object it is assumed that the
*		two 3D objects are known to be valid.
* \param	gObj		First object.
* \param	dObj		Second object.
*/
static WlzErrorNum WlzGreyIncValuesInDomain3D(WlzObject *gObj, WlzObject *dObj)
{
  WlzPlaneDomain *gPD,
  		 *dPD;
  WlzVoxelValues *gVV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gPD = gObj->domain.p;
  gVV = gObj->values.vox;
  dPD = dObj->domain.p;
  if((dPD->plane1 < gPD->plane1) || (dPD->lastpl > gPD->lastpl))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    int		p,
		p0,
		p1;

    p0 = ALG_MAX(dPD->plane1, gPD->plane1);
    p1 = ALG_MIN(dPD->lastpl, gPD->lastpl);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(p = p0; p <= p1; ++p)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	WlzDomain   *doms;
	WlzValues    *vals;
	WlzErrorNum  errNum2D = WLZ_ERR_NONE;

	doms = dPD->domains + p - dPD->plane1;
	vals = gVV->values + p - gPD->plane1;;
	if(((*doms).core != NULL) && ((*vals).core != NULL))
	{
	  WlzObject    *obj2D = NULL;

	  if((obj2D = WlzAssignObject(
		      WlzMakeMain(WLZ_2D_DOMAINOBJ, *doms, *vals, NULL, NULL,
		                  &errNum2D), NULL)) != NULL)
	  {
	    errNum2D = WlzGreyIncValues2D(obj2D);
	    (void )WlzFreeObj(obj2D);
	  }
	}
#ifdef _OPENMP
#pragma omp critical
	{
#endif
	  if((errNum == WLZ_ERR_NONE) && (errNum2D != WLZ_ERR_NONE))
	  {
	    errNum = errNum2D;
	  }
#ifdef _OPENMP
	}
#endif
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArithmetic
* \brief	Increments all valus of the firstobjct which are within
* 		the domain of the second object. The domain of the first
* 		object must cover that of the second.
*		Because this is a static object it is assumed that the
*		two 2D objects are known to be valid.
* \param	gObj		First object.
* \param	dObj		Second object.
*/
static WlzErrorNum WlzGreyIncValuesInDomain2D(WlzObject *gObj, WlzObject *dObj)
{
  WlzObject	*tObj;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tObj = WlzAssignObject(
         WlzMakeMain(WLZ_2D_DOMAINOBJ, dObj->domain, gObj->values, NULL, NULL,
	             &errNum), NULL);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGreyIncValues2D(tObj);
    (void )WlzFreeObj(tObj);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArithmetic
* \brief	Increments all values within the given object.
* \param	obj		Given object.
*/
WlzErrorNum 	WlzGreyIncValues2D(WlzObject *obj)
{
  WlzGreyWSpace gWSp;
  WlzIntervalWSpace iWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  errNum = WlzInitGreyScan(obj, &iWSp, &gWSp);
  while((errNum == WLZ_ERR_NONE) &&
        ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
  {
    int	     i,
    	     len;
    WlzGreyP gP;
    
    gP = gWSp.u_grintptr;
    len = iWSp.rgtpos - iWSp.lftpos + 1;
    switch(gWSp.pixeltype)
    {
      case WLZ_GREY_INT:
        for(i = 0; i < len; ++i)
	{
	  *(gP.inp)++ += 1;
	}
	break;
      case WLZ_GREY_SHORT:
        for(i = 0; i < len; ++i)
	{
	  *(gP.shp)++ += 1;
	}
	break;
      case WLZ_GREY_UBYTE:
        for(i = 0; i < len; ++i)
	{
	  *(gP.ubp)++ += 1;
	}
	break;
      case WLZ_GREY_FLOAT:
        for(i = 0; i < len; ++i)
	{
	  *(gP.flp)++ += 1.0f;
	}
	break;
      case WLZ_GREY_DOUBLE:
        for(i = 0; i < len; ++i)
	{
	  *(gP.dbp)++ += 1.0;
	}
	break;
      default:
        break;
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}

/*! 
* \return       Object with transformed grey-values.
* \ingroup      WlzArithmetic
* \brief        Apply a binary operation (add subtract etc) to
*               each pixel value in the given object. The operand value
*               is in <tt>pval</tt>.
* \param    o1	Input object
* \param    pval	Pixel value for binary operation.
* \param    op		Opertor
* \param    dstErr	Error return.
* \par      Source:
*                WlzScalarArithmeticOp.c
*/
WlzObject *WlzScalarBinaryOp2(
  WlzObject	*o1,
  WlzPixelV	pval,
  WlzBinaryOperatorType	op,
  WlzErrorNum	*dstErr)
{
  WlzObject	*obj=NULL, *tmp3;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzValues	values;
  WlzPixelV	old_bckgrnd, new_bckgrnd;
  WlzGreyType	new_grey_type;
  int		p;

  /* check object pointers */
  if( (o1 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* check object types - WLZ_EMPTY_OBJ is legal */
  if( errNum == WLZ_ERR_NONE ){
    switch( o1->type ){

    case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
    case WLZ_3D_DOMAINOBJ: /* FALLTHROUGH */
    case WLZ_TRANS_OBJ:
      break;

    case WLZ_EMPTY_OBJ:
      obj = WlzMakeEmpty(&errNum);
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    }
  }

  /* check domains and valuetables */
  if( (errNum == WLZ_ERR_NONE) && (obj == NULL) ){
    switch( o1->type ){

    case WLZ_2D_DOMAINOBJ:
      if( (o1->domain.core == NULL) ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( (o1->values.core == NULL) ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( (o1->domain.core == NULL) ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( (o1->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN) ){
	errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
	break;
      }
      if( (o1->values.core == NULL) ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      if( (o1->values.vox->type != WLZ_VOXELVALUETABLE_GREY) ){
	errNum = WLZ_ERR_VOXELVALUES_TYPE;
	break;
      }
      break;

    case WLZ_TRANS_OBJ:
      if( (o1->domain.core == NULL) ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( (o1->values.core == NULL) ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      if((values.obj = WlzScalarBinaryOp2(o1->values.obj, pval,
					 op, &errNum)) != NULL){
	obj = WlzMakeMain(WLZ_TRANS_OBJ, o1->domain, values,
			  NULL, NULL, &errNum);
	break;
      }
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* set up temporary object */
  if((errNum == WLZ_ERR_NONE) && (obj == NULL)){
    switch( o1->type ){

    case WLZ_2D_DOMAINOBJ:
      values.core = NULL;
      if( (tmp3 = WlzMakeMain(WLZ_2D_DOMAINOBJ, o1->domain,
			      values, NULL, NULL, &errNum)) == NULL ){
	break;
      }
      old_bckgrnd = WlzGetBackground(o1, NULL);
      switch( WlzGreyTableTypeToGreyType(o1->values.core->type, NULL) ){
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	new_grey_type = WlzGreyTableTypeToGreyType(o1->values.core->type,
						   NULL);
	new_bckgrnd = old_bckgrnd;
	break;	
      case WLZ_GREY_UBYTE:
	new_grey_type = WLZ_GREY_SHORT;
	new_bckgrnd.type = WLZ_GREY_SHORT;
	new_bckgrnd.v.shv = old_bckgrnd.v.ubv;
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
      if(errNum == WLZ_ERR_NONE){
	values.v = WlzNewValueTb(tmp3,
				 WlzGreyTableType(WLZ_GREY_TAB_RAGR,
						  new_grey_type, NULL),
				 new_bckgrnd, &errNum);
	tmp3->values = WlzAssignValues( values, NULL );
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      values.core = NULL;
      if((tmp3 = WlzMakeMain(WLZ_3D_DOMAINOBJ, o1->domain, values,
			     NULL, NULL, &errNum)) == NULL){
	break;
      }

      /* now a new destination voxeltable */
      old_bckgrnd = WlzGetBackground(o1, NULL);
      switch( old_bckgrnd.type ){
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	new_grey_type = old_bckgrnd.type;
	new_bckgrnd = old_bckgrnd;
	break;	
      case WLZ_GREY_UBYTE:
	new_grey_type = WLZ_GREY_SHORT;
	new_bckgrnd.type = WLZ_GREY_SHORT;
	new_bckgrnd.v.shv = old_bckgrnd.v.ubv;
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
      }
      if(errNum == WLZ_ERR_NONE){
	values.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					 tmp3->domain.p->plane1,
					 tmp3->domain.p->lastpl,
					 new_bckgrnd, NULL, &errNum);
	if( values.vox == NULL ){break;}

	for(p=tmp3->domain.p->plane1; p <= tmp3->domain.p->lastpl; p++){
	  WlzValues	values2d;

	  /* currently test for NULL domain to imply WLZ_EMPTY_DOMAIN */
	  if( tmp3->domain.p->domains[p-tmp3->domain.p->plane1].core
	     == NULL ){
	    values2d.core = NULL;
	  }
	  else {
	    WlzObject	*tmp2d;
	    values2d.core = NULL;
	    if((tmp2d = WlzMakeMain(WLZ_2D_DOMAINOBJ,
		             tmp3->domain.p->domains[p-tmp3->domain.p->plane1],
		             values2d, NULL, NULL, &errNum)) != NULL){
	      values2d.v = WlzNewValueTb(
		tmp2d, WlzGreyTableType(WLZ_GREY_TAB_RAGR, new_grey_type,
					NULL), new_bckgrnd, &errNum);
	      WlzFreeObj( tmp2d );
	    }
	  }
	  values.vox->values[p-tmp3->domain.p->plane1] = 
	    WlzAssignValues(values2d, NULL);
	}
	tmp3->values = WlzAssignValues(values, NULL);
      }

      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* apply operation and free space */
  if((errNum == WLZ_ERR_NONE) && (obj == NULL)){
    if((errNum = WlzScalarBinaryOp(o1, pval, tmp3, op)) != WLZ_ERR_NONE){
      WlzFreeObj( tmp3 );
    }
    else {
      obj = tmp3;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return obj;
}


/* function:     WlzScalarAdd    */
/*! 
* \ingroup      WlzArithmetic
* \brief        Add a value to each pixel of an object.
*
* \return       Object with value added to each pixel.
* \param    o1	Input object
* \param    pval	Value to be added.
* \param    dstErr	Error return.
* \par      Source:
*                WlzScalarArithmeticOp.c
*/
WlzObject *WlzScalarAdd(
  WlzObject	*o1,
  WlzPixelV	pval,
  WlzErrorNum	*dstErr)
{
  return WlzScalarBinaryOp2(o1, pval, WLZ_BO_ADD, dstErr);
}

/* function:     WlzScalarSubtract    */
/*! 
* \ingroup      WlzArithmetic
* \brief        Subtract a value from each pixel of an object.
*
* \return       Object with value subtacted from each pixel.
* \param    o1	Input object
* \param    pval	Value to be subtracted.
* \param    dstErr	Error return.
* \par      Source:
*                WlzScalarArithmeticOp.c
*/
WlzObject *WlzScalarSubtract(
  WlzObject	*o1,
  WlzPixelV	pval,
  WlzErrorNum	*dstErr)
{
  return WlzScalarBinaryOp2(o1, pval, WLZ_BO_SUBTRACT, dstErr);
}

/* function:     WlzScalarMultiply    */
/*! 
* \ingroup      WlzArithmetic
* \brief        Multiply each pixel of an object.
*
* \return       Object with multiplied pixel values.
* \param    o1	Input object
* \param    pval	Multiplication factor.
* \param    dstErr	Error return.
* \par      Source:
*                WlzScalarArithmeticOp.c
*/
WlzObject *WlzScalarMultiply(
  WlzObject	*o1,
  WlzPixelV	pval,
  WlzErrorNum	*dstErr)
{
  return WlzScalarBinaryOp2(o1, pval, WLZ_BO_MULTIPLY, dstErr);
}

/* function:     WlzScalarDivide    */
/*! 
* \ingroup      WlzArithmetic
* \brief        Divide each pixel of an object.
*
* \return       Object with each pixel divided.
* \param    o1	Input object
* \param    pval	Division value.
* \param    dstErr	Error return.
* \par      Source:
*                WlzScalarArithmeticOp.c
*/
WlzObject *WlzScalarDivide(
  WlzObject	*o1,
  WlzPixelV	pval,
  WlzErrorNum	*dstErr)
{
  return WlzScalarBinaryOp2(o1, pval, WLZ_BO_DIVIDE, dstErr);
}

/*!
* \return	New woolz object or NULL on error.
* \ingroup	WlzArithmetic
* \brief	Scales the values of the given Woolz object so that
*               \f$v_{new} = m v_{given} + a.\f$
* \param	iObj			Given object.
* \param	m			Value to multiply object values by.
* \param	a			Value to add to product.
* \param	rGType			Required grey type for returned object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzScalarMulAdd(WlzObject *iObj, WlzPixelV m, WlzPixelV a,
				WlzGreyType rGType, WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(iObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(iObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(iObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        rObj = WlzScalarMulAdd2D(iObj, m, a, rGType, &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
        rObj = WlzScalarMulAdd3D(iObj, m, a, rGType, &errNum);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New woolz object or NULL on error.
* \ingroup	WlzArithmetic
* \brief	Scales the values of the given 2D Woolz object so that
* 		\f$v_{new} = m v_{given} + a.\f$ The input object is known
* 		to be a valid 2D domain object with grey values.
* \param	iObj			Given object.
* \param	m			Value to multiply object values by.
* \param	a			Value to add to product.
* \param	rGType			Required grey type for returned object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzScalarMulAdd2D(WlzObject *iObj, WlzPixelV m, WlzPixelV a,
				WlzGreyType rGType, WlzErrorNum *dstErr)
{
  WlzValues	rValues;
  WlzObjectType rVType;
  WlzPixelV	bgdV;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rValues.core = NULL;
  bgdV = WlzGetBackground(iObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    rVType = WlzGreyTableTypeToTableType(iObj->values.v->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rVType = WlzGreyTableType(rVType, rGType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rValues.v = WlzNewValueTb(iObj, rVType, bgdV, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, iObj->domain, rValues,
    		       iObj->plist, iObj->assoc, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(rGType)
    {
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_RGBA:  /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	WlzValueConvertPixel(&m, m, WLZ_GREY_DOUBLE);
	WlzValueConvertPixel(&a, a, WLZ_GREY_DOUBLE);
	errNum = WlzScalarMulAddSet2D(rObj, iObj, m.v.dbv, a.v.dbv);
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj == NULL)
    {
      (void )WlzFreeValueTb(rValues.v);
    }
    else
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New woolz object or NULL on error.
* \ingroup	WlzArithmetic
* \brief	Scales the values of the given 3D Woolz object so that
* 		\f$v_{new} = m v_{given} + a.\f$ The input object is known
* 		to be a valid 3D domain object with grey values.
* \param	iObj			Given object.
* \param	m			Value to multiply object values by.
* \param	a			Value to add to product.
* \param	rGType			Required grey type for returned object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzScalarMulAdd3D(WlzObject *iObj, WlzPixelV m, WlzPixelV a,
				WlzGreyType rGType, WlzErrorNum *dstErr)
{
  WlzValues	rValues;
  WlzObjectType rVType;
  WlzPixelV	bgdV;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rValues.core = NULL;
  bgdV = WlzGetBackground(iObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    rVType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, rGType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rValues.vox = WlzNewValuesVox(iObj, rVType, bgdV, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idP,
    		plMin,
    		plMax;
    WlzPlaneDomain *iPDom;
    WlzVoxelValues *iVox,
    		*rVox;

    iPDom = iObj->domain.p;
    iVox = iObj->values.vox;
    rVox = rValues.vox;
    plMin = iPDom->plane1;
    plMax = iPDom->lastpl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idP = plMin; idP <= plMax; ++idP)
    {
      WlzErrorNum errNum2D = WLZ_ERR_NONE;

      if(errNum == WLZ_ERR_NONE)
      {
	int	      idO;
	WlzDomain *iDom2D;
	WlzValues *iVal2D,
		  *rVal2D;

	idO = idP - iPDom->plane1;
	iDom2D = iPDom->domains + idO;
	iVal2D = iVox->values + idO;
	rVal2D = rVox->values + idO;
	if(((*iDom2D).core != NULL) && ((*iVal2D).core != NULL) &&
	   ((*rVal2D).core != NULL))
	{
	  WlzObject *iObj2D = NULL,
		    *rObj2D = NULL;
	  
	  if((iObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *iDom2D, *iVal2D,
				   NULL, NULL, &errNum2D)) != NULL)
	  {
	    rObj2D = WlzScalarMulAdd2D(iObj2D, m, a, rGType, &errNum2D);
	  }
	  if(errNum2D == WLZ_ERR_NONE)
	  {
	    *rVal2D = WlzAssignValues(rObj2D->values, NULL);
	  }
	  (void )WlzFreeObj(iObj2D);
	  (void )WlzFreeObj(rObj2D);
	}
      }
#ifdef _OPENMP
#pragma omp critical
      {
#endif
	if(errNum2D != WLZ_ERR_NONE)
	{
	  errNum = errNum2D;
	}
#ifdef _OPENMP
      }
#endif
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, iObj->domain, rValues,
    		       iObj->plist, iObj->assoc, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj == NULL)
    {
      (void )WlzFreeVoxelValueTb(rValues.vox);
    }
    else
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
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
* \ingroup	WlzArithmetic
* \brief	Sets the values of the return object from the input object
* 		using simple linear scaling, see WlzScalarMulAdd(). The
* 		objects are known to be 2D, have the same domain.
* \param	rObj
* \param	iObj
* \param	m
* \param	a
*/
static WlzErrorNum WlzScalarMulAddSet2D(WlzObject *rObj, WlzObject *iObj,
				     double m, double a)
{
  int		bufLen;
  WlzGreyWSpace iGWSp,
  		rGWSp;
  WlzIntervalWSpace iIWSp,
  		    rIWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bufLen = iObj->domain.i->lastkl - iObj->domain.i->kol1 + 1;
  if((bufLen != iObj->domain.i->lastkl - iObj->domain.i->kol1 + 1) ||
     (bufLen < 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if(bufLen > 0)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(iObj, &iIWSp, &iGWSp);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(rObj, &rIWSp, &rGWSp);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      double	*buf = NULL;

      if((buf = AlcMalloc(sizeof(double) * bufLen)) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	while((errNum = WlzNextGreyInterval(&iIWSp)) == WLZ_ERR_NONE)
	{
	  int	t,
	  	idN,
	  	itvLen;
	  double f;


	  itvLen = iIWSp.colrmn;
	  (void )WlzNextGreyInterval(&rIWSp);
	  switch(iGWSp.pixeltype)
	  {
	    case WLZ_GREY_INT:
	      WlzValueCopyIntToDouble(buf, iGWSp.u_grintptr.inp, itvLen);
	      break;
	    case WLZ_GREY_SHORT:
	      WlzValueCopyShortToDouble(buf, iGWSp.u_grintptr.shp, itvLen);
	      break;
	    case WLZ_GREY_UBYTE:
	      WlzValueCopyUByteToDouble(buf, iGWSp.u_grintptr.ubp, itvLen);
	      break;
	    case WLZ_GREY_FLOAT:
	      WlzValueCopyFloatToDouble(buf, iGWSp.u_grintptr.flp, itvLen);
	      break;
	    case WLZ_GREY_DOUBLE:
	      WlzValueCopyDoubleToDouble(buf, iGWSp.u_grintptr.dbp, itvLen);
	      break;
	    case WLZ_GREY_RGBA:
	      WlzValueCopyRGBAToDouble(buf, iGWSp.u_grintptr.rgbp, itvLen);
	      break;
	    default:
	      break;
	  }
	  switch(rGWSp.pixeltype)
	  {
	    case WLZ_GREY_UBYTE:
	      for(idN = 0; idN < itvLen; ++idN)
	      {
		f = (buf[idN] * m) + a;
		f = WLZ_CLAMP(f, 0, 255);
		rGWSp.u_grintptr.ubp[idN] = WLZ_NINT(f);
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      for(idN = 0; idN < itvLen; ++idN)
	      {
		f = (buf[idN] * m) + a;
		f = WLZ_CLAMP(f, SHRT_MIN, SHRT_MAX);
		rGWSp.u_grintptr.shp[idN] = WLZ_NINT(f);
	      }
	      break;
	    case WLZ_GREY_INT:
	      for(idN = 0; idN < itvLen; ++idN)
	      {
		f = (buf[idN] * m) + a;
		f = WLZ_CLAMP(f, INT_MIN, INT_MAX);
		rGWSp.u_grintptr.inp[idN] = WLZ_NINT(f);
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      for(idN = 0; idN < itvLen; ++idN)
	      {
		WlzUInt	u;

		f = (buf[idN] * m) + a;
		f = WLZ_CLAMP(f, 0, 255);
		t = WLZ_NINT(f);
                WLZ_RGBA_RGBA_SET(u, t, t, t, 255);
		rGWSp.u_grintptr.inp[idN] = u;
	      }
	    case WLZ_GREY_FLOAT:
	      for(idN = 0; idN < itvLen; ++idN)
	      {
		double	t;

		t = (buf[idN] * m) + a;
		rGWSp.u_grintptr.flp[idN] = WLZ_CLAMP(t, -(FLT_MAX), FLT_MAX);
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      for(idN = 0; idN < itvLen; ++idN)
	      {
		rGWSp.u_grintptr.dbp[idN] = (buf[idN] * m) + a;
	      }
	      break;
	    default:
	      break;
	  }
	}
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
      AlcFree(buf);
    }
  }
  return(errNum);
}
