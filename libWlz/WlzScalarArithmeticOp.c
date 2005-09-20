#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzScalarArithmeticOp.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Functions which apply scalar arithmetic operations
* 		to domain objects.
* \ingroup	WlzArithmetic
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>


/* function:     WlzScalarBinaryOp2    */
/*! 
* \ingroup      WlzArithmetic
* \brief        Apply a binary operation (add subtract etc) to
 each pixel value in the given object. The operand value is in
 <tt>pval</tt>.
*
* \return       Object with transformed grey-values.
* \param    o1	Input object
* \param    pval	Pixel value for binary operation.
* \param    op	Opertor
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

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
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
      if( values.obj = WlzScalarBinaryOp2(o1->values.obj, pval,
					 op, &errNum) ){
	obj = WlzMakeMain(WLZ_TRANS_OBJ, o1->domain, values,
			  NULL, NULL, &errNum);
	break;
      }
      break;
    }
  }

  /* set up temporary object */
  if( (errNum == WLZ_ERR_NONE) && (obj == NULL) ){
    switch( o1->type ){

    case WLZ_2D_DOMAINOBJ:
      values.core = NULL;
      if( (tmp3 = WlzMakeMain(WLZ_2D_DOMAINOBJ, o1->domain,
			      values, NULL, NULL, &errNum)) == NULL ){
	break;
      }
      old_bckgrnd = WlzGetBackground(o1, NULL);
      switch( WlzGreyTableTypeToGreyType(o1->values.core->type, NULL) ){
      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_FLOAT:
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
      }
      values.v = WlzNewValueTb(tmp3,
			       WlzGreyTableType(WLZ_GREY_TAB_RAGR,
						new_grey_type, NULL),
			       new_bckgrnd, &errNum);
      tmp3->values = WlzAssignValues( values, NULL );
      break;

    case WLZ_3D_DOMAINOBJ:
      values.core = NULL;
      if( (tmp3 = WlzMakeMain(WLZ_3D_DOMAINOBJ, o1->domain, values,
			      NULL, NULL, &errNum)) == NULL ){
	break;
      }

      /* now a new destination voxeltable */
      old_bckgrnd = WlzGetBackground(o1, NULL);
      switch( old_bckgrnd.type ){
      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
	new_grey_type = old_bckgrnd.type;
	new_bckgrnd = old_bckgrnd;
	break;	
      case WLZ_GREY_UBYTE:
	new_grey_type = WLZ_GREY_SHORT;
	new_bckgrnd.type = WLZ_GREY_SHORT;
	new_bckgrnd.v.shv = old_bckgrnd.v.ubv;
	break;
      }
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
	  if( tmp2d = WlzMakeMain
	     (WLZ_2D_DOMAINOBJ,
	      tmp3->domain.p->domains[p-tmp3->domain.p->plane1],
	      values2d, NULL, NULL, &errNum) ){
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

      break;

    }
  }

  /* apply operation and free space */
  if( (errNum == WLZ_ERR_NONE) && (obj == NULL) ){
    if( (errNum = WlzScalarBinaryOp(o1, pval, tmp3, op)) != WLZ_ERR_NONE ){
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
