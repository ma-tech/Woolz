#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzSepTrans_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzSepTrans.c
* \author       Richard Baldock
* \date         May 2003
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
* \brief	Execute a separable transform on a 2D domain object.
* \ingroup	WlzValuesFilters
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/* function:     WlzSepTrans    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        Perform 2D seperable transform on a 2D grey-level image.
* 		It is the users responsibility to ensure that the grey-value
* 		types are appropriate.
* 		Now extended to include compound objects.
*
* \return       Pointer to transformed object
* \param    obj	Input object pointer
* \param    x_fun	Convolution function to be applied in the x-direction (along the rows).
* \param    x_params	Parameter pointer to be passed to the x-function.
* \param    y_fun	Convolution function to be applied in the y-direction (down the columns).
* \param    y_params	Parameter pointer to be passed to the y-function.
* \param    dstErr	error return.
* \par      Source:
*                WlzSepTrans.c
*/
WlzObject *WlzSepTrans(
  WlzObject		*obj,
  WlzIntervalConvFunc	x_fun,
  void			*x_params,
  WlzIntervalConvFunc	y_fun,
  void			*y_params,
  WlzErrorNum		*dstErr)
{
  WlzIntervalWSpace	iwspace;
  WlzGreyWSpace		gwspace;
  WlzSepTransWSpace	stwspc;
  WlzValues		values;
  WlzObject		*obj1, *obj2;
  int			i, width, height;
  WlzCompoundArray	*cobj1, *cobj2;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object pointers and type */
  obj2 = NULL;
  stwspc.outbuf.p.dbp = NULL;
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      /* check domain and valuetable */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      break;

    default:
    case WLZ_3D_DOMAINOBJ:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_TRANS_OBJ:
      if((obj1 = WlzSepTrans(obj->values.obj,
			     x_fun, x_params, y_fun, y_params,
			     &errNum)) != NULL){
	values.obj = obj1;
	return WlzMakeMain(obj->type, obj->domain, values,
			   NULL, obj, dstErr);
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      cobj1 = (WlzCompoundArray *) obj;
      if((cobj2 = WlzMakeCompoundArray(cobj1->type, 1, cobj1->n, NULL,
				       cobj1->otype, &errNum)) != NULL){
	/* transform each object, ignore type errors */
	for(i=0; i < cobj1->n; i++){
	  cobj2->o[i] =
	    WlzAssignObject(WlzSepTrans(cobj1->o[i],
					x_fun, x_params, y_fun, y_params,
					&errNum), NULL);
	}
	return (WlzObject *) cobj2;
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);
    }
  }

  /* make space for a calculation buffer - assume worst case of doubles */
  if( errNum == WLZ_ERR_NONE ){
    width  = obj->domain.i->lastkl - obj->domain.i->kol1 + 1;
    height = obj->domain.i->lastln - obj->domain.i->line1 + 1;
    stwspc.inbuf.type = WlzGreyTableTypeToGreyType(obj->values.core->type,
						   NULL);
    stwspc.bckgrnd = WlzGetBackground(obj, NULL);

    switch( stwspc.inbuf.type ){
    case WLZ_GREY_INT:
    case WLZ_GREY_SHORT:
    case WLZ_GREY_UBYTE:
      stwspc.outbuf.type = WLZ_GREY_INT;
      break;
    default:
      stwspc.outbuf.type = stwspc.inbuf.type;
      break;
    }

    if( (stwspc.outbuf.p.dbp = (double *) 
	 AlcMalloc(sizeof(double)*(width>height?width:height))) == NULL){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  /* tranpose the object - interchange x & y coordinates */
  if( errNum == WLZ_ERR_NONE ){
    obj1 = WlzTransposeObj(obj, &errNum);
  }

  /* perform the y convolution */
  if((errNum == WLZ_ERR_NONE) &&
     ((errNum = WlzInitGreyScan(obj1, &iwspace, &gwspace)) == WLZ_ERR_NONE))
  {
    while( (errNum = WlzNextGreyInterval(&iwspace)) == WLZ_ERR_NONE ){
      stwspc.inbuf.p.inp = gwspace.u_grintptr.inp;
      stwspc.len = iwspace.rgtpos - iwspace.lftpos + 1;
      if((errNum = (*y_fun)(&stwspc, y_params)) == WLZ_ERR_NONE){
	switch(  stwspc.inbuf.type ){
	case WLZ_GREY_INT:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.inp++)
	    *stwspc.inbuf.p.inp = stwspc.outbuf.p.inp[i];
	  break;
	case WLZ_GREY_SHORT:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.shp++)
	    *stwspc.inbuf.p.shp = stwspc.outbuf.p.inp[i];
	  break;
	case WLZ_GREY_UBYTE:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.ubp++)
	    *stwspc.inbuf.p.ubp = stwspc.outbuf.p.inp[i];
	  break;
	case WLZ_GREY_FLOAT:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.flp++)
	    *stwspc.inbuf.p.flp = stwspc.outbuf.p.flp[i];
	  break;
	case WLZ_GREY_DOUBLE:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.dbp++)
	    *stwspc.inbuf.p.dbp = stwspc.outbuf.p.dbp[i];
	  break;
	case WLZ_GREY_RGBA:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.rgbp++)
	    *stwspc.inbuf.p.rgbp = stwspc.outbuf.p.rgbp[i];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
	}
      }
      else {
	break; /* break from the while loop */
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  /* rotate back and free temporary object */
  if( errNum == WLZ_ERR_NONE ){
    obj2 = WlzTransposeObj(obj1, &errNum);
    WlzFreeObj(obj1);
  }

  /* perform x convolution */
  if((errNum == WLZ_ERR_NONE) &&
     ((errNum = WlzInitGreyScan(obj2, &iwspace, &gwspace)) == WLZ_ERR_NONE))
  {
    while((errNum = WlzNextGreyInterval(&iwspace)) == WLZ_ERR_NONE){
      stwspc.inbuf.p.inp = gwspace.u_grintptr.inp;
      stwspc.len = iwspace.rgtpos - iwspace.lftpos + 1;
      if( (errNum = (*x_fun)(&stwspc, x_params)) == WLZ_ERR_NONE ){
	switch(  gwspace.pixeltype ){
	case WLZ_GREY_INT:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.inp++)
	    *stwspc.inbuf.p.inp = stwspc.outbuf.p.inp[i];
	  break;
	case WLZ_GREY_SHORT:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.shp++)
	    *stwspc.inbuf.p.shp = stwspc.outbuf.p.inp[i];
	  break;
	case WLZ_GREY_UBYTE:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.ubp++)
	    *stwspc.inbuf.p.ubp = stwspc.outbuf.p.inp[i];
	  break;
	case WLZ_GREY_FLOAT:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.flp++)
	    *stwspc.inbuf.p.flp = stwspc.outbuf.p.flp[i];
	  break;
	case WLZ_GREY_DOUBLE:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.dbp++)
	    *stwspc.inbuf.p.dbp = stwspc.outbuf.p.dbp[i];
	  break;
	case WLZ_GREY_RGBA:
	  for(i=0; i < stwspc.len; i++, stwspc.inbuf.p.rgbp++)
	    *stwspc.inbuf.p.rgbp = stwspc.outbuf.p.rgbp[i];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
	}
      }
      else {
	break; /* break from the while loop */
      }
    }
  }
  if( errNum == WLZ_ERR_EOO ){
    errNum = WLZ_ERR_NONE;
  }

  /* free the buffer and return transformed object */
  if( stwspc.outbuf.p.dbp ){
    AlcFree((void *) stwspc.outbuf.p.dbp);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return obj2;
}

