#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzSepTrans.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Executes a separable image transform on a 2D Woolz
*		domain object. It's up to the user to make sure the
*		resultant is a legal value for the grey-value type.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzSepTrans						*
*   Date       : Thu Jul 17 19:58:46 1997				*
*************************************************************************
*   Synopsis   :Perform separable transform on a grey-level image	*
*   Returns    :WlzObject *: transformed object				*
*   Parameters :WlzObject *obj: object for transform			*
*		WlzIntervalConvFunc	x_fun: convolution function to	*
*			be applied the the intervals - x-direction	*
*		void		*x_params: parameter pointer to be	*
*			passed to the convolution function		*
*		WlzIntervalConvFunc	y_fun: convolution function for	*
*			the y-direction					*
*		
*   Global refs:None.							*
************************************************************************/

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
      if( obj1 = WlzSepTrans(obj->values.obj,
			     x_fun, x_params, y_fun, y_params, &errNum) ){
	values.obj = obj1;
	return WlzMakeMain(obj->type, obj->domain, values,
			   NULL, obj, dstErr);
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
      (*y_fun)(&stwspc, y_params);
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
    while( (errNum = WlzNextGreyInterval(&iwspace)) == WLZ_ERR_NONE ){
      stwspc.inbuf.p.inp = gwspace.u_grintptr.inp;
      stwspc.len = iwspace.rgtpos - iwspace.lftpos + 1;
      (*x_fun)(&stwspc, x_params);
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

