#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzGauss.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Mon May 26 17:37:52 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzValuesFilters
* \brief        Apply a Gaussian filter to a 2D grey-level image.
 Uses WlzSepTrans and for colour images can only do smoothing correctly
 (i.e. derivative zero).
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

#define AFACTOR	100

/* function:     WlzGauss2    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        Gaussian filter of grey-level 2D woolz object. x- and
 y-coordinate width parameters and derivative degree can be independently
 specified. For derivative zero, i.e. Gaussian smoothing, the filter is
 normalised. Derivatives are derivative of the normalised filter. RGB
 data will only return values for smoothing, higher derivatives are not
 implemented. The width parameter is the full-width half-height of the
 Gaussian distribution.
*
* \return       Pointer to transformed object
* \param    obj	Input object
* \param    wx	x-direction width parameter
* \param    wy	y-direction width parameter
* \param    x_deriv	x-direction derivative
* \param    y_deriv	y-direction derivative
* \param    wlzErr	error return
* \par      Source:
*                WlzGauss.c
*/
WlzObject *WlzGauss2(
  WlzObject	*obj,
  double	wx,
  double	wy,
  int		x_deriv,
  int		y_deriv,
  WlzErrorNum	*wlzErr)
{
  WlzObject		*newobj;
  Wlz1DConvMask		x_params, y_params;
  float 		alpha, sum;
  int 			i, n, value;
  WlzErrorNum		errNum = WLZ_ERR_UNSPECIFIED;
    
  /* check object, don't need to check type etc. because WlzSepTrans
     does it */
  if( obj == NULL )
  {
    errNum = WLZ_ERR_OBJECT_NULL;
    if(wlzErr)
    {
      *wlzErr = errNum;
    }
    return NULL;
  }

  /* do need to check for rgb grey type */
  if( WlzGreyTypeFromObj(obj, &errNum) == WLZ_GREY_RGBA ){
    if( (x_deriv != 0) || (y_deriv != 0) ){
      errNum = WLZ_ERR_UNIMPLEMENTED;
    }
  }
  if( errNum != WLZ_ERR_NONE ){
    if(wlzErr)
    {
      *wlzErr = errNum;
    }
    return NULL;
  }    
    
  alpha = (float) 4.0 * log( (double) 2.0 );
    
  /* set up x function parameters */
  x_params.mask_size = (((int) wx * 4)/2)*2 + 1;
  if( (x_params.mask_values = (int *)
       AlcMalloc(sizeof(int) * x_params.mask_size)) == NULL){
    errNum = WLZ_ERR_MEM_ALLOC;
    if(wlzErr)
    {
      *wlzErr = errNum;
    }
    return NULL;
  }
  n = x_params.mask_size / 2;
    
  switch( x_deriv ){

  case 0:
    for(i=0, sum = -AFACTOR; i <= n; i++){
      value = (int) (AFACTOR * exp(((double) -alpha*i*i/wx/wx)));
      *(x_params.mask_values+n-i) = value;
      *(x_params.mask_values+n+i) = value;
      sum += 2 * value;
    }
    x_params.norm_factor = sum;
    break;

  case 1:
    *(x_params.mask_values+n) = 0.0;
    for(i=1, sum = 0; i <= n; i++){
      value = (int) AFACTOR * i * exp(((double) -alpha*i*i/wx/wx));
      *(x_params.mask_values+n-i) = value;
      *(x_params.mask_values+n+i) = -value;
      sum += value;
    }
/*    sum *= -wx / 2 / sqrt( log( (double) 2 ) / WLZ_M_PI );*/
    if( n > 0 )
      x_params.norm_factor = sum;
    else
      x_params.norm_factor = 1;
    break;

  case 2:
    for(i=0; i <= n; i++){
      value = (int) AFACTOR * (alpha * i*i / wx/wx -1) *
	exp(((double) -alpha*i*i/wx/wx));
      *(x_params.mask_values+n-i) = value;
      *(x_params.mask_values+n+i) = value;
    }
    x_params.norm_factor = *(x_params.mask_values+n) * wx*wx*wx / 4 / alpha
      / sqrt( log( (double) 2 ) / WLZ_M_PI );
    break;

  default:
    AlcFree((void *) x_params.mask_values);
    errNum = WLZ_ERR_PARAM_DATA;
    if(wlzErr)
    {
      *wlzErr = errNum;
    }
    return NULL;

  }
    
  /* set up y function parameters */
  y_params.mask_size = (((int) wy * 4)/2)*2 + 1;
  if( (y_params.mask_values = (int *)
       AlcMalloc(sizeof(int) * y_params.mask_size)) == NULL){
    AlcFree((void *) x_params.mask_values);
    errNum = WLZ_ERR_MEM_ALLOC;
    if(wlzErr)
    {
      *wlzErr = errNum;
    }
    return NULL;
  }
  n = y_params.mask_size / 2;
    
  switch( y_deriv ){

  case 0:
    for(i=0, sum = -AFACTOR; i <= n; i++){
      value = (int) AFACTOR * exp(((double) -alpha*i*i/wy/wy));
      *(y_params.mask_values+n-i) = value;
      *(y_params.mask_values+n+i) = value;
      sum += 2 * value;
    }
    y_params.norm_factor = sum;
    break;

  case 1:
    *(y_params.mask_values+n) = 0.0;
    for(i=1, sum = 0; i <= n; i++){
      value = (int) AFACTOR * i * exp(((double) -alpha*i*i/wy/wy));
      *(y_params.mask_values+n-i) = value;
      *(y_params.mask_values+n+i) = -value;
      sum += value;
    }
/*    sum *= -wy / 2 / sqrt( log( (double) 2 ) /WLZ_M_PI );*/
    if( n > 0 )
      y_params.norm_factor = sum;
    else
      y_params.norm_factor = 1;
    break;

  case 2:
    for(i=0; i <= n; i++){
      value = (int) AFACTOR * (alpha * i*i / wy/wy -1) *
	exp(((double) -alpha*i*i/wy/wy));
      *(y_params.mask_values+n-i) = value;
      *(y_params.mask_values+n+i) = value;
    }
    y_params.norm_factor = *(y_params.mask_values+n) * wy*wy*wy / 4 / alpha
      / sqrt( log( (double) 2 ) / WLZ_M_PI );
    break;

  default:
    AlcFree((void *) x_params.mask_values);
    AlcFree((void *) y_params.mask_values);
    errNum = WLZ_ERR_PARAM_DATA;
    if(wlzErr)
    {
      *wlzErr = errNum;
    }
    return NULL;

  }
    
  newobj = WlzSepTrans(obj,
		       Wlz1DConv, (void *) &x_params,
		       Wlz1DConv, (void *) &y_params,
		       &errNum);
  AlcFree((void *) x_params.mask_values);
  AlcFree((void *) y_params.mask_values);
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  return(newobj);
}
  

/* function:     Wlz1DConv    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        Perform a 1D convolution on a 1D array of data. Typically
 used to pass to WlzSepTrans(). The params variable is a 1D convolution mask.
*
* \return       Woolz error
* \param    stwspc	Separable transfom work space
* \param    params	parameters passed from top-level calling function,
 unchanged by WlzSepTrans()
* \par      Source:
*                WlzGauss.c
*/  
WlzErrorNum Wlz1DConv(
  WlzSepTransWSpace	*stwspc,
  void			*params)
{
  Wlz1DConvMask	*convParams = (Wlz1DConvMask *) params;
  int 		i, j, n, *mask, factor, length;
  int		intSum;
  double	dblSum;
  WlzGreyP	inbuf, outbuf;
  UINT		red, green, blue;
    
  /* set some local parameters */
  n = convParams->mask_size / 2;
  mask = convParams->mask_values + n;
  factor = convParams->norm_factor;
  inbuf = stwspc->inbuf.p;
  outbuf = stwspc->outbuf.p;
  length = stwspc->len;

  /* calculate the new value  - use int for UBYTE, short and int
     double otherwise, separate rgb values each use UINT */
  switch( stwspc->inbuf.type ){

  case WLZ_GREY_INT:

    /* first convolve up to the half-width of the mask */
    for(i=0; (i < n) && (i < length); i++){
      intSum = *inbuf.inp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.inp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	intSum += inbuf.inp[-((i>j)?j:i)] * mask[-j];
      }
      inbuf.inp++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }

    /* now the central portion */
    while( i < (length-n-1) ){
      intSum = *inbuf.inp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.inp[j] * mask[j];
	intSum += inbuf.inp[-j] * mask[-j];
      }
      inbuf.inp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }   

    /* now the last bit within a half-width of the end */
    while( i < length ){
      intSum = *inbuf.inp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.inp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	intSum += inbuf.inp[-j] * mask[-j];
      }
      inbuf.inp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }      
    break;

  case WLZ_GREY_SHORT:

    /* first convolve up to the half-width of the mask */
    for(i=0; (i < n) && (i < length); i++){
      intSum = *inbuf.shp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.shp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	intSum += inbuf.shp[-((i>j)?j:i)] * mask[-j];
      }
      inbuf.shp++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }

    /* now the central portion */
    while( i < (length-n-1) ){
      intSum = *inbuf.shp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.shp[j] * mask[j];
	intSum += inbuf.shp[-j] * mask[-j];
      }
      inbuf.shp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }   

    /* now the last bit within a half-width of the end */
    while( i < length ){
      intSum = *inbuf.shp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.shp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	intSum += inbuf.shp[-j] * mask[-j];
      }
      inbuf.shp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }      
    break;

  case WLZ_GREY_UBYTE:

    /* first convolve up to the half-width of the mask */
    for(i=0; (i < n) && (i < length); i++){
      intSum = *inbuf.ubp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.ubp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	intSum += inbuf.ubp[-((i>j)?j:i)] * mask[-j];
      }
      inbuf.ubp++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }

    /* now the central portion */
    while( i < (length-n-1) ){
      intSum = *inbuf.ubp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.ubp[j] * mask[j];
	intSum += inbuf.ubp[-j] * mask[-j];
      }
      inbuf.ubp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }   

    /* now the last bit within a half-width of the end */
    while( i < length ){
      intSum = *inbuf.ubp * mask[0];
      for(j=1; j <= n; j++){
	intSum += inbuf.ubp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	intSum += inbuf.ubp[-j] * mask[-j];
      }
      inbuf.ubp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = intSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = intSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) intSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = intSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = intSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(intSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }      
    break;

  case WLZ_GREY_FLOAT:

    /* first convolve up to the half-width of the mask */
    for(i=0; (i < n) && (i < length); i++){
      dblSum = *inbuf.flp * mask[0];
      for(j=1; j <= n; j++){
	dblSum += inbuf.flp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	dblSum += inbuf.flp[-((i>j)?j:i)] * mask[-j];
      }
      inbuf.flp++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = dblSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = dblSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) dblSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = dblSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = dblSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(dblSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }

    /* now the central portion */
    while( i < (length-n-1) ){
      dblSum = *inbuf.flp * mask[0];
      for(j=1; j <= n; j++){
	dblSum += inbuf.flp[j] * mask[j];
	dblSum += inbuf.flp[-j] * mask[-j];
      }
      inbuf.flp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = dblSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = dblSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) dblSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = dblSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = dblSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(dblSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }   

    /* now the last bit within a half-width of the end */
    while( i < length ){
      dblSum = *inbuf.flp * mask[0];
      for(j=1; j <= n; j++){
	dblSum += inbuf.flp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	dblSum += inbuf.flp[-j] * mask[-j];
      }
      inbuf.flp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = dblSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = dblSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) dblSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = dblSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = dblSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(dblSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }      
    break;

  case WLZ_GREY_DOUBLE:

    /* first convolve up to the half-width of the mask */
    for(i=0; (i < n) && (i < length); i++){
      dblSum = *inbuf.dbp * mask[0];
      for(j=1; j <= n; j++){
	dblSum += inbuf.dbp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	dblSum += inbuf.dbp[-((i>j)?j:i)] * mask[-j];
      }
      inbuf.dbp++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = dblSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = dblSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) dblSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = dblSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = dblSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(dblSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }

    /* now the central portion */
    while( i < (length-n-1) ){
      dblSum = *inbuf.dbp * mask[0];
      for(j=1; j <= n; j++){
	dblSum += inbuf.dbp[j] * mask[j];
	dblSum += inbuf.dbp[-j] * mask[-j];
      }
      inbuf.dbp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = dblSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = dblSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) dblSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = dblSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = dblSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(dblSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }   

    /* now the last bit within a half-width of the end */
    while( i < length ){
      dblSum = *inbuf.dbp * mask[0];
      for(j=1; j <= n; j++){
	dblSum += inbuf.dbp[(j+i)>(length-1)?-i+length-1:j] * mask[j];
	dblSum += inbuf.dbp[-j] * mask[-j];
      }
      inbuf.dbp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = dblSum/factor;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = dblSum/factor;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) dblSum/factor;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = dblSum/factor;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = dblSum/factor;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(dblSum/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, red, red, 255);
	*outbuf.rgbp++;
	break;
      }
    }      
    break;

  case WLZ_GREY_RGBA:

    /* first convolve up to the half-width of the mask */
    for(i=0; (i < n) && (i < length); i++){
      red = WLZ_RGBA_RED_GET(*inbuf.rgbp) * mask[0];
      green = WLZ_RGBA_GREEN_GET(*inbuf.rgbp) * mask[0];
      blue = WLZ_RGBA_BLUE_GET(*inbuf.rgbp) * mask[0];
      for(j=1; j <= n; j++){
	red += WLZ_RGBA_RED_GET(inbuf.rgbp[(j+i)>(length-1)?-i+length-1:j])
	  * mask[j];
	red += WLZ_RGBA_RED_GET(inbuf.rgbp[-((i>j)?j:i)]) * mask[-j];
	green += WLZ_RGBA_GREEN_GET(inbuf.rgbp[(j+i)>(length-1)?-i+length-1:j])
	  * mask[j];
	green += WLZ_RGBA_GREEN_GET(inbuf.rgbp[-((i>j)?j:i)]) * mask[-j];
	blue += WLZ_RGBA_BLUE_GET(inbuf.rgbp[(j+i)>(length-1)?-i+length-1:j])
	  * mask[j];
	blue += WLZ_RGBA_BLUE_GET(inbuf.rgbp[-((i>j)?j:i)]) * mask[-j];
      }
      inbuf.rgbp++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(red/factor, 0, 255);
	green = WLZ_CLAMP(green/factor, 0, 255);
	blue = WLZ_CLAMP(blue/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, green, blue, 255);
	*outbuf.rgbp++;
	break;
      }
    }

    /* now the central portion */
    while( i < (length-n-1) ){
      red = WLZ_RGBA_RED_GET(*inbuf.rgbp) * mask[0];
      green = WLZ_RGBA_GREEN_GET(*inbuf.rgbp) * mask[0];
      blue = WLZ_RGBA_BLUE_GET(*inbuf.rgbp) * mask[0];
      for(j=1; j <= n; j++){
	red += WLZ_RGBA_RED_GET(inbuf.rgbp[j]) * mask[j];
	red += WLZ_RGBA_RED_GET(inbuf.rgbp[-j]) * mask[-j];
	green += WLZ_RGBA_GREEN_GET(inbuf.rgbp[j]) * mask[j];
	green += WLZ_RGBA_GREEN_GET(inbuf.rgbp[-j]) * mask[-j];
	blue += WLZ_RGBA_BLUE_GET(inbuf.rgbp[j]) * mask[j];
	blue += WLZ_RGBA_BLUE_GET(inbuf.rgbp[-j]) * mask[-j];
      }
      inbuf.rgbp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(red/factor, 0, 255);
	green = WLZ_CLAMP(green/factor, 0, 255);
	blue = WLZ_CLAMP(blue/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, green, blue, 255);
	*outbuf.rgbp++;
	break;
      }
    }   

    /* now the last bit within a half-width of the end */
    while( i < length ){
      red = WLZ_RGBA_RED_GET(*inbuf.rgbp) * mask[0];
      green = WLZ_RGBA_GREEN_GET(*inbuf.rgbp) * mask[0];
      blue = WLZ_RGBA_BLUE_GET(*inbuf.rgbp) * mask[0];
      for(j=1; j <= n; j++){
	red += WLZ_RGBA_RED_GET(inbuf.rgbp[(j+i)>(length-1)?-i+length-1:j])
	  * mask[j];
	red += WLZ_RGBA_RED_GET(inbuf.rgbp[-j]) * mask[-j];
	green += WLZ_RGBA_GREEN_GET(inbuf.rgbp[(j+i)>(length-1)?-i+length-1:j])
	  * mask[j];
	green += WLZ_RGBA_GREEN_GET(inbuf.rgbp[-j]) * mask[-j];
	blue += WLZ_RGBA_BLUE_GET(inbuf.rgbp[(j+i)>(length-1)?-i+length-1:j])
	  * mask[j];
	blue += WLZ_RGBA_BLUE_GET(inbuf.rgbp[-j]) * mask[-j];
      }
      inbuf.rgbp++;
      i++;
      switch( stwspc->outbuf.type ){
      case WLZ_GREY_INT:
	*outbuf.inp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_SHORT:
	*outbuf.shp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_UBYTE:
	*outbuf.ubp++ = (UBYTE) (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_FLOAT:
	*outbuf.flp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_DOUBLE:
	*outbuf.dbp++ = (red+green+blue)/factor/3.0;
	break;
      case WLZ_GREY_RGBA:
	red = WLZ_CLAMP(red/factor, 0, 255);
	green = WLZ_CLAMP(green/factor, 0, 255);
	blue = WLZ_CLAMP(blue/factor, 0, 255);
	WLZ_RGBA_RGBA_SET(*outbuf.rgbp, red, green, blue, 255);
	*outbuf.rgbp++;
	break;
      }
    }      
    break;

  }

  return WLZ_ERR_NONE;
}
