#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzThreshold_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzThreshold.c
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
* \brief	Thresholds 2D or 3D domain objects with values.
* \ingroup	WlzThreshold
*/

#include <stdlib.h>
#include <Wlz.h>

static WlzObject *WlzThreshold3d(WlzObject	*obj,
				 WlzPixelV	threshV,
				 WlzThresholdType highlow,
				 WlzErrorNum	*dstErr);

/*!
* \ingroup	WlzThreshold
* \brief	Computes the number of intervals for WlzThreshold()
* 		from integral grey values or floating point if not
* 		comparing equality. The following parameters are
* 		required:
*  		<ul>
*  		  <li>N</li> Number of intervals to be computed: nints.
*  		  <li>NL1</li> Minimum line coordinate: nl1.
*  		  <li>NLL</li> Maximum line coordinate: nll.
*  		  <li>NK1</li> Minimum column coordinate nk1.
*  		  <li>NKL</li> Maximum column coordinate nkl.
*  		  <li>G</li>   The grey pointer: g.
*  		  <li>IWS</li> The interval workspace: iwsp.
*  		  <li>TV</li>  The threshold value: thresh_i.
*  		  <li>P</li>   The grey pointer union member: eg ubp.
*  		  <li>OP</li>  The comparison operator: < for WLZ_THRESH_LOW,
*  		  	       >= for WLZ_THRESH_HIGH and == for
*  		  	       WLZ_THRESH_EQUAL.
*		  <li>K</li>   Column index: colno.
*                 <li>OV</li>   Over interval: over.
*  		</ul>
*/
#define WLZ_THRESH_ADD_ITV_1(N,NL1,NLL,NK1,NKL,G,IWS,TV,P,OP,K,OV) \
{ \
  for((K)=((IWS).lftpos);(K)<=((IWS).rgtpos);++(K)) \
  { \
    if(*((G).P)OP(TV)) \
    { \
      if(!(OV)) \
      { \
	(OV)=1; \
	if(((IWS).linpos)<(NL1)) \
	{ \
	  (NL1)=((IWS).linpos); \
	} \
	if(((IWS).linpos)>(NLL)) \
	{ \
	  (NLL)=((IWS).linpos); \
	} \
	if((K)<(NK1)) \
	{ \
	  (NK1)=(K); \
	} \
      } \
    } \
    else \
    { \
      if(OV) \
      { \
	if((K)>(NKL)) \
	{ \
	  (NKL)=(K); \
	} \
	(OV)=0; \
	++(N); \
      } \
    } \
    ++((G).P); \
  } \
}

/*!
* \ingroup	WlzThreshold
* \brief	Computes the number of intervals for WlzThreshold()
* 		from floating point grey values using an equality
* 		operator. The following parameters are
* 		required:
*  		<ul>
*  		  <li>N</li> Number of intervals to be computed: nints.
*  		  <li>NL1</li> Minimum line coordinate: nl1.
*  		  <li>NLL</li> Maximum line coordinate: nll.
*  		  <li>NK1</li> Minimum column coordinate nk1.
*  		  <li>NKL</li> Maximum column coordinate nkl.
*  		  <li>G</li>   The grey pointer: g.
*  		  <li>IWS</li> The interval workspace: iwsp.
*  		  <li>TV</li>  The threshold value: thresh_i.
*  		  <li>P</li>   The grey pointer union member: eg ubp.
*  		  <li>E</li>   The epsilon value for floating point
*  		  	       comparison, should be > 0.0.
*		  <li>K</li>   Column index: colno.
*                 <li>OV</li>   Over interval: over.
*  		</ul>
*/
#define WLZ_THRESH_ADD_ITV_FE_1(N,NL1,NLL,NK1,NKL,G,IWS,TV,P,E,K,OV) \
{ \
  for((K)=((IWS).lftpos);(K)<=((IWS).rgtpos);++(K)) \
  { \
    if((*((G).P)<((TV)-(E)))||(*((G).P)>((TV)+(E)))) \
    { \
      if(OV) \
      { \
	if((K)>(NKL)) \
	{ \
	  (NKL)=(K); \
	} \
	(OV)=0; \
	++(N); \
      } \
    } \
    else \
    { \
      if(!(OV)) \
      { \
	(OV)=1; \
	if(((IWS).linpos)<(NL1)) \
	{ \
	  (NL1)=((IWS).linpos); \
	} \
	if(((IWS).linpos)>(NLL)) \
	{ \
	  (NLL)=((IWS).linpos); \
	} \
	if((K)<(NK1)) \
	{ \
	  (NK1)=(K); \
	} \
      } \
    } \
    ++((G).P); \
  } \
}

/*!
* \ingroup	WlzThreshold
* \brief	Computes the number of intervals for WlzThreshold()
* 		from RGBA grey values using a modulus operator to
* 		compute a scalar grey value from the RGBA value. The
* 		following parameters are required:
*  		<ul>
*  		  <li>N</li> Number of intervals to be computed: nints.
*  		  <li>NL1</li> Minimum line coordinate: nl1.
*  		  <li>NLL</li> Maximum line coordinate: nll.
*  		  <li>NK1</li> Minimum column coordinate nk1.
*  		  <li>NKL</li> Maximum column coordinate nkl.
*  		  <li>G</li>   The grey pointer: g.
*  		  <li>IWS</li> The interval workspace: iwsp.
*  		  <li>TV</li>  The threshold value: thresh_i.
*  		  <li>OP</li>  The comparison operator: < for WLZ_THRESH_LOW,
*  		  	       >= for WLZ_THRESH_HIGH and == for
*  		  	       WLZ_THRESH_EQUAL.
*		  <li>K</li>   Column index: colno.
*                 <li>OV</li>   Over interval: over.
*  		</ul>
*/
#define WLZ_THRESH_ADD_ITV_RGB_1(N,NL1,NLL,NK1,NKL,G,IWS,TV,OP,K,OV) \
{ \
  for((K)=((IWS).lftpos);(K)<=((IWS).rgtpos);++(K)) \
  { \
    if(WLZ_RGBA_MODULUS_2(*((G).rgbp))OP(TV)) \
    { \
      if(!(OV)) \
      { \
	(OV)=1; \
	if(((IWS).linpos)<(NL1)) \
	{ \
	  (NL1)=((IWS).linpos); \
	} \
	if(((IWS).linpos)>(NLL)) \
	{ \
	  (NLL)=((IWS).linpos); \
	} \
	if((K)<(NK1)) \
	{ \
	  (NK1)=(K); \
	} \
      } \
    } \
    else \
    { \
      if(OV) \
      { \
	if((K)>(NKL)) \
	{ \
	  (NKL)=(K); \
	} \
	(OV)=0; \
	++(N); \
      } \
    } \
    ++((G).rgbp); \
  } \
}

/*!
* \ingroup	WlzThreshold
* \brief	Constructs the intervals for WlzThreshold()
* 		from integral grey values or floating point if not
* 		comparing equality. The following parameters are
* 		required:
*  		<ul>
*  		  <li>N</li> Number of intervals to be computed: nints.
*  		  <li>NK1</li> Minimum column coordinate nk1.
*  		  <li>G</li>   The grey pointer: g.
*  		  <li>ITV</li> The interval pointer.
*  		  <li>IWS</li> The interval workspace: iwsp.
*  		  <li>TV</li>  The threshold value: thresh_i.
*  		  <li>P</li>   The grey pointer union member: eg ubp.
*  		  <li>OP</li>  The comparison operator: < for WLZ_THRESH_LOW,
*  		  	       >= for WLZ_THRESH_HIGH and == for
*  		  	       WLZ_THRESH_EQUAL.
*		  <li>K</li>   Column index: colno.
*                 <li>OV</li>   Over interval: over.
*  		</ul>
*/
#define WLZ_THRESH_ADD_ITV_2(N,NK1,G,ITV,IWS,TV,P,OP,K,OV) \
{ \
  for((K)=((IWS).lftpos);(K)<=((IWS).rgtpos);++(K)) \
  { \
    if((*(G).P)OP(TV)) \
    { \
      if(!(OV)) \
      { \
	(OV)=1; \
	((ITV)->ileft)=(K)-(NK1); \
      } \
    } \
    else \
    { \
      if(OV) \
      { \
	(OV)=0; \
	((ITV)->iright)=(K)-(NK1)-1; \
	++(N); \
	++(ITV); \
      } \
    } \
    ++((G).P); \
  } \
}

/*!
* \ingroup	WlzThreshold
* \brief	Constructs the intervals for WlzThreshold()
* 		from floating point grey values using an equality
* 		operator. The following parameters are required:
*  		<ul>
*  		  <li>N</li> Number of intervals to be computed: nints.
*  		  <li>NK1</li> Minimum column coordinate nk1.
*  		  <li>G</li>   The grey pointer: g.
*  		  <li>ITV</li> The interval pointer.
*  		  <li>IWS</li> The interval workspace: iwsp.
*  		  <li>TV</li>  The threshold value: thresh_i.
*  		  <li>P</li>   The grey pointer union member: eg ubp.
*  		  <li>E</li>   The epsilon value for floating point
*  		  	       comparison, should be > 0.0.
*		  <li>K</li>   Column index: colno.
*                 <li>OV</li>   Over interval: over.
*  		</ul>
*/
#define WLZ_THRESH_ADD_ITV_FE_2(N,NK1,G,ITV,IWS,TV,P,E,K,OV) \
{ \
  for((K)=((IWS).lftpos);(K)<=((IWS).rgtpos);++(K)) \
  { \
    if((*((G).P)<((TV)-(E)))||(*((G).P)>((TV)+(E)))) \
    { \
      if(OV) \
      { \
	(OV)=0; \
	((ITV)->iright)=(K)-(NK1)-1; \
	++(N); \
	++(ITV); \
      } \
    } \
    else \
    { \
      if(!(OV)) \
      { \
	(OV)=1; \
	((ITV)->ileft)=(K)-(NK1); \
      } \
    } \
    ++((G).P); \
  } \
}

/*!
* \ingroup	WlzThreshold
* \brief	Constructs the intervals for WlzThreshold()
* 		from RGBA grey values using a modulus operator to
* 		compute a scalar grey value from the RGBA value. The
* 		following parameters are required:
*  		<ul>
*  		  <li>N</li> Number of intervals to be computed: nints.
*  		  <li>NK1</li> Minimum column coordinate nk1.
*  		  <li>G</li>   The grey pointer: g.
*  		  <li>ITV</li> The interval pointer.
*  		  <li>IWS</li> The interval workspace: iwsp.
*  		  <li>TV</li>  The threshold value: thresh_i.
*  		  <li>P</li>   The grey pointer union member: eg ubp.
*  		  <li>OP</li>  The comparison operator: < for WLZ_THRESH_LOW,
*  		  	       >= for WLZ_THRESH_HIGH and == for
*  		  	       WLZ_THRESH_EQUAL.
*		  <li>K</li>   Column index: colno.
*                 <li>OV</li>   Over interval: over.
*  		</ul>
*/
#define WLZ_THRESH_ADD_ITV_RGB_2(N,NK1,G,ITV,IWS,TV,P,OP,K,OV) \
{ \
  for((K)=((IWS).lftpos);(K)<=((IWS).rgtpos);++(K)) \
  { \
    if((WLZ_RGBA_MODULUS_2(*(G).P))OP(TV)) \
    { \
      if(!(OV)) \
      { \
	(OV)=1; \
	((ITV)->ileft)=(K)-(NK1); \
      } \
    } \
    else \
    { \
      if(OV) \
      { \
	(OV)=0; \
	((ITV)->iright)=(K)-(NK1)-1; \
	++(N); \
	++(ITV); \
      } \
    } \
    ++((G).P); \
  } \
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzThreshold
* \brief	Thresholds a woolz grey-level object, 2D or 3D.
* \param	obj			Object to be thresholded.
* \param	threshV			Threshold pixel value.
* \param	highlow			Mode parameter with possible values:
*					<ul>
*					<li> WLZ_THRESH_LOW - thresholded
*					object is of values < given value.
*					</li>
*					<li> WLZ_THRESH_HIGH - thresholded
*					object is of values >= given value.
*					</li>
*					<li> WLZ_THRESH_EQUAL - thresholded
*					object is of values == given value.
*					</li>
*					</ul>
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzObject *WlzThreshold(WlzObject	*obj,
			WlzPixelV	threshV,
			WlzThresholdType highlow,
			WlzErrorNum	*dstErr)
{
  WlzObject		*nobj = NULL;
  WlzIntervalDomain	*idom = NULL;
  WlzGreyP		g;
  int			colno, nints;
  int			over;
  int			nl1,nll,nk1,nkl;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzInterval		*itvl = NULL, *jtvl = NULL;
  int			thresh_i;
  float			thresh_f;
  double		thresh_d;
  WlzDomain		domain;
  WlzValues		values;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  const float		eps_f = 1.0e-6;
  const double		eps_d = 1.0e-12;

  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      /* check object 2D domain and valuetable */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
	break;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzThreshold3d(obj, threshV, highlow, dstErr);

    case WLZ_TRANS_OBJ:
      if((nobj = WlzThreshold(obj->values.obj, threshV, highlow,
			      &errNum)) != NULL){
	values.obj = nobj;
	return WlzMakeMain(obj->type, obj->domain, values,
			   NULL, obj, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }

  /* check the highlow flag */
  if( errNum == WLZ_ERR_NONE ){
    switch( highlow ){

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;

    case WLZ_THRESH_LOW:
    case WLZ_THRESH_HIGH:
    case WLZ_THRESH_EQUAL:
      break;
    }
  }
    
  /* get the threshold value - this does not need to be the same
     as the valuetable pixel type */
  if( errNum == WLZ_ERR_NONE ){
    switch( threshV.type ){
    case WLZ_GREY_INT:
      thresh_i = threshV.v.inv;
      thresh_f = (float )thresh_i;
      thresh_d = thresh_i;
      break;
    case WLZ_GREY_SHORT:
      thresh_i = (int )(threshV.v.shv);
      thresh_f = (float )thresh_i;
      thresh_d = thresh_i;
      break;
    case WLZ_GREY_UBYTE:
      thresh_i = (int )(threshV.v.ubv);
      thresh_f = (float )thresh_i;
      thresh_d = thresh_i;
      break;
    case WLZ_GREY_FLOAT:
      thresh_f = threshV.v.flv;
      thresh_d = thresh_f;
      thresh_i = (int )thresh_f;
      break;
    case WLZ_GREY_DOUBLE:
      thresh_d = threshV.v.dbv;
      thresh_f = (float )thresh_d;
      thresh_i = (int )thresh_d;
      break;
    case WLZ_GREY_RGBA:
      thresh_d = WLZ_RGBA_MODULUS(threshV.v.rgbv);
      thresh_f = (float )thresh_d;
      thresh_i = (int )thresh_d;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
    }
  }
    
  /*
   * first pass - find line and column bounds of thresholded
   * object and number of intervals.
   */
  if( errNum == WLZ_ERR_NONE ){
    idom = obj->domain.i;
    nl1 = idom->lastln;
    nll = idom->line1;
    nk1 = idom->lastkl;
    nkl = idom->kol1;
    (void) WlzInitGreyScan(obj, &iwsp, &gwsp);
    nints = 0;
    if( gwsp.pixeltype == WLZ_GREY_RGBA ){
      thresh_i *= thresh_i;
    }
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
      g = gwsp.u_grintptr;
      over = 0;
      switch(gwsp.pixeltype){
	case WLZ_GREY_INT:
	  switch(highlow) {
	    case WLZ_THRESH_LOW:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,inp,
	                           <,colno,over);
	      break;
	    case WLZ_THRESH_HIGH:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,inp,
	                           >=,colno,over);
	      break;
	    case WLZ_THRESH_EQUAL:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,inp,
	                           ==,colno,over);
	      break;
	  }
	  break;
	case WLZ_GREY_SHORT:
	  switch(highlow) {
	    case WLZ_THRESH_LOW:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,shp,
	                           <,colno,over);
	      break;
	    case WLZ_THRESH_HIGH:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,shp,
	                           >=,colno,over);
	      break;
	    case WLZ_THRESH_EQUAL:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,shp,
	                           ==,colno,over);
	      break;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  switch(highlow) {
	    case WLZ_THRESH_LOW:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,ubp,
	                           <,colno,over);
	      break;
	    case WLZ_THRESH_HIGH:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,ubp,
	                           >=,colno,over);
	      break;
	    case WLZ_THRESH_EQUAL:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,ubp,
	                           ==,colno,over);
	      break;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  switch(highlow) {
	    case WLZ_THRESH_LOW:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,flp,
	                           <,colno,over);
	      break;
	    case WLZ_THRESH_HIGH:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,flp,
	                           >=,colno,over);
	      break;
	    case WLZ_THRESH_EQUAL:
	      WLZ_THRESH_ADD_ITV_FE_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,flp,
				      eps_f,colno,over);
	      break;
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  switch(highlow) {
	    case WLZ_THRESH_LOW:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,dbp,
	                           <,colno,over);
	      break;
	    case WLZ_THRESH_HIGH:
	      WLZ_THRESH_ADD_ITV_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,dbp,
	                           >=,colno,over);
	      break;
	    case WLZ_THRESH_EQUAL:
	      WLZ_THRESH_ADD_ITV_FE_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,dbp,
				      eps_d,colno,over);
	      break;
	  }
	case WLZ_GREY_RGBA:
	  switch(highlow) {
	    case WLZ_THRESH_LOW:
	      WLZ_THRESH_ADD_ITV_RGB_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,
	                               <,colno,over);
	      break;
	    case WLZ_THRESH_HIGH:
	      WLZ_THRESH_ADD_ITV_RGB_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,
	                               >=,colno,over);
	      break;
	    case WLZ_THRESH_EQUAL:
	      WLZ_THRESH_ADD_ITV_RGB_1(nints,nl1,nll,nk1,nkl,g,iwsp,thresh_i,
	                               ==,colno,over);
	      break;
	  }
	  break;
	default:
	  break;
      }

      if (over == 1) {
	if (colno > nkl)
	  nkl = colno;
	over = 0;
	nints++;
      }
    }
    nkl--;	/* since we have looked at points beyond interval ends */
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  /* domain structure */
  if( errNum == WLZ_ERR_NONE ){
    if( nints > 0 ){
      if((idom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				       nl1, nll, nk1, nkl, &errNum)) != NULL){
	if( (itvl = (WlzInterval *)
	     AlcMalloc(nints * sizeof(WlzInterval))) == NULL ){
	  errNum = WLZ_ERR_MEM_ALLOC;
	  WlzFreeIntervalDomain(idom);
	}
	else {
	  idom->freeptr = AlcFreeStackPush(idom->freeptr, (void *)itvl, NULL);
	}
      }

      /*
       * second pass - construct intervals
       */
      if( errNum == WLZ_ERR_NONE ){
	errNum = WlzInitGreyScan(obj, &iwsp, &gwsp);
	nints = 0;
	jtvl = itvl;
      }

      /* find thresholded endpoints */
      if( errNum == WLZ_ERR_NONE ){
	while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	  if( iwsp.linpos < nl1 || iwsp.linpos > nll ){
	    continue;
	  }
	  g = gwsp.u_grintptr;
	  over = 0;
	  switch(gwsp.pixeltype){
	    case WLZ_GREY_INT:
	      switch(highlow) {
		case WLZ_THRESH_LOW:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,inp,
                                       <,colno,over);
		  break;
		case WLZ_THRESH_HIGH:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,inp,
                                       >=,colno,over);
		  break;
		case WLZ_THRESH_EQUAL:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,inp,
                                       ==,colno,over);
		  break;
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      switch(highlow) {
		case WLZ_THRESH_LOW:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,shp,
                                       <,colno,over);
		  break;
		case WLZ_THRESH_HIGH:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,shp,
                                       >=,colno,over);
		  break;
		case WLZ_THRESH_EQUAL:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,shp,
                                       ==,colno,over);
		  break;
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      switch(highlow) {
		case WLZ_THRESH_LOW:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,ubp,
                                       <,colno,over);
		  break;
		case WLZ_THRESH_HIGH:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,ubp,
                                       >=,colno,over);
		  break;
		case WLZ_THRESH_EQUAL:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,ubp,
                                       ==,colno,over);
		  break;
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      switch(highlow) {
		case WLZ_THRESH_LOW:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,flp,
                                       <,colno,over);
		  break;
		case WLZ_THRESH_HIGH:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,flp,
                                       >=,colno,over);
		  break;
		case WLZ_THRESH_EQUAL:
		  WLZ_THRESH_ADD_ITV_FE_2(nints,nk1,g,itvl,iwsp,thresh_i,flp,
					  eps_f,colno,over);
		  break;
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      switch(highlow) {
		case WLZ_THRESH_LOW:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,dbp,
                                       <,colno,over);
		  break;
		case WLZ_THRESH_HIGH:
		  WLZ_THRESH_ADD_ITV_2(nints,nk1,g,itvl,iwsp,thresh_i,dbp,
                                       >=,colno,over);
		  break;
		case WLZ_THRESH_EQUAL:
		  WLZ_THRESH_ADD_ITV_FE_2(nints,nk1,g,itvl,iwsp,thresh_i,dbp,
					  eps_d,colno,over);
		  break;
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      switch(highlow) {
		case WLZ_THRESH_LOW:
		  WLZ_THRESH_ADD_ITV_RGB_2(nints,nk1,g,itvl,iwsp,thresh_i,dbp,
		                           <,colno,over);
		  break;
		case WLZ_THRESH_HIGH:
		  WLZ_THRESH_ADD_ITV_RGB_2(nints,nk1,g,itvl,iwsp,thresh_i,dbp,
		                           >=,colno,over);
		  break;
		case WLZ_THRESH_EQUAL:
		  WLZ_THRESH_ADD_ITV_RGB_2(nints,nk1,g,itvl,iwsp,thresh_i,dbp,
					   ==,colno,over);
		  break;
	      }
	      break;
	    default:
	      break;
	  }
	  if (over == 1) {
	    over = 0;
	    itvl->iright = colno - nk1 - 1;
	    nints++;
	    itvl++;
	  }
	  /*
	   * end of line ?
	   */
	  if (iwsp.intrmn == 0) {
	    WlzMakeInterval(iwsp.linpos, idom, nints, jtvl);
	    jtvl = itvl;
	    nints = 0;
	  }
	}
	if( errNum == WLZ_ERR_EOO ){
	  errNum = WLZ_ERR_NONE;
	}
      }
    } else {
      /* no thresholded points - make a dummy domain anyway */
      return WlzMakeEmpty(dstErr);
    }
  }

  /* main object */
  if( errNum == WLZ_ERR_NONE ){
    domain.i = idom;
    nobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, obj->values,
		       obj->plist, obj, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(nobj);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzThreshold
* \brief	Private function used to threshold 3D objects.
* \param	obj			Object to be thresholded.
* \param	threshV			Threshold pixel value.
* \param	highlow			Mode parameter with possible values:
*					<ul>
*					<li> WLZ_THRESH_HIGH - thresholded
*					object is of values >= threshold value.
*					</li>
*					<li> WLZ_THRESH_LOW - thresholded
*					object is of values < threshold value.
*					</li>
*					<li> WLZ_THRESH_EQUAL - thresholded
*					object is of values == given value.
*					</li>
*					</ul>
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
static WlzObject *WlzThreshold3d(WlzObject	*obj,
				 WlzPixelV	threshV,
				 WlzThresholdType highlow,
				 WlzErrorNum	*dstErr)
{
  /*----LOCAL AUTOMATIC VARIABLES-----*/
  WlzObject		*obj1;
  WlzPlaneDomain	*pdom, *npdom;
  WlzVoxelValues	*voxtab, *nvoxtab;
  WlzDomain		domain;
  WlzValues		vals;
  int			i, p, nplanes;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* no need to check the object pointer or type because this procedure
     can only be accessed via WlzThreshold. The domain and valuetable
     must be checked however */
  obj1 = NULL;
  if( obj->domain.p == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->values.vox == NULL ){
    errNum = WLZ_ERR_VALUES_NULL;
  }

  /* check types */
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->domain.p->type ){

    case WLZ_PLANEDOMAIN_DOMAIN:
      break;

    default:
      errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
      break;
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->values.vox->type ){

    case WLZ_VOXELVALUETABLE_GREY:
      break;

    default:
      errNum = WLZ_ERR_VOXELVALUES_TYPE;
      break;
    }
  }

  /* make new planedomain and voxelvaluetable */
  if( errNum == WLZ_ERR_NONE ){
    pdom = obj->domain.p;
    voxtab = obj->values.vox;
    npdom = WlzMakePlaneDomain(pdom->type,
			       pdom->plane1, pdom->lastpl,
			       pdom->line1, pdom->lastln,
			       pdom->kol1, pdom->lastkl, &errNum);
  }
    
  if((errNum == WLZ_ERR_NONE) &&
     ((nvoxtab = WlzMakeVoxelValueTb(voxtab->type, voxtab->plane1,
				     voxtab->lastpl, voxtab->bckgrnd,
				     NULL, &errNum)) == NULL) ){
    WlzFreePlaneDomain(npdom);
  }

  if( errNum == WLZ_ERR_NONE ){
    /* copy voxel_sizes */
    for(i=0; i < 3; i++){
      npdom->voxel_size[i] = pdom->voxel_size[i];
    }
    /* threshold each plane */
    nplanes = pdom->lastpl - pdom->plane1 + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(p = 0; p < nplanes; ++p){
      if(errNum == WLZ_ERR_NONE){
  	WlzDomain	*domains,
			*ndomains;
  	WlzValues	*values,
			*nvalues;
	WlzObject	*gObj2D = NULL,
			*tObj2D = NULL;
	WlzErrorNum 	errNum2D = WLZ_ERR_NONE;

	domains = pdom->domains + p;
	values = voxtab->values + p;
	ndomains = npdom->domains + p;
	nvalues = nvoxtab->values + p;
	if(((*domains).core == NULL) || ((*values).core == NULL)){
	  (*ndomains).core = NULL;
	  (*nvalues).core = NULL;
	}
	else if((gObj2D = WlzAssignObject(
	                  WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				      NULL, NULL, &errNum2D), NULL)) != NULL){

	  if( gObj2D->domain.i != NULL ){
	    if((tObj2D = WlzThreshold(gObj2D, threshV, highlow,
	                              &errNum2D)) != NULL){
	      if(tObj2D->type == WLZ_2D_DOMAINOBJ){
		*ndomains = WlzAssignDomain(tObj2D->domain, NULL);
		*nvalues = WlzAssignValues(tObj2D->values, NULL);
	      }
	      else {
		(*ndomains).core = NULL;
		(*nvalues).core = NULL;
	      }
	      (void )WlzFreeObj(tObj2D);
	    }
	  } else {
	    (*ndomains).core = NULL;
	    (*nvalues).core = NULL;
	  }
	  (void )WlzFreeObj(gObj2D);
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
    if(errNum != WLZ_ERR_NONE){
      WlzFreePlaneDomain(npdom);
      WlzFreeVoxelValueTb(nvoxtab);
    }
  }
  /* standardise the plane domain */
  if((errNum == WLZ_ERR_NONE) &&
     ((errNum = WlzStandardPlaneDomain(npdom, nvoxtab)) != WLZ_ERR_NONE) ){
    WlzFreePlaneDomain( npdom );
    WlzFreeVoxelValueTb( nvoxtab );
  }

  /* return a new object */
  if( errNum == WLZ_ERR_NONE ){
    domain.p = npdom;
    vals.vox = nvoxtab;
    if((obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, vals,
			   NULL, obj, &errNum)) != NULL){
      /*	nvoxtab->original = obj1; */
      nvoxtab->original_table = WlzAssignValues(obj->values, NULL);
    }
    else {
      WlzFreePlaneDomain( npdom );
      WlzFreeVoxelValueTb( nvoxtab );
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return obj1;
}
