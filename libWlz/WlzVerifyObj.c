#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzVerifyObj.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for verifying Woolz objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

/************************************************************************
*   Function   : WlzVerifyObject					*
*   Date       : Fri Nov 29 09:32:54 1996				*
*************************************************************************
*   Synopsis   :verify the object data, fix if possible and fix!=0	*
*		currently only the domains of 2D and 3D domain objects	*
*		can be verified. This must be extended to all objects.	*
*   Returns    :WlzErrorNum: all errors are reported but the returned	*
*		error will depend on whether the error could be fixed	*
*		if fix == 0 then the error is not fixed and returned	*
*   Parameters :WlzObject	*obj: object pointer to be checked	*
*		int 		fix: if != 0 then attempt to fix the	*
*			error.						*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum 
WlzVerifyObject(WlzObject	*obj, 
		int 		fix)
{
  /* local variables */
  WlzObject	*tmpobj;
  WlzErrorNum	wlzerrno=WLZ_ERR_NONE;
  WlzDomain	*domains=NULL;
  WlzValues	*values =NULL, vals;
  int		p, nplanes;

  /* check pointer */
  if( obj == NULL ){
    return(WLZ_ERR_OBJECT_NULL);
  }

  /* check domains */
  switch( obj->type ){

  case WLZ_2D_DOMAINOBJ:
    wlzerrno = WlzVerifyIntervalDomain(obj->domain, fix);
    break;

  case WLZ_3D_DOMAINOBJ:
    switch( obj->domain.p->type ){

    case WLZ_PLANEDOMAIN_DOMAIN:
      nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
      domains = obj->domain.p->domains;
      if( obj->values.core ){
	values = obj->values.vox->values;
      }
      vals.core = NULL;
      for(p=0; p < nplanes; p++){
	if( values ){
	  tmpobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[p], values[p],
			       NULL, NULL, NULL);
	}
	else {
	  tmpobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains[p], vals,
			       NULL, NULL, NULL);
	}
	if( WlzVerifyObject(tmpobj, fix) != WLZ_ERR_NONE ){
	  wlzerrno = WLZ_ERR_PLANEDOMAIN_DATA;
	}
	if( tmpobj ){
	  WlzFreeObj( tmpobj );
	}
      }
      break;

    default:
      break;

    }
    break;

  default:
    break;

  }

  return( wlzerrno );
}



/************************************************************************
*   Function   : WlzVerifyIntervalDomain				*
*   Returns    :WlzErrorNum: one of WLZ_ERR_INTERVALDOMAIN_NULL, WLZ_ERR_LINE_DATA*
*		WLZ_ERR_COLUMN_DATA, WLZ_ERR_NONE and errors from 	*
*		WlzVerifyIntervalLine.					*
*   Parameters :WlzDomain	dom: woolz domain union			*
*		int 		fix: if != zero attempt to fix		*
*   Date       : Mon Oct 14 15:47:58 1996				*
*   Synopsis   :							*
*	verify_intervaldomain - check domain parameters are valid, fix
*	if possible for if "fix" > 0
************************************************************************/

WlzErrorNum 
WlzVerifyIntervalDomain(WlzDomain	dom,
			int 		fix)
{
  /* local variables */
  WlzIntervalLine 	*intvlines;
  WlzIntervalDomain 	*idom = dom.i;
  int 		i;
  WlzErrorNum	wlzerrno = WLZ_ERR_NONE;

  /* check pointer */
  if( idom == NULL ){
    return(WLZ_ERR_INTERVALDOMAIN_NULL);
  }

  /* check line and column values */
  if( idom->line1 > idom->lastln ){
    if( fix ){
      idom->lastln = idom->line1;
    } else {
      return(WLZ_ERR_LINE_DATA);
    }
  }
  if( idom->kol1 > idom->lastkl ){
    if( fix ){
      idom->lastkl = idom->kol1;
    } else {
      return(WLZ_ERR_COLUMN_DATA);
    }
  }

  /* check intervals if type 1 domain */
  switch( idom->type ){

  case WLZ_INTERVALDOMAIN_INTVL:
    intvlines = idom->intvlines;
    for(i=idom->line1; i <= idom->lastln; i++, intvlines++){
      if( (wlzerrno = WlzVerifyIntervalLine(intvlines, fix)) != WLZ_ERR_NONE){
	return( wlzerrno );
      }
    }
    break;

  default:
    break;

  }

  /* if we have got this far then all ok or fixed */
  return( wlzerrno );
}

/************************************************************************
*   Function   : WlzVerifyIntervalLine					*
*   Returns    :WlzErrorNum: one of WLZ_ERR_NONE, WLZ_ERR_INTERVALLINE_NULL,	*
*		WLZ_ERR_INTERVAL_ADJACENT or errors returned by		*
*		WlzVerifyInterval.					*
*   Parameters :WlzIntervalLine *intvline: interval line structure	*
*		int fix: if != 0 then attempt to fix.			*
*   Date       : Mon Oct 14 15:49:46 1996				*
*   Synopsis   :							*
************************************************************************/

WlzErrorNum 
WlzVerifyIntervalLine(WlzIntervalLine *intvline,
		      int fix)
{
  /* local variables */
  WlzInterval	*intvs;
  int		i, j;
  WlzErrorNum	wlzerrno = WLZ_ERR_NONE;

  /* check pointer */
  if( intvline == NULL ){
    return(WLZ_ERR_INTERVALLINE_NULL);
  }

  /* check number of intervals */
  if( intvline->nintvs < 0 ){
    if( fix ){
      intvline->nintvs = 0;
    } else {
      return(WLZ_ERR_INTERVAL_NUMBER);
    }
  }

  /* check each interval */
  intvs = intvline->intvs;
  for(i=0; i < intvline->nintvs; i++){
    if( (wlzerrno = WlzVerifyInterval(intvs+i, fix)) != 0 ){
      return(wlzerrno);
    }
    if( i ){
      /* check for touching or overlap */
      if( (intvs[i].ileft - intvs[i-1].iright) < 2 ){
	if( fix ){
	  /* remove the current interval */
	  intvs[i-1].iright = intvs[i].iright;
	  for(j=i; j < (intvline->nintvs - 1); j++){
	    intvs[j] = intvs[j+1];
	  }
	  intvline->nintvs--;
	  i--;
	} else {
	  return(WLZ_ERR_INTERVAL_ADJACENT);
	}
      }
    }
  }

  return( wlzerrno );
}

/************************************************************************
*   Function   : WlzVerifyInterval					*
*   Returns    :WlzErrorNum: one of WLZ_ERR_INTERVAL_NULL, WLZ_ERR_INTERVAL_DATA,	*
*		WLZ_ERR_INTERVAL_BOUND, WLZ_ERR_NONE			*
*   Parameters :WlzInterval	*intv: interval pointer to test		*
*		int fix: if != 0 then attempt to fix			*
*   Date       : Mon Oct 14 15:50:54 1996				*
*   Synopsis   :							*
************************************************************************/

WlzErrorNum
WlzVerifyInterval(WlzInterval	*intv,
		  int 		fix)
{
  WlzErrorNum	wlzerrno=WLZ_ERR_NONE;

  /* check pointer */
  if( intv == NULL ){
    return(WLZ_ERR_INTERVAL_NULL);
  }

  /* check interval lower bound */
  if( intv->ileft < 0 ){
    if( fix ){
      intv->ileft = 0;
    } else {
      return(WLZ_ERR_INTERVAL_BOUND);
    }
  }

  /* check interval */
  if( intv->ileft > intv->iright ){
    if( fix ){
      intv->ileft = intv->iright;
    } else {
      return(WLZ_ERR_INTERVAL_DATA);
    }
  }

  return( wlzerrno );
}
