#pragma ident "MRC HGU $Id$"
/************************************************************************
*   Copyright  :   1994 Medical Research Council, UK.                   *
*                  All rights reserved.                                 *
*************************************************************************
*   Address    :   MRC Human Genetics Unit,                             *
*                  Western General Hospital,                            *
*                  Edinburgh, EH4 2XU, UK.                              *
*************************************************************************
*   Project    :   Woolz Library					*
*   File       :   WlzEmpty.c						*
*************************************************************************
* This module has been copied from the original woolz library and       *
* modified for the public domain distribution. The original authors of  *
* the code and the original file headers and comments are in the        *
* HISTORY file.                                                         *
*************************************************************************
*   Author Name :  Richard Baldock					*
*   Author Login:  richard@hgu.mrc.ac.uk				*
*   Date        :  Mon Oct 30 17:34:21 2000				*
*   $Revision$							*
*   $Name$								*
*   Synopsis    : 							*
*************************************************************************
*   Maintenance :  date - name - comments (Last changes at the top)	*
************************************************************************/

#include <Wlz.h>

 
int WlzIsEmpty(WlzObject *obj, WlzErrorNum *wlzErr)
{
  WlzObject		tmpobj;
  WlzDomain		domain;
  WlzPlaneDomain	*pldom;
  int			i, p, emptyFlg = 0;
  WlzErrorNum		errNum = WLZ_ERR_NONE; 

  /* check the object */
  if( obj == NULL )
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else {
	switch( obj->domain.i->type ){
	case WLZ_INTERVALDOMAIN_INTVL:
	  for(i=obj->domain.i->line1; i <=  obj->domain.i->lastln; i++){
	    if( obj->domain.i->intvlines[i-obj->domain.i->line1].nintvs > 0 ){
	      emptyFlg = 0;
	      break;
	    }
	    else {
	      emptyFlg = 1;
	    }
	  }
	  break;

	default:
	  break;
	}
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else {
	/* scan through planes */
	pldom = obj->domain.p;
	tmpobj.type = WLZ_2D_DOMAINOBJ;
	tmpobj.values.core = NULL;
	tmpobj.plist = NULL;
	tmpobj.assoc = NULL;
	for(p=0; (p <= (pldom->lastpl - pldom->plane1)) &&
	      (errNum == WLZ_ERR_NONE); p++){
	  domain = pldom->domains[p];
	  if( domain.core ){
	    tmpobj.domain = domain;
	    if( WlzIsEmpty(&tmpobj, &errNum) ){
	      emptyFlg = 1;
	    }
	    else {
	      emptyFlg = 0;
	      break;
	    }
	  }
	  else {
	    emptyFlg = 1;
	  }
	}
      }
      break;

    case WLZ_EMPTY_OBJ:
      emptyFlg = 1;
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if(wlzErr)
  {
    *wlzErr = errNum;
  }

  return emptyFlg;
}
