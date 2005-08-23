#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzEmpty.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Sep 24 08:03:21 2003
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
* \ingroup      WlzDomainOps
* \brief        Convenience procedure to check empty status of an object.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>


/* function:     WlzIsEmpty    */
/*! 
* \ingroup      WlzDomainOps
* \brief        Convenience procedure to check if an object is empty.
This include objects with zero area or volume.
*
* \return       Zero if non-empty, 1 if empty.
* \param    obj	Input object
* \param    wlzErr	Error return
* \par      Source:
*                WlzEmpty.c
*/ 
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
