#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzVerifyObj.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed Sep 24 17:15:06 2003
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
* \ingroup      WlzError
* \brief        Functions for verifying Woolz objects.
*               
* \todo         Extend verification beyond just 2D and 3D domain
 objects.
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>

/* function:     WlzVerifyObject    */
/*! 
* \ingroup      WlzError
* \brief        Verify the object data, fix if possible if paramter
<tt> fix != 0</tt>. Curretnly only the domains of 2D and 3D objects
 can be verified.
*
* \return       Error detected. If fix != 0 then the error may be fixed.
* \param    obj	Input object to be verified
* \param    fix	If fix != 0 then attempt to fix the error detected.
* \par      Source:
*                WlzVerifyObj.c
*/
WlzErrorNum WlzVerifyObject(
  WlzObject	*obj, 
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
	if( domains[p].core == NULL ){
	  continue;
	}
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

/* function:     WlzVerifyIntervalDomain    */
/*! 
* \ingroup      WlzError
* \brief        Verify the input interval domain, check domain parameters,
 intervals etc. Fix errors if requested.
*
* \return       Error detected in the woolz interval domain. Possible values
 are <tt>WLZ_ERR_INTERVALDOMAIN_NULL, WLZ_ERR_LINE_DATA, WLZ_ERR_COLUMN_DATA,
 WLZ_ERR_NONE </tt> and errors from WlzVerifyIntervalLine().
* \param    dom	Input interval domain to be checked.
* \param    fix	Attempt to fix errors if <tt>fix != 0</tt>
* \par      Source:
*                WlzVerifyObj.c
*/
WlzErrorNum  WlzVerifyIntervalDomain(
  WlzDomain	dom,
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

/* function:     WlzVerifyIntervalLine    */
/*! 
* \ingroup      WlzError
* \brief        Detect errors in an interval line structure, fix if
 requested.
*
* \return       Detected error, one of <tt>WLZ_ERR_NONE, WLZ_ERR_INTERVALLINE_NULL, WLZ_ERR_INTERVAL_ADJACENT</tt> or errors returned by WlzVerifyInterval().
* \param    intvline	Interval line to be checked
* \param    fix	Fix error if parameter <tt>fix != 0</tt>
* \par      Source:
*                WlzVerifyObj.c
*/
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

/* function:     WlzVerifyInterval    */
/*! 
* \ingroup      WlzError
* \brief        Check an interval structure, fix if requested.
*
* \return       Error detected in interval, one of <tt>WLZ_ERR_INTERVAL_NULL, WLZ_ERR_INTERVAL_DATA, WLZ_ERR_INTERVAL_BOUND, WLZ_ERR_NONE</tt>.
* \param    intv	Interval structure to be tested.
* \param    fix	Fix error if parameter <tt>fix != 0</tt>.
* \par      Source:
*                WlzVerifyObj.c
*/
WlzErrorNum WlzVerifyInterval(
  WlzInterval	*intv,
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
