#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyScan.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz object grey value scanning functions.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <limits.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzInitGreyScan					*
*   Returns    :							*
*   Parameters :							*
*   Date       : Mon Oct 14 12:19:30 1996				*
*   Synopsis   :							*
* initialise interval and grey scanning of standard object in standard  *
* direction                                                             *
************************************************************************/

WlzErrorNum 
WlzInitGreyScan(WlzObject	*obj,	/* raster and tranpl of 0 provided */
		WlzIntervalWSpace	*iwsp,
		WlzGreyWSpace	*gwsp)
{
  return( WlzInitGreyRasterScan(obj, iwsp, gwsp, WLZ_RASTERDIR_ILIC, 0) );
}


/************************************************************************
*   Function   : WlzInitGreyRasterScan					*
*   Returns    :							*
*   Parameters :							*
*   Date       : Mon Oct 14 12:23:59 1996				*
*   Synopsis   :							*
* as initgreyscan, but with choice of raster direction and transplanting*
************************************************************************/

WlzErrorNum 
WlzInitGreyRasterScan(WlzObject 	*obj,
		      WlzIntervalWSpace *iwsp,
		      WlzGreyWSpace 	*gwsp,
		      WlzRasterDir 	raster,
		      int 		tranpl)
{
  WlzErrorNum	wlzerrno;

  if( (wlzerrno = WlzInitRasterScan(obj,iwsp,raster)) != WLZ_ERR_NONE ){
    return( wlzerrno );
  }
  return( WlzInitGreyWSpace(obj,iwsp,gwsp,tranpl) );
}

/************************************************************************
*   Function   : WlzInitGreyWSpace					*
*   Returns    :							*
*   Parameters :							*
*   Date       : Mon Oct 14 12:24:22 1996				*
*   Synopsis   :							*
* attach grey workspace to interval scanning workspace                  *
************************************************************************/

WlzErrorNum 
WlzInitGreyWSpace(WlzObject 		*obj,
		  WlzIntervalWSpace 	*iwsp,
		  WlzGreyWSpace 	*gwsp,
		  int 			tranpl)
{
  WlzValues 	vdmn;
  WlzErrorNum	wlzerrno = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    wlzerrno = WLZ_ERR_OBJECT_NULL;
  }
  else if((iwsp == NULL) || (gwsp == NULL))
  {
    wlzerrno = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    iwsp->gryptr = gwsp;
    /* Set up gwsp */
    gwsp->gtable = obj->values;
    vdmn = gwsp->gtable;
    gwsp->pixeltype = WlzGreyTableTypeToGreyType(vdmn.v->type, &wlzerrno);
  }
  if(wlzerrno == WLZ_ERR_NONE)
  {
    gwsp->gdomaintype = WlzGreyTableTypeToTableType(vdmn.v->type, &wlzerrno);
  }
  if(wlzerrno == WLZ_ERR_NONE)
  {
    switch (gwsp->gdomaintype) {

      case WLZ_GREY_TAB_RAGR:
	/* synchronise grey table line pointer to interval domain
	   - pointer when incremented points to first line */
	gwsp->gline = vdmn.v->vtblines + iwsp->linpos - vdmn.v->line1;
	break;

      case WLZ_GREY_TAB_RECT:
	/* rectangular grey table does not have valueline pointers */
	break;

      case WLZ_GREY_TAB_INTL:
	/* synchronise grey table line pointer to interval domain
	   - pointer when incremented points to first line */
	gwsp->gline = (WlzValueLine *) (vdmn.i->vil + iwsp->linpos
	    - vdmn.i->line1);
	break;

      default:
	wlzerrno = WLZ_ERR_VALUES_TYPE;
	break;

    }
    gwsp->intptr = iwsp;
    gwsp->tranpl = tranpl;
  }
  return( wlzerrno );
}



/************************************************************************
*   Function   : WlzNextGreyInterval					*
*   Returns    :							*
*   Parameters :							*
*   Date       : Mon Oct 14 12:26:36 1996				*
*   Synopsis   :							*
* obtain next interval and its grey table.
* if tranpl and gvio=0 and scan already started (check line number)
*	then must collect previous grey interval
* if tranpl and gvio=1 then must output this interval
*	(unless scan finished).
************************************************************************/

WlzErrorNum 
WlzNextGreyInterval(WlzIntervalWSpace *iwsp)
{
  WlzGreyWSpace *gwsp = iwsp->gryptr;
  WlzErrorNum	wlzerrno = WLZ_ERR_NONE;

  if(gwsp->tranpl && (gwsp->gvio == 0) &&
      (iwsp->linpos != (iwsp->linbot - iwsp->lineraster)))
  {
    wlzerrno = WlzGreyInterval(iwsp);
  }
  if(wlzerrno == WLZ_ERR_NONE)
  {
    wlzerrno = WlzNextInterval(iwsp);
  }
  if((wlzerrno == WLZ_ERR_NONE) && (!gwsp->tranpl || (gwsp->gvio == 1)))
  {
    wlzerrno = WlzGreyInterval(iwsp);
  }
  return(wlzerrno);
}


/************************************************************************
*   Function   : WlzGreyInterval					*
*   Returns    :							*
*   Parameters :							*
*   Date       : Mon Oct 14 12:28:05 1996				*
*   Synopsis   :							*
* handle grey table for an interval.
* This must follow a call to "nextinterval"
************************************************************************/

WlzErrorNum 
WlzGreyInterval(WlzIntervalWSpace *iwsp)
{
  WlzGreyP  g;
  WlzValues vdmn;
  WlzGreyWSpace *gwsp = iwsp->gryptr;
  int i;
  WlzValueIntervalLine *vil;
  WlzValueLine *val;

  vdmn = gwsp->gtable;
  switch( gwsp->gdomaintype ){

  case WLZ_GREY_TAB_RAGR:
    if (iwsp->nwlpos){
      /* new line signal, move to
	 next line of greytable */
      gwsp->gline += iwsp->nwlpos*iwsp->lineraster;
    }

    /* pointer to grey values */
    switch( gwsp->pixeltype ){

    case WLZ_GREY_INT:
      g.inp = (gwsp->gline->values.inp + iwsp->lftpos
	       - vdmn.v->kol1 - gwsp->gline->vkol1);
      break;

    case WLZ_GREY_SHORT:
      g.shp = (gwsp->gline->values.shp + iwsp->lftpos
	       - vdmn.v->kol1 - gwsp->gline->vkol1);
      break;

    case WLZ_GREY_UBYTE:
      g.ubp = (gwsp->gline->values.ubp + iwsp->lftpos
	       - vdmn.v->kol1 - gwsp->gline->vkol1);
      break;

    case WLZ_GREY_FLOAT:
      g.flp = (gwsp->gline->values.flp + iwsp->lftpos
	       - vdmn.v->kol1 - gwsp->gline->vkol1);
      break;

    case WLZ_GREY_DOUBLE:
      g.dbp = (gwsp->gline->values.dbp + iwsp->lftpos
	       - vdmn.v->kol1 - gwsp->gline->vkol1);
      break;

    default:
      return( WLZ_ERR_GREY_TYPE );
      /* break */

    }
    break;

  case WLZ_GREY_TAB_RECT:
    /* pointer to grey values */
    switch (gwsp->pixeltype) {

    case WLZ_GREY_INT:
      g.inp = (vdmn.r->values.inp + (iwsp->linpos - vdmn.r->line1)
	       * vdmn.r->width + iwsp->lftpos - vdmn.r->kol1);
      break;

    case WLZ_GREY_SHORT:
      g.shp = (vdmn.r->values.shp + (iwsp->linpos - vdmn.r->line1)
	       * vdmn.r->width + iwsp->lftpos - vdmn.r->kol1);
      break;

    case WLZ_GREY_UBYTE:
      g.ubp = (vdmn.r->values.ubp + (iwsp->linpos - vdmn.r->line1)
	       * vdmn.r->width + iwsp->lftpos - vdmn.r->kol1);
      break;

    case WLZ_GREY_FLOAT:
      g.flp = (vdmn.r->values.flp + (iwsp->linpos - vdmn.r->line1)
	       * vdmn.r->width + iwsp->lftpos - vdmn.r->kol1);
      break;

    case WLZ_GREY_DOUBLE:
      g.dbp = (vdmn.r->values.dbp + (iwsp->linpos - vdmn.r->line1)
	       * vdmn.r->width + iwsp->lftpos - vdmn.r->kol1);
      break;

    default:
      return( WLZ_ERR_GREY_TYPE );
      /* break */

    }
    break;

  case WLZ_GREY_TAB_INTL:
    if (iwsp->nwlpos){				
      /* new line signal, move to
	 next line of greytable */
      gwsp->gline =
	(WlzValueLine *) (((WlzValueIntervalLine *) gwsp->gline)
			  + iwsp->nwlpos*iwsp->lineraster);
    }

    vil = (WlzValueIntervalLine *) gwsp->gline;
    val = vil->vtbint;
    for (i=0; i<vil->nintvs; i++) {
      if( (vdmn.i->kol1 + val->vkol1 <= iwsp->lftpos)
	  && (vdmn.i->kol1 + val->vlastkl >= iwsp->rgtpos) ) {
	switch (gwsp->pixeltype) {

	case WLZ_GREY_INT:
	  g.inp = val->values.inp + iwsp->lftpos - vdmn.i->kol1 - val->vkol1;
	  break;

	case WLZ_GREY_SHORT:
	  g.shp = val->values.shp + iwsp->lftpos - vdmn.i->kol1 - val->vkol1;
	  break;

	case WLZ_GREY_UBYTE:
	  g.ubp = val->values.ubp + iwsp->lftpos - vdmn.i->kol1 - val->vkol1;
	  break;

	case WLZ_GREY_FLOAT:
	  g.flp = val->values.flp + iwsp->lftpos - vdmn.i->kol1 - val->vkol1;
	  break;

	case WLZ_GREY_DOUBLE:
	  g.dbp = val->values.dbp + iwsp->lftpos - vdmn.i->kol1 - val->vkol1;
	  break;

	default:
	  return( WLZ_ERR_GREY_TYPE );
	  /* break */

	}
	goto G_OK;
      }
      val++;
    }

    return( WLZ_ERR_VALUES_DATA );

  G_OK:
    break;

  default:
    return( WLZ_ERR_VALUES_TYPE );

  }


  if (gwsp->tranpl != 0) {	/* transplant mode */
    /* 
     * MOD 12/1/92 - in transplant mode, we will always copy to or from
     * an integer array.  This change is made so that routines like
     * seqpar can more easily deal with differing pixel types.
     */
    switch (gwsp->gvio) {

    case 0 :		/* input to object */
      switch(gwsp->pixeltype ) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp->colrmn; i++)
	  *g.inp++ = gwsp->u_grintptr.inp[i];
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp->colrmn; i++)
	  *g.shp++ = CLAMP (gwsp->u_grintptr.inp[i], SHRT_MIN, SHRT_MAX);
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp->colrmn; i++)
	  *g.ubp++ = CLAMP (gwsp->u_grintptr.inp[i], 0, 255);
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp->colrmn; i++)
	  *g.flp++ = gwsp->u_grintptr.inp[i];
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp->colrmn; i++)
	  *g.dbp++ = gwsp->u_grintptr.inp[i];
	break;

      default:
	/* can't get here because grey-type already checked */
	break;

      }
      break;

    case 1 :		/* output from object */
      switch(gwsp->pixeltype ) {

      case WLZ_GREY_INT:
	for (i=0; i<iwsp->colrmn; i++)
	  gwsp->u_grintptr.inp[i] = *g.inp++;
	break;

      case WLZ_GREY_SHORT:
	for (i=0; i<iwsp->colrmn; i++)
	  gwsp->u_grintptr.inp[i] = *g.shp++;
	break;

      case WLZ_GREY_UBYTE:
	for (i=0; i<iwsp->colrmn; i++)
	  gwsp->u_grintptr.inp[i] = *g.ubp++;
	break;

      case WLZ_GREY_FLOAT:
	for (i=0; i<iwsp->colrmn; i++)
	  gwsp->u_grintptr.inp[i] = *g.flp++;
	break;

      case WLZ_GREY_DOUBLE:
	for (i=0; i<iwsp->colrmn; i++)
	  gwsp->u_grintptr.inp[i] = *g.dbp++;
	break;

      default:
	/* can't get here because grey-type already checked */
	break;

      }
      break;

    default:
      return( WLZ_ERR_INT_DATA );

    }
  }
  else {
    gwsp->u_grintptr = g;	/* direct mode */
  }

  return( WLZ_ERR_NONE );
}
