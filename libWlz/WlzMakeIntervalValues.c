#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMakeIntervalValues.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Sep 26 11:29:23 2003
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
* \ingroup      WlzAllocation
* \brief        Makes an interval values table.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdlib.h>
#include <Wlz.h>

/* function:     WlzMakeIntervalValues    */
/*! 
* \ingroup      WlzAllocation
* \brief        Make in interval values table to match the input object
The table will have linkcount set to zero.
*
* \return       New interval values table.
* \param    type	Required table type.
* \param    obj	Input object.
* \param    bckgrnd	Values table background value.
* \param    dstErr	Error return.
* \par Constraints
 The woolz object type must resolve to WLZ_GREY_TAB_INTL using
 WlzGreyTableTypeToTableType()or an error will be reported and NULL
 returned.  For historical reasons
 the type encodes the table type as well as the grey-value type.
 
* \par      Source:
*                WlzMakeIntervalValues.c
*/
WlzIntervalValues *
WlzMakeIntervalValues(WlzObjectType	type,
		      WlzObject 	*obj,
		      WlzPixelV		bckgrnd,
		      WlzErrorNum	*dstErr)
{
  WlzValues 		v;
  WlzValueIntervalLine 	*vil;
  WlzValueLine 		*val;
  WlzGreyP 		g;
  WlzIntervalWSpace 	iwsp;
  WlzIntervalDomain 	*idom;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the values table type */
  v.i = NULL;
  (void) WlzGreyTableTypeToTableType(type, &errNum);
  if( errNum == WLZ_ERR_NONE ){
    (void) WlzGreyTableTypeToGreyType(type, &errNum);
  }

  /* check the object */
  if( (errNum == WLZ_ERR_NONE) && (obj == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      idom = obj->domain.i;
      break;

    case WLZ_EMPTY_OBJ:
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }

  /*
   * allocate space for basic, per line and per interval structures
   */
  if((errNum == WLZ_ERR_NONE) &&
     ((v.i = (WlzIntervalValues *) 
      AlcCalloc(sizeof(WlzIntervalValues)
		+ (idom->lastln - idom->line1 + 1)
		* sizeof(WlzValueIntervalLine) +
		WlzIntervalCount(idom, NULL) * sizeof(WlzValueLine), 1))
     == NULL) ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  if( errNum == WLZ_ERR_NONE ){
    vil = (WlzValueIntervalLine *) (v.i + 1);
    val = (WlzValueLine *) (vil + idom->lastln - idom->line1 + 1);
    v.i->bckgrnd = bckgrnd;
    switch( WlzGreyTableTypeToGreyType(type, NULL) ){

    case WLZ_GREY_INT:
      g.inp = (int *) AlcCalloc(WlzArea(obj, NULL), sizeof(int));
      break;

    case WLZ_GREY_SHORT:
      g.shp = (short *) AlcCalloc(WlzArea(obj, NULL), sizeof(short));
      break;

    case WLZ_GREY_UBYTE:
      g.ubp = (UBYTE *) AlcCalloc(WlzArea(obj, NULL), sizeof(UBYTE));
      break;

    case WLZ_GREY_FLOAT:
      g.flp = (float *) AlcCalloc(WlzArea(obj, NULL), sizeof(float));
      break;

    case WLZ_GREY_DOUBLE:
      g.dbp = (double *) AlcCalloc(WlzArea(obj, NULL), sizeof(double));
      break;

    case WLZ_GREY_RGBA:
      g.rgbp = (UINT *) AlcCalloc(WlzArea(obj, NULL), sizeof(UINT));
      break;

    default:
      WlzFreeValues( v );
      v.i = NULL;
      errNum = WLZ_ERR_GREY_TYPE;
      break;
    }
  }
  if( (errNum == WLZ_ERR_NONE) && (g.inp == NULL) ){
    WlzFreeValues( v );
    v.i = NULL;
    errNum = WLZ_ERR_MEM_ALLOC;
  }

  /*
   * fill in structure values and initialise scanning
   */
  if( errNum == WLZ_ERR_NONE ){
    v.i->type = type;
    v.i->freeptr = AlcFreeStackPush(v.i->freeptr, (void *)g.inp, NULL);
    v.i->line1 = idom->line1;
    v.i->lastln = idom->lastln;
    v.i->kol1 = idom->kol1;
    v.i->width = idom->lastkl - idom->kol1 + 1;
    v.i->vil = vil;
    v.i->linkcount = 0;
    v.i->original_table.core = NULL;
    vil--;
    if((errNum = WlzInitRasterScan(obj, &iwsp,
    				   WLZ_RASTERDIR_ILIC)) != WLZ_ERR_NONE ){
      WlzFreeValues( v );
      v.i = NULL;
    }
  }

  /*
   * fill in the line and interval structures,
     note grey-type already checked.
   */
  if( errNum == WLZ_ERR_NONE ){
    while( (errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
      if (iwsp.nwlpos != 0) {
	vil += iwsp.nwlpos;
	vil->vtbint = val;
      }
      vil->nintvs++;
      val->vkol1 = iwsp.lftpos - v.i->kol1;
      val->vlastkl = iwsp.rgtpos - v.i->kol1;
      switch( WlzGreyTableTypeToGreyType(type, NULL) ){

      case WLZ_GREY_INT:
	val->values.inp = g.inp;
	g.inp += iwsp.colrmn;
	break;

      case WLZ_GREY_SHORT:
	val->values.shp = g.shp;
	g.shp += iwsp.colrmn;
	break;

      case WLZ_GREY_UBYTE:
	val->values.ubp = g.ubp;
	g.ubp += iwsp.colrmn;
	break;

      case WLZ_GREY_FLOAT:
	val->values.flp = g.flp;
	g.flp += iwsp.colrmn;
	break;

      case WLZ_GREY_DOUBLE:
	val->values.dbp = g.dbp;
	g.dbp += iwsp.colrmn;
	break;

      case WLZ_GREY_RGBA:
	val->values.rgbp = g.rgbp;
	g.rgbp += iwsp.colrmn;
	break;

      }
      val++;
    }
    switch( errNum ){
    case WLZ_ERR_NONE:
    case WLZ_ERR_EOO:
      errNum = WLZ_ERR_NONE;
      break;
    default:
      WlzFreeValues( v );
      v.i = NULL;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(v.i);
}
 
