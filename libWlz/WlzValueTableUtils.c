#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzValueTableUtils.c
* Date:         March 1999
* Author:       Bill Hill, Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for computing value amd value table types.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzGreyTableType					*
* Returns:	WlzObjectType:		Type of grey table.		*
* Purpose:	Computes a grey table type from table and grey types.	*
* Global refs:	-							*
* Parameters:	WlzObjectType tableType: The basic table type.		*
*		WlzGreyType greyType:	The grey type.			*
*		WlzErrorNum *wlzErr:	Destination error pointer, may	*
*					BE NULL.			*
************************************************************************/
WlzObjectType	WlzGreyTableType(WlzObjectType tableType,
			         WlzGreyType greyType,
				 WlzErrorNum *wlzErr)
{
  WlzObjectType gTabType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(tableType)
  {
    case WLZ_GREY_TAB_RAGR:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_RAGR_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_RAGR_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_RAGR_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_RAGR_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_RAGR_DOUBLE;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_RECT:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_RECT_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_RECT_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_RECT_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_RECT_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_RECT_DOUBLE;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_INTL:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_INTL_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_INTL_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_INTL_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_INTL_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_INTL_DOUBLE;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    default:
      errNum = WLZ_ERR_VALUES_TYPE;
      break;
  }
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  return(gTabType);
}

/************************************************************************
* Function:	WlzGreyTableTypeToGreyType				*
* Returns:	WlzGreyType:		Type of grey value.		*
* Purpose:	Computes the type of grey from a grey table type.	*
* Global refs:	-							*
* Parameters:	WlzObjectType tableType: The basic table type.		*
*		WlzErrorNum *wlzErr:	Destination error pointer, may	*
*					BE NULL.			*
************************************************************************/
WlzGreyType	WlzGreyTableTypeToGreyType(WlzObjectType gTabType,
				           WlzErrorNum *wlzErr)
{
  WlzGreyType	greyType = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gTabType)
  {
    case WLZ_VALUETABLE_RAGR_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_INT:
      greyType = WLZ_GREY_INT;
      break;
    case WLZ_VALUETABLE_RAGR_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_SHORT:
      greyType = WLZ_GREY_SHORT;
      break;
    case WLZ_VALUETABLE_RAGR_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_UBYTE:
      greyType = WLZ_GREY_UBYTE;
      break;
    case WLZ_VALUETABLE_RAGR_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_FLOAT:
      greyType = WLZ_GREY_FLOAT;
      break;
    case WLZ_VALUETABLE_RAGR_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_DOUBLE:
      greyType = WLZ_GREY_DOUBLE;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  return(greyType);
}

/************************************************************************
* Function:	WlzGreyTableTypeToTableType				*
* Returns:	WlzGreyType:		Type of grey value.		*
* Purpose:	Computes the type of table from a  grey table type.	*
* Global refs:	-							*
* Parameters:	WlzObjectType tableType: The basic table type.		*
*		WlzErrorNum *wlzErr:	Destination error pointer, may	*
*					BE NULL.			*
************************************************************************/
WlzObjectType WlzGreyTableTypeToTableType(WlzObjectType gTabType,
				 	  WlzErrorNum *wlzErr)
{
  WlzObjectType	tableType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gTabType)
  {
    case WLZ_VALUETABLE_RAGR_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_DOUBLE:
      tableType = WLZ_GREY_TAB_RAGR;
      break;
    case WLZ_VALUETABLE_RECT_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_DOUBLE:
      tableType = WLZ_GREY_TAB_RECT;
      break;
    case WLZ_VALUETABLE_INTL_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_DOUBLE:
      tableType = WLZ_GREY_TAB_INTL;
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
  }
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  return(tableType);
}
