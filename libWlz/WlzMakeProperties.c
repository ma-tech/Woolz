#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzMakeProperties.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Creates and free's a simple property list.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

/************************************************************************
*   Function   : WlzMakeCoreProperties					*
*   Returns    :WlzSimpleProperty *: pointer to property structure	*
*   Parameters :int size: size of data space required			*
*   Date       : Mon Oct 14 15:43:24 1996				*
*   Synopsis   :Allocate space for a WlzSimpleProperty which is a	*
*		structure with size bytes allocated for data		*
************************************************************************/

WlzSimpleProperty *
WlzMakeSimpleProperty(int size, WlzErrorNum *dstErr)
{
  WlzSimpleProperty	*p=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  if( (p = (WlzSimpleProperty *) AlcCalloc(1,size)) == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else {
    p->type = WLZ_PROPERTY_SIMPLE;
    p->linkcount = 0;
    p->size = size;

    if( (p->prop = AlcMalloc(size)) == NULL ){
      AlcFree( p );
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( p );
}

/************************************************************************
*   Function   : WlzFreeProperty					*
*   Date       : Mon Oct 21 19:58:10 1996				*
*************************************************************************
*   Synopsis   :Free space allocated for a WlzSimpleProperty		*
*   Returns    :WlzErrorNum: WLZ_ERR_NONE on success			*
*   Parameters :WlzSimpleProperty *prop: property structure to be freed	*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzFreeProperty(WlzSimpleProperty *prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (prop == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(prop->linkcount), &errNum) ){
    if( prop->prop != NULL ){
      AlcFree( prop->prop );
    }
    AlcFree( (void *) prop );
  }

  return errNum;
}

