#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzMakeCompound.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Makes Woolz compound objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

/************************************************************************
*   Function   : WlzMakeCompoundArray					*
*   Returns    :WlzCompoundArray *: compound array object		*
*   Parameters :WlzObjectType	type: compound object type:		*
*		WLZ_COMPOUND_ARR_1 or WLZ_COMPOUND_ARR_2		*
*		int 	mode: see below					*
*		int 	n: number of objects				*
*		WlzObject **ol: input object list			*
*		WlzObjectType otype: woolz object type if objects are	*
*			type checked.					*
*   Date       : Mon Oct 14 15:24:58 1996				*
*   Synopsis   :							*
* allocate a struct compounda.  Various mode-switched behaviour:
* mode==1:	Allocate empty array space for n objects.
* mode==2:	Array is input parameter ol; link.
* mode==3:	Array is input parameter ol; allocate space and copy,
*		incrementing object linkcounts.
************************************************************************/

WlzCompoundArray *
WlzMakeCompoundArray(WlzObjectType	type,
		     int 		mode,
		     int 		n,
		     WlzObject 		**ol,
		     WlzObjectType	otype,
		     WlzErrorNum	*dstErr)
{
  WlzCompoundArray	*co=NULL;
  int			i;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  
  if (type != WLZ_COMPOUND_ARR_1 && type != WLZ_COMPOUND_ARR_2) {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  /*
   * If appropriate, check objects have suitable type
   */
  if ((errNum == WLZ_ERR_NONE) && (type == WLZ_COMPOUND_ARR_1)) {
    switch(mode) {

    case 1:
      break;

    case 2:
    case 3:
      for (i=0; i<n; i++) {
	if (ol[i] != NULL && ol[i]->type != otype) {
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
	}
      }
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* Allocate space and fill if appropriate
     Note no default case because the mode has been checked above */
  if( errNum == WLZ_ERR_NONE ){
    switch (mode) {

    case 1:
      if( (co = (WlzCompoundArray *)
	   AlcCalloc(1,sizeof(WlzCompoundArray)
		     +n*sizeof(WlzObject *))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      co->o = (WlzObject **) (co+1);
      break;

    case 2:
      if( (co = (WlzCompoundArray *)
	   AlcCalloc(1,sizeof(WlzCompoundArray))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      co->o = ol;
      break;

    case 3:
      if( (co = (WlzCompoundArray *)
	   AlcCalloc(1,sizeof(WlzCompoundArray)
		     +n*sizeof(WlzObject *))) == NULL ){
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
      co->o = (WlzObject **) (co+1);
      for (i=0; i<n; i++) {
	co->o[i] = WlzAssignObject(ol[i], NULL);
      }
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    co->n = n;
    co->type = type;
    if (type == WLZ_COMPOUND_ARR_1){
      co->otype = otype;
    }
    co->linkcount = 0;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(co);
}
