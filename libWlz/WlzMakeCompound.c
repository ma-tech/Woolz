#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMakeCompound.c
* \author       Richard Baldock
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Makes Woolz compound objects.
* \ingroup	WlzAllocation
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>

/*!
* \return	New Woolz compound array object.
* \ingroup	WlzAllocation
* \brief	Makes a new Woolz compound array object.
* \param	type			Type of compound array which should
*					be either WLZ_COMPOUND_ARR_1 or
*					WLZ_COMPOUND_ARR_2.
* \param	mode			Action to be performed, which may have
*					the values:
*					<ul>
*					<li>
*					1: Allocate an empty array with
*					space for n objects.
*					</li>
*					<li>
* 					2: Array members are given in the
*					input object list, objects are linked.
*					</li>
*					<li>
* 					3: Array members are given in the
*					input object list, objects copied.
*					</li>
*					</ul>
* \param	n			Number of objects.
* \param	ol			Input object list.
* \param	otype			Woolz object type if objects are
*					type checked.
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
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
