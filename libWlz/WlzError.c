#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzError.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Error related functions for the Woolz library.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 28-01-00 bill	Add ALG_ERR_CONVERGENCE.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzErrorFromAlg
* Returns:	WlzErrorNum:		Woolz error code.
* Purpose:	Converts an Alg error code to a Woolz error code.
* Global refs:	-
* Parameters:	AlgError algErr:	Given Alg error code.
************************************************************************/
WlzErrorNum	WlzErrorFromAlg(AlgError algErr)
{
  WlzErrorNum	wlzErr = WLZ_ERR_ALG;

  switch(algErr)
  {
    case ALG_ERR_NONE:
      wlzErr = WLZ_ERR_NONE;
      break;
    case ALG_ERR_FUNC:
      wlzErr = WLZ_ERR_PARAM_DATA;
      break;
    case ALG_ERR_MALLOC:
      wlzErr = WLZ_ERR_MEM_ALLOC;
      break;
    case ALG_ERR_SINGULAR:
      wlzErr = WLZ_ERR_ALG_SINGULAR;
      break;
    case ALG_ERR_HOMOGENEOUS:
      wlzErr = WLZ_ERR_ALG_HOMOGENEOUS;
      break;
    case ALG_ERR_CONVERGENCE:
      wlzErr = WLZ_ERR_ALG_CONVERGENCE;
      break;
    case ALG_ERR_DIVZERO:
      wlzErr = WLZ_ERR_ALG;
      break;
  }
  return(wlzErr);
}
